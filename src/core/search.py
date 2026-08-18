"""Search orchestrator.

Coordinates `SearchQuery` execution across one or more database clients on
worker threads, consults the disk cache, and pushes per-source results onto a
queue that the Tk main loop drains with `after()` polling.

Design notes:
    * One worker thread per (client, query) — independent network calls don't
      block each other.
    * Per-source results are returned as soon as ready (the UI fills in
      incrementally) — we don't wait for the slowest client.
    * Errors from one client never abort the others; they're packaged as
      SearchResult(status=ERROR) and reported in the status bar.
"""

from __future__ import annotations

import queue
import threading
import time
from typing import Callable, Optional

from .cache import DiskCache
from .models import SearchQuery, SearchResult, SearchStatus, ProteinVariant
from ..databases import AVAILABLE_CLIENTS, BAIT_SOURCES, DatabaseClient


class SearchOrchestrator:
    def __init__(
        self,
        cache: Optional[DiskCache] = None,
        email: str = "",
        api_key: str = "",
        timeout_s: int = 30,
    ) -> None:
        self.cache = cache
        self.email = email
        self.api_key = api_key
        self.timeout_s = timeout_s
        self._queue: "queue.Queue[SearchResult]" = queue.Queue()
        self._workers: list[threading.Thread] = []
        self._cancel = threading.Event()

    @property
    def queue(self) -> "queue.Queue[SearchResult]":
        return self._queue

    def cancel(self) -> None:
        """Signal cooperative cancellation — workers check between sources."""
        self._cancel.set()

    def reset(self) -> None:
        self._cancel.clear()
        self._workers = [t for t in self._workers if t.is_alive()]

    def start(self, query: SearchQuery, on_done: Optional[Callable[[], None]] = None) -> int:
        """Spawn one worker per enabled source. Returns the number of workers.

        Bait clients (BLAST, Foldseek) are sequenced *after* the sequence-DB
        workers finish, because they need bait sequences harvested from
        those results. The wait is bounded by the slowest sequence client.
        """
        self.reset()
        sources = [s for s in query.sources if s in AVAILABLE_CLIENTS]
        if not sources:
            return 0
        sequence_first = [s for s in sources if s not in BAIT_SOURCES]
        bait_sources = [s for s in sources if s in BAIT_SOURCES]

        total = len(sequence_first) + len(bait_sources)
        outstanding = {"n": total}
        lock = threading.Lock()
        seq_done = threading.Event()
        seq_remaining = {"n": len(sequence_first)}
        seq_lock = threading.Lock()
        if not sequence_first:
            seq_done.set()

        def finished_one() -> None:
            with lock:
                outstanding["n"] -= 1
                if outstanding["n"] == 0 and on_done:
                    on_done()

        def finished_seq() -> None:
            with seq_lock:
                seq_remaining["n"] -= 1
                if seq_remaining["n"] == 0:
                    seq_done.set()
            finished_one()

        for src in sequence_first:
            client_cls = AVAILABLE_CLIENTS[src]
            t = threading.Thread(
                target=self._run_one,
                args=(client_cls, query, finished_seq),
                daemon=True,
                name=f"search-{src}",
            )
            t.start()
            self._workers.append(t)

        for src in bait_sources:
            t = threading.Thread(
                target=self._run_bait,
                args=(src, query, seq_done, finished_one),
                daemon=True,
                name=f"search-{src}",
            )
            t.start()
            self._workers.append(t)
        return total

    # ---- internals -------------------------------------------------------

    def _run_bait(self, source: str, query: SearchQuery, seq_done: threading.Event, done: Callable[[], None]) -> None:
        """Bait clients (BLAST / Foldseek) need a query sequence. Wait for
        the sequence clients to finish, harvest the longest sequence per
        gene from their cached results, then run one job per gene.
        """
        start = time.time()
        try:
            cached = self._maybe_cached(source, query)
            if cached is not None:
                self._queue.put(SearchResult(
                    source, query, SearchStatus.OK if cached else SearchStatus.EMPTY,
                    variants=cached, message="(cached)", elapsed_s=time.time() - start,
                ))
                return
            seq_done.wait(timeout=600)
            if self._cancel.is_set():
                self._queue.put(SearchResult(source, query, SearchStatus.ERROR, message="cancelled"))
                return
            baits = self._gather_bait_sequences(query)
            if not baits:
                self._queue.put(SearchResult(
                    source, query, SearchStatus.EMPTY,
                    message="no bait sequences — re-search with 'Fetch sequences' checked",
                    elapsed_s=time.time() - start,
                ))
                return
            client = AVAILABLE_CLIENTS[source](
                email=self.email, api_key=self.api_key,
                timeout_s=self.timeout_s, bait_sequences=baits,
            )
            variants = client.search(query)
            self._maybe_store(source, query, variants)
            status = SearchStatus.OK if variants else SearchStatus.EMPTY
            self._queue.put(SearchResult(
                source=source, query=query, status=status, variants=variants,
                elapsed_s=time.time() - start,
            ))
        except Exception as e:
            self._queue.put(SearchResult(
                source=source, query=query, status=SearchStatus.ERROR,
                message=str(e), elapsed_s=time.time() - start,
            ))
        finally:
            done()

    def _gather_bait_sequences(self, query: SearchQuery) -> dict[str, str]:
        """Pull the longest sequence found per gene_symbol from the cache."""
        baits: dict[str, str] = {}
        if not self.cache:
            return baits
        for src in ("UniProt", "NCBI", "Ensembl"):
            variants = self.cache.get(query.cache_key(src)) or []
            for v in variants:
                if not v.sequence:
                    continue
                cur = baits.get(v.gene_symbol)
                if cur is None or len(v.sequence) > len(cur):
                    baits[v.gene_symbol] = v.sequence
        return baits

    def _run_one(self, client_cls: type[DatabaseClient], query: SearchQuery, done: Callable[[], None]) -> None:
        source = client_cls.name
        start = time.time()
        try:
            if self._cancel.is_set():
                self._queue.put(SearchResult(source, query, SearchStatus.ERROR, message="cancelled"))
                return

            variants = self._maybe_cached(source, query)
            cached = variants is not None
            if not cached:
                client = client_cls(email=self.email, api_key=self.api_key, timeout_s=self.timeout_s)
                variants = client.search(query)
                self._maybe_store(source, query, variants)

            status = SearchStatus.OK if variants else SearchStatus.EMPTY
            note = "(cached)" if cached else ""
            self._queue.put(SearchResult(
                source=source,
                query=query,
                status=status,
                variants=variants,
                message=note,
                elapsed_s=time.time() - start,
            ))
        except Exception as e:  # client raised — report but don't crash others
            self._queue.put(SearchResult(
                source=source,
                query=query,
                status=SearchStatus.ERROR,
                message=str(e),
                elapsed_s=time.time() - start,
            ))
        finally:
            done()

    def _maybe_cached(self, source: str, query: SearchQuery) -> Optional[list[ProteinVariant]]:
        if not self.cache:
            return None
        return self.cache.get(query.cache_key(source))

    def _maybe_store(self, source: str, query: SearchQuery, variants: list[ProteinVariant]) -> None:
        if not self.cache:
            return
        self.cache.put(query.cache_key(source), variants)
