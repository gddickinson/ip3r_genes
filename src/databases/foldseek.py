"""Foldseek web-service client — structure-based remote-homology search.

Foldseek (https://search.foldseek.com) lets you search a query sequence (or
structure) against AlphaFold DB / PDB / ESM Atlas using structural alignment,
which can find homologs at <20% sequence identity ("twilight zone") that
BLAST misses. The structure-first question this family poses — whether
the *fold* is conserved across distant lineages even when sequence is not —
is exactly the regime where Foldseek beats BLAST.

This client works differently from the sequence-DB clients because Foldseek
runs as an asynchronous *job*:

    1. POST sequence(s) to `/api/ticket/msa` → get a ticket id
    2. Poll `/api/ticket/{id}` until status == "COMPLETE"
    3. GET `/api/result/{id}/0` → JSON list of structural hits

The search is opt-in (slower than sequence searches; consumes a public-server
job slot). The client therefore requires `query.include_sequence=True` AND a
prior sequence-DB result to use as the bait — it picks the longest sequence
returned for each gene symbol and runs one Foldseek search per gene.

For a `SearchQuery` alone (no bait sequence), this client returns [] rather
than fail-loudly — the orchestrator surfaces it as EMPTY status so the user
sees "no Foldseek bait available; re-search with Fetch sequences first".
"""

from __future__ import annotations

import time
from typing import Any, Optional

import requests

from ..core.models import ProteinVariant, SearchQuery
from .base import DatabaseClient


FOLDSEEK_BASE = "https://search.foldseek.com/api"
# AlphaFold DB UniProt subset + PDB are the most useful default targets.
DEFAULT_DATABASES = ["afdb50", "pdb100"]


class FoldseekClient(DatabaseClient):
    name = "Foldseek"

    def __init__(
        self,
        email: str = "",
        api_key: str = "",
        timeout_s: int = 60,
        bait_sequences: Optional[dict[str, str]] = None,
        max_hits: int = 50,
        databases: Optional[list[str]] = None,
        max_wait_s: int = 180,
    ) -> None:
        super().__init__(email=email, api_key=api_key, timeout_s=timeout_s)
        self.bait_sequences = bait_sequences or {}
        self.max_hits = max_hits
        self.databases = databases or DEFAULT_DATABASES
        self.max_wait_s = max_wait_s

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        if not self.bait_sequences:
            return []
        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        for symbol in query.gene_symbols:
            bait = self.bait_sequences.get(symbol)
            if not bait:
                continue
            ticket = self._submit(bait)
            if not ticket:
                continue
            if not self._wait_for(ticket):
                continue
            hits = self._fetch_results(ticket)
            for h in hits[: self.max_hits]:
                v = self._hit_to_variant(h, symbol)
                if v and v.key() not in seen:
                    seen.add(v.key())
                    variants.append(v)
        return variants

    # ---- internals -------------------------------------------------------

    def _submit(self, sequence: str) -> Optional[str]:
        fasta = f">query\n{sequence}\n"
        try:
            r = requests.post(
                f"{FOLDSEEK_BASE}/ticket/msa",
                data={
                    "q": fasta,
                    "database[]": self.databases,
                    "mode": "3diaa",
                },
                timeout=self.timeout_s,
            )
            r.raise_for_status()
            data = r.json()
        except (requests.RequestException, ValueError):
            return None
        return data.get("id")

    def _wait_for(self, ticket_id: str) -> bool:
        deadline = time.time() + self.max_wait_s
        while time.time() < deadline:
            try:
                r = requests.get(f"{FOLDSEEK_BASE}/ticket/{ticket_id}", timeout=self.timeout_s)
                r.raise_for_status()
                status = r.json().get("status", "")
            except (requests.RequestException, ValueError):
                return False
            if status == "COMPLETE":
                return True
            if status == "ERROR":
                return False
            time.sleep(2)
        return False

    def _fetch_results(self, ticket_id: str) -> list[dict[str, Any]]:
        try:
            r = requests.get(f"{FOLDSEEK_BASE}/result/{ticket_id}/0", timeout=self.timeout_s)
            r.raise_for_status()
            data = r.json()
        except (requests.RequestException, ValueError):
            return []
        results = data.get("results", [])
        if not results:
            return []
        alignments = results[0].get("alignments", [])
        flat: list[dict[str, Any]] = []
        for db_block in alignments:
            if isinstance(db_block, list):
                flat.extend(db_block)
            elif isinstance(db_block, dict):
                flat.append(db_block)
        return flat

    def _hit_to_variant(self, hit: dict[str, Any], symbol: str) -> Optional[ProteinVariant]:
        target = hit.get("target") or hit.get("targetId") or ""
        if not target:
            return None
        prob = hit.get("prob") or hit.get("probability")
        seq_id = hit.get("seqId") or hit.get("pident")
        tm = hit.get("tmScore") or hit.get("tm_score")
        bits = []
        if prob is not None:
            bits.append(f"prob={prob}")
        if tm is not None:
            bits.append(f"TM={tm}")
        if seq_id is not None:
            bits.append(f"seqId={seq_id}")
        desc = "Foldseek hit · " + " ".join(bits) if bits else "Foldseek structural hit"
        # Foldseek often returns "<DB>_<accession>" — try to make a useful URL.
        url = ""
        if target.startswith("AF-"):
            acc = target.split("-")[1] if "-" in target else target
            url = f"https://alphafold.ebi.ac.uk/entry/{acc}"
        else:
            url = f"https://www.rcsb.org/structure/{target.split('_')[0]}"
        try:
            length = int(hit.get("tLen", 0) or 0) or None
        except (TypeError, ValueError):
            length = None
        return ProteinVariant(
            source=self.name,
            accession=target,
            gene_symbol=symbol,
            species=str(hit.get("taxName", "") or ""),
            length_aa=length,
            description=desc,
            gene_id=symbol,
            url=url,
            raw={k: str(v) for k, v in hit.items()},
        )
