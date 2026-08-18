"""Ensembl REST client (rest.ensembl.org).

Two query paths depending on whether a species was given:

  * With a species → `/lookup/symbol/{species}/{gene}` resolves the gene id,
    then `/lookup/id/{gene_id}?expand=1` enumerates transcripts. Each
    transcript that has a `Translation` becomes one ProteinVariant.
    `/xrefs/symbol/` is kept as a fallback — see `_symbol_to_ids`.

  * Without a species → fall back to the cross-species symbol-based search
    by trying a curated list of common reference species (see species.py),
    so that the user gets ortholog-style results without having to know
    Ensembl's per-species URL.
"""

from __future__ import annotations

import time
from typing import Any

import requests

from ..core.models import ProteinVariant, SearchQuery
from ..utils.species import DEFAULT_SPECIES_PANEL, ensembl_species_slug
from .base import DatabaseClient


BASE = "https://rest.ensembl.org"
HEADERS = {"Accept": "application/json"}


class EnsemblClient(DatabaseClient):
    name = "Ensembl"

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        species_list = [query.species] if query.species else DEFAULT_SPECIES_PANEL
        for symbol in query.gene_symbols:
            for species in species_list:
                slug = ensembl_species_slug(species)
                try:
                    gene_ids = self._symbol_to_ids(slug, symbol)
                except requests.RequestException as e:
                    if query.species:
                        raise RuntimeError(f"Ensembl symbol lookup failed: {e}") from e
                    continue  # silent skip when scanning the panel
                for gid in gene_ids[: query.max_results_per_source]:
                    try:
                        gene = self._lookup_gene(gid)
                    except requests.RequestException:
                        continue
                    for v in self._gene_to_variants(gene, symbol, query.include_sequence):
                        if v.key() not in seen:
                            seen.add(v.key())
                            variants.append(v)
                time.sleep(0.1)
        return variants

    # ---- internals -------------------------------------------------------

    def _get(self, path: str, attempts: int = 3) -> Any:
        """GET with retries.

        `rest.ensembl.org` is intermittently slow rather than down — a symbol
        lookup that normally answers in under a second was measured at 55 s on
        2026-08-18, and it also returns transient 500/503s under load. One slow
        response should not cost the session an entire source, so read timeouts
        and 5xx are retried with a backoff before the error is reported.
        """
        last: Exception | None = None
        for attempt in range(attempts):
            try:
                r = requests.get(BASE + path, headers=HEADERS,
                                 timeout=self.timeout_s)
            except requests.Timeout as exc:
                last = exc
            else:
                if r.status_code == 404:
                    return None
                if r.status_code < 500:
                    r.raise_for_status()
                    return r.json()
                last = requests.HTTPError(f"{r.status_code} from Ensembl",
                                          response=r)
            time.sleep(1.5 * (attempt + 1))
        raise last if last else RuntimeError("Ensembl request failed")

    def _symbol_to_ids(self, species_slug: str, symbol: str) -> list[str]:
        """Gene id(s) for a symbol, via `lookup/symbol` with an `xrefs` fallback.

        S0 (2026-08-18) measured a per-species fault at Ensembl:
        `/xrefs/symbol/homo_sapiens/{symbol}` stalls indefinitely — no
        response and no error, for `BRCA2` as well as `ITPR1` — while the same
        endpoint answers in 0.6 s for `danio_rerio` and
        `/lookup/symbol/homo_sapiens/{symbol}` answers normally. A retry budget
        cannot fix that, because there is nothing to retry against; the human
        ITPR preset lost Ensembl entirely
        (`results/s0_baseline/ensembl_endpoint_probe.tsv`).

        So `lookup/symbol` is now the primary path. It returns the one gene
        Ensembl considers canonical for the symbol, which is what this client
        wants anyway. `xrefs/symbol` — which can return several ids — is kept
        as a fallback for the cases where `lookup/symbol` finds nothing, so no
        recall is lost when the outage clears.
        """
        gene = self._get(f"/lookup/symbol/{species_slug}/{symbol}")
        if gene and gene.get("id"):
            return [gene["id"]]

        # Only a clean 404 (`_get` → None) reaches the fallback. A transport
        # error is deliberately allowed to propagate instead: retrying it
        # against `xrefs/symbol` would walk straight into the stalling
        # endpoint and burn `attempts × timeout_s` per gene per species for
        # nothing.
        data = self._get(f"/xrefs/symbol/{species_slug}/{symbol}?object_type=gene")
        if not data:
            return []
        return [entry["id"] for entry in data if entry.get("type") == "gene"]

    def _lookup_gene(self, gene_id: str) -> dict[str, Any]:
        data = self._get(f"/lookup/id/{gene_id}?expand=1")
        return data or {}

    def _gene_to_variants(self, gene: dict[str, Any], symbol: str, include_seq: bool) -> list[ProteinVariant]:
        out: list[ProteinVariant] = []
        if not gene:
            return out
        species = (gene.get("species") or "").replace("_", " ").capitalize()
        gene_id = gene.get("id", "")
        for tx in gene.get("Transcript", []) or []:
            trans = tx.get("Translation")
            if not trans:
                continue
            acc = trans.get("id", "")
            if not acc:
                continue
            length = trans.get("length")
            try:
                length = int(length) if length is not None else None
            except (TypeError, ValueError):
                length = None
            v = ProteinVariant(
                source=self.name,
                accession=acc,
                gene_symbol=symbol,
                species=species,
                length_aa=length,
                description=(tx.get("display_name") or tx.get("biotype") or ""),
                transcript_id=tx.get("id", ""),
                gene_id=gene_id,
                url=f"https://www.ensembl.org/{gene.get('species','')}/Gene/Summary?g={gene_id}",
                raw={"gene": gene_id, "transcript": tx.get("id"), "biotype": tx.get("biotype", "")},
            )
            if include_seq:
                v.sequence = self._fetch_sequence(acc)
            out.append(v)
        return out

    def _fetch_sequence(self, protein_id: str) -> str:
        try:
            r = requests.get(
                f"{BASE}/sequence/id/{protein_id}?type=protein",
                headers={"Accept": "text/plain"},
                timeout=self.timeout_s,
            )
            if r.status_code == 200:
                return r.text.strip()
        except requests.RequestException:
            pass
        return ""
