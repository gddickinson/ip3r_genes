"""NCBI Entrez client via Biopython.

Query strategy:
    1. esearch in the `protein` database for `<gene>[Gene Name] AND <organism>[Organism]`
       (organism clause dropped if no species was given).
    2. esummary on the resulting protein UIDs to get accession, length, title,
       organism, and the linked gene id.
    3. Optionally efetch FASTA for sequence content if include_sequence is set.

We deliberately use the `protein` index rather than going gene → protein via
elink, because the protein index already enumerates every isoform / variant
record, which is what the user actually wants.
"""

from __future__ import annotations

import time
from typing import Any

from Bio import Entrez

from ..core.models import ProteinVariant, SearchQuery
from .base import DatabaseClient


class NCBIClient(DatabaseClient):
    name = "NCBI"

    def __init__(self, email: str = "", api_key: str = "", timeout_s: int = 30) -> None:
        super().__init__(email=email, api_key=api_key, timeout_s=timeout_s)
        Entrez.email = email or "anonymous@example.com"
        if api_key:
            Entrez.api_key = api_key

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        for symbol in query.gene_symbols:
            term = f"{symbol}[Gene Name]"
            if query.species:
                term += f" AND \"{query.species}\"[Organism]"
            elif query.taxon_id:
                term += f" AND txid{query.taxon_id}[Organism:exp]"
            try:
                uids = self._esearch(term, query.max_results_per_source)
            except Exception as e:
                raise RuntimeError(f"NCBI esearch failed: {e}") from e
            if not uids:
                continue
            for summary in self._esummary_chunks(uids):
                variant = self._summary_to_variant(summary, symbol)
                if variant and variant.key() not in seen:
                    seen.add(variant.key())
                    variants.append(variant)
            if query.include_sequence:
                self._fill_sequences(variants)
            time.sleep(0.34)  # be polite — NCBI allows 3 req/s without a key
        return variants

    # ---- internals -------------------------------------------------------

    def _esearch(self, term: str, retmax: int) -> list[str]:
        with Entrez.esearch(db="protein", term=term, retmax=retmax) as h:
            data = Entrez.read(h)
        return list(data.get("IdList", []))

    def _esummary_chunks(self, uids: list[str], chunk: int = 100) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        for i in range(0, len(uids), chunk):
            batch = uids[i : i + chunk]
            with Entrez.esummary(db="protein", id=",".join(batch)) as h:
                summaries = Entrez.read(h)
            out.extend(summaries)
            time.sleep(0.34)
        return out

    def _summary_to_variant(self, s: dict[str, Any], symbol: str) -> ProteinVariant | None:
        accession = str(s.get("AccessionVersion") or s.get("Caption") or "")
        if not accession:
            return None
        try:
            length = int(s.get("Length", 0)) or None
        except (TypeError, ValueError):
            length = None
        taxid = s.get("TaxId")
        try:
            taxid = int(taxid) if taxid is not None else None
        except (TypeError, ValueError):
            taxid = None
        return ProteinVariant(
            source=self.name,
            accession=accession,
            gene_symbol=symbol,
            species=str(s.get("Organism", "") or ""),
            taxon_id=taxid,
            length_aa=length,
            description=str(s.get("Title", "") or ""),
            url=f"https://www.ncbi.nlm.nih.gov/protein/{accession}",
            raw={k: str(v) for k, v in s.items()},
        )

    def _fill_sequences(self, variants: list[ProteinVariant]) -> None:
        ids = [v.accession for v in variants if not v.sequence]
        if not ids:
            return
        try:
            with Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text") as h:
                fasta = h.read()
        except Exception:
            return
        by_acc = {v.accession: v for v in variants}
        cur_acc, cur_seq_parts = None, []
        for line in fasta.splitlines():
            if line.startswith(">"):
                if cur_acc and cur_acc in by_acc:
                    by_acc[cur_acc].sequence = "".join(cur_seq_parts)
                header = line[1:].split()[0]
                cur_acc = header
                cur_seq_parts = []
            else:
                cur_seq_parts.append(line.strip())
        if cur_acc and cur_acc in by_acc:
            by_acc[cur_acc].sequence = "".join(cur_seq_parts)
