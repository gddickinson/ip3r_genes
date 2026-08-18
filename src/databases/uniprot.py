"""UniProt REST client (rest.uniprot.org).

UniProt's search API takes a Lucene-style query string. We build:
    gene:<symbol> [ AND organism_id:<taxon> | AND organism_name:"<sp>" ]
and request the JSON view with only the fields we render.

UniProtKB returns one entry per *protein*; isoforms inside an entry appear in
`comments` of type "ALTERNATIVE PRODUCTS". We surface the canonical entry plus
each named isoform as separate ProteinVariants so length differences are
visible in the results table.
"""

from __future__ import annotations

from typing import Any

import requests

from ..core.models import ProteinVariant, SearchQuery
from .base import DatabaseClient


SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,gene_names,organism_name,organism_id,length,protein_name,sequence,cc_alternative_products"


class UniProtClient(DatabaseClient):
    name = "UniProt"

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        for symbol in query.gene_symbols:
            parts = [f"gene:{symbol}"]
            if query.taxon_id:
                parts.append(f"organism_id:{query.taxon_id}")
            elif query.species:
                parts.append(f'organism_name:"{query.species}"')
            params = {
                "query": " AND ".join(parts),
                "format": "json",
                "size": str(min(query.max_results_per_source, 500)),
                "fields": FIELDS,
            }
            try:
                r = requests.get(SEARCH_URL, params=params, timeout=self.timeout_s)
                r.raise_for_status()
                data = r.json()
            except requests.RequestException as e:
                raise RuntimeError(f"UniProt search failed: {e}") from e

            for entry in data.get("results", []):
                for v in self._entry_to_variants(entry, symbol, query.include_sequence):
                    if v.key() not in seen:
                        seen.add(v.key())
                        variants.append(v)
        return variants

    # ---- internals -------------------------------------------------------

    def _entry_to_variants(self, entry: dict[str, Any], symbol: str, include_seq: bool) -> list[ProteinVariant]:
        acc = entry.get("primaryAccession", "")
        if not acc:
            return []
        org = entry.get("organism", {}) or {}
        species = org.get("scientificName", "")
        taxon = org.get("taxonId")
        try:
            taxon = int(taxon) if taxon is not None else None
        except (TypeError, ValueError):
            taxon = None
        seq_info = entry.get("sequence", {}) or {}
        length = seq_info.get("length")
        seq_value = seq_info.get("value", "") if include_seq else ""
        prot = (entry.get("proteinDescription", {}) or {}).get("recommendedName", {}) or {}
        desc = (prot.get("fullName", {}) or {}).get("value", "") or entry.get("uniProtkbId", "")

        canonical = ProteinVariant(
            source=self.name,
            accession=acc,
            gene_symbol=symbol,
            species=species,
            taxon_id=taxon,
            length_aa=length,
            description=desc,
            gene_id=symbol,
            sequence=seq_value,
            url=f"https://www.uniprot.org/uniprotkb/{acc}",
            raw={"uniProtkbId": entry.get("uniProtkbId", "")},
        )
        out: list[ProteinVariant] = [canonical]

        for comment in entry.get("comments", []) or []:
            if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
                continue
            for iso in comment.get("isoforms", []) or []:
                ids = iso.get("isoformIds", []) or []
                if not ids:
                    continue
                iso_acc = ids[0]
                if iso_acc == acc:
                    continue
                out.append(ProteinVariant(
                    source=self.name,
                    accession=iso_acc,
                    gene_symbol=symbol,
                    species=species,
                    taxon_id=taxon,
                    length_aa=None,
                    description=f"isoform {iso.get('name', {}).get('value', '')} of {acc}",
                    gene_id=symbol,
                    url=f"https://www.uniprot.org/uniprotkb/{iso_acc}",
                    raw={"parent": acc, "isoform": iso},
                ))
        return out
