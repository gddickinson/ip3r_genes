"""AlphaFold Protein Structure Database client.

The AlphaFold DB (https://alphafold.ebi.ac.uk) hosts predicted 3D structures
for >200M proteins, keyed by UniProt accession. There is no "search by gene
symbol" endpoint, so this client takes a different shape from the sequence
clients: it pivots off UniProt — first resolves `gene` → UniProt accessions
via the UniProt search API, then fetches the AlphaFold prediction metadata
for each.

What the user gets:
    * One ProteinVariant per UniProt accession that has an AlphaFold model.
    * `description` field includes mean pLDDT (folding confidence) and the
      direct model URL. The user can paste that URL into PyMOL / ChimeraX or
      open it in a browser to see the 3D structure.
    * `raw` payload preserves the full prediction record (model URLs in PDB,
      mmCIF, BCIF; PAE JSON URL; coverage info) for the detail pane.

Why this is a "cutting-edge" complement to BLAST: AlphaFold DB lets users
find proteins whose *fold* is similar to the query, mediated through
structure-aware downstream tools (Foldseek — see foldseek.py). When the
sequence is in the twilight zone (<20% identity), structure is conserved
where sequence is not: a divergent candidate that still folds into the
IP3R architecture — the beta-trefoil IP3-binding core over a six-TM pore,
assembled as a fourfold-symmetric tetramer — is a family member however
modest its identity to ITPR1/2/3.
"""

from __future__ import annotations

import time
from typing import Any

import requests

from ..core.models import ProteinVariant, SearchQuery
from .base import DatabaseClient
from .uniprot import SEARCH_URL as UNIPROT_SEARCH_URL


ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"


class AlphaFoldClient(DatabaseClient):
    name = "AlphaFold"

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        for symbol in query.gene_symbols:
            accessions = self._uniprot_accessions(symbol, query)
            for acc in accessions[: query.max_results_per_source]:
                pred = self._fetch_prediction(acc)
                if not pred:
                    continue
                v = self._prediction_to_variant(pred, symbol, acc)
                if v and v.key() not in seen:
                    seen.add(v.key())
                    variants.append(v)
                time.sleep(0.05)  # AlphaFold API is generous but be polite
        return variants

    # ---- internals -------------------------------------------------------

    def _uniprot_accessions(self, symbol: str, query: SearchQuery) -> list[str]:
        parts = [f"gene:{symbol}"]
        if query.taxon_id:
            parts.append(f"organism_id:{query.taxon_id}")
        elif query.species:
            parts.append(f'organism_name:"{query.species}"')
        params = {
            "query": " AND ".join(parts),
            "format": "json",
            "size": str(min(query.max_results_per_source, 100)),
            "fields": "accession",
        }
        try:
            r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=self.timeout_s)
            r.raise_for_status()
            data = r.json()
        except requests.RequestException as e:
            raise RuntimeError(f"AlphaFold (via UniProt) lookup failed: {e}") from e
        return [e.get("primaryAccession", "") for e in data.get("results", []) if e.get("primaryAccession")]

    def _fetch_prediction(self, accession: str) -> dict[str, Any] | None:
        try:
            r = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=self.timeout_s)
        except requests.RequestException:
            return None
        if r.status_code == 404:
            return None
        if r.status_code != 200:
            return None
        try:
            data = r.json()
        except ValueError:
            return None
        if isinstance(data, list):
            return data[0] if data else None
        return data

    def _prediction_to_variant(self, pred: dict[str, Any], symbol: str, acc: str) -> ProteinVariant | None:
        model_url = pred.get("pdbUrl") or pred.get("cifUrl") or ""
        plddt = pred.get("globalMetricValue")  # mean pLDDT 0-100
        plddt_txt = f"pLDDT={plddt:.1f}" if isinstance(plddt, (int, float)) else ""
        organism = pred.get("organismScientificName", "") or pred.get("organism", "")
        taxon = pred.get("taxId")
        try:
            taxon = int(taxon) if taxon is not None else None
        except (TypeError, ValueError):
            taxon = None
        try:
            length = int(pred.get("uniprotEnd", 0) or 0)
        except (TypeError, ValueError):
            length = 0
        desc_bits = [pred.get("uniprotDescription", "") or "", plddt_txt]
        desc = " · ".join(b for b in desc_bits if b)
        return ProteinVariant(
            source=self.name,
            accession=acc,
            gene_symbol=symbol,
            species=organism,
            taxon_id=taxon,
            length_aa=length or None,
            description=desc or "AlphaFold predicted structure",
            gene_id=symbol,
            url=model_url or f"https://alphafold.ebi.ac.uk/entry/{acc}",
            raw={
                "modelEntityId": pred.get("modelEntityId", ""),
                "pdbUrl": pred.get("pdbUrl", ""),
                "cifUrl": pred.get("cifUrl", ""),
                "bcifUrl": pred.get("bcifUrl", ""),
                "paeImageUrl": pred.get("paeImageUrl", ""),
                "paeDocUrl": pred.get("paeDocUrl", ""),
                "globalMetricValue": plddt,
            },
        )
