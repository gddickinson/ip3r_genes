"""Structural evidence for a candidate.

Reuses the existing AlphaFold and Foldseek clients but packages them into a
single `StructuralEvidence` object investigation can consume:

    * AlphaFold model URL + global pLDDT + (when available) per-residue pLDDT
      summary stats over each TM region
    * Foldseek hit list against the canonical family panel (opt-in — Foldseek
      is async and slow; default off, --foldseek to enable)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests

from ..databases.foldseek import FoldseekClient


ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"


@dataclass
class StructuralEvidence:
    accession: str
    alphafold_available: bool = False
    alphafold_url: str = ""
    cif_url: str = ""
    pae_doc_url: str = ""
    pae_image_url: str = ""
    global_plddt: Optional[float] = None
    plddt_in_tm: Optional[float] = None     # mean pLDDT averaged over TM regions
    foldseek_hits: list[dict] = field(default_factory=list)


def fetch_alphafold(accession: str, timeout_s: int = 30) -> Optional[dict]:
    if not accession:
        return None
    try:
        r = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=timeout_s)
    except requests.RequestException:
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


def fetch_per_residue_plddt(accession: str, timeout_s: int = 30) -> list[float]:
    """Pull the per-residue pLDDT vector from AlphaFold's PDB file.

    pLDDT is encoded in the B-factor column of the PDB ATOM records (one
    value per residue, repeated for each atom). We parse just CA atoms.
    """
    pdb_url = f"{ALPHAFOLD_API}/{accession}"
    pred = fetch_alphafold(accession, timeout_s)
    if not pred:
        return []
    url = pred.get("pdbUrl", "")
    if not url:
        return []
    try:
        r = requests.get(url, timeout=timeout_s)
        if r.status_code != 200:
            return []
    except requests.RequestException:
        return []
    plddt: list[float] = []
    seen_residues: set[int] = set()
    for line in r.text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            res_seq = int(line[22:26])
            b_factor = float(line[60:66])
        except (ValueError, IndexError):
            continue
        if res_seq in seen_residues:
            continue
        seen_residues.add(res_seq)
        plddt.append(b_factor)
    return plddt


def gather_structural_evidence(
    accession: str,
    transmembrane_regions: Optional[list[tuple[int, int]]] = None,
    foldseek_bait_sequence: str = "",
    run_foldseek: bool = False,
    timeout_s: int = 30,
) -> StructuralEvidence:
    ev = StructuralEvidence(accession=accession)
    pred = fetch_alphafold(accession, timeout_s)
    if pred:
        ev.alphafold_available = True
        ev.alphafold_url = pred.get("pdbUrl", "")
        ev.cif_url = pred.get("cifUrl", "")
        ev.pae_doc_url = pred.get("paeDocUrl", "")
        ev.pae_image_url = pred.get("paeImageUrl", "")
        gv = pred.get("globalMetricValue")
        if isinstance(gv, (int, float)):
            ev.global_plddt = float(gv)

        if transmembrane_regions:
            per_res = fetch_per_residue_plddt(accession, timeout_s)
            if per_res:
                total, n = 0.0, 0
                for s, e in transmembrane_regions:
                    s = max(1, s); e = min(len(per_res), e)
                    if s <= e:
                        total += sum(per_res[s - 1:e])
                        n += (e - s + 1)
                if n > 0:
                    ev.plddt_in_tm = total / n

    if run_foldseek and foldseek_bait_sequence:
        client = FoldseekClient(
            bait_sequences={"query": foldseek_bait_sequence},
            max_hits=20,
            max_wait_s=240,
            timeout_s=timeout_s,
        )
        from ..core.models import SearchQuery
        hits = client.search(SearchQuery(gene_symbols=["query"]))
        for h in hits:
            ev.foldseek_hits.append({
                "target": h.accession,
                "description": h.description,
                "url": h.url,
                **{k: v for k, v in h.raw.items()},
            })

    return ev
