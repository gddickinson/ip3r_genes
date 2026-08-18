"""Full domain architecture for a candidate via InterPro.

The discovery scorer only checked for two specific Pfam IDs. Investigation
asks the wider question: *what are ALL the domains in this protein, and is
the family-of-interest domain the dominant one?*

The signal we want:
    * coverage of the family-defining IPR/Pfam over the whole protein
    * how many TM helices are predicted (IP3Rs: 6 per subunit)
    * what other domain families are co-present (FN3, PKD, etc.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests


INTERPRO_ALL = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot"


@dataclass
class DomainAnnotation:
    source: str          # 'pfam' | 'interpro' | 'panther' | 'cdd' | ...
    accession: str       # e.g. PF12166 / IPR027272
    name: str
    type: str            # 'family' | 'domain' | 'homologous_superfamily' | ...
    locations: list[tuple[int, int]] = field(default_factory=list)

    def coverage(self, protein_length: int) -> float:
        if not self.locations or protein_length <= 0:
            return 0.0
        covered = sum(max(0, e - s + 1) for s, e in self.locations)
        return min(1.0, covered / protein_length)


@dataclass
class DomainArchitecture:
    accession: str
    protein_length: int = 0
    annotations: list[DomainAnnotation] = field(default_factory=list)
    n_transmembrane: int = 0       # filled in by uniprot_detail; cached here
    family_signatures: list[str] = field(default_factory=list)   # e.g. ["IPR027272", "PF12166", ...]

    def find_family_coverage(self, family_keywords: list[str]) -> float:
        """Best coverage of any annotation whose name contains a keyword."""
        kw = [k.lower() for k in family_keywords]
        best = 0.0
        for a in self.annotations:
            if any(k in a.name.lower() for k in kw):
                best = max(best, a.coverage(self.protein_length))
        return best

    def has_family_signature(self, family_keywords: list[str]) -> bool:
        kw = [k.lower() for k in family_keywords]
        return any(any(k in a.name.lower() for k in kw) for a in self.annotations)


def fetch_domain_architecture(accession: str, timeout_s: int = 30) -> Optional[DomainArchitecture]:
    if not accession:
        return None
    try:
        r = requests.get(f"{INTERPRO_ALL}/{accession}/", timeout=timeout_s)
    except requests.RequestException:
        return None
    if r.status_code != 200:
        return None
    try:
        data = r.json()
    except ValueError:
        return None

    arch = DomainArchitecture(accession=accession)
    protein_len = 0
    for entry in data.get("results", []):
        meta = entry.get("metadata", {}) or {}
        proteins = entry.get("proteins", []) or []
        locs: list[tuple[int, int]] = []
        for p in proteins:
            try:
                protein_len = max(protein_len, int(p.get("protein_length", 0) or 0))
            except (TypeError, ValueError):
                pass
            for loc in p.get("entry_protein_locations", []) or []:
                for frag in loc.get("fragments", []) or []:
                    try:
                        locs.append((int(frag["start"]), int(frag["end"])))
                    except (KeyError, TypeError, ValueError):
                        continue
        ann = DomainAnnotation(
            source=str(meta.get("source_database") or ""),
            accession=str(meta.get("accession") or ""),
            name=str(meta.get("name") or ""),
            type=str(meta.get("type") or ""),
            locations=locs,
        )
        arch.annotations.append(ann)
    arch.protein_length = protein_len
    arch.annotations.sort(key=lambda a: (-len(a.locations), a.source))
    return arch
