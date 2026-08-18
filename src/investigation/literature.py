"""PubMed search for a candidate accession + organism.

Two complementary queries:
    1. PubMed search "{accession} OR {uniprot_id}" — direct accession hits
    2. PubMed search "{organism} AND <family keyword>" — relevant
       publications on the organism (broader; catches papers describing an
       IP3-receptor-like channel without citing the specific accession)

Both use Biopython's Bio.Entrez (already a dependency).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

from Bio import Entrez

from ..utils.family import LITERATURE_KEYWORD


@dataclass
class LiteratureEvidence:
    accession: str
    direct_pmids: list[str] = field(default_factory=list)
    organism_family_pmids: list[str] = field(default_factory=list)
    direct_titles: list[tuple[str, str]] = field(default_factory=list)         # (pmid, title)
    organism_family_titles: list[tuple[str, str]] = field(default_factory=list)
    notes: str = ""


def search_literature(
    accession: str,
    uniprot_id: str = "",
    organism: str = "",
    email: str = "",
    max_titles: int = 5,
    timeout_s: int = 30,
) -> LiteratureEvidence:
    if email:
        Entrez.email = email
    ev = LiteratureEvidence(accession=accession)

    # Direct accession search
    terms = [accession] + ([uniprot_id] if uniprot_id else [])
    direct_q = " OR ".join(f'"{t}"' for t in terms)
    try:
        with Entrez.esearch(db="pubmed", term=direct_q, retmax=20) as h:
            data = Entrez.read(h)
        ev.direct_pmids = list(data.get("IdList", []))
    except Exception as e:
        ev.notes += f"direct search failed: {e}; "

    time.sleep(0.34)
    if organism:
        org_q = (f'"{organism}"[Title/Abstract] AND '
                 f'"{LITERATURE_KEYWORD}"[Title/Abstract]')
        try:
            with Entrez.esearch(db="pubmed", term=org_q, retmax=20) as h:
                data = Entrez.read(h)
            ev.organism_family_pmids = list(data.get("IdList", []))
        except Exception as e:
            ev.notes += f"organism search failed: {e}; "

    # Pull titles for the top results
    def _titles(pmids: list[str]) -> list[tuple[str, str]]:
        if not pmids:
            return []
        try:
            with Entrez.esummary(db="pubmed", id=",".join(pmids[:max_titles])) as h:
                rec = Entrez.read(h)
        except Exception:
            return []
        out = []
        for r in rec:
            out.append((str(r.get("Id", "")), str(r.get("Title", ""))))
        return out

    ev.direct_titles = _titles(ev.direct_pmids)
    time.sleep(0.34)
    ev.organism_family_titles = _titles(ev.organism_family_pmids)
    return ev
