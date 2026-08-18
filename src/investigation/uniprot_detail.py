"""Detailed UniProt entry — cross-refs, features, GO terms, references.

The standard search client returns just enough to build a ProteinVariant.
Investigation pulls the full UniProtKB JSON to surface:
    * transmembrane / topology features (count, positions)
    * GO terms (especially channel / ion-transport related)
    * existing PDB / Pfam / PANTHER / OMA / OrthoDB cross-refs
    * literature references (publication ids → see literature.py for PubMed)
    * protein evidence level (experimental vs predicted)
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

import requests


UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb"


@dataclass
class UniProtDetail:
    accession: str
    uniprot_id: str = ""
    reviewed: bool = False
    sequence_length: int = 0
    organism: str = ""
    taxon_id: Optional[int] = None
    protein_name: str = ""
    gene_names: list[str] = field(default_factory=list)
    n_transmembrane: int = 0
    transmembrane_regions: list[tuple[int, int]] = field(default_factory=list)
    feature_counts: dict[str, int] = field(default_factory=dict)
    go_terms: list[tuple[str, str, str]] = field(default_factory=list)   # (id, name, aspect)
    xref_counts: dict[str, int] = field(default_factory=dict)
    pdb_ids: list[str] = field(default_factory=list)
    pubmed_ids: list[str] = field(default_factory=list)
    evidence_level: str = ""
    raw_keywords: list[str] = field(default_factory=list)

    def channel_related_go_terms(self) -> list[tuple[str, str, str]]:
        keys = ("channel", "ion transport", "mechanosensitive", "calcium ion", "cation")
        return [g for g in self.go_terms if any(k in g[1].lower() for k in keys)]


def fetch_uniprot_detail(accession: str, timeout_s: int = 30) -> Optional[UniProtDetail]:
    if not accession:
        return None
    try:
        r = requests.get(f"{UNIPROT_JSON}/{accession}.json", timeout=timeout_s)
    except requests.RequestException:
        return None
    if r.status_code != 200:
        return None
    try:
        data = r.json()
    except ValueError:
        return None

    organism = data.get("organism", {}) or {}
    seq = data.get("sequence", {}) or {}
    desc = (data.get("proteinDescription", {}) or {}).get("recommendedName", {}) or {}
    prot_name = (desc.get("fullName", {}) or {}).get("value", "") or data.get("uniProtkbId", "")

    detail = UniProtDetail(
        accession=data.get("primaryAccession", accession),
        uniprot_id=data.get("uniProtkbId", ""),
        reviewed=(data.get("entryType", "").lower().startswith("uniprotkb reviewed")
                  or "swiss-prot" in data.get("entryType", "").lower()),
        sequence_length=int(seq.get("length", 0) or 0),
        organism=str(organism.get("scientificName", "") or ""),
        taxon_id=(int(organism.get("taxonId")) if organism.get("taxonId") is not None else None),
        protein_name=str(prot_name),
        evidence_level=str(data.get("proteinExistence", "") or ""),
    )

    for g in data.get("genes", []) or []:
        gn = (g.get("geneName") or {}).get("value")
        if gn:
            detail.gene_names.append(str(gn))
        for n in g.get("orderedLocusNames", []) or []:
            v = n.get("value") if isinstance(n, dict) else n
            if v:
                detail.gene_names.append(str(v))

    feats = data.get("features", []) or []
    detail.feature_counts = dict(Counter(f.get("type", "?") for f in feats))
    for f in feats:
        if f.get("type") == "Transmembrane":
            loc = f.get("location", {}) or {}
            try:
                start = int((loc.get("start") or {}).get("value", 0))
                end = int((loc.get("end") or {}).get("value", 0))
            except (TypeError, ValueError):
                continue
            if start and end:
                detail.transmembrane_regions.append((start, end))
    detail.n_transmembrane = len(detail.transmembrane_regions)

    xrefs = data.get("uniProtKBCrossReferences", []) or []
    detail.xref_counts = dict(Counter(x.get("database", "?") for x in xrefs))
    for x in xrefs:
        db = x.get("database")
        if db == "PDB":
            detail.pdb_ids.append(x.get("id", ""))
        elif db == "GO":
            # x['properties'] is a list of {key,value}; aspect is "GoTerm" like "F:..."
            term_name = ""
            aspect = ""
            for p in x.get("properties", []) or []:
                if p.get("key") == "GoTerm":
                    val = p.get("value", "")
                    if ":" in val:
                        aspect, term_name = val.split(":", 1)
            detail.go_terms.append((x.get("id", ""), term_name.strip(), aspect.strip()))

    refs = data.get("references", []) or []
    for ref in refs:
        cit = ref.get("citation", {}) or {}
        for x in cit.get("citationCrossReferences", []) or []:
            if x.get("database") == "PubMed":
                detail.pubmed_ids.append(x.get("id", ""))

    detail.raw_keywords = [k.get("name", "") for k in data.get("keywords", []) or []]
    return detail
