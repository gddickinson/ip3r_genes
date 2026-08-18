"""Lite synteny check via EBI Proteins API and NCBI Gene cross-references.

True synteny (à la MCScanX / SyntenyMap) needs whole-genome ortholog tables
that we don't have at hand. The lite version asks a much smaller question:

    For this UniProt accession, can we resolve a genomic location, and
    if so what are the immediately flanking genes on the same contig?

The flanking-gene set is the input the user takes to a separate synteny
tool (or just inspects manually). For oomycetes / arthropods etc. where
NCBI Gene neighbors aren't always populated, this often returns just the
genome assembly link without flanking genes — that's still useful to log.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests


PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"


@dataclass
class SyntenyContext:
    accession: str
    organism: str = ""
    genome_assembly: str = ""
    chromosome: str = ""
    location: str = ""
    flanking_genes: list[str] = field(default_factory=list)
    ncbi_gene_id: str = ""
    notes: str = ""


def fetch_synteny_context(accession: str, timeout_s: int = 30) -> Optional[SyntenyContext]:
    if not accession:
        return None
    try:
        r = requests.get(
            f"{PROTEINS_API}/coordinates/{accession}",
            headers={"Accept": "application/json"},
            timeout=timeout_s,
        )
    except requests.RequestException:
        return None
    if r.status_code != 200:
        # No coordinates available. Still record the gap as a SyntenyContext.
        return SyntenyContext(
            accession=accession,
            notes=f"EBI proteins coordinates endpoint returned {r.status_code}; "
                  f"genomic coordinates not available for this entry.",
        )
    try:
        data = r.json()
    except ValueError:
        return None
    if isinstance(data, list):
        data = data[0] if data else {}

    organism = (data.get("taxid", {}) or {}).get("commonName", "") or ""
    chrom = data.get("genomicLocation", {}).get("chromosome", "") if isinstance(data.get("genomicLocation"), dict) else ""

    ctx = SyntenyContext(
        accession=accession,
        organism=organism,
        chromosome=chrom or "",
        location="",
        genome_assembly=(data.get("genomicLocation", {}) or {}).get("assembly", "") if isinstance(data.get("genomicLocation"), dict) else "",
    )
    # The EBI coordinates schema can vary; we surface what we got and let
    # the case-file render decide how to display it.
    gn_locs = data.get("gnCoordinate", [])
    if isinstance(gn_locs, list) and gn_locs:
        first = gn_locs[0]
        ctx.genome_assembly = first.get("genomicLocation", {}).get("assemblyName", ctx.genome_assembly)
        ex = first.get("genomicLocation", {}).get("exon", [])
        if ex:
            start = ex[0].get("genomeLocation", {}).get("begin", {}).get("position")
            end   = ex[-1].get("genomeLocation", {}).get("end", {}).get("position")
            if start and end:
                ctx.location = f"{start}-{end}"
    return ctx
