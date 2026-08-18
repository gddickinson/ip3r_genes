"""Deep-dive investigation of a single candidate accession.

After the discovery scorer flags candidates (a fourth ITPR? an ITPR in a
lineage reported to lack one?), this
module gathers convergent evidence across independent lines of analysis to
let the user (or follow-up scripts) decide whether a candidate is:

    * a real novel paralog of the family,
    * a known protein in a different family that happens to share one domain,
    * a misannotated protein (the family it really belongs to is wrong on
      record) — the gap between 'scores well' and 'is a real gene'.

The seven lines of evidence (see roadmap §15):

    1. Full domain architecture        — domain_arch.py
    2. Detailed UniProt cross-refs     — uniprot_detail.py
    3. Predicted structure + pLDDT     — structure.py
    4. Structural homology (Foldseek)  — structure.py (opt-in, slow)
    5. Phylogenetic placement vs panel — phylo_context.py
    6. Genomic context / synteny       — synteny_lite.py
    7. Literature evidence             — literature.py

A pipeline orchestrator gathers all seven, then synthesis.py renders a
"case file" with a verdict and suggested follow-up experiments.
"""

from .pipeline import (
    Investigator,
    InvestigationResult,
    InvestigationOptions,
    write_investigation,
)
from .synthesis import render_case_file

__all__ = [
    "Investigator",
    "InvestigationResult",
    "InvestigationOptions",
    "write_investigation",
    "render_case_file",
]
