"""Novel-paralog discovery pipeline.

Glues together the analysis module (sequence outlier detection) with
structure / domain / orthology cross-checks to score and rank candidates
for being a *novel paralog* of the query gene family — e.g. a hypothetical
a fourth vertebrate ITPR.

Contents:
    candidates.py — composite scoring + ranked report
"""

from .candidates import (
    DiscoveryConfig,
    DiscoveryReport,
    Candidate,
    build_signature_set,
    discover_novel_paralogs,
    write_discovery,
)
from .domain_scan import run_domain_scan
from .exhaustive import run_exhaustive_hunt

__all__ = [
    "DiscoveryConfig", "DiscoveryReport", "Candidate", "build_signature_set",
    "discover_novel_paralogs", "write_discovery", "run_domain_scan",
    "run_exhaustive_hunt",
]
