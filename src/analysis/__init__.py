"""Sequence-analysis subsystem.

Operates on a list of ProteinVariant objects (each carrying a sequence).
Everything is independent of GUI / network — same code path serves the
GUI "Analyze" tab and the headless CLI `--analyze` flag.

Contents:
    alignment.py — pairwise alignment + simple progressive MSA
    distance.py  — identity / distance matrix on a variant set
    tree.py      — neighbor-joining tree → Newick + ASCII rendering
    clusters.py  — walk the tree to identify groups; flag novel candidates
    mutations.py — substitution / indel calling vs a chosen reference
    evolution.py — per-column conservation, variable-region detection

The orchestrator entry point is `analyse(variants, …) -> AnalysisResult`.
"""

from .alignment import progressive_msa, pairwise_align, AlignmentRow
from .distance import identity_matrix, DistanceTable
from .tree import build_nj_tree, tree_to_newick, tree_to_ascii
from .clusters import cluster_variants, ClusterReport
from .mutations import call_mutations, MutationCall
from .evolution import conservation_per_column, variable_regions
from .pipeline import analyse, AnalysisResult

__all__ = [
    "progressive_msa", "pairwise_align", "AlignmentRow",
    "identity_matrix", "DistanceTable",
    "build_nj_tree", "tree_to_newick", "tree_to_ascii",
    "cluster_variants", "ClusterReport",
    "call_mutations", "MutationCall",
    "conservation_per_column", "variable_regions",
    "analyse", "AnalysisResult",
]
