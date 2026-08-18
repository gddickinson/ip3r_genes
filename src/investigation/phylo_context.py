"""Place a candidate in the phylogenetic context of a reference panel.

Given the candidate sequence and a panel of canonical-family sequences
(ITPR1 = Q14643, ITPR2 = Q14571, ITPR3 = Q14573), build an MSA + NJ tree and
report:

    * the candidate's nearest panel member (by identity)
    * the candidate's bootstrap-style "topology classification":
        - "within paralog X" (clusters inside an existing panel branch)
        - "sister to paralog X" (sits as a sister branch to a paralog)
        - "outgroup" (deeper than the family's MRCA — surprising / suspect)
    * a Newick string for the panel + candidate so the user can open it in
      FigTree / iTOL
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests

from ..analysis.alignment import progressive_msa
from ..analysis.distance import identity_matrix
from ..analysis.tree import build_nj_tree, tree_to_ascii, tree_to_newick


UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb"


@dataclass
class PhyloContext:
    candidate_accession: str
    panel_labels: list[str] = field(default_factory=list)
    nearest_panel_label: str = ""
    nearest_panel_identity: float = 0.0
    pairwise_identities: dict[str, float] = field(default_factory=dict)
    newick: str = ""
    ascii_tree: str = ""
    topology_class: str = ""    # 'within' | 'sister' | 'outgroup' | 'unknown'

    def summary(self) -> str:
        lines = [f"Phylogenetic placement of {self.candidate_accession}:"]
        if self.nearest_panel_label:
            lines.append(f"  Nearest panel member: {self.nearest_panel_label} "
                         f"({self.nearest_panel_identity:.0%} identity)")
        if self.topology_class:
            lines.append(f"  Topology class: {self.topology_class}")
        if self.pairwise_identities:
            lines.append("  Pairwise identities to panel:")
            for k, v in sorted(self.pairwise_identities.items(), key=lambda x: -x[1]):
                lines.append(f"    {k}: {v:.0%}")
        return "\n".join(lines)


def _fetch_fasta(accession: str, timeout_s: int = 20) -> str:
    try:
        r = requests.get(f"{UNIPROT_FASTA}/{accession}.fasta", timeout=timeout_s)
        if r.status_code != 200:
            return ""
    except requests.RequestException:
        return ""
    return "".join(line.strip() for line in r.text.splitlines() if not line.startswith(">") and line.strip())


def build_phylo_context(
    candidate_accession: str,
    candidate_sequence: str,
    panel: list[tuple[str, str]],
    timeout_s: int = 20,
) -> PhyloContext:
    """`panel` is [(label, sequence_or_accession)]. If sequence is empty, we
    fetch it from UniProt by accession.
    """
    ctx = PhyloContext(candidate_accession=candidate_accession)

    # Resolve panel sequences
    rows: list[tuple[str, str]] = []
    for label, value in panel:
        seq = value
        if not seq or not any(c.isalpha() for c in seq[:10]):
            # treat as accession
            seq = _fetch_fasta(value, timeout_s)
        if seq:
            rows.append((label, seq))
    if not candidate_sequence:
        candidate_sequence = _fetch_fasta(candidate_accession, timeout_s)
    if not candidate_sequence:
        return ctx

    cand_label = f"CANDIDATE_{candidate_accession}"
    all_rows = rows + [(cand_label, candidate_sequence)]

    # MSA + identity matrix
    msa = progressive_msa(all_rows, use_mafft=False)
    if not msa:
        return ctx
    table = identity_matrix(msa)
    ctx.panel_labels = [l for l, _ in rows]

    # pairwise identities candidate ↔ each panel member
    cand_idx = next((i for i, lbl in enumerate(table.labels) if lbl == cand_label), -1)
    if cand_idx < 0:
        return ctx
    best = (None, 0.0)
    for j, lbl in enumerate(table.labels):
        if j == cand_idx:
            continue
        identity = table.identity[cand_idx][j]
        ctx.pairwise_identities[lbl] = identity
        if identity > best[1]:
            best = (lbl, identity)
    if best[0]:
        ctx.nearest_panel_label = best[0]
        ctx.nearest_panel_identity = best[1]

    # NJ tree
    tree = build_nj_tree(table)
    ctx.newick = tree_to_newick(tree)
    ctx.ascii_tree = tree_to_ascii(tree)

    # Topology class (heuristic): we look at the candidate's sister-clade
    # composition. If the sister is exactly one panel member → 'within';
    # if a clade of multiple panel members → 'sister'; if a single non-panel
    # leaf or an empty sister → 'outgroup'.
    ctx.topology_class = _classify_topology(tree, cand_label, set(ctx.panel_labels))
    return ctx


def _classify_topology(tree, candidate_label: str, panel_set: set[str]) -> str:
    """Walk the tree to find candidate's nearest internal node and inspect
    the sister clade's leaves.
    """
    target = None
    for term in tree.get_terminals():
        if term.name == candidate_label:
            target = term
            break
    if target is None:
        return "unknown"
    # Find parent path
    path = tree.get_path(target)
    if len(path) < 2:
        return "unknown"
    parent = path[-2]  # the parent clade
    sister_leaves: list[str] = []
    for child in parent.clades:
        if child is target:
            continue
        for term in child.get_terminals():
            if term.name:
                sister_leaves.append(term.name)
    if not sister_leaves:
        return "outgroup"
    panel_leaves = [l for l in sister_leaves if l in panel_set]
    if len(panel_leaves) == 1 and len(sister_leaves) == 1:
        return f"within {panel_leaves[0]}"
    if panel_leaves:
        return f"sister to {{{', '.join(sorted(panel_leaves))}}}"
    return "outgroup"
