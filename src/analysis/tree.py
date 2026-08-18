"""Phylogenetic tree construction via Bio.Phylo.

Input: a `DistanceTable`. Output: a Bio.Phylo Tree object plus convenience
renderers — Newick string (for FigTree / iTOL / Dendroscope) and an ASCII
tree (for the GUI's text-only analysis pane and for terminal CLI output).

Default constructor is Neighbor-Joining (the same first-pass tool Dong et
standard route to a first placement). UPGMA is offered as an alternative
when a clock-like model is appropriate.
"""

from __future__ import annotations

import io
from typing import Literal

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from .distance import DistanceTable


def _to_biopython_matrix(table: DistanceTable) -> DistanceMatrix:
    return DistanceMatrix(names=list(table.labels), matrix=table.as_lower_triangle())


def build_nj_tree(table: DistanceTable, method: Literal["nj", "upgma"] = "nj") -> Tree:
    if len(table.labels) < 3:
        # Degenerate — Bio.Phylo NJ needs ≥3 taxa. Build a trivial pair tree.
        return _trivial_tree(table.labels)
    matrix = _to_biopython_matrix(table)
    constructor = DistanceTreeConstructor()
    if method == "upgma":
        return constructor.upgma(matrix)
    return constructor.nj(matrix)


def _trivial_tree(labels: list[str]) -> Tree:
    from Bio.Phylo.BaseTree import Clade
    if not labels:
        return Tree(root=Clade(name="empty"))
    if len(labels) == 1:
        return Tree(root=Clade(name=labels[0]))
    root = Clade(branch_length=0.0)
    root.clades = [Clade(branch_length=0.5, name=l) for l in labels]
    return Tree(root=root)


def tree_to_newick(tree: Tree) -> str:
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()


def tree_to_ascii(tree: Tree, width: int = 80) -> str:
    buf = io.StringIO()
    Phylo.draw_ascii(tree, file=buf, column_width=width)
    return buf.getvalue()
