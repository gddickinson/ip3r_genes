"""Pairwise identity / distance computation.

`identity_matrix` consumes an MSA (list of AlignmentRow with equal-length
strings) and returns a `DistanceTable` where:
    * `identity[i][j]` ∈ [0, 1] — fraction of aligned positions where rows i
      and j have the same residue (gaps in either count as mismatch).
    * `distance[i][j]` = 1 - identity[i][j]

With `covered_only=True`, identity is computed over the mutually covered
columns only (positions where *both* rows have a residue). A fragment
aligned to a full-length protein is then scored on the region it actually
covers instead of being dragged toward 0 by the missing residues — the
fragment-aware mode S6 uses for the census identity matrix.

This is the canonical input for Bio.Phylo's DistanceTreeConstructor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .alignment import AlignmentRow


@dataclass
class DistanceTable:
    labels: list[str]
    identity: list[list[float]] = field(default_factory=list)
    distance: list[list[float]] = field(default_factory=list)

    def as_lower_triangle(self) -> list[list[float]]:
        """Return the distance matrix as a lower-triangular list-of-lists
        — the shape Bio.Phylo.TreeConstruction.DistanceMatrix expects.
        """
        n = len(self.labels)
        out: list[list[float]] = []
        for i in range(n):
            row = [self.distance[i][j] for j in range(i + 1)]
            out.append(row)
        return out


def identity_matrix(
    rows: Iterable[AlignmentRow], covered_only: bool = False,
) -> DistanceTable:
    rows = list(rows)
    labels = [r.label for r in rows]
    n = len(rows)
    width = max((len(r.aligned) for r in rows), default=0)
    identity = [[0.0] * n for _ in range(n)]
    distance = [[0.0] * n for _ in range(n)]
    for i in range(n):
        identity[i][i] = 1.0
        for j in range(i + 1, n):
            same = compared = 0
            ai, aj = rows[i].aligned, rows[j].aligned
            for k in range(width):
                ci = ai[k] if k < len(ai) else "-"
                cj = aj[k] if k < len(aj) else "-"
                if ci == "-" and cj == "-":
                    continue
                if covered_only and (ci == "-" or cj == "-"):
                    continue
                compared += 1
                if ci == cj:
                    same += 1
            iden = (same / compared) if compared else 0.0
            identity[i][j] = identity[j][i] = iden
            d = 1.0 - iden
            distance[i][j] = distance[j][i] = d
    return DistanceTable(labels=labels, identity=identity, distance=distance)
