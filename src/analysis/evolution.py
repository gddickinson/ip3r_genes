"""Per-column conservation + variable-region detection on an MSA.

Conservation score: 1 - Shannon entropy normalized to the residue alphabet,
treating gaps as a 21st character. A score of 1.0 is fully conserved; 0.0
is maximally variable.

Variable regions: contiguous windows of columns where the average
conservation falls below a threshold — these are candidates for
functional-divergence sites (e.g. the IP3-binding core, which is where
ITPR1/2/3 differ in ligand affinity despite an identical pore).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from .alignment import AlignmentRow


_AA_PLUS_GAP = 21  # 20 standard amino acids + gap


@dataclass
class VariableRegion:
    start: int    # 1-based column
    end: int      # 1-based column, inclusive
    mean_conservation: float

    def width(self) -> int:
        return self.end - self.start + 1


def conservation_per_column(rows: Iterable[AlignmentRow]) -> list[float]:
    rows = list(rows)
    if not rows:
        return []
    width = max(len(r.aligned) for r in rows)
    n = len(rows)
    out: list[float] = []
    max_entropy = math.log(min(n, _AA_PLUS_GAP))
    for k in range(width):
        counts: dict[str, int] = {}
        for r in rows:
            c = r.aligned[k] if k < len(r.aligned) else "-"
            counts[c] = counts.get(c, 0) + 1
        h = 0.0
        for v in counts.values():
            p = v / n
            if p > 0:
                h -= p * math.log(p)
        cons = 1.0 - (h / max_entropy if max_entropy > 0 else 0.0)
        out.append(max(0.0, min(1.0, cons)))
    return out


def variable_regions(
    conservation: list[float],
    threshold: float = 0.4,
    min_width: int = 5,
) -> list[VariableRegion]:
    regions: list[VariableRegion] = []
    n = len(conservation)
    i = 0
    while i < n:
        if conservation[i] < threshold:
            j = i
            while j < n and conservation[j] < threshold:
                j += 1
            if j - i >= min_width:
                mean = sum(conservation[i:j]) / (j - i)
                regions.append(VariableRegion(start=i + 1, end=j, mean_conservation=mean))
            i = j
        else:
            i += 1
    return regions
