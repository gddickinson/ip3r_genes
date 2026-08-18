"""Cluster variants into groups; flag novel-group candidates.

Two complementary outputs:

1. **Hierarchical clusters** derived by cutting the NJ tree at a chosen
   distance threshold. Each cluster is interpreted as a putative group of
   close homologs — useful for separating ITPR1-like vs ITPR2-like vs
   ITPR3-like sequences, or for spotting un-annotated paralogs.

2. **Novel-candidate flags** — variants whose nearest neighbor in the
   distance matrix is below a user-set identity threshold (default 50 %).
   These are the "odd ones out" that may be (a) misannotated, (b) members
   of a yet-unrecognized subfamily, or (c) true distant homologs that
   warrant manual review.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .distance import DistanceTable


@dataclass
class ClusterReport:
    threshold: float
    clusters: list[list[str]] = field(default_factory=list)
    novel_candidates: list[tuple[str, str, float]] = field(default_factory=list)
    # ^ (label, nearest_other_label, best_identity_seen)

    def summary(self) -> str:
        lines = [
            f"Clustering threshold: identity ≥ {self.threshold:.0%}",
            f"Found {len(self.clusters)} cluster(s).",
        ]
        for i, c in enumerate(self.clusters, 1):
            lines.append(f"  Cluster {i} ({len(c)} member{'s' if len(c) != 1 else ''}):")
            for m in c:
                lines.append(f"    • {m}")
        if self.novel_candidates:
            lines.append("")
            lines.append(f"Novel-group candidates (nearest neighbor < {self.threshold:.0%} identity):")
            for label, neighbor, ident in self.novel_candidates:
                lines.append(f"  • {label}  → closest: {neighbor} ({ident:.0%})")
        return "\n".join(lines)


def cluster_variants(table: DistanceTable, identity_threshold: float = 0.5) -> ClusterReport:
    """Single-link clustering by identity threshold.

    Two variants are grouped together if their pairwise identity ≥
    `identity_threshold`. This is intentionally simple — it captures
    "these belong to the same protein family" without needing a Newick
    parser. For deeper analysis (HOGs, time-slice grouping), use the OMA
    client.
    """
    n = len(table.labels)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            if table.identity[i][j] >= identity_threshold:
                union(i, j)

    groups: dict[int, list[str]] = {}
    for i, label in enumerate(table.labels):
        groups.setdefault(find(i), []).append(label)
    clusters = sorted(groups.values(), key=lambda g: -len(g))

    novel: list[tuple[str, str, float]] = []
    for i, label in enumerate(table.labels):
        best_j, best_id = -1, 0.0
        for j in range(n):
            if j == i:
                continue
            if table.identity[i][j] > best_id:
                best_id = table.identity[i][j]
                best_j = j
        if best_j >= 0 and best_id < identity_threshold:
            novel.append((label, table.labels[best_j], best_id))

    return ClusterReport(
        threshold=identity_threshold,
        clusters=clusters,
        novel_candidates=novel,
    )
