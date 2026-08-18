"""Family presence/absence matrix + lineage-retention mining.

The core comparative-genomics signal for a paralog's fate: a gene that
is *present* in some lineages and *absent* in others is either a novel
family member with lineage-specific retention (e.g. a paralog kept in
monotreme/xenarthran, lost in most mammals) or an annotation gap worth
chasing. This module turns the current results table into that view:

    * rows    = family members (gene symbols; unnamed models grouped as
                homolog clusters when a distance matrix is supplied)
    * columns = species in the result set
    * cells   = highest-quality evidence found (protein count)

`retention_calls()` then classifies each row: core (present everywhere),
lineage-specific retention (present in ≥2 species but absent in ≥2), or
singleton. Pure post-processing — no network calls.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from ..core.models import ProteinVariant
from .distance import DistanceTable


@dataclass
class PresenceReport:
    members: list[str]                       # row labels
    species: list[str]                       # column labels
    counts: dict[tuple[str, str], int]       # (member, species) -> n variants
    calls: dict[str, str] = field(default_factory=dict)   # member -> verdict
    notes: list[str] = field(default_factory=list)

    def text_matrix(self) -> str:
        w = max((len(m) for m in self.members), default=8) + 2
        sp_short = [s.split()[0][:10] for s in self.species]
        colw = max(10, max((len(s) for s in sp_short), default=10) + 1)
        lines = [" " * w + "".join(f"{s:>{colw}}" for s in sp_short)]
        for m in self.members:
            row = "".join(
                f"{self.counts.get((m, s), 0) or '·':>{colw}}" for s in self.species
            )
            lines.append(f"{m:<{w}}" + row)
        return "\n".join(lines)

    def summary(self) -> str:
        lines = [
            f"Presence/absence matrix — {len(self.members)} family member(s) × "
            f"{len(self.species)} species",
            "",
            self.text_matrix(),
            "",
            "Retention calls:",
        ]
        for m in self.members:
            lines.append(f"  {m:24s} {self.calls.get(m, '?')}")
        lines += [""] + [f"Note: {n}" for n in self.notes]
        return "\n".join(lines)


def _member_label(v: ProteinVariant, cluster_of: dict[str, str]) -> str:
    """Group variants into family members. Named genes group by root symbol
    (itpr1a/itpr1b → ITPR1); unnamed models group by homolog cluster
    when available, else by their own id."""
    sym = v.gene_symbol
    if sym and not sym.startswith("ENS") and ":" not in sym:
        root = sym.upper().rstrip("ABX")
        return root if len(root) >= 3 else sym.upper()
    return cluster_of.get(f"{v.source}|{v.accession}|{v.gene_symbol}", sym)


def build_presence_matrix(
    variants: list[ProteinVariant],
    distances: Optional[DistanceTable] = None,
    analysis_label_for: Optional[dict[str, str]] = None,
    label_species: Optional[dict[str, str]] = None,
    cluster_identity: float = 0.35,
) -> PresenceReport:
    """Build the member × species matrix.

    When `distances` (+ label maps from the analysis step) are given,
    unnamed variants that are mutually ≥ `cluster_identity` identical are
    merged into one 'homolog cluster' row — so the armadillo/platypus/frog
    unnamed orthologs become a single family member instead of
    three singletons.
    """
    # Single-linkage clustering of unnamed variants via the distance table.
    cluster_of: dict[str, str] = {}
    if distances and analysis_label_for:
        unnamed = [v for v in variants
                   if v.gene_symbol.startswith("ENS") or ":" in v.gene_symbol]
        keyed = [(f"{v.source}|{v.accession}|{v.gene_symbol}", v) for v in unnamed]
        keyed = [(k, v) for k, v in keyed
                 if analysis_label_for.get(k) in (distances.labels or [])]
        parent: dict[str, str] = {k: k for k, _ in keyed}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i, (ka, _va) in enumerate(keyed):
            for kb, _vb in keyed[i + 1:]:
                la, lb = analysis_label_for[ka], analysis_label_for[kb]
                ia, ib = distances.labels.index(la), distances.labels.index(lb)
                if distances.identity[ia][ib] >= cluster_identity:
                    parent[find(ka)] = find(kb)
        groups: dict[str, list[str]] = defaultdict(list)
        for k, _v in keyed:
            groups[find(k)].append(k)
        for gi, (_root, members) in enumerate(sorted(groups.items()), 1):
            name = f"homolog_cluster_{gi}" if len(members) > 1 else None
            for k in members:
                cluster_of[k] = name or k.split("|")[1]

    counts: dict[tuple[str, str], int] = defaultdict(int)
    species_set: set[str] = set()
    members_set: set[str] = set()
    for v in variants:
        if not v.species:
            continue
        m = _member_label(v, cluster_of)
        counts[(m, v.species)] += 1
        species_set.add(v.species)
        members_set.add(m)

    species = sorted(species_set)
    members = sorted(members_set)
    report = PresenceReport(members=members, species=species, counts=dict(counts))

    # Retention calls
    n_sp = len(species)
    for m in members:
        present = [s for s in species if counts.get((m, s), 0)]
        absent = [s for s in species if not counts.get((m, s), 0)]
        if len(present) == n_sp and n_sp > 1:
            report.calls[m] = f"CORE — present in all {n_sp} species"
        elif len(present) >= 2 and len(absent) >= 2:
            report.calls[m] = (
                f"LINEAGE-SPECIFIC RETENTION — kept in {len(present)} "
                f"({', '.join(p.split()[0] for p in present[:4])}"
                + ("…" if len(present) > 4 else "") + f"), absent in {len(absent)} "
                f"({', '.join(a.split()[0] for a in absent[:4])}"
                + ("…" if len(absent) > 4 else "") + ") — lineage-specific retention"
            )
        elif len(present) == 1:
            report.calls[m] = (f"singleton — only {present[0]}; novel duplicate, "
                              "annotation gap in other genomes, or artifact")
        else:
            report.calls[m] = f"present in {len(present)}/{n_sp} species"

    if any(c.startswith("LINEAGE-SPECIFIC") for c in report.calls.values()):
        report.notes.append(
            "Lineage-specific retention is the signature of a paralog with a "
            "fate: verify absences with a BLAST bait search before "
            "claiming loss (they may be unannotated, not absent).")
    report.notes.append(
        "Counts are records in the current result table — a '·' can mean "
        "true absence, an unannotated genome, or simply a species the search "
        "didn't cover.")
    return report
