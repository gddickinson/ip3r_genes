"""The tandem-duplication test, and the control that proves it can fire.

miniprot aligns each bait independently, so if a genome encodes the same part of
a protein twice the *same* bait aligns twice, at two disjoint genomic places.
That geometry is the signature, and it is read here off the sweep's own
alignments (S19's per-accession cache, so nothing is re-aligned).

**The geometry alone is not the test, and this is where the family bites.** An
ITPR2 bait aligns at the ITPR1 and ITPR3 genes as well as its own — the three
paralogues are 61-68 % identical — so "the same bait aligns twice at disjoint
positions" fires in essentially every vertebrate genome and measures paralogy,
not duplication. Run that way the detector reached a specificity of 0.16. So the
pair has to be inside the cell's **own** loci: the alignments are clustered and
attributed exactly as the sweep does it (`cluster_loci` → `filter_loci` →
`s5_classify.cell_loci`, which is where D14 lives), and a hit at another
paralogue's gene therefore belongs to that paralogue's locus and is not offered
here.

Three classes, because the distance is the biology:

* **`within_locus`** — both alignments inside one S5 locus cluster, so part of
  the protein is encoded twice inside one gene's span: either an internal
  duplication, or two neighbouring genes the 10 kb clustering merged.
* **`tandem`** — two loci of the same cell, same contig, within
  `TANDEM_GAP_BP`: two genes side by side.
* **`dispersed`** — two loci further apart, or on different contigs.

**The control is the point.** A detector that never fires and a family with no
duplicates look identical, so the detector is scored as a classifier of a copy
count it never sees: S16's committed `n_copies` per genome × cell, decided by
coverage, identity and aligned length and by no pairwise geometry at all.
Sensitivity and specificity are written out whatever they are, the 3R teleost
groups are checked separately as the one place two copies are expected for a
known reason, and `s21_test_arch` constructs a duplicate pair and a single
alignment and requires the rule to fire on the first and not the second.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as S16                                             # noqa: E402
import s19_panel as PANEL                                         # noqa: E402
import s21_lib as L                                               # noqa: E402
import s5_classify as CL                                          # noqa: E402
import s5_sweep_lib as SW                                         # noqa: E402

#: An alignment covering less than this fraction of its bait is not evidence of
#: a copy: a shared module hit at 5 % coverage aligns twice all over a genome.
MIN_ALN_COV = 0.30

#: How much of the shorter of two query spans must overlap for the two
#: alignments to be "the same part of the protein". Half is deliberately
#: generous: a real duplicate pair usually overlaps almost completely, and a
#: tighter bar would make the detector's failures look like absences.
Q_OVERLAP = 0.50

#: Same-contig separation under which a pair is called tandem rather than
#: dispersed. The widest ITPR gene span measured in the sweep is under 1 Mb, so
#: 1 Mb keeps a pair of neighbouring genes together without merging two genes
#: that are megabases apart.
TANDEM_GAP_BP = 1_000_000

PAIR_COLS = ("accession", "organism", "vclass", "cell", "bait", "clade",
             "family", "class", "contig_a", "start_a", "end_a", "strand_a",
             "cov_a", "contig_b", "start_b", "end_b", "strand_b", "cov_b",
             "q_overlap_frac", "genomic_gap_bp", "same_contig", "same_strand",
             "within_one_locus")

CELL_COLS = ("accession", "organism", "vclass", "cell", "detector_fires",
             "n_pairs", "n_within_locus", "n_tandem", "n_dispersed",
             "s16_copies", "s16_same_contig_copies", "truth_duplicated",
             "agreement")


def _q_overlap(a, b) -> float:
    lo = max(a.q_start, b.q_start)
    hi = min(a.q_end, b.q_end)
    if hi < lo:
        return 0.0
    shorter = min(a.q_end - a.q_start + 1, b.q_end - b.q_start + 1)
    return (hi - lo + 1) / shorter if shorter else 0.0


def _disjoint(a, b) -> bool:
    return a.contig != b.contig or a.end < b.start or b.end < a.start


def _gap(a, b) -> int:
    if a.contig != b.contig:
        return -1
    return max(0, max(a.start, b.start) - min(a.end, b.end) - 1)


def pair_class(a, b, same_locus: bool) -> str:
    if same_locus:
        return "within_locus"
    if a.contig == b.contig and _gap(a, b) <= TANDEM_GAP_BP:
        return "tandem"
    return "dispersed"


def cell_alignments(alns: list) -> dict[str, list[tuple]]:
    """`cell -> [(locus_index, alignment)]`, by the sweep's own attribution.

    The clustering, the identity floor and `cell_loci` are the sweep's, imported
    unchanged: a locus won by the RyR baits is never offered to an ITPR cell,
    and only inside the winning family is the paralogue question asked (D14).
    """
    loci = SW.filter_loci(SW.cluster_loci(alns))
    out: dict[str, list[tuple]] = {}
    for cell in L.CELLS:
        primary, secondary = CL.cell_loci(loci, cell)
        mine = primary + secondary
        for li, locus in enumerate(mine):
            for a in locus.alns:
                out.setdefault(cell, []).append((li, a))
    return out


def detect(alns: list, min_cov: float = MIN_ALN_COV,
           q_overlap: float = Q_OVERLAP) -> list[tuple]:
    """Every duplicate-geometry pair inside one cell's own loci.

    Returns `(cell, bait, aln_a, aln_b, q_overlap, same_locus)`.
    """
    out = []
    for cell, items in sorted(cell_alignments(alns).items()):
        by_bait: dict[str, list[tuple]] = {}
        for li, a in items:
            if a.coverage >= min_cov:
                by_bait.setdefault(a.bait, []).append((li, a))
        for bait, group in sorted(by_bait.items()):
            group = sorted(group, key=lambda x: (x[1].contig, x[1].start))
            for i, (la, a) in enumerate(group):
                for lb, b in group[i + 1:]:
                    if not _disjoint(a, b):
                        continue
                    ov = _q_overlap(a, b)
                    if ov < q_overlap:
                        continue
                    out.append((cell, bait, a, b, ov, la == lb))
    return out


def run(on_progress=None) -> tuple[list[dict], list[dict]]:
    """Detector pairs and the per-cell control, over the whole sweep."""
    meta = {r["id"]: {"clade": r["clade"], "family": r["family"],
                      "length": int(float(r["length"]))}
            for r in L.read_tsv(L.BAIT_MANIFEST)}
    copies = {(r["accession"], r["cell"]): r
              for r in L.read_tsv(L.RESULTS / "duplication" / "copy_number.tsv")}
    summaries = L.load_summaries()
    accs = sorted(summaries)
    prows: list[dict] = []
    fires: dict[tuple[str, str], dict] = {}
    for i, acc in enumerate(accs, 1):
        s = summaries[acc]
        alns = PANEL.genome_alignments(acc, meta)
        for cell, bait, a, b, ov, same_locus in detect(alns):
            cls = pair_class(a, b, same_locus)
            prows.append({
                "accession": acc, "organism": s.get("organism", ""),
                "vclass": s.get("vclass", ""), "cell": cell,
                "bait": bait.split("|")[0], "clade": a.clade,
                "family": a.family, "class": cls,
                "contig_a": a.contig, "start_a": a.start, "end_a": a.end,
                "strand_a": a.strand, "cov_a": round(a.coverage, 4),
                "contig_b": b.contig, "start_b": b.start, "end_b": b.end,
                "strand_b": b.strand, "cov_b": round(b.coverage, 4),
                "q_overlap_frac": round(ov, 4), "genomic_gap_bp": _gap(a, b),
                "same_contig": int(a.contig == b.contig),
                "same_strand": int(a.strand == b.strand),
                "within_one_locus": int(same_locus)})
            d = fires.setdefault((acc, cell),
                                 {"n_pairs": 0, "n_within_locus": 0,
                                  "n_tandem": 0, "n_dispersed": 0})
            d["n_pairs"] += 1
            d[f"n_{cls}"] += 1
        if on_progress and (i % 50 == 0 or i == len(accs)):
            on_progress(i, len(accs))
    crows = []
    for acc in accs:
        s = summaries[acc]
        for cell in L.CELLS:
            cp = copies.get((acc, cell))
            if cp is None:
                continue
            n_copies = int(float(cp.get("n_copies") or 0))
            d = fires.get((acc, cell), {"n_pairs": 0, "n_within_locus": 0,
                                        "n_tandem": 0, "n_dispersed": 0})
            # A `within_locus` pair is a duplication inside one gene's span, not
            # a second gene, so it is deliberately not counted as the detector
            # firing on copy number — the control would otherwise score the
            # detector against a claim it is not making.
            fired = int(d["n_tandem"] + d["n_dispersed"] > 0)
            truth = int(n_copies >= 2)
            crows.append({
                "accession": acc, "organism": s.get("organism", ""),
                "vclass": s.get("vclass", ""), "cell": cell,
                "detector_fires": fired, **d, "s16_copies": n_copies,
                "s16_same_contig_copies": int(float(
                    cp.get("same_contig_copies") or 0)),
                "truth_duplicated": truth,
                "agreement": "tp" if fired and truth else
                             "fn" if truth else
                             "fp" if fired else "tn"})
    return prows, crows


def control_summary(crows: list[dict]) -> list[dict]:
    """Sensitivity and specificity of the detector against S16's copy call."""
    out = []
    for cell in ["all"] + list(L.CELLS):
        sub = [r for r in crows if cell == "all" or r["cell"] == cell]
        tp = sum(1 for r in sub if r["agreement"] == "tp")
        fn = sum(1 for r in sub if r["agreement"] == "fn")
        fp = sum(1 for r in sub if r["agreement"] == "fp")
        tn = sum(1 for r in sub if r["agreement"] == "tn")
        out.append({
            "cell": cell, "n_cells": len(sub), "tp": tp, "fn": fn, "fp": fp,
            "tn": tn,
            "sensitivity": round(tp / (tp + fn), 4) if tp + fn else float("nan"),
            "specificity": round(tn / (tn + fp), 4) if tn + fp else float("nan"),
            "n_truth_duplicated": tp + fn,
            "can_fire": int(tp + fp > 0)})
    return out


def class_summary(prows: list[dict]) -> list[dict]:
    """Pair counts by class and cell, plus how many genomes each reaches."""
    out = []
    for cls in ("within_locus", "tandem", "dispersed"):
        for cell in ["all"] + list(L.CELLS):
            sub = [r for r in prows if r["class"] == cls
                   and (cell == "all" or r["cell"] == cell)]
            if not sub:
                continue
            out.append({
                "class": cls, "cell": cell, "n_pairs": len(sub),
                "n_genomes": len({r["accession"] for r in sub}),
                "median_q_overlap": round(L.quantile(
                    [r["q_overlap_frac"] for r in sub], 0.5), 3),
                "median_gap_bp": round(L.quantile(
                    [float(r["genomic_gap_bp"]) for r in sub], 0.5), 1),
                "frac_same_strand": round(
                    sum(r["same_strand"] for r in sub) / len(sub), 4)})
    return out


def teleost_check(crows: list[dict]) -> list[dict]:
    """The 3R control: does the detector fire where S16 placed a 3R duplicate?

    S16's `teleost_copies.tsv` is the independent statement — three groups
    defined by which whole-genome duplications a lineage has been through — so
    the detector agreeing with it there is a check on the detector in the one
    place two copies are expected for a known reason.
    """
    tel = L.read_tsv(L.TELEOST)
    idx = {(r["accession"], r["cell"]): r for r in crows}
    out = []
    for t in tel:
        for cell in L.PARALOGS:
            n = int(float(t.get(f"{cell}_copies") or 0))
            c = idx.get((t["accession"], cell))
            if c is None:
                continue
            out.append({
                "accession": t["accession"], "organism": t["organism"],
                "group": t.get("group", ""), "cell": cell,
                "s16_copies": n, "detector_fires": c["detector_fires"],
                "n_tandem": c["n_tandem"], "n_dispersed": c["n_dispersed"],
                "expected": int(n >= 2),
                "agrees": int(c["detector_fires"] == int(n >= 2))})
    return out


def bh(pvals: list[float]) -> list[float]:
    return S16.benjamini_hochberg(pvals)
