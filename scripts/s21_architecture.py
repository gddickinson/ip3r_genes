"""The architecture itself: exon counts, coding length, span, intron sizes.

Four scope rules, each a positive test naming the number it fired on, because
an exon count is only architecture when the model that produced it is a whole
gene in an assembly that can hold one:

* **A1** the locus is one of S16's committed gene copies (`is_copy`) — its
  coverage clears 0.50, its identity the sweep's own floor and its alignment
  500 residues. Read back, never recomputed.
* **A2** the contig spans the gene (D4). Without it an exon count measures
  contig lengths: S19 measured 15.2 % of the sweep's cells as false negatives
  and every one is an assembly, so a `fragment` cell has fewer exons because
  the assembly is broken and not because the gene is.
* **A3** the model covers at least `COV_ARCH` of its bait. A model covering
  half its bait has half the exons; S16's copy bar is 0.50 and that is too
  loose for a count. The bar is swept in `sensitivity()` and the count barely
  moves, which is the reason it can be quoted.
* **A4** the bait reaches a usable frame (`s21_frame.MIN_BAIT_ON_REF`), which
  is what makes a boundary comparable between paralogues. A locus failing only
  A4 keeps its own counts and is dropped from the column-frame tests alone.

Cross-paralogue comparisons are **paired within genome** (D16): intron size,
assembly quality and annotation completeness all scale with the assembly, so an
unpaired comparison measures the assemblies. Every pair is a two-sided exact
sign test with ties dropped and counted (S8's rule) and the family of tests is
BH-corrected together (S9's rule).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s15_lib as S15                                             # noqa: E402
import s16_lib as S16                                             # noqa: E402
import s21_blocks as B                                            # noqa: E402
import s21_frame as FR                                            # noqa: E402
import s21_lib as L                                               # noqa: E402

#: A3's bar. Swept over `COV_BARS` in `sensitivity()`.
COV_ARCH = 0.90
COV_BARS = [0.50, 0.70, 0.80, 0.90, 0.95, 0.98]

#: The measurements every comparison is run on, and whether more is bigger.
METRICS = ("n_exons", "cds_bp", "span_bp", "median_intron_bp",
           "max_intron_bp", "mean_exon_bp", "total_intron_bp")

ARCH_COLS = ("accession", "organism", "vclass", "vorder", "cell", "locus_idx",
             "mp_id", "bait", "bait_paralog", "frame_via_cell", "reference",
             "strand", "contig", "start", "end", "coverage", "identity",
             "q_start", "q_end", "residues", "n_blocks", "n_exons", "n_merges", "n_introns",
             "frame_steps", "cds_bp", "span_bp", "total_intron_bp",
             "median_intron_bp", "max_intron_bp", "min_intron_bp_seen",
             "mean_exon_bp", "frame_balance_bp", "canonical_junctions",
             "minor_junctions", "non_canonical_junctions",
             "frac_canonical", "in_scope", "excluded_by")

INTRON_COLS = ("accession", "cell", "locus_idx", "intron_index", "length",
               "splice_class", "frame_ok", "col_left", "col_right",
               "intron_phase", "q_left", "q_right")


def scope_exclusion(row: dict, frames_ok: set[tuple[str, str]],
                    cov: float = COV_ARCH) -> str:
    """Why this locus is not an architecture — "" when it is one."""
    if row["cell"] not in L.CELLS:
        return f"A0_cell_{row['cell'] or 'unassigned'}"
    if not row["is_copy"]:
        return "A1_not_a_gene_copy"
    if not row["contig_spans_gene"]:
        return "A2_contig_cannot_span_the_gene"
    if row["coverage"] < cov:
        return f"A3_coverage_below_{cov}"
    if (row["bait"], row["cell"]) not in frames_ok:
        return "A4_no_usable_frame"
    return ""


def usable_frames(frames: dict) -> set[tuple[str, str]]:
    return {(r["bait"], r["cell"]) for r in FR.frame_rows(frames)
            if r["usable"]}


def intron_phase(gff_phase: int) -> int:
    """Classical intron phase from the following exon's GFF phase.

    GFF phase is the number of bases to drop from the start of the feature to
    reach a codon boundary, so a following exon at phase 0 means the intron
    falls between codons (phase 0), at phase 1 means two bases of the
    interrupted codon lie upstream (phase 2), and at phase 2 means one does
    (phase 1). Getting this inverted makes phases 1 and 2 swap, which changes
    every shared-intron call without changing any count.
    """
    return (3 - gff_phase) % 3


def build(models: list[dict], pairs: list[dict], frames: dict,
          min_intron: int) -> tuple[list[dict], list[dict]]:
    """The architecture rows and the intron rows, from the cached measurements.

    The block pairs are re-read rather than re-derived, so the splice class on a
    junction here is the same byte the calibration was made on.
    """
    ok = usable_frames(frames)
    by_locus: dict[tuple[str, str, int], list[dict]] = {}
    for p in pairs:
        by_locus.setdefault((p["accession"], p["cell"], p["locus_idx"]),
                            []).append(p)
    labels = L.bait_labels()
    arows: list[dict] = []
    irows: list[dict] = []
    for m in models:
        key = (m["accession"], m["cell"], m["locus_idx"])
        junc = sorted(by_locus.get(key, []), key=lambda p: p["pair_index"])
        cell = m["cell"]
        ref = L.FRAME_ACC.get(cell, "")
        meta = labels.get(m["bait"]) or {}
        via_cell = int(meta.get("family") == "ITPR" and not meta.get("paralog"))
        canon = sum(1 for p in junc if p["splice_class"] == "canonical")
        minor = sum(1 for p in junc if p["splice_class"] == "minor")
        noncan = sum(1 for p in junc if p["splice_class"] == "non_canonical")
        # The exon merge: with the calibrated floor at the smallest spliceable
        # gap, no junction is below it, so `n_merges` is 0 by measurement.
        merged = [p for p in junc if p["gap_bp"] < min_intron]
        n_exons = m["n_blocks"] - len(merged)
        ilen = [p["gap_bp"] for p in junc if p["gap_bp"] >= min_intron]
        excl = scope_exclusion(m, ok)
        arows.append({
            "accession": m["accession"], "organism": m["organism"],
            "vclass": m["vclass"], "vorder": m["vorder"], "cell": cell,
            "locus_idx": m["locus_idx"], "mp_id": m["mp_id"],
            "bait": m["bait"], "bait_paralog": m["bait_paralog"],
            "frame_via_cell": via_cell, "reference": ref,
            "strand": m["strand"], "contig": m["contig"], "start": m["start"],
            "end": m["end"], "coverage": m["coverage"],
            "identity": m["identity"], "q_start": m["q_start"],
            "q_end": m["q_end"], "residues": m["residues"],
            "n_blocks": m["n_blocks"], "n_exons": n_exons,
            "n_merges": len(merged), "n_introns": len(ilen),
            "frame_steps": m["frame_steps"], "cds_bp": m["cds_bp"],
            "span_bp": m["span_bp"], "total_intron_bp": sum(ilen),
            "median_intron_bp": round(S15.median([float(x) for x in ilen]), 1)
            if ilen else 0.0,
            "max_intron_bp": max(ilen) if ilen else 0,
            "min_intron_bp_seen": min(ilen) if ilen else 0,
            "mean_exon_bp": round(m["cds_bp"] / n_exons, 2) if n_exons else 0.0,
            "frame_balance_bp": m["frame_balance_bp"],
            "canonical_junctions": canon, "minor_junctions": minor,
            "non_canonical_junctions": noncan,
            "frac_canonical": round((canon + minor) / len(junc), 4)
            if junc else 0.0,
            "in_scope": int(not excl), "excluded_by": excl})
        if excl:
            continue
        for p in junc:
            if p["gap_bp"] < min_intron:
                continue
            irows.append({
                "accession": m["accession"], "cell": cell,
                "locus_idx": m["locus_idx"], "intron_index": p["pair_index"],
                "length": p["gap_bp"], "splice_class": p["splice_class"],
                "frame_ok": p["frame_ok"],
                "col_left": FR.column_of(frames, m["bait"], cell,
                                         p["prev_end_q"]),
                "col_right": FR.column_of(frames, m["bait"], cell,
                                          p["next_start_q"]),
                "intron_phase": intron_phase(p["next_phase_gff"]),
                "q_left": p["prev_end_q"], "q_right": p["next_start_q"]})
    return arows, irows


# --------------------------------------------------------------------------
# summaries
# --------------------------------------------------------------------------
def summarise(arows: list[dict], keys: tuple[str, ...]) -> list[dict]:
    """Median and spread of every metric, grouped."""
    groups: dict[tuple, list[dict]] = {}
    for r in arows:
        if not r["in_scope"]:
            continue
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    out = []
    for gk, rows in sorted(groups.items(), key=lambda kv: [str(x) for x in kv[0]]):
        d = {k: v for k, v in zip(keys, gk)}
        d["n_loci"] = len(rows)
        d["n_genomes"] = len({r["accession"] for r in rows})
        for met in METRICS:
            xs = [float(r[met]) for r in rows]
            d[f"{met}_median"] = round(S15.median(xs), 1)
            d[f"{met}_p10"] = round(L.quantile(xs, 0.10), 1)
            d[f"{met}_p90"] = round(L.quantile(xs, 0.90), 1)
            d[f"{met}_min"] = round(min(xs), 1)
            d[f"{met}_max"] = round(max(xs), 1)
        out.append(d)
    return out


def paired(arows: list[dict], stratum: str = "") -> list[dict]:
    """D16: every cross-paralogue comparison, paired within genome.

    One locus per (genome, cell) — the best-covered — because a genome carrying
    two ITPR1 copies would otherwise contribute two differences to one pair and
    weight itself twice.
    """
    best: dict[tuple[str, str], dict] = {}
    for r in arows:
        if not r["in_scope"]:
            continue
        if stratum and r["vclass"] != stratum:
            continue
        k = (r["accession"], r["cell"])
        if k not in best or r["coverage"] > best[k]["coverage"]:
            best[k] = r
    cells = list(L.PARALOGS) + [L.CONTROL_CELL]
    rows = []
    for i, a in enumerate(cells):
        for b in cells[i + 1:]:
            for met in METRICS:
                diffs = []
                for acc in {k[0] for k in best}:
                    ra, rb = best.get((acc, a)), best.get((acc, b))
                    if ra and rb:
                        diffs.append(float(ra[met]) - float(rb[met]))
                st = S15.sign_test(diffs)
                rows.append({
                    "stratum": stratum or "all", "cell_a": a, "cell_b": b,
                    "metric": met, "n_pairs": st["n"] + st["n_ties"],
                    "n_a_greater": st["n_pos"], "n_b_greater": st["n_neg"],
                    "n_ties": st["n_ties"],
                    "median_diff": round(S15.median(diffs), 1) if diffs else
                    float("nan"), "p": st["p"]})
    qs = S16.benjamini_hochberg([r["p"] if r["p"] == r["p"] else 1.0
                                 for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = round(q, 6)
    return rows


def sensitivity(models: list[dict], pairs: list[dict], frames: dict,
                min_intron: int) -> list[dict]:
    """Every count recomputed at each coverage bar — A3 measured, not assumed."""
    ok = usable_frames(frames)
    by_locus: dict[tuple[str, str, int], list[dict]] = {}
    for p in pairs:
        by_locus.setdefault((p["accession"], p["cell"], p["locus_idx"]),
                            []).append(p)
    rows = []
    for bar in COV_BARS:
        sel = [m for m in models if not scope_exclusion(m, ok, bar)]
        for cell in L.CELLS:
            sub = [m for m in sel if m["cell"] == cell]
            if not sub:
                continue
            ex = [float(m["n_blocks"]
                        - sum(1 for p in by_locus.get(
                            (m["accession"], m["cell"], m["locus_idx"]), [])
                            if p["gap_bp"] < min_intron)) for m in sub]
            rows.append({
                "cov_bar": bar, "cell": cell, "n_loci": len(sub),
                "n_genomes": len({m["accession"] for m in sub}),
                "n_exons_median": round(S15.median(ex), 1),
                "n_exons_p10": round(L.quantile(ex, 0.10), 1),
                "n_exons_p90": round(L.quantile(ex, 0.90), 1),
                "cds_bp_median": round(S15.median(
                    [float(m["cds_bp"]) for m in sub]), 1),
                "is_operating_point": int(bar == COV_ARCH)})
    return rows
