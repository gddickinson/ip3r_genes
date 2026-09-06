"""s23_classify.py — copy number, not paralog cells.

S5's ledger has one row per genome x paralog: three named ITPR cells and a RyR
control, and every cell answers "is *this* paralog here". That model is a
statement about vertebrates. ITPR1/2/3 are a 2R product and S5b found the trio
absent below the cyclostomes, so outside the vertebrates the cells do not
exist and asking which of three is present is asking a question with no
referent.

The question here is **how many ITPRs this genome has**. That changes what a
locus is for, and three things follow.

**1. Loci are counted, not assigned.** Every locus won by the ITPR baits (D14,
positive test on alignment score — `s5_sweep_lib.Locus.family`) is graded by
how much of its bait it covers: `full`, `fragment` or `scrap`. The genome's
copy number is its count of `full` loci.

**2. Counting is deliberately conservative, because a split gene inflates it.**
Two loci a few kb apart on one strand, each covering a *different* part of the
bait, are far more likely one gene the aligner broke than two genes that
happen to be adjacent and complementary. S5b measured what happens when this
is ignored: before `MIN_LOCUS_IDENTITY` existed, *Lissotriton* ITPR3 read as
10 loci. Status was never wrong there — the best locus always won — but copy
number *is* the result in this task, so the error S5 could tolerate is the one
this task cannot. `merge_split_loci` folds them, records that it did, and the
ledger carries both the raw and the merged count so the merge is auditable
rather than invisible.

**3. Both controls are cells, and they mean different things.** The RyR cell
is a positive control only where a RyR is expected — inside Metazoa. The
**MIR cell is a control in every genome** (bait rule B2): a locus for the
MIR-domain sharer proves the search reached this assembly, in the clades where
the receptor's absence is the claim and RyR's absence is not evidence of
anything. And it is simultaneously a negative control for D14 — a MIR locus
called ITPR would be the sharpest possible failure of the family call, so the
ledger reports MIR-won loci explicitly rather than dropping them.

Statuses per genome (the brief's four, plus the annotation split S5 uses):
  found_annotated      >=1 full locus, overlapping an annotated gene
  found_unannotated    >=1 full locus, annotation exists but no gene there
  found_no_annotation  >=1 full locus, the assembly carries no annotation
  fragment_only        loci, but none full
  assembly_gap         fragment-only, and the best locus abuts a contig edge
                       or an N-run — the absence is the assembly's, not the
                       genome's
  no_locus             nothing the ITPR baits win
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s5_classify import _annotate_locus, _locus_flags, name_family  # noqa: E402
from s5_sweep_lib import COV_FOUND, Locus                           # noqa: E402

#: Locus grades by bait coverage. `COV_FOUND` (0.70) is S5's own "found" bar,
#: reused so a full locus here and a found cell there mean the same thing.
COV_FULL = COV_FOUND
#: Below this a locus is a `scrap` — the shared channel module matching
#: something, not a gene model. Set at the coverage S5b measured its junk
#: population topping out at (30 %), one level of evidence below `fragment`.
COV_SCRAP = 0.30

#: Merge rule (2). Two same-strand loci within this distance are candidates
#: for being one split gene. Ten times `s5_sweep_lib.LOCUS_GAP`, because the
#: clustering gap is what already failed to join them.
MERGE_GAP_BP = 100_000
#: ...and they are only merged if their bait query spans overlap by no more
#: than this fraction of the shorter one. Complementary halves of a bait are
#: one gene; two loci each covering the same half are two genes.
MERGE_MAX_QUERY_OVERLAP = 0.20

#: The three roles a bait can carry, matching `s23_bait_spec`.
FAMILY_ITPR = "ITPR"
FAMILY_RYR = "RYR"
CONTROL_MIR = "MIR"


# ------------------------------------------------------------------ merging

def _query_overlap_frac(a, b) -> float:
    lo = max(a.q_start, b.q_start)
    hi = min(a.q_end, b.q_end)
    if hi < lo:
        return 0.0
    shorter = min(a.q_end - a.q_start + 1, b.q_end - b.q_start + 1)
    return (hi - lo + 1) / shorter if shorter > 0 else 1.0


def merge_split_loci(loci: list[Locus], gap: int = MERGE_GAP_BP,
                     max_overlap: float = MERGE_MAX_QUERY_OVERLAP
                     ) -> tuple[list[Locus], list[dict]]:
    """Fold neighbouring loci that look like one split gene.

    Returns (merged loci, merge records). A merge record names both loci and
    the two numbers the decision was made on, so the ledger can show its
    working — the alternative is a copy number nobody can check.
    """
    by_key: dict[tuple, list[Locus]] = {}
    for L in loci:
        by_key.setdefault((L.contig, L.strand), []).append(L)
    out: list[Locus] = []
    merges: list[dict] = []
    for (contig, strand), group in by_key.items():
        group.sort(key=lambda L: L.start)
        cur = None
        for L in group:
            if cur is None:
                cur = Locus(contig, strand, L.start, L.end, list(L.alns))
                continue
            dist = L.start - cur.end
            ov = _query_overlap_frac(cur.best, L.best)
            if 0 <= dist <= gap and ov <= max_overlap:
                merges.append({
                    "contig": contig, "strand": strand,
                    "kept_start": cur.start, "kept_end": cur.end,
                    "merged_start": L.start, "merged_end": L.end,
                    "gap_bp": dist, "query_overlap": round(ov, 3),
                    "kept_bait": cur.best.bait, "merged_bait": L.best.bait,
                    "why": (f"{dist:,} bp apart on one strand and their bait "
                            f"spans overlap {ov:.0%} — complementary parts of "
                            "one gene, counted once")})
                cur.end = max(cur.end, L.end)
                cur.alns.extend(L.alns)
            else:
                out.append(cur)
                cur = Locus(contig, strand, L.start, L.end, list(L.alns))
        if cur is not None:
            out.append(cur)
    return out, merges


# ------------------------------------------------------------------ grading

def grade(coverage: float) -> str:
    if coverage >= COV_FULL:
        return "full"
    if coverage >= COV_SCRAP:
        return "fragment"
    return "scrap"


def locus_row(L: Locus, family: str) -> dict:
    a = L.best_of_family(family) or L.best
    return {
        "contig": L.contig, "strand": L.strand, "start": L.start, "end": L.end,
        "family": L.family, "bait": a.bait, "band": a.clade,
        "identity": round(a.identity, 4), "coverage": round(a.coverage, 4),
        "span_coverage": round(a.span_coverage, 4), "score": a.score,
        "aligned_aa": a.aligned_aa, "bait_len": a.bait_len,
        "q_span": [a.q_start, a.q_end], "frameshifts": a.frameshifts,
        "stop_codons": a.stop_codons, "mp_id": a.mp_id,
        "n_baits": len(L.alns), "grade": grade(a.coverage),
        "family_margin": L.family_margin(),
    }


def classify_genome(loci: list[Locus], gene_index, seqlens: dict, fetch_fn,
                    merge: bool = True) -> dict:
    """One genome's copy-number row, its control cells and its loci.

    `fetch_fn(contig, start, end) -> str` supplies sequence for the N-run and
    contig-edge flags, exactly as in S5.
    """
    itpr_raw = [L for L in loci if L.family == FAMILY_ITPR]
    ryr = [L for L in loci if L.family == FAMILY_RYR]
    mir = [L for L in loci if L.family == CONTROL_MIR]

    merges: list[dict] = []
    itpr = itpr_raw
    if merge:
        itpr, merges = merge_split_loci(itpr_raw)

    rows = []
    for L in sorted(itpr, key=lambda L: -L.best.score):
        d = locus_row(L, FAMILY_ITPR)
        region = fetch_fn(L.contig, L.start - 2000, L.end + 2000)
        d.update(_locus_flags(L, seqlens, region))
        _annotate_locus(d, gene_index, (L.best_of_family(FAMILY_ITPR)
                                        or L.best).cds_blocks, "ITPR1")
        # `_annotate_locus` takes a *cell* class; outside the vertebrates there
        # is none, so the paralog-match field it writes is meaningless here and
        # is replaced by the family-level question this task actually asks.
        d["annot_names_family"] = name_family(
            (d.get("annot_gene") or {}).get("name", "")) == "ITPR"
        d.pop("annot_paralog_matches", None)
        rows.append(d)

    full = [d for d in rows if d["grade"] == "full"]
    frag = [d for d in rows if d["grade"] == "fragment"]
    scrap = [d for d in rows if d["grade"] == "scrap"]

    if full:
        if gene_index is None:
            status = "found_no_annotation"
        else:
            status = ("found_annotated" if full[0].get("annot_gene")
                      else "found_unannotated")
    elif frag or scrap:
        best = (frag or scrap)[0]
        status = ("assembly_gap" if best.get("contig_edge") or best.get("n_gap")
                  else "fragment_only")
    else:
        status = "no_locus"

    return {
        "status": status,
        "n_full": len(full), "n_fragment": len(frag), "n_scrap": len(scrap),
        "n_loci_raw": len(itpr_raw), "n_loci_merged": len(itpr),
        "n_merges": len(merges),
        "best_coverage": rows[0]["coverage"] if rows else 0.0,
        "best_identity": rows[0]["identity"] if rows else 0.0,
        "best_family_margin": rows[0]["family_margin"] if rows else 0.0,
        "annotated_full": sum(1 for d in full if d.get("annot_gene")),
        "family_named_full": sum(1 for d in full if d["annot_names_family"]),
        # controls
        "ryr_loci": len(ryr),
        "ryr_full": sum(1 for L in ryr
                        if grade((L.best_of_family(FAMILY_RYR)
                                  or L.best).coverage) == "full"),
        "mir_loci": len(mir),
        "mir_full": sum(1 for L in mir
                        if grade((L.best_of_family(CONTROL_MIR)
                                  or L.best).coverage) == "full"),
        "loci": rows,
        "control_loci": [locus_row(L, L.family) for L in ryr + mir],
        "merges": merges,
    }


def control_verdict(row: dict, expects_ryr: bool,
                    has_control_bait: bool = True) -> tuple[str, str]:
    """Did the search demonstrably reach this assembly? (bait rule B2)

    Returns (verdict, why), in this order:

      `controlled_by_target`  the ITPR baits found a full locus. A control
                              exists to make a *negative* interpretable; a
                              genome where the family itself was found has
                              proved the search reached it, and asking for a
                              second proof would report *Dictyostelium* —
                              whose iplA this sweep recovers and matches to
                              its own annotation — as uncontrolled.
      `controlled`            a MIR-domain locus was found.
      `controlled_ryr_only`   no MIR locus, but a RyR one where RyR is
                              expected.
      `no_control_bait`       the panel carries no control bait for this
                              genome's clade, so nothing could have fired.
                              Distinct from `uncontrolled`: one is a silent
                              control, the other is an absent one, and only
                              the first is evidence about the assembly.
      `uncontrolled`          a control was available and none fired. No
                              absence claim may rest on this genome.

    S5b could write "the RyR control fired in all 309" only because it had a
    control that was *expected* in all 309. This is the same sentence made
    sayable outside the vertebrates.
    """
    if row.get("n_full", 0) > 0:
        return "controlled_by_target", (
            f"{row['n_full']} full ITPR locus/loci found — the search reached "
            "this assembly, so no separate proof-of-search is needed")
    mir_ok = row["mir_loci"] > 0
    ryr_ok = row["ryr_loci"] > 0
    if mir_ok:
        return "controlled", (
            f"{row['mir_loci']} MIR-domain locus/loci found"
            + (f" and {row['ryr_loci']} RyR" if ryr_ok else ""))
    if expects_ryr and ryr_ok:
        return "controlled_ryr_only", (
            f"no MIR locus but {row['ryr_loci']} RyR locus/loci — the search "
            "reached the assembly")
    if not has_control_bait and not expects_ryr:
        return "no_control_bait", (
            "the panel carries no MIR control bait for this genome's clade "
            "(controls are built for the clades carrying a negative claim), "
            "so no control could have fired here")
    return "uncontrolled", (
        "neither the MIR-domain control nor "
        + ("the RyR control" if expects_ryr else "an expected RyR")
        + " was found; no absence claim may rest on this genome")
