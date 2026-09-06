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
**control cell fires in every genome** (bait rule B2): a locus won by a
control bait proves the search reached this assembly, in the clades where the
receptor's absence is the claim and RyR's absence is not evidence of anything.
Which protein family that control is drawn from is measured per clade
(`s23_control_select`), and the MIR-domain sharer stays in the panel whatever
the measurement says, because it is the only control inside the family's own
signature set: a control locus called ITPR would be the sharpest possible
failure of the family call, so the ledger reports control-won loci explicitly
rather than dropping them.

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

import s23_calibration as cal                                       # noqa: E402
from s5_classify import _annotate_locus, _locus_flags, name_family  # noqa: E402
from s5_sweep_lib import COV_FOUND, Aln, Locus                      # noqa: E402

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

#: Kingdom-level buckets, so a bait's origin and a genome's group can be
#: compared. The two vocabularies differ by history — the bait's `group` is
#: the S20 reference-proteome DB it came out of, the genome's is S4's
#: manifest column — and the cross-kingdom control tier is a comparison
#: between them, so it needs one bucket both map into.
KINGDOM_GROUP = {
    "metazoa": "metazoa", "metazoa_nonvert": "metazoa",
    "fungi": "fungi", "viridiplantae": "viridiplantae",
    "protista_other": "protist", "sar": "protist",
    "amoebozoa": "protist", "discoba": "protist", "other": "protist",
}


def kingdom_group(g: str) -> str:
    return KINGDOM_GROUP.get((g or "").strip(), (g or "").strip())


#: The three roles a bait can carry, matching `s23_bait_spec`.
FAMILY_ITPR = "ITPR"
FAMILY_RYR = "RYR"
CONTROL_ROLE = "CONTROL"
CONTROL_MIR = CONTROL_ROLE


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


# ------------------------------------------------------------ copy counting

def cds_footprint_bp(a: Aln) -> int:
    """The genomic bases the alignment's exons actually occupy."""
    return sum(e - s + 1 for s, e in (a.cds_blocks or [])) or (a.end - a.start + 1)


def _overlap_frac(a: Aln, b: Aln) -> float:
    ov = min(a.end, b.end) - max(a.start, b.start) + 1
    if ov <= 0:
        return 0.0
    shorter = min(a.end - a.start + 1, b.end - b.start + 1)
    return ov / shorter if shorter > 0 else 1.0


def full_copies(loci: list[Locus], family: str,
                cov_full: float = None,
                max_overlap: float = None) -> list[Aln]:
    """The distinct complete gene models the sweep found, genome-wide.

    **Why this is not a count of loci.** A locus is a cluster of alignments
    chained at `LOCUS_GAP`, and outside the vertebrates `-G` is 650 kb for the
    metazoa, so a cluster is a large object: in the S23a pilot *Drosophila*'s
    single 22 kb Itpr gene sat inside a locus spanning 297 kb — 13x the gene —
    because 26 alignments from other baits chained across it. The status call
    was right (the best alignment covered the gene's CDS completely), but two
    real genes inside one such chain would have been counted **once**, and
    copy number is what this task delivers.

    So a copy is defined on the evidence rather than on the cluster: a genomic
    interval carrying an alignment that covers a complete bait. Two alignments
    that overlap are the same copy — several baits hitting one gene is the
    normal case, not two genes — and two that do not are two copies, however
    the clustering happened to group them. Greedy by score, so the copy is
    anchored on the best evidence for it.
    """
    if cov_full is None:
        cov_full = COV_FULL
    if max_overlap is None:
        max_overlap = cal.copy_max_overlap()[0]
    cands = [a for L in loci for a in L.alns
             if a.family == family and a.coverage >= cov_full]
    cands.sort(key=lambda a: (-a.score, a.contig, a.start))
    kept: list[Aln] = []
    for a in cands:
        if any(a.contig == b.contig and _overlap_frac(a, b) > max_overlap
               for b in kept):
            continue
        kept.append(a)
    return sorted(kept, key=lambda a: (a.contig, a.start))


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


def _describe(L: Locus, family: str, gene_index) -> dict:
    """A locus row with its annotation attached, but no sequence fetched."""
    d = locus_row(L, family)
    best = L.best_of_family(family) or L.best
    _annotate_locus(d, gene_index, best.cds_blocks, "ITPR1")
    # `_annotate_locus` takes a *cell* class; outside the vertebrates there is
    # none, so the paralog-match field it writes is meaningless here and is
    # replaced by the family-level question this task actually asks.
    d["annot_names_family"] = name_family(
        (d.get("annot_gene") or {}).get("name", "")) == "ITPR"
    d.pop("annot_paralog_matches", None)
    d["span_bp"] = L.end - L.start + 1
    d["cds_footprint_bp"] = cds_footprint_bp(best)
    d["span_inflation"] = round(d["span_bp"] / max(1, d["cds_footprint_bp"]), 2)
    return d


def cross_group_support(loci: list[Locus], family: str, own_group: str,
                        bait_meta: dict | None) -> list[str]:
    """Which *other* kingdom-level groups' baits also align across these loci.

    The second tier of the control, and it costs nothing to collect: the panel
    carries a control bait per control clade and searches all of them in every
    genome, so a fungal genome is already being probed with plant, protist and
    metazoan control baits alongside its own.

    That matters because the two tiers answer different questions. A control
    bait from the genome's own clade firing shows the **assembly is
    searchable**. A bait from another kingdom firing on the same gene shows
    the method **crosses kingdom-level divergence in this assembly** — which
    is the thing an absence claim actually needs, since any ITPR here would
    have to be found by baits at least that far away.
    """
    if not bait_meta or not own_group:
        return []
    mine = kingdom_group(own_group)
    groups = set()
    for L in loci:
        for a in L.alns:
            if a.family != family or a.coverage < COV_SCRAP:
                continue
            g = kingdom_group((bait_meta.get(a.bait) or {}).get("group", ""))
            if g and g != mine:
                groups.add(g)
    return sorted(groups)


def classify_genome(loci: list[Locus], gene_index, seqlens: dict, fetch_fn,
                    merge: bool = True, call_floor: float | None = None,
                    bait_meta: dict | None = None, own_group: str = ""
                    ) -> dict:
    """One genome's copy-number row, its control cells and its loci.

    `fetch_fn(contig, start, end) -> str` supplies sequence for the N-run and
    contig-edge flags, exactly as in S5.

    `loci` arrives filtered only at the **recording** floor
    (`s23_calibration.RECORD_MIN_IDENTITY`). The **call** floor is applied
    here, and the clusters it excludes are kept in the summary as
    `below_floor_loci` — they are the population the floor is calibrated
    against, and a sweep that filtered them out at search time could not
    measure its own threshold.
    """
    if call_floor is None:
        call_floor, _ = cal.call_min_identity()
    itpr_all = [L for L in loci if L.family == FAMILY_ITPR]
    itpr_raw = [L for L in itpr_all if L.best.identity >= call_floor]
    below = [L for L in itpr_all if L.best.identity < call_floor]
    # The control cells are held to the *same* floor as the family, and to the
    # same grade bar. A control that can be satisfied by evidence too weak to
    # be called a gene is not a control: it would let a genome be declared
    # `controlled` — and its silence about the receptor read as biology — on a
    # chained fragment of a mannosyltransferase. Recording at
    # `RECORD_MIN_IDENTITY` made this reachable, so it is closed here.
    ryr_all = [L for L in loci if L.family == FAMILY_RYR]
    ctl_all = [L for L in loci if L.family == CONTROL_ROLE]

    def _called(pool: list[Locus], family: str) -> list[Locus]:
        return [L for L in pool
                if L.best.identity >= call_floor
                and grade((L.best_of_family(family) or L.best).coverage)
                != "scrap"]

    ryr = _called(ryr_all, FAMILY_RYR)
    ctl = _called(ctl_all, CONTROL_ROLE)

    merges: list[dict] = []
    itpr = itpr_raw
    if merge:
        itpr, merges = merge_split_loci(itpr_raw)

    rows = []
    for L in sorted(itpr, key=lambda L: -L.best.score):
        d = _describe(L, FAMILY_ITPR, gene_index)
        region = fetch_fn(L.contig, L.start - 2000, L.end + 2000)
        d.update(_locus_flags(L, seqlens, region))
        rows.append(d)

    below_rows = [_describe(L, FAMILY_ITPR, gene_index)
                  for L in sorted(below, key=lambda L: -L.best.score)[:40]]

    copies = full_copies(itpr, FAMILY_ITPR)
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
        "call_min_identity": round(call_floor, 4),
        # `n_full` is the copy number: distinct non-overlapping complete gene
        # models (`full_copies`), not the count of full-graded clusters, which
        # is kept beside it as `n_full_loci` so the difference is visible.
        "n_full": len(copies), "n_full_loci": len(full),
        "n_fragment": len(frag), "n_scrap": len(scrap),
        "n_below_floor": len(below),
        "best_below_floor_identity": (round(max(
            (L.best.identity for L in below), default=0.0), 4)),
        "copies": [{"contig": a.contig, "start": a.start, "end": a.end,
                    "strand": a.strand, "bait": a.bait,
                    "identity": round(a.identity, 4),
                    "coverage": round(a.coverage, 4), "score": a.score,
                    "cds_footprint_bp": cds_footprint_bp(a)} for a in copies],
        "below_floor_loci": below_rows,
        "n_loci_raw": len(itpr_raw), "n_loci_merged": len(itpr),
        "n_merges": len(merges),
        "best_coverage": rows[0]["coverage"] if rows else 0.0,
        "best_identity": rows[0]["identity"] if rows else 0.0,
        "best_family_margin": rows[0]["family_margin"] if rows else 0.0,
        "annotated_full": sum(1 for d in full if d.get("annot_gene")),
        "family_named_full": sum(1 for d in full if d["annot_names_family"]),
        # controls
        "ryr_loci_recorded": len(ryr_all),
        "control_loci_recorded": len(ctl_all),
        "ryr_loci": len(ryr),
        "ryr_full": sum(1 for L in ryr
                        if grade((L.best_of_family(FAMILY_RYR)
                                  or L.best).coverage) == "full"),
        "control_loci": len(ctl),
        "control_full": sum(1 for L in ctl
                            if grade((L.best_of_family(CONTROL_ROLE)
                                      or L.best).coverage) == "full"),
        "control_cross_groups": cross_group_support(
            ctl, CONTROL_ROLE, own_group, bait_meta),
        "loci": rows,
        "control_locus_rows": [locus_row(L, L.family)
                               for L in ryr + ctl],
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
      `controlled_cross_kingdom`
                              a control locus was found, **and** a control
                              bait from another kingdom aligns across it. The
                              strongest form: the search demonstrably crosses
                              kingdom-level divergence in this assembly, which
                              is the distance any receptor here would have to
                              be found across.
      `controlled`            a control locus was found, from this genome's
                              own clade only.
      `controlled_ryr_only`   no control locus, but a RyR one where RyR is
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
    ctl_ok = row["control_loci"] > 0
    ryr_ok = row["ryr_loci"] > 0
    cross = row.get("control_cross_groups") or []
    if ctl_ok and cross:
        return "controlled_cross_kingdom", (
            f"{row['control_loci']} control locus/loci found, and baits from "
            f"{', '.join(cross)} align across them as well — the search "
            "crosses kingdom-level divergence in this assembly"
            + (f"; {row['ryr_loci']} RyR locus/loci" if ryr_ok else ""))
    if ctl_ok:
        return "controlled", (
            f"{row['control_loci']} control locus/loci found, but only from baits "
            "in this genome's own group — the assembly is searchable; how far "
            "the search reaches is not shown"
            + (f"; {row['ryr_loci']} RyR locus/loci" if ryr_ok else ""))
    if expects_ryr and ryr_ok:
        return "controlled_ryr_only", (
            f"no control locus but {row['ryr_loci']} RyR locus/loci — the "
            "search "
            "reached the assembly")
    if not has_control_bait and not expects_ryr:
        return "no_control_bait", (
            "the panel carries no control bait for this genome's clade "
            "(controls are built for the clades carrying a negative claim), "
            "so no control could have fired here")
    return "uncontrolled", (
        "neither the proof-of-search control nor "
        + ("the RyR control" if expects_ryr else "an expected RyR")
        + " was found; no absence claim may rest on this genome")
