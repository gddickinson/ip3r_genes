"""The second half of S21's negative controls, split to stay inside the budget.

T10-T20: the shared-intron test, the duplication detector and the annotation
half. The constructed-input helpers and the pass/fail bookkeeping live in
`s21_test_arch`, which this module imports and which calls `run()` here, so the
two halves cannot record a failure differently (the `s3_report.py` /
`s3_report_d10.py` split, applied to a test suite).
"""

from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_annot as AN                                            # noqa: E402
import s21_calibrate_intron as CAL                                # noqa: E402
import s21_introns as I                                           # noqa: E402
import s21_tandem as TD                                           # noqa: E402
from s21_test_arch import _Aln, _Region, check                    # noqa: E402


# --------------------------------------------------------------------------
# T10-T13 — the shared-intron test
# --------------------------------------------------------------------------
def t10_shared_test_reaches_both_extremes() -> None:
    """Identical intron sets must match fully; disjoint sets not at all."""
    cols = set(range(1, 1001))
    pos = {(10 * i, i % 3) for i in range(1, 51)}
    other = {(10 * i + 3, i % 3) for i in range(1, 51)}
    same = I.shared_pair(pos, pos, cols, 0)
    none = I.shared_pair(pos, other, cols, 0)
    check("T10 shared-intron test reaches both extremes",
          same["observed"] == len(pos) and same["p_analytic"] < 1e-6
          and none["observed"] == 0 and none["p_analytic"] == 1.0
          and I._poisson_binomial_tail([0.3] * 60, 0) == 1.0,
          f"same {same['observed']}/{same['p_analytic']:.3g}, "
          f"none {none['observed']}/{none['p_analytic']:.3g}")


def t11_poisson_binomial_matches_brute_force() -> None:
    """The exact tail must equal an enumeration on a case small enough to enumerate."""
    probs = [0.1, 0.25, 0.5, 0.75]
    for k in range(len(probs) + 1):
        brute = 0.0
        for bits in itertools.product((0, 1), repeat=len(probs)):
            if sum(bits) >= k:
                pr = 1.0
                for b, p in zip(bits, probs):
                    pr *= p if b else 1 - p
                brute += pr
        got = I._poisson_binomial_tail(probs, k)
        if abs(got - brute) > 1e-12:
            check("T11 Poisson-binomial tail is exact", False,
                  f"k={k}: {got} vs {brute}")
            return
    check("T11 Poisson-binomial tail is exact", True)


def t12_both_nulls_agree() -> None:
    """The seeded permutation and the analytic null must land in the same place."""
    cols = set(range(1, 2001))
    pos_a = {(7 * i, i % 3) for i in range(1, 61)}
    pos_b = {(7 * i, i % 3) for i in range(1, 31)} | \
            {(7 * i + 2, i % 3) for i in range(31, 61)}
    rng = random.Random(1)
    d = I.shared_pair(pos_a, pos_b, cols, 0, rng)
    close = abs(d["expected_perm"] - d["expected_analytic"]) < 0.5
    check("T12 permutation and analytic nulls agree", close,
          f"perm {d['expected_perm']} vs analytic {d['expected_analytic']}")


def t13_phase_must_match_for_a_shared_intron() -> None:
    """Same column, different phase, is not one ancestral intron."""
    cols = set(range(1, 1001))
    a = {(100, 0), (200, 1)}
    b = {(100, 1), (200, 2)}
    d = I.shared_pair(a, b, cols, 0)
    check("T13 phase must match for a shared intron", d["observed"] == 0,
          f"observed {d['observed']}")


# --------------------------------------------------------------------------
# T14-T16 — the tandem detector
# --------------------------------------------------------------------------
def t14_detector_fires_and_declines() -> None:
    bait = "Q14643|ITPR1|Homo_sapiens|ITPR1"
    one = [_Aln("c1", 1000, 9000, "+", bait, 1, 900)]
    two = one + [_Aln("c1", 3_000_000, 3_009_000, "+", bait, 1, 900)]
    got_one = TD.detect(one)
    got_two = TD.detect(two)
    check("T14 detector fires on a duplicate pair and declines a single gene",
          not got_one and len(got_two) == 1 and got_two[0][0] == "ITPR1",
          f"single {len(got_one)} pairs, duplicate {len(got_two)} pairs")


def t15_within_locus_class_is_reachable() -> None:
    """`within_locus` is 0 over the sweep; the branch must still be able to fire."""
    bait = "Q14643|ITPR1|Homo_sapiens|ITPR1"
    alns = [_Aln("c1", 1000, 9000, "+", bait, 1, 900),
            _Aln("c1", 14000, 22000, "+", bait, 1, 900)]
    got = TD.detect(alns)
    cls = [TD.pair_class(a, b, same) for _c, _b, a, b, _o, same in got]
    check("T15 within_locus duplication is reachable",
          got and cls == ["within_locus"], f"classes {cls}")


def t16_paralogy_does_not_read_as_duplication() -> None:
    """A bait's hit at another paralogue's gene must not pair inside its own cell.

    The confound that made the first version of this detector useless: the ITPR
    paralogues are 61-68 % identical, so every bait aligns at all three genes.
    Here an ITPR1 bait aligns weakly at a locus the ITPR2 baits win, and the
    pair must not be offered to the ITPR1 cell.
    """
    b1 = "Q14643|ITPR1|Homo_sapiens|ITPR1"
    b2 = "Q14571|ITPR2|Homo_sapiens|ITPR2"
    alns = [
        _Aln("c1", 1000, 9000, "+", b1, 1, 900, score=5000.0, clade="ITPR1"),
        _Aln("c2", 1000, 9000, "+", b2, 1, 900, score=5000.0, clade="ITPR2",
             bait_len=1000),
        _Aln("c2", 1000, 9000, "+", b1, 1, 900, score=1200.0, clade="ITPR1"),
    ]
    got = TD.detect(alns)
    itpr1_pairs = [g for g in got if g[0] == "ITPR1"]
    check("T16 a hit at another paralogue's gene is not a duplication",
          not itpr1_pairs, f"{len(itpr1_pairs)} ITPR1 pairs from paralogy")


# --------------------------------------------------------------------------
# T17-T19 — the annotation half
# --------------------------------------------------------------------------
def t17_terminus_classes_each_reachable() -> None:
    exons = [{"start": 1000, "end": 1099}, {"start": 2000, "end": 2099},
             {"start": 3000, "end": 3099}]
    seq = "A" * 4000
    region = _Region(seq, start=1)
    arch = {"accession": "X", "organism": "o", "vclass": "v", "cell": "ITPR1",
            "locus_idx": 0, "state": "fragmentary"}
    models = [{"gene_id": "g1", "symbol": "s", "biotype": "protein_coding",
               "low": 1000, "high": 2050, "n_blocks": 2},
              {"gene_id": "g2", "symbol": "s", "biotype": "protein_coding",
               "low": 3000, "high": 3099, "n_blocks": 1}]
    rows = AN.terminus_rows(arch, region, models, exons, "+")
    verdicts = {(r["gene_id"], r["terminus"]): (r["vs_alignment"],
                                                r["is_gene_terminus"])
                for r in rows}
    ok = (verdicts[("g1", "acceptor_low")][1] == 1
          and verdicts[("g1", "donor_high")][0] == "mid_exon"
          and verdicts[("g2", "acceptor_low")][0] == "at_exon_boundary"
          and verdicts[("g2", "donor_high")][1] == 1)
    check("T17 gene termini excluded, mid-exon and boundary both reachable",
          ok, str(verdicts))


def t18_intron_terminus_is_its_own_class() -> None:
    exons = [{"start": 1000, "end": 1099}, {"start": 2000, "end": 2099}]
    region = _Region("A" * 3000, start=1)
    arch = {"accession": "X", "organism": "o", "vclass": "v", "cell": "ITPR1",
            "locus_idx": 0, "state": "split"}
    models = [{"gene_id": "g1", "symbol": "s", "biotype": "protein_coding",
               "low": 1000, "high": 1500, "n_blocks": 1},
              {"gene_id": "g2", "symbol": "s", "biotype": "protein_coding",
               "low": 1600, "high": 2099, "n_blocks": 1}]
    rows = AN.terminus_rows(arch, region, models, exons, "+")
    classes = {r["vs_alignment"] for r in rows if not r["is_gene_terminus"]}
    v = AN.locus_verdict(arch, rows)
    check("T18 a terminus inside an intron is its own class",
          classes == {"in_intron"}
          and v["verdict"] == "structure_disagreement",
          f"classes {classes}, verdict {v['verdict']}")


def t19_concordance_can_report_disagreement() -> None:
    exons = [{"start": 1000, "end": 1099}, {"start": 2000, "end": 2099},
             {"start": 3000, "end": 3099}]
    arch = {"accession": "X", "organism": "o", "vclass": "v", "cell": "ITPR1",
            "locus_idx": 0, "state": ""}
    agree = AN.concordance_row(arch, [(1000, 1099), (2000, 2099),
                                      (3000, 3099)], exons, "RefSeq")
    off = AN.concordance_row(arch, [(1000, 1099), (2010, 2089),
                                    (3000, 3099)], exons, "RefSeq")
    check("T19 concordance can report disagreement",
          agree["frac_exact"] == 1.0 and off["frac_exact"] < 1.0,
          f"agree {agree['frac_exact']}, offset {off['frac_exact']}")


# --------------------------------------------------------------------------
# T20-T21 — the calibration, and the suite itself
# --------------------------------------------------------------------------
def t20_calibration_refuses_and_can_separate() -> None:
    """It must refuse with no artefact population and place a bar with one."""
    spliceable = [{"gap_bp": 500, "splice_class": "canonical"}] * 5000
    refused = CAL.calibrate(spliceable)
    mixed = spliceable + [{"gap_bp": 2, "splice_class": "non_canonical"}] * 400
    placed = CAL.calibrate(mixed)
    check("T20 intron calibration refuses without a second population",
          refused["separated"] == 0
          and refused["floor_bp"] == refused["smallest_gap_bp"]
          and placed["separated"] == 1 and placed["floor_bp"] > 2,
          f"refused {refused['separated']}/{refused['floor_bp']}, "
          f"placed {placed['separated']}/{placed.get('floor_bp')}")


def run() -> None:
    t10_shared_test_reaches_both_extremes()
    t11_poisson_binomial_matches_brute_force()
    t12_both_nulls_agree()
    t13_phase_must_match_for_a_shared_intron()
    t14_detector_fires_and_declines()
    t15_within_locus_class_is_reachable()
    t16_paralogy_does_not_read_as_duplication()
    t17_terminus_classes_each_reachable()
    t18_intron_terminus_is_its_own_class()
    t19_concordance_can_report_disagreement()
    t20_calibration_refuses_and_can_separate()
