"""Constructed negative controls for the S16 rules, run before anything is
written (the pattern of `s5_bait_screen.self_test()`, `s15_test_loss.py`).

Every rule in a duplication analysis returns a plausible number when it is
wrong, and two of S16's failure modes are silent:

* a **split model counted as a duplicate** manufactures the very thing the
  task is testing — and the merge rule fired on 0 of 2,146 loci in this
  sweep, so a merge that could not fire at all would look identical;
* a **circular paralogon test** — leaving the family's own genes in the
  windows, or building the permutation null from a different link set than
  the test — produces a beautiful p-value out of the thing being explained.

So the tests here are mostly tests on *refusal* and on *reachability*: the
rule must fire on a constructed case, and must not fire on its near-miss.

Run: `python3 scripts/s16_test_dup.py` (exit 0 = all pass).
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_paralogon as P                                      # noqa: E402
import s16_teleost as T                                        # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + (f" — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def _loc(acc="A", cell="ITPR1", contig="c1", start=1000, end=2000,
         q=(1, 2700), cov=1.0, ident=0.9, aa=2700, strand="+", lesions=0,
         bait="X|ITPR1|Sp"):
    return {
        "accession": acc, "organism": "Test sp", "vclass": "Actinopteri",
        "vorder": "Testiformes", "assembly_level": "Chromosome",
        "contig_n50": 50_000_000, "contig_spans_gene": True, "cell": cell,
        "cell_status": "found_annotated", "slot": "cell", "locus_idx": 0,
        "locus_clade": cell, "locus_family": "ITPR", "assigned_via": "primary",
        "contig": contig, "start": start, "end": end, "strand": strand,
        "bait": bait, "identity": ident, "coverage": cov, "aligned_aa": aa,
        "q_start": q[0], "q_end": q[1], "bait_len": 2700, "frameshifts": 0,
        "stop_codons": lesions, "lesions": lesions,
        "lesion_density": 1000.0 * lesions / aa if aa else 0.0,
        "family_margin": 1.0, "paralog_margin": 1.0, "contig_edge": False,
        "longest_n_run": 0, "annot_gene": "", "annot_frac_cds": 0.0,
        "annot_paralog_matches": False, "merged_models": 1, "merged_from": "",
    }


# ------------------------------------------------------------ the merge

def t1_split_model_is_one_copy() -> None:
    """T1 — two alignments over overlapping halves of one bait, close
    together on one strand, are one gene. The halves are each built *above*
    the copy bar on purpose: that is the only configuration in which the
    merge changes the count, and a test whose halves fall below the bar
    anyway would pass whether the merge worked or not."""
    a = _loc(start=1_000, end=60_000, q=(1, 1560), cov=0.58, aa=1560)
    b = _loc(start=70_000, end=130_000, q=(1400, 2700), cov=0.58, aa=1301)
    check("T1 both halves clear the copy bar on their own",
          L.is_copy(a) and L.is_copy(b))
    merged, records = L.merge_split_models([a, b])
    check("T1 split model merges to one copy",
          len(merged) == 1 and len(records) == 1,
          f"{len(merged)} loci, {len(records)} merge records")
    check("T1 merged record carries both alignments' spans",
          bool(merged) and merged[0]["q_start"] == 1
          and merged[0]["q_end"] == 2700
          and merged[0]["merged_models"] == 2)
    check("T1 unmerged the same pair would have been two copies",
          sum(1 for r in (a, b) if L.is_copy(r)) == 2
          and sum(1 for r in merged if L.is_copy(r)) == 1)


def t2_tandem_duplicates_do_not_merge() -> None:
    """T2 — the other half of T1. Two alignments each covering the *whole*
    bait are two genes however close they are; a merge rule that folded them
    would delete real duplications."""
    a = _loc(start=1_000, end=60_000, q=(1, 2700))
    b = _loc(start=70_000, end=130_000, q=(1, 2700))
    merged, records = L.merge_split_models([a, b])
    check("T2 two complete models 10 kb apart stay two copies",
          len(merged) == 2 and not records, f"{len(merged)} loci")


def t3_merge_respects_strand_and_contig() -> None:
    """T3 — complementary spans on opposite strands, or on two contigs, are
    not one gene the aligner split."""
    a = _loc(start=1_000, end=60_000, q=(1, 1300), cov=0.48, aa=1300)
    b = dict(_loc(start=70_000, end=130_000, q=(1310, 2700), cov=0.51,
                  aa=1390), strand="-")
    c = dict(_loc(contig="c2", start=1_000, end=60_000, q=(1310, 2700),
                  cov=0.51, aa=1390))
    check("T3 opposite strands do not merge",
          len(L.merge_split_models([a, b])[0]) == 2)
    check("T3 different contigs do not merge",
          len(L.merge_split_models([a, c])[0]) == 2)


def t4_merge_is_order_invariant() -> None:
    """T4 — the merge must not depend on the order the sweep listed the loci
    in, or a copy count becomes a fact about a JSON file."""
    recs = [_loc(start=1_000, end=60_000, q=(1, 1300), cov=0.48, aa=1300),
            _loc(start=70_000, end=130_000, q=(1310, 2700), cov=0.51, aa=1390),
            _loc(contig="c2", start=5_000, end=70_000, q=(1, 2700))]
    a = L.merge_split_models(recs)[0]
    b = L.merge_split_models(list(reversed(recs)))[0]
    key = lambda rs: [(r["contig"], r["start"], r["end"]) for r in rs]  # noqa: E731
    check("T4 merge is order-invariant", key(a) == key(b))


def t5_copy_floors_fire() -> None:
    """T5 — each copy floor rejects its own violation and nothing else."""
    cases = [
        ("identity_below_sweep_floor", _loc(ident=0.20)),
        ("coverage_below_full", _loc(cov=0.30)),
        ("too_short", _loc(cov=1.0, aa=200)),
        ("cell_unassigned", _loc(cell="unassigned")),
        ("", _loc()),
    ]
    for expected, rec in cases:
        got = L.copy_exclusion(rec)
        check(f"T5 copy rule returns {expected or '(a copy)'}",
              got == expected, f"got {got or '(a copy)'}")


# --------------------------------------------------------- the 2R test

def t6_family_genes_leave_their_own_windows() -> None:
    """T6 — the circularity guard. An ITPR or RyR gene left inside a window
    would let the family's own paralogy count as evidence that its blocks
    are paralogous, which is the thing being explained."""
    by = {"c1": [("c1", i * 1000, i * 1000 + 500, "+", f"g{i}",
                  sym, "protein_coding")
                 for i, sym in enumerate(
                     ["AAA", "ITPR1", "BBB", "RYR2", "CCC", "ITPR3P",
                      "DDD", "EEE", "FFF", "GGG"])]}
    win = P.window_genes(by, "c1", 4000, 4500, 10)
    syms = [g[5] for g in win]
    check("T6 no ITPR/RyR symbol survives in a window",
          not any(s.upper().startswith(P.FAMILY_PREFIXES) for s in syms),
          ";".join(syms))
    check("T6 the non-family neighbours do survive",
          {"AAA", "BBB", "CCC", "DDD"} <= set(syms))


def t7_dating_filters_links() -> None:
    """T7 — the duplication-node filter must actually discriminate. It is
    what separates a 2R ohnolog pair from an older duplication whose copies
    happen to sit in these blocks, and it changed this task's answer."""
    pairs = {"a": [("x", "Vertebrata"), ("y", "Opisthokonta")]}
    allp = P.links_between(["a"], {"x", "y"}, pairs)
    dated = P.links_between(["a"], {"x", "y"}, pairs,
                            set(P.LEVELS_2R_WINDOW))
    check("T7 undated test sees both links", len(allp) == 2)
    check("T7 dated test sees only the vertebrate one",
          len(dated) == 1 and dated[0][1] == "x")
    check("T7 the null is built from the same link set",
          P.paralogs_of(["a"], pairs, set(P.LEVELS_2R_WINDOW)) == {"x"}
          and P.paralogs_of(["a"], pairs, None) == {"x", "y"})


def t8_null_uses_real_windows() -> None:
    """T8 — D17. The null must be drawn from the real gene order, and a
    window's paralogs must be findable in it, or the test cannot reject."""
    rng = random.Random(1)
    windows = [["p1", "z1", "z2"], ["z3", "z4", "z5"]] * 50
    null = P.permutation_null({"p1"}, 3, windows, rng, n=2000)
    check("T8 null is non-degenerate and bounded",
          0 < sum(null) / len(null) < 1.0,
          f"mean {sum(null) / len(null):.3f}")
    empty = P.permutation_null(set(), 3, windows, rng, n=200)
    check("T8 a window with no paralogs scores 0", set(empty) == {0})


def t9_family_root_collapses_arrays() -> None:
    """T9 — a tandem array is one duplication, not a paralogon."""
    fams = {P.family_root(s) for s in
            ("ZNF763", "ZNF799", "ZNF433", "SLC22A10")}
    check("T9 nine zinc fingers are one family",
          fams == {"ZNF", "SLC"}, ";".join(sorted(fams)))


def t10_bh_is_monotone_and_bounded() -> None:
    q = L.benjamini_hochberg([0.001, 0.02, 0.5, 1.0])
    check("T10 BH is monotone, bounded and >= p",
          q == sorted(q) and all(0 <= x <= 1 for x in q)
          and all(a >= b for a, b in zip(q, [0.001, 0.02, 0.5, 1.0])),
          f"{[round(x, 4) for x in q]}")
    check("T10 BH of an empty family is empty",
          L.benjamini_hochberg([]) == [])


# --------------------------------------------------------- the 3R test

def t11_groups_are_assigned_by_rule() -> None:
    cases = [
        ({"vclass": "Actinopteri", "vorder": "Semionotiformes",
          "organism": "Lepisosteus oculatus"}, "pre_3R_outgroup"),
        ({"vclass": "Cladistia", "vorder": "Polypteriformes",
          "organism": "Erpetoichthys calabaricus"}, "pre_3R_outgroup"),
        ({"vclass": "Actinopteri", "vorder": "Salmoniformes",
          "organism": "Salmo salar"}, "extra_wgd"),
        ({"vclass": "Actinopteri", "vorder": "Cypriniformes",
          "organism": "Cyprinus carpio"}, "extra_wgd"),
        ({"vclass": "Actinopteri", "vorder": "Cypriniformes",
          "organism": "Danio rerio"}, "teleost"),
        ({"vclass": "Mammalia", "vorder": "Primates",
          "organism": "Homo sapiens"}, "non_actinopterygian"),
    ]
    for row, expected in cases:
        got = T.group_of(row)
        check(f"T11 {row['organism']} -> {expected}", got == expected,
              f"got {got}")


def t12_cross_anchor_agreement_can_be_chance() -> None:
    """T12 — the reachability half. Under independent lineage-specific
    duplications the anchors carry no shared information and agreement must
    land near 0.5; a statistic that returns 1.0 on random blocks would make
    the 3R result unfalsifiable."""
    rng = random.Random(7)
    vocab = [f"S{i}" for i in range(400)]
    sets = []
    for g in range(30):
        for c in range(2):
            keys = frozenset(rng.sample(vocab, 12))
            sets.append({"accession": f"G{g}", "organism": f"Sp {g}",
                         "vclass": "Actinopteri", "vorder": f"Ord{g}",
                         "cell": T.FOCAL, "group": "teleost",
                         "contig": f"c{c}", "start": 1, "bait": "b|X|Y",
                         "keys": keys, "label": f"G{g}#{c}"})
    r3, _ = T.block_consistency(sets, log=lambda *a: None)
    frac = float(dict((r[0], r[1]) for r in r3[1:])["fraction_agreeing"])
    check("T12 random blocks agree at about chance", 0.4 <= frac <= 0.75,
          f"fraction_agreeing {frac}")


def t13_cross_anchor_agreement_finds_a_real_partition() -> None:
    """T13 — the positive half of T12. Two genuinely shared blocks must be
    recovered even when each genome's own copy order is scrambled."""
    rng = random.Random(11)
    block_a = [f"A{i}" for i in range(20)]
    block_b = [f"B{i}" for i in range(20)]
    sets = []
    for g in range(30):
        ka = frozenset(rng.sample(block_a, 10))
        kb = frozenset(rng.sample(block_b, 10))
        pair = [ka, kb] if g % 2 == 0 else [kb, ka]
        for c, keys in enumerate(pair):
            sets.append({"accession": f"G{g}", "organism": f"Sp {g}",
                         "vclass": "Actinopteri", "vorder": f"Ord{g}",
                         "cell": T.FOCAL, "group": "teleost",
                         "contig": f"c{c}", "start": 1, "bait": "b|X|Y",
                         "keys": keys, "label": f"G{g}#{c}"})
    r3, _ = T.block_consistency(sets, log=lambda *a: None)
    d = dict((r[0], r[1]) for r in r3[1:])
    check("T13 a real shared partition is recovered",
          float(d["fraction_agreeing"]) >= 0.95
          and float(d["p_binomial_two_sided"]) < 1e-6,
          f"{d['fraction_agreeing']} p={d['p_binomial_two_sided']}")


def t14_binomial_is_exact() -> None:
    check("T14 binomial tail is exact",
          abs(T.binom_p(10, 10) - 1 / 1024) < 1e-12
          and abs(T.binom_p(0, 10) - 1.0) < 1e-12)


def main() -> int:
    print("S16 negative controls")
    for fn in (t1_split_model_is_one_copy, t2_tandem_duplicates_do_not_merge,
               t3_merge_respects_strand_and_contig, t4_merge_is_order_invariant,
               t5_copy_floors_fire, t6_family_genes_leave_their_own_windows,
               t7_dating_filters_links, t8_null_uses_real_windows,
               t9_family_root_collapses_arrays, t10_bh_is_monotone_and_bounded,
               t11_groups_are_assigned_by_rule,
               t12_cross_anchor_agreement_can_be_chance,
               t13_cross_anchor_agreement_finds_a_real_partition,
               t14_binomial_is_exact):
        fn()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        return 1
    print("\nall S16 negative controls pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
