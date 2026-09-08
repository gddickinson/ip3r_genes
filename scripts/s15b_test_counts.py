"""s15b_test_counts.py — constructed negative controls for S15b's rules,
run before anything is written (`s5_bait_screen.self_test()`'s pattern).

These are checks on *refusal* and on *reachability*, because both failure
modes are silent here.  A Dollo routine that cannot find a loss returns
zero, which is the answer this task expects, so the count would look right
while measuring nothing — T1-T4 put losses on known edges and require them
found, merged and bounded correctly.  A recoder that quietly disagrees
with S15a would make the whole sensitivity matrix a comparison against
something that never ran — T5 requires it to reproduce all 927 committed
states exactly.  And a model fitter that fits an invariant character
returns its own starting point dressed as a rate — T9/T10 require the
refusal and require it to be based on a measured monotone profile.

Each test names the rule it is a control for.  Mutation testing is run by
`--mutate`: three deliberate breakages, each of which must be caught.
"""

from __future__ import annotations

import math
import sys

import s15b_coding as coding
import s15b_dollo as dollo
import s15b_fossils as fossils
import s15b_lib as lib
import s15b_mk as mk

#: a small balanced tree, and one with a polytomy, both written out so a
#: reader can check the expected answers by eye
BALANCED = "(((a,b)AB,(c,d)CD)ABCD,((e,f)EF,(g,h)GH)EFGH)root;"
POLYTOMY = "((a,b,c,d)P,(e,f)EF)root;"


def _t(name: str, rule: str, ok: bool, detail: str = "") -> dict:
    return dict(name=name, rule=rule, passed=bool(ok), detail=detail)


def _char(tree_tips, present, absent, unknown=()) -> dict:
    ch = {}
    for t in tree_tips:
        ch[t] = 1 if t in present else 0 if t in absent else None
    for t in unknown:
        ch[t] = None
    return ch


def run(verbose: bool = False) -> list[dict]:
    out: list[dict] = []

    # ---------------------------------------------------------- the parser
    root = lib.parse_newick(BALANCED)
    tips = [t.name for t in lib.tips(root)]
    out.append(_t("T0 newick round-trip", "parse_newick keeps labels",
                  lib.to_newick(root) == BALANCED and
                  sorted(tips) == list("abcdefgh"),
                  f"{len(tips)} tips"))

    unary = lib.parse_newick("((x)Sp1,(y)Sp2)root;")
    out.append(_t("T0b unary nodes kept", "the sweep's one-tip-per-assembly "
                  "design puts a species node above every tip",
                  [n.name for n in lib.preorder(unary)
                   if len(n.children) == 1] == ["Sp1", "Sp2"]))

    # ------------------------------------------------- Dollo reachability
    ch = _char(tips, present=set(tips) - {"a"}, absent={"a"})
    r = dollo.dollo(root, ch)
    out.append(_t("T1 a single loss is found", "the count must be reachable",
                  r["n_max"] == 1 and r["n_min"] == 1 and
                  r["losses"][0]["node"] == "a",
                  f"n_max={r['n_max']} node={r['losses'][0]['node']
                                             if r['losses'] else None}"))

    ch = _char(tips, present=set(tips) - {"a", "b"}, absent={"a", "b"})
    r = dollo.dollo(root, ch)
    out.append(_t("T2 sister losses merge into their parent edge",
                  "a clade that lost the gene is one event, not two",
                  r["n_max"] == 1 and r["losses"][0]["node"] == "AB",
                  f"n_max={r['n_max']} node={r['losses'][0]['node']
                                             if r['losses'] else None}"))

    ch = _char(tips, present={"a", "b", "c", "d"},
               absent={"e", "f", "g", "h"})
    r = dollo.dollo(root, ch)
    out.append(_t("T3 unpinned, a clade-wide absence is ancestral, not a loss",
                  "Dollo gains once at the MRCA of the present tips; this "
                  "is the reading ITPR2/ITPR3 need, because S13 places "
                  "their duplication inside the tree",
                  r["gain_node"] == "ABCD" and r["n_max"] == 0,
                  f"gain={r['gain_node']} n_max={r['n_max']}"))

    rp = dollo.dollo(root, ch, pin_gain_at_root=True)
    out.append(_t("T3b pinned at the root, the same absence is one loss",
                  "the reading the family character needs, licensed by "
                  "S20/S23's eukaryote-wide range",
                  rp["gain_node"] == "root" and rp["n_max"] == 1 and
                  rp["losses"][0]["node"] == "EFGH",
                  f"gain={rp['gain_node']} n_max={rp['n_max']}"))

    poly = lib.parse_newick(POLYTOMY)
    ptips = [t.name for t in lib.tips(poly)]
    ch = _char(ptips, present={"c", "d", "e", "f"}, absent={"a", "b"})
    r = dollo.dollo(poly, ch)
    out.append(_t("T4 a polytomy bounds the count rather than fixing it",
                  "two losses under an unresolved node may be one event",
                  r["n_max"] == 2 and r["n_min"] == 1 and
                  all(x["parent_is_polytomy"] for x in r["losses"]),
                  f"n_max={r['n_max']} n_min={r['n_min']}"))

    ch = _char(ptips, present={"e", "f"}, absent=set(),
               unknown={"a", "b", "c", "d"})
    r = dollo.dollo(poly, ch, pin_gain_at_root=True)
    out.append(_t("T4b an all-unknown subtree supports no loss",
                  "`?` is not absence — asked with the gain pinned, so the "
                  "subtree is inside the reconstruction and could have "
                  "carried one",
                  r["n_max"] == 0, f"n_max={r['n_max']}"))

    # --------------------------------------------- the recoder is S15a's
    rows = lib.read_tsv(lib.MATRIX)
    re_ = coding.recode_all(rows, **coding.BASE_SETTING)
    mism = [(r["accession"], r["cell"], r["state"], d["state"])
            for r, d in zip(rows, re_) if r["state"] != d["state"]]
    out.append(_t("T5 the recoder reproduces S15a exactly at its own setting",
                  "the sensitivity matrix is a comparison against S15a and "
                  "must start from S15a",
                  not mism, f"{len(mism)} mismatches of {len(rows)}"))

    # the ladder must be ordered: presence can only fall as evidence tightens
    counts = []
    for ev in coding.EVIDENCE_NAMES:
        rr = coding.recode_all(rows, ev, True, True)
        counts.append(sum(1 for x in rr
                          if x["state"] in coding.PRESENT_STATES))
    out.append(_t("T6 the evidence ladder is monotone",
                  "a rung named after what it refuses must not admit more",
                  all(counts[i] >= counts[i + 1]
                      for i in range(len(counts) - 1)),
                  " >= ".join(str(c) for c in counts)))

    # R5 and R6 can only ever *add* absences when switched off
    on = sum(1 for x in coding.recode_all(rows, "no_recon", True, True)
             if x["state"] == "absent")
    off = sum(1 for x in coding.recode_all(rows, "no_recon", False, False)
              if x["state"] == "absent")
    out.append(_t("T7 turning D45/D4 off can only manufacture losses",
                  "the rules are protective; removing one cannot rescue a "
                  "cell", off >= on, f"{on} -> {off}"))

    # a constructed cell must still be able to reach `absent`
    probe = dict(ledger_status="absent", control_ok="1", best_coverage="0",
                 recon_coverage="0.0", n_spare_itpr_loci="0",
                 contig_spans_gene="1")
    out.append(_t("T8 the loss state is still reachable after recoding",
                  "S15a T8's property, re-asserted on this module",
                  coding.recode(probe, **coding.BASE_SETTING)["state"]
                  == "absent"))

    # --------------------------------------------------------- the Mk half
    small = lib.parse_newick(BALANCED)
    lib.apply_lengths(small, lib.branch_lengths(small, "unit"))
    inv = _char(tips, present=set(tips), absent=set())
    fits = [mk.fit(small, inv, m) for m in mk.MODELS]
    out.append(_t("T9 an invariant character is refused, not fitted",
                  "a fitted rate on an invariant character is the "
                  "optimiser's starting point",
                  all(f["fitted"] == 0 and f["reason"] for f in fits),
                  ";".join(str(f["fitted"]) for f in fits)))

    prof = mk.profile(small, inv, "irreversible")
    out.append(_t("T10 the refusal rests on a measured profile",
                  "monotone is reported, not asserted",
                  prof["shape"] == "monotone decreasing" and
                  prof["at_boundary"] == 1, prof["shape"]))

    var = _char(tips, present=set(tips) - {"a", "b"}, absent={"a", "b"})
    f = mk.fit(small, var, "irreversible")
    out.append(_t("T11 a varying character *is* fitted",
                  "the refusal must be about the character, not the model",
                  f["fitted"] == 1 and isinstance(f["rate"], float)
                  and f["rate"] > 0, f"rate={f['rate']}"))

    # the likelihood must be a likelihood: monotone in the right direction
    ll_lo = mk.loglik(small, var, "irreversible", (1e-6,))
    ll_hi = mk.loglik(small, var, "irreversible", (f["rate"],))
    out.append(_t("T12 the fitted optimum beats the boundary",
                  "a fit that cannot beat rate zero is not a fit",
                  ll_hi > ll_lo, f"{ll_lo:.3f} -> {ll_hi:.3f}"))

    # pruning against a hand-checkable case: one tip, one branch
    two = lib.parse_newick("(a,b)root;")
    lib.apply_lengths(two, lib.branch_lengths(two, "unit"))
    q = 0.3
    got = mk.loglik(two, {"a": 1, "b": 0}, "irreversible", (q,))
    want = math.log(math.exp(-q) * (1 - math.exp(-q)))
    out.append(_t("T13 pruning matches the closed form on a two-tip tree",
                  "the rescaling must not change the likelihood",
                  abs(got - want) < 1e-9, f"{got:.9f} vs {want:.9f}"))

    # ---------------------------------------------------------- the fossils
    live = fossils.dead_locus(
        dict(verdict="elevated_lesions", coverage="1.0", stop_codons="0"),
        "present_single_locus")
    out.append(_t("T14 a complete, stop-free locus is not a fossil",
                  "S10's one-sided rule: zero stops falsifies decay",
                  live["is_fossil"] == 0, live["reading"]))

    dead = fossils.dead_locus(
        dict(verdict="elevated_lesions", coverage="0.3", stop_codons="7"),
        "absent")
    out.append(_t("T15 a fossil is still findable",
                  "a screen that cannot fire is not a screen",
                  dead["is_fossil"] == 1, dead["reading"]))

    ctrl = fossils.dead_locus(
        dict(verdict="elevated_lesions", coverage="1.0", stop_codons="0"),
        "not_in_matrix")
    out.append(_t("T16 the RyR control cannot be a fossil by absence of a "
                  "state", "the control is not an ITPR cell and has no row "
                  "in the character matrix",
                  ctrl["is_fossil"] == 0 and
                  ctrl["state_reading_available"] == 0, ctrl["reading"]))

    # BH must be a real correction
    qs = lib.bh([0.001, 0.02, 0.5, 0.9])
    out.append(_t("T17 BH is monotone and never shrinks a p-value",
                  "one family, one correction (S9's rule)",
                  all(q >= p for q, p in zip(qs, [0.001, 0.02, 0.5, 0.9]))
                  and all(qs[i] <= qs[i + 1] + 1e-12
                          for i in range(len(qs) - 1)),
                  ",".join(f"{q:.4f}" for q in qs)))

    if verbose:
        for r in out:
            print(f"  [{'ok ' if r['passed'] else 'FAIL'}] {r['name']}"
                  + (f" — {r['detail']}" if r["detail"] else ""))
    return out


def self_test(verbose: bool = False) -> dict:
    res = run(verbose=verbose)
    failed = [r for r in res if not r["passed"]]
    return dict(n_tests=len(res), n_passed=len(res) - len(failed),
                n_failed=len(failed),
                failures=[r["name"] for r in failed],
                status="fail" if failed else "pass",
                mutations_caught=[])


# --------------------------------------------------------- mutation testing

def mutate(verbose: bool = False) -> list[str]:
    """Break three rules deliberately; each must be caught by the suite."""
    caught = []

    orig = dollo.dollo

    def loose(root, char, pin_gain_at_root=False):
        """Count any subtree holding an absent tip, not only a lost one."""
        res = orig(root, char, pin_gain_at_root)
        counts = dollo._subtree_counts(root, char)
        gain = (root if pin_gain_at_root
                else dollo.mrca_of_present(root, char)[0])
        if gain is None:
            return res
        bad = [c for c in gain.children if counts[id(c)][1] > 0]
        res = dict(res)
        res["n_max"] = len(bad)
        res["losses"] = [dict(node=c.name, parent=gain.name, parent_degree=2,
                              n_tips_lost=counts[id(c)][1], n_tips_unknown=0,
                              tips="", siblings_lost_under_parent=1,
                              parent_is_polytomy=0) for c in bad]
        return res
    dollo.dollo = loose
    if [r for r in run() if r["name"].startswith("T1") and not r["passed"]]:
        caught.append("counted any subtree holding an absent tip as a loss "
                      "-> T1 fired")
    dollo.dollo = orig

    orig_levels = coding.EVIDENCE_LEVELS
    orig_names = coding.EVIDENCE_NAMES
    # make the rung named after refusing everything the most permissive one
    coding.EVIDENCE_LEVELS = orig_levels[:-1] + (
        ("full_locus_only", 0.0, ("R1", "R2", "R3")),)
    coding.EVIDENCE_NAMES = tuple(e[0] for e in coding.EVIDENCE_LEVELS)
    if [r for r in run() if r["name"].startswith("T6") and not r["passed"]]:
        caught.append("made the strictest rung admit more than the loosest "
                      "-> T6 fired")
    coding.EVIDENCE_LEVELS = orig_levels
    coding.EVIDENCE_NAMES = orig_names

    orig_inv = mk.is_invariant
    mk.is_invariant = lambda char: False
    if [r for r in run() if r["name"].startswith("T9") and not r["passed"]]:
        caught.append("let the fitter fit an invariant character -> T9 fired")
    mk.is_invariant = orig_inv

    if verbose:
        for c in caught:
            print(f"  caught: {c}")
    return caught


if __name__ == "__main__":
    st = self_test(verbose=True)
    if "--mutate" in sys.argv:
        st["mutations_caught"] = mutate(verbose=True)
    print(f"\n{st['n_passed']}/{st['n_tests']} passed — {st['status']}")
    sys.exit(0 if st["status"] == "pass" else 1)
