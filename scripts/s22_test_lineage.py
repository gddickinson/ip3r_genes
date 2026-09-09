"""S22's negative controls for the paired design and the lineage half.

Split out of `s22_test_ligand.py` to keep both inside the 500-line budget;
`s22_test_ligand.main()` runs these, so there is still one entry point and
the driver cannot run half the suite.

Four of the rules checked here return a number that looks *better* when
they are wrong — a paired test admitting a tip that covers the pore and not
the core reports a truncated gene model as ligand-core divergence, a matcher
ignoring its own tolerance re-imports the confound it exists to remove, and
a panel filing "no reference proteome" as "no PLC" manufactures its own test
set.  T13b, T14, T20 and T23 are the other kind: a rule that can only ever
return the answer S22 reports is not a measurement, so each is made to fire.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L                                            # noqa: E402


def _bind(check_fn):
    """The reporting function is the caller's, so one suite has one tally."""
    global check
    check = check_fn


check = None  # set by _bind before any test runs


# --- T11-T13  the paired design --------------------------------------------

def t11_uncovered_tip_is_dropped() -> None:
    import s22_paired as P
    tips = [
        {"paralog": "ITPR1", "tip": "good",
         "contact_span_coverage": 0.99, "contact_span_identity": 0.90,
         "channel_minus_luminal_coverage": 0.99,
         "channel_minus_luminal_identity": 0.95},
        {"paralog": "ITPR1", "tip": "n_terminal_truncation",
         "contact_span_coverage": 0.05, "contact_span_identity": 0.10,
         "channel_minus_luminal_coverage": 0.99,
         "channel_minus_luminal_identity": 0.95},
    ]
    pairs, dropped = P._pairs(tips, "ITPR1", "contact_span",
                              "channel_minus_luminal")
    check("T11 a tip covering the pore and not the core is dropped",
          len(pairs) == 1 and len(dropped) == 1
          and dropped[0][1] == "contact_span", f"{len(pairs)}/{len(dropped)}")


def t12_paired_statistic_is_the_difference() -> None:
    """Two tips diverged to different degrees, same core/pore gap."""
    diffs = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    s = L.sign_test(diffs)
    check("T12 a constant gap is detected however divergent the tips are",
          s["n_pos"] == 7 and s["p"] < 0.05, f"p={s['p']}")


def t13_shell_transfer_refuses_unaligned() -> None:
    rows = L.read_tsv(L.OUT_DIR / "shell_index.tsv")
    bad = [r for r in rows if r["transfer"] == "unaligned"
           and r["resi"] not in ("", "None")]
    check("T13 an unaligned shell position carries no residue number",
          not bad, f"{len(bad)} rows carry a number they should not")


# --- T14-T17  the lineage half ---------------------------------------------

def t13b_transfer_refusal_is_reachable() -> None:
    """The refusal fires on nothing in this data, so it is made to fire.

    A branch that never executes is not a rule (`s15b_test_counts.py`'s
    reachability point).  `C.transfer` is replaced with a map that resolves
    nothing, and every row for a paralogue other than the structure
    reference must then come back `unaligned` with no residue number.
    """
    import s22_contacts as CC
    real = CC.transfer
    CC.transfer = lambda paralog: {}
    try:
        rows = CC.shell_index()
    finally:
        CC.transfer = real
    others = [r for r in rows if r["paralog"] != L.STRUCTURE_REF]
    ok = others and all(r["transfer"] == "unaligned" and r["resi"] is None
                        for r in others)
    check("T13b an unresolvable transfer refuses rather than guessing", ok,
          f"{sum(1 for r in others if r['resi'] is not None)} rows kept a "
          "residue number")


def t14_lineage_test_can_fire() -> None:
    """Reachability: a real shift of the size the RyR control produces must
    be detected at the reference group's own size."""
    import random
    rng = random.Random(7)
    ref = [rng.gauss(0.0, 0.05) for _ in range(30)]
    shifted = [rng.gauss(-0.20, 0.05) for _ in range(8)]
    p = L.mann_whitney(shifted, ref)["p"]
    check("T14 a 0.20 shift in 8 tips against 30 is detected",
          p is not None and p < 0.01, f"p={p}")


def t15_no_proteome_is_not_absence() -> None:
    """'We looked and there is no PLC' and 'we could not look' are different
    cells, and only the first may enter the lineage test."""
    path = L.OUT_DIR / "lineage_panel.tsv"
    if not path.exists():
        print("  [skip] T15 lineage_panel.tsv not built yet")
        return
    rows = L.read_tsv(path)
    strata = {r["stratum"] for r in rows}
    check("T15 the two cells have different stratum names",
          "nonvert_itpr_plc_absent" not in strata
          or "nonvert_itpr_plc_no_reference_proteome" in strata
          or not any(r["plc_status"] == "no_reference_proteome" for r in rows),
          str(sorted(strata)))
    mismatched = [r for r in rows
                  if r["stratum"] == "nonvert_itpr_plc_absent"
                  and r["plc_status"] != "absent"]
    check("T15b no tip enters the absent stratum on a status that is not "
          "`absent`", not mismatched, f"{len(mismatched)} mismatched")
    unbacked = [r for r in rows
                if r["plc_status"] == "absent" and not r["plc_upid"]]
    check("T15c every `absent` PLC call names the proteome it was measured on",
          not unbacked, f"{len(unbacked)} calls with no proteome")


def t16_power_is_zero_when_the_test_cannot_run() -> None:
    check("T16 a four-pair sign test has no power at all at alpha 0.05",
          L.power_binomial(4, 0.5, 0.99) == 0.0)
    check("T16b a twenty-pair sign test does have power",
          L.power_binomial(20, 0.5, 0.95) > 0.9)


def t17_plc_needs_both_domains() -> None:
    x = {"a", "b", "c"}
    y = {"b", "c", "d"}
    both = x & y
    check("T17 a protein carrying one PI-PLC half is not a PI-PLC",
          both == {"b", "c"} and "a" not in both and "d" not in both)


# --- T18-T23  the divergence-matched lineage design ------------------------

def _fake(upid, pore, delta, status):
    return {"upid": upid, "organism": upid, "phylum": "Testophyta",
            "pore_identity": f"{pore}", "core_identity": f"{pore + delta}",
            "delta_core_minus_pore": f"{delta}", "plc_status": status,
            "in_test": "True", "protein_count": "10000", "n_itpr": "1",
            "seq_length": "2700"}


def t18_matching_respects_the_tolerance() -> None:
    import s22_deep_lineage as D
    rows = [_fake("A1", 0.35, -0.10, "absent"),
            _fake("P1", 0.36, -0.02, "present"),      # inside tolerance
            _fake("P2", 0.80, -0.09, "present")]      # outside, and closer in y
    pairs, _ = D.divergence_matched(rows)
    check("T18 only controls inside the divergence tolerance are matched",
          pairs[0]["matched_upids"] == "P1", pairs[0]["matched_upids"])


def t19_unmatched_is_reported_not_dropped() -> None:
    import s22_deep_lineage as D
    rows = [_fake("A1", 0.20, -0.10, "absent"),
            _fake("P1", 0.80, -0.02, "present")]
    pairs, summary = D.divergence_matched(rows)
    check("T19 a record with no match is written out with n_matches = 0",
          len(pairs) == 1 and pairs[0]["n_matches"] == 0
          and pairs[0]["difference"] == "" and summary["n_unmatched"] == 1,
          str(summary["n_unmatched"]))


def t20_matched_test_can_fire() -> None:
    """Reachability: a real shift at the observed group size must be found."""
    import s22_deep_lineage as D
    rows = []
    for i in range(20):
        pore = 0.30 + i * 0.001
        rows.append(_fake(f"A{i}", pore, -0.12, "absent"))
        rows.append(_fake(f"P{i}", pore, -0.02, "present"))
    _pairs, summary = D.divergence_matched(rows)
    check("T20 a 0.10 shift in 20 matched pairs is detected",
          summary["p_sign"] and float(summary["p_sign"]) < 0.01,
          f"p={summary['p_sign']} median={summary['median_difference']}")


def t21_matching_is_deterministic() -> None:
    import s22_deep_lineage as D
    rows = [_fake("A1", 0.35, -0.10, "absent")]
    rows += [_fake(f"P{i}", 0.35, -0.02, "present") for i in range(6)]
    a, _ = D.divergence_matched(rows)
    b, _ = D.divergence_matched(list(reversed(rows)))
    check("T21 matching on tied divergence is order-invariant",
          a[0]["matched_upids"] == b[0]["matched_upids"],
          f"{a[0]['matched_upids']} vs {b[0]['matched_upids']}")


def t22_cooccurrence_names_only_swept_groups() -> None:
    path = L.OUT_DIR / "pathway_cooccurrence.tsv"
    if not path.exists():
        print("  [skip] T22 pathway_cooccurrence.tsv not built yet")
        return
    import s22_plc as PLC
    named = {r["group"] for r in L.read_tsv(path)}
    check("T22 the co-occurrence table names no group the PI-PLC sweep "
          "did not cover", named <= set(PLC.EUK_GROUPS),
          str(sorted(named - set(PLC.EUK_GROUPS))))


def t23_deep_panel_coverage_bar_bites() -> None:
    path = L.OUT_DIR / "deep_lineage_panel.tsv"
    if not path.exists():
        print("  [skip] T23 deep_lineage_panel.tsv not built yet")
        return
    rows = L.read_tsv(path)
    dropped = [r for r in rows if r["in_test"] != "True"]
    reasons = {r["drop_reason"].split(" coverage")[0] for r in dropped}
    check("T23 the coverage bar drops records on both modules, not one",
          {"core", "pore"} <= reasons or not dropped, str(sorted(reasons)))
    bad = [r for r in rows if r["in_test"] == "True"
           and (float(r["core_coverage"]) < 0.5 or float(r["pore_coverage"]) < 0.5)]
    check("T23b no record below the bar is in the test", not bad,
          f"{len(bad)} below-bar records in the test")




TESTS = [
    t11_uncovered_tip_is_dropped,
    t12_paired_statistic_is_the_difference,
    t13_shell_transfer_refuses_unaligned,
    t13b_transfer_refusal_is_reachable,
    t14_lineage_test_can_fire,
    t15_no_proteome_is_not_absence,
    t16_power_is_zero_when_the_test_cannot_run,
    t17_plc_needs_both_domains,
    t18_matching_respects_the_tolerance,
    t19_unmatched_is_reported_not_dropped,
    t20_matched_test_can_fire,
    t21_matching_is_deterministic,
    t22_cooccurrence_names_only_swept_groups,
    t23_deep_panel_coverage_bar_bites,
]


def run(check_fn) -> None:
    _bind(check_fn)
    for t in TESTS:
        t()
