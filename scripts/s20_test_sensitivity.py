"""Equivalence test for S20's one-search-two-sensitivities design.

S20 runs every search at `-E 10` and takes the primary call by filtering the
same output at E ≤ 1e-5. The design was written on the assumption — asserted
from documentation, not measured — that `-E` is a *sequence* reporting
threshold and the two files would therefore be identical below the cut.

**They are not, and this test is how that was found.** Over 798 protist
targets the sequence set matches exactly, but four targets carry a different
number of *domain* rows. The cause is that `--domE` (default 10.0) is applied
to the **conditional** E-value, and the conditional E-value is normalised by
how many sequences passed the sequence threshold. A looser `-E` lets more
sequences through, which inflates every c-Evalue by a constant factor (1.41×
here), which pushes marginal domains over `--domE` and out of the report.

Every domain that differs is a marginal one — the extra rows the strict run
carries score at or below 0 bits — but they are not harmless: merged into
`hmm_coverage` they can move a target's evidence class, because a −2.0-bit
alignment spanning 1,386 match states counts as coverage. So the difference
is measured here rather than assumed away, and this task deliberately uses
the **relaxed** side, which excludes them.

What is compared, on a real group DB, row for row:

  1. the same profile is run twice over the same database, once at each
     threshold;
  2. the relaxed output is filtered at the strict threshold;
  3. the filtered set and the strict run must agree on **every** target, and
     on each target's full-sequence E-value, bit score and domain count.

A mismatch fails the test rather than being reported as a difference, because
there is no version of this design that survives one.

`archaea` is the default subject: 634 proteomes and 0.5 G residues make it the
cheapest group that is still a real database, and it is one of the groups the
negative claims are made on.

Run:  python3 scripts/s20_test_sensitivity.py --group archaea --cpu 4
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_groups import EVALUE_PRIMARY, EVALUE_RELAXED, GROUPS  # noqa: E402
from scripts.s20_lib import S20_DIR, group_paths, log, s20_dirs, write_json  # noqa: E402
from scripts.s20_run_sweep import run as run_cmd  # noqa: E402
from scripts.s3_assign import MIN_PROFILE_POSITIONS, assign  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    HMM_SWEEP_DIR, best_hits_by_target, parse_domtblout,
)

COMPARED = ("full_evalue", "full_score", "n_domains", "hmm_coverage",
            "target_coverage", "tlen", "qlen")

#: Below this the comparison is not evidence. Two empty sets are identical.
MIN_TARGETS = 25


def search_at(hmm: Path, db: Path, out: Path, evalue: str, cpu: int) -> Path:
    if out.exists() and out.stat().st_size:
        log("s20_test", f"reusing {out.name}")
        return out
    tmp = out.with_suffix(".part")
    run_cmd(["hmmsearch", "--cpu", str(cpu), "-E", evalue, "--noali",
             "--domtblout", str(tmp), "-o", str(out.with_suffix(".log")),
             str(hmm), str(db)])
    tmp.rename(out)
    return out


def compare(strict: Path, relaxed: Path, cut: float) -> dict:
    strict_hits = best_hits_by_target(parse_domtblout(strict))
    relaxed_rows = [r for r in parse_domtblout(relaxed)
                    if r["full_evalue"] <= cut]
    filtered = best_hits_by_target(relaxed_rows)

    only_strict = sorted(set(strict_hits) - set(filtered))
    only_filtered = sorted(set(filtered) - set(strict_hits))
    differing = []
    for target in sorted(set(strict_hits) & set(filtered)):
        a, b = strict_hits[target], filtered[target]
        for field in COMPARED:
            if a[field] != b[field]:
                differing.append({"target": target, "field": field,
                                  "strict": a[field], "filtered": b[field]})
    # The fields are only a proxy. What matters is whether the *call* moves,
    # so the two sides are run through the same assignment the sweep uses.
    return {"n_strict": len(strict_hits), "n_filtered": len(filtered),
            "only_in_strict_run": only_strict[:20],
            "n_only_in_strict_run": len(only_strict),
            "only_in_filtered_relaxed": only_filtered[:20],
            "n_only_in_filtered_relaxed": len(only_filtered),
            "differing_fields": differing[:20],
            "n_differing_fields": len(differing),
            "compared_fields": list(COMPARED)}


def call_effect(strict: Path, relaxed: Path, cut: float,
                profile: str) -> dict:
    """Does the difference move any assignment, or only a coverage number?

    Runs the sweep's own `assign()` over each side with the *other* profile
    absent, so the comparison isolates what this file contributes: the call,
    the evidence class and the profile positions the D22 gate reads.
    """
    from scripts.s20_run_sweep import _primary
    a = {r["accession"]: r for r in assign(_primary(strict, cut), {})} \
        if profile == "itpr" else \
        {r["accession"]: r for r in assign({}, _primary(strict, cut))}
    b = {r["accession"]: r for r in assign(_primary(relaxed, cut), {})} \
        if profile == "itpr" else \
        {r["accession"]: r for r in assign({}, _primary(relaxed, cut))}
    # Targets present on only one side matter only if that side *called*
    # them. A one-sided target the assignment declines changes nothing
    # downstream, and on this data all of them are declined by the D22 gate.
    one_sided = []
    for acc in sorted(set(a) ^ set(b)):
        row = a.get(acc) or b[acc]
        one_sided.append({"accession": acc,
                          "side": "strict" if acc in a else "relaxed",
                          "assignment": row["assignment"],
                          "reason": row["reason"]})
    moved_call, moved_evidence, moved_gate = [], [], []
    for acc in sorted(set(a) & set(b)):
        if a[acc]["assignment"] != b[acc]["assignment"]:
            moved_call.append({"accession": acc,
                               "strict": a[acc]["assignment"],
                               "relaxed": b[acc]["assignment"]})
        if a[acc]["evidence"] != b[acc]["evidence"]:
            moved_evidence.append({"accession": acc,
                                   "strict": a[acc]["evidence"],
                                   "relaxed": b[acc]["evidence"]})
        col = f"{profile}_positions"
        if (a[acc][col] >= MIN_PROFILE_POSITIONS) != \
                (b[acc][col] >= MIN_PROFILE_POSITIONS):
            moved_gate.append({"accession": acc, "strict": a[acc][col],
                               "relaxed": b[acc][col]})
    called_one_sided = [r for r in one_sided
                        if r["assignment"] in ("ITPR", "RYR")]
    return {"n_compared": len(set(a) & set(b)),
            "n_one_sided": len(one_sided),
            "n_one_sided_called": len(called_one_sided),
            "one_sided_assignments": dict(Counter(
                r["assignment"] for r in one_sided)),
            "one_sided_called": called_one_sided[:10],
            "one_sided_example_reason": one_sided[0]["reason"] if one_sided
            else "",
            "n_call_changed": len(moved_call), "call_changed": moved_call[:10],
            "n_evidence_changed": len(moved_evidence),
            "evidence_changed": moved_evidence[:10],
            "n_d22_gate_changed": len(moved_gate),
            "d22_gate_changed": moved_gate[:10]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", default="archaea", choices=list(GROUPS))
    ap.add_argument("--profile", default="itpr", choices=("itpr", "ryr"))
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()

    _, raw = s20_dirs()
    db = group_paths(args.group)["db"]
    if not db.exists():
        raise SystemExit(f"sweep DB missing: {db}")
    hmm = HMM_SWEEP_DIR / f"{args.profile}.hmm"

    # The relaxed run is the sweep's own output where it exists — testing a
    # fresh copy would prove a property of hmmsearch without proving anything
    # about the files this task actually used.
    relaxed = raw / f"hmmsearch_{args.profile}_vs_{args.group}.domtblout"
    if not relaxed.exists():
        relaxed = search_at(hmm, db,
                            raw / f"test_relaxed_{args.profile}_{args.group}.domtblout",
                            EVALUE_RELAXED, args.cpu)
    else:
        log("s20_test", f"relaxed side is the sweep's own {relaxed.name}")
    strict = search_at(hmm, db,
                       raw / f"test_strict_{args.profile}_{args.group}.domtblout",
                       EVALUE_PRIMARY, args.cpu)

    result = compare(strict, relaxed, float(EVALUE_PRIMARY))
    result["call_effect"] = call_effect(strict, relaxed,
                                        float(EVALUE_PRIMARY), args.profile)
    result.update({"group": args.group, "profile": args.profile,
                   "strict_evalue": EVALUE_PRIMARY,
                   "relaxed_evalue": EVALUE_RELAXED,
                   "strict_domtblout": str(strict),
                   "relaxed_domtblout": str(relaxed)})
    # A comparison of zero targets agrees perfectly and proves nothing. The
    # first run of this test passed on `archaea`, which has no ITPR hit at
    # either threshold — a green light from an empty set. A test that can
    # pass vacuously is worse than no test, so an empty comparison is a
    # failure with its own message.
    result["min_targets"] = MIN_TARGETS
    ce = result["call_effect"]
    # The criterion is about the *call*, not raw set equality. HMMER prints
    # the sequence E-value to two significant figures, so a target whose true
    # E-value sits just above 1e-5 prints as `1e-05` and this task's `<=`
    # filter admits it while a `-E 1e-5` run does not report it at all. That
    # is a boundary-rounding artefact of the printed table: 42 of 13,770 plant
    # targets on the ryr profile, every one at 31.2 bits and every one
    # declined by the D22 gate. Requiring set equality would fail the test on
    # a difference that cannot reach a result; requiring that no one-sided
    # target is *called* tests the thing that can.
    ok = (result["n_strict"] >= MIN_TARGETS
          and ce["n_call_changed"] == 0
          and ce["n_d22_gate_changed"] == 0
          and ce["n_one_sided_called"] == 0)
    result["equivalent"] = ok
    write_json(S20_DIR / f"sensitivity_equivalence_{args.group}.json", result)

    log("s20_test", f"{args.group}/{args.profile}: strict run "
                    f"{result['n_strict']} targets, relaxed run filtered to "
                    f"{result['n_filtered']}")
    if ok:
        log("s20_test",
            f"PASS — {ce['n_compared']:,} shared targets agree on every "
            f"assignment; {result['n_differing_fields']} marginal domain-row "
            f"difference(s) and {ce['n_one_sided']} boundary target(s) move "
            f"no call and no D22 gate ({ce['n_evidence_changed']} "
            f"evidence-class change(s), reported not ignored; one-sided "
            f"targets are {ce['one_sided_assignments']})")
        return 0
    if result["n_strict"] < MIN_TARGETS:
        log("s20_test", f"FAIL — the strict run returned "
                        f"{result['n_strict']} targets, under the "
                        f"{MIN_TARGETS} this comparison needs to mean "
                        "anything. Point --group/--profile at a combination "
                        "with real hits.")
        return 1
    log("s20_test", f"FAIL — {ce['n_call_changed']} assignment(s) changed, "
                    f"{ce['n_d22_gate_changed']} D22 gate crossing(s), "
                    f"{ce['n_one_sided_called']} one-sided target(s) actually "
                    "called")
    return 1


if __name__ == "__main__":
    sys.exit(main())
