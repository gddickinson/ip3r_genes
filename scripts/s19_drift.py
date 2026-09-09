"""S19 — iterative-search drift, recoded from the raw logs.

jackhmmer builds each round's profile from the previous round's inclusion
list, so one wrong inclusion is enough to start a walk out of the family.
This project made seven such runs and D10's kill criterion was written to
catch the walk. What S19 adds is not a fourth opinion but a **measurement of
the criterion itself**: each rule re-evaluated per round on that round's own
included list, so the reader can see which rule could have fired, when, and
what the run's model was made of by then.

The result the sessions kept meeting is D10b, and it is general rather than
about this family. Off-family accretion *dilutes* the sister-family share, so
K1 — the rule written to catch drift — moves the **wrong way** while a run
drifts. `s20_protista_other` grew from 721 to 26,148 targets while its sister
share fell from 0.28 % to 0.16 %, and K2 never fired either because no single
round grew by more than 3.18x. Only the round ceiling stopped it, and a
ceiling is a budget, not a diagnosis.

So this module measures a **replacement** signal on the same logs: the
off-family share and its rise, which is what actually moves when a run
drifts. It is reported beside K1's own trace, not substituted for it, and it
is offered as a proposal rather than applied — changing the criterion would
change a committed verdict, which is a decision and not a methods result.

Nothing is re-run. Everything is parsed from the archived logs and the
committed convergence tables.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import s3_kill as K                                             # noqa: E402
import s19_lib as S                                             # noqa: E402

#: The single-pass baseline each run should be read against — the same
#: database searched once, without iteration.
def hmmsearch_baseline(group: str) -> tuple[int, int]:
    """(targets scored, targets called family) by one hmmsearch pass."""
    path = (S.HMM_DIR / "hmmsearch_assignments.tsv" if group == "vertebrata"
            else S.S20_DIR / f"assignments_{group}.tsv")
    if not path.exists():
        return (0, 0)
    rows = S.read_tsv(path)
    return (len(rows), sum(1 for r in rows
                           if r.get("assignment") in ("ITPR", "RYR")))


def _completed(path: Path) -> bool:
    """HMMER writes `[ok]` as the last line of a run that finished itself.

    Kept because it separates a run that stopped because it ran out of rounds
    from one that was interrupted, and the two are not the same evidence.
    """
    try:
        with open(path, "rb") as fh:
            fh.seek(max(0, path.stat().st_size - 400))
            return b"[ok]" in fh.read()
    except OSError:
        return False


def runs_table(log=S.log) -> tuple[list[list], list[list], dict]:
    verdicts = S.jackhmmer_verdicts()
    rows, round_rows = [], []
    summary: dict = {}
    for tag, group, _path in S.JACK_RUNS:
        conv = S.convergence_rounds(tag)
        if not conv:
            continue
        v = verdicts.get(tag, {})
        scored, called = hmmsearch_baseline(group)
        counts = [int(S.fnum(r["new_targets"], float, 0)) for r in conv]
        included = [int(S.fnum(r["n_included"], float, 0)) for r in conv]
        acc_rounds = int(S.fnum(v.get("accepted_rounds"), float, len(conv)))
        final = included[-1] if included else 0
        accepted = int(S.fnum(v.get("accepted_targets"), float, 0))
        stats = _run_stats(tag)
        # A run whose session did not commit a wall clock says so rather
        # than reading zero: the four S20 runs carry no timing at all, and a
        # 0.0 h in a cost table is a claim that they were free.
        secs = S.fnum(stats.get("elapsed_s"), float, 0.0)
        rows.append([tag, v.get("accession", ""), group, scored, called,
                     len(conv), v.get("verdict", ""), v.get("rule") or "-",
                     acc_rounds, accepted, final,
                     round(secs / 3600, 2) if secs else "",
                     "measured" if secs else "not recorded by the session",
                     int(_completed(S.jackhmmer_log(tag))),
                     round(final / max(1, called), 2),
                     round(accepted / max(1, called), 2),
                     v.get("reason", "")])
        for r in conv:
            round_rows.append([
                tag, group, int(S.fnum(r["round"], float, 0)),
                int(S.fnum(r["new_targets"], float, 0)),
                int(S.fnum(r["n_included"], float, 0)),
                int(S.fnum(r["n_own"], float, 0)),
                int(S.fnum(r["n_sister"], float, 0)),
                int(S.fnum(r["n_offfamily"], float, 0)),
                S.fnum(r["sister_frac"], float, 0.0),
                S.fnum(r["sister_rise"], float, 0.0),
                S.fnum(r["growth"], float, 0.0) if r["growth"] else "",
                _offfamily_frac(r), _offfamily_rise(conv, r),
                int(int(S.fnum(r["round"], float, 0)) <= acc_rounds)])
        summary[tag] = {
            "group": group, "verdict": v.get("verdict", ""),
            "rule": v.get("rule", ""), "rounds": len(conv),
            "accepted_rounds": acc_rounds, "accepted_targets": accepted,
            "final_targets": final, "new_targets": counts,
            "hmmsearch_called": called,
        }
    S.write_tsv(S.out_dir() / "jackhmmer_runs.tsv",
                ["run", "seed", "database", "hmmsearch_targets_scored",
                 "hmmsearch_called_family", "rounds_run", "verdict",
                 "kill_rule", "accepted_rounds", "accepted_targets",
                 "final_targets", "elapsed_h", "timing_source",
                 "log_completed",
                 "final_over_hmmsearch_called",
                 "accepted_over_hmmsearch_called", "reason"], rows)
    S.write_tsv(S.out_dir() / "jackhmmer_rounds.tsv",
                ["run", "database", "round", "new_targets", "n_included",
                 "n_own", "n_sister", "n_offfamily", "sister_frac",
                 "sister_rise", "growth", "offfamily_frac", "offfamily_rise",
                 "accepted"], round_rows)
    n_clean = sum(1 for r in rows if r[6] == "clean")
    log(f"jackhmmer_runs: {len(rows)} runs, {n_clean} clean, "
        f"{len(rows) - n_clean} killed")
    return rows, round_rows, summary


def _run_stats(tag: str) -> dict:
    for path in (S.HMM_DIR / f"sweep_stats_jackhmmer_{tag}.json",
                 S.S20_DIR / "sweep_stats.json"):
        blob = S.read_json(path)
        runs = blob.get("jackhmmer") or {}
        if tag in runs:
            return runs[tag]
        for key, val in runs.items():
            if isinstance(val, dict) and val.get("tag") == tag:
                return val
    return {}


def _offfamily_frac(r: dict) -> float:
    n = S.fnum(r["n_included"], float, 0.0)
    return round(S.fnum(r["n_offfamily"], float, 0.0) / n, 4) if n else 0.0


def _offfamily_rise(conv: list[dict], r: dict) -> float:
    idx = conv.index(r)
    if idx == 0:
        return 0.0
    return round(_offfamily_frac(r) - _offfamily_frac(conv[idx - 1]), 4)


# ---------------------------------------------------- the criterion, measured

#: The proposed replacement for K1, stated here so the report can quote a
#: threshold rather than describe one. Not applied — see the module docstring.
MAX_OFFFAMILY_RISE = K.MAX_SISTER_RISE          # the same 0.10, other axis


def criterion_trace(round_rows: list[list], log=S.log,
                    write: bool = True) -> list[list]:
    """When each rule could first have fired, per run.

    K1 as written (sister-family share rising by more than 0.10 in a round),
    K2 (a round including more than 10x the previous one), K3 (the round
    ceiling), and the proposed off-family counterpart of K1 on the same
    rounds. Every rule is evaluated on every run, including the runs it did
    not fire on, because a criterion that is only shown where it worked is
    not a measurement of a criterion.
    """
    by_run: dict[str, list[list]] = {}
    for r in round_rows:
        by_run.setdefault(r[0], []).append(r)
    rows = []
    for tag, rs in by_run.items():
        rs.sort(key=lambda r: r[2])
        group = rs[0][1]
        first = {"K1": None, "K2": None, "K3": None, "OFF": None}
        for r in rs:
            rnd, sister_rise, growth = r[2], r[9], r[10]
            off_rise = r[12]
            if first["K1"] is None and sister_rise > K.MAX_SISTER_RISE:
                first["K1"] = rnd
            if (first["K2"] is None and growth != ""
                    and float(growth) > K.MAX_GROWTH):
                first["K2"] = rnd
            if first["K3"] is None and rnd >= K.MAX_ITER:
                first["K3"] = rnd
            if first["OFF"] is None and off_rise > MAX_OFFFAMILY_RISE:
                first["OFF"] = rnd
        sister = [r[8] for r in rs]
        off = [r[11] for r in rs]
        # "never" rather than an empty cell: a blank here reads as missing
        # data, and the whole point of the table is that K1 never fires.
        rows.append([tag, group, len(rs),
                     first["K1"] or "never", first["K2"] or "never",
                     first["K3"] or "never", first["OFF"] or "never",
                     round(sister[0], 4), round(sister[-1], 4),
                     round(sister[-1] - sister[0], 4),
                     round(off[0], 4), round(off[-1], 4),
                     round(off[-1] - off[0], 4),
                     int(sister[-1] < sister[0])])
    rows.sort(key=lambda r: (r[1], r[0]))
    if not write:
        # The negative controls call this on constructed rounds. Writing
        # then would overwrite the committed table with two synthetic rows —
        # which is exactly what happened on the first build, and the report
        # read 2 runs where there are 7.
        return rows
    S.write_tsv(S.out_dir() / "kill_criterion_trace.tsv",
                ["run", "database", "rounds", "first_round_K1",
                 "first_round_K2", "first_round_K3",
                 "first_round_offfamily_rise", "sister_frac_first",
                 "sister_frac_last", "sister_frac_change",
                 "offfamily_frac_first", "offfamily_frac_last",
                 "offfamily_frac_change", "sister_share_fell"], rows)
    fell = sum(1 for r in rows if r[13])
    log(f"kill_criterion_trace: the sister share *fell* over the run in "
        f"{fell}/{len(rows)} runs — D10b, measured")
    return rows


#: A finished model with more than this share of its included targets scored
#: by neither profile is a model that walked out of the family. The bar is
#: stated rather than fitted, and it does not need fitting: the seven runs
#: separate into 0.23-0.34 and 0.81-0.99 with nothing between.
DRIFTED_OFFFAMILY_FRAC = 0.50


def criterion_validation(trace: list[list], log=S.log) -> list[list]:
    """Score each kill rule as a classifier against the outcome it exists for.

    The outcome is measured on the **finished model**, not taken from any
    session's prose: the share of a run's final included set that neither
    profile scores at all. A run that ends 97 % off-family walked out of the
    family whatever its verdict said; one that ends 23 % off-family did not.
    The seven runs separate cleanly at `DRIFTED_OFFFAMILY_FRAC`, so the
    labels are not a judgement call.

    Then each rule is scored on whether it fires on those runs and not on the
    others. **The proposed off-family rule is a proposal, not a change**: it
    is the same 0.10 threshold K1 already uses, moved to the axis that
    actually moves, and it is validated here on seven runs after the fact.
    Applying it would overturn committed verdicts, which is a decision for a
    session that owns those tables and not a methods result.
    """
    rows = []
    labels = {r[0]: int(r[11] >= DRIFTED_OFFFAMILY_FRAC) for r in trace}
    fired = {
        "K1 (sister-family share rises > 0.10 in a round)":
            {r[0]: int(r[3] != "never") for r in trace},
        "K2 (a round includes > 10x the previous)":
            {r[0]: int(r[4] != "never") for r in trace},
        "K3 (the 10-round ceiling)":
            {r[0]: int(r[5] != "never") for r in trace},
        "proposed: off-family share rises > 0.10 in a round":
            {r[0]: int(r[6] != "never") for r in trace},
    }
    for rule, calls in fired.items():
        tp = sum(1 for k, v in calls.items() if v and labels[k])
        fp = sum(1 for k, v in calls.items() if v and not labels[k])
        fn = sum(1 for k, v in calls.items() if not v and labels[k])
        tn = sum(1 for k, v in calls.items() if not v and not labels[k])
        rows.append([rule, len(calls), tp, fp, fn, tn,
                     round(tp / max(1, tp + fn), 4),
                     round(tn / max(1, tn + fp), 4),
                     ";".join(sorted(k for k, v in calls.items()
                                     if v and not labels[k])) or "-",
                     ";".join(sorted(k for k, v in calls.items()
                                     if not v and labels[k])) or "-"])
    S.write_tsv(S.out_dir() / "kill_criterion_validation.tsv",
                ["rule", "n_runs", "true_positive", "false_positive",
                 "false_negative", "true_negative", "sensitivity",
                 "specificity", "fired_on_a_run_that_did_not_drift",
                 "missed_a_run_that_did"], rows)
    S.write_tsv(S.out_dir() / "drift_outcome.tsv",
                ["run", "database", "offfamily_frac_final", "drifted",
                 "bar"],
                [[r[0], r[1], r[11], labels[r[0]], DRIFTED_OFFFAMILY_FRAC]
                 for r in trace])
    log(f"kill_criterion_validation: {sum(labels.values())}/{len(labels)} "
        f"runs drifted by the finished-model measure")
    return rows


def completeness_effect(summary: dict, log=S.log) -> list[list]:
    """What iteration bought, per database, in family records.

    The one number a completeness argument can rest on: does an iterated
    model hold family records the single pass did not? It is asked of the
    **accepted** rounds only, and separately of the whole run, so a reader
    can see that the extra rounds of a drifting run bought non-family
    targets and nothing else.
    """
    rows = []
    for group in S.ALL_GROUPS:
        _scored, called = hmmsearch_baseline(group)
        runs = [t for t, s in summary.items() if s["group"] == group]
        if not runs:
            rows.append([group, called, 0, "", "", "", "no run"])
            continue
        accepted = S.jackhmmer_targets(group)
        hmm_called = _called_set(group)
        family_in_accepted = accepted & hmm_called
        rows.append([group, called, len(runs), len(accepted),
                     len(family_in_accepted), len(accepted - hmm_called),
                     ";".join(sorted(runs))])
    S.write_tsv(S.out_dir() / "iteration_yield.tsv",
                ["database", "hmmsearch_called_family", "n_runs",
                 "accepted_targets_union", "of_which_called_family",
                 "of_which_not_called_family", "runs"], rows)
    log(f"iteration_yield: {len(rows)} databases")
    return rows


def _called_set(group: str) -> set[str]:
    path = (S.HMM_DIR / "hmmsearch_assignments.tsv" if group == "vertebrata"
            else S.S20_DIR / f"assignments_{group}.tsv")
    return {S.acc_key(r["accession"]) for r in S.read_tsv(path)
            if r.get("assignment") in ("ITPR", "RYR")}


def seed_effect(summary: dict, log=S.log) -> list[list]:
    """Do differently-seeded runs over one database reach the same family set?

    The three vertebrate runs were seeded from a human paralog, a fly gene
    and an amoebozoan gene — as unlike each other as this family allows — and
    S3 reported the family content identical across all three. Recomputed
    here on the accepted rounds, which is the set a completeness claim may
    use, rather than on the final models.
    """
    import s3_hmm_lib as H
    verdicts = S.jackhmmer_verdicts()
    sets: dict[str, set[str]] = {}
    for tag, group, _p in S.JACK_RUNS:
        if group != "vertebrata":
            continue
        path = S.jackhmmer_log(tag)
        if not path.exists():
            continue
        parsed = H.parse_jackhmmer_log(path)
        acc_rounds = int(S.fnum(verdicts.get(tag, {}).get("accepted_rounds"),
                                float, len(parsed["rounds"])))
        sets[tag] = {S.header_accession(">" + a)
                     for a in K.accepted_targets(parsed["rounds"], acc_rounds)}
    called = _called_set("vertebrata")
    rows = []
    for tag, s in sets.items():
        others = set().union(*[v for k, v in sets.items() if k != tag]) \
            if len(sets) > 1 else set()
        fam = s & called
        other_fam = others & called
        rows.append([tag, verdicts.get(tag, {}).get("accession", ""), len(s),
                     len(fam), len(s - called), len(fam - other_fam),
                     len(s - others), len(called - s)])
    if sets:
        shared = set.intersection(*[s & called for s in sets.values()])
        rows.append(["ALL THREE (intersection)", "", "", len(shared), "", "",
                     "", len(called - set.union(*sets.values()))])
    S.write_tsv(S.out_dir() / "seed_effect.tsv",
                ["run", "seed", "accepted_targets", "of_which_family",
                 "of_which_not_family", "family_unique_to_this_seed",
                 "targets_unique_to_this_seed",
                 "hmmsearch_family_not_in_this_run"], rows)
    log(f"seed_effect: {len(sets)} vertebrate runs compared")
    return rows


def run(log=S.log) -> dict:
    rows, round_rows, summary = runs_table(log)
    trace = criterion_trace(round_rows, log)
    validation = criterion_validation(trace, log)
    completeness_effect(summary, log)
    seed_effect(summary, log)
    out = {
        "n_runs": len(rows),
        "clean": sum(1 for r in rows if r[6] == "clean"),
        "killed": sum(1 for r in rows if r[6] == "killed"),
        "by_rule": {rule: sum(1 for r in rows if r[7] == rule)
                    for rule in ("K1", "K2", "K3")},
        "sister_share_fell_in": sum(1 for r in trace if r[13]),
        "k1_never_fired_in": sum(1 for r in trace if r[3] == "never"),
        "offfamily_rule_would_fire_in": sum(1 for r in trace
                                            if r[6] != "never"),
        "drifted_by_finished_model": sum(
            1 for r in trace if r[11] >= DRIFTED_OFFFAMILY_FRAC),
        "rule_scores": {r[0]: {"sensitivity": r[6], "specificity": r[7]}
                        for r in validation},
        "runs": summary,
    }
    S.write_json(S.out_dir() / "drift_summary.json", out)
    return out


if __name__ == "__main__":
    run()
