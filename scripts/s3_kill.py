"""S3 — D10's coded kill criterion for iterative searches.

D10: "Iterative searches need a coded kill criterion. jackhmmer runs that
diverge must be killed by a rule recorded in code and reported, not by eye."

For this family divergence has a specific, measurable shape. An ITPR-seeded
jackhmmer run drifts by walking across the shared domain architecture into
the ryanodine receptors — the two families share every diagnostic Pfam, so
a model that loosens by one round can swallow 6,807 RyRs and report a
beautifully converged, entirely wrong result. That is the failure this
criterion exists to catch, and it is checkable per round because each
round's included set is in the log.

Three rules, evaluated on every round in order. The first that fires ends
the run's *usable* rounds; everything from that round on is reported but
excluded from the census merge.

  K1  Sister-family *drift*. The share of the model's included targets
      assigned to the opposite family rises more than MAX_SISTER_RISE above
      its round-1 value.

      The first version of this rule tested the sister share against a flat
      5 % ceiling, and it fired on every run at round 1 — because a single
      human ITPR1 sequence searched at E <= 1e-5 already returns 32 %
      ryanodine receptors before any iteration has happened. That is not
      contamination; ITPR and RyR are genuine homologues sharing the whole
      pore, and a search sensitive enough to find deep family members is
      necessarily sensitive enough to find its sister family. The *level*
      of sister content is a fact about the two families' shared ancestry,
      so only the *change* in it can be attributed to iterating the model.
      Round 1 is therefore the baseline and K1 watches the rise.

      Measured on the *cumulative* included set — the alignment the next
      round's model is actually built from — not on the round's increment,
      because a model is shaped by everything in it, not by what arrived
      last.
  K2  Explosive growth. The included set grows to more than MAX_GROWTH × its
      previous size, once past a floor that keeps round 2's normal jump from
      tripping it. A model that suddenly rests on ten times as much sequence
      has stopped being a model of this family.
  K3  Ceiling. The run reached MAX_ITER rounds *without* converging, so no
      claim of completeness can rest on it. A run that converges on its
      final allowed round has converged, and K3 does not apply to it.

A run that trips none of them and ends in CONVERGED is what the completeness
argument is allowed to use.
"""

from __future__ import annotations

# K1 — how far the sister-family share may rise above its round-1 baseline,
# in percentage points of the included set. A run that begins with the
# family's natural sister content and holds it has not drifted; one whose
# sister share climbs has. Ten points is a third of the observed baseline,
# so a model would have to take on a substantial new body of RyR sequence
# to trip it.
MAX_SISTER_RISE = 0.10
MIN_ROUND_FOR_K1 = 2     # round 1 *is* the baseline
MAX_GROWTH = 10.0        # K2
MIN_PREV_FOR_K2 = 20     # K2 floor: don't trip on 1 → 15
MAX_ITER = 10            # K3


def evaluate(rounds: list[dict], sister_of: dict[str, str],
             own_family: str, converged: bool = False) -> dict:
    """Apply K1–K3 to a jackhmmer run's parsed rounds.

    `rounds` are s3_hmm_lib.parse_jackhmmer_log() dicts, each with an
    `included` list of target names (cumulative, as jackhmmer reports it).
    `sister_of` maps a target name to the family the profiles assigned it ("ITPR" / "RYR" / "unassigned").
    `own_family` is the seed's family, and `converged` is jackhmmer's own
    verdict — a run that converges on its last allowed round has converged,
    and K3 must not punish it for using every round it was given.

    Returns {"verdict", "accepted_rounds", "killed_at", "rule", "reason",
    "per_round": [...]} — per_round carries the numbers each rule was
    evaluated on, so the criterion is reported and not just applied.
    """
    sister = "RYR" if own_family == "ITPR" else "ITPR"
    per_round: list[dict] = []
    killed_at = None
    rule = reason = ""
    prev_included = None
    baseline: float | None = None

    for rd in rounds:
        included = rd.get("included") or []
        n = len(included)
        calls = [sister_of.get(t, "unassigned") for t in included]
        n_sister = sum(1 for c in calls if c == sister)
        n_own = sum(1 for c in calls if c == own_family)
        # Everything the profiles place in neither family. Reported, not
        # used as a kill rule: this bucket also holds real family members
        # whose records are too fragmentary to clear D22's span gate, so
        # counting it as contamination would punish the run for finding
        # exactly what a sensitive iterative search is supposed to find.
        n_off = n - n_sister - n_own
        frac = round(n_sister / n, 4) if n else 0.0
        if baseline is None and n:
            baseline = frac
        rise = round(frac - baseline, 4) if baseline is not None else 0.0
        growth = round(n / prev_included, 2) if prev_included else None
        per_round.append({
            "round": rd["round"], "new_targets": rd["new_targets"],
            "n_included": n, "n_own": n_own, "n_sister": n_sister,
            "n_offfamily": n_off, "sister_frac": frac,
            "sister_rise": rise, "growth": growth,
        })
        if killed_at is None:
            if rd["round"] >= MIN_ROUND_FOR_K1 and rise > MAX_SISTER_RISE:
                killed_at, rule = rd["round"], "K1"
                reason = (f"after round {rd['round']} the {sister} share of "
                          f"the model was {frac:.1%} ({n_sister}/{n}), "
                          f"{rise:+.1%} on its round-1 baseline of "
                          f"{baseline:.1%} — over the "
                          f"{MAX_SISTER_RISE:.0%}-point drift limit")
            elif (prev_included and prev_included >= MIN_PREV_FOR_K2
                  and growth and growth > MAX_GROWTH):
                killed_at, rule = rd["round"], "K2"
                reason = (f"round {rd['round']} grew the included set to "
                          f"{n} targets, {growth}× the previous round's "
                          f"{prev_included} (limit {MAX_GROWTH}×)")
        prev_included = n or prev_included

    if killed_at is None and not converged and len(rounds) >= MAX_ITER:
        killed_at, rule = MAX_ITER, "K3"
        reason = (f"reached the {MAX_ITER}-round ceiling without converging; "
                  "no completeness claim may rest on this run")

    accepted = (killed_at - 1) if killed_at else len(rounds)
    return {
        "verdict": "killed" if killed_at else "clean",
        "accepted_rounds": accepted,
        "killed_at": killed_at,
        "rule": rule,
        "reason": reason or "no kill rule fired",
        "per_round": per_round,
    }


def accepted_targets(rounds: list[dict], accepted_rounds: int) -> set[str]:
    """Target names included as of the last accepted round — the hits a
    killed run is still allowed to contribute. The per-round lists are
    cumulative, so this is a union only to tolerate a round whose table was
    truncated."""
    out: set[str] = set()
    for rd in rounds:
        if rd["round"] > accepted_rounds:
            break
        out.update(rd.get("included") or [])
    return out
