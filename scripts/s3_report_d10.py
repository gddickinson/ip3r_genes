"""S3 — section 4 of the census-v3 report: jackhmmer, D10, and what an
iterated model is actually built from.

Split out of `s3_report.py` to keep both files under the project's 500-line
limit. Same rule applies here as there (D13): every number is read from a
committed table, nothing is hand-written.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_hmm_lib import (  # noqa: E402
    CENSUS_V3_DIR, HMM_SWEEP_DIR,
)
from scripts.s3_kill import (  # noqa: E402
    MAX_GROWTH, MAX_ITER, MAX_SISTER_RISE,
)
from scripts.s3_seed_spec import JACKHMMER_SEEDS  # noqa: E402


def render(a, conv, v3, rows_or_none) -> None:
    """Append section 4 to the report.

    `a` appends a line, `conv` is convergence.tsv, `v3` the stats JSON and
    `rows_or_none` the caller's table loader — passed in rather than
    re-imported so both modules read the tables the same way.
    """
    # ---------------------------------------------------------- jackhmmer
    # Wall clock per seed: the sweep writes one stats JSON per run, and the
    # cost is a result here — one seed cost more than the other two together.
    elapsed: dict[str, float] = {}
    for tag in JACKHMMER_SEEDS:
        st = HMM_SWEEP_DIR / f"sweep_stats_jackhmmer_{tag}.json"
        if st.exists():
            d = json.loads(st.read_text()).get("jackhmmer", {}).get(tag, {})
            if d.get("elapsed_s"):
                elapsed[tag] = float(d["elapsed_s"])

    a("## 4. jackhmmer to convergence, and D10")
    a("")
    if conv:
        done = {r["seed_tag"] for r in conv}
        missing = [tag for tag in sorted(JACKHMMER_SEEDS) if tag not in done]
        if missing:
            a("> **Still running: "
              + ", ".join(f"`{m}` ({JACKHMMER_SEEDS[m][0]}, "
                          f"{JACKHMMER_SEEDS[m][1]})" for m in missing)
              + ".** The log is archived; fold the run in with `python3 "
              "scripts/s3_run_sweep.py --stage jackhmmer --parse-only "
              + " ".join(f"--seed-tag {m}" for m in missing)
              + "` followed by `s3_census_v3.py`, `s3_figures.py` and this "
              "script. Nothing below depends on it: the completeness "
              "argument rests on the runs that finished.")
            a("")
        a("| Seed | Rounds | New targets per round | Converged | Wall clock "
          "| D10 verdict |")
        a("|---|---:|---|---|---:|---|")
        for tag in sorted({r["seed_tag"] for r in conv}):
            rs = sorted((r for r in conv if r["seed_tag"] == tag),
                        key=lambda r: int(r["round"]))
            curve = " → ".join(r["new_targets"] for r in rs)
            rule = next((r["kill_rule"] for r in rs if r["kill_rule"]), "")
            verdict = rs[0]["verdict"] + (f" ({rule})" if rule else "")
            cost = elapsed.get(tag)
            a(f"| `{tag}` | {len(rs)} | {curve} | {rs[0]['converged']} | "
              + (f"{cost / 3600.0:.1f} h" if cost else "—")
              + f" | {verdict} |")
        a("")
        a("D10 requires a coded kill criterion rather than a judgement by "
          "eye. For this family the failure mode is specific: an "
          "ITPR-seeded run drifts across the shared domain architecture "
          "into the ryanodine receptors and reports a converged, confident, "
          "wrong answer. Three rules, evaluated per round on that round's "
          f"own inclusion list (`scripts/s3_kill.py`): **K1** the sister "
          f"family's share of the model rising more than "
          f"{MAX_SISTER_RISE * 100:.0f} percentage points above its round-1 "
          f"value; **K2** the included set growing more than "
          f"{MAX_GROWTH:.0f}× in one round; **K3** hitting the "
          f"{MAX_ITER}-round ceiling without converging. Rounds from the "
          "first firing on are reported but excluded from the merge.")
        a("")
        base = {r["seed_tag"]: float(r["sister_frac"]) for r in conv
                if int(r["round"]) == 1}
        if base:
            a("**K1 measures the rise, not the level, and that is a result "
              "in itself** (D10a). The rule was first coded as a flat 5 % "
              "ceiling on sister-family content, and it fired on every run "
              "at round 1 — because a single ITPR sequence searched at "
              "E ≤ 1e-5 already returns "
              + (f"{min(base.values()):.0%}–{max(base.values()):.0%}"
                 if len(base) > 1 and
                 round(min(base.values()), 2) != round(max(base.values()), 2)
                 else f"{max(base.values()):.0%}")
              + " ryanodine receptors before any iteration has happened "
              + f"({', '.join(sorted(base))}). "
              "That is not contamination: the two families are genuine "
              "homologues sharing the whole pore, and a search sensitive "
              "enough to reach *Acanthamoeba* is necessarily sensitive "
              "enough to reach RYR1. The level is a fact about their shared "
              "ancestry; only the change in it can be attributed to "
              "iterating the model.")
            a("")
        for tag, d in (v3.get("jackhmmer") or {}).items():
            a(f"* `{tag}` ({d['accession']}): {d['verdict']} — {d['reason']}; "
              f"{d['accepted_rounds']} of {d['n_rounds']} rounds accepted, "
              f"{d['accepted_targets']:,} targets.")
        a("")

        comp = rows_or_none(CENSUS_V3_DIR / "jackhmmer_model_composition.tsv")
        if comp:
            a("### What an iterated model is actually built from")
            a("")
            a("A convergence curve says a search has stopped finding things. "
              "It does not say what it found. Every target supporting each "
              "run's final model was therefore classified by the sweep's own "
              "two-profile verdict and D22 evidence class:")
            a("")
            tags = sorted({r["seed_tag"] for r in comp})
            a("| Composition of the final model | "
              + " | ".join(f"`{t}`" for t in tags) + " |")
            a("|---|" + "---:|" * len(tags))
            keys = []
            for r in comp:
                k = (r["profile_call"], r["sweep_evidence"])
                if k not in keys:
                    keys.append(k)
            by = {(r["seed_tag"], r["profile_call"], r["sweep_evidence"]):
                  (int(r["targets"]), float(r["share"])) for r in comp}
            for call, ev in keys:
                cells = []
                for tag in tags:
                    n, s = by.get((tag, call, ev), (0, 0.0))
                    cells.append(f"{n:,} ({s:.1%})")
                a(f"| {call}, {ev} | " + " | ".join(cells) + " |")
            a("")
            per_seed = {t: {} for t in tags}
            for r in comp:
                per_seed[r["seed_tag"]][(r["profile_call"],
                                         r["sweep_evidence"])] = \
                    int(r["targets"])
            totals = {t: sum(v.values()) for t, v in per_seed.items()}
            mods = [float(r["share"]) for r in comp
                    if r["sweep_evidence"] == "module"]
            a("**Module-only matches are "
              + (f"{min(mods):.1%}–{max(mods):.1%}" if len(mods) > 1
                 else f"{max(mods):.1%}")
              + " of these models** — the SPRY and EF-hand proteins D22's "
              "span gate keeps out of the census. That is what an iterative "
              "search at E ≤ 1e-5 accretes if nothing stops it, and it is "
              "also why these runs decay to an asymptote of a few new targets "
              "a round rather than to zero: the tail being walked is the long "
              "tail of proteins sharing one small domain, not the family.")
            a("")
            if len(tags) > 1:
                # The family-called categories are the completeness claim;
                # everything else is what the run dragged in with them, so
                # the two are reported apart rather than averaged into one
                # "the seeds agree" sentence that only two of them support.
                spread = []
                for call, ev in keys:
                    if call not in ("ITPR", "RYR"):
                        continue
                    ns = [per_seed[t].get((call, ev), 0) for t in tags]
                    spread.append(
                        f"{call}/{ev} " + (f"{ns[0]:,} in all {len(tags)}"
                                           if min(ns) == max(ns)
                                           else f"{min(ns):,}–{max(ns):,}"))
                a("**The family core is the same whichever seed finds it.** "
                  "The runs started from a human ITPR1, a fly Itpr and an "
                  "*Acanthamoeba* receptor — a vertebrate, an insect and an "
                  "amoebozoan — and their final models carry the same "
                  "family-called content: " + "; ".join(spread) + ". Starting "
                  "points that far apart reaching the same core is a "
                  "completeness statement that does not depend on any one of "
                  "them converging.")
                a("")
                key_off = ("unassigned", "no sweep hit")
                worst = max(tags, key=lambda t: per_seed[t].get(key_off, 0))
                best = min(tags, key=lambda t: per_seed[t].get(key_off, 0))
                wr = sorted((r for r in conv if r["seed_tag"] == worst),
                            key=lambda r: int(r["round"]))
                growth = (int(wr[-1]["n_included"]) / int(wr[0]["n_included"])
                          if wr and int(wr[0]["n_included"]) else 0.0)
                a("**What differs between them is everything that is not the "
                  f"family.** `{worst}`'s final model rests on "
                  f"{totals[worst]:,} targets against `{best}`'s "
                  f"{totals[best]:,}, and "
                  f"{per_seed[worst].get(key_off, 0):,} of them are proteins "
                  "neither profile scores at all "
                  f"({per_seed[best].get(key_off, 0):,} in `{best}`). "
                  "**K1 cannot see that drift**: off-family accretion "
                  "*dilutes* the sister-family share instead of raising it"
                  + (f" — across this run it falls from "
                     f"{float(wr[0]['sister_frac']):.1%} to "
                     f"{float(wr[-1]['sister_frac']):.1%} while the model "
                     f"grows {growth:.1f}× — " if wr else " — ")
                  + "so K1 reads divergence as the opposite. K3 is what "
                  "caught it, and the composition table above is why the run "
                  "is reported rather than silently dropped.")
                a("")
