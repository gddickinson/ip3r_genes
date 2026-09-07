"""S7 step 1 — the substitution model, chosen in two stages.

The brief's command is `-m MFP`: ModelFinder's exhaustive scan. Measured
on *this* alignment and *this* machine it does not finish. 11 of up to
1,232 protein models in 11 min 12 s at 8 threads — ~59 s a model,
projecting **~20 hours** — and IQ-TREE's own thread benchmark says 4
threads, not 8, so more cores do not buy the run back.

The cost is the free-rate models. LG, LG+I, LG+G4 and LG+I+G4 took ~20 s
between them; LG+R2 … LG+R8 took 90–150 s each. They are also not
optional here: on this alignment LG+R5 beats LG+I+G4 by **682 BIC
units**, so a scan that dropped `+R` to save the time would return a
measurably worse model.

So the selection is greedy, in two stages, and both stages' tables are
committed:

  A  every matrix, rate heterogeneity fixed at +G4 / +I+G4 (the fast
     ones) -> the best exchangeability matrix by BIC
  B  that matrix alone, every rate-heterogeneity model IQ-TREE offers,
     +F and not -> the final model

**The limitation this buys, stated rather than buried**: the matrix that
wins under +G4 need not be the matrix that would win under +R5. Stage
A's full table is committed so the margin it won by can be read off, and
`model_selection.tsv` carries every model either stage scored.

Run:  python3 scripts/s7_model.py
Out:  results/phylogeny/model_matrix.*, model_rate.*, model_selection.tsv,
      model_choice.json
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.s6_lib import write_tsv                          # noqa: E402
from scripts.s7_lib import PHYLO_DIR                          # noqa: E402

INPUT = PHYLO_DIR / "input.fasta"
STAGE_A = PHYLO_DIR / "model_matrix"
STAGE_B = PHYLO_DIR / "model_rate"
CHOICE = PHYLO_DIR / "model_choice.json"
TABLE = PHYLO_DIR / "model_selection.tsv"

THREADS = 4          # IQ-TREE's own benchmark on this alignment (see S7 report)
SEED = 42

#: Stage A's rate models — the two that are cheap enough to run over
#: every matrix. Stage B replaces this with the full list.
STAGE_A_RATES = "G4,I+G4"

#: The measured projection that made the exhaustive scan unusable, kept
#: here so the report renders it rather than repeating a number by hand.
EXHAUSTIVE_PROBE = {
    "models_scored": 11,
    "of_up_to": 1232,
    "seconds": 672,
    "threads": 8,
    "projection_hours": round(1232 * (672 / 11) / 3600, 1),
    "why_freerate_matters": {
        "LG+I+G4_bic": 433732.851,
        "LG+R5_bic": 433050.518,
        "bic_gain": round(433732.851 - 433050.518, 3),
    },
}

ROW = re.compile(r"^\s*(\d+)\s+(\S+)\s+(-?[\d.]+)\s+(\d+)\s+"
                 r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")


def parse_models(log: Path, stage: str) -> list[dict]:
    out = []
    for line in log.read_text().splitlines():
        m = ROW.match(line)
        if not m:
            continue
        out.append({"stage": stage, "rank": int(m.group(1)),
                    "model": m.group(2), "lnL": m.group(3),
                    "df": m.group(4), "AIC": m.group(5),
                    "AICc": m.group(6), "BIC": m.group(7)})
    return out


def best_by_bic(rows: list[dict]) -> dict:
    return min(rows, key=lambda r: float(r["BIC"]))


def reported_best(prefix: Path) -> str:
    m = re.search(r"Best-fit model.*?: (\S+)",
                  Path(str(prefix) + ".iqtree").read_text())
    return m.group(1) if m else ""


def run(prefix: Path, extra: list[str]) -> float:
    cmd = ["iqtree2", "-s", str(INPUT), "--prefix", str(prefix),
           "-m", "MF", "-T", str(THREADS), "-seed", str(SEED)] + extra
    print("running:", " ".join(cmd), flush=True)
    t0 = time.time()
    rc = subprocess.run(cmd).returncode
    dt = time.time() - t0
    if rc:
        raise SystemExit(f"iqtree2 model selection failed (rc={rc})")
    print(f"  done in {dt:.0f} s", flush=True)
    return dt


def matrix_of(model: str) -> str:
    """`LG+F+R5` -> `LG` — the exchangeability matrix, without its
    frequency or rate terms."""
    return model.split("+")[0]


def main() -> int:
    if not INPUT.exists():
        raise SystemExit("run `python3 scripts/s7_run.py prep` first")
    rows: list[dict] = []
    timings: dict[str, float] = {}

    log_a = Path(str(STAGE_A) + ".log")
    if not log_a.exists():
        timings["stage_a_s"] = run(STAGE_A, ["--mrate", STAGE_A_RATES])
    rows_a = parse_models(log_a, "A_matrix")
    best_a = best_by_bic(rows_a)
    matrix = matrix_of(best_a["model"])
    print(f"stage A: {len(rows_a)} models, best {best_a['model']} "
          f"(BIC {best_a['BIC']}) -> matrix {matrix}")

    log_b = Path(str(STAGE_B) + ".log")
    if not log_b.exists():
        timings["stage_b_s"] = run(STAGE_B, ["--mset", matrix])
    rows_b = parse_models(log_b, "B_rate")
    best_b = best_by_bic(rows_b)
    final = reported_best(STAGE_B) or best_b["model"]
    print(f"stage B: {len(rows_b)} models, best {final} "
          f"(BIC {best_b['BIC']})")

    rows = rows_a + rows_b
    write_tsv(TABLE, ["stage", "rank", "model", "lnL", "df", "AIC",
                      "AICc", "BIC"], rows)

    # Stage A's margin — how much better the chosen matrix was than the
    # runner-up *matrix*, which is what the greedy step actually bets on.
    by_matrix: dict[str, float] = {}
    for r in rows_a:
        m = matrix_of(r["model"])
        by_matrix[m] = min(by_matrix.get(m, float("inf")), float(r["BIC"]))
    order = sorted(by_matrix.items(), key=lambda t: t[1])
    runner = order[1] if len(order) > 1 else ("", float("nan"))

    choice = {
        "final_model": final,
        "matrix": matrix,
        "stage_a": {"models_scored": len(rows_a),
                    "rates_tested": STAGE_A_RATES,
                    "matrices_scored": len(by_matrix),
                    "best_model": best_a["model"],
                    "best_bic": float(best_a["BIC"]),
                    "runner_up_matrix": runner[0],
                    "runner_up_bic": runner[1],
                    "bic_margin": round(runner[1] - order[0][1], 3)},
        "stage_b": {"models_scored": len(rows_b),
                    "best_model": final,
                    "best_bic": float(best_b["BIC"]),
                    "bic_gain_over_stage_a": round(
                        float(best_a["BIC"]) - float(best_b["BIC"]), 3)},
        "threads": THREADS,
        "seed": SEED,
        "exhaustive_probe": EXHAUSTIVE_PROBE,
        "runtimes_s": {k: round(v, 1) for k, v in timings.items()},
        "why_two_stage": ("the exhaustive -m MFP scan was measured at "
                          f"{EXHAUSTIVE_PROBE['projection_hours']} h on "
                          "this alignment; see the module docstring"),
    }
    CHOICE.write_text(json.dumps(choice, indent=1) + "\n")
    print(json.dumps(choice, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
