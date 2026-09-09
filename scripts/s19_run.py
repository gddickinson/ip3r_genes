"""S19 driver — methods results: what the search was worth.

Stages, in order (later ones read earlier ones' tables):

  contribution  per-channel contribution, head-to-head inside one database,
                per-gene recovery, and cost against yield
  contiguity    the false-negative rate on a family with no losses, and the
                floor scan that calibrates D4's a-priori contiguity bar
  panel         the exact bait-panel ablation, re-clustered from the 309
                retained miniprot GFFs, plus the bait-identity design rule
  drift         iterative-search drift recoded from the raw jackhmmer logs,
                with each kill rule scored against a measured outcome
  inference     where the inference methods ran out (S9/S13/S15/S17 limits)
  tables        methods_stats.json — rules, parameters, self-test, SHA-256
  figures       the four figures, from committed tables only
  report        report.md, rendered purely from the committed tables

Everything runs from committed artefacts; nothing downloads. The one job
that touches bulk data is the accession universe — one `grep '^>'` pass per
reference-proteome FASTA — and it is cached under the data root, so a rerun
is offline and deterministic.

`s19_test_methods.py` runs **before anything is written** and the driver
refuses to continue if it fails. A failed stage does not stop the ones after
it that do not depend on it, and the report marks an unfinished section
*not run yet*.

Usage:

    python scripts/s19_run.py                 # every stage
    python scripts/s19_run.py --only panel
    python scripts/s19_run.py --from drift
    python scripts/s19_run.py --list
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402

STAGES = ("contribution", "contiguity", "panel", "drift", "inference",
          "tables", "figures", "report")

#: Stages that must not run on a stale table written by an earlier failure.
#: `tables` and `report` read everything, so they run last by construction.
MODULE = {
    "contribution": "s19_contribution", "contiguity": "s19_contiguity",
    "panel": "s19_panel", "drift": "s19_drift", "inference": "s19_inference",
    "tables": "s19_tables", "figures": "s19_figures", "report": "s19_report",
}


def self_test() -> bool:
    path = Path(__file__).with_name("s19_test_methods.py")
    proc = subprocess.run([sys.executable, str(path)], text=True)
    return proc.returncode == 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", action="append", choices=STAGES,
                    help="run just this stage (repeatable)")
    ap.add_argument("--from", dest="start", choices=STAGES,
                    help="run this stage and everything after it")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    ap.add_argument("--skip-self-test", action="store_true",
                    help="for debugging a stage the self-test depends on")
    args = ap.parse_args()

    if args.list:
        for s in STAGES:
            print(s)
        return 0

    todo = list(args.only) if args.only else list(STAGES)
    if args.start:
        todo = list(STAGES[STAGES.index(args.start):])
    S.out_dir()

    # The suite runs twice. Before any stage writes, only the checks on
    # constructed inputs can run — the four that read a committed table
    # report as skipped. The `tables` stage re-runs the whole suite once
    # those tables exist and records the result in `methods_stats.json`, so
    # a check that could not run before is not silently never run.
    if not args.skip_self_test and not self_test():
        S.log("negative controls FAILED — refusing to write")
        return 1

    failed = []
    for name in todo:
        t0 = time.time()
        S.log(f"=== stage {name}")
        try:
            module = __import__(MODULE[name])
            module.run()
            S.log(f"=== {name} done in {time.time() - t0:.1f}s")
        except Exception:                                    # noqa: BLE001
            traceback.print_exc()
            failed.append(name)
            S.log(f"=== {name} FAILED after {time.time() - t0:.1f}s")
    status = S.read_json(S.OUT_DIR / "methods_stats.json").get("self_test", {})
    if status:
        S.log(f"negative controls after the run: {status.get('n_checks')} "
              f"checks, {status.get('n_failed')} failed")
    if failed:
        S.log(f"S19 stages failed: {', '.join(failed)}")
        return 1
    S.log(f"S19 stages complete: {', '.join(todo)} -> {S.OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
