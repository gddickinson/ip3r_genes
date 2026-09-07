"""S9b driver — the whole model stage, in order, unattended.

    codeml  ->  relax  ->  tables  ->  report  ->  figures

Every stage is resumable, so this can be started, interrupted and started
again: codeml reuses any job whose output is complete, RELAX reuses any
paralog whose json is on disk, and the last three stages are pure
re-derivations from the committed tables.

    python scripts/s9_run.py                  # everything still outstanding
    python scripts/s9_run.py --list
    python scripts/s9_run.py --from tables
    python scripts/s9_run.py --only codeml --workers 7

A stage that fails does **not** stop the ones after it that do not depend
on it — `tables`, `report` and `figures` render whatever has landed, and
the report marks an unfinished section *not run yet* rather than omitting
it. That is deliberate for an overnight run: waking up to a partial S9 that
says which parts are partial is worth more than waking up to nothing.

The exit code is non-zero if any stage failed, so a wrapper can tell.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = Path(__file__).resolve().parent
PY = sys.executable

STAGES = ("codeml", "relax", "tables", "report", "figures")


def stage_cmd(stage: str, workers: int, relax_workers: int) -> list[str]:
    if stage == "codeml":
        return [PY, str(SCRIPTS / "s9_codeml.py"), "--workers", str(workers),
                "--timeout-h", "16"]
    if stage == "relax":
        return [PY, str(SCRIPTS / "s9_relax.py"),
                "--workers", str(relax_workers), "--timeout-h", "16"]
    return [PY, str(SCRIPTS / f"s9_{stage}.py")]


def run_stage(stage: str, workers: int, relax_workers: int) -> tuple[bool, float]:
    cmd = stage_cmd(stage, workers, relax_workers)
    print(f"\n{'=' * 70}\n== {stage}  ({time.strftime('%Y-%m-%d %H:%M:%S')})\n"
          f"== {' '.join(cmd)}\n{'=' * 70}", flush=True)
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=ROOT)
    secs = time.time() - t0
    ok = proc.returncode == 0
    print(f"== {stage}: {'ok' if ok else f'FAILED rc={proc.returncode}'} "
          f"in {secs / 60:.1f} min", flush=True)
    return ok, secs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7,
                    help="concurrent codeml processes")
    ap.add_argument("--relax-workers", type=int, default=6,
                    help="threads HyPhy may use per RELAX run")
    ap.add_argument("--only", action="append", default=[],
                    choices=STAGES)
    ap.add_argument("--from", dest="start", default="", choices=("",) + STAGES)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for s in STAGES:
            print(s)
        return 0

    todo = list(STAGES)
    if args.start:
        todo = todo[STAGES.index(args.start):]
    if args.only:
        todo = [s for s in todo if s in args.only]

    failed: list[str] = []
    t0 = time.time()
    for stage in todo:
        ok, _ = run_stage(stage, args.workers, args.relax_workers)
        if not ok:
            failed.append(stage)
    print(f"\n{'=' * 70}\nS9b finished in {(time.time() - t0) / 3600:.2f} h  "
          f"({time.strftime('%Y-%m-%d %H:%M:%S')})", flush=True)
    if failed:
        print("FAILED stages: " + ", ".join(failed))
    else:
        print("all stages ok -> results/selection/report.md")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
