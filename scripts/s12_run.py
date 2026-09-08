"""The S12 driver — ordered stages, each resumable, none blocking the rest.

    refs -> crossmap -> runs -> quantify -> atlas -> tables
         -> figures -> report

`--only` / `--from` / `--list` select stages.  A failed stage does **not**
stop the ones after it that do not depend on it: `tables`, `figures` and
`report` render whatever has landed and the report marks an unfinished
section *not run yet*, because waking up to a partial S12 that says which
parts are partial is worth more than waking up to nothing (S9's rule).

The self-test runs **before anything is written** and a failure refuses the
run (`s10_run.py`'s rule).

Run:  python scripts/s12_run.py [--from quantify] [--spots 4000000]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

HERE = Path(__file__).resolve().parent
PY = sys.executable

STAGES = ["refs", "crossmap", "runs", "quantify", "atlas", "tables",
          "figures", "report"]


def cmd(stage: str, args) -> list[str]:
    base = [PY, str(HERE / f"s12_{stage}.py")]
    if stage in ("refs", "crossmap", "runs", "quantify", "atlas") \
            and args.species:
        base += ["--species", args.species]
    if stage == "quantify":
        base += ["--spots", str(args.spots), "--threads", str(args.threads)]
    if stage in ("refs", "crossmap"):
        base += ["--threads", str(args.threads)]
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--from", dest="start", default="")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--species", default="")
    ap.add_argument("--spots", type=int, default=4_000_000)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--skip-self-test", action="store_true")
    args = ap.parse_args()

    if args.list:
        for s in STAGES:
            print(s)
        return 0

    if not args.skip_self_test:
        import s12_test_expression
        s12_test_expression.self_test()

    todo = STAGES
    if args.start:
        if args.start not in STAGES:
            raise SystemExit(f"unknown stage {args.start}")
        todo = STAGES[STAGES.index(args.start):]
    if args.only:
        todo = [s for s in STAGES if s in args.only]

    failed: list[str] = []
    for stage in todo:
        print(f"\n=== [s12] {stage} " + "=" * (56 - len(stage)), flush=True)
        t0 = time.time()
        rc = subprocess.run(cmd(stage, args)).returncode
        print(f"--- {stage}: rc={rc} in {time.time() - t0:.0f}s", flush=True)
        if rc != 0:
            failed.append(stage)
    if failed:
        print(f"\nfailed stages: {', '.join(failed)}")
        return 1
    print("\nall stages ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
