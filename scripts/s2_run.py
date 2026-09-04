"""S2 driver — the census-v2 stages in order.

    python scripts/s2_run.py --list
    python scripts/s2_run.py                 # everything
    python scripts/s2_run.py --from call     # resume at a stage
    python scripts/s2_run.py --only report

Stages are separate processes rather than imported functions so that a
stage that needs the network (`interpro`, `uniprot`) and one that does not
(`call`, `figures`, `report`) can be run independently — every stage after
`sequences` reads only committed tables, so the analysis half reruns
offline in seconds.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STAGES = [
    ("interpro", "enumerate the seed signatures to exhaustion (network)"),
    ("uniprot", "architecture + lineage sweep (network)"),
    ("call", "the positive ITPR/RYR call on every record"),
    ("sequences", "seeded-space FASTA, representatives, core panel (network)"),
    ("verify", "sequence-level check of the call (MAFFT, ~6 min)"),
    ("delta", "census v1 → v2 delta, with a verdict per missing accession"),
    ("figures", "the five figures"),
    ("report", "render report.md from the committed tables"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", action="append")
    ap.add_argument("--from", dest="start")
    args = ap.parse_args()

    names = [s for s, _ in STAGES]
    if args.list:
        for s, why in STAGES:
            print(f"  {s:<10s} {why}")
        return 0
    run = args.only or names
    if args.start:
        if args.start not in names:
            raise SystemExit(f"unknown stage {args.start!r}; see --list")
        run = names[names.index(args.start):]

    for stage in run:
        if stage not in names:
            raise SystemExit(f"unknown stage {stage!r}; see --list")
        script = ROOT / "scripts" / f"s2_{stage}.py"
        print(f"\n=== s2:{stage} ===", flush=True)
        t0 = time.time()
        rc = subprocess.run([sys.executable, str(script)], cwd=ROOT).returncode
        print(f"=== s2:{stage} rc={rc} in {time.time() - t0:.0f}s", flush=True)
        if rc != 0:
            return rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
