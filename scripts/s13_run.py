"""The S13 driver — ordered stages, each resumable, none blocking the rest.

    python3 scripts/s13_run.py --list
    python3 scripts/s13_run.py                    # everything
    python3 scripts/s13_run.py --only reconcile
    python3 scripts/s13_run.py --from figures

Runs `s13_test_recon.py` **before writing anything** and refuses to continue
if it fails: the reconciliation's every rule returns a plausible number when it
is wrong, and two of them already did.

A failed stage does not stop the ones after it that do not depend on it, and
the report marks an unfinished section *not run yet* — waking up to a partial
S13 that says which parts are partial is worth more than waking up to nothing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

RECON = ROOT / "results" / "reconciliation"
PY = sys.executable
# The figures need matplotlib, which lives in the recorded conda env (D18).
FIG_PY = "/opt/anaconda3/envs/piezo1/bin/python"

STAGES = [
    ("species_tree", [PY, "scripts/s13_species_tree.py"],
     "the curated calibrated species tree (D15) — an input, validated first"),
    ("reconcile", [PY, "scripts/s13_reconcile.py"],
     "the topology x cyclostome x support matrix, plus the rooting check"),
    ("losses", [PY, "scripts/s13_losses.py"],
     "every implied loss asked of the S5 genome ledger"),
    ("cyclostome", [PY, "scripts/s13_cyclostome.py"],
     "the six tips the deep placement rests on: S8's call and the "
     "long-branch check"),
    ("stats", None, "SHA-256 of every table, parameters, self-test status"),
    ("figures", [FIG_PY, "scripts/s13_figures.py"], "the four figures"),
    ("report", [PY, "scripts/s13_report.py"], "report.md, rendered from the "
     "committed tables only (D13)"),
]

TABLES = ["species_tree.nwk", "species_tree_calibrations.tsv",
          "reconciliation_events.tsv", "duplication_placement.tsv",
          "reconciliation_summary.tsv", "rooting_check.tsv",
          "paralog_presence.tsv", "loss_verdicts.tsv",
          "corroborated_losses.tsv", "cyclostome_loci.tsv",
          "branch_lengths.tsv"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_stats(selftest_ok: bool, problems: list) -> None:
    import s13_lib as L
    import s13_losses as LO
    import s13_reconcile as R
    stats = {
        "task": "S13",
        "written_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "self_test": {"passed": selftest_ok, "n_groups": 0,
                      "problems": problems},
        "parameters": {
            "support_bar": {"sh_alrt": L.MIN_ALRT, "ufboot": L.MIN_UFBOOT,
                            "source": "S7's own 'well supported' bar"},
            "topologies": {k: v[1] for k, v in R.TOPOLOGIES.items()},
            "variants": list(R.VARIANTS),
            "duplication_rule": "non-binary LCA (Vernot et al. 2008); losses "
                                "per Zmasek & Eddy 2001",
            "ledger_status_groups": {
                "present": sorted(LO.PRESENT), "absent": sorted(LO.ABSENT),
                "remnant": sorted(LO.REMNANT),
                "undecidable": sorted(LO.UNDECIDABLE)},
        },
        "sha256": {},
    }
    import s13_test_recon as T
    stats["self_test"]["n_groups"] = len(T.TESTS)
    for name in TABLES:
        p = RECON / name
        if p.exists():
            stats["sha256"][name] = sha256(p)
    (RECON / "reconciliation_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n")
    print(f"wrote {RECON/'reconciliation_stats.json'} "
          f"({len(stats['sha256'])} tables)")


def main() -> int:
    ap = argparse.ArgumentParser()
    names = [s[0] for s in STAGES]
    ap.add_argument("--only", action="append", choices=names)
    ap.add_argument("--from", dest="start", choices=names)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--skip-tests", action="store_true",
                    help="do not run the negative controls first (never use "
                         "this for a build whose output is committed)")
    args = ap.parse_args()

    if args.list:
        for n, _, why in STAGES:
            print(f"{n:14s} {why}")
        return 0

    problems: list = []
    if not args.skip_tests:
        import s13_test_recon as T
        print(f"S13 negative controls ({len(T.TESTS)} groups)")
        problems = T.self_test(verbose=False)
        for p in problems:
            print(f"  ! {p}")
        if problems:
            print("refusing to write: the self-tests failed")
            return 1
        print("  all passed\n")

    todo = names
    if args.only:
        todo = [n for n in names if n in args.only]
    elif args.start:
        todo = names[names.index(args.start):]

    RECON.mkdir(parents=True, exist_ok=True)
    failed = []
    for name, cmd, why in STAGES:
        if name not in todo:
            continue
        print(f"--- {name}: {why}")
        t0 = time.time()
        if cmd is None:
            write_stats(not problems, problems)
            rc = 0
        else:
            rc = subprocess.run(cmd, cwd=ROOT).returncode
        print(f"--- {name}: rc={rc} in {time.time()-t0:.1f}s\n")
        if rc:
            failed.append(name)
    if failed:
        print(f"failed stages: {', '.join(failed)}")
        return 1
    print("S13 complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
