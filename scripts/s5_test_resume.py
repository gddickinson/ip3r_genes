"""s5_test_resume.py — the two properties a resumable sweep has to have.

The S5 driver is resumable: it skips any genome carrying a `.sweep.done`
marker and reuses any existing non-empty `miniprot.gff`. That is what makes a
309-genome sweep survivable across interruptions, and it is only safe if two
things hold.

**1. A failed or killed run leaves nothing reusable.** miniprot's output used
to be written straight to `miniprot.gff`, so a run killed mid-write left a
truncated file that the next run silently accepted as complete — and a GFF cut
at a line boundary parses perfectly, so the resulting genome would carry a
clean-looking ledger row that is simply short of loci. There is no way to
detect it afterwards. `run_miniprot` now writes to a `.partial` sibling and
renames on success; this test proves a failure leaves neither file.

**2. Re-running a genome reproduces its result.** Otherwise a `--redo`, a
resumed shard and a first run are three different experiments, and the ledger
depends on which one happened to write each row. `--full` extends this through
miniprot itself rather than only the parse and classification, which is the
stronger claim and the one D24 cares about (S3 found MAFFT is not reproducible
at `--thread -1`; the same question has to be asked of every tool here).

Fields excluded from the comparison are the ones that *must* differ between
runs — wall-clock timings and the run timestamp — and they are listed
explicitly rather than filtered by pattern, so a field that starts varying
cannot hide behind a wildcard.

Usage:
  python scripts/s5_test_resume.py                    # both, cached miniprot
  python scripts/s5_test_resume.py --full             # re-run miniprot too
  python scripts/s5_test_resume.py --accession GCF_...
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

for extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(extra).is_dir() and extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root              # noqa: E402
import s5_sweep_lib as lib                                 # noqa: E402

SWEEP = get_data_root() / "genome_sweep"

#: Fields that legitimately differ between two runs of the same genome.
VOLATILE = {"timings", "total_s", "run_at"}


def strip_volatile(summary: dict) -> dict:
    out = copy.deepcopy(summary)
    for k in VOLATILE:
        out.pop(k, None)
    return out


def diff_paths(a, b, path: str = "") -> list[str]:
    """Every leaf where two summaries disagree, addressed by path."""
    if type(a) is not type(b):
        return [f"{path}: type {type(a).__name__} vs {type(b).__name__}"]
    if isinstance(a, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a:
                out.append(f"{path}.{k}: missing in first")
            elif k not in b:
                out.append(f"{path}.{k}: missing in second")
            else:
                out += diff_paths(a[k], b[k], f"{path}.{k}")
        return out
    if isinstance(a, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff_paths(x, y, f"{path}[{i}]")
        return out
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def test_atomic_write() -> list[str]:
    """A failing miniprot must leave neither the target nor a partial."""
    work = SWEEP / "_resume_test"
    work.mkdir(parents=True, exist_ok=True)
    target = work / "miniprot.gff"
    target.unlink(missing_ok=True)
    real_run = subprocess.run
    failures: list[str] = []

    def fake_run(cmd, *a, **kw):
        # Let miniprot "write" some output, then fail — the exact shape of a
        # killed run: a non-empty file and a non-zero exit.
        if kw.get("stdout") is not None and hasattr(kw["stdout"], "write"):
            kw["stdout"].write("##gff-version 3\nNC_1\tmp\tmRNA\t1\t9\t.\t+\t.\tID=MP1\n")
        return subprocess.CompletedProcess(cmd, 1, "", "simulated failure")

    subprocess.run = fake_run                        # noqa: S001 — test double
    try:
        try:
            lib.run_miniprot(Path("/dev/null"), Path("/dev/null"), target)
        except RuntimeError:
            pass                                     # expected
        else:
            failures.append("run_miniprot did not raise on a non-zero exit")
    finally:
        subprocess.run = real_run

    if target.exists():
        failures.append(f"a failed run left {target.name} "
                        f"({target.stat().st_size} bytes) — the next run "
                        "would reuse it as complete")
    leftovers = list(work.glob("*.partial"))
    if leftovers:
        failures.append(f"a failed run left {len(leftovers)} .partial file(s)")
    for f in work.glob("*"):
        f.unlink()
    work.rmdir()
    return failures


def test_idempotent(accession: str, full: bool) -> list[str]:
    """Re-running a swept genome must reproduce its summary."""
    d = SWEEP / accession
    summary = d / "summary.json"
    if not (d / ".sweep.done").exists() or not summary.exists():
        return [f"{accession} is not a completed sweep — nothing to compare"]
    before = json.loads(summary.read_text())
    if full:
        (d / "miniprot.gff").unlink(missing_ok=True)
    proc = subprocess.run(
        [sys.executable, "scripts/s5_run_sweep.py", "--only", accession,
         "--redo", "--threads", "2"],
        cwd=PROJECT_ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        return [f"re-run failed: {proc.stderr.strip()[:300]}"]
    after = json.loads(summary.read_text())
    diffs = diff_paths(strip_volatile(before), strip_volatile(after))
    return [f"summary differs at {d}" for d in diffs[:12]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", default="GCF_901000725.3",
                    help="a completed sweep to re-run (default: Takifugu)")
    ap.add_argument("--full", action="store_true",
                    help="delete the cached GFF so miniprot re-runs too")
    args = ap.parse_args()

    results: dict[str, list[str]] = {}
    print("1. a failed run leaves nothing reusable")
    results["atomic_write"] = test_atomic_write()
    print("   " + ("PASS" if not results["atomic_write"] else "FAIL"))

    print(f"2. re-running {args.accession} reproduces its summary"
          + (" (miniprot included)" if args.full else " (cached miniprot)"))
    results["idempotent"] = test_idempotent(args.accession, args.full)
    print("   " + ("PASS" if not results["idempotent"] else "FAIL"))

    failed = {k: v for k, v in results.items() if v}
    if failed:
        print("\nFAILURES:")
        for name, msgs in failed.items():
            for m in msgs:
                print(f"  {name}: {m}")
        raise SystemExit(1)
    print(f"\nPASS: both resume properties hold"
          + (" through miniprot" if args.full else ""))


if __name__ == "__main__":
    main()
