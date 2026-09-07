"""S9 step 3 — run the codeml suite (PAML) defined in `s9_jobs.py`.

Every job is resumable — an `mlc` that already parses is reused, and a
`runmode -2` job is reused when its distance matrices are on disk — so the
suite can be run in stages and interrupted without losing work. Each job
has its own directory under `results/selection/codeml/<name>/`, which is
also how several codeml processes run at once: codeml writes fixed
filenames into its working directory and two jobs sharing one would
overwrite each other's output.

    python scripts/s9_codeml.py --list
    python scripts/s9_codeml.py --only 'm0_|pair_|two_ratio|bs_'
    python scripts/s9_codeml.py --workers 6

Results are collected into `codeml_results.json`, which the report and the
figures render from (D13). Nothing downstream re-parses an `mlc`.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import OUT_DIR  # noqa: E402
from s9_codeml_lib import owner_alive, run_job  # noqa: E402
from s9_jobs import CODEML_DIR, build_jobs  # noqa: E402

LIVE = ROOT / "results" / "session_live.json"
RESULTS_JSON = OUT_DIR / "codeml_results.json"

#: Rough relative cost per model, measured on this alignment. Used only to
#: order the queue longest-first: with a fixed worker count the makespan of
#: a mixed queue is set by whether the long jobs start early, and the
#: whole-tree branch models (57 tips) cost several times a per-paralog site
#: model (13-19 tips). Nothing about the results depends on it.
MODEL_COST = {"pair": 1, "m0": 4, "m1a": 8, "m2a": 20, "m7": 35, "m8": 45,
              "two": 60, "bs": 90}


def cost_rank(job) -> int:
    """Estimated cost of one job — bigger runs first."""
    model = job.name.split("_")[0]
    base = MODEL_COST.get(model, 10)
    # the whole-tree jobs read the full alignment, the paralog jobs a subset
    if job.seqfile.name == "codon_trimmed.phy":
        base *= 3
    return base


def live(stage: str, done: int, total: int, note: str = "") -> None:
    LIVE.write_text(json.dumps({
        "task": "S9", "stage": stage,
        "steps": [{"label": stage, "done": done, "total": total}],
        "note": note, "ts": time.strftime("%H:%M:%S")}))


def job_complete(job) -> bool:
    """Output on disk that a *finished* run leaves.

    An existing `mlc` is not enough: codeml writes it at the start and
    appends, so a job that is still running — or one that was killed —
    leaves a file that exists and holds no likelihood.
    """
    if job.pairwise:
        m = job.workdir / "2ML.dS"
        return m.exists() and m.stat().st_size > 0
    mlc = job.mlc()
    return mlc.exists() and "lnL(ntime" in mlc.read_text()


def load_results() -> dict:
    if RESULTS_JSON.exists():
        try:
            return json.loads(RESULTS_JSON.read_text())
        except Exception:  # noqa: BLE001
            pass
    return {}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--skip-slow", action="store_true",
                    help="omit the M7/M8 site models")
    ap.add_argument("--only", default="",
                    help="regex over job names")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--timeout-h", type=float, default=10.0)
    args = ap.parse_args()

    CODEML_DIR.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(args.skip_slow)
    if args.only:
        pat = re.compile(args.only)
        jobs = [j for j in jobs if pat.search(j.name)]
    jobs.sort(key=cost_rank, reverse=True)
    if args.list:
        n_done = 0
        for j in jobs:
            state = "    "
            if job_complete(j):
                state, n_done = "done", n_done + 1
            elif owner_alive(j.workdir / "RUNNING"):
                state = "run "
            print(f"{state}  {j.name:24} {j.description}")
        print(f"\n{n_done} of {len(jobs)} complete")
        return

    print(f"{len(jobs)} codeml jobs -> {CODEML_DIR}")
    results = load_results()
    done = 0

    def work(job):
        # "Was it already on disk?" has to be asked *before* run_job, which
        # returns a reused result indistinguishable from a fresh one. A
        # resumed job that recorded 0.0 minutes would quietly turn the
        # report's CPU-time figure into a fiction.
        reused = job_complete(job)
        t0 = time.time()
        res = run_job(job, timeout_s=int(args.timeout_h * 3600))
        return job, res, (None if reused else time.time() - t0)

    # as_completed, not map: map yields in submission order, so one slow
    # site model would hold back every result behind it — including from
    # the json, which is what a resumed run reads.
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(work, j) for j in jobs]
        for fut in as_completed(futures):
            job, res, secs = fut.result()
            done += 1
            live("codeml", done, len(jobs), job.name)
            if res is None:
                status = "FAILED/TIMEOUT"
            elif job.pairwise:
                status = "pairwise matrices written"
            else:
                status = (f"lnL={res.lnL:.3f} np={res.np} "
                          f"w={[round(x, 4) for x in res.omegas]}")
            print(f"[{done}/{len(jobs)}] {job.name:24} "
                  f"{'reused' if secs is None else f'{secs / 60:7.1f} min'}  "
                  f"{status}", flush=True)
            if res:
                prev = results.get(job.name, {})
                row = {
                    "lnL": res.lnL, "np": res.np, "kappa": res.kappa,
                    "omegas": res.omegas, "tree_length": res.tree_length,
                    "site_classes": res.site_classes,
                    "description": job.description,
                    "settings": {k: str(v) for k, v in job.settings.items()},
                    "minutes": (prev.get("minutes") if secs is None
                                else round(secs / 60, 2))}
                results[job.name] = row
                # Per-job, beside its own output: two drivers can run at
                # once (the site models and the whole-tree branch models
                # are hours apart in cost), and a single shared json would
                # let the second clobber the first's entries.
                (job.workdir / "result.json").write_text(
                    json.dumps(row, indent=2, sort_keys=True))
    RESULTS_JSON.write_text(json.dumps(results, indent=2, sort_keys=True))
    live("codeml", len(jobs), len(jobs), f"{len(results)} results")
    print(f"\n{len(results)} results -> {RESULTS_JSON}")


if __name__ == "__main__":
    main()
