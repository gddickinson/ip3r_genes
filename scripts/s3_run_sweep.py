"""S3 steps 2–3 — the two-profile hmmsearch sweep and the jackhmmer
convergence runs.

Against a concatenated reference-proteome DB built by s3_fetch_proteomes.py:

  1. `hmmsearch --domtblout` with **both** profiles at E ≤ 1e-5, then
     best-profile assignment with a bit-score margin (s3_assign.py). Running
     both is the point: a one-profile sweep cannot tell an ITPR from the
     ryanodine receptor sitting in the same result set (D14).
  2. `jackhmmer` from three single-sequence seeds — human ITPR1, the fly
     Itpr and the *Acanthamoeba* ITPR — iterated to convergence under a
     10-round ceiling, with per-round inclusion counts *and* per-round
     included target lists parsed from the log so D10's kill criterion
     (s3_kill.py) is evaluated on evidence rather than by eye.

Reporting and inclusion thresholds are both E ≤ 1e-5, so the iteration
criterion matches the claimed sweep threshold.

Raw outputs (domtblout + logs) → <data_root>/hmmer/. Summaries → the
committed results/hmm_sweep/.

Run:  python3 scripts/s3_run_sweep.py --stage hmmsearch
      python3 scripts/s3_run_sweep.py --stage jackhmmer --seed-tag itpr1_human
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_assign import (  # noqa: E402
    ASSIGN_FIELDS, assign, load_profile_hits,
)
from scripts.s3_hmm_lib import (  # noqa: E402
    HMM_SWEEP_DIR, parse_jackhmmer_log, read_fasta, write_fasta, write_tsv,
)
from scripts.s3_seed_spec import JACKHMMER_SEEDS  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

PROFILE_NAMES = ("itpr", "ryr")


def log(msg: str) -> None:
    print(f"[s3_sweep] {msg}", flush=True)


def run(cmd: list[str]) -> float:
    log("$ " + " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log(f"  rc={proc.returncode} in {dt:.0f}s")
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-3000:])
        raise SystemExit(f"command failed: {' '.join(cmd[:2])}")
    return dt


def live(task: str, steps: list[tuple[str, bool]]) -> None:
    """Dashboard live panel (session protocol)."""
    (PROJECT_ROOT / "results").mkdir(exist_ok=True)
    (PROJECT_ROOT / "results" / "session_live.json").write_text(json.dumps({
        "task": task,
        "steps": [{"label": lbl, "done": done} for lbl, done in steps],
    }, indent=1))


# ---------------------------------------------------------------- hmmsearch
def run_hmmsearch(db: Path, raw_dir: Path, cpu: int, evalue: str) -> dict:
    """Both profiles over the sweep DB, then the margin assignment."""
    stats: dict = {"profiles": {}}
    domtbls: dict[str, Path] = {}
    for name in PROFILE_NAMES:
        hmm = HMM_SWEEP_DIR / f"{name}.hmm"
        if not hmm.exists():
            raise SystemExit(f"profile missing: {hmm} — run s3_build_seed.py")
        domtbl = raw_dir / f"hmmsearch_{name}_vs_vertref.domtblout"
        logf = raw_dir / f"hmmsearch_{name}_vs_vertref.log"
        dt = run(["hmmsearch", "--cpu", str(cpu), "-E", evalue, "--noali",
                  "--domtblout", str(domtbl), "-o", str(logf),
                  str(hmm), str(db)])
        domtbls[name] = domtbl
        hits = load_profile_hits(domtbl)
        stats["profiles"][name] = {"targets": len(hits),
                                   "elapsed_s": round(dt, 1),
                                   "domtblout": str(domtbl)}
        log(f"{name}.hmm: {len(hits)} targets at E<={evalue}")
        live("S3", [("hmmsearch itpr.hmm", True),
                    ("hmmsearch ryr.hmm", name == "ryr"),
                    ("jackhmmer ×3", False), ("census v3", False)])

    rows = assign(load_profile_hits(domtbls["itpr"]),
                  load_profile_hits(domtbls["ryr"]))
    write_tsv(HMM_SWEEP_DIR / "hmmsearch_assignments.tsv", ASSIGN_FIELDS, rows)
    calls = {c: sum(1 for r in rows if r["assignment"] == c)
             for c in ("ITPR", "RYR", "unassigned")}
    stats["assignment"] = {"targets": len(rows), "calls": calls}
    log(f"assignment: {len(rows)} targets → {calls}")
    return stats


# ---------------------------------------------------------------- jackhmmer
def parse_jackhmmer(tag: str, accession: str, raw_dir: Path,
                    elapsed_s: float | None = None) -> dict:
    """Turn an archived jackhmmer log into the committed convergence outputs.

    Separate from running the search so the parse can be redone offline
    against the archived log — the same discipline S2 applied to InterPro's
    raw pages. `--parse-only` re-runs just this, which matters because the
    per-round *included target lists* D10 is evaluated on come from the log
    and nothing else.
    """
    logf = raw_dir / f"jackhmmer_{tag}.log"
    if not logf.exists():
        raise SystemExit(f"no archived log for {tag}: {logf}")
    conv = parse_jackhmmer_log(logf)
    rounds = conv["rounds"]
    log(f"jackhmmer {tag}: {len(rounds)} rounds, "
        f"converged={conv['converged']}, "
        f"trailing incomplete round={conv['trailing_incomplete_round']}")

    write_tsv(HMM_SWEEP_DIR / f"jackhmmer_convergence_{tag}.tsv",
              ["seed_tag", "accession", "round", "new_targets",
               "targets_in_msa", "n_included", "converged"],
              [{"seed_tag": tag, "accession": accession, "round": r["round"],
                "new_targets": r["new_targets"],
                "targets_in_msa": r["targets_in_msa"] if r["targets_in_msa"]
                is not None else "",
                "n_included": len(r["included"]),
                "converged": conv["converged"]} for r in rounds])
    # The per-round included sets are bulky, so they are archived beside the
    # raw log rather than committed.
    (raw_dir / f"jackhmmer_{tag}_rounds.json").write_text(json.dumps(
        {"tag": tag, "accession": accession, "converged": conv["converged"],
         "trailing_incomplete_round": conv["trailing_incomplete_round"],
         "rounds": rounds}, indent=1))
    out = {"tag": tag, "accession": accession, "rounds": len(rounds),
           "converged": conv["converged"],
           "domtblout": str(raw_dir / f"jackhmmer_{tag}.domtblout"),
           "log": str(logf),
           "rounds_json": str(raw_dir / f"jackhmmer_{tag}_rounds.json"),
           "new_targets": [r["new_targets"] for r in rounds]}
    if elapsed_s is not None:
        out["elapsed_s"] = round(elapsed_s, 1)
    else:
        # A re-parse must not erase how long the search took: that number
        # can only be measured while it runs.
        prior = HMM_SWEEP_DIR / f"sweep_stats_jackhmmer_{tag}.json"
        if prior.exists():
            was = json.loads(prior.read_text()).get(
                "jackhmmer", {}).get(tag, {}).get("elapsed_s")
            if was is not None:
                out["elapsed_s"] = was
    return out


def run_jackhmmer(tag: str, accession: str, db: Path, raw_dir: Path,
                  cpu: int, evalue: str, max_iter: int) -> dict:
    """One jackhmmer run from a single-sequence seed, to convergence."""
    seeds = read_fasta(HMM_SWEEP_DIR / "itpr_seed.faa")
    match = [(h, s) for h, s in seeds.items() if h.startswith(accession + "|")]
    if not match:
        raise SystemExit(f"jackhmmer seed {accession} not in itpr_seed.faa")
    header, seq = match[0]
    seed_fa = raw_dir / f"jackhmmer_{tag}_seed.faa"
    write_fasta(seed_fa, [(header, seq)])

    domtbl = raw_dir / f"jackhmmer_{tag}.domtblout"
    logf = raw_dir / f"jackhmmer_{tag}.log"
    dt = run(["jackhmmer", "--cpu", str(cpu), "-N", str(max_iter),
              "-E", evalue, "--incE", evalue, "--noali",
              "--domtblout", str(domtbl), "-o", str(logf),
              str(seed_fa), str(db)])
    return parse_jackhmmer(tag, accession, raw_dir, elapsed_s=dt)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path,
                    default=require_data_root() / "proteomes" /
                    "vertebrata_refprot.fasta")
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--max-iter", type=int, default=10)
    ap.add_argument("--stage", choices=["all", "hmmsearch", "jackhmmer"],
                    default="all")
    ap.add_argument("--parse-only", action="store_true",
                    help="re-derive the jackhmmer convergence tables from "
                         "the archived logs without searching again")
    ap.add_argument("--seed-tag", action="append",
                    choices=sorted(JACKHMMER_SEEDS),
                    help="jackhmmer seed(s) to run (default: all); separate "
                         "processes may run one tag each in parallel — every "
                         "output file is per-tag")
    args = ap.parse_args()
    if not args.parse_only and not args.db.exists():
        raise SystemExit(f"sweep DB missing: {args.db} "
                         "— run s3_fetch_proteomes.py")

    raw_dir = require_data_root() / "hmmer"
    raw_dir.mkdir(parents=True, exist_ok=True)
    stats: dict = {"db": str(args.db), "evalue": args.evalue,
                   "cpu": args.cpu, "max_iter": args.max_iter}

    if args.stage in ("all", "hmmsearch"):
        stats["hmmsearch"] = run_hmmsearch(
            args.db, raw_dir, args.cpu, args.evalue)

    tags = args.seed_tag or sorted(JACKHMMER_SEEDS)
    if args.stage in ("all", "jackhmmer"):
        stats["jackhmmer"] = {
            tag: (parse_jackhmmer(tag, JACKHMMER_SEEDS[tag][0], raw_dir)
                  if args.parse_only else
                  run_jackhmmer(tag, JACKHMMER_SEEDS[tag][0], args.db, raw_dir,
                                args.cpu, args.evalue, args.max_iter))
            for tag in tags}

    suffix = args.stage if args.stage != "jackhmmer" \
        else "jackhmmer_" + "_".join(tags)
    out = HMM_SWEEP_DIR / f"sweep_stats_{suffix}.json"
    out.write_text(json.dumps(stats, indent=1))
    log(f"stage '{args.stage}' complete → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
