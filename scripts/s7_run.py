"""S7 driver — IQ-TREE 2 ML phylogeny over S6's trimmed alignment.

Subcommands (stdlib-only; iqtree2 must be on PATH):

  python3 scripts/s7_run.py prep    # write input.fasta from the MSA
  python3 scripts/s7_run.py ml      # main search: 1000 UFBoot + 1000
                                    # SH-aLRT under the model chosen by
                                    # `s7_model.py` (resumable)
  python3 scripts/s7_run.py bnni    # the same search with --bnni, the
                                    # robustness check on UFBoot under
                                    # model violation (a deep, saturated
                                    # alignment is exactly the case the
                                    # flag exists for)
  python3 scripts/s7_run.py au      # constrained ITPR1/2/3 sister
                                    # topologies + AU test

**Threads are pinned, not AUTO (D24).** The brief's command says
`-T AUTO`; `-T AUTO` picks its thread count from the machine's current
load, and IQ-TREE is only reproducible at a *fixed* seed and a *fixed*
thread count. S6 made the same call for MAFFT, where `--thread -1` gave
two different profiles from the same seeds. Everything downstream of
this tree is built from it, so the thread count is recorded in
`tree_stats.json` beside the SHA-256 of the input and the treefile.

Inputs   results/msa_v2/trimmed.fasta + representatives.tsv
Outputs  results/phylogeny/ (input.fasta, itpr_ml.*, tree_stats.json,
                             au/*)
"""

from __future__ import annotations

import atexit
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.s6_lib import read_fasta, write_fasta, write_live  # noqa: E402
from scripts.s7_lib import (  # noqa: E402
    MSA_DIR, PHYLO_DIR, group_labels, load_groups, parse_newick,
)
from scripts.s7_model import CHOICE  # noqa: E402

TRIMMED = MSA_DIR / "trimmed.fasta"
INPUT = PHYLO_DIR / "input.fasta"
PREFIX = PHYLO_DIR / "itpr_ml"
TREEFILE = Path(str(PREFIX) + ".treefile")
IQREPORT = Path(str(PREFIX) + ".iqtree")
BNNI_PREFIX = PHYLO_DIR / "itpr_ml_bnni"
STATS = PHYLO_DIR / "tree_stats.json"
AU_DIR = PHYLO_DIR / "au"

THREADS = 4          # pinned — IQ-TREE's own benchmark (D24)
#: The three constrained searches are independent, so they run at once.
#: Their thread count is *lower and separately recorded* rather than
#: reused: 3 x 3 fits this machine's 10 cores, and D24 needs the number
#: each search actually ran at, not the number the main search used.
CONSTRAINT_THREADS = 3
SEED = 42
UFBOOT = 1000
ALRT = 1000
AU_REPLICATES = 10000


def claim(prefix: Path) -> None:
    """Refuse a `--prefix` a live process already owns.

    Two IQ-TREE runs sharing one prefix interleave their output in
    silence: each writes the same `.ckp.gz`, `.treefile` and `.log`, and
    whichever finishes last wins some files and loses others. It is not
    hypothetical — a previous session's `au` was interrupted, its three
    constrained searches were orphaned and kept running for 50 minutes,
    a fresh `au` started three more on the same prefixes, and the AU
    test was handed whichever `.treefile` happened to be on disk. The
    only visible symptom was a likelihood that looked wrong.

    A stale lock (the owner is gone) is taken over rather than treated
    as an error, so an interrupted run does not need cleaning by hand.
    """
    lock = Path(str(prefix) + ".lock")
    if lock.exists():
        try:
            pid = int(lock.read_text().split()[0])
        except (ValueError, IndexError):
            pid = 0
        alive = False
        if pid:
            try:
                os.kill(pid, 0)          # signal 0 — existence check only
                alive = True
            except ProcessLookupError:
                alive = False            # the owner is gone: stale lock
            except PermissionError:
                alive = True             # exists, owned by another user
        if alive:
            raise SystemExit(
                f"{prefix.name}: pid {pid} is already writing this prefix. "
                f"Wait for it, or kill it and delete {lock}.")
        print(f"note: taking over a stale lock from pid {pid} ({lock.name})")
    lock.write_text(f"{os.getpid()}\n")
    atexit.register(lambda: lock.unlink(missing_ok=True))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prep_input() -> tuple[int, int]:
    PHYLO_DIR.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(TRIMMED)
    write_fasta(seqs, INPUT)
    n, cols = len(seqs), len(next(iter(seqs.values())))
    print(f"input: {n} seqs x {cols} cols -> {INPUT}")
    return n, cols


def _steps() -> list[tuple[str, bool]]:
    """The live panel, read off disk rather than off the caller.

    Each stage used to pass its own booleans, so `run_au` announced the
    `--bnni` check as done whether or not it had ever run — and the
    dashboard showed a completed step for a search that was still going.
    A progress panel that reports what a function believes rather than
    what exists is worse than none.
    """
    return [
        ("prepare input", INPUT.exists()),
        ("model selection (2 stages)", CHOICE.exists()),
        ("IQ-TREE ML search", IQREPORT.exists()),
        ("UFBoot --bnni check",
         Path(str(BNNI_PREFIX) + ".iqtree").exists()),
        ("AU test", (PHYLO_DIR / "au_test.tsv").exists()),
        ("report + figures", (PHYLO_DIR / "report.md").exists()
         and (PHYLO_DIR / "figures" / "sister_au.png").exists()),
    ]


def chosen_model() -> str:
    if not CHOICE.exists():
        raise SystemExit("run `python3 scripts/s7_model.py` first")
    return json.loads(CHOICE.read_text())["final_model"]


def run_ml() -> int:
    n, cols = prep_input()
    # Resume on the .iqtree report, NOT the .treefile: IQ-TREE writes a
    # treefile at checkpoints during the run, so a treefile can exist
    # without support values and without the search having finished. The
    # .iqtree report is written once, at the end.
    if IQREPORT.exists():
        print(f"resume: {IQREPORT} already exists — skipping search")
        return 0
    claim(PREFIX)
    cmd = ["iqtree2", "-s", str(INPUT), "--prefix", str(PREFIX),
           "-m", chosen_model(), "-B", str(UFBOOT), "-alrt", str(ALRT),
           "-T", str(THREADS), "-seed", str(SEED)]
    print("running:", " ".join(cmd))
    write_live("S7", _steps())
    t0 = time.time()
    rc = subprocess.run(cmd).returncode
    dt = time.time() - t0
    print(f"iqtree2 rc={rc} in {dt:.0f} s")
    if rc:
        return rc
    record_stats(cmd, dt, n, cols)
    write_live("S7", _steps())
    return 0


def run_bnni() -> int:
    """Re-run with --bnni: UFBoot's own guard against model violation.

    Not a second answer — a check. If the ML topology and the clade
    supports the claims rest on move under `--bnni`, the supports are
    the artefact the flag exists to expose, and the report says so.
    """
    # Resume on the report, not the treefile (see `run_au`).
    rep = Path(str(BNNI_PREFIX) + ".iqtree")
    cmd, dt = None, None
    if not rep.exists():
        claim(BNNI_PREFIX)
        cmd = ["iqtree2", "-s", str(INPUT), "--prefix", str(BNNI_PREFIX),
               "-m", best_model(), "-B", str(UFBOOT), "-alrt", str(ALRT),
               "--bnni", "-T", str(THREADS), "-seed", str(SEED)]
        print("running:", " ".join(cmd))
        t0 = time.time()
        rc = subprocess.run(cmd).returncode
        dt = time.time() - t0
        print(f"iqtree2 --bnni rc={rc} in {dt:.0f} s")
        if rc:
            return rc
    if rep.exists():
        _record_bnni_stats(cmd, dt)
    write_live("S7", _steps())
    return 0


def logged_command(prefix: Path) -> str:
    """The command IQ-TREE itself recorded in its log.

    Read back rather than reconstructed, so a run resumed in a later
    invocation records the command that actually produced the tree —
    including a run whose driver was a different version of this file.
    """
    log = Path(str(prefix) + ".log")
    if not log.exists():
        return ""
    for line in log.read_text().splitlines()[:40]:
        if line.startswith("Command:"):
            return line.split(":", 1)[1].strip()
    return ""


def _record_bnni_stats(cmd: list[str] | None,
                       runtime_s: float | None) -> None:
    """The check's own settings, beside the search's (D24).

    A robustness check whose command and log-likelihood are not written
    down cannot be re-run to the same answer, and the whole point of it
    is that its answer is comparable with the main run's.
    """
    if not STATS.exists():
        return
    rep = Path(str(BNNI_PREFIX) + ".iqtree")
    d = json.loads(STATS.read_text())
    d["bnni_check"] = {
        "command": " ".join(cmd) if cmd else logged_command(BNNI_PREFIX),
        "threads": THREADS,
        "seed": SEED,
        "runtime_s": (round(runtime_s, 1) if runtime_s is not None
                      else d.get("bnni_check", {}).get("runtime_s", "")),
        "model": best_model(rep) if rep.exists() else "",
        "log_likelihood": log_likelihood(rep) if rep.exists() else None,
        "sha256_treefile": sha256(Path(str(BNNI_PREFIX) + ".treefile")),
        "why": ("UFBoot is optimistic under model violation; --bnni is "
                "its own guard. A check, not a second answer — the "
                "reported tree is the main search's. Claim-by-claim "
                "comparison: bnni_comparison.tsv"),
    }
    STATS.write_text(json.dumps(d, indent=1) + "\n")
    print(f"recorded the --bnni run in {STATS.name}")


def best_model(report: Path = IQREPORT) -> str:
    """The model the tree was actually built under.

    Read back out of the .iqtree file rather than from `model_choice.json`,
    so a treefile and the model quoted beside it cannot come from two
    different runs.
    """
    m = re.search(r"Model of substitution: (\S+)", report.read_text())
    return m.group(1) if m else chosen_model()


def log_likelihood(report: Path = IQREPORT) -> float | None:
    m = re.search(r"Log-likelihood of the tree: (-?[\d.]+)",
                  report.read_text())
    return float(m.group(1)) if m else None


def record_stats(cmd: list[str], runtime_s: float, n: int, cols: int) -> None:
    ver = subprocess.run(["iqtree2", "--version"], capture_output=True,
                         text=True).stdout.split("\n")[0].strip()
    stats = {
        "iqtree": {
            "version": ver,
            "command": " ".join(cmd),
            "threads": THREADS,
            "thread_note": ("D24 — pinned; -T AUTO chooses from machine "
                            "load and is not reproducible"),
            "seed": SEED,
            "runtime_s": round(runtime_s, 1),
            "n_sequences": n,
            "columns": cols,
            "model": best_model(),
            "model_selected_by": "scripts/s7_model.py (two-stage, greedy)",
            "log_likelihood": log_likelihood(),
            "ufboot_replicates": UFBOOT,
            "alrt_replicates": ALRT,
            "sha256_in": sha256(INPUT),
            "sha256_treefile": sha256(TREEFILE),
        },
    }
    STATS.write_text(json.dumps(stats, indent=1) + "\n")
    print(f"wrote {STATS}")


# ------------------------------------------------------------------ AU test

def run_au() -> int:
    # imported here, not at module load, so `ml` can start before the
    # constraint logic (which needs the ML tree) exists.
    from scripts.s7_constraints import au_constraints, corrected_membership
    AU_DIR.mkdir(parents=True, exist_ok=True)
    model = best_model()
    groups = load_groups()
    tree = parse_newick(TREEFILE.read_text())
    all_leaves = tree.leaf_names()
    outgroup = group_labels(groups, "RYR") & all_leaves
    core, _audit = corrected_membership(groups, tree)

    # Build every constraint file first, then run the three searches
    # together — they are independent, and three at CONSTRAINT_THREADS
    # fits the machine that one at THREADS does not fill.
    jobs: list[tuple[str, Path, subprocess.Popen | None]] = []
    for name, nwk in au_constraints(core, outgroup, all_leaves).items():
        cpath = AU_DIR / f"{name}.constraint.nwk"
        cpath.write_text(nwk + "\n")
        prefix = AU_DIR / name
        tf = Path(str(prefix) + ".treefile")
        # Resume on the .iqtree report, not the .treefile — the same
        # trap `run_ml` documents: IQ-TREE writes a treefile at
        # checkpoints, so a constrained search killed part-way leaves a
        # treefile of an unfinished topology and the AU test would be
        # run on it without a symptom.
        done = Path(str(prefix) + ".iqtree")
        proc = None
        if not done.exists():
            claim(prefix)
            cmd = ["iqtree2", "-s", str(INPUT), "--prefix", str(prefix),
                   "-m", model, "-g", str(cpath),
                   "-T", str(CONSTRAINT_THREADS), "-seed", str(SEED)]
            print("running:", " ".join(cmd), flush=True)
            proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL)
        jobs.append((name, tf, proc))

    trees = [TREEFILE.read_text().strip()]
    names = ["ML"]
    t0 = time.time()
    for name, tf, proc in jobs:
        if proc is not None and proc.wait():
            return proc.returncode
        if not tf.exists() or not Path(str(AU_DIR / name) + ".iqtree").exists():
            raise SystemExit(f"constrained search {name} did not finish")
        trees.append(tf.read_text().strip())
        names.append(name)
    print(f"three constrained searches done in {time.time() - t0:.0f} s")

    (AU_DIR / "candidates.nwk").write_text("\n".join(trees) + "\n")
    (AU_DIR / "candidates.names").write_text("\n".join(names) + "\n")
    claim(AU_DIR / "autest")
    cmd = ["iqtree2", "-s", str(INPUT), "--prefix", str(AU_DIR / "autest"),
           "-m", model, "-z", str(AU_DIR / "candidates.nwk"),
           "-n", "0", "-zb", str(AU_REPLICATES), "-au",
           "-T", str(THREADS), "-seed", str(SEED)]
    print("running:", " ".join(cmd))
    rc = subprocess.run(cmd).returncode
    if rc == 0:
        write_au_table(names)
        _record_au_stats(model)
        _compare_best_tree()
        write_live("S7", _steps())
    return rc


def _compare_best_tree() -> None:
    """If a constrained search beat the free one, say which claims move.

    A constraint that costs negative likelihood means the reported tree
    is not the global optimum (D34). Which tree that is comes out of
    `au_test.tsv` rather than being named here, so a rerun on different
    data compares whatever actually won.
    """
    import csv
    tsv = PHYLO_DIR / "au_test.tsv"
    if not tsv.exists():
        return
    with open(tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    try:
        best = min(rows, key=lambda r: float(r["deltaL"]))
        ml = next(r for r in rows if r["tree"] == "ML")
    except (ValueError, KeyError, StopIteration):
        return
    if best["tree"] == "ML" or float(ml["deltaL"]) <= 0:
        print("the unconstrained tree is the best of the set — "
              "no alternative-tree comparison needed")
        return
    tree = AU_DIR / f"{best['tree']}.treefile"
    print(f"note: {best['tree']} beat the unconstrained search by "
          f"{float(ml['deltaL']):.1f} logL — comparing the claim set "
          f"against it (D34)")
    subprocess.run([sys.executable, str(Path(__file__).with_name(
        "s7_bnni.py")), "--tree", str(tree), "--label",
        f"{best['tree']} (best tree found)", "--out",
        "best_tree_comparison.tsv"])


def _record_au_stats(model: str) -> None:
    """The AU run's own settings, beside the search's (D24)."""
    if not STATS.exists():
        return
    d = json.loads(STATS.read_text())
    d["au_test"] = {
        "model": model,
        "constrained_search_threads": CONSTRAINT_THREADS,
        "constrained_searches_run_in_parallel": 3,
        "au_replicates": AU_REPLICATES,
        "seed": SEED,
        "hypotheses": "H1_12 / H2_13 / H3_23 — see scripts/s7_constraints.py",
    }
    STATS.write_text(json.dumps(d, indent=1) + "\n")


def write_au_table(names: list[str]) -> None:
    """Parse IQ-TREE's `-au` block into `au_test.tsv` (D13: the report
    is rendered from this table, never from the .iqtree file)."""
    from scripts.s7_constraints import HYPOTHESES
    text = Path(str(AU_DIR / "autest") + ".iqtree").read_text()
    block = text.split("USER TREES", 1)[-1]
    rows = []
    for line in block.splitlines():
        m = re.match(r"\s*(\d+)\s+(-?[\d.]+)\s+([\d.e+-]+)\s+(.*)", line)
        if not m:
            continue
        idx = int(m.group(1))
        nums = re.findall(r"[\d.]+(?:e-?\d+)?", m.group(4))
        # columns: bp-RELL p-KH p-SH c-ELW p-AU (each followed by +/-)
        rows.append({
            "idx": idx,
            "tree": names[idx - 1] if idx <= len(names) else f"T{idx}",
            "hypothesis": HYPOTHESES.get(
                names[idx - 1] if idx <= len(names) else "", ""),
            "logL": m.group(2),
            "deltaL": m.group(3),
            "bp_RELL": nums[0] if len(nums) > 0 else "",
            "p_KH": nums[1] if len(nums) > 1 else "",
            "p_SH": nums[2] if len(nums) > 2 else "",
            "c_ELW": nums[3] if len(nums) > 3 else "",
            "p_AU": nums[4] if len(nums) > 4 else "",
        })
    cols = ["idx", "tree", "hypothesis", "logL", "deltaL", "bp_RELL",
            "p_KH", "p_SH", "c_ELW", "p_AU"]
    out = AU_DIR.parent / "au_test.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {out} ({len(rows)} trees)")


def main() -> int:
    sub = sys.argv[1] if len(sys.argv) > 1 else "ml"
    if sub == "prep":
        prep_input()
        return 0
    if sub == "ml":
        return run_ml()
    if sub == "bnni":
        return run_bnni()
    if sub == "au":
        return run_au()
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
