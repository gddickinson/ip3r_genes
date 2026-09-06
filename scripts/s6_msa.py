"""S6 step 2 — MAFFT L-INS-i + trimAl over the representative set.

Reads  results/msa_v2/representatives.{fasta,tsv}   (s6_select_reps.py)
Writes results/msa_v2/aln.fasta         full L-INS-i alignment
       results/msa_v2/trimmed.fasta     trimAl -automated1
       results/msa_v2/align_stats.json  commands, versions, runtimes, SHA-256
       results/msa_v2/identity_covered.tsv    fragment-aware
       results/msa_v2/identity_classic.tsv    gaps count as mismatch
       results/msa_v2/conservation.tsv        per-column, trimmed MSA
       results/msa_v2/coverage.tsv            per-sequence, trimmed MSA

**Single-threaded, by decision (D24).** MAFFT L-INS-i at `--thread -1` is
not reproducible on this machine: S3 aligned the same 22 RyR seeds twice and
got 8,510 and 8,468 columns, because the iterative refinement stage combines
partial results in whatever order the threads finish. D24 says anything a
committed artefact is built from is aligned single-threaded, and names S6
explicitly. The whole downstream chain — the tree, the AU test, the codon
alignment, the conservation profile — is built from this file, so it is
pinned to `--thread 1` and the SHA-256 of every input and output is recorded
in `align_stats.json`, so a rebuild that drifts is visible in the committed
data rather than only in a column count.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s6_msa.py
      ... --align-only     just MAFFT (hours; run it in the background)
      ... --stats-only     everything after MAFFT, from a finished aln
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from s6_lib import MSA_DIR, read_fasta, read_tsv, write_live  # noqa: E402

REPS_FASTA = MSA_DIR / "representatives.fasta"
REPS_TSV = MSA_DIR / "representatives.tsv"
ALN = MSA_DIR / "aln.fasta"
TRIMMED = MSA_DIR / "trimmed.fasta"
STATS = MSA_DIR / "align_stats.json"

#: The brief allows `--auto` "if size forces it", and the choice must be
#: recorded. L-INS-i is O(N^2) full pairwise DP before the progressive
#: stage; at this family's ~2,800 aa subunit the practical ceiling on one
#: thread is around here; the set S6 selects is inside it.
LINSI_MAX = 200


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(cmd: list[str], stream: str = "stderr") -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
        text = p.stderr if stream == "stderr" else p.stdout
        return next((l.strip() for l in text.splitlines() if l.strip()), "?")
    except OSError:
        return "not found"


def run_logged(cmd: list[str], stdout_path: Path | None,
               log_path: Path) -> float:
    t0 = time.time()
    with open(log_path, "w") as log:
        if stdout_path:
            with open(stdout_path, "w") as out:
                subprocess.run(cmd, stdout=out, stderr=log, check=True)
        else:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                           check=True)
    return time.time() - t0


def gap_pct(seqs: dict[str, str]) -> float:
    total = sum(len(s) for s in seqs.values())
    gaps = sum(s.count("-") for s in seqs.values())
    return 100.0 * gaps / total if total else 0.0


def ragged(seqs: dict[str, str]) -> list[str]:
    """S1's lesson: a silent MAFFT failure degrades to a star alignment."""
    lens = {len(s) for s in seqs.values()}
    return [] if len(lens) <= 1 else sorted(lens)


def write_matrix(path: Path, labels: list[str], m) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label"] + labels)
        for lab, row in zip(labels, m):
            w.writerow([lab] + [f"{v:.4f}" for v in row])


# ------------------------------------------------------------------ stages

def align(stats: dict) -> dict:
    reps = read_fasta(REPS_FASTA)
    n = len(reps)
    use_linsi = n <= LINSI_MAX
    cmd = (["mafft", "--localpair", "--maxiterate", "1000"] if use_linsi
           else ["mafft", "--auto"])
    cmd += ["--thread", "1", "--anysymbol", str(REPS_FASTA)]
    print(f"MAFFT {'L-INS-i' if use_linsi else '--auto'} on {n} sequences, "
          f"{sum(len(s) for s in reps.values()):,} residues, single-threaded "
          f"(D24). This takes hours.", flush=True)
    t = run_logged(cmd, ALN, MSA_DIR / "mafft.log")
    aln = read_fasta(ALN)
    if len(aln) != n:
        raise SystemExit(f"MAFFT returned {len(aln)} of {n} sequences")
    # The alignment must be of *this* input, not of whatever was on disk
    # when the run started. A long MAFFT and an edited selector overlap
    # easily, and the failure is silent: `align_stats.json` would record
    # the SHA-256 of a file the alignment was not built from, which is
    # exactly the drift D24 exists to make visible.
    drift = [k for k, v in aln.items()
             if v.replace("-", "").upper() != reps.get(k, "").upper()]
    if drift or set(aln) != set(reps):
        raise SystemExit(
            f"aln.fasta is not an alignment of the current "
            f"representatives.fasta: {len(drift)} sequences differ, "
            f"{len(set(aln) ^ set(reps))} labels differ. Re-run MAFFT.")
    bad = ragged(aln)
    if bad:
        raise SystemExit(f"alignment is ragged, lengths {bad} — MAFFT failed "
                         f"silently (see mafft.log)")
    stats["mafft"] = {
        "aligner": "L-INS-i" if use_linsi else "--auto",
        "why": (f"n={n} <= {LINSI_MAX}, the brief's L-INS-i applies"
                if use_linsi else
                f"n={n} > {LINSI_MAX}, size forced --auto"),
        "command": " ".join(cmd), "version": tool_version(["mafft", "--version"]),
        "threads": 1, "thread_note": "D24 — pinned; --thread -1 is not "
                                     "reproducible for this pipeline",
        "runtime_s": round(t, 1), "n_sequences": n,
        "columns": len(next(iter(aln.values()))),
        "gap_pct": round(gap_pct(aln), 2),
        "sha256_in": sha256(REPS_FASTA), "sha256_out": sha256(ALN),
    }
    print(f"  aligned: {stats['mafft']['columns']} columns, "
          f"{stats['mafft']['gap_pct']}% gaps, {t/60:.1f} min", flush=True)
    return stats


def trim(stats: dict) -> dict:
    cmd = ["trimal", "-in", str(ALN), "-out", str(TRIMMED), "-automated1"]
    t = run_logged(cmd, None, MSA_DIR / "trimal.log")
    trimmed = {k.split()[0]: v for k, v in read_fasta(TRIMMED).items()}
    bad = ragged(trimmed)
    if bad:
        raise SystemExit(f"trimmed alignment is ragged, lengths {bad}")
    aln_cols = stats.get("mafft", {}).get(
        "columns") or len(next(iter(read_fasta(ALN).values())))
    stats["trimal"] = {
        "command": " ".join(cmd),
        "version": tool_version(["trimal", "--version"], "stdout"),
        "runtime_s": round(t, 1), "columns_in": aln_cols,
        "columns_kept": len(next(iter(trimmed.values()))),
        "pct_kept": round(100 * len(next(iter(trimmed.values()))) / aln_cols, 1),
        "gap_pct": round(gap_pct(trimmed), 2), "sha256_out": sha256(TRIMMED),
    }
    print(f"  trimAl: {stats['trimal']['columns_kept']} of {aln_cols} columns "
          f"kept ({stats['trimal']['pct_kept']}%)", flush=True)
    return stats


def matrices(stats: dict) -> dict:
    from src.analysis.alignment import AlignmentRow
    from src.analysis.distance import identity_matrix
    from src.analysis.evolution import conservation_per_column
    import s6_rep_spec as spec

    trimmed = {k.split()[0]: v for k, v in read_fasta(TRIMMED).items()}
    meta = {r["label"]: r for r in read_tsv(REPS_TSV)}
    order = {g: i for i, g in enumerate(spec.GROUP_ORDER)}
    labels = sorted(trimmed, key=lambda l: (
        order.get(meta.get(l, {}).get("group", ""), 99),
        meta.get(l, {}).get("species", l)))
    rows = [AlignmentRow(label=l, aligned=trimmed[l]) for l in labels]

    cov = identity_matrix(rows, covered_only=True)
    cls = identity_matrix(rows)
    write_matrix(MSA_DIR / "identity_covered.tsv", labels, cov.identity)
    write_matrix(MSA_DIR / "identity_classic.tsv", labels, cls.identity)

    cons = conservation_per_column(rows)
    with open(MSA_DIR / "conservation.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["column", "conservation"])
        for i, c in enumerate(cons, 1):
            w.writerow([i, f"{c:.4f}"])

    L = len(next(iter(trimmed.values())))
    with open(MSA_DIR / "coverage.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "group", "species", "coverage", "ungapped_aa"])
        for l in labels:
            ung = L - trimmed[l].count("-")
            w.writerow([l, meta.get(l, {}).get("group", ""),
                        meta.get(l, {}).get("species", ""),
                        f"{ung / L:.4f}", ung])
    stats["matrices"] = {
        "n_labels": len(labels), "trimmed_columns": L,
        "mean_conservation": round(sum(cons) / len(cons), 4),
        "columns_over_0.9": sum(1 for c in cons if c >= 0.9),
    }
    stats["sites"] = site_classes(trimmed, labels)
    print(f"  matrices: {len(labels)} labels, mean conservation "
          f"{stats['matrices']['mean_conservation']}", flush=True)
    return stats


def site_classes(trimmed: dict[str, str], labels: list[str]) -> dict:
    """Constant / variable-uninformative / parsimony-informative columns.

    The number S7 actually has to work with. trimAl `-automated1` is a
    heuristic, and a trimming that kept plenty of columns but cut the
    informative ones would look fine in a percentage and starve the tree;
    this counts what is left. A column is parsimony-informative when at
    least two different residues each appear at least twice — the standard
    definition, and the one IQ-TREE reports, so the two can be compared.
    """
    from collections import Counter
    L = len(next(iter(trimmed.values())))
    const = uninf = inf = allgap = 0
    for k in range(L):
        col = Counter(trimmed[l][k] for l in labels
                      if trimmed[l][k] not in "-Xx")
        if not col:
            allgap += 1
        elif len(col) == 1:
            const += 1
        elif sum(1 for c in col.values() if c >= 2) >= 2:
            inf += 1
        else:
            uninf += 1
    return {"constant": const, "variable_uninformative": uninf,
            "parsimony_informative": inf, "all_gap_or_X": allgap,
            "pct_informative": round(100 * inf / L, 1)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--align-only", action="store_true")
    ap.add_argument("--stats-only", action="store_true")
    args = ap.parse_args()
    if not REPS_FASTA.exists():
        raise SystemExit("run scripts/s6_select_reps.py first")
    stats = json.loads(STATS.read_text()) if STATS.exists() else {}

    if not args.stats_only:
        write_live("S6", [("select representatives", True),
                          ("resolve sequences", True), ("MAFFT L-INS-i", False),
                          ("trimAl", False), ("matrices + figures", False)])
        stats = align(stats)
        STATS.write_text(json.dumps(stats, indent=1))
        if args.align_only:
            return 0
    write_live("S6", [("select representatives", True),
                      ("resolve sequences", True), ("MAFFT L-INS-i", True),
                      ("trimAl", False), ("matrices + figures", False)])
    stats = trim(stats)
    STATS.write_text(json.dumps(stats, indent=1))
    write_live("S6", [("select representatives", True),
                      ("resolve sequences", True), ("MAFFT L-INS-i", True),
                      ("trimAl", True), ("matrices + figures", False)])
    stats = matrices(stats)
    STATS.write_text(json.dumps(stats, indent=1))
    write_live("S6", [("select representatives", True),
                      ("resolve sequences", True), ("MAFFT L-INS-i", True),
                      ("trimAl", True), ("matrices + figures", True)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
