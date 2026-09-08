"""S12's cross-mapping control — measured, not asserted.

The reference is a **closed set**: a read from a gene not in it cannot be
assigned away, so if the three ITPR paralogs were similar enough at the
nucleotide level for a 100 nt read to cross between them, reads from one
would inflate another and the whole task would be measuring the wrong gene.
The PIEZO project's S12 argued this away in a caveat ("well under 70 %
nucleotide identity, far below what hisat2 end-to-end will align across").
That is an argument, and the ITPR paralogs are more similar to each other
than the PIEZO ones are, so here it is measured instead.

The test is direct: **tile every reference sequence with synthetic reads of
the library's own read length and map them back with the same aligner
settings**.  A read from `ITPR2|x` must be assigned to `ITPR2|x`.  Anything
else is this reference's cross-mapping rate, and it is written to a table
whether it is zero or not.

Reads are tiled, not sampled, so the test is exhaustive over positions
rather than a random subset — and deterministic, so a rerun that changes
the number is a change in the reference or the aligner, never in a seed.

Output (committed):
    results/expression/crossmap_control.tsv

Run:  python scripts/s12_crossmap.py [--read-len 100] [--step 10]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s12_lib import DATA, OUT_DIR, live, read_fasta, tool_bin, write_tsv  # noqa: E402
from s12_quantify import MIN_MAPQ  # noqa: E402

COLS = ["species", "source_seq", "role", "reads_simulated", "aligned",
        "assigned_self", "assigned_other", "unaligned", "below_mapq",
        "worst_other_seq", "worst_other_reads", "self_rate", "cross_rate"]


def tile(seq: str, read_len: int, step: int) -> list[str]:
    return [seq[i:i + read_len]
            for i in range(0, max(1, len(seq) - read_len + 1), step)]


def run_species(short: str, read_len: int, step: int, threads: int,
                log=print) -> list[dict]:
    fa = DATA / "refs" / f"{short}.fasta"
    idx = DATA / "index" / short
    if not fa.exists() or not Path(str(idx) + ".1.ht2").exists():
        log(f"  !! no reference or index for {short}")
        return []
    refs = read_fasta(fa)
    rows: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "reads.fa"
        origin: dict[str, str] = {}
        with open(q, "w") as fh:
            for name, seq in refs.items():
                if name.startswith("decoy_"):
                    continue
                for i, r in enumerate(tile(seq, read_len, step)):
                    rid = f"{len(origin)}"
                    origin[rid] = name
                    fh.write(f">{rid}\n{r}\n")
        if not origin:
            return []
        p = subprocess.run(
            [tool_bin("hisat2"), "-x", str(idx), "-f", "-U", str(q),
             "-p", str(threads), "--no-spliced-alignment", "--no-unal"],
            capture_output=True, text=True)
        if p.returncode != 0:
            log(f"  !! hisat2 failed for {short}: {p.stderr[:200]}")
            return []

    per: dict[str, dict] = {}
    for name in {v for v in origin.values()}:
        per[name] = {"aligned": 0, "self": 0, "other": 0, "lowq": 0,
                     "others": {}}
    for line in p.stdout.split("\n"):
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        flag = int(f[1])
        if flag & 0x100 or flag & 0x800 or flag & 0x4:
            continue
        src = origin.get(f[0])
        if src is None:
            continue
        rec = per[src]
        if int(f[4]) < MIN_MAPQ:
            rec["lowq"] += 1
            continue
        rec["aligned"] += 1
        if f[2] == src:
            rec["self"] += 1
        else:
            rec["other"] += 1
            rec["others"][f[2]] = rec["others"].get(f[2], 0) + 1

    totals: dict[str, int] = {}
    for rid, name in origin.items():
        totals[name] = totals.get(name, 0) + 1
    roles = _roles(short)
    for name, rec in sorted(per.items()):
        n = totals[name]
        worst = max(rec["others"].items(), key=lambda kv: kv[1],
                    default=("", 0))
        rows.append({
            "species": short, "source_seq": name,
            "role": roles.get(name, ""), "reads_simulated": n,
            "aligned": rec["aligned"], "assigned_self": rec["self"],
            "assigned_other": rec["other"],
            "unaligned": n - rec["aligned"] - rec["lowq"],
            "below_mapq": rec["lowq"],
            "worst_other_seq": worst[0], "worst_other_reads": worst[1],
            "self_rate": round(rec["self"] / n, 5) if n else 0,
            "cross_rate": round(rec["other"] / n, 5) if n else 0})
    return rows


def _roles(short: str) -> dict[str, str]:
    import csv
    path = OUT_DIR / "reference_table.tsv"
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] == short:
                out[r["seq"]] = r["role"]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="")
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()

    import s12_panel
    panel, _ = s12_panel.build_panel()
    if args.species:
        want = set(args.species.split(","))
        panel = [s for s in panel if s.short in want]

    rows: list[dict] = []
    for i, sp in enumerate(panel, 1):
        print(f"[crossmap] {sp.name}", flush=True)
        live("crossmap", i - 1, len(panel), sp.name)
        r = run_species(sp.short, args.read_len, args.step, args.threads)
        for x in r:
            print(f"    {x['source_seq']:24s} {x['assigned_self']:>6}/"
                  f"{x['reads_simulated']:<6} self, "
                  f"{x['assigned_other']} elsewhere"
                  + (f" (worst: {x['worst_other_seq']} "
                     f"x{x['worst_other_reads']})"
                     if x["assigned_other"] else ""))
        rows += r
    write_tsv(OUT_DIR / "crossmap_control.tsv", COLS, rows)
    live("crossmap", len(panel), len(panel), "done")
    worst = max((r["cross_rate"] for r in rows), default=0)
    print(f"\nwrote crossmap_control.tsv ({len(rows)} rows); "
          f"worst cross-mapping rate {worst:.5f}")


if __name__ == "__main__":
    main()
