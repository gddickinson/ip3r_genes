"""S12 step 3 — stream each SRA run and count reads per reference and junction.

Per run, `fastq-dump` streams the first N spots straight into `hisat2`
against the small per-species reference; the SAM is parsed in-process and
only counts are kept.  Nothing but counts touches the disk, so the whole
panel costs no meaningful storage.

Three things are counted for every reference sequence:

  * `reads` — primary alignments at or above the MAPQ floor;
  * `junction_reads` — alignments whose *contiguous* aligned block covers
    at least `--anchor` nt on **both** sides of an exon junction.  Genomic
    DNA carryover cannot produce these, so they are the evidence that a
    transcript exists rather than a locus;
  * `covered_bases` — the union of reference positions covered, so reads
    piled on one repeat cannot pass for expression.

And, separately, **each junction is counted on its own**.  That is the
measurement S10 handed here: the omission cases carry 55 junctions apiece
that no annotated gene model spans, and "this gene is transcribed" is a
weaker claim than "these particular junctions, which the annotation does
not have, are spliced in a real library".

`--no-spliced-alignment` is deliberate.  The reference is already spliced —
it is CDS, not genome — so a read crossing a junction is *contiguous* here.
Leaving hisat2's spliced mode on would let it open a gap inside the CDS and
call an intron that does not exist, and those alignments would then be
counted as spanning a junction they had actually skipped.

Resumable: a run whose counts file exists is skipped.

Outputs (bulk, under the data root):
    <data_root>/expression/counts/<species>.<run>.tsv
    <data_root>/expression/counts/<species>.<run>.junctions.tsv
    <data_root>/expression/counts/<species>.<run>.coverage.tsv

Run:  python scripts/s12_quantify.py [--species ...] [--spots 4000000]
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s12_lib import DATA, OUT_DIR, live, tool_bin  # noqa: E402

COUNTS = DATA / "counts"
MIN_MAPQ = 10            # hisat2 gives ~1 to multi-mappers; 10 keeps unique
DEFAULT_SPOTS = 4_000_000
DEFAULT_ANCHOR = 8
COVERAGE_BINS = 200


def load_junctions() -> dict[tuple[str, str], list[dict]]:
    path = OUT_DIR / "junctions.tsv"
    out: dict[tuple[str, str], list[dict]] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out.setdefault((row["species"], row["seq"]), []).append({
                "index": int(row["junction_index"]),
                "offset": int(row["cds_offset"]),
                "annotated": int(row["annotated"]),
                "class": row["class"]})
    return out


def ref_meta(species: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    path = OUT_DIR / "reference_table.tsv"
    if not path.exists():
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["species"] == species:
                out[row["seq"]] = {"length": int(row["length"] or 0),
                                   "role": row["role"]}
    return out


def aligned_blocks(pos: int, cigar: str) -> list[tuple[int, int]]:
    """Reference intervals (1-based, inclusive) covered by the alignment."""
    blocks: list[tuple[int, int]] = []
    ref, num = pos, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num or 0)
        num = ""
        if ch in "M=X":
            blocks.append((ref, ref + n - 1))
            ref += n
        elif ch in "DN":
            ref += n
        # I, S, H, P consume no reference
    return blocks


def spanned(blocks: list[tuple[int, int]], offsets: list[int], anchor: int
            ) -> list[int]:
    """Which junction offsets a *contiguous* block reads through.

    The anchor is required on both sides, and on the same block.  Requiring
    it on one side only admits an alignment that has run a few bases past
    the junction and stopped, which is not evidence of splicing — S10 hit
    exactly that and its junction rule grew its second half from it.
    """
    hit: list[int] = []
    for lo, hi in blocks:
        for j in offsets:
            if lo <= j - anchor + 1 and hi >= j + anchor:
                hit.append(j)
    return hit


def quantify_run(species: str, run: str, spots: int, threads: int,
                 anchor: int, junctions: dict, log=print,
                 force: bool = False) -> dict | None:
    idx = DATA / "index" / species
    out_path = COUNTS / f"{species}.{run}.tsv"
    if out_path.exists() and not force:
        return {"run": run, "species": species, "skipped": True}
    if not Path(str(idx) + ".1.ht2").exists():
        log(f"  !! no hisat2 index for {species} — run s12_refs.py first")
        return None

    COUNTS.mkdir(parents=True, exist_ok=True)
    meta = ref_meta(species)
    t0 = time.time()

    fq = subprocess.Popen(
        [tool_bin("fastq-dump"), "-X", str(spots), "--split-spot",
         "--skip-technical", "-Z", run],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    # hisat2's summary carries the library size, which is the denominator of
    # any per-million figure.  It has to be the reads actually streamed, not
    # the number requested: a run can hold fewer spots than --spots asks for
    # and normalising by the request would deflate every rate in it.
    summary = COUNTS / f".{species}.{run}.summary"
    ht = subprocess.Popen(
        [tool_bin("hisat2"), "-x", str(idx), "-U", "-", "-p", str(threads),
         "--no-spliced-alignment", "--no-unal",
         "--summary-file", str(summary)],
        stdin=fq.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=1)
    if fq.stdout:
        fq.stdout.close()

    counts: dict[str, int] = {}
    junc_reads: dict[str, int] = {}
    per_junction: dict[tuple[str, int], int] = {}
    covered: dict[str, set] = {}
    n_aligned = n_lowq = 0

    assert ht.stdout is not None
    for line in ht.stdout:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        flag = int(f[1])
        if flag & 0x100 or flag & 0x800 or flag & 0x4:
            continue
        if int(f[4]) < MIN_MAPQ:
            n_lowq += 1
            continue
        ref, pos, cigar = f[2], int(f[3]), f[5]
        n_aligned += 1
        counts[ref] = counts.get(ref, 0) + 1
        blocks = aligned_blocks(pos, cigar)
        cov = covered.setdefault(ref, set())
        for lo, hi in blocks:
            cov.update(range(lo, hi + 1))
        js = junctions.get((species, ref))
        if not js:
            continue
        hits = spanned(blocks, [j["offset"] for j in js], anchor)
        if hits:
            junc_reads[ref] = junc_reads.get(ref, 0) + 1
            for off in set(hits):
                per_junction[(ref, off)] = per_junction.get((ref, off), 0) + 1

    ht.wait()
    fq.wait()
    if ht.returncode not in (0, None):
        err = (ht.stderr.read() if ht.stderr else "")[:300]
        log(f"  !! hisat2 exit {ht.returncode} for {run}: {err}")
        summary.unlink(missing_ok=True)
        return None

    n_reads = _library_size(summary)
    summary.unlink(missing_ok=True)
    if n_reads == 0:
        log(f"  !! {run}: fastq-dump streamed no reads (run withdrawn, or "
            f"the SRA toolkit is not configured) — not written")
        return None

    _write_coverage(COUNTS / f"{species}.{run}.coverage.tsv", covered, meta)
    _write_junctions(COUNTS / f"{species}.{run}.junctions.tsv", species,
                     junctions, per_junction)
    with open(out_path, "w") as fh:
        fh.write(f"# library_reads\t{n_reads}\n")
        fh.write("seq\treads\tjunction_reads\tcovered_bases\tref_len\n")
        for ref in sorted(set(counts) | set(meta)):
            fh.write(f"{ref}\t{counts.get(ref, 0)}\t{junc_reads.get(ref, 0)}"
                     f"\t{len(covered.get(ref, ()))}"
                     f"\t{meta.get(ref, {}).get('length', 0)}\n")
    log(f"  {run}: {n_reads:,} reads in, {n_aligned} aligned "
        f"(MAPQ>={MIN_MAPQ}), {sum(junc_reads.values())} junction-spanning "
        f"[{time.time() - t0:.0f}s]")
    return {"run": run, "species": species, "aligned": n_aligned,
            "low_mapq": n_lowq, "library_reads": n_reads,
            "spots_requested": spots, "seconds": round(time.time() - t0),
            "skipped": False}


def _library_size(summary: Path) -> int:
    if not summary.exists():
        return 0
    for line in summary.read_text().split("\n"):
        s = line.strip()
        if s.endswith("reads; of these:"):
            try:
                return int(s.split()[0])
            except ValueError:
                return 0
    return 0


def _write_coverage(path: Path, covered: dict, meta: dict) -> None:
    """Binned coverage for the real references — where along the ORF the
    reads fall, which is the direct form of the S10 question."""
    with open(path, "w") as fh:
        fh.write("seq\tn_bins\tbin_width\tcoverage\n")
        for ref, positions in sorted(covered.items()):
            if meta.get(ref, {}).get("role") in ("decoy", "housekeeping"):
                continue
            length = meta.get(ref, {}).get("length", 0)
            if not length:
                continue
            width = max(1, length // COVERAGE_BINS)
            bins = [0] * COVERAGE_BINS
            for p in positions:
                bins[min(COVERAGE_BINS - 1, (p - 1) // width)] += 1
            fh.write(f"{ref}\t{COVERAGE_BINS}\t{width}\t"
                     f"{','.join(str(b) for b in bins)}\n")


def _write_junctions(path: Path, species: str, junctions: dict,
                     per_junction: dict) -> None:
    """Every junction of every reference, whether or not a read crossed it.

    Junctions with no reads are written as zeros rather than omitted: the
    denominator — how many of a gene's junctions were *asked* — is half the
    result, and a table of only the ones that fired cannot supply it.
    """
    with open(path, "w") as fh:
        fh.write("seq\tjunction_index\tcds_offset\tannotated\tclass\treads\n")
        for (sp, ref), js in sorted(junctions.items()):
            if sp != species:
                continue
            for j in js:
                fh.write(f"{ref}\t{j['index']}\t{j['offset']}"
                         f"\t{j['annotated']}\t{j['class']}"
                         f"\t{per_junction.get((ref, j['offset']), 0)}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="")
    ap.add_argument("--spots", type=int, default=DEFAULT_SPOTS)
    ap.add_argument("--anchor", type=int, default=DEFAULT_ANCHOR)
    ap.add_argument("--threads", type=int,
                    default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--runs", default="")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    runs_path = OUT_DIR / "runs_selected.tsv"
    if not runs_path.exists():
        raise SystemExit("run scripts/s12_runs.py first")
    with open(runs_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if args.species:
        want = set(args.species.split(","))
        rows = [r for r in rows if r["species"] in want]
    if args.runs:
        want_runs = set(args.runs.split(","))
        rows = [r for r in rows if r["run"] in want_runs]
    if args.limit:
        rows = rows[:args.limit]

    junctions = load_junctions()
    ok = failed = skipped = 0
    for i, r in enumerate(rows, 1):
        print(f"[{i}/{len(rows)}] {r['species']} {r['run']} ({r['tissue']})",
              flush=True)
        live("quantify", i - 1, len(rows),
             f"{r['species']} {r['run']} ({r['tissue']})")
        m = quantify_run(r["species"], r["run"], args.spots, args.threads,
                         args.anchor, junctions, force=args.force)
        if m is None:
            failed += 1
        elif m.get("skipped"):
            skipped += 1
        else:
            ok += 1
    live("quantify", len(rows), len(rows), "done")
    print(f"\n{ok} quantified, {skipped} already done, {failed} failed")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
