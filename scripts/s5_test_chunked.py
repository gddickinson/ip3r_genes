"""s5_test_chunked.py — prove the chunked miniprot path returns the same loci.

The two largest assemblies in the manifest (*Protopterus annectens* 40 Gbp,
*Lissotriton helveticus* 23 Gbp) cannot be indexed by miniprot in one pass, so
they take `run_miniprot_chunked`: split into whole-contig chunks, align each,
concatenate the GFFs. That path is ported and, until this test, had never been
executed in this project — and it will run unattended on the two genomes least
likely to be checked by hand.

Two specific ways it can be silently wrong:

  **Coordinates.** The design claims chunks hold complete contigs so every
  reported coordinate stays valid without translation. If a contig were ever
  split, coordinates in later chunks would be offset and the loci would land
  in the wrong place — with nothing in the output to say so.

  **Alignment IDs.** miniprot numbers alignments from `MP000001` *per run*, so
  a concatenated multi-chunk GFF repeats every ID. CDS rows are bound to their
  parent by ID, so a naive parse would attach one alignment's exons to a
  different alignment's mRNA — silently inflating coverage and changing calls.

So the test is equivalence, not smoke: run a genome we already have both ways
and require the same loci, the same baits, the same coverage and the same cell
statuses, with boundaries matching to `BOUNDARY_TOLERANCE_BP`. It also measures
miniprot's peak RSS per Gbp — the number that says whether a production chunk
fits in this machine's memory before *Protopterus* finds out.

What it found on first run, on *Takifugu rubripes*:

  * the two paths agree on all 10 loci and all 4 cell statuses, with a worst
    boundary difference of 9 bp;
  * the ID collision is **real** — 167 mRNA rows carry only 60 distinct raw
    IDs across 6 chunks — so the guard in `parse_miniprot_gff` is load-bearing
    rather than defensive;
  * miniprot wants 6.3-8.0 GB of RSS per Gbp of reference, which put the
    ported 4 Gbp chunk at ~25 GB on a 34 GB machine. `CHUNK_BP` is now 2.5
    Gbp (~20 GB), because the giants run unattended.

Usage:
  python scripts/s5_test_chunked.py                      # smallest swept genome
  python scripts/s5_test_chunked.py --accession GCF_...  # a specific one
  python scripts/s5_test_chunked.py --chunk-mbp 60       # force more chunks
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

for extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(extra).is_dir() and extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root              # noqa: E402
from fetch_genomes import fna_path, read_manifest          # noqa: E402
import s5_sweep_lib as lib                                 # noqa: E402
import s5_classify as clf                                  # noqa: E402
import s5_genome_io as gio                                 # noqa: E402

ALL_CELLS = (*lib.CLASSES, lib.CONTROL_CLASS)


def peak_rss_gb() -> float:
    """Peak RSS of children **so far** — a high-water mark, not a gauge.

    The first version of this test reported `after - before` per phase, which
    is always <= 0 for the second phase however much memory it used, and duly
    reported 0.00 GB for the chunked run. It is a running maximum: only the
    absolute value after each phase means anything, and only the phase that
    set a new high can be attributed the increase.
    """
    raw = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return raw / 1e9 if sys.platform == "darwin" else raw * 1024 / 1e9


#: Locus boundaries may differ slightly between the two paths without either
#: being wrong: miniprot extends a terminal exon using flanking reference
#: context, and at a chunk edge there is less of it. Measured on Takifugu the
#: differences were 3 and 9 bp on 2 of 10 loci, with identical bait, coverage,
#: family and cell status. A genuine coordinate-offset bug — a contig split
#: across chunks — shifts loci by thousands to millions of bp, so a tolerance
#: this tight still catches the failure this test exists for.
BOUNDARY_TOLERANCE_BP = 50


def loci_key(loci: list) -> dict:
    """Identity of each locus, excluding the boundaries that may jitter."""
    return {(L.contig, L.strand, L.best.bait, L.family): (L.start, L.end,
                                                          round(L.best.coverage, 4))
            for L in loci}


def compare_loci(a: dict, b: dict) -> tuple[list[str], int]:
    """(failures, worst boundary delta seen)."""
    failures: list[str] = []
    worst = 0
    only_a, only_b = set(a) - set(b), set(b) - set(a)
    if only_a or only_b:
        failures.append(f"{len(only_a)} loci only unchunked, {len(only_b)} "
                        f"only chunked (same contig/strand/bait/family key)")
        for k in list(only_a)[:3]:
            failures.append(f"      unchunked-only {k} {a[k]}")
        for k in list(only_b)[:3]:
            failures.append(f"      chunked-only   {k} {b[k]}")
    for k in set(a) & set(b):
        (s1, e1, c1), (s2, e2, c2) = a[k], b[k]
        worst = max(worst, abs(s1 - s2), abs(e1 - e2))
        if abs(s1 - s2) > BOUNDARY_TOLERANCE_BP or \
           abs(e1 - e2) > BOUNDARY_TOLERANCE_BP:
            failures.append(
                f"      {k[0]}{k[1]} moved {abs(s1 - s2)}/{abs(e1 - e2)} bp "
                f"(> {BOUNDARY_TOLERANCE_BP}): {s1}-{e1} vs {s2}-{e2}")
        if c1 != c2:
            failures.append(f"      {k[0]}{k[1]} coverage {c1} vs {c2}")
    return failures, worst


def cell_signature(loci: list, fna: Path, idx: dict, seqlens: dict) -> dict:
    def fetch(c, s, e):
        return gio.fetch_region(fna, idx, c, s, e)
    return {c: clf.classify_class(loci, c, None, seqlens, fetch)["status"]
            for c in ALL_CELLS}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession")
    ap.add_argument("--chunk-mbp", type=int, default=80,
                    help="chunk size in Mbp — small enough to force several")
    ap.add_argument("--threads", type=int, default=2,
                    help="kept low: the sweep may be running alongside")
    args = ap.parse_args()

    manifest = {r["accession"]: r for r in read_manifest()}
    if args.accession:
        acc = args.accession
    else:                       # smallest genome already on disk
        have = [a for a in manifest if fna_path(a) is not None]
        if not have:
            raise SystemExit("no genome on disk — run the sweep first")
        acc = min(have, key=lambda a: int(manifest[a]["total_length_bp"]))
    fna = fna_path(acc)
    if fna is None:
        raise SystemExit(f"{acc} is not on disk")
    row = manifest[acc]
    genome_bp = int(row["total_length_bp"])
    print(f"chunked-path equivalence test on {acc} — {row['organism']}, "
          f"{genome_bp / 1e9:.2f} Gbp")

    baits = get_data_root() / "genome_sweep" / "baits.faa"
    if not baits.exists():
        seqs, _ = lib.load_baits()
        baits.parent.mkdir(parents=True, exist_ok=True)
        lib.write_fasta(seqs, baits)
    _, bait_meta = lib.load_baits()
    max_intron = lib.max_intron_for(genome_bp)

    work = get_data_root() / "genome_sweep" / "_chunk_test"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    idx = gio.read_fai(gio.build_fai(fna))
    seqlens = {k: v[0] for k, v in idx.items()}

    try:
        # --- whole-genome reference
        t0 = time.time()
        whole = work / "whole.gff"
        lib.run_miniprot(fna, baits, whole, threads=args.threads,
                         max_intron=max_intron)
        t_whole, rss_whole = time.time() - t0, peak_rss_gb()
        loci_a = lib.cluster_loci(lib.parse_miniprot_gff(whole, bait_meta))
        cells_a = cell_signature(loci_a, fna, idx, seqlens)
        print(f"  unchunked: {len(loci_a)} loci in {t_whole:.0f}s, "
              f"peak child RSS {rss_whole:.2f} GB "
              f"({rss_whole / (genome_bp / 1e9):.2f} GB per Gbp)")

        # --- chunked
        t0 = time.time()
        chunked = work / "chunked.gff"
        n_chunks = lib.run_miniprot_chunked(
            fna, baits, chunked, work / "chunks", threads=args.threads,
            max_intron=max_intron, chunk_bp=args.chunk_mbp * 1_000_000)
        t_chunk, rss_after = time.time() - t0, peak_rss_gb()
        loci_b = lib.cluster_loci(lib.parse_miniprot_gff(chunked, bait_meta))
        cells_b = cell_signature(loci_b, fna, idx, seqlens)
        plan = json.loads((work / "chunks" / "chunk_plan.json").read_text())
        print(f"  chunked  : {len(loci_b)} loci in {t_chunk:.0f}s over "
              f"{n_chunks} chunks; high-water RSS now {rss_after:.2f} GB "
              + ("(a new peak — the chunked path used more)"
                 if rss_after > rss_whole + 0.01
                 else "(unchanged, so each chunk stayed under the whole-genome peak)"))
        print(f"             chunk sizes (Gbp): "
              + ", ".join(f"{b / 1e9:.2f}" for b in plan["chunk_bp_actual"]))

        # --- the assertions
        failures, worst_delta = compare_loci(loci_key(loci_a), loci_key(loci_b))
        if len(loci_a) != len(loci_b):
            failures.insert(0, f"locus count differs: {len(loci_a)} unchunked "
                               f"vs {len(loci_b)} chunked")
        if cells_a != cells_b:
            failures.append(f"cell statuses differ:\n      unchunked {cells_a}"
                            f"\n      chunked   {cells_b}")
        # every chunk must hold whole contigs — the coordinate claim
        contigs_by_chunk = plan.get("chunk_contigs", [])
        if sum(contigs_by_chunk) != len(idx):
            failures.append(
                f"chunks hold {sum(contigs_by_chunk)} contigs but the genome "
                f"has {len(idx)} — a contig was split or dropped")
        # the mp_id collision guard: a multi-chunk GFF must reuse raw IDs
        raw_ids = [ln.split("ID=")[1].split(";")[0]
                   for ln in chunked.read_text().splitlines()
                   if "\tmRNA\t" in ln and "ID=" in ln]
        if n_chunks > 1 and len(raw_ids) == len(set(raw_ids)):
            print("  note: no duplicate mRNA IDs across chunks in this "
                  "genome, so the collision guard was not exercised here")
        elif n_chunks > 1:
            print(f"  collision guard exercised: {len(raw_ids)} mRNA rows, "
                  f"{len(set(raw_ids))} distinct raw IDs")

        per_gbp = rss_whole / (genome_bp / 1e9)
        proj = per_gbp * lib.CHUNK_BP / 1e9
        print(f"\n  worst boundary difference: {worst_delta} bp "
              f"(tolerance {BOUNDARY_TOLERANCE_BP})")
        print(f"  memory: {per_gbp:.2f} GB RSS per Gbp -> a "
              f"{lib.CHUNK_BP / 1e9:.1f} Gbp production chunk projects to "
              f"~{proj:.1f} GB")

        if failures:
            print("\nFAILED:")
            for f in failures:
                print(f"  - {f}")
            raise SystemExit(1)
        print(f"\nPASS: {len(loci_a)} loci and all {len(ALL_CELLS)} cell "
              "statuses identical between the two paths")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
