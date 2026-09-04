"""s5_run_sweep.py — S5 driver: per-genome miniprot sweep + tblastn rescue.

For each manifest genome: fetch (resumable, md5-verified), download the GFF3
annotation if the assembly has one, miniprot the committed bait panel
(`--gff --trans`), cluster alignments into loci, classify ITPR1/ITPR2/ITPR3
**and the RyR control**, and for cells with no locus run a tblastn rescue
whose hits are only counted OUTSIDE already-found loci — otherwise
cross-paralog similarity would make "absent" impossible to call.

The RyR cell is the internal positive control (session brief, D14). A genome
where it comes back `absent` or `no_locus` has an assembly or pipeline
problem, and the driver says so on the line rather than leaving it to be
noticed in the ledger.

Outputs per genome -> <data_root>/genome_sweep/<acc>/
  miniprot.gff          raw miniprot output (##PAF/##STA + GFF3)
  genes_slim.tsv        annotated-gene intervals (if the assembly is annotated)
  tblastn_<class>.tsv   rescue output (only when rescue ran)
  novel_models.faa      translations of unannotated / fragmentary loci
  summary.json          the per-genome ledger record
  .sweep.done / .sweep.failed   resume markers

Usage:
  python scripts/s5_run_sweep.py --only GCF_000002315.6
  python scripts/s5_run_sweep.py --limit 5 --threads 8 --delete-after
  python scripts/s5_run_sweep.py --retry-failed --redo
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

for extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(extra).is_dir() and extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root            # noqa: E402
from fetch_genomes import (read_manifest, fetch_one, purge_one,  # noqa: E402
                           fna_path)
import s5_genome_io as gio                                # noqa: E402
import s5_sweep_lib as lib                                # noqa: E402
import s5_classify as clf                                 # noqa: E402
import s5_rescue as rescue_lib                            # noqa: E402

SWEEP_VERSION = "s5.1"
ALL_CELLS = (*lib.CLASSES, lib.CONTROL_CLASS)

#: Statuses whose gene models are worth carrying into census v4.
NOVEL_STATUSES = {"found_unannotated", "found_no_annotation",
                  "fragment", "assembly_gap"}


def sweep_root() -> Path:
    p = get_data_root() / "genome_sweep"
    p.mkdir(parents=True, exist_ok=True)
    return p


def sweep_dir(acc: str) -> Path:
    p = sweep_root() / acc
    p.mkdir(parents=True, exist_ok=True)
    return p


def ensure_bait_files() -> tuple[Path, dict, dict]:
    """Write the shared bait FASTA + per-cell rescue FASTAs once."""
    root = sweep_root()
    seqs, meta = lib.load_baits()

    def write_if_changed(path: Path, content: dict[str, str]) -> None:
        tmp = path.with_suffix(path.suffix + ".tmp")
        lib.write_fasta(content, tmp)
        if path.exists() and path.read_bytes() == tmp.read_bytes():
            tmp.unlink()
            return
        tmp.replace(path)          # panel changed -> summaries go stale

    baits = root / "baits.faa"
    write_if_changed(baits, seqs)
    rescue = {}
    for cell, ids in rescue_lib.rescue_baits(meta).items():
        p = root / f"rescue_{cell}.faa"
        write_if_changed(p, {i: seqs[i] for i in ids if i in seqs})
        rescue[cell] = p
    return baits, rescue, meta


def filter_hsps_outside(hsps: list[dict], known: list[tuple], pad: int = 5000):
    out = []
    for h in hsps:
        lo, hi = sorted((h["sstart"], h["send"]))
        if any(c == h["contig"] and lo <= e + pad and hi >= s - pad
               for c, s, e in known):
            continue
        out.append(h)
    return out


def load_gene_index(row: dict, acc: str, gdir: Path, out: Path):
    """The assembly's annotated genes, downloaded once and cached slim."""
    if row.get("annotated") != "Y":
        return None
    gff_gz = gio.ensure_annotation(acc, gdir, quiet=True)
    if not gff_gz:
        return None
    slim = out / "genes_slim.tsv"
    if not slim.exists():
        genes = gio.slim_genes(gff_gz)
        with open(slim, "w") as fh:
            for g in genes:
                fh.write("\t".join(map(str, g)) + "\n")
    else:
        genes = [(f[0], int(f[1]), int(f[2]), f[3], f[4], f[5], f[6])
                 for f in (ln.split("\t") for ln in
                           slim.read_text().splitlines())]
    return gio.GeneIndex(genes)


def run_alignment(row: dict, fna: Path, baits: Path, out: Path,
                  threads: int) -> tuple[list, list, int, int, str]:
    """miniprot (chunked if the genome is too big to index) -> loci."""
    gff = out / "miniprot.gff"
    genome_bp = int(row.get("total_length_bp") or 0)
    max_intron = lib.max_intron_for(genome_bp)
    prev = {}
    if (out / "summary.json").exists():
        try:
            prev = json.loads((out / "summary.json").read_text())
        except json.JSONDecodeError:
            prev = {}
    baits_version = hashlib.md5(baits.read_bytes()).hexdigest()[:12]
    stale = (prev.get("max_intron") != max_intron
             or prev.get("baits_version") != baits_version)
    n_chunks = 1
    if not gff.exists() or gff.stat().st_size == 0 or stale:
        if genome_bp > lib.CHUNK_BP:          # too big for one miniprot index
            n_chunks = lib.run_miniprot_chunked(
                fna, baits, gff, out / "_chunks", threads, max_intron)
            shutil.rmtree(out / "_chunks", ignore_errors=True)
        else:
            lib.run_miniprot(fna, baits, gff, threads=threads,
                             max_intron=max_intron)
    else:
        n_chunks = prev.get("miniprot_chunks", 1)
    alns = lib.parse_miniprot_gff(gff, lib.load_baits()[1])
    return alns, lib.cluster_loci(alns), max_intron, n_chunks, baits_version


def run_rescue(cells: dict, fna: Path, out: Path, rescue: dict, baits: Path,
               bait_meta: dict, gene_index, known: list, fetch_fn,
               threads: int) -> None:
    need = [c for c in ALL_CELLS if cells[c]["status"] == "no_locus"]
    if not need:
        return
    db = out / "blastdb" / "genome"
    fresh = [c for c in need
             if not (out / f"tblastn_{c}.tsv").exists()
             or (out / f"tblastn_{c}.tsv").stat().st_size == 0]
    if fresh:                       # skip the db build when every tsv is cached
        rescue_lib.run_makeblastdb(fna, db)
    for c in need:
        tsv = out / f"tblastn_{c}.tsv"
        if c in fresh:
            rescue_lib.run_tblastn(rescue[c], db, tsv, threads=threads)
        hsps = filter_hsps_outside(rescue_lib.parse_tblastn(tsv), known)
        summ = rescue_lib.summarize_tblastn(hsps)
        regions = rescue_lib.cluster_hsps(hsps)[:rescue_lib.MAX_REGIONS]
        for r in regions:
            rescue_lib.attribute_region(fetch_fn, r, baits,
                                        out / "_regions", bait_meta)
        rescue_lib.annotate_regions(regions, gene_index)
        own, ambiguous, elsewhere = rescue_lib.split_by_margin(
            regions, c, bait_meta)
        cells[c]["rescue"] = summ or {"n_hsps": 0}
        cells[c]["rescue_regions"] = own
        cells[c]["ambiguous_regions"] = ambiguous
        cells[c]["attributed_elsewhere"] = elsewhere
        cells[c]["status"] = ("tblastn_trace" if own else
                              "tblastn_trace_ambiguous" if ambiguous
                              else "absent")
    shutil.rmtree(out / "_regions", ignore_errors=True)
    if fresh:
        shutil.rmtree(out / "blastdb", ignore_errors=True)   # ~1 GB/genome


def collect_novel(acc: str, alns: list, cells: dict, others: list) -> dict:
    """Gene models worth carrying into census v4.

    The rule is one sentence: keep every locus whose annotation does not
    already name it correctly. That is wider than "unannotated", and
    deliberately so — the pilot's first genome found the two teleost ITPR1
    3R duplicates sitting side by side, one annotated `itpr1b` and one left
    as an unnamed `LOC101074739`. The unnamed one is a real ITPR1 gene that
    no name-based census can reach, so it has to travel into census v4; a
    rule keying on `annot_gene` alone would have dropped it, because there
    *is* a gene model there. A locus named for the wrong paralog is kept for
    the same reason — that is the correction list S18 is after.
    """
    tr_by_id = {a.mp_id: a.translation for a in alns if a.translation}
    novel: dict[str, str] = {}

    def add(d: dict, label: str) -> None:
        tr = tr_by_id.get(d["mp_id"], "")
        if not tr:
            return
        novel[f"{acc}|{label}|{d['contig']}:{d['start']}-{d['end']}"
              f"{d['strand']}|cov{d['coverage']:.2f}|id{d['identity']:.2f}"
              f"|len{len(tr)}"] = tr

    for c in ALL_CELLS:
        cell = cells[c]
        for d in cell.get("loci", []):
            named_right = bool(d.get("annot_gene")) and bool(
                d.get("annot_paralog_matches"))
            if cell["status"] in NOVEL_STATUSES or not named_right:
                add(d, c)
    for d in others:
        add(d, d["clade"])
    return novel


def process_genome(row: dict, baits: Path, rescue: dict, bait_meta: dict,
                   threads: int, no_rescue: bool) -> dict:
    acc = row["accession"]
    out = sweep_dir(acc)
    t0 = time.time()
    timings: dict[str, float] = {}

    gdir = fetch_one(acc, quiet=True)
    fna = fna_path(acc)
    if fna is None:
        raise RuntimeError("no .fna after fetch")
    timings["fetch_s"] = round(time.time() - t0, 1)

    t = time.time()
    gene_index = load_gene_index(row, acc, gdir, out)
    timings["annotation_s"] = round(time.time() - t, 1)

    t = time.time()
    idx = gio.read_fai(gio.build_fai(fna))
    seqlens = {k: v[0] for k, v in idx.items()}
    timings["fai_s"] = round(time.time() - t, 1)

    t = time.time()
    alns, loci, max_intron, n_chunks, baits_version = run_alignment(
        row, fna, baits, out, threads)
    timings["miniprot_s"] = round(time.time() - t, 1)

    def fetch_fn(c, s, e):
        return gio.fetch_region(fna, idx, c, s, e)

    cells = {c: clf.classify_class(loci, c, gene_index, seqlens, fetch_fn)
             for c in ALL_CELLS}
    others = clf.other_loci(loci, gene_index, seqlens, fetch_fn)
    known = [(L.contig, L.start, L.end) for L in loci]

    if not no_rescue:
        t = time.time()
        run_rescue(cells, fna, out, rescue, baits, bait_meta, gene_index,
                   known, fetch_fn, threads)
        timings["rescue_s"] = round(time.time() - t, 1)

    novel = collect_novel(acc, alns, cells, others)
    if novel:
        lib.write_fasta(novel, out / "novel_models.faa")

    meta = gio.read_assembly_meta(gdir)
    control = cells[lib.CONTROL_CLASS]
    summary = {
        "accession": acc, "organism": row["organism"],
        "vclass": row.get("vclass", ""), "vorder": row.get("vorder", ""),
        "reasons": row.get("reasons", ""), "annotated": row.get("annotated"),
        "assembly_level": meta.get("assembly_level") or row.get("level"),
        "has_annotation_index": gene_index is not None,
        "sweep_version": SWEEP_VERSION,
        "max_intron": max_intron,
        "max_intron_capped": lib.intron_rule_capped(
            int(row.get("total_length_bp") or 0)),
        "miniprot_chunks": n_chunks,
        "baits_version": baits_version, "n_baits": len(bait_meta),
        "genome_bp": int(row.get("total_length_bp") or 0),
        "contig_n50": int(row.get("contig_n50") or 0),
        "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_alignments": len(alns), "n_loci": len(loci),
        "cells": cells, "other_loci": others,
        "control_ok": control["status"] not in ("absent", "no_locus"),
        "control_status": control["status"],
        "control_n_loci": control["n_loci"],
        "novel_model_count": len(novel),
        "timings": timings,
        "total_s": round(time.time() - t0, 1),
    }
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    return summary


def one_line(summary: dict) -> str:
    bits = []
    for c in lib.CLASSES:
        cell = summary["cells"][c]
        extra = f"x{cell['n_loci']}" if cell["n_loci"] > 1 else ""
        bits.append(f"{c}:{cell['status']}{extra}")
    ctrl = summary["cells"][lib.CONTROL_CLASS]
    bits.append(f"[control RYR:{ctrl['status']}x{ctrl['n_loci']}]"
                + ("" if summary["control_ok"] else "  ** CONTROL FAILED **"))
    return "  ".join(bits) + (f"  (+{len(summary['other_loci'])} other)"
                              if summary["other_loci"] else "")


def write_live(pending: list, done: int, current: str) -> None:
    """Dashboard live panel."""
    path = PROJECT_ROOT / "results" / "session_live.json"
    try:
        path.write_text(json.dumps({
            "task": "S5", "workers": 1,
            "steps": [{"label": f"{r['accession']} {r['organism']}",
                       "done": i < done} for i, r in enumerate(pending)],
            "current": current}, indent=1))
    except OSError:
        pass


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="comma-separated accessions")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-fasta-gb", type=float, default=12.0,
                    help="skip genomes above this FASTA estimate "
                         "(the giants get their own chunked run)")
    ap.add_argument("--smallest-first", action="store_true",
                    help="pure size order (default: margin species first)")
    ap.add_argument("--no-rescue", action="store_true")
    ap.add_argument("--redo", action="store_true", help="rerun even if .sweep.done")
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--delete-after", action="store_true",
                    help="purge the genome dir after a successful sweep")
    ap.add_argument("--no-prefetch", action="store_true")
    ap.add_argument("--shard", help="i/N — every Nth pending genome, so "
                                    "parallel workers stay disjoint")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = read_manifest()
    if args.only:
        keep = {a.strip() for a in args.only.split(",")}
        rows = [r for r in rows if r["accession"] in keep]
    else:
        rows = [r for r in rows if float(r["fasta_gb"]) <= args.max_fasta_gb]

    def is_margin(r):
        return not r["reasons"].startswith("order_rep")
    if args.smallest_first:
        rows.sort(key=lambda r: float(r["fasta_gb"]))
    else:
        rows.sort(key=lambda r: (0 if is_margin(r) else 1, float(r["fasta_gb"])))

    pending = []
    for r in rows:
        d = sweep_root() / r["accession"]
        if (d / ".sweep.done").exists() and not args.redo:
            continue
        if (d / ".sweep.failed").exists() and not (args.retry_failed or args.redo):
            continue
        pending.append(r)
    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        pending = [r for k, r in enumerate(pending) if k % n == i]
        print(f"shard {i}/{n}: ", end="")
    if args.limit:
        pending = pending[:args.limit]
    print(f"{len(pending)} genomes to sweep "
          f"({sum(float(r['fasta_gb']) for r in pending):.0f} GB FASTA est.)")
    if args.dry_run:
        for r in pending:
            print(f"  {r['accession']}  {r['organism']}  ~{r['fasta_gb']} GB  "
                  f"-G {lib.max_intron_for(int(r['total_length_bp'] or 0)):,}  "
                  f"[{r['reasons']}]")
        return

    baits, rescue, bait_meta = ensure_bait_files()
    prefetch_thread = None

    def prefetch(acc: str):
        try:
            fetch_one(acc, quiet=True)
        except Exception as e:                            # noqa: BLE001
            print(f"  (prefetch {acc} failed: {e})")

    n_ok = n_fail = n_control_fail = 0
    for i, row in enumerate(pending):
        acc = row["accession"]
        if prefetch_thread:
            prefetch_thread.join()
        if not args.no_prefetch and i + 1 < len(pending):
            prefetch_thread = threading.Thread(
                target=prefetch, args=(pending[i + 1]["accession"],), daemon=True)
            prefetch_thread.start()
        else:
            prefetch_thread = None
        print(f"[{i + 1}/{len(pending)}] {acc} {row['organism']} "
              f"(~{row['fasta_gb']} GB) ...", flush=True)
        write_live(pending, i, f"{acc} {row['organism']}")
        out = sweep_dir(acc)
        try:
            summary = process_genome(row, baits, rescue, bait_meta,
                                     args.threads, args.no_rescue)
            (out / ".sweep.done").write_text(
                f"{SWEEP_VERSION} {datetime.now(timezone.utc).isoformat()}\n")
            (out / ".sweep.failed").unlink(missing_ok=True)
            n_ok += 1
            n_control_fail += 0 if summary["control_ok"] else 1
            print(f"  {one_line(summary)}  [{summary['total_s']:.0f}s]",
                  flush=True)
            if args.delete_after:
                purge_one(acc)
        except Exception as e:                            # noqa: BLE001
            n_fail += 1
            (out / ".sweep.failed").write_text(
                f"{datetime.now(timezone.utc).isoformat()}\n{e}\n"
                + traceback.format_exc())
            print(f"  ! FAILED: {e}", flush=True)
    write_live(pending, len(pending), "")
    print(f"done: {n_ok} swept, {n_fail} failed, "
          f"{len(pending) - n_ok - n_fail} skipped")
    if n_control_fail:
        print(f"  ** {n_control_fail} genome(s) failed the RyR positive "
              "control — assembly or pipeline, not biology; investigate "
              "before their ITPR cells are read")


if __name__ == "__main__":
    main()
