"""s23_run_sweep.py — S23 driver: per-genome copy-number sweep + rescue.

For each manifest genome: fetch (resumable, md5-verified), download the GFF3 if
the assembly has one, miniprot the committed S23 panel, cluster alignments into
loci, and classify **copy number** rather than paralog cells
(`s23_classify.py`). Where the ITPR baits win nothing, a tblastn rescue asks
whether any trace is there at all.

**The control logic is where this differs from S5 and it is not a detail.**
S5 can say "a genome with no RyR locus is broken" because every vertebrate has
three RyRs. Outside Metazoa that is false — S20 found architecture-level RyR in
2 of 6,928 non-metazoan proteomes — so this sweep carries the **MIR-domain
sharer** as its universal positive control (bait rule B2) and RyR only where a
RyR is expected. A genome that fails its control is not a negative result; it
is excluded from every absence claim, and `s23_classify.control_verdict` is
what decides which.

Outputs per genome -> <data_root>/s23_sweep/<acc>/
  miniprot.gff          raw miniprot output
  genes_slim.tsv        annotated-gene intervals (if the assembly is annotated)
  tblastn_itpr.tsv      rescue output (only when rescue ran)
  novel_models.faa      translations of unannotated / fragmentary ITPR loci
  summary.json          the per-genome ledger record
  .sweep.done / .sweep.failed   resume markers

Usage:
  python3 scripts/s23_run_sweep.py --only GCF_000001215.4
  python3 scripts/s23_run_sweep.py --reasons anchor --threads 8
  python3 scripts/s23_run_sweep.py --limit 10 --smallest-first --delete-after
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

for _extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(_extra).is_dir() and _extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = _extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root                  # noqa: E402
from fetch_genomes import fetch_one, fna_path, purge_one       # noqa: E402
import s23_bait_spec as spec                                   # noqa: E402
import s23_calibration as cal                                  # noqa: E402
import s23_classify as clf                                     # noqa: E402
import s5_genome_io as gio                                     # noqa: E402
import s5_rescue as rescue_lib                                 # noqa: E402
import s5_sweep_lib as lib                                     # noqa: E402

SWEEP_VERSION = "s23.1"
MANIFEST = PROJECT_ROOT / "results" / "s23_scope" / "genome_manifest_s23.tsv"
BAITS_DIR = PROJECT_ROOT / "results" / "s23_baits"
#: How many ITPR baits a rescue query carries, spread across clade bands.
RESCUE_BAITS = 6


def sweep_root() -> Path:
    d = get_data_root() / "s23_sweep"
    d.mkdir(parents=True, exist_ok=True)
    return d


def sweep_dir(acc: str) -> Path:
    d = sweep_root() / acc
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_manifest() -> list[dict]:
    import csv
    with open(MANIFEST) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ensure_bait_files() -> tuple[Path, Path, dict, dict]:
    """(panel faa, rescue faa, id->seq, id->meta) written under the data root.

    The rescue query is a *spread* of ITPR baits across clade bands, derived
    from the committed manifest rather than listed, so it cannot drift from
    the panel it is supposed to represent.
    """
    seqs, meta = lib.load_baits(BAITS_DIR)
    work = sweep_root() / "_baits"
    work.mkdir(parents=True, exist_ok=True)
    panel = work / "baits.faa"
    lib.write_fasta(seqs, panel)

    picked: list[str] = []
    for band in spec.BAND_ORDER:
        hit = sorted(bid for bid, m in meta.items()
                     if m["family"] == spec.FAMILY_ITPR and m["band"] == band)
        if hit:
            picked.append(hit[0])
        if len(picked) >= RESCUE_BAITS:
            break
    rescue = work / "rescue_itpr.faa"
    lib.write_fasta({b: seqs[b] for b in picked}, rescue)
    return panel, rescue, seqs, meta


def load_gene_index(row: dict, acc: str, gdir: Path, out: Path):
    """The assembly's annotated genes, or None if it carries no gene set."""
    if row.get("annotated") != "Y":
        return None
    slim = out / "genes_slim.tsv"
    if not slim.exists() or slim.stat().st_size == 0:
        gff = gio.ensure_annotation(acc, gdir, quiet=True)
        if gff is None:
            return None
        genes = gio.slim_genes(gff)
        with open(slim, "w") as fh:
            for g in genes:
                fh.write("\t".join(str(x) for x in g) + "\n")
    genes = []
    with open(slim) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7:
                genes.append((f[0], int(f[1]), int(f[2]), f[3], f[4], f[5],
                              f[6]))
    return gio.GeneIndex(genes) if genes else None


def run_alignment(row: dict, fna: Path, panel: Path, out: Path, meta: dict,
                  threads: int) -> tuple[list, list, int, int, str]:
    gff = out / "miniprot.gff"
    group = row.get("group") or "other"
    max_intron = cal.max_intron_for(group)
    genome_bp = int(row.get("total_length_bp") or 0)
    prev = {}
    if (out / "summary.json").exists():
        try:
            prev = json.loads((out / "summary.json").read_text())
        except json.JSONDecodeError:
            prev = {}
    baits_version = hashlib.md5(panel.read_bytes()).hexdigest()[:12]
    stale = (prev.get("max_intron") != max_intron
             or prev.get("baits_version") != baits_version)
    n_chunks = 1
    if not gff.exists() or gff.stat().st_size == 0 or stale:
        if genome_bp > lib.CHUNK_BP:
            n_chunks = lib.run_miniprot_chunked(
                fna, panel, gff, out / "_chunks", threads, max_intron)
            shutil.rmtree(out / "_chunks", ignore_errors=True)
        else:
            lib.run_miniprot(fna, panel, gff, threads=threads,
                             max_intron=max_intron)
    else:
        n_chunks = prev.get("miniprot_chunks", 1)
    alns = lib.parse_miniprot_gff(gff, meta)
    loci = lib.filter_loci(lib.cluster_loci(alns))
    return alns, loci, max_intron, n_chunks, baits_version


def run_rescue(fna: Path, out: Path, rescue_faa: Path, panel: Path,
               meta: dict, gene_index, known: list, fetch_fn,
               threads: int) -> dict:
    """tblastn for a genome the ITPR baits won nothing in."""
    db = out / "blastdb" / "genome"
    tsv = out / "tblastn_itpr.tsv"
    fresh = not tsv.exists() or tsv.stat().st_size == 0
    if fresh:
        rescue_lib.run_makeblastdb(fna, db)
        rescue_lib.run_tblastn(rescue_faa, db, tsv, threads=threads)
    hsps = [h for h in rescue_lib.parse_tblastn(tsv)
            if not any(h["contig"] == c and h["start"] <= e + 5000
                       and h["end"] >= s - 5000 for c, s, e in known)]
    summ = rescue_lib.summarize_tblastn(hsps) or {"n_hsps": 0}
    regions = rescue_lib.cluster_hsps(hsps)[:rescue_lib.MAX_REGIONS]
    for r in regions:
        rescue_lib.attribute_region(fetch_fn, r, panel, out / "_regions", meta)
    rescue_lib.annotate_regions(regions, gene_index)
    # Family-level only. S5 had to attribute a rescue region to one of three
    # paralogs; here the question is whether anything of the family is there,
    # so the paralog margin (`ATTRIBUTION_REL_MARGIN`) never enters.
    itpr_regions = [r for r in regions
                    if rescue_lib.region_family(r, meta)[0] == spec.FAMILY_ITPR]
    shutil.rmtree(out / "_regions", ignore_errors=True)
    if fresh:
        shutil.rmtree(out / "blastdb", ignore_errors=True)
    return {"summary": summ, "n_regions": len(regions),
            "itpr_regions": itpr_regions[:6],
            "n_itpr_regions": len(itpr_regions)}


def collect_novel(acc: str, row: dict) -> dict:
    """Translations of ITPR loci that no annotated family gene covers."""
    out = {}
    for i, d in enumerate(row["loci"], 1):
        if d["grade"] == "scrap" or d.get("annot_names_family"):
            continue
        seq = d.pop("_translation", "")
        if seq:
            out[f"{acc}|locus{i}|{d['contig']}:{d['start']}-{d['end']}|"
                f"{d['grade']}"] = seq
    return out


def process_genome(row: dict, panel: Path, rescue_faa: Path, meta: dict,
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

    idx = gio.read_fai(gio.build_fai(fna))
    seqlens = {k: v[0] for k, v in idx.items()}

    t = time.time()
    alns, loci, max_intron, n_chunks, baits_version = run_alignment(
        row, fna, panel, out, meta, threads)
    timings["miniprot_s"] = round(time.time() - t, 1)

    def fetch_fn(c, s, e):
        return gio.fetch_region(fna, idx, c, s, e)

    # keep each locus's translation for the novel-model FASTA before the
    # classifier's rows lose the Aln objects
    by_pos = {(a.contig, a.start): a.translation for a in alns}
    cn = clf.classify_genome(loci, gene_index, seqlens, fetch_fn)
    for d in cn["loci"]:
        d["_translation"] = by_pos.get((d["contig"], d["start"]), "")

    group = row.get("group") or "other"
    expects_ryr = group == "metazoa"
    # Does the panel carry a control bait this genome's clade could match?
    # Controls are built per *absence* clade (bait rule B2), so most genomes
    # have none, and "no control fired" there means nothing.
    clades = {row.get("class", ""), row.get("phylum", "")}
    has_ctl = any(m.get("family") == spec.CONTROL_MIR
                  and m.get("band") in clades for m in meta.values())
    verdict, why = clf.control_verdict(cn, expects_ryr, has_ctl)

    rescue = None
    if cn["status"] == "no_locus" and not no_rescue:
        t = time.time()
        known = [(L.contig, L.start, L.end) for L in loci]
        rescue = run_rescue(fna, out, rescue_faa, panel, meta, gene_index,
                            known, fetch_fn, threads)
        timings["rescue_s"] = round(time.time() - t, 1)
        if rescue["n_itpr_regions"]:
            cn["status"] = "tblastn_trace"

    novel = collect_novel(acc, cn)
    if novel:
        lib.write_fasta(novel, out / "novel_models.faa")
    for d in cn["loci"]:
        d.pop("_translation", None)

    holds, bar, bar_why = cal.spans_a_gene(row.get("contig_n50", 0), group)
    meta_asm = gio.read_assembly_meta(gdir)
    summary = {
        "accession": acc, "organism": row["organism"], "group": group,
        "kingdom": row.get("kingdom", ""), "phylum": row.get("phylum", ""),
        "class": row.get("class", ""), "order": row.get("order", ""),
        "reasons": row.get("reasons", ""), "annotated": row.get("annotated"),
        "assembly_level": meta_asm.get("assembly_level") or row.get("level"),
        "has_annotation_index": gene_index is not None,
        "sweep_version": SWEEP_VERSION, "max_intron": max_intron,
        "miniprot_chunks": n_chunks, "baits_version": baits_version,
        "n_baits": len(meta),
        "genome_bp": int(row.get("total_length_bp") or 0),
        "contig_n50": int(row.get("contig_n50") or 0),
        "spans_gene": holds, "contiguity_bar_bp": bar, "bar_source": bar_why,
        "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_alignments": len(alns), "n_loci_all_families": len(loci),
        "min_locus_identity": lib.MIN_LOCUS_IDENTITY,
        "expects_ryr": expects_ryr, "has_control_bait": has_ctl,
        "control_verdict": verdict, "control_why": why,
        "rescue": rescue,
        "novel_model_count": len(novel),
        "timings": timings, "total_s": round(time.time() - t0, 1),
        **{k: v for k, v in cn.items() if k != "loci"},
        "loci": cn["loci"],
    }
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    return summary


def one_line(s: dict) -> str:
    bits = [f"{s['status']}", f"n_full={s['n_full']}"]
    if s["n_fragment"] or s["n_scrap"]:
        bits.append(f"(+{s['n_fragment']}frag/{s['n_scrap']}scrap)")
    if s["n_merges"]:
        bits.append(f"[{s['n_merges']} split-merge]")
    bits.append(f"MIR:{s['mir_loci']}")
    if s["expects_ryr"]:
        bits.append(f"RYR:{s['ryr_loci']}")
    if s["control_verdict"] == "uncontrolled":
        bits.append("** UNCONTROLLED **")
    if not s["spans_gene"]:
        bits.append("[below contiguity bar]")
    return "  ".join(bits)


def write_live(done: int, total: int, current: str) -> None:
    path = PROJECT_ROOT / "results" / "session_live.json"
    try:
        path.write_text(json.dumps({
            "task": "S23a", "workers": 1,
            "steps": [{"label": f"genomes swept {done}/{total}",
                       "done": done >= total},
                      {"label": f"current: {current}", "done": False}],
        }, indent=1))
    except OSError:
        pass


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--reasons", action="append", default=[],
                    help="only genomes carrying this scope rule (G1-G5)")
    ap.add_argument("--groups", action="append", default=[])
    ap.add_argument("--limit", type=int)
    ap.add_argument("--smallest-first", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--no-rescue", action="store_true")
    ap.add_argument("--delete-after", action="store_true")
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()

    rows = read_manifest()
    if args.only:
        rows = [r for r in rows if r["accession"] in set(args.only)]
    if args.reasons:
        want = set(args.reasons)
        rows = [r for r in rows
                if want & set(r["reasons"].split(";"))]
    if args.groups:
        rows = [r for r in rows if r.get("group") in set(args.groups)]
    if args.smallest_first:
        rows.sort(key=lambda r: int(r.get("total_length_bp") or 0))
    if args.limit:
        rows = rows[:args.limit]

    panel, rescue_faa, _seqs, meta = ensure_bait_files()
    print(f"{len(rows)} genome(s); panel {len(meta)} baits\n{cal.summary()}\n")

    done = failed = 0
    for i, row in enumerate(rows, 1):
        acc = row["accession"]
        out = sweep_dir(acc)
        write_live(i - 1, len(rows), row["organism"])
        if (out / ".sweep.done").exists() and not args.redo:
            print(f"[{i}/{len(rows)}] {acc} {row['organism'][:34]:34s} cached")
            done += 1
            continue
        print(f"[{i}/{len(rows)}] {acc} {row['organism'][:34]:34s} "
              f"{int(row.get('total_length_bp') or 0)/1e6:.0f} Mbp … ",
              end="", flush=True)
        try:
            s = process_genome(row, panel, rescue_faa, meta, args.threads,
                               args.no_rescue)
            (out / ".sweep.done").write_text(
                datetime.now(timezone.utc).isoformat())
            (out / ".sweep.failed").unlink(missing_ok=True)
            print(f"{s['total_s']:.0f}s  {one_line(s)}")
            done += 1
        except Exception as exc:                       # noqa: BLE001
            (out / ".sweep.failed").write_text(traceback.format_exc())
            print(f"FAILED: {exc}")
            failed += 1
        finally:
            if args.delete_after:
                purge_one(acc)
    write_live(done, len(rows), "complete")
    print(f"\n{done} done, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
