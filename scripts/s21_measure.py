"""The measurement pass: every locus's blocks, and the genome under them.

One pass per genome. For each locus the sweep placed, this reads the model's
CDS blocks out of the retained `miniprot.gff` and the genome under the model's
own span out of the retained FASTA, and records **every consecutive block
pair** with the gap between them in gene order, the splice dinucleotides at
that gap's two edges and whether the reading frame carries over.

The pair table is where the two things S21 cannot assume come from.

* **What an intron is.** A gap below the real minimum intron length is an indel
  the aligner walked around, not an intron, and merging those is the difference
  between a 58-exon gene and a 70-exon one. The threshold is calibrated in
  `s21_calibrate_intron` against the splice dinucleotides, which are evidence
  the alignment score did not produce.
* **Whether a boundary is trustworthy.** A junction whose flanking
  dinucleotides are `GT..AG` is a splice site the genome agrees with; one that
  is not may still be an intron the aligner placed a few bases off, and a
  boundary claim resting on it is weaker. Both are recorded per junction rather
  than aggregated away.

The pair table is bulk (about 150,000 rows over the sweep) and is cached under
the data root; every committed table is derived from it, so the whole task
re-runs offline and deterministically once this pass has been made.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_blocks as B                                            # noqa: E402
import s21_lib as L                                               # noqa: E402

PAIR_COLS = ("accession", "cell", "locus_idx", "mp_id", "bait", "strand",
             "pair_index", "prev_end_q", "next_start_q", "gap_bp",
             "intron_start", "intron_end", "donor", "acceptor",
             "splice_class", "frame_ok", "next_phase_gff", "prev_len_bp",
             "next_len_bp")

MODEL_COLS = ("accession", "organism", "vclass", "vorder", "cell", "locus_idx",
              "mp_id", "bait", "bait_paralog", "bait_len", "strand", "contig",
              "start", "end", "identity", "coverage", "aligned_aa",
              "frameshifts", "stop_codons", "n_blocks", "frame_steps",
              "span_bp", "cds_bp", "q_start", "q_end", "residues",
              "frame_balance_bp", "stop_codon_bp", "contig_spans_gene",
              "is_copy", "integrity", "lesion_density", "cell_status")


def cache_dir() -> Path:
    d = L.data_root() / "s21"
    d.mkdir(parents=True, exist_ok=True)
    return d


def measure_genome(acc: str, summary: dict, loci: list[dict],
                   labels: dict) -> tuple[list[dict], list[dict]]:
    """Model rows and block-pair rows for one genome.

    Every locus the sweep recorded is measured, whatever its status: a
    fragment's exon count is not architecture, but the *reason* it is not is a
    number S21 has to be able to print, and excluding it here would leave the
    scope table with no denominator.
    """
    gff = L.sweep_dir(acc) / "miniprot.gff"
    fna = L.genome_fna(acc)
    if not gff.exists():
        return [], []
    by_mp = {r["mp_id"]: r for r in loci}
    models = B.read_models(gff, set(by_mp))
    mrows: list[dict] = []
    prows: list[dict] = []
    for mp_id, rec in by_mp.items():
        m = models.get(mp_id)
        if m is None or not m["blocks"]:
            continue
        arch = B.architecture(m, 1)          # no merging: raw blocks
        bait_acc = rec["bait"].split("|")[0]
        mrows.append({
            "accession": acc, "organism": summary.get("organism", ""),
            "vclass": summary.get("vclass", ""),
            "vorder": summary.get("vorder", ""),
            "cell": rec["cell"], "locus_idx": rec["locus_idx"], "mp_id": mp_id,
            "bait": bait_acc,
            "bait_paralog": (labels.get(bait_acc) or {}).get("paralog", ""),
            "bait_len": rec["bait_len"], "strand": m["strand"],
            "contig": m["contig"], "start": m["start"], "end": m["end"],
            "identity": rec["identity"], "coverage": rec["coverage"],
            "aligned_aa": rec["aligned_aa"], "frameshifts": m["frameshifts"],
            "stop_codons": m["stop_codons"], "n_blocks": arch["n_blocks"],
            "frame_steps": arch["frame_steps"], "span_bp": arch["span_bp"],
            "cds_bp": arch["cds_bp"], "q_start": arch["q_start"],
            "q_end": arch["q_end"], "residues": arch["residues"],
            "frame_balance_bp": arch["frame_balance_bp"],
            "stop_codon_bp": arch["stop_codon_bp"],
            "contig_spans_gene": int(rec["contig_spans_gene"]),
            "is_copy": int(rec["is_copy"]), "integrity": rec.get("integrity", ""),
            "lesion_density": rec.get("lesion_density", 0.0),
            "cell_status": rec.get("cell_status", "")})
        region = None
        if fna is not None:
            region = L.LocusRegion(acc, fna, m["contig"], m["start"], m["end"])
        bl = m["blocks"]
        for i, (p, n) in enumerate(zip(bl, bl[1:])):
            g = B.gap_bp(p, n, m["strand"])
            lo, hi = B.intron_interval(p, n, m["strand"])
            donor = acceptor = ""
            if region is not None and hi >= lo:
                donor, acceptor = L.splice_pair(region, (lo, hi), m["strand"])
            prows.append({
                "accession": acc, "cell": rec["cell"],
                "locus_idx": rec["locus_idx"], "mp_id": mp_id,
                "bait": bait_acc, "strand": m["strand"], "pair_index": i,
                "prev_end_q": p["q_end"], "next_start_q": n["q_start"],
                "gap_bp": g, "intron_start": lo, "intron_end": hi,
                "donor": donor, "acceptor": acceptor,
                "splice_class": L.splice_class(donor, acceptor),
                "frame_ok": int(B.next_phase(p) == n["phase"]),
                "next_phase_gff": n["phase"],
                "prev_len_bp": p["end"] - p["start"] + 1,
                "next_len_bp": n["end"] - n["start"] + 1})
    return mrows, prows


def run(on_progress=None) -> tuple[Path, Path]:
    """Measure every genome, writing the two cache tables under the data root."""
    labels = L.bait_labels()
    summaries = L.load_summaries()
    loci = L.loci_with_mp(summaries)
    by_acc: dict[str, list[dict]] = {}
    for r in loci:
        by_acc.setdefault(r["accession"], []).append(r)
    accs = sorted(by_acc)
    mrows: list[dict] = []
    prows: list[dict] = []
    t0 = time.time()
    for i, acc in enumerate(accs, 1):
        summary = summaries.get(acc) or {}
        m, p = measure_genome(acc, summary, by_acc[acc], labels)
        mrows += m
        prows += p
        if on_progress and (i % 20 == 0 or i == len(accs)):
            on_progress(i, len(accs), time.time() - t0)
    mp = cache_dir() / "model_rows.tsv"
    pp = cache_dir() / "block_pairs.tsv"
    L.write_tsv(mp, mrows, list(MODEL_COLS))
    L.write_tsv(pp, prows, list(PAIR_COLS))
    (cache_dir() / "measure_stats.json").write_text(json.dumps({
        "n_genomes": len(accs), "n_models": len(mrows), "n_pairs": len(prows),
        "seconds": round(time.time() - t0, 1)}, indent=1))
    return mp, pp


def load_models() -> list[dict]:
    rows = L.read_tsv(cache_dir() / "model_rows.tsv")
    for r in rows:
        for k in ("locus_idx", "bait_len", "start", "end", "aligned_aa",
                  "frameshifts", "stop_codons", "n_blocks", "frame_steps",
                  "span_bp", "cds_bp", "q_start", "q_end", "residues",
                  "frame_balance_bp", "stop_codon_bp", "contig_spans_gene",
                  "is_copy"):
            r[k] = int(float(r.get(k) or 0))
        for k in ("identity", "coverage", "lesion_density"):
            r[k] = float(r.get(k) or 0.0)
    return rows


def load_pairs() -> list[dict]:
    rows = L.read_tsv(cache_dir() / "block_pairs.tsv")
    for r in rows:
        for k in ("locus_idx", "pair_index", "prev_end_q", "next_start_q",
                  "gap_bp", "intron_start", "intron_end", "frame_ok",
                  "next_phase_gff", "prev_len_bp", "next_len_bp"):
            r[k] = int(float(r.get(k) or 0))
    return rows


if __name__ == "__main__":
    def prog(i, n, s):
        print(f"  [measure] {i}/{n} genomes  {s:.0f}s", flush=True)
    m, p = run(prog)
    print(m, p)
