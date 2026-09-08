"""S12 step 5 — the independent cross-check, and the species that can answer it.

The read evidence and the deposit evidence are different instruments.  S10
ran the deposit test on its two cases and could not answer it: *Nibea
albiflora* has 43 mRNA records and *Dissostichus eleginoides* 10, so a
search returning nothing showed the species has almost no transcript
deposits rather than that the gene is untranscribed.

S12's panel adds a third species, and it changes the picture: **NCBI holds
37,166 mRNA records for *Dissostichus mawsoni***.  There the same test has a
denominator that can answer, and it asks something the reads cannot — does
an *independently submitted, independently assembled* transcript cross the
junctions the annotation does not model?

Two controls, both necessary:

  * **the annotated junctions of the same genes**, probed identically.  If
    the deposits recover junctions the annotation models but none it does
    not, that is a fact about the deposit's own gene-model bias; if they
    recover both, the unannotated junctions are ordinary transcript.
  * **the locus's own genomic sequence**, S10's control unchanged.  A probe
    is 180 nt of contiguous *spliced* sequence, so nothing genomic can span
    its junction — many hits, zero spans, or the spanning rule is broken.

Outputs (committed):
    results/expression/atlas_resources.tsv
    results/expression/atlas_probes.tsv
    results/expression/atlas_summary.tsv

Run:  python scripts/s12_atlas.py [--species ...] [--retmax 40000]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s12_panel  # noqa: E402
from s10_evidence import Region, classify_introns  # noqa: E402
from s10_gff import read_annotation_window, read_miniprot_model  # noqa: E402
from s10_probes import (  # noqa: E402
    build_probes, fetch_species_nucleotides, genomic_control, resource_counts,
    search_probes,
)
from s12_lib import (  # noqa: E402
    CACHE, GENOMES, OUT_DIR, SWEEP, live, write_tsv,
)
from s12_refs import genome_files, locus_models  # noqa: E402

MANIFEST = Path(__file__).resolve().parents[1] / "results" / "genome_manifest.tsv"
ANNOTATED_CLASS = ("within_one_model",)
INFORMATIVE = ("between_models", "model_to_gap", "unannotated")


def taxids() -> dict[str, int]:
    out: dict[str, int] = {}
    if not MANIFEST.exists():
        return out
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out[r["accession"]] = int(r["taxid"])
            except (KeyError, ValueError):
                continue
    return out


def probe_locus(sp, locus, rec, log=print) -> tuple[list[dict], dict]:
    """Probe one locus's junctions against the species deposits + control."""
    fna, gff = genome_files(sp.accession)
    sweep_gff = SWEEP / sp.accession / "miniprot.gff"
    model = read_miniprot_model(sweep_gff, rec["mp_id"])
    if not model.get("cds"):
        return [], {}
    lo = min(b["start"] for b in model["cds"]) - 200
    hi = max(b["end"] for b in model["cds"]) + 200
    region = Region(fna, model["contig"], lo, hi)
    ann = read_annotation_window(gff, model["contig"], lo, hi) if gff else {}
    introns = classify_introns(model, ann.get("genes", []), region)

    cache = CACHE / "atlas" / sp.accession / locus.cell
    archive = CACHE / "atlas" / sp.accession
    probes = (build_probes(region, model, introns, classes=INFORMATIVE)
              + build_probes(region, model, introns, classes=ANNOTATED_CLASS))
    if not probes:
        return [], {}

    fetched = fetch_species_nucleotides(TAXID[sp.accession], archive,
                                        retmax=RETMAX)
    mrna = fetched.get("mrna", {})
    subject = mrna.get("path")
    deposits = search_probes(probes, subject, cache, "mrna") if subject \
        else {"status": "unavailable", "by_probe": {}}
    control = search_probes(probes, genomic_control(region, cache), cache,
                            "genomic")

    rows = []
    for p in probes:
        d = deposits["by_probe"].get(p["junction_index"], {})
        c = control["by_probe"].get(p["junction_index"], {})
        rows.append({
            "species": sp.short, "organism": sp.name,
            "accession": sp.accession, "cell": locus.cell, "role": locus.role,
            "junction_index": p["junction_index"],
            "intron_class": p["intron_class"],
            "annotated": int(p["intron_class"] in ANNOTATED_CLASS),
            "splice_class": p["splice_class"], "probe_len": p["probe_len"],
            "deposit_status": deposits["status"],
            "deposit_hits": d.get("n_hits", 0),
            "deposit_spanning": d.get("n_spanning", 0),
            "deposit_spanning_anchor_only": d.get("n_spanning_anchor_only", 0),
            "deposit_best_identity": d.get("best_identity", ""),
            "deposit_best_subject": d.get("best_subject", ""),
            "control_hits": c.get("n_hits", 0),
            "control_spanning": c.get("n_spanning", 0)})
    meta = {"n_records": mrna.get("n_records", 0),
            "n_total": mrna.get("n_total", 0),
            "status": mrna.get("status", "unavailable"),
            **_record_lengths(subject)}
    return rows, meta


def _record_lengths(fasta) -> dict:
    """Median and maximum deposited record length.

    The deposit test's answer here is a property of what was deposited, and
    length is the property that decides it: a 2,700-residue channel needs
    8 kb of transcript, and a library of 550 nt single-pass reads cannot
    hold one however many records it has. The number therefore belongs in
    a committed table rather than in the report's prose (D13).
    """
    if not fasta or not Path(fasta).exists():
        return {"median_record_nt": "", "max_record_nt": ""}
    lens, n = [], 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    if not lens:
        return {"median_record_nt": "", "max_record_nt": ""}
    lens.sort()
    m = len(lens)
    med = lens[m // 2] if m % 2 else (lens[m // 2 - 1] + lens[m // 2]) // 2
    return {"median_record_nt": med, "max_record_nt": lens[-1]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="")
    ap.add_argument("--retmax", type=int, default=40000)
    args = ap.parse_args()

    global TAXID, RETMAX
    TAXID, RETMAX = taxids(), args.retmax

    panel, _ = s12_panel.build_panel()
    if args.species:
        want = set(args.species.split(","))
        panel = [s for s in panel if s.short in want]

    resources: list[dict] = []
    probe_rows: list[dict] = []

    for i, sp in enumerate(panel, 1):
        tax = TAXID.get(sp.accession, 0)
        print(f"[atlas] {sp.name} (taxid {tax})", flush=True)
        live("atlas", i - 1, len(panel), sp.name)
        counts = resource_counts(tax, CACHE / "atlas" / sp.accession /
                                 "resources.json") if tax else {}
        loci = locus_models(sp.accession)
        fetched_meta: dict = {}
        for locus in sp.loci:
            key = (locus.cell, locus.contig, locus.start, locus.end,
                   locus.strand)
            rec = loci.get(key)
            if rec is None:
                continue
            rows, meta = probe_locus(sp, locus, rec)
            probe_rows += rows
            fetched_meta = meta or fetched_meta
            span = sum(r["deposit_spanning"] for r in rows)
            print(f"    {locus.cell:6s} {len(rows):>3} probes  "
                  f"{sum(r['deposit_hits'] for r in rows):>5} deposit hits  "
                  f"{span:>4} spanning  "
                  f"(control {sum(r['control_hits'] for r in rows)} hits, "
                  f"{sum(r['control_spanning'] for r in rows)} spanning)")
        resources.append({
            "species": sp.short, "organism": sp.name, "taxid": tax,
            "nuccore_total": counts.get("nuccore_total", ""),
            "nuccore_mrna": counts.get("nuccore_mrna", ""),
            "nuccore_tsa": counts.get("nuccore_tsa", ""),
            "sra_rnaseq": counts.get("sra_rnaseq", ""),
            "mrna_fetched": fetched_meta.get("n_records", 0),
            "mrna_available": fetched_meta.get("n_total", 0),
            "median_record_nt": fetched_meta.get("median_record_nt", ""),
            "max_record_nt": fetched_meta.get("max_record_nt", ""),
            "fetch_status": fetched_meta.get("status", "unavailable")})

    summary = _summarise(probe_rows)
    write_tsv(OUT_DIR / "atlas_resources.tsv", RESOURCE_COLS, resources)
    write_tsv(OUT_DIR / "atlas_probes.tsv", PROBE_COLS, probe_rows)
    write_tsv(OUT_DIR / "atlas_summary.tsv", SUMMARY_COLS, summary)
    live("atlas", len(panel), len(panel), "done")
    print(f"\nwrote atlas_probes.tsv ({len(probe_rows)} probes), "
          f"atlas_summary.tsv ({len(summary)} rows)")


def _summarise(rows: list[dict]) -> list[dict]:
    acc: dict[tuple, dict] = {}
    for r in rows:
        key = (r["species"], r["cell"], r["annotated"])
        a = acc.setdefault(key, {
            "species": r["species"], "organism": r["organism"],
            "cell": r["cell"], "role": r["role"], "annotated": r["annotated"],
            "probes": 0, "deposit_hits": 0, "probes_with_hit": 0,
            "deposit_spanning": 0, "probes_spanned": 0,
            "control_hits": 0, "control_spanning": 0,
            "deposit_status": r["deposit_status"]})
        a["probes"] += 1
        a["deposit_hits"] += r["deposit_hits"]
        a["probes_with_hit"] += 1 if r["deposit_hits"] else 0
        a["deposit_spanning"] += r["deposit_spanning"]
        a["probes_spanned"] += 1 if r["deposit_spanning"] else 0
        a["control_hits"] += r["control_hits"]
        a["control_spanning"] += r["control_spanning"]
    return sorted(acc.values(),
                  key=lambda r: (r["species"], r["cell"], -r["annotated"]))


RESOURCE_COLS = ["species", "organism", "taxid", "nuccore_total",
                 "nuccore_mrna", "nuccore_tsa", "sra_rnaseq", "mrna_fetched",
                 "mrna_available", "median_record_nt", "max_record_nt",
                 "fetch_status"]
PROBE_COLS = ["species", "organism", "accession", "cell", "role",
              "junction_index", "intron_class", "annotated", "splice_class",
              "probe_len", "deposit_status", "deposit_hits",
              "deposit_spanning", "deposit_spanning_anchor_only",
              "deposit_best_identity", "deposit_best_subject", "control_hits",
              "control_spanning"]
SUMMARY_COLS = ["species", "organism", "cell", "role", "annotated", "probes",
                "probes_with_hit", "deposit_hits", "probes_spanned",
                "deposit_spanning", "control_hits", "control_spanning",
                "deposit_status"]


if __name__ == "__main__":
    main()
