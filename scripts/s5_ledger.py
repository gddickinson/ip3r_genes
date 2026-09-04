"""s5_ledger.py — assemble the per-genome ledger from the sweep summaries.

Reads every `summary.json` under `<data_root>/genome_sweep/` and writes the
task's deliverable tables. Nothing is recomputed here: the summaries carry the
calls, this arranges them (D13 — the report renders from the data, so the two
cannot drift).

Outputs -> results/genome_ledger/
  genome_ledger.tsv        one row per genome x class cell, with evidence paths
  genome_ledger_wide.tsv   one row per genome, the four cells as columns
  ledger_status_counts.tsv status x class tally, and the same split by vclass
  control_failures.tsv     genomes where the RyR positive control did not fire
  locus_margins.tsv        every found locus's family and paralog margins
  rescue_regions.tsv       every tblastn rescue region, its attribution and
                           the genes it overlaps — fragment-level evidence,
                           which is what the attribution margin acts on
  ledger_stats.json        totals, coverage of the manifest, run timings

Usage:
  python scripts/s5_ledger.py
  python scripts/s5_ledger.py --pilot     # label the run as a pilot subset
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from src.utils.data_root import get_data_root            # noqa: E402
import s5_sweep_lib as lib                                # noqa: E402
import s5_classify as clf                                # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "genome_ledger"
MANIFEST = PROJECT_ROOT / "results" / "genome_manifest.tsv"
ALL_CELLS = (*lib.CLASSES, lib.CONTROL_CLASS)

#: Statuses in the order a reader should meet them: best evidence first,
#: then the two kinds of "not found", then the two kinds of absence.
STATUS_ORDER = ["found_annotated", "found_unannotated", "found_no_annotation",
                "fragment", "assembly_gap", "tblastn_trace",
                "tblastn_trace_ambiguous", "absent", "no_locus"]


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def load_summaries() -> list[dict]:
    root = get_data_root() / "genome_sweep"
    out = []
    for d in sorted(root.glob("*/")):
        path = d / "summary.json"
        if not (d / ".sweep.done").exists() or not path.exists():
            continue
        try:
            out.append(json.loads(path.read_text()))
        except json.JSONDecodeError:
            print(f"  ! unreadable summary: {path}")
    return out


def cell_row(s: dict, cell_class: str, evidence_dir: str) -> dict:
    cell = s["cells"][cell_class]
    best = (cell.get("loci") or [{}])[0]
    ag = best.get("annot_gene") or {}
    rescue = cell.get("rescue") or {}
    regions = cell.get("rescue_regions") or []
    return {
        "accession": s["accession"], "organism": s["organism"],
        "vclass": s.get("vclass", ""), "vorder": s.get("vorder", ""),
        "class": cell_class,
        "is_control": int(cell_class == lib.CONTROL_CLASS),
        "status": cell["status"],
        "n_loci": cell["n_loci"], "n_primary": cell.get("n_primary", 0),
        "n_secondary": cell.get("n_secondary", 0),
        "best_coverage": cell.get("best_coverage", ""),
        "best_identity": cell.get("best_identity", ""),
        "family_margin": cell.get("best_family_margin", ""),
        "paralog_margin": cell.get("best_paralog_margin", ""),
        "assigned_via": cell.get("best_assigned_via", ""),
        "contig": best.get("contig", ""), "start": best.get("start", ""),
        "end": best.get("end", ""), "strand": best.get("strand", ""),
        "bait": best.get("bait", ""),
        "frameshifts": best.get("frameshifts", ""),
        "stop_codons": best.get("stop_codons", ""),
        "contig_edge": int(bool(best.get("contig_edge"))),
        "longest_n_run": best.get("longest_n_run", ""),
        "annot_gene": ag.get("name", ""),
        "annot_paralog": cell.get("annot_paralog", ""),
        "annot_paralog_matches": int(bool(best.get("annot_paralog_matches"))),
        "rescue_hsps": rescue.get("n_hsps", ""),
        "rescue_best_evalue": rescue.get("best_evalue", ""),
        "rescue_regions": len(regions),
        "rescue_ambiguous": len(cell.get("ambiguous_regions") or []),
        "rescue_elsewhere": len(cell.get("attributed_elsewhere") or []),
        "annotated_assembly": s.get("annotated", ""),
        "assembly_level": s.get("assembly_level", ""),
        "genome_bp": s.get("genome_bp", ""),
        "contig_n50": s.get("contig_n50", ""),
        "max_intron": s.get("max_intron", ""),
        "max_intron_capped": int(bool(s.get("max_intron_capped"))),
        "contig_spans_gene": int(lib.spans_a_gene(int(s.get("contig_n50") or 0))),
        "reasons": s.get("reasons", ""),
        "evidence": f"{evidence_dir}/{s['accession']}/summary.json",
    }


def margin_rows(s: dict) -> list[dict]:
    """Every found locus's margins, with whether the annotation agrees.

    The rows where `annot_paralog` is non-empty are loci whose paralog
    identity is established independently of the alignment, so they are what
    `s5_rescue.ATTRIBUTION_REL_MARGIN` can be calibrated against.
    """
    out = []
    for cell_class in ALL_CELLS:
        cell = s["cells"][cell_class]
        for d in cell.get("loci", []):
            out.append({
                "accession": s["accession"], "organism": s["organism"],
                "vclass": s.get("vclass", ""), "class": cell_class,
                "status": cell["status"],
                "contig": d.get("contig", ""), "start": d.get("start", ""),
                "coverage": d.get("coverage", ""),
                "identity": d.get("identity", ""),
                "family_margin": d.get("family_margin", ""),
                "paralog_margin": d.get("paralog_margin", ""),
                "assigned_via": d.get("assigned_via", ""),
                "bait": d.get("bait", ""),
                "annot_gene": (d.get("annot_gene") or {}).get("name", ""),
                "annot_paralog": d.get("annot_paralog", ""),
                "annot_agrees": int(bool(d.get("annot_paralog_matches"))),
            })
    return out


def region_rows(s: dict) -> list[dict]:
    """Every rescue region, with the attribution it received and why.

    S5a could only calibrate the attribution margin on *complete* loci, which
    bounds a fragment's separation from above without measuring it. These are
    the fragments themselves. A region overlapping a gene the assembly names
    for a paralog carries an identity established independently of the bait
    scores — the same trick, one evidence level down, and the calibration
    S5a deferred to here.
    """
    out = []
    for cell_class in ALL_CELLS:
        cell = s["cells"][cell_class]
        for kind in ("rescue_regions", "ambiguous_regions",
                     "attributed_elsewhere"):
            for r in cell.get(kind, []):
                genes = r.get("genes") or []
                named = ""
                for g in genes:
                    p = clf.name_paralog(g.get("name", ""))
                    if p:
                        named = p
                        break
                out.append({
                    "accession": s["accession"], "organism": s["organism"],
                    "vclass": s.get("vclass", ""),
                    "cell": cell_class, "cell_status": cell["status"],
                    "region_kind": kind,
                    "contig": r.get("contig", ""), "start": r.get("start", ""),
                    "end": r.get("end", ""),
                    "n_hsps": r.get("n_hsps", ""),
                    "aligned_aa": r.get("aligned_aa", ""),
                    "best_evalue": r.get("best_evalue", ""),
                    "best_pident": r.get("best_pident", ""),
                    "assigned_clade": r.get("assigned_clade", ""),
                    "attribution_rel_margin": r.get("attribution_rel_margin", ""),
                    "family_call": r.get("family_call", ""),
                    "family_rel_margin": r.get("family_rel_margin", ""),
                    "clade_bits": ";".join(
                        f"{k}={v}" for k, v in (r.get("clade_bits") or {}).items()),
                    "overlapping_genes": ";".join(
                        g.get("name", "") for g in genes[:4]),
                    "annot_paralog": named,
                    "annot_agrees": int(bool(named)
                                        and named == r.get("assigned_clade")),
                })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", action="store_true",
                    help="label the run a pilot subset of the manifest")
    args = ap.parse_args()

    summaries = load_summaries()
    if not summaries:
        raise SystemExit("no completed sweeps found — run scripts/s5_run_sweep.py")
    evidence_dir = str(get_data_root() / "genome_sweep")

    rows, margins, wide, regions = [], [], [], []
    for s in summaries:
        for cell_class in ALL_CELLS:
            rows.append(cell_row(s, cell_class, evidence_dir))
        margins.extend(margin_rows(s))
        regions.extend(region_rows(s))
        w = {"accession": s["accession"], "organism": s["organism"],
             "vclass": s.get("vclass", ""), "vorder": s.get("vorder", ""),
             "annotated": s.get("annotated", ""),
             "assembly_level": s.get("assembly_level", ""),
             "genome_bp": s.get("genome_bp", ""),
             "contig_n50": s.get("contig_n50", ""),
             "contig_spans_gene": int(lib.spans_a_gene(
                 int(s.get("contig_n50") or 0))),
             "reasons": s.get("reasons", ""),
             "n_loci": s.get("n_loci", 0),
             "other_loci": len(s.get("other_loci") or []),
             "novel_models": s.get("novel_model_count", 0),
             "control_ok": int(bool(s.get("control_ok"))),
             "sweep_s": s.get("total_s", "")}
        for cell_class in ALL_CELLS:
            cell = s["cells"][cell_class]
            w[cell_class] = cell["status"]
            w[f"{cell_class}_n"] = cell["n_loci"]
        n_found = sum(1 for c in lib.CLASSES
                      if s["cells"][c]["status"].startswith("found"))
        w["itpr_found"] = n_found
        wide.append(w)

    cell_cols = ["accession", "organism", "vclass", "vorder", "class",
                 "is_control", "status", "n_loci", "n_primary", "n_secondary",
                 "best_coverage", "best_identity", "family_margin",
                 "paralog_margin", "assigned_via", "contig_spans_gene",
                 "contig", "start", "end",
                 "strand", "bait", "frameshifts", "stop_codons", "contig_edge",
                 "longest_n_run", "annot_gene", "annot_paralog",
                 "annot_paralog_matches", "rescue_hsps", "rescue_best_evalue",
                 "rescue_regions", "rescue_ambiguous", "rescue_elsewhere",
                 "annotated_assembly", "assembly_level", "genome_bp",
                 "contig_n50", "max_intron", "max_intron_capped", "reasons",
                 "evidence"]
    write_tsv(OUT_DIR / "genome_ledger.tsv", cell_cols,
              sorted(rows, key=lambda r: (r["accession"], r["class"])))

    wide_cols = (["accession", "organism", "vclass", "vorder", "annotated",
                  "assembly_level", "genome_bp", "contig_n50",
                  "contig_spans_gene", "reasons"]
                 + [c for cell in ALL_CELLS for c in (cell, f"{cell}_n")]
                 + ["itpr_found", "n_loci", "other_loci", "novel_models",
                    "control_ok", "sweep_s"])
    write_tsv(OUT_DIR / "genome_ledger_wide.tsv", wide_cols,
              sorted(wide, key=lambda r: r["accession"]))

    write_tsv(OUT_DIR / "locus_margins.tsv",
              ["accession", "organism", "vclass", "class", "status", "contig",
               "start", "coverage", "identity", "family_margin",
               "paralog_margin", "assigned_via", "bait", "annot_gene",
               "annot_paralog", "annot_agrees"], margins)

    write_tsv(OUT_DIR / "rescue_regions.tsv",
              ["accession", "organism", "vclass", "cell", "cell_status",
               "region_kind", "contig", "start", "end", "n_hsps",
               "aligned_aa", "best_evalue", "best_pident", "assigned_clade",
               "attribution_rel_margin", "family_call", "family_rel_margin",
               "clade_bits", "overlapping_genes", "annot_paralog",
               "annot_agrees"], regions)

    by_status: dict[tuple, int] = Counter()
    by_class_vclass: dict[tuple, int] = Counter()
    for r in rows:
        by_status[(r["class"], r["status"])] += 1
        by_class_vclass[(r["vclass"], r["class"], r["status"])] += 1
    tally = [{"class": c, "status": st, "n": n}
             for (c, st), n in sorted(by_status.items(),
                                      key=lambda kv: (ALL_CELLS.index(kv[0][0]),
                                                      STATUS_ORDER.index(kv[0][1])
                                                      if kv[0][1] in STATUS_ORDER
                                                      else 99))]
    tally += [{"vclass": v, "class": c, "status": st, "n": n}
              for (v, c, st), n in sorted(by_class_vclass.items())]
    write_tsv(OUT_DIR / "ledger_status_counts.tsv",
              ["vclass", "class", "status", "n"], tally)

    failures = [w for w in wide if not w["control_ok"]]
    write_tsv(OUT_DIR / "control_failures.tsv",
              ["accession", "organism", "vclass", "annotated",
               "assembly_level", "contig_n50", "RYR", "RYR_n", "reasons"],
              failures)

    manifest = read_tsv(MANIFEST) if MANIFEST.exists() else []
    swept = {s["accession"] for s in summaries}
    stats = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "is_pilot": bool(args.pilot),
        "genomes_swept": len(summaries),
        "genomes_in_manifest": len(manifest),
        "manifest_remaining": len([m for m in manifest
                                   if m["accession"] not in swept]),
        "cells": len(rows), "loci_recorded": len(margins),
        "rescue_regions": len(regions),
        "control_failures": len(failures),
        "itpr_span_stats": lib.itpr_span_stats(),
        "genomes_below_contiguity_bar": sum(
            1 for w in wide if not w["contig_spans_gene"]),
        "novel_models": sum(w["novel_models"] for w in wide),
        "status_counts": {f"{c}:{st}": n for (c, st), n in by_status.items()},
        "itpr_found_per_genome": dict(Counter(w["itpr_found"] for w in wide)),
        "median_sweep_s": sorted(s["total_s"] for s in summaries)[
            len(summaries) // 2],
        "total_sweep_s": round(sum(s["total_s"] for s in summaries), 1),
        "total_genome_bp": sum(int(s.get("genome_bp") or 0) for s in summaries),
    }
    (OUT_DIR / "ledger_stats.json").write_text(json.dumps(stats, indent=1) + "\n")

    label = "pilot ledger" if args.pilot else "genome ledger"
    print(f"{label}: {len(summaries)} genomes, {len(rows)} cells, "
          f"{len(margins)} loci")
    for cell_class in ALL_CELLS:
        counts = {st: n for (c, st), n in by_status.items() if c == cell_class}
        ordered = sorted(counts.items(),
                         key=lambda kv: STATUS_ORDER.index(kv[0])
                         if kv[0] in STATUS_ORDER else 99)
        tag = " (control)" if cell_class == lib.CONTROL_CLASS else ""
        print(f"  {cell_class}{tag}: "
              + ", ".join(f"{st} {n}" for st, n in ordered))
    below = [w for w in wide if not w["contig_spans_gene"]]
    if below:
        span = lib.itpr_span_stats()["median"]
        print(f"  {len(below)} genome(s) below D4's contiguity bar "
              f"(contig N50 < the median measured ITPR span, {span:,} bp) — "
              "no absence claim may rest on them: "
              + ", ".join(w["accession"] for w in below[:5]))
    if failures:
        print(f"  ** {len(failures)} RyR control failure(s): "
              + ", ".join(f["accession"] for f in failures[:5]))
    print(f"wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
