"""S10 step 1 — rank every recovered locus by annotation loss, pick two cases.

Reads the sweep's per-genome `summary.json` (S5's own record of what it found
and what the annotation said about it) plus `genes_slim.tsv` (the annotation's
gene table, archived beside it) and writes three tables:

* `case_ranking.tsv` — every locus S5 recovered, with its loss, its failure
  mode, whether it passed E1–E5 and which rule it failed on. This is the
  denominator: it is what makes "the two worst" a statement about a
  distribution rather than about two loci someone liked.
* `annotation_depth.tsv` — the E5 control, per genome: how many annotated
  genes are longer than each locus, and the annotation's own gene-length
  distribution. A genome excluded by E5 is visible here with the number that
  excluded it.
* `cases.tsv` — the two selected cases with everything the later steps need.

Nothing here recomputes a call. `annotation_loss` comes from the `frac_cds`
values S5 measured at sweep time, so the ranking cannot disagree with the
ledger it is derived from (D13).
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_case_spec as spec                                    # noqa: E402
from s10_case_spec import (FRAG_MIN_FRAC_CDS, control_strength,  # noqa: E402
                           eligible, mode_of, select_cases, severity_key)
from s5_calibration import spans_a_gene                          # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "annotation_bugs"
CONTROL_CLASS = "RYR"


def _data_root() -> Path:
    sys.path.insert(0, str(ROOT))
    from src.utils.data_root import require_data_root
    return require_data_root()


def gene_length_profile(genes_slim: Path) -> dict:
    """Sorted protein-coding gene spans, and the summary E5 is read off.

    Only `protein_coding` genes count. A pseudogene span would inflate the
    control with exactly the models this task is disputing — *Nibea
    albiflora*'s longest family models are all pseudogenes, and crediting the
    annotation with being able to build a 53 kb gene on the strength of one it
    declined to translate would be circular.
    """
    spans: list[int] = []
    with open(genes_slim) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[6] != "protein_coding":
                continue
            spans.append(int(f[2]) - int(f[1]) + 1)
    spans.sort()
    n = len(spans)
    return {"spans": spans, "n_protein_coding": n,
            "median_span": spans[n // 2] if n else 0,
            "p90_span": spans[int(n * 0.9)] if n else 0,
            "max_span": spans[-1] if n else 0}


def locus_records(sweep_root: Path) -> tuple[list[dict], list[dict]]:
    """One record per recovered ITPR locus, plus the per-genome E5 control."""
    records, depth = [], []
    for summary in sorted(sweep_root.glob("*/summary.json")):
        d = json.loads(summary.read_text())
        acc = d["accession"]
        slim = summary.parent / "genes_slim.tsv"
        prof = (gene_length_profile(slim) if slim.exists() else
                {"spans": [], "n_protein_coding": 0, "median_span": 0,
                 "p90_span": 0, "max_span": 0})
        spans = prof["spans"]
        for cell_name, cell in d["cells"].items():
            if cell_name == CONTROL_CLASS:
                continue
            loci = cell.get("loci") or []
            if not loci:
                continue
            best = loci[0]
            span = best["end"] - best["start"] + 1
            og = best.get("overlapping_genes") or []
            frags = [g for g in og if g.get("frac_cds", 0) > FRAG_MIN_FRAC_CDS]
            loss = round(1.0 - max([g["frac_cds"] for g in og], default=0.0), 4)
            n_longer = len(spans) - bisect.bisect_left(spans, span)
            rec = {
                "accession": acc, "organism": d["organism"],
                "vclass": d.get("vclass", ""), "vorder": d.get("vorder", ""),
                "cell": cell_name, "cell_status": cell["status"],
                "contig": best["contig"], "start": best["start"],
                "end": best["end"], "strand": best["strand"],
                "span": span, "mp_id": best["mp_id"], "bait": best["bait"],
                "coverage": best["coverage"], "identity": best["identity"],
                "frameshifts": best.get("frameshifts", 0),
                "stop_codons": best.get("stop_codons", 0),
                "cds_bp": best.get("cds_bp", 0),
                "contig_edge": bool(best.get("contig_edge")),
                "n_gap": bool(best.get("n_gap")),
                # D4's bar, taken from S5's own function rather than
                # re-expressed here, so E3 and the ledger cannot disagree.
                "contig_spans_gene": spans_a_gene(int(d.get("contig_n50") or 0)),
                "annotated": d.get("annotated") == "Y",
                "has_annotation_index": bool(d.get("has_annotation_index")),
                "contig_n50": int(d.get("contig_n50") or 0),
                "assembly_level": d.get("assembly_level", ""),
                "n_genes_longer": n_longer,
                "annot_max_gene_span": prof["max_span"],
                "annot_median_gene_span": prof["median_span"],
                "loss": loss, "n_fragments": len(frags),
                "best_frac_cds": round(1.0 - loss, 4),
                "sum_frac_cds": round(sum(g["frac_cds"] for g in frags), 4),
                "fragment_names": ";".join(g["name"] for g in sorted(
                    frags, key=lambda x: -x["frac_cds"])),
                "control_strength": control_strength(d["cells"], cell_name),
                "reasons": d.get("reasons", ""),
            }
            rec["headroom"] = round(rec["contig_n50"] / max(1, span), 2)
            ok, failed = eligible(rec)
            rec["eligible"], rec["failed_rule"] = ok, failed
            rec["mode"] = mode_of(rec["n_fragments"], loss) if ok else ""
            records.append(rec)
        depth.append({
            "accession": acc, "organism": d["organism"],
            "annotated": d.get("annotated", ""),
            "n_protein_coding": prof["n_protein_coding"],
            "median_gene_span": prof["median_span"],
            "p90_gene_span": prof["p90_span"], "max_gene_span": prof["max_span"],
            "n_itpr_loci": sum(1 for k, c in d["cells"].items()
                               if k != CONTROL_CLASS and c.get("loci")),
            "n_found_annotated": sum(
                1 for c in d["cells"].values()
                if c.get("status") == "found_annotated"),
        })
    records.sort(key=severity_key)
    depth.sort(key=lambda r: r["accession"])
    return records, depth


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


RANK_COLS = ["accession", "organism", "vclass", "vorder", "cell", "cell_status",
             "contig", "start", "end", "strand", "span", "mp_id", "bait",
             "coverage", "identity", "frameshifts", "stop_codons", "cds_bp",
             "contig_spans_gene", "contig_n50", "headroom", "assembly_level",
             "n_genes_longer", "annot_max_gene_span", "annot_median_gene_span",
             "loss", "best_frac_cds", "sum_frac_cds", "n_fragments",
             "fragment_names", "control_strength", "eligible", "failed_rule",
             "mode", "reasons"]

DEPTH_COLS = ["accession", "organism", "annotated", "n_protein_coding",
              "median_gene_span", "p90_gene_span", "max_gene_span",
              "n_itpr_loci", "n_found_annotated"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    sweep_root = _data_root() / "genome_sweep"
    records, depth = locus_records(sweep_root)
    cases = select_cases(records)

    write_tsv(args.out / "case_ranking.tsv", records, RANK_COLS)
    write_tsv(args.out / "annotation_depth.tsv", depth, DEPTH_COLS)
    write_tsv(args.out / "cases.tsv", cases,
              ["case_id", "selected_as", "selection_rank_in_mode"] + RANK_COLS)

    elig = [r for r in records if r["eligible"]]
    stats = {
        "n_loci_recovered": len(records),
        "n_eligible": len(elig),
        "n_excluded_by_rule": {
            rule: sum(1 for r in records if r["failed_rule"] == rule)
            for rule, _ in spec.ELIGIBILITY},
        "median_loss_eligible": (
            sorted(r["loss"] for r in elig)[len(elig) // 2] if elig else None),
        "n_loss_zero": sum(1 for r in elig if r["loss"] == 0.0),
        "n_by_mode": {m: sum(1 for r in elig if r["mode"] == m)
                      for m in spec.FAILURE_MODES},
        "n_failures": sum(1 for r in elig if r["mode"]),
        "cases": [{k: c[k] for k in
                   ("case_id", "selected_as", "accession", "organism", "cell",
                    "loss", "n_fragments", "control_strength", "headroom",
                    "n_genes_longer")} for c in cases],
    }
    (args.out / "selection_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n")

    print(f"[s10] {len(records)} recovered loci, {len(elig)} eligible, "
          f"{stats['n_failures']} failures")
    for rule, _ in spec.ELIGIBILITY:
        n = stats["n_excluded_by_rule"][rule]
        if n:
            print(f"       excluded by {rule}: {n}")
    for c in cases:
        print(f"  {c['case_id']}: {c['organism']} {c['cell']} "
              f"({c['selected_as']}) loss={c['loss']} "
              f"frags={c['n_fragments']} control={c['control_strength']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
