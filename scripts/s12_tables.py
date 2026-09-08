"""S12 step 4 — the committed tables, written once from the per-run counts.

Everything the report and the figures render comes from here (D13); nothing
downstream re-parses a counts file or re-reads a SAM.

The detection rule is the brief's, and all three parts are needed:

  * **>= MIN_JUNCTION_READS junction-spanning reads.**  The load-bearing
    one.  A junction read crosses a splice point, so genomic-DNA carryover
    in the library cannot produce it.
  * **>= MIN_READS reads.**  A single junction read is a single read.
  * **more reads than this run's own decoy.**  The decoy is the same
    sequence reversed — identical length and base composition, no
    homology — so it measures what this reference collects by accident in
    this library, not in some other one.

Two things are computed that the PIEZO port did not have, and both come
from S10's hand-off.  `junction_support.tsv` scores **each junction
separately**, split by whether the annotation models it, because "the gene
is transcribed" and "the junctions no annotated model spans are spliced"
are different claims and only the second answers S10.  And
`annotation_gap_coverage.tsv` measures how much of each locus's read
coverage falls where no annotated CDS block lies — the read-level form of
the annotation-loss number S10 ranked on.

Outputs (committed, under results/expression/):
    run_metrics.tsv  expression_by_run.tsv  expression_by_locus.tsv
    expression_by_tissue.tsv  junction_support.tsv
    annotation_gap_coverage.tsv  coverage_profile.tsv  expression_stats.json
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s12_lib import (  # noqa: E402
    DATA, OUT_DIR, sha256, tool_version, write_tsv,
)
from s12_quantify import COUNTS, DEFAULT_ANCHOR, MIN_MAPQ  # noqa: E402

MIN_JUNCTION_READS = 2
MIN_READS = 5


def _crossmap_worst() -> float | str:
    """The worst measured cross-mapping rate, from `crossmap_control.tsv`.

    Recorded in the stats file because the reference is a closed set and
    that is the number which says whether the closure is safe. Absent
    rather than assumed when the control has not been run.
    """
    rows = _rows(OUT_DIR / "crossmap_control.tsv")
    if not rows:
        return "not measured"
    try:
        return max(float(r.get("cross_rate") or 0) for r in rows)
    except ValueError:
        return "not measured"


def _self_test_status() -> str:
    """Whether the negative controls pass, recorded in the stats file.

    A stats file that does not say the self-test ran leaves a reader unable
    to tell a passing build from one where the checks were skipped.
    """
    import s12_test_expression
    try:
        return "pass" if s12_test_expression.main() == 0 else "fail"
    except Exception as exc:                          # noqa: BLE001
        return f"error:{type(exc).__name__}"


def _rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_counts(species: str, run: str) -> tuple[int, dict[str, dict]]:
    path = COUNTS / f"{species}.{run}.tsv"
    if not path.exists():
        return 0, {}
    lib = 0
    out: dict[str, dict] = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("# library_reads"):
                lib = int(line.split("\t")[1])
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            r = dict(zip(header, f))
            out[r["seq"]] = {"reads": int(r["reads"]),
                             "junction_reads": int(r["junction_reads"]),
                             "covered_bases": int(r["covered_bases"]),
                             "ref_len": int(r["ref_len"])}
    return lib, out


def detect(reads: int, junction_reads: int, decoy: int) -> tuple[int, str]:
    """The three-part call, with the part that failed named."""
    fails = []
    if junction_reads < MIN_JUNCTION_READS:
        fails.append(f"junction_reads<{MIN_JUNCTION_READS}")
    if reads < MIN_READS:
        fails.append(f"reads<{MIN_READS}")
    if reads <= decoy:
        fails.append("not_above_decoy")
    return (0, ";".join(fails)) if fails else (1, "")


# ---------------------------------------------------------------------------

def build(log=print) -> dict:
    refs = _rows(OUT_DIR / "reference_table.tsv")
    if not refs:
        raise SystemExit("no reference_table.tsv — run s12_refs.py first")
    runs = _rows(OUT_DIR / "runs_selected.tsv")
    ref_by = {(r["species"], r["seq"]): r for r in refs}
    real = [r for r in refs if r["role"] != "decoy"]

    metrics: list[dict] = []
    by_run: list[dict] = []
    junc_total: dict[tuple, dict] = {}
    cov_total: dict[tuple, list[int]] = {}
    cov_meta: dict[tuple, tuple[int, int]] = {}
    missing = 0

    for r in runs:
        sp, run = r["species"], r["run"]
        lib, counts = load_counts(sp, run)
        if not counts:
            missing += 1
            continue
        metrics.append({"species": sp, "organism": r["organism"], "run": run,
                        "tissue": r["tissue"], "tissue_from": r["tissue_from"],
                        "study": r["study"], "library_reads": lib,
                        "spots": r["spots"], "read_len": r["read_len"],
                        "layout": r["layout"],
                        "aligned": sum(c["reads"] for c in counts.values())})
        for ref in real:
            if ref["species"] != sp:
                continue
            seq = ref["seq"]
            c = counts.get(seq, {"reads": 0, "junction_reads": 0,
                                 "covered_bases": 0,
                                 "ref_len": int(ref["length"] or 0)})
            d = counts.get(f"decoy_{seq}", {"reads": 0, "covered_bases": 0})
            called, why = detect(c["reads"], c["junction_reads"], d["reads"])
            by_run.append({
                "species": sp, "organism": r["organism"], "run": run,
                "tissue": r["tissue"], "seq": seq, "role": ref["role"],
                "cell": ref.get("cell", ""),
                "failure_mode": ref.get("failure_mode", ""),
                "reads": c["reads"], "junction_reads": c["junction_reads"],
                "covered_bases": c["covered_bases"],
                "ref_len": c["ref_len"], "decoy_reads": d["reads"],
                "decoy_covered_bases": d.get("covered_bases", 0),
                "library_reads": lib,
                "reads_per_million": round(1e6 * c["reads"] / lib, 3) if lib else 0,
                "detected": called, "failed_criteria": why})

        for j in _rows(COUNTS / f"{sp}.{run}.junctions.tsv"):
            key = (sp, j["seq"], int(j["junction_index"]))
            rec = junc_total.setdefault(key, {
                "species": sp, "seq": j["seq"],
                "junction_index": int(j["junction_index"]),
                "cds_offset": int(j["cds_offset"]),
                "annotated": int(j["annotated"]), "class": j["class"],
                "reads": 0, "runs_with_reads": 0, "runs_tested": 0})
            n = int(j["reads"])
            rec["reads"] += n
            rec["runs_tested"] += 1
            rec["runs_with_reads"] += 1 if n else 0

        for c in _rows(COUNTS / f"{sp}.{run}.coverage.tsv"):
            key = (sp, c["seq"])
            bins = [int(x) for x in c["coverage"].split(",")]
            acc = cov_total.setdefault(key, [0] * len(bins))
            for i, v in enumerate(bins):
                acc[i] += v
            cov_meta[key] = (int(c["n_bins"]), int(c["bin_width"]))

    if missing:
        log(f"  !! {missing} selected run(s) have no counts file yet")

    junction_rows = sorted(junc_total.values(),
                           key=lambda r: (r["species"], r["seq"],
                                          r["junction_index"]))
    locus_rows = _by_locus(by_run, junction_rows, ref_by)
    tissue_rows = _by_tissue(by_run)
    gap_rows = _gap_coverage(cov_total, cov_meta, ref_by)
    cov_rows = [{"species": k[0], "seq": k[1], "n_bins": cov_meta[k][0],
                 "bin_width": cov_meta[k][1],
                 "coverage": ",".join(str(x) for x in v)}
                for k, v in sorted(cov_total.items())]

    written = {}
    for name, cols, rows in (
            ("run_metrics.tsv", METRIC_COLS, metrics),
            ("expression_by_run.tsv", RUN_COLS, by_run),
            ("expression_by_locus.tsv", LOCUS_COLS, locus_rows),
            ("expression_by_tissue.tsv", TISSUE_COLS, tissue_rows),
            ("junction_support.tsv", JUNC_COLS, junction_rows),
            ("annotation_gap_coverage.tsv", GAP_COLS, gap_rows),
            ("coverage_profile.tsv", COV_COLS, cov_rows)):
        p = write_tsv(OUT_DIR / name, cols, rows)
        written[name] = {"rows": len(rows), "sha256": sha256(p)}
        log(f"  {name}: {len(rows)} rows")

    # Hash *every* committed table in the directory, not only the ones this
    # module writes: the references, the run selection, the cross-mapping
    # control and the deposit cross-check are written by other stages, and
    # a provenance record that covers half the directory is worse than
    # none because it looks complete.
    for p in sorted(OUT_DIR.glob("*.tsv")):
        written.setdefault(p.name, {
            "rows": max(0, sum(1 for _ in open(p)) - 1),
            "sha256": sha256(p)})

    stats = {
        "task": "S12", "written_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "self_test": _self_test_status(),
        "tools": {name: tool_version(name)
                  for name in ("hisat2", "fastq-dump", "miniprot", "blastn")},
        "parameters": {
            "min_junction_reads": MIN_JUNCTION_READS, "min_reads": MIN_READS,
            "junction_anchor_nt": DEFAULT_ANCHOR, "min_mapq": MIN_MAPQ,
            "decoy": "reversed CDS (length- and composition-matched)"},
        "crossmap_worst_rate": _crossmap_worst(),
        "runs_quantified": len(metrics),
        "runs_selected": len(runs),
        "runs_missing_counts": missing,
        "tables": written,
    }
    (OUT_DIR / "expression_stats.json").write_text(json.dumps(stats, indent=1))
    return stats


def _by_locus(by_run, junction_rows, ref_by) -> list[dict]:
    """One row per species x reference sequence, pooled over runs."""
    acc: dict[tuple, dict] = {}
    for r in by_run:
        key = (r["species"], r["seq"])
        a = acc.setdefault(key, {
            "species": r["species"], "organism": r["organism"],
            "seq": r["seq"], "role": r["role"], "cell": r["cell"],
            "failure_mode": r["failure_mode"], "runs": 0, "runs_detected": 0,
            "reads": 0, "junction_reads": 0, "decoy_reads": 0,
            "decoy_covered_bases": 0,
            "library_reads": 0, "tissues_detected": set(),
            "max_covered_bases": 0, "ref_len": r["ref_len"]})
        a["runs"] += 1
        a["runs_detected"] += r["detected"]
        a["reads"] += r["reads"]
        a["junction_reads"] += r["junction_reads"]
        a["decoy_reads"] += r["decoy_reads"]
        a["decoy_covered_bases"] = max(a["decoy_covered_bases"],
                                       r["decoy_covered_bases"])
        a["library_reads"] += r["library_reads"]
        a["max_covered_bases"] = max(a["max_covered_bases"], r["covered_bases"])
        if r["detected"] and r["tissue"] != "unknown":
            a["tissues_detected"].add(r["tissue"])

    junc: dict[tuple, dict] = {}
    for j in junction_rows:
        k = (j["species"], j["seq"])
        d = junc.setdefault(k, {"a_tot": 0, "a_hit": 0, "u_tot": 0, "u_hit": 0})
        side = "a" if j["annotated"] else "u"
        d[f"{side}_tot"] += 1
        d[f"{side}_hit"] += 1 if j["reads"] else 0

    out = []
    for key, a in sorted(acc.items()):
        j = junc.get(key, {"a_tot": 0, "a_hit": 0, "u_tot": 0, "u_hit": 0})
        ref = ref_by.get(key, {})
        out.append({
            **{k: v for k, v in a.items() if k != "tissues_detected"},
            "tissues_detected": ";".join(sorted(a["tissues_detected"])),
            "n_tissues_detected": len(a["tissues_detected"]),
            "reads_per_million": round(
                1e6 * a["reads"] / a["library_reads"], 3)
            if a["library_reads"] else 0,
            "coverage_frac": round(a["max_covered_bases"] / a["ref_len"], 4)
            if a["ref_len"] else 0,
            "junctions_annotated": j["a_tot"],
            "junctions_annotated_covered": j["a_hit"],
            "junctions_unannotated": j["u_tot"],
            "junctions_unannotated_covered": j["u_hit"],
            "annotation_loss": ref.get("annotation_loss", ""),
            "cell_status": ref.get("cell_status", ""),
            "annot_gene": ref.get("annot_gene", ""),
            "annot_biotype": ref.get("annot_biotype", ""),
            "annot_frac_cds": ref.get("annot_frac_cds", ""),
            "genome_copies": ref.get("genome_copies", "")})
    return out


def _by_tissue(by_run) -> list[dict]:
    acc: dict[tuple, dict] = {}
    for r in by_run:
        key = (r["species"], r["seq"], r["tissue"])
        a = acc.setdefault(key, {
            "species": r["species"], "organism": r["organism"],
            "seq": r["seq"], "role": r["role"], "cell": r["cell"],
            "tissue": r["tissue"], "runs": 0, "runs_detected": 0, "reads": 0,
            "junction_reads": 0, "decoy_reads": 0, "library_reads": 0})
        a["runs"] += 1
        a["runs_detected"] += r["detected"]
        for f in ("reads", "junction_reads", "decoy_reads", "library_reads"):
            a[f] += r[f]
    out = []
    for a in acc.values():
        called, why = detect(a["reads"], a["junction_reads"], a["decoy_reads"])
        out.append({**a, "detected": called, "failed_criteria": why,
                    "reads_per_million": round(
                        1e6 * a["reads"] / a["library_reads"], 3)
                    if a["library_reads"] else 0})
    return sorted(out, key=lambda r: (r["species"], r["seq"], r["tissue"]))


def _gap_coverage(cov_total, cov_meta, ref_by) -> list[dict]:
    """How much read coverage falls where no annotated CDS block lies.

    Bins are assigned to the annotated or unannotated side by the exon they
    contain, read from the junction table's own exon classification: a bin
    between two junctions the annotation models is annotated territory.
    This is the read-level form of S10's annotation-loss measurement, and it
    is the one number in S12 that can be compared directly against it.
    """
    junc_by_seq: dict[tuple, list[dict]] = {}
    for j in _rows(OUT_DIR / "junctions.tsv"):
        junc_by_seq.setdefault((j["species"], j["seq"]), []).append(j)
    out = []
    for key, bins in sorted(cov_total.items()):
        js = sorted(junc_by_seq.get(key, []), key=lambda x: int(x["cds_offset"]))
        if not js:
            continue
        n_bins, width = cov_meta[key]
        ann_edges = [(int(a["cds_offset"]), int(b["cds_offset"]))
                     for a, b in zip(js, js[1:])
                     if int(a["annotated"]) and int(b["annotated"])]
        ann_reads = un_reads = 0
        ann_bins = un_bins = ann_hit = un_hit = 0
        for i, v in enumerate(bins):
            lo, hi = i * width + 1, (i + 1) * width
            inside = any(s < lo and hi <= e for s, e in ann_edges)
            if inside:
                ann_bins += 1
                ann_reads += v
                ann_hit += 1 if v else 0
            else:
                un_bins += 1
                un_reads += v
                un_hit += 1 if v else 0
        ref = ref_by.get(key, {})
        total = ann_reads + un_reads
        out.append({
            "species": key[0], "organism": ref.get("organism", ""),
            "seq": key[1], "cell": ref.get("cell", ""),
            "role": ref.get("role", ""),
            "annotation_loss": ref.get("annotation_loss", ""),
            "bins_annotated": ann_bins, "bins_unannotated": un_bins,
            "bins_annotated_with_reads": ann_hit,
            "bins_unannotated_with_reads": un_hit,
            "covered_annotated": ann_reads, "covered_unannotated": un_reads,
            "frac_coverage_unannotated": round(un_reads / total, 4)
            if total else 0})
    return out


METRIC_COLS = ["species", "organism", "run", "tissue", "tissue_from", "study",
               "library_reads", "spots", "read_len", "layout", "aligned"]
RUN_COLS = ["species", "organism", "run", "tissue", "seq", "role", "cell",
            "failure_mode", "reads", "junction_reads", "covered_bases",
            "ref_len", "decoy_reads", "decoy_covered_bases", "library_reads",
            "reads_per_million", "detected", "failed_criteria"]
LOCUS_COLS = ["species", "organism", "seq", "role", "cell", "failure_mode",
              "cell_status", "annot_gene", "annot_biotype", "annot_frac_cds",
              "genome_copies", "annotation_loss", "runs",
              "runs_detected", "reads", "junction_reads", "decoy_reads",
              "decoy_covered_bases", "reads_per_million",
              "max_covered_bases", "ref_len",
              "coverage_frac", "junctions_annotated",
              "junctions_annotated_covered", "junctions_unannotated",
              "junctions_unannotated_covered", "n_tissues_detected",
              "tissues_detected", "library_reads"]
TISSUE_COLS = ["species", "organism", "seq", "role", "cell", "tissue", "runs",
               "runs_detected", "reads", "junction_reads", "decoy_reads",
               "reads_per_million", "detected", "failed_criteria",
               "library_reads"]
JUNC_COLS = ["species", "seq", "junction_index", "cds_offset", "annotated",
             "class", "reads", "runs_with_reads", "runs_tested"]
GAP_COLS = ["species", "organism", "seq", "cell", "role", "annotation_loss",
            "bins_annotated", "bins_unannotated", "bins_annotated_with_reads",
            "bins_unannotated_with_reads", "covered_annotated",
            "covered_unannotated", "frac_coverage_unannotated"]
COV_COLS = ["species", "seq", "n_bins", "bin_width", "coverage"]


if __name__ == "__main__":
    s = build()
    print(f"\n{s['runs_quantified']}/{s['runs_selected']} runs -> "
          f"{OUT_DIR/'expression_stats.json'}")
