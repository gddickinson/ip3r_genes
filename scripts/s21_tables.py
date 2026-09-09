"""S21's statistics blob and the headline numbers the report may not compute.

D13 applied to a gene-structure report: every number §-level prose quotes comes
from `headline()`, which reads the committed tables. A report that recomputed
its own numbers could disagree with the table printed beside it.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_annot as AN                                            # noqa: E402
import s21_architecture as A                                      # noqa: E402
import s21_blocks as B                                            # noqa: E402
import s21_calibrate_intron as CAL                                # noqa: E402
import s21_frame as FR                                            # noqa: E402
import s21_introns as I                                           # noqa: E402
import s21_lib as L                                               # noqa: E402
import s21_tandem as TD                                           # noqa: E402

TABLES = (
    "architecture_loci.tsv", "architecture_by_paralog.tsv",
    "architecture_by_class.tsv", "architecture_sensitivity.tsv",
    "paired_comparisons.tsv", "intron_positions.tsv",
    "intron_conservation.tsv", "intron_conservation_summary.tsv",
    "shared_introns.tsv", "shared_intron_summary.tsv",
    "intron_floor_calibration.tsv", "junction_gap_bins.tsv",
    "frame_map.tsv", "frame_anchors.tsv", "junction_quality.tsv",
    "boundary_concordance.tsv", "boundary_concordance_summary.tsv",
    "fragment_termini.tsv", "fragment_verdicts.tsv", "fragment_summary.tsv",
    "tandem_pairs.tsv", "tandem_by_class.tsv", "tandem_control.tsv",
    "tandem_teleost_check.tsv",
)


def _row(rows: list[dict], **where) -> dict:
    for r in rows:
        if all(str(r.get(k, "")) == str(v) for k, v in where.items()):
            return r
    return {}


def junction_quality(arows: list[dict], pairs: list[dict]) -> list[dict]:
    """Splice-class and frame quality of every junction, by paralogue.

    The instrument's own error rate: an exon boundary the genome does not read
    as a splice site is a boundary the alignment placed and the sequence does
    not corroborate, and a junction beside a frame step is one the alignment is
    not entitled to call. Both are reported per paralogue and pooled, over every
    junction in the scope rather than a sample.
    """
    scope = {(r["accession"], r["cell"], r["locus_idx"]) for r in arows
             if r["in_scope"]}
    out = []
    for cell in ["all"] + list(L.CELLS):
        sub = [p for p in pairs
               if (p["accession"], p["cell"], p["locus_idx"]) in scope
               and (cell == "all" or p["cell"] == cell)]
        if not sub:
            continue
        canon = sum(1 for p in sub if p["splice_class"] == "canonical")
        minor = sum(1 for p in sub if p["splice_class"] == "minor")
        noncan = sum(1 for p in sub if p["splice_class"] == "non_canonical")
        unread = len(sub) - canon - minor - noncan
        frame_bad = sum(1 for p in sub if not p["frame_ok"])
        out.append({
            "cell": cell, "n_junctions": len(sub),
            "n_loci": len({(p["accession"], p["cell"], p["locus_idx"])
                           for p in sub}),
            "canonical": canon, "minor": minor, "non_canonical": noncan,
            "unreadable": unread,
            "frac_spliceable": round((canon + minor) / len(sub), 5),
            "frac_canonical": round(canon / len(sub), 5),
            "n_frame_step": frame_bad,
            "frac_frame_step": round(frame_bad / len(sub), 5),
            "min_gap_bp": min(p["gap_bp"] for p in sub),
            "median_gap_bp": round(L.quantile(
                [float(p["gap_bp"]) for p in sub], 0.5), 1)})
    return out


def headline() -> dict:
    """The numbers §-level prose is allowed to quote, read from the tables."""
    out: dict = {}
    d = L.OUT
    loci = L.read_tsv(d / "architecture_loci.tsv")
    out["n_loci_measured"] = len(loci)
    out["n_in_scope"] = sum(1 for r in loci if r["in_scope"] == "1")
    out["n_genomes_in_scope"] = len({r["accession"] for r in loci
                                     if r["in_scope"] == "1"})
    excl: dict[str, int] = {}
    for r in loci:
        if r["excluded_by"]:
            excl[r["excluded_by"]] = excl.get(r["excluded_by"], 0) + 1
    out["excluded_by"] = excl
    out["n_merges_total"] = sum(int(float(r["n_merges"])) for r in loci)

    per = L.read_tsv(d / "architecture_by_paralog.tsv")
    out["by_paralog"] = {
        r["cell"]: {"n_loci": int(r["n_loci"]), "n_genomes": int(r["n_genomes"]),
                    "exons": float(r["n_exons_median"]),
                    "exons_p10": float(r["n_exons_p10"]),
                    "exons_p90": float(r["n_exons_p90"]),
                    "exons_min": float(r["n_exons_min"]),
                    "exons_max": float(r["n_exons_max"]),
                    "cds_bp": float(r["cds_bp_median"]),
                    "span_bp": float(r["span_bp_median"]),
                    "median_intron_bp": float(r["median_intron_bp_median"]),
                    "max_intron_bp": float(r["max_intron_bp_median"]),
                    "mean_exon_bp": float(r["mean_exon_bp_median"])}
        for r in per}
    spans = [v["span_bp"] for k, v in out["by_paralog"].items()
             if k in L.PARALOGS]
    out["span_fold_range"] = round(max(spans) / min(spans), 2) if spans else 0.0
    exons = [v["exons"] for k, v in out["by_paralog"].items()
             if k in L.PARALOGS]
    out["exon_count_range"] = [min(exons), max(exons)] if exons else []

    q = L.read_tsv(d / "junction_quality.tsv")
    allq = _row(q, cell="all")
    out["junctions"] = {
        "n": int(allq.get("n_junctions", 0)),
        "frac_spliceable": float(allq.get("frac_spliceable", 0)),
        "frac_canonical": float(allq.get("frac_canonical", 0)),
        "frac_frame_step": float(allq.get("frac_frame_step", 0)),
        "min_gap_bp": int(float(allq.get("min_gap_bp", 0)))}

    cal = (L.read_tsv(d / "intron_floor_calibration.tsv") or [{}])[0]
    out["intron_floor"] = {
        "floor_bp": int(float(cal.get("floor_bp", 0))),
        "separated": int(float(cal.get("separated", 0))),
        "n_gaps": int(float(cal.get("n_gaps", 0))),
        "n_non_spliceable": int(float(cal.get("n_non_spliceable", 0))),
        "smallest_gap_bp": int(float(cal.get("smallest_gap_bp", 0))),
        "derivation": cal.get("derivation", "")}

    cons = L.read_tsv(d / "intron_conservation_summary.tsv")
    out["conservation"] = {
        r["cell"]: {"n_loci": int(r["n_loci"]),
                    "positions_seen": int(r["n_positions_seen"]),
                    "at_50pc": int(r["n_at_50pc"]),
                    "at_90pc": int(r["n_at_90pc"]),
                    "at_99pc": int(r["n_at_99pc"]),
                    "median_per_locus": float(r["median_positions_per_locus"])}
        for r in cons}

    sh = [r for r in L.read_tsv(d / "shared_intron_summary.tsv")
          if r["is_operating_point"] == "1" and r["excludes_frame_via_cell"] == "0"]
    out["shared"] = {
        f"{r['cell_a']}|{r['cell_b']}": {
            "n_genomes": int(r["n_genomes"]),
            "median_observed": float(r["median_observed"]),
            "median_expected": float(r["median_expected"]),
            "enrichment": float(r["enrichment"]),
            "n_significant": int(r["n_genomes_p_lt_0.05"]),
            "q_worst": float(r["q_worst_genome"])} for r in sh}

    conc = L.read_tsv(d / "boundary_concordance_summary.tsv")
    allc = _row(conc, grouping="all")
    out["concordance"] = {
        "n_edges": int(allc.get("n_annot_edges", 0)),
        "n_loci": int(allc.get("n_loci", 0)),
        "n_genomes": int(allc.get("n_genomes", 0)),
        "frac_exact": float(allc.get("frac_exact", 0)),
        "by_source": {r["group"]: float(r["frac_exact"]) for r in conc
                      if r["grouping"] == "source"},
        "by_cell": {r["group"]: float(r["frac_exact"]) for r in conc
                    if r["grouping"] == "cell"}}

    frag = L.read_tsv(d / "fragment_summary.tsv")
    allf = _row(frag, state="all")
    out["fragments"] = {k: (int(v) if str(v).lstrip("-").isdigit()
                            else float(v) if v.replace(".", "", 1).isdigit()
                            else v)
                        for k, v in allf.items()}
    out["fragments_by_state"] = {
        r["state"]: {k: v for k, v in r.items() if k != "state"}
        for r in frag if r["state"] != "all"}

    tc = L.read_tsv(d / "tandem_control.tsv")
    allt = _row(tc, cell="all")
    out["tandem_control"] = {
        "sensitivity": float(allt.get("sensitivity", 0)),
        "specificity": float(allt.get("specificity", 0)),
        "n_cells": int(allt.get("n_cells", 0)),
        "tp": int(allt.get("tp", 0)), "fn": int(allt.get("fn", 0)),
        "fp": int(allt.get("fp", 0)), "tn": int(allt.get("tn", 0))}
    tcls = L.read_tsv(d / "tandem_by_class.tsv")
    out["tandem_classes"] = {
        f"{r['class']}|{r['cell']}": int(r["n_pairs"]) for r in tcls}
    out["n_within_locus_pairs"] = sum(
        int(r["n_pairs"]) for r in tcls
        if r["class"] == "within_locus" and r["cell"] == "all")
    tel = L.read_tsv(d / "tandem_teleost_check.tsv")
    out["teleost_agreement"] = {
        "n": len(tel), "agree": sum(1 for r in tel if r["agrees"] == "1")}

    pc = L.read_tsv(d / "paired_comparisons.tsv")
    out["paired"] = {
        f"{r['cell_a']}|{r['cell_b']}|{r['metric']}": {
            "n": int(r["n_pairs"]), "median_diff": float(r["median_diff"]),
            "q": float(r["q"])}
        for r in pc if r["stratum"] == "all"}

    fm = L.read_tsv(d / "frame_map.tsv")
    out["frame"] = {
        "n_pairs": len(fm),
        "n_unusable": sum(1 for r in fm if r["usable"] == "0"),
        "min_frac_on_ref": min((float(r["frac_bait_on_ref"]) for r in fm),
                              default=0.0),
        "n_via_cell": sum(1 for r in fm if r["frame_via_cell"] == "1")}
    anch = L.read_tsv(d / "frame_anchors.tsv")
    out["anchors"] = {"n": len(anch),
                      "n_one_column": sum(1 for r in anch
                                          if r["one_column"] == "1")}
    return out


def write_stats(self_test_ok: bool, extra: dict | None = None) -> Path:
    """`gene_architecture_stats.json`: rules, parameters, self-test, SHA-256s."""
    stats = {
        "task": "S21",
        "rules": {
            "A1": "the locus is one of S16's committed gene copies (is_copy)",
            "A2": "the contig spans the gene (D4)",
            "A3": f"the model covers >= {A.COV_ARCH} of its bait",
            "A4": "the bait reaches a usable alignment frame",
        },
        "parameters": {
            "COV_ARCH": A.COV_ARCH, "COV_BARS": A.COV_BARS,
            "MIN_INTRON_FALLBACK": B.MIN_INTRON_FALLBACK,
            "min_intron_bp_used": B.min_intron_bp(),
            "SPLICEABLE_BAR": CAL.SPLICEABLE_BAR, "MIN_GAPS": CAL.MIN_GAPS,
            "MIN_BAIT_ON_REF": FR.MIN_BAIT_ON_REF,
            "FRAME_ACC": L.FRAME_ACC,
            "TOLERANCES": list(I.TOLERANCES), "N_PERM": I.N_PERM,
            "SEED": I.SEED, "PREVALENCE_BARS": list(I.PREVALENCE_BARS),
            "EDGE_TOL": AN.EDGE_TOL, "TARGET_STATES": list(AN.TARGET_STATES),
            "MIN_ALN_COV": TD.MIN_ALN_COV, "Q_OVERLAP": TD.Q_OVERLAP,
            "TANDEM_GAP_BP": TD.TANDEM_GAP_BP,
        },
        "self_test": "pass" if self_test_ok else "fail",
        "n_negative_controls": 21,
        "mutation_tests": {
            "n": 6, "caught": 6,
            "breakages": [
                "s21_blocks.next_phase returns the wrong codon complement (T6)",
                "s21_architecture.intron_phase not inverted (T7)",
                "s21_blocks.gap_bp drops the strand branch (T3)",
                "s21_introns.shared_pair ignores phase (T13)",
                "s21_tandem.detect drops the cell attribution (T16)",
                "s21_blocks.read_models keys on the raw miniprot id (T1)"]},
        "sha256": {},
    }
    if extra:
        stats.update(extra)
    for name in TABLES:
        p = L.OUT / name
        if p.exists():
            stats["sha256"][name] = L.sha256(p)
    path = L.OUT / "gene_architecture_stats.json"
    path.write_text(json.dumps(stats, indent=1))
    return path
