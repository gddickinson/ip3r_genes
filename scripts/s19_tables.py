"""S19's statistics blob and the headline numbers the report may not compute.

D13 applied to a methods report: every number §-level prose quotes comes from
`headline()`, which reads the committed tables. A report that recomputed its
own numbers could disagree with the table beside it.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_drift as DR                                          # noqa: E402
import s19_floor as CT                                          # noqa: E402
import s19_lib as S                                             # noqa: E402
import s19_panel as PN                                          # noqa: E402

TABLES = (
    "census_growth.tsv", "method_sets.tsv", "method_pairwise.tsv",
    "head_to_head.tsv", "gene_recovery.tsv", "gene_recovery_by_cell.tsv",
    "gene_recovery_by_scope.tsv", "db_status_combined.tsv",
    "method_cost.tsv", "contiguity_cells.tsv", "contiguity_bins.tsv",
    "contiguity_by_stratum.tsv", "contiguity_tests.tsv",
    "absence_floor.tsv", "floor_clade_composition.tsv",
    "residual_neighbourhood.tsv", "panel_cells.tsv", "panel_recall.tsv",
    "panel_changes.tsv", "panel_by_group.tsv", "panel_unfilled_slots.tsv",
    "panel_single_bait.tsv", "panel_identity_recall.tsv",
    "panel_min_identity.tsv", "panel_validation.tsv", "jackhmmer_runs.tsv",
    "jackhmmer_rounds.tsv", "kill_criterion_trace.tsv",
    "kill_criterion_validation.tsv", "drift_outcome.tsv",
    "iteration_yield.tsv", "seed_effect.tsv",
    "reconciliation_loss_audit.tsv", "synteny_power.tsv",
    "codon_model_power.tsv", "saturation.tsv", "tree_resolution.tsv",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _row(rows: list[dict], **where) -> dict:
    for r in rows:
        if all(str(r.get(k, "")) == str(v) for k, v in where.items()):
            return r
    return {}


def headline() -> dict:
    """The numbers §-level prose is allowed to quote, read from the tables."""
    out: dict = {}
    d = S.out_dir()

    floor = S.read_tsv(d / "absence_floor.tsv")
    cells = S.read_tsv(d / "contiguity_cells.tsv")
    for series in ("itpr_present", "ryr_sister"):
        sub = [c for c in cells if c["control"] == series]
        fn = sum(1 for c in sub if c["false_negative"] == "1")
        lo, hi = S.wilson(fn, len(sub))
        at_bar = _row(floor, series=series, is_d4_bar="1")
        out[series] = {
            "n_cells": len(sub), "n_false_negative": fn,
            "rate": round(fn / max(1, len(sub)), 4),
            "wilson": [round(lo, 4), round(hi, 4)],
            "rate_at_d4_bar": S.fnum(at_bar.get("fn_rate"), float, 0.0),
            "genomes_at_d4_bar": int(S.fnum(at_bar.get("n_genomes_kept"),
                                            float, 0)),
        }
    out["d4_bar"] = S.contiguity_bar()
    out["n_genomes"] = len({c["accession"] for c in cells})
    out["not_counted_as_control"] = sum(1 for c in cells if not c["control"])

    tests = S.read_tsv(d / "contiguity_tests.tsv")
    out["contiguity_effect"] = {
        r["series"]: r["detail"] for r in tests
        if r["test"] == "logit_found~log10(contig_n50)"}
    out["itpr_vs_ryr_p"] = _row(
        tests, test="fisher_itpr_vs_ryr_false_negative").get("p", "")

    recall = S.read_tsv(d / "panel_recall.tsv")
    out["panel"] = {r["panel"]: {"n_baits": int(r["n_baits"]),
                                 "recall": float(r["itpr_recall"]),
                                 "found": int(r["n_itpr_found"]),
                                 "delta": int(r["delta_vs_full"])}
                    for r in recall}
    out["panel_validation_mismatches"] = len(
        S.read_tsv(d / "panel_validation.tsv"))
    out["panel_cells_validated"] = out["n_genomes"] * len(S.CELLS)
    ident = S.read_tsv(d / "panel_identity_recall.tsv")
    out["single_bait"] = [
        {"bin": r["identity_bin"], "n": int(r["n_measurements"]),
         "recall": float(r["recall"]), "lo": float(r["wilson_lo"])}
        for r in ident]

    h2h = S.read_tsv(d / "head_to_head.tsv")
    out["head_to_head_gene_scale"] = {
        r["database"]: {"pfam": int(r["n_pfam"]), "hmm": int(r["n_hmm"]),
                        "shared": int(r["n_shared"]),
                        "hmm_only": int(r["n_hmm_only"]),
                        "pfam_only": int(r["n_pfam_only"])}
        for r in h2h if r["length_band"].startswith(">=2000") and r["n_pfam"]}
    out["head_to_head_tail"] = {
        r["database"]: {"pfam": int(r["n_pfam"]), "hmm": int(r["n_hmm"]),
                        "hmm_only": int(r["n_hmm_only"]),
                        "pfam_only": int(r["n_pfam_only"])}
        for r in h2h if r["length_band"].startswith("<1000") and r["n_pfam"]}

    rec = S.read_tsv(d / "gene_recovery.tsv")
    genes = [r for r in rec if r["gene_present"] == "1"]
    by_channel: dict[str, int] = {}
    for r in genes:
        by_channel[r["recovery_channel"]] = by_channel.get(
            r["recovery_channel"], 0) + 1
    out["gene_recovery"] = {
        "n_genes": len(genes), "by_channel": by_channel,
        "invisible": sum(v for k, v in by_channel.items()
                         if k.startswith("genome_only")),
        "by_scope": {r["scope"]: float(r["frac"])
                     for r in S.read_tsv(d / "gene_recovery_by_scope.tsv")},
        "by_cell": {r["cell"]: float(r["frac"])
                    for r in S.read_tsv(d / "gene_recovery_by_cell.tsv")},
    }

    val = S.read_tsv(d / "kill_criterion_validation.tsv")
    out["kill_rules"] = {r["rule"]: {"sensitivity": float(r["sensitivity"]),
                                     "specificity": float(r["specificity"]),
                                     "tp": int(r["true_positive"]),
                                     "fp": int(r["false_positive"])}
                         for r in val}
    trace = S.read_tsv(d / "kill_criterion_trace.tsv")
    out["drift"] = {
        "n_runs": len(trace),
        "k1_never": sum(1 for r in trace if r["first_round_K1"] == "never"),
        "sister_fell": sum(1 for r in trace if r["sister_share_fell"] == "1"),
        "drifted": sum(1 for r in S.read_tsv(d / "drift_outcome.tsv")
                       if r["drifted"] == "1"),
    }

    out["inference"] = S.read_json(d / "inference_summary.json")
    growth = S.read_tsv(d / "census_growth.tsv")
    out["census"] = [{"v": r["census"], "n": int(r["n_records"]),
                      "added": int(r["n_added"]),
                      "channel": r["channel_added"]} for r in growth]
    cost = S.read_tsv(d / "method_cost.tsv")
    out["recorded_hours"] = round(sum(S.fnum(r["elapsed_h"], float, 0.0)
                                      for r in cost), 2)
    out["channels_without_timing"] = [
        r["method"] for r in cost if not r["elapsed_h"]]
    return out


def write_stats(log=S.log) -> Path:
    d = S.out_dir()
    payload = {
        "task": "S19 — methods results",
        "rules": {
            "control_series": (
                "a cell is a control iff S15a states its gene present "
                "(923 ITPR cells) or it is the RyR sister cell (309); "
                "`paralog_unassignable` is in neither"),
            "false_negative": ("a control cell the S5 ledger did not grade "
                               "`found_*`"),
            "d4_bar_bp": S.contiguity_bar(),
            "fn_targets": list(CT.FN_TARGETS),
            "floor_grid": CT.FLOOR_GRID,
            "neighbourhood_present": CT.NEIGHBOURHOOD_PRESENT,
            "panel_cov_found": PN.COV_FOUND,
            "drifted_offfamily_frac": DR.DRIFTED_OFFFAMILY_FRAC,
            "proposed_offfamily_rise": DR.MAX_OFFFAMILY_RISE,
        },
        "self_test": _self_test_status(),
        "headline": headline(),
        "sha256": {name: sha256(d / name) for name in TABLES
                   if (d / name).exists()},
        "missing_tables": [n for n in TABLES if not (d / n).exists()],
    }
    path = S.write_json(d / "methods_stats.json", payload)
    log(f"methods_stats.json: {len(payload['sha256'])} tables checksummed, "
        f"{len(payload['missing_tables'])} missing")
    return path


def _self_test_status() -> dict:
    import subprocess
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).with_name("s19_test_methods.py"))],
        capture_output=True, text=True)
    lines = [x for x in proc.stdout.splitlines() if x.strip().startswith(
        ("ok ", "FAIL"))]
    return {"exit_code": proc.returncode, "n_checks": len(lines),
            "n_failed": sum(1 for x in lines if x.strip().startswith("FAIL")),
            "passed": proc.returncode == 0}


def run(log=S.log) -> dict:
    write_stats(log)
    return headline()


if __name__ == "__main__":
    run()
