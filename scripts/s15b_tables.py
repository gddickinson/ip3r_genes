"""s15b_tables.py — every committed S15b table, written once.

The report and the figures render from these files and never recompute
(D13).  `headline()` is the set of numbers the report is not allowed to
compute for itself, and `write_stats()` records the parameters, the
self-test status and the SHA-256 of every table beside them, so a rebuild
that drifts is visible in the data (S6/S8/S9/S10's discipline).
"""

from __future__ import annotations

import json
import time

import s15b_coding as coding
import s15b_lib as lib

TABLES = {
    "dollo_counts": [
        "coding", "character", "evidence", "use_r5", "use_r6", "is_base",
        "n_tips_scored", "n_present", "n_absent", "n_undecided",
        "n_unknown_tips", "dollo_losses_max", "dollo_losses_min",
        "gain_node", "gain_rule", "mrca_of_present", "resolution_note"],
    "loss_edges": [
        "coding", "character", "evidence", "use_r5", "use_r6", "node",
        "parent", "parent_degree", "parent_is_polytomy",
        "siblings_lost_under_parent", "n_tips_lost", "n_tips_unknown",
        "tips"],
    "sensitivity_matrix": [
        "coding", "character", "evidence", "use_r5", "use_r6",
        "branch_lengths", "is_base", "n_tips_scored", "n_present",
        "n_absent", "n_undecided", "dollo_losses_max", "dollo_losses_min",
        "gain_node", "gain_rule", "mk_fitted", "mk_best_model", "mk_loss_rate",
        "mk_irrev_rate", "mk_irrev_loglik", "mk_refusal"],
    "manufactured_losses": [
        "accession", "organism", "vclass", "cell", "ledger_status",
        "base_state", "base_rule", "recon_coverage", "recon_n_contigs",
        "n_spare_itpr_loci", "contig_spans_gene", "n_settings_absent",
        "n_settings", "loosest_evidence", "needs_r5_off", "needs_r6_off",
        "released_by", "base_reason"],
    "mk_fits": [
        "coding", "character", "evidence", "use_r5", "use_r6",
        "branch_lengths", "n_absent", "model", "fitted", "n_params", "rate",
        "rate_gain", "loglik", "aic", "profile_shape", "profile_argmax",
        "reason", "first_seen_note"],
    "mk_profile": [
        "branch_lengths", "model", "axis", "rate", "loglik", "shape",
        "argmax", "at_boundary"],
    "resolution_bound": [
        "metric", "value", "note"],
    "fossil_loci": [
        "accession", "organism", "vclass", "cell", "contig", "coverage",
        "identity", "aligned_aa", "frameshifts", "stop_codons",
        "lesion_density", "verdict", "cell_state", "elevated", "incomplete",
        "not_live", "state_reading_available", "is_fossil", "reading"],
    "fossil_by_cell": [
        "cell", "n_scored", "n_above_bar", "n_fossils", "frac_above_bar",
        "median_identity", "median_identity_above_bar"],
    "lesion_by_class": [
        "cell", "vclass", "matched", "n_genomes", "n", "n_pos", "n_neg",
        "n_ties", "median_diff", "direction", "p", "q_bh", "underpowered"],
    "lesion_class_controls": [
        "cell", "vclass", "direction", "n", "q_bh", "n_above_bar",
        "pos_above", "neg_above", "p_above", "n_below_bar", "pos_below",
        "neg_below", "p_below", "median_identity", "median_others_identity",
        "sibling_of", "note"],
}


def headline(ctx: dict) -> dict:
    """The numbers the report reads rather than derives."""
    dc = ctx["dollo_counts"]
    base = [r for r in dc if r["is_base"]]
    fam_base = next(r for r in base if r["coding"] == "family")
    grid = ctx["sensitivity_matrix"]
    fam = [r for r in grid if r["coding"] == "family"]
    par = [r for r in grid if r["coding"] == "paralog"]
    settings = {(r["evidence"], r["use_r5"], r["use_r6"]) for r in grid}
    fam_bad = {(r["evidence"], r["use_r5"], r["use_r6"]) for r in fam
               if r["dollo_losses_max"] > 0}
    par_bad = {(r["evidence"], r["use_r5"], r["use_r6"]) for r in par
               if r["dollo_losses_max"] > 0}
    bl_spread = {}
    for r in grid:
        key = (r["coding"], r["character"], r["evidence"], r["use_r5"],
               r["use_r6"])
        bl_spread.setdefault(key, set()).add(r["dollo_losses_max"])
    fos = ctx["fossil_stats"]
    cls = [t for t in ctx["lesion_by_class"] if t["p"] != ""]
    sig = [t for t in cls if float(t.get("q_bh") or 1.0) <= 0.05]
    return dict(
        n_settings=len(settings),
        dollo_family_base=fam_base["dollo_losses_max"],
        dollo_paralog_base=max(r["dollo_losses_max"] for r in base
                               if r["coding"] == "paralog"),
        n_settings_family_manufactures=len(fam_bad),
        n_settings_paralog_manufactures=len(par_bad),
        max_family_losses=max(r["dollo_losses_max"] for r in fam),
        max_paralog_losses=max(r["dollo_losses_max"] for r in par),
        max_paralog_losses_min=max(r["dollo_losses_min"] for r in par),
        n_cells_ever_absent=len(ctx["manufactured_losses"]),
        n_cells=len({(r["accession"], r["cell"])
                     for r in ctx["character_matrix"]}),
        branch_length_changes_dollo=int(any(len(v) > 1
                                            for v in bl_spread.values())),
        gain_node_itpr1=next(r["gain_node"] for r in base
                             if r["character"] == "ITPR1"),
        gain_node_itpr2=next(r["gain_node"] for r in base
                             if r["character"] == "ITPR2"),
        gain_node_itpr3=next(r["gain_node"] for r in base
                             if r["character"] == "ITPR3"),
        gain_node_family=fam_base["gain_node"],
        # read off the grid, not off mk_fits: the fits are memoised on the
        # character vector, so mk_fits records the setting each *distinct*
        # character was first fitted under and cannot be filtered by setting
        n_mk_fitted_base=sum(1 for r in grid
                             if r["evidence"] == coding.BASE_SETTING[
                                 "evidence"]
                             and r["use_r5"] and r["use_r6"]
                             and r["mk_fitted"]),
        n_mk_rows_base=sum(1 for r in grid
                           if r["evidence"] == coding.BASE_SETTING["evidence"]
                           and r["use_r5"] and r["use_r6"]),
        n_mk_fits_total=sum(1 for r in ctx["mk_fits"] if r["fitted"]),
        n_mk_refusals_total=sum(1 for r in ctx["mk_fits"]
                                if not r["fitted"]),
        n_scored_loci=fos["n_scored"], n_above_lesion_bar=fos["n_above_bar"],
        n_fossils_itpr=fos["n_fossils_itpr"],
        n_above_bar_full_coverage=fos["n_above_bar_full_coverage"],
        n_class_tests=len(cls), n_class_significant=len(sig),
        strongest_class=(min(sig, key=lambda t: float(t["q_bh"]))["vclass"]
                         if sig else ""),
        strongest_class_cell=(min(sig, key=lambda t: float(t["q_bh"]))["cell"]
                              if sig else ""),
        strongest_class_q=(min(float(t["q_bh"]) for t in sig)
                           if sig else ""),
    )


def write_all(ctx: dict) -> dict[str, str]:
    lib.OUT.mkdir(parents=True, exist_ok=True)
    written = {}
    for name, cols in TABLES.items():
        rows = ctx.get(name)
        if rows is None:
            continue
        path = lib.OUT / f"{name}.tsv"
        lib.write_tsv(path, rows, cols)
        written[f"{name}.tsv"] = lib.sha256(path)
    return written


def write_stats(ctx: dict, written: dict[str, str], self_test: dict) -> None:
    payload = dict(
        task="S15b",
        written_at=time.strftime("%Y-%m-%d %H:%M:%S"),
        inputs=dict(
            character_matrix=str(lib.MATRIX.relative_to(lib.PROJECT)),
            species_tree=str(lib.NEWICK.relative_to(lib.PROJECT)),
            integrity_loci=str(lib.INTEGRITY_LOCI.relative_to(lib.PROJECT)),
            integrity_pairs=str(lib.INTEGRITY_PAIRS.relative_to(lib.PROJECT)),
            calibrations=str(lib.CALIBRATIONS.relative_to(lib.PROJECT))),
        parameters=dict(
            evidence_levels=list(coding.EVIDENCE_NAMES),
            recon_bar_calibrated=coding.BAR_CALIBRATED,
            recon_bar_gap=[coding.BAR_GAP_LO, coding.BAR_GAP_HI],
            base_setting=coding.BASE_SETTING,
            branch_length_schemes=list(lib.BL_SCHEMES),
            mk_models=list(ctx.get("mk_models", ())),
            lesion_bar=ctx["fossil_stats"]["lesion_bar"],
            class_min_n=ctx.get("class_min_n"),
        ),
        base_reproduction=ctx.get("base_reproduction"),
        polytomy_profile=ctx.get("polytomy_profile"),
        fossil=ctx["fossil_stats"],
        headline=ctx["headline"],
        self_test=self_test,
        tables=written)
    (lib.OUT / "loss_counts_stats.json").write_text(
        json.dumps(payload, indent=1, default=str), encoding="utf-8")


def resolution_bound_rows(profile: dict, dollo_rows: list[dict]) -> list[dict]:
    """The bound any count on this tree is read against, as a table."""
    with_losses = [r for r in dollo_rows if r["dollo_losses_max"] > 0]
    gap = [r for r in with_losses
           if r["dollo_losses_max"] != r["dollo_losses_min"]]
    return [
        dict(metric="internal_nodes", value=profile["n_internal"],
             note="including the unary species nodes the sweep's "
                  "one-tip-per-assembly design creates"),
        dict(metric="unary_nodes", value=profile["n_unary"],
             note="one per species represented by a single assembly"),
        dict(metric="polytomies", value=profile["n_polytomies"],
             note="internal nodes of degree 3 or more"),
        dict(metric="max_polytomy_degree", value=profile["max_degree"],
             note="the largest unresolved node in the tree"),
        dict(metric="edges_below_polytomies",
             value=profile["n_edges_below_polytomies"],
             note="edges whose sibling relationships the tree does not "
                  "resolve; a loss on one of these is bounded, not counted"),
        dict(metric="settings_with_any_loss", value=len(with_losses),
             note="rows of dollo_counts.tsv placing at least one loss"),
        dict(metric="settings_where_bound_is_wide", value=len(gap),
             note="rows where the maximum and minimum independent counts "
                  "differ, i.e. where a polytomy carries more than one loss"),
    ]
