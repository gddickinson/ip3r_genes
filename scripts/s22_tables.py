"""S22's stats file and `headline()` — the numbers the report may not compute.

D13 in this task's own terms: the report renders from the committed tables,
and the handful of numbers a reader meets in the first paragraph are read
out here once so the report cannot arrive at a different value for them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402
import s22_shells as SH  # noqa: E402
import s22_contacts as C  # noqa: E402
import s22_paired as P  # noqa: E402
import s22_plc as PLC  # noqa: E402
import s22_deep_lineage as DL  # noqa: E402

TABLES = [
    "module_map.tsv", "ligand_shells.tsv", "shell_agreement.tsv",
    "shell_index.tsv", "tip_divergence.tsv", "module_contrast.tsv",
    "module_contrast_dropped.tsv", "column_contrast.tsv",
    "paralog_contrast.tsv", "contact_test.tsv", "shell_constraint.tsv",
    "shell_trend.tsv", "omega_by_module.tsv", "omega_module_test.tsv",
    "omega_by_shell.tsv", "omega_contacts.tsv", "plc_profiles.tsv",
    "plc_repertoire.tsv", "pathway_by_proteome.tsv",
    "pathway_cooccurrence.tsv", "plc_absent_taxa.tsv", "lineage_panel.tsv",
    "lineage_test.tsv", "lineage_power.tsv", "deep_lineage_panel.tsv",
    "deep_lineage_test.tsv", "deep_lineage_covariates.tsv",
    "deep_lineage_matched.tsv", "deep_lineage_matched_summary.tsv",
    "deep_lineage_matched_power.tsv", "deep_lineage_ryr_control.tsv",
    "deep_lineage_power.tsv",
]


def _num(rows, **where):
    for r in rows:
        if all(str(r.get(k, "")) == str(v) for k, v in where.items()):
            return r
    return {}


def headline() -> dict:
    """The numbers §1 and the abstract are allowed to state."""
    out: dict = {}
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")
    out["core_span"] = {r["paralog"]: f"{r['start']}-{r['end']}"
                        for r in mods if r["module"] == "ligand_core"
                        and r["is_primary"] == "True"}
    out["pore_span"] = {r["paralog"]: f"{r['start']}-{r['end']}"
                        for r in mods if r["module"] == "pore_module"
                        and r["is_primary"] == "True"}

    agree = L.read_tsv(L.OUT_DIR / "shell_agreement.tsv")
    out["n_structures"] = len(agree)
    out["s0_contacts_recovered"] = all(
        r["n_s0_contacts_recovered"] == r["n_s0_contacts"] for r in agree)
    shells = L.read_tsv(L.OUT_DIR / "ligand_shells.tsv")
    out["n_pocket_residues"] = len(shells)
    out["n_consensus_contacts"] = sum(
        1 for r in shells
        if int(r["n_structures_contact"]) >= C.CONSENSUS_MAJORITY)
    out["extra_contacts_vs_s0"] = sorted(
        int(r["resi"]) for r in shells
        if int(r["n_structures_contact"]) >= C.CONSENSUS_MAJORITY
        and r["is_s0_contact"] != "True")

    con = L.read_tsv(L.OUT_DIR / "module_contrast.tsv")
    out["module_primary"] = {
        r["paralog"]: {"n_tips": r["n_tips"],
                       "core": r["mean_core_identity"],
                       "pore": r["mean_pore_identity"],
                       "difference": r["mean_difference"],
                       "q": r["q_wilcoxon"], "direction": r["direction"]}
        for r in con if r["is_primary"] == "True"}
    out["module_loop_in"] = {
        r["paralog"]: {"difference": r["mean_difference"],
                       "q": r["q_wilcoxon"], "direction": r["direction"]}
        for r in con if r["core_definition"] == "contact_span"
        and r["pore_definition"] == "channel_all"}

    ct = L.read_tsv(L.OUT_DIR / "contact_test.tsv")
    out["contacts"] = {
        (r["paralog"], r["background"]): {"contact": r["mean_contact_jsd"],
                                          "background": r["mean_background_jsd"],
                                          "q": r["q_permutation"]}
        for r in ct if r["contact_set"] == "s0_contact" and r["layer"] == "deep"}
    tr = L.read_tsv(L.OUT_DIR / "shell_trend.tsv")
    out["trend"] = {r["paralog"]: (r["rho_distance_vs_jsd"],
                                   r["p_distance_vs_jsd"]) for r in tr}

    co = L.read_tsv(L.OUT_DIR / "pathway_cooccurrence.tsv")
    out["cooccurrence"] = {r["group"]: {k: r[k] for k in
                                        ("n_proteomes", "n_itpr", "n_plc",
                                         "n_both", "n_itpr_no_plc")}
                           for r in co}
    out["n_itpr_no_plc"] = sum(int(r["n_itpr_no_plc"]) for r in co)
    lt = L.read_tsv(L.OUT_DIR / "lineage_test.tsv")
    out["lineage"] = [r for r in lt if r["p_mannwhitney"]]
    pw = L.read_tsv(L.OUT_DIR / "lineage_power.tsv")
    out["power"] = pw

    deep = L.OUT_DIR / "deep_lineage_matched_summary.tsv"
    if deep.exists():
        out["matched"] = L.read_tsv(deep)[0]
        out["matched_power"] = L.read_tsv(
            L.OUT_DIR / "deep_lineage_matched_power.tsv")
        out["deep_strata"] = L.read_tsv(L.OUT_DIR / "deep_lineage_test.tsv")
        out["deep_covariates"] = L.read_tsv(
            L.OUT_DIR / "deep_lineage_covariates.tsv")
        rc = L.read_tsv(L.OUT_DIR / "deep_lineage_ryr_control.tsv")
        out["ryr_control"] = next((r for r in rc if r["tip"] == "SUMMARY"), {})
        panel = L.read_tsv(L.OUT_DIR / "deep_lineage_panel.tsv")
        out["deep_panel_n"] = len(panel)
        out["deep_panel_in_test"] = sum(1 for r in panel
                                        if r["in_test"] == "True")
    return out


def stats() -> dict:
    present = [t for t in TABLES if (L.OUT_DIR / t).exists()]
    return {
        "task": "S22",
        "rules": {
            "ligand_core_primary": "minimal contiguous span containing every "
                                   "IP3 contact measured on the structure",
            "pore_module_primary": "PF00520 channel span less the luminal "
                                   "loop located by geometry in S17",
            "sensitivity_definitions": ["ibc_literature (ITPR1 224-604, R05)",
                                        "channel_all (PF00520 as drawn)"],
            "pairing": "one core and one pore number per orthologue; a tip "
                       f"below {P.MIN_MODULE_COVERAGE} coverage of either "
                       "module is dropped and counted",
            "contact_backgrounds": ["rest_of_core", "rest_of_pocket"],
            "plc_definition": "a protein carrying both PF00387 and PF00388",
            "lineage_control": "the ryanodine receptors — same pore, no IP3",
            "divergence_matching": f"each PLC-absent record matched to up to "
                                   f"{DL.MAX_MATCHES} PLC-present records "
                                   f"within {DL.MATCH_TOLERANCE} pore "
                                   "identity of its own",
        },
        "parameters": {
            "search_radius_A": SH.SEARCH_RADIUS_A,
            "shell_edges_A": [[lo, hi, n] for lo, hi, n in SH.SHELL_EDGES],
            "structures": [SH.S0_STRUCTURE] + list(SH.REPLICATES),
            "consensus_majority": C.CONSENSUS_MAJORITY,
            "permutation_iters": C.PERM_ITERS,
            "primary_layer": L.PRIMARY_LAYER,
            "plc_pfam": PLC.PLC_PFAM,
            "evalue_primary": PLC.EVALUE_PRIMARY,
            "evalue_sweep": PLC.EVALUE_SWEEP,
            "groups_swept_for_plc": list(PLC.EUK_GROUPS),
        },
        "modules": {r["paralog"] + "/" + r["definition"]:
                    {"start": r["start"], "end": r["end"],
                     "n_residues": r["n_residues"], "check": r["check"]}
                    for r in L.read_tsv(L.OUT_DIR / "module_map.tsv")},
        "tables": {t: {"sha256": L.sha256(L.OUT_DIR / t),
                       "rows": max(0, sum(1 for _ in open(L.OUT_DIR / t)) - 1)}
                   for t in present},
        "tables_missing": [t for t in TABLES if t not in present],
    }


def main() -> int:
    self_test = __import__("s22_test_ligand").main()
    s = stats()
    s["self_test"] = "passed" if self_test == 0 else "FAILED"
    (L.OUT_DIR / "ligand_site_stats.json").write_text(json.dumps(s, indent=2))
    L.log(f"wrote ligand_site_stats.json — {len(s['tables'])} tables, "
          f"{len(s['tables_missing'])} missing, self-test {s['self_test']}")
    return 0 if self_test == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
