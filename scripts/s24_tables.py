"""S24's committed tables and `supplementary_stats.json`.

`supp_figure_stats.tsv` is the brief's "their stats table": one row per
supplementary figure with the file it was drawn from, the numbers it draws,
and the SHA-256 of every png and pdf, so a rebuild that drifts is visible in
the data (S6/S8/S9/S15's discipline).

`headline()` is the set of numbers the report is not allowed to compute for
itself (D13).
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s24_lib as L                                            # noqa: E402

FIGURE_SOURCES = {
    1: ["msa_v2/aln.fasta", "msa_v2/trimmed.fasta",
        "msa_v2/column_map.tsv", "msa_v2/representatives.tsv"],
    2: ["constraint/variants.tsv",
        "constraint/paralog_variant_positions.tsv",
        "constraint/constraint_ITPR1_Q14643.tsv",
        "constraint/constraint_ITPR2_Q14571.tsv",
        "constraint/constraint_ITPR3_Q14573.tsv"],
    3: ["constraint/ortholog_shape.tsv", "constraint/layer_sizes.tsv",
        "constraint/orthologs_manifest.tsv"],
    4: ["selection/codon_aln.fasta", "selection/codon_trimmed.fasta",
        "selection/tip_codes.tsv", "selection/cds_status.tsv"],
    5: ["constraint/painted/constraint_reference_ITPR1_7LHF.pdb",
        "constraint/painted/constraint_reference_ITPR2_9YKK.pdb",
        "constraint/painted/constraint_reference_ITPR3_8TKG.pdb",
        "constraint/painted/selection_reference_ITPR3_8TKG.pdb"],
    6: ["constraint/painted/constraint_reference_ITPR2_9YKK.pdb",
        "constraint/painted/constraint_reference_ITPR3_8TKG.pdb",
        "constraint/variant_by_element.tsv", "constraint/variants.tsv"],
}

FIGURE_STEM = {
    1: "SuppFig1_representative_alignment",
    2: "SuppFig2_labelled_positions",
    3: "SuppFig3_paralog_alignments",
    4: "SuppFig4_codon_alignment",
    5: "SuppFig5_constraint_on_channel",
    6: "SuppFig6_variants_on_structure",
}

FIGURE_CAPTION = {
    1: "The representative alignment, and the columns the tree actually saw",
    2: "The ligand core and the pore module at residue resolution across the "
       "three paralogues",
    3: "The within-paralogue alignments the constraint map is computed on",
    4: "The trimmed codon alignment behind every omega estimate",
    5: "The constraint map painted onto the channel",
    6: "Every labelled variant, and the per-element enrichment test",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write(stats: dict) -> dict:
    rows = []
    for num in sorted(FIGURE_STEM):
        stem = FIGURE_STEM[num]
        key = f"supp_fig_{num}"
        drawn = stats.get(key, {})
        for ext in ("png", "pdf"):
            path = L.FIG_DIR / f"{stem}.{ext}"
            rows.append({
                "number": num, "figure": stem, "format": ext,
                "caption": FIGURE_CAPTION[num],
                "sources": ";".join(FIGURE_SOURCES[num]),
                "bytes": path.stat().st_size if path.exists() else "",
                "sha256": sha256(path) if path.exists() else "",
                "drawn_from": json.dumps(drawn, sort_keys=True)
                if ext == "png" else "",
            })
    L.write_tsv(L.OUT_DIR / "supp_figure_stats.tsv", rows,
                ["number", "figure", "format", "caption", "sources", "bytes",
                 "sha256", "drawn_from"])
    return {"rows": len(rows),
            "figures": len(FIGURE_STEM),
            "missing": [r["figure"] for r in rows if not r["sha256"]]}


def headline(stats: dict) -> dict:
    """The numbers the report may quote, computed once, here."""
    g = stats.get("guards", {})
    f1 = stats.get("supp_fig_1", {})
    f2 = stats.get("supp_fig_2", {})
    f3 = stats.get("supp_fig_3", {})
    f4 = stats.get("supp_fig_4", {})
    f6 = stats.get("supp_fig_6", {})
    audit = stats.get("audit", {})
    return {
        "columns_checked": g.get("columns_checked"),
        "sequences_checked": g.get("sequences"),
        "variant_residues_checked": g.get("variant_residues_checked"),
        "aligned_partners_checked": g.get("aligned_partners_checked"),
        "aln_columns": f1.get("input_columns"),
        "trimmed_columns": f1.get("kept_columns"),
        "median_occupancy_kept": f1.get("median_occupancy_kept"),
        "median_occupancy_cut": f1.get("median_occupancy_cut"),
        "pathogenic_positions_drawn": f2.get("positions_drawn"),
        "pathogenic_conserved_all_three":
            f2.get("conserved_across_all_three"),
        "orthologs_screened": f3.get("sequences_screened"),
        "orthologs_dropped": f3.get("sequences_dropped"),
        "codons_in": f4.get("codons_in"),
        "codons_kept": f4.get("codons_kept"),
        "structures_eligible_for_variants": sorted(
            f6.get("structures_drawn", {})),
        "structures_refused": [c["file"] for c in
                               f6.get("numbering_checks", [])
                               if not c["carries_human_numbering"]],
        # Flat scalars, because the manuscript's claims ledger addresses a
        # JSON value by a dotted key path and cannot index a list.
        "itpr1_afdb_model_residues": next(
            (c["residues"] for c in f6.get("numbering_checks", [])
             if c["file"] == "constraint_model_Q14643.pdb"), None),
        "itpr1_reference_agree": next(
            (c["agree"] for c in f6.get("numbering_checks", [])
             if c["file"] == "constraint_reference_ITPR1_7LHF.pdb"), None),
        "element_tests": f6.get("element_tests"),
        "element_tests_significant": f6.get("significant_after_bh"),
        "figures_audited": audit.get("figures_audited"),
        "audit_findings": audit.get("findings"),
        "audit_findings_by_status": audit.get("findings_by_status"),
    }
