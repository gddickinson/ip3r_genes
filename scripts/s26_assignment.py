"""S26 — which paper each results entry is primary in, and why (the data).

One row per entry under `results/`, or per file group for a directory that
P5s splits. Each row gives the paper the entry is primary in, the rule that
placed it there, what the entry is, the paper it came closest to going to
instead, and the rule that kept it out of that paper. Paper `-` means the
entry is excluded, and the reason column says why it is not a result of any
paper. `s26_assign.py` enforces the table; this module only declares it.

The rules are in `papers/paper_rules.md`, committed before this table.
"""

from __future__ import annotations

#: (entry, file glob inside it or "*", paper, rule, what it is,
#:  nearest other paper, why not there)
ASSIGNMENT: list[tuple[str, str, str, str, str, str, str]] = [
    # ---------------------------------------------------------- range
    ("results/benchmark_controls", "*", "range", "P2",
     "S1: recall and specificity of the family call on labelled panels, "
     "with the ryanodine receptors as the sharp decoy",
     "origin", "P2: it is the negative control of the family call, and the "
     "family call is the range paper's instrument"),
    ("results/census_v2", "*", "range", "P1",
     "S2: the uncapped InterPro enumeration with a positive ITPR/RYR call "
     "on every record", "archive", "P1: it answers where the family is, "
     "not how well it is recorded"),
    ("results/census_v3", "*", "range", "P1",
     "S3: the two-instrument census and the jackhmmer completeness argument",
     "archive", "P1"),
    ("results/hmm_sweep", "*", "range", "P2",
     "S3: the ITPR and RyR profiles and the vertebrate proteome sweep",
     "retention", "P2: the profile margin is the range paper's positive "
     "family test"),
    ("results/s20_sweep", "*", "range", "P1",
     "S20: the non-vertebrate proteome sweep, 6,928 reference proteomes",
     "archive", "P1"),
    ("results/s23_baits", "*", "range", "P2",
     "S23: the non-vertebrate bait panel and its per-clade controls",
     "retention", "P2: the MIR-sharer control is what makes an absence a "
     "measurement"),
    ("results/s23_scope", "*", "range", "P1",
     "S23: the 194-genome non-vertebrate sweep and the controlled absences",
     "archive", "P1"),
    ("results/census_v5", "*", "range", "P3",
     "S20: census v4 plus the non-vertebrate sweep", "archive",
     "P3: the census is the range paper's denominator"),
    ("results/census_v6", "*", "range", "P3",
     "S23: census v5 plus the non-vertebrate gene models", "origin",
     "P3: the representatives are drawn from it but it is not a tree result"),
    ("results/structures", "*", "range", "P2",
     "S11: the fold of what the census calls an IP3 receptor, against the "
     "cryo-EM IP3R and RyR references, and AlphaFold DB coverage",
     "constraint", "P2: a structural family call is a third, independent "
     "control of the call the range paper makes"),
    ("results/methods", "jackhmmer_*", "range", "P5s",
     "S19: the jackhmmer runs re-derived from their raw logs", "retention",
     "P2: D10's kill criterion guards the range paper's completeness claim"),
    ("results/methods", "kill_criterion_*", "range", "P5s",
     "S19: D10's kill criterion scored as a classifier", "retention", "P2"),
    ("results/methods", "drift_*", "range", "P5s",
     "S19: the drift outcome measured on the finished models", "retention",
     "P2"),
    ("results/methods", "iteration_yield.tsv", "range", "P5s",
     "S19: what each jackhmmer round added", "retention", "P2"),
    ("results/methods", "seed_effect.tsv", "range", "P5s",
     "S19: the effect of the jackhmmer seed", "retention", "P2"),
    ("results/methods", "figures/fig_s19_drift.*", "range", "P5s",
     "S19: the drift figure", "retention", "P2"),
    # ---------------------------------------------------------- origin
    ("results/msa_v2", "*", "origin", "P1",
     "S6: the representative alignment and its identity matrices",
     "range", "P1: it is built to place the paralogues, not to find them"),
    ("results/phylogeny", "*", "origin", "P1",
     "S7: the ML tree and the AU test of the three sister hypotheses",
     "range", "P1"),
    ("results/synteny", "*", "origin", "P1",
     "S8: flanking-gene neighbourhoods and the consensus paralog caller",
     "retention", "P1: its caller is cited by retention, but its result is "
     "the paralogon"),
    ("results/reconciliation", "*", "origin", "P1",
     "S13: gene-tree/species-tree reconciliation, dated duplications",
     "retention", "P1: its loss audit is a step in dating, not a count"),
    ("results/duplication", "*", "origin", "P1",
     "S16: 2R ohnology, the 3R teleost copies and the copy-number landscape",
     "retention", "P1"),
    ("results/gene_architecture", "*", "origin", "P1",
     "S21: exon/intron architecture and the shared-ancestral-intron test",
     "archive", "P1: its headline is the shared intron core; the "
     "fragment-terminus half is cited by the archive paper"),
    ("results/methods", "reconciliation_loss_audit.tsv", "origin", "P5s",
     "S19: where the reconciliation's loss inference ran out", "retention",
     "P2: the limit of the origin paper's own dating"),
    ("results/methods", "synteny_power.tsv", "origin", "P5s",
     "S19: synteny reach against synteny accuracy", "retention", "P2"),
    ("results/methods", "tree_resolution.tsv", "origin", "P5s",
     "S19: the size of the AU test's confidence set", "range", "P2"),
    ("results/methods", "inference_summary.json", "origin", "P5s",
     "S19: the inference-limit summary", "constraint",
     "P2: three of its five limits are the origin paper's"),
    ("results/supplementary", "figures/SuppFig1_*", "origin", "P5s",
     "S24: the representative alignment and the columns the tree saw",
     "constraint", "P5: it draws msa_v2, which is primary here"),
    # ---------------------------------------------------------- retention
    ("results/genome_manifest.tsv", "*", "retention", "P3",
     "S4: the declared 309-genome vertebrate denominator", "archive",
     "P3: both papers use it and each restates it; it was built for the "
     "count"),
    ("results/genome_manifest_notes.md", "*", "retention", "P3",
     "S4: the scope document rendered from the manifest", "archive", "P3"),
    ("results/s5_baits", "*", "retention", "P2",
     "S5: the 38-bait panel with the RyR presence control", "range",
     "P2: the RyR control is what makes a missing cell searchable"),
    ("results/genome_ledger", "*", "retention", "P1",
     "S5: the 309-genome ledger, one row per genome and paralogue",
     "archive", "P1"),
    ("results/loss_dynamics", "*", "retention", "P1",
     "S15a: the character matrix and the loss instrument", "origin", "P1"),
    ("results/loss_counts", "*", "retention", "P1",
     "S15b: the Dollo counts and the sensitivity matrix", "origin", "P1"),
    ("results/methods", "contiguity_*", "retention", "P5s",
     "S19: the ledger's false-negative rate against contiguity", "archive",
     "P2: a retention count is only a measurement beside its miss rate"),
    ("results/methods", "absence_floor.tsv", "retention", "P5s",
     "S19: the contiguity floor calibration", "archive", "P2"),
    ("results/methods", "floor_clade_composition.tsv", "retention", "P5s",
     "S19: what each floor costs, by clade", "archive", "P2"),
    ("results/methods", "residual_neighbourhood.tsv", "retention", "P5s",
     "S19: misses whose whole region is gone", "archive", "P2"),
    ("results/methods", "panel_*", "retention", "P5s",
     "S19: the exact bait-panel ablation", "range",
     "P2: the panel is the retention sweep's instrument"),
    ("results/methods", "figures/fig_s19_contiguity.*", "retention", "P5s",
     "S19: the contiguity figure", "archive", "P2"),
    ("results/methods", "figures/fig_s19_panel.*", "retention", "P5s",
     "S19: the panel ablation figure", "range", "P2"),
    ("results/methods", "methods_stats.json", "retention", "P5s",
     "S19: the task's stats record", "archive",
     "P5s: the record goes with the half the thesis kept it beside"),
    ("results/methods", "report.md", "retention", "P5s",
     "S19: the task's rendered report", "archive", "P5s"),
    # ---------------------------------------------------------- archive
    ("results/census_v4", "*", "archive", "P1",
     "S5b: census v3 plus the genome models, with how a database holds each",
     "retention", "P1: its result is the db-status of each locus"),
    ("results/annotation_bugs", "*", "archive", "P1",
     "S10: two annotation failures proved at the exon", "retention", "P1"),
    ("results/expression", "*", "archive", "P2",
     "S12: RNA-seq junction evidence at the loci the annotation loses",
     "retention", "P2: the transcription evidence is the control that the "
     "lost loci are genes"),
    ("results/annotation_audit", "*", "archive", "P1",
     "S18: the locus and protein-record audit, and 297 corrections",
     "retention", "P1"),
    ("results/methods", "census_growth.tsv", "archive", "P5s",
     "S19: what each search channel added to the census", "range",
     "P1: the question is what a database-only search would have missed"),
    ("results/methods", "method_*", "archive", "P5s",
     "S19: channel sets, pairwise overlap and cost", "range", "P1"),
    ("results/methods", "head_to_head.tsv", "archive", "P5s",
     "S19: the sweep against the database records", "range", "P1"),
    ("results/methods", "gene_recovery*.tsv", "archive", "P5s",
     "S19: per-gene recovery by channel, cell and scope", "retention", "P1"),
    ("results/methods", "db_status_combined.tsv", "archive", "P5s",
     "S19: how databases hold every demonstrated gene", "retention", "P1"),
    ("results/methods", "contribution_summary.json", "archive", "P5s",
     "S19: the contribution summary", "range", "P1"),
    ("results/methods", "figures/fig_s19_contribution.*", "archive", "P5s",
     "S19: the contribution figure", "range", "P1"),
    # ---------------------------------------------------------- constraint
    ("results/selection", "*", "constraint", "P1",
     "S9: codon models across the vertebrate family", "origin",
     "P1: its result is purifying selection, not the duplication's timing"),
    ("results/constraint", "*", "constraint", "P1",
     "S17: per-residue constraint and the variant classifier", "ligand",
     "P1: the ligand paper asks one module's question of these layers"),
    ("results/methods", "codon_model_power.tsv", "constraint", "P5s",
     "S19: the power of the codon models", "origin", "P2"),
    ("results/methods", "saturation.tsv", "constraint", "P5s",
     "S19: pairwise synonymous saturation", "origin", "P2"),
    ("results/supplementary", "figures/SuppFig3_*", "constraint", "P5s",
     "S24: the within-paralogue alignments", "ligand", "P5"),
    ("results/supplementary", "figures/SuppFig4_*", "constraint", "P5s",
     "S24: the trimmed codon alignment", "origin", "P5"),
    ("results/supplementary", "figures/SuppFig5_*", "constraint", "P5s",
     "S24: the constraint map painted on the channel", "ligand", "P5"),
    ("results/supplementary", "figures/SuppFig6_*", "constraint", "P5s",
     "S24: every labelled variant on the structure", "ligand", "P5"),
    # ---------------------------------------------------------- ligand
    ("results/ligand_site", "*", "ligand", "P1",
     "S22: the ligand core against the pore, the shells and the PLC panel",
     "constraint", "P1: one module's question, with its own controls (P2) "
     "and four figures (P4)"),
    ("results/supplementary", "figures/SuppFig2_*", "ligand", "P5s",
     "S24: the ligand core and the pore at residue resolution",
     "constraint", "P5"),
    # ---------------------------------------------------------- excluded
    ("results/s0_baseline", "*", "-", "excluded",
     "S0: the literature baseline and the reference table every paper "
     "resolves its citations against", "-",
     "an input all six papers read and none claims as a result"),
    ("results/toolchain_manifest.txt", "*", "-", "excluded",
     "S1: the tool versions every task shells out to", "-",
     "deposited with every paper; not a result"),
    ("results/session_live.json", "*", "-", "excluded",
     "the dashboard's progress panel", "-", "session state"),
    ("results/2026-08-18_091322_itpr1_itpr2_itpr3", "*", "-", "excluded",
     "S0 smoke test of the search application", "-",
     "an application test, not an analysis"),
    ("results/2026-08-18_092259_itpr1a_itpr1b_itpr2", "*", "-", "excluded",
     "S0 smoke test of the search application", "-",
     "an application test, not an analysis"),
    ("results/2026-08-18_123301_itpr1_itpr2_itpr3", "*", "-", "excluded",
     "S0 smoke test of the search application", "-",
     "an application test, not an analysis"),
    ("results/2026-08-18_123618_s0_smoke_human", "*", "-", "excluded",
     "S0 smoke test of the search application", "-",
     "an application test, not an analysis"),
    ("results/supplementary", "figure_*", "-", "excluded",
     "S24: the audit of the single manuscript's figures", "-",
     "each paper builds its own figure audit"),
    ("results/supplementary", "supp*", "-", "excluded",
     "S24: the manuscript's supplementary-figure record", "-",
     "each paper records its own figure manifest"),
    ("results/supplementary", "report.md", "-", "excluded",
     "S24: the rendered report of the manuscript's figure audit", "-",
     "each paper builds its own figure audit"),
]

#: The thesis chapter each paper corresponds to, for the cross-check against
#: `thesis/chapter_assignment.tsv`. None: the chapter's entries are excluded
#: from the series (the baseline and the methods chapter's infrastructure).
THESIS_CHAPTER_PAPER: dict[int, str | None] = {
    2: None, 3: "range", 4: "retention", 5: "range", 6: "origin",
    7: "origin", 8: "origin", 9: "retention", 10: "constraint",
    11: "constraint", 12: "ligand", 13: "archive", 14: None,
}

#: Every place this assignment departs from the thesis's, with the reason.
#: A departure that is not declared here fails the build, and so does a
#: declared departure that no longer departs.
DEPARTURES: dict[str, str] = {
    "results/methods":
        "The thesis keeps S19 whole in Chapter 4 beside the search it "
        "measures (T4). Here P2 forces a split: the false-negative rate is "
        "the retention paper's central control, the drift audit is the "
        "range paper's, the channel contribution is the archive paper's "
        "result, and the inference limits bound the origin and constraint "
        "papers' own claims. This is the departure the thesis's Appendix "
        "E.2 predicted.",
    "results/census_v4":
        "The thesis places census v4 with the vertebrate sweep that made "
        "it. Its result, how a database holds each recovered locus, is the "
        "archive paper's question, and the retention paper never uses it.",
    "results/structures":
        "The thesis places S11 with the channel. In the series its "
        "load-bearing result is the structural family call, a third "
        "instrument confirming the range paper's ITPR/RyR separation; the "
        "constraint paper paints S17's own tables and does not need it.",
    "results/supplementary":
        "The thesis places S24 with its methods chapter. Each supplementary "
        "figure draws one paper's data, so P5 sends each to that paper and "
        "the manuscript's own figure audit is excluded.",
}
