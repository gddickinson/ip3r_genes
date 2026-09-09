"""S25 — every figure the thesis places, and where it comes from.

One row per figure: the chapter it belongs to, the slug the text places it
by, the committed source stem it is copied from (relative to the repository
root, without extension), and a one-line stub for `figure_manifest.tsv`.

Three rules are enforced from this table by `s25_figures.py`:

* **T7** — a figure is *copied*, never re-plotted, and must exist as both png
  (what GitHub renders) and pdf (what the typeset document uses). A missing
  either way is a build error.
* **A slug is used once.** Two results directories both commit a figure
  called `profile_separation`, and two more both commit `copy_number`; the
  slug here is the thesis's own name for a figure and is unique, so a
  chapter cannot silently place its neighbour's panel.
* **Nothing committed is quietly left out.** Every png under a
  `results/*/figures/` directory that the assignment places in a chapter must
  appear here or in `UNPLACED`. That is rule T6 again: a figure that exists
  and is not in the thesis has to be a decision rather than an oversight.

Numbering is per chapter and in order of first placement in the text, so a
chapter can be rewritten without renumbering anything anywhere else.
"""

from __future__ import annotations

R = "results"
D = "docs/figures"

#: (chapter, slug, source stem, stub caption)
FIGURES: list[tuple[int, str, str, str]] = [
    # ------------------------------------------------ 1. the receptor
    (1, "signal_hierarchy", f"{D}/signal_hierarchy",
     "Where the receptor sits in the phosphoinositide pathway"),
    (1, "channel_structure", f"{D}/channel_structure",
     "The channel, measured on PDB 6DQN"),
    (1, "gating_logic", f"{D}/gating_logic",
     "What opens the channel and what closes it"),
    (1, "domain_architecture", f"{D}/domain_architecture",
     "The domain architecture the two families share"),
    (1, "discovery_timeline", f"{D}/discovery_timeline",
     "Fifty years of the receptor, as the literature records it"),
    # ------------------------------------------------ 2. the baseline
    (2, "review_gene_architecture", f"{D}/gene_architecture",
     "Exon counts and genomic spans for the three human paralogues"),
    (2, "review_taxonomic_range", f"{D}/taxonomic_range",
     "What the databases said the family's range was"),
    (2, "review_family_separation", f"{D}/family_separation",
     "Identity between and within the two families"),
    (2, "review_conservation", f"{D}/conservation_profile",
     "Per-column conservation on human ITPR1 numbering"),
    (2, "review_alignment_windows", f"{D}/alignment_windows",
     "The measured functional sites, in alignment windows"),
    (2, "review_regulation", f"{D}/regulation_map",
     "The curated regulator map"),
    (2, "review_disease", f"{D}/disease_map",
     "The curated disease-variant map, by how far each source localises"),
    # ---------------------------------- 3. two families, one architecture
    (3, "census_space", f"{R}/census_v2/figures/census_space",
     "The enumerated search space, by signature"),
    (3, "census_lineage", f"{R}/census_v2/figures/census_lineage",
     "The census by lineage, with both family calls"),
    (3, "census_lengths", f"{R}/census_v2/figures/census_lengths",
     "Length distributions of the two families"),
    (3, "census_margin", f"{R}/census_v2/figures/census_margin",
     "The labelled-bait identity margin, inside and outside the no-call band"),
    (3, "profile_separation_vert",
     f"{R}/census_v3/figures/profile_separation",
     "Profile separation with the no-call band drawn"),
    (3, "instrument_agreement",
     f"{R}/census_v3/figures/instrument_agreement",
     "Which instrument calls each record"),
    (3, "jackhmmer_convergence",
     f"{R}/census_v3/figures/jackhmmer_convergence",
     "Convergence against sister-family contamination"),
    (3, "proteome_copy_number",
     f"{R}/census_v3/figures/proteome_copy_number",
     "Copy number per proteome against annotation depth"),
    # ------------------------------------- 4. the vertebrate sweep, S19
    (4, "ledger_by_class", f"{R}/genome_ledger/figures/ledger_by_class",
     "The three-paralogue ledger across 309 genomes, by class"),
    (4, "ledger_status", f"{R}/genome_ledger/figures/ledger_status",
     "Cell status across the sweep"),
    (4, "contiguity_confound",
     f"{R}/genome_ledger/figures/contiguity_confound",
     "Recovery against assembly contiguity, binned on the bar"),
    (4, "ledger_copy_number", f"{R}/genome_ledger/figures/copy_number",
     "Copy number per genome across the vertebrate sweep"),
    (4, "s19_contiguity", f"{R}/methods/figures/fig_s19_contiguity",
     "The measured false-negative rate against contiguity"),
    (4, "s19_panel", f"{R}/methods/figures/fig_s19_panel",
     "What each bait was worth: the panel ablation"),
    (4, "s19_contribution", f"{R}/methods/figures/fig_s19_contribution",
     "What each search channel contributed"),
    (4, "s19_drift", f"{R}/methods/figures/fig_s19_drift",
     "The kill criterion scored as a classifier of drift"),
    # --------------------------------------- 5. outside the vertebrates
    (5, "range_by_phylum", f"{R}/s20_sweep/figures/range_by_phylum",
     "The family across the eukaryotes, as a fraction of swept proteomes"),
    (5, "profile_separation_euk",
     f"{R}/s20_sweep/figures/profile_separation",
     "The same separation outside the vertebrates"),
    (5, "plant_fungal_chase", f"{R}/s20_sweep/figures/plant_fungal_chase",
     "Every plant and fungal record chased individually"),
    (5, "jackhmmer_s20", f"{R}/s20_sweep/figures/jackhmmer_s20",
     "Per-group convergence and its sister-family trace"),
    (5, "nonvert_copy_number", f"{R}/s23_scope/figures/copy_number",
     "Copy number outside the vertebrates"),
    (5, "absence_at_genome", f"{R}/s23_scope/figures/absence_at_genome",
     "Each absence claim beside the controlled genomes behind it"),
    (5, "identity_floor", f"{R}/s23_scope/figures/identity_floor",
     "The measured identity floor and the populations it separates"),
    (5, "span_inflation", f"{R}/s23_scope/figures/span_inflation",
     "Locus span against coding footprint, by group"),
    # ------------------------------------------ 6. the alignment and tree
    (6, "msa_conservation", f"{R}/msa_v2/figures/msa_conservation",
     "Per-column conservation with the architecture mapped through the "
     "alignment"),
    (6, "msa_identity_heatmap",
     f"{R}/msa_v2/figures/msa_identity_heatmap",
     "All-pairs identity across the representative set"),
    (6, "msa_group_identity", f"{R}/msa_v2/figures/msa_group_identity",
     "Identity within and between groups"),
    (6, "msa_coverage", f"{R}/msa_v2/figures/msa_coverage",
     "Per-sequence coverage of the alignment"),
    (6, "supp_representative_alignment",
     f"{R}/supplementary/figures/SuppFig1_representative_alignment",
     "The alignment, and the columns the tree actually saw"),
    (6, "tree_ml_rooted", f"{R}/phylogeny/figures/tree_ml_rooted",
     "The rooted maximum-likelihood phylogram"),
    (6, "sister_au", f"{R}/phylogeny/figures/sister_au",
     "The AU test on the three sister hypotheses"),
    (6, "support_profile", f"{R}/phylogeny/figures/support_profile",
     "Node support as a joint condition, not two marginals"),
    (6, "paralog_placement", f"{R}/phylogeny/figures/paralog_placement",
     "Where the tree places each tip against its census label"),
    # ----------------------------------------------- 7. where they came from
    (7, "synteny_paralogon", f"{R}/synteny/figures/synteny_paralogon",
     "The three human neighbourhoods with their shared ohnologue families"),
    (7, "synteny_pair_classes",
     f"{R}/synteny/figures/synteny_pair_classes",
     "Every pair class against its matched random-window control"),
    (7, "synteny_caller", f"{R}/synteny/figures/synteny_caller",
     "The consensus paralogue caller, calibrated and swept"),
    (7, "synteny_clade_decay", f"{R}/synteny/figures/synteny_clade_decay",
     "How far a neighbourhood travels"),
    (7, "s16_copy_number", f"{R}/duplication/figures/s16_copy_number",
     "Copy number grouped by whole-genome-duplication status"),
    (7, "s16_paralogon", f"{R}/duplication/figures/s16_paralogon",
     "The 2R test against its real-window permutation null"),
    (7, "s16_dcs", f"{R}/duplication/figures/s16_dcs",
     "Double-conserved synteny for the teleost copies"),
    (7, "s16_blocks", f"{R}/duplication/figures/s16_blocks",
     "Cross-anchor block identity, and whether the anchors agree"),
    (7, "recon_dated_backbone",
     f"{R}/reconciliation/figures/recon_dated_backbone",
     "The dated backbone with each calibration's spread as a band"),
    (7, "recon_matrix", f"{R}/reconciliation/figures/recon_matrix",
     "The reconciliation matrix over topologies and variants"),
    (7, "recon_losses", f"{R}/reconciliation/figures/recon_losses",
     "The loss audit against the genome ledger"),
    (7, "recon_cyclostome",
     f"{R}/reconciliation/figures/recon_cyclostome",
     "The six cyclostome loci, and the long-branch check behind them"),
    # ------------------------------------------------------- 8. the gene
    (8, "architecture_by_paralog",
     f"{R}/gene_architecture/figures/architecture_by_paralog",
     "Exon count and genomic span on separate axes"),
    (8, "intron_positions",
     f"{R}/gene_architecture/figures/intron_positions",
     "Intron positions by alignment column and phase"),
    (8, "junction_quality",
     f"{R}/gene_architecture/figures/junction_quality",
     "The aligner's own junction error rate, drawn on the complement"),
    (8, "fragments_and_duplicates",
     f"{R}/gene_architecture/figures/fragments_and_duplicates",
     "Where a fragmentary annotation stops, and the duplication detector"),
    # ------------------------------------------------- 9. retention
    (9, "character_matrix",
     f"{R}/loss_dynamics/figures/s15_character_matrix",
     "The character matrix, ordered by contiguity with the bar drawn"),
    (9, "reconstruction", f"{R}/loss_dynamics/figures/s15_reconstruction",
     "Reference coverage reassembled across contigs"),
    (9, "integrity", f"{R}/loss_dynamics/figures/s15_integrity",
     "Lesion density, and the confounder that predicts it"),
    (9, "synteny_reach", f"{R}/loss_dynamics/figures/s15_synteny_reach",
     "What a trace region gives the synteny caller to work with"),
    (9, "reconstruction_bar",
     f"{R}/loss_counts/figures/reconstruction_bar",
     "The calibrated bar with the gap's two edges drawn"),
    (9, "sensitivity_matrix",
     f"{R}/loss_counts/figures/sensitivity_matrix",
     "The loss count across four axes of sensitivity"),
    (9, "mk_profile", f"{R}/loss_counts/figures/mk_profile",
     "The likelihood along a rate grid, for every model and axis"),
    (9, "lesion_strata", f"{R}/loss_counts/figures/lesion_strata",
     "Lesion density stratified by class, with its contiguity control"),
    # ------------------------------------------------- 10. selection
    (10, "omega_by_paralog",
     f"{R}/selection/figures/s9_omega_by_paralog",
     "Per-paralogue omega on a logarithmic axis"),
    (10, "branch_contrast", f"{R}/selection/figures/s9_branch_contrast",
     "The two-ratio and branch-site contrasts"),
    (10, "dnds_saturation", f"{R}/selection/figures/s9_dnds_saturation",
     "Pairwise saturation, drawn log-log with the neutral diagonal"),
    (10, "bs_restarts", f"{R}/selection/figures/s9_bs_restarts",
     "Branch-site fits across four initial omega"),
    (10, "supp_codon_alignment",
     f"{R}/supplementary/figures/SuppFig4_codon_alignment",
     "The trimmed codon alignment behind every estimate"),
    # ------------------------------------------------- 11. the channel
    (11, "s11_panel", f"{R}/structures/figures/s11_panel",
     "The structural panel: references, states and negative controls"),
    (11, "s11_tm_calibration",
     f"{R}/structures/figures/s11_tm_calibration",
     "TM-align scores with both published bars drawn"),
    (11, "s11_plddt_domains",
     f"{R}/structures/figures/s11_plddt_domains",
     "Per-domain confidence, against an outside-domains contrast"),
    (11, "s11_afdb_coverage",
     f"{R}/structures/figures/s11_afdb_coverage",
     "AlphaFold coverage against record length"),
    (11, "s17_channel_profile",
     f"{R}/constraint/figures/s17_channel_profile",
     "Constraint along the channel, binned within element boundaries"),
    (11, "s17_elements", f"{R}/constraint/figures/s17_elements",
     "Constraint by element, with the composition-free metric beside it"),
    (11, "s17_functional_sites",
     f"{R}/constraint/figures/s17_functional_sites",
     "The measured functional residues against their own elements"),
    (11, "s17_variant_classifier",
     f"{R}/constraint/figures/s17_variant_classifier",
     "Four conservation layers as variant classifiers"),
    (11, "supp_paralog_alignments",
     f"{R}/supplementary/figures/SuppFig3_paralog_alignments",
     "The within-paralogue alignments the map is computed on"),
    (11, "supp_constraint_on_channel",
     f"{R}/supplementary/figures/SuppFig5_constraint_on_channel",
     "The constraint map painted onto the channel"),
    (11, "supp_variants_on_structure",
     f"{R}/supplementary/figures/SuppFig6_variants_on_structure",
     "Every labelled variant, and the per-element enrichment test"),
    # ---------------------------------------------- 12. the ligand site
    (12, "s22_modules", f"{R}/ligand_site/figures/s22_fig1_modules",
     "The core against the pore, under both pore definitions"),
    (12, "s22_shells", f"{R}/ligand_site/figures/s22_fig2_shells",
     "Constraint against distance from the ligand"),
    (12, "s22_omega", f"{R}/ligand_site/figures/s22_fig3_omega",
     "Per-site rates by module and by shell"),
    (12, "s22_lineage", f"{R}/ligand_site/figures/s22_fig4_lineage",
     "The lineage test beside the shift it could have detected"),
    (12, "supp_labelled_positions",
     f"{R}/supplementary/figures/SuppFig2_labelled_positions",
     "The two modules at residue resolution in all three paralogues"),
    # -------------------------------------------------- 13. the archive
    (13, "s18_by_source",
     f"{R}/annotation_audit/figures/s18_fig1_by_source",
     "Locus state by annotation source, raw and above the contiguity bar"),
    (13, "s18_calibration",
     f"{R}/annotation_audit/figures/s18_fig2_calibration",
     "The completeness bar drawn inside the distribution it sits in"),
    (13, "s18_family_vs_control",
     f"{R}/annotation_audit/figures/s18_fig3_family_vs_control",
     "This family against its sister in the same assemblies"),
    (13, "s18_protein_side",
     f"{R}/annotation_audit/figures/s18_fig4_protein_side",
     "The protein records: what the name claims against what the sequence "
     "is"),
    (13, "exon_tracks", f"{R}/annotation_bugs/figures/exon_tracks",
     "Two loci at true genomic width, exons never widened"),
    (13, "annotation_loss", f"{R}/annotation_bugs/figures/annotation_loss",
     "Annotation loss across every recovered locus"),
    (13, "fragment_tiling",
     f"{R}/annotation_bugs/figures/fragment_tiling",
     "The annotated proteins tiled back onto the genome's own loci"),
    (13, "case_validation", f"{R}/annotation_bugs/figures/validation",
     "The five independent checks on each case"),
    (13, "s12_junctions", f"{R}/expression/figures/s12_junctions",
     "Every junction, crossed or not, with the decoy floor drawn"),
    (13, "s12_detection", f"{R}/expression/figures/s12_detection",
     "Detection per locus against its own reversed decoy"),
    (13, "s12_gap_coverage", f"{R}/expression/figures/s12_gap_coverage",
     "Read coverage over the sequence the annotation loses"),
    (13, "s12_instruments", f"{R}/expression/figures/s12_instruments",
     "The deposit cross-check and its genomic negative control"),
    # -------------------------------------------------- 14. methods
    (14, "census_growth", f"{R}/census_v2/figures/census_growth",
     "How the census grew across its editions, and on what evidence"),
]

#: Committed figures deliberately not placed, with the reason.
#:
#: The only eight are the app's own bundle plots from the four S0 smoke-test
#: runs. They are written automatically by `src/utils/results_writer.py` for
#: any search bundle, are drawn from whatever that one search returned rather
#: than from a committed analysis table, and say nothing the census figures do
#: not say better. Rule T7 would admit them — they are committed beside their
#: data — but they are not results, and the smoke test's result is that the
#: three databases answered at all.
UNPLACED: dict[str, str] = {
    f"{R}/{run}/figures/{stem}":
        "app bundle plot from an S0 smoke-test search, not an analysis figure"
    for run in ("2026-08-18_091322_itpr1_itpr2_itpr3",
                "2026-08-18_092259_itpr1a_itpr1b_itpr2",
                "2026-08-18_123301_itpr1_itpr2_itpr3",
                "2026-08-18_123618_s0_smoke_human")
    for stem in ("length_histogram", "sources_bar")
}
