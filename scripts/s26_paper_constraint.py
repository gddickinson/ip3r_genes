"""S26 paper configuration: Constraint along the receptor.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/constraint/`.
"""

TITLE = ('Purifying selection on the IP₃ receptor is strongest at the gate '
         'and filter and weakest in a luminal loop beside them')
SHORT_TITLE = 'Constraint along the receptor'
QUESTION = 'Where along the IP₃ receptor does purifying selection act?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('omega', 'main', 'selection/figures/s9_omega_by_paralog', 'Purifying selection in all three paralogues'),
    ('branches', 'main', 'selection/figures/s9_branch_contrast', 'No branch escapes it'),
    ('elements', 'main', 'constraint/figures/s17_elements', 'Constraint by element'),
    ('profile', 'main', 'constraint/figures/s17_channel_profile', 'Constraint along the channel'),
    ('sites', 'main', 'constraint/figures/s17_functional_sites', 'The measured functional sites'),
    ('classifier', 'main', 'constraint/figures/s17_variant_classifier', 'Constraint as a variant classifier'),
    ('codon_checks', 'ed', ['selection/figures/s9_dnds_saturation', 'selection/figures/s9_bs_restarts'], 'Saturation and the branch-site restarts'),
    ('paralog_alignments', 'supp', 'supplementary/figures/SuppFig3_paralog_alignments', 'The within-paralogue alignments'),
    ('codon_alignment', 'supp', 'supplementary/figures/SuppFig4_codon_alignment', 'The trimmed codon alignment'),
    ('constraint_on_channel', 'supp', 'supplementary/figures/SuppFig5_constraint_on_channel', 'The constraint map on the channel'),
    ('variants_on_structure', 'supp', 'supplementary/figures/SuppFig6_variants_on_structure', 'Every labelled variant on the structure'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('the composition-free metric and the curated-subset recomputation', 'results/constraint/metric_controls.tsv'),
    ('the curated residue-level variants the harvest must recover', 'results/constraint/curated_variant_control.tsv'),
    ('the whole-protein and shared-position contrasts', 'results/constraint/variant_constraint_test.tsv'),
    ('the transcript numbering check', 'results/constraint/clinvar_transcripts.tsv'),
    ('the branch-site restarts', 'results/selection/bs_restarts.tsv'),
    ('synonymous saturation', 'results/methods/saturation.tsv'),
    ('the power of the codon models', 'results/methods/codon_model_power.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = ["C01", "T53", "T54", "C115", "C116", "C117",
                            "C118", "CO22", "CO23", "CO24", "C127"]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = (
    "Read alone, this paper claims where purifying selection acts on the "
    "vertebrate IP3 receptor, measured on two independent axes it builds "
    "itself. On the rate axis, all three paralogues are held far below "
    "neutrality (one-ratio omega 0.02 to 0.04), ITPR1 most tightly and with "
    "selection intensified relative to its sisters, and no site is "
    "positively selected. On the conservation axis, computed on 249 to 265 "
    "orthologues per paralogue, the gate and the selectivity filter are the "
    "most conserved elements, the gate is identical between the three human "
    "copies, and a geometrically located luminal loop inside the pore domain "
    "is the least conserved sequence in the receptor. Conservation "
    "classifies ClinVar missense positions, best on the family-wide layer, "
    "and the per-residue tables are a stratification resource for 1,546 "
    "uncertain variants. Its codon sequences, orthologue sets, controls and "
    "variant harvest are all its own; the companion papers supply only the "
    "tree it conditions on and the genome sweep its orthologues are drawn "
    "from, both of which it declares and uses as inputs."
)
STANDALONE_CLAIMS: list[str] = ["C100", "C101", "C102", "C105", "C119",
                                 "C121", "CO57", "C129", "C131", "C127",
                                 "C128", "C135"]
