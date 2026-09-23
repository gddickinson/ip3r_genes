"""S26 paper configuration: Retention across 309 vertebrate genomes.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/retention/`.
"""

TITLE = 'No vertebrate has lost an IP₃ receptor paralogue'
SHORT_TITLE = 'Retention across 309 vertebrate genomes'
QUESTION = 'Has any vertebrate lineage lost one of its three IP₃ receptor genes?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('census', 'main', 'genome_ledger/figures/ledger_by_class', 'The three-paralogue census across 309 genomes'),
    ('false_negatives', 'main', 'methods/figures/fig_s19_contiguity', 'The miss rate against assembly contiguity'),
    ('panel', 'main', 'methods/figures/fig_s19_panel', 'The bait panel ablated'),
    ('matrix', 'main', 'loss_dynamics/figures/s15_character_matrix', 'The character matrix'),
    ('reconstruction', 'main', 'loss_dynamics/figures/s15_reconstruction', 'Reassembling a shattered gene'),
    ('sensitivity', 'main', 'loss_counts/figures/sensitivity_matrix', 'What it takes to manufacture a loss'),
    ('status', 'ed', 'genome_ledger/figures/ledger_status', 'The ledger status of every cell'),
    ('confound', 'ed', 'genome_ledger/figures/contiguity_confound', 'Recovery against assembly contiguity'),
    ('synteny_reach', 'ed', 'loss_dynamics/figures/s15_synteny_reach', 'Why synteny cannot reach the undecided cells'),
    ('reconstruction_bar', 'ed', 'loss_counts/figures/reconstruction_bar', 'The reconstruction bar moved across its gap'),
    ('mk_profile', 'ed', 'loss_counts/figures/mk_profile', 'The rate models refused'),
    ('integrity', 'ed', 'loss_dynamics/figures/s15_integrity', 'Reading-frame integrity and its confounders'),
    ('lesion_strata', 'ed', 'loss_counts/figures/lesion_strata', 'The lesion signal by class'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('the ryanodine receptor presence control in every genome', 'results/genome_ledger/control_failures.tsv'),
    ('the reconstruction bar calibrated against decoy regions', 'results/loss_dynamics/recon_calibration.tsv'),
    ('the false-negative rate on genes independently known present', 'results/methods/contiguity_cells.tsv'),
    ('the bait-panel ablation that reproduces the ledger', 'results/methods/panel_recall.tsv'),
    ('every loss the settings can manufacture', 'results/loss_counts/manufactured_losses.tsv'),
    ('the identity-matched within-genome lesion test', 'results/loss_dynamics/integrity_tests.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = ["C01", "C26", "C85", "RE03", "RE01"]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = """
If no other paper in the series is published, this one still claims that
across a declared scope of 309 vertebrate genome assemblies, searched with a
positive control that fired in every genome, no genome x paralogue cell of
927 reaches a state that licenses a loss and Dollo parsimony places zero
losses on the IP3 receptor family character; that the search behind that
zero misses 140 of 923 known-present cells overall but only 5 of 563 in
assemblies contiguous enough to hold the gene, measured on its own output;
that a family-level loss can be manufactured in only 2 of 32 analytical
settings; and that the full bait-panel ablation reproduces the ledger cell
for cell. None of this needs the origin paper's tree. Paralogue labels come
from the bait panel's own labelled sequences, the cell calls use no synteny
or gene tree, and the family-level count, the paper's strongest claim, uses
no paralogue label at all.
"""
STANDALONE_CLAIMS: list[str] = ["C28", "C85", "C86", "C91", "C35",
                                 "RE25", "RE24", "C92", "C45"]
