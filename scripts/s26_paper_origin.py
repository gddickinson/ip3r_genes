"""S26 paper configuration: Where ITPR1, ITPR2 and ITPR3 came from.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/origin/`.
"""

TITLE = 'The three vertebrate IP₃ receptors arose in two separate early vertebrate duplications'
SHORT_TITLE = 'Where ITPR1, ITPR2 and ITPR3 came from'
QUESTION = 'How did one ancestral gene become the three vertebrate IP₃ receptors?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('tree', 'main', 'phylogeny/figures/tree_ml_rooted', 'The family tree rooted on the ryanodine receptors'),
    ('sister', 'main', 'phylogeny/figures/sister_au', 'The sister-pair test'),
    ('paralogon', 'main', 'duplication/figures/s16_paralogon', 'The 2R paralogon test'),
    ('dated', 'main', 'reconciliation/figures/recon_dated_backbone', 'Where the duplications sit on the dated species tree'),
    ('neighbourhood', 'main', 'synteny/figures/synteny_paralogon', "Each paralogue's genomic neighbourhood"),
    ('introns', 'main', 'gene_architecture/figures/intron_positions', 'The shared intron core'),
    ('teleost', 'main', 'duplication/figures/s16_dcs', 'The teleost genome duplication'),
    ('alignment', 'ed', ['msa_v2/figures/msa_conservation', 'msa_v2/figures/msa_coverage'], 'The representative alignment'),
    ('identity', 'ed', 'msa_v2/figures/msa_identity_heatmap', 'All-pairs identity across the representatives'),
    ('group_identity', 'ed', 'msa_v2/figures/msa_group_identity', 'Identity within and between groups'),
    ('support', 'ed', ['phylogeny/figures/support_profile', 'phylogeny/figures/paralog_placement'], 'Node support and paralogue placement'),
    ('synteny', 'ed', ['synteny/figures/synteny_pair_classes', 'synteny/figures/synteny_caller', 'synteny/figures/synteny_clade_decay'], 'Neighbourhood conservation against its null'),
    ('reconciliation', 'ed', ['reconciliation/figures/recon_matrix', 'reconciliation/figures/recon_cyclostome', 'reconciliation/figures/recon_losses'], 'The reconciliation matrix and the cyclostome loci'),
    ('duplication', 'ed', ['duplication/figures/s16_copy_number', 'duplication/figures/s16_blocks'], 'Copy number and the ancestral block'),
    ('architecture', 'ed', ['gene_architecture/figures/architecture_by_paralog', 'gene_architecture/figures/junction_quality', 'gene_architecture/figures/fragments_and_duplicates'], "The gene's architecture and its instrument"),
    ('representative_alignment', 'supp', 'supplementary/figures/SuppFig1_representative_alignment', 'The alignment and the columns the tree saw'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('the AU test of all three sister hypotheses', 'results/phylogeny/au_test.tsv'),
    ('the --bnni model-violation check', 'results/phylogeny/bnni_comparison.tsv'),
    ('the matched random-window synteny null', 'results/synteny/pair_stats.tsv'),
    ('the real-window paralogon permutation null with the ryanodine receptors as positive control', 'results/duplication/paralogon_test.tsv'),
    ('the 309-genome replication null', 'results/duplication/replication_null.tsv'),
    ('every rooting of the vertebrate subtree scored', 'results/reconciliation/rooting_check.tsv'),
    ('the shared-intron null with the ryanodine receptors as control', 'results/gene_architecture/shared_intron_summary.tsv'),
    ("the AU test's confidence set", 'results/methods/tree_resolution.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = ["C51", "OR01", "OR02", "OR03", "OR04", "OR05",
                            "C181", "C182"]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = """
Read alone, this paper claims that the three vertebrate IP₃ receptor genes
are ohnologues made by two duplications on two different branches of the
early vertebrate tree: *ITPR1* split from the ancestor of *ITPR2* and *ITPR3*
on the vertebrate stem, and *ITPR2* from *ITPR3* on the gnathostome stem. It
rests on its own tree and topology test (the *ITPR2* plus *ITPR3* pairing is
the only one retained), its own dated reconciliation, its own neighbourhood
paralogy test replicated across 309 genomes against a real-window null, and
its own intron-position test, each with the ryanodine receptors run through
the same instrument as the control. The teleost result, that one ancestral
duplication doubled *ITPR1* alone, also needs nothing from the other papers.
It uses the companion range paper only for the census it samples its
representatives from, and states that sample in its own methods.
"""
STANDALONE_CLAIMS: list[str] = ["C60", "OR06", "OR07", "C82", "C80",
                                 "C67", "C69", "C70", "C207", "C76",
                                 "C78", "C218"]
