"""S26 paper configuration: How the family is recorded.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/archive/`.
"""

TITLE = 'Four in five vertebrate IP₃ receptor genes cannot be reached from any protein database'
SHORT_TITLE = 'How the family is recorded'
QUESTION = 'How completely do the public databases record the IP₃ receptor genes that vertebrate genomes carry?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('by_source', 'main', 'annotation_audit/figures/s18_fig1_by_source', 'Annotation state by annotation source'),
    ('family_control', 'main', 'annotation_audit/figures/s18_fig3_family_vs_control', 'The family against its sister family'),
    ('exon_tracks', 'main', 'annotation_bugs/figures/exon_tracks', 'Two annotation failures at the exon'),
    ('junctions', 'main', 'expression/figures/s12_junctions', 'The lost genes are spliced'),
    ('contribution', 'main', 'methods/figures/fig_s19_contribution', 'What a database-only search would have missed'),
    ('protein_side', 'main', 'annotation_audit/figures/s18_fig4_protein_side', 'The protein records'),
    ('calibration', 'ed', 'annotation_audit/figures/s18_fig2_calibration', 'The completeness bar'),
    ('cases', 'ed', ['annotation_bugs/figures/annotation_loss', 'annotation_bugs/figures/fragment_tiling'], 'How the cases were chosen, and what the annotation delivers at them'),
    ('case_checks', 'ed', 'annotation_bugs/figures/validation', 'Three checks on each case'),
    ('expression', 'ed', ['expression/figures/s12_detection', 'expression/figures/s12_instruments'], 'Detection against decoys, and reads against deposits'),
    ('gap_coverage', 'ed', 'expression/figures/s12_gap_coverage', 'Read coverage outside the annotation'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('the ryanodine receptors audited in the same assemblies', 'results/annotation_audit/family_vs_control.tsv'),
    ('the completeness bar validated', 'results/annotation_audit/complete_calibration.tsv'),
    ('the per-genome annotation-depth control', 'results/annotation_bugs/annotation_depth.tsv'),
    ('exon-boundary concordance with an independent pipeline', 'results/annotation_bugs/boundary_concordance.tsv'),
    ('the reversed-decoy mapping floor', 'results/expression/expression_by_locus.tsv'),
    ('the cross-mapping control between paralogues', 'results/expression/crossmap_control.tsv'),
    ('the head-to-head restricted to one database', 'results/methods/head_to_head.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = [
    "C01", "AR33", "AR34", "AR27", "AR28", "AR30", "AR31", "C157", "AR32",
    "C04", "AR03", "AR01",
]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = """
If no other paper in the series is published, this one still claims that
the public record of a well-studied vertebrate gene family is far worse than
the family: 744 of the 923 IP3 receptor genes a genome sweep demonstrates
in 309 vertebrate assemblies are reachable by no protein-database search.
The record fails at the rate a size-matched sister family fails in the same
assemblies, so the failure belongs to the archives. It is decided by
whether a curated RefSeq or a submitter GenBank annotation delivered the
gene, by assembly contiguity, and by missing gene symbols, not by the gene.
Two failures are proved at the exon and by RNA-seq junction reads, and 297
coordinate-level corrections are released. The genome sweep is used here
only as a source of loci; its own miss rate, reported by the retention
paper, can only make these failure rates conservative.
"""
STANDALONE_CLAIMS: list[str] = [
    "AR01", "AR02", "AR36", "AR50", "AR52", "AR57", "AR59", "C160", "C161",
    "C163", "C166", "C173",
]
