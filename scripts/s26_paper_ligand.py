"""S26 paper configuration: The ligand site.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/ligand/`.
"""

TITLE = ("The IP₃ receptor's binding core is not more constrained than its "
         "pore, and the answer depends on where the pore ends")
SHORT_TITLE = 'The ligand site'
QUESTION = 'Is the IP₃-binding core under different constraint from the pore it gates?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('modules', 'main', 'ligand_site/figures/s22_fig1_modules', 'The ligand core against the pore'),
    ('shells', 'main', 'ligand_site/figures/s22_fig2_shells', 'Constraint against distance to the ligand'),
    ('omega', 'main', 'ligand_site/figures/s22_fig3_omega', 'The same contrast on the rate axis'),
    ('lineage', 'main', 'ligand_site/figures/s22_fig4_lineage', 'Lineages without the upstream enzyme'),
    ('labelled_positions', 'supp', 'supplementary/figures/SuppFig2_labelled_positions', 'The core and the pore at residue resolution'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('the ten measured contacts recovered in six depositions', 'results/ligand_site/shell_agreement.tsv'),
    ("the ryanodine receptors as the lineage test's positive control", 'results/ligand_site/deep_lineage_ryr_control.tsv'),
    ('the contrast under all four module definitions', 'results/ligand_site/module_contrast.tsv'),
    ("the lineage test's power", 'results/ligand_site/lineage_power.tsv'),
    ("the vertebrate proteomes as the enzyme search's positive control", 'results/ligand_site/plc_repertoire.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = ["C239", "C240", "C241", "C222", "C246",
                            "LI29"]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = (
    "Read alone, this paper claims three things about the one module of the "
    "IP3 receptor that its sister family does not share functionally. "
    "First, paired inside each of roughly 250 vertebrate orthologues per "
    "paralogue, the IP3-binding core is less conserved than the pore module "
    "in ITPR1 and ITPR3 and level in ITPR2, and this ordering reverses in "
    "all three paralogues when the luminal loop is counted as pore, so the "
    "comparison is meaningful only with that boundary declared. Second, "
    "measured as all-atom distance in six IP3-bound structures, the ten "
    "contact residues are more constrained than the rest of the core but "
    "not than the rest of a pocket about 15 angstrom across. Third, 64 "
    "eukaryotic proteomes carry the receptor without the phospholipase C "
    "that makes IP3, and matched on divergence their ligand core is not "
    "relaxed, against a ryanodine receptor control the test detects. Every "
    "one of these rests on tables this paper computes; it takes only the "
    "per-residue conservation layers and orthologue alignments of the "
    "constraint paper as inputs, and it re-declares their scope in its own "
    "methods."
)
STANDALONE_CLAIMS: list[str] = ["C233", "C236", "C242", "C243", "C244",
                                 "C225", "LI18", "C247", "C252"]
