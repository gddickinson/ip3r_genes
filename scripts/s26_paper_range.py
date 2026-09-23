"""S26 paper configuration: The family's range across the eukaryotes.

The paper's title, the one question it answers (P1), the figures it places
(P4, P5), the controls it measures itself (P2), the claims that declare its
scope (P3), and what it claims if no other paper in the series is ever
published (P6). Enforced by `s26_rules.py` and `s26_figures.py`; the prose is
in `papers/range/`.
"""

TITLE = 'The IP₃ receptor is ancestrally eukaryotic and has been lost repeatedly outside the animals'
SHORT_TITLE = "The family's range across the eukaryotes"
QUESTION = 'In which eukaryotic lineages does an IP₃ receptor gene exist?'

#: (slug, kind, source stem(s) under results/, caption stub). Numbered by
#: first mention in the body, never here.
FIGURES = [
    ('range', 'main', 's20_sweep/figures/range_by_phylum', 'The family across the eukaryotes, by proteome'),
    ('separation', 'main', 'census_v3/figures/profile_separation', 'Two profiles separate the IP₃ receptors from the ryanodine receptors'),
    ('absence', 'main', 's23_scope/figures/absence_at_genome', 'Absences confirmed in controlled genomes'),
    ('chase', 'main', 's20_sweep/figures/plant_fungal_chase', 'Every plant and fungal record chased to a verdict'),
    ('copies', 'main', 's23_scope/figures/copy_number', 'Copy number outside the vertebrates'),
    ('fold', 'main', 'structures/figures/s11_tm_calibration', 'The fold separates the two families'),
    ('enumeration', 'ed', ['census_v2/figures/census_space', 'census_v2/figures/census_lineage'], "The uncapped enumeration of the family's signatures"),
    ('margin', 'ed', 'census_v2/figures/census_margin', 'The labelled-bait margin across the census'),
    ('census_growth', 'ed', ['census_v2/figures/census_growth', 'census_v2/figures/census_lengths'], 'How the census grew, and the length of what it holds'),
    ('instruments', 'ed', ['census_v3/figures/instrument_agreement', 'census_v3/figures/proteome_copy_number'], 'Two instruments on one census'),
    ('completeness', 'ed', ['census_v3/figures/jackhmmer_convergence', 's20_sweep/figures/jackhmmer_s20', 'methods/figures/fig_s19_drift'], 'Iterative search to convergence, and the drift guard'),
    ('separation_outside', 'ed', 's20_sweep/figures/profile_separation', 'Profile separation outside the vertebrates'),
    ('genome_instrument', 'ed', ['s23_scope/figures/identity_floor', 's23_scope/figures/span_inflation'], "The genomic sweep's locus floor and span"),
    ('structures', 'ed', ['structures/figures/s11_panel', 'structures/figures/s11_plddt_domains'], 'The structure panel and model confidence'),
    ('afdb', 'ed', 'structures/figures/s11_afdb_coverage', 'AlphaFold DB coverage of the census'),
]

#: (what it controls, the committed table it is measured in). P2: each table
#: must be primary in this paper.
CONTROLS = [
    ('ryanodine receptor decoys: specificity of the family call', 'results/benchmark_controls/negative_controls.tsv'),
    ('labelled positives: recall of the family call', 'results/benchmark_controls/positive_controls.tsv'),
    ('the architecture rule scored against gene symbols it never sees', 'results/census_v2/rule_audit.tsv'),
    ('profile assignment calibrated against the architecture call', 'results/hmm_sweep/calibration_summary.json'),
    ('the per-clade search control that makes a genome absence a measurement', 'results/s23_scope/control_ledger.tsv'),
    ('the structural family call against both reference families', 'results/structures/tm_vs_reference.tsv'),
    ('the jackhmmer drift criterion scored as a classifier', 'results/methods/kill_criterion_validation.tsv'),
]

#: Ledger claims that declare the paper's denominators (P3).
SCOPE_CLAIMS: list[str] = [
    "C04",   # 763 vertebrate reference proteomes
    "C05",   # 14,414,821 vertebrate proteins
    "C03",   # 6,928 non-vertebrate reference proteomes
    "C06",   # 63,144,898 non-vertebrate proteins
    "C02",   # 194 non-vertebrate genome assemblies
    "RA44",  # 100.4 Gbp
    "T13",   # the enumerated space, 15,421 proteins
]

#: P6: what this paper claims if none of the others is ever published, and
#: the ledger claims that answer rests on (each primary in this paper).
STANDALONE = """
Read alone, this paper claims that the IP3 receptor is an ancestral
eukaryotic gene that has been lost independently in land plants, in Dikarya
and in several other fungal, protist and parasitic lineages, and that it is
absent from every prokaryote. Each part rests on evidence measured inside the
paper: a family call benchmarked on ryanodine receptor decoys and audited
against gene symbols it never reads; a sweep of 7,691 reference proteomes
with every negative made at two sensitivities and beside a within-kingdom
presence; all 99 plant and fungal records chased to a verdict; and 35
clade-level absences confirmed in genome assemblies that each carry a
positive control chosen for their clade by measurement. It needs no other
paper in the series, cites none, and is the one the others take their
family call from.
"""
STANDALONE_CLAIMS: list[str] = [
    "RA03", "RA04",          # the call rejects every RyR decoy
    "RA06", "RA07",          # the architecture rule audited on symbols
    "C10", "C11", "C12", "C13",  # presence, and the plant/Dikarya absences
    "RA28",                  # no archaeal proteome
    "RA32", "RA33",          # each negative beside a positive control
    "C21", "C22", "C23", "C24",  # the chase: real genes, no contamination
    "C15", "C16", "C17",     # 35 absences hold, no genome uncontrolled
]
