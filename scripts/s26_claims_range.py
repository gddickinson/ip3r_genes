"""S26 ledger: the load-bearing numbers of the range paper.

`CARRY` names claims already declared by the manuscript's ledger (`C…`), the
thesis's (`T…`) or another paper, by identifier; `CLAIMS` declares numbers no
earlier ledger has, with the `RA` prefix `s26_claims.PREFIX` gives this
paper. Every value must appear in the paper's own text.

Three manuscript claims about the same cells are *not* carried: C137 (the
AlphaFold coverage, which the manuscript ledger holds as the fraction 0.209
while the text states a percentage) and C141-C143 (pLDDT medians held at two
decimals while the text states one). RA67-RA70 check the same cells at the
precision the paper actually prints.
"""

S20 = "results/s20_sweep/report.md"
S23 = "results/s23_scope/report.md"
S11 = "results/structures/report.md"
PFV = "results/s20_sweep/plant_fungal_verdicts.tsv"

CARRY: list[str] = [
    # scope
    "C02", "C03", "C04", "C05", "C06", "T13", "T14",
    # the census
    "C07", "C08", "C09", "C10",
    # absences
    "C11", "C12", "C13", "C14", "C15", "C16", "C17",
    # copy number
    "C18", "C19", "C20",
    # the plant and fungal chase
    "C21", "C22", "C23", "C24",
    # structures
    "C138", "C139", "C140",
    # the enumeration and the architecture rule
    "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22", "T23", "T24",
    # the profiles and census v3
    "T25", "T26", "T27", "T29", "T30", "T31", "T32", "T33", "T34", "T35",
    "T36", "T37", "T38",
]


def _c(i, claim, source, op, expect, **kw):
    return dict(id=f"RA{i:02d}", claim=claim, source=source, op=op,
                expect=str(expect), **kw)


BC = "results/benchmark_controls/summary.json"
CAL = "results/hmm_sweep/calibration_summary.json"
SEED = "results/methods/seed_effect.tsv"
DRIFT = "results/methods/drift_outcome.tsv"
KILL = "results/methods/kill_criterion_validation.tsv"
MAN = "results/s23_scope/manifest_stats.json"
CHOICE = "results/s23_baits/control_profile_choice.tsv"
BAITS = "results/s23_baits/bait_build_stats.json"
LSTAT = "results/s23_scope/ledger_stats.json"
ABS = "results/s23_scope/absence_at_genome.tsv"
SPAN = "results/s23_baits/span_calibration.json"
COPY = "results/s23_scope/copy_number_ledger.tsv"
PLDDT = "results/structures/plddt_summary.tsv"

CLAIMS: list[dict] = [
    # ------------------------------------------------- the benchmark
    _c(1, "recall of the family call is 24 of 25", BC, "json", 24,
       key="recall_hit"),
    _c(2, "25 positive controls", BC, "json", 25, key="recall_n"),
    _c(3, "specificity is 31 of 31", BC, "json", 31,
       key="specificity_hit"),
    _c(4, "all 6 ryanodine receptor decoys are rejected", BC, "json", 6,
       key="ryr_pass"),
    _c(5, "Drosophila Itpr fell 0.008 short of the breadth identity",
       "results/benchmark_controls/report.md", "grep", "short by 0.008"),
    # ------------------------------------------------- the architecture audit
    _c(6, "the architecture rule agrees with 2,479 of 2,479 IP3R symbols",
       "results/census_v2/rule_audit.tsv", "cell", 2479,
       where={"test": "architecture call vs symbol ITPR"}, column="agree"),
    _c(7, "and with 2,492 of 2,492 RyR symbols",
       "results/census_v2/rule_audit.tsv", "cell", 2492,
       where={"test": "architecture call vs symbol RYR"}, column="agree"),
    # ------------------------------------------------- profile calibration
    _c(8, "with seeds removed the profiles agree on 11,875 records", CAL,
       "json", 11875, key="seeds_excluded.agree"),
    _c(9, "and disagree on 1", CAL, "json", 1,
       key="seeds_excluded.disagree"),
    _c(10, "the profiles call 2,074 architecture-unassigned records ITPR",
       CAL, "json", 2074, key="confusion_s2_vs_profile.unassigned → ITPR"),
    _c(11, "and 240 RyR", CAL, "json", 240,
       key="confusion_s2_vs_profile.unassigned → RYR"),
    # ------------------------------------------------- iteration and drift
    _c(12, "three seeds intersect on 4,785 family records", SEED, "cell",
       4785, where={"run": "ALL THREE (intersection)"},
       column="of_which_family"),
    _c(13, "and differ by at most 162", SEED, "cell", 162,
       where={"run": "itpr_acanthamoeba"},
       column="family_unique_to_this_seed"),
    _c(14, "the least-drifted drifted run sits at an off-family share of "
       "0.81", DRIFT, "cell", "0.81", where={"run": "itpr_acanthamoeba"},
       column="offfamily_frac_final", tol=0.005),
    _c(15, "the most off-family clean run sits at 0.34", DRIFT, "cell",
       "0.34", where={"run": "itpr1_human"}, column="offfamily_frac_final",
       tol=0.005),
    _c(16, "the sister-share kill rule fires on none of the drifted runs",
       KILL, "cell", 0, where={"rule__startswith": "K1"},
       column="true_positive"),
    _c(17, "the round ceiling catches all 3 drifted runs", KILL, "cell", 3,
       where={"rule__startswith": "K3"}, column="true_positive"),
    _c(18, "the off-family rule catches 3 of 3 drifted runs", KILL, "cell",
       3, where={"rule__startswith": "proposed"}, column="true_positive"),
    _c(19, "and clears 4 of 4 clean ones", KILL, "cell", 4,
       where={"rule__startswith": "proposed"}, column="true_negative"),
    # ------------------------------------------------- the eukaryote sweep
    _c(20, "28,137 targets were scored by at least one profile", S20,
       "grep", "28,137 targets were scored by at least one profile"),
    _c(21, "2,769 by both above the floor, 1 inside the no-call band", S20,
       "grep", "of those **1/2,769 (0%) fall inside"),
    _c(22, "45 of the 135 clades swept carry a call", S20, "grep",
       "**45 clades of the 135 swept carry an ITPR call**"),
    _c(23, "292 of 311 arthropod proteomes", S20, "grep",
       "| Arthropoda | 311 | 292/311"),
    _c(24, "106 of 110 nematode proteomes", S20, "grep",
       "| Nematoda | 110 | 106/110"),
    _c(25, "15 of 48 Chlorophyta proteomes", S20, "grep",
       "| Chlorophyta | 48 | 15/48"),
    _c(26, "18 of 34 Mucoromycota proteomes", S20, "grep",
       "| Mucoromycota | 34 | 18/34"),
    _c(27, "6 of 16 Chytridiomycota proteomes", S20, "grep",
       "| Chytridiomycota | 16 | 6/16"),
    _c(28, "0 of 634 archaeal proteomes", S20, "grep",
       "| Archaea | 634 | 0/634"),
    _c(29, "3,537 bacterial proteomes, one per genus", S20, "grep",
       "| bacteria_genus | 3,537 |"),
    _c(30, "a non-seed full-length RyR in Salpingoeca rosetta", S20, "grep",
       "| F2UC37 | *Salpingoeca rosetta* | 5340 | 1387.2 | 59% | 65% | no |"),
    _c(31, "architecture-level RyR in 2 non-vertebrate proteomes", S20,
       "grep", "The 2 architecture-level RYR call(s)"),
    # ------------------------------------------------- the controls
    _c(32, "MIR returns 633 matches in land plants and 4,376 in Dikarya",
       S20, "grep", "returns 633 substantial matches in land plants, 4,376 "
       "substantial matches in Dikarya"),
    _c(33, "PF08709: none in land plants against 26 in Chlorophyta", S20,
       "grep", "0 substantial PF08709 match(es), against 26 in"),
    _c(34, "PF08709: none in Dikarya against 16 in Mucoromycota", S20,
       "grep", "0 substantial PF08709 match(es), against 16 in"),
    _c(35, "iteration added 9 targets the single pass never reported", S20,
       "grep", "9 target(s) entered an accepted model that the single "
       "profile pass never reported"),
    _c(36, "the two in an empty lineage are mannosyltransferases", S20,
       "grep", "**2 of them sit in a lineage this task calls empty**"),
    # ------------------------------------------------- the chase
    _c(37, "64 plant records were chased", PFV, "count", 64,
       where={"kingdom": "Viridiplantae"}),
    _c(38, "22 plant records are real genes", PFV, "count", 22,
       where={"kingdom": "Viridiplantae", "verdict": "real_gene"}),
    _c(39, "25 fungal records are real genes", PFV, "count", 25,
       where={"kingdom": "Fungi", "verdict": "real_gene"}),
    _c(40, "the lowest plant cross-kingdom identity is 19.9 %", PFV, "cell",
       "19.9", where={"accession": "A0AAE0BCU2"},
       column="out_kingdom_pident"),
    _c(41, "the highest surviving plant identity is 39.9 %", PFV, "cell",
       "39.9", where={"accession": "A0AAE0G204"},
       column="out_kingdom_pident"),
    _c(42, "the highest surviving fungal identity is 33.5 %", PFV, "cell",
       "33.5", where={"accession": "A0A507FK04"},
       column="out_kingdom_pident"),
    _c(43, "no record of the 99 exceeds 45.8 %", PFV, "cell", "45.8",
       where={"accession": "A0AAE0LLW0"}, column="out_kingdom_pident"),
    # ------------------------------------------------- the genome scope
    _c(44, "the genome scope is 100.4 Gbp", MAN, "json", "100.4",
       key="total_gbp"),
    _c(45, "drawn from 24,596 eukaryotic reference assemblies", MAN, "json",
       24596, key="assemblies_in_pool"),
    _c(46, "less 6,216 vertebrate ones", MAN, "json", 6216,
       key="vertebrate_excluded"),
    _c(47, "three absences no summary had named", S23, "grep",
       "Bacillariophyta 0/16, Rhodophyta 0/12, and **Cestoda 0/11"),
    _c(48, "myosin head in 36 of 36 Aconoidasida proteomes", CHOICE, "cell",
       36, where={"clade": "Aconoidasida"}, column="with_hit"),
    _c(49, "and 23 of 23 Conoidasida", CHOICE, "cell", 23,
       where={"clade": "Conoidasida"}, column="with_hit"),
    _c(50, "SMC N-terminal domain in 12 of 12 red algal proteomes", CHOICE,
       "cell", 12, where={"clade": "Rhodophyta"}, column="with_hit"),
    _c(51, "the panel carries 108 baits", BAITS, "json", 108, key="baits"),
    _c(52, "37 of them IP3 receptors", BAITS, "json", 37, key="itpr"),
    _c(53, "16 ryanodine receptors", BAITS, "json", 16, key="ryr"),
    _c(54, "55 clade controls", BAITS, "json", 55, key="mir_controls"),
    _c(55, "all 56 control clades are filled", BAITS, "json", 56,
       key="control_clades"),
    _c(56, "115 genomes are controlled across a kingdom boundary", LSTAT,
       "json", 115, key="control_verdicts.controlled_cross_kingdom"),
    _c(57, "77 are controlled by the receptor itself", LSTAT, "json", 77,
       key="control_verdicts.controlled_by_target"),
    _c(58, "Ascomycota: none in 31 controlled assemblies", ABS, "cell", 31,
       where={"clade": "Ascomycota"}, column="genomes_controlled"),
    _c(59, "Streptophyta: none in 25", ABS, "cell", 25,
       where={"clade": "Streptophyta"}, column="genomes_controlled"),
    _c(60, "Basidiomycota: none in 17", ABS, "cell", 17,
       where={"clade": "Basidiomycota"}, column="genomes_controlled"),
    _c(61, "two land-plant assemblies carry only tblastn traces", ABS,
       "cell", 2, where={"clade": "Streptophyta"},
       column="genomes_with_trace_only"),
    _c(62, "sixteen gene spans were measured", SPAN, "json", 16,
       key="n_spans"),
    _c(63, "the shortest span is 3,739 bp", SPAN, "json", 3739,
       key="min_span_bp"),
    _c(64, "the longest 324,840 bp", SPAN, "json", 324840,
       key="max_span_bp"),
    _c(65, "the median 19,435 bp", SPAN, "json", 19435,
       key="median_span_bp"),
    # ------------------------------------------------- copy number
    _c(66, "Stentor coeruleus carries 13", COPY, "cell", 13,
       where={"organism": "Stentor coeruleus"}, column="n_full"),
    _c(71, "Dysidea avara carries 8", COPY, "cell", 8,
       where={"organism": "Dysidea avara"}, column="n_full"),
    _c(72, "Cymbomonas tetramitiformis carries three gene models", COPY,
       "cell", 3, where={"organism": "Cymbomonas tetramitiformis"},
       column="n_full"),
    # ------------------------------------------------- structures
    _c(67, "the IP3-binding core is the best-modelled domain at 83.9",
       PLDDT, "cell", "83.9", where={"pfam": "PF08709"},
       column="median_plddt", tol=0.05),
    _c(68, "the pore the worst of the named domains at 71.0", PLDDT, "cell",
       "71.0", where={"pfam": "PF00520"}, column="median_plddt",
       tol=0.05),
    _c(69, "against 69.5 outside annotated domains", PLDDT, "cell", "69.5",
       where={"pfam": "-"}, column="median_plddt", tol=0.05),
    _c(70, "AlphaFold DB models 20.9 % of the census's protein records",
       S11, "grep", "**20.9 % of the census's UniProt-shaped ITPR records "
       "have a usable predicted model**"),
    _c(73, "none of 81 control pairs reaches the same-fold bar", S11,
       "grep", "**0 of 81** control pairs reach the 0.5 same-fold bar"),
    _c(74, "same paralogue, different state: median 0.78; different IP3Rs "
       "0.43", S11, "grep", "score a median **0.78**; two different IP3 "
       "receptors score **0.43**"),
    _c(75, "2 of 7 deep models reach the same-fold bar", S11, "grep",
       "**2 of 7** deep models reach the 0.5 same-fold bar"),
]
