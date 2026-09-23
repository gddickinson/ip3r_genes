"""S26 ledger: the load-bearing numbers of the retention paper.

`CARRY` names claims already declared by the manuscript's ledger (`C…`), the
thesis's (`T…`) or another paper, by identifier; `CLAIMS` declares numbers no
earlier ledger has, with the `RE` prefix `s26_claims.PREFIX` gives this
paper. Every value must appear in the paper's own text.
"""

CARRY: list[str] = [
    # scope and the sweep
    "C01", "C25", "C26", "C27", "C28", "C29", "C30", "C31",
    # the false-negative rate and the contiguity bar
    "C35", "C36", "C37", "C38", "C39", "C41", "C42", "C43", "C44",
    # the panel ablation
    "C45", "C46", "C47", "C48", "C49", "C50",
    # states, parsimony and the sensitivity grid
    "C85", "C86", "C87", "C88", "C89", "C90", "C91", "C92", "C93",
    "C94", "C95", "C96",
    # reading-frame integrity
    "C97", "C98", "C99",
]

_LD = "results/loss_dynamics/loss_dynamics_stats.json"
_LC = "results/loss_counts/loss_counts_stats.json"
_MS = "results/methods/methods_stats.json"
_CT = "results/methods/contiguity_tests.tsv"
_CS = "results/methods/contiguity_by_stratum.tsv"
_AF = "results/methods/absence_floor.tsv"
_RC = "results/loss_dynamics/recon_calibration.tsv"
_SC = "results/loss_dynamics/state_counts.tsv"
_BB = "results/s5_baits/bait_build_stats.json"
_IC = "results/s5_baits/intron_calibration.tsv"
_MC = "results/genome_ledger/margin_calibration.json"
_IT = "results/loss_dynamics/integrity_tests.tsv"
_LB = "results/loss_counts/lesion_by_class.tsv"
_LK = "results/loss_counts/lesion_class_controls.tsv"

CLAIMS: list[dict] = [
    # ------------------------------------------------------------ scope
    dict(id="RE01", claim="161 order representatives enter the scope",
         source="results/genome_manifest.tsv", op="count",
         where={"reasons__startswith": "order_rep"}, expect="161"),
    dict(id="RE02", claim="169 margin species enter the scope",
         source="results/genome_manifest_notes.md", op="grep",
         expect="the union is **169 species**"),
    dict(id="RE03", claim="the bait panel holds 38 baits",
         source=_BB, op="json", key="n_baits", expect="38"),
    dict(id="RE04", claim="30 of them are IP3 receptors",
         source=_BB, op="json", key="n_itpr", expect="30"),
    dict(id="RE05", claim="8 are ryanodine receptors",
         source=_BB, op="json", key="n_ryr", expect="8"),
    dict(id="RE06", claim="the panel totals 121,294 residues",
         source=_BB, op="json", key="total_residues", expect="121294"),
    dict(id="RE07", claim="the widest IP3R intron is 152,216 bp (human ITPR1)",
         source=_IC, op="cell", column="max_intron_any_transcript",
         where={"species": "homo_sapiens", "symbol": "ITPR1"},
         expect="152216"),
    dict(id="RE08", claim="the widest RyR intron is 227,927 bp (human RYR2)",
         source=_IC, op="cell", column="max_intron_any_transcript",
         where={"species": "homo_sapiens", "symbol": "RYR2"},
         expect="227927"),
    dict(id="RE09", claim="the inherited attribution margin was 0.333",
         source=_MC, op="json", key="inherited_threshold", expect="0.333"),
    dict(id="RE10", claim="the margin in force is 0.22",
         source=_MC, op="json", key="threshold_in_force", expect="0.22"),
    dict(id="RE11", claim="542 annotation-established loci were measured",
         source=_MC, op="json", key="n_contested", expect="542"),
    dict(id="RE12", claim="368 of them fall below the inherited margin",
         source=_MC, op="json", key="n_contested_below_inherited",
         expect="368"),
    dict(id="RE13", claim="all 53 annotated rescue fragments agree",
         source=_MC, op="json", key="fragment_n_correct", expect="53"),
    dict(id="RE14", claim="the contiguity bar is the median span of 28 genes",
         source="results/genome_ledger/ledger_stats.json", op="json",
         key="itpr_span_stats.n", expect="28"),
    # ------------------------------------------------ false negatives
    dict(id="RE15", claim="a found cell's median contig N50 is 3,396,515 bp",
         source=_CT, op="grep", expect="median 3,396,515 vs 23,460"),
    dict(id="RE16", claim="the odds rise 8.10-fold per tenfold N50",
         source="results/methods/report.md", op="grep", expect="8.10"),
    dict(id="RE17", claim="and 19.99-fold in the RyR series",
         source=_CT, op="grep", expect="OR per 10x=19.99"),
    dict(id="RE18", claim="chromosome-level assemblies miss 3 IP3R cells",
         source=_CS, op="cell", column="n_false_negative",
         where={"stratum": "assembly_level", "value": "Chromosome",
                "series": "itpr_present"}, expect="3"),
    dict(id="RE19", claim="of 512", source=_CS, op="cell", column="n_cells",
         where={"stratum": "assembly_level", "value": "Chromosome",
                "series": "itpr_present"}, expect="512"),
    dict(id="RE20", claim="and 0 of 172 RyR cells", source=_CS, op="cell",
         column="n_cells", where={"stratum": "assembly_level",
                                  "value": "Chromosome",
                                  "series": "ryr_sister"}, expect="172"),
    dict(id="RE21", claim="below the bar ITPR1 misses 0.3917", source=_CT,
         op="cell", column="effect",
         where={"test": "false_negative_rate_below_bar", "series": "ITPR1"},
         expect="0.3917"),
    dict(id="RE22", claim="ITPR2 misses 0.4333", source=_CT, op="cell",
         column="effect",
         where={"test": "false_negative_rate_below_bar", "series": "ITPR2"},
         expect="0.4333"),
    dict(id="RE23", claim="ITPR3 misses 0.3", source=_CT, op="cell",
         column="effect",
         where={"test": "false_negative_rate_below_bar", "series": "ITPR3"},
         expect="0.3000"),
    dict(id="RE24", claim="above the bar the IP3R series has 563 cells",
         source=_AF, op="cell", column="n_cells",
         where={"contig_n50_floor": "142212", "series": "itpr_present"},
         expect="563"),
    dict(id="RE25", claim="of which 5 are missed", source=_AF, op="cell",
         column="n_false_negative",
         where={"contig_n50_floor": "142212", "series": "itpr_present"},
         expect="5"),
    dict(id="RE26", claim="120 genomes fall below the bar", source=_LD,
         op="json", key="headline.n_genomes_below_bar", expect="120"),
    # ---------------------------------------------------------- panel
    dict(id="RE27", claim="dropping ITPR2's baits costs 240 cells",
         source=_MS, op="json", key="headline.panel.drop_ITPR2_baits.delta",
         expect="-240"),
    dict(id="RE28", claim="19 panels were simulated, the full one and "
         "18 ablations", source="results/methods/panel_recall.tsv",
         op="nunique", column="panel", expect="19"),
    dict(id="RE29", claim="2,179 call changes across all panels",
         source="results/methods/panel_changes.tsv", op="count",
         expect="2179"),
    dict(id="RE30", claim="68 of them are gains",
         source="results/methods/panel_changes.tsv", op="count",
         where={"direction": "gain"}, expect="68"),
    # --------------------------------------------------------- states
    dict(id="RE31", claim="783 cells are a single full-coverage locus",
         source=_SC, op="cell", column="n",
         where={"scope": "all", "state": "present_single_locus"},
         expect="783"),
    dict(id="RE32", claim="90 are truncated by the assembly", source=_SC,
         op="cell", column="n",
         where={"scope": "all", "state": "present_truncated"}, expect="90"),
    dict(id="RE33", claim="7 are partial", source=_SC, op="cell", column="n",
         where={"scope": "all", "state": "present_partial"}, expect="7"),
    dict(id="RE34", claim="11 of 120 genomes below the bar fall short of 2.5 "
         "gene-equivalents", source=_LD, op="json",
         key="headline.n_below_2_5_below_bar", expect="11"),
    # -------------------------------------------------- reconstruction
    dict(id="RE35", claim="44 positive regions calibrate the bar", source=_RC,
         op="cell", column="n_positive", where={"metric": "coverage"},
         expect="44"),
    dict(id="RE36", claim="9 negative regions", source=_RC, op="cell",
         column="n_negative", where={"metric": "coverage"}, expect="9"),
    dict(id="RE37", claim="20 co-shattered regions", source=_RC, op="cell",
         column="n_co_trace", where={"metric": "coverage"}, expect="20"),
    dict(id="RE38", claim="positive median coverage 0.7948", source=_RC,
         op="cell", column="positive_median", where={"metric": "coverage"},
         expect="0.7948"),
    dict(id="RE39", claim="negative median 0.0303", source=_RC, op="cell",
         column="negative_median", where={"metric": "coverage"},
         expect="0.0303"),
    dict(id="RE40", claim="co-shattered median 0.6308", source=_RC, op="cell",
         column="co_trace_median", where={"metric": "coverage"},
         expect="0.6308"),
    dict(id="RE41", claim="the gap's upper edge is 0.1737", source=_RC,
         op="cell", column="gap_hi", where={"metric": "coverage"},
         expect="0.1737"),
    dict(id="RE58", claim="the gap's lower edge is 0.1", source=_RC,
         op="cell", column="gap_lo", where={"metric": "coverage"},
         expect="0.1"),
    dict(id="RE42", claim="the operating point is 0.13685", source=_RC,
         op="cell", column="operating_point", where={"metric": "coverage"},
         expect="0.13685"),
    dict(id="RE43", claim="432 rescue regions were offered to synteny",
         source=_LD, op="json", key="headline.synteny_regions",
         expect="432"),
    dict(id="RE44", claim="273 sit on a contig with no annotated gene",
         source=_LD, op="json", key="headline.synteny_no_neighbourhood",
         expect="273"),
    dict(id="RE45", claim="8 reach the caller's floor", source=_LD, op="json",
         key="headline.synteny_reached", expect="8"),
    dict(id="RE46", claim="the median region's contig extends 39.936 kb",
         source="results/loss_dynamics/synteny_reach_summary.tsv", op="cell",
         column="median_contig_extent_kb", where={"window": "fixed10"},
         expect="39.936"),
    # ---------------------------------------------- tree, parsimony, Mk
    dict(id="RE47", claim="the taxonomy tree has 468 internal nodes",
         source="results/loss_dynamics/tree_polytomies.tsv", op="count",
         expect="468"),
    dict(id="RE48", claim="49 of them are polytomies",
         source="results/loss_dynamics/tree_polytomies.tsv", op="count",
         where={"is_polytomy": "1"}, expect="49"),
    dict(id="RE49", claim="the family coding's worst case is 1 genome",
         source=_LC, op="json", key="headline.max_family_losses",
         expect="1"),
    dict(id="RE50", claim="twelve Mk combinations are refused at the base "
         "setting", source=_LC, op="json", key="headline.n_mk_rows_base",
         expect="12"),
    # ------------------------------------------------------- integrity
    dict(id="RE51", claim="ITPR3 exceeds its matched siblings in 39 genomes",
         source=_IT, op="cell", column="n_pos",
         where={"cell": "ITPR3", "matched": "1"}, expect="39"),
    dict(id="RE52", claim="to 14", source=_IT, op="cell", column="n_neg",
         where={"cell": "ITPR3", "matched": "1"}, expect="14"),
    dict(id="RE53", claim="7 family loci fire a fossil reading", source=_LC,
         op="json", key="headline.n_fossils_itpr", expect="7"),
    dict(id="RE54", claim="the bird ITPR3 excess is 25 genomes", source=_LB,
         op="cell", column="n_pos",
         where={"cell": "ITPR3", "vclass": "Aves"}, expect="25"),
    dict(id="RE55", claim="to 2", source=_LB, op="cell", column="n_neg",
         where={"cell": "ITPR3", "vclass": "Aves"}, expect="2"),
    dict(id="RE56", claim="21 of the 27 bird pairs lie below the bar",
         source=_LK, op="cell", column="n_below_bar",
         where={"cell": "ITPR3", "vclass": "Aves"}, expect="21"),
    dict(id="RE57", claim="6 pairs lie above it", source=_LK, op="cell",
         column="n_above_bar", where={"cell": "ITPR3", "vclass": "Aves"},
         expect="6"),
]
