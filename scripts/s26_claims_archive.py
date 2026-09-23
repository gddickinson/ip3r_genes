"""S26 ledger: the load-bearing numbers of the archive paper.

`CARRY` names claims already declared by the manuscript's ledger (`C…`), the
thesis's (`T…`) or another paper, by identifier; `CLAIMS` declares numbers no
earlier ledger has, with the `AR` prefix `s26_claims.PREFIX` gives this
paper. Every value must appear in the paper's own text.

The paper's headline, that four in five demonstrated IP₃ receptor genes are
unreachable, is declared here from `results/methods/gene_recovery.tsv`
rather than carried from the manuscript's C168/C169. Those read
`methods_stats.json`, which P5s assigns to the retention paper, and they
count the ryanodine receptor control cell in with the family (940 of
1,232). This paper states the family (744 of 923) and the control (196 of
309) separately.
"""

A = "results/annotation_audit/"
B = "results/annotation_bugs/"
E = "results/expression/"
M = "results/methods/"
ITPR = "ITPR1,ITPR2,ITPR3"
ALL = "all scorable loci"
D4 = "loci on a contig that spans the gene (D4)"

CARRY: list[str] = [
    "C01", "C04", "C43",                      # scope: 309 genomes, 763 proteomes, bar
    "C32", "C33", "C34",                      # census v4 db-status
    "C148", "C149", "C150",                   # non-coding-only excess
    "C156", "C157", "C158", "C159", "C160", "C161", "C162",
    "C163", "C164", "C165",                   # corrections
    "C166", "C167",                           # S10 eligibility
    "C173", "C174", "C175",                   # S12
]


def _c(i, claim, source, op, expect, **kw):
    return dict(id=f"AR{i:02d}", claim=claim, source=source, op=op,
                expect=str(expect), **kw)


def _cell(i, claim, source, column, where, expect, **kw):
    return _c(i, claim, source, "cell", expect, column=column, where=where,
              **kw)


GR = M + "gene_recovery.tsv"
BYCELL = M + "gene_recovery_by_cell.tsv"
SCOPE = M + "gene_recovery_by_scope.tsv"
FVC = A + "family_vs_control.tsv"
SRC = A + "by_source.tsv"
STATS = A + "annotation_audit_stats.json"
V4 = "results/census_v4/census_v4_stats.json"
RANK = B + "case_ranking.tsv"
BLOCKS = B + "block_accounting.tsv"
RF = B + "reading_frame.tsv"
PROBE = B + "probe_summary.tsv"
AUDIT = B + "assembly_audit.tsv"
CORR = A + "correction_summary.tsv"

CLAIMS: list[dict] = [
    # ---------------------------------------------------------- reachability
    _c(1, "923 IP3 receptor genes are demonstrated in the DNA", GR, "count",
       923, where={"gene_present": "1", "cell__in": ITPR}),
    _c(2, "744 of them have no protein record that reaches them", GR, "count",
       744, where={"cell__in": ITPR,
                   "recovery_channel__startswith": "genome_only"}),
    _c(3, "1,232 genome-by-cell genes are demonstrated, both families", GR,
       "count", 1232, where={"gene_present": "1"}),
    _cell(4, "257 ITPR1 genes are unreachable", BYCELL,
          "n_protein_db_invisible", {"cell": "ITPR1"}, 257),
    _cell(5, "of 309 ITPR1 genes", BYCELL, "n_genes_present",
          {"cell": "ITPR1"}, 309),
    _cell(6, "260 ITPR2 genes are unreachable", BYCELL,
          "n_protein_db_invisible", {"cell": "ITPR2"}, 260),
    _cell(7, "of 307 ITPR2 genes", BYCELL, "n_genes_present",
          {"cell": "ITPR2"}, 307),
    _cell(8, "227 ITPR3 genes are unreachable", BYCELL,
          "n_protein_db_invisible", {"cell": "ITPR3"}, 227),
    _cell(9, "of 307 ITPR3 genes", BYCELL, "n_genes_present",
          {"cell": "ITPR3"}, 307),
    _cell(10, "196 ryanodine receptor control cells are unreachable", BYCELL,
          "n_protein_db_invisible", {"cell": "RYR"}, 196),
    _cell(11, "of 309 control cells", BYCELL, "n_genes_present",
          {"cell": "RYR"}, 309),
    _c(12, "289 unreachable IP3R genes are in species without a reference "
       "proteome", GR, "count", 289,
       where={"cell__in": ITPR,
              "recovery_channel": "genome_only:no_reference_proteome"}),
    _c(13, "186 in species with only fragmentary records", GR, "count", 186,
       where={"cell__in": ITPR, "recovery_channel":
              "genome_only:records_exist_but_none_full_length"}),
    _c(14, "254 where records exist and none resolves to the paralogue", GR,
       "count", 254, where={"cell__in": ITPR, "recovery_channel":
                            "genome_only:no_record_resolves_to_this_paralog"}),
    _c(15, "15 in species with no family record at all", GR, "count", 15,
       where={"cell__in": ITPR, "recovery_channel":
              "genome_only:no_family_record_for_species"}),
    _c(16, "179 IP3R genes are held as DNA and as a resolving record", GR,
       "count", 179, where={"cell__in": ITPR, "recovery_channel":
                            "protein_database_and_genome"}),
    _cell(17, "order representatives lose 413 genes", SCOPE,
          "n_protein_db_invisible", {"scope": "order representative"}, 413),
    _cell(18, "of 556", SCOPE, "n_genes_present",
          {"scope": "order representative"}, 556),
    _cell(19, "margin species lose 465 genes", SCOPE,
          "n_protein_db_invisible", {"scope": "margin species"}, 465),
    _cell(20, "of 592", SCOPE, "n_genes_present",
          {"scope": "margin species"}, 592),
    _cell(21, "Pfam enumeration returns 3,135 gene-scale vertebrate records",
          M + "head_to_head.tsv", "n_pfam",
          {"database": "vertebrata", "length_band": ">=2000 aa (gene-scale)"},
          3135),
    _cell(22, "the profile pair returns 3,136", M + "head_to_head.tsv",
          "n_hmm", {"database": "vertebrata",
                    "length_band": ">=2000 aa (gene-scale)"}, 3136),
    _cell(23, "3,135 are shared", M + "head_to_head.tsv", "n_shared",
          {"database": "vertebrata", "length_band": ">=2000 aa (gene-scale)"},
          3135),
    _c(24, "the sweep's models come from 224 genomes", V4, "json", 224,
       key="genomes_contributing"),
    _c(25, "488 of the models are IP3 receptors", V4, "json", 488,
       key="itpr_models"),
    _c(26, "3 lie inside a gene named for another paralogue", V4, "json", 3,
       key="itpr_by_db_status.annotated_other_paralog"),
    # ---------------------------------------------------------- scope
    _c(27, "2,144 gene-scale loci are audited", STATS, "json", 2144,
       key="headline.n_loci"),
    _c(28, "1,059 of them are IP3 receptor loci", STATS, "json", 1059,
       key="headline.n_itpr_loci"),
    _c(29, "255 loci sit in assemblies with no gene set", STATS, "json", 255,
       key="headline.n_loci_no_gene_set"),
    _c(30, "932 IP3 receptor loci are scorable", STATS, "json", 932,
       key="headline.n_scorable_itpr_loci"),
    _c(31, "942 ryanodine receptor loci are scorable", STATS, "json", 942,
       key="headline.n_control_loci"),
    _c(32, "8,306 of the protein records are vertebrate", STATS, "json",
       8306, key="headline.n_protein_vertebrate"),
    _c(33, "168 assemblies carry a RefSeq annotation",
       "results/genome_manifest.tsv", "count", 168,
       where={"annotation_source": "RefSeq"}),
    _c(34, "141 carry a GenBank annotation or none",
       "results/genome_manifest.tsv", "count", 141,
       where={"annotation_source": "GenBank"}),
    _c(35, "the bait panel holds 38 proteins",
       "results/s5_baits/bait_manifest.tsv", "count", 38),
    # ---------------------------------------------------------- locus states
    _cell(36, "689 IP3R loci are complete", FVC, "n_itpr_failing",
          {"scope": ALL, "measure": "complete"}, 689),
    _cell(37, "243 are not", FVC, "n_itpr_failing",
          {"scope": ALL, "measure": "any annotation failure"}, 243),
    _cell(38, "149 are fragmentary", FVC, "n_itpr_failing",
          {"scope": ALL, "measure": "fragmentary"}, 149),
    _cell(39, "14 are split", FVC, "n_itpr_failing",
          {"scope": ALL, "measure": "split"}, 14),
    _cell(40, "38 are unannotated", FVC, "n_itpr_failing",
          {"scope": ALL, "measure": "unannotated"}, 38),
    _cell(41, "the bar is validated on 1,077 loci",
          A + "complete_calibration.tsv", "n", {"usable": "1"}, 1077),
    _cell(42, "14 of them fall below the bar",
          A + "complete_calibration.tsv", "n_below_bar", {"usable": "1"}, 14),
    _cell(43, "712 loci are complete at a bar of 0.3",
          A + "state_sensitivity.tsv", "n_complete",
          {"bar": "0.3", "min_piece": "0.05"}, 712),
    _cell(44, "635 at a bar of 0.95", A + "state_sensitivity.tsv",
          "n_complete", {"bar": "0.95", "min_piece": "0.05"}, 635),
    _cell(45, "on good contigs 38 IP3R loci fail", FVC, "n_itpr_failing",
          {"scope": D4, "measure": "any annotation failure"}, 38),
    _cell(46, "of 568", FVC, "n_itpr", {"scope": D4,
                                         "measure": "any annotation failure"},
          568),
    _cell(47, "243 of 433 bird loci are complete", A + "state_by_class.tsv",
          "n_complete", {"vclass": "Aves"}, 243),
    _cell(48, "of 433 bird loci", A + "state_by_class.tsv", "n_loci",
          {"vclass": "Aves"}, 433),
    _cell(49, "124 bird loci are fragmentary", A + "state_by_class.tsv",
          "n_fragmentary", {"vclass": "Aves"}, 124),
    # ---------------------------------------------------------- the control
    _cell(50, "734 ryanodine receptor loci are complete", FVC,
          "n_ryr_failing", {"scope": ALL, "measure": "complete"}, 734),
    _cell(51, "208 fail", FVC, "n_ryr_failing",
          {"scope": ALL, "measure": "any annotation failure"}, 208),
    _cell(52, "the family-control difference has q = 0.13", FVC, "q_bh",
          {"scope": ALL, "measure": "any annotation failure"}, "0.13",
          tol=0.005),
    _cell(53, "on good contigs 35 control loci fail", FVC, "n_ryr_failing",
          {"scope": D4, "measure": "any annotation failure"}, 35),
    _cell(54, "of 632", FVC, "n_ryr", {"scope": D4,
                                        "measure": "any annotation failure"},
          632),
    _cell(55, "and the difference there has q = 0.82", FVC, "q_bh",
          {"scope": D4, "measure": "any annotation failure"}, "0.82",
          tol=0.005),
    _cell(56, "the non-coding-only excess has q = 0.006", FVC, "q_bh",
          {"scope": ALL, "measure": "noncoding"}, "0.006", tol=0.0005),
    # ---------------------------------------------------------- the archive
    _cell(57, "1,161 RefSeq loci are complete", SRC, "n_refseq",
          {"scope": ALL, "state": "complete"}, 1161),
    _cell(58, "of 1,175", SRC, "n_refseq_total",
          {"scope": ALL, "state": "complete"}, 1175),
    _cell(59, "262 GenBank loci are complete", SRC, "n_genbank",
          {"scope": ALL, "state": "complete"}, 262),
    _cell(60, "of 699", SRC, "n_genbank_total",
          {"scope": ALL, "state": "complete"}, 699),
    _cell(61, "263 GenBank loci are fragmentary", SRC, "n_genbank",
          {"scope": ALL, "state": "fragmentary"}, 263),
    _cell(62, "92 GenBank loci are unannotated", SRC, "n_genbank",
          {"scope": ALL, "state": "unannotated"}, 92),
    _cell(63, "10 RefSeq loci are unannotated", SRC, "n_refseq",
          {"scope": ALL, "state": "unannotated"}, 10),
    _cell(64, "on good contigs 1,009 RefSeq loci are complete", SRC,
          "n_refseq", {"scope": D4, "state": "complete"}, 1009),
    _cell(65, "of 1,015", SRC, "n_refseq_total",
          {"scope": D4, "state": "complete"}, 1015),
    _cell(66, "against 118 GenBank loci", SRC, "n_genbank",
          {"scope": D4, "state": "complete"}, 118),
    _cell(67, "of 185", SRC, "n_genbank_total",
          {"scope": D4, "state": "complete"}, 185),
    # ---------------------------------------------------------- proteins
    _c(68, "5 records carry a name from the sister family", STATS, "json", 5,
       key="headline.n_name_wrong_family"),
    _c(69, "74 name a paralogue the panel cannot call", STATS, "json", 74,
       key="headline.n_paralog_not_callable"),
    _c(70, "66 are named for the superfamily", STATS, "json", 66,
       key="headline.n_name_family_ambiguous"),
    _cell(71, "9,841 records are called family by both instruments",
          A + "pfam_recall.tsv", "n_records",
          {"axis": "which instruments called it", "bucket": "both"}, 9841),
    _cell(72, "1,554 by one instrument only", A + "pfam_recall.tsv",
          "n_records", {"axis": "which instruments called it",
                        "bucket": "v4+s20"}, 1554),
    # ---------------------------------------------------------- cases
    _c(73, "the sweep recovered 880 IP3 receptor loci", RANK, "count", 880),
    _c(74, "the annotation delivers 375 of the eligible loci", RANK, "count",
       375, where={"eligible": "True", "mode": ""}),
    _c(75, "360 eligible loci have zero annotation loss", RANK, "count", 360,
       where={"eligible": "True", "loss": "0.0"}),
    _c(76, "the fifth rule removes 6 loci", RANK, "count", 6,
       where={"failed_rule": "E5"}),
    _cell(77, "case A has 56 coding exons", BLOCKS, "n_aligned_exons",
          {"case_id": "case_a"}, 56),
    _cell(78, "totalling 8,021 bp", BLOCKS, "aligned_cds_bp",
          {"case_id": "case_a"}, 8021),
    _cell(79, "case A's contig N50 is 84.24 times the locus", B + "cases.tsv",
          "headroom", {"case_id": "case_a"}, "84.24"),
    _cell(80, "case A splices into 2,673 codons", RF, "cds_codons",
          {"case_id": "case_a", "cell": "ITPR2"}, 2673),
    _cell(81, "against 14.2 stops expected", RF, "expected_stops",
          {"case_id": "case_a", "cell": "ITPR2"}, "14.2"),
    _c(82, "case A has 54 GT-AG introns", B + "case_introns.tsv", "count",
       54, where={"case_id": "case_a", "splice_class": "canonical"}),
    _c(83, "and 55 introns in all", B + "case_introns.tsv", "count", 55,
       where={"case_id": "case_a"}),
    _cell(84, "52 of case A's boundaries are shared by a majority", PROBE,
          "n_shared_with_majority", {"case_id": "case_a"}, 52),
    _cell(85, "of 35 reference genomes", PROBE, "n_reference_genomes",
          {"case_id": "case_a"}, 35),
    _cell(86, "SSPN flanks the paralogue in 197 swept vertebrates",
          B + "flank_consensus_check.tsv", "consensus_species",
          {"case_id": "case_a", "symbol": "SSPN"}, 197),
    _cell(87, "of 215", B + "flank_consensus_check.tsv",
          "consensus_species_total", {"case_id": "case_a", "symbol": "SSPN"},
          215),
    _cell(88, "case A's annotation files 7,380 genes as pseudogenes", AUDIT,
          "annot_pseudogene", {"case_id": "case_a"}, 7380),
    _cell(89, "of 23,345", AUDIT, "annot_genes_total", {"case_id": "case_a"},
          23345),
    _cell(90, "case B has 60 coding exons", BLOCKS, "n_aligned_exons",
          {"case_id": "case_b"}, 60),
    _cell(91, "totalling 8,064 bp", BLOCKS, "aligned_cds_bp",
          {"case_id": "case_b"}, 8064),
    _cell(92, "covered by 3 gene models", BLOCKS, "n_annotated_gene_models",
          {"case_id": "case_b"}, 3),
    _cell(93, "delivering 4,664 bp", BLOCKS, "annotated_cds_bp_on_gene",
          {"case_id": "case_b"}, 4664),
    _cell(94, "leaving 3,400 bp unmodelled", BLOCKS, "uncovered_bp",
          {"case_id": "case_b"}, 3400),
    _cell(95, "in 28 blocks", BLOCKS, "n_uncovered_blocks",
          {"case_id": "case_b"}, 28),
    _cell(96, "case B gives 2,687 codons", RF, "cds_codons",
          {"case_id": "case_b", "cell": "ITPR3"}, 2687),
    _cell(97, "against 11.4 stops expected", RF, "expected_stops",
          {"case_id": "case_b", "cell": "ITPR3"}, "11.4"),
    _c(98, "all 59 of case B's introns are GT-AG", B + "case_introns.tsv",
       "count", 59, where={"case_id": "case_b",
                           "splice_class": "canonical"}),
    _cell(99, "59 boundaries shared by a majority", PROBE,
          "n_shared_with_majority", {"case_id": "case_b"}, 59),
    _cell(100, "of 40 reference genomes", PROBE, "n_reference_genomes",
          {"case_id": "case_b"}, 40),
]

CLAIMS += [
    _cell(101, "case A's species has 43 transcript records", PROBE,
          "n_total_mrna", {"case_id": "case_a"}, 43),
    _cell(102, "case B's species has 10", PROBE, "n_total_mrna",
          {"case_id": "case_b"}, 10),
    _cell(103, "the genomic control makes 112 hits at case A", PROBE,
          "hits_genomic_control", {"case_id": "case_a"}, 112),
    _cell(104, "and 50 at case B", PROBE, "hits_genomic_control",
          {"case_id": "case_b"}, 50),
    # ---------------------------------------------------------- expression
    _c(105, "536,000,000 reads were streamed", E + "run_metrics.tsv", "sum",
       536000000, column="library_reads"),
    _c(106, "from 16 studies", E + "run_metrics.tsv", "nunique", 16,
       column="study"),
    _c(107, "all 7 lost loci are in the expression panel",
       E + "expression_by_locus.tsv", "count", 7, where={"role": "failed"}),
    _c(108, "the decoys collected 42 reads", E + "expression_by_run.tsv",
       "sum", 42, column="decoy_reads"),
    _c(109, "over 469 run-by-locus comparisons", E + "expression_by_run.tsv",
       "count", 469),
    _c(110, "12,500 synthetic reads tile the references",
       E + "crossmap_control.tsv", "sum", 12500, column="reads_simulated"),
    # ---------------------------------------------------------- corrections
    _c(111, "20 corrections concern unannotated genes", CORR, "sum", 20,
       column="n", where={"cls": "C1_unannotated"}),
    _c(112, "54 concern genes held only by a non-coding feature", CORR, "sum",
       54, column="n", where={"cls": "C2_noncoding_only"}),
    _c(113, "217 concern incomplete models", CORR, "sum", 217, column="n",
       where={"cls": "C3_incomplete"}),
    _c(114, "1 concerns a wrong paralogue name", CORR, "sum", 1, column="n",
       where={"cls": "C5_wrong_paralog_name"}),
    _c(115, "5 concern protein records named for the wrong family", CORR,
       "sum", 5, column="n", where={"cls": "C6_protein_wrong_family"}),
]
