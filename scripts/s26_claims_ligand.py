"""S26 ledger: the load-bearing numbers of the ligand-site paper.

`CARRY` names claims already declared by the manuscript's ledger (`C…`), the
thesis's (`T…`) or another paper, by identifier; `CLAIMS` declares numbers no
earlier ledger has, with the `LI` prefix `s26_claims.PREFIX` gives this
paper. Every value must appear in the paper's own text.

Several manuscript claims check a value the paper prints rounded (a q-value
of 0.01725 is written 0.017). Those are declared here again with a tolerance
rather than carried, because the appearance rule needs the printed value; the
carried and the re-declared rows read the same cell of the same table.
"""

_T = "results/ligand_site/"

CARRY: list[str] = [
    "C222", "C223", "C225", "C231", "C232", "C232b", "C232c",
    "C233", "C234", "C236", "C237", "C239", "C240", "C241",
    "C242", "C243", "C244", "C246", "C246b", "C246c", "C247",
    "C248c", "C248d", "C249", "C250", "C251", "C252", "C253",
    "C254", "C255", "C257", "T58",
]


def _c(i, claim, source, op, expect, **kw):
    return dict(id=f"LI{i:02d}", claim=claim, source=_T + source, op=op,
                expect=expect, **kw)


CLAIMS: list[dict] = [
    # ---------------------------------------------------------- modules
    _c(1, "the primary ITPR1 ligand core starts at residue 265",
       "module_map.tsv", "cell", "265", column="start",
       where={"paralog": "ITPR1", "definition": "contact_span"}),
    _c(2, "and ends at residue 569", "module_map.tsv", "cell", "569",
       column="end", where={"paralog": "ITPR1", "definition": "contact_span"}),
    _c(3, "the primary ITPR1 pore module starts at residue 2363",
       "module_map.tsv", "cell", "2363", column="start",
       where={"paralog": "ITPR1", "definition": "channel_minus_luminal"}),
    _c(4, "and ends at residue 2608", "module_map.tsv", "cell", "2608",
       column="end",
       where={"paralog": "ITPR1", "definition": "channel_minus_luminal"}),
    _c(5, "ITPR1 drops three orthologues on the coverage rule",
       "module_contrast.tsv", "cell", "3", column="n_dropped",
       where={"paralog": "ITPR1", "is_primary": "True"}),
    _c(6, "ITPR2 shows no core-pore difference (q = 0.13)",
       "module_contrast.tsv", "cell", "0.13", column="q_wilcoxon", tol=0.005,
       where={"paralog": "ITPR2", "is_primary": "True"}),
    _c(7, "with the published core, ITPR2 tips to a core lead of 0.0049",
       "module_contrast.tsv", "cell", "0.0049", column="mean_difference",
       where={"paralog": "ITPR2", "core_definition": "ibc_literature",
              "pore_definition": "channel_minus_luminal"}),
    # ---------------------------------------------------------- pocket
    _c(8, "8TKG adds two contacts to the ten", "shell_agreement.tsv", "cell",
       "12", column="n_contact_le_4.5A", where={"pdb_id": "8TKG"}),
    _c(9, "Ala276 is a contact in four of the six depositions",
       "shell_agreement.tsv", "count", "4",
       where={"extra_contacts__in": "276;411,276"}),
    _c(10, "the pocket holds 125 residues", "shell_trend.tsv", "cell", "125",
       column="n_residues", where={"paralog": "ITPR1"}),
    _c(11, "the contact shell holds twelve residues", "shell_constraint.tsv",
       "cell", "12", column="n_residues",
       where={"paralog": "ITPR1", "shell": "contact"}),
    _c(12, "the outermost shell holds 59 residues", "shell_constraint.tsv",
       "cell", "59", column="n_residues",
       where={"paralog": "ITPR1", "shell": "fourth"}),
    _c(13, "the ITPR1 whole-protein mean conservation is 0.7366",
       "shell_constraint.tsv", "cell", "0.7366",
       column="whole_protein_mean_jsd",
       where={"paralog": "ITPR1", "shell": "contact"}),
    _c(14, "ITPR1 contacts score 0.813", "contact_test.tsv", "cell", "0.813",
       column="mean_contact_jsd", tol=0.0005,
       where={"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR1", "layer": "deep"}),
    _c(15, "against 0.744 for the rest of the core", "contact_test.tsv",
       "cell", "0.744", column="mean_background_jsd", tol=0.0005,
       where={"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR1", "layer": "deep"}),
    _c(16, "contacts beat the rest of the core in ITPR2 (q = 0.017)",
       "contact_test.tsv", "cell", "0.017", column="q_permutation",
       tol=0.0005,
       where={"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR2", "layer": "deep"}),
    _c(17, "and in ITPR3 (q = 0.041)", "contact_test.tsv", "cell", "0.041",
       column="q_permutation", tol=0.0005,
       where={"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR3", "layer": "deep"}),
    _c(18, "against the rest of the pocket the smallest q is 0.12",
       "contact_test.tsv", "cell", "0.12", column="q_permutation",
       tol=0.005,
       where={"contact_set": "s0_contact", "background": "rest_of_pocket",
              "paralog": "ITPR2", "layer": "deep"}),
    _c(19, "the permutation ran 200,000 iterations", "contact_test.tsv",
       "cell", "200000", column="iters",
       where={"contact_set": "s0_contact", "background": "rest_of_core",
              "paralog": "ITPR1", "layer": "deep"}),
    _c(20, "conservation falls with distance in ITPR2 (rho = -0.436)",
       "shell_trend.tsv", "cell", "-0.436", column="rho_distance_vs_jsd",
       tol=0.0005, where={"paralog": "ITPR2"}),
    # ---------------------------------------------------------- rate axis
    _c(21, "the ITPR1 core's mean beta is 0.0251", "omega_module_test.tsv",
       "cell", "0.0251", column="mean_core_beta",
       where={"paralog": "ITPR1", "is_primary": "True"}),
    _c(22, "against 0.0446 for the pore", "omega_module_test.tsv", "cell",
       "0.0446", column="mean_pore_beta",
       where={"paralog": "ITPR1", "is_primary": "True"}),
    _c(23, "the ITPR3 core evolves faster at q = 0.058",
       "omega_module_test.tsv", "cell", "0.058", column="q_mannwhitney",
       tol=0.0005, where={"paralog": "ITPR3", "is_primary": "True"}),
    # ---------------------------------------------------------- enzyme
    _c(24, "22 fungal proteomes carry an ITPR and no PI-PLC",
       "pathway_cooccurrence.tsv", "cell", "22", column="n_itpr_no_plc",
       where={"group": "fungi"}),
    _c(25, "28 other protist proteomes do", "pathway_cooccurrence.tsv",
       "cell", "28", column="n_itpr_no_plc",
       where={"group": "protista_other"}),
    _c(26, "10 non-vertebrate animal proteomes do",
       "pathway_cooccurrence.tsv", "cell", "10", column="n_itpr_no_plc",
       where={"group": "metazoa_nonvert"}),
    _c(27, "4 green-plant proteomes do", "pathway_cooccurrence.tsv", "cell",
       "4", column="n_itpr_no_plc", where={"group": "viridiplantae"}),
    _c(28, "the 2,764 non-vertebrate eukaryotic proteomes of the panel",
       "plc_repertoire.tsv", "count", "2764",
       where={"group__in": "fungi,metazoa_nonvert,protista_other,"
                            "viridiplantae"}),
    # ---------------------------------------------------------- lineage
    _c(29, "all 662 non-vertebrate receptors from reference proteomes",
       "deep_lineage_panel.tsv", "count", "662"),
    _c(30, "the RyR control shifts by -0.0996 on the representative alignment",
       "lineage_test.tsv", "cell", "-0.0996", column="shift",
       where={"test": "ryr_control vs every ITPR tip"}),
    _c(31, "on six RyR tips", "lineage_test.tsv", "cell", "6",
       column="n_test", where={"test": "ryr_control vs every ITPR tip"}),
    _c(32, "all 15 Mucoromycota receptors lack the enzyme",
       "deep_lineage_test.tsv", "cell", "15", column="n_plc_absent",
       where={"stratum": "Mucoromycota"}),
    _c(33, "only three phylum strata are testable", "deep_lineage_test.tsv",
       "count", "3", where={"testable": "True", "note": ""}),
    _c(34, "matched within 0.03 pore identity",
       "deep_lineage_matched_summary.tsv", "cell", "0.03",
       column="match_tolerance"),
    _c(35, "17 matched pairs fall the other way",
       "deep_lineage_matched_summary.tsv", "cell", "17",
       column="n_core_less_relaxed"),
    _c(36, "sign test p = 0.87", "deep_lineage_matched_summary.tsv", "cell",
       "0.87", column="p_sign", tol=0.005),
    _c(37, "power 0.99 at a shift of 0.05", "deep_lineage_matched_power.tsv",
       "cell", "0.99", column="power", tol=0.005, where={"shift": "0.05"}),
    _c(38, "power 0.59 at a shift of 0.03", "deep_lineage_matched_power.tsv",
       "cell", "0.59", column="power", tol=0.005, where={"shift": "0.03"}),
]
