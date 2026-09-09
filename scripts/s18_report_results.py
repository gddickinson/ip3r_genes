"""S18's results half — the two records, the correction list and the priors.

Split out of `s18_report.py` to keep both under 500 lines, and taking the
caller's loader and formatters so the two halves cannot read the tables
differently (the `s3_report.py` / `s3_report_d10.py` pattern).

The section order is the order a reader needs it in: what the annotation does
with the gene, who produced that annotation, whether it is this family or
vertebrate annotation in general, what the protein records say, what the
proteomes that returned nothing turn out to hold — and only then the
corrections, because a correction list read before its evidence is a list of
assertions.

The priors and the hand-off live in `s18_report_priors.py`, a third piece of
the same split, because this half would otherwise pass the 500-line budget.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_report_priors as RP                                    # noqa: E402

FAILS = ("unannotated", "noncoding", "fragmentary", "split")


def _row(rows, **kw):
    for r in rows:
        if all(str(r.get(k)) == str(v) for k, v in kw.items()):
            return r
    return {}


def results(load, load_json, num, pct, share, pfmt, table, MISSING
            ) -> list[str]:
    stats = load_json("annotation_audit_stats.json")
    H = stats.get("headline", {})
    L: list[str] = []
    L += _s5_states(load, num, pct, share, table, MISSING, H)
    L += _s6_source(load, num, pct, share, pfmt, table, MISSING)
    L += _s7_control(load, num, pct, pfmt, table, MISSING, H)
    L += _s8_protein(load, num, pct, share, table, MISSING, H)
    L += _s9_zero(load, num, table, MISSING, H)
    L += _s10_corrections(load, num, table, MISSING, H)
    L += RP.priors(load, num, pct, H)
    L += RP.handoff(pct, H)
    return L


# --------------------------------------------------------------------------
def _s5_states(load, num, pct, share, table, MISSING, H) -> list[str]:
    cells = load("state_by_cell.tsv")
    if not cells:
        return ["## 5. What the annotation does with the gene", "", MISSING]
    order = ["ITPR1", "ITPR2", "ITPR3", "RYR"]
    cells = sorted([r for r in cells if r["cell"] in order],
                   key=lambda r: order.index(r["cell"]))
    cls = load("state_by_class.tsv")
    ranked = sorted([r for r in cls if int(r["n_loci"]) >= 20],
                    key=lambda r: -float(r["frac_complete"]))
    # Top three and bottom three, deduplicated: with six qualifying classes
    # the two slices overlap and the table prints Actinopteri twice.
    worst, seen = [], set()
    for r in ranked[:3] + ranked[-3:]:
        if r["vclass"] not in seen:
            seen.add(r["vclass"])
            worst.append(r)
    L = [
        "## 5. What the annotation does with the gene",
        "",
        f"Of the {num(H.get('n_scorable_itpr_loci'))} ITPR loci in assemblies "
        f"that ship a gene set, **{pct(H.get('frac_itpr_complete'))}** are "
        f"delivered as one complete annotated model and "
        f"**{pct(H.get('frac_itpr_any_failure'))}** are not. The failures "
        f"are not evenly shaped: {num(H.get('n_itpr_unannotated'))} loci have "
        f"no same-strand annotated feature at all and "
        f"{num(H.get('n_itpr_noncoding'))} are held only by a non-coding one, "
        f"so the database serves no protein for either.",
        "",
        "*(Figure `s18_fig3_family_vs_control` panel A.)*",
        "",
    ]
    L += table(cells, [
        ("cell", "cell"), ("loci", "n_loci"),
        ("complete", "n_complete"), ("split", "n_split"),
        ("fragmentary", "n_fragmentary"), ("non-coding only", "n_noncoding"),
        ("unannotated", "n_unannotated"),
        ("complete (share)", "frac_complete"),
    ])
    L += ["",
          "The single largest determinant is not the paralog but the "
          f"assembly: above D4's contiguity bar the ITPR failure rate falls "
          f"from {pct(H.get('frac_itpr_any_failure'))} to "
          f"**{pct(H.get('frac_itpr_failure_above_d4'))}**. Two thirds of "
          "what looks like an annotation problem is a contig too short to "
          "hold a 2,700-residue gene.",
          ""]
    if worst:
        L += [f"### 5.1 By vertebrate class", "",
              f"The best and worst three of the {len(ranked)} classes with at "
              f"least 20 scorable loci:", ""]
        L += table(worst, [
            ("class", "vclass"), ("loci", "n_loci"),
            ("complete", "frac_complete"),
            ("fragmentary", "frac_fragmentary"),
            ("unannotated", "frac_unannotated"),
        ])
        L += [""]
    return L


def _s6_source(load, num, pct, share, pfmt, table, MISSING) -> list[str]:
    src = load("by_source.tsv")
    if not src:
        return ["## 6. Who produced the annotation (D9)", "", MISSING]
    all_scope = [r for r in src if r["scope"] == "all scorable loci"]
    d4 = [r for r in src if r["scope"].startswith("loci on a contig")]
    a_comp = _row(all_scope, state="complete")
    d_comp = _row(d4, state="complete")
    L = [
        "## 6. Who produced the annotation (D9)",
        "",
        "D9 says a RefSeq gene set and a submitter-deposited GenBank one are "
        "not comparable evidence. They are not close.",
        "",
        f"Over every scorable locus, **{pct(a_comp.get('frac_refseq'))}** of "
        f"the {num(a_comp.get('n_refseq_total'))} loci in RefSeq (`GCF_`) "
        f"annotations are complete against **{pct(a_comp.get('frac_genbank'))}"
        f"** of the {num(a_comp.get('n_genbank_total'))} in submitter "
        f"GenBank (`GCA_`) ones. Every state differs at "
        f"q ≤ {pfmt(max(float(r['q_bh']) for r in all_scope))}.",
        "",
        "*(Figure `s18_fig1_by_source`.)*",
        "",
    ]
    L += table(all_scope, [
        ("state", "state"), ("RefSeq n", "n_refseq"),
        ("RefSeq", "frac_refseq"), ("GenBank n", "n_genbank"),
        ("GenBank", "frac_genbank"), ("q (BH)", "q_bh"),
    ])
    L += ["",
          "### 6.1 The confounder, controlled",
          "",
          "GenBank assemblies are less contiguous, and a locus on a contig "
          "too short to hold the gene cannot be annotated completely by "
          "anybody — so the contrast is run again over the loci that clear "
          "D4's bar. It narrows and does not close: "
          f"**{pct(d_comp.get('frac_refseq'))}** against "
          f"**{pct(d_comp.get('frac_genbank'))}** on "
          f"{num(d_comp.get('n_refseq_total'))} and "
          f"{num(d_comp.get('n_genbank_total'))} loci "
          f"(q = {pfmt(d_comp.get('q_bh'))}). About "
          f"{share(float(d_comp.get('frac_genbank', 0)) - float(a_comp.get('frac_genbank', 0)), 1.0 - float(a_comp.get('frac_genbank', 1)))}"
          " of the gap between the two archives is assembly quality; the "
          "rest is the gene set.",
          ""]
    L += table(d4, [
        ("state", "state"), ("RefSeq n", "n_refseq"),
        ("RefSeq", "frac_refseq"), ("GenBank n", "n_genbank"),
        ("GenBank", "frac_genbank"), ("q (BH)", "q_bh"),
    ])
    L += [""]
    return L


def _s7_control(load, num, pct, pfmt, table, MISSING, H) -> list[str]:
    fvc = load("family_vs_control.tsv")
    if not fvc:
        return ["## 7. Is it this family, or is it vertebrate annotation?",
                "", MISSING]
    a = _row(fvc, scope="all scorable loci", measure="any annotation failure")
    d = [r for r in fvc if r["scope"].startswith("loci on a contig")
         and r["measure"] == "any annotation failure"][0]
    nc = _row(fvc, scope="all scorable loci", measure="noncoding")
    L = [
        "## 7. Is it this family, or is it vertebrate annotation?",
        "",
        "The audit's premise is that a ~2,700-residue, ~58-exon gene with a "
        "sister family sharing every diagnostic domain should be badly "
        "recorded. The control says otherwise, and that is the result.",
        "",
        f"The ITPR cells fail on **{pct(a['frac_itpr_failing'])}** of "
        f"{num(a['n_itpr'])} loci; the ryanodine receptors, in the same "
        f"assemblies through the same pipelines, on "
        f"**{pct(a['frac_ryr_failing'])}** of {num(a['n_ryr'])} "
        f"(q = {pfmt(a['q_bh'])}). Above D4's bar the two are "
        f"{pct(d['frac_itpr_failing'])} and {pct(d['frac_ryr_failing'])} "
        f"(q = {pfmt(d['q_bh'])}). **No overall difference survives "
        f"correction.** How this family is recorded is how vertebrate genes "
        f"of this size are recorded.",
        "",
        "*(Figure `s18_fig3_family_vs_control` panel B.)*",
        "",
    ]
    L += table([r for r in fvc if r["scope"] == "all scorable loci"], [
        ("measure", "measure"), ("ITPR n", "n_itpr_failing"),
        ("ITPR", "frac_itpr_failing"), ("RyR n", "n_ryr_failing"),
        ("RyR", "frac_ryr_failing"), ("q (BH)", "q_bh"),
    ])
    ncd = [r for r in fvc if r["scope"].startswith("loci on a contig")
           and r["measure"] == "noncoding"][0]
    L += ["",
          "One state does separate, and it is the family-specific one: an "
          f"ITPR locus is **{float(nc['frac_itpr_failing']) / max(1e-9, float(nc['frac_ryr_failing'])):.1f}×** "
          f"more likely than a RyR locus to be held only by a non-coding "
          f"feature — {nc['n_itpr_failing']} loci against "
          f"{nc['n_ryr_failing']} (q = {pfmt(nc['q_bh'])}). That is the "
          "*Podiceps* failure at scale: a gene the annotation identifies "
          "correctly, names correctly, and files as non-coding, so no "
          "protein record is ever created and no name-based search can reach "
          "it.",
          "",
          "**And it does not survive the contiguity control.** Above D4's bar "
          f"the same comparison is {ncd['n_itpr_failing']} against "
          f"{ncd['n_ryr_failing']} on {ncd['n_itpr']} and {ncd['n_ryr']} loci "
          f"(q = {pfmt(ncd['q_bh'])}) — the effect is confined to assemblies "
          "too broken to carry the gene, which is where a submitter has most "
          "reason to file a model as non-coding in the first place. The "
          "excess is real in the record set and this measurement cannot say "
          "it is a fact about the family rather than about the assemblies the "
          "family's loci happen to sit in.",
          ""]
    return L


def _s8_protein(load, num, pct, share, table, MISSING, H) -> list[str]:
    pv = load("protein_name_verdicts.tsv")
    pf = load("pfam_recall.tsv")
    if not pv:
        return ["## 8. What the protein databases call these records", "",
                MISSING]
    n = float(H.get("n_protein_records") or 0)
    unusable = (float(H.get("n_symbol_placeholder") or 0)
                + float(H.get("n_symbol_absent") or 0))
    L = [
        "## 8. What the protein databases call these records",
        "",
        f"Each of the {num(n)} full-length family protein records was scored "
        f"against the committed 38-bait panel by blastp: best bit score in "
        f"each family, assigned to the winner only when it beats the loser by "
        f"more than D7's relative margin, then the paralog inside the winning "
        f"family — and the record's own gene symbol and protein name read "
        f"through the *same* verdict rule the genome half uses.",
        "",
        f"**The family call is not in dispute.** The sequence disagrees with "
        f"the census call on {num(H.get('n_protein_seq_family_disagrees'))} "
        f"of {num(n)} records, and the databases name the sister family at "
        f"only {num(H.get('n_name_wrong_family'))} — every one of them a "
        f"non-vertebrate record whose best family score is under 200 bits, "
        f"i.e. where the sequence barely separates the families either. "
        f"D14's hazard, which S0 measured at 49 % of zebrafish PF08709 "
        f"records, does not appear as a *naming* error in the full-length "
        f"record set.",
        "",
        f"**The paralog call is not in dispute either.** "
        f"{num(H.get('n_symbol_wrong_paralog'))} of "
        f"{num(H.get('n_protein_vertebrate'))} vertebrate records carry a "
        f"symbol naming a paralog the panel assigns elsewhere, and a further "
        f"{num(H.get('n_paralog_not_callable'))} claim RYR3 — which this "
        f"panel cannot call, because S5's slot table records no RYR3 bait, "
        f"so those are reported as outside the instrument's reach rather "
        f"than as errors.",
        "",
        f"**What is in dispute is whether the records are findable.** "
        f"{num(H.get('n_symbol_placeholder'))} records carry a placeholder "
        f"gene symbol and {num(H.get('n_symbol_absent'))} carry none at all: "
        f"**{share(unusable, n)} of the family's full-length protein records "
        f"have no usable gene symbol.** On the protein-name side "
        f"{num(H.get('n_name_family_ambiguous'))} are named for the "
        f"superfamily — \"RyR/IP3R Homology associated domain-containing "
        f"protein\" — which names the family and its sister together and "
        f"therefore separates neither.",
        "",
        "*(Figure `s18_fig4_protein_side` panel A.)*",
        "",
    ]
    if pf:
        both = _row(pf, axis="which instruments called it", bucket="both")
        v4 = _row(pf, axis="which instruments called it", bucket="v4+s20")
        L += [
            "### 8.1 Pfam recall — would a signature query have found them?",
            "",
            "The signature that *names* this family is PF08709. A reader "
            "looking for the IP3 receptors would query it, so the recall "
            "question is what that query returns out of the records the "
            "project's own two instruments call family.",
            "",
        ]
        L += table([r for r in pf if r["axis"] == "census call"]
                   + [r for r in pf if r["axis"] == "which instruments called it"],
                   [("axis", "axis"), ("bucket", "bucket"),
                    ("records", "n_records"),
                    ("carry PF08709", "frac_with_naming_pfam"),
                    ("complete architecture", "frac_complete_architecture")])
        L += ["",
              f"A PF08709 query reaches {pct(both.get('frac_with_naming_pfam'))} "
              f"of the records both instruments agree on, and "
              f"{pct(v4.get('frac_with_naming_pfam'))} of those only one "
              f"instrument found. The signature is a good but not complete "
              f"index of its own family, and the deficit is concentrated "
              f"outside the vertebrates.",
              ""]
        grp = [r for r in pf if r["axis"] == "taxonomic group"]
        L += table(sorted(grp, key=lambda r: float(r["frac_with_naming_pfam"])),
                   [("group", "bucket"), ("records", "n_records"),
                    ("carry PF08709", "frac_with_naming_pfam"),
                    ("complete architecture", "frac_complete_architecture")])
        L += [""]
    return L


def _s9_zero(load, num, table, MISSING, H) -> list[str]:
    z = load("zero_hit_proteomes.tsv")
    if not z:
        return ["## 9. The proteomes that returned nothing", "", MISSING]
    L = [
        "## 9. The proteomes that returned nothing",
        "",
        f"S3 swept 764 vertebrate reference proteomes with both profiles and "
        f"{num(H.get('n_zero_hit_proteomes'))} came back with no family hit "
        f"at all. On its own that is uninterpretable — an ITPR-shaped hole in "
        f"a proteome is either a gene the species lacks or a gene its gene "
        f"caller did not find — so each was resolved against an assembly of "
        f"its own species, which the S5 sweep searched with an instrument "
        f"that owes the gene caller nothing.",
        "",
        f"**All {num(H.get('n_zero_hit_gene_caller'))} are gene-caller "
        f"failures.** Every one of the 15 species has a genome in the S4 "
        f"scope, and in every one the genomic sweep recovers at least one "
        f"ITPR cell at over half the bait's length while the proteome holds "
        f"none. Not one is `genome_also_empty`, and not one is "
        f"`undecidable_no_genome` — the verdict that exists so a species with "
        f"no assembly in scope could not be reported as an absence.",
        "",
        "*(Figure `s18_fig4_protein_side` panel B.)*",
        "",
    ]
    L += table(sorted(z, key=lambda r: -int(r["proteome_proteins"])), [
        ("proteome", "proteome_id"), ("species", "species"),
        ("class", "vclass"), ("proteins", "proteome_proteins"),
        ("assembly", "genome_accession"), ("ITPR loci", "n_itpr_loci"),
        ("cells recovered", "n_recovered"), ("verdict", "verdict"),
    ])
    L += [""]
    return L


def _s10_corrections(load, num, table, MISSING, H) -> list[str]:
    summ = load("correction_summary.tsv")
    if not summ:
        return ["## 10. The correction list", "", MISSING]
    corr = load("corrections.tsv")
    high = [r for r in corr if r["priority"] == "high" and r["vetoed"] != "1"]
    L = [
        "## 10. The correction list",
        "",
        f"{num(H.get('n_corrections'))} correction items, "
        f"{num(H.get('n_corrections_high'))} of them `high` priority, "
        f"{num(H.get('n_corrections_vetoed'))} withheld under **D6**. Each "
        f"row carries the assembly, the coordinates, the current state, the "
        f"current name, the proposal, the archived evidence file a curator "
        f"can open, and the numbers the class fired on. → "
        f"`results/annotation_audit/corrections.tsv`",
        "",
        "**Priority is a rule, not an impression.** `high` needs the gene "
        "demonstrably present *and* the assembly demonstrably able to carry "
        "it *and* the reading frame intact; a partial recovery or an unscored "
        "ORF is `medium`; anything below D4's contiguity bar is `low`, "
        "because there the annotation's silence may be the assembly's fault "
        "and not the annotator's.",
        "",
        "**D6 is applied as a column, not a filter.** If S15's ORF screen "
        "scored a locus lesion-rich, the audit does not tell RefSeq to "
        "resurrect it — but the row is still written, flagged, with its "
        "reason, because a locus this audit declined to correct is evidence "
        "about the audit.",
        "",
    ]
    L += table(summ, [
        ("class", "cls"), ("priority", "priority"), ("items", "n"),
        ("vetoed (D6)", "n_vetoed"), ("RefSeq", "n_refseq"),
        ("GenBank", "n_genbank"),
    ])
    if high:
        L += ["", "### 10.1 The first ten `high`-priority items", "",
              "Sorted as the table is written; the full list is the file.", ""]
        L += table(high[:10], [
            ("assembly", "assembly"), ("species", "organism"),
            ("cell", "cell"), ("locus", "contig"),
            ("state", "current_state"), ("current name", "current_name"),
            ("proposal", "proposal"),
        ])
    L += [""]
    return L


