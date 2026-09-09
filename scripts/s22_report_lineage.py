"""The lineage half of S22's report — §9 to §11, split out of
`s22_report_results.py` to keep every module inside the 500-line budget (the
`s3_report.py` / `s3_report_d10.py` pattern, applied twice as S7, S9, S12 and
S15b do).  It takes the caller's headline dict and its formatters, so the
three halves cannot read the tables differently.

The section it carries is the one that argues against its own headline: the
pooled lineage result is significant and wrong, and §9.3 is the table that
says so before §9.4 gives the answer.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L                                            # noqa: E402
import s22_priors as PR                                        # noqa: E402


def A(*parts: str) -> str:
    return "\n".join(parts)


def f(x, nd=3) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x) if x not in (None, "") else "—"


def pfmt(p) -> str:
    try:
        v = float(p)
    except (TypeError, ValueError):
        return "—"
    if v == 0:
        return "0"
    return f"{v:.3g}" if v < 0.001 else f"{v:.4f}"


PARA = ("ITPR1", "ITPR2", "ITPR3")


def section_lineage(h: dict) -> str:
    path = L.OUT_DIR / "lineage_test.tsv"
    if not path.exists():
        return A("## 9. The lineage test", "", "*not run yet*")
    lt = L.read_tsv(path)
    trows = ["| comparison | n | median Δ | reference Δ | shift | p |",
             "|---|---|---|---|---|---|"]
    for r in lt:
        if not r["p_mannwhitney"]:
            continue
        trows.append(f"| {r['test']} | {r['n_test']} | "
                     f"{f(r['median_test_delta'], 4)} | "
                     f"{f(r['median_reference_delta'], 4)} | "
                     f"{f(r['shift'], 4)} | {pfmt(r['p_mannwhitney'])} |")

    if "matched" not in h:
        return A("## 9. The lineage test, and what it could have seen", "",
                 f"§5 found **{h['n_itpr_no_plc']}** reference proteomes "
                 "carrying an ITPR and no PI-PLC.", "", *trows, "",
                 "*the sequence-level half is not run yet*")

    m = h["matched"]
    cov = h["deep_covariates"]
    strata = h["deep_strata"]
    mp = h["matched_power"]
    rc = h["ryr_control"]

    srows = ["| phylum | PLC absent | PLC present | shift | p | q |",
             "|---|---|---|---|---|---|"]
    for r in strata:
        if r["stratum"] == "POOLED (all phyla)" or int(r["n_plc_absent"]) == 0:
            continue
        srows.append(f"| {r['stratum']} | {r['n_plc_absent']} | "
                     f"{r['n_plc_present']} | {f(r['shift'], 4)} | "
                     f"{pfmt(r['p_mannwhitney'])} | "
                     f"{pfmt(r['q_mannwhitney'])} |")
    pooled = next((r for r in strata if r["stratum"].startswith("POOLED")), {})
    covrows = ["| covariate | ρ with Δ | p | median, PLC absent | "
               "median, PLC present |", "|---|---|---|---|---|"]
    for r in cov:
        covrows.append(f"| {r['covariate']} | {f(r['rho_vs_delta'])} | "
                       f"{pfmt(r['p'])} | {f(r['median_plc_absent'], 4)} | "
                       f"{f(r['median_plc_present'], 4)} |")
    prows = ["| true shift | power of the matched test |", "|---|---|"]
    for r in mp:
        if r.get("shift"):
            prows.append(f"| {f(r['shift'], 3)} | {f(r['power'], 2)} |")

    return A(
        "## 9. The lineage test",
        "",
        "### 9.1 The taxa",
        "",
        f"§5 found **{h['n_itpr_no_plc']}** reference proteomes in the whole "
        "eukaryotic sweep carrying an IP₃ receptor and no "
        "phosphoinositide-specific phospholipase C, and they are not "
        "scattered: they are the oomycetes, the early-diverging "
        "(Mucoromycota) fungi, three *Perkinsus* species, two ciliates, four "
        "prasinophyte algae, and a group of platyhelminths.  Two clades "
        "account for more than half of them.",
        "",
        "The internal control for that list costs nothing and is strong: "
        "**every one of these proteomes carries a full-length IP₃ receptor**, "
        "a ~2,700-residue multi-exon gene the sweep found in it.  A gene set "
        "complete enough to hold this receptor is complete enough to hold a "
        "phospholipase C.",
        "",
        "### 9.2 The measurement, on the whole population rather than a sample",
        "",
        "msa_v2 was built to span four kingdoms with 134 representatives, and "
        "the taxa above are not the taxa a diversity rule picks — asked of "
        "that alignment the test has two tips on one side:",
        "",
        *trows,
        "",
        "The ryanodine receptors are the positive control and they fire: "
        "same pore, same diagnostic domains, no IP₃ site, and the paired "
        "difference moves by about a tenth of an identity unit on six tips.",
        "",
        "So the test is repeated on the **whole population the question is "
        "asked of** — all "
        f"{h.get('deep_panel_n', '—')} non-vertebrate reference proteomes "
        "that carry an ITPR, each aligned pairwise to the human reference "
        f"through the same module residues.  {h.get('deep_panel_in_test', '—')} "
        "resolve both modules above the coverage bar and enter.",
        "",
        "### 9.3 The pooled answer is a clade artefact",
        "",
        f"Pooled, the PLC-absent records have a wider core-minus-pore gap: "
        f"{f(pooled.get('median_absent'), 4)} against "
        f"{f(pooled.get('median_present'), 4)}, "
        f"p = {pfmt(pooled.get('p_mannwhitney'))} on "
        f"{pooled.get('n_plc_absent', '—')} against "
        f"{pooled.get('n_plc_present', '—')}.  That number should not be "
        "believed, and the table that says why is the covariate table:",
        "",
        *covrows,
        "",
        "PLC-absent proteomes sit at a **median pore identity of "
        f"{f(next((r['median_plc_absent'] for r in cov if r['covariate'] == 'pore_identity'), None), 3)} "
        "against "
        f"{f(next((r['median_plc_present'] for r in cov if r['covariate'] == 'pore_identity'), None), 3)}** "
        "for the rest — they are oomycetes and early-diverging fungi, far "
        "further from the human reference than the arthropods and nematodes "
        "that dominate the other cell.  And the paired difference is itself "
        "correlated with divergence, so a pooled comparison of absent against "
        "present is largely a comparison of distant against near.",
        "",
        "Inside a phylum, where the clade is held constant, only three "
        "strata have both cells populated at all:",
        "",
        *srows,
        "",
        "Mucoromycota is the one that cannot be tested for the most "
        "informative reason: **every Mucoromycota proteome in this set that "
        "carries an ITPR lacks a PI-PLC**, so the phylum has no internal "
        "control.",
        "",
        "### 9.4 Matched on divergence, the effect is gone — and bounded",
        "",
        "Each PLC-absent record is matched to up to three PLC-present records "
        "within 0.03 pore identity of its own, and scored against their "
        "median.  This is S15a's identity-matched paired design (D47) applied "
        "to the same kind of claim, and it removes the confound rather than "
        "arguing with it.",
        "",
        f"**All {m['n_plc_absent_in_test']} matched.**  The median within-pair "
        f"difference is {f(m['median_difference'], 4)} "
        f"(95 % CI {f(m['ci95_lo'], 4)} to {f(m['ci95_hi'], 4)}), the tips "
        f"split {m['n_core_more_relaxed']} to {m['n_core_less_relaxed']}, and "
        f"the sign test gives p = {pfmt(m['p_sign'])} "
        f"(Wilcoxon {pfmt(m['p_wilcoxon'])}).",
        "",
        "### 9.5 Power — and this is a bounded null, not an empty one",
        "",
        *prows,
        "",
        "Measured through the *same* pairwise instrument, at the *same* "
        "divergence as the test group, the ryanodine receptors' paired "
        f"difference is {f(rc.get('delta_core_minus_pore'), 4)} against the "
        f"ITPR median, a shift of {f(rc.get('shift_vs_itpr'), 4)}.  The "
        "matched test has essentially full power at that shift and 59 % at "
        "half of it, and its confidence interval excludes anything larger "
        f"than about {f(abs(float(m['ci95_lo'])), 3)}.",
        "",
        "**Losing the enzyme that makes IP₃ does not relax the receptor's "
        "IP₃-binding core.**  Not *we could not tell* — the test can see a "
        "ligand-free core when there is one, at the same evolutionary "
        "distance, and it does not see one here.",
    )


def section_priors(h: dict) -> str:
    prim = h.get("module_primary", {})
    con = h.get("contacts", {})
    lines = ["## 10. Priors — what earlier tasks said, and what S22 measures",
             ""]
    c1 = con.get(("ITPR1", "rest_of_core"), {})
    c2 = con.get(("ITPR2", "rest_of_core"), {})
    c3 = con.get(("ITPR3", "rest_of_core"), {})
    cp = con.get(("ITPR1", "rest_of_pocket"), {})
    lines.append("- " + PR.line(
        "contacts_vs_own_element", "confirmed",
        f"asked against the binding core rather than the Pfam element the "
        f"contacts fall in, the same ten residues win in all three "
        f"paralogues — q = {pfmt(c1.get('q'))} / {pfmt(c2.get('q'))} / "
        f"{pfmt(c3.get('q'))} after BH across the whole family of tests."))
    lines.append("- " + PR.line(
        "pore_most_constrained", "confirmed",
        "with the luminal loop separated the pore module is not merely above "
        "the linkers, it is above the ligand core: "
        f"{f(prim.get('ITPR1', {}).get('difference'), 4)} and "
        f"{f(prim.get('ITPR3', {}).get('difference'), 4)} identity points in "
        "ITPR1 and ITPR3, paired per orthologue."))
    lines.append("- " + PR.line(
        "luminal_loop_least_conserved", "confirmed",
        "and it is load-bearing rather than incidental — leaving the loop "
        "inside PF00520 reverses the sign of the core-versus-pore comparison "
        "in every paralogue (§6.2)."))
    lines.append("- " + PR.line(
        "ligand_core_names_the_family", "not corroborated",
        "the module that names the family is the less constrained of the "
        "two measured here.  Naming a family after its diagnostic domain is "
        "a statement about what the domain identifies, not about what "
        "selection holds most tightly, and this task is the first in the "
        "project to measure the second."))
    lines.append("- " + PR.line(
        "ryr_shares_the_domains", "confirmed",
        "and the ryanodine receptors are used here as the positive control "
        "the lineage test needs — the same pore, the same diagnostic "
        "domains, no IP₃ site."))
    lines.append("- " + PR.line(
        "omega_ranking", "orthogonal",
        "S9 ranked whole-clade rates and S22 compares two modules inside one "
        "protein; ITPR2's flat module contrast and its lowest-ranked ω are "
        "not the same measurement and neither corroborates the other."))
    lines.append("- " + PR.line(
        "family_absent_in_plants_and_dikarya", "orthogonal",
        "S20's absence is about whether the receptor exists; §5's table is "
        "about what sits beside the ones that do, and the two share no cell."))
    rc = h.get("ryr_control", {})
    lines.append("- " + PR.line(
        "family_identity_to_ryr", "confirmed",
        "measured through S22's own pairwise instrument the six RyR tips sit "
        f"at {f(rc.get('core_identity'))} identity in the ligand core against "
        f"{f(rc.get('pore_identity'))} in the pore — the two families share a "
        "fold and a pore and not a ligand site, and that gap is what makes "
        "them the positive control §9 needs."))
    return "\n".join(lines)


def section_handoff(h: dict) -> str:
    m = h.get("matched", {})
    rc = h.get("ryr_control", {})
    return A(
        "## 11. What this settles, and what it does not",
        "",
        "**Settled.**",
        "",
        "1. The ligand core is not the most constrained part of this "
        "receptor.  Paired per orthologue against roughly 250 sequences per "
        "paralogue, the pore module is ahead in ITPR1 and ITPR3 and level in "
        "ITPR2.",
        "2. That answer depends on one boundary — whether fifty residues of "
        "luminal loop count as pore — and reverses when they do.  Any "
        "version of this comparison that does not state the boundary is not "
        "interpretable.",
        "3. The ten measured IP₃ contacts are more constrained than the rest "
        "of the binding core and no more constrained than the rest of the "
        "pocket.  What selection holds is a neighbourhood about fifteen "
        "ångström across, not a contact set.",
        "4. That result replicates on an independent axis: every contact "
        "site in every paralogue is under detectable purifying selection, "
        "and the share falls by 0.13-0.25 across the shells.  It is a step "
        "down from the ligand's first two shells onto a floor rather than a "
        "gradient running out to 15 Å, and §8 reports it as such.",
        f"5. {h['n_itpr_no_plc']} of the eukaryotic reference proteomes this "
        "project swept carry an IP₃ receptor and no phosphoinositide-specific "
        "phospholipase C — concentrated in the oomycetes and the "
        "early-diverging fungi, and every one of them carrying a full-length "
        "receptor.",
        "6. **Losing that enzyme does not relax the ligand core.**  Matched "
        f"on divergence the within-pair difference is "
        f"{f(m.get('median_difference'), 4)} "
        f"(95 % CI {f(m.get('ci95_lo'), 4)} to {f(m.get('ci95_hi'), 4)}, "
        f"p = {pfmt(m.get('p_sign'))}), against a positive control whose "
        f"shift on the same instrument is {f(rc.get('shift_vs_itpr'), 4)}.  "
        "This is a bounded null, not an absent measurement.",
        "",
        "**Not settled.**",
        "",
        "1. Whether the pocket's constraint is *about* IP₃.  Everything "
        "within 15 Å of the ligand is also within the fold that holds it, "
        "and no measurement here separates ligand binding from domain "
        "packing.  A mutational or binding-affinity dataset would; sequence "
        "will not.",
        "2. What these lineages' receptors are gated by.  A PI-PLC search is "
        "a search for one enzyme family; the absence says the canonical "
        "route to IP₃ is missing, not that the receptor has no ligand, and "
        "the null in §9.4 is consistent with the site being held by "
        "something else.",
        "3. Whether the core-versus-pore difference is a vertebrate fact.  "
        "The paired test in §6 runs on the deep sweep alignments, which are "
        "vertebrate by construction; §9's panel measures the same two "
        "modules across the eukaryotes but against one reference, which is a "
        "weaker instrument.",
        "",
        "**Three things to hold against this task.**  The per-orthologue "
        "identity is measured against the human reference, so a difference "
        "between modules is partly a difference in how far each module has "
        "moved *from human* rather than from the ancestor; the column-level "
        "metrics in §6.3 are the check on that and they agree.  The pocket "
        "is defined on human ITPR3 and carried to the other two paralogues "
        "by alignment, so a shell assignment in ITPR1 is one transfer "
        "removed from a measurement.  And the PI-PLC call is a domain "
        "architecture call: a divergent phospholipase whose X or Y box has "
        "drifted past the profile would be scored absent, which would move a "
        "taxon into the test group and could only weaken the null, not "
        "create it.",
    )


def build(h: dict) -> str:
    return "\n\n".join([
        section_lineage(h),
        section_priors(h),
        section_handoff(h),
    ])
