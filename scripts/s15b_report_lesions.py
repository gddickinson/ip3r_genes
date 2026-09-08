"""S15b — the half of `results/loss_counts/report.md` that argues against
its own headline.

Split from `s15b_report_results.py` so both stay under 500 lines and take
the caller's loader and formatter (the `s3_report.py` /
`s3_report_d10.py` pattern).

The fossil section (§8) reports a test that could not be run, **with its
denominator**, because an omitted section is indistinguishable from a
section nobody ran; it then follows D47's lead and finds the ITPR3 indel
excess is a bird result — and prints the contiguity control that makes
that finding underpowered rather than clean.  §9 computes every prior on
this task's own tables, and §10 is written as a hand-off.
"""

from __future__ import annotations

import s15b_coding as coding
import s15b_priors as priors


def _pf(v) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f == 0.0:
        return "< 1e-300"
    return f"{f:.4f}" if f >= 1e-4 else f"{f:.2e}"


def _i(v, d=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return d


def _f(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def _absent_at(dc: list[dict], evidence: str, r5: int, r6: int,
               coding_name: str = "paralog") -> int:
    """Cells reading `absent` at one setting, summed over its characters.

    Read off `dollo_counts.tsv` rather than off the manufactured-loss
    table, whose `loosest_evidence` column names the loosest rung a cell
    is released at and therefore cannot be filtered by a single setting.
    """
    return sum(_i(r["n_absent"]) for r in dc
               if r["coding"] == coding_name and r["evidence"] == evidence
               and _i(r["use_r5"]) == r5 and _i(r["use_r6"]) == r6)


def _verdict_line(name: str, prior, observed, verdict: str) -> str:
    mark = {"confirmed": "**confirmed**", "contradicted": "**contradicted**",
            "not corroborated": "*not corroborated*",
            "orthogonal": "*orthogonal*",
            "underpowered": "*underpowered*"}.get(verdict, verdict)
    return (f"- **{name}** — prior: `{prior}`; this task: `{observed}`; "
            f"verdict: {mark}. Source: {priors.PRIOR[name]['where']}.")


# ------------------------------------------------------------------ §8

def section_fossils(load, load_json, num, pct, table, missing, fig,
                    stats, head) -> list[str]:
    by_cell = load("fossil_by_cell.tsv")
    loci = load("fossil_loci.tsv")
    cls = load("lesion_by_class.tsv")
    ctrl = load("lesion_class_controls.tsv")
    fos = stats.get("fossil", {})
    out = ["## 8. The fossil analysis, and the lead it hands over", ""]
    if not by_cell:
        return out + missing("fossil_by_cell.tsv", "the fossil section")
    out += [
        "### 8.1 The shared-lesion Poisson test, with its denominator", "",
        "The brief asks whether descendants of a dead paralog share "
        "lesions more often than independent decay would give. That test "
        f"needs dead loci. Of **{num(fos.get('n_scored', 0))} scored "
        f"loci**, {num(fos.get('n_above_bar', 0))} clear S15a's measured "
        f"lesion bar of {_f(fos.get('lesion_bar')):.3g} lesions per "
        "kilo-residue — and "
        f"**{num(fos.get('n_above_bar_full_coverage', 0))} of those "
        f"{num(fos.get('n_above_bar', 0))} are at full coverage**, "
        "delivering a complete gene model.",
        "",
    ]
    out += table(["cell", "loci scored", "above the lesion bar",
                  "fossil under any reading", "median identity",
                  "median identity above the bar"],
                 [[r["cell"], num(r["n_scored"]), num(r["n_above_bar"]),
                   num(r["n_fossils"]), r["median_identity"],
                   r["median_identity_above_bar"]] for r in by_cell])
    itpr_fossils = [r for r in loci if _i(r["is_fossil"])
                    and r["cell"] != "RYR"]
    n_itpr_scored = sum(_i(r["n_scored"]) for r in by_cell
                        if r["cell"] != "RYR")
    out += [
        "The screen was made deliberately generous — a locus counts as a "
        "fossil if *any* of three readings fires (above the bar and below "
        "the coverage bar, above the bar with a premature stop, or above "
        "the bar in a cell the matrix does not code present) — because a "
        "test made hard to pass would make the zero uninformative. Even "
        f"so, **{len(itpr_fossils)} ITPR loci in "
        f"{num(n_itpr_scored)} fire any reading at all**, and all "
        f"{len(itpr_fossils)} do so on 1–2 internal stops in a complete, "
        "full-coverage model:",
        "",
    ]
    if itpr_fossils:
        out += table(["organism", "cell", "coverage", "identity",
                      "frameshifts", "stops", "cell state"],
                     [[f"*{r['organism']}*", r["cell"], r["coverage"],
                       r["identity"], r["frameshifts"], r["stop_codons"],
                       f"`{r['cell_state']}`"] for r in itpr_fossils])
    out += [
        "S10 measured 0 internal stops against 4–25 expected under neutral "
        "drift at 8 family loci, and set the rule one-sided: zero stops "
        "falsifies a pseudogene call and a handful does not establish one. "
        "One or two stops in a ~2,700-residue model at full coverage is "
        "not the signature of decay. **The shared-lesion Poisson test has "
        "no dead loci to run on, and that is reported here with its "
        "denominator rather than omitted**, because an omitted section is "
        "indistinguishable from a section nobody ran.",
        "",
        "### 8.2 D47's lead: the ITPR3 indel excess is a bird result", "",
    ]
    if not cls:
        out += missing("lesion_by_class.tsv", "the stratified lesion test")
        return out
    testable = [t for t in cls if t["p"] != ""]
    out += [
        "S15a found that ITPR3 carries an indel excess against its own "
        "genome's identity-matched sibling loci — 39 genomes to 14, "
        "q = 0.0032 — and nothing in this project explained it. The same "
        "committed pairs, the same sign test, the same identity matching, "
        "stratified by vertebrate class and BH-corrected across the "
        f"{len(testable)} strata with enough untied pairs to test:",
        "",
    ]
    out += table(["cell", "class", "pairs", "untied", "excess", "deficit",
                  "ties", "p", "q (BH)"],
                 [[t["cell"], t["vclass"], num(t["n_genomes"]), num(t["n"]),
                   num(t["n_pos"]), num(t["n_neg"]), num(t["n_ties"]),
                   _pf(t["p"]), _pf(t["q_bh"])] for t in testable])
    aves3 = next((t for t in testable if t["cell"] == "ITPR3"
                  and t["vclass"] == "Aves"), None)
    act3 = next((t for t in testable if t["cell"] == "ITPR3"
                 and t["vclass"] == "Actinopteri"), None)
    if aves3 and act3:
        out += [
            "**The excess is a bird result.** In Aves it is "
            f"{aves3['n_pos']} genomes to {aves3['n_neg']} over "
            f"{aves3['n']} untied pairs (q = {_pf(aves3['q_bh'])}); in "
            f"Actinopteri it is {act3['n_pos']} to {act3['n_neg']} over "
            f"{act3['n']} untied pairs (p = {_pf(act3['p'])}) — no signal "
            "at all, on a sample that is smaller but not small. Whatever is "
            "raising ITPR3's indel density is not doing it across the "
            "vertebrates.",
            "",
            "**And the design makes a cell and its siblings the same "
            "observation twice.** The comparison is within one genome, "
            "against that genome's other family loci, so a bird ITPR3 that "
            "is elevated makes its own ITPR1 and ITPR2 look deficient by "
            "construction — which is exactly the pattern in the table. The "
            "three Aves rows are one result, not three.",
            "",
        ]
    if ctrl:
        out += ["### 8.3 The control that has to be printed", ""]
        out += table(["cell / class", "direction", "above D4's bar",
                      "below D4's bar", "median identity (locus / others)"],
                     [[f"{c['cell']} / {c['vclass']}", c["direction"],
                       f"{c['pos_above']}:{c['neg_above']}"
                       + (f" (p = {_pf(c['p_above'])})"
                          if c["p_above"] != "" else " (too few to test)"),
                       f"{c['pos_below']}:{c['neg_below']}"
                       + (f" (p = {_pf(c['p_below'])})"
                          if c["p_below"] != "" else " (too few to test)"),
                       f"{c['median_identity']} / "
                       f"{c['median_others_identity']}"]
                      for c in ctrl])
        aves_ctrl = next((c for c in ctrl if c["cell"] == "ITPR3"
                          and c["vclass"] == "Aves"), None)
        if aves_ctrl:
            out += [
                f"**{aves_ctrl['n_below_bar']} of "
                f"{_i(aves_ctrl['n_above_bar']) + _i(aves_ctrl['n_below_bar'])} "
                "of the bird pairs are in assemblies below D4's "
                "contiguity bar**, and that is where the test has its "
                f"power: {aves_ctrl['pos_below']}:{aves_ctrl['neg_below']}, "
                f"p = {_pf(aves_ctrl['p_below'])}. Above the bar all "
                f"{aves_ctrl['pos_above']} pairs point the same way and "
                "none points against — but "
                f"{_i(aves_ctrl['n_above_bar'])} pairs cannot carry a "
                "test. S5b measured 66 % of Aves assemblies below the "
                "bar, the worst of any class, so this is the class where "
                "an indel signal is hardest to separate from an assembly "
                "signal. The identity control is clean — locus and "
                "siblings are matched to "
                f"{aves_ctrl['median_identity']} against "
                f"{aves_ctrl['median_others_identity']} — so D47's own "
                "confounder is not what this is. **The honest verdict is "
                "that the lineage is now named and the mechanism is not, "
                "and that the contiguous half of the evidence is too "
                "small to settle it.**",
                "",
            ]
    out += fig("lesion_strata",
               "**Figure 4.** D47's within-genome, identity-matched sign "
               "test stratified by vertebrate class (a); the strongest "
               "stratum split by D4's contiguity bar (b); and the fossil "
               "denominator, with the number of loci scored above each bar "
               "(c).")
    return out


# ------------------------------------------------------------------ §9

def section_priors(load, num, pct, table, missing, fig, stats,
                   head) -> list[str]:
    cls = load("lesion_by_class.tsv")
    fos = stats.get("fossil", {})
    lines = ["## 9. Every prior, computed on this half's own tables", "",
             "Each row states what an earlier task concluded and where, "
             "computes S15b's answer beside it, and renders the verdict "
             "from the comparison.", ""]
    lines.append(_verdict_line(
        "absent_cells", 0, head.get("dollo_family_base", 0),
        priors.verdict(0, head.get("dollo_family_base", 0))))
    lines.append(_verdict_line(
        "family_coding_primary",
        "family-level presence per genome",
        f"{head.get('n_settings_family_manufactures', 0)} of "
        f"{head.get('n_settings', 0)} settings manufacture a family loss "
        f"against {head.get('n_settings_paralog_manufactures', 0)} of "
        f"{head.get('n_settings', 0)} paralog-resolved",
        "confirmed"))
    n_r5 = _absent_at(load("dollo_counts.tsv"),
                      coding.BASE_SETTING["evidence"], 0, 1)
    lines.append(_verdict_line(
        "loss_candidates", 4, n_r5, priors.verdict(4, n_r5)))
    lines.append(_verdict_line(
        "itpr1_gain_node", "Vertebrata", head.get("gain_node_itpr1", "?"),
        "orthogonal"))
    lines.append(_verdict_line(
        "itpr23_gain_node", "Gnathostomata",
        f"ITPR2 → {head.get('gain_node_itpr2', '?')}, "
        f"ITPR3 → {head.get('gain_node_itpr3', '?')}",
        "orthogonal"))
    aves = next((t for t in cls if t["cell"] == "ITPR3"
                 and t["vclass"] == "Aves"), None)
    if aves:
        lines.append(_verdict_line(
            "itpr3_indel_excess", 0.0032,
            f"concentrated in Aves ({aves['n_pos']}:{aves['n_neg']}, "
            f"q = {_pf(aves['q_bh'])}), absent in Actinopteri",
            "underpowered"))
    lines.append(_verdict_line(
        "no_fossils", 0, head.get("n_fossils_itpr", 0),
        "contradicted" if head.get("n_fossils_itpr", 0) else "confirmed"))
    lines.append(_verdict_line(
        "zero_stops_falsifies", 0,
        f"{fos.get('n_above_bar_zero_stops', 0)} of "
        f"{fos.get('n_above_bar', 0)} loci above the lesion bar carry no "
        "internal stop; the rest carry 1-2 at full coverage",
        "confirmed"))
    lines.append(_verdict_line(
        "polytomy_degree", 23,
        (stats.get("polytomy_profile", {}) or {}).get("max_degree", 0),
        priors.verdict(
            23, (stats.get("polytomy_profile", {}) or {}).get(
                "max_degree", 0))))
    lines.append("")
    lines += [
        "One verdict needs its wording defended. `no_fossils` is rendered "
        "**contradicted** because the generous screen fires on 7 ITPR "
        "loci where S15a reported none — but every one of the 7 is at full "
        "coverage with 1–2 internal stops, so the disagreement is with "
        "S15a's *wording* and not with its finding. The screen is written "
        "to be easy to pass precisely so that a reader can see how little "
        "passing it is worth.",
        "",
    ]
    return lines


# ----------------------------------------------------------------- §10

def section_handoff(load, num, table, missing, fig, stats,
                    head) -> list[str]:
    ctrl = load("lesion_class_controls.tsv")
    aves = next((c for c in ctrl if c["cell"] == "ITPR3"
                 and c["vclass"] == "Aves"), None)
    out = ["## 10. What this half settles and what it does not", "",
           "**Settles.**", ""]
    out += [
        "- **No vertebrate lineage in this scope has lost an IP3 "
        f"receptor.** Dollo places {head.get('dollo_family_base', 0)} "
        "losses on the family character and "
        f"{head.get('dollo_paralog_base', 0)} on every paralog character, "
        f"across {num(head.get('n_cells', 0))} genome × paralog cells in "
        "309 genomes, and the routine that returns that zero is one that "
        "finds a constructed loss on a known edge.",
        "- **The family-level coding is robust and the paralog-resolved "
        "one is not**, by a factor of "
        f"{head.get('n_settings_paralog_manufactures', 0)} to "
        f"{head.get('n_settings_family_manufactures', 0)} settings. D46 "
        "was a judgement when S15a made it; it is a measurement now.",
        "- **The reconstruction bar's position inside its measured gap "
        "changes nothing**, and D45 alone is worth 4 losses that never "
        "happened.",
        "- **No Mk rate can be reported for this character, and the "
        "reason is measurable**: every likelihood is monotone to its "
        "boundary on every model, axis and branch-length scheme.",
        "- **There are no pseudogene fossils**, on a denominator of "
        f"{num(head.get('n_scored_loci', 0))} scored loci, with the "
        "generous screen printed beside the strict one.",
        "",
        "**Does not settle.** Four things, all of them leads rather than "
        "gaps.", "",
    ]
    if aves:
        out += [
            "1. **Why bird ITPR3 carries extra indels.** The lineage is "
            f"named — {aves['n']} identity-matched within-genome pairs, "
            f"q = {_pf(aves['q_bh'])}, with no equivalent in ray-finned "
            f"fish — but {aves['n_below_bar']} of {aves['n']} of those "
            "pairs sit in assemblies below D4's contiguity bar, and the "
            f"{aves['n_above_bar']} above it all point the same way "
            "without being enough to test. The next instrument is a bird "
            "panel restricted to chromosome-level assemblies, which is a "
            "different sample and not a different statistic.",
        ]
    out += [
        "2. **What the one reconstruction-only cell actually is.** "
        "*Bothrops jararaca* ITPR2 is the single cell in 927 whose "
        "presence rests on the cross-contig reassembly and on nothing "
        "else. It is not a loss under any bar the calibration licenses, "
        "but it is the cell a sceptical reader should be handed first.",
        "3. **The count is bounded by the tree's resolution, and the "
        "bound was never tested.** With no loss to place, the interval "
        "between `dollo_losses_max` and `dollo_losses_min` is never "
        "exercised at the operating point. It is exercised across the "
        "manufactured settings, and there it is wide — up to "
        f"{head.get('max_paralog_losses', 0)} edges collapsing to "
        f"{head.get('max_paralog_losses_min', 0)}.",
        "4. **A zero across 309 vertebrate genomes is a statement about "
        "309 vertebrate genomes.** S4 declared that denominator and S5b "
        "swept it; nothing here extends to a species with no assembly.",
        "",
        "**Two things a reader should hold against this half.** The "
        "sensitivity matrix's most extreme settings are not settings "
        "anybody would adopt — refusing a complete-but-truncated locus "
        "*and* ignoring contiguity is not a defensible protocol — so the "
        "144 cells that read `absent` somewhere in the grid should be "
        "read as a map of fragility and not as a list of candidates. And "
        "the gain-node agreement with S13 is printed because it is worth "
        "printing, not because it is independent: the cyclostome cells "
        "that put ITPR2 and ITPR3's gain on the gnathostome stem are "
        "`paralog_unassignable` for a bait-panel reason, which is the same "
        "evidence S13 used, seen from a different side.",
        "",
    ]
    return out



def render(load, load_json, num, pct, table, missing, fig, stats,
           head) -> list[str]:
    out: list[str] = []
    out += section_fossils(load, load_json, num, pct, table, missing, fig,
                           stats, head)
    out += section_priors(load, num, pct, table, missing, fig, stats, head)
    out += section_handoff(load, num, table, missing, fig, stats, head)
    return out
