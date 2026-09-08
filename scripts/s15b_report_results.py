"""S15b — the counts half of `results/loss_counts/report.md`.

Split from `s15b_report.py`, and split again from
`s15b_report_lesions.py`, so all three stay under 500 lines and take the
caller's loader and formatter — the `s3_report.py` / `s3_report_d10.py`
pattern applied twice, as S7, S9 and S12 do.

This half carries the count itself (§5), the sensitivity matrix that is
this task's deliverable (§6) and the Mk section's refusal (§7).  The half
that argues against the result — the fossil denominator, D47's stratified
lead with its contiguity control, the priors table and the hand-off —
lives in `s15b_report_lesions.py`.
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


# ------------------------------------------------------------------ §5

def section_count(load, num, pct, table, missing, fig, head) -> list[str]:
    dc = load("dollo_counts.tsv")
    out = ["## 5. The count", ""]
    if not dc:
        return out + missing("dollo_counts.tsv", "the Dollo count")
    base = [r for r in dc if _i(r["is_base"])]
    out += [
        "Under S15a's own setting, on both codings:",
        "",
    ]
    out += table(
        ["coding", "character", "present", "absent", "undecided",
         "gain placed at", "Dollo losses (max / min)"],
        [[r["coding"], f"`{r['character']}`", num(r["n_present"]),
          num(r["n_absent"]), num(r["n_undecided"]), r["gain_node"],
          f"{r['dollo_losses_max']} / {r['dollo_losses_min']}"]
         for r in base])
    out += [
        "**No vertebrate lineage in this scope has lost an IP3 receptor.** "
        "The primary, family-level count is "
        f"{head.get('dollo_family_base', 0)}, and so is every "
        "paralog-resolved one. The count is zero because no cell reaches "
        "the state that would license it, not because the routine cannot "
        "return anything else: `s15b_test_counts.py` T1 puts a loss on a "
        "known edge and requires it found, T2 requires two sister losses "
        "merged into their parent, and T8 requires a constructed "
        "contiguous, controlled, spare-free empty cell to come back "
        "`absent`.",
        "",
        "### 5.1 The gain nodes, which nobody asked this instrument for",
        "",
        "Dollo places the single gain at the MRCA of the tips that carry "
        "the character, and for the paralog characters that node is not "
        "pinned. The three answers are "
        f"**ITPR1 → {head.get('gain_node_itpr1', '?')}**, "
        f"**ITPR2 → {head.get('gain_node_itpr2', '?')}**, "
        f"**ITPR3 → {head.get('gain_node_itpr3', '?')}** — which is "
        "exactly where S13 placed the two duplications, recovered here by "
        "a method that reads no gene tree, no alignment and no "
        "reconciliation. It is a consistency check and not a second "
        "result: the reason ITPR2 and ITPR3 have no cyclostome tip is "
        "that S15a coded those cells `paralog_unassignable` for want of a "
        "cyclostome-labelled bait, and the reason S13 placed the "
        "duplication there is a reconciliation. The agreement is worth "
        "printing; it is not independent evidence.",
        "",
    ]
    return out


# ------------------------------------------------------------------ §6

def section_sensitivity(load, num, pct, table, missing, fig,
                        head) -> list[str]:
    dc = load("dollo_counts.tsv")
    ml = load("manufactured_losses.tsv")
    out = ["## 6. The sensitivity matrix — the deliverable", ""]
    if not dc:
        return out + missing("dollo_counts.tsv", "the sensitivity matrix")
    out += [
        "With no loss to place, what is worth reporting is the shape of "
        "the zero: how far the rules must move before a loss appears, "
        "which rule has to move, and what each move buys. "
        f"{head.get('n_settings', 0)} settings were walked "
        f"({len(coding.EVIDENCE_NAMES)} evidence rungs × D45 on/off × D4 "
        "on/off), on both codings, under all three branch-length schemes.",
        "",
    ]
    out += fig("sensitivity_matrix",
               "**Figure 1.** Dollo losses manufactured by every setting of "
               "the rule chain, for the primary family-level coding (a) and "
               "the paralog-resolved one (b). The violet ring marks S15a's "
               "own operating point. A zero is drawn as an explicit zero "
               "and never as an empty cell, because an empty cell reads as "
               "*not measured*.")
    fam_rows, par_rows = [], []
    for ev in coding.EVIDENCE_NAMES:
        line_f, line_p = [f"`{ev}`"], [f"`{ev}`"]
        for r5, r6 in ((1, 1), (0, 1), (1, 0), (0, 0)):
            g = [r for r in dc if r["evidence"] == ev
                 and _i(r["use_r5"]) == r5 and _i(r["use_r6"]) == r6]
            line_f.append(str(max((_i(r["dollo_losses_max"]) for r in g
                                   if r["coding"] == "family"), default=0)))
            par = [r for r in g if r["coding"] == "paralog"]
            mx = max((_i(r["dollo_losses_max"]) for r in par), default=0)
            mn = max((_i(r["dollo_losses_min"]) for r in par), default=0)
            line_p.append(f"{mx} / {mn}" if mx else "0")
        fam_rows.append(line_f)
        par_rows.append(line_p)
    cols = ["evidence rung", "both on (S15a)", "D45 off", "D4 off",
            "both off"]
    out += ["**Family-level coding (primary, D46).**", ""]
    out += table(cols, fam_rows)
    out += ["**Paralog-resolved coding (max / min independent, worst of "
            "ITPR1/2/3).**", ""]
    out += table(cols, par_rows)
    out += [
        "Read across the two tables, the result is an asymmetry rather "
        "than a number. The **primary coding manufactures a loss in "
        f"{head.get('n_settings_family_manufactures', 0)} of "
        f"{head.get('n_settings', 0)} settings**, and its worst case is "
        f"{head.get('max_family_losses', 0)} genome in 309 — reached only "
        "by refusing everything except a complete locus *and* ignoring "
        "D4's contiguity bar at the same time. The **paralog-resolved "
        "coding manufactures one in "
        f"{head.get('n_settings_paralog_manufactures', 0)} of "
        f"{head.get('n_settings', 0)}**, up to "
        f"{head.get('max_paralog_losses', 0)} loss edges "
        f"({head.get('max_paralog_losses_min', 0)} independent under the "
        "tree's own resolution). That is D46 measured rather than "
        "asserted: a per-paralog absence is fragile to every one of these "
        "knobs and a family-level absence is not.",
        "",
        "### 6.1 The three rungs that change nothing", "",
        "Moving the reconstruction bar across the **whole gap the "
        f"calibration measured** — `{coding.BAR_GAP_LO}` (the decoy's "
        f"maximum) to `{coding.BAR_GAP_HI}` (the lowest candidate) — "
        "manufactures no loss on either coding. The bar's position inside "
        "its own uncertainty is not what any result here rests on.",
        "",
        "### 6.2 What each knob is worth", "",
    ]
    if ml:
        by_release: dict[str, list[dict]] = {}
        for r in ml:
            by_release.setdefault(r["released_by"], []).append(r)
        rows = []
        for key, g in sorted(by_release.items(),
                             key=lambda kv: -len(kv[1]))[:8]:
            rows.append([key, num(len(g)),
                         ", ".join(sorted({x["cell"] for x in g})),
                         ", ".join(sorted({x["vclass"] for x in g}))[:60]])
        out += table(["the loosest setting that releases the cell",
                      "cells", "paralog cells", "classes"], rows)
        base_ev = coding.BASE_SETTING["evidence"]
        n_r5_only = _absent_at(dc, base_ev, 0, 1)
        no_knob = [r for r in ml if not _i(r["needs_r5_off"])
                   and not _i(r["needs_r6_off"])]
        out += [
            f"**{num(head.get('n_cells_ever_absent', 0))} of "
            f"{num(head.get('n_cells', 0))} cells read `absent` under at "
            "least one of the 32 settings**, and every one of them names "
            "in its row what held it at S15a's operating point.",
            "",
            "Turning **D45 off alone** — refusing to treat a genome's "
            "spare, unassignable family loci as evidence — manufactures "
            f"{n_r5_only} losses immediately, at every evidence rung, "
            "and they are the cyclostome ITPR2/ITPR3 cells in *Petromyzon "
            "marinus* and *Myxine glutinosa*. Both genomes carry three "
            "ITPR loci apiece. That is what D45 is worth.",
            "",
        ]
        if no_knob:
            top = no_knob[0]
            out += [
                "**Tightening the evidence alone**, with both protective "
                f"rules left on, manufactures {len(no_knob)} loss in the "
                f"whole scope: *{top['organism']}* {top['cell']} "
                f"({top['accession']}), whose reference reassembles at "
                f"{_f(top['recon_coverage']):.3f} across "
                f"{top['recon_n_contigs']} contigs in an assembly "
                "contiguous enough to carry the gene, with no spare family "
                "locus to explain it. It reads `absent` only if a "
                "cross-contig reassembly at "
                f"{_f(top['recon_coverage']):.0%} of the reference is "
                "refused as evidence — six times the bar the calibration "
                "measured. It is not a loss; it is the single cell in 927 "
                "whose presence rests on the reconstruction instrument "
                "alone.",
                "",
            ]
    out += fig("reconstruction_bar",
               "**Figure 2.** The 43 cells R4 places, by what would catch "
               "each if the reconstruction bar rose (a), and the losses "
               "that rise manufactures under each rule setting (b). The "
               "calibrated bar and both edges of its measured gap are "
               "drawn, not stated.")
    return out


# ------------------------------------------------------------------ §7

def section_mk(load, load_json, num, pct, table, missing, fig,
               stats, head) -> list[str]:
    prof = load("mk_profile.tsv")
    fits = load("mk_fits.tsv")
    out = ["## 7. Mk models — stated, not fitted", ""]
    if not prof:
        return out + missing("mk_profile.tsv", "the Mk section")
    out += [
        "The brief asks for ER, ARD and an irreversible model with the "
        "gain rate pinned to zero. On the primary character all three are "
        "**refused**, and the refusal is measured rather than asserted.",
        "",
        "An invariant character contains no transition to estimate. "
        "Profiling each likelihood along a log-spaced rate grid over "
        "eight orders of magnitude gives, at every one of the "
        f"{len(set((r['branch_lengths'], r['model'], r.get('axis', '')) for r in prof))} "
        "model × branch-length × axis combinations, a monotone curve with "
        "its maximum on the grid's boundary:",
        "",
    ]
    seen, rows = set(), []
    for r in prof:
        key = (r["model"], r.get("axis", "loss"), r["branch_lengths"])
        if key in seen:
            continue
        seen.add(key)
        rows.append([f"`{r['model']}`", r.get("axis", "loss"),
                     r["branch_lengths"], r["shape"],
                     f"{_f(r['argmax']):.3g}",
                     "yes" if _i(r["at_boundary"]) else "no"])
    out += table(["model", "axis profiled", "branch lengths", "shape",
                  "arg max", "at the boundary"], rows)
    out += [
        "The loss axis falls and ARD's **gain** axis rises, to the edge of "
        "the grid in both directions — which is what unidentifiability "
        "looks like when it is drawn rather than argued. ARD is profiled "
        "on its gain axis and not on its diagonal, because the diagonal is "
        "ER by construction and would put the same curve on the figure "
        "twice under two names.",
        "",
        "So no rate is reported for the primary character. A fitter run "
        "on it would return its own starting point, and all "
        f"{head.get('n_mk_rows_base', 0) - head.get('n_mk_fitted_base', 0)} "
        f"of {head.get('n_mk_rows_base', 0)} model fits at S15a's "
        "operating point — every character, under every branch-length "
        "scheme — are refusals for exactly that reason.",
        "",
        "### 7.1 The informative version", "",
        "Where the sensitivity matrix *does* produce variation the fits "
        "are real, and they are reported for what they are: "
        f"{head.get('n_mk_fits_total', 0)} fits against "
        f"{head.get('n_mk_refusals_total', 0)} refusals, one row per "
        "*distinct* character (the fits are memoised on the character "
        "vector, so a rung that changes no cell costs nothing). "
        "The rate tracks the number of losses the setting manufactured "
        "and the units the branch lengths are in, which is the whole "
        "point — it is a property of the filter, not of the family.",
        "",
    ]
    fitted = [r for r in fits if _i(r["fitted"])]
    if fitted:
        worst = max(fitted, key=lambda r: _i(r["n_absent"]))
        same = [r for r in fitted
                if r["character"] == worst["character"]
                and r["evidence"] == worst["evidence"]
                and r["use_r5"] == worst["use_r5"]
                and r["use_r6"] == worst["use_r6"]
                and r["model"] == "irreversible"]
        rows = [[r["branch_lengths"], f"`{r['model']}`", num(r["n_absent"]),
                 f"{_f(r['rate']):.4g}", f"{_f(r['loglik']):.2f}",
                 f"{_f(r['aic']):.2f}"]
                for r in sorted(same, key=lambda r: r["branch_lengths"])]
        if rows:
            out += [
                f"The most extreme cell of the matrix — `{worst['coding']}` "
                f"coding, character `{worst['character']}`, evidence "
                f"`{worst['evidence']}`, both protective rules off, "
                f"{worst['n_absent']} cells made absent — under the "
                "irreversible model:",
                "",
            ]
            out += table(["branch lengths", "model", "cells absent",
                          "loss rate", "lnL", "AIC"], rows)
            rates = [_f(r["rate"]) for r in same if _f(r["rate"]) > 0]
            spread = (max(rates) / min(rates)) if rates else 0.0
            out += [
                f"The three rates span a factor of {spread:,.0f} across the "
                "three schemes while describing the same character. That "
                "is the branch-length axis doing the only thing it can do "
                "in this task: set the units a rate is quoted in.",
                "",
            ]
    bl = head.get("branch_length_changes_dollo", 0)
    out += [
        "### 7.2 The branch-length axis, reported rather than omitted", "",
        f"Across all {head.get('n_settings', 0)} settings and all three "
        "schemes, the number of settings at which a branch-length scheme "
        f"changes a Dollo count is **{bl}**. Parsimony counts edges and "
        "does not read a length, so this could not have come out any "
        "other way — but the brief names branch lengths as an axis, and a "
        "matrix that quietly dropped one axis would be "
        "indistinguishable from one that had tested it.",
        "",
    ]
    out += fig("mk_profile",
               "**Figure 3.** The primary character's likelihood profile "
               "under all three models and all three branch-length schemes "
               "(a) — monotone to the boundary, which is the argument for "
               "refusing to fit it — and every rate that *can* be fitted, "
               "against the number of losses its setting manufactured (b).")
    return out



def render(load, load_json, num, pct, table, missing, fig, stats,
           head) -> list[str]:
    """The counts half, then the half that argues against it."""
    import s15b_report_lesions as lesions
    out: list[str] = []
    out += section_count(load, num, pct, table, missing, fig, head)
    out += section_sensitivity(load, num, pct, table, missing, fig, head)
    out += section_mk(load, load_json, num, pct, table, missing, fig,
                      stats, head)
    out += lesions.render(load, load_json, num, pct, table, missing, fig,
                          stats, head)
    return out
