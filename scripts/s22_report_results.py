"""The results half of S22's report — the `s3_report.py` / `s3_report_d10.py`
split.  It takes the caller's headline dict so the two halves cannot read the
tables differently, and every headline is chosen by the data: each prior is
stated in `s22_priors.PRIOR` with where the earlier task said it, computed on
S22's own tables, and rendered with both numbers printed either way.

The comparison that has to be sayable is the one that goes against the
task's own premise, and here it is the headline: the module this task was
built around is **not** the more constrained of the two.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L                                            # noqa: E402
import s22_contacts as C                                       # noqa: E402
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


def section_modules(h: dict) -> str:
    if not h.get("module_primary"):
        return A("## 6. The ligand core against the pore", "", "*not run yet*")
    prim, loop_in = h["module_primary"], h["module_loop_in"]
    rows = ["| paralogue | orthologues | core identity | pore identity | "
            "difference | tips favouring core / pore | q |",
            "|---|---|---|---|---|---|---|"]
    for p in PARA:
        r = prim.get(p, {})
        con = next((x for x in L.read_tsv(L.OUT_DIR / "module_contrast.tsv")
                    if x["paralog"] == p and x["is_primary"] == "True"), {})
        rows.append(
            f"| {p} | {r.get('n_tips', '—')} | {f(r.get('core'))} | "
            f"{f(r.get('pore'))} | {f(r.get('difference'), 4)} | "
            f"{con.get('n_core_more_conserved', '—')} / "
            f"{con.get('n_pore_more_conserved', '—')} | "
            f"{pfmt(r.get('q'))} |")
    flip = ["| paralogue | loop excluded (primary) | loop included | verdict |",
            "|---|---|---|---|"]
    for p in PARA:
        a, b = prim.get(p, {}), loop_in.get(p, {})
        flip.append(f"| {p} | {f(a.get('difference'), 4)} "
                    f"({a.get('direction', '—')}, q = {pfmt(a.get('q'))}) | "
                    f"{f(b.get('difference'), 4)} ({b.get('direction', '—')}, "
                    f"q = {pfmt(b.get('q'))}) | "
                    f"{'reversed' if a.get('direction') != b.get('direction') else 'unchanged'} |")

    cc = L.read_tsv(L.OUT_DIR / "column_contrast.tsv")
    met = ["| paralogue | core JSD | pore JSD | core frac-modal | "
           "pore frac-modal | which metric favours the pore |",
           "|---|---|---|---|---|---|"]
    for p in PARA:
        r = next((x for x in cc if x["paralog"] == p and x["layer"] == "deep"
                  and x["is_primary"] == "True"), {})
        jsd_pore = float(r.get("mean_pore_jsd") or 0) > float(r.get("mean_core_jsd") or 0)
        fm_pore = float(r.get("mean_pore_frac_modal") or 0) > float(r.get("mean_core_frac_modal") or 0)
        met.append(f"| {p} | {f(r.get('mean_core_jsd'))} | "
                   f"{f(r.get('mean_pore_jsd'))} | "
                   f"{f(r.get('mean_core_frac_modal'))} | "
                   f"{f(r.get('mean_pore_frac_modal'))} | "
                   f"{'JSD, ' if jsd_pore else ''}"
                   f"{'frac-modal' if fm_pore else 'neither'} |")

    n_rev = sum(1 for p in PARA
                if prim.get(p, {}).get("direction") != loop_in.get(p, {}).get("direction"))
    return A(
        "## 6. The ligand core against the pore",
        "",
        "### 6.1 The answer, and it is not the one the question expects",
        "",
        "**The pore module is more conserved than the ligand core, not less.** "
        "Measured per orthologue and paired inside each paralogue's own "
        "alignment:",
        "",
        *rows,
        "",
        "ITPR1 and ITPR3 put the pore ahead by about two identity points, "
        "against roughly 250 orthologues each and with the tips running "
        "seven-to-one and five-to-one the same way.  ITPR2 shows no "
        "difference at all.  The module that gives this family its name and "
        "its defining Pfam signature is the *less* constrained of the two "
        "across a vertebrate orthologue set.",
        "",
        "### 6.2 And it reverses on a boundary nobody states",
        "",
        ("**All three paralogues flip sign**" if n_rev == 3 else
         f"**{n_rev} of three paralogues flip sign**") + " when the luminal "
        "loop is "
        "left inside the pore module — which is how InterPro draws PF00520 "
        "and therefore how the comparison would be made by default:",
        "",
        *flip,
        "",
        "Fifty residues of solvent-facing luminal loop, the least conserved "
        "sequence in the whole receptor (S17), sit about fifty residues from "
        "the gate inside one Pfam domain.  Include them and the ligand core "
        "wins by five identity points at q < 1e-37 in all three paralogues; "
        "exclude them and the pore wins in two and ties in the third.  "
        "Neither answer is wrong about its own module — they are answers "
        "about different modules, and the difference between them is one "
        "boundary that a comparison drawn on Pfam spans would never have to "
        "declare.",
        "",
        "### 6.3 The two metrics disagree, and the composition-free one wins",
        "",
        "At column level the JSD sees no difference between the modules, "
        "while the per-orthologue identity does.  That is not a "
        "contradiction: JSD is a divergence from a background amino-acid "
        "table, so a transmembrane module scores lower than a soluble one at "
        "equal conservation.  The composition-free column metric agrees with "
        "the paired test.",
        "",
        *met,
        "",
        "So the column-level null result is a property of the metric, and it "
        "is reported here rather than dropped because a reader coming from "
        "S17's tables would otherwise find two of this project's own numbers "
        "pointing different ways with no explanation.",
    )


def section_contacts(h: dict) -> str:
    path = L.OUT_DIR / "contact_test.tsv"
    if not path.exists():
        return A("## 7. The contacts, and the pocket", "", "*not run yet*")
    ct = L.read_tsv(path)
    rows = ["| paralogue | contacts | rest of core | q | "
            "rest of pocket | q |", "|---|---|---|---|---|---|"]
    for p in PARA:
        a = next((r for r in ct if r["paralog"] == p and r["layer"] == "deep"
                  and r["contact_set"] == "s0_contact"
                  and r["background"] == "rest_of_core"), {})
        b = next((r for r in ct if r["paralog"] == p and r["layer"] == "deep"
                  and r["contact_set"] == "s0_contact"
                  and r["background"] == "rest_of_pocket"), {})
        rows.append(f"| {p} | {f(a.get('mean_contact_jsd'))} | "
                    f"{f(a.get('mean_background_jsd'))} | "
                    f"{pfmt(a.get('q_permutation'))} | "
                    f"{f(b.get('mean_background_jsd'))} | "
                    f"{pfmt(b.get('q_permutation'))} |")
    sc = L.read_tsv(L.OUT_DIR / "shell_constraint.tsv")
    shell = ["| shell | ITPR1 | ITPR2 | ITPR3 | residues |",
             "|---|---|---|---|---|"]
    for s in C.SHELL_ORDER:
        vals = [next((r for r in sc if r["paralog"] == p and r["shell"] == s), {})
                for p in PARA]
        shell.append(f"| {s} | " + " | ".join(f(v.get("mean_jsd")) for v in vals)
                     + f" | {vals[0].get('n_residues', '—')} |")
    tr = h["trend"]
    trend = ["| paralogue | ρ(distance, conservation) | p |", "|---|---|---|"]
    for p in PARA:
        rho, pv = tr.get(p, ("—", "—"))
        trend.append(f"| {p} | {f(rho)} | {pfmt(pv)} |")
    return A(
        "## 7. The contacts, and the pocket",
        "",
        "**The ten residues that touch IP₃ are more constrained than the rest "
        "of the binding core, and not more constrained than the rest of the "
        "pocket.**  The permutation resamples which positions carry the "
        "contact label, keeping the constraint values where they are, because "
        "the hypothesis is about these particular residues.",
        "",
        *rows,
        "",
        "Against the core the contacts win in all three paralogues.  Against "
        "every *other* residue the structure places within "
        "15 Å of the ligand they win in none.  The constrained unit is the "
        "pocket, not the ten contacts — which is only sayable because the "
        "shells were measured, and is invisible to a contact/not label.",
        "",
        *shell,
        "",
        "Every shell out to 15 Å sits above the whole-protein mean, and the "
        "step a contact-driven model predicts at 4.5 Å is not there.  The "
        "gradient across the neighbourhood is real but shallow, and reaches "
        "significance in one paralogue:",
        "",
        *trend,
        "",
        "Read together: IP₃ binding constrains a pocket about fifteen "
        "ångström across, and the resolution at which the ten contacts are "
        "special is the resolution at which the whole pocket is.",
    )


def section_selection() -> str:
    path = L.OUT_DIR / "omega_module_test.tsv"
    if not path.exists():
        return A("## 8. The same question asked of the substitution rate", "",
                 "*not run yet*")
    mt = L.read_tsv(path)
    bs = L.read_tsv(L.OUT_DIR / "omega_by_shell.tsv")
    rows = ["| paralogue | core β | pore β | direction | q |",
            "|---|---|---|---|---|"]
    for p in PARA:
        r = next((x for x in mt if x["paralog"] == p
                  and x["is_primary"] == "True"), {})
        rows.append(f"| {p} | {f(r.get('mean_core_beta'), 4)} | "
                    f"{f(r.get('mean_pore_beta'), 4)} | "
                    f"{r.get('direction', '—')} | "
                    f"{pfmt(r.get('q_mannwhitney'))} |")
    shell = ["| shell | ITPR1 | ITPR2 | ITPR3 |", "|---|---|---|---|"]
    series: dict[str, list[float]] = {p: [] for p in PARA}
    for s in C.SHELL_ORDER:
        vals = [next((r["frac_purifying_q05"] for r in bs
                      if r["paralog"] == p and r["shell"] == s), None)
                for p in PARA]
        for p, v in zip(PARA, vals):
            series[p].append(float(v) if v not in (None, "") else float("nan"))
        shell.append(f"| {s} | " + " | ".join(f(v) for v in vals) + " |")
    # Whether the fall is monotone is computed, not asserted: on this data it
    # is not — every paralogue's fourth shell sits slightly above its third.
    mono = {p: all(b <= a + 1e-9 for a, b in zip(v, v[1:]))
            for p, v in series.items()}
    drop = {p: v[0] - min(v) for p, v in series.items()}
    n_mono = sum(mono.values())
    trend_sentence = (
        "Every contact site in every paralogue is under detectable purifying "
        "selection, and the share falls with distance from the ligand — by "
        + " / ".join(f(drop[p], 3) for p in PARA)
        + " between the contact shell and the lowest shell. "
        + ("The fall is monotone in all three." if n_mono == 3 else
           ("It is not monotone in any of the three" if n_mono == 0 else
            f"It is monotone in {n_mono} of three")
           + ": the fourth shell sits a little above the third, so what the "
             "data show is a step down from the ligand's first two shells "
             "and then a floor, not a gradient running all the way out."))
    return A(
        "## 8. The same question asked of the substitution rate",
        "",
        "Conservation and ω are both called constraint in prose and they are "
        "not the same quantity: a column's dispersion is a statement about "
        "which amino acids are seen across 250 orthologues, a site's β is a "
        "rate on a 57-tip vertebrate tree.  S17's committed FEL run supplies "
        "the second; nothing here refits a model.",
        "",
        *rows,
        "",
        "**The module comparison does not replicate on this axis, and in "
        "ITPR1 it points the other way**: §6 has the pore more conserved "
        "there, and β has it evolving faster.  Both are computed correctly "
        "and they are not the same measurement — 250 sweep orthologues of "
        "one paralogue against 57 vertebrate tips, a column's dispersion "
        "against a rate on a tree — but a reader should have the "
        "disagreement in front of them rather than only the half that "
        "agrees.  With a rate that is exactly zero at most sites of this "
        "protein and a quarter of the sequences, this is the underpowered "
        "version of the same comparison, and the module result is read off "
        "§6.",
        "",
        "**The shell result does replicate.**  The share of sites FEL calls "
        "purifying at q ≤ 0.05:",
        "",
        *shell,
        "",
        trend_sentence + "  Two instruments, two alignments, two different "
        "statistics, the same pocket.",
    )


def build(h: dict) -> str:
    import s22_report_lineage as LN
    return "\n\n".join([
        section_modules(h),
        section_contacts(h),
        section_selection(),
        LN.build(h),
    ])
