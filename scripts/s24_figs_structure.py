"""Supplementary Figures 5-6: the constraint map on the channel, and the
labelled variants with their enrichment test.

Three decisions worth naming.

* **The trace is drawn from the file S17 painted, and coloured from its
  B-factor column.** Nothing here re-computes a conservation value or
  re-maps a residue: if the picture disagreed with `constraint_*.tsv` it
  would be a bug in S17's painting step, which is the point of drawing it.
  Unscored residues are **grey, never the low end of the scale** — S17
  wrote them as −1 precisely so that "we could not score this" and "this is
  the least conserved part of the receptor" cannot look the same, and the
  luminal loop is both.

* **A structure may carry human variant positions only if it passes an
  identity test, and two of the four candidates fail.** A ClinVar position
  is a number in human canonical numbering; a deposition numbers its own
  construct. `numbering_check` requires *every* residue the structure shares
  with the human per-residue table to carry the same amino acid. Human
  *ITPR2* (9YKK) and *ITPR3* (8TKG) pass at 2,168/2,168 and 2,210/2,210.
  The *ITPR1* reference is a rat structure (736/2,300) and the AlphaFold DB
  model of human *ITPR1* is the 2,695-residue Q14643-4 isoform (479/2,695),
  so neither carries *ITPR1*'s 55 pathogenic positions and the figure says
  so rather than drawing them on coordinates that are not theirs.

* **The enrichment panel plots the odds ratio with its Fisher p, not the
  raw count.** Eleven pathogenic positions in RIH_N is a statement about
  RIH_N only against the share of the protein RIH_N is, and the elements
  here run from 2 residues to 618.
"""

from __future__ import annotations

import math
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.lines import Line2D                            # noqa: E402

import figstyle as F                                           # noqa: E402
import s24_lib as L                                            # noqa: E402

PAINTED = L.CONSTRAINT_DIR / "painted"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

#: The cryo-EM reference S11 chose for each paralogue, and the human model
#: AlphaFold DB serves for it. Both are candidates for carrying variant
#: positions and both are put through the same test.
REFERENCE_FILE = {
    "ITPR1": "constraint_reference_ITPR1_7LHF.pdb",
    "ITPR2": "constraint_reference_ITPR2_9YKK.pdb",
    "ITPR3": "constraint_reference_ITPR3_8TKG.pdb",
}
MODEL_FILE = {"ITPR1": "constraint_model_Q14643.pdb",
              "ITPR3": "constraint_model_Q14573.pdb"}
SELECTION_FILE = {"ITPR3": "selection_reference_ITPR3_8TKG.pdb"}

UNSCORED = "#d6d5cf"


def numbering_check(path: Path, paralog: str) -> dict:
    """Does this structure carry the human canonical numbering?

    Positive test: every residue the file and the human per-residue table
    share must carry the same amino acid. Reports the fraction rather than a
    bare pass/fail, because the two that fail fail for different reasons and
    both are worth printing.
    """
    atoms = L.read_ca_pdb(path)
    index = L.residue_index(paralog)
    agree = sum(1 for a in atoms
                if index.get(a["resi"], {}).get("aa") == THREE_TO_ONE.get(a["aa"]))
    return {"file": path.name, "paralog": paralog, "residues": len(atoms),
            "agree": agree, "frac_agree": round(agree / len(atoms), 4),
            "carries_human_numbering": agree == len(atoms)}


def _trace(ax, atoms, title, letter, cmap, vmin, vmax, note="", flat=None):
    u, v = L.principal_axes(atoms)
    if flat is not None:
        # Supplementary Fig. 5 already paints the constraint map; here the
        # trace is context for the markers and a second heat map would both
        # repeat it and bury them.
        ax.plot(v, u, color=F.GRID, lw=0.35, zorder=1)
        ax.scatter(v, u, s=3.2, color=flat, lw=0, zorder=2)
        ax.set_aspect("equal")
        ax.axis("off")
        F.panel(ax, letter, title, pad=2.0)
        if note:
            ax.annotate(note, xy=(0.5, -0.01), xycoords="axes fraction",
                        ha="center", va="top", fontsize=F.FS_NOTE,
                        color=F.MUTED)
        return None
    scored = [(a, x, y) for a, x, y in zip(atoms, v, u) if a["b"] >= 0]
    blank = [(a, x, y) for a, x, y in zip(atoms, v, u) if a["b"] < 0]
    ax.plot(v, u, color=F.GRID, lw=0.35, zorder=1)
    if blank:
        ax.scatter([b[1] for b in blank], [b[2] for b in blank], s=3.2,
                   color=UNSCORED, lw=0, zorder=2)
    sc = ax.scatter([s[1] for s in scored], [s[2] for s in scored], s=3.2,
                    c=[s[0]["b"] for s in scored], cmap=cmap, vmin=vmin,
                    vmax=vmax, lw=0, zorder=3)
    ax.set_aspect("equal")
    ax.axis("off")
    F.panel(ax, letter, title, pad=2.0)
    if note:
        ax.annotate(note, xy=(0.5, -0.02), xycoords="axes fraction",
                    ha="center", va="top", fontsize=F.FS_NOTE, color=F.MUTED)
    return sc


# ------------------------------------------------------------------ fig 5

def fig_constraint_on_channel(out: Path, stats: dict) -> None:
    fig = plt.figure(figsize=(F.W_FULL, 7.6))
    gs = fig.add_gridspec(2, 2, hspace=0.20, wspace=0.02,
                          left=0.055, right=0.985, top=0.955, bottom=0.135)
    axes = gs.subplots()
    rec = {}
    letters = iter("abcd")
    sc = None
    for paralog, ax in zip(L.PARALOGS, axes.ravel()):
        atoms = L.read_ca_pdb(PAINTED / REFERENCE_FILE[paralog])
        scored = [a["b"] for a in atoms if a["b"] >= 0]
        pdb = REFERENCE_FILE[paralog].split("_")[-1].replace(".pdb", "")
        sc = _trace(ax, atoms,
                    f"{paralog}  ({pdb}, {len(atoms):,} resolved residues)",
                    next(letters), "viridis", 40, 95,
                    note=f"median deep-layer constraint "
                         f"{st.median(scored):.0f}")
        rec[paralog] = {"file": REFERENCE_FILE[paralog],
                        "residues": len(atoms), "scored": len(scored),
                        "median_b": round(st.median(scored), 2)}

    ax = axes.ravel()[3]
    atoms = L.read_ca_pdb(PAINTED / SELECTION_FILE["ITPR3"])
    scored = [a["b"] for a in atoms if a["b"] >= 0]
    ceiling = sum(1 for b in scored if b >= 99.999) / len(scored)
    sel = _trace(ax, atoms,
                 "ITPR3  (8TKG) painted with the selection layer instead",
                 next(letters), "cividis", 60, 100,
                 note=f"{100 * ceiling:.0f} % of scored residues sit at the "
                      f"ceiling: no non-synonymous\nsubstitution anywhere in "
                      f"the tree. {sum(1 for a in atoms if a['b'] < 0)} "
                      f"residues unscored")
    rec["ITPR3_selection"] = {"file": SELECTION_FILE["ITPR3"],
                              "residues": len(atoms), "scored": len(scored),
                              "median_b": round(st.median(scored), 2)}

    cax = fig.add_axes([0.10, 0.082, 0.30, 0.011])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("panels a–c: deep-layer constraint\n"
                 "(Jensen–Shannon divergence × 100)", fontsize=F.FS_NOTE)
    cb.ax.tick_params(labelsize=F.FS_TICK)
    cax2 = fig.add_axes([0.58, 0.082, 0.30, 0.011])
    cb2 = fig.colorbar(sel, cax=cax2, orientation="horizontal")
    cb2.set_label("panel d: selection layer\n(1 − dN rate) × 100",
                  fontsize=F.FS_NOTE)
    cb2.ax.tick_params(labelsize=F.FS_TICK)
    fig.legend(handles=[Line2D([], [], marker="o", ms=3, lw=0,
                               color=UNSCORED,
                               label="residue with no score (written −1 by "
                                     "S17, never 0)")],
               loc="lower center", bbox_to_anchor=(0.5, -0.004),
               fontsize=F.FS_TICK, frameon=False)
    F.save(fig, out / "SuppFig5_constraint_on_channel")
    plt.close(fig)
    stats["supp_fig_5"] = rec


# ------------------------------------------------------------------ fig 6

def _bh(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = min(prev, pvals[i] * n / (n - rank + 1))
        q[i] = val
        prev = val
    return q


def fig_variants_on_structure(out: Path, stats: dict) -> None:
    variants = L.read_tsv(L.CONSTRAINT_DIR / "variants.tsv")
    by_element = L.read_tsv(L.CONSTRAINT_DIR / "variant_by_element.tsv")

    checks = [numbering_check(PAINTED / p, g)
              for g, p in REFERENCE_FILE.items()]
    checks += [numbering_check(PAINTED / p, g) for g, p in MODEL_FILE.items()]
    # Where a paralogue has more than one structure that passes, the
    # experimental one is drawn: the predicted model of the same protein is
    # the same coordinates argued for rather than measured.
    eligible: dict[str, dict] = {}
    for c in checks:
        if not c["carries_human_numbering"]:
            continue
        cur = eligible.get(c["paralog"])
        if cur is None or (cur["file"].startswith("constraint_model")
                           and c["file"].startswith("constraint_reference")):
            eligible[c["paralog"]] = c

    fig = plt.figure(figsize=(F.W_FULL - 0.55, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[3.1, 1.3, 2.4], hspace=0.62,
                          wspace=0.06, left=0.13, right=0.98, top=0.95,
                          bottom=0.20)

    drawn = {}
    for k, paralog in enumerate(sorted(eligible)):
        ax = fig.add_subplot(gs[0, k])
        chk = eligible[paralog]
        pdb = chk["file"].replace(".pdb", "").split("_")[-1]
        atoms = L.read_ca_pdb(PAINTED / chk["file"])
        # A light ramp: the trace is context here, the markers are the
        # subject, and a full-contrast heat map would bury them.
        _trace(ax, atoms, f"{paralog}  ({pdb})", "ab"[k], None, 0, 1,
               flat=F.GRID)
        u, v = L.principal_axes(atoms)
        at = {a["resi"]: (x, y) for a, x, y in zip(atoms, v, u)}
        marks = {"P/LP": 0, "B/LB": 0}
        seen: set[tuple[str, int]] = set()
        for var in variants:
            if var["gene"] != paralog or var["class_bucket"] not in marks:
                continue
            key = (var["class_bucket"], int(var["resi"]))
            if key in seen or int(var["resi"]) not in at:
                continue
            seen.add(key)
            x, y = at[int(var["resi"])]
            ax.plot(x, y, marker="o", ms=4.2, mew=0.5,
                    markeredgecolor=F.SURFACE, zorder=6,
                    color=F.CLINICAL["pathogenic"]
                    if var["class_bucket"] == "P/LP"
                    else F.CLINICAL["benign"])
            marks[var["class_bucket"]] += 1
        ax.annotate(f"{marks['P/LP']} pathogenic, {marks['B/LB']} benign",
                    xy=(0.5, -0.01), xycoords="axes fraction", ha="center",
                    va="top", fontsize=F.FS_NOTE, color=F.MUTED)
        drawn[paralog] = marks

    ax = fig.add_subplot(gs[1, :])
    checks.sort(key=lambda c: -c["frac_agree"])
    ys = list(range(len(checks)))[::-1]
    ax.barh(ys, [c["frac_agree"] for c in checks], height=0.6,
            color=[F.QUALITY["complete"] if c["carries_human_numbering"]
                   else F.QUALITY["unannotated"] for c in checks])
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{c['paralog']}  "
                        f"{c['file'].replace('constraint_', '')[:-4]}"
                        for c in checks], fontsize=F.FS_NOTE)
    for y, c in zip(ys, checks):
        ax.annotate(f"{c['agree']:,}/{c['residues']:,}",
                    (c["frac_agree"], y), xytext=(3, 0),
                    textcoords="offset points", va="center",
                    fontsize=F.FS_NOTE, color=F.MUTED)
    ax.set_xlim(0, 1.2)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("residues whose amino acid matches the human per-residue "
                  "table")
    F.despine(ax)
    F.hgrid(ax, "x")
    F.panel(ax, "c", "which structures may carry a human variant position")

    ax = fig.add_subplot(gs[2, :])
    rows = [r for r in by_element if r["class_bucket"] == "P/LP"
            and L.f(r["odds_ratio"]) is not None]
    rows.sort(key=lambda r: (L.PARALOGS.index(r["gene"]),
                             L.ELEMENT_ORDER.index(r["element"])
                             if r["element"] in L.ELEMENT_ORDER else 99))
    qs = _bh([L.f(r["p_fisher"], 1.0) for r in rows])
    xs = list(range(len(rows)))
    # An element whose 2x2 table has a zero cell has an infinite odds ratio.
    # Plotted as it comes it lands off the axes and vanishes, which is how a
    # significant result disappears from its own figure: it is drawn at the
    # ceiling with a triangle instead, and the key says why.
    finite = [L.f(r["odds_ratio"], 0.0) for r in rows
              if math.isfinite(L.f(r["odds_ratio"], 0.0))]
    ceiling = max(finite) * 2.2
    n_inf = 0
    for x, r, q in zip(xs, rows, qs):
        raw = L.f(r["odds_ratio"], 0.0)
        capped = not math.isfinite(raw)
        odds = ceiling if capped else max(raw, 0.02)
        n_inf += int(capped)
        ax.plot([x, x], [1.0, odds], color=F.GRID, lw=0.8, zorder=1)
        ax.plot(x, odds, marker="^" if capped else "o",
                ms=5.0 if capped else (4.5 if q < 0.05 else 3.0),
                color=F.PARALOG[r["gene"]], mew=0.8, zorder=3,
                markeredgecolor=F.INK if q < 0.05 else "none")
        ax.annotate(r["n_positions"], (x, odds), xytext=(0, 6),
                    textcoords="offset points", ha="center",
                    fontsize=F.FS_NOTE, color=F.MUTED)
    ax.set_ylim(min(finite) / 2.2, ceiling * 2.4)
    # A rule between genes, because the same element name appears once per
    # gene and without it the three runs read as one series.
    for x in range(1, len(rows)):
        if rows[x]["gene"] != rows[x - 1]["gene"]:
            ax.axvline(x - 0.5, color=F.GRID, lw=0.8)
    ax.axhline(1.0, color=F.INK, lw=0.7, ls="--")
    ax.set_yscale("log")
    ax.set_xlim(-0.7, len(rows) - 0.3)
    ax.set_xticks(xs)
    ax.set_xticklabels([L.ELEMENT_LABEL.get(r["element"], r["element"])
                        for r in rows], rotation=40, ha="right",
                       fontsize=F.FS_NOTE)
    for x, r in zip(xs, rows):
        ax.get_xticklabels()[x].set_color(F.PARALOG[r["gene"]])
    ax.set_ylabel("odds ratio against\nthe rest of the protein")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "d", "where the pathogenic variants are, per element")
    ax.annotate(f"ringed = q < 0.05 after Benjamini–Hochberg over all "
                f"{len(rows)} element tests; the number above each point\n"
                f"is its count of pathogenic positions; a triangle is an "
                f"infinite odds ratio (a zero cell), drawn at the ceiling.\n"
                f"Element labels are coloured by gene.",
                xy=(0.0, -0.60), xycoords="axes fraction", va="top",
                fontsize=F.FS_NOTE, color=F.MUTED)

    fig.legend(handles=[
        Line2D([], [], marker="o", ms=3.4, lw=0,
               color=F.CLINICAL["pathogenic"], label="pathogenic / likely"),
        Line2D([], [], marker="o", ms=3.4, lw=0,
               color=F.CLINICAL["benign"], label="benign / likely")],
        loc="upper right", bbox_to_anchor=(0.985, 0.99),
        fontsize=F.FS_TICK, frameon=False)
    F.save(fig, out / "SuppFig6_variants_on_structure")
    plt.close(fig)

    stats["supp_fig_6"] = {
        "numbering_checks": checks,
        "structures_drawn": {p: eligible[p]["file"] for p in eligible},
        "positions_drawn": drawn,
        "element_tests": len(rows),
        "significant_after_bh": sum(1 for q in qs if q < 0.05),
        "infinite_odds_ratios": n_inf,
    }


FIGURES = {
    "supp5": fig_constraint_on_channel,
    "supp6": fig_variants_on_structure,
}
