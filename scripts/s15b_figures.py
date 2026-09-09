"""s15b_figures.py — the four S15b figures, from committed S15b tables only
(D13, D19).

Fig 1  **The deliverable.**  The sensitivity matrix as it is: the evidence
       ladder down, the two protective rules across, the Dollo count in
       each cell, once for the primary family coding and once for the
       paralog-resolved one.  S15a's own operating point is ringed.  The
       figure exists so a reader can see *where* in the grid a loss first
       appears and how far that is from where the task stands.

Fig 2  **The reconstruction bar and the one cell it decides.**  Every
       reassembled cell's coverage against the calibrated bar and its
       measured gap, coloured by what would catch the cell if the bar rose
       — beside the cumulative count of losses the rising bar manufactures
       under each rule setting.  The bar is *drawn*, with both gap edges,
       not stated in a caption.

Fig 3  **Why no Mk rate is reported.**  The likelihood profile of the
       primary character under all three models and all three
       branch-length schemes: monotone to the boundary, which is the whole
       argument for refusing to fit it.  ARD is drawn on its **gain**
       axis, not its diagonal — the diagonal is ER by construction and
       would put the same curve on the figure twice — and it is the one
       curve that rises, which is what unidentifiability looks like.  Beside it, the rates that *are*
       fitted on the manufactured settings, plotted against the number of
       losses the setting manufactured — because that is what they are a
       property of.

Fig 4  **The lesion result, stratified.**  D47's ITPR3 indel excess split
       by vertebrate class, with its contiguity control beside it.  The
       within-genome paired design makes a cell and its siblings the same
       comparison read from two sides, so the sign counts are drawn as
       opposed bars rather than as independent series.

Two rules the panels follow.  A count of zero is drawn as an explicit
zero and never as an empty cell, because an empty cell reads as *not
measured* and the zero is the result.  And the manufactured counts are on
a **linear** axis with the base setting ringed, since the claim is about
distance from the operating point and a log axis would flatter it.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

import figstyle as fs
import s15b_coding as coding
import s15b_lib as lib

RULE_COLS = [(1, 1), (0, 1), (1, 0), (0, 0)]
RULE_LABEL = {(1, 1): "both on\n(S15a)", (0, 1): "D45 off",
              (1, 0): "D4 off", (0, 0): "both off"}
CATCH_COLOUR = {"R5": fs.STATUS["tblastn_trace_ambiguous"],
                "R6": fs.STATUS["tblastn_trace"],
                "R7": fs.STATUS["absent"]}
CATCH_LABEL = {"R5": "D45: the genome has spare family loci",
               "R6": "D4: the assembly cannot span the gene",
               "R7": "nothing — the cell reads absent"}
#: the same three, short enough to be axis ticks at 6.6 pt
CATCH_TICK = {"R5": "spare loci (D45)", "R6": "contiguity (D4)",
              "R7": "nothing: absent"}


def _t(name):
    return lib.read_tsv(lib.OUT / f"{name}.tsv")


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


# ------------------------------------------------------ fig 1: the matrix

def fig_sensitivity() -> None:
    rows = _t("dollo_counts")
    ladder = list(coding.EVIDENCE_NAMES)
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.5),
                             constrained_layout=True)
    # One colour scale across both panels. Scaled per panel, the single
    # manufactured loss in the family-level grid is drawn in the same
    # maroon as the 45 in the paralogue-resolved one, and the panel whose
    # whole point is that it manufactures almost nothing reads as the
    # panel that manufactures most (found by the S24 figure audit).
    vmax = max(
        (_i(r["dollo_losses_max"]) for r in rows
         if r["coding"] in ("family", "paralog")
         and r["evidence"] in ladder), default=1) or 1
    for ax, coding_name, title in zip(
            axes, ("family", "paralog"),
            ("family-level presence per genome (primary, D46)",
             "paralog-resolved, worst of ITPR1/2/3")):
        grid = []
        for ev in ladder:
            line = []
            for r5, r6 in RULE_COLS:
                cells = [r for r in rows if r["coding"] == coding_name
                         and r["evidence"] == ev
                         and _i(r["use_r5"]) == r5 and _i(r["use_r6"]) == r6]
                line.append(max((_i(c["dollo_losses_max"]) for c in cells),
                                default=0))
            grid.append(line)
        ax.imshow(grid, cmap="Reds", vmin=0, vmax=vmax, aspect="auto")
        for i, line in enumerate(grid):
            for j, v in enumerate(line):
                ax.text(j, i, str(v), ha="center", va="center",
                        fontsize=fs.FS_NOTE,
                        color=fs.SURFACE if v > 0.55 * vmax else fs.INK)
        base_i = ladder.index(coding.BASE_SETTING["evidence"])
        ax.add_patch(Rectangle((-0.5, base_i - 0.5), 1, 1, fill=False,
                               edgecolor=fs.ACCENT, linewidth=1.6, zorder=5))
        ax.set_xticks(range(len(RULE_COLS)))
        ax.set_xticklabels([RULE_LABEL[k] for k in RULE_COLS],
                           fontsize=fs.FS_TICK)
        ax.set_yticks(range(len(ladder)))
        ax.set_yticklabels(ladder, fontsize=fs.FS_TICK)
        if coding_name == "family":
            ax.set_ylabel("evidence accepted as presence\n"
                          "(loosest at the top)", fontsize=fs.FS_LABEL)
        fs.panel(ax, "a" if coding_name == "family" else "b", title)
        ax.tick_params(length=0)
    fig.suptitle("Dollo losses manufactured by each setting of the rule "
                 "chain (309 vertebrate genomes);\nthe violet ring is "
                 "S15a's own operating point", fontsize=fs.FS_SUPTITLE)
    fs.save(fig, lib.FIGS / "sensitivity_matrix")
    plt.close(fig)


# --------------------------------------------- fig 2: the bar and its cell

def fig_recon_bar() -> None:
    matrix = lib.read_tsv(lib.MATRIX)
    frag = [r for r in matrix if r["state"] == "present_fragmented"]
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.9),
                             constrained_layout=True,
                             gridspec_kw=dict(width_ratios=(1.45, 1.0)))

    ax = axes[0]
    # one row per catch class, and the within-row offset is the point's
    # rank in its own row: no hash, no RNG, so the figure is reproducible
    # (D24 applied to a jitter).
    rows_by_catch = {"R5": [], "R6": [], "R7": []}
    for r in sorted(frag, key=lambda x: _f(x["recon_coverage"])):
        catch = ("R5" if _i(r["n_spare_itpr_loci"]) > 0 else
                 "R6" if not _i(r["contig_spans_gene"]) else "R7")
        rows_by_catch[catch].append(r)
    for row, (catch, group) in enumerate(rows_by_catch.items()):
        for k, r in enumerate(group):
            ax.plot(_f(r["recon_coverage"]), row + 0.06 * (k % 5) - 0.12,
                    "o", ms=3.4, mfc=CATCH_COLOUR[catch], mec=fs.SURFACE,
                    mew=0.4, zorder=3)
        ax.text(1.0, row, f" {len(group)}", fontsize=fs.FS_NOTE,
                color=fs.MUTED, va="center")
    ax.axvspan(coding.BAR_GAP_LO, coding.BAR_GAP_HI, color=fs.HILITE,
               zorder=0)
    for x, lab in ((coding.BAR_GAP_LO, "decoy max"),
                   (coding.BAR_GAP_HI, "lowest candidate")):
        ax.axvline(x, color=fs.FAINT, lw=0.7, ls=":", zorder=1)
    ax.axvline(coding.BAR_CALIBRATED, color=fs.ACCENT, lw=1.1, zorder=2)
    ax.text(coding.BAR_CALIBRATED, 2.62, " calibrated bar, in its "
            "measured gap", fontsize=fs.FS_NOTE, color=fs.ACCENT,
            ha="left", va="top")
    ax.set_xlim(0, 1.06)
    ax.set_ylim(-0.55, 2.7)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels([CATCH_TICK[k] for k in ("R5", "R6", "R7")],
                       fontsize=fs.FS_TICK)
    ax.set_xlabel("reference reassembled across contigs (fraction)",
                  fontsize=fs.FS_LABEL)
    ax.tick_params(axis="y", length=0)
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "a", f"the {len(frag)} reassembled cells, by what catches "
                      f"each if the bar rises")

    ax = axes[1]
    bars = [b / 100 for b in range(0, 101, 2)]
    for (r5, r6), style in zip(RULE_COLS, ("-", "--", "-.", ":")):
        ys = []
        for b in bars:
            n = 0
            for r in matrix:
                if _f(r["recon_coverage"]) >= b:
                    continue
                if r["ledger_status"].startswith("found") or \
                        r["ledger_status"] in ("assembly_gap", "fragment"):
                    continue
                if r5 and _i(r["n_spare_itpr_loci"]) > 0:
                    continue
                if r6 and not _i(r["contig_spans_gene"]):
                    continue
                n += 1
            ys.append(n)
        ax.plot(bars, ys, style, color=fs.INK if (r5 and r6) else fs.MUTED,
                lw=1.4 if (r5 and r6) else 0.9, label=
                RULE_LABEL[(r5, r6)].replace("\n", " "))
    ax.axvline(coding.BAR_CALIBRATED, color=fs.ACCENT, lw=1.1)
    ax.set_xlabel("reconstruction bar", fontsize=fs.FS_LABEL)
    ax.set_ylabel("cells reading absent", fontsize=fs.FS_LABEL)
    ax.set_xlim(0, 1.0)
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="upper left",
              handlelength=1.6)
    fs.panel(ax, "b", "and the losses that makes")
    fs.save(fig, lib.FIGS / "reconstruction_bar")
    plt.close(fig)


# ------------------------------------------------------- fig 3: the Mk half

def fig_mk() -> None:
    prof = _t("mk_profile")
    fits = [r for r in _t("mk_fits") if _i(r["fitted"])]
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.8),
                             constrained_layout=True)

    ax = axes[0]
    styles = {"unit": "-", "ultrametric": "--", "calibrated": ":"}
    colours = {"ER (loss)": fs.PARALOG["ITPR1"],
               "ARD (gain)": fs.PARALOG["ITPR2"],
               "irreversible (loss)": fs.PARALOG["ITPR3"]}
    series = {"ER (loss)": ("ER", "loss"), "ARD (gain)": ("ARD", "gain"),
              "irreversible (loss)": ("irreversible", "loss")}
    for scheme in lib.BL_SCHEMES:
        for label, (model, axis) in series.items():
            g = [r for r in prof if r["branch_lengths"] == scheme
                 and r["model"] == model and r.get("axis", "loss") == axis]
            if not g:
                continue
            g.sort(key=lambda r: _f(r["rate"]))
            best = max(_f(r["loglik"]) for r in g)
            ax.plot([_f(r["rate"]) for r in g],
                    [_f(r["loglik"]) - best for r in g],
                    styles[scheme], color=colours[label], lw=1.0)
    ax.set_xscale("log")
    # the three models differ by four orders of magnitude in how fast the
    # likelihood falls; on a linear axis only the steepest is visible and
    # the other two read as flat, which is the opposite of the claim
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("transition rate", fontsize=fs.FS_LABEL)
    ax.set_ylabel("log likelihood, relative to its own maximum",
                  fontsize=fs.FS_LABEL)
    fs.despine(ax)
    fs.hgrid(ax)
    handles = [Line2D([], [], color=colours[m], lw=1.2, label=m)
               for m in series]
    handles += [Line2D([], [], color=fs.MUTED, ls=styles[s], lw=1.0, label=s)
                for s in lib.BL_SCHEMES]
    ax.legend(handles=handles, fontsize=fs.FS_NOTE, frameon=False, ncol=2,
              loc="lower left", handlelength=1.6)
    fs.panel(ax, "a", "the primary character: monotone to the boundary")

    ax = axes[1]
    if fits:
        for scheme, mark in zip(lib.BL_SCHEMES, ("o", "s", "^")):
            g = [r for r in fits if r["branch_lengths"] == scheme
                 and r["model"] == "irreversible"]
            ax.plot([_i(r["n_absent"]) for r in g],
                    [_f(r["rate"]) for r in g], mark, ms=3.0,
                    mfc="none", mew=0.8,
                    color=fs.INK if scheme == "unit" else fs.MUTED,
                    label=scheme, ls="none")
        ax.set_yscale("log")
    else:                                                # pragma: no cover
        ax.text(0.5, 0.5, "no setting produced a fittable character",
                ha="center", va="center", fontsize=fs.FS_NOTE,
                color=fs.MUTED, transform=ax.transAxes)
    ax.set_xlabel("cells the setting made absent", fontsize=fs.FS_LABEL)
    ax.set_ylabel("fitted loss rate (irreversible)", fontsize=fs.FS_LABEL)
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right",
              handlelength=1.0)
    fs.panel(ax, "b", "every rate that can be fitted is a rate of the filter")
    fs.save(fig, lib.FIGS / "mk_profile")
    plt.close(fig)


# ---------------------------------------------------- fig 4: D47's lead

def fig_lesions() -> None:
    tests = [t for t in _t("lesion_by_class") if not _i(t["underpowered"])]
    ctrl = _t("lesion_class_controls")
    fossil = _t("fossil_by_cell")
    fig, axes = plt.subplots(1, 3, figsize=(fs.W_FULL, 2.8),
                             constrained_layout=True,
                             gridspec_kw=dict(width_ratios=(1.35, 1.0, 1.0)))

    ax = axes[0]
    tests.sort(key=lambda t: (t["vclass"], t["cell"]))
    labels, y = [], []
    for i, t in enumerate(tests):
        col = fs.PARALOG.get(t["cell"], fs.GROUP["RYR"])
        ax.barh(i, _i(t["n_pos"]), color=col, height=0.7)
        ax.barh(i, -_i(t["n_neg"]), color=col, height=0.7, alpha=0.42)
        q = _f(t["q_bh"], 1.0)
        if q <= 0.05:
            ax.text(_i(t["n_pos"]) + 1.2, i, f"q={lib.pfmt(q)}",
                    fontsize=fs.FS_NOTE, va="center", color=fs.INK)
        labels.append(f"{t['cell']} · {t['vclass']}")
        y.append(i)
    ax.axvline(0, color=fs.INK, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("genomes with a lesion deficit  ←→  excess",
                  fontsize=fs.FS_LABEL)
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "a", "within-genome, identity-matched, by class")

    ax = axes[1]
    strong = sorted(ctrl, key=lambda c: _f(c["q_bh"], 1.0))[:1]
    if strong:
        c = strong[0]
        groups = [("above D4's bar", _i(c["pos_above"]), _i(c["neg_above"])),
                  ("below D4's bar", _i(c["pos_below"]), _i(c["neg_below"]))]
        for i, (lab, pos, neg) in enumerate(groups):
            ax.bar(i - 0.17, pos, width=0.32,
                   color=fs.PARALOG.get(c["cell"], fs.INK), label=None)
            ax.bar(i + 0.17, neg, width=0.32,
                   color=fs.PARALOG.get(c["cell"], fs.INK), alpha=0.42)
            ax.text(i, max(pos, neg) + 0.6, f"n={pos + neg}",
                    ha="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([g[0] for g in groups], fontsize=fs.FS_TICK)
        ax.set_ylabel("genomes", fontsize=fs.FS_LABEL)
        ax.text(0.5, -0.30, "solid: genomes with an excess.  pale: with a "
                "deficit", transform=ax.transAxes, ha="center",
                fontsize=fs.FS_NOTE, color=fs.MUTED)
        fs.panel(ax, "b", f"{c['cell']}/{c['vclass']} by contiguity")
    fs.despine(ax)
    fs.hgrid(ax)

    ax = axes[2]
    cells = [f for f in fossil if _i(f["n_scored"])]
    xs = range(len(cells))
    ax.bar(xs, [_i(f["n_above_bar"]) for f in cells], width=0.62,
           color=[fs.PARALOG.get(f["cell"], fs.GROUP["RYR"]) for f in cells])
    ax.bar(xs, [_i(f["n_fossils"]) for f in cells], width=0.62,
           color=fs.STATUS["absent"])
    for i, f in enumerate(cells):
        ax.text(i, _i(f["n_above_bar"]) + 0.6, str(_i(f["n_scored"])),
                ha="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([f["cell"] for f in cells], fontsize=fs.FS_TICK)
    ax.set_ylabel("loci above the lesion bar\n(number above the bar: loci "
                  "scored)", fontsize=fs.FS_LABEL)
    ax.legend(handles=[Patch(facecolor=fs.STATUS["absent"],
                             label="fossil under any reading")],
              fontsize=fs.FS_NOTE, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.17), handlelength=1.0)
    ax.set_ylim(0, max(_i(f["n_above_bar"]) for f in cells) * 1.25)
    fs.despine(ax)
    fs.hgrid(ax)
    fs.panel(ax, "c", "no dead loci to read a lesion off")
    fs.save(fig, lib.FIGS / "lesion_strata")
    plt.close(fig)


FIGURES = {
    "sensitivity_matrix": fig_sensitivity,
    "reconstruction_bar": fig_recon_bar,
    "mk_profile": fig_mk,
    "lesion_strata": fig_lesions,
}


def main(only: list[str] | None = None) -> None:
    fs.use()
    lib.FIGS.mkdir(parents=True, exist_ok=True)
    for name, fn in FIGURES.items():
        if only and name not in only:
            continue
        fn()
        print(f"  [fig] {name}")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] or None)
