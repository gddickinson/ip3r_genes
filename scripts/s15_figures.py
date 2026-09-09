"""s15_figures.py — the four S15a figures, from committed S15 tables only
(D13, D19).

Fig 1  **The character matrix.**  Every genome x paralog cell, ordered by
       assembly contiguity, with D4's bar drawn.  The state vocabulary is
       the figure's whole point: the panel exists so a reader can see that
       the red the S5 ledger showed is gone, and see *where* it went.

Fig 2  **The reconstruction, and what it is measured against.**  The three
       populations `s15_calibrate_recon.py` reads the bar off, with the
       operating point drawn and the gap shaded.  A calibration figure
       that asked to be believed would not be one, so both edges of the
       gap are marks and not a caption.

Fig 3  **Why synteny could not answer.**  Not a bar of call rates — the
       claim is about *availability*, so the panel is the joint
       distribution of what a trace region has to work with: genes on its
       own contig against informative flank keys, with the caller's
       four-key floor drawn, beside S8's own accuracy on the same axis.
       Putting accuracy and reach on one axis is what makes "accurate and
       unavailable" a readable sentence.

Fig 4  **The ORF screen and its confounders.**  Lesion density against
       each of the three things that could produce it without a gene
       being dead, and the paired within-genome test that removes two of
       them.  The identity panel is the one that matters and it is drawn
       first, because contiguity is the confounder everyone expects and
       identity is the one that turned out to be real.

Two rules the panels follow.  Density axes are **logarithmic with a zero
band**, because 72 % of intact loci carry no lesion at all and a linear
axis puts the entire calibration population on one pixel.  And the state
colours come from `figstyle.STATUS`, so a cell that is red here is red
for the same reason it was red in S5's ledger figures.
"""

from __future__ import annotations

import collections

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import figstyle as fs
import s15_lib as lib
import s15_states as st
import s5_calibration as s5cal

#: state -> colour, mapped onto the ledger's own evidence scale
STATE_COLOUR = {
    "present_single_locus": fs.STATUS["found_annotated"],
    "present_truncated": fs.STATUS["found_no_annotation"],
    "present_partial": fs.STATUS["fragment"],
    "present_fragmented": fs.STATUS["assembly_gap"],
    "paralog_unassignable": fs.STATUS["tblastn_trace_ambiguous"],
    "undecidable_contiguity": fs.STATUS["tblastn_trace"],
    "no_control": "#8a897f",
    "absent": fs.STATUS["absent"],
}
STATE_LABEL = {
    "present_single_locus": "one locus",
    "present_truncated": "truncated by the assembly",
    "present_partial": "partial locus",
    "present_fragmented": "reassembled across contigs",
    "paralog_unassignable": "paralog unassignable",
    "undecidable_contiguity": "undecidable (contiguity)",
    "no_control": "no positive control",
    "absent": "absent",
}
STATE_ORDER = ["present_single_locus", "present_truncated", "present_partial",
               "present_fragmented", "paralog_unassignable",
               "undecidable_contiguity", "no_control", "absent"]


def _t(name):
    return lib.read_tsv(lib.OUT / name)


def _f(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def _p(v: float) -> str:
    """A p-value that underflowed to 0.0 is not 'p = 0'."""
    if v != v:
        return "p n/a"
    return "p < 1e-300" if v <= 0.0 else f"p = {v:.0e}"


# ------------------------------------------------------------------ fig 1

def fig_matrix() -> None:
    matrix = _t("character_matrix.tsv")
    copies = _t("implied_copies.tsv")
    order = sorted(copies, key=lambda r: _f(r["contig_n50"]))
    pos = {r["accession"]: i for i, r in enumerate(order)}
    by_cell: dict[str, dict[str, str]] = collections.defaultdict(dict)
    for r in matrix:
        by_cell[r["cell"]][r["accession"]] = r["state"]

    fig = plt.figure(figsize=(fs.W_FULL, 4.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.75, 0.52, 1.35], hspace=0.30)
    ax = fig.add_subplot(gs[0])
    # ITPR1 on top: the rows read in the family's own order (PARALOG_ORDER),
    # which is the order every other figure in the project uses
    for row, cell in enumerate(reversed(lib.ITPR_CELLS)):
        for acc, i in pos.items():
            state = by_cell[cell].get(acc)
            if not state:
                continue
            ax.bar(i, 0.82, bottom=row + 0.09, width=1.0,
                   color=STATE_COLOUR.get(state, "#cccccc"), linewidth=0)
    n_below = sum(1 for r in order
                  if not s5cal.spans_a_gene(int(_f(r["contig_n50"]))))
    ax.axvline(n_below - 0.5, color=fs.INK, linewidth=0.9, zorder=5)
    # the note sits inside the panel: above the bars it collides with the
    # panel title, and a title a reader has to disentangle is not a title
    ax.annotate(f"D4 contiguity bar\n{n_below} genomes to the left cannot\n"
                f"hold this gene on one contig",
                xy=(n_below - 0.5, 0.02), xycoords=("data", "axes fraction"),
                xytext=(-5, 0), textcoords="offset points",
                ha="right", va="bottom", fontsize=fs.FS_NOTE, color=fs.INK,
                bbox=dict(facecolor=fs.SURFACE, edgecolor="none", pad=1.6,
                          alpha=0.93))
    ax.set_yticks([i + 0.5 for i in range(3)])
    ax.set_yticklabels(list(reversed(lib.ITPR_CELLS)))
    ax.set_ylim(0, 3.0)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_xlabel("309 vertebrate genomes, ordered by contig N50")
    fs.despine(ax, keep=("left",))
    ax.set_xticks([])
    fs.panel(ax, "a", "every genome × paralog cell, and the state its "
                      "evidence supports")

    # the legend gets its own strip: at 6.4 pt eight entries cannot share a
    # row with a panel title without one of them being unreadable
    lax = fig.add_subplot(gs[1])
    lax.axis("off")
    present = [s_ for s_ in STATE_ORDER
               if any(s_ == v for d in by_cell.values() for v in d.values())]
    absent_states = [s_ for s_ in STATE_ORDER if s_ not in present]
    lax.legend(handles=[Patch(facecolor=STATE_COLOUR[s_],
                              label=STATE_LABEL[s_]) for s_ in present],
               loc="upper center", ncol=3, fontsize=fs.FS_NOTE,
               title=("states reached; not reached anywhere: "
                      + ", ".join(STATE_LABEL[s_] for s_ in absent_states)),
               title_fontproperties={"size": fs.FS_NOTE, "weight": "bold"})

    ax = fig.add_subplot(gs[2])
    above = [r for r in copies if r["contig_spans_gene"] == "1"]
    below = [r for r in copies if r["contig_spans_gene"] != "1"]
    # bins centred on multiples of 0.25 so the three-paralog line falls in
    # the middle of its bar and not on the edge between two
    bins = [i * 0.25 - 0.125 for i in range(0, 22)]
    ax.hist([[min(_f(r["implied_copies"]), 5.0) for r in below],
             [min(_f(r["implied_copies"]), 5.0) for r in above]],
            bins=bins, stacked=True, linewidth=0,
            color=[fs.STATUS["assembly_gap"], fs.STATUS["found_annotated"]],
            label=[f"below D4's bar (n = {len(below)})",
                   f"above D4's bar (n = {len(above)})"])
    ax.axvline(3.0, color=fs.INK, linewidth=0.9, linestyle=(0, (3, 2)))
    ax.annotate("three paralogs", xy=(3.0, 0.97), xycoords=("data",
                                                            "axes fraction"),
                xytext=(4, 0), textcoords="offset points",
                fontsize=fs.FS_NOTE, color=fs.INK, va="top")
    ax.set_xlabel("ITPR gene-equivalents the assembly holds "
                  "(placed loci + reassembly; ≥ 5 pooled)")
    ax.set_ylabel("genomes")
    ax.set_xlim(-0.2, 5.2)
    fs.hgrid(ax)
    fs.despine(ax)
    ax.legend(loc="upper left", fontsize=fs.FS_NOTE)
    n_short = sum(1 for r in above if _f(r["implied_copies"]) < 2.5)
    fs.panel(ax, "b", f"{n_short} of {len(above)} genomes above the bar fall "
                      f"short of three")
    fs.save(fig, lib.FIGS / "s15_character_matrix")
    plt.close(fig)


# ------------------------------------------------------------------ fig 2

def fig_reconstruction() -> None:
    rows = [r for r in _t("reconstruction.tsv") if r["bait"]]
    cal = _t("recon_calibration.tsv")
    cov = next((r for r in cal if r["metric"] == "coverage"), {})
    scopes = ["own_clade", "co_trace", "decoy_accounted"]
    label = {"own_clade": "the cell's own clade\n(the candidate gene)",
             "co_trace": "a paralog that is itself\nshattered here",
             "decoy_accounted": "a paralog already placed\nat a locus "
                                "(the decoy)"}
    colour = {"own_clade": fs.STATUS["found_annotated"],
              "co_trace": fs.STATUS["assembly_gap"],
              "decoy_accounted": fs.STATUS["absent"]}

    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.5),
                             gridspec_kw={"width_ratios": [1.15, 1.0],
                                          "wspace": 0.28})
    ax = axes[0]
    lo, hi = _f(cov.get("gap_lo")), _f(cov.get("gap_hi"))
    op = _f(cov.get("operating_point"))
    ax.axhspan(lo, hi, color=fs.HILITE, zorder=0)
    ax.axhline(op, color=fs.INK, linewidth=0.9, zorder=1)
    rng = __import__("random").Random(15)
    for i, sc in enumerate(scopes):
        vals = [_f(r["coverage"]) for r in rows if r["scope"] == sc]
        ax.scatter([i + rng.uniform(-0.17, 0.17) for _ in vals], vals,
                   s=11, color=colour[sc], alpha=0.85, linewidth=0, zorder=3)
        ax.annotate(f"n = {len(vals)}", xy=(i, 1.03), ha="center",
                    fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.annotate(f"bar {op:.3f}\n(gap {lo:.3f}–{hi:.3f})",
                xy=(0.995, op), xycoords=("axes fraction", "data"),
                xytext=(0, 4), textcoords="offset points", ha="right",
                va="bottom", fontsize=fs.FS_NOTE, color=fs.INK)
    ax.set_xticks(range(3))
    ax.set_xticklabels([label[s] for s in scopes], fontsize=fs.FS_NOTE)
    ax.set_ylim(-0.03, 1.10)
    ax.set_ylabel("fraction of the reference reassembled\noutside every "
                  "placed locus")
    fs.hgrid(ax)
    fs.despine(ax)
    fs.panel(ax, "a", "the bar is read off an accounted-for gene")

    ax = axes[1]
    own = [r for r in rows if r["scope"] == "own_clade"]
    ax.scatter([int(r["n_contigs"]) for r in own],
               [_f(r["coverage"]) for r in own], s=13,
               color=fs.STATUS["found_annotated"], linewidth=0, alpha=0.85)
    ax.axhline(op, color=fs.INK, linewidth=0.9)
    ax.set_xlabel("contigs the gene is spread over")
    ax.set_ylabel("fraction of the reference recovered")
    ax.set_ylim(-0.03, 1.05)
    fs.hgrid(ax)
    fs.despine(ax)
    fs.panel(ax, "b", "every undecided cell is one gene in pieces")
    fs.save(fig, lib.FIGS / "s15_reconstruction")
    plt.close(fig)


# ------------------------------------------------------------------ fig 3

def fig_synteny() -> None:
    reach = [r for r in _t("synteny_reach.tsv")
             if r["window"] == "informative10"]
    acc = _t("caller_accuracy_by_keys.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.4),
                             gridspec_kw={"wspace": 0.3})

    ax = axes[0]
    counts = collections.Counter(
        (min(int(r["n_genes_on_contig"]), 6), min(int(r["n_keys"]), 6))
        for r in reach)
    xs = [k[0] for k in counts]
    ys = [k[1] for k in counts]
    ns = [counts[k] for k in counts]
    ax.scatter(xs, ys, s=[8 + 2.4 * n for n in ns],
               color=fs.STATUS["tblastn_trace"], linewidth=0, alpha=0.9)
    for (x, y), n in counts.items():
        if n >= 20:
            ax.annotate(str(n), xy=(x, y), ha="center", va="center",
                        fontsize=fs.FS_NOTE, color=fs.SURFACE)
    ax.axhline(3.5, color=fs.INK, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.annotate("the caller needs 4 keys", xy=(6, 3.5), xytext=(0, 3),
                textcoords="offset points", ha="right",
                fontsize=fs.FS_NOTE, color=fs.INK)
    ax.set_xlabel("coding genes on the region's own contig (≥ 6 pooled)")
    ax.set_ylabel("informative flank keys (≥ 6 pooled)")
    fs.hgrid(ax)
    fs.despine(ax)
    n0 = sum(1 for r in reach if int(r["n_genes_on_contig"]) == 0)
    fs.panel(ax, "a", f"{n0} of {len(reach)} regions: no gene on the contig")

    ax = axes[1]
    labels, rates, accs = [], [], []
    for r in acc:
        lo, hi = int(r["keys_lo"]), int(r["keys_hi"])
        labels.append(f"{lo}" if lo == hi else
                      (f"{lo}+" if hi >= 90 else f"{lo}–{hi}"))
        rates.append(_f(r["call_rate"]))
        accs.append(_f(r["accuracy"]) if r["accuracy"] not in ("", "nan")
                    else 0.0)
    x = range(len(labels))
    ax.bar([i - 0.19 for i in x], rates, width=0.36, linewidth=0,
           color=fs.STATUS["found_no_annotation"], label="call rate")
    ax.bar([i + 0.19 for i in x], accs, width=0.36, linewidth=0,
           color=fs.STATUS["found_annotated"], label="accuracy when called")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_xlabel("informative flank keys available")
    ax.set_ylim(0, 1.08)
    fs.hgrid(ax)
    fs.despine(ax)
    # One label per contiguous run of empty bins, centred on the run. Two
    # adjacent zero bins each carrying the same two-line note overlap into
    # unreadable text, which is what the S24 figure audit found here.
    runs, start = [], None
    for i, r in enumerate(rates + [1.0]):
        if r == 0.0 and start is None:
            start = i
        elif r != 0.0 and start is not None:
            runs.append((start, i - 1))
            start = None
    for lo, hi in runs:
        ax.annotate("below the\ncaller's floor", xy=((lo + hi) / 2, 0.03),
                    ha="center", va="bottom", fontsize=fs.FS_NOTE,
                    color=fs.MUTED)
        if hi > lo:
            ax.plot([lo - 0.34, hi + 0.34], [0.015, 0.015], color=fs.FAINT,
                    lw=0.7)
    ax.legend(loc="upper left", fontsize=fs.FS_NOTE)
    fs.panel(ax, "b", "accurate wherever called: reach failed, not accuracy")
    fs.save(fig, lib.FIGS / "s15_synteny_reach")
    plt.close(fig)


# ------------------------------------------------------------------ fig 4

def fig_integrity() -> None:
    loci = [r for r in _t("integrity_loci.tsv") if r["scored"] == "1"]
    tests = _t("integrity_tests.tsv")
    cov = _t("integrity_covariates.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(fs.W_FULL, 2.35),
                             gridspec_kw={"wspace": 0.34,
                                          "width_ratios": [1, 1, 0.95]})

    def dens(r):
        return max(_f(r["lesion_density"]), 0.04)      # the zero band

    # the RyR control is 886 of 1,760 loci, so it is drawn first: plotted
    # last it hides the three paralogs the comparison is about
    order = ([r for r in loci if r["cell"] == "RYR"]
             + [r for r in loci if r["cell"] != "RYR"])
    for ax, (xkey, xlab, logx, title) in zip(
            axes[:2],
            [("identity", "identity of the locus to its bait", False,
              "the alignment"),
             ("contig_n50", "contig N50 (bp)", True, "the assembly")]):
        colours = [fs.GROUP.get(r["cell"], fs.FAINT) for r in order]
        ax.scatter([_f(r[xkey]) for r in order], [dens(r) for r in order],
                   s=5, c=colours, linewidth=0, alpha=0.6)
        ax.set_yscale("log")
        if logx:
            ax.set_xscale("log")
        row = next((c for c in cov if c["covariate"] == xkey
                    and c["subset"] == "all"), {})
        ax.set_xlabel(xlab)
        ax.set_ylabel("disabling lesions per kilo-residue\n(0 shown at 0.04)")
        fs.hgrid(ax)
        fs.despine(ax)
        fs.panel(ax, "a" if xkey == "identity" else "b", title)
        # the statistic goes inside: three panels across a 6.7 in text
        # block cannot carry it in the title without the titles colliding
        ax.annotate(f"ρ = {_f(row.get('rho')):+.2f}\n"
                    f"{_p(_f(row.get('p')))}",
                    xy=(0.97, 0.03), xycoords="axes fraction", ha="right",
                    va="bottom", fontsize=fs.FS_NOTE, color=fs.INK)
    axes[0].legend(handles=[Patch(facecolor=fs.GROUP[c], label=c)
                            for c in lib.ALL_CELLS],
                   loc="upper right", ncol=2, fontsize=fs.FS_NOTE)

    ax = axes[2]
    matched = [t for t in tests if t["matched"] == "1"]
    cells = [t["cell"] for t in matched]
    pos = [int(t["n_pos"]) for t in matched]
    neg = [int(t["n_neg"]) for t in matched]
    y = range(len(cells))
    ax.barh([i + 0.0 for i in y], [-n for n in neg], height=0.6, linewidth=0,
            color=fs.STATUS["found_annotated"], label="fewer than siblings")
    ax.barh([i + 0.0 for i in y], pos, height=0.6, linewidth=0,
            color=fs.STATUS["absent"], label="more than siblings")
    for i, t in enumerate(matched):
        q = _f(t["q_bh"])
        ax.annotate(f"q = {q:.3f}" if q >= 0.001 else "q < 0.001",
                    xy=(pos[i], i), xytext=(4, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=fs.FS_NOTE,
                    color=fs.INK if q < 0.05 else fs.MUTED)
    ax.set_xlim(-max(neg) * 1.12, max(pos) * 1.62)
    ax.axvline(0, color=fs.INK, linewidth=0.7)
    ax.set_yticks(list(y))
    ax.set_yticklabels(cells)
    ax.set_xlabel("genomes (identity-matched\nwithin-genome pairs)")
    fs.despine(ax, keep=("bottom",))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
              fontsize=fs.FS_NOTE)
    fs.panel(ax, "c", "within one genome")
    fs.save(fig, lib.FIGS / "s15_integrity")
    plt.close(fig)


FIGURES = {"matrix": fig_matrix, "reconstruction": fig_reconstruction,
           "synteny": fig_synteny, "integrity": fig_integrity}


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
