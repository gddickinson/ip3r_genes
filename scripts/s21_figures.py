"""The four S21 figures, from committed S21 tables only (D13, D19).

Four decisions worth naming.

* The architecture figure puts **exon count and genomic span on separate axes**
  and never on one pair of bars. The whole result is that one is conserved and
  the other is not, and a shared scale on which 58 exons and 244,000 bp both
  appear would hide it.
* Span is on a **logarithmic** axis and exon count is not: spans run over two
  orders of magnitude across the scope and counts over a factor of 1.5.
* The shared-intron figure plots **observed against expected per genome**, with
  the identity line drawn, rather than a bar of enrichments. The claim is that
  every genome sits far off the diagonal in the same direction, and a mean
  enrichment cannot show that a single genome does.
* The junction figure is drawn on the **complement** — the fraction of junctions
  the genome does *not* read as a splice pair — because at 99.9 % agreement a
  bar chart of the agreement is four full bars and says nothing. An exact zero is
  written as `0`, never left as an empty space (S15b's rule).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import figstyle as FS                                             # noqa: E402
import s21_lib as L                                               # noqa: E402

FS.use()
import matplotlib.pyplot as plt                                   # noqa: E402

CELLS = list(L.PARALOGS) + [L.CONTROL_CELL]


def _colour(cell: str) -> str:
    return FS.GROUP.get(cell, FS.GROUP["RYR"])


def _label(cell: str) -> str:
    return "RyR" if cell == "RYR" else cell


def _f(row: dict, key: str, default=0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------------
def fig_architecture() -> None:
    """Exon count conserved, genomic span not — and the sensitivity beside it."""
    per = {r["cell"]: r for r in L.read_tsv(L.OUT / "architecture_by_paralog.tsv")}
    loci = L.read_tsv(L.OUT / "architecture_loci.tsv")
    sens = L.read_tsv(L.OUT / "architecture_sensitivity.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FS.W_FULL, 2.5))

    ax = axes[0]
    xs = range(len(CELLS))
    for i, cell in enumerate(CELLS):
        r = per.get(cell)
        if not r:
            continue
        med = _f(r, "n_exons_median")
        lo, hi = _f(r, "n_exons_p10"), _f(r, "n_exons_p90")
        ax.bar(i, med, color=_colour(cell), width=0.62)
        ax.errorbar(i, med, yerr=[[med - lo], [hi - med]], fmt="none",
                    ecolor=FS.INK, elinewidth=0.8, capsize=2.0)
        ax.annotate(f"{med:.0f}", (i, hi), xytext=(0, 3),
                    textcoords="offset points", ha="center",
                    fontsize=FS.FS_TICK, color=FS.INK)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([_label(c) for c in CELLS])
    ax.set_ylabel("coding exons per gene")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "a", "exons: 57-58 against 104")

    ax = axes[1]
    for i, cell in enumerate(CELLS):
        vals = [_f(r, "span_bp") / 1000.0 for r in loci
                if r["cell"] == cell and r["in_scope"] == "1"]
        if not vals:
            continue
        parts = ax.violinplot([vals], positions=[i], widths=0.7,
                              showextrema=False, showmedians=True)
        for body in parts["bodies"]:
            body.set_facecolor(_colour(cell))
            body.set_alpha(0.75)
            body.set_edgecolor("none")
        parts["cmedians"].set_color(FS.INK)
        parts["cmedians"].set_linewidth(0.9)
    ax.set_yscale("log")
    ax.set_xticks(list(range(len(CELLS))))
    ax.set_xticklabels([_label(c) for c in CELLS])
    ax.set_ylabel("genomic span (kb)")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "b", "span: a 4.2-fold spread")

    ax = axes[2]
    # ITPR1 and ITPR3 sit on the same value at every bar, so the series are
    # given distinct markers: with one marker the upper line simply hides the
    # other and the figure reads as a missing paralogue.
    markers = {"ITPR1": "o", "ITPR2": "s", "ITPR3": "^", "RYR": "D"}
    for cell in CELLS:
        rows = sorted((r for r in sens if r["cell"] == cell),
                      key=lambda r: _f(r, "cov_bar"))
        if not rows:
            continue
        ax.plot([_f(r, "cov_bar") for r in rows],
                [_f(r, "n_exons_median") for r in rows],
                marker=markers.get(cell, "o"), markersize=3.0,
                markerfacecolor="none", markeredgewidth=0.9, linewidth=1.0,
                color=_colour(cell), label=_label(cell))
    ax.set_xlabel("bait-coverage bar")
    ax.set_ylabel("median exon count")
    ax.legend(frameon=False, fontsize=FS.FS_TICK, loc="center right")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "c", "flat across the bar")
    fig.tight_layout()
    FS.save(fig, L.FIGS / "architecture_by_paralog")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_introns() -> None:
    """Intron positions: how many are conserved, and how many are shared."""
    cons = L.read_tsv(L.OUT / "intron_conservation.tsv")
    summ = {r["cell"]: r for r in
            L.read_tsv(L.OUT / "intron_conservation_summary.tsv")}
    shared = [r for r in L.read_tsv(L.OUT / "shared_introns.tsv")
              if r["is_operating_point"] == "1"]
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))

    ax = axes[0]
    for cell in CELLS:
        vals = sorted((_f(r, "prevalence") for r in cons if r["cell"] == cell),
                      reverse=True)
        if not vals:
            continue
        n = int(_f(summ.get(cell, {}), "n_at_90pc"))
        ax.plot(range(1, len(vals) + 1), vals, linewidth=1.2,
                color=_colour(cell),
                label=f"{_label(cell)} ({n} at 90 %)")
    ax.axhline(0.90, color=FS.MUTED, linewidth=0.7, linestyle=(0, (3, 2)))
    ax.annotate("90 % of loci", xy=(0.98, 0.90), xycoords=("axes fraction",
                                                           "data"),
                xytext=(0, 3), textcoords="offset points", ha="right",
                fontsize=FS.FS_TICK, color=FS.MUTED)
    ax.set_xlabel("intron position, ranked by prevalence")
    ax.set_ylabel("fraction of loci carrying it")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=FS.FS_TICK, loc="lower left")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "a", "a universal core, per paralogue")

    ax = axes[1]
    pair_colour = {("ITPR1", "ITPR2"): FS.PARALOG["ITPR1"],
                   ("ITPR1", "ITPR3"): FS.PARALOG["ITPR3"],
                   ("ITPR2", "ITPR3"): FS.PARALOG["ITPR2"]}
    for (a, b), col in pair_colour.items():
        pts = [(max(_f(r, "expected_analytic"), 0.16), _f(r, "observed"))
               for r in shared if r["cell_a"] == a and r["cell_b"] == b]
        if pts:
            ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=5.0,
                       color=col, alpha=0.7, linewidths=0,
                       label=f"{a} vs {b}")
    ctrl = [(max(_f(r, "expected_analytic"), 0.16), _f(r, "observed"))
            for r in shared if r["cell_b"] == "RYR"]
    if ctrl:
        ax.scatter([p[0] for p in ctrl], [p[1] for p in ctrl], s=5.0,
                   color=FS.GROUP["RYR"], alpha=0.7, linewidths=0,
                   label="ITPR vs RyR (control)")
    lim = [0.15, 3.0]
    ax.plot(lim, lim, color=FS.MUTED, linewidth=0.7, linestyle=(0, (3, 2)))
    ax.annotate("chance", xy=(1.6, 1.6), xytext=(3, -9),
                textcoords="offset points", fontsize=FS.FS_TICK,
                color=FS.MUTED)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlim(*lim)
    ax.set_ylim(0, 90)
    ax.set_xlabel("intron positions expected by chance")
    ax.set_ylabel("positions shared")
    ax.legend(frameon=False, fontsize=FS.FS_TICK, loc="center left",
              bbox_to_anchor=(0.0, 0.45))
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "b", "one genome per point")
    fig.tight_layout()
    FS.save(fig, L.FIGS / "intron_positions")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_junctions() -> None:
    """What the genome says about the boundaries the aligner placed."""
    q = {r["cell"]: r for r in L.read_tsv(L.OUT / "junction_quality.tsv")}
    bins = L.read_tsv(L.OUT / "junction_gap_bins.tsv")
    conc = L.read_tsv(L.OUT / "boundary_concordance_summary.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FS.W_FULL, 2.4))

    ax = axes[0]
    for i, cell in enumerate(CELLS):
        r = q.get(cell, {})
        pct = 100.0 * (1.0 - _f(r, "frac_spliceable"))
        ax.bar(i, pct, color=_colour(cell), width=0.62)
        ax.annotate(f"{pct:.2f}" if pct else "0", (i, pct), xytext=(0, 3),
                    textcoords="offset points", ha="center",
                    fontsize=FS.FS_TICK, color=FS.INK)
    ax.set_xticks(list(range(len(CELLS))))
    ax.set_xticklabels([_label(c) for c in CELLS])
    ax.set_ylabel("junctions not spliceable (%)")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "a", "the instrument's error rate")

    ax = axes[1]
    xs, ys, ns = [], [], []
    for r in bins:
        lo = r["gap_lo"] or "0"
        hi = r["gap_hi"] or ""
        xs.append(f"{lo}-{hi}" if hi else f">{lo}")
        ys.append(100.0 * _f(r, "frac_spliceable"))
        ns.append(int(_f(r, "n")))
    ax.bar(range(len(xs)), ys, color=FS.BLUES[3], width=0.7)
    # The bins differ by four orders of magnitude in size — the 21-30 bp bin
    # holds seven gaps and the >10 kb bin holds 7,032 — so each bar carries the
    # count it rests on. Without it the one bin below 100 % reads as a real
    # population rather than as one non-canonical junction in seven.
    for i, (y, n) in enumerate(zip(ys, ns)):
        ax.annotate(f"{n}", (i, y), xytext=(0, 2), textcoords="offset points",
                    ha="center", fontsize=FS.FS_TICK - 2, color=FS.MUTED,
                    rotation=90)
    ax.set_xticks(list(range(len(xs))))
    ax.set_xticklabels(xs, rotation=90, fontsize=FS.FS_TICK - 1)
    ax.set_ylim(0, 128)
    ax.set_xlabel("gap between consecutive blocks (bp)")
    ax.set_ylabel("spliceable (%)")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "b", "no sub-intron population")

    ax = axes[2]
    rows = [r for r in conc if r["grouping"] == "cell"] + \
           [r for r in conc if r["grouping"] == "source"]
    labels, vals, cols = [], [], []
    for r in rows:
        labels.append(_label(r["group"]))
        vals.append(100.0 * _f(r, "frac_exact"))
        cols.append(_colour(r["group"]) if r["grouping"] == "cell"
                    else FS.MUTED)
    ax.barh(range(len(labels)), vals, color=cols, height=0.66)
    ax.set_yticks(list(range(len(labels))))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("edges on an exon boundary (%)")
    FS.hgrid(ax, axis="x")
    FS.despine(ax)
    FS.panel(ax, "c", "an independent pipeline agrees")
    fig.tight_layout()
    FS.save(fig, L.FIGS / "junction_quality")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_fragments() -> None:
    """Where a fragmentary annotation stops, and whether anything splices there."""
    summ = {r["state"]: r for r in L.read_tsv(L.OUT / "fragment_summary.tsv")}
    verd = L.read_tsv(L.OUT / "fragment_verdicts.tsv")
    ctrl = {r["cell"]: r for r in L.read_tsv(L.OUT / "tandem_control.tsv")}
    fig, axes = plt.subplots(1, 3, figsize=(FS.W_FULL, 2.4))

    order = ["annotation_failure", "structure_disagreement", "mixed",
             "broken_at_real_junctions", "no_internal_terminus"]
    colour = {"annotation_failure": FS.QUALITY["unannotated"],
              "structure_disagreement": FS.QUALITY["noncoding"],
              "mixed": FS.QUALITY["fragmentary"],
              "broken_at_real_junctions": FS.QUALITY["complete"],
              "no_internal_terminus": FS.QUALITY["split"]}
    ax = axes[0]
    states = [s for s in ("split", "fragmentary") if s in summ]
    drawn: list[str] = []
    for i, state in enumerate(states):
        left = 0.0
        total = int(_f(summ[state], "n_loci")) or 1
        for v in order:
            n = int(_f(summ[state], v))
            if not n:
                continue
            if v not in drawn:
                drawn.append(v)
            ax.barh(i, 100.0 * n / total, left=left, height=0.62,
                    color=colour[v])
            left += 100.0 * n / total
    ax.set_yticks(list(range(len(states))))
    ax.set_yticklabels([f"{s}\n(n={int(_f(summ[s], 'n_loci'))})"
                        for s in states])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("loci (%)")
    # The legend is built from the segments actually drawn, in both bars: with
    # per-bar labels the `mixed` segment appeared in the fragmentary bar and in
    # no legend, because the split bar has none of it.
    handles = [plt.Rectangle((0, 0), 1, 1, color=colour[v]) for v in drawn]
    ax.legend(handles, [v.replace("_", " ") for v in drawn], frameon=False,
              fontsize=FS.FS_TICK - 1, loc="upper center",
              bbox_to_anchor=(0.5, -0.42), ncol=2)
    FS.hgrid(ax, axis="x")
    FS.despine(ax)
    FS.panel(ax, "a", "verdict per locus")

    ax = axes[1]
    cats = [("termini_at_exon_boundary", "at a junction"),
            ("termini_mid_exon", "inside an exon"),
            ("termini_in_intron", "inside an intron")]
    row = summ.get("all", {})
    total = int(_f(row, "internal_termini")) or 1
    cols = [FS.QUALITY["complete"], FS.QUALITY["unannotated"],
            FS.QUALITY["noncoding"]]
    for i, ((key, lab), col) in enumerate(zip(cats, cols)):
        n = int(_f(row, key))
        ax.bar(i, 100.0 * n / total, color=col, width=0.62)
        ax.annotate(f"{n}", (i, 100.0 * n / total), xytext=(0, 3),
                    textcoords="offset points", ha="center",
                    fontsize=FS.FS_TICK, color=FS.INK)
    ax.set_xticks(list(range(len(cats))))
    ax.set_xticklabels([c[1] for c in cats], rotation=20, ha="right")
    ax.set_ylabel(f"of {total} internal termini (%)")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "b", "where the annotation stops")

    ax = axes[2]
    cells = [c for c in CELLS if c in ctrl]
    for i, cell in enumerate(cells):
        sens = _f(ctrl[cell], "sensitivity") * 100.0
        spec = _f(ctrl[cell], "specificity") * 100.0
        ax.bar(i - 0.17, sens, width=0.32, color=_colour(cell))
        ax.bar(i + 0.17, spec, width=0.32, color=_colour(cell), alpha=0.45)
    ax.set_xticks(list(range(len(cells))))
    ax.set_xticklabels([_label(c) for c in cells])
    ax.set_ylim(0, 105)
    ax.set_ylabel("%")
    ax.annotate("solid: sensitivity\npale: specificity", xy=(0.02, 0.06),
                xycoords="axes fraction", fontsize=FS.FS_TICK,
                color=FS.MUTED)
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "c", "duplication detector vs S16's copy call")
    fig.tight_layout()
    FS.save(fig, L.FIGS / "fragments_and_duplicates")
    plt.close(fig)


def main() -> int:
    L.FIGS.mkdir(parents=True, exist_ok=True)
    for fn in (fig_architecture, fig_introns, fig_junctions, fig_fragments):
        print(f"  [fig] {fn.__name__}")
        fn()
    return 0


if __name__ == "__main__":
    sys.exit(main())
