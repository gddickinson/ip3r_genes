"""S2 figures — rendered only from the committed census tables (D13, D19).

Four panels, each answering one of the task's questions:

  census_space      what the enumeration actually contains, by lineage
  census_lengths    the two families as length distributions, and the
                    fragment pile that neither of them explains
  census_margin     the sequence-level check on the architecture call
  census_growth     what a domain search adds over a search by name

No network, no re-derivation: every number is read back out of
`results/census_v2/*.tsv`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402

import figstyle                                     # noqa: E402
from scripts.s2_lib import OUT_DIR, read_tsv        # noqa: E402

FIG_DIR = OUT_DIR / "figures"

CALL_COLOUR = {
    "ITPR": figstyle.PARALOG["ITPR1"],
    "RYR": figstyle.GROUP["RYR"],
    "unassigned": "#a9a79e",
}
CALL_ORDER = ["ITPR", "RYR", "unassigned"]


def provenance(fig, note: str) -> None:
    """Corner tag: what kind of number this figure is made of (D19)."""
    from matplotlib.transforms import offset_copy
    tr = offset_copy(fig.transFigure, fig=fig, x=0, y=-7, units="points")
    fig.text(0.995, 0.0, note, ha="right", va="top", transform=tr,
             fontsize=5.6, color="#8a897f")


def census() -> list[dict]:
    return read_tsv(OUT_DIR / "census_v2.tsv")


# ------------------------------------------------------------------ fig 1
def fig_space() -> None:
    """Size and composition, kept on separate axes.

    Counts span four orders of magnitude, so they need a log axis — and a
    stacked bar on a log axis is unreadable, because its segments no longer
    add up. Size is therefore one panel and composition another, the second
    one linear, where a proportion means what it looks like.
    """
    figstyle.use()
    rows = census()
    groups: dict[str, dict[str, int]] = {}
    for r in rows:
        g = groups.setdefault(r["group"] or "unclassified",
                              {c: 0 for c in CALL_ORDER})
        g[r["call"]] = g.get(r["call"], 0) + 1
    order = sorted(groups, key=lambda g: -sum(groups[g].values()))
    totals = [sum(groups[g].values()) for g in order]

    fig, axes = plt.subplots(1, 2, figsize=(figstyle.W_FULL, 3.0),
                             gridspec_kw={"width_ratios": [1.0, 1.25]})
    y = list(range(len(order)))

    ax = axes[0]
    ax.barh(y, totals, height=0.66, color="#52514e", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.7, max(totals) * 3.2)
    for i, n in enumerate(totals):
        ax.text(n * 1.35, i, f"{n:,}", va="center", fontsize=6.2,
                color="#52514e")
    ax.set_xlabel("records (log scale)")
    figstyle.panel(ax, "a", "size of the search space")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")

    ax = axes[1]
    left = [0.0] * len(order)
    for call in CALL_ORDER:
        vals = [100 * groups[g][call] / sum(groups[g].values()) for g in order]
        ax.barh(y, vals, left=left, height=0.66, color=CALL_COLOUR[call],
                label=call, zorder=3)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("composition (% of the group's records)")
    figstyle.panel(ax, "b", "what the records are")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), fontsize=7)

    provenance(fig, "computed — census_v2.tsv")
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    figstyle.save(fig, FIG_DIR / "census_space")
    plt.close(fig)


# ------------------------------------------------------------------ fig 2
def fig_lengths() -> None:
    figstyle.use()
    rows = census()
    fig, axes = plt.subplots(1, 2, figsize=(figstyle.W_FULL, 2.7))

    ax = axes[0]
    # Every record is binned; the axis is clipped at 6,000 aa so the two
    # peaks are resolvable. The legend counts are therefore the full counts,
    # not what happens to be inside the view.
    top = max(int(r["length"]) for r in rows)
    bins = range(0, top + 100, 100)
    over = 0
    for call in CALL_ORDER:
        vals = [int(r["length"]) for r in rows if r["call"] == call]
        over += sum(1 for v in vals if v > 6000)
        ax.hist(vals, bins=bins, color=CALL_COLOUR[call], alpha=0.82,
                label=f"{call} (n={len(vals):,})", zorder=3)
    ax.set_xlim(0, 6000)
    ax.axvspan(2000, 3600, color="#2a78d6", alpha=0.07, zorder=1)
    ax.text(2800, ax.get_ylim()[1] * 0.94, "size band", ha="center",
            fontsize=6, color="#52514e")
    ax.set_xlabel("length (aa)")
    ax.set_ylabel("records")
    ax.legend(frameon=False, fontsize=6.4)
    figstyle.panel(ax, "a", "length was not used to call, and separates "
                            "them anyway")
    if over:
        ax.text(0.985, 0.80, f"{over} records > 6,000 aa\nbinned, off-axis",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.6,
                color="#8a897f")
    figstyle.despine(ax)
    figstyle.hgrid(ax)

    ax = axes[1]
    una = sorted(int(r["length"]) for r in rows if r["call"] == "unassigned")
    called = sorted(int(r["length"]) for r in rows if r["call"] != "unassigned")
    for vals, colour, lab in ((called, "#52514e", "called"),
                              (una, CALL_COLOUR["unassigned"], "unassigned")):
        if not vals:
            continue
        cum = [100 * (i + 1) / len(vals) for i in range(len(vals))]
        ax.plot(vals, cum, color=colour, lw=1.4, label=lab, zorder=3)
    ax.axvline(2000, color="#b3261e", lw=0.8, ls="--", zorder=2)
    ax.set_xscale("log")
    ax.set_xlabel("length (aa, log)")
    ax.set_ylabel("cumulative % of records")
    ax.legend(frameon=False, fontsize=6.4, loc="upper left")
    figstyle.panel(ax, "b", "the unassigned pile is a fragment pile")
    figstyle.despine(ax)
    figstyle.hgrid(ax)

    provenance(fig, "computed — census_v2.tsv")
    fig.tight_layout()
    figstyle.save(fig, FIG_DIR / "census_lengths")
    plt.close(fig)


# ------------------------------------------------------------------ fig 3
def fig_margin() -> None:
    """The sequence check, with D7's no-call band drawn rather than implied.

    Without the band the four sign flips near zero read as failures of the
    architecture call. They are not: the test declines to separate two
    families at that distance, which is why the band exists.
    """
    figstyle.use()
    path = OUT_DIR / "sequence_check.tsv"
    meta_path = OUT_DIR / "sequence_check_meta.json"
    if not path.exists() or not meta_path.exists():
        print("  [s2] no sequence_check.tsv — skipping census_margin")
        return
    meta = json.loads(meta_path.read_text())
    band = meta["no_call_band"]
    ag = meta["agreement"]["full_alignment"]
    rows = [r for r in read_tsv(path) if r["metric"] == "full_alignment"]
    rows.sort(key=lambda r: float(r["margin"]))

    fig, ax = plt.subplots(figsize=(figstyle.W_FULL, 5.6))
    ys = list(range(len(rows)))
    ax.axvspan(-band, band, color="#d6d5cf", alpha=0.55, zorder=1)
    ax.scatter([float(r["margin"]) for r in rows], ys,
               c=[CALL_COLOUR.get(r["arch_call"], "#a9a79e") for r in rows],
               s=18, zorder=3, linewidths=0)
    ax.axvline(0, color="#52514e", lw=0.9, zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r['phylum'] or r['group']}  {r['accession']}"
                        for r in rows], fontsize=5.2)
    ax.set_ylim(-1, len(rows))
    ax.set_xlabel("identity to nearest ITPR bait "
                  "\u2212 identity to nearest RyR bait")
    ax.text(0, len(rows) - 0.2, f"no-call band  \u00b1{band:.2f}",
            ha="center", va="bottom", fontsize=6, color="#52514e")
    ax.set_title(f"architecture call vs sequence: "
                 f"{ag['decisive_agree']}/{ag['decisive_n']} agree where the "
                 f"margin decides; {ag['band_n']} rows inside the band",
                 fontsize=7.6, loc="left", pad=12)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=4,
                          color=CALL_COLOUR[c], label=f"architecture: {c}")
               for c in ("ITPR", "RYR")]
    ax.legend(handles=handles, frameon=False, fontsize=6.6, loc="lower right")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")
    provenance(fig, "computed — sequence_check.tsv")
    fig.tight_layout()
    figstyle.save(fig, FIG_DIR / "census_margin")
    plt.close(fig)


# ------------------------------------------------------------------ fig 4
def fig_growth() -> None:
    """What a domain search adds over a search that knows the family's name.

    The interesting quantity is the *fraction* invisible to a name search,
    so that is the axis; the absolute size of each group is annotated
    rather than plotted, since it is already fig. 1a.
    """
    figstyle.use()
    rows = read_tsv(OUT_DIR / "delta_by_group.tsv")
    rows = [r for r in rows if int(r["v2_records"]) > 0]
    rows.sort(key=lambda r: -int(r["v2_records"]))
    fig, ax = plt.subplots(figsize=(figstyle.W_FULL, 2.8))
    y = list(range(len(rows)))
    pct = [float(r["pct_new"]) for r in rows]
    ax.barh(y, [100 - p for p in pct], height=0.66, color="#184f95", zorder=3,
            label="already found by name")
    ax.barh(y, pct, left=[100 - p for p in pct], height=0.66, color="#86b6ef",
            zorder=3, label="added by the domain enumeration")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['group']}  (n={int(r['v2_records']):,})"
                        for r in rows], fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of the group's census-v2 records")
    # The caveat travels with the figure: census v1 is the app's own search
    # output, mostly a human smoke test, so "added" is near-total outside
    # the vertebrates by construction. Read alone, this panel would
    # overstate what the enumeration proved.
    ax.set_title("census v1 is the app's search bundles, not a family-wide "
                 "harvest — outside\nthe vertebrates it barely searched, so "
                 "these are the searches' scope, not their failure",
                 fontsize=6.4, loc="left", color="#52514e", pad=8)
    for i, (r, p_new) in enumerate(zip(rows, pct)):
        # Keep the annotation legible whichever segment it lands on: white
        # inside the dark bar, grey outside it when that bar is too short.
        inside = (100 - p_new) >= 14
        ax.text(2 if inside else (100 - p_new) + 2, i,
                f"{int(r['in_v1'])} of {int(r['v2_records']):,}",
                va="center", fontsize=5.8,
                color="#ffffff" if inside else "#52514e")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), fontsize=7)
    provenance(fig, "computed — delta_by_group.tsv")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    figstyle.save(fig, FIG_DIR / "census_growth")
    plt.close(fig)


# ------------------------------------------------------------------ fig 5
def fig_lineage() -> None:
    """Where the family is, outside the animals — by phylum, in taxa.

    Records are the wrong unit here: one well-sequenced alga can contribute
    a dozen. The bar is therefore *taxa*, split into those with at least one
    ITPR call and those that put records in the search space without one —
    which is what makes the Streptophyta and Dikarya zeros legible as zeros
    rather than as absent rows.
    """
    figstyle.use()
    rows = [r for r in read_tsv(OUT_DIR / "lineage_calls.tsv")
            if r["group"] not in ("Vertebrata", "Metazoa (non-vertebrate)")]
    rows.sort(key=lambda r: (r["group"], -int(r["itpr_taxa"]),
                             -int(r["taxa"])))
    fig, ax = plt.subplots(figsize=(figstyle.W_FULL, 4.4))
    y = list(range(len(rows)))
    hit = [int(r["itpr_taxa"]) for r in rows]
    rest = [int(r["taxa"]) - int(r["itpr_taxa"]) for r in rows]
    ax.barh(y, hit, height=0.64, color=CALL_COLOUR["ITPR"], zorder=3,
            label="taxa with an ITPR call")
    ax.barh(y, rest, left=hit, height=0.64, color="#d6d5cf", zorder=3,
            label="taxa in the search space with none")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['phylum']}   ({r['group']})" for r in rows],
                       fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("taxa")
    for i, r in enumerate(rows):
        if int(r["itpr_taxa"]) == 0:
            ax.text(int(r["taxa"]) + 0.6, i, "none", va="center",
                    fontsize=5.8, color="#b3261e")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), fontsize=7)
    provenance(fig, "computed — lineage_calls.tsv")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    figstyle.save(fig, FIG_DIR / "census_lineage")
    plt.close(fig)


FIGURES = {"census_space": fig_space, "census_lengths": fig_lengths,
           "census_margin": fig_margin, "census_growth": fig_growth,
           "census_lineage": fig_lineage}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(k)
        return 0
    for name in (args.only or list(FIGURES)):
        print(f"[s2] {name}")
        FIGURES[name]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
