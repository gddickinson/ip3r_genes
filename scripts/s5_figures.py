"""s5_figures.py — the S5 ledger figures, from the committed tables only.

D13/D19 applied to figures: every panel is drawn from `results/genome_ledger/`
and `results/census_v4/`, never from a sweep summary or a recomputed number,
so a figure cannot disagree with the report beside it.

Four figures:

  ledger_status        what the sweep found, per paralog, with the RyR
                       positive control drawn beside them rather than
                       described — the control's own bar is the reader's
                       check that a red ITPR bar means biology
  ledger_by_class      the same, split by vertebrate class: where the
                       evidence is thin is a taxonomic fact, not a global one
  contiguity_confound  the result S5a found and this tests at scale — status
                       against assembly contiguity, with D4's bar drawn, and
                       the per-paralog recovery either side of it
  copy_number          loci per genome per paralog, which is where the
                       teleost 3R duplicates and any real expansion show up

Usage:
  python scripts/s5_figures.py [--only ledger_status]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402

import figstyle as fs                                        # noqa: E402
import s5_sweep_lib as lib                                   # noqa: E402

LEDGER = PROJECT_ROOT / "results" / "genome_ledger"
FIGDIR = LEDGER / "figures"
ALL_CELLS = (*lib.CLASSES, lib.CONTROL_CLASS)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def _statuses_present(rows: list[dict]) -> list[str]:
    seen = {r["status"] for r in rows}
    ordered = [s for s in fs.STATUS_ORDER if s in seen]
    return ordered + sorted(seen - set(ordered))


def _stack(ax, categories: list[str], counts: dict, statuses: list[str],
           horizontal: bool = False) -> None:
    base = [0.0] * len(categories)
    for st in statuses:
        vals = [counts.get((c, st), 0) for c in categories]
        colour = fs.STATUS.get(st, "#8a897f")
        if horizontal:
            ax.barh(categories, vals, left=base, color=colour,
                    label=fs.STATUS_LABEL.get(st, st), height=0.72)
        else:
            ax.bar(categories, vals, bottom=base, color=colour,
                   label=fs.STATUS_LABEL.get(st, st), width=0.72)
        base = [b + v for b, v in zip(base, vals)]


def fig_ledger_status(ledger: list[dict]) -> None:
    statuses = _statuses_present(ledger)
    counts = Counter((r["class"], r["status"]) for r in ledger)
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 3.0))
    _stack(ax, list(ALL_CELLS), counts, statuses)
    n_genomes = len({r["accession"] for r in ledger})
    ax.set_ylabel(f"genomes (n = {n_genomes})")
    ax.set_title("What the genomic sweep found, per paralog", loc="left")
    ax.axvline(2.5, color="#52514e", lw=0.8, ls=":")
    ax.text(3.0, ax.get_ylim()[1] * 0.97, "positive\ncontrol", ha="center",
            va="top", fontsize=6.5, color="#52514e")
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=6.2, frameon=False, ncol=2, loc="center left",
              bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "ledger_status"))
    plt.close(fig)


def fig_ledger_by_class(ledger: list[dict]) -> None:
    statuses = _statuses_present(ledger)
    classes = [c for c, _ in Counter(
        r["vclass"] for r in ledger if r["vclass"]).most_common()]
    classes = [c for c in classes
               if sum(1 for r in ledger if r["vclass"] == c) >= 8]
    fig, axes = plt.subplots(1, len(lib.CLASSES),
                             figsize=(fs.W_FULL, 0.34 * len(classes) + 1.4),
                             sharey=True)
    for ax, cell in zip(axes, lib.CLASSES):
        rows = [r for r in ledger if r["class"] == cell]
        counts = Counter((r["vclass"], r["status"]) for r in rows)
        # fraction, so classes of very different size stay comparable
        totals = {c: sum(counts.get((c, s), 0) for s in statuses)
                  for c in classes}
        frac = {(c, s): (counts.get((c, s), 0) / totals[c] if totals[c] else 0)
                for c in classes for s in statuses}
        _stack(ax, classes, frac, statuses, horizontal=True)
        ax.set_xlim(0, 1)
        ax.set_title(cell, loc="left", fontsize=8)
        ax.set_xlabel("fraction of genomes")
        fs.despine(ax, keep=("bottom",))
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(len(classes)))
    axes[0].set_yticklabels([f"{c} ({sum(1 for r in ledger if r['vclass'] == c and r['class'] == lib.CLASSES[0])})"
                             for c in classes], fontsize=6.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=fs.STATUS.get(s, "#8a897f"))
               for s in statuses]
    fig.legend(handles, [fs.STATUS_LABEL.get(s, s) for s in statuses],
               fontsize=6.2, frameon=False, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Evidence for each paralog, by vertebrate class",
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fs.save(fig, str(FIGDIR / "ledger_by_class"))
    plt.close(fig)


def fig_contiguity(ledger: list[dict]) -> None:
    """S5a's confounder, tested at scale: recovery against contiguity.

    Binned, not smoothed. A sliding window averages the x of its members, so
    a window straddling D4's bar reports a contiguity no genome in it has and
    draws the curve *through* the threshold the panel exists to show. Fixed
    bins with the bar on a bin edge cannot do that, and carrying n per bin
    keeps the reader from reading a 2-genome bin as a rate.
    """
    bar = lib.itpr_span_stats()["median"]
    edges = [0, 50_000, bar, 1_000_000, 10_000_000, float("inf")]
    labels = ["<50 kb", f"50 kb-\n{bar / 1000:.0f} kb",
              f"{bar / 1000:.0f} kb-\n1 Mb", "1-10 Mb", ">10 Mb"]

    def bin_of(n50: int) -> int:
        for i in range(len(edges) - 1):
            if edges[i] <= n50 < edges[i + 1]:
                return i
        return len(edges) - 2

    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.9))

    ax = axes[0]
    width = 0.8 / len(lib.CLASSES)
    n_per_bin = Counter(bin_of(int(r["contig_n50"] or 0))
                        for r in ledger if r["class"] == lib.CLASSES[0])
    for k, cell in enumerate(lib.CLASSES):
        rows = [r for r in ledger if r["class"] == cell and r["contig_n50"]]
        tot: Counter = Counter()
        hit: Counter = Counter()
        for r in rows:
            b = bin_of(int(r["contig_n50"]))
            tot[b] += 1
            hit[b] += r["status"].startswith("found")
        xs = [i + (k - 1) * width for i in range(len(labels))]
        ys = [hit[i] / tot[i] if tot[i] else 0.0 for i in range(len(labels))]
        ax.bar(xs, ys, width=width, color=fs.PARALOG[cell], label=cell)
    for i in range(len(labels)):
        if n_per_bin.get(i):
            ax.text(i, 1.04, f"n={n_per_bin[i]}", ha="center", fontsize=5.8,
                    color="#52514e")
    ax.axvline(1.5, color="#b3261e", lw=1.1, ls="--")
    ax.text(1.42, 0.60, "D4 bar — below this a\ncontig cannot hold the gene",
            fontsize=6.0, color="#b3261e", va="center", ha="right")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=6.3)
    ax.set_xlabel("contig N50")
    ax.set_ylabel("fraction found")
    ax.set_ylim(0, 1.12)
    fs.panel(ax, "a", "Recovery tracks assembly contiguity")
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=6.2, frameon=False, loc="upper center", ncol=3,
              bbox_to_anchor=(0.5, -0.22))

    ax = axes[1]
    spans = {}
    for cell in lib.CLASSES:
        rows = [r for r in ledger if r["class"] == cell]
        above = [r for r in rows if r["contig_spans_gene"] == "1"]
        below = [r for r in rows if r["contig_spans_gene"] != "1"]
        spans[cell] = (
            sum(1 for r in above if r["status"].startswith("found")) / max(1, len(above)),
            sum(1 for r in below if r["status"].startswith("found")) / max(1, len(below)),
            len(above), len(below))
    x = list(range(len(lib.CLASSES)))
    ax.bar([i - 0.19 for i in x], [spans[c][0] for c in lib.CLASSES],
           width=0.36, color=[fs.PARALOG[c] for c in lib.CLASSES])
    ax.bar([i + 0.19 for i in x], [spans[c][1] for c in lib.CLASSES],
           width=0.36, color=[fs.PARALOG[c] for c in lib.CLASSES],
           alpha=0.42, hatch="///", edgecolor="white")
    for i, c in enumerate(lib.CLASSES):
        ax.text(i - 0.19, spans[c][0] + 0.02, f"{spans[c][0]:.0%}",
                ha="center", fontsize=6.0)
        ax.text(i + 0.19, spans[c][1] + 0.02, f"{spans[c][1]:.0%}",
                ha="center", fontsize=6.0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n{lib.itpr_span_stats() and ''}".strip()
                        for c in lib.CLASSES])
    ax.set_xticklabels(list(lib.CLASSES))
    ax.set_ylabel("fraction found")
    ax.set_ylim(0, 1.18)
    n_above, n_below = spans[lib.CLASSES[0]][2], spans[lib.CLASSES[0]][3]
    fs.panel(ax, "b", f"Either side of the bar "
                      f"({n_above} vs {n_below} genomes)")
    fs.despine(ax)
    fs.hgrid(ax)
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#52514e"),
               plt.Rectangle((0, 0), 1, 1, fc="#52514e", alpha=0.42,
                             hatch="///", ec="white")]
    ax.legend(handles, ["contig spans a gene", "contig does not"],
              fontsize=6.2, frameon=False, loc="upper center", ncol=2,
              bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "contiguity_confound"))
    plt.close(fig)


def fig_copy_number(ledger: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.6))
    found = [r for r in ledger if r["status"].startswith("found")]
    maxn = max((int(r["n_loci"]) for r in found), default=1)
    bins = list(range(1, min(maxn, 8) + 1))
    width = 0.8 / len(ALL_CELLS)
    for k, cell in enumerate(ALL_CELLS):
        rows = [r for r in found if r["class"] == cell]
        tally = Counter(min(int(r["n_loci"]), 8) for r in rows)
        total = max(1, len(rows))
        colour = (fs.PARALOG[cell] if cell in fs.PARALOG
                  else fs.GROUP["RYR"])
        ax.bar([b + (k - 1.5) * width for b in bins],
               [tally.get(b, 0) / total for b in bins],
               width=width, color=colour,
               label=cell + (" (control)" if cell == lib.CONTROL_CLASS else ""))
    ax.set_xticks(bins)
    ax.set_xticklabels([str(b) if b < 8 else "8+" for b in bins])
    ax.set_xlabel("loci found in the genome")
    ax.set_ylabel("fraction of genomes where found")
    ax.set_title("Copy number per genome", loc="left")
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=6.5, frameon=False, ncol=4)
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "copy_number"))
    plt.close(fig)


FIGURES = {
    "ledger_status": fig_ledger_status,
    "ledger_by_class": fig_ledger_by_class,
    "contiguity_confound": fig_contiguity,
    "copy_number": fig_copy_number,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(f"  {k}")
        return

    fs.use()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    ledger = read_tsv(LEDGER / "genome_ledger.tsv")
    for name in (args.only or list(FIGURES)):
        FIGURES[name](ledger)
        print(f"  drew {name}")
    print(f"wrote {FIGDIR}")


if __name__ == "__main__":
    main()
