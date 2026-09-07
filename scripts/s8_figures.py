"""The four S8 figures, from committed S8 tables only (D13, D19).

  1. `synteny_pair_classes`  mean Jaccard per pair class beside its own
                             matched random-window control -- the figure the
                             whole task rests on, because a Jaccard means
                             nothing without the null it beat.
  2. `synteny_paralogon`     the three human neighbourhoods drawn as gene
                             tracks, with the ohnologous flank families the
                             root-key test recovered linked between them.
  3. `synteny_caller`        the consensus caller's two distributions --
                             annotation-confirmed loci against random
                             windows -- and the threshold sweep that chose
                             its operating point.
  4. `synteny_clade_decay`   within-paralog Jaccard split into same-class and
                             cross-class pairs, which is where the ITPR3
                             neighbourhood's decay lives.

Nothing here recomputes a Jaccard, a consensus or a call: every number is
read out of `pair_stats.tsv`, `paralogon_shared.tsv`, `flanks.tsv`,
`caller_calibration.tsv`, `caller_null.tsv` and `caller_frac_sweep.tsv`.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s8_figures.py [--only slug]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.patches import Patch, Rectangle              # noqa: E402

import figstyle as fs                                        # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / "synteny"
FIGS = OUT / "figures"
HUMAN = "GCF_000001405.40"
CLASSES = ("ITPR1", "ITPR2", "ITPR3", "RYR")
COLOUR = {**fs.PARALOG, "RYR": fs.GROUP["RYR"]}


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stats(window="fixed10", key="relaxed"):
    return [r for r in read_tsv(OUT / "pair_stats.tsv")
            if r["window"] == window and r["key"] == key]


# --------------------------------------------------- 1. pair classes

def fig_pair_classes():
    rows = {r["pair_class"]: r for r in stats() if r["stratum"] == "all"}
    order = [f"within_{c}" for c in CLASSES] + [
        k for k in rows if k.startswith("cross_")]
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 3.0))
    y = range(len(order))
    obs = [float(rows[k]["mean_j"]) for k in order]
    ctlv = [float(rows[k]["control_mean_j"]) for k in order]
    cols = [COLOUR.get(k.replace("within_", ""), "#8a897f") for k in order]
    ax.barh(list(y), obs, height=0.62, color=cols, zorder=3)
    ax.scatter(ctlv, list(y), s=22, marker="D", color=fs.INK, zorder=4,
               label="matched random-window control")
    for i, k in enumerate(order):
        r = rows[k]
        n = int(r["n_pairs"])
        txt = f"  {float(r['mean_j']):.3f}   n={n:,}"
        if float(r["mean_j"]) > float(r["control_mean_j"]):
            txt += f"   {float(r['ratio']):.0f}x null"
        ax.text(max(obs) * 0.012 + float(r["mean_j"]), i, txt, va="center",
                fontsize=fs.FS_TICK, color=fs.INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([k.replace("within_", "within ")
                        .replace("cross_", "").replace("_vs_", " vs ")
                        for k in order], fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, max(obs) * 1.45)
    ax.set_xlabel("mean Jaccard of flanking-gene symbol sets")
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    ax.legend(loc="lower right", frameon=False, fontsize=fs.FS_TICK)
    fs.panel(ax, "", "Flanking-gene sharing, one locus per species, "
                     "+/-10 coding genes")
    fig.tight_layout()
    return fs.save(fig, FIGS / "synteny_pair_classes")


# ----------------------------------------------------- 2. paralogon

def _human_track(cell: str) -> list[tuple[str, str, int]]:
    """(symbol, root, signed rank) for the human locus's fixed10 flanks."""
    out = []
    for r in read_tsv(OUT / "flanks.tsv"):
        if (r["accession"] != HUMAN or r["window"] != "fixed10"
                or r["cell"] != cell):
            continue
        rank = int(r["rank"]) * (-1 if r["side"] == "up" else 1)
        out.append((r["symbol"], r["root_key"], rank))
    return sorted(out, key=lambda t: t[2])


def fig_paralogon():
    shared = [r for r in read_tsv(OUT / "paralogon_shared.tsv")
              if r["pair"].endswith("@fixed10")]
    roots = {r["key"]: r for r in shared}
    cells = ("ITPR1", "ITPR2", "ITPR3")
    fig, axes = plt.subplots(2, 1, figsize=(fs.W_FULL, 4.2),
                             gridspec_kw={"height_ratios": [2.1, 1.0]})
    ax = axes[0]
    pos = {c: 2 - i for i, c in enumerate(cells)}
    xy = {}
    for c in cells:
        track = _human_track(c)
        yv = pos[c]
        ax.plot([-10.6, 10.6], [yv, yv], color=fs.GRID, lw=0.8, zorder=1)
        for sym, root, rank in track:
            hit = root in roots
            ax.add_patch(Rectangle((rank - 0.42, yv - 0.16), 0.84, 0.32,
                                   facecolor=COLOUR[c] if hit else "#e2e0da",
                                   edgecolor=fs.INK if hit else "none",
                                   lw=0.6, zorder=3))
            if hit:
                xy[(c, root)] = (rank, yv)
                # the connectors run downward from the top track, so a label
                # above a lower box lands on the line that reaches it
                above = c == cells[0]
                ax.text(rank, yv + (0.28 if above else -0.28), sym,
                        ha="center", va="bottom" if above else "top",
                        fontsize=fs.FS_TICK - 1, color=fs.INK, zorder=5,
                        bbox=dict(boxstyle="square,pad=0.12", fc="white",
                                  ec="none"))
        ax.add_patch(Rectangle((-0.42, yv - 0.22), 0.84, 0.44,
                               facecolor=fs.INK, edgecolor="none", zorder=3))
        ax.text(-11.2, yv, c, ha="right", va="center", fontsize=fs.FS_TICK,
                color=COLOUR[c], fontweight="bold")
    for root in roots:
        pts = [xy[k] for k in xy if k[1] == root]
        if len(pts) == 2:
            (x1, y1), (x2, y2) = pts
            ax.annotate("", xy=(x2, y2 - 0.18 if y2 > y1 else y2 + 0.18),
                        xytext=(x1, y1 - 0.18 if y1 > y2 else y1 + 0.18),
                        arrowprops=dict(arrowstyle="-", color=fs.INK,
                                        lw=0.9, linestyle=(0, (3, 2))),
                        zorder=2)
    ax.set_xlim(-11.6, 11.4)
    ax.set_ylim(-0.55, 2.75)
    ax.set_yticks([])
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.set_xticklabels(["-10", "-5", "gene", "+5", "+10"])
    ax.set_xlabel("flanking coding genes, human")
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "a", "Human ITPR neighbourhoods; shaded = a gene family "
                      "shared with another ITPR locus")

    ax = axes[1]
    if shared:
        rows_ = sorted(shared, key=lambda r: -float(r["enrichment"]))
        y, h, seen = list(range(len(rows_))), 0.26, []
        for i, r in enumerate(rows_):
            a, b = r["pair"].split("@")[0].split("_vs_")
            seen += [a, b]
            # each bar wears its own paralog's colour, so the reader does not
            # have to remember which side of the pair is "first"
            ax.barh(i + h, float(r["frac_a"]), height=h, color=COLOUR[a],
                    zorder=3)
            ax.barh(i, float(r["frac_b"]), height=h, color=COLOUR[b], zorder=3)
            ax.barh(i - h, float(r["background"]), height=h, color="#a9a79e",
                    zorder=3)
            ax.text(float(r["frac_a"]) + 0.01, i + h, a, va="center",
                    fontsize=fs.FS_TICK - 1, color=fs.INK)
            ax.text(float(r["frac_b"]) + 0.01, i, b, va="center",
                    fontsize=fs.FS_TICK - 1, color=fs.INK)
            ax.text(float(r["background"]) + 0.01, i - h,
                    f"random windows ({float(r['background']):.3f})",
                    va="center", fontsize=fs.FS_TICK - 1, color=fs.INK)
        ax.set_yticks(y)
        ax.set_yticklabels([r["key"] for r in rows_], fontsize=fs.FS_TICK)
        ax.invert_yaxis()
        ax.set_xlabel("fraction of species whose locus carries the family")
        ax.set_xlim(0, 1.0)
    else:
        ax.text(0.5, 0.5, "no shared flank family passed the prevalence bar",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    fs.panel(ax, "b", "Prevalence of each shared family (root key) against "
                      "its random-window background")
    fig.tight_layout()
    return fs.save(fig, FIGS / "synteny_paralogon")


# -------------------------------------------------------- 3. caller

def fig_caller():
    calib = read_tsv(OUT / "caller_calibration.tsv")
    null = read_tsv(OUT / "caller_null.tsv")
    sweep = read_tsv(OUT / "caller_frac_sweep.tsv")
    meta = json.loads((OUT / "synteny_stats.json").read_text())
    picked = meta["caller"]["consensus_frac"]

    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.6))
    ax = axes[0]
    top = max([int(r["best_score"]) for r in calib + null] + [1])
    bins = list(range(top + 2))
    cc = Counter(int(r["best_score"]) for r in calib)
    nc = Counter(int(r["best_score"]) for r in null)
    n_c, n_n = max(len(calib), 1), max(len(null), 1)
    ax.bar([b - 0.2 for b in bins], [cc.get(b, 0) / n_c for b in bins],
           width=0.4, color=fs.PARALOG["ITPR1"], zorder=3,
           label=f"annotation-confirmed loci (n={len(calib)})")
    ax.bar([b + 0.2 for b in bins], [nc.get(b, 0) / n_n for b in bins],
           width=0.4, color="#a9a79e", zorder=3,
           label=f"random windows (n={len(null)})")
    max_null = max([int(r["best_score"]) for r in null] + [0])
    ax.axvline(max_null + 0.5, color=fs.INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
    ax.text(max_null + 0.65, 0.62, "highest score any\nrandom window reached",
            fontsize=fs.FS_TICK - 1, color=fs.INK, va="top")
    ax.set_xlabel("keys shared with the best paralog consensus")
    ax.set_ylabel("fraction of loci")
    ax.set_xticks(bins)
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(loc="upper center", frameon=False, fontsize=fs.FS_TICK - 1)
    fs.panel(ax, "a", "Consensus overlap")

    ax = axes[1]
    x = [float(r["consensus_frac"]) for r in sweep]
    ax.plot(x, [float(r["call_rate"]) for r in sweep], "-o", ms=3.5,
            color=fs.PARALOG["ITPR1"], zorder=3, label="call rate")
    ax.plot(x, [float(r["false_call_rate"]) for r in sweep], "-o", ms=3.5,
            color=fs.PARALOG["ITPR2"], zorder=3, label="false-call rate")
    ax.plot(x, [float(r["margin"]) for r in sweep], "-o", ms=3.5,
            color=fs.INK, zorder=3, label="difference (what is maximised)")
    ax.axvline(picked, color=fs.GRID, lw=6, zorder=1)
    ax.text(picked, 0.94, f"chosen\n{picked:g}", ha="center", va="top",
            fontsize=fs.FS_TICK - 1, color=fs.INK)
    ax.set_xlabel("consensus threshold (fraction of species)")
    ax.set_ylabel("rate")
    ax.set_ylim(0, 1.0)
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(loc="center left", frameon=False, fontsize=fs.FS_TICK - 1)
    fs.panel(ax, "b", "Operating point, measured")
    fig.tight_layout()
    return fs.save(fig, FIGS / "synteny_caller")


# --------------------------------------------------- 4. clade decay

def fig_clade_decay():
    rows = stats()
    by = defaultdict(dict)
    for r in rows:
        if r["pair_class"].startswith("within_"):
            by[r["pair_class"].replace("within_", "")][r["stratum"]] = r
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.7))
    w = 0.36
    xs = list(range(len(CLASSES)))
    for k, (stratum, alpha, lab) in enumerate(
            (("same_vclass", 1.0, "same vertebrate class"),
             ("cross_vclass", 0.45, "different classes"))):
        vals = [float(by[c][stratum]["mean_j"]) for c in CLASSES]
        ax.bar([x + (k - 0.5) * w for x in xs], vals, width=w,
               color=[COLOUR[c] for c in CLASSES], alpha=alpha, zorder=3,
               label=lab, edgecolor="none")
        for x, c, v in zip(xs, CLASSES, vals):
            ax.text(x + (k - 0.5) * w, v + 0.008, f"{v:.2f}", ha="center",
                    va="bottom", fontsize=fs.FS_TICK - 1, color=fs.INK)
    ctl = max(float(by[c]["all"]["control_mean_j"]) for c in CLASSES)
    ax.axhline(ctl, color=fs.INK, lw=0.9, ls=(0, (3, 2)), zorder=4)
    ax.text(len(CLASSES) - 0.5, ctl + 0.006, "random-window control",
            ha="right", va="bottom", fontsize=fs.FS_TICK - 1, color=fs.INK)
    ax.set_xticks(xs)
    ax.set_xticklabels(list(CLASSES))
    ax.set_ylabel("mean Jaccard")
    handles = [Patch(facecolor="#8a897f", alpha=1.0,
                     label="same vertebrate class"),
               Patch(facecolor="#8a897f", alpha=0.45,
                     label="different classes")]
    ax.legend(handles=handles, loc="upper right", frameon=False,
              fontsize=fs.FS_TICK)
    fs.despine(ax)
    fs.hgrid(ax)
    fs.panel(ax, "", "Neighbourhood conservation within and across "
                     "vertebrate classes")
    fig.tight_layout()
    return fs.save(fig, FIGS / "synteny_clade_decay")


FIGURES = {
    "synteny_pair_classes": fig_pair_classes,
    "synteny_paralogon": fig_paralogon,
    "synteny_caller": fig_caller,
    "synteny_clade_decay": fig_clade_decay,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(k)
        return
    fs.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for slug in (args.only or list(FIGURES)):
        paths = FIGURES[slug]()
        print(f"[s8-fig] {slug} -> {', '.join(p.name for p in paths)}")


if __name__ == "__main__":
    main()
