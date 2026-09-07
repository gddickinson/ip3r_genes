"""The four S9 figures, from committed S9 tables only (D13, D19).

  1. `s9_omega_by_paralog`  per-paralog one-ratio ω with the curated-CDS
                            sensitivity estimate beside it, on a log axis
                            with neutrality *drawn*. At ω ≈ 0.03 a linear
                            axis puts every bar on the floor and the one
                            number a reader wants — how far below 1 — is
                            the one the picture hides.
  2. `s9_dnds_saturation`   pairwise dN against dS per paralog, with the
                            neutral diagonal and the saturation bar. This
                            is the panel that qualifies every ratio in the
                            report, so it is drawn rather than described.
  3. `s9_branch_contrast`   two-ratio background vs foreground ω per
                            paralog clade, and RELAX's k beside it — two
                            different questions, so two panels, not two
                            y-scales on one.
  4. `s9_bs_restarts`       every branch-site restart against its own null.
                            A point below the identity line is a local
                            optimum, which is why model A is restarted by
                            construction here.

Nothing recomputes a likelihood, a ratio or a test: every value is read out
of `omega_table.tsv`, `pairwise_dnds.tsv`, `lrt_table.tsv`,
`bs_restarts.tsv` and `relax_table.tsv`.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s9_figures.py [--only slug]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.lines import Line2D                          # noqa: E402
from matplotlib.patches import Patch                          # noqa: E402

import figstyle as fs                                        # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results" / "selection"
FIGS = OUT / "figures"
PARALOGS = fs.PARALOG_ORDER


def load(name: str) -> list[dict]:
    p = OUT / name
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def by_job(rows: list[dict]) -> dict[str, dict]:
    return {r["job"]: r for r in rows}


# ---- 1. one-ratio ω --------------------------------------------------------

def fig_omega() -> None:
    jobs = by_job(load("omega_table.tsv"))
    full = [f(jobs.get(f"m0_{p}", {}).get("omega")) for p in PARALOGS]
    cur = [f(jobs.get(f"m0_{p}_curated", {}).get("omega")) for p in PARALOGS]
    if not any(v is not None for v in full):
        print("  skip s9_omega_by_paralog — no one-ratio results yet")
        return
    fig, ax = plt.subplots(figsize=(fs.W_HALF * 1.4, 2.7))
    w = 0.36
    for i, p in enumerate(PARALOGS):
        c = fs.PARALOG[p]
        if full[i] is not None:
            ax.bar(i - w / 2, full[i], w, color=c, zorder=3)
        if cur[i] is not None:
            ax.bar(i + w / 2, cur[i], w, color=c, alpha=0.4, zorder=3,
                   hatch="///", edgecolor="white", linewidth=0)
    ax.axhline(1.0, color=fs.INK, lw=1.0, ls="--", zorder=4)
    # Inside the axes, left-aligned: at x = n - 0.5 with ha="right" the
    # label lands on the frame and is clipped away.
    ax.annotate("neutrality (ω = 1)", xy=(-0.42, 1.0), xytext=(0, 3),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=fs.FS_TICK, color=fs.INK)
    ax.set_yscale("log")
    lo = min(v for v in full + cur if v is not None)
    ax.set_ylim(lo / 2.2, 2.6)
    ax.set_xlim(-0.55, len(PARALOGS) - 0.45)
    ax.set_xticks(range(len(PARALOGS)))
    ax.set_xticklabels(PARALOGS)
    ax.set_ylabel("one-ratio ω  (dN/dS)")
    ax.legend(handles=[
        Patch(facecolor=fs.MUTED, edgecolor="none", label="all CDS"),
        Patch(facecolor=fs.MUTED, alpha=0.4, edgecolor="white",
              hatch="///", label="curated CDS only")],
        frameon=False, fontsize=fs.FS_TICK, loc="upper right",
        handlelength=1.4)
    fs.despine(ax)
    fs.hgrid(ax)
    fs.panel(ax, "", "Purifying selection on every paralog")
    fig.tight_layout()
    fs.save(fig, FIGS / "s9_omega_by_paralog")
    plt.close(fig)


# ---- 2. saturation ---------------------------------------------------------

def fig_saturation() -> None:
    pairs = load("pairwise_dnds.tsv")
    if not pairs:
        print("  skip s9_dnds_saturation — no pairwise results yet")
        return
    bar = 1.5
    # Both axes log. On a log-x / linear-y plot the neutral diagonal is not
    # a line at all, and the first draft's dashed "ω = 1" ran off the panel
    # in the first pixel — leaving the *saturation* bar as the only line on
    # the figure, which a reader would read as neutrality. Log-log puts
    # every constant-ω contour back on a straight line.
    fig, axes = plt.subplots(1, len(PARALOGS), figsize=(fs.W_FULL, 2.6),
                             sharex=True, sharey=True)
    xs_lim = (0.02, 200.0)
    for ax, p in zip(axes, PARALOGS):
        rows = [r for r in pairs if r["set"] == p]
        xs = [f(r["dS"], 0.0) for r in rows]
        ys = [f(r["dN"], 0.0) for r in rows]
        ax.scatter(xs, ys, s=6, color=fs.PARALOG[p], alpha=0.5,
                   linewidths=0, zorder=3)
        for om, style, lab in ((1.0, "--", "ω = 1"), (0.05, ":", "ω = 0.05")):
            ax.plot(xs_lim, [xs_lim[0] * om, xs_lim[1] * om],
                    color=fs.INK if om == 1.0 else fs.MUTED,
                    lw=0.9, ls=style, zorder=4)
            del lab
        ax.axvspan(bar, xs_lim[1], color=fs.MUTED, alpha=0.10, zorder=1,
                   linewidth=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*xs_lim)
        ax.set_ylim(0.002, 3.0)
        n_sat = sum(int(r["saturated"]) for r in rows)
        fs.panel(ax, "", f"{p}  ({100 * n_sat / max(1, len(rows)):.0f} % "
                         f"past dS = {bar})")
        ax.set_xlabel("dS")
        fs.despine(ax)
        fs.hgrid(ax)
    axes[0].set_ylabel("dN")
    # The ω contours are labelled in a legend, not on the lines: at these
    # limits both exit through the top edge and an in-line label lands off
    # the panel, which is how the first draft lost them.
    axes[0].legend(handles=[
        Line2D([], [], color=fs.INK, lw=0.9, ls="--", label="ω = 1"),
        Line2D([], [], color=fs.MUTED, lw=0.9, ls=":", label="ω = 0.05")],
        frameon=False, fontsize=fs.FS_TICK, loc="lower right",
        handlelength=1.8)
    axes[-1].annotate("shaded: dS past\nthe saturation bar",
                      xy=(0.97, 0.06), xycoords="axes fraction", ha="right",
                      fontsize=fs.FS_TICK, color=fs.MUTED)
    fig.tight_layout()
    fs.save(fig, FIGS / "s9_dnds_saturation")
    plt.close(fig)


# ---- 3. branch contrast + RELAX -------------------------------------------

def fig_branch_contrast() -> None:
    jobs = by_job(load("omega_table.tsv"))
    relax = {r["paralog"]: r for r in load("relax_table.tsv")}
    have_tr = any(f"two_ratio_{p}" in jobs for p in PARALOGS)
    if not have_tr and not relax:
        print("  skip s9_branch_contrast — no branch or RELAX results yet")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.6))
    w = 0.36
    for i, p in enumerate(PARALOGS):
        oms = (jobs.get(f"two_ratio_{p}", {}).get("omegas") or "").split(";")
        bg = f(oms[0]) if oms and oms[0] else None
        fg = f(oms[1]) if len(oms) > 1 else None
        if bg is not None:
            ax1.bar(i - w / 2, bg, w, color=fs.MUTED, zorder=3)
        if fg is not None:
            ax1.bar(i + w / 2, fg, w, color=fs.PARALOG[p], zorder=3)
    ax1.set_xticks(range(len(PARALOGS)))
    ax1.set_xticklabels(PARALOGS)
    ax1.set_ylabel("ω")
    ax1.legend(handles=[
        Line2D([], [], marker="s", ls="", color=fs.MUTED, label="background"),
        Line2D([], [], marker="s", ls="", color=fs.INK, label="foreground clade")],
        frameon=False, fontsize=fs.FS_TICK, loc="upper left")
    fs.despine(ax1)
    fs.hgrid(ax1)
    fs.panel(ax1, "A", "Two-ratio: each clade against the rest")

    ks = [f(relax.get(p, {}).get("k")) for p in PARALOGS]
    for i, p in enumerate(PARALOGS):
        if ks[i] is not None:
            ax2.bar(i, ks[i], 0.55, color=fs.PARALOG[p], zorder=3)
    ax2.axhline(1.0, color=fs.INK, lw=1.0, ls="--", zorder=4)
    ax2.annotate("k = 1 (no change)", xy=(len(PARALOGS) - 0.6, 1.0),
                 xytext=(0, 4), textcoords="offset points", ha="right",
                 va="bottom", fontsize=fs.FS_TICK, color=fs.INK)
    ax2.set_xticks(range(len(PARALOGS)))
    ax2.set_xticklabels(PARALOGS)
    ax2.set_ylabel("RELAX k")
    if not any(k is not None for k in ks):
        ax2.text(0.5, 0.5, "RELAX not run yet", transform=ax2.transAxes,
                 ha="center", va="center", color=fs.MUTED,
                 fontsize=fs.FS_TICK)
    fs.despine(ax2)
    fs.hgrid(ax2)
    fs.panel(ax2, "B", "Relaxed (k < 1) or intensified (k > 1)")
    fig.tight_layout()
    fs.save(fig, FIGS / "s9_branch_contrast")
    plt.close(fig)


# ---- 4. branch-site restarts ----------------------------------------------

def fig_bs_restarts() -> None:
    rows = load("bs_restarts.tsv")
    if not rows:
        print("  skip s9_bs_restarts — no branch-site results yet")
        return
    fig, ax = plt.subplots(figsize=(fs.W_HALF * 1.5, 2.7))
    for p in PARALOGS:
        sub = [r for r in rows if r["paralog"] == p]
        for r in sub:
            null = f(r["lnL_null"])
            alt = f(r["lnL"])
            if null is None or alt is None:
                continue
            ax.scatter(alt - null, f(r["initial_omega"]),
                       s=52 if r["is_best"] == "1" else 26,
                       color=fs.PARALOG[p],
                       edgecolor=fs.INK if r["is_best"] == "1" else "none",
                       linewidths=0.8, zorder=3)
    ax.axvline(0.0, color=fs.INK, lw=1.0, ls="--", zorder=4)
    ax.set_xlabel("lnL(alternative) − lnL(its own null)")
    ax.set_ylabel("initial ω of the restart")
    ax.annotate("left of the line = local optimum,\nnot a result",
                xy=(0.02, 0.06), xycoords="axes fraction",
                fontsize=fs.FS_TICK, color=fs.MUTED)
    ax.legend(handles=fs.paralog_handles()
              + [Line2D([], [], marker="o", ls="", color="none",
                        markeredgecolor=fs.INK, label="best restart")],
              frameon=False, fontsize=fs.FS_TICK, loc="upper right")
    fs.despine(ax)
    fs.hgrid(ax)
    fs.panel(ax, "", "Branch-site model A: every restart against its null")
    fig.tight_layout()
    fs.save(fig, FIGS / "s9_bs_restarts")
    plt.close(fig)


FIGURES = {"s9_omega_by_paralog": fig_omega,
           "s9_dnds_saturation": fig_saturation,
           "s9_branch_contrast": fig_branch_contrast,
           "s9_bs_restarts": fig_bs_restarts}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(k)
        return
    fs.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for slug, fn in FIGURES.items():
        if args.only and slug not in args.only:
            continue
        print(f"  {slug}")
        fn()


if __name__ == "__main__":
    main()
