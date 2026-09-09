"""The four S19 figures, from committed S19 tables only (D13, D19).

Three rules this task's figures inherit and one it adds.

* **Binned, never smoothed** (`s5_figures.py`): the contiguity panel's bins
  are placed *on* D4's bar rather than across it, because a window
  straddling the threshold reports an error rate no genome in it has, and
  the threshold is what the panel exists to show.
* **The bar is drawn, not captioned** (`s15_figures.py`): D4's contiguity
  bar and TM-align-style published thresholds belong on the axis.
* **A zero is drawn as a zero** (`s15b_figures.py`): an empty bar reads as
  *not measured*, and several of these zeros are the result.
* And the one S19 needs: the ablation panel is drawn as **delta from the
  full panel**, not as absolute recall, because every panel scores between
  0.58 and 0.85 and on an absolute axis the whole result is one pixel.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                 # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

import figstyle as F                                            # noqa: E402
import s19_drift as DR                                          # noqa: E402
import s19_lib as S                                             # noqa: E402

FIG_DIR = S.OUT_DIR / "figures"
SERIES_COLOUR = {"itpr_present": F.PARALOG["ITPR1"],
                 "ryr_sister": F.GROUP["RYR"]}
SERIES_LABEL = {"itpr_present": "ITPR cells (S15: gene present)",
                "ryr_sister": "RyR sister cell (present in every vertebrate)"}


def _num(v, default=0.0):
    return S.fnum(v, float, default)


# ------------------------------------------------------- 1. the false negative

def fig_contiguity(log=S.log) -> None:
    """The measurement the task rests on, and the threshold it calibrates."""
    F.use()
    bins = S.read_tsv(S.OUT_DIR / "contiguity_bins.tsv")
    floor = S.read_tsv(S.OUT_DIR / "absence_floor.tsv")
    bar = S.contiguity_bar()

    fig, axes = plt.subplots(1, 3, figsize=(F.W_FULL, 2.5))

    ax = axes[0]
    labels = []
    for r in bins:
        if r["bin"] not in labels:
            labels.append(r["bin"])
    x = range(len(labels))
    for series in ("itpr_present", "ryr_sister"):
        ys, los, his = [], [], []
        for lab in labels:
            row = next((r for r in bins
                        if r["bin"] == lab and r["series"] == series), None)
            ys.append(_num(row["rate"]) if row else 0.0)
            los.append(_num(row["wilson_lo"]) if row else 0.0)
            his.append(_num(row["wilson_hi"]) if row else 0.0)
        ax.errorbar(list(x), ys,
                    yerr=[[y - l for y, l in zip(ys, los)],
                          [h - y for y, h in zip(ys, his)]],
                    marker="o", ms=3.2, lw=1.2, capsize=2,
                    color=SERIES_COLOUR[series], label=SERIES_LABEL[series])
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("false-negative rate")
    ax.set_xlabel("contig N50")
    ax.axvline(1.5, color=F.MUTED, lw=0.8, ls="--")
    ax.set_ylim(-0.03, 0.72)
    ax.annotate(f"D4 bar, {bar / 1000:.0f} kb", xy=(1.5, 0.34),
                xytext=(4, 0), textcoords="offset points",
                fontsize=6, color=F.MUTED, va="center")
    ax.legend(loc="upper right", frameon=False, fontsize=5.5,
              handlelength=1.2, borderpad=0.2)
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "a", "misses track the assembly")

    ax = axes[1]
    for series in ("itpr_present", "ryr_sister"):
        mine = sorted((r for r in floor if r["series"] == series),
                      key=lambda r: _num(r["contig_n50_floor"]))
        xs = [max(1e3, _num(r["contig_n50_floor"])) for r in mine]
        ax.plot(xs, [_num(r["fn_rate"]) for r in mine], marker="o", ms=2.6,
                lw=1.2, color=SERIES_COLOUR[series],
                label=SERIES_LABEL[series].split(" (")[0])
    ax.set_ylim(-0.006, 0.175)
    ax.axvline(bar, color=F.MUTED, lw=0.8, ls="--")
    ax.axhline(0.05, color="#b3261e", lw=0.7, ls=":")
    ax.annotate("5 %", xy=(1e3, 0.05), xytext=(1, 2),
                textcoords="offset points", fontsize=6, color="#b3261e")
    ax.set_xscale("log")
    ax.set_xlabel("contig N50 floor applied (bp)")
    ax.set_ylabel("residual false-negative rate")
    ax.legend(loc="upper right", frameon=False, fontsize=5.5,
              handlelength=1.2, borderpad=0.2)
    ax.annotate(f"D4 bar", xy=(bar, 0.16), xytext=(4, 0),
                textcoords="offset points", fontsize=6, color=F.MUTED,
                va="center")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "b", "the a-priori bar, calibrated")

    ax = axes[2]
    mine = sorted((r for r in floor if r["series"] == "itpr_present"),
                  key=lambda r: _num(r["contig_n50_floor"]))
    xs = [max(1e3, _num(r["contig_n50_floor"])) for r in mine]
    ax.plot(xs, [_num(r["frac_genomes_kept"]) for r in mine], marker="o",
            ms=2.6, lw=1.2, color=F.MUTED)
    ax.axvline(bar, color=F.MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("contig N50 floor applied (bp)")
    ax.set_ylabel("fraction of the 309 genomes retained")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "c", "what the floor costs")

    fig.tight_layout()
    F.save(fig, FIG_DIR / "fig_s19_contiguity")
    plt.close(fig)
    log("figure: fig_s19_contiguity")


# --------------------------------------------------------------- 2. the panel

def fig_panel(log=S.log) -> None:
    """Breadth against paralog coverage, and the design rule behind it."""
    F.use()
    recall = S.read_tsv(S.OUT_DIR / "panel_recall.tsv")
    ident = S.read_tsv(S.OUT_DIR / "panel_identity_recall.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 2.9),
                             gridspec_kw={"width_ratios": [1.55, 1]})

    ax = axes[0]
    rows = [r for r in recall if r["panel"] != "full38"]
    rows.sort(key=lambda r: _num(r["delta_vs_full"]))
    names = [r["panel"] for r in rows]
    deltas = [_num(r["delta_vs_full"]) for r in rows]
    colours = ["#b3261e" if d < -10 else "#2a78d6" if d >= 0 else "#a9a79e"
               for d in deltas]
    y = range(len(names))
    ax.barh(list(y), deltas, color=colours, height=0.72)
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{n}  ({r['n_baits']})" for n, r in
                        zip(names, rows)], fontsize=6)
    ax.axvline(0, color=F.INK, lw=0.8)
    # Symmetric log: the ablations that matter differ by 1-2 cells and the
    # one that removes every labelled bait differs by 783. On a linear axis
    # the whole breadth result — which is the finding — is a single pixel;
    # the linear window below 10 keeps those deltas readable while the
    # 783-cell bar stays on the same axis rather than being cropped away.
    ax.set_xscale("symlog", linthresh=10, linscale=0.6)
    ax.set_xlim(-1400, 60)
    ax.set_xticks([-1000, -100, -10, 0, 10])
    ax.set_xticklabels(["-1000", "-100", "-10", "0", "+10"])
    ax.set_xlabel("ITPR cells found, change from the 38-bait panel "
                  "(symmetric log, linear below 10)")
    # The longest bar's label goes *inside* it: outside, it lands on the
    # tick label of its own row.
    for yy, d in zip(y, deltas):
        inside = d < -400
        ax.annotate(f"{d:+.0f}", xy=(d, yy),
                    xytext=(4 if inside or d >= 0 else -4, 0),
                    textcoords="offset points", fontsize=5.5, va="center",
                    ha="left" if inside or d >= 0 else "right",
                    color="white" if inside else F.MUTED)
    F.despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    F.panel(ax, "a", "breadth is nearly free; paralog coverage is not")

    ax = axes[1]
    rows = [r for r in ident if r["identity_bin"] != "no alignment"]
    x = range(len(rows))
    ax.bar(list(x), [_num(r["recall"]) for r in rows], color=F.PARALOG["ITPR3"],
           width=0.7)
    ax.errorbar(list(x), [_num(r["recall"]) for r in rows],
                yerr=[[_num(r["recall"]) - _num(r["wilson_lo"]) for r in rows],
                      [_num(r["wilson_hi"]) - _num(r["recall"])
                       for r in rows]],
                fmt="none", ecolor=F.INK, lw=0.8, capsize=2)
    for xx, r in zip(x, rows):
        ax.annotate(f"n={r['n_measurements']}", xy=(xx, 0.02),
                    fontsize=5.5, ha="center", color="white", rotation=90,
                    va="bottom")
    ax.set_xticks(list(x))
    ax.set_xticklabels([r["identity_bin"] for r in rows], rotation=45,
                       ha="right", fontsize=6)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("recall of one bait, alone")
    ax.set_xlabel("bait-to-target identity")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "b", "one bait suffices at any identity above 0.5")

    fig.tight_layout()
    F.save(fig, FIG_DIR / "fig_s19_panel")
    plt.close(fig)
    log("figure: fig_s19_panel")


# ------------------------------------------------------- 3. per-method yield

def fig_contribution(log=S.log) -> None:
    """Where a profile HMM earns its place, and what no protein DB holds."""
    F.use()
    h2h = S.read_tsv(S.OUT_DIR / "head_to_head.tsv")
    rec = S.read_tsv(S.OUT_DIR / "gene_recovery_by_cell.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 2.7),
                             gridspec_kw={"width_ratios": [1.35, 1]})

    ax = axes[0]
    dbs = [d for d in ("vertebrata", "metazoa_nonvert", "protista_other",
                       "viridiplantae", "fungi")
           if any(r["database"] == d and r["n_pfam"] for r in h2h)]
    bands = ["<1000 aa (profile tail)", "1000-1999 aa (partial)",
             ">=2000 aa (gene-scale)"]
    band_colour = (F.BLUES[1], F.BLUES[3], F.BLUES[5])
    # The **disagreement rate** — the share of the union of the two methods
    # that only one of them returned — rather than stacked counts on a log
    # axis. Stacking on a log axis misreads proportions by construction, and
    # the claim here is a proportion: at gene scale the two methods return
    # the same set, and everything the profile adds is in the tail.
    width = 0.26
    for i, band in enumerate(bands):
        fr, ns = [], []
        for db in dbs:
            r = next((x for x in h2h if x["database"] == db
                      and x["length_band"] == band), None)
            union = _num(r["n_union"]) if r else 0.0
            only = (_num(r["n_hmm_only"]) + _num(r["n_pfam_only"])) if r else 0.0
            fr.append(only / union if union else 0.0)
            ns.append(int(union))
        xs = [j + (i - 1) * width for j in range(len(dbs))]
        ax.bar(xs, fr, width, color=band_colour[i],
               label=band.split(" (")[0])
        for xx, v, n in zip(xs, fr, ns):
            # An exact zero is written as `0`, never left as an absent bar:
            # at gene scale the zero is the result.
            ax.annotate("0" if v == 0 else f"{v:.2f}", xy=(xx, v),
                        xytext=(0, 2.5), textcoords="offset points",
                        fontsize=5, ha="center", color=F.MUTED, rotation=90)
    ax.set_xticks(list(range(len(dbs))))
    ax.set_xticklabels(dbs, rotation=18, ha="right", fontsize=6)
    ax.set_ylim(0, 1.42)
    ax.set_ylabel("share of the union only one method found")
    ax.legend(loc="upper center", frameon=False, fontsize=5.5, ncol=3,
              handlelength=1.1, borderpad=0.2, columnspacing=1.0,
              bbox_to_anchor=(0.5, 1.03))
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "a", "the profile's gain sits in the tail")

    ax = axes[1]
    cells = [r["cell"] for r in rec]
    fr = [_num(r["frac"]) for r in rec]
    lo = [_num(r["wilson_lo"]) for r in rec]
    hi = [_num(r["wilson_hi"]) for r in rec]
    cols = [F.PARALOG.get(c, F.GROUP["RYR"]) for c in cells]
    x = range(len(cells))
    ax.bar(list(x), fr, color=cols, width=0.62)
    ax.errorbar(list(x), fr,
                yerr=[[a - b for a, b in zip(fr, lo)],
                      [a - b for a, b in zip(hi, fr)]],
                fmt="none", ecolor=F.INK, lw=0.8, capsize=2)
    ax.set_xticks(list(x))
    ax.set_xticklabels(cells)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("genes no protein database holds")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "b", "three quarters of the genes are DNA only")

    fig.tight_layout()
    F.save(fig, FIG_DIR / "fig_s19_contribution")
    plt.close(fig)
    log("figure: fig_s19_contribution")


# --------------------------------------------------------------- 4. the drift

def fig_drift(log=S.log) -> None:
    """The rule written to catch drift, against the axis that moves."""
    F.use()
    rounds = S.read_tsv(S.OUT_DIR / "jackhmmer_rounds.tsv")
    val = S.read_tsv(S.OUT_DIR / "kill_criterion_validation.tsv")
    outcome = {r["run"]: r["drifted"] == "1"
               for r in S.read_tsv(S.OUT_DIR / "drift_outcome.tsv")}

    fig, axes = plt.subplots(1, 3, figsize=(F.W_FULL, 2.5))
    runs = []
    for r in rounds:
        if r["run"] not in runs:
            runs.append(r["run"])

    for ax, col, title, letter in (
            (axes[0], "sister_frac",
             "K1's axis: the sister-family share", "a"),
            (axes[1], "offfamily_frac",
             "the axis that actually moves", "b")):
        for run in runs:
            mine = sorted((r for r in rounds if r["run"] == run),
                          key=lambda r: int(r["round"]))
            drifted = outcome.get(run, False)
            ax.plot([int(r["round"]) for r in mine],
                    [_num(r[col]) for r in mine],
                    marker="o", ms=2.4, lw=1.1,
                    color="#b3261e" if drifted else F.PARALOG["ITPR1"],
                    alpha=0.95 if drifted else 0.55)
        ax.set_xlabel("jackhmmer round")
        ax.set_ylabel("share of the included set")
        ax.set_ylim(-0.02, 1.02)
        F.despine(ax)
        F.hgrid(ax)
        F.panel(ax, letter, title)
    axes[1].axhline(DR.DRIFTED_OFFFAMILY_FRAC, color=F.MUTED, lw=0.8, ls="--")
    axes[1].annotate("drifted", xy=(1, DR.DRIFTED_OFFFAMILY_FRAC),
                     xytext=(0, 3), textcoords="offset points", fontsize=6,
                     color=F.MUTED)
    axes[0].plot([], [], color="#b3261e", lw=1.1,
                 label="run that drifted (3)")
    axes[0].plot([], [], color=F.PARALOG["ITPR1"], lw=1.1, alpha=0.55,
                 label="run that did not (4)")
    axes[0].legend(loc="upper right", frameon=False, fontsize=6)

    ax = axes[2]
    # The committed rule names are full sentences, which is right in a table
    # and impossible on an axis; the short forms are a display map here and
    # the table keeps the sentence.
    short = {"K1": "K1\nsister", "K2": "K2\ngrowth", "K3": "K3\nceiling"}
    names = [short.get(r["rule"].split(" (")[0], "proposed\noff-family")
             for r in val]
    x = range(len(names))
    ax.bar([i - 0.18 for i in x], [_num(r["sensitivity"]) for r in val],
           0.34, color=F.PARALOG["ITPR1"], label="sensitivity")
    ax.bar([i + 0.18 for i in x], [_num(r["specificity"]) for r in val],
           0.34, color=F.PARALOG["ITPR3"], label="specificity")
    for i, r in enumerate(val):
        for off, key in ((-0.18, "sensitivity"), (0.18, "specificity")):
            v = _num(r[key])
            ax.annotate("0" if v == 0 else f"{v:.2f}", xy=(i + off, v),
                        xytext=(0, 2), textcoords="offset points",
                        fontsize=5.5, ha="center", color=F.MUTED)
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, fontsize=5.5)
    ax.set_ylim(0, 1.42)
    ax.legend(loc="upper center", frameon=False, fontsize=6, ncol=2,
              bbox_to_anchor=(0.5, 1.02), columnspacing=1.0,
              handlelength=1.2)
    ax.set_ylabel("over 7 runs")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "c", "each kill rule against the outcome")

    fig.tight_layout()
    F.save(fig, FIG_DIR / "fig_s19_drift")
    plt.close(fig)
    log("figure: fig_s19_drift")


def run(log=S.log) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig_contiguity(log)
    fig_panel(log)
    fig_contribution(log)
    fig_drift(log)


if __name__ == "__main__":
    run()
