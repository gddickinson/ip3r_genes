"""The four S11 figures, from committed S11 tables only (D13, D19).

Two decisions worth stating.

**The calibration figure draws the bars, it does not describe them.** TM-align
has two published interpretation bars — 0.17, below which a pair is
indistinguishable from two random structures, and 0.50, above which they
share a fold — and the whole point of a calibration panel is that a reader
can see where the negative controls fall relative to them. A panel that
printed the numbers in a caption instead would be asking to be believed.

**The confidence figure is per domain, never per model.** A mean pLDDT over
a 2,700-residue multi-domain channel averages a well-predicted β-trefoil
with hundreds of residues of linker; the resulting single number is high
enough to look reassuring and says nothing about the part any claim rests
on. The domain rows are ordered along the subunit, so the figure reads
N-terminus to pore.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.lines import Line2D                            # noqa: E402
from matplotlib.patches import Patch                           # noqa: E402

import figstyle as fs                                          # noqa: E402
import s11_tables as tb                                        # noqa: E402
from s11_tmalign import TM_FOLD, TM_RANDOM                     # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "structures"
FIGS = OUT / "figures"

#: The AFDB outcome scale, drawn from the project's validated evidence ramp
#: rather than from fresh colours: a usable model is the best evidence, an
#: isoform model the worst thing that still looks like evidence.
AFDB_STATUS = [
    ("full", "modelled, full length", fs.STATUS["found_annotated"]),
    ("partial", "modelled, partial", fs.STATUS["found_no_annotation"]),
    ("isoform", "isoform model only", fs.STATUS["tblastn_trace_ambiguous"]),
    ("absent", "no model", fs.STATUS["absent"]),
]


def _load(name: str) -> list[dict]:
    path = OUT / name
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _f(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _short(label: str) -> str:
    """Domain labels for an axis: drop the review's LaTeX and the
    parenthetical gloss. Truncating instead — as the first version did —
    turned "IP$_3$-binding core (β-trefoil)" into "IP3-binding core (b-trefoi",
    which reads as a typo rather than as a shortened label."""
    out = label.replace("$_3$", "3")
    if "(" in out:
        out = out.split("(")[0].strip()
    return out


#: Controls take the *lighter* neutral. `fs.MUTED` and
#: `fs.GROUP["invert_metazoa"]` are the same hex, so using MUTED for the
#: controls put two identical swatches in the legend under two different
#: names — "ITPR, no paralog" and "negative control" — which is worse than
#: no legend.
CONTROL_COLOUR = fs.GROUP["plant"]


def _role_colour(row: dict) -> str:
    if row.get("role") == "control":
        return CONTROL_COLOUR
    para = (row.get("paralog") or "").upper()
    if para in fs.PARALOG:
        return fs.PARALOG[para]
    if row.get("call") == "RYR" or para.startswith("RYR"):
        return fs.GROUP["RYR"]
    return fs.GROUP["invert_metazoa"]


# --------------------------------------------------------------------------
# Figure 1 — what AlphaFold DB holds
# --------------------------------------------------------------------------

def fig_afdb_coverage() -> None:
    rows = [r for r in _load("afdb_coverage_summary.tsv")
            if r["scope"] == "census" and r["group"] != "ALL"]
    if not rows:
        return
    rows.sort(key=lambda r: -int(r["n_applicable"] or 0))
    groups = [r["group"] for r in rows]
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(fs.W_FULL, 3.2), width_ratios=[1.35, 1.0])

    left = [0.0] * len(rows)
    for key, label, colour in AFDB_STATUS:
        vals = []
        for i, r in enumerate(rows):
            denom = int(r["n_applicable"] or 0) or 1
            vals.append(100.0 * int(r[key] or 0) / denom)
        ax.barh(range(len(rows)), vals, left=left, color=colour,
                edgecolor="white", linewidth=0.4, height=0.72, label=label)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{g}  ({r['n_applicable']})"
                        for g, r in zip(groups, rows)])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of UniProt-shaped census records")
    fs.panel(ax, "a", "AlphaFold DB coverage of the ITPR census")
    fs.despine(ax)
    # Legend below the axes: at eight groups the bars run the full width,
    # so an inset legend sits on top of the shortest group's data.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), frameon=False,
              fontsize=6, ncol=2, handlelength=1.2, columnspacing=1.0)

    # Panel b — the length story. AFDB's monomer pipeline stops short of a
    # full-length receptor, so coverage is not uniform across the family;
    # it is a function of how long the record is.
    census = _load("afdb_coverage_census.tsv")
    modelled, missing = [], []
    for r in census:
        if r["afdb_status"] == "not_applicable":
            continue
        ln = _f(r["length"])
        if not ln:
            continue
        (modelled if _f(r["model_coverage"]) >= 0.95 else missing).append(ln)
    bins = list(range(0, 4001, 200))
    ax2.hist([missing, modelled], bins=bins, stacked=True,
             color=[fs.STATUS["absent"], fs.STATUS["found_annotated"]],
             edgecolor="white", linewidth=0.3,
             label=["no usable model", "usable model"])
    ax2.set_xlabel("record length (aa)")
    ax2.set_ylabel("census records")
    ax2.set_xlim(0, 4000)
    fs.panel(ax2, "b", "coverage against record length")
    fs.despine(ax2)
    fs.hgrid(ax2)
    ax2.legend(loc="upper right", frameon=False, fontsize=6)
    fs.save(fig, FIGS / "s11_afdb_coverage")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 2 — the calibration figure the brief asks for
# --------------------------------------------------------------------------

def fig_tm_calibration() -> None:
    rows = _load("tm_vs_reference.tsv")
    if not rows:
        return
    # Every panel member is plotted, references included. `_best_against`
    # already excludes a structure's pair with itself, so a reference's
    # position is its score against the *other* references and is not
    # circular; dropping the IP3R references while keeping the state panel
    # — which is also experimental references — was an inconsistency that
    # left the ITPR1 and ITPR2 anchors off their own calibration figure.
    rows.sort(key=lambda r: (r["role"] != "control", -_f(r["best_itpr"])))
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(fs.W_FULL, 3.4), width_ratios=[1.0, 1.15])

    xs = [_f(r["best_itpr"]) for r in rows]
    ys = [_f(r["best_ryr"]) for r in rows]
    for r, x, y in zip(rows, xs, ys):
        marker = "s" if r["role"] == "control" else (
            "^" if r["role"] in ("reference", "state_panel") else "o")
        called = r.get("structural_call") in ("ITPR", "RYR")
        colour = _role_colour(r)
        # Open markers are the structures the fold test **declined**. The
        # gate is the whole point of the panel, so a reader has to be able
        # to see which side of it each structure fell — otherwise the three
        # negative controls and the six partial models are indistinguishable
        # from the calls, which is the confusion the gate exists to prevent.
        ax.scatter(x, y, s=28, marker=marker,
                   facecolor=colour if called else "none",
                   edgecolor="white" if called else colour,
                   linewidth=0.5 if called else 1.0, zorder=3)
    lim = 0.95
    ax.plot([0, lim], [0, lim], color=fs.MUTED, lw=0.7, ls=":", zorder=1)
    for bar, label in ((TM_RANDOM, "random"), (TM_FOLD, "same fold")):
        ax.axvline(bar, color=fs.MUTED, lw=0.7, ls="--", zorder=1)
        ax.axhline(bar, color=fs.MUTED, lw=0.7, ls="--", zorder=1)
        # Bar labels at the *foot* of their line: the upper-left corner is
        # the only empty region left once the calls cluster top-right, so
        # that is where the legend has to go.
        ax.text(bar + 0.012, 0.012, label, fontsize=6, color=fs.MUTED,
                va="bottom", ha="left", rotation=90)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("best TM-score vs an IP$_3$R reference")
    ax.set_ylabel("best TM-score vs a RyR reference")
    fs.panel(ax, "a", "the family call, structurally")
    fs.despine(ax)
    ax.legend(handles=[
        Patch(facecolor=fs.PARALOG["ITPR1"], label="ITPR1"),
        Patch(facecolor=fs.PARALOG["ITPR2"], label="ITPR2"),
        Patch(facecolor=fs.PARALOG["ITPR3"], label="ITPR3"),
        Patch(facecolor=fs.GROUP["invert_metazoa"], label="ITPR, no paralog"),
        Patch(facecolor=fs.GROUP["RYR"], label="RyR reference"),
        Patch(facecolor=CONTROL_COLOUR, label="negative control"),
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor=fs.INK, markersize=4.5,
               label="declined (below the fold bar)"),
    ], loc="upper left", frameon=False, fontsize=6, handlelength=1.1,
        labelspacing=0.35, borderaxespad=0.9)

    # Panel b — every score as a strip, by what the pair is. This is where
    # the controls earn their place: the floor is drawn, not asserted. The
    # classification comes from `s11_tables.pair_class` so the figure and
    # the report cannot disagree about what a pair is; the line break is a
    # display concern and is applied here only.
    classes = tb.pair_scores(_load("tm_scores.tsv"),
                             _load("structure_manifest.tsv"))
    order = [k for k in tb.PAIR_CLASSES if classes.get(k)]
    labels = [k.replace(", ", ",\n").replace(" vs ", " vs\n") for k in order]
    for i, key in enumerate(order):
        vals = classes[key]
        colour = (CONTROL_COLOUR if "control" in key else
                  fs.GROUP["RYR"] if key.startswith("RyR vs") else
                  fs.PARALOG["ITPR1"])
        ax2.scatter([i] * len(vals), vals, s=14, color=colour, alpha=0.55,
                    edgecolor="none", zorder=3)
        med = sorted(vals)[len(vals) // 2]
        ax2.plot([i - 0.28, i + 0.28], [med, med], color=fs.INK, lw=1.4,
                 zorder=4)
    for bar, label in ((TM_RANDOM, "random"), (TM_FOLD, "same fold")):
        ax2.axhline(bar, color=fs.MUTED, lw=0.7, ls="--", zorder=1)
        # Left-aligned at the axis edge: the right-hand end is where the
        # control strip sits, and the label landed on top of its points.
        ax2.text(-0.45, bar + 0.012, label, fontsize=6, color=fs.MUTED,
                 ha="left", va="bottom")
    ax2.set_xlim(-0.55, len(order) - 0.45)
    ax2.set_xticks(range(len(order)))
    ax2.set_xticklabels([l.replace("IP3R", "IP$_3$R") for l in labels],
                        fontsize=6)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("TM-score (the smaller normalisation)")
    fs.panel(ax2, "b", "what a TM-score means on this panel")
    fs.despine(ax2)
    fs.hgrid(ax2)
    fs.save(fig, FIGS / "s11_tm_calibration")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3 — per-domain confidence
# --------------------------------------------------------------------------

def fig_plddt() -> None:
    rows = [r for r in _load("plddt_domains.tsv") if r["placed"] == "1"]
    if not rows:
        return
    summary = _load("plddt_summary.tsv")
    labels = [s["label"] for s in summary]
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.6))
    by_label: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_label[r["label"]].append(r)
    for i, label in enumerate(labels):
        vals = by_label.get(label, [])
        for r in vals:
            ax.scatter(i + 0.0, _f(r["mean_plddt"]), s=20,
                       color=_role_colour(r), alpha=0.7, edgecolor="none",
                       zorder=3)
        if vals:
            ms = sorted(_f(r["mean_plddt"]) for r in vals)
            ax.plot([i - 0.3, i + 0.3], [ms[len(ms) // 2]] * 2,
                    color=fs.INK, lw=1.4, zorder=4)
    for bar, label in ((70, "confident"), (90, "very high")):
        ax.axhline(bar, color=fs.MUTED, lw=0.7, ls="--", zorder=1)
        ax.text(len(labels) - 0.4, bar, f"{label} ", fontsize=6,
                color=fs.MUTED, ha="right", va="bottom")
    ax.text(-0.45, 50.6, "axis starts at AlphaFold's 'very low' bound (50)",
            fontsize=5.5, color=fs.MUTED, va="bottom")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([_short(l) for l in labels], rotation=18, ha="right",
                       fontsize=6.5)
    ax.set_ylabel("mean pLDDT")
    # The axis starts at 50 because that is AlphaFold's own band boundary
    # (below 50 is "very low"), not at a number chosen to fill the panel —
    # and nothing on this panel is below it. A 0-100 axis puts every point
    # in the top third and makes the 70 and 90 bands unreadable, which are
    # the only things a pLDDT is read against.
    ax.set_ylim(50, 100)
    fs.panel(ax, "", "AlphaFold confidence, per domain, along the subunit")
    fs.despine(ax)
    fs.hgrid(ax)
    fs.save(fig, FIGS / "s11_plddt_domains")
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 4 — the panel itself
# --------------------------------------------------------------------------

def fig_panel() -> None:
    man = [r for r in _load("structure_manifest.tsv")]
    if not man:
        return
    man.sort(key=lambda r: (r["role"], -_f(r["resolved_residues"])))
    fig, ax = plt.subplots(figsize=(fs.W_FULL, max(2.6, 0.14 * len(man) + 1)))
    ys = range(len(man))
    for y, r in zip(ys, man):
        n = _f(r["resolved_residues"])
        expected = _f(r["expected_length"])
        if expected:
            ax.barh(y, expected, color="#e8e6df", height=0.7, zorder=1)
        ax.barh(y, n, color=_role_colour(r), height=0.7, zorder=2,
                alpha=0.9 if r["status"] == "ok" else 0.35)
    ax.set_yticks(list(ys))
    ax.set_yticklabels(
        [f"{r['id'][:30]}  [{r['role']}]" for r in man], fontsize=5.4)
    ax.invert_yaxis()
    ax.set_xlabel("residues — full bar: sequence length, filled: resolved / modelled")
    fs.panel(ax, "", "the S11 structure panel")
    fs.despine(ax)
    fs.save(fig, FIGS / "s11_panel")
    plt.close(fig)


FIGURES = {
    "afdb_coverage": fig_afdb_coverage,
    "tm_calibration": fig_tm_calibration,
    "plddt_domains": fig_plddt,
    "panel": fig_panel,
}


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args(argv)
    if args.list:
        for k in sorted(FIGURES):
            print(k)
        return 0
    fs.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for name in (args.only or sorted(FIGURES)):
        FIGURES[name]()
        print(f"[s11] figure {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
