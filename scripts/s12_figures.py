"""The four S12 figures, from committed S12 tables only (D13, D19).

Three decisions worth stating.

**The junction figure is per junction, not per gene.**  A bar saying "this
gene has junction-spanning reads" is the claim S10 could already almost
make.  The claim S12 exists to make is that the *particular* junctions no
annotated model spans are spliced, so the figure plots them individually
along the transcript, with the annotated ones drawn beside them in the same
panel as the internal control.

**The decoy floor is drawn, not stated.**  A detection rule whose third
clause is "more reads than its own decoy" has to show where the decoy sits;
a caption saying "decoys were zero" asks to be believed.

**Coverage is drawn on a log count axis with a zero-safe floor.**  These
libraries differ ~40x in depth and the interesting comparison is between a
locus and its own decoy within a run, not between runs; on a linear axis
every low-depth run collapses onto the baseline and the panel becomes a
picture of sequencing depth.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import figstyle as F  # noqa: E402
from s12_lib import OUT_DIR  # noqa: E402

FIGS = OUT_DIR / "figures"
CELL_ORDER = ["ITPR1", "ITPR2", "ITPR3", "RYR"]
ROLE_MARK = {"failed": "o", "control_paralog": "s", "control_family": "^",
             "housekeeping": "D"}


def rows(name: str) -> list[dict]:
    p = OUT_DIR / name
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _f(r, k, d=0.0):
    try:
        return float(r.get(k) or d)
    except ValueError:
        return d


def _i(r, k, d=0):
    try:
        return int(float(r.get(k) or d))
    except ValueError:
        return d


def _species_order(rs):
    return sorted({r["species"] for r in rs})


def _colour(cell):
    return F.GROUP.get(cell, F.PARALOG.get(cell, F.FAINT))


def _abbrev(organism: str, cell: str = "") -> str:
    """`Nibea albiflora` + ITPR2 -> `N. albiflora ITPR2`, safely.

    Defensive about a missing or one-word name: these labels are built from
    whichever table the panel reads, and a panel rendered while an upstream
    stage is still running would otherwise die on an empty string rather
    than draw what it has.
    """
    parts = (organism or "").split()
    name = (f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2
            else (parts[0] if parts else "?"))
    return f"{name} {cell}".strip()


# ---------------------------------------------------------------------------
# Fig 1 — the detection result, locus by locus, against its own decoy
# ---------------------------------------------------------------------------

def fig_detection(log=print) -> None:
    loci = [r for r in rows("expression_by_locus.tsv")
            if r["role"] != "housekeeping"]
    if not loci:
        log("  fig1: no expression_by_locus.tsv")
        return
    order = _species_order(loci)
    fig, axes = plt.subplots(1, len(order), figsize=(F.W_FULL, 2.6),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, sp in zip(axes, order):
        sub = sorted([r for r in loci if r["species"] == sp],
                     key=lambda r: CELL_ORDER.index(r["cell"])
                     if r["cell"] in CELL_ORDER else 9)
        x = np.arange(len(sub))
        reads = [max(_i(r, "reads"), 0) for r in sub]
        decoy = [max(_i(r, "decoy_reads"), 0) for r in sub]
        ax.bar(x - 0.19, [max(v, 0.4) for v in reads], width=0.36,
               color=[_colour(r["cell"]) for r in sub], zorder=3)
        ax.bar(x + 0.19, [max(v, 0.4) for v in decoy], width=0.36,
               color=F.GRID, edgecolor=F.FAINT, linewidth=0.5, zorder=3,
               label="reversed decoy")
        # The "annotation loses this gene" mark is a *marker*, not a text
        # glyph: figstyle.save() refuses a figure needing a glyph the
        # document font lacks, and it refused this panel's first version
        # over a multiplication sign.
        fail_x = [xi for xi, r in zip(x, sub) if r["role"] == "failed"]
        if fail_x:
            ax.plot(fail_x, [0.36] * len(fail_x), marker="x", lw=0,
                    markersize=3.6, markeredgewidth=1.0, color=F.ACCENT,
                    zorder=5, clip_on=False)
        ax.set_yscale("log")
        ax.set_ylim(0.3, max(4, max(reads + [1]) * 3))
        ax.set_xticks(x)
        ax.set_xticklabels([r["cell"] for r in sub], fontsize=F.FS_TICK)
        ax.set_title(sub[0]["organism"], fontsize=F.FS_TITLE,
                     style="italic", color=F.INK)
        F.despine(ax)
        F.hgrid(ax)
    axes[0].set_ylabel("reads (pooled over runs)", fontsize=F.FS_LABEL)
    # A log axis cannot draw a zero, so every bar is floored at 0.4 — well
    # below the y-limit's first tick.  The legend has to say so, or a decoy
    # that collected nothing reads as a decoy that collected a little.
    handles = [plt.Rectangle((0, 0), 1, 1, color=F.FAINT,
                             label="recovered locus"),
               plt.Rectangle((0, 0), 1, 1, facecolor=F.GRID,
                             edgecolor=F.FAINT, linewidth=0.5,
                             label="its reversed decoy (bars at the axis "
                                   "floor are zero)"),
               plt.Line2D([], [], color=F.ACCENT, lw=0, marker="x",
                          markersize=3.6, markeredgewidth=1.0,
                          label="the annotation loses this gene")]
    axes[0].legend(handles=handles, fontsize=F.FS_NOTE, frameon=False,
                   ncol=3, loc="lower left", bbox_to_anchor=(0, 1.02))
    fig.suptitle("Reads on each recovered locus, against its own "
                 "composition-matched decoy", fontsize=F.FS_SUPTITLE,
                 y=1.12)
    fig.tight_layout()
    F.save(fig, FIGS / "s12_detection")
    plt.close(fig)
    log("  fig1: s12_detection")


# ---------------------------------------------------------------------------
# Fig 2 — junction-by-junction support along each failed transcript
# ---------------------------------------------------------------------------

def fig_junctions(log=print) -> None:
    js = rows("junction_support.tsv")
    refs = {(r["species"], r["seq"]): r
            for r in rows("expression_by_locus.tsv")}
    failed = [k for k, r in refs.items() if r["role"] == "failed"]
    if not js or not failed:
        log("  fig2: no junction_support.tsv")
        return
    failed.sort(key=lambda k: (k[0], k[1]))
    n = len(failed)
    fig, axes = plt.subplots(n, 1, figsize=(F.W_FULL, 0.62 * n + 0.9),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, key in zip(axes, failed):
        sub = sorted([r for r in js
                      if (r["species"], r["seq"]) == key],
                     key=lambda r: _i(r, "cds_offset"))
        if not sub:
            continue
        xs = [_i(r, "cds_offset") for r in sub]
        ys = [_i(r, "reads") for r in sub]
        ann = [_i(r, "annotated") for r in sub]
        ax.vlines(xs, 0, [max(y, 0) for y in ys], linewidth=0.8,
                  color=[F.BLUES[4] if a else F.ACCENT for a in ann],
                  zorder=3)
        ax.scatter([x for x, y in zip(xs, ys) if y == 0],
                   [0] * sum(1 for y in ys if y == 0), s=5, marker="x",
                   color=F.FAINT, zorder=4, linewidths=0.6)
        ref = refs[key]
        ax.set_ylabel(_abbrev(ref.get("organism", ""), ref["cell"]),
                      fontsize=F.FS_TICK, rotation=0, ha="right", va="center")
        ax.set_yscale("symlog", linthresh=1)
        ax.tick_params(labelsize=F.FS_TICK)
        F.despine(ax)
        F.hgrid(ax)
    axes[-1].set_xlabel("position in the spliced coding sequence (nt)",
                        fontsize=F.FS_LABEL)
    handles = [plt.Line2D([], [], color=F.ACCENT, lw=1.4,
                          label="junction the annotation does not model"),
               plt.Line2D([], [], color=F.BLUES[4], lw=1.4,
                          label="junction the annotation models"),
               plt.Line2D([], [], color=F.FAINT, lw=0, marker="x",
                          markersize=3, label="no read crossed it")]
    axes[0].legend(handles=handles, fontsize=F.FS_NOTE, frameon=False,
                   ncol=3, loc="lower left", bbox_to_anchor=(0, 1.05))
    fig.tight_layout()
    F.save(fig, FIGS / "s12_junctions")
    plt.close(fig)
    log("  fig2: s12_junctions")


# ---------------------------------------------------------------------------
# Fig 3 — where the reads fall relative to what the annotation delivers
# ---------------------------------------------------------------------------

def fig_gap_coverage(log=print) -> None:
    gap = [r for r in rows("annotation_gap_coverage.tsv")
           if r["role"] in ("failed", "control_paralog", "control_family")]
    if not gap:
        log("  fig3: no annotation_gap_coverage.tsv")
        return
    gap.sort(key=lambda r: (-_f(r, "frac_coverage_unannotated"), r["species"]))
    y = np.arange(len(gap))
    fig, ax = plt.subplots(figsize=(F.W_FULL, 0.26 * len(gap) + 1.0))
    ax.barh(y, [_f(r, "frac_coverage_unannotated") for r in gap],
            color=[_colour(r["cell"]) for r in gap], zorder=3, height=0.66)
    for yi, r in zip(y, gap):
        loss = r.get("annotation_loss", "")
        if loss not in ("", None):
            ax.plot(_f(r, "annotation_loss"), yi, marker="|", markersize=7,
                    color=F.INK, zorder=5, markeredgewidth=1.1)
    ax.set_yticks(y)
    ax.set_yticklabels([_abbrev(r.get("organism", ""), r["cell"])
                        for r in gap], fontsize=F.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("fraction of read coverage falling outside any annotated "
                  "coding block", fontsize=F.FS_LABEL)
    ax.tick_params(labelsize=F.FS_TICK)
    F.despine(ax, keep=("bottom",))
    ax.grid(axis="x", color=F.GRID, linewidth=0.5, zorder=0)
    ax.plot([], [], marker="|", color=F.INK, lw=0, markersize=7,
            markeredgewidth=1.1, label="S10 annotation loss (coding bp)")
    ax.legend(fontsize=F.FS_NOTE, frameon=False, loc="lower right")
    # Descriptive, not concluding. Whether the two measures agree is a
    # result the report computes and renders a verdict on; a figure title
    # that states the conclusion would say it whatever the bars did.
    fig.suptitle("Read coverage outside the annotation, against the coding "
                 "footprint S10 measured", fontsize=F.FS_SUPTITLE, y=0.99)
    fig.tight_layout()
    F.save(fig, FIGS / "s12_gap_coverage")
    plt.close(fig)
    log("  fig3: s12_gap_coverage")


# ---------------------------------------------------------------------------
# Fig 4 — the two instruments side by side: reads vs deposits
# ---------------------------------------------------------------------------

def fig_instruments(log=print) -> None:
    summ = rows("atlas_summary.tsv")
    loci = {(r["species"], r["cell"]): r
            for r in rows("expression_by_locus.tsv")}
    if not summ:
        log("  fig4: no atlas_summary.tsv")
        return
    organism = {r["species"]: r["organism"] for r in summ}
    keys = sorted({(r["species"], r["cell"]) for r in summ})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(F.W_FULL, 2.5))
    y = np.arange(len(keys))
    labels = []
    reads_frac, dep_frac = [], []
    for k in keys:
        loc = loci.get(k, {})
        tot = _i(loc, "junctions_unannotated")
        hit = _i(loc, "junctions_unannotated_covered")
        reads_frac.append(hit / tot if tot else 0.0)
        s = [r for r in summ if (r["species"], r["cell"]) == k
             and _i(r, "annotated") == 0]
        pt = sum(_i(r, "probes") for r in s)
        ps = sum(_i(r, "probes_spanned") for r in s)
        dep_frac.append(ps / pt if pt else 0.0)
        labels.append(_abbrev(organism.get(k[0], k[0]), k[1]))
    for ax, vals, title in ((ax1, reads_frac,
                             "streamed RNA-seq reads"),
                            (ax2, dep_frac,
                             "submitted transcript records")):
        ax.barh(y, vals, color=[_colour(k[1]) for k in keys], height=0.66,
                zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(labels if ax is ax1 else [], fontsize=F.FS_TICK)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.02)
        ax.set_xlabel("unannotated junctions with evidence",
                      fontsize=F.FS_LABEL)
        ax.set_title(title, fontsize=F.FS_TITLE)
        ax.tick_params(labelsize=F.FS_TICK)
        F.despine(ax, keep=("bottom",))
        ax.grid(axis="x", color=F.GRID, linewidth=0.5, zorder=0)
    fig.tight_layout()
    F.save(fig, FIGS / "s12_instruments")
    plt.close(fig)
    log("  fig4: s12_instruments")


FIGURES = {"detection": fig_detection, "junctions": fig_junctions,
           "gap": fig_gap_coverage, "instruments": fig_instruments}


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    args = ap.parse_args()
    F.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for name, fn in FIGURES.items():
        if args.only and name not in args.only:
            continue
        fn()


if __name__ == "__main__":
    main()
