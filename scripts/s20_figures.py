"""The S20 figures — rendered only from committed S20 tables (D13, D19).

Four, one per thing the task establishes:

  range_by_phylum      where in the eukaryotic tree the family is, as a
                       fraction of *swept proteomes* rather than of records,
                       so a well-sequenced phylum cannot outvote a sparse one
  profile_separation   D14 outside the vertebrates: both profiles' scores on
                       every S20 target, with the no-call band drawn
  plant_fungal_chase   the verdicts, and the cross-kingdom identities the
                       contamination test rests on, with its threshold drawn
  jackhmmer_s20        per-group convergence beside the sister-family trace
                       K1 is evaluated on

Colour follows the project's rule: categorical colour is spent on family
identity (ITPR blue, RyR the accent violet that marks the sister family
everywhere in this project) and nothing else. Taxonomy is carried by row
order and labels, magnitude by the single blue ramp.

Run:  python3 scripts/s20_figures.py            # all four
      python3 scripts/s20_figures.py --only range_by_phylum
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as fs  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    CENSUS_V5_DIR, FIG_DIR, S20_DIR, log, read_json, read_tsv,
)
from scripts.s20_verdicts import CONTAMINANT_PIDENT, OUTLIER_PIDENT  # noqa: E402
from scripts.s3_assign import MIN_SCORE, REL_MARGIN  # noqa: E402

ITPR_C = fs.PARALOG["ITPR1"]
RYR_C = fs.GROUP["RYR"]
NONE_C = fs.FAINT

VERDICT_C = {
    "real_gene": fs.STATUS["found_annotated"],
    "fragment": fs.STATUS["fragment"],
    "no_genome_backing": fs.STATUS["assembly_gap"],
    "cross_kingdom_outlier": fs.STATUS["tblastn_trace_ambiguous"],
    "module_only": fs.STATUS["tblastn_trace"],
    "contaminant_suspect": fs.STATUS["absent"],
    "unresolved": fs.GRID,
}
VERDICT_ORDER = ["real_gene", "fragment", "no_genome_backing",
                 "cross_kingdom_outlier", "module_only",
                 "contaminant_suspect", "unresolved"]


def _rows(name: str, base: Path = S20_DIR) -> list[dict]:
    path = base / name
    return read_tsv(path) if path.exists() else []


# ------------------------------------------------------ 1. range by phylum
def range_by_phylum(stem: Path) -> None:
    """Presence as a fraction of *swept proteomes*, per clade.

    Per proteome rather than per record, so a well-sequenced phylum cannot
    outvote a sparse one — 300 arthropod proteomes and 3 placozoan ones each
    get a rate, not a count.

    The prokaryote groups collapse to one row apiece. They contribute 38
    bacterial and 5 archaeal phyla, every one of them zero, and spending
    forty rows on a single negative claim would bury the eukaryote result
    the figure exists to show. The row label carries the phylum count so
    nothing is hidden by the collapse.
    """
    presence = _rows("proteome_presence.tsv")
    taxa = {int(r["taxid"]): r for r in _rows("taxonomy.tsv")}
    if not presence:
        raise SystemExit("no proteome_presence.tsv — run s20_run_sweep.py")

    COLLAPSE = ("Bacteria", "Archaea")
    agg: dict[tuple[str, str], list[int]] = {}
    clades: dict[str, set] = {}
    for r in presence:
        tx = taxa.get(int(r["taxid"]), {})
        group = tx.get("group", "unclassified") or "unclassified"
        clade = tx.get("clade") or tx.get("phylum") or "unclassified"
        clades.setdefault(group, set()).add(clade)
        key = (group, "" if group in COLLAPSE else clade)
        cell = agg.setdefault(key, [0, 0])
        cell[0] += 1
        cell[1] += 1 if r["itpr_status"] == "present" else 0

    keep = {k: v for k, v in agg.items()
            if v[0] >= 3 or k[0] in COLLAPSE}
    rate = {}
    for (g, _), v in keep.items():
        tot = rate.setdefault(g, [0, 0])
        tot[0] += v[0]
        tot[1] += v[1]
    order = sorted(keep.items(),
                   key=lambda kv: (-rate[kv[0][0]][1] / rate[kv[0][0]][0],
                                   kv[0][0], -kv[1][1] / kv[1][0], -kv[1][0]))

    fig, ax = plt.subplots(figsize=(fs.W_FULL,
                                    min(fs.H_MAX, 1.0 + 0.16 * len(order))))
    ys = list(range(len(order)))
    fracs = [v[1] / v[0] for _, v in order]
    ax.barh(ys, fracs, height=0.7,
            color=[fs.BLUES[-1] if f > 0 else fs.STATUS["absent"]
                   for f in fracs],
            edgecolor="none")
    labels = []
    for (g, clade), _ in order:
        labels.append(clade if clade else
                      f"{g} — all {len(clades[g])} phyla")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100 %"], fontsize=fs.FS_TICK)
    ax.set_xlabel("reference proteomes carrying an ITPR call",
                  fontsize=fs.FS_LABEL)
    for y, (_, v) in zip(ys, order):
        ax.text(v[1] / v[0] + 0.012, y, f"{v[1]}/{v[0]}", va="center",
                fontsize=fs.FS_NOTE, color=fs.MUTED)
    last = None
    for y, ((g, _), _) in zip(ys, order):
        if g != last:
            ax.text(1.0, y - 0.62, g, ha="right", va="bottom",
                    fontsize=fs.FS_NOTE, color=fs.ACCENT, style="italic")
            last = g
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    fs.panel(ax, "", "Where the family is, per swept proteome")
    fig.tight_layout()
    fs.save(fig, stem)
    plt.close(fig)


# ------------------------------------------------- 2. profile separation
def profile_separation(stem: Path) -> None:
    pts = []
    for path in sorted(S20_DIR.glob("assignments_*.tsv")):
        for r in read_tsv(path):
            try:
                pts.append((float(r["itpr_score"] or 0),
                            float(r["ryr_score"] or 0), r["assignment"]))
            except ValueError:
                continue
    if not pts:
        raise SystemExit("no assignments_*.tsv — run s20_run_sweep.py")

    fig, ax = plt.subplots(figsize=(fs.W_HALF * 1.55, fs.W_HALF * 1.35))
    for call, colour, z in (("unassigned", NONE_C, 1), ("RYR", RYR_C, 2),
                            ("ITPR", ITPR_C, 3)):
        xs = [p[0] for p in pts if p[2] == call]
        ys = [p[1] for p in pts if p[2] == call]
        ax.scatter(xs, ys, s=4, c=colour, alpha=0.55, linewidths=0,
                   zorder=z, label=f"{call} ({len(xs)})")
    lim = max([p[0] for p in pts] + [p[1] for p in pts]) * 1.15
    # The no-call band is where the assignment rule declines to speak:
    # |win - lose| / win <= REL_MARGIN, i.e. x(1-m) <= y <= x/(1-m).
    edge = [MIN_SCORE * 0.5, lim]
    ax.fill_between(edge, [e * (1 - REL_MARGIN) for e in edge],
                    [e / (1 - REL_MARGIN) for e in edge],
                    color=fs.HILITE, zorder=0, linewidth=0)
    ax.plot(edge, edge, color=fs.GRID, lw=0.6, zorder=0)
    ax.axvline(MIN_SCORE, color=fs.GRID, lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.axhline(MIN_SCORE, color=fs.GRID, lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.set_xscale("symlog", linthresh=MIN_SCORE)
    ax.set_yscale("symlog", linthresh=MIN_SCORE)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("itpr.hmm bit score", fontsize=fs.FS_LABEL)
    ax.set_ylabel("ryr.hmm bit score", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right",
              markerscale=2.0)
    ax.text(0.03, 0.97,
            f"shaded: the {REL_MARGIN:.0%} no-call band\n"
            f"dashed: the {MIN_SCORE:.0f}-bit floor\n"
            f"0 on an axis = that profile did not score it",
            transform=ax.transAxes, fontsize=fs.FS_NOTE, va="top",
            color=fs.MUTED, linespacing=1.4)
    fs.despine(ax)
    fs.panel(ax, "", "D14 outside the vertebrates")
    fig.tight_layout()
    fs.save(fig, stem)
    plt.close(fig)


# ------------------------------------------------- 3. the plant/fungal chase
def plant_fungal_chase(stem: Path) -> None:
    rows = _rows("plant_fungal_verdicts.tsv")
    if not rows:
        raise SystemExit("no plant_fungal_verdicts.tsv — run s20_verdicts.py")

    classes: list[str] = []
    for r in rows:
        key = f"{r['kingdom']} / {r['lineage_class'].replace('_', ' ')}"
        if key not in classes:
            classes.append(key)
    classes.sort()

    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.6),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax = axes[0]
    left = [0.0] * len(classes)
    for v in VERDICT_ORDER:
        widths = [sum(1 for r in rows
                      if f"{r['kingdom']} / {r['lineage_class'].replace('_', ' ')}" == c
                      and r["verdict"] == v) for c in classes]
        if not any(widths):
            continue
        ax.barh(range(len(classes)), widths, left=left, height=0.68,
                color=VERDICT_C.get(v, fs.GRID), edgecolor="none",
                label=v.replace("_", " "))
        left = [a + b for a, b in zip(left, widths)]
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("records chased", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, ncol=2,
              loc="lower right", handlelength=1.0)
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    fs.panel(ax, "a", "Verdict per record")

    ax = axes[1]
    vals = [(float(r["out_kingdom_pident"]), r["verdict"]) for r in rows
            if r["out_kingdom_pident"] not in ("", None)]
    for i, (p, v) in enumerate(sorted(vals)):
        ax.scatter(p, i, s=9, c=VERDICT_C.get(v, fs.GRID), linewidths=0)
    ax.axvline(CONTAMINANT_PIDENT, color=fs.STATUS["absent"], lw=0.8)
    ax.axvline(OUTLIER_PIDENT, color=fs.FAINT, lw=0.8, ls=(0, (2, 2)))
    ax.text(CONTAMINANT_PIDENT - 1.5, len(vals) * 0.5,
            f"contamination call\n{CONTAMINANT_PIDENT:.0f} %", ha="right",
            fontsize=fs.FS_NOTE, color=fs.STATUS["absent"])
    ax.set_xlabel("identity to the nearest outside-kingdom protein (%)",
                  fontsize=fs.FS_LABEL)
    ax.set_ylabel("records, sorted", fontsize=fs.FS_LABEL)
    ax.set_yticks([])
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "b", "The contamination test")
    fig.tight_layout()
    fs.save(fig, stem)
    plt.close(fig)


# ----------------------------------------------------- 4. jackhmmer per group
def jackhmmer_s20(stem: Path) -> None:
    """Per-group convergence beside the drift K1 watches — and the drift it
    does not.

    Panel b plots the sister-family share (what K1 evaluates) *and* the
    off-family share on the same axes, because that contrast is the finding.
    S3 recorded it as D10b: off-family accretion dilutes the sister share
    instead of raising it, so a model can grow thirty-fold into proteins of
    neither family with K1 reading a flat zero throughout. A panel showing
    only the sister trace would be a flat line at 0 and would look like a
    clean run.
    """
    rows = _rows("jackhmmer_convergence_s20.tsv")
    if not rows:
        raise SystemExit("no jackhmmer_convergence_s20.tsv — "
                         "run s20_jackhmmer.py")
    groups = sorted({r["group"] for r in rows})
    verdicts = ((read_json(S20_DIR / "jackhmmer_verdicts_s20.json") or {})
                .get("verdicts", {}))
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.6))
    shades = [fs.BLUES[-1], fs.BLUES[-3], fs.BLUES[-5], fs.MUTED]
    for i, g in enumerate(groups):
        sub = sorted((r for r in rows if r["group"] == g),
                     key=lambda r: int(r["round"]))
        xs = [int(r["round"]) for r in sub]
        colour = shades[i % len(shades)]
        n_inc = [int(r["n_included"] or 0) for r in sub]
        axes[0].plot(xs, n_inc, marker="o", ms=2.6, lw=1.0,
                     color=colour, label=g)
        # Mark the round the rule actually fired on. `kill_rule` in the
        # convergence table is a run-level field repeated on every row, so
        # reading it per row would mark round 1 for a run killed at round 3;
        # `killed_at` is the round, and it lives in the verdicts JSON.
        v = verdicts.get(g, {})
        at = v.get("killed_at")
        if at:
            here = [r for r in sub if int(r["round"]) == int(at)]
            if here:
                y = int(here[0]["n_included"] or 0)
                axes[0].scatter([at], [y], s=44, facecolors="none",
                                edgecolors=fs.STATUS["absent"], linewidths=1.0,
                                zorder=5)
                axes[0].annotate(v.get("rule", ""), (at, y),
                                 textcoords="offset points", xytext=(5, -11),
                                 fontsize=fs.FS_NOTE,
                                 color=fs.STATUS["absent"])
        axes[1].plot(xs, [float(r["sister_frac"] or 0) for r in sub],
                     marker="o", ms=2.6, lw=1.1, color=colour,
                     label=f"{g} — sister family")
        off = [int(r["n_offfamily"] or 0) / max(int(r["n_included"] or 1), 1)
               for r in sub]
        axes[1].plot(xs, off, marker="^", ms=2.6, lw=1.0, ls=(0, (3, 2)),
                     color=colour, label=f"{g} — neither family")

    axes[0].set_yscale("log")
    axes[0].set_xlabel("jackhmmer round", fontsize=fs.FS_LABEL)
    axes[0].set_ylabel("targets in the model", fontsize=fs.FS_LABEL)
    axes[0].legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(axes[0]); fs.hgrid(axes[0])
    fs.panel(axes[0], "a", "Convergence, and where D10 killed it")

    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_xlabel("jackhmmer round", fontsize=fs.FS_LABEL)
    axes[1].set_ylabel("share of the model", fontsize=fs.FS_LABEL)
    axes[1].legend(fontsize=fs.FS_NOTE, frameon=False, loc="center right")
    axes[1].text(0.03, 0.96, "K1 watches the solid line;\nthe dashed one is "
                 "what it cannot see", transform=axes[1].transAxes,
                 fontsize=fs.FS_NOTE, va="top", color=fs.MUTED,
                 linespacing=1.4)
    fs.despine(axes[1]); fs.hgrid(axes[1])
    fs.panel(axes[1], "b", "Drift, sister and off-family")
    fig.tight_layout()
    fs.save(fig, stem)
    plt.close(fig)


FIGURES = {
    "range_by_phylum": range_by_phylum,
    "profile_separation": profile_separation,
    "plant_fungal_chase": plant_fungal_chase,
    "jackhmmer_s20": jackhmmer_s20,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(" ", k)
        return 0
    fs.use()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rc = 0
    for name in (args.only or list(FIGURES)):
        try:
            FIGURES[name](FIG_DIR / name)
            log("s20_figures", f"{name} ✓")
        except SystemExit as exc:
            log("s20_figures", f"{name} skipped: {exc}")
            rc = 1
    _ = CENSUS_V5_DIR
    return rc


if __name__ == "__main__":
    sys.exit(main())
