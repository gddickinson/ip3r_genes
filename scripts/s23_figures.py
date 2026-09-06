"""s23_figures.py — the S23 figures, from the committed tables only.

D13/D19 applied to figures: every panel is drawn from `results/s23_scope/`,
never from a sweep summary or a recomputed number, so a figure cannot disagree
with the report beside it.

Four figures, and each one draws the thing the text would otherwise assert:

  copy_number        how many ITPR genes a genome has, outside the
                     vertebrates. The vertebrate trio is drawn as a reference
                     line rather than as a category, because ITPR1/2/3 are a
                     2R product and the question here is copy number.
  absence_at_genome  the negative claims, with the **controlled** genome count
                     drawn beside the swept one. An absence bar whose control
                     bar is short is a claim standing on nothing, and that has
                     to be visible rather than stated in a footnote.
  identity_floor     the threshold S23b measured: annotation-confirmed loci
                     against annotation-contradicted ones, with the floor
                     drawn where the measurement put it. If the two
                     populations touch, the figure shows them touching.
  span_inflation     locus span against CDS footprint, per group — the
                     measurement that says why copy number is counted on
                     alignments and not on clusters.

Usage:
  python3 scripts/s23_figures.py [--only copy_number]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402

import figstyle as fs                                        # noqa: E402

OUT = PROJECT_ROOT / "results" / "s23_scope"
FIGDIR = OUT / "figures"

#: Kingdom-level groups in reading order, with the palette they borrow. The
#: colours come from `figstyle.GROUP`, which already encodes plant / protist /
#: fungi / invertebrate, so an S23 figure and an S6 alignment figure name the
#: same lineage in the same colour.
GROUPS = [("metazoa", "invert_metazoa"), ("viridiplantae", "plant"),
          ("fungi", "fungi"), ("sar", "protist"), ("amoebozoa", "protist"),
          ("discoba", "protist"), ("other", "protist")]
GROUP_COLOUR = {g: fs.GROUP[k] for g, k in GROUPS}
GROUP_LABEL = {"metazoa": "Metazoa (non-vert.)", "viridiplantae": "Viridiplantae",
               "fungi": "Fungi", "sar": "SAR", "amoebozoa": "Amoebozoa",
               "discoba": "Discoba", "other": "other eukaryotes"}


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def _int(v, default: int = 0) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _float(v, default: float = 0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ------------------------------------------------------------------ figures

def fig_copy_number(ledger: list[dict]) -> None:
    """Copy number per genome, by group — controlled genomes only."""
    order = [g for g, _ in GROUPS]
    rows = [r for r in ledger
            if (r.get("control_verdict") or "").startswith("controlled")]
    by_group: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        by_group[r.get("group") or "other"][_int(r.get("n_full"))] += 1
    groups = [g for g in order if by_group.get(g)]
    max_n = max((max(c) for c in by_group.values() if c), default=0)
    bins = list(range(0, min(max_n, 6) + 1))

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.9))
    width = 0.8 / max(1, len(groups))
    for i, g in enumerate(groups):
        c = by_group[g]
        vals = [sum(v for k, v in c.items() if k == b) if b < bins[-1]
                else sum(v for k, v in c.items() if k >= b) for b in bins]
        xs = [b + (i - (len(groups) - 1) / 2) * width for b in bins]
        ax.bar(xs, vals, width=width * 0.92, color=GROUP_COLOUR[g],
               label=f"{GROUP_LABEL[g]} (n={sum(c.values())})")
    ax.set_xticks(bins)
    ax.set_xticklabels([str(b) if b < bins[-1] else f"{b}+" for b in bins])
    ax.set_xlabel("complete ITPR gene models per genome")
    ax.set_ylabel(f"genomes (n = {len(rows)})")
    ax.set_title("Copy number outside the vertebrates", loc="left")
    ax.axvline(3, color="#52514e", lw=0.8, ls=":")
    ax.text(3.06, ax.get_ylim()[1] * 0.95, "vertebrate\nparalog count",
            fontsize=6.2, va="top", color="#52514e")
    fs.despine(ax)
    fs.hgrid(ax)
    ax.legend(fontsize=6.2, frameon=False, loc="upper right")
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "copy_number"))
    plt.close(fig)


def fig_absence(absences: list[dict]) -> None:
    """Every absence clade: swept, controlled, and what was found."""
    rows = [r for r in absences if _int(r.get("genomes_swept"))]
    rows.sort(key=lambda r: -_int(r.get("swept_proteomes")))
    rows = rows[:22]
    labels = [f"{r['clade']} ({r['swept_proteomes']})" for r in rows]
    y = list(range(len(rows)))[::-1]

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 0.24 * len(rows) + 1.3))
    swept = [_int(r["genomes_swept"]) for r in rows]
    ctl = [_int(r["genomes_controlled"]) for r in rows]
    found = [_int(r["genomes_with_full_itpr"]) for r in rows]
    ax.barh(y, swept, color="#e7e6e0", height=0.72, label="genomes searched")
    ax.barh(y, ctl, color="#a9a79e", height=0.72,
            label="control fired (claim admissible)")
    ax.barh(y, found, color=fs.STATUS["found_annotated"], height=0.5,
            label="ITPR found")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.2)
    ax.set_xlabel("genomes (clade's swept proteome count in brackets)")
    ax.set_title("The absence claims, taken to assembly level", loc="left")
    fs.despine(ax, keep=("left", "bottom"))
    fs.hgrid(ax, axis="x")
    ax.legend(fontsize=6.2, frameon=False, loc="lower right")
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "absence_at_genome"))
    plt.close(fig)


def fig_identity_floor(loci: list[dict], cal: dict | None) -> None:
    """The measured floor: confirmed loci against contradicted ones."""
    pops = {
        "confirmed": [_float(r["identity"]) for r in loci
                      if r["evidence"] == "confirmed"],
        "contradicted": [_float(r["identity"]) for r in loci
                         if r["evidence"] in ("contradicted", "sister")],
        "unnamed": [_float(r["identity"]) for r in loci
                    if r["evidence"] in ("unnamed", "no_annotation")],
    }
    colours = {"confirmed": fs.STATUS["found_annotated"],
               "contradicted": fs.STATUS["absent"],
               "unnamed": "#a9a79e"}
    labels = {"confirmed": "the assembly's own annotation names a family gene",
              "contradicted": "…names a different gene (or a RyR)",
              "unnamed": "no informative annotation"}

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.6))
    for i, (k, vals) in enumerate(pops.items()):
        if not vals:
            continue
        ax.scatter(vals, [i + 0.0] * len(vals), s=11, alpha=0.55,
                   color=colours[k], edgecolors="none",
                   label=f"{labels[k]} (n={len(vals)})")
    if cal:
        floor = float(cal.get("call_min_identity", 0))
        rec = float(cal.get("record_min_identity", 0))
        ax.axvline(floor, color="#52514e", lw=1.0)
        ax.text(floor + 0.004, len(pops) - 0.35,
                f"call floor {floor:.2f}", fontsize=6.4, color="#52514e")
        ax.axvline(rec, color="#a9a79e", lw=0.8, ls=":")
        ax.text(rec + 0.004, -0.42, f"recorded from {rec:.2f}", fontsize=6.0,
                color="#8a897f")
    ax.set_yticks(range(len(pops)))
    ax.set_yticklabels(["confirmed", "contradicted", "no evidence"],
                       fontsize=6.6)
    ax.set_ylim(-0.6, len(pops) - 0.4)
    ax.set_xlabel("locus identity (best alignment)")
    ax.set_title("Where the locus identity floor sits, and what it separates",
                 loc="left")
    fs.despine(ax)
    ax.legend(fontsize=6.0, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, -0.32), ncol=1)
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "identity_floor"))
    plt.close(fig)


def fig_span_inflation(spans: list[dict]) -> None:
    """Locus span against CDS footprint — why copies are counted on alignments."""
    by_group: dict[str, list] = defaultdict(list)
    for r in spans:
        v = _float(r.get("span_inflation"))
        if v > 0:
            by_group[r.get("group") or "other"].append(v)
    groups = [g for g, _ in GROUPS if by_group.get(g)]

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.6))
    for i, g in enumerate(groups):
        vals = by_group[g]
        ax.scatter(vals, [i] * len(vals), s=13, alpha=0.5,
                   color=GROUP_COLOUR[g], edgecolors="none")
    ax.axvline(1.0, color="#52514e", lw=0.8, ls=":")
    ax.text(1.05, len(groups) - 0.4, "span = CDS footprint", fontsize=6.2,
            color="#52514e")
    ax.set_xscale("log")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([GROUP_LABEL[g] for g in groups], fontsize=6.6)
    ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.set_xlabel("locus span ÷ aligned CDS footprint (log scale)")
    ax.set_title("A locus is much bigger than the gene inside it", loc="left")
    fs.despine(ax)
    fig.tight_layout()
    fs.save(fig, str(FIGDIR / "span_inflation"))
    plt.close(fig)


FIGURES = {
    "copy_number": ("copy_number_ledger.tsv", fig_copy_number),
    "absence_at_genome": ("absence_at_genome.tsv", fig_absence),
    "identity_floor": ("locus_identity.tsv", fig_identity_floor),
    "span_inflation": ("locus_span.tsv", fig_span_inflation),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k, (src, _) in FIGURES.items():
            print(f"  {k:20s} <- {src}")
        return 0
    fs.use()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    cal_path = OUT / "locus_calibration.json"
    cal = json.loads(cal_path.read_text()) if cal_path.exists() else None
    want = args.only or list(FIGURES)
    made = 0
    for name in want:
        src, fn = FIGURES[name]
        path = OUT / src
        if not path.exists():
            print(f"  ! {name}: {src} not built yet")
            continue
        rows = read_tsv(path)
        fn(rows, cal) if name == "identity_floor" else fn(rows)
        print(f"  {name} <- {src} ({len(rows)} rows)")
        made += 1
    print(f"{made} figure(s) -> {FIGDIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
