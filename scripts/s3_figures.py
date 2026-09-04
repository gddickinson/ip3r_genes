"""The S3 figures — rendered only from the committed S3 tables (D13, D19).

Four figures, each answering one question the session had to settle:

  profile_separation   Do two profiles actually separate ITPR from RyR, and
                       does the no-call band sit where a band should sit?
  jackhmmer_convergence  Does iterating find anything the profile missed,
                       and did any run diverge under D10?
  instrument_agreement Which instrument calls each record, and what does the
                       second one add over the architecture rule alone?
  proteome_copy_number How many ITPRs does a vertebrate reference proteome
                       carry, and which carry none?

No figure module opens a database, reads a structure or recomputes a call.

Run:  python3 scripts/s3_figures.py [--only <slug>]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import figstyle  # noqa: E402
from scripts.s3_assign import REL_MARGIN  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    CENSUS_V3_DIR, HMM_SWEEP_DIR, read_tsv,
)

FIG_DIR = CENSUS_V3_DIR / "figures"
CALL_COLOUR = {"ITPR": figstyle.PARALOG["ITPR1"], "RYR": figstyle.GROUP["RYR"],
               "unassigned": figstyle.FAINT, "conflict": "#b3261e"}


def _fnum(row: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return default


# ------------------------------------------------------- 1. profile separation
def fig_profile_separation():
    import matplotlib.pyplot as plt

    rows = read_tsv(HMM_SWEEP_DIR / "calibration_assignments.tsv")
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(figstyle.W_FULL, 2.9),
        gridspec_kw={"width_ratios": [1.25, 1.0]})

    for call in ("unassigned", "RYR", "ITPR"):
        pts = [(_fnum(r, "itpr_score"), _fnum(r, "ryr_score"))
               for r in rows if r["s2_arch_call"] == call]
        if not pts:
            continue
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=3.0,
                   linewidths=0, alpha=0.55, color=CALL_COLOUR[call],
                   label=f"{call} ({len(pts):,})", rasterized=True)
    lim = max(max(_fnum(r, "itpr_score"), _fnum(r, "ryr_score"))
              for r in rows) * 1.05
    # The no-call band: |win - lose| / win < REL_MARGIN, i.e. the wedge
    # around the diagonal. Drawn, not described.
    ax.fill_between([0, lim], [0, lim * (1 - REL_MARGIN)],
                    [0, lim / (1 - REL_MARGIN)],
                    color=figstyle.HILITE, zorder=0)
    ax.plot([0, lim], [0, lim], color=figstyle.GRID, lw=0.6, zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("itpr.hmm bit score")
    ax.set_ylabel("ryr.hmm bit score")
    # Lower right is the only empty quadrant: the two families occupy the
    # two axes and the band between them is empty, which is the result.
    leg = ax.legend(loc="lower right", frameon=False, markerscale=3.5,
                    handletextpad=0.3, borderpad=0.2, labelspacing=0.35,
                    title="census v2 architecture call")
    leg.get_title().set_fontsize(figstyle.FS_TICK)
    figstyle.panel(ax, "a", "every record scored by both profiles")
    figstyle.despine(ax)
    ax.annotate(f"{REL_MARGIN:.0%} no-call band —\nempty", ha="center",
                xy=(lim * 0.50, lim * 0.50), xytext=(lim * 0.60, lim * 0.80),
                fontsize=figstyle.FS_TICK, color=figstyle.MUTED,
                arrowprops=dict(arrowstyle="-", lw=0.5,
                                color=figstyle.FAINT))

    # Panel b — the margin's distribution, which is what the call rests on.
    for call in ("ITPR", "RYR"):
        vals = [_fnum(r, "rel_margin") for r in rows
                if r["assignment"] == call]
        ax2.hist(vals, bins=40, range=(0, 1), histtype="stepfilled",
                 alpha=0.75, color=CALL_COLOUR[call], label=call, lw=0)
    ax2.axvline(REL_MARGIN, color=figstyle.INK, lw=0.8, ls=(0, (3, 2)))
    ax2.annotate(f"D7 band\n{REL_MARGIN:.0%}", xy=(REL_MARGIN, 0.92),
                 xycoords=("data", "axes fraction"),
                 xytext=(4, 0), textcoords="offset points",
                 fontsize=figstyle.FS_TICK, color=figstyle.MUTED, va="top")
    ax2.set_yscale("log")
    ax2.set_xlabel("relative bit-score margin  (win − lose) / win")
    ax2.set_ylabel("records")
    ax2.legend(frameon=False, loc="upper center")
    figstyle.panel(ax2, "b", "how far apart the two profiles land")
    figstyle.despine(ax2)
    figstyle.hgrid(ax2)
    fig.tight_layout(pad=0.5, w_pad=1.4)
    return fig, "profile_separation"


# ------------------------------------------------------ 2. jackhmmer curve
def fig_jackhmmer_convergence():
    import matplotlib.pyplot as plt

    rows = read_tsv(CENSUS_V3_DIR / "convergence.tsv")
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(figstyle.W_FULL, 2.6))
    tags = sorted({r["seed_tag"] for r in rows})
    colours = [figstyle.PARALOG["ITPR1"], figstyle.PARALOG["ITPR2"],
               figstyle.PARALOG["ITPR3"]]
    for tag, colour in zip(tags, colours):
        pts = sorted((int(r["round"]), int(r["new_targets"]))
                     for r in rows if r["seed_tag"] == tag)
        killed = {r["kill_rule"] for r in rows if r["seed_tag"] == tag} - {""}
        label = tag + (f"  [D10 {'/'.join(killed)}]" if killed else "")
        ax.plot([p[0] for p in pts], [max(p[1], 0.4) for p in pts],
                marker="o", ms=3.2, lw=1.1, color=colour, label=label)
        sis = sorted((int(r["round"]), float(r["sister_rise"]),
                      float(r["sister_frac"])) for r in rows
                     if r["seed_tag"] == tag)
        ax2.plot([s[0] for s in sis], [s[1] * 100 for s in sis],
                 marker="o", ms=3.2, lw=1.1, color=colour,
                 label=f"{tag}  (baseline {sis[0][2]:.0%})")

    ax.set_yscale("log")
    ax.set_xlabel("jackhmmer iteration")
    ax.set_ylabel("new targets included")
    ax.legend(frameon=False, fontsize=figstyle.FS_TICK)
    figstyle.panel(ax, "a", "the completeness curve")
    figstyle.despine(ax)
    figstyle.hgrid(ax)

    from scripts.s3_kill import MAX_SISTER_RISE
    ax2.axhline(MAX_SISTER_RISE * 100, color="#b3261e", lw=0.9, ls=(0, (3, 2)))
    ax2.annotate(f"D10 K1 limit  +{MAX_SISTER_RISE * 100:.0f} points",
                 xy=(0.98, MAX_SISTER_RISE * 100),
                 xycoords=("axes fraction", "data"),
                 xytext=(0, 3), textcoords="offset points", ha="right",
                 fontsize=figstyle.FS_TICK, color="#b3261e")
    ax2.axhline(0, color=figstyle.GRID, lw=0.6)
    # Scale to the data, not to the rule: a run that drifts off-family drives
    # the sister share *down* (the accretion dilutes it), and a fixed window
    # around K1's ceiling would clip exactly the trace worth seeing.
    lo = min(float(r["sister_rise"]) for r in rows) * 100
    ax2.set_ylim(min(lo - 3, -6), MAX_SISTER_RISE * 100 + 3)
    ax2.set_xlabel("jackhmmer iteration")
    ax2.set_ylabel("RyR share of the model,\nchange from round 1 (points)")
    ax2.legend(frameon=False, fontsize=figstyle.FS_TICK,
               loc="upper left")
    figstyle.panel(ax2, "b", "did iterating drift it into the sister family?")
    figstyle.despine(ax2)
    figstyle.hgrid(ax2)
    fig.tight_layout(pad=0.5, w_pad=1.6)
    return fig, "jackhmmer_convergence"


# -------------------------------------------------- 3. instrument agreement
def fig_instrument_agreement():
    import matplotlib.pyplot as plt

    rows = read_tsv(CENSUS_V3_DIR / "census_v3.tsv")
    order = ["both", "architecture", "profile", "neither"]
    label = {"both": "both instruments", "architecture": "architecture only",
             "profile": "profile only", "neither": "neither"}
    calls = ["ITPR", "RYR", "unassigned", "conflict"]
    counts = Counter((r["instruments"], r["call"]) for r in rows)

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(figstyle.W_FULL, 2.7),
        gridspec_kw={"width_ratios": [1.0, 1.0]})
    left = [0.0] * len(order)
    for call in calls:
        vals = [counts.get((inst, call), 0) for inst in order]
        if not any(vals):
            continue
        ax.barh(range(len(order)), vals, left=left, height=0.62,
                color=CALL_COLOUR[call], label=call, linewidth=0)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([label[o] for o in order])
    ax.invert_yaxis()
    ax.set_xlabel("census v3 records")
    ax.legend(frameon=False, ncol=2, fontsize=figstyle.FS_TICK)
    for i, tot in enumerate(left):
        if tot:
            ax.annotate(f"{int(tot):,}", xy=(tot, i), xytext=(3, 0),
                        textcoords="offset points", va="center",
                        fontsize=figstyle.FS_TICK, color=figstyle.MUTED)
    ax.set_xlim(0, max(left) * 1.16)
    figstyle.panel(ax, "a", "which instrument calls each record")
    figstyle.despine(ax)
    figstyle.hgrid(ax, axis="x")

    # Panel b — what the profile changed, transition by transition.
    changes = read_tsv(CENSUS_V3_DIR / "call_changes.tsv")
    trans = Counter(f"{r['v2_call']} → {r['call']}" for r in changes)
    items = trans.most_common(6)[::-1]
    if items:
        colours = [CALL_COLOUR.get(k.split(" → ")[1], figstyle.FAINT)
                   for k, _ in items]
        ax2.barh(range(len(items)), [v for _, v in items], height=0.62,
                 color=colours, linewidth=0)
        ax2.set_yticks(range(len(items)))
        ax2.set_yticklabels([k for k, _ in items])
        for i, (_, v) in enumerate(items):
            ax2.annotate(f"{v:,}", xy=(v, i), xytext=(3, 0),
                         textcoords="offset points", va="center",
                         fontsize=figstyle.FS_TICK, color=figstyle.MUTED)
        ax2.set_xlim(0, max(v for _, v in items) * 1.20)
    ax2.set_xlabel("records whose call changed from census v2")
    figstyle.panel(ax2, "b", "what the second instrument moved")
    figstyle.despine(ax2)
    figstyle.hgrid(ax2, axis="x")
    fig.tight_layout(pad=0.5, w_pad=1.4)
    return fig, "instrument_agreement"


# ------------------------------------------------- 4. per-proteome copy number
def fig_proteome_copy_number():
    import matplotlib.pyplot as plt

    rows = read_tsv(CENSUS_V3_DIR / "census_v3.tsv")
    manifest = read_tsv(HMM_SWEEP_DIR / "proteome_manifest.tsv")
    # A handful of reference proteomes share a taxon id (different strains
    # or assemblies of one organism), and the census keys on taxon, so the
    # unit here is the taxon, not the proteome file. Saying so matters:
    # 763 proteomes cover 758 distinct taxa.
    sizes: dict[int, int] = {}
    for m in manifest:
        tid = int(m["Organism Id"])
        sizes[tid] = sizes.get(tid, 0) + int(m["Protein count"])
    per_taxon: Counter[int] = Counter()
    for r in rows:
        if r["call"] == "ITPR" and r["taxon_id"] and int(r["taxon_id"]) in sizes:
            per_taxon[int(r["taxon_id"])] += 1

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(figstyle.W_FULL, 2.5))
    counts = [per_taxon.get(t, 0) for t in sizes]
    hist = Counter(counts)
    top = 20                       # the tail is long and thin; fold it in
    ks = list(range(0, top + 1))
    vals = [hist.get(k, 0) for k in ks] + [sum(v for k, v in hist.items()
                                               if k > top)]
    # A gap before the folded tail bar, so its tick cannot collide with 20.
    xs = list(range(top + 1)) + [top + 2]
    ax.bar(xs, vals, width=0.78, linewidth=0,
           color=[("#b3261e" if k == 0 else figstyle.PARALOG["ITPR1"])
                  for k in ks] + [figstyle.FAINT])
    ax.set_xticks(list(range(0, top + 1, 5)) + [top + 2])
    ax.set_xticklabels([str(k) for k in range(0, top + 1, 5)] + [f">{top}"])
    for x, v in zip(xs, vals):
        if v >= 15:
            ax.annotate(f"{v}", xy=(x, v), xytext=(0, 2),
                        textcoords="offset points", ha="center",
                        fontsize=figstyle.FS_TICK, color=figstyle.MUTED)
    ax.annotate(f"{hist.get(0, 0)} taxa with none", xy=(0, hist.get(0, 0)),
                xytext=(0.10, 0.80), textcoords="axes fraction",
                fontsize=figstyle.FS_TICK, color="#b3261e",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#b3261e"))
    ax.set_xlabel("ITPR records in the proteome")
    ax.set_ylabel("taxa")
    ax.set_ylim(0, max(vals) * 1.18)
    figstyle.panel(ax, "a", f"{len(sizes)} vertebrate taxa "
                            f"({len(manifest)} reference proteomes)")
    figstyle.despine(ax)
    figstyle.hgrid(ax)

    # Panel b — the record count is a count of *annotations*, so it has to be
    # read against how deeply the proteome is annotated at all. A taxon with
    # 40 ITPR records is a well-covered genome, not a family expansion; a
    # taxon with none is, on this evidence, a thin gene set.
    xs2 = [sizes[t] for t in sizes]
    ys2 = [per_taxon.get(t, 0) for t in sizes]
    ax2.scatter(xs2, ys2, s=5.0, linewidths=0, alpha=0.5,
                color=figstyle.PARALOG["ITPR1"], rasterized=True)
    zero = [(sizes[t], 0) for t in sizes if per_taxon.get(t, 0) == 0]
    if zero:
        ax2.scatter([p[0] for p in zero], [p[1] for p in zero], s=11.0,
                    linewidths=0, color="#b3261e",
                    label=f"no ITPR record ({len(zero)})")
        med = sorted(p[0] for p in zero)[len(zero) // 2]
        ax2.legend(frameon=False, loc="upper left",
                   fontsize=figstyle.FS_TICK, markerscale=1.4)
        ax2.annotate(f"median gene set of those 15:\n{med:,} proteins",
                     xy=(med, 0), xytext=(-14, 30),
                     textcoords="offset points",
                     ha="right", fontsize=figstyle.FS_TICK, color="#b3261e",
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="#b3261e"))
    ax2.set_xscale("log")
    ax2.set_xlabel("proteins in the reference proteome")
    ax2.set_ylabel("ITPR records")
    figstyle.panel(ax2, "b", "copy number tracks annotation depth")
    figstyle.despine(ax2)
    figstyle.hgrid(ax2)
    fig.tight_layout(pad=0.5, w_pad=1.4)
    return fig, "proteome_copy_number"


FIGURES = {
    "profile_separation": fig_profile_separation,
    "jackhmmer_convergence": fig_jackhmmer_convergence,
    "instrument_agreement": fig_instrument_agreement,
    "proteome_copy_number": fig_proteome_copy_number,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for slug in FIGURES:
            print(" ", slug)
        return 0

    figstyle.use()
    import matplotlib.pyplot as plt
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for slug in (args.only or list(FIGURES)):
        fig, stem = FIGURES[slug]()
        paths = figstyle.save(fig, FIG_DIR / stem)
        plt.close(fig)
        print(f"[s3_figs] {stem} → {', '.join(p.name for p in paths)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
