"""The four S6 figures, from committed S6 tables only (D13, D19).

  1. `msa_identity_heatmap`   pairwise identity over mutually covered columns,
                              ordered by group, with the group strip beside it.
  2. `msa_conservation`       per-column conservation along the trimmed MSA
                              with human ITPR1's domain architecture mapped
                              **through the alignment** onto it.
  3. `msa_coverage`           per-sequence coverage of the trimmed MSA by
                              group — which tips are fragments.
  4. `msa_group_identity`     mean between-group identity, which is the
                              alignment's own preview of S7's sister
                              question, and its answer's error bar.

Nothing here recomputes an alignment or an identity: every number is read
out of `identity_covered.tsv`, `conservation.tsv`, `coverage.tsv` and
`representatives.tsv`, plus S0's committed InterPro domain coordinates for
figure 2. A figure that recomputed would be a second measurement wearing the
first one's caption.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s6_figures.py [--only slug]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.patches import Patch                         # noqa: E402

import figstyle as fs                                        # noqa: E402
import s6_rep_spec as spec                                   # noqa: E402
from s6_lib import MSA_DIR, read_fasta, read_tsv             # noqa: E402

FIGS = MSA_DIR / "figures"
DOMAINS = (Path(__file__).resolve().parents[1] / "results" / "s0_baseline" /
           "review_figures" / "domain_coords.tsv")
HUMAN_ITPR1 = "Q14643"


# ------------------------------------------------------------------- data

def load_identity(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    labels = rows[0][1:]
    m = [[float(v) for v in r[1:]] for r in rows[1:]]
    return labels, m


def load_meta() -> dict[str, dict]:
    return {r["label"]: r for r in read_tsv(MSA_DIR / "representatives.tsv")}


def present_groups(meta: dict[str, dict], labels) -> list[str]:
    have = {meta[l]["group"] for l in labels if l in meta}
    return [g for g in spec.GROUP_ORDER if g in have]


#: Tick-length names. `GROUP_LABEL` is written for a legend, where
#: "vertebrates, paralog unassigned" reads correctly; as an axis tick it
#: wraps into the panel below and pushes "ryanodine receptors" off the
#: right edge.
SHORT = {"vertebrate_basal": "vert.\nbasal", "invert_metazoa": "inverts",
         "RYR": "RyR", "ITPR1": "ITPR1", "ITPR2": "ITPR2", "ITPR3": "ITPR3",
         "plant": "plants", "protist": "protists", "fungi": "fungi"}


# --------------------------------------------------------------- figure 1

def fig_identity_heatmap() -> Path:
    labels, m = load_identity(MSA_DIR / "identity_covered.tsv")
    meta = load_meta()
    groups = present_groups(meta, labels)
    fig, (strip, ax) = plt.subplots(
        1, 2, figsize=(fs.W_FULL, 5.6),
        gridspec_kw={"width_ratios": [0.3, 12], "wspace": 0.03})
    im = ax.imshow(m, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Pairwise identity over mutually covered columns, "
                 f"{len(labels)} representatives", pad=6)
    for i, l in enumerate(labels):
        strip.barh(i, 1, height=1.0,
                   color=fs.GROUP.get(meta.get(l, {}).get("group", ""),
                                      "#cccccc"))
    strip.set_ylim(len(labels) - 0.5, -0.5)
    strip.set_xticks([])
    strip.set_yticks([])
    for s in strip.spines.values():
        s.set_visible(False)
    # Legend under the panel, not inside it: over a 134 x 134 heatmap an
    # inset legend sits on top of the data it is explaining.
    ax.legend([Patch(facecolor=fs.GROUP[g]) for g in groups],
              [fs.GROUP_LABEL[g] for g in groups], fontsize=6,
              loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=5,
              frameon=False, handlelength=1.0, columnspacing=1.2)
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("identity (covered columns)")
    return fs.save(fig, FIGS / "msa_identity_heatmap")[0]


# --------------------------------------------------------------- figure 2

def human_domain_columns() -> list[tuple[str, int, int, str]]:
    """Human ITPR1's Pfam spans, in *trimmed-MSA column* coordinates.

    The domain table is in residue coordinates on Q14643. Mapping it onto
    the alignment means walking that row of the trimmed MSA and counting
    ungapped positions — so a domain lands where the alignment actually put
    it, not where a proportional guess would.
    """
    if not DOMAINS.exists():
        return []
    trimmed = {k.split()[0]: v for k, v in
               read_fasta(MSA_DIR / "trimmed.fasta").items()}
    meta = load_meta()
    row = next((l for l in trimmed
                if meta.get(l, {}).get("accession") == HUMAN_ITPR1), None)
    if row is None:
        return []
    seq = trimmed[row]
    res_to_col: dict[int, int] = {}
    n = 0
    for col, ch in enumerate(seq):
        if ch != "-":
            n += 1
            res_to_col.setdefault(n, col)
    out = []
    for d in read_tsv(DOMAINS):
        if d["accession"] != HUMAN_ITPR1:
            continue
        cols = [res_to_col[r] for r in range(int(d["start"]), int(d["end"]) + 1)
                if r in res_to_col]
        if cols:
            out.append((d["label"], min(cols), max(cols), d["shared_class"]))
    return out


def fig_conservation() -> Path:
    cons = [float(r["conservation"]) for r in
            read_tsv(MSA_DIR / "conservation.tsv")]
    win = 25
    roll = [sum(cons[max(0, i - win // 2):i + win // 2 + 1])
            / len(cons[max(0, i - win // 2):i + win // 2 + 1])
            for i in range(len(cons))]
    doms = human_domain_columns()
    fig, (ax, dax) = plt.subplots(
        2, 1, figsize=(fs.W_FULL, 2.9), sharex=True,
        gridspec_kw={"height_ratios": [5, 1], "hspace": 0.12})
    x = range(1, len(cons) + 1)
    ax.plot(x, roll, color=fs.PARALOG["ITPR1"], lw=1.0)
    ax.fill_between(x, roll, color=fs.PARALOG["ITPR1"], alpha=0.15)
    ax.axhline(sum(cons) / len(cons), color="#52514e", lw=0.6, ls="--")
    ax.set_ylabel(f"conservation\n(rolling {win})")
    ax.set_ylim(0, 1)
    ax.set_xlim(1, len(cons))
    fs.despine(ax)
    fs.hgrid(ax)
    ax.set_title("Per-column conservation of the trimmed alignment, with "
                 "human ITPR1's domains mapped through it", pad=5)
    palette = {"shared": "#4a3aa7", "generic": "#a9a79e", "ryr_only": "#eb6834"}
    seen = []
    for label, c0, c1, cls in doms:
        w = c1 - c0 + 1
        dax.barh(0, w, left=c0 + 1, height=0.55,
                 color=palette.get(cls, "#8a897f"))
        # A 107 aa domain is ~40 columns wide here; its name does not fit
        # inside the bar, and drawn anyway it overprints its neighbours.
        if w >= 0.055 * len(cons):
            dax.text((c0 + c1) / 2 + 1, 0, label, ha="center", va="center",
                     fontsize=5, color="white")
        else:
            dax.text((c0 + c1) / 2 + 1, -0.62, label, ha="center", va="top",
                     fontsize=4.4, color="#3a3936", rotation=0)
        if cls not in seen:
            seen.append(cls)
    names = {"shared": "shared with RyR", "generic": "generic pore",
             "ryr_only": "RyR-only"}
    dax.legend([Patch(facecolor=palette.get(c, "#8a897f")) for c in seen],
               [names.get(c, c) for c in seen], fontsize=5, ncol=len(seen),
               loc="upper right", bbox_to_anchor=(1.0, 2.9), frameon=False,
               handlelength=1.0)
    dax.set_ylim(-1.4, 0.5)
    dax.set_yticks([])
    dax.set_xlabel("trimmed MSA column")
    for s in dax.spines.values():
        s.set_visible(False)
    return fs.save(fig, FIGS / "msa_conservation")[0]


# --------------------------------------------------------------- figure 3

def fig_coverage() -> Path:
    rows = read_tsv(MSA_DIR / "coverage.tsv")
    groups = [g for g in spec.GROUP_ORDER
              if any(r["group"] == g for r in rows)]
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.7))
    for i, g in enumerate(groups):
        vals = sorted(float(r["coverage"]) for r in rows if r["group"] == g)
        # Offsets alternate about the centre. Assigned monotonically with
        # rank instead, the column becomes a diagonal and reads as a trend
        # that is only the sort order.
        step = 0.30 / max((len(vals) - 1) // 2, 1)
        xs = [i + ((j + 1) // 2) * step * (1 if j % 2 else -1)
              for j in range(len(vals))]
        ax.scatter(xs, vals, s=11, color=fs.GROUP[g], zorder=3,
                   edgecolors="white", linewidths=0.3)
        med = vals[len(vals) // 2]
        ax.plot([i - 0.36, i + 0.36], [med, med], color="#3a3936", lw=1.1,
                zorder=4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([SHORT.get(g, fs.GROUP_LABEL[g]) for g in groups],
                       fontsize=6)
    ax.set_ylabel("coverage of trimmed MSA")
    ax.set_ylim(0, 1.02)
    ax.set_title("Per-sequence coverage of the trimmed alignment "
                 "(bar = group median)", pad=5)
    fs.despine(ax)
    fs.hgrid(ax)
    return fs.save(fig, FIGS / "msa_coverage")[0]


# --------------------------------------------------------------- figure 4

def group_identity_table() -> tuple[list[str], dict, dict]:
    labels, m = load_identity(MSA_DIR / "identity_covered.tsv")
    meta = load_meta()
    idx = defaultdict(list)
    for i, l in enumerate(labels):
        idx[meta.get(l, {}).get("group", "?")].append(i)
    groups = [g for g in spec.GROUP_ORDER if g in idx]
    means, spread = {}, {}
    for a in groups:
        for b in groups:
            vals = [m[i][j] for i in idx[a] for j in idx[b] if i != j]
            if vals:
                vals.sort()
                means[(a, b)] = sum(vals) / len(vals)
                spread[(a, b)] = (vals[len(vals) // 4],
                                  vals[3 * len(vals) // 4])
    return groups, means, spread


def fig_group_identity() -> Path:
    groups, means, spread = group_identity_table()
    trio = [g for g in spec.PARALOGS if g in groups]
    pairs = [(trio[i], trio[j]) for i in range(len(trio))
             for j in range(i + 1, len(trio))]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.6),
                                  gridspec_kw={"width_ratios": [1, 1.5]})
    for k, (a, b) in enumerate(pairs):
        lo, hi = spread[(a, b)]
        mu = means[(a, b)]
        ax.plot([k, k], [lo, hi], color="#52514e", lw=1.0, zorder=2)
        ax.scatter([k], [mu], s=42, color=fs.GROUP[a], zorder=3,
                   edgecolors=fs.GROUP[b], linewidths=1.6)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([f"{a} × {b}" for a, b in pairs], fontsize=6)
    ax.set_ylabel("mean identity (covered)")
    ax.set_title("Between-paralog identity", pad=5)
    ax.set_xlim(-0.5, len(pairs) - 0.5)
    fs.despine(ax)
    fs.hgrid(ax)

    order = [g for g in groups]
    m = [[means.get((a, b), float("nan")) for b in order] for a in order]
    im = ax2.imshow(m, cmap="viridis", vmin=0, vmax=1)
    ax2.set_xticks(range(len(order)))
    ax2.set_yticks(range(len(order)))
    ticks = [SHORT.get(g, fs.GROUP_LABEL[g]).replace("\n", " ")
             for g in order]
    ax2.set_xticklabels(ticks, rotation=45, ha="right", fontsize=6)
    ax2.set_yticklabels(ticks, fontsize=6)
    ax2.set_title("Mean identity between groups", pad=5)
    for i in range(len(order)):
        for j in range(len(order)):
            v = m[i][j]
            if v == v:
                ax2.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=4.6,
                         color="white" if v < 0.6 else "#1b1b18")
    plt.colorbar(im, ax=ax2, fraction=0.04, pad=0.02)
    return fs.save(fig, FIGS / "msa_group_identity")[0]


FIGURES = {
    "msa_identity_heatmap": fig_identity_heatmap,
    "msa_conservation": fig_conservation,
    "msa_coverage": fig_coverage,
    "msa_group_identity": fig_group_identity,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(k)
        return 0
    fs.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for slug in (args.only or list(FIGURES)):
        p = FIGURES[slug]()
        print(f"  {slug} -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
