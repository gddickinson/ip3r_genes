"""S7 figure 1 — the rooted ML phylogram.

Reads  results/phylogeny/rooted.nwk        (s7_report.py)
       results/phylogeny/membership_audit.tsv
       results/msa_v2/representatives.tsv
Writes results/phylogeny/figures/tree_ml_rooted.{png,pdf}

Design. Branches in neutral ink; identity carried by a coloured dot at
each tip plus a labelled box behind each paralog clade — only the three
vertebrate paralogs and the RyR outgroup wear a hue, because scattered
tip dots in nine colours cannot clear the all-pairs CVD floor and the
tree does not need them to (position, the boxes and the species name
already say which is which). Support appears twice: a filled dot on
every branch clearing IQ-TREE's own joint threshold (UFBoot >= 95 and
SH-aLRT >= 80), and explicit `aLRT/UFBoot` numbers on the deep nodes the
paper's claims rest on.

**The reassigned tips are read from `membership_audit.tsv`, not listed
here**, so the ring markers on the figure and the constraint sets the AU
test was run on cannot drift apart (D13 applied to a figure).

Run:  python3 scripts/s7_figure.py    (matplotlib needed)
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from scripts.s7_lib import (  # noqa: E402
    PHYLO_DIR, Node, load_groups, parse_newick, support_of,
)

ROOTED = PHYLO_DIR / "rooted.nwk"
AUDIT = PHYLO_DIR / "membership_audit.tsv"
FIGS = PHYLO_DIR / "figures"
CORE = ("ITPR1", "ITPR2", "ITPR3")
BOXED = CORE + ("RYR",)
MIN_BRACKET = 5          # tips before a non-core run earns a bracket label


def read_audit() -> tuple[dict[str, str], set[str]]:
    """membership_audit.tsv -> (relabelled tip -> new group, unconstrained)."""
    relabel: dict[str, str] = {}
    free: set[str] = set()
    if not AUDIT.exists():
        return relabel, free
    with open(AUDIT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["rule"] == "reassigned":
                relabel[r["label"]] = r["assigned_to"]
            elif r["rule"] == "unconstrained":
                free.add(r["label"])
    return relabel, free


def ladderize(n: Node) -> int:
    if n.is_leaf:
        return 1
    sizes = [ladderize(c) for c in n.children]
    order = sorted(range(len(sizes)), key=lambda i: sizes[i])
    n.children = [n.children[i] for i in order]
    return sum(sizes)


def layout(root: Node):
    """x = distance from root, y = tip order; internals midway."""
    xs: dict[int, float] = {id(root): 0.0}
    ys: dict[int, float] = {}
    tip_y = [0]

    def down(n: Node) -> None:
        for c in n.children:
            xs[id(c)] = xs[id(n)] + max(c.length, 0.0)
            down(c)

    def up(n: Node) -> float:
        if n.is_leaf:
            ys[id(n)] = tip_y[0]
            tip_y[0] += 1
            return ys[id(n)]
        kid_ys = [up(c) for c in n.children]
        ys[id(n)] = (min(kid_ys) + max(kid_ys)) / 2.0
        return ys[id(n)]

    down(root)
    up(root)
    return xs, ys


def tip_text(label: str, groups: dict, species_count: Counter) -> str:
    r = groups.get(label, {})
    sp = (r.get("species", label) or label).split(" (")[0]
    if species_count[r.get("species", label)] > 1:
        acc = (r.get("accession", "") or "").split("|")[0]
        return f"{sp}  {acc}"
    return sp


def contiguous_runs(tips, label_of):
    """Maximal runs of the same group in tip order -> (group, i0, i1)."""
    runs, start = [], 0
    for i in range(1, len(tips) + 1):
        if i == len(tips) or label_of(tips[i]) != label_of(tips[start]):
            runs.append((label_of(tips[start]), start, i - 1))
            start = i
    return runs


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    import figstyle as fs
    fs.use()

    groups = load_groups()
    relabel, _free = read_audit()
    tree = parse_newick(ROOTED.read_text())
    ladderize(tree)
    xs, ys = layout(tree)
    tips = tree.leaves()
    n = len(tips)
    species_count = Counter(groups.get(t.name, {}).get("species", t.name)
                            for t in tips)

    def grp_of(t: Node) -> str:
        return relabel.get(t.name, groups.get(t.name, {}).get("group", ""))

    # Deep nodes carrying the claims: each boxed group's largest clade
    # under tree-corrected membership, plus its two ancestors (the
    # paralog splits the sister question is about).
    parent: dict[int, Node] = {}
    for nd in tree.walk():
        for c in nd.children:
            parent[id(c)] = nd
    corrected: dict[str, set[str]] = {g: set() for g in BOXED}
    for t in tips:
        g = grp_of(t)
        if g in corrected:
            corrected[g].add(t.name)
    annotate: set[int] = set()
    clade_node: dict[str, Node] = {}
    for g in BOXED:
        mem = corrected[g]
        best, best_n = None, 0
        for nd in tree.walk():
            if nd.is_leaf:
                continue
            ls = nd.leaf_names()
            if ls <= mem and len(ls) > best_n:
                best, best_n = nd, len(ls)
        if best is None:
            continue
        clade_node[g] = best
        annotate.add(id(best))
        up = best
        for _ in range(2):
            if id(up) in parent:
                up = parent[id(up)]
                annotate.add(id(up))

    fig_h = min(fs.H_MAX, n * 0.062 + 1.15)
    fig, ax = plt.subplots(figsize=(fs.W_FULL, fig_h))
    xmax = max(xs.values())
    x_tip_text = xmax * 1.012

    for nd in tree.walk():
        x0 = xs[id(nd)]
        if not nd.is_leaf:
            kid_ys = [ys[id(c)] for c in nd.children]
            ax.plot([x0, x0], [min(kid_ys), max(kid_ys)],
                    color=fs.INK, lw=0.5, solid_capstyle="round", zorder=2)
        for c in nd.children:
            ax.plot([x0, xs[id(c)]], [ys[id(c)], ys[id(c)]],
                    color=fs.INK, lw=0.5, solid_capstyle="round", zorder=2)

    for nd in tree.walk():
        if nd.is_leaf or nd is tree:
            continue
        alrt, uf = support_of(nd.name)
        if uf is not None and uf >= 95 and (alrt is None or alrt >= 80):
            ax.plot(xs[id(nd)], ys[id(nd)], "o", ms=1.8, mew=0,
                    color=fs.INK, zorder=4)
        if id(nd) in annotate and uf is not None and id(nd) not in {
                id(v) for v in clade_node.values()}:
            lab = f"{alrt:.0f}/{uf:.0f}" if alrt is not None else f"{uf:.0f}"
            ax.annotate(lab, (xs[id(nd)], ys[id(nd)]),
                        xytext=(-2.5, 2.0), textcoords="offset points",
                        ha="right", va="bottom", fontsize=fs.FS_NOTE,
                        color=fs.INK, zorder=6,
                        bbox=dict(boxstyle="square,pad=0.1", fc=fs.SURFACE,
                                  ec="none", alpha=0.9))

    for t in tips:
        g = grp_of(t)
        colour = fs.GROUP.get(g, fs.FAINT) if g in BOXED else fs.FAINT
        ax.plot(xs[id(t)], ys[id(t)], "o", ms=2.0, mew=0, color=colour,
                zorder=4)
        txt = tip_text(t.name, groups, species_count)
        if t.name in relabel:
            ax.plot(xs[id(t)], ys[id(t)], "o", ms=4.6, mfc="none",
                    mec=fs.INK, mew=0.6, zorder=5)
            txt += f"  → {relabel[t.name]}"
        ax.text(x_tip_text, ys[id(t)], txt, fontsize=3.9, va="center",
                color=fs.INK, fontstyle="italic")

    # Boxes go beyond the longest tip label; the label widths are only
    # known once drawn and depend on the x-limits, so iterate to a fixed
    # point rather than guessing an offset.
    xlim_r = xmax * 1.55
    x_box_r = x_box_lab = xmax * 1.3
    for _ in range(6):
        ax.set_xlim(-xmax * 0.01, xlim_r)
        fig.canvas.draw()
        inv = ax.transData.inverted()
        rend = fig.canvas.get_renderer()
        right = max((inv.transform((txt.get_window_extent(rend).x1, 0))[0]
                     for txt in ax.texts), default=xmax * 1.3)
        per_px = (xlim_r + xmax * 0.01) / max(ax.get_window_extent().width, 1)
        x_box_r = right + per_px * 0.06 * fig.dpi
        x_box_lab = x_box_r + per_px * 0.035 * fig.dpi
        new_r = x_box_lab + per_px * 0.40 * fig.dpi
        if abs(new_r - xlim_r) < xmax * 0.004:
            xlim_r = new_r
            break
        xlim_r = new_r
    ax.set_xlim(-xmax * 0.01, xlim_r)

    boxed_labels: list = []
    for g in BOXED:
        nd = clade_node.get(g)
        if nd is None:
            continue
        kid_ys = [ys[id(l)] for l in nd.leaves()]
        y0, y1 = min(kid_ys) - 0.45, max(kid_ys) + 0.45
        x0 = xs[id(parent.get(id(nd), nd))] - xmax * 0.004
        colour = fs.GROUP[g]
        for lw, fc, alpha, z in ((0.8, colour, 0.085, 0), (0.9, "none", 1, 3)):
            ax.add_patch(FancyBboxPatch(
                (x0, y0), x_box_r - x0, y1 - y0,
                boxstyle="round,pad=0,rounding_size=" + str(xmax * 0.012),
                linewidth=lw, edgecolor=colour, facecolor=fc, alpha=alpha,
                zorder=z,
                mutation_aspect=(x_box_r - x0) / max(y1 - y0, 1e-9)))
        alrt, uf = support_of(nd.name)
        sup = (f"{alrt:.0f}/{uf:.0f}"
               if (alrt is not None and uf is not None) else "")
        bt = ax.text(x_box_lab, (y0 + y1) / 2,
                     f"{g}\nn = {len(kid_ys)}" + (f"\n{sup}" if sup else ""),
                     rotation=90, ha="left", va="center",
                     fontsize=fs.FS_LABEL, fontweight="bold", color=colour,
                     zorder=6, multialignment="center", linespacing=1.2)
        # kept so the group labels below can step around them: a boxed
        # clade's label and a neighbouring group's occupy the same
        # column, and a short run next to a paralog clade overprints it.
        boxed_labels.append(bt)

    # A run earns a label if it is long enough — *or* if it is the
    # largest run its group has anywhere on the tree. Without the second
    # clause a group that never forms a run of MIN_BRACKET is invisible:
    # the plants (runs of 2 and 1) and the fungi (a run of 2) were
    # unlabelled, so two of the nine groups could not be found in the
    # figure at all, and the family's range across the kingdoms is one
    # of the results this project exists to report. The clause costs at
    # most one extra label per group, so a group cannot be missing and
    # the panel cannot fill with singletons either.
    # Every non-boxed group is marked, and marked **once**. A bracket
    # goes on each run of 2+ tips so the reader can see the group is
    # scattered; the *text* goes only on the group's largest run, and
    # names the group's total. Labelling every qualifying run instead
    # put three "protists" labels, a "plants" and a "5 of 19" at the
    # same x within a few tips of each other, and they overprinted into
    # mush — a legend that has to be decoded is not a legend. Where two
    # groups' labels would still overlap in y (rotated text is tall in
    # data units, and a 2-tip run is not), the later one steps outward
    # by one column, measured from the drawn text rather than guessed.
    runs = [r for r in contiguous_runs(tips, grp_of) if r[0] not in BOXED]
    total: dict[str, int] = {}
    for grp, i0, i1 in runs:
        total[grp] = total.get(grp, 0) + (i1 - i0 + 1)
    label_run: dict[str, tuple] = {}
    for r in runs:
        cur = label_run.get(r[0])
        if cur is None or (r[2] - r[1]) > (cur[2] - cur[1]):
            label_run[r[0]] = r
    # Occupied space is a *rectangle*, not a point: a boxed clade's
    # label is three rotated lines ("ITPR1 / n = 13 / 93/95") and so is
    # wider in x than one step. Recording it at a single x let the
    # neighbouring group label step once and still land on top of it.
    rend0 = fig.canvas.get_renderer()

    def _bbox(t):
        bb = t.get_window_extent(rend0).transformed(ax.transData.inverted())
        return (min(bb.y0, bb.y1), max(bb.y0, bb.y1),
                min(bb.x0, bb.x1), max(bb.x0, bb.x1))

    placed: list[tuple[float, float, float, float]] = [
        _bbox(t) for t in boxed_labels]
    for grp, i0, i1 in runs:
        n_run = i1 - i0 + 1
        if n_run < 2 and label_run[grp] != (grp, i0, i1):
            continue
        y0, y1 = ys[id(tips[i0])] - 0.4, ys[id(tips[i1])] + 0.4
        if n_run >= 2:
            ax.plot([x_box_r, x_box_r], [y0, y1], color=fs.FAINT, lw=1.0,
                    solid_capstyle="butt", zorder=3)
        if label_run[grp] != (grp, i0, i1):
            continue
        text = f"{fs.GROUP_LABEL.get(grp, grp)}\nn = {total[grp]}"
        yc = (y0 + y1) / 2
        t = ax.text(x_box_lab, yc, text, rotation=90, ha="left",
                    va="center", fontsize=fs.FS_NOTE, color=fs.MUTED,
                    zorder=6, multialignment="center", linespacing=1.2)
        ly0, ly1, lx0, lx1 = _bbox(t)
        gap, width = xmax * 0.012, max(lx1 - lx0, xmax * 0.02)
        x = x_box_lab
        for _ in range(8):
            hit = [q for q in placed
                   if not (ly1 < q[0] or ly0 > q[1])          # y overlap
                   and not (x + width < q[2] or x > q[3])]     # and x
            if not hit:
                break
            x = max(q[3] for q in hit) + gap
        t.set_x(x)
        placed.append((ly0, ly1, x, x + width))

    bar = round(xmax * 0.1, 1) or round(xmax * 0.1, 2)
    ax.plot([0.0, bar], [n + 1.5, n + 1.5], color=fs.INK, lw=1.1)
    ax.text(bar / 2, n + 2.6, f"{bar} substitutions/site", ha="center",
            va="top", fontsize=fs.FS_NOTE, color=fs.MUTED)
    # The ring entry appears only when a ring does. On this tree the
    # relabelling rule fires on nothing, and a key that names a marker
    # the figure does not carry asserts a correction that was never
    # made — the drift D13 exists to prevent, in a legend.
    key = [(0.22, 1.8, dict(mew=0, color=fs.INK),
            "UFBoot ≥ 95 and SH-aLRT ≥ 80")]
    if relabel:
        key.append((0.56, 4.4, dict(mfc="none", mec=fs.INK, mew=0.6),
                    "census label overturned by the tree"))
    for xf, ms, kw, label in key:
        ax.plot([xf], [n + 1.5], "o", ms=ms,
                transform=ax.get_yaxis_transform(), clip_on=False, **kw)
        ax.annotate(label, xy=(xf, n + 1.5),
                    xycoords=ax.get_yaxis_transform(),
                    xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=fs.FS_NOTE, color=fs.MUTED,
                    annotation_clip=False)

    ax.set_ylim(n + 3.0, -1.4)
    ax.axis("off")
    ax.set_title(f"Maximum-likelihood phylogeny of the ITPR family "
                 f"({n} proteins, rooted on the ryanodine receptors)",
                 fontsize=fs.FS_SUPTITLE, fontweight="bold", color=fs.INK,
                 pad=6, loc="left")

    FIGS.mkdir(parents=True, exist_ok=True)
    fs.save(fig, FIGS / "tree_ml_rooted")
    plt.close(fig)
    print(f"figure: {FIGS}/tree_ml_rooted.png (+.pdf), {n} tips")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
