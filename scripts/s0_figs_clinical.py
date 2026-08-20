"""Figures for §1 and §9 — the chronology, and where disease sits.

`discovery_timeline`  four decades, from a soluble messenger to a variant set
`disease_map`         the clinical panel placed on the three paralogues

Both render from `milestones.tsv` / `disease_sites.tsv`, whose citation keys
are validated against the bibliography by `s0_figdata_curated.py`. Variant
positions are drawn at the resolution the cited source gives — see that
script's docstring for what `point`, `domain` and `gene` promise.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

import figstyle
import s0_fig_lib as lib

LANE = {"found": ("#2a78d6", "what it is"),
        "works": ("#184f95", "how it works"),
        "disease": ("#b3261e", "what goes wrong")}
PARALOGS = ["ITPR1", "ITPR2", "ITPR3"]


def discovery_timeline():
    """Two columns of milestones over the field's own publication record.

    An earlier version put every milestone on a true-year spine and nudged
    colliding labels apart; three papers in 1989 made the leader lines cross
    and the displacement was more visible than the chronology. The dates are
    printed instead, and the *to-scale* temporal information is given properly
    by panel b — the year distribution of the 137 references the review is
    built from, which is real data rather than a decoration.
    """
    figstyle.use()
    rows = lib.load_tsv(lib.DATA / "milestones.tsv")
    rows.sort(key=lambda r: (int(r["year"]), r["ref_id"]))
    refs = lib.load_tsv(lib.S0 / "references.tsv")

    fig = plt.figure(figsize=(lib.W_REVIEW, 3.55))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.30], hspace=0.34)
    ax = fig.add_subplot(gs[0])
    axh = fig.add_subplot(gs[1])

    half = (len(rows) + 1) // 2
    cols = [rows[:half], rows[half:]]
    for ci, col in enumerate(cols):
        x0 = ci * 0.52
        for ri, r in enumerate(col):
            y = half - 1 - ri
            colour, _ = LANE[r["lane"]]
            ax.plot([x0], [y], marker="o", ms=3.4, mfc=colour, mec="#ffffff",
                    mew=0.5, zorder=3)
            ax.text(x0 + 0.022, y, r["year"], ha="left", va="center",
                    fontsize=figstyle.FS_NOTE - 0.4, color=figstyle.MUTED)
            ax.text(x0 + 0.072, y, r["label"], ha="left", va="center",
                    fontsize=figstyle.FS_NOTE - 0.2, color=figstyle.INK)
        ax.plot([x0, x0], [-0.4, half - 0.6], color=figstyle.GRID, lw=1.2,
                zorder=1)
    ax.set_xlim(-0.03, 1.04)
    ax.set_ylim(-0.75, half - 0.35)
    ax.set_xticks([])
    ax.set_yticks([])
    figstyle.despine(ax, keep=())
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3.4,
                              color=c, label=lab)
                       for c, lab in LANE.values()],
              loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.10),
              fontsize=figstyle.FS_NOTE)
    figstyle.panel(ax, "a", f"{rows[0]['year']}–{rows[-1]['year']}, "
                            f"{len(rows)} milestones")

    # ---- b: the bibliography's own year distribution, to scale
    years = [int(r["year"]) for r in refs if r["year"].isdigit()]
    lo, hi = min(years), max(years)
    edges = np.arange(lo, hi + 2)
    counts = np.histogram(years, bins=edges)[0]
    axh.bar(edges[:-1], counts, width=0.86, color="#9ec5f4", zorder=3)
    top = max(counts) * 1.30
    for r in rows:
        axh.plot([int(r["year"])], [top * 0.93], marker="v", ms=2.8,
                 color=LANE[r["lane"]][0], zorder=4)
    axh.set_xlim(lo - 1.2, hi + 1.2)
    axh.set_ylim(0, top)
    axh.set_ylabel("refs")
    axh.set_xlabel("year of publication")
    figstyle.hgrid(axh)
    figstyle.panel(axh, "b", f"the {len(years)} references this review cites, "
                             f"and where the milestones fall")
    lib.provenance(fig, "curated", "years read from references.tsv")
    return lib.save(fig, "discovery_timeline")


# ------------------------------------------------------------- disease map

def disease_map():
    figstyle.use()
    sites = lib.load_tsv(lib.DATA / "disease_sites.tsv")
    doms = lib.load_tsv(lib.DATA / "domain_coords.tsv")
    meta = lib.load_json(lib.DATA / "structure_meta.json")
    by_sym: dict[str, list] = {}
    length: dict[str, int] = {}
    for d in doms:
        by_sym.setdefault(d["symbol"], []).append(d)
        length[d["symbol"]] = int(d["length_aa"])

    fig, ax = plt.subplots(figsize=(lib.W_REVIEW, 3.65))
    row_h = 3.0
    for gi, sym in enumerate(PARALOGS):
        base = (len(PARALOGS) - 1 - gi) * row_h
        lib.domain_track(ax, by_sym[sym], base, length[sym], height=0.42,
                         label_min_aa=10 ** 9)
        ax.text(-70, base, sym, ha="right", va="center",
                fontsize=figstyle.FS_LABEL, color=figstyle.INK)
        rows = [s for s in sites if s["gene"] == sym]
        above = below = 0
        for s in rows:
            colour = lib.MECHANISM_COLOUR[s["mechanism"]]
            if s["kind"] == "point":
                x = int(s["where"])
                above += 1
                ytxt = base + 0.42 + 0.34 * above
                ax.plot([x, x], [base + 0.22, ytxt - 0.06], color=colour,
                        lw=0.6)
                ax.plot([x], [base + 0.22], marker="v", ms=3.4, color=colour)
                ax.text(x + 40, ytxt, s["label"], ha="left", va="center",
                        fontsize=figstyle.FS_NOTE - 0.5, color=colour)
            elif s["kind"] == "domain":
                hits = [d for d in by_sym[sym] if d["pfam"] == s["where"]]
                st, en = int(hits[0]["start"]), int(hits[0]["end"])
                above += 1
                ytxt = base + 0.42 + 0.34 * above
                ax.add_patch(Rectangle((st, base - 0.30), en - st, 0.60,
                                       facecolor="none", edgecolor=colour,
                                       linewidth=0.9, zorder=6))
                ax.plot([(st + en) / 2, (st + en) / 2],
                        [base + 0.30, ytxt - 0.06], color=colour, lw=0.6)
                ax.text((st + en) / 2 + 40, ytxt, s["label"], ha="left",
                        va="center", fontsize=figstyle.FS_NOTE - 0.5,
                        color=colour)
            else:
                below += 1
                yb = base - 0.52 - 0.40 * below
                ax.plot([1, length[sym]], [yb] * 2, color=colour, lw=2.0,
                        solid_capstyle="butt", alpha=0.55)
                ax.text(length[sym] + 60, yb, s["label"], ha="left",
                        va="center", fontsize=figstyle.FS_NOTE - 0.5,
                        color=colour)

    # The recurrent multisystem variant sits just past the measured gate —
    # a coincidence worth showing, since §9.3's mechanism is the tetramer.
    gate = sorted(int("".join(c for c in r if c.isdigit()))
                  for r in meta["gate_lining_residues"])
    var = int([s["where"] for s in sites
               if s["kind"] == "point" and s["where"] == "2524"][0])
    ax.annotate(f"{var - gate[-1]} residues past the gate measured in "
                f"{meta['pdb_id']}",
                xy=(var, -0.24), xytext=(2745, -0.30), ha="left",
                va="center",
                fontsize=figstyle.FS_NOTE - 0.5, color=figstyle.MUTED,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=figstyle.GRID))

    ax.set_xlim(-620, 4450)
    ax.set_ylim(-2.35, 2 * row_h + 3.45)
    ax.set_yticks([])
    ax.set_xticks(range(0, 2801, 500))
    ax.set_xlabel("residue")
    figstyle.despine(ax, keep=("bottom",))
    ax.legend(handles=[Patch(facecolor=c, label=k)
                       for k, c in lib.MECHANISM_COLOUR.items()]
                      + [Line2D([], [], marker="v", ls="none", ms=3.4,
                                color=figstyle.MUTED,
                                label="residue named by the source"),
                         Patch(facecolor="none", edgecolor=figstyle.MUTED,
                               label="localised to a domain only"),
                         Line2D([], [], lw=2.0, color=figstyle.MUTED,
                                alpha=0.55, label="not localised")],
              loc="upper center", ncol=4, fontsize=figstyle.FS_NOTE - 0.4,
              bbox_to_anchor=(0.5, 1.02), columnspacing=1.4)
    lib.provenance(fig, "curated", "positions at the resolution the sources give")
    return lib.save(fig, "disease_map")
