"""The orientation figure: the pathway, the channel and the sequence.

`receptor_overview` is the figure a reader meets first. It answers three
questions a census of this family cannot be read without: what makes the
ligand, what the channel looks like, and what the protein's sequence actually
is at the two positions the rest of the document keeps returning to.

Its three panels do not share a provenance, so each says on the panel where it
comes from rather than leaving one corner tag to cover all three:

  a  **schematic.** The phosphoinositide pathway as the literature describes
     it. No data, not to scale.
  b  **measured.** Every dimension is read from `structure_meta.json`, which
     `s0_figdata_structure.py` measured on PDB 6DQN, and the domain track from
     `domain_coords.tsv`, which is InterPro's own per-accession coordinates.
     Nothing here is drawn to a remembered number.
  c  **computed.** Residues from `align_blocks.tsv`, the committed MAFFT
     alignment windows, anchored on the sites panel b measures rather than on
     a residue list from a paper.

The panels are deliberately complementary to the figures that follow rather
than a summary of them: Figure 1.2 is the same structure as a measurement with
its pore profile, Figure 1.5 is the domain architecture set against the sister
family, and Chapter 2's window figure is the sequence comparison across both
families. This one exists so that a reader has the object in mind before any
of those are argued.
"""

from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.patches import (Circle, FancyArrowPatch, FancyBboxPatch,
                                Polygon, Rectangle)

import figstyle
import s0_fig_lib as lib

#: Panel c shows the family only. The sister-family comparison is Chapter 2's
#: figure, and repeating it here would make this an argument rather than an
#: orientation.
SEQ_ROWS = ["Q14643", "Q14571", "Q14573", "P29993"]
SEQ_LABEL = {"Q14643": "ITPR1 human", "Q14571": "ITPR2 human",
             "Q14573": "ITPR3 human", "P29993": "Itpr fly"}
#: (window, heading, colour) — the colour ties each window to the domain it
#: sits in on panel b's track.
SEQ_WINDOWS = [("IP$_3$ contact 1", "in the IP$_3$-binding core",
                lib.DOMAIN_COLOUR["PF08709"]),
               ("selectivity filter", "in the pore",
                lib.DOMAIN_COLOUR["PF00520"])]

MEMBRANE = "#dcd9d0"
ER_FILL = "#eef3fb"
CYTO = "#f6f5f1"
CHANNEL = "#2a78d6"


def _arrow(ax, xy0, xy1, colour=figstyle.MUTED, lw=0.7, style="-|>",
           rad=0.0):
    ax.add_patch(FancyArrowPatch(
        xy0, xy1, arrowstyle=style, mutation_scale=6, lw=lw, color=colour,
        shrinkA=1.5, shrinkB=1.5,
        connectionstyle=f"arc3,rad={rad}"))


def _pathway(ax) -> None:
    """Panel a — what makes the ligand, and what the ligand opens."""
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    ax.add_patch(FancyBboxPatch(
        (3, 4), 94, 85, boxstyle="round,pad=0,rounding_size=6",
        facecolor=CYTO, edgecolor="#d5d3cc", lw=0.6))
    ax.add_patch(Rectangle((3, 84), 94, 5, facecolor=MEMBRANE,
                           edgecolor="none"))
    ax.text(5, 60, "cytosol", fontsize=figstyle.FS_NOTE - 0.6,
            color=figstyle.MUTED)

    # the plasma-membrane receptor and the enzyme it activates
    ax.add_patch(Rectangle((10, 82.5), 7, 8, facecolor="#8a897f",
                           edgecolor="none"))
    ax.text(13.5, 95, "agonist", ha="center",
            fontsize=figstyle.FS_NOTE - 0.6, color=figstyle.MUTED)
    _arrow(ax, (13.5, 94), (13.5, 91))
    ax.text(4, 71, "cell-surface\nreceptor",
            fontsize=figstyle.FS_NOTE - 0.8, color=figstyle.MUTED)
    _arrow(ax, (18, 83.5), (27.5, 77), rad=-0.15)
    ax.add_patch(FancyBboxPatch(
        (27, 70), 13, 7, boxstyle="round,pad=0,rounding_size=2",
        facecolor="#ffffff", edgecolor="#8a897f", lw=0.6))
    ax.text(33.5, 73.5, "PLC", ha="center", va="center",
            fontsize=figstyle.FS_NOTE, color=figstyle.INK)

    # PIP2 is cleaved into DAG, which stays in the membrane, and IP3, which
    # leaves it. Only the second is what this family responds to.
    ax.text(50, 86.5, "PIP$_2$", ha="center", va="center",
            fontsize=figstyle.FS_NOTE - 0.4, color=figstyle.MUTED)
    _arrow(ax, (40, 74), (47, 84), rad=0.18)
    ax.text(66, 86.5, "DAG", ha="center", va="center",
            fontsize=figstyle.FS_NOTE - 0.4, color=figstyle.MUTED)
    _arrow(ax, (55, 86.5), (61, 86.5))
    ax.text(50, 68, "IP$_3$", ha="center", va="center",
            fontsize=figstyle.FS_NOTE + 1.0, color=CHANNEL)
    _arrow(ax, (50, 83), (50, 71.5), colour=CHANNEL)
    _arrow(ax, (50, 64.5), (54, 51), colour=CHANNEL, rad=0.22)

    # the store, and the channel in its membrane
    ax.add_patch(FancyBboxPatch(
        (12, 15), 62, 29, boxstyle="round,pad=0,rounding_size=8",
        facecolor=ER_FILL, edgecolor="#c3d4ec", lw=0.7))
    ax.text(19, 24, "endoplasmic reticulum\nCa$^{2+}$ store",
            fontsize=figstyle.FS_NOTE - 0.6, color="#4d6a95", va="center")
    ax.add_patch(Rectangle((50, 40.5), 16, 6.5, facecolor=CHANNEL,
                           edgecolor="none"))
    for x in (52.5, 56.5, 60.5, 64.5):
        ax.add_patch(Rectangle((x - 1.2, 40.5), 0.8, 6.5,
                               facecolor="#ffffff", edgecolor="none"))
    ax.text(47, 43.5, "IP$_3$\nreceptor", ha="right", va="center",
            fontsize=figstyle.FS_NOTE, color=CHANNEL)

    _arrow(ax, (62, 47.5), (66, 60), colour="#b3261e", rad=-0.2)
    ax.text(68, 63, "Ca$^{2+}$", fontsize=figstyle.FS_NOTE, color="#b3261e")
    for x, y in ((60, 57), (71, 57), (75, 66), (64, 66)):
        ax.add_patch(Circle((x, y), 1.0, facecolor="#b3261e",
                            edgecolor="none", alpha=0.7))
    ax.text(96, 44, "effectors:\nmitochondria,\nkinases,\ntranscription",
            ha="right", va="top", fontsize=figstyle.FS_NOTE - 0.8,
            color=figstyle.MUTED)
    _arrow(ax, (78, 60), (84, 49), rad=-0.15)

    ax.annotate("the gene for this channel is what this thesis counts",
                xy=(58, 40), xytext=(9, 8), fontsize=figstyle.FS_NOTE - 0.4,
                color=CHANNEL, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=CHANNEL,
                                connectionstyle="arc3,rad=-0.25"))
    figstyle.panel(ax, "a", "the pathway that makes the ligand (schematic)")


def _channel(ax, meta: dict) -> None:
    """Panel b, upper — the receptor at the dimensions 6DQN was measured at.

    Every coordinate is read from the measurements file. The one drawing
    decision is the outline of the cytosolic mass, which is a shape rather
    than a measurement; its *extent* is the measured overall height and
    diameter, so the figure is to scale even where the silhouette is not.
    """
    z0, z1 = (float(v) for v in meta["tm_span_z_A"])
    gate_z = float(meta["gate_z_A"])
    filt_z = float(meta["filter_z_A"])
    lig_z = float(meta["ip3_sites"]["A"][2])
    lig_r = float(meta["ip3_radial_offset_from_axis_A"])
    height = float(meta["overall_height_A"])
    half = float(meta["overall_max_diameter_A"]) / 2
    top = z0 + height

    ax.set_xlim(-half - 26, half + 38)
    ax.set_ylim(z0 - 26, top + 20)
    ax.axis("off")

    ax.add_patch(Polygon(
        [(-half, top), (half, top), (half * 0.60, z1 + 22),
         (half * 0.26, z1), (-half * 0.26, z1), (-half * 0.60, z1 + 22)],
        closed=True, facecolor="#e8eef8", edgecolor="#bcd0ea", lw=0.7))
    ax.text(0, top - 14, "cytosolic solenoid", ha="center", va="center",
            fontsize=figstyle.FS_NOTE, color="#4d6a95")

    ax.add_patch(Rectangle((-half - 26, z0), 2 * half + 64, z1 - z0,
                           facecolor=MEMBRANE, edgecolor="none"))
    ax.text(-half - 22, (z0 + z1) / 2, "ER membrane", va="center",
            fontsize=figstyle.FS_NOTE - 0.8, color=figstyle.MUTED)
    ax.add_patch(Rectangle((-11, z0 - 8), 22, (z1 - z0) + 16,
                           facecolor="#ffffff", edgecolor="#bcd0ea", lw=0.7))

    for z, colour, name in ((gate_z, "#b3261e", "gate"),
                            (filt_z, "#c08a2b", "selectivity filter")):
        ax.plot([-11, 11], [z, z], color=colour, lw=1.8,
                solid_capstyle="butt")
        ax.annotate(name, xy=(11, z), xytext=(34, z),
                    fontsize=figstyle.FS_NOTE - 0.4, color=colour,
                    va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=colour))

    for sign in (-1, 1):
        ax.add_patch(Circle((sign * lig_r, lig_z), 7.5, facecolor=CHANNEL,
                            edgecolor="none", alpha=0.85))
    ax.annotate("IP$_3$ site", xy=(lig_r, lig_z),
                xytext=(half * 0.80, lig_z - 20),
                fontsize=figstyle.FS_NOTE - 0.4, color=CHANNEL, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=CHANNEL))

    # the two numbers that make the coupling problem concrete
    ax.plot([-lig_r, 0], [lig_z, gate_z], color=figstyle.MUTED, lw=0.6,
            ls=":")
    ax.text(-half * 0.82, (lig_z + gate_z) / 2,
            f"{float(meta['ip3_to_gate_mean_A']):.0f} Å\nthrough space",
            ha="right", va="center", fontsize=figstyle.FS_NOTE - 0.6,
            color=figstyle.MUTED)
    ax.annotate("", xy=(half + 6, lig_z), xytext=(half + 6, gate_z),
                arrowprops=dict(arrowstyle="<->", lw=0.6,
                                color=figstyle.MUTED))
    ax.text(half + 10, (lig_z + gate_z) / 2,
            f"{float(meta['ip3_axial_rise_above_gate_A']):.0f} Å\non the axis",
            fontsize=figstyle.FS_NOTE - 0.6, color=figstyle.MUTED,
            va="center")
    ax.text(0, z0 - 17, "ER lumen", ha="center",
            fontsize=figstyle.FS_NOTE - 0.8, color=figstyle.MUTED)
    figstyle.panel(ax, "b", f"the channel, measured on {meta['pdb_id']}")


def _track(ax, domains: list[dict]) -> None:
    """Panel b, lower — the same subunit as a linear architecture."""
    length = int(domains[0]["length_aa"])
    ax.set_xlim(-40, length + 40)
    ax.set_ylim(-0.55, 0.75)
    lib.domain_track(ax, domains, 0.0, length, height=0.42, label_min_aa=260)
    ax.set_yticks([])
    ax.set_xticks([1, 1000, 2000, length])
    ax.set_xticklabels(["1", "1,000", "2,000", f"{length:,}"])
    ax.tick_params(labelsize=figstyle.FS_TICK - 0.4, length=2)
    figstyle.despine(ax, keep=("bottom",))
    ax.set_xlabel("residue of one subunit (human ITPR3)",
                  fontsize=figstyle.FS_LABEL - 0.6)


def _sequence(ax, rows: list[dict], window: str, heading: str,
              colour: str, first: bool) -> None:
    """Panel c — the residues themselves at one anchored window."""
    cells = {(r["accession"], int(r["column"])): r for r in rows
             if r["window"] == window}
    cols = sorted({int(r["column"]) for r in rows if r["window"] == window})
    for yi, acc in enumerate(SEQ_ROWS):
        y = len(SEQ_ROWS) - 1 - yi
        for xi, col in enumerate(cols):
            cell = cells[(acc, col)]
            res = cell["residue"]
            same = len({cells[(a, col)]["residue"] for a in SEQ_ROWS}) == 1
            ax.add_patch(Rectangle((xi - 0.5, y - 0.5), 1, 1,
                                   facecolor=colour if same else "#f2f1ec",
                                   edgecolor="#ffffff", linewidth=0.4))
            if cell["is_marked"] == "1":
                ax.add_patch(Rectangle((xi - 0.5, y - 0.5), 1, 1,
                                       facecolor="none", edgecolor="#b3261e",
                                       linewidth=0.8, zorder=4))
            ax.text(xi, y, res, ha="center", va="center", zorder=5,
                    fontsize=figstyle.FS_NOTE - 1.0, family="monospace",
                    color="#ffffff" if same else figstyle.MUTED)
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_ylim(-0.6, len(SEQ_ROWS) - 0.4)
    ax.set_xticks([])
    ax.set_yticks(range(len(SEQ_ROWS)))
    if first:
        ax.set_yticklabels([SEQ_LABEL[a] for a in reversed(SEQ_ROWS)],
                           fontsize=figstyle.FS_TICK - 0.4)
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=0)
    figstyle.despine(ax, keep=())
    lo = cells[(SEQ_ROWS[0], cols[0])]["itpr1_pos"]
    hi = cells[(SEQ_ROWS[0], cols[-1])]["itpr1_pos"]
    ax.set_xlabel(f"ITPR1 {lo}–{hi}", fontsize=figstyle.FS_NOTE - 0.6,
                  labelpad=2, color=figstyle.MUTED)
    ax.set_title(f"{window}, {heading}", fontsize=figstyle.FS_NOTE + 0.2,
                 color=colour, pad=3)


def receptor_overview():
    figstyle.use()
    meta = lib.load_json(lib.DATA / "structure_meta.json")
    coords = [d for d in lib.load_tsv(lib.DATA / "domain_coords.tsv")
              if d["accession"] == "Q14573"]
    blocks = lib.load_tsv(lib.DATA / "align_blocks.tsv")

    fig = plt.figure(figsize=(lib.W_REVIEW, 4.45))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.60)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.02], wspace=0.20)
    right = top[1].subgridspec(2, 1, height_ratios=[1.0, 0.30], hspace=0.42)
    bottom = outer[1].subgridspec(1, 2, width_ratios=[13, 15], wspace=0.16)

    _pathway(fig.add_subplot(top[0]))
    _channel(fig.add_subplot(right[0]), meta)
    _track(fig.add_subplot(right[1]), coords)

    first_ax = None
    for i, (window, heading, colour) in enumerate(SEQ_WINDOWS):
        ax = fig.add_subplot(bottom[i])
        _sequence(ax, blocks, window, heading, colour, first=(i == 0))
        first_ax = first_ax or ax
    # The panel letter goes on the figure rather than through `panel()`,
    # which would overwrite the window heading each grid needs.
    pos = first_ax.get_position()
    fig.text(0.012, pos.y1 + 0.075, "c", fontsize=figstyle.FS_LETTER,
             fontweight="bold", va="bottom", color=figstyle.INK)
    fig.text(0.045, pos.y1 + 0.075,
             "the sequence at those two positions, in this family only",
             fontsize=figstyle.FS_TITLE, va="bottom", color=figstyle.INK)

    lib.provenance(fig, "schematic",
                   "panel a only; b measured on 6DQN, c from the "
                   "committed alignment")
    return lib.save(fig, "receptor_overview")
