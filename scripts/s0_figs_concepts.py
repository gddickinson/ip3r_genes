"""Mechanism figures for §3, §4 and §5.

`gating_logic`     the co-agonist bell, sequential binding, all-four-sites
`regulation_map`   every input in §4, grouped by mechanism and signed
`signal_hierarchy` blip -> puff -> wave, and frequency as the encoded variable

The two curve figures are **schematic**: they draw the shape of a published
relationship with plausible constants, not digitised data, and are tagged that
way on the canvas. `regulation_map` is **curated** — its content is
`regulators.tsv`, whose every row is checked against the bibliography.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch

import figstyle
import s0_fig_lib as lib

EFFECT_COLOUR = {"activates": "#0f7d3d", "inhibits": "#b3261e",
                 "biphasic": "#c08a2b", "context": "#8a897f"}
EFFECT_LABEL = {"activates": "increases release",
                "inhibits": "decreases release",
                "biphasic": "biphasic",
                "context": "sign depends on paralogue / cell type"}


def gating_logic():
    figstyle.use()
    fig, axes = plt.subplots(1, 3, figsize=(lib.W_REVIEW, 2.05),
                             gridspec_kw={"width_ratios": [1.0, 0.72, 0.95],
                                          "wspace": 0.44})
    ax, ax2, ax3 = axes

    # ---- a: the bell. Activation and inhibition as two opposed Hill terms.
    ca = np.logspace(-2, 2, 400)
    act = ca ** 2 / (0.2 ** 2 + ca ** 2)
    inh = 3.0 ** 3 / (3.0 ** 3 + ca ** 3)
    po = act * inh
    po = po / po.max()
    ax.plot(ca, po, color="#184f95", lw=1.2)
    ax.fill_between(ca, 0, po, color="#cde2fb", lw=0)
    peak = ca[po.argmax()]
    ax.axvline(peak, color=figstyle.MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.annotate("activation", xy=(0.030, 0.72), color="#184f95",
                fontsize=figstyle.FS_NOTE, ha="center")
    ax.annotate("inhibition", xy=(38, 0.72), color="#b3261e",
                fontsize=figstyle.FS_NOTE, ha="center")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("cytosolic [Ca$^{2+}$] (µM)")
    ax.set_ylabel("open probability\n(normalised)")
    figstyle.hgrid(ax)
    figstyle.panel(ax, "a", "Ca$^{2+}$ is co-agonist and inhibitor")

    # ---- b: the order of binding is the safety catch. Stacked, not in a
    # row: at this panel width a horizontal chain gives boxes taller than
    # they are wide, and "R·IP3·Ca2+" does not fit inside one.
    states = [("R", "#e5e4df"), ("R·IP$_3$", "#9ec5f4"),
              ("R·IP$_3$·Ca$^{2+}$", "#2a78d6"), ("open", "#184f95")]
    for i, (name, colour) in enumerate(states):
        y = len(states) - 1 - i
        ax2.add_patch(FancyBboxPatch(
            (0.06, y - 0.17), 0.88, 0.34,
            boxstyle="round,pad=0.012,rounding_size=0.06",
            facecolor=colour, edgecolor="none"))
        ax2.text(0.50, y, name, ha="center", va="center",
                 fontsize=figstyle.FS_NOTE, zorder=4,
                 color="#ffffff" if i >= 2 else figstyle.INK)
        if i:
            ax2.annotate("", xy=(0.50, y + 0.19), xytext=(0.50, y + 0.81),
                         arrowprops=dict(arrowstyle="->", lw=0.7,
                                         color=figstyle.MUTED))
    for i, lab in enumerate(("+ IP$_3$", "+ Ca$^{2+}$")):
        ax2.text(0.56, len(states) - 1.5 - i, lab, ha="left", va="center",
                 fontsize=figstyle.FS_NOTE - 0.3, color=figstyle.MUTED)
    ax2.text(0.50, -0.62, "Ca$^{2+}$ cannot act first — the ordering\n"
                          "is what stops spontaneous release",
             ha="center", va="center", fontsize=figstyle.FS_NOTE - 0.4,
             color=figstyle.MUTED, linespacing=1.3)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-0.95, 3.35)
    ax2.set_xticks([])
    ax2.set_yticks([])
    figstyle.despine(ax2, keep=())
    figstyle.panel(ax2, "b", "the two ligands bind in order")

    # ---- c: why "all four" makes a switch out of a graded input
    p = np.linspace(0, 1, 300)
    ax3.plot(p, 1 - (1 - p) ** 4, color="#a9a79e", lw=1.1,
             label="any one site")
    ax3.plot(p, p ** 4, color="#184f95", lw=1.4, label="all four sites")
    ax3.plot(p, p, color=figstyle.GRID, lw=0.8, ls=(0, (3, 2)))
    ax3.set_xlabel("occupancy per site")
    ax3.set_ylabel("channel activated")
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1.02)
    ax3.legend(loc="upper left", fontsize=figstyle.FS_NOTE - 0.4)
    figstyle.hgrid(ax3)
    figstyle.panel(ax3, "c", "the requirement is a threshold")
    lib.provenance(fig, "schematic", "shapes after the cited work; "
                                     "constants illustrative")
    return lib.save(fig, "gating_logic")


# ------------------------------------------------------------- regulation

def regulation_map():
    figstyle.use()
    rows = lib.load_tsv(lib.DATA / "regulators.tsv")
    groups: dict[str, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["mechanism"], []).append(r)

    per_line, lines = 4, []
    for mech, items in groups.items():
        for i in range(0, len(items), per_line):
            lines.append((mech if i == 0 else "", items[i:i + per_line],
                          items[i]["section"]))

    fig, ax = plt.subplots(figsize=(lib.W_REVIEW,
                                    0.44 * len(lines) + 0.62))
    for li, (mech, items, sec) in enumerate(lines):
        y = len(lines) - 1 - li
        if mech:
            ax.text(-0.06, y, mech, ha="right", va="center",
                    fontsize=figstyle.FS_LABEL, color=figstyle.INK)
            ax.text(-0.06, y - 0.30, f"§{sec}", ha="right", va="center",
                    fontsize=figstyle.FS_NOTE - 0.6, color=figstyle.MUTED)
        for ci, it in enumerate(items):
            colour = EFFECT_COLOUR[it["effect"]]
            ax.add_patch(FancyBboxPatch(
                (ci * 1.02 + 0.03, y - 0.26), 0.96, 0.52,
                boxstyle="round,pad=0.01,rounding_size=0.09",
                facecolor=colour, alpha=0.13, edgecolor=colour,
                linewidth=0.8))
            ax.text(ci * 1.02 + 0.51, y, it["name"], ha="center", va="center",
                    fontsize=figstyle.FS_NOTE - 0.5, color=figstyle.INK)
    ax.set_xlim(-1.30, per_line * 1.02 + 0.10)
    ax.set_ylim(-0.62, len(lines) - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    figstyle.despine(ax, keep=())
    ax.legend(handles=[Patch(facecolor=EFFECT_COLOUR[k], alpha=0.35,
                             edgecolor=EFFECT_COLOUR[k], label=EFFECT_LABEL[k])
                       for k in ("activates", "inhibits", "biphasic",
                                 "context")],
              loc="lower center", ncol=2, bbox_to_anchor=(0.44, -0.03),
              fontsize=figstyle.FS_NOTE)
    lib.provenance(fig, "curated", "every entry cited in §3–§4")
    return lib.save(fig, "regulation_map")


# --------------------------------------------------------------- hierarchy

def signal_hierarchy():
    figstyle.use()
    fig = plt.figure(figsize=(lib.W_REVIEW, 2.05))
    outer = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.80, 1.0], wspace=0.46)
    # The blip is ~10× smaller than the wave. Sharing one y axis would draw it
    # as a flat line, so each event gets its own band and the amplitudes are
    # stated instead of implied.
    stack = outer[0].subgridspec(3, 1, hspace=0.22)
    ax2 = fig.add_subplot(outer[1])
    ax3 = fig.add_subplot(outer[2])
    rng = np.random.default_rng(11)
    t = np.linspace(0, 10, 1200)

    def kick(centres, amp, tau, width=0.0):
        y = np.zeros_like(t)
        for c in centres:
            m = t >= c
            y[m] += amp * np.exp(-(t[m] - c) / tau)
            if width:
                y += amp * 0.35 * np.exp(-((t - c) / width) ** 2)
        return y

    rows = [("wave — the cell", kick([3.0], 9.0, 2.6, 0.9), "#184f95", "10×"),
            ("puff — a cluster", kick([2.4, 6.6], 3.0, 0.55, 0.12),
             "#2a78d6", "3×"),
            ("blip — one channel", kick([2.1, 5.4, 7.8], 0.9, 0.18, 0.05),
             "#a9a79e", "1×")]
    for i, (name, y, colour, rel) in enumerate(rows):
        axi = fig.add_subplot(stack[i])
        axi.plot(t, y + rng.normal(0, 0.02 * max(y), y.size), color=colour,
                 lw=0.8)
        axi.set_xlim(0, 10)
        axi.set_ylim(-0.12 * max(y), 1.35 * max(y))
        axi.set_yticks([])
        axi.text(0.02, 0.96, name, transform=axi.transAxes, va="top",
                 fontsize=figstyle.FS_NOTE - 0.3, color=colour)
        axi.text(0.99, 0.96, rel, transform=axi.transAxes, va="top",
                 ha="right", fontsize=figstyle.FS_NOTE - 0.5,
                 color=figstyle.MUTED)
        figstyle.despine(axi, keep=("bottom",))
        if i < 2:
            axi.tick_params(labelbottom=False)
        else:
            axi.set_xlabel("time (s)")
        if i == 0:
            figstyle.panel(axi, "a", "release is quantal, and recruits")

    # ---- b: the same three events, placed by scale
    pts = [("blip", 0.3, 0.02, "#a9a79e"), ("puff", 3.0, 0.6, "#2a78d6"),
           ("wave", 60.0, 20.0, "#184f95")]
    for name, space, time, colour in pts:
        ax2.scatter([space], [time], s=46, color=colour, zorder=3)
        ax2.annotate(name, xy=(space, time), xytext=(space * 1.9, time * 1.5),
                     fontsize=figstyle.FS_NOTE, color=colour)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.1, 400)
    ax2.set_ylim(0.005, 200)
    ax2.set_xlabel("spatial extent (µm)")
    ax2.set_ylabel("duration (s)")
    figstyle.hgrid(ax2, axis="both")
    figstyle.panel(ax2, "b", "four orders of magnitude")

    # ---- c: frequency, not amplitude, carries the dose
    for i, (dose, freq, colour) in enumerate([("low", 0.28, "#9ec5f4"),
                                              ("mid", 0.62, "#3987e5"),
                                              ("high", 1.15, "#184f95")]):
        spikes = np.arange(1.0, 10.0, 1.0 / freq)
        y = kick(spikes, 3.0, 0.35, 0.10)
        ax3.plot(t, y + i * 4.6, color=colour, lw=0.8)
        ax3.text(10.3, i * 4.6 + 1.4, f"{dose} agonist",
                 fontsize=figstyle.FS_NOTE, color=colour, va="center")
    ax3.set_xlim(0, 10)
    ax3.set_xlabel("time (s)")
    ax3.set_ylabel("[Ca$^{2+}$]")
    ax3.set_yticks([])
    figstyle.despine(ax3, keep=("bottom",))
    figstyle.panel(ax3, "c", "frequency encodes the dose")
    lib.provenance(fig, "schematic", "illustrative traces, not recordings")
    return lib.save(fig, "signal_hierarchy")
