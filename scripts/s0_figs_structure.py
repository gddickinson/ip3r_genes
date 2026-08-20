"""Figures for §2 — the subunit's domains, and the channel it builds.

`domain_architecture`  every signature that defines an IP3R is also in a RyR
`channel_structure`    the tetramer, its pore and the distance IP3 acts across

Both render from committed tables only (`domain_coords.tsv`, `structure_*`).
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch, Rectangle

import figstyle
import s0_fig_lib as lib

ORDER = ["ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "RYR3"]


def domain_architecture():
    figstyle.use()
    doms = lib.load_tsv(lib.DATA / "domain_coords.tsv")
    by_sym: dict[str, list] = {}
    length = {}
    for d in doms:
        by_sym.setdefault(d["symbol"], []).append(d)
        length[d["symbol"]] = int(d["length_aa"])

    fig = plt.figure(figsize=(lib.W_REVIEW, 3.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.52)
    ax = fig.add_subplot(gs[0])

    ymap = {s: len(ORDER) - 1 - i for i, s in enumerate(ORDER)}
    for sym in ORDER:
        y = ymap[sym]
        # No in-bar text: at this scale a 200 aa domain is 4 mm wide and any
        # label overflows it. Panel b is the key, in the same colours.
        lib.domain_track(ax, by_sym[sym], y, length[sym], label_min_aa=10**9)
        ax.text(-140, y, sym, ha="right", va="center",
                fontsize=figstyle.FS_LABEL, color=figstyle.INK)
        ax.text(length[sym] + 90, y, f"{length[sym]:,} aa", ha="left",
                va="center", fontsize=figstyle.FS_NOTE, color=figstyle.MUTED)
    ax.set_xlim(-780, 6050)
    ax.set_ylim(-0.72, len(ORDER) - 0.10)
    ax.set_yticks([])
    ax.set_xticks(range(0, 5001, 1000))
    ax.set_xlabel("residue")
    lib.domain_legend(ax, loc="upper right", bbox_to_anchor=(1.005, 1.06),
                      fontsize=figstyle.FS_NOTE)
    figstyle.despine(ax, keep=("bottom",))
    figstyle.panel(ax, "a", "drawn to a common scale from InterPro coordinates")

    # ---- b: the same table, counted. Columns ordered as they occur along the
    # chain — shared signatures first, then the pore, then the RyR-only set.
    counts: dict[tuple, int] = {}
    first_at: dict[str, float] = {}
    klass: dict[str, str] = {}
    for d in doms:
        key = (d["symbol"], d["pfam"])
        counts[key] = counts.get(key, 0) + 1
        klass[d["pfam"]] = d["shared_class"]
        rel = int(d["start"]) / int(d["length_aa"])
        first_at[d["pfam"]] = min(first_at.get(d["pfam"], 9.9), rel)
    rank = {"shared": 0, "generic": 1, "ryr_only": 2}
    pfams = sorted(klass, key=lambda p: (rank.get(klass[p], 3), first_at[p]))

    ax2 = fig.add_subplot(gs[1])
    for xi, pf in enumerate(pfams):
        for yi, sym in enumerate(ORDER):
            n = counts.get((sym, pf), 0)
            y = len(ORDER) - 1 - yi
            if n:
                ax2.add_patch(Rectangle((xi - 0.42, y - 0.40), 0.84, 0.80,
                                        facecolor=lib.DOMAIN_COLOUR.get(pf),
                                        linewidth=0))
                ax2.text(xi, y, str(n), ha="center", va="center",
                         fontsize=figstyle.FS_NOTE, color="#ffffff")
            else:
                ax2.plot([xi], [y], marker=".", ms=1.6, color=figstyle.GRID)
    ax2.set_xticks(range(len(pfams)))
    ax2.set_xticklabels([f"{lib.DOMAIN_SHORT.get(p, p)}\n{p}" for p in pfams],
                        fontsize=figstyle.FS_NOTE - 0.6)
    ax2.set_yticks(range(len(ORDER)))
    ax2.set_yticklabels(list(reversed(ORDER)), fontsize=figstyle.FS_TICK)
    ax2.set_xlim(-0.6, len(pfams) - 0.4)
    ax2.set_ylim(-0.6, len(ORDER) - 0.4)
    ax2.tick_params(length=0)
    figstyle.despine(ax2, keep=())
    figstyle.panel(ax2, "b", "copies per subunit — blank means the signature "
                             "is absent")
    lib.provenance(fig, "measured", "InterPro / Pfam, re-derived")
    return lib.save(fig, "domain_architecture")


# ------------------------------------------------------------------ fig 2

def _domain_of(pos: int, doms: list[dict]) -> str | None:
    for d in doms:
        if int(d["start"]) <= pos <= int(d["end"]):
            return d["pfam"]
    return None


def _rot(xy: np.ndarray, deg: float) -> np.ndarray:
    t = np.deg2rad(deg)
    r = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return xy @ r.T


def channel_structure():
    figstyle.use()
    meta = lib.load_json(lib.DATA / "structure_meta.json")
    ca = lib.load_tsv(lib.DATA / "structure_ca.tsv")
    pore = lib.load_tsv(lib.DATA / "structure_pore.tsv")
    doms = [d for d in lib.load_tsv(lib.DATA / "domain_coords.tsv")
            if d["accession"] == meta["uniprot"]]

    res = np.array([int(r["resseq"]) for r in ca])
    xyz = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in ca])
    colours = [lib.DOMAIN_COLOUR.get(_domain_of(int(p), doms), "#cfcec8")
               for p in res]
    tm_lo, tm_hi = meta["tm_span_z_A"]
    ip3 = np.array(meta["ip3_sites"]["A"])

    fig = plt.figure(figsize=(lib.W_REVIEW, 4.05))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.42, 1.0],
                          height_ratios=[1.0, 0.82], wspace=0.30, hspace=0.42)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 1])

    # ---- a: side view, four subunits reconstructed by C4 rotation
    for ax, proj in ((axA, "side"), (axB, "top")):
        for turn in (0, 90, 180, 270):
            p = _rot(xyz[:, :2], turn)
            pts = (np.column_stack([p[:, 0], xyz[:, 2]]) if proj == "side"
                   else p)
            segs, cols = [], []
            for i in range(len(pts) - 1):
                if res[i + 1] - res[i] != 1:
                    continue
                segs.append([pts[i], pts[i + 1]])
                cols.append(colours[i])
            ax.add_collection(LineCollection(
                segs, colors=cols, linewidths=0.45,
                alpha=0.95 if proj == "side" else 0.75))

    axA.axhspan(tm_lo, tm_hi, color="#f2f1ec", zorder=0)
    axA.text(116, (tm_lo + tm_hi) / 2, "membrane", rotation=90, ha="center",
             va="center", fontsize=figstyle.FS_NOTE, color=figstyle.MUTED)
    for turn in (0, 90, 180, 270):
        p = _rot(ip3[None, :2], turn)[0]
        axA.plot([p[0]], [ip3[2]], marker="o", ms=4.2,
                 mfc="#b3261e", mec="#ffffff", mew=0.5, zorder=6)
    gate_z = meta["gate_z_A"]
    axA.plot([0], [gate_z], marker="_", ms=7, color="#14140f", zorder=6)

    # The review quotes "roughly 100 Å from the gate". That is the rise along
    # the axis; the site is 62 Å off it, so the through-space separation is
    # larger. The measuring bar sits in the margin and only the hypotenuse is
    # drawn on the particle, so neither obscures the trace.
    px = _rot(ip3[None, :2], 0)[0][0]
    bar = -122.0
    axA.annotate("", xy=(bar, ip3[2]), xytext=(bar, gate_z),
                 arrowprops=dict(arrowstyle="<->", lw=0.7, color="#b3261e",
                                 shrinkA=0, shrinkB=0))
    for zz in (ip3[2], gate_z):
        axA.plot([bar - 3.5, bar + 3.5], [zz, zz], lw=0.7, color="#b3261e")
    axA.text(bar - 6, (ip3[2] + gate_z) / 2,
             f"{meta['ip3_axial_rise_above_gate_A']:.0f} Å along the axis",
             rotation=90, ha="center", va="center",
             fontsize=figstyle.FS_NOTE, color="#b3261e")
    axA.plot([px, 0], [ip3[2], gate_z], color="#b3261e", lw=0.9, zorder=5)
    axA.text(px * 0.45, (ip3[2] + gate_z) / 2 - 2,
             f"{meta['ip3_to_gate_A']['A']:.0f} Å\nthrough space",
             fontsize=figstyle.FS_NOTE, color="#b3261e", ha="center",
             va="center", linespacing=1.2, zorder=7,
             bbox=dict(boxstyle="round,pad=0.18", fc="#ffffff", ec="none",
                       alpha=0.88))
    axA.annotate("IP$_3$ site", xy=(px, ip3[2]), xytext=(-104, 62),
                 fontsize=figstyle.FS_NOTE, color="#b3261e", ha="center",
                 arrowprops=dict(arrowstyle="-", lw=0.5, color="#b3261e"))
    axA.annotate("gate", xy=(0, gate_z), xytext=(34, gate_z + 3),
                 fontsize=figstyle.FS_NOTE, color=figstyle.INK,
                 arrowprops=dict(arrowstyle="-", lw=0.5,
                                 color=figstyle.MUTED))
    axA.set_xlim(-135, 135)
    axA.set_ylim(-125, 80)
    axA.set_xlabel("Å from the four-fold axis")
    axA.set_ylabel("Å along the four-fold axis")
    figstyle.panel(axA, "a", f"{meta['pdb_id']} — "
                   + meta["description"].replace("IP3", "IP$_3$"))
    axA.set_aspect("equal")

    axB.set_xlim(-135, 135)
    axB.set_ylim(-135, 135)
    axB.set_aspect("equal")
    axB.set_xticks([])
    axB.set_yticks([])
    figstyle.despine(axB, keep=())
    figstyle.panel(axB, "b", "viewed down the axis")

    # ---- c: the pore
    pz = np.array([float(r["z_along_axis"]) for r in pore])
    pr = np.array([float(r["min_heavy_atom_radius"]) for r in pore])
    axC.fill_between(pz, 0, pr, color="#dbe9fb", zorder=1)
    axC.plot(pz, pr, color="#184f95", lw=1.0, zorder=2)
    for key, lab, ha, dx in (("filter", "selectivity filter", "right", -2.5),
                             ("gate", "gate", "left", 2.5)):
        z, r = meta[f"{key}_z_A"], meta[f"{key}_min_radius_A"]
        axC.plot([z], [r], marker="v", ms=3.6, color="#b3261e", zorder=4)
        names = ", ".join(n.capitalize() for n in meta[f"{key}_lining_residues"])
        axC.annotate(f"{lab}\n{names}\n{r:.1f} Å", xy=(z + dx, r + 1.0),
                     ha=ha, va="bottom",
                     fontsize=figstyle.FS_NOTE - 0.5, color=figstyle.INK,
                     linespacing=1.25)
    axC.axvspan(tm_lo, tm_hi, color="#f2f1ec", zorder=0)
    axC.set_xlim(pz.min(), pz.max())
    axC.set_ylim(0, 17.5)
    axC.set_xlabel("Å along the four-fold axis")
    axC.set_ylabel("min. atom–axis\ndistance (Å)")
    figstyle.hgrid(axC)
    figstyle.panel(axC, "c", "the permeation path")

    handles = [Patch(facecolor=lib.DOMAIN_COLOUR[p],
                     label=lib.DOMAIN_SHORT[p])
               for p in ("PF08709", "PF02815", "PF01365", "PF08454", "PF00520")]
    handles.append(Patch(facecolor="#cfcec8", label="not in a Pfam domain"))
    axA.legend(handles=handles, loc="upper center", ncol=3, fontsize=5.9,
               bbox_to_anchor=(0.5, -0.13), columnspacing=0.8)
    lib.provenance(fig, "measured",
                   f"PDB {meta['pdb_id']}, {meta['resolution_A']} Å")
    return lib.save(fig, "channel_structure")
