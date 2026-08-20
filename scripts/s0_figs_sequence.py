"""Figures drawn from the alignments — §3.1, §6, §7.1, §9.

`conservation_profile`  where the family is constrained, along human IP3R1
`alignment_windows`     the pore aligns with RyR residue for residue; the
                        IP3-binding site does not

Both render from `align_*.tsv`, which `s0_figdata_align.py` writes from the
committed control panels.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch, Rectangle

import figstyle
import s0_fig_lib as lib

#: Row order for the window figure: the ITPR grade, then the sister family.
ROW_ORDER = ["Q14643", "Q14571", "Q14573", "P29993", "Q9Y0A1",
             "P21817", "Q92736", "Q15413"]
ITPR_ROWS = ROW_ORDER[:5]


def _smooth(y: np.ndarray, win: int) -> np.ndarray:
    k = np.ones(win) / win
    return np.convolve(y, k, mode="same")


def conservation_profile():
    figstyle.use()
    cons = lib.load_tsv(lib.DATA / "align_conservation.tsv")
    doms = [d for d in lib.load_tsv(lib.DATA / "domain_coords.tsv")
            if d["symbol"] == "ITPR1"]
    blocks = lib.load_tsv(lib.DATA / "align_blocks.tsv")
    ident = lib.load_tsv(lib.DATA / "align_identity.tsv")
    ameta = lib.load_json(lib.DATA / "align_meta.json")

    pos = np.array([int(r["itpr1_pos"]) for r in cons])
    val = np.array([float(r["conservation"]) for r in cons])
    gap = np.array([float(r["gap_fraction"]) for r in cons])
    length = int(doms[0]["length_aa"])

    # Nested grids: the domain track must sit tight under the profile it
    # annotates, while the lower block needs room for its own titles.
    fig = plt.figure(figsize=(lib.W_REVIEW, 4.75))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.30], hspace=0.34)
    top = outer[0].subgridspec(2, 1, height_ratios=[1.0, 0.24], hspace=0.08)
    bot = outer[1].subgridspec(1, 2, width_ratios=[1.0, 0.70], wspace=0.46)
    ax = fig.add_subplot(top[0])
    axd = fig.add_subplot(top[1], sharex=ax)
    axm = fig.add_subplot(bot[0])
    axl = fig.add_subplot(bot[1])

    # ---- a: constraint along the chain
    ax.fill_between(pos, 0, _smooth(val, 25), color="#cde2fb", lw=0)
    ax.plot(pos, _smooth(val, 25), color="#184f95", lw=0.7)
    ax.plot(pos, _smooth(gap, 25), color="#eb6834", lw=0.6, alpha=0.85)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("conservation")
    ax.set_xlim(0, length)
    ax.tick_params(labelbottom=False)
    figstyle.hgrid(ax)

    # the sites measured on the structure, carried onto ITPR1 numbering
    marks: dict[str, list[int]] = {}
    for b in blocks:
        if b["is_marked"] == "1" and b["accession"] == "Q14643" \
                and b["itpr1_pos"]:
            marks.setdefault(b["kind"], []).append(int(b["itpr1_pos"]))
    from matplotlib.lines import Line2D
    site_handles = []
    for kind, colour, lab in (("ligand", "#b3261e", "IP$_3$ contact"),
                              ("pore", "#c08a2b", "filter / gate")):
        xs = sorted(set(marks.get(kind, [])))
        ax.plot(xs, [0.972] * len(xs), marker="v", ls="none", ms=3.0,
                color=colour, clip_on=False)
        site_handles.append(Line2D([], [], marker="v", ls="none", ms=3.0,
                                   color=colour, label=f"{lab} (measured)"))
    ax.legend(handles=[
        Patch(facecolor="#184f95", label="conservation (25-aa mean)"),
        Patch(facecolor="#eb6834", label="gap fraction")] + site_handles,
        loc="lower left", ncol=2, bbox_to_anchor=(0.0, -0.04),
        fontsize=figstyle.FS_NOTE)
    figstyle.panel(ax, "a", "25 IP$_3$R sequences aligned; human IP$_3$R1 "
                            "numbering")

    # ---- domain track under the profile, same x axis
    lib.domain_track(axd, doms, 0.0, length, height=0.55, label_min_aa=10**9)
    axd.set_ylim(-0.55, 0.55)
    axd.set_yticks([])
    axd.set_xlabel("residue in human IP$_3$R1")
    figstyle.despine(axd, keep=("bottom",))

    # ---- c: all-pairs identity across both families
    labels = {r["a"]: r["label_a"] for r in ident if r["set"] == "superfamily"}
    labels.update({r["b"]: r["label_b"] for r in ident
                   if r["set"] == "superfamily"})
    piv = {(r["a"], r["b"]): float(r["identity_covered"])
           for r in ident if r["set"] == "superfamily"}
    n = len(ROW_ORDER)
    m = np.full((n, n), np.nan)
    for i, a in enumerate(ROW_ORDER):
        for j, b in enumerate(ROW_ORDER):
            v = piv.get((a, b), piv.get((b, a)))
            if v is not None:
                m[i, j] = v * 100
    im = axm.imshow(m, cmap="Blues", vmin=0, vmax=100)
    for i in range(n):
        for j in range(n):
            if not np.isnan(m[i, j]):
                axm.text(j, i, f"{m[i, j]:.0f}", ha="center", va="center",
                         fontsize=figstyle.FS_NOTE - 1.2,
                         color="#ffffff" if m[i, j] > 55 else figstyle.INK)
    axm.set_xticks(range(n))
    axm.set_yticks(range(n))
    short = [labels[a].replace("H. sapiens", "human")
             .replace("D. melanogaster", "fly")
             .replace("C. elegans", "worm") for a in ROW_ORDER]
    axm.set_xticklabels(short, rotation=90, fontsize=figstyle.FS_NOTE - 0.8)
    axm.set_yticklabels(short, fontsize=figstyle.FS_NOTE - 0.8)
    axm.tick_params(length=0)
    axm.axhline(4.5, color=figstyle.INK, lw=0.7)
    axm.axvline(4.5, color=figstyle.INK, lw=0.7)
    cb = fig.colorbar(im, ax=axm, fraction=0.040, pad=0.02, shrink=0.80)
    cb.set_label("% identity", fontsize=figstyle.FS_NOTE)
    cb.ax.tick_params(labelsize=figstyle.FS_NOTE - 1)
    cb.outline.set_visible(False)
    figstyle.panel(axm, "b", "identity over mutually covered columns")

    # ---- d: what the paralogue block actually spans
    hum = [("Q14643", "Q14571"), ("Q14643", "Q14573"), ("Q14571", "Q14573")]
    cross = [(a, b) for a in ROW_ORDER[:3] for b in ROW_ORDER[5:]]
    grades = [(a, b) for a in ROW_ORDER[:3] for b in ROW_ORDER[3:5]]
    groups = [("human IP$_3$R\npairs", hum, "#2a78d6"),
              ("IP$_3$R vs the\ninvertebrate grade", grades, "#8a897f"),
              ("IP$_3$R vs RyR", cross, "#4a3aa7")]
    for xi, (name, pairs, colour) in enumerate(groups):
        vals = [100 * (piv.get(p) or piv.get((p[1], p[0]))) for p in pairs]
        axl.scatter([xi + 0.06 * (k - len(vals) / 2) for k in range(len(vals))],
                    vals, s=9, color=colour, zorder=3)
        axl.plot([xi - 0.28, xi + 0.28], [np.mean(vals)] * 2, color=colour,
                 lw=1.2, zorder=2)
        axl.text(xi, max(vals) + 4.0, f"{np.mean(vals):.0f}%", ha="center",
                 fontsize=figstyle.FS_NOTE, color=colour)
    axl.set_xticks(range(len(groups)))
    axl.set_xticklabels([g[0] for g in groups],
                        fontsize=figstyle.FS_NOTE - 0.4)
    axl.set_ylim(0, 100)
    axl.set_ylabel("% identity")
    axl.set_xlim(-0.6, len(groups) - 0.4)
    figstyle.hgrid(axl)
    figstyle.panel(axl, "c", "the three comparisons")
    lib.provenance(fig, "computed",
                   f"MAFFT {ameta['mafft'][0]['version']} on the committed "
                   f"control panel")
    return lib.save(fig, "conservation_profile")


# --------------------------------------------------------------- windows

def alignment_windows():
    figstyle.use()
    rows = lib.load_tsv(lib.DATA / "align_blocks.tsv")
    ameta = lib.load_json(lib.DATA / "align_meta.json")
    wins: list[str] = []
    for r in rows:
        if r["window"] not in wins:
            wins.append(r["window"])
    cells = {(r["window"], r["accession"], int(r["column"])): r for r in rows}
    cols = {w: sorted({int(r["column"]) for r in rows if r["window"] == w})
            for w in wins}
    label = {r["accession"]: r["label"] for r in rows}

    ligand = [w for w in wins if cells[(w, ROW_ORDER[0], cols[w][0])]["kind"]
              == "ligand"]
    pore = [w for w in wins if w not in ligand]

    fig = plt.figure(figsize=(lib.W_REVIEW, 3.05))
    gs = fig.add_gridspec(
        2, max(len(ligand), len(pore)), hspace=0.46, wspace=0.30,
        height_ratios=[1, 1])
    axes = []
    for r, group in enumerate((ligand, pore)):
        for c, w in enumerate(group):
            axes.append((fig.add_subplot(gs[r, c]), w, c == 0))

    for ax, w, first in axes:
        cc = cols[w]
        # Shade a cell where it matches the consensus of the IP3R rows, so the
        # RyR rows go pale exactly where the two families diverge.
        consensus = {}
        for col in cc:
            got = [cells[(w, a, col)]["residue"] for a in ITPR_ROWS]
            got = [g for g in got if g not in "-."]
            consensus[col] = max(set(got), key=got.count) if got else None
        for yi, acc in enumerate(ROW_ORDER):
            y = len(ROW_ORDER) - 1 - yi
            for xi, col in enumerate(cc):
                cell = cells[(w, acc, col)]
                res = cell["residue"]
                hit = res == consensus[col] and res not in "-."
                ax.add_patch(Rectangle(
                    (xi - 0.5, y - 0.5), 1, 1,
                    facecolor="#2a78d6" if hit else "#f2f1ec",
                    edgecolor="#ffffff", linewidth=0.4))
                if cell["is_marked"] == "1":
                    ax.add_patch(Rectangle(
                        (xi - 0.5, y - 0.5), 1, 1, facecolor="none",
                        edgecolor="#b3261e", linewidth=0.7, zorder=4))
                ax.text(xi, y, res, ha="center", va="center",
                        fontsize=figstyle.FS_NOTE - 1.0, zorder=5,
                        color="#ffffff" if hit else figstyle.MUTED,
                        family="monospace")
        ax.set_xlim(-0.5, len(cc) - 0.5)
        ax.set_ylim(-0.7, len(ROW_ORDER) - 0.3)
        ax.set_xticks([])
        ax.set_yticks(range(len(ROW_ORDER)) if first else [])
        if first:
            ax.set_yticklabels(
                [label[a].replace("H. sapiens", "human")
                 .replace("D. melanogaster", "fly")
                 .replace("C. elegans", "worm")
                 for a in reversed(ROW_ORDER)],
                fontsize=figstyle.FS_NOTE - 0.8)
        ax.tick_params(length=0)
        figstyle.despine(ax, keep=())
        first_pos = min(int(cells[(w, ROW_ORDER[0], c)]["itpr3_pos"])
                        for c in cc)
        ax.set_title(f"{w}\n(IP$_3$R3 {first_pos}–{first_pos + len(cc) - 1})",
                     fontsize=figstyle.FS_NOTE, loc="center", pad=3,
                     color=figstyle.INK, linespacing=1.25)

    fig.legend(handles=[
        Patch(facecolor="#2a78d6", label="matches the IP$_3$R consensus"),
        Patch(facecolor="#f2f1ec", label="differs"),
        Patch(facecolor="#ffffff", edgecolor="#b3261e", linewidth=1.0,
              label="measured site (IP$_3$ contact, filter or gate)")],
        loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.045),
        fontsize=figstyle.FS_NOTE)
    lib.provenance(fig, "computed",
                   f"MAFFT {ameta['mafft'][1]['version']}, "
                   f"windows from PDB 6DQN")
    return lib.save(fig, "alignment_windows")
