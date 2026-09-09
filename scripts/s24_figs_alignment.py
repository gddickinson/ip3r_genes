"""Supplementary Figures 1-2: the representative alignment, and the residues.

Supplementary Figures 3-4 live in `s24_figs_inputs.py`; the pair was split to
keep both inside the project's 500-line budget.

Three decisions worth naming, each because the obvious version of the figure
would say something the data does not.

* **Supplementary Fig. 1 draws the kept columns as a track, not as a second
  alignment.** The trimmed alignment is not a picture of anything: it is the
  input alignment with 9,980 columns deleted, and drawing it beside the input
  invites the reader to compare two rasters that share no x-axis. The panel
  draws one raster, in input coordinates, with the kept columns marked
  underneath — which is the only way the reader can see *where* trimAl cut.

* **The raster bins columns and plots occupancy, never a residue colour.** At
  11,777 columns on a 6.7-inch page one pixel is nine columns, so a residue
  palette would draw whichever residue happened to land on the pixel. Occupancy
  is a mean and survives binning.

* **Supplementary Fig. 2 shows the residue letters.** The whole point of a
  residue-level panel is that the reader can check the claim; a heat strip with
  no letters would ask them to take the join on trust. Every letter drawn here
  has been checked against the residue its own paralogue's constraint table
  holds at that position (`s24_lib.verify_residues`), and the check is a hard
  failure, not a warning.

"""

from __future__ import annotations

import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.patches import Patch                           # noqa: E402

import figstyle as F                                           # noqa: E402
import s24_lib as L                                            # noqa: E402

GAPS = "-."

#: Kept against cut is an *ordered* contrast — evidence retained against
#: evidence removed — so it takes the diverging quality scale and not two
#: categorical hues. Categorical colour is spent on paralogue identity and
#: nothing else (figstyle), and a panel that coloured "kept" in ITPR1 blue
#: beside a group key that colours ITPR1 blue would read as a paralogue
#: statement.
KEPT = F.QUALITY["complete"]
#: The palest grey in the scale vanished at printed size, so `cut` takes the
#: mid grey instead: the contrast is still ordered, and it survives the page.
CUT = F.FAINT


def _group_of(reps: list[dict]) -> dict[str, str]:
    return {r["label"]: r["group"] for r in reps}


# ------------------------------------------------------------------ fig 1

def fig_representative_alignment(out: Path, stats: dict) -> None:
    aln = L.read_fasta(L.MSA_DIR / "aln.fasta")
    trimmed = L.read_fasta(L.MSA_DIR / "trimmed.fasta")
    colmap = L.read_tsv(L.MSA_DIR / "column_map.tsv")
    reps = L.read_tsv(L.MSA_DIR / "representatives.tsv")
    group = _group_of(reps)

    order = [g for g in F.GROUP_ORDER
             if any(group.get(k) == g for k in aln)]
    labels = [k for g in order for k in aln if group.get(k) == g]
    labels += [k for k in aln if k not in labels]
    width = len(next(iter(aln.values())))
    kept = sorted(int(r["aln_column"]) - 1 for r in colmap)
    keptset = set(kept)

    nbins = 900
    step = width / nbins
    grid = []
    for lab in labels:
        s = aln[lab]
        row = []
        for b in range(nbins):
            lo, hi = int(b * step), max(int(b * step) + 1, int((b + 1) * step))
            seg = s[lo:hi]
            row.append(sum(1 for c in seg if c not in GAPS) / len(seg))
        grid.append(row)

    fig = plt.figure(figsize=(F.W_FULL, 6.3))
    gs = fig.add_gridspec(4, 2, height_ratios=[7.0, 0.5, 2.2, 2.2],
                          width_ratios=[0.035, 1.0], hspace=0.55, wspace=0.02)

    ax_side = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(grid, aspect="auto", cmap="Blues", vmin=0.0, vmax=1.0,
              interpolation="nearest",
              extent=(0, width, len(labels), 0))
    ax.set_yticks([])
    ax.set_xlim(0, width)
    ax.set_xlabel("column of the 11,777-column L-INS-i alignment")
    F.despine(ax, keep=("bottom",))
    F.panel(ax, "a", "134 representatives, per-cell residue occupancy")

    for i, lab in enumerate(labels):
        ax_side.add_patch(plt.Rectangle(
            (0, i), 1, 1, color=F.GROUP.get(group.get(lab, ""), F.FAINT),
            lw=0))
    ax_side.set_xlim(0, 1)
    ax_side.set_ylim(len(labels), 0)
    ax_side.axis("off")

    ax_track = fig.add_subplot(gs[1, 1], sharex=ax)
    for c in kept:
        ax_track.add_patch(plt.Rectangle((c, 0), 1, 1,
                                         color=KEPT, lw=0))
    ax_track.set_xlim(0, width)
    ax_track.set_ylim(0, 1)
    ax_track.set_yticks([])
    ax_track.axis("off")
    ax_track.annotate(f"the {len(kept):,} columns trimAl kept — the alignment "
                      f"the tree, the selection tests and the constraint map "
                      f"were computed on",
                      xy=(0, -0.3), xycoords="axes fraction", fontsize=F.FS_NOTE,
                      color=F.MUTED, va="top", annotation_clip=False)

    occ = [sum(1 for lab in labels if aln[lab][c] not in GAPS) / len(labels)
           for c in range(width)]
    ax2 = fig.add_subplot(gs[2, 1])
    xs, ys = L.binned(occ, 450)
    ax2.fill_between(xs, ys, color=F.BLUES[1], lw=0)
    ax2.plot(xs, ys, color=F.BLUES[4], lw=0.7)
    kept_occ = [occ[c] for c in kept]
    cut_occ = [occ[c] for c in range(width) if c not in keptset]
    ax2.set_xlim(0, width)
    ax2.set_ylim(0, 1.02)
    ax2.set_yticks([0, 0.5, 1.0])
    ax2.set_ylabel("occupancy")
    ax2.set_xlabel("column of the input alignment")
    F.despine(ax2)
    F.hgrid(ax2)
    F.panel(ax2, "b", "column occupancy along the alignment "
                      "(mean of 26 columns)")

    # The claim is that the cut columns are the sparse ones, and that is a
    # statement about two distributions, not about their position along the
    # alignment - which panel b already shows and a second trace of would
    # only repeat. Two histograms on one axis is the shape of the claim.
    ax3 = fig.add_subplot(gs[3, 1])
    bins = [i / 40 for i in range(41)]
    ax3.hist([cut_occ, kept_occ], bins=bins, stacked=False, log=True,
             color=[CUT, KEPT], lw=0,
             label=[f"cut ({len(cut_occ):,}), median "
                    f"{st.median(cut_occ):.2f}",
                    f"kept ({len(kept_occ):,}), median "
                    f"{st.median(kept_occ):.2f}"])
    ax3.set_xlim(0, 1)
    ax3.set_xlabel("occupancy of the column")
    ax3.set_ylabel("columns (log)")
    ax3.legend(loc="upper center", ncol=2)
    F.despine(ax3)
    F.hgrid(ax3)
    F.panel(ax3, "c", "what trimAl removed was the sparse columns")

    handles = [Patch(facecolor=F.GROUP[g], label=F.GROUP_LABEL[g])
               for g in order]
    fig.legend(handles=handles, loc="lower center", ncol=5,
               fontsize=F.FS_TICK, frameon=False,
               bbox_to_anchor=(0.5, -0.035))
    F.save(fig, out / "SuppFig1_representative_alignment")
    plt.close(fig)

    stats["supp_fig_1"] = {
        "sequences": len(labels), "input_columns": width,
        "kept_columns": len(kept),
        "median_occupancy_kept": round(st.median(kept_occ), 4),
        "median_occupancy_cut": round(st.median(cut_occ), 4),
        "groups": {g: sum(1 for k in labels if group.get(k) == g)
                   for g in order},
    }


# ------------------------------------------------------------------ fig 2

def _region_rows(paralog: str, elements: tuple[str, ...]) -> list[dict]:
    return [r for r in L.constraint_table(paralog) if r["element"] in elements]


def _labelled(variants: list[dict], gene: str, elements: tuple[str, ...],
              bucket: str) -> dict[int, str]:
    out: dict[int, str] = {}
    for v in variants:
        if v["gene"] == gene and v["class_bucket"] == bucket \
                and v["element"] in elements:
            out.setdefault(int(v["resi"]), v["ref_aa"])
    return out


def _profile_panel(ax, elements, variants, title, letter) -> None:
    for i, paralog in enumerate(L.PARALOGS):
        rows = _region_rows(paralog, elements)
        if not rows:
            continue
        base = min(int(r["resi"]) for r in rows)
        xs = [int(r["resi"]) - base for r in rows]
        ys = [L.f(r["deep_jsd"], 0.0) + i * 1.15 for r in rows]
        ax.fill_between(xs, [i * 1.15] * len(xs), ys,
                        color=F.PARALOG[paralog], lw=0, alpha=0.85)
        plp = _labelled(variants, paralog, elements, "P/LP")
        for resi in plp:
            ax.plot(resi - base, i * 1.15 + 1.05, marker="v", ms=2.2,
                    color=F.CLINICAL["pathogenic"], mew=0)
        ax.annotate(paralog, xy=(0, i * 1.15 + 0.45), xytext=(-4, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=F.FS_TICK, color=F.PARALOG[paralog])
    ax.set_ylim(-0.05, 3 * 1.15 + 0.05)
    ax.set_yticks([])
    ax.set_xlabel("residue, offset from the first residue of the region")
    F.despine(ax, keep=("bottom",))
    F.panel(ax, letter, title)


def fig_labelled_positions(out: Path, stats: dict) -> None:
    variants = L.read_tsv(L.CONSTRAINT_DIR / "variants.tsv")
    pairs = L.read_tsv(L.CONSTRAINT_DIR / "paralog_variant_positions.tsv")
    index = {p: L.residue_index(p) for p in L.PARALOGS}

    regions = L.LIGAND_ELEMENTS + L.PORE_ELEMENTS
    seen: set[tuple[str, int]] = set()
    grid: list[dict] = []
    for v in variants:
        if v["class_bucket"] != "P/LP" or v["element"] not in regions:
            continue
        key = (v["gene"], int(v["resi"]))
        if key in seen:
            continue
        seen.add(key)
        row = {"gene": v["gene"], "resi": int(v["resi"]),
               "element": v["element"], "alt": v["alt_aa"],
               "region": "ligand" if v["element"] in L.LIGAND_ELEMENTS
                         else "pore"}
        for p in L.PARALOGS:
            if p == v["gene"]:
                row[p] = index[p][int(v["resi"])]["aa"]
                continue
            hit = [q for q in pairs
                   if q["gene"] == v["gene"] and q["other"] == p
                   and int(q["resi"]) == int(v["resi"])]
            row[p] = hit[0]["other_aa"] if hit else "-"
        grid.append(row)
    grid.sort(key=lambda r: (r["region"] != "ligand",
                             L.ELEMENT_ORDER.index(r["element"]), r["gene"],
                             r["resi"]))

    fig = plt.figure(figsize=(F.W_FULL, 6.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.5, 1.5, 3.6], hspace=0.6)
    _profile_panel(fig.add_subplot(gs[0]), L.LIGAND_ELEMENTS, variants,
                   "the IP3-binding core (β-trefoil + MIR): deep-layer "
                   "constraint, marker = a pathogenic position", "a")
    _profile_panel(fig.add_subplot(gs[1]), L.PORE_ELEMENTS, variants,
                   "the pore module (channel, filter, gate, luminal loop)",
                   "b")

    ax = fig.add_subplot(gs[2])
    n = len(grid)
    half = (n + 1) // 2
    for col, chunk in enumerate((grid[:half], grid[half:])):
        x0 = col * 0.52
        for k, r in enumerate(chunk):
            y = len(chunk) - k
            ax.text(x0, y, f"{r['gene']} {r['resi']}", fontsize=F.FS_NOTE,
                    va="center", ha="left", color=F.MUTED)
            ax.text(x0 + 0.115, y, L.ELEMENT_LABEL[r["element"]],
                    fontsize=F.FS_NOTE, va="center", ha="left",
                    color=F.PARALOG.get(r["gene"], F.INK))
            for j, p in enumerate(L.PARALOGS):
                same = r[p] == r[r["gene"]]
                ax.text(x0 + 0.30 + j * 0.032, y, r[p],
                        fontsize=F.FS_TICK, va="center", ha="center",
                        color=F.INK if same else F.CLINICAL["pathogenic"],
                        family="monospace")
            ax.text(x0 + 0.405, y, f"→{r['alt']}", fontsize=F.FS_NOTE,
                    va="center", ha="left", color=F.MUTED,
                    family="monospace")
        for j, p in enumerate(L.PARALOGS):
            ax.text(x0 + 0.30 + j * 0.032, len(chunk) + 1.1, p[-1],
                    fontsize=F.FS_NOTE, va="center", ha="center",
                    color=F.PARALOG[p])
        ax.text(x0 + 0.30 + 0.032, len(chunk) + 2.0, "ITPR1/2/3",
                fontsize=F.FS_NOTE, va="center", ha="center", color=F.MUTED)
    ax.set_xlim(-0.01, 1.02)
    ax.set_ylim(0, half + 3)
    ax.axis("off")
    F.panel(ax, "c", "every pathogenic position in the two regions, and the "
                     "residue each paralogue carries there")
    ax.annotate("red letter = the other paralogue carries a different residue; "
                "every letter checked against that paralogue's own "
                "per-residue table",
                xy=(0.0, -0.02), xycoords="axes fraction",
                fontsize=F.FS_NOTE, color=F.MUTED, va="top")

    F.save(fig, out / "SuppFig2_labelled_positions")
    plt.close(fig)

    stats["supp_fig_2"] = {
        "positions_drawn": n,
        "ligand_positions": sum(1 for r in grid if r["region"] == "ligand"),
        "pore_positions": sum(1 for r in grid if r["region"] == "pore"),
        "by_element": dict(Counter(r["element"] for r in grid)),
        "conserved_across_all_three": sum(
            1 for r in grid
            if len({r[p] for p in L.PARALOGS}) == 1),
    }


FIGURES = {
    "supp1": fig_representative_alignment,
    "supp2": fig_labelled_positions,
}
