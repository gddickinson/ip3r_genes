"""Figures for §6.1, §7.1 and §7.5 — genes, family separation, taxonomic range.

`gene_architecture`  the 6.5-fold span asymmetry that became open question Q7
`family_separation`  IP3R and RyR separate by a measured margin, not by name
`taxonomic_range`    the metazoan family that holds plant and fungal records

All three render from tables committed by S0 and S1.
"""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import figstyle
import s0_fig_lib as lib

PARALOGS = ["ITPR1", "ITPR2", "ITPR3"]


def gene_architecture():
    figstyle.use()
    rows = {r["symbol"]: r for r in lib.load_tsv(lib.S0 / "gene_structure.tsv")}

    fig, axes = plt.subplots(1, 2, figsize=(lib.W_REVIEW, 2.35),
                             gridspec_kw={"width_ratios": [1.0, 0.58],
                                          "wspace": 0.34})
    ax, ax2 = axes

    # ---- a: genomic span, with everything that is *not* different moved
    # into the row label, so nothing has to fit inside a 76 kb bar
    y = np.arange(len(PARALOGS))[::-1]
    spans = [int(rows[s]["genomic_span_bp"]) for s in PARALOGS]
    for yi, sym, span in zip(y, PARALOGS, spans):
        ax.barh(yi, span / 1000, height=0.52,
                color=figstyle.PARALOG[sym], zorder=3)
        ax.text(span / 1000 + 12, yi, f"{span/1000:,.0f} kb", va="center",
                ha="left", fontsize=figstyle.FS_NOTE, color=figstyle.INK)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"$\\it{{{s}}}$\n{rows[s]['n_exons_canonical']} exons · "
         f"{int(rows[s]['canonical_protein_aa']):,} aa" for s in PARALOGS],
        fontsize=figstyle.FS_NOTE, linespacing=1.4)
    ax.set_xlim(0, 700)
    ax.set_xticks(range(0, 501, 100))
    ax.set_xlabel("genomic span (kb)")
    ax.tick_params(axis="y", length=0)
    figstyle.hgrid(ax, axis="x")
    ratio = max(spans) / min(spans)
    ax.annotate("", xy=(646, 2.0), xytext=(646, 0.0),
                arrowprops=dict(arrowstyle="<->", lw=0.7,
                                color=figstyle.MUTED))
    ax.text(638, 1.0, f"{ratio:.1f}×", fontsize=figstyle.FS_LABEL,
            color=figstyle.INK, ha="right", va="center", fontweight="bold")
    figstyle.panel(ax, "a", "same protein, same exon count, different gene")

    # ---- b: the asymmetry as one number per gene
    dens = [int(rows[s]["genomic_span_bp"]) /
            int(rows[s]["canonical_protein_aa"]) for s in PARALOGS]
    ax2.bar(range(3), dens, width=0.56,
            color=[figstyle.PARALOG[s] for s in PARALOGS], zorder=3)
    for i, d in enumerate(dens):
        ax2.text(i, d + 4, f"{d:.0f}", ha="center",
                 fontsize=figstyle.FS_NOTE, color=figstyle.INK)
    ax2.set_xticks(range(3))
    ax2.set_xticklabels([f"$\\it{{{s}}}$" for s in PARALOGS],
                        fontsize=figstyle.FS_TICK)
    ax2.set_ylabel("bp of locus per\nresidue of protein")
    ax2.set_ylim(0, max(dens) * 1.22)
    figstyle.hgrid(ax2)
    figstyle.panel(ax2, "b", "packing")
    lib.provenance(fig, "measured", "Ensembl release 15.12")
    return lib.save(fig, "gene_architecture")


# ------------------------------------------------------------- separation

TRUTH_STYLE = {"ITPR": ("#2a78d6", "IP$_3$R (true positives)"),
               "RYR": ("#4a3aa7", "ryanodine receptors"),
               "other": ("#a9a79e", "other decoys")}


def family_separation():
    figstyle.use()
    bait = lib.load_tsv(lib.S1 / "bait_margin.tsv")
    zf = lib.load_tsv(lib.S0 / "zebrafish_itpr.tsv")

    fig = plt.figure(figsize=(lib.W_REVIEW, 2.35))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.86, 0.80], wspace=0.44)
    ax, ax2, ax3 = (fig.add_subplot(gs[i]) for i in range(3))

    # ---- a: the two identities, against each other
    for truth, (colour, lab) in TRUTH_STYLE.items():
        sel = [r for r in bait if r["truth"] == truth]
        ax.scatter([float(r["id_to_itpr_bait_cov"]) for r in sel],
                   [float(r["id_to_ryr_bait_cov"]) for r in sel],
                   s=11, color=colour, label=lab, zorder=3,
                   edgecolor="#ffffff", linewidth=0.3)
    lim = (0.10, 1.02)
    ax.plot(lim, lim, color=figstyle.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=2)
    ax.text(0.74, 0.60, "equally close\nto both baits", rotation=39,
            fontsize=figstyle.FS_NOTE - 0.6, color=figstyle.MUTED,
            ha="center", va="center")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("identity to the IP$_3$R bait")
    ax.set_ylabel("identity to the RyR bait")
    ax.set_aspect("equal")
    figstyle.hgrid(ax, axis="both")
    figstyle.panel(ax, "a", "every control, both baits")

    # ---- b: the margin, and the band inside which it is not a call
    ax2.axhspan(-0.10, 0.10, color="#f2f1ec", zorder=0)
    ax2.text(-0.48, 0.13, "no-call band (D7)",
             fontsize=figstyle.FS_NOTE - 0.6, color=figstyle.MUTED,
             va="bottom", ha="left")
    rng = np.random.default_rng(4)
    for xi, truth in enumerate(("ITPR", "RYR", "other")):
        colour, _ = TRUTH_STYLE[truth]
        vals = [float(r["margin_cov"]) for r in bait if r["truth"] == truth]
        ax2.scatter(xi + rng.uniform(-0.17, 0.17, len(vals)), vals, s=10,
                    color=colour, zorder=3, edgecolor="#ffffff",
                    linewidth=0.3)
    ax2.axhline(0, color=figstyle.MUTED, lw=0.6)
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(["IP$_3$R", "RyR", "other"],
                        fontsize=figstyle.FS_TICK)
    ax2.set_xlim(-0.55, 2.55)
    ax2.set_ylabel("margin  (IP$_3$R − RyR identity)")
    figstyle.hgrid(ax2)
    figstyle.panel(ax2, "b", "no overlap")

    # ---- c: what one domain query in one genome actually returns
    def klass(r):
        name = r["protein_name"].lower()
        if r["gene"].lower().startswith("itpr"):
            return "IP$_3$R"
        if "ryanodine" in name or r["gene"].lower().startswith("ryr"):
            return ("unnamed RyR locus"
                    if r["gene"].lower().startswith("loc") else "RyR")
        return "other"
    order = ["IP$_3$R", "RyR", "unnamed RyR locus", "other"]
    colour = {"IP$_3$R": "#2a78d6", "RyR": "#4a3aa7",
              "unnamed RyR locus": "#8d84d3", "other": "#a9a79e"}
    counts = {k: 0 for k in order}
    for r in zf:
        counts[klass(r)] += 1
    bottom = 0
    for k in order:
        if not counts[k]:
            continue
        ax3.bar(0, counts[k], bottom=bottom, width=0.62, color=colour[k],
                zorder=3)
        ax3.text(0.40, bottom + counts[k] / 2,
                 f"{k} — {counts[k]}", va="center", ha="left",
                 fontsize=figstyle.FS_NOTE - 0.4, color=figstyle.INK)
        bottom += counts[k]
    ax3.set_xlim(-0.45, 2.3)
    ax3.set_ylim(0, bottom * 1.04)
    ax3.set_xticks([])
    ax3.set_ylabel("protein records")
    figstyle.hgrid(ax3)
    figstyle.panel(ax3, "c", "PF08709 in zebrafish")
    fig.legend(handles=[Patch(facecolor=c, label=l)
                        for c, l in TRUTH_STYLE.values()],
               loc="lower center", ncol=3, bbox_to_anchor=(0.36, -0.10),
               fontsize=figstyle.FS_NOTE)
    lib.provenance(fig, "measured", "S1 control benchmark; InterPro")
    return lib.save(fig, "family_separation")


# ---------------------------------------------------------------- range

def taxonomic_range():
    figstyle.use()
    rows = lib.load_tsv(lib.S0 / "interpro_taxonomy_counts.tsv")
    doms = lib.load_tsv(lib.DATA / "domain_coords.tsv")
    kingdoms = ["Metazoa", "SAR (stramenopiles/alveolates/rhizaria)",
                "Discoba (incl. kinetoplastids)", "Fungi", "Viridiplantae",
                "Amoebozoa", "Bacteria", "Archaea"]
    by = {r["taxon"]: r for r in rows}

    fig, axes = plt.subplots(1, 2, figsize=(lib.W_REVIEW, 2.45),
                             gridspec_kw={"width_ratios": [1.0, 0.78],
                                          "wspace": 0.30})
    ax, ax2 = axes

    # ---- a: where the defining signature is found
    flagged = {"Fungi", "Viridiplantae"}
    y = np.arange(len(kingdoms))[::-1]
    for yi, tax in zip(y, kingdoms):
        n = int(by[tax]["n_PF08709"])
        if n:
            ax.barh(yi, n, height=0.58,
                    color="#b3261e" if tax in flagged else "#2a78d6", zorder=3)
        ax.text(max(n, 0.55) * 1.3, yi, f"{n:,}", va="center", ha="left",
                fontsize=figstyle.FS_NOTE,
                color="#b3261e" if n == 0 else figstyle.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([t.split(" (")[0] for t in kingdoms],
                       fontsize=figstyle.FS_TICK)
    ax.set_xscale("log")
    ax.set_xlim(0.5, 90_000)
    ax.set_xlabel("proteins carrying PF08709 (log scale)")
    ax.tick_params(axis="y", length=0)
    figstyle.hgrid(ax, axis="x")
    at, sc = (int(by["Arabidopsis thaliana"]["n_PF08709"]),
              int(by["Saccharomyces cerevisiae"]["n_PF08709"]))
    ax.text(70_000, 5.45, f"yet $\\it{{A.\\ thaliana}}$ carries {at}\n"
                          f"and $\\it{{S.\\ cerevisiae}}$ {sc}",
            fontsize=figstyle.FS_NOTE, color="#b3261e", ha="right",
            va="center", linespacing=1.35)
    ax.legend(handles=[Patch(facecolor="#b3261e",
                             label="lineages textbooks say lack the family")],
              loc="lower right", fontsize=figstyle.FS_NOTE,
              bbox_to_anchor=(1.0, 0.02))
    figstyle.panel(ax, "a", "the family is overwhelmingly metazoan")

    # ---- b: and the signature does not even find every real receptor
    sigs = ["PF08709", "PF02815", "PF01365", "PF08454", "PF00520"]
    prots = [("Q14643", "IP$_3$R1\nhuman"), ("P29993", "Itpr\nfly"),
             ("Q9Y0A1", "itr-1\nworm"), ("Q9NA13", "iplA\n$\\it{Dicty.}$")]
    have = {(d["accession"], d["pfam"]) for d in doms}
    for xi, (acc, _) in enumerate(prots):
        for yi, pf in enumerate(sigs):
            yy = len(sigs) - 1 - yi
            present = (acc, pf) in have
            ax2.add_patch(plt.Rectangle(
                (xi - 0.40, yy - 0.38), 0.80, 0.76,
                facecolor=lib.DOMAIN_COLOUR[pf] if present else "#ffffff",
                edgecolor="none" if present else "#b3261e",
                linewidth=0 if present else 0.8, zorder=3))
            if not present:
                ax2.plot([xi], [yy], marker="x", ms=3.6, mew=0.9,
                         color="#b3261e", zorder=4)
    ax2.set_xticks(range(len(prots)))
    ax2.set_xticklabels([p[1] for p in prots], fontsize=figstyle.FS_NOTE,
                        linespacing=1.35)
    ax2.set_yticks(range(len(sigs)))
    ax2.set_yticklabels([f"{lib.DOMAIN_SHORT[p]}\n{p}"
                         for p in reversed(sigs)],
                        fontsize=figstyle.FS_NOTE - 0.6, linespacing=1.25)
    ax2.set_xlim(-0.6, len(prots) - 0.4)
    ax2.set_ylim(-0.6, len(sigs) - 0.4)
    ax2.tick_params(length=0)
    figstyle.despine(ax2, keep=())
    figstyle.panel(ax2, "b", "signatures carried, per receptor")
    lib.provenance(fig, "measured", "InterPro, 2026-08-18")
    return lib.save(fig, "taxonomic_range")
