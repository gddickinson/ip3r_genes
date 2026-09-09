"""The four S13 figures, from committed S13 tables only (D13, D19).

  1. `recon_dated_backbone`  the dated species tree with every duplication the
                             reconciliation places drawn *on the branch it maps
                             to*, on a real time axis. The brief's dated
                             backbone. Ages come from
                             `species_tree_calibrations.tsv`, placements from
                             `duplication_placement.tsv`; nothing is positioned
                             by eye.
  2. `recon_matrix`          the topology x variant matrix — where the deepest
                             paralog duplication maps in every cell, with the
                             event counts. The figure exists because agreement
                             across the matrix is the robustness claim, and a
                             sentence saying "it never moves" is not checkable.
  3. `recon_losses`          the loss audit: what the 53 implied losses turn
                             into once each is asked of the S5 genome ledger.
  4. `recon_cyclostome`      the six tips the deep placement rests on — their
                             root-to-tip distances against all 57 (the
                             long-branch check) and S8's independent flank call
                             per locus against its own null.

The dated backbone is drawn on a **linear time axis with the calibration
spread shown as a band**, not as a cladogram with ages written on it: the
whole point of a placement result is which interval the event falls in, and a
cladogram cannot show an interval.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s13_figures.py [--only slug]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.patches import Patch                         # noqa: E402

import figstyle as fs                                        # noqa: E402
from s13_lib import RECON_DIR, read_tsv                      # noqa: E402
import s13_species_tree as ST                                # noqa: E402

FIGS = RECON_DIR / "figures"

VERDICT_COLOUR = {
    "sampled_in_gene_tree": fs.STATUS["found_annotated"],
    "sampling_artefact": fs.STATUS["found_no_annotation"],
    "paralog_unassignable": fs.STATUS["assembly_gap"],
    "undecidable": fs.STATUS["tblastn_trace_ambiguous"],
    "loss_with_remnant": fs.STATUS["tblastn_trace"],
    "corroborated_loss": fs.STATUS["absent"],
    "no_genome_in_manifest": "#efeee8",
}
VERDICT_LABEL = {
    "sampled_in_gene_tree": "in the gene tree",
    "sampling_artefact": "in the genome, unsampled",
    "paralog_unassignable": "loci present, paralog unassignable",
    "undecidable": "undecidable",
    "loss_with_remnant": "loss with a remnant",
    "corroborated_loss": "corroborated loss",
    "no_genome_in_manifest": "no genome in scope",
}


# ------------------------------------------------------------------ fig 1

def _layout(node):
    """Species-tree layout: name -> (age Ma, y), tips evenly spaced."""
    pos, y = {}, [0.0]

    def rec(n):
        if isinstance(n, str):
            pos[n] = (0.0, y[0])
            y[0] += 1.0
            return y[0] - 1.0
        name, kids = n
        ys = [rec(k) for k in kids]
        pos[name] = (ST.CALIBRATIONS[name][0], sum(ys) / len(ys))
        return pos[name][1]

    rec(node)
    return pos, y[0]


def fig_dated_backbone():
    cal = {r["node"]: r for r in read_tsv(RECON_DIR / "species_tree_calibrations.tsv")}
    place = [r for r in read_tsv(RECON_DIR / "duplication_placement.tsv")
             if r["topology"] == "ml"]
    pos, n_tips = _layout(ST.TOPOLOGY)

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 5.0))
    # branches
    def draw(n):
        if isinstance(n, str):
            return
        name, kids = n
        x0, y0 = pos[name]
        ys = []
        for k in kids:
            kname = k if isinstance(k, str) else k[0]
            x1, y1 = pos[kname]
            ax.plot([x0, x1], [y1, y1], color=fs.MUTED, lw=0.8, zorder=1)
            ys.append(y1)
            draw(k)
        ax.plot([x0, x0], [min(ys), max(ys)], color=fs.MUTED, lw=0.8, zorder=1)
    draw(ST.TOPOLOGY)

    for tip in ST.walk_tips(ST.TOPOLOGY):
        _, ty = pos[tip]
        ax.text(-6, ty, tip.replace("_", " "), fontsize=4.6, va="center",
                ha="left", color=fs.INK, style="italic")

    # calibration spread on the nodes a placement uses
    used = {r["species_node"] for r in place}
    for name in used:
        if name not in pos:
            continue
        x, ny = pos[name]
        row = cal[name]
        ax.plot([float(row["age_lo"]), float(row["age_hi"])], [ny, ny],
                color=fs.PARALOG["ITPR1"], lw=3.0, alpha=0.20,
                solid_capstyle="butt", zorder=2)

    # the duplications, on the branch each maps to. One marker per distinct
    # (variant, node, paralog set): the matrix repeats the same placement in
    # several cells and drawing each would make agreement look like weight.
    variant_colour = {"with_cyclostome": fs.PARALOG["ITPR2"],
                      "without_cyclostome": fs.PARALOG["ITPR3"],
                      "support_collapsed": fs.MUTED}
    variant_label = {"with_cyclostome": "with the cyclostome loci",
                     "without_cyclostome": "with them dropped",
                     "support_collapsed": "with unsupported nodes collapsed"}
    seen_marks, seen, drawn = set(), {}, set()
    for r in place:
        key = (r["variant"], r["species_node"], r["paralog_groups"])
        if key in seen_marks:
            continue
        seen_marks.add(key)
        name = r["species_node"]
        if name not in pos:
            continue
        x, ny = pos[name]
        off = seen.get(name, 0)
        seen[name] = off + 1
        drawn.add(r["variant"])
        ax.plot([x], [ny + 0.55 + 0.55 * off], marker="s", ms=4.5,
                color=variant_colour[r["variant"]], mec="white", mew=0.5,
                zorder=4)
        ax.annotate(r["paralog_groups"].replace("+unplaced", "+cyclostome")
                    .replace("ITPR", ""),
                    xy=(x, ny + 0.55 + 0.55 * off), xytext=(4, 0),
                    textcoords="offset points", fontsize=4.6,
                    va="center", color=fs.INK)

    ax.set_xlim(620, -95)
    ax.set_ylim(-1.2, n_tips + 2.4)
    ax.set_yticks([])
    ax.set_xlabel("million years before present")
    ax.set_xticks([600, 500, 400, 300, 200, 100, 0])
    fs.despine(ax, keep=("bottom",))
    fs.hgrid(ax, axis="x")
    ax.legend(handles=[Patch(facecolor=variant_colour[v], label=variant_label[v])
                       for v in variant_colour if v in drawn]
              + [Patch(facecolor=fs.PARALOG["ITPR1"], alpha=0.20,
                       label="published spread of the node's age")],
              loc="upper left", frameon=False, fontsize=5.4,
              title="paralog duplication placed", title_fontsize=5.4)
    fs.panel(ax, "", "Where the reconciliation puts each duplication, on the "
                     "dated species tree")
    fig.tight_layout()
    return fs.save(fig, FIGS / "recon_dated_backbone")


# ------------------------------------------------------------------ fig 2

def fig_matrix():
    rows = [r for r in read_tsv(RECON_DIR / "reconciliation_summary.tsv")]
    topos = list(dict.fromkeys(r["topology"] for r in rows))
    variants = list(dict.fromkeys(r["variant"] for r in rows))
    cell = {(r["topology"], r["variant"]): r for r in rows}
    node_colour = {"Vertebrata": fs.PARALOG["ITPR2"],
                   "Gnathostomata": fs.PARALOG["ITPR3"]}

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.5),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    for i, t in enumerate(topos):
        for j, v in enumerate(variants):
            r = cell.get((t, v))
            if r is None or r.get("note"):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=5.2,
                        color=fs.MUTED)
                continue
            node = r["deepest_paralog_dup_node"]
            ax.add_patch(plt.Rectangle((j - .45, i - .42), .9, .84,
                                       facecolor=node_colour.get(node, "#ddd"),
                                       alpha=.30, lw=0))
            ax.text(j, i + .10, node, ha="center", va="center", fontsize=5.4,
                    color=fs.INK)
            ax.text(j, i - .19,
                    f"{r['n_duplications']} dup / {r['implied_losses']} loss",
                    ha="center", va="center", fontsize=4.8, color=fs.MUTED)
    ax.set_xticks(range(len(variants)))
    ax.set_xticklabels([v.replace("_", "\n") for v in variants], fontsize=5.4)
    ax.set_yticks(range(len(topos)))
    ax.set_yticklabels(topos, fontsize=5.6)
    ax.set_xlim(-.6, len(variants) - .4)
    ax.set_ylim(-.6, len(topos) - .4)
    ax.invert_yaxis()
    fs.despine(ax, keep=())
    fs.panel(ax, "a", "Deepest paralog duplication, per cell")

    # b — the rooting check
    rk = read_tsv(RECON_DIR / "rooting_check.tsv")
    out = [r for r in rk if r["edge_index"] == "-1"][0]
    others = [r for r in rk if r["edge_index"] != "-1"]
    tot = [int(r["total_events"]) for r in others]
    ax2.hist(tot, bins=range(min(tot), max(tot) + 2), color=fs.STATUS["fragment"],
             edgecolor="white", linewidth=.4)
    ax2.axvline(int(out["total_events"]), color=fs.PARALOG["ITPR1"], lw=1.2)
    ax2.annotate("outgroup rooting", xy=(int(out["total_events"]), 0),
                 xytext=(3, 26), textcoords="offset points", fontsize=5.0,
                 color=fs.PARALOG["ITPR1"], rotation=90)
    ax2.axvline(min(tot), color=fs.PARALOG["ITPR2"], lw=1.2, linestyle=(0, (3, 2)))
    ax2.annotate("minimum-event rooting", xy=(min(tot), 0), xytext=(3, 26),
                 textcoords="offset points", fontsize=5.0,
                 color=fs.PARALOG["ITPR2"], rotation=90)
    ax2.set_xlabel("duplications + implied losses")
    ax2.set_ylabel("rootings")
    fs.despine(ax2)
    fs.hgrid(ax2)
    fs.panel(ax2, "b", "Every rooting of the vertebrate subtree")
    fig.tight_layout()
    return fs.save(fig, FIGS / "recon_matrix")


# ------------------------------------------------------------------ fig 3

def fig_losses():
    rows = read_tsv(RECON_DIR / "loss_verdicts.tsv")
    paralogs = [p for p in ("ITPR1", "ITPR2", "ITPR3")]
    order = [v for v in VERDICT_LABEL if any(r["verdict"] == v for r in rows)]
    counts = {(r["paralog"], r["verdict"]): int(r["n_species"]) for r in rows}

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.3),
                                  gridspec_kw={"width_ratios": [1.5, 1]})
    left = [0.0] * len(paralogs)
    for v in order:
        vals = [counts.get((p, v), 0) for p in paralogs]
        ax.barh(range(len(paralogs)), vals, left=left, height=.6,
                color=VERDICT_COLOUR[v], edgecolor="white", linewidth=.4,
                label=VERDICT_LABEL[v])
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(range(len(paralogs)))
    ax.set_yticklabels(paralogs)
    ax.invert_yaxis()
    ax.set_xlabel("species in the sampled species tree")
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2,
              frameon=False, fontsize=5.0)
    fs.panel(ax, "a", "What each implied loss turns out to be")

    # b — implied losses per matrix cell, so the raw number is visible beside
    #     the number that survives the audit
    summ = [r for r in read_tsv(RECON_DIR / "reconciliation_summary.tsv")
            if not r.get("note")]
    labels = [f"{r['topology']}\n{r['variant'].replace('_', ' ')}" for r in summ]
    vals = [int(r["implied_losses"]) for r in summ]
    real = counts.get(("all", "corroborated_loss"), 0)
    ax2.bar(range(len(vals)), vals, color=fs.STATUS["fragment"],
            edgecolor="white", linewidth=.4)
    ax2.axhline(real, color=fs.STATUS["absent"], lw=1.2)
    ax2.annotate(f"corroborated by the genome sweep: {real}",
                 xy=(0, real), xytext=(0, 4), textcoords="offset points",
                 fontsize=5.0, color=fs.STATUS["absent"])
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, fontsize=3.8, rotation=90)
    ax2.set_ylabel("implied losses")
    fs.despine(ax2)
    fs.hgrid(ax2)
    fs.panel(ax2, "b", "Implied against corroborated")
    fig.tight_layout()
    return fs.save(fig, FIGS / "recon_losses")


# ------------------------------------------------------------------ fig 4

def fig_cyclostome():
    bl = read_tsv(RECON_DIR / "branch_lengths.tsv")
    cy = read_tsv(RECON_DIR / "cyclostome_loci.tsv")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.4))
    vals = [float(r["root_to_tip"]) for r in bl]
    ax.hist(vals, bins=18, color=fs.STATUS["fragment"], edgecolor="white",
            linewidth=.4)
    for r in bl:
        if r["is_cyclostome"]:
            ax.axvline(float(r["root_to_tip"]), color=fs.GROUP["vertebrate_basal"],
                       lw=1.0)
    ax.set_xlabel("root-to-tip distance (substitutions/site)")
    ax.set_ylabel("vertebrate tips")
    ax.legend(handles=[Patch(facecolor=fs.GROUP["vertebrate_basal"],
                             label="the six cyclostome loci")],
              loc="upper right", frameon=False, fontsize=5.2)
    fs.despine(ax)
    fs.hgrid(ax)
    fs.panel(ax, "a", "Are the cyclostome tips long branches?")

    # b — S8's independent call per locus, against its own null
    labs, cols, texts = [], [], []
    for r in cy:
        sp = r["species"].split()[0][:4]
        labs.append(f"{sp} {r['s8_annot_gene'][-6:] or '?'}")
        call = r["s8_paralog_call"]
        cols.append(fs.PARALOG.get(call, fs.STATUS["assembly_gap"]))
        texts.append(f"{call or 'no call'} ({r['s8_null_verdict'] or 'n/a'})")
    ax2.barh(range(len(labs)), [1] * len(labs), height=.6, color=cols,
             edgecolor="white", linewidth=.4)
    for i, t in enumerate(texts):
        ax2.text(0.03, i, t, va="center", fontsize=5.0, color="white"
                 if cols[i] != fs.STATUS["assembly_gap"] else fs.INK)
    for i, r in enumerate(cy):
        ax2.text(1.03, i, f"pair {r['pair_support'] or '-'}", va="center",
                 fontsize=4.8, color=fs.MUTED)
    ax2.text(1.03, len(cy) - 0.35, "SH-aLRT/UFBoot of the tree's\nown "
             "orthology pair for this locus",
             va="top", fontsize=4.6, color=fs.MUTED)
    ax2.set_yticks(range(len(labs)))
    ax2.set_yticklabels(labs, fontsize=5.0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 1.6)
    ax2.set_xticks([])
    fs.despine(ax2, keep=())
    # The title is computed, not asserted: two of the six loci made no call
    # at all, and "every one inside its null" reported them as agreeing
    # (found by the S24 figure audit).
    n_called = sum(1 for r in cy if r["s8_null_verdict"] == "within_null")
    n_nocall = sum(1 for r in cy if not r["s8_paralog_call"]
                   or r["s8_paralog_call"] == "no_call")
    fs.panel(ax2, "b", f"S8's flank call per locus — {n_called} of "
                       f"{len(cy)} inside its null, {n_nocall} no call")
    fig.tight_layout()
    return fs.save(fig, FIGS / "recon_cyclostome")


FIGURES = {
    "recon_dated_backbone": fig_dated_backbone,
    "recon_matrix": fig_matrix,
    "recon_losses": fig_losses,
    "recon_cyclostome": fig_cyclostome,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for k in FIGURES:
            print(k)
        return
    fs.use()
    FIGS.mkdir(parents=True, exist_ok=True)
    for slug in (args.only or list(FIGURES)):
        paths = FIGURES[slug]()
        print(f"[s13-fig] {slug} -> {', '.join(p.name for p in paths)}")


if __name__ == "__main__":
    main()
