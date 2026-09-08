"""The four S16 figures, from committed S16 tables only (D13, D19).

Four decisions worth naming:

* **The copy-number figure is grouped by 3R status, not by taxonomy.** The
  result is a contrast between three groups defined by which whole-genome
  duplications a lineage has been through, and a per-class bar would bury it
  inside Actinopteri.
* **The RyR control is drawn beside the ITPRs everywhere it exists**, in the
  accent violet `figstyle.GROUP["RYR"]` reserves for it. A panel showing only
  the family under test cannot show that the instrument works.
* **The 2R panel draws its null**, as a line across the bars, rather than
  quoting it in a caption: 80 % of genomes carrying a paralog link means
  nothing until 2.6 % is visible on the same axis.
* **The DCS panel plots the two copies against each other**, so the
  disjointness is a geometric fact on the figure — points off both axes
  would mean both copies carry the same ancestral symbols, which is what a
  tandem duplication looks like.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402

import figstyle as FS                                          # noqa: E402
import s16_lib as L                                            # noqa: E402
import s16_tables as TB                                        # noqa: E402
import s16_teleost as T                                        # noqa: E402

FIG_DIR = L.OUT_DIR / "figures"
GROUP_LABEL = {"pre_3R_outgroup": "pre-3R ray-fins\n(bichir, gar, bowfin)",
               "teleost": "teleosts\n(3R)",
               "extra_wgd": "extra WGD\n(sturgeon, salmon)"}
GROUP_ORDER = ["pre_3R_outgroup", "teleost", "extra_wgd"]


def _rows(name: str) -> list[dict]:
    p = L.OUT_DIR / name
    return L.read_tsv(p) if p.exists() else []


# --------------------------------------------------- 1. copy-number landscape

def fig_copy_number(out: Path) -> None:
    tel = _rows("teleost_copies.tsv")
    clade = _rows("copy_number_by_clade.tsv")
    if not tel or not clade:
        return
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))

    ax = axes[0]
    cells = list(L.PARALOGS) + [L.CONTROL_CELL]
    width = 0.2
    for k, cell in enumerate(cells):
        col = (FS.PARALOG.get(cell) or FS.GROUP["RYR"])
        xs, ys = [], []
        for i, grp in enumerate(GROUP_ORDER):
            sub = [r for r in tel
                   if r["group"] == grp and r["contig_spans_gene"] == "1"]
            if not sub:
                continue
            xs.append(i + (k - 1.5) * width)
            ys.append(sum(int(r[f"{cell}_copies"]) for r in sub) / len(sub))
        ax.bar(xs, ys, width=width, color=col, label=cell,
               edgecolor="none")
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER], fontsize=6)
    ax.set_ylabel("mean gene copies per genome")
    FS.hgrid(ax)
    FS.despine(ax)
    ax.legend(fontsize=6, frameon=False, ncol=2, loc="upper left")
    FS.panel(ax, "a", "3R doubled ITPR1 alone — and all three RyRs")

    ax = axes[1]
    order = [r for r in clade if r["cell"] == "ITPR1"
             and int(r["n_above_bar"]) >= 4]
    order.sort(key=lambda r: -float(r["frac_multi_above_bar"]))
    labels = [r["vclass"] for r in order]
    for k, cell in enumerate(cells):
        col = (FS.PARALOG.get(cell) or FS.GROUP["RYR"])
        vals = []
        for lab in labels:
            m = [r for r in clade if r["vclass"] == lab and r["cell"] == cell]
            vals.append(float(m[0]["frac_multi_above_bar"]) if m else 0.0)
        ax.bar([i + (k - 1.5) * width for i in range(len(labels))], vals,
               width=width, color=col, edgecolor="none")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=6)
    ax.set_ylabel("genomes with >1 copy")
    ax.set_ylim(0, 1.05)
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "b", "only the ray-fins carry a second ITPR")
    fig.tight_layout()
    FS.save(fig, out / "s16_copy_number")
    plt.close(fig)


# ------------------------------------------------------------- 2. the 2R test

def fig_paralogon(out: Path) -> None:
    summary = _rows("paralogon_species_summary.tsv")
    tests = _rows("paralogon_test.tsv")
    if not summary or not tests:
        return
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))

    ax = axes[0]
    pairs = [r["pair"] for r in summary]
    frac = [float(r["frac_with_link"]) for r in summary]
    frac2 = [float(r["frac_with_2R_link"]) for r in summary]
    null = float(summary[0]["null_frac"])
    x = range(len(pairs))
    ax.bar([i - 0.19 for i in x], frac, width=0.36, color=FS.BLUES[2],
           label="any paralog link", edgecolor="none")
    ax.bar([i + 0.19 for i in x], frac2, width=0.36, color=FS.BLUES[5],
           label="2R-dated link", edgecolor="none")
    ax.axhline(null, color=FS.CLINICAL["pathogenic"], linewidth=1.0,
               linestyle="--")
    ax.annotate(f"matched random windows, {null:.1%}",
                xy=(len(pairs) - 0.55, null), xytext=(0, 4),
                textcoords="offset points", fontsize=6,
                color=FS.CLINICAL["pathogenic"], ha="right")
    ax.set_xticks(list(x))
    ax.set_xticklabels([p.replace("_vs_", "\nvs ") for p in pairs], fontsize=6)
    ax.set_ylabel("genomes with a link")
    ax.set_ylim(0, 1.0)
    FS.hgrid(ax)
    FS.despine(ax)
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    FS.panel(ax, "a", "replicated in 309 genomes, against its null")

    ax = axes[1]
    prim = [r for r in tests if int(r["window_n"]) == TB.PRIMARY_WINDOW
            and not r["pair_class"].startswith("ITPR_vs_RYR")]
    keyed = {(r["level_set"], r["group_a"], r["group_b"]): r for r in prim}
    labs = sorted({(r["group_a"], r["group_b"], r["pair_class"])
                   for r in prim}, key=lambda t: (t[2], t[0], t[1]))
    y = range(len(labs))
    for off, level_set, shade, lab in (
            (-0.19, "all_levels", FS.BLUES[2], "any paralog link"),
            (0.19, TB.PRIMARY_LEVEL_SET, FS.BLUES[5], "2R-dated link")):
        vals = [int(keyed[(level_set, a, b)]["observed_links"])
                if (level_set, a, b) in keyed else 0 for a, b, _ in labs]
        ax.barh([i + off for i in y], vals, height=0.36, color=shade,
                edgecolor="none", label=lab)
    nul = [float(keyed[("all_levels", a, b)]["null_mean"])
           if ("all_levels", a, b) in keyed else 0.0 for a, b, _ in labs]
    ax.plot(nul, list(y), marker="D", linestyle="none", markersize=3.0,
            color=FS.CLINICAL["pathogenic"], label="null (mean, undated)")
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{a}–{b}" + (" (control)" if c.startswith("RYR")
                                      else "") for a, b, c in labs],
                       fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("paralog links in the human windows")
    ax.set_xticks(range(0, 5))
    FS.hgrid(ax, axis="x")
    FS.despine(ax)
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    FS.panel(ax, "b", "the human windows, links and null")
    fig.tight_layout()
    FS.save(fig, out / "s16_paralogon")
    plt.close(fig)


# ------------------------------------------------------- 3. double-conserved

def fig_dcs(out: Path) -> None:
    dcs = _rows("dcs_test.tsv")
    if not dcs:
        return
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))

    ax = axes[0]
    xa = [int(r["shared_tetrapod_a"]) for r in dcs]
    xb = [int(r["shared_tetrapod_b"]) for r in dcs]
    disj = [r["tetrapod_partition_disjoint"] == "1" for r in dcs]
    seen: collections.Counter = collections.Counter(zip(xa, xb))
    for (a, b), n in seen.items():
        ok = any(d for d, x, y in zip(disj, xa, xb) if (x, y) == (a, b))
        ax.plot([a], [b], marker="o", linestyle="none",
                markersize=2.4 + 1.5 * (n ** 0.5),
                color=FS.PARALOG["ITPR1"] if ok else FS.QUALITY["fragmentary"],
                alpha=0.85)
    lim = max(xa + xb) + 1
    ax.plot([0, lim], [0, 0], color=FS.GRID, linewidth=0.6)
    ax.plot([0, 0], [0, lim], color=FS.GRID, linewidth=0.6)
    ax.set_xlim(-0.4, lim)
    ax.set_ylim(-0.4, lim)
    ax.set_xlabel("ancestral-block symbols kept by copy A")
    ax.set_ylabel("kept by copy B")
    FS.hgrid(ax)
    FS.despine(ax)
    FS.panel(ax, "a", "both copies keep part of one block")

    ax = axes[1]
    cats = ["both keep\nancestral symbols", "and the two sets\nare disjoint"]
    tet = [sum(1 for r in dcs if r["both_share_tetrapod"] == "1"),
           sum(1 for r in dcs if r["tetrapod_partition_disjoint"] == "1")]
    fish = [sum(1 for r in dcs if r["both_share_fish"] == "1"),
            sum(1 for r in dcs if r["fish_partition_disjoint"] == "1")]
    x = range(len(cats))
    ax.bar([i - 0.19 for i in x], tet, width=0.36, color=FS.BLUES[2],
           label="vs the tetrapod block", edgecolor="none")
    ax.bar([i + 0.19 for i in x], fish, width=0.36, color=FS.BLUES[5],
           label="vs the pre-3R ray-finned block", edgecolor="none")
    ax.axhline(len(dcs), color=FS.MUTED, linewidth=0.8, linestyle=":")
    ax.annotate(f"all {len(dcs)} two-copy genomes", xy=(1.4, len(dcs)),
                xytext=(0, 3), textcoords="offset points", fontsize=6,
                color=FS.MUTED, ha="right")
    ax.set_xticks(list(x))
    ax.set_xticklabels(cats, fontsize=6)
    ax.set_ylabel("genomes")
    ax.set_ylim(0, len(dcs) * 1.18)
    FS.hgrid(ax)
    FS.despine(ax)
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    FS.panel(ax, "b", "the 3R signature, two references")
    fig.tight_layout()
    FS.save(fig, out / "s16_dcs")
    plt.close(fig)


# ------------------------------------------- 4. blocks, dispersion, sensitivity

def fig_blocks(out: Path) -> None:
    r3 = {r["metric"]: r["value"] for r in _rows("r3_summary.tsv")}
    assign = _rows("block_assignments.tsv")
    sens = _rows("copy_sensitivity.tsv")
    if not assign or not sens:
        return
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))

    ax = axes[0]
    by_anchor: dict[str, list[float]] = collections.defaultdict(list)
    for r in assign:
        by_anchor[r["anchor_organism"]].append(float(r["margin"]))
    names = sorted(by_anchor, key=lambda k: -len(by_anchor[k]))
    ax.boxplot([by_anchor[n] for n in names], vert=False, widths=0.55,
               patch_artist=True,
               boxprops={"facecolor": FS.BLUES[1], "edgecolor": FS.INK,
                         "linewidth": 0.6},
               medianprops={"color": FS.INK, "linewidth": 0.9},
               whiskerprops={"color": FS.INK, "linewidth": 0.6},
               capprops={"color": FS.INK, "linewidth": 0.6},
               flierprops={"marker": "o", "markersize": 1.6,
                           "markerfacecolor": FS.MUTED,
                           "markeredgecolor": "none"})
    ax.set_yticklabels([f"$\\it{{{n.replace(' ', chr(92) + ' ')}}}$"
                        for n in names], fontsize=6)
    ax.set_xlabel("assignment margin (Jaccard difference)")
    frac = r3.get("fraction_agreeing", "")
    FS.hgrid(ax, axis="x")
    FS.despine(ax)
    FS.panel(ax, "a", f"every anchor agrees ({frac} of calls)")

    ax = axes[1]
    for cell in list(L.PARALOGS) + [L.CONTROL_CELL]:
        rows = [r for r in sens
                if r["cell"] == cell and r["scope"] == "above_contiguity_bar"]
        rows.sort(key=lambda r: float(r["cov_bar"]))
        ax.plot([float(r["cov_bar"]) for r in rows],
                [int(r["genomes_with_2plus"]) for r in rows],
                marker="o", markersize=2.6, linewidth=1.0,
                color=FS.PARALOG.get(cell) or FS.GROUP["RYR"], label=cell)
    ax.axvline(L.COV_FULL, color=FS.MUTED, linewidth=0.8, linestyle=":")
    ax.annotate("operating point", xy=(L.COV_FULL, 0), xytext=(3, 6),
                textcoords="offset points", fontsize=6, color=FS.MUTED)
    ax.set_xlabel("coverage bar a locus must clear to count as a copy")
    ax.set_ylabel("genomes with >1 copy")
    FS.hgrid(ax)
    FS.despine(ax)
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    FS.panel(ax, "b", "the counts do not move with the bar")
    fig.tight_layout()
    FS.save(fig, out / "s16_blocks")
    plt.close(fig)


def run(out_dir: Path | None = None, log=print) -> None:
    FS.use()
    out = out_dir or FIG_DIR
    out.mkdir(parents=True, exist_ok=True)
    for fn in (fig_copy_number, fig_paralogon, fig_dcs, fig_blocks):
        fn(out)
        log(f"[figures] {fn.__name__}")


def main() -> None:
    run()


if __name__ == "__main__":
    main()
