"""The four S10 figures, from committed S10 tables only (D13, D19).

The exon-track figure is the one the brief asks for and the one that has to be
right. Two decisions in it:

* **The two cases are drawn on their own coordinate axes, not a shared one.**
  The loci are 52 kb and 68 kb, on different chromosomes of different species;
  a shared axis would say nothing and would compress one of them.
* **Exons are drawn at true genomic width and never widened to be visible.** A
  56-exon gene over 52 kb has exons averaging 143 bp — 0.3 % of the locus — so
  most of them are thinner than a line. Fattening them for legibility would
  draw a gene whose coding fraction looks like 40 % when it is 15 %, and the
  figure's whole argument is about which parts of the locus carry coding
  sequence. Legibility comes from a minimum *stroke* on a zero-width mark
  instead, which does not distort the extent.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.patches import Patch, Rectangle                # noqa: E402

import figstyle as fs                                          # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "annotation_bugs"
FIGS = OUT / "figures"

#: What the annotation holds at an exon. Ordered worst-to-best so the legend
#: reads as a scale, and drawn from the project's validated evidence ramp
#: rather than from fresh colours.
EXON_COLOUR = {
    "in_annotated_cds": fs.QUALITY["high"] if "high" in fs.QUALITY
    else fs.BLUES[4],
    "in_annotated_noncoding": fs.FAINT,
    "intergenic": fs.CLINICAL["pathogenic"],
}
EXON_LABEL = {
    "in_annotated_cds": "in an annotated CDS",
    "in_annotated_noncoding": "inside a gene, not coding",
    "intergenic": "no annotated gene",
}


def load(name: str) -> list[dict]:
    path = OUT / name
    if not path.exists():
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))


def _by_case(rows: list[dict]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for r in rows:
        out.setdefault(r["case_id"], []).append(r)
    return out


# --------------------------------------------------------------------------
def fig_exon_tracks() -> None:
    """Figure 1 — the exon track: what the annotation has where the gene is."""
    cases = load("cases.tsv")
    exons = _by_case(load("case_exons.tsv"))
    models = _by_case(load("annotated_models.tsv"))
    introns = _by_case(load("case_introns.tsv"))
    neigh = _by_case(load("locus_neighbourhood.tsv"))
    if not cases or not exons:
        return
    fs.use()
    fig, axes = plt.subplots(len(cases), 1, figsize=(fs.W_FULL, 4.4))
    axes = [axes] if len(cases) == 1 else list(axes)
    for ax, case, letter in zip(axes, cases, "ABCD"):
        cid = case["case_id"]
        ex = exons.get(cid, [])
        lo = min(int(e["start"]) for e in ex)
        hi = max(int(e["end"]) for e in ex)
        kb = 1000.0

        # track 1 — the recovered gene model
        for e in ex:
            s, w = int(e["start"]), int(e["end"]) - int(e["start"]) + 1
            ax.add_patch(Rectangle(((s - lo) / kb, 0.62), w / kb, 0.22,
                                   facecolor=EXON_COLOUR[e["class"]],
                                   edgecolor=EXON_COLOUR[e["class"]],
                                   linewidth=0.35))
        ax.hlines(0.73, 0, (hi - lo) / kb, color=fs.MUTED, linewidth=0.4,
                  zorder=0)

        # track 2 — the annotated gene models over the same coordinates
        for m in models.get(cid, []):
            s, e2 = int(m["start"]), int(m["end"])
            colour = fs.ACCENT if m["pseudo"] == "1" else fs.BLUES[3]
            ax.add_patch(Rectangle(((s - lo) / kb, 0.28),
                                   (e2 - s + 1) / kb, 0.18,
                                   facecolor=colour, edgecolor=colour,
                                   linewidth=0.35))
            ax.text(((s + e2) / 2 - lo) / kb, 0.24,
                    m["name"] + (" (pseudogene)" if m["pseudo"] == "1" else ""),
                    ha="center", va="top", fontsize=fs.FS_NOTE - 0.6,
                    color=fs.MUTED)
        if not models.get(cid):
            ax.text(((hi - lo) / 2) / kb, 0.37,
                    "no annotated gene model on any exon",
                    ha="center", va="center", fontsize=fs.FS_NOTE,
                    color=fs.CLINICAL["pathogenic"], style="italic")

        # the junctions the annotation broke the gene at
        for i in introns.get(cid, []):
            if i["class"] != "between_models":
                continue
            x = ((int(i["start"]) + int(i["end"])) / 2 - lo) / kb
            # Stops above the model track: drawn down to the labels it
            # strikes through the gene names, which is where the reader looks
            # to see *which* two models the split is between.
            ax.vlines(x, 0.46, 0.90, color=fs.INK, linewidth=0.7,
                      linestyle=(0, (2, 1.4)))
            ax.text(x, 0.93, "split", ha="center", va="bottom",
                    fontsize=fs.FS_NOTE - 0.4, color=fs.INK)

        # The flanking genes, named — the evidence that the annotation was
        # working right up to the gap. Only the *nearest* on each side is
        # drawn: the table holds five per side and they sit within a few kb of
        # each other, so drawing them all overprints one label on the next and
        # the panel loses the one thing this track is for.
        for side in ("upstream", "downstream"):
            near = sorted((n for n in neigh.get(cid, [])
                           if n["side"] == side and int(n["distance_bp"]) >= 0),
                          key=lambda n: int(n["distance_bp"]))
            if not near:
                continue
            n = near[0]
            d = int(n["distance_bp"])
            x = 0.0 if side == "upstream" else (hi - lo) / kb
            ax.annotate(f"{n['name']}  ({d / kb:.0f} kb)",
                        xy=(x, 0.08),
                        xytext=(-4 if side == "upstream" else 4, 0),
                        textcoords="offset points",
                        ha="right" if side == "upstream" else "left",
                        va="center", fontsize=fs.FS_NOTE - 0.4,
                        color=fs.MUTED)
        ax.hlines(0.08, 0, (hi - lo) / kb, color=fs.GRID,
                  linewidth=0.5, zorder=0)

        ax.set_xlim(-((hi - lo) / kb) * 0.22, ((hi - lo) / kb) * 1.22)
        ax.set_ylim(0.0, 1.05)
        ax.set_yticks([])
        ax.set_xlabel(f"kb along {case['contig']} "
                      f"(locus starts at {lo:,})", fontsize=fs.FS_TICK)
        fs.despine(ax, keep=("bottom",))
        n_ex = len(ex)
        n_int = sum(1 for e in ex if e["class"] == "intergenic")
        fs.panel(ax, letter,
                 f"{case['organism']} {case['cell']} — {n_int} of {n_ex} "
                 f"exons have no annotated gene ({case['selected_as']})")
    axes[0].legend(handles=[Patch(facecolor=EXON_COLOUR[k], label=EXON_LABEL[k])
                            for k in EXON_COLOUR]
                   + [Patch(facecolor=fs.BLUES[3], label="annotated model"),
                      Patch(facecolor=fs.ACCENT, label="annotated pseudogene")],
                   loc="upper center", bbox_to_anchor=(0.5, 1.42), ncol=5,
                   fontsize=fs.FS_TICK - 0.3)
    fig.tight_layout(h_pad=2.6)
    fs.save(fig, FIGS / "exon_tracks")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_annotation_loss() -> None:
    """Figure 2 — the denominator: how rare these failures are."""
    rank = load("case_ranking.tsv")
    if not rank:
        return
    elig = [r for r in rank if r["eligible"] == "True"]
    fs.use()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.W_FULL, 2.5))

    losses = sorted(float(r["loss"]) for r in elig)
    bins = [i / 20 for i in range(21)]
    ax1.hist(losses, bins=bins, color=fs.BLUES[3], edgecolor=fs.SURFACE,
             linewidth=0.4)
    ax1.set_yscale("log")
    ax1.set_xlabel("annotation loss  (1 − best single model's share of the CDS)")
    ax1.set_ylabel("eligible loci")
    ax1.axvline(0.5, color=fs.INK, linewidth=0.7, linestyle=(0, (2, 1.4)))
    ax1.text(0.51, ax1.get_ylim()[1] * 0.5, "failure", fontsize=fs.FS_NOTE,
             color=fs.INK, ha="left", va="top")
    fs.hgrid(ax1)
    fs.panel(ax1, "A", f"{len(elig)} loci the annotation could have got right")

    # the E5 control, drawn as the thing it is: locus span against the longest
    # gene the same annotation builds anywhere.
    for r in rank:
        if r["failed_rule"] not in ("", "E5"):
            continue
        x = int(r["span"]) / 1000
        y = max(1, int(r["annot_max_gene_span"] or 1)) / 1000
        excluded = r["failed_rule"] == "E5"
        ax2.scatter(x, y, s=6 if not excluded else 14,
                    facecolor=fs.CLINICAL["pathogenic"] if excluded
                    else fs.FAINT,
                    edgecolor="none", alpha=0.9 if excluded else 0.35,
                    zorder=3 if excluded else 1)
    lim = [10, 3000]
    ax2.plot(lim, lim, color=fs.INK, linewidth=0.6, linestyle=(0, (2, 1.4)))
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(*lim)
    ax2.set_ylim(*lim)
    ax2.set_xlabel("locus span (kb)")
    ax2.set_ylabel("longest gene the same\nannotation builds (kb)")
    ax2.legend(handles=[
        Patch(facecolor=fs.FAINT, label="passes E5"),
        Patch(facecolor=fs.CLINICAL["pathogenic"], label="excluded by E5")],
        loc="lower right", fontsize=fs.FS_TICK - 0.3)
    fs.hgrid(ax2, axis="both")
    fs.panel(ax2, "B", "E5: can this annotation build a gene this long?")
    fig.tight_layout()
    fs.save(fig, FIGS / "annotation_loss")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_tiling() -> None:
    """Figure 3 — which residues of the protein each annotated model delivers."""
    cases = load("cases.tsv")
    tiles = _by_case(load("fragment_tiling.tsv"))
    if not cases or not tiles:
        return
    fs.use()
    fig, axes = plt.subplots(len(cases), 1, figsize=(fs.W_FULL, 3.2))
    axes = [axes] if len(cases) == 1 else list(axes)
    for ax, case, letter in zip(axes, cases, "ABCD"):
        cid = case["case_id"]
        rows = [r for r in tiles.get(cid, [])
                if r["tiles"] == "1" and r["at_own_locus"] == "1"]
        want = f"|{case['cell']}|{case['contig']}:{case['start']}-{case['end']}"
        mine = [r for r in rows if want in r["best_locus"]]
        length = max([int(r["subject_aa"]) for r in mine], default=0)
        if not length:
            length = int(case.get("cds_bp", 0) or 3) // 3
        ax.add_patch(Rectangle((1, 0.55), length, 0.16, facecolor=fs.GRID,
                               edgecolor="none"))
        ax.text(length + length * 0.01, 0.63, f"{length:,} aa",
                va="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
        # Labels alternate between two rows. A 66-residue model and a
        # 415-residue one starting 52 residues apart put their labels on top of
        # each other at any single height, and the pair that overlaps is
        # exactly the pair a reader needs to tell apart.
        for i, r in enumerate(sorted(mine, key=lambda x: int(x["s_start"]))):
            s, e = int(r["s_start"]), int(r["s_end"])
            colour = fs.ACCENT if r["pseudo"] == "1" else fs.BLUES[4]
            ax.add_patch(Rectangle((s, 0.55), e - s + 1, 0.16,
                                   facecolor=colour, edgecolor=fs.SURFACE,
                                   linewidth=0.5))
            ax.text((s + e) / 2, 0.50 - 0.16 * (i % 2),
                    f"{r['name']}\n{s}–{e}", ha="center", va="top",
                    fontsize=fs.FS_NOTE - 0.8, color=fs.MUTED)
        if not mine:
            # An empty bar states a fact but not the interesting one. Where the
            # genome's family-named models went instead is the result.
            elsewhere = [r for r in tiles.get(cid, [])
                         if r.get("naming") == "name_mismatch"]
            note = "no annotated model delivers any residue of this gene"
            if elsewhere:
                r = elsewhere[0]
                loc = r["best_locus"].split("|")
                note += (f"\nthe genome's `{r['name']}`-named model "
                         f"({int(r['protein_aa']):,} aa) delivers residues "
                         f"{r['s_start']}–{r['s_end']} of the "
                         f"{loc[1] if len(loc) > 1 else '?'} locus instead")
            ax.text(length / 2, 0.36, note.replace("`", ""), ha="center",
                    va="top", fontsize=fs.FS_NOTE, color=fs.MUTED,
                    style="italic")
        covered = sum(int(r["s_end"]) - int(r["s_start"]) + 1 for r in mine)
        ax.set_xlim(-length * 0.02, length * 1.12)
        ax.set_ylim(0.10, 0.95)
        ax.set_yticks([])
        ax.set_xlabel("residue of the recovered protein", fontsize=fs.FS_TICK)
        fs.despine(ax, keep=("bottom",))
        fs.panel(ax, letter, f"{case['organism']} {case['cell']} — annotated "
                 f"models deliver {covered:,} of {length:,} residues "
                 f"({covered / max(1, length):.0%})")
    fig.tight_layout(h_pad=2.2)
    fs.save(fig, FIGS / "fragment_tiling")
    plt.close(fig)


# --------------------------------------------------------------------------
def fig_validation() -> None:
    """Figure 4 — the checks on S10's own evidence."""
    conc = _by_case(load("boundary_concordance.tsv"))
    introns = _by_case(load("case_introns.tsv"))
    orf = load("reading_frame.tsv")
    cases = load("cases.tsv")
    if not cases:
        return
    fs.use()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(fs.W_FULL, 2.3))

    # Both curves saturate at 1.0 over most of their length, so plain markers
    # of equal size hide whichever is drawn first. Sizes and fills differ, and
    # the later curve is drawn unfilled, so the overlap reads as overlap.
    for i, case in enumerate(cases):
        rows = conc.get(case["case_id"], [])
        if not rows:
            continue
        fracs = sorted(float(r["frac_refs"]) for r in rows)
        colour = fs.PARALOG.get(case["cell"], fs.BLUES[3])
        ax1.plot(range(1, len(fracs) + 1), fracs,
                 marker="o" if i == 0 else "s",
                 markersize=3.4 if i == 0 else 2.2,
                 markerfacecolor=colour if i == 0 else "none",
                 markeredgecolor=colour, markeredgewidth=0.7,
                 linewidth=0.8, color=colour,
                 label=f"{case['organism'].split()[0][0]}. "
                       f"{case['organism'].split()[1]} {case['cell']}")
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("exon boundary (ranked)")
    ax1.set_ylabel("fraction of reference\ngenomes sharing it")
    ax1.legend(loc="lower right", fontsize=fs.FS_TICK - 0.6)
    fs.hgrid(ax1)
    # "independently annotated", not "RefSeq": the reference set is whichever
    # swept genomes corroborate the alignment, and it is a mixture of RefSeq
    # and submitter annotations. Naming one of them would overstate it.
    fs.panel(ax1, "A", "Boundaries vs annotated genomes")

    labels, values, colours = [], [], []
    for case in cases:
        counts = Counter(r["splice_class"] for r in introns.get(case["case_id"], []))
        for klass in ("canonical", "minor", "non_canonical", "unreadable"):
            if not counts.get(klass):
                continue
            labels.append(f"{case['cell']}\n{klass}")
            values.append(counts[klass])
            colours.append(fs.BLUES[4] if klass == "canonical" else fs.FAINT)
    ax2.bar(range(len(values)), values, color=colours, width=0.68)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=fs.FS_TICK - 1.2)
    ax2.set_ylabel("introns")
    fs.hgrid(ax2)
    fs.panel(ax2, "B", "Splice-site dinucleotides")

    rows = [r for r in orf if r.get("verdict") in
            ("open_reading_frame", "disrupted")]
    if rows:
        xs = [float(r["expected_stops"] or 0) for r in rows]
        ys = [int(r["internal_stops"]) for r in rows]
        cols = [fs.ACCENT if r["is_case_locus"] == "1" else fs.FAINT
                for r in rows]
        ax3.scatter(xs, ys, c=cols, s=16, zorder=3)
        top = max(xs + [1]) * 1.15
        ax3.plot([0, top], [0, top], color=fs.INK, linewidth=0.6,
                 linestyle=(0, (2, 1.4)))
        ax3.set_xlim(0, top)
        ax3.set_ylim(-0.6, max(ys + [1]) + 1.0)
        ax3.legend(handles=[Patch(facecolor=fs.ACCENT, label="case locus"),
                            Patch(facecolor=fs.FAINT,
                                  label="other family locus,\nsame genome")],
                   loc="upper left", fontsize=fs.FS_TICK - 0.6)
    ax3.set_xlabel("stops expected if neutral")
    ax3.set_ylabel("internal stops observed")
    fs.hgrid(ax3)
    fs.panel(ax3, "C", "Is the reading frame open?")
    if rows:
        ax3.annotate("no locus carries a single\ninternal stop",
                     xy=(0.96, 0.52), xycoords="axes fraction",
                     fontsize=fs.FS_NOTE, color=fs.MUTED, ha="right")
    fig.tight_layout()
    fs.save(fig, FIGS / "validation")
    plt.close(fig)


def main() -> int:
    FIGS.mkdir(parents=True, exist_ok=True)
    for fn in (fig_exon_tracks, fig_annotation_loss, fig_tiling, fig_validation):
        try:
            fn()
            print(f"[s10] figure {fn.__name__}")
        except Exception as exc:                      # noqa: BLE001
            print(f"[s10] !! {fn.__name__}: {type(exc).__name__}: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
