"""Supplementary Figures 3-4: the inputs the constraint map and every omega
estimate were computed on.

Split out of `s24_figs_alignment.py` to keep both inside the project's
500-line budget; `binned` moved to `s24_lib` at the same time so the two
halves cannot bin a profile differently and draw curves that are not
comparable.

Two decisions worth naming.

* **Supplementary Fig. 3's occupancy panel is per residue of the human
  reference, not per alignment column.** The three deep alignments have three
  different widths and no shared coordinate; the reference each was built
  around is the only axis all three can be read on.

* **Supplementary Fig. 4 plots codon columns, not nucleotide columns.** The
  codon alignment is 9,759 nucleotides wide and every claim in S9 is about
  codons; an axis in nucleotides would make a three-fold difference between
  the two look like a property of the data.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402

import figstyle as F                                           # noqa: E402
import s24_lib as L                                            # noqa: E402

GAPS = "-."

# ------------------------------------------------------------------ fig 3

def fig_paralog_alignments(out: Path, stats: dict) -> None:
    shape = L.read_tsv(L.CONSTRAINT_DIR / "ortholog_shape.tsv")
    layers = L.read_tsv(L.CONSTRAINT_DIR / "layer_sizes.tsv")
    manifest = L.read_tsv(L.CONSTRAINT_DIR / "orthologs_manifest.tsv")

    fig, axes = plt.subplots(3, 1, figsize=(F.W_FULL, 6.2))

    ax = axes[0]
    for paralog in L.PARALOGS:
        rows = L.constraint_table(paralog)
        xs = [int(r["resi"]) for r in rows]
        ys = [L.f(r["deep_occupancy"], 0.0) for r in rows]
        bx, by = L.binned(ys, 400)
        ax.plot([xs[0] + b for b in bx], by, lw=0.8,
                color=F.PARALOG[paralog], label=paralog)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("residue of the human reference (each paralogue's own "
                  "UniProt numbering)")
    ax.set_ylabel("occupancy")
    ax.legend(loc="lower right", ncol=3)
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "a", "how many of the 249–265 orthologues cover each residue")

    ax = axes[1]
    for i, paralog in enumerate(L.PARALOGS):
        rows = [r for r in shape if r["paralog"] == paralog]
        vals = sorted(L.f(r["frac_in_ref_columns"], 0.0) for r in rows)
        ys = [i + (k - len(vals) / 2) / max(len(vals), 1) * 0.55
              for k in range(len(vals))]
        drop = [(v, y) for v, y, r in zip(vals, ys, sorted(
            rows, key=lambda q: L.f(q["frac_in_ref_columns"], 0.0)))
            if r["kept"] != "True"]
        ax.scatter(vals, ys, s=4, color=F.PARALOG[paralog], lw=0, alpha=0.55)
        if drop:
            ax.scatter([d[0] for d in drop], [d[1] for d in drop], s=22,
                       facecolor="none", edgecolor=F.CLINICAL["pathogenic"],
                       lw=0.9, zorder=4)
        bar = {L.f(r["bar"]) for r in rows}
        for b in bar:
            ax.plot([b, b], [i - 0.42, i + 0.42], color=F.INK, lw=0.9)
    ax.set_yticks(range(3))
    ax.set_yticklabels(L.PARALOGS)
    ax.set_xlabel("fraction of the sequence's own residues that land in a "
                  "reference column")
    F.despine(ax)
    F.hgrid(ax, "x")
    F.panel(ax, "b", "the shape screen bait coverage cannot do; the ring is "
                     "the one sequence it dropped")

    ax = axes[2]
    order = ["deep", "shallow", "vert", "family"]
    xs, hs, cs, texts = [], [], [], []
    pos = 0
    for lname in order:
        for r in layers:
            if r["layer"] != lname:
                continue
            xs.append(pos)
            hs.append(int(r["n_sequences"]))
            cs.append(F.PARALOG.get(r["paralog"], F.MUTED))
            texts.append(f"{lname}\n{r['paralog']}")
            pos += 1
        pos += 0.6
    ax.bar(xs, hs, color=cs, width=0.8)
    for x, h in zip(xs, hs):
        ax.annotate(str(h), (x, h), xytext=(0, 2),
                    textcoords="offset points", ha="center",
                    fontsize=F.FS_NOTE, color=F.MUTED)
    ax.set_xticks(xs)
    ax.set_xticklabels(texts, fontsize=F.FS_NOTE)
    ax.set_ylabel("sequences")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "c", "the four conservation layers, and the depth control "
                     "beside each deep set")

    fig.tight_layout()
    F.save(fig, out / "SuppFig3_paralog_alignments")
    plt.close(fig)

    dropped = [r for r in shape if r["kept"] != "True"]
    stats["supp_fig_3"] = {
        "sequences_screened": len(shape),
        "sequences_dropped": len(dropped),
        "dropped": [r["name"] for r in dropped],
        "layer_sizes": {f"{r['layer']}:{r['paralog']}": int(r["n_sequences"])
                        for r in layers},
        "manifest_rows": len(manifest),
    }


# ------------------------------------------------------------------ fig 4

def fig_codon_alignment(out: Path, stats: dict) -> None:
    aln = L.read_fasta(L.SELECTION_DIR / "codon_aln.fasta")
    trimmed = L.read_fasta(L.SELECTION_DIR / "codon_trimmed.fasta")
    codes = {r["code"]: r for r in L.read_tsv(L.SELECTION_DIR /
                                              "tip_codes.tsv")}
    status = {r["label"]: r for r in L.read_tsv(L.SELECTION_DIR /
                                                "cds_status.tsv")}
    ncod = len(next(iter(aln.values()))) // 3
    ntrim = len(next(iter(trimmed.values()))) // 3

    occ = []
    for c in range(ncod):
        cod = [s[3 * c:3 * c + 3] for s in aln.values()]
        occ.append(sum(1 for x in cod if x and x[0] not in GAPS) / len(cod))

    fig, axes = plt.subplots(3, 1, figsize=(F.W_FULL, 6.0))

    ax = axes[0]
    xs, ys = L.binned(occ, 500)
    ax.fill_between(xs, ys, color=F.BLUES[1], lw=0)
    ax.plot(xs, ys, color=F.BLUES[4], lw=0.7)
    ax.set_xlim(0, ncod)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("occupancy")
    ax.set_xlabel(f"codon of the {ncod:,}-codon PAL2NAL alignment")
    F.despine(ax)
    F.hgrid(ax)
    # The note goes in the title, not in the panel: the fill reaches the top
    # of the axes almost everywhere, so any in-plot annotation lands on the
    # data it is describing.
    F.panel(ax, "a", f"the codon alignment, per-codon occupancy over 57 tips; "
                     f"trimAl kept {ntrim:,} of {ncod:,} codons "
                     f"({100 * ntrim / ncod:.0f} %)")

    ax = axes[1]
    sets = ["ITPR1", "ITPR2", "ITPR3", "background"]
    pos = 0
    ticks, ticklabels = [], []
    for s in sets:
        members = [c for c in trimmed if codes.get(c, {}).get("set") == s]
        members.sort(key=lambda c: codes[c]["species"])
        cov = [sum(1 for k in range(ntrim)
                   if trimmed[c][3 * k] not in GAPS) / ntrim for c in members]
        ax.bar(range(pos, pos + len(members)), cov, width=0.85,
               color=F.PARALOG.get(s, F.MUTED))
        ticks.append(pos + len(members) / 2 - 0.5)
        ticklabels.append(f"{s}\n(n = {len(members)})")
        pos += len(members) + 2
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, fontsize=F.FS_NOTE)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("coverage of the\ntrimmed codon alignment")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "b", "per-tip coverage, by the selection set the tree put the "
                     "tip in")

    ax = axes[2]
    rows = [status[codes[c]["label"]] for c in trimmed
            if codes[c]["label"] in status]
    by_route = defaultdict(list)
    for r in rows:
        by_route[r["route"]].append(int(r["n_masked"] or 0))
    xs, hs, cs = [], [], []
    for i, (route, vals) in enumerate(sorted(by_route.items())):
        jitter = [i + (k - len(vals) / 2) / max(len(vals), 1) * 0.6
                  for k in range(len(vals))]
        ax.scatter(jitter, sorted(vals), s=9, lw=0,
                   color=F.BLUES[4] if route == "uniprot" else F.BLUES[2],
                   label=f"{route} (n = {len(vals)})")
        xs.append(i)
    ax.set_xticks(range(len(by_route)))
    ax.set_xticklabels(sorted(by_route), fontsize=F.FS_TICK)
    ax.set_ylabel("codons masked to NNN")
    ax.legend(loc="upper left", ncol=2)
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "c", "what validation masked, by the route the coding "
                     "sequence came from")

    fig.tight_layout()
    F.save(fig, out / "SuppFig4_codon_alignment")
    plt.close(fig)

    stats["supp_fig_4"] = {
        "tips": len(trimmed), "codons_in": ncod, "codons_kept": ntrim,
        "pct_kept": round(100 * ntrim / ncod, 1),
        "routes": {k: len(v) for k, v in by_route.items()},
        "codons_masked_total": sum(sum(v) for v in by_route.values()),
        "sets": dict(Counter(codes[c]["set"] for c in trimmed)),
    }


FIGURES = {
    "supp3": fig_paralog_alignments,
    "supp4": fig_codon_alignment,
}
