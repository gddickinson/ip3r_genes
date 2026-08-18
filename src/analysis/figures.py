"""Publication-style figures for the analysis report.

All figures are optional — when matplotlib isn't installed, the report
writer falls back to ASCII / no-image content. The figures we render:

    identity_heatmap.png  — pairwise % identity heatmap
    tree.png              — neighbor-joining tree (Bio.Phylo + matplotlib)
    conservation.png      — per-column conservation plot, with variable
                            regions shaded
    length_histogram.png  — histogram of protein lengths, coloured by gene
    sources_bar.png       — variants per source × gene
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Optional

try:  # matplotlib optional
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    HAVE_MPL = True
except ImportError:  # pragma: no cover
    HAVE_MPL = False

from ..core.models import ProteinVariant
from .distance import DistanceTable
from .evolution import VariableRegion


def identity_heatmap(
    table: DistanceTable, out_path: Path,
    title: str = "Pairwise sequence identity",
) -> Optional[Path]:
    if not HAVE_MPL or not table.labels:
        return None
    n = len(table.labels)
    fig_size = max(6, min(16, 0.4 * n + 4))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(table.identity, cmap="viridis", vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(table.labels, rotation=90, fontsize=7)
    ax.set_yticklabels(table.labels, fontsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Identity")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def tree_figure(
    tree, out_path: Path, title: str = "Neighbor-joining tree",
    highlight_labels: Optional[set[str]] = None,
    highlight_color: str = "#c0392b",
) -> Optional[Path]:
    """Render the tree; tip labels in `highlight_labels` are drawn bold in
    `highlight_color` (used to flag novel candidates)."""
    if not HAVE_MPL or tree is None:
        return None
    try:
        from Bio import Phylo
    except ImportError:
        return None
    highlight_labels = highlight_labels or set()
    n_terminals = max(1, len(tree.get_terminals()))
    fig, ax = plt.subplots(figsize=(11, max(4, 0.3 * n_terminals + 2)))
    Phylo.draw(tree, axes=ax, do_show=False, branch_labels=None)
    for text in ax.texts:
        if text.get_text().strip() in highlight_labels:
            text.set_color(highlight_color)
            text.set_fontweight("bold")
    if highlight_labels:
        title += "  ·  novel candidates in red"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def conservation_plot(
    conservation: list[float],
    regions: list[VariableRegion],
    out_path: Path,
    title: str = "Per-column conservation (Shannon-derived)",
) -> Optional[Path]:
    if not HAVE_MPL or not conservation:
        return None
    fig, ax = plt.subplots(figsize=(12, 3.5))
    xs = list(range(1, len(conservation) + 1))
    ax.plot(xs, conservation, linewidth=0.8, color="#2c3e50")
    ax.fill_between(xs, conservation, 0, alpha=0.2, color="#2c3e50")
    for r in regions:
        ax.add_patch(Rectangle(
            (r.start - 0.5, 0), r.width(), 1.0,
            color="#e74c3c", alpha=0.15,
        ))
    ax.set_xlim(1, len(conservation))
    ax.set_ylim(0, 1)
    ax.set_xlabel("MSA column")
    ax.set_ylabel("Conservation (1 = identical)")
    ax.set_title(f"{title}  ·  {len(regions)} variable region(s) highlighted")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def length_histogram(variants: list[ProteinVariant], out_path: Path) -> Optional[Path]:
    if not HAVE_MPL or not variants:
        return None
    by_gene: dict[str, list[int]] = defaultdict(list)
    for v in variants:
        if v.length_aa:
            by_gene[v.gene_symbol.upper()].append(v.length_aa)
    if not by_gene:
        return None
    fig, ax = plt.subplots(figsize=(9, 4))
    bins_max = max(max(vals) for vals in by_gene.values())
    bins = list(range(0, bins_max + 100, 100))
    for gene, lens in sorted(by_gene.items()):
        ax.hist(lens, bins=bins, alpha=0.6, label=f"{gene} (n={len(lens)})")
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Variants")
    ax.set_title("Protein-length distribution per gene")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def sources_bar(variants: list[ProteinVariant], out_path: Path) -> Optional[Path]:
    if not HAVE_MPL or not variants:
        return None
    counts: dict[tuple[str, str], int] = defaultdict(int)
    for v in variants:
        counts[(v.gene_symbol.upper(), v.source)] += 1
    if not counts:
        return None
    genes = sorted({g for g, _ in counts})
    sources = sorted({s for _, s in counts})
    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(genes) + 3), 4))
    width = 0.8 / max(1, len(sources))
    for i, src in enumerate(sources):
        vals = [counts.get((g, src), 0) for g in genes]
        xs = [j + i * width for j in range(len(genes))]
        ax.bar(xs, vals, width=width, label=src)
    ax.set_xticks([j + width * (len(sources) - 1) / 2 for j in range(len(genes))])
    ax.set_xticklabels(genes)
    ax.set_ylabel("Variants")
    ax.set_title("Variants per gene × source")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
