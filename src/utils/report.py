"""Render a publication-style report for a results bundle.

Produces:
    report.md   — GitHub-flavoured markdown with embedded figure references
    report.html — standalone HTML (no external CSS) for sharing

The markdown is the canonical version; HTML is generated from it with a
minimal renderer so users don't need pandoc / markdown-cli installed.
"""

from __future__ import annotations

import html
from datetime import datetime
from pathlib import Path
from typing import Optional

from ..analysis.figures import (
    HAVE_MPL,
    conservation_plot,
    identity_heatmap,
    length_histogram,
    sources_bar,
    tree_figure,
)
from ..analysis.pipeline import AnalysisResult
from ..analysis.tree import build_nj_tree
from ..core.models import ProteinVariant, SearchQuery


def write_report(
    out_dir: Path,
    variants: list[ProteinVariant],
    query: SearchQuery,
    analysis: Optional[AnalysisResult] = None,
    discovery_summary: str = "",
    notes: str = "",
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    fig_paths: dict[str, Optional[Path]] = {}
    if HAVE_MPL:
        fig_paths["sources_bar"] = sources_bar(variants, figures_dir / "sources_bar.png")
        fig_paths["lengths"] = length_histogram(variants, figures_dir / "length_histogram.png")
        if analysis:
            if analysis.distances:
                fig_paths["heatmap"] = identity_heatmap(
                    analysis.distances, figures_dir / "identity_heatmap.png",
                )
            if analysis.newick:
                tree_obj = build_nj_tree(analysis.distances) if analysis.distances else None
                if tree_obj is not None:
                    fig_paths["tree"] = tree_figure(tree_obj, figures_dir / "tree.png")
            if analysis.conservation:
                fig_paths["conservation"] = conservation_plot(
                    analysis.conservation,
                    analysis.variable_regions,
                    figures_dir / "conservation.png",
                )

    md = _render_markdown(variants, query, analysis, discovery_summary, notes, fig_paths)
    md_path = out_dir / "report.md"
    md_path.write_text(md)

    html_path = out_dir / "report.html"
    html_path.write_text(_md_to_html(md, title=f"Protein Variant Report — {', '.join(query.gene_symbols)}"))

    return {
        "report_md": str(md_path),
        "report_html": str(html_path),
        "figures": {k: (str(p) if p else "") for k, p in fig_paths.items()},
    }


def _render_markdown(
    variants: list[ProteinVariant],
    query: SearchQuery,
    analysis: Optional[AnalysisResult],
    discovery_summary: str,
    notes: str,
    fig_paths: dict[str, Optional[Path]],
) -> str:
    L: list[str] = []
    L.append(f"# Protein Variant Report — {', '.join(query.gene_symbols)}")
    L.append("")
    L.append(f"_Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}_")
    L.append("")
    L.append("## 1. Query")
    L.append("")
    L.append(f"- **Genes searched**: `{', '.join(query.gene_symbols)}`")
    L.append(f"- **Species filter**: `{query.species or '(any)'}`")
    L.append(f"- **Sources queried**: `{', '.join(query.sources)}`")
    L.append(f"- **Sequences fetched**: `{'yes' if query.include_sequence else 'no'}`")
    L.append(f"- **Total variants returned**: **{len(variants)}**")
    if notes:
        L.append("")
        L.append(f"_Notes:_ {notes}")

    # 2. Per-source overview
    L.append("")
    L.append("## 2. Per-source overview")
    L.append("")
    if fig_paths.get("sources_bar"):
        L.append(f"![Variants per gene × source](figures/sources_bar.png)")
        L.append("")
    if fig_paths.get("lengths"):
        L.append(f"![Protein length distribution](figures/length_histogram.png)")
        L.append("")

    # 3. Cross-database discrepancies
    L.append("## 3. Cross-database length agreement")
    L.append("")
    L.append("Disagreement between sources on the maximum protein length flags "
             "incomplete or fragmented annotation, which is how a real gene ends "
             "up recorded as two short models or none at all.")
    L.append("")
    by_gene_src: dict[str, dict[str, int]] = {}
    for v in variants:
        if not v.length_aa:
            continue
        by_gene_src.setdefault(v.gene_symbol.upper(), {}).setdefault(v.source, 0)
        by_gene_src[v.gene_symbol.upper()][v.source] = max(
            by_gene_src[v.gene_symbol.upper()][v.source], v.length_aa,
        )
    if by_gene_src:
        L.append("| Gene | Per-source max length | Δ | Flag |")
        L.append("|------|-----------------------|---|------|")
        for gene, src_max in sorted(by_gene_src.items()):
            if len(src_max) < 2:
                continue
            vals = list(src_max.values())
            delta = max(vals) - min(vals)
            flag = "⚠️" if delta >= 100 else "ok"
            cells = ", ".join(f"{s}={n}" for s, n in sorted(src_max.items()))
            L.append(f"| {gene} | {cells} | {delta} | {flag} |")

    # 4. Sequence analysis
    if analysis:
        L.append("")
        L.append("## 4. Sequence analysis")
        L.append("")
        L.append(f"- **Variants with sequences**: {analysis.n_analyzed} / {analysis.n_input}")
        if analysis.distances and analysis.distances.labels:
            n = len(analysis.distances.labels)
            mean_id = sum(
                analysis.distances.identity[i][j]
                for i in range(n) for j in range(i + 1, n)
            ) / max(1, n * (n - 1) // 2)
            L.append(f"- **Mean pairwise identity**: {mean_id:.1%}")
        if analysis.clusters:
            L.append(f"- **Identity-based clusters**: {len(analysis.clusters.clusters)}")
            L.append(f"- **Novel-group candidates** (nearest < threshold): "
                     f"{len(analysis.clusters.novel_candidates)}")
        if analysis.variable_regions:
            L.append(f"- **Variable regions detected**: {len(analysis.variable_regions)}")
        L.append("")
        if fig_paths.get("heatmap"):
            L.append("### Pairwise identity heatmap")
            L.append("")
            L.append("![Identity heatmap](figures/identity_heatmap.png)")
            L.append("")
        if fig_paths.get("tree"):
            L.append("### Neighbor-joining tree")
            L.append("")
            L.append("![NJ tree](figures/tree.png)")
            L.append("")
            L.append("Newick form (paste into FigTree / iTOL):")
            L.append("```")
            L.append(analysis.newick or "(empty)")
            L.append("```")
            L.append("")
        if fig_paths.get("conservation"):
            L.append("### Per-column conservation")
            L.append("")
            L.append("![Conservation](figures/conservation.png)")
            L.append("")
        if analysis.clusters and analysis.clusters.clusters:
            L.append("### Identity-based clusters")
            L.append("")
            L.append(f"_(single-link clustering at ≥ {analysis.clusters.threshold:.0%} identity)_")
            L.append("")
            for i, cluster in enumerate(analysis.clusters.clusters, 1):
                L.append(f"**Cluster {i}** — {len(cluster)} member(s):")
                for m in cluster:
                    L.append(f"  - `{m}`")
                L.append("")
        if analysis.clusters and analysis.clusters.novel_candidates:
            L.append("### Outlier / novel-group flags")
            L.append("")
            L.append("| Variant | Nearest in set | Identity |")
            L.append("|---------|----------------|----------|")
            for label, neighbour, ident in analysis.clusters.novel_candidates:
                L.append(f"| `{label}` | `{neighbour}` | {ident:.0%} |")
            L.append("")
        if analysis.variable_regions:
            L.append("### Variable regions")
            L.append("")
            L.append("| Cols | Width | Mean conservation |")
            L.append("|------|-------|-------------------|")
            for r in analysis.variable_regions[:20]:
                L.append(f"| {r.start}–{r.end} | {r.width()} | {r.mean_conservation:.2f} |")
            L.append("")
        if analysis.mutations:
            L.append(f"### Mutations vs reference `{analysis.reference_label}`")
            L.append("")
            L.append(f"_({len(analysis.mutations)} total positions; first 25 shown)_")
            L.append("")
            L.append("| Query | Position | Change |")
            L.append("|-------|----------|--------|")
            for m in analysis.mutations[:25]:
                L.append(f"| `{m.query_label}` | {m.ref_pos} | `{m.short()}` |")
            L.append("")

    # 5. Discovery
    if discovery_summary:
        L.append("## 5. Novel-paralog discovery")
        L.append("")
        L.append("```")
        L.append(discovery_summary.rstrip())
        L.append("```")
        L.append("")

    # 6. Methods (so the report is self-contained for a colleague)
    L.append("## 6. Methods")
    L.append("")
    L.append(
        "Variants were aggregated from NCBI Entrez (`Bio.Entrez`), Ensembl REST "
        "and UniProt REST. Predicted 3D structures were pulled from the "
        "AlphaFold Protein Structure Database; structural homologs were "
        "(optionally) searched via Foldseek. Pairwise alignments use BLOSUM62 "
        "with open / extend gap penalties of −10 / −0.5. The multiple sequence "
        "alignment is a star alignment seeded by the longest sequence (or "
        "MAFFT --auto when available). Pairwise identity excludes columns where "
        "both rows are gaps. Neighbor-joining trees are constructed via "
        "Bio.Phylo.TreeConstruction with the 1 − identity distance matrix. "
        "Per-column conservation is 1 − Shannon entropy normalized to "
        "min(n_rows, 21). Novel-paralog candidates are scored on six independent "
        "criteria (size, domain signature, structural fold, sequence outlier, "
        "taxonomic breadth, cluster exclusion); see `discovery/candidates.tsv` "
        "for full evidence per candidate."
    )
    L.append("")
    L.append("---")
    L.append("")
    L.append("_Report generated by the Protein Variant Finder._")
    return "\n".join(L) + "\n"


def _md_to_html(md: str, title: str) -> str:
    """Tiny markdown → HTML renderer (no external deps).

    Handles: headings, paragraphs, code fences, inline code, images, links,
    tables, lists, bold/italic. Good enough for our report layout.
    """
    lines = md.splitlines()
    out: list[str] = [
        "<!doctype html>", '<html lang="en"><head><meta charset="utf-8">',
        f"<title>{html.escape(title)}</title>",
        "<style>",
        "body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;",
        "  max-width:980px;margin:2em auto;padding:0 1.5em;color:#222;line-height:1.5;}",
        "h1{border-bottom:2px solid #2c3e50;padding-bottom:0.3em;}",
        "h2{margin-top:1.8em;color:#2c3e50;border-bottom:1px solid #ddd;padding-bottom:0.2em;}",
        "h3{margin-top:1.4em;color:#34495e;}",
        "code,pre{font-family:'SF Mono',Menlo,Consolas,monospace;font-size:0.92em;}",
        "code{background:#f4f6f8;padding:0.1em 0.3em;border-radius:3px;}",
        "pre{background:#f4f6f8;padding:1em;border-radius:6px;overflow-x:auto;}",
        "table{border-collapse:collapse;margin:0.5em 0 1em 0;}",
        "th,td{border:1px solid #d0d7de;padding:0.4em 0.7em;text-align:left;}",
        "th{background:#f4f6f8;}",
        "img{max-width:100%;height:auto;border:1px solid #eee;border-radius:4px;}",
        "blockquote{border-left:4px solid #2c3e50;margin:0;padding-left:1em;color:#555;}",
        "</style></head><body>",
    ]
    in_pre = False
    in_table = False
    in_list = False
    for raw in lines:
        line = raw.rstrip()
        if line.startswith("```"):
            if in_pre:
                out.append("</pre>")
                in_pre = False
            else:
                out.append("<pre>")
                in_pre = True
            continue
        if in_pre:
            out.append(html.escape(line))
            continue
        if line.startswith("# "):
            out.append(f"<h1>{html.escape(line[2:])}</h1>"); continue
        if line.startswith("## "):
            out.append(f"<h2>{html.escape(line[3:])}</h2>"); continue
        if line.startswith("### "):
            out.append(f"<h3>{html.escape(line[4:])}</h3>"); continue
        if line.startswith("|") and "|" in line[1:]:
            # table — including separator row
            if "---" in line:
                continue
            cells = [c.strip() for c in line.strip("|").split("|")]
            if not in_table:
                out.append("<table>")
                in_table = True
                out.append("<tr>" + "".join(f"<th>{_inline(c)}</th>" for c in cells) + "</tr>")
            else:
                out.append("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in cells) + "</tr>")
            continue
        else:
            if in_table:
                out.append("</table>")
                in_table = False
        if line.startswith("- "):
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{_inline(line[2:])}</li>")
            continue
        else:
            if in_list:
                out.append("</ul>")
                in_list = False
        if not line:
            out.append("")
            continue
        if line.startswith("!["):
            # image: ![alt](path)
            try:
                alt = line[2:line.index("]")]
                src = line[line.index("(") + 1:line.index(")")]
                out.append(f'<p><img alt="{html.escape(alt)}" src="{html.escape(src)}"></p>')
                continue
            except ValueError:
                pass
        if line.startswith("---"):
            out.append("<hr>"); continue
        out.append(f"<p>{_inline(line)}</p>")
    if in_pre:
        out.append("</pre>")
    if in_table:
        out.append("</table>")
    if in_list:
        out.append("</ul>")
    out.append("</body></html>")
    return "\n".join(out)


def _inline(text: str) -> str:
    s = html.escape(text)
    # bold
    while "**" in s:
        a = s.find("**"); b = s.find("**", a + 2)
        if b == -1:
            break
        s = s[:a] + f"<strong>{s[a+2:b]}</strong>" + s[b+2:]
    # inline code
    import re
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    # links [text](url)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    return s
