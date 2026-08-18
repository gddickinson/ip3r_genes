"""Bundle a search result set into a timestamped results folder.

Layout produced:

    results/2026-08-18_191234_ip3r/
        results.csv         — full table, one row per ProteinVariant
        results.json        — same data, full fidelity
        results.fasta       — sequences (only variants that have them)
        metadata.json       — query, sources, timing, environment
        summary.md          — human-readable report with per-source counts,
                              length-distribution table, and Dong-style
                              cross-DB discrepancy section

The writer is invoked from two places:
    * the GUI File menu ("Save results bundle…")
    * the CLI `python run.py --headless --preset … --save-results`
"""

from __future__ import annotations

import json
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

from ..core.models import ProteinVariant, SearchQuery, SearchResult
from .exporters import write_csv, write_fasta, write_json


def _slugify(text: str) -> str:
    out = "".join(c if c.isalnum() else "_" for c in text.strip())
    return out.strip("_").lower() or "results"


def make_bundle_dir(root: Path, label: str = "results") -> Path:
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    dir_path = root / f"{stamp}_{_slugify(label)}"
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def write_bundle(
    out_dir: Path,
    variants: list[ProteinVariant],
    query: SearchQuery,
    per_source_results: Iterable[SearchResult] = (),
    notes: str = "",
) -> dict:
    """Write all bundle files. Returns a small dict summarizing what was written."""
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path  = out_dir / "results.csv"
    json_path = out_dir / "results.json"
    fasta_path = out_dir / "results.fasta"
    meta_path = out_dir / "metadata.json"
    summary_path = out_dir / "summary.md"

    n_csv   = write_csv(variants, csv_path)
    n_json  = write_json(variants, json_path)
    n_fasta = write_fasta(variants, fasta_path)

    per_source = [
        {
            "source": r.source,
            "status": r.status.value,
            "n_variants": len(r.variants),
            "elapsed_s": round(r.elapsed_s, 2),
            "message": r.message,
        }
        for r in per_source_results
    ]
    metadata = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "query": {
            "gene_symbols": query.gene_symbols,
            "species": query.species,
            "taxon_id": query.taxon_id,
            "sources": query.sources,
            "max_results_per_source": query.max_results_per_source,
            "include_sequence": query.include_sequence,
        },
        "per_source": per_source,
        "totals": {
            "variants": len(variants),
            "csv_rows": n_csv,
            "json_rows": n_json,
            "fasta_records": n_fasta,
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "notes": notes,
    }
    meta_path.write_text(json.dumps(metadata, indent=2))
    summary_path.write_text(_render_summary(variants, query, per_source, notes))

    return {
        "dir": str(out_dir),
        "csv": str(csv_path) if n_csv else "",
        "json": str(json_path) if n_json else "",
        "fasta": str(fasta_path) if n_fasta else "",
        "metadata": str(meta_path),
        "summary": str(summary_path),
        "n_variants": len(variants),
    }


def _render_summary(
    variants: list[ProteinVariant],
    query: SearchQuery,
    per_source: list[dict],
    notes: str,
) -> str:
    lines: list[str] = []
    lines.append(f"# Search Results — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("## Query")
    lines.append("")
    lines.append(f"- **Genes**: {', '.join(query.gene_symbols)}")
    lines.append(f"- **Species**: {query.species or '(any)'}")
    lines.append(f"- **Sources**: {', '.join(query.sources)}")
    lines.append(f"- **Max per source**: {query.max_results_per_source}")
    lines.append(f"- **Sequences fetched**: {'yes' if query.include_sequence else 'no'}")
    if notes:
        lines.append("")
        lines.append(f"_Notes_: {notes}")

    lines.append("")
    lines.append("## Per-source results")
    lines.append("")
    lines.append("| Source | Status | Variants | Elapsed (s) | Message |")
    lines.append("|--------|--------|----------|-------------|---------|")
    for s in per_source:
        lines.append(f"| {s['source']} | {s['status']} | {s['n_variants']} | {s['elapsed_s']} | {s['message']} |")

    lines.append("")
    lines.append(f"**Total variants**: {len(variants)}")

    # Per-gene × per-source table with length range — the Dong-style discrepancy view
    lines.append("")
    lines.append("## Variants by gene × source (length range in aa)")
    lines.append("")
    by_gs: dict[tuple[str, str], list[ProteinVariant]] = defaultdict(list)
    for v in variants:
        by_gs[(v.gene_symbol.upper(), v.source)].append(v)
    if by_gs:
        lines.append("| Gene | Source | Hits | Length range | Example accession |")
        lines.append("|------|--------|------|--------------|--------------------|")
        for (gene, source), vs in sorted(by_gs.items()):
            lens = [v.length_aa for v in vs if v.length_aa]
            rng = f"{min(lens)} – {max(lens)}" if lens else "n/a"
            lines.append(f"| {gene} | {source} | {len(vs)} | {rng} | {vs[0].accession} |")
    else:
        lines.append("_(no variants)_")

    # Spot Dong-style cross-source discrepancies in protein length
    lines.append("")
    lines.append("## Cross-database length discrepancies")
    lines.append("")
    lines.append(
        "Same gene reported with significantly different max protein lengths "
        "across sources — flags incomplete annotation, which in a ~2,700-residue "
        "gene is the commonest way a real paralog goes uncounted."
    )
    lines.append("")
    by_gene: dict[str, dict[str, int]] = defaultdict(dict)
    for (gene, source), vs in by_gs.items():
        lens = [v.length_aa for v in vs if v.length_aa]
        if lens:
            by_gene[gene][source] = max(lens)
    flagged = 0
    lines.append("| Gene | Per-source max length | Δ | Flag |")
    lines.append("|------|----------------------|---|------|")
    for gene, src_max in sorted(by_gene.items()):
        if len(src_max) < 2:
            continue
        vals = list(src_max.values())
        delta = max(vals) - min(vals)
        flag = "⚠️ check annotations" if delta >= 100 else "ok"
        if delta >= 100:
            flagged += 1
        cells = ", ".join(f"{s}={n}" for s, n in sorted(src_max.items()))
        lines.append(f"| {gene} | {cells} | {delta} | {flag} |")
    if flagged == 0:
        lines.append("")
        lines.append("_No significant length discrepancies (Δ ≥ 100 aa) detected._")

    # Genes that returned zero hits — biologically meaningful (a real loss)
    lines.append("")
    lines.append("## Genes with zero hits")
    lines.append("")
    requested = {g.upper() for g in query.gene_symbols}
    found = {v.gene_symbol.upper() for v in variants}
    missing = sorted(requested - found)
    if missing:
        lines.append("These query genes returned no hits from any source:")
        lines.append("")
        for g in missing:
            lines.append(f"- **{g}** — absent or unannotated in the searched species/databases")
        lines.append("")
        lines.append(
            "_Note_: absence may be biologically real — but at this stage it is only an absence from a *database*. Confirm against raw assembly sequence before calling a loss."
        )
    else:
        lines.append("_All requested genes returned at least one variant._")

    lines.append("")
    return "\n".join(lines) + "\n"
