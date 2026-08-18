"""End-to-end analysis pipeline.

`analyse(variants, …)` is the single entry point both the GUI Analyze tab
and the CLI `--analyze` flag use. It:

    1. Filters variants to those with a real sequence (skips short isoforms
       and DB hits that came back metadata-only).
    2. Builds the MSA (star or MAFFT).
    3. Computes the identity matrix.
    4. Builds a neighbor-joining tree.
    5. Clusters at the user-chosen identity threshold; flags novel
       candidates.
    6. (Optional) Calls mutations against a chosen reference.
    7. Scores conservation per column; locates variable regions.

The whole result is returned as `AnalysisResult` (also written to disk in
the results bundle as `analysis/`).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..core.models import ProteinVariant
from .alignment import AlignmentRow, progressive_msa
from .clusters import ClusterReport, cluster_variants
from .distance import DistanceTable, identity_matrix
from .evolution import VariableRegion, conservation_per_column, variable_regions
from .mutations import MutationCall, call_mutations
from .tree import build_nj_tree, tree_to_ascii, tree_to_newick


@dataclass
class AnalysisResult:
    n_input: int
    n_analyzed: int                       # variants that had a sequence
    msa: list[AlignmentRow] = field(default_factory=list)
    distances: Optional[DistanceTable] = None
    newick: str = ""
    ascii_tree: str = ""
    clusters: Optional[ClusterReport] = None
    mutations: list[MutationCall] = field(default_factory=list)
    reference_label: str = ""
    conservation: list[float] = field(default_factory=list)
    variable_regions: list[VariableRegion] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)

    def text_summary(self) -> str:
        lines = [
            f"Analysis summary: {self.n_analyzed}/{self.n_input} variants had sequences.",
        ]
        if self.skipped:
            lines.append(f"  Skipped (no sequence): {', '.join(self.skipped[:5])}"
                         + ("…" if len(self.skipped) > 5 else ""))
        if self.distances and self.distances.labels:
            n = len(self.distances.labels)
            mean_id = (
                sum(self.distances.identity[i][j]
                    for i in range(n) for j in range(i + 1, n))
                / max(1, n * (n - 1) // 2)
            )
            lines.append(f"  Mean pairwise identity: {mean_id:.1%}")
        if self.clusters:
            lines.append("")
            lines.append(self.clusters.summary())
        if self.variable_regions:
            lines.append("")
            lines.append(f"Variable regions (length ≥ 5, mean conservation < threshold): "
                         f"{len(self.variable_regions)}")
            for r in self.variable_regions[:10]:
                lines.append(f"  • cols {r.start}-{r.end} (width {r.width()}, "
                             f"conservation {r.mean_conservation:.2f})")
        if self.mutations:
            lines.append("")
            lines.append(f"Mutation calls vs reference '{self.reference_label}': "
                         f"{len(self.mutations)} positions")
            for m in self.mutations[:15]:
                lines.append(f"  • {m.query_label}: {m.short()}")
            if len(self.mutations) > 15:
                lines.append(f"  …and {len(self.mutations) - 15} more.")
        if self.ascii_tree:
            lines.append("")
            lines.append("Neighbor-joining tree (ASCII):")
            lines.append(self.ascii_tree.rstrip())
        return "\n".join(lines)


def _label_for(v: ProteinVariant) -> str:
    """Concise, unique label used in alignments and trees."""
    species = (v.species or "?").split()[0][:8]
    return f"{v.source[:3]}_{v.gene_symbol}_{species}_{v.accession[:14]}"


def _pick_reference(variants: list[ProteinVariant], requested: str = "") -> Optional[ProteinVariant]:
    if not variants:
        return None
    if requested:
        for v in variants:
            if v.accession == requested or _label_for(v) == requested:
                return v
    # Default: longest sequence — usually the canonical isoform.
    return max(variants, key=lambda v: len(v.sequence or ""))


def analyse(
    variants: Iterable[ProteinVariant],
    identity_threshold: float = 0.5,
    conservation_threshold: float = 0.6,
    reference_accession: str = "",
    use_mafft: bool = False,
    max_variants: int = 60,
) -> AnalysisResult:
    variants = list(variants)
    with_seq = [v for v in variants if v.sequence]
    skipped = [_label_for(v) for v in variants if not v.sequence]

    if not with_seq:
        return AnalysisResult(n_input=len(variants), n_analyzed=0, skipped=skipped)

    if len(with_seq) > max_variants:
        # Trim to the longest N — short isoform fragments are noise for
        # tree-building.
        with_seq = sorted(with_seq, key=lambda v: -len(v.sequence))[:max_variants]

    # _label_for is concise but not guaranteed unique (truncation can
    # collide across isoforms of the same protein). Append a suffix when
    # needed so downstream Bio.Phylo doesn't reject the matrix.
    raw_labels = [_label_for(v) for v in with_seq]
    seen: dict[str, int] = {}
    unique_labels: list[str] = []
    for lbl in raw_labels:
        n = seen.get(lbl, 0)
        seen[lbl] = n + 1
        unique_labels.append(lbl if n == 0 else f"{lbl}_{n+1}")
    labelled = list(zip(unique_labels, [v.sequence for v in with_seq]))
    label_by_variant = {id(v): unique_labels[i] for i, v in enumerate(with_seq)}
    rows = progressive_msa(labelled, use_mafft=use_mafft)
    distances = identity_matrix(rows)
    clusters = cluster_variants(distances, identity_threshold=identity_threshold)
    tree = build_nj_tree(distances, method="nj")
    newick = tree_to_newick(tree)
    ascii_tree = tree_to_ascii(tree)
    conservation = conservation_per_column(rows)
    var_regions = variable_regions(conservation, threshold=conservation_threshold)

    reference = _pick_reference(with_seq, reference_accession)
    mutations: list[MutationCall] = []
    ref_label_used = ""
    if reference:
        ref_label_used = label_by_variant.get(id(reference), _label_for(reference))
        for v in with_seq:
            if v is reference:
                continue
            v_label = label_by_variant.get(id(v), _label_for(v))
            mutations.extend(call_mutations(ref_label_used, reference.sequence, v_label, v.sequence))

    return AnalysisResult(
        n_input=len(variants),
        n_analyzed=len(with_seq),
        msa=rows,
        distances=distances,
        newick=newick,
        ascii_tree=ascii_tree,
        clusters=clusters,
        mutations=mutations,
        reference_label=ref_label_used,
        conservation=conservation,
        variable_regions=var_regions,
        skipped=skipped,
    )


def write_analysis(out_dir: Path, result: AnalysisResult) -> dict:
    """Persist analysis artefacts inside a results bundle's `analysis/` subdir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}

    if result.msa:
        msa_path = out_dir / "alignment.fasta"
        msa_path.write_text("".join(f">{r.label}\n{r.aligned}\n" for r in result.msa))
        written["alignment_fasta"] = str(msa_path)

    if result.newick:
        nwk_path = out_dir / "tree.newick"
        nwk_path.write_text(result.newick + "\n")
        written["tree_newick"] = str(nwk_path)
        ascii_path = out_dir / "tree.txt"
        ascii_path.write_text(result.ascii_tree)
        written["tree_ascii"] = str(ascii_path)

    if result.distances:
        m_path = out_dir / "identity_matrix.tsv"
        with m_path.open("w") as f:
            f.write("\t" + "\t".join(result.distances.labels) + "\n")
            for i, label in enumerate(result.distances.labels):
                cells = "\t".join(f"{result.distances.identity[i][j]:.3f}"
                                  for j in range(len(result.distances.labels)))
                f.write(f"{label}\t{cells}\n")
        written["identity_matrix"] = str(m_path)

    if result.mutations:
        mut_path = out_dir / "mutations.tsv"
        with mut_path.open("w") as f:
            f.write("query\tposition\tref_aa\tquery_aa\tkind\tshort\n")
            for m in result.mutations:
                f.write(f"{m.query_label}\t{m.ref_pos}\t{m.ref_aa}\t{m.query_aa}\t{m.kind}\t{m.short()}\n")
        written["mutations"] = str(mut_path)

    summary_path = out_dir / "summary.txt"
    summary_path.write_text(result.text_summary())
    written["summary"] = str(summary_path)

    if result.conservation:
        cons_path = out_dir / "conservation.tsv"
        with cons_path.open("w") as f:
            f.write("column\tconservation\n")
            for i, c in enumerate(result.conservation, 1):
                f.write(f"{i}\t{c:.3f}\n")
        written["conservation"] = str(cons_path)

    manifest = {
        "n_input": result.n_input,
        "n_analyzed": result.n_analyzed,
        "reference_label": result.reference_label,
        "n_mutations": len(result.mutations),
        "n_clusters": len(result.clusters.clusters) if result.clusters else 0,
        "n_novel_candidates": len(result.clusters.novel_candidates) if result.clusters else 0,
        "n_variable_regions": len(result.variable_regions),
        "files": written,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    written["manifest"] = str(out_dir / "manifest.json")
    return written
