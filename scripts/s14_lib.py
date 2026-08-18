"""S14 shared helpers — manuscript assembly, figures, deposit manifest, claims.

Used by `s14_figures.py`, `s14_deposit.py`, `s14_claims.py` and the driver
`s14_assemble.py`. Stdlib only: the manuscript package must build on a bare
Python so a reviewer can regenerate it without the analysis toolchain.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

# ---------------------------------------------------------------- page size
#: Printed width of a full-width figure, in inches. A4 at 2.0 cm margins gives
#: a 17.0 cm text block. `figstyle` imports these so figures are drawn at the
#: size they are placed at.
W_FULL = 6.7
W_HALF = 3.25
#: Tallest a figure may be before the page height limits it.
H_MAX = 8.6

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
MS = ROOT / "manuscript"
MS_FIGS = MS / "figures"

# ---------------------------------------------------------------- manuscript

#: Section files, in the order they are stitched into `manuscript.md`.
#: The manuscript does not exist until S14 writes it; `s14_assemble.py`
#: exits non-zero while any of these is missing, which is the intended
#: behaviour for a package that is not yet written.
SECTION_ORDER = [
    "00_frontmatter.md",
    "01_main.md",
    "02_results_range.md",       # the family across the eukaryotes
    "03_results_census.md",      # the vertebrate three-paralog census
    "04_results_origin.md",      # where ITPR1/2/3 came from (2R, 3R)
    "05_results_fates.md",       # retention, loss, decay
    "06_results_machine.md",     # constraint, structure, clinical variants
    "07_results_annotation.md",  # how badly the family is recorded
    "08_discussion.md",
    "09_methods_search.md",
    "10_methods_analysis.md",
    "11_figure_legends.md",
    "12_extended_data.md",
    "13_supplementary.md",
    "14_data_availability.md",
    "15_references.md",
]

SECTION_FRONTMATTER = "00_frontmatter.md"
SECTION_FIGURE_LEGENDS = "11_figure_legends.md"
SECTION_EXTENDED_DATA = "12_extended_data.md"
SECTION_SUPPLEMENTARY = "13_supplementary.md"

# ------------------------------------------------------------------ figures

#: (number, slug, source path relative to results/, one-line caption stub).
#: This is the *plan*: the source paths are where the roadmap's tasks write
#: their figures, and `s14_figures.py` fails loudly on any that is missing
#: rather than quietly shipping an incomplete set. Renumber freely — the
#: number here is the publication number, and stale files in
#: `manuscript/figures/` are deleted on every build.
MAIN_FIGURES = [
    (1, "eukaryote_range", "hmm_sweep/s20/s20_presence_figure",
     "The ITPR family across the eukaryotes"),
    (2, "ledger_by_class", "figures/genome_ledger_by_class",
     "The three-paralog census across the vertebrate genome scope"),
    (3, "phylogeny", "phylogeny/figures/tree_ml_rooted",
     "One family, three vertebrate paralogs, rooted on the ryanodine "
     "receptors"),
    (4, "duplication", "duplication/figures/duplication_history",
     "Where the three vertebrate paralogs came from"),
    (5, "loss_dynamics", "loss_dynamics/figures/loss_dynamics",
     "Three paralogs, three fates"),
    (6, "constraint", "constraint/figures/constraint",
     "One machine: constraint, structure and the clinical variants"),
    (7, "annotation_audit", "annotation_audit/figures/annotation_audit",
     "How the family is recorded"),
]

#: Extended Data figures. Several bundle more than one source panel file.
EXTENDED_FIGURES = [
    (1, "invertebrates", ["figures/s23_copy_number"],
     "ITPR copy number across the invertebrate and protist genomes"),
    (2, "synteny", ["synteny/figures/synteny_tracks",
                    "synteny/figures/jaccard_heatmap",
                    "synteny/figures/pairclass_box"],
     "Each paralog has its own genomic neighbourhood"),
    (3, "reconciliation", ["reconciliation/figures/reconciliation"],
     "Dating the duplications that made ITPR1, ITPR2 and ITPR3"),
    (4, "selection", ["selection/figures/selection"],
     "Selection across the three paralogs"),
    (5, "structures", ["structures/figures/structures"],
     "Predicted and experimental structures across the family"),
    (6, "architecture", ["architecture/figures/architecture"],
     "Gene architecture and the shared intron map"),
    (7, "annotation_bugs", ["annotation_validation/figures/annotation_bugs"],
     "Annotation errors validated at the molecular level"),
    (8, "ledger_heatmap", ["figures/genome_ledger_heatmap"],
     "The full per-genome evidence matrix"),
    (9, "expression", ["expression/figures/expression"],
     "Transcription and splicing evidence"),
    (10, "alignment", ["msa_v2/figures/conservation_trimmed",
                       "msa_v2/figures/identity_heatmap_covered"],
     "The representative alignment"),
    (11, "methods", ["methods/figures/methods"],
     "What the search methods were worth"),
]

#: Supplementary figures — the alignments and structures the main figures
#: rest on. Nothing in them is re-aligned or re-rendered from new data, so a
#: supplementary panel that disagreed with its main figure would be a bug.
SUPPLEMENTARY_FIGURES = [
    (1, "alignment_phylogeny", ["alignments/figures/supp_aln_msa_v2"],
     "The alignment behind the phylogeny"),
    (2, "alignment_pore", ["alignments/figures/supp_aln_pore"],
     "The pore module and the IP3-binding core, residue by residue"),
    (3, "alignment_constraint", ["alignments/figures/supp_aln_deep"],
     "The alignments the constraint map is computed on"),
    (4, "alignment_codon", ["alignments/figures/supp_aln_codon"],
     "The codon alignment behind every omega estimate"),
    (5, "structure_constraint", ["alignments/figures/supp_struct_constraint"],
     "The constraint map on the channel"),
    (6, "structure_findings", ["alignments/figures/supp_struct_findings"],
     "What the structures show"),
]

# ------------------------------------------------------------------ deposit

#: Directories deposited wholesale (relative to results/).
DEPOSIT_DIRS = [
    "alignments", "annotation_audit", "annotation_validation", "architecture",
    "benchmark_controls", "census_v2", "census_v3", "census_v4", "census_v5",
    "constraint", "duplication", "expression", "figures", "hmm_sweep",
    "loss_dynamics", "methods", "msa_v2", "phylogeny", "reconciliation",
    "s23_baits", "s5_baits", "selection", "structures", "synteny",
]

#: Individual files at the results/ root that are deposited.
DEPOSIT_ROOT_GLOBS = ["*.tsv", "*.csv", "*.md", "*.txt", "*.faa", "*.json"]

#: Extensions never deposited (working files, caches, binaries).
DEPOSIT_SKIP_SUFFIXES = {".pyc", ".log"}

#: Exact filenames never deposited.
DEPOSIT_SKIP_NAMES = {".DS_Store", "Thumbs.db",
                      ".pdf_build.md", ".pdf_preamble.tex"}

#: Bulk data classes deliberately excluded, with how to regenerate each.
#: Filled in by the tasks that create each class — an entry here without a
#: working regeneration command is worse than no entry, so add one only when
#: the command has been run.
BULK_EXCLUSIONS: list[tuple[str, str, str, str, str]] = [
    # (what, size, source, manifest committed in the repo, command)
]


def deposit_group(rel_parts: tuple[str, ...]) -> str:
    """Group label for the deposit manifest: one row group per source tree."""
    if rel_parts[0] == "results" and len(rel_parts) > 1:
        return f"results/{rel_parts[1]}" if len(rel_parts) > 2 else "results"
    return rel_parts[0] if len(rel_parts) > 1 else "root"


def png_size_inches(path: Path) -> tuple[float, float]:
    """(width, height) in inches from a PNG's IHDR and pHYs chunks.

    Stdlib only, so the manuscript package still builds without matplotlib.
    Falls back to 400 dpi (the project's savefig default) when a file carries
    no physical-size chunk.
    """
    import struct

    with path.open("rb") as fh:
        if fh.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"{path} is not a PNG")
        px_w = px_h = 0
        dpi = 400.0
        while True:
            head = fh.read(8)
            if len(head) < 8:
                break
            length, kind = struct.unpack(">I4s", head)
            data = fh.read(length)
            fh.read(4)                                   # CRC
            if kind == b"IHDR":
                px_w, px_h = struct.unpack(">II", data[:8])
            elif kind == b"pHYs" and len(data) >= 9 and data[8] == 1:
                ppm_x = struct.unpack(">I", data[:4])[0]  # pixels per metre
                if ppm_x:
                    dpi = ppm_x * 0.0254
            elif kind in (b"IDAT", b"IEND"):
                break
    if not px_w:
        raise ValueError(f"{path}: no IHDR")
    return px_w / dpi, px_h / dpi


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_json(path: Path):
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def size_str(n_bytes: int) -> str:
    """Compact human size (1 decimal above KB)."""
    x = float(n_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if x < 1024 or unit == "GB":
            return f"{x:.0f} {unit}" if unit == "B" else f"{x:.1f} {unit}"
        x /= 1024.0
    return f"{x:.1f} GB"
