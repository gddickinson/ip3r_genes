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
    "04b_results_architecture.md",  # the gene itself: exons, introns (S21)
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
#: Every entry names a figure a completed task committed. `s14_figures.py`
#: fails loudly on any that is missing rather than quietly shipping an
#: incomplete set. Renumber freely - the number here is the publication
#: number, and stale files in `manuscript/figures/` are deleted on every
#: build.
MAIN_FIGURES = [
    (1, "eukaryote_range", "s20_sweep/figures/range_by_phylum",
     "The ITPR family across the eukaryotes"),
    (2, "vertebrate_census", "genome_ledger/figures/ledger_by_class",
     "The three-paralog census across 309 vertebrate genomes"),
    (3, "phylogeny", "phylogeny/figures/tree_ml_rooted",
     "One family, three vertebrate paralogs, rooted on the ryanodine "
     "receptors"),
    (4, "paralogon", "duplication/figures/s16_paralogon",
     "Where the three vertebrate paralogs came from"),
    (5, "retention", "loss_dynamics/figures/s15_character_matrix",
     "Three paralogs, retained in every vertebrate genome searched"),
    (6, "constraint", "constraint/figures/s17_channel_profile",
     "The pore-forming half: what cannot change, and what can"),
    (7, "annotation", "annotation_audit/figures/s18_fig1_by_source",
     "How the family is recorded, and by whom"),
]

#: Extended Data figures. Several bundle more than one source panel file;
#: panels are lettered a, b, c ... in the order listed here.
EXTENDED_FIGURES = [
    (1, "outside_vertebrates",
     ["s23_scope/figures/copy_number",
      "s23_scope/figures/absence_at_genome",
      "s20_sweep/figures/plant_fungal_chase"],
     "Copy number and controlled absence outside the vertebrates"),
    (2, "vertebrate_sweep",
     ["genome_ledger/figures/ledger_status",
      "genome_ledger/figures/contiguity_confound",
      "genome_ledger/figures/copy_number"],
     "The vertebrate genomic sweep and its contiguity confounder"),
    (3, "methods",
     ["methods/figures/fig_s19_contiguity",
      "methods/figures/fig_s19_panel",
      "methods/figures/fig_s19_contribution",
      "methods/figures/fig_s19_drift"],
     "What the search methods were worth"),
    (4, "alignment",
     ["msa_v2/figures/msa_conservation",
      "msa_v2/figures/msa_identity_heatmap",
      "msa_v2/figures/msa_coverage"],
     "The representative alignment every downstream result stands on"),
    (5, "phylogeny_detail",
     ["phylogeny/figures/sister_au",
      "phylogeny/figures/support_profile",
      "phylogeny/figures/paralog_placement"],
     "The sister test, node support and paralog placement"),
    (6, "synteny",
     ["synteny/figures/synteny_paralogon",
      "synteny/figures/synteny_pair_classes",
      "synteny/figures/synteny_caller",
      "synteny/figures/synteny_clade_decay"],
     "Each paralog has its own genomic neighbourhood"),
    (7, "reconciliation",
     ["reconciliation/figures/recon_dated_backbone",
      "reconciliation/figures/recon_matrix",
      "reconciliation/figures/recon_losses",
      "reconciliation/figures/recon_cyclostome"],
     "Dating the duplications that made ITPR1, ITPR2 and ITPR3"),
    (8, "duplication",
     ["duplication/figures/s16_copy_number",
      "duplication/figures/s16_dcs",
      "duplication/figures/s16_blocks"],
     "The teleost genome duplication and the ancestral block"),
    (9, "gene_architecture",
     ["gene_architecture/figures/architecture_by_paralog",
      "gene_architecture/figures/intron_positions",
      "gene_architecture/figures/junction_quality",
      "gene_architecture/figures/fragments_and_duplicates"],
     "The gene: exon structure, shared introns, and where a fragmentary "
     "annotation stops"),
    (10, "retention_evidence",
     ["loss_dynamics/figures/s15_reconstruction",
      "loss_dynamics/figures/s15_integrity",
      "loss_dynamics/figures/s15_synteny_reach",
      "loss_counts/figures/sensitivity_matrix"],
     "The loss instrument, and what it takes to manufacture a loss"),
    (11, "selection",
     ["selection/figures/s9_omega_by_paralog",
      "selection/figures/s9_branch_contrast",
      "selection/figures/s9_dnds_saturation",
      "selection/figures/s9_bs_restarts"],
     "Selection across the three paralogs"),
    (12, "constraint_detail",
     ["constraint/figures/s17_elements",
      "constraint/figures/s17_functional_sites",
      "constraint/figures/s17_variant_classifier"],
     "Constraint by element, at the ligand site, and as a variant "
     "classifier"),
    (13, "ligand_site",
     ["ligand_site/figures/s22_fig1_modules",
      "ligand_site/figures/s22_fig2_shells",
      "ligand_site/figures/s22_fig3_omega",
      "ligand_site/figures/s22_fig4_lineage"],
     "The ligand core against the pore, and the lineages that lost the "
     "enzyme upstream"),
    (14, "structures",
     ["structures/figures/s11_panel",
      "structures/figures/s11_tm_calibration",
      "structures/figures/s11_plddt_domains",
      "structures/figures/s11_afdb_coverage"],
     "Predicted and experimental structures across the family"),
    (15, "annotation_quality",
     ["annotation_audit/figures/s18_fig2_calibration",
      "annotation_audit/figures/s18_fig3_family_vs_control",
      "annotation_audit/figures/s18_fig4_protein_side"],
     "The annotation audit against its own controls"),
    (16, "annotation_bugs",
     ["annotation_bugs/figures/exon_tracks",
      "annotation_bugs/figures/annotation_loss",
      "annotation_bugs/figures/validation",
      "expression/figures/s12_junctions"],
     "Two annotation failures validated at the exon and by RNA-seq"),
]

#: Supplementary figures - the alignments and structures the main figures
#: rest on. S24 drew these and filled the list; before it ran the list was
#: deliberately empty, because `s14_figures.py` exits non-zero on a missing
#: figure and a placeholder entry for a figure nobody had drawn would have
#: made the S14a build fail for a reason that was not a defect.
#:
#: One source stem each: unlike the Extended Data set, every supplementary
#: figure here is a single multi-panel file, so no letter suffix is appended.
SUPPLEMENTARY_FIGURES: list[tuple[int, str, list[str], str]] = [
    (1, "representative_alignment",
     ["supplementary/figures/SuppFig1_representative_alignment"],
     "The representative alignment, and the columns the tree actually saw"),
    (2, "labelled_positions",
     ["supplementary/figures/SuppFig2_labelled_positions"],
     "The ligand core and the pore module at residue resolution across the "
     "three paralogues"),
    (3, "paralog_alignments",
     ["supplementary/figures/SuppFig3_paralog_alignments"],
     "The within-paralogue alignments the constraint map is computed on"),
    (4, "codon_alignment",
     ["supplementary/figures/SuppFig4_codon_alignment"],
     "The trimmed codon alignment behind every omega estimate"),
    (5, "constraint_on_channel",
     ["supplementary/figures/SuppFig5_constraint_on_channel"],
     "The constraint map painted onto the channel"),
    (6, "variants_on_structure",
     ["supplementary/figures/SuppFig6_variants_on_structure"],
     "Every labelled variant, and the per-element enrichment test"),
]

# ------------------------------------------------------------------ deposit

#: Directories deposited wholesale (relative to results/): one per completed
#: ledger task, plus the five census editions. A directory named here that
#: does not exist is reported by `s14_deposit.py` rather than skipped, so a
#: renamed results tree cannot silently shrink the deposit.
DEPOSIT_DIRS = [
    "annotation_audit", "annotation_bugs", "benchmark_controls",
    "census_v2", "census_v3", "census_v4", "census_v5", "census_v6",
    "constraint", "duplication", "expression", "gene_architecture",
    "genome_ledger", "hmm_sweep", "ligand_site", "loss_counts",
    "loss_dynamics", "methods", "msa_v2", "phylogeny",
    "reconciliation", "s0_baseline", "s20_sweep", "s23_baits", "s23_scope",
    "s5_baits", "selection", "structures", "supplementary", "synteny",
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
    # Every command below was run to produce the data it regenerates, and
    # every one takes its scope from the committed manifest beside it, so a
    # reader reproduces the same set rather than a similar one. They need the
    # analysis conda environment (see results/toolchain_manifest.txt); the
    # manuscript package itself builds on a bare Python.
    ("Genome assemblies and their annotations, 503 genomes",
     "624 GB", "NCBI Datasets",
     "results/genome_manifest.tsv, results/s23_scope/genome_manifest_s23.tsv",
     "python scripts/fetch_genomes.py"),
    ("UniProt reference proteomes, 7,691 proteomes, and the concatenated "
     "search databases built from them",
     "60 GB", "UniProt release FTP",
     "results/hmm_sweep/proteome_manifest.tsv, "
     "results/s20_sweep/proteome_manifest_*.tsv",
     "python scripts/s3_fetch_proteomes.py --group vertebrata && "
     "python scripts/s20_fetch.py"),
    ("Per-genome vertebrate sweep output: miniprot GFF, rescue output and "
     "summary.json for each of the 309 genomes",
     "939 MB", "generated",
     "results/genome_ledger/genome_ledger.tsv",
     "python scripts/s5_run_sweep.py"),
    ("Per-genome non-vertebrate sweep output for each of the 194 genomes",
     "726 MB", "generated",
     "results/s23_scope/copy_number_ledger.tsv",
     "python scripts/s23_run_sweep.py"),
    ("HMMER output: hmmsearch domain tables and jackhmmer logs for every "
     "profile sweep and convergence run",
     "1.3 GB", "generated",
     "results/hmm_sweep/sweep_stats_*.json, "
     "results/s20_sweep/jackhmmer_convergence_s20.tsv",
     "python scripts/s3_run_sweep.py"),
    ("Structure files: RCSB mmCIF downloads, AlphaFold DB models, the "
     "CA-trace panel and the TM-align pair cache",
     "4.0 GB", "RCSB PDB and AlphaFold DB",
     "results/structures/structure_manifest.tsv",
     "python scripts/s11_run.py --only panel"),
    ("Archived API responses: InterPro pages, UniProt records, NCBI "
     "datasets dumps and taxonomy lookups",
     "1.9 GB", "generated (archived so every parse re-runs offline)",
     "each task's stats json records what it fetched",
     "re-run the owning task; every fetch is cached and resumable"),
    ("RNA-seq: streamed SRA reads, HISAT2 indexes and per-run counts",
     "15 MB committed-side; reads are streamed and not retained",
     "NCBI SRA",
     "results/expression/runs_selected.tsv",
     "python scripts/s12_run.py --only quantify"),
    ("BLAST databases built from the swept proteomes and genomes",
     "77 MB", "generated",
     "results/s5_baits/baits.faa, results/s23_baits/baits.faa",
     "rebuilt automatically by the sweep drivers"),
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
