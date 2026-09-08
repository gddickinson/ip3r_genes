"""S10's committed tables — written once, read by the report and the figures.

Nothing downstream re-parses a GFF, re-runs blastp or re-reads the genome
(D13). If a number is in the report it is in one of these files, and if it is
in a figure it is in the same file the report read it from.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

EXON_COLS = ["case_id", "accession", "cell", "exon_index", "start", "end",
             "length", "q_start", "q_end", "aln_identity", "phase", "class",
             "annot_gene", "annot_gene_id", "annot_biotype", "annot_pseudo",
             "annot_cds_overlap_bp", "annot_cds_overlap_frac"]

INTRON_COLS = ["case_id", "accession", "cell", "intron_index", "start", "end",
               "length", "class", "left_models", "right_models", "donor",
               "acceptor", "splice_class"]

MODEL_COLS = ["case_id", "accession", "cell", "gene_id", "name", "locus_tag",
              "biotype", "pseudo", "start", "end", "strand", "span",
              "n_cds_blocks", "cds_bp", "n_transcripts", "protein_ids",
              "overlap_aligned_cds_bp", "frac_aligned_cds",
              "n_aligned_exons_covered", "q_start", "q_end", "description"]

BLOCK_COLS = ["case_id", "accession", "cell", "n_aligned_exons",
              "n_aligned_cds_blocks", "aligned_cds_bp",
              "n_annotated_gene_models", "n_annotated_cds_blocks",
              "annotated_cds_bp_on_gene", "frac_aligned_cds_annotated",
              "frac_aligned_cds_translated", "n_uncovered_blocks",
              "uncovered_bp", "n_untranslated_blocks", "untranslated_bp"]

NEIGHBOUR_COLS = ["case_id", "accession", "cell", "side", "gene_id", "name",
                  "biotype", "pseudo", "start", "end", "strand",
                  "distance_bp", "description"]


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def stamp(rows: list[dict], case: dict) -> list[dict]:
    """Put the case identity on every row of every table.

    Both cases' rows live in one file per table. A per-case file would be
    tidier to write and worse to read: a reader comparing the two cases would
    have to join, and the report's own tables would not be checkable against a
    single source.
    """
    ident = {"case_id": case.get("case_id", ""),
             "accession": case.get("accession", ""),
             "cell": case.get("cell", "")}
    return [{**ident, **r} for r in rows]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_stats(out: Path, payload: dict, tables: list[Path]) -> None:
    """`annotation_bugs_stats.json` — the run's parameters and table hashes.

    The SHA-256 of every committed table goes in, as in S6, S8 and S9, so a
    re-run that drifts is visible in the data rather than only in a diff.
    """
    payload = dict(payload)
    payload["tables"] = {
        p.name: {"sha256": sha256(p), "bytes": p.stat().st_size}
        for p in sorted(tables) if p.exists()}
    (out / "annotation_bugs_stats.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n")
