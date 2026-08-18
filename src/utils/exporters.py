"""Export ProteinVariant lists to FASTA / CSV / JSON."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from ..core.models import ProteinVariant


def _wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def write_fasta(variants: list[ProteinVariant], path: str | Path) -> int:
    """Write variants as multi-FASTA. Returns the number of records written.

    Variants without a sequence are skipped (we don't fabricate placeholder
    sequences — exporter is for downstream tools that need real residues).
    """
    p = Path(path)
    written = 0
    with p.open("w") as f:
        for v in variants:
            if not v.sequence:
                continue
            header = f">{v.source}|{v.accession}|{v.gene_symbol}|{v.species}|len={v.length_aa or len(v.sequence)}"
            f.write(header + "\n")
            f.write(_wrap(v.sequence) + "\n")
            written += 1
    return written


def write_csv(variants: list[ProteinVariant], path: str | Path) -> int:
    p = Path(path)
    if not variants:
        p.write_text("")
        return 0
    cols = ["source", "gene_symbol", "accession", "species", "taxon_id",
            "length_aa", "transcript_id", "gene_id", "description", "url"]
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for v in variants:
            row = {c: getattr(v, c, "") for c in cols}
            w.writerow(row)
    return len(variants)


def write_json(variants: list[ProteinVariant], path: str | Path) -> int:
    p = Path(path)
    payload = [asdict(v) for v in variants]
    p.write_text(json.dumps(payload, indent=2, default=str))
    return len(payload)
