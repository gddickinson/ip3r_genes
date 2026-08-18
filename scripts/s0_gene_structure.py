#!/usr/bin/env python3
"""S0 step 4 (continued) — turn the gene-architecture ``[lit]`` claims into ``[db]``.

`docs/ip3r_background.md` asserts from textbook knowledge that each human ITPR
is "a ~58-60 exon gene spanning hundreds of kb", and gives a cytogenetic locus
per paralog. Both are checkable against live databases, so S0 checks them
rather than citing them:

  * exon count + genomic span   <- Ensembl REST ``lookup/symbol?expand=1``
  * cytogenetic band            <- Ensembl REST ``overlap/region/.../band``

Writes ``results/s0_baseline/gene_structure.tsv``.

Note on Ensembl (S0 smoke test): ``/xrefs/symbol/homo_sapiens/...`` — the
endpoint ``src/databases/ensembl.py`` resolves symbols with — stalls
indefinitely, while ``/lookup/symbol/homo_sapiens/...`` answers normally. This
script deliberately uses ``lookup/symbol``, which is also the fix recorded for
the client in the S0 session log.

Usage:  python scripts/s0_gene_structure.py [--out results/s0_baseline]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

ENSEMBL = "https://rest.ensembl.org"
TIMEOUT_S = 120
RETRIES = 3

GENES = [
    ("ITPR1", "homo_sapiens", "3p26.1"),      # third field: the [lit] claim
    ("ITPR2", "homo_sapiens", "12p11.23"),
    ("ITPR3", "homo_sapiens", "6p21.31"),
]


def _get(path: str, params: dict | None = None):
    params = dict(params or {})
    params["content-type"] = "application/json"
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.get(f"{ENSEMBL}{path}", params=params,
                             timeout=TIMEOUT_S,
                             headers={"Accept": "application/json"})
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except Exception as exc:                       # noqa: BLE001
            last = exc
            print(f"      retry {attempt + 1}/{RETRIES}: {exc}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"Ensembl failed after {RETRIES} tries: {path} ({last})")


def cytoband(species: str, chrom: str, start: int, end: int) -> str:
    """Cytogenetic band(s) the gene's span overlaps, e.g. '3p26.1'."""
    payload = _get(f"/overlap/region/{species}/{chrom}:{start}-{end}",
                   {"feature": "band"})
    if not payload:
        return ""
    bands = sorted({b.get("id", "") for b in payload if b.get("id")})
    arm = chrom
    return ";".join(f"{arm}{b}" for b in bands)


def gene_row(symbol: str, species: str, lit_locus: str) -> list:
    print(f"  {symbol} ...")
    gene = _get(f"/lookup/symbol/{species}/{symbol}", {"expand": "1"})
    if gene is None:
        raise RuntimeError(f"Ensembl has no gene {symbol} in {species}")
    chrom = gene["seq_region_name"]
    start, end = int(gene["start"]), int(gene["end"])
    span = end - start + 1

    transcripts = gene.get("Transcript", [])
    canonical = None
    for t in transcripts:
        if t.get("is_canonical"):
            canonical = t
            break
    if canonical is None and transcripts:
        # fall back to the protein-coding transcript with the most exons
        coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
        canonical = max(coding or transcripts,
                        key=lambda t: len(t.get("Exon", [])))

    n_exon_canonical = len(canonical.get("Exon", [])) if canonical else 0
    coding_tx = [t for t in transcripts if t.get("biotype") == "protein_coding"]
    exon_counts = [len(t.get("Exon", [])) for t in coding_tx] or [0]
    protein_len = ""
    if canonical and canonical.get("Translation"):
        protein_len = canonical["Translation"].get("length", "")

    band = cytoband(species, chrom, start, end)
    print(f"      {gene['id']}  chr{chrom}:{start:,}-{end:,}  span={span:,} bp  "
          f"exons(canonical)={n_exon_canonical}  band={band}  "
          f"lit_locus={lit_locus}")
    return [symbol, gene["id"], chrom, start, end, span,
            canonical.get("id", "") if canonical else "",
            n_exon_canonical, min(exon_counts), max(exon_counts),
            len(coding_tx), len(transcripts), protein_len, band, lit_locus,
            int(band.startswith(lit_locus.split(".")[0])) if band else 0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="results/s0_baseline")
    args = ap.parse_args()
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)

    print("S0 gene structure from Ensembl (lookup/symbol, expand=1)")
    rows = [gene_row(*g) for g in GENES]

    header = ["symbol", "ensembl_gene", "chromosome", "start", "end",
              "genomic_span_bp", "canonical_transcript", "n_exons_canonical",
              "n_exons_min_coding", "n_exons_max_coding", "n_coding_transcripts",
              "n_transcripts", "canonical_protein_aa", "cytoband_db",
              "cytoband_lit", "band_arm_matches_lit"]
    path = out / "gene_structure.tsv"
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")
    print(f"  wrote {path.relative_to(ROOT)}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
