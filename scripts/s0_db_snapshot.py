#!/usr/bin/env python3
"""S0 step 4 — re-derive every ``[db]`` number in ``docs/ip3r_background.md``.

The background document carries a database snapshot taken by the setup
session. Nothing may be quoted from it until it has been re-derived by this
project's own run, so this script re-queries every one of those numbers and
writes them to committed tables:

    results/s0_baseline/interpro_signature_counts.tsv   section 4, table 1
    results/s0_baseline/interpro_taxonomy_counts.tsv    section 4, table 2
    results/s0_baseline/reference_proteins.tsv          sections 2 + 3
    results/s0_baseline/pfam_architecture.tsv           section 1 + section 3
    results/s0_baseline/zebrafish_itpr.tsv              section 2
    results/s0_baseline/db_snapshot_meta.json           date, endpoints, versions

`scripts/s0_report.py` renders the S0 report purely from these files, so the
report and the data cannot drift (Decisions D13).

Usage:  python scripts/s0_db_snapshot.py [--out results/s0_baseline]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
SLEEP_S = 0.34          # EBI politeness
TIMEOUT_S = 60
RETRIES = 4

# ---------------------------------------------------------------- what to ask

# Section 4, table 1: protein count per Pfam signature.
SIGNATURES = [
    ("PF08709", "Ins145_P3_rec", "the IP3-binding core (beta-trefoil)"),
    ("PF01365", "RYDR_ITPR (RIH)", "RyR and IP3R homology domain"),
    ("PF08454", "RIH_assoc", "RIH-associated"),
    ("PF02815", "MIR", "also in O-mannosyltransferases - not family-specific"),
    ("PF00520", "Ion_trans", "generic voltage-gated-channel pore domain"),
]

# Section 4, table 2: taxonomic distribution. Counted for the two signatures
# the background document tabulates.
TAXA = [
    ("Metazoa", 33208),
    ("SAR (stramenopiles/alveolates/rhizaria)", 2698737),
    ("Discoba (incl. kinetoplastids)", 2611352),
    ("Fungi", 4751),
    ("Viridiplantae", 33090),
    ("Amoebozoa", 554915),
    ("Bacteria", 2),
    ("Archaea", 2157),
    ("Arabidopsis thaliana", 3702),
    ("Saccharomyces cerevisiae", 559292),
    ("Paramecium tetraurelia", 5888),
    ("Trypanosoma brucei", 5691),
    ("Chlamydomonas reinhardtii", 3055),
    ("Dictyostelium discoideum", 44689),
]
TAXONOMY_SIGNATURES = ["PF08709", "PF01365"]

# Sections 2 and 3: the reference panel and the sister family.
REFERENCE_PROTEINS = [
    ("ITPR1", "Q14643", "ITPR"),
    ("ITPR2", "Q14571", "ITPR"),
    ("ITPR3", "Q14573", "ITPR"),
    ("RYR1", "P21817", "RYR"),
    ("RYR2", "Q92736", "RYR"),
    ("RYR3", "Q15413", "RYR"),
]

# Section 3: the RyR-specific Pfams the background document names, plus the
# four it says RYR1 shares with ITPR. Every one is checked per protein.
ARCHITECTURE_PFAMS = [
    ("PF02815", "MIR"),
    ("PF08709", "Ins145_P3_rec"),
    ("PF01365", "RYDR_ITPR"),
    ("PF08454", "RIH_assoc"),
    ("PF00520", "Ion_trans"),
    ("PF02026", "RyR"),
    ("PF06459", "RyR_TM4-6"),
    ("PF21119", "RyR_jsol"),
    ("PF00622", "SPRY"),
]

ZEBRAFISH_QUERY = "taxonomy_id:7955 AND xref:pfam-PF08709"


# ------------------------------------------------------------------ transport

def _get(url: str, params: dict | None = None) -> dict | None:
    """GET with retries. Returns None on a 204/404 (InterPro's 'no members')."""
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT_S,
                             headers={"Accept": "application/json"})
            if r.status_code in (204, 404):
                return None
            r.raise_for_status()
            time.sleep(SLEEP_S)
            return r.json()
        except Exception as exc:                       # noqa: BLE001
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET failed after {RETRIES} tries: {url} ({last})")


def interpro_count(pfam: str, taxon_id: int | None = None) -> int:
    """Protein count for a Pfam signature, optionally restricted to a taxon.

    InterPro reports the size of a filtered set in the ``count`` field of any
    listing page, so a page_size=1 request is the cheapest way to read it.
    A signature with no members in a taxon returns 204/404, i.e. a true zero.
    """
    url = f"{INTERPRO}/protein/uniprot/entry/pfam/{pfam}"
    if taxon_id is not None:
        url += f"/taxonomy/uniprot/{taxon_id}"
    payload = _get(url, {"page_size": 1})
    if payload is None:
        return 0
    return int(payload.get("count", 0))


def uniprot_entry(acc: str) -> dict:
    payload = _get(f"{UNIPROT}/{acc}.json")
    if payload is None:
        raise RuntimeError(f"UniProt returned nothing for {acc}")
    return payload


def uniprot_search(query: str, fields: str, size: int = 100) -> list[dict]:
    payload = _get(f"{UNIPROT}/search",
                   {"query": query, "fields": fields, "size": size,
                    "format": "json"})
    return (payload or {}).get("results", [])


# -------------------------------------------------------------------- helpers

def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if c is None else str(c) for c in row) + "\n")
    print(f"  wrote {path.relative_to(ROOT)}  ({len(rows)} rows)")


def entry_pfams(entry: dict) -> dict[str, int]:
    """Pfam accession -> number of cross-references, from a UniProt entry."""
    counts: dict[str, int] = {}
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "Pfam":
            continue
        acc = xref.get("id", "")
        n = 1
        for prop in xref.get("properties", []):
            if prop.get("key") == "MatchStatus":
                try:
                    n = int(prop.get("value", "1"))
                except ValueError:
                    n = 1
        counts[acc] = counts.get(acc, 0) + n
    return counts


def entry_gene(entry: dict) -> str:
    genes = entry.get("genes") or []
    if genes and genes[0].get("geneName"):
        return genes[0]["geneName"].get("value", "")
    return ""


# ----------------------------------------------------------------- collectors

def collect_signature_counts(out: Path) -> list[list]:
    print("[1/5] InterPro protein counts per Pfam signature")
    rows = []
    for pfam, name, note in SIGNATURES:
        n = interpro_count(pfam)
        print(f"       {pfam} {name:<18} {n:>8,}")
        rows.append([pfam, name, n, note])
    write_tsv(out / "interpro_signature_counts.tsv",
              ["pfam", "name", "n_proteins", "note"], rows)
    return rows


def collect_taxonomy_counts(out: Path) -> list[list]:
    print("[2/5] InterPro taxonomic distribution")
    rows = []
    for label, taxid in TAXA:
        counts = {}
        for pfam in TAXONOMY_SIGNATURES:
            counts[pfam] = interpro_count(pfam, taxid)
        print(f"       {label:<42} " +
              "  ".join(f"{p}={counts[p]:>6,}" for p in TAXONOMY_SIGNATURES))
        rows.append([label, taxid] + [counts[p] for p in TAXONOMY_SIGNATURES])
    write_tsv(out / "interpro_taxonomy_counts.tsv",
              ["taxon", "taxon_id"] + [f"n_{p}" for p in TAXONOMY_SIGNATURES],
              rows)
    return rows


def collect_reference_proteins(out: Path) -> tuple[list[list], list[list]]:
    print("[3/5] UniProt reference panel + sister family")
    ref_rows, arch_rows = [], []
    for symbol, acc, family in REFERENCE_PROTEINS:
        entry = uniprot_entry(acc)
        length = entry.get("sequence", {}).get("length")
        gene = entry_gene(entry)
        name = (entry.get("proteinDescription", {})
                     .get("recommendedName", {})
                     .get("fullName", {}).get("value", ""))
        reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
        pfams = entry_pfams(entry)
        print(f"       {symbol:<6} {acc}  {length:>5} aa  gene={gene}  "
              f"pfam={len(pfams)}")
        ref_rows.append([symbol, acc, family, gene, length,
                         "reviewed" if reviewed else "unreviewed", name])
        for pfam, pfname in ARCHITECTURE_PFAMS:
            arch_rows.append([symbol, acc, family, pfam, pfname,
                              int(pfam in pfams), pfams.get(pfam, 0)])
    write_tsv(out / "reference_proteins.tsv",
              ["symbol", "accession", "family", "uniprot_gene", "length_aa",
               "status", "protein_name"], ref_rows)
    write_tsv(out / "pfam_architecture.tsv",
              ["symbol", "accession", "family", "pfam", "pfam_name",
               "present", "n_matches"], arch_rows)
    return ref_rows, arch_rows


def collect_zebrafish(out: Path) -> list[list]:
    print("[4/5] Zebrafish PF08709 carriers")
    results = uniprot_search(
        ZEBRAFISH_QUERY,
        "accession,id,gene_names,protein_name,length,reviewed", size=200)
    rows = []
    for entry in results:
        rows.append([
            entry.get("primaryAccession", ""),
            entry_gene(entry),
            entry.get("sequence", {}).get("length"),
            "reviewed" if entry.get("entryType", "").startswith(
                "UniProtKB reviewed") else "unreviewed",
            (entry.get("proteinDescription", {})
                  .get("recommendedName", {})
                  .get("fullName", {}).get("value", "")
             or entry.get("proteinDescription", {})
                  .get("submissionNames", [{}])[0]
                  .get("fullName", {}).get("value", "")),
        ])
    rows.sort(key=lambda r: (r[1] or "zzz", -(r[2] or 0)))
    for r in rows[:25]:
        print(f"       {r[0]:<12} {str(r[1]):<12} {str(r[2]):>6} aa  {r[3]}")
    write_tsv(out / "zebrafish_itpr.tsv",
              ["accession", "gene", "length_aa", "status", "protein_name"],
              rows)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="results/s0_baseline")
    args = ap.parse_args()
    out = (ROOT / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"S0 database snapshot -> {out.relative_to(ROOT)}")
    sig = collect_signature_counts(out)
    tax = collect_taxonomy_counts(out)
    ref, arch = collect_reference_proteins(out)
    zeb = collect_zebrafish(out)

    print("[5/5] metadata")
    meta = {
        "task": "S0",
        "derived_on": date.today().isoformat(),
        "endpoints": {
            "interpro_signature_count":
                f"{INTERPRO}/protein/uniprot/entry/pfam/{{PF}}/?page_size=1",
            "interpro_taxonomy_count":
                f"{INTERPRO}/protein/uniprot/entry/pfam/{{PF}}"
                f"/taxonomy/uniprot/{{taxid}}/?page_size=1",
            "uniprot_entry": f"{UNIPROT}/{{acc}}.json",
            "uniprot_search": f"{UNIPROT}/search?query={{q}}",
        },
        "zebrafish_query": ZEBRAFISH_QUERY,
        "row_counts": {
            "interpro_signature_counts.tsv": len(sig),
            "interpro_taxonomy_counts.tsv": len(tax),
            "reference_proteins.tsv": len(ref),
            "pfam_architecture.tsv": len(arch),
            "zebrafish_itpr.tsv": len(zeb),
        },
    }
    (out / "db_snapshot_meta.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"  wrote {(out / 'db_snapshot_meta.json').relative_to(ROOT)}")
    print("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
