#!/usr/bin/env python3
"""Pfam domain **coordinates** for the reference and sister panels (review figures).

`results/s0_baseline/pfam_architecture.tsv` records which family signatures a
protein carries and how many copies. The figures need something it does not
hold: *where* each copy sits along the chain, so the subunit can be drawn to
scale and the ITPR / RyR comparison of §7.1 can be shown rather than asserted.

Source is the InterPro API, one query per accession:

    entry/pfam/protein/uniprot/{acc}  ->  entry_protein_locations[].fragments[]

Output (committed, small):

    results/s0_baseline/review_figures/domain_coords.tsv

Re-run whenever the reference panel in `src/utils/family.py` changes.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "results" / "s0_baseline" / "review_figures"
OUT_TSV = OUT_DIR / "domain_coords.tsv"
REF_TSV = ROOT / "results" / "s0_baseline" / "reference_proteins.tsv"
PANEL_FASTA = ROOT / "results" / "benchmark_controls" / "panel_positives.fasta"

#: The non-vertebrate single-*itpr* grade. Included because the signature that
#: defines this family does not find all of it: *Dictyostelium* iplA is a
#: characterised IP3 receptor that carries neither PF08709 nor the pore domain,
#: which is a fact about the records the census is built from (§7.5, Q1/Q4).
GRADE_PANEL = ["P29993", "Q9Y0A1", "Q9NA13"]

API = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/?page_size=100"

#: Pfam short names are terse; these are the review's own words for them so a
#: figure legend reads like §2.1 rather than like a database dump.
PFAM_LABEL = {
    "PF08709": "IP$_3$-binding core (β-trefoil)",
    "PF02815": "MIR",
    "PF01365": "RIH",
    "PF08454": "RIH-associated",
    "PF00520": "Ion_trans pore",
    "PF02026": "RyR repeat",
    "PF06459": "RyR TM4-6",
    "PF21119": "RyR jsol",
    "PF00622": "SPRY",
    "PF13499": "EF-hand pair",
}
#: Which family each signature is diagnostic for, as §2.1/§7.1 use the terms.
PFAM_CLASS = {
    "PF08709": "shared", "PF02815": "shared", "PF01365": "shared",
    "PF08454": "shared", "PF00520": "generic",
    "PF02026": "ryr_only", "PF06459": "ryr_only", "PF21119": "ryr_only",
    "PF00622": "ryr_only", "PF13499": "ryr_only",
}


def fetch(acc: str, retries: int = 4) -> list[dict]:
    """Every Pfam entry hit on one accession, with fragment coordinates."""
    url = API.format(acc=acc)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "ip3r-review-figures"})
            with urllib.request.urlopen(req, timeout=90) as fh:
                return json.load(fh).get("results", [])
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == retries - 1:
                raise SystemExit(f"InterPro failed for {acc}: {exc}")
            time.sleep(2 * (attempt + 1))
    return []


def read_panel() -> list[dict]:
    lines = REF_TSV.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t")))
            for ln in lines[1:] if ln.strip()]
    return rows + read_grade()


def read_grade() -> list[dict]:
    """The invertebrate / non-metazoan grade, from the committed FASTA."""
    seqs, key, head, buf = {}, None, "", []
    for ln in PANEL_FASTA.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            if key:
                seqs[key] = (head, "".join(buf))
            head = ln[1:]
            key = head.split("|")[0]
            buf = []
        else:
            buf.append(ln.strip())
    if key:
        seqs[key] = (head, "".join(buf))
    out = []
    for acc in GRADE_PANEL:
        if acc not in seqs:
            raise SystemExit(f"{acc} not in {PANEL_FASTA.name}")
        head, seq = seqs[acc]
        out.append({"symbol": head.split("|")[1], "accession": acc,
                    "family": "ITPR_grade", "length_aa": str(len(seq))})
    return out


def main() -> int:
    panel = read_panel()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for prot in panel:
        acc, sym = prot["accession"], prot["symbol"]
        entries = fetch(acc)
        n = 0
        for res in entries:
            pfam = res["metadata"]["accession"]
            for prot_hit in res["proteins"]:
                for loc in prot_hit["entry_protein_locations"]:
                    for frag in loc["fragments"]:
                        rows.append({
                            "symbol": sym,
                            "accession": acc,
                            "family": prot["family"],
                            "length_aa": prot["length_aa"],
                            "pfam": pfam,
                            "pfam_name": res["metadata"]["name"],
                            "label": PFAM_LABEL.get(pfam, pfam),
                            "shared_class": PFAM_CLASS.get(pfam, "other"),
                            "start": frag["start"],
                            "end": frag["end"],
                        })
                        n += 1
        print(f"  {sym:6s} {acc:8s} {n:2d} domain copies")
        time.sleep(0.34)                    # InterPro politeness

    rows.sort(key=lambda r: (r["family"], r["symbol"], r["start"]))
    cols = ["symbol", "accession", "family", "length_aa", "pfam", "pfam_name",
            "label", "shared_class", "start", "end"]
    with OUT_TSV.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {OUT_TSV.relative_to(ROOT)}  ({len(rows)} domain copies)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
