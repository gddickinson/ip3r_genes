"""S2 step 2 — the architecture sweep that makes the ITPR/RYR call cheap.

The call needs each record's **full Pfam architecture**, not just the seed
signature it was found by. Asking InterPro protein-by-protein would be
~15,000 requests against an API that was answering one request in twelve
when this was written. UniProt answers the same question in one streamed
request, because `xref_pfam` is a returnable field: the query

    (xref:pfam-PF08709) OR (xref:pfam-PF01365) OR (xref:pfam-PF08454)

returns every protein in the seeded search space with its complete Pfam
list, its length, its gene names, its organism and its full lineage. The
architecture column is UniProt's copy of the same Pfam matches InterPro
serves, so this is a second view of one underlying dataset — which is
exactly why `s2_call.py` reconciles the two lists instead of trusting
either alone.

The raw TSV (~14 MB with full lineages) is archived under the data root,
and so is the parsed intermediate — every column of it is carried into
`census_v2.tsv`, so committing both would put the same 4 MB in the repo
twice. What the repo keeps from this step is the group tally.
"""

from __future__ import annotations

import argparse
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.data_root import require_data_root          # noqa: E402
from scripts.s2_lib import (                               # noqa: E402
    OUT_DIR, ROOT, SEED_PFAMS, fetch, live, write_tsv,
)

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = ("accession,reviewed,protein_name,gene_names,organism_name,"
          "organism_id,length,xref_pfam,fragment,protein_existence,lineage")

#: Coarse taxonomic buckets, tested in order — the first clade in a record's
#: lineage that matches wins. Vertebrata is split out of Metazoa because the
#: three-paralog cell model is a vertebrate statement; everything below it
#: is where the census's open questions live.
GROUPS = [
    ("Vertebrata", "Vertebrata"),
    ("Metazoa", "Metazoa (non-vertebrate)"),
    ("Fungi", "Fungi"),
    ("Viridiplantae", "Viridiplantae"),
    ("Sar", "SAR"),
    ("Discoba", "Discoba"),
    ("Amoebozoa", "Amoebozoa"),
    ("Eukaryota", "Eukaryota (other)"),
    ("Bacteria", "Bacteria"),
    ("Archaea", "Archaea"),
    ("Viruses", "Viruses"),
]

RANKS = ["kingdom", "phylum", "class", "order", "family"]

COLS = ["accession", "reviewed", "gene", "protein_name", "species",
        "taxon_id", "length", "fragment", "protein_existence", "pfams",
        "n_pfam", "group"] + RANKS


def raw_path() -> Path:
    d = require_data_root(ROOT) / "raw_api" / "uniprot"
    d.mkdir(parents=True, exist_ok=True)
    return d / "s2_seed_sweep.tsv"


def download(dest: Path) -> int:
    query = " OR ".join(f"(xref:pfam-{p})" for p in SEED_PFAMS)
    url = (f"{STREAM}?query={urllib.parse.quote(query)}"
           f"&format=tsv&fields={FIELDS}")
    status, body, _ = fetch(url, timeout_s=300, tries=6)
    if status != 200:
        raise RuntimeError(f"UniProt stream returned HTTP {status}")
    dest.write_bytes(body)
    return len(body.decode("utf-8").splitlines()) - 1


def parse_lineage(lineage: str) -> tuple[str, dict[str, str]]:
    """UniProt's ranked lineage string → (coarse group, rank columns).

    Entries look like `Chordata (phylum), Mammalia (class), …`. Ranks are
    read where UniProt states one; the coarse group is decided by clade
    membership, which is stated for every record whether or not the ranks
    are.
    """
    names, ranked = [], {}
    for part in lineage.split(", "):
        part = part.strip()
        if not part:
            continue
        if part.endswith(")") and " (" in part:
            name, rank = part.rsplit(" (", 1)
            rank = rank[:-1]
        else:
            name, rank = part, ""
        names.append(name)
        if rank in RANKS and rank not in ranked:
            ranked[rank] = name
    group = "unclassified"
    nameset = set(names)
    for clade, label in GROUPS:
        if clade in nameset:
            group = label
            break
    return group, ranked


def parse(src: Path) -> list[dict]:
    rows: list[dict] = []
    with src.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue

            def col(name: str) -> str:
                return f[idx[name]].strip()

            pfams = [p for p in col("Pfam").split(";") if p]
            group, ranked = parse_lineage(col("Taxonomic lineage"))
            genes = col("Gene Names").split()
            rows.append({
                "accession": col("Entry"),
                "reviewed": int(col("Reviewed") == "reviewed"),
                "gene": genes[0] if genes else "",
                "protein_name": col("Protein names")[:120],
                "species": col("Organism"),
                "taxon_id": col("Organism (ID)"),
                "length": col("Length"),
                "fragment": col("Fragment"),
                "protein_existence": col("Protein existence"),
                "pfams": ";".join(sorted(pfams)),
                "n_pfam": len(pfams),
                "group": group,
                **{r: ranked.get(r, "") for r in RANKS},
            })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parse-only", action="store_true",
                    help="reparse the archived TSV, no network")
    args = ap.parse_args()

    dest = raw_path()
    live("S2", [("uniprot sweep", False), ("parse", False)])
    if not args.parse_only:
        n = download(dest)
        print(f"[s2] UniProt stream: {n} records → {dest}", flush=True)
    if not dest.exists():
        raise SystemExit(f"no archived sweep at {dest}; run without "
                         f"--parse-only")
    live("S2", [("uniprot sweep", True), ("parse", False)])

    rows = parse(dest)
    write_tsv(dest.parent / "uniprot_records.tsv", rows, COLS)
    groups: dict[str, int] = {}
    for r in rows:
        groups[r["group"]] = groups.get(r["group"], 0) + 1
    write_tsv(OUT_DIR / "uniprot_groups.tsv",
              [{"group": g, "records": n} for g, n in
               sorted(groups.items(), key=lambda kv: -kv[1])],
              ["group", "records"])
    live("S2", [("uniprot sweep", True), ("parse", True)])
    print(f"[s2] intermediate → {dest.parent / 'uniprot_records.tsv'}")
    print(f"[s2] parsed {len(rows)} records into "
          f"{len(groups)} taxonomic groups")
    for g, n in sorted(groups.items(), key=lambda kv: -kv[1]):
        print(f"     {n:6d}  {g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
