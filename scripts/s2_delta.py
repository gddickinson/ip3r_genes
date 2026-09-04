"""S2 step 6 — what the enumeration adds over the app's own searches.

Census v1 is what the application returns when it is pointed at the family
by name: the committed search bundles under `results/2026-*/` plus S1's
control panels. Those are *keyword and accession* searches — they find what
is already labelled ITPR. The enumeration is a *domain* search: it finds
what carries the signature, whatever it is called.

The delta is therefore the answer to "how much of this family is invisible
to a search that knows its name?", and it is reported by taxonomic group,
because the answer is not uniform — a well-annotated vertebrate is found by
either route, and a protist is not.

Every v1 accession that census v2 does **not** contain is listed
individually in `delta_missing_from_v2.tsv`, each with a verdict, because
that direction is the only place a systematic hole in the enumeration
would show up. Three of the four verdicts are structural rather than
faults: an accession outside UniProt's namespace was never in the searched
space, a `-2` isoform is represented by its canonical parent, and S1's
decoy panel is *supposed* to be absent. The fourth — a UniProt protein
with no seed signature annotated — is the real limit of a domain census,
and those are resolved individually against UniProt's own Pfam list rather
than assumed.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s2_lib import (                                    # noqa: E402
    OUT_DIR, ROOT, SEED_PFAMS, fetch, read_tsv, write_tsv,
)

#: UniProt accession syntax (UniProtKB accession format, both patterns).
UNIPROT_ACC = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")

BUNDLES = sorted((ROOT / "results").glob("2026-*/results.csv"))
PANELS = [ROOT / "results" / "benchmark_controls" / "panel_positives.json",
          ROOT / "results" / "benchmark_controls" / "panel_decoys.json"]


def load_v1() -> dict[str, dict]:
    """Every accession the app-level searches and the S1 panels produced."""
    out: dict[str, dict] = {}
    for path in BUNDLES:
        with path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                acc = (row.get("accession") or "").strip()
                if not acc:
                    continue
                out.setdefault(acc, {
                    "accession": acc, "gene": row.get("gene_symbol", ""),
                    "species": row.get("species", ""),
                    "length": row.get("length_aa", ""),
                    "origin": f"bundle:{path.parent.name}",
                    "source": row.get("source", "")})
    for path in PANELS:
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        records = data if isinstance(data, list) else data.get("records", [])
        for r in records:
            acc = (r.get("accession") or "").strip()
            if acc:
                out.setdefault(acc, {
                    "accession": acc, "gene": r.get("gene", ""),
                    "species": r.get("species", ""),
                    "length": r.get("length", ""),
                    "origin": f"panel:{path.stem}", "source": "S1"})
    return out


def resolve(accessions: list[str]) -> dict[str, dict]:
    """Ask UniProt what Pfam signatures these proteins actually carry.

    The point is to distinguish "the enumeration missed it" from "UniProt
    annotates no seed signature on it", which are different failures with
    different owners.
    """
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = ("https://rest.uniprot.org/uniprotkb/stream?query="
           + urllib.parse.quote(query)
           + "&format=tsv&fields=accession,length,xref_pfam,protein_name")
    status, body, _ = fetch(url, timeout_s=120, tries=5)
    if status != 200:
        return {}
    out = {}
    lines = body.decode("utf-8").splitlines()
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < 3:
            continue
        out[f[0]] = {"length": f[1],
                     "pfams": ";".join(sorted(p for p in f[2].split(";") if p))}
    return out


def classify(v1: dict[str, dict], only_v1: list[str],
             v2: dict[str, dict], offline: bool) -> list[dict]:
    rows, unresolved = [], []
    for acc in only_v1:
        r = dict(v1[acc])
        parent = acc.split("-")[0]
        if "panel_decoys" in r["origin"]:
            r["verdict"] = "s1_decoy (correctly absent)"
        elif not UNIPROT_ACC.match(parent):
            r["verdict"] = "outside UniProt namespace"
        elif "-" in acc and parent in v2:
            r["verdict"] = f"isoform of censused {parent}"
        else:
            r["verdict"] = "unresolved"
            unresolved.append(acc)
        r["pfams"] = ""
        rows.append(r)
    if unresolved and not offline:
        found = resolve(unresolved)
        for r in rows:
            if r["verdict"] != "unresolved":
                continue
            info = found.get(r["accession"])
            if info is None:
                r["verdict"] = "not in UniProt"
                continue
            r["pfams"] = info["pfams"]
            r["length"] = r["length"] or info["length"]
            carries = [p for p in SEED_PFAMS if p in info["pfams"]]
            r["verdict"] = ("ENUMERATION HOLE — carries "
                            + "+".join(carries)) if carries else (
                "no seed signature annotated")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--offline", action="store_true",
                    help="skip the UniProt resolution of unexplained rows")
    args = ap.parse_args()

    v2 = {r["accession"]: r for r in read_tsv(OUT_DIR / "census_v2.tsv")}
    v1 = load_v1()
    shared = set(v1) & set(v2)
    only_v1 = sorted(set(v1) - set(v2))
    only_v2 = sorted(set(v2) - set(v1))

    per_group: dict[str, dict] = {}
    for acc, r in v2.items():
        g = per_group.setdefault(r["group"] or "unclassified", {
            "group": r["group"] or "unclassified", "v2_records": 0,
            "v2_itpr": 0, "in_v1": 0, "new_in_v2": 0, "v2_species": set()})
        g["v2_records"] += 1
        g["v2_itpr"] += int(r["call"] == "ITPR")
        g["v2_species"].add(r["taxon_id"])
        if acc in v1:
            g["in_v1"] += 1
        else:
            g["new_in_v2"] += 1
    rows = []
    for g in per_group.values():
        g["v2_species"] = len(g["v2_species"])
        g["pct_new"] = (f"{100 * g['new_in_v2'] / g['v2_records']:.1f}"
                        if g["v2_records"] else "")
        rows.append(g)
    rows.sort(key=lambda r: -r["v2_records"])
    write_tsv(OUT_DIR / "delta_by_group.tsv", rows,
              ["group", "v2_records", "v2_itpr", "v2_species", "in_v1",
               "new_in_v2", "pct_new"])

    missing = classify(v1, only_v1, v2, args.offline)
    write_tsv(OUT_DIR / "delta_missing_from_v2.tsv", missing,
              ["accession", "verdict", "gene", "species", "length", "pfams",
               "origin", "source"])
    verdicts: dict[str, int] = {}
    for r in missing:
        key = r["verdict"].split(" — ")[0]
        if key.startswith("isoform of censused"):
            key = "isoform of a censused canonical"
        verdicts[key] = verdicts.get(key, 0) + 1
    write_tsv(OUT_DIR / "delta_missing_verdicts.tsv",
              [{"verdict": k, "accessions": n}
               for k, n in sorted(verdicts.items(), key=lambda kv: -kv[1])],
              ["verdict", "accessions"])

    summary = {
        "v1_accessions": len(v1), "v2_accessions": len(v2),
        "shared": len(shared), "only_v1": len(only_v1),
        "only_v2": len(only_v2),
        "bundles": [p.parent.name for p in BUNDLES],
        "growth_factor": round(len(v2) / len(v1), 1) if v1 else None,
    }
    (OUT_DIR / "delta_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"[s2] census v1 {len(v1)} → v2 {len(v2)} "
          f"(×{summary['growth_factor']}); {len(only_v1)} v1 accessions "
          f"absent from v2:")
    for k, n in sorted(verdicts.items(), key=lambda kv: -kv[1]):
        print(f"     {n:5d}  {k}")
    for r in rows:
        print(f"     {r['group']:<26s} {r['v2_records']:6d} records, "
              f"{r['pct_new']:>5s} % new")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
