"""S20 step 2 — taxid → ranked lineage for every proteome and every hit.

The range result *is* a taxonomy statement, so every row in this task has to
carry where it sits, and the placement has to come from one source rather
than from the organism string. UniProt's taxonomy API gives the full ranked
lineage per taxid; this fetches it in batches, archives every raw response
under the data root, and re-parses the committed table out of the archive so
the lineage is reproducible offline.

The coarse `group` bucket is **imported from `s2_uniprot.GROUPS`, not
re-declared** — S2 bucketed the InterPro census with that table and census
v5 has to stack on the same one, or "Viridiplantae 0/15 in S2" and
"Viridiplantae N in S20" would be counts of different things.

Run:  python3 scripts/s20_taxa.py --groups           # every swept proteome
      python3 scripts/s20_taxa.py --taxids 3702,4932
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from scripts.s2_uniprot import GROUPS as COARSE_GROUPS  # noqa: E402
from scripts.s20_groups import ALL_GROUPS  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    S20_DIR, log, read_tsv, require_data_root, swept_manifest, write_tsv,
)

#: The **exact-id** endpoint, not `taxonomy/search`. UniProt's search field
#: `tax_id:` is hierarchical — `tax_id:3702` matches *Arabidopsis thaliana*
#: and every strain and subspecies under it — so a 100-term OR query can
#: match far more than 100 records, the page caps at `size`, and some of the
#: taxids actually asked for fall off the end and come back as blank rows.
#: That is exactly what happened on the first run of this module: a batch
#: returned 100 results, none of which were the *Arabidopsis* and
#: *Chlamydomonas* it had been asked for. `taxonIds/` returns one record per
#: id and nothing else, and the fetch verifies the set it got back.
TAXON_IDS = "https://rest.uniprot.org/taxonomy/taxonIds"

#: **Below** the endpoint's 25-record page cap, measured rather than assumed:
#: asking for 50 or 100 ids returns 200 OK with 25 records and no indication
#: that the rest were dropped. With the batch under the cap, a short page can
#: only mean a taxid UniProt genuinely has no record for, so "missing" and
#: "truncated" cannot be confused — which is the distinction the first
#: version of this module got wrong in the other direction.
BATCH = 20
PAGE_CAP = 25      # measured: 50 or 100 ids come back as 200 OK with 25 rows
RETRIES = 4

RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]
TAXA_FIELDS = (["taxid", "scientific_name", "rank", "group"] + RANKS
               + ["clade", "clade_from", "lineage"])

#: Lineage names that are containers rather than groups anyone reports a
#: result by — skipped when falling back to an unranked clade.
_NOT_A_CLADE = {"cellular organisms", "Eukaryota", "Bacteria", "Archaea",
                "Opisthokonta", "Sar", "SAR", "Amorphea", "Diaphoretickes"}


def derive_clade(ranked: dict[str, str], names: list[str]) -> tuple[str, str]:
    """The deepest group this taxon can be *reported under*, and where it came
    from.

    UniProt gives no `phylum` rank to some of the lineages this task most
    needs to name — the choanoflagellates and filastereans that sit at the
    animal transition, the cryptophytes, the apusozoans. Reporting them as
    "unclassified" would put the family's most interesting neighbours in the
    same bucket as a parse failure. So the label falls back phylum → class →
    the first informative unranked lineage name, and records which, so the
    column can never be mistaken for a phylum assignment.
    """
    if ranked.get("phylum"):
        return ranked["phylum"], "phylum"
    if ranked.get("class"):
        return ranked["class"], "class"
    for name in names:
        if name and name not in _NOT_A_CLADE:
            return name, "lineage"
    return "unclassified", "none"

TAXA_TSV = S20_DIR / "taxonomy.tsv"


def archive_dir() -> Path:
    d = require_data_root() / "raw_api" / "uniprot_taxonomy"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _batch_path(taxids: list[int]) -> Path:
    return archive_dir() / f"taxonids_{taxids[0]}_{len(taxids)}.json"


def _results(payload) -> list[dict]:
    return payload["results"] if isinstance(payload, dict) else payload


def fetch_batch(taxids: list[int]) -> list[dict]:
    """One exact-id taxonomy batch, archived verbatim and reused if present.

    Verifies what came back against what was asked for: a batch that returns
    the right *number* of records but not the right *ones* is the failure
    this endpoint was chosen to avoid, and silence about it is how it went
    unnoticed the first time.
    """
    dest = _batch_path(taxids)
    if dest.exists():
        try:
            return _results(json.loads(dest.read_text()))
        except (json.JSONDecodeError, KeyError):
            dest.unlink()
    url = f"{TAXON_IDS}/{','.join(str(t) for t in taxids)}"
    params = {"format": "json",
              "fields": "id,scientific_name,rank,lineage"}
    last = ""
    for attempt in range(1, RETRIES + 1):
        try:
            r = requests.get(url, params=params, timeout=180)
            if r.status_code == 200:
                rows = _results(r.json())
                got = {int(x["taxonId"]) for x in rows if x.get("taxonId")}
                extra = got - set(taxids)
                if extra:
                    raise RuntimeError(
                        f"taxonomy batch returned {len(extra)} taxid(s) that "
                        f"were not requested: {sorted(extra)[:5]}")
                if len(rows) >= PAGE_CAP:
                    raise RuntimeError(
                        f"taxonomy batch returned {len(rows)} records for "
                        f"{len(taxids)} ids — at the endpoint's {PAGE_CAP}-row "
                        "page cap, so the page may be truncated; lower BATCH")
                dest.write_text(r.text)
                return rows
            last = f"HTTP {r.status_code}"
        except (requests.RequestException, json.JSONDecodeError) as exc:
            last = type(exc).__name__
        time.sleep(3 * attempt)
    raise RuntimeError(f"taxonomy batch failed ({last}): {taxids[:3]}…")


def parse_record(rec: dict) -> dict:
    """One taxonomy record → the committed row, ranks and coarse group."""
    lineage = rec.get("lineage", []) or []
    names = [e.get("scientificName", "") for e in lineage]
    ranked: dict[str, str] = {}
    for entry in lineage:
        rank = entry.get("rank", "")
        if rank in RANKS and rank not in ranked:
            ranked[rank] = entry.get("scientificName", "")
    # A record's own rank belongs in its own column too — a proteome whose
    # taxid *is* a genus has no genus in its ancestry.
    own_rank = rec.get("rank", "")
    if own_rank in RANKS and own_rank not in ranked:
        ranked[own_rank] = rec.get("scientificName", "")

    nameset = set(names) | {rec.get("scientificName", "")}
    group = "unclassified"
    for clade, label in COARSE_GROUPS:
        if clade in nameset:
            group = label
            break
    clade, clade_from = derive_clade(ranked, names)
    row = {"taxid": rec.get("taxonId"),
           "scientific_name": rec.get("scientificName", ""),
           "rank": own_rank, "group": group,
           "clade": clade, "clade_from": clade_from,
           "lineage": "; ".join(n for n in names if n)}
    row.update({r: ranked.get(r, "") for r in RANKS})
    return row


def resolve(taxids, existing: dict[int, dict] | None = None) -> dict[int, dict]:
    """taxid → lineage row for every requested taxid, batched and archived."""
    out = dict(existing or {})
    todo = sorted({int(t) for t in taxids} - set(out))
    if not todo:
        return out
    log("s20_taxa", f"resolving {len(todo)} taxids "
                    f"({len(out)} already cached)")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        for rec in fetch_batch(chunk):
            row = parse_record(rec)
            if row["taxid"] is not None:
                out[int(row["taxid"])] = row
        if (i // BATCH) % 10 == 0 or i + BATCH >= len(todo):
            log("s20_taxa", f"  {min(i + BATCH, len(todo))}/{len(todo)}")
    missing = [t for t in todo if t not in out]
    if missing:
        # A taxid UniProt's taxonomy search does not return is recorded as
        # such rather than dropped: an empty lineage row is visible in the
        # table, a missing row is not.
        log("s20_taxa", f"{len(missing)} taxids returned no taxonomy record")
        for t in missing:
            out[t] = {"taxid": t, "scientific_name": "", "rank": "",
                      "group": "unclassified", "clade": "unclassified",
                      "clade_from": "none", "lineage": "",
                      **{r: "" for r in RANKS}}
    return out


def load_table() -> dict[int, dict]:
    if not TAXA_TSV.exists():
        return {}
    return {int(r["taxid"]): r for r in read_tsv(TAXA_TSV)}


def save_table(table: dict[int, dict]) -> None:
    S20_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(TAXA_TSV, TAXA_FIELDS,
              [table[t] for t in sorted(table)])
    log("s20_taxa", f"{len(table)} taxa → {TAXA_TSV}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", action="store_true",
                    help="resolve every swept proteome's taxid")
    ap.add_argument("--census", action="store_true",
                    help="also resolve every taxid in census v4, so v5 can "
                         "carry lineage columns on every row (D8) rather than "
                         "only on the rows this task added")
    ap.add_argument("--taxids", help="comma-separated extra taxids")
    ap.add_argument("--reparse", action="store_true",
                    help="rebuild the table from the archived responses "
                         "(no network) — use after changing the parse")
    args = ap.parse_args()

    want: set[int] = set()
    if args.groups:
        want |= {r["taxid"] for r in swept_manifest(ALL_GROUPS)}
    if args.census:
        census = (PROJECT_ROOT / "results" / "census_v4" / "census_v4.tsv")
        want |= {int(r["taxon_id"]) for r in read_tsv(census)
                 if (r.get("taxon_id") or "").isdigit()}
        # S5b's genomic gene models carry a genome accession and no taxid,
        # so without the manifest's taxids 488 real ITPR records would sit in
        # census v5 with no lineage at all — and D8 says every row carries one.
        manifest = PROJECT_ROOT / "results" / "genome_manifest.tsv"
        if manifest.exists():
            want |= {int(r["taxid"]) for r in read_tsv(manifest)
                     if (r.get("taxid") or "").isdigit()}
    if args.taxids:
        want |= {int(t) for t in args.taxids.split(",") if t.strip()}
    if not want:
        ap.error("give --groups and/or --taxids")

    table = resolve(want, {} if args.reparse else load_table())
    save_table(table)
    by_group: dict[str, int] = {}
    for t in want:
        g = table[t]["group"]
        by_group[g] = by_group.get(g, 0) + 1
    for g, n in sorted(by_group.items(), key=lambda kv: -kv[1]):
        log("s20_taxa", f"  {g:26s} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
