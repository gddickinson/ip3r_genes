"""S20 step 1 (data) — fetch the six groups' reference proteomes.

Reuses S3's download and concatenation code unchanged, so the non-vertebrate
sweep DB is built exactly as the vertebrate one was. Two things S20 adds:

**Sampling.** `bacteria_genus` is 17,981 proteomes reduced to one per genus
by the rule in `s20_groups.py`. The manifest records every available
proteome's genus and how many the pick stood for, so the sample is auditable
from the committed table.

**A 404 is an exclusion, not an abort.** S3 aborts on any failed download,
which is right for a claim whose denominator must be complete. Here the
manifest is UniProt's current *proteome list* and the files come from the
current *release FTP tree*, and the two are not always in step: a proteome
listed but not published 404s permanently, and retrying cannot fix it. So a
404 is recorded in `proteome_unavailable_<group>.tsv` with the organism and
protein count it took out of the denominator, and the manifest marks that
row `unavailable` instead of `swept`. Every other failure still aborts,
because those are network problems a re-run does fix.

Run:  python3 scripts/s20_fetch.py --group fungi
      python3 scripts/s20_fetch.py --all
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from scripts.s20_groups import ALL_GROUPS, GROUPS, apply_sample  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    MANIFEST_FIELDS, group_paths, live, log, write_json, write_tsv,
)
from scripts.s3_fetch_proteomes import (  # noqa: E402
    STREAM, concatenate, download_all, spot_check,
)

FTP_ROOT = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/reference_proteomes")


def fetch_proteome_list(group: str, cache: Path) -> list[dict]:
    """UniProt's reference-proteome list for a group, cached verbatim."""
    if cache.exists():
        text = cache.read_text()
        log("s20_fetch", f"reusing cached proteome list {cache.name}")
    else:
        r = requests.get(STREAM.format(query=GROUPS[group]["query"]),
                         timeout=600)
        r.raise_for_status()
        text = r.text
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(text)
    rows = []
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        upid, organism, taxid, pcount = line.split("\t")
        rows.append({"upid": upid, "organism": organism,
                     "taxid": int(taxid), "protein_count": int(pcount)})
    return rows


def fetch_group(group: str, ftp_root: str, manifest_only: bool) -> int:
    spec = GROUPS[group]
    paths = group_paths(group)
    raw_list = paths["manifest"].with_name(f"_raw_proteome_list_{group}.tsv")

    available = fetch_proteome_list(group, raw_list)
    selected, sample_stats = apply_sample(group, available)
    log("s20_fetch", f"{group}: {len(available)} reference proteomes → "
                     f"{len(selected)} selected "
                     f"({sample_stats['proteins_selected']:,} proteins) "
                     f"[{sample_stats['rule']}]")
    if manifest_only:
        write_json(paths["stats"].with_name(f"sample_{group}.json"),
                   sample_stats)
        return 0

    paths["dir"].mkdir(parents=True, exist_ok=True)
    dl = download_all(selected, paths["dir"], f"{ftp_root}/{spec['domain']}")

    hard = [f for f in dl["failures"] if not f.endswith("failed:404")]
    gone = {f.split()[0] for f in dl["failures"] if f.endswith("failed:404")}
    if hard:
        log("s20_fetch", f"ABORT: {len(hard)} non-404 failures — re-run to "
                         f"resume ({hard[:5]})")
        return 1

    by_upid = {r["upid"]: r for r in selected}
    if gone:
        write_tsv(paths["unavailable"],
                  ["upid", "organism", "taxid", "protein_count", "reason"],
                  [dict(by_upid[u], reason="404 in current release FTP tree")
                   for u in sorted(gone)])
        log("s20_fetch", f"{len(gone)} listed proteomes are not published in "
                         f"the current release → {paths['unavailable'].name}")

    swept = [r for r in selected if r["upid"] not in gone]
    write_tsv(paths["manifest"], MANIFEST_FIELDS,
              [{**r, "group": group,
                "genus": r.get("genus", ""),
                "genus_n_available": r.get("genus_n_available", ""),
                "status": "unavailable" if r["upid"] in gone else "swept"}
               for r in selected])

    stats = concatenate(swept, paths["dir"], paths["db"])
    stats.update({"group": group, "query": spec["query"],
                  "sample": sample_stats, "download": dl,
                  "n_unavailable": len(gone),
                  "n_swept_proteomes": len(swept),
                  "first_header_parsed": spot_check(paths["db"])})
    write_json(paths["stats"], stats)
    log("s20_fetch", f"{group} DB ready: {stats['n_seqs']:,} seqs, "
                     f"{stats['n_residues']:,} residues, "
                     f"{stats['db_bytes'] / 1e9:.1f} GB")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", action="append", choices=list(GROUPS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--ftp-root", default=FTP_ROOT)
    args = ap.parse_args()

    groups = list(ALL_GROUPS) if args.all else (args.group or [])
    if not groups:
        ap.error("give --group (repeatable) or --all")

    for i, g in enumerate(groups):
        live([(f"fetch {name}", name in groups[:i]) for name in groups]
             + [("hmmsearch", False), ("verdicts", False),
                ("jackhmmer", False), ("census v5", False)])
        rc = fetch_group(g, args.ftp_root, args.manifest_only)
        if rc:
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
