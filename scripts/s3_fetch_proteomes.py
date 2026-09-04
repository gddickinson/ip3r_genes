"""S3 step 2 — download UniProt reference proteomes for one taxonomic group.

The sweep database (the *denominator* of the HMM completeness claim) is
every UniProt **reference** proteome in the group: one canonical (per-gene)
FASTA per proteome, from the per-proteome files of the current UniProt
release. Isoform / "additional" sets are deliberately excluded — they add
isoforms of genes already present, not new genes, and the census counts
genes.

S3 uses `vertebrata`; S20 reuses the same code path and manifest format for
the non-vertebrate eukaryotes, so the two sweeps stay comparable.

Outputs (per group; `vertebrata` keeps the unsuffixed names):
  results/hmm_sweep/proteome_manifest[_<group>].tsv    upid, organism, taxid
  <data_root>/proteomes/<group>/UPxxx_taxid.fasta.gz   one file per proteome
  <data_root>/proteomes/<group>_refprot.fasta          concatenated sweep DB
  results/hmm_sweep/proteome_db_stats[_<group>].json   totals + failures

Resumable: gzip-valid files are skipped, so re-run after a failure.

Run:  python3 scripts/s3_fetch_proteomes.py --group vertebrata
      python3 scripts/s3_fetch_proteomes.py --list-groups
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from scripts.s3_hmm_lib import HMM_SWEEP_DIR, parse_uniprot_header  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

STREAM = ("https://rest.uniprot.org/proteomes/stream?query={query}"
          "&format=tsv&fields=upid,organism,organism_id,protein_count")
FTP_ROOT = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/reference_proteomes")
N_WORKERS = 8
RETRIES = 3

# Taxon groups. `vertebrata` plus the four non-vertebrate eukaryote groups
# partition Eukaryota, so no proteome is fetched or swept twice.
GROUPS: dict[str, dict] = {
    "vertebrata": {
        "query": "taxonomy_id:7742+AND+reference:true",
        "domain": "Eukaryota",
        "note": "S3's sweep set — where all three ITPR paralogs live",
    },
    "metazoa_nonvert": {
        "query": "taxonomy_id:33208+NOT+taxonomy_id:7742+AND+reference:true",
        "domain": "Eukaryota",
        "note": "invertebrates — the single-Itpr grade (S20)",
    },
    "fungi": {
        "query": "taxonomy_id:4751+AND+reference:true",
        "domain": "Eukaryota",
        "note": "S2 found calls only in early-diverging phyla, none in Dikarya",
    },
    "viridiplantae": {
        "query": "taxonomy_id:33090+AND+reference:true",
        "domain": "Eukaryota",
        "note": "S2: every call is Chlorophyta, Streptophyta 0/15 (S20)",
    },
    "protista_other": {
        "query": ("taxonomy_id:2759+NOT+taxonomy_id:33208"
                  "+NOT+taxonomy_id:33090+NOT+taxonomy_id:4751"
                  "+AND+reference:true"),
        "domain": "Eukaryota",
        "note": "amoebozoa, SAR, discoba — where iplA and the deep grade sit",
    },
}


def group_paths(group: str) -> tuple[Path, Path, Path, Path]:
    """(manifest, proteome dir, concatenated DB, stats) for a group."""
    suffix = "" if group == "vertebrata" else f"_{group}"
    root = require_data_root() / "proteomes"
    return (HMM_SWEEP_DIR / f"proteome_manifest{suffix}.tsv",
            root / group,
            root / f"{group}_refprot.fasta",
            HMM_SWEEP_DIR / f"proteome_db_stats{suffix}.json")


def log(msg: str) -> None:
    print(f"[s3_proteomes] {msg}", flush=True)


def fetch_manifest(query: str, manifest: Path) -> list[dict]:
    """Download (or reuse) the reference-proteome list for a query."""
    HMM_SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    if manifest.exists():
        text = manifest.read_text()
        log(f"reusing committed manifest {manifest}")
    else:
        r = requests.get(STREAM.format(query=query), timeout=180)
        r.raise_for_status()
        text = r.text
        manifest.write_text(text)
        log(f"manifest fetched → {manifest}")
    rows = []
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        upid, organism, taxid, pcount = line.split("\t")
        rows.append({"upid": upid, "organism": organism,
                     "taxid": int(taxid), "protein_count": int(pcount)})
    log(f"{len(rows)} reference proteomes, "
        f"{sum(r['protein_count'] for r in rows):,} proteins (UniProtKB count)")
    return rows


def gzip_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except OSError:
        return False


def download_one(row: dict, dest_dir: Path, ftp_base: str) -> tuple[str, str]:
    """Returns (upid, status) where status ∈ ok / cached / failed:<why>."""
    upid, taxid = row["upid"], row["taxid"]
    dest = dest_dir / f"{upid}_{taxid}.fasta.gz"
    if gzip_ok(dest):
        return upid, "cached"
    url = f"{ftp_base}/{upid}/{upid}_{taxid}.fasta.gz"
    last = ""
    for attempt in range(1, RETRIES + 1):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                if r.status_code == 404:
                    return upid, "failed:404"
                r.raise_for_status()
                tmp = dest.with_suffix(".part")
                with tmp.open("wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                tmp.rename(dest)
            if gzip_ok(dest):
                return upid, "ok"
            last = "gzip-invalid"
            dest.unlink(missing_ok=True)
        except requests.RequestException as exc:
            last = type(exc).__name__
        time.sleep(2 * attempt)
    return upid, f"failed:{last}"


def download_all(rows: list[dict], dest_dir: Path, ftp_base: str) -> dict:
    dest_dir.mkdir(parents=True, exist_ok=True)
    counts = {"ok": 0, "cached": 0, "failed": 0}
    failures: list[str] = []
    t0 = time.time()
    with ThreadPoolExecutor(N_WORKERS) as pool:
        futs = {pool.submit(download_one, row, dest_dir, ftp_base): row
                for row in rows}
        for i, fut in enumerate(as_completed(futs), 1):
            upid, status = fut.result()
            if status.startswith("failed"):
                counts["failed"] += 1
                failures.append(f"{upid} {status}")
                log(f"  {upid}: {status}")
            else:
                counts[status] += 1
            if i % 50 == 0 or i == len(rows):
                log(f"  {i}/{len(rows)} ({counts['ok']} new, "
                    f"{counts['cached']} cached, {counts['failed']} failed, "
                    f"{time.time() - t0:.0f}s)")
    return {"counts": counts, "failures": failures,
            "elapsed_s": round(time.time() - t0, 1)}


def concatenate(rows: list[dict], dest_dir: Path, db_path: Path) -> dict:
    """Stream-decompress every proteome into one plain-FASTA sweep DB."""
    n_seqs = n_res = 0
    t0 = time.time()
    with db_path.open("w") as out:
        for i, row in enumerate(rows, 1):
            src = dest_dir / f"{row['upid']}_{row['taxid']}.fasta.gz"
            with gzip.open(src, "rt") as f:
                for line in f:
                    if line.startswith(">"):
                        n_seqs += 1
                    else:
                        n_res += len(line.strip())
                    out.write(line)
            if i % 100 == 0 or i == len(rows):
                log(f"  concat {i}/{len(rows)} ({n_seqs:,} seqs)")
    return {"db_path": str(db_path), "n_proteomes": len(rows),
            "n_seqs": n_seqs, "n_residues": n_res,
            "db_bytes": db_path.stat().st_size,
            "concat_elapsed_s": round(time.time() - t0, 1)}


def spot_check(db_path: Path) -> dict:
    """Parse the first header as a format sanity check."""
    with db_path.open() as f:
        first = f.readline().strip()
    return parse_uniprot_header(first[1:]) if first.startswith(">") else {}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", default="vertebrata")
    ap.add_argument("--query", help="override the group's UniProt query")
    ap.add_argument("--domain", help="FTP domain dir (Eukaryota/Archaea/Bacteria)")
    ap.add_argument("--list-groups", action="store_true")
    ap.add_argument("--manifest-only", action="store_true",
                    help="fetch the proteome list and stop (counts + size check)")
    ap.add_argument("--no-concat", action="store_true")
    ap.add_argument("--ftp-root", default=FTP_ROOT,
                    help="reference-proteomes root; the EBI mirror "
                         "(https://ftp.ebi.ac.uk/pub/databases/uniprot/"
                         "current_release/knowledgebase/reference_proteomes) "
                         "works when ftp.uniprot.org is unreachable")
    args = ap.parse_args()

    if args.list_groups:
        for name, g in GROUPS.items():
            print(f"  {name:18s} {g['domain']:10s} {g['note']}")
            print(f"  {'':18s} {g['query']}")
        return 0

    spec = GROUPS.get(args.group, {})
    query = args.query or spec.get("query")
    domain = args.domain or spec.get("domain", "Eukaryota")
    if not query:
        ap.error(f"unknown group {args.group!r} and no --query given "
                 f"(known: {', '.join(GROUPS)})")

    manifest, dest_dir, db_path, stats_path = group_paths(args.group)
    log(f"group={args.group} domain={domain}")
    rows = fetch_manifest(query, manifest)
    if args.manifest_only:
        return 0

    dl = download_all(rows, dest_dir, f"{args.ftp_root}/{domain}")
    if dl["counts"]["failed"]:
        suffix = "" if args.group == "vertebrata" else f"_{args.group}"
        (HMM_SWEEP_DIR / f"proteome_fetch_failures{suffix}.txt").write_text(
            "\n".join(dl["failures"]) + "\n")
        log(f"ABORT before concat: {dl['counts']['failed']} downloads failed "
            "— re-run to resume")
        return 1
    if args.no_concat:
        log(f"downloads complete ({dl['counts']}); concat skipped")
        return 0

    stats = concatenate(rows, dest_dir, db_path)
    stats.update({"group": args.group, "query": query, "download": dl,
                  "first_header_parsed": spot_check(db_path)})
    stats_path.write_text(json.dumps(stats, indent=1))
    log(f"DB ready: {stats['n_seqs']:,} seqs, {stats['n_residues']:,} "
        f"residues, {stats['db_bytes'] / 1e9:.1f} GB → {db_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
