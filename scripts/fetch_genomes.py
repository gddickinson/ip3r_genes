"""fetch_genomes.py — download manifest genomes into <data_root>/genomes/.

Resumable: a completed genome leaves a `.done` marker and is skipped on
rerun; an interrupted download is simply re-fetched, because nothing is
installed until it has verified. Every extracted file is checked against the
`md5sum.txt` NCBI packs inside the zip — a truncated download is a silent
corruption otherwise, and S5 would search a truncated genome and report an
absence.

For operation within a storage budget, `--search-cmd '… {fna} …'
--delete-after-search` runs a command per genome and purges the FASTA
afterwards, so the peak footprint is one genome rather than the manifest's
553 GB. The S5 driver can also import `fetch_one()` / `purge_one()`.

Usage:
  python3 scripts/fetch_genomes.py                              # all pending
  python3 scripts/fetch_genomes.py --only GCF_…,GCF_…
  python3 scripts/fetch_genomes.py --limit 2 --smallest-first   # smoke test
  python3 scripts/fetch_genomes.py --dry-run
  python3 scripts/fetch_genomes.py --verify                     # re-check .done
  python3 scripts/fetch_genomes.py --purge GCF_…
  python3 scripts/fetch_genomes.py \\
      --search-cmd 'tblastn -query bait.faa -subject {fna}' --delete-after-search
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.data_root import require_data_root  # noqa: E402
from scripts.s4_manifest_lib import datasets_bin  # noqa: E402

MANIFEST = PROJECT_ROOT / "results" / "genome_manifest.tsv"


def read_manifest(path: Path = MANIFEST) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genomes_root() -> Path:
    return require_data_root() / "genomes"


def genome_dir(accession: str) -> Path:
    return genomes_root() / accession


def is_done(accession: str) -> bool:
    return (genome_dir(accession) / ".done").exists()


def fna_path(accession: str) -> Path | None:
    hits = sorted(genome_dir(accession).glob("*.fna"))
    return hits[0] if hits else None


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_extract(zf: zipfile.ZipFile, dest: Path) -> None:
    for name in zf.namelist():
        if name.startswith(("/", "..")) or ".." in Path(name).parts:
            raise RuntimeError(f"unsafe path in zip: {name}")
    zf.extractall(dest)


def _verify_md5s(staging: Path) -> int:
    """Check every file listed in NCBI's md5sum.txt; return count verified."""
    md5file = staging / "md5sum.txt"
    if not md5file.exists():
        raise RuntimeError("zip is missing md5sum.txt")
    n = 0
    for line in md5file.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        want, rel = parts
        target = staging / rel
        if not target.exists():
            raise RuntimeError(f"md5sum.txt names a missing file: {rel}")
        got = _md5(target)
        if got != want:
            raise RuntimeError(f"md5 mismatch for {rel}: {got} != {want}")
        n += 1
    if n == 0:
        raise RuntimeError("md5sum.txt contained no checksums")
    return n


def fetch_one(accession: str, quiet: bool = False, attempts: int = 4) -> Path:
    """Download, verify and install one genome. Returns the genome dir.

    Download → extract → verify are retried *together* with backoff: NCBI's
    endpoint intermittently drops streams, and a truncated zip surfaces as an
    md5 failure rather than as a download error, so retrying only the
    transfer would install corrupt sequence.
    """
    dest = genome_dir(accession)
    if is_done(accession):
        return dest
    zips = genomes_root() / "_zips"
    staging = genomes_root() / "_staging" / accession
    zips.mkdir(parents=True, exist_ok=True)
    zip_path = zips / f"{accession}.zip"

    last_err = ""
    n_verified = 0
    for attempt in range(1, attempts + 1):
        if staging.exists():
            shutil.rmtree(staging)
        zip_path.unlink(missing_ok=True)
        proc = subprocess.run(
            [datasets_bin(), "download", "genome", "accession", accession,
             "--include", "genome,seq-report", "--filename", str(zip_path),
             "--no-progressbar"], capture_output=True, text=True)
        if proc.returncode != 0:
            last_err = f"datasets download failed: {proc.stderr.strip()[:200]}"
        else:
            try:
                with zipfile.ZipFile(zip_path) as zf:
                    _safe_extract(zf, staging)
                n_verified = _verify_md5s(staging)
                break
            except (zipfile.BadZipFile, RuntimeError) as exc:
                last_err = str(exc)[:200]
        if attempt == attempts:
            raise RuntimeError(
                f"{accession}: {last_err} (after {attempts} attempts)")
        if not quiet:
            print(f"  retry {attempt}/{attempts - 1} for {accession}: "
                  f"{last_err}")
        time.sleep(5 * attempt)

    payload = staging / "ncbi_dataset" / "data"
    src_dir = payload / accession
    if not src_dir.is_dir():
        raise RuntimeError(f"unexpected zip layout for {accession}")
    if dest.exists():
        shutil.rmtree(dest)
    shutil.move(str(src_dir), str(dest))
    for extra in payload.glob("*.jsonl"):   # assembly_data_report etc.
        shutil.move(str(extra), str(dest / extra.name))
    # Keep the checksums beside the installed copy so --verify can re-check
    # a genome months later without re-fetching. Rewritten to basenames and
    # to *only the files that were kept* (the zip also checksums
    # dataset_catalog.json, which the install drops): --verify can then treat
    # a name it cannot find as a failure instead of skipping it, which is the
    # only way a deleted or truncated .fna gets caught.
    kept = []
    for line in (staging / "md5sum.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and (dest / Path(parts[1]).name).exists():
            kept.append(f"{parts[0]}  {Path(parts[1]).name}")
    (dest / "md5sum.txt").write_text("\n".join(kept) + "\n")

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    (dest / ".done").write_text(
        f"accession={accession}\nfetched={stamp}\n"
        f"md5_verified_files={n_verified}\n")
    zip_path.unlink(missing_ok=True)
    # Only this accession's staging dir — wiping the shared _staging root
    # would delete a concurrent fetch's files mid-verification.
    shutil.rmtree(staging, ignore_errors=True)
    if not quiet:
        size_gb = sum(p.stat().st_size for p in dest.rglob("*")) / 1e9
        print(f"  {accession}: {size_gb:.2f} GB on disk, "
              f"{n_verified} files md5-verified")
    return dest


def verify_one(accession: str) -> tuple[bool, str]:
    """Re-check an installed genome against its own kept md5sum.txt."""
    dest = genome_dir(accession)
    if not is_done(accession):
        return False, "not fetched"
    md5file = dest / "md5sum.txt"
    if not md5file.exists():
        return False, "no md5sum.txt kept (fetched before --verify existed)"
    n = 0
    for line in md5file.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        want, name = parts
        target = dest / name
        if not target.exists():
            return False, f"missing file: {name}"
        if _md5(target) != want:
            return False, f"md5 mismatch: {name}"
        n += 1
    return (n > 0), (f"{n} files re-verified" if n else "nothing to verify")


def purge_one(accession: str) -> None:
    """Delete one genome dir (guarded to stay inside data_root/genomes)."""
    dest = genome_dir(accession).resolve()
    root = genomes_root().resolve()
    if root not in dest.parents:
        raise RuntimeError(f"refusing to purge outside genomes root: {dest}")
    if dest.exists():
        shutil.rmtree(dest)
        print(f"  {accession}: purged")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--only", help="comma-separated accessions")
    ap.add_argument("--limit", type=int, help="stop after N genomes")
    ap.add_argument("--max-gb", type=float,
                    help="stop before the cumulative FASTA estimate exceeds this")
    ap.add_argument("--smallest-first", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="re-check installed genomes' md5s, then exit")
    ap.add_argument("--purge", help="comma-separated accessions to delete")
    ap.add_argument("--search-cmd",
                    help="shell command per fetched genome; {fna} {acc} {dir} expand")
    ap.add_argument("--delete-after-search", action="store_true",
                    help="purge each genome after its --search-cmd succeeds")
    args = ap.parse_args()

    if args.purge:
        for acc in args.purge.split(","):
            purge_one(acc.strip())
        return 0
    if args.delete_after_search and not args.search_cmd:
        ap.error("--delete-after-search requires --search-cmd")

    rows = read_manifest(args.manifest)
    if args.only:
        keep = {a.strip() for a in args.only.split(",")}
        rows = [r for r in rows if r["accession"] in keep]
        missing = keep - {r["accession"] for r in rows}
        if missing:
            ap.error(f"not in manifest: {sorted(missing)}")

    if args.verify:
        bad = 0
        for r in rows:
            if not is_done(r["accession"]):
                continue
            ok, msg = verify_one(r["accession"])
            print(f"  {'ok ' if ok else 'BAD'} {r['accession']} {msg}")
            bad += 0 if ok else 1
        print(f"{bad} genome(s) failed verification")
        return 1 if bad else 0

    if args.smallest_first:
        rows.sort(key=lambda r: int(r["total_length_bp"]))

    pending = [r for r in rows if not is_done(r["accession"])]
    print(f"{len(rows)} selected, {len(rows) - len(pending)} already done, "
          f"{len(pending)} to fetch")
    done_count, cum_gb = 0, 0.0
    for r in pending:
        acc, est = r["accession"], float(r["fasta_gb"])
        if args.limit and done_count >= args.limit:
            break
        if args.max_gb and cum_gb + est > args.max_gb:
            print(f"  stopping: {acc} ({est:.2f} GB) would exceed --max-gb")
            break
        print(f"[{done_count + 1}] {acc} {r['organism']} (~{est:.2f} GB)")
        if args.dry_run:
            done_count, cum_gb = done_count + 1, cum_gb + est
            continue
        fetch_one(acc)
        done_count, cum_gb = done_count + 1, cum_gb + est
        if args.search_cmd:
            cmd = args.search_cmd.format(fna=fna_path(acc), acc=acc,
                                         dir=genome_dir(acc))
            print(f"  running: {cmd}")
            rc = subprocess.run(cmd, shell=True).returncode
            if rc != 0:
                print(f"  ! search command failed (rc={rc}); genome kept")
                continue
            if args.delete_after_search:
                purge_one(acc)
    print(f"Fetched {done_count} genomes (~{cum_gb:.1f} GB est. FASTA)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
