"""S20 shared helpers: paths, group manifests, lineage, TSV plumbing.

Everything structural that more than one S20 script needs. The heavy lifting
is reused from S3 — `s3_hmm_lib` for FASTA/domtblout/TSV, `s3_assign` for the
margin call — because S20's whole point is that the *same* instrument is
pointed at a different part of the tree. A second implementation of the call
would make the range result incomparable with S3's.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_hmm_lib import HMM_SWEEP_DIR, read_tsv, write_tsv  # noqa: E402,F401
from src.utils.data_root import require_data_root  # noqa: E402

S20_DIR = PROJECT_ROOT / "results" / "s20_sweep"
CENSUS_V5_DIR = PROJECT_ROOT / "results" / "census_v5"
FIG_DIR = S20_DIR / "figures"

MANIFEST_FIELDS = ["upid", "organism", "taxid", "protein_count",
                   "group", "genus", "genus_n_available", "status"]


def log(tag: str, msg: str) -> None:
    print(f"[{tag}] {msg}", flush=True)


def s20_dirs() -> tuple[Path, Path]:
    """(committed results dir, bulk raw dir under the data root)."""
    S20_DIR.mkdir(parents=True, exist_ok=True)
    raw = require_data_root() / "hmmer" / "s20"
    raw.mkdir(parents=True, exist_ok=True)
    return S20_DIR, raw


def group_paths(group: str) -> dict[str, Path]:
    """Every path one group's data lives at."""
    root = require_data_root() / "proteomes"
    return {
        "manifest": S20_DIR / f"proteome_manifest_{group}.tsv",
        "dir": root / group,
        "db": root / f"{group}_refprot.fasta",
        "stats": S20_DIR / f"proteome_db_stats_{group}.json",
        "unavailable": S20_DIR / f"proteome_unavailable_{group}.tsv",
    }


def load_group_manifest(group: str) -> list[dict]:
    """The swept proteome list for a group, ints restored."""
    rows = read_tsv(group_paths(group)["manifest"])
    for r in rows:
        r["taxid"] = int(r["taxid"])
        r["protein_count"] = int(r["protein_count"])
        r["genus_n_available"] = int(r["genus_n_available"] or 0)
    return rows


def swept_manifest(groups) -> list[dict]:
    """Every proteome actually swept, across groups (status == 'swept')."""
    out = []
    for g in groups:
        if group_paths(g)["manifest"].exists():
            out += [r for r in load_group_manifest(g) if r["status"] == "swept"]
    return out


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=str) + "\n")


def read_json(path: Path):
    return json.loads(path.read_text()) if path.exists() else None


def live(steps: list[tuple[str, bool]], extra: dict | None = None) -> None:
    """Dashboard live panel (session protocol)."""
    payload = {"task": "S20",
               "steps": [{"label": l, "done": d} for l, d in steps]}
    if extra:
        payload.update(extra)
    (PROJECT_ROOT / "results").mkdir(exist_ok=True)
    (PROJECT_ROOT / "results" / "session_live.json").write_text(
        json.dumps(payload, indent=1))


# ------------------------------------------- exact hit → proteome attribution
def _scan_headers(args) -> tuple[str, list[str]]:
    """Header-only pass over one proteome gz → (upid, accessions it holds).

    Module-level and picklable so a process pool can run one per worker;
    gzip decompression does not release the GIL usefully, so threads would
    not help here.
    """
    import gzip
    path, wanted = args
    upid = path.name.split("_")[0]
    found: list[str] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[:1] != ">":
                continue
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) >= 3 else line[1:].split()[0]
            if acc in wanted:
                found.append(acc)
    return upid, found


def accession_to_upid(group: str, wanted: set[str]) -> dict[str, str]:
    """acc → the reference proteome the sequence actually came from.

    The sweep runs against one concatenated DB, which loses which proteome
    each sequence belongs to, and the obvious substitute — matching on the
    header's taxid — is wrong: 3 protist, 18 plant and 27 fungal taxids carry
    **two** reference proteomes each, so a taxid key credits both proteomes
    with either one's hits and inflates the presence denominator's numerator.
    So the attribution is measured, by a header-only pass over each proteome
    file, and cached under the data root because it is expensive and fixed.
    """
    from concurrent.futures import ProcessPoolExecutor
    cache = (require_data_root() / "hmmer" / "s20" /
             f"acc2upid_{group}.tsv")
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        cached = {r["accession"]: r["upid"] for r in read_tsv(cache)}
        if wanted <= set(cached):
            return {a: cached[a] for a in wanted}
    files = sorted(group_paths(group)["dir"].glob("*.fasta.gz"))
    log("s20_lib", f"{group}: attributing {len(wanted)} hits across "
                   f"{len(files)} proteome files")
    out: dict[str, str] = {}
    with ProcessPoolExecutor() as pool:
        for upid, accs in pool.map(_scan_headers,
                                   ((p, wanted) for p in files),
                                   chunksize=8):
            for a in accs:
                out.setdefault(a, upid)
    write_tsv(cache, ["accession", "upid"],
              [{"accession": a, "upid": u} for a, u in sorted(out.items())])
    log("s20_lib", f"{group}: {len(out)}/{len(wanted)} hits attributed")
    return out
