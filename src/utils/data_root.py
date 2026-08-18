"""Bulk-data root resolution.

Publication-grade work needs genome downloads, proteome sets, BLAST
databases, and structure files — tens to hundreds of GB that must not live
in the git repo. The active bulk-storage location is whatever path is in
`data_root.txt` at the project root; sessions read it through
`get_data_root()` and never hard-code storage paths.

To move bulk storage (e.g. onto a different external drive): edit
`data_root.txt` to the new location and re-run any session — the standard
subdirectories are created on demand.

If the configured path is unavailable (drive unplugged) `get_data_root()`
falls back to `<project>/data` and says so, rather than dying mid-session.
That fallback is fine for a few MB of cached API pages and fatal for a
400 GB genome download onto a laptop disk, so anything bulk calls
`require_data_root()` instead, which raises. The session protocol runs
`python -m src.utils.data_root --require` at the start of every session
for exactly this reason.
"""

from __future__ import annotations

from pathlib import Path

SUBDIRS = [
    "genomes",          # assembly FASTAs for the tblastn/miniprot sweep
    "proteomes",        # UniProt reference proteomes for hmmsearch
    "blast_db",         # formatted BLAST/DIAMOND databases
    "hmmer",            # profile HMMs + search outputs
    "structures",       # predicted/experimental structure files
    "raw_api",          # raw API dumps (InterPro pages, Compara JSON)
    "results_archive",  # bundles too large for the repo
    "sra",              # RNA-seq evidence downloads
]


def get_data_root(project_root: Path | None = None, create: bool = True) -> Path:
    """Resolve the bulk-data root; optionally ensure subdirs exist."""
    project_root = project_root or Path(__file__).resolve().parents[2]
    cfg = project_root / "data_root.txt"
    root = project_root / "data"
    if cfg.exists():
        configured = Path(cfg.read_text().strip()).expanduser()
        parent_ok = configured.exists() or configured.parent.exists()
        if parent_ok:
            root = configured
        else:
            print(f"[data_root] configured path unavailable ({configured}) — "
                  f"falling back to {root}")
    if create:
        try:
            for sub in SUBDIRS:
                (root / sub).mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise RuntimeError(
                f"data root {root} is not writable ({e}) — edit data_root.txt"
            ) from e
    return root


def describe(project_root: Path | None = None) -> str:
    root = get_data_root(project_root, create=False)
    exists = root.exists()
    return f"data root: {root}  ({'available' if exists else 'NOT AVAILABLE'})"


def require_data_root(project_root: Path | None = None) -> Path:
    """The bulk-data root, or a hard failure.

    Use this before any download or sweep. The point is that an unplugged
    drive stops the session at the top rather than after 200 GB has landed
    on the internal disk.
    """
    project_root = project_root or Path(__file__).resolve().parents[2]
    cfg = project_root / "data_root.txt"
    if not cfg.exists():
        raise RuntimeError(
            f"no data_root.txt at {project_root} — bulk storage is not "
            f"configured; see PUBLICATION_ROADMAP.md 'Storage'")
    configured = Path(cfg.read_text().strip()).expanduser()
    if not (configured.exists() or configured.parent.exists()):
        raise RuntimeError(
            f"configured data root is unavailable: {configured}\n"
            f"Attach the drive (or edit data_root.txt) before running "
            f"anything that downloads or sweeps bulk data.")
    return get_data_root(project_root)


def free_bytes(path: Path) -> int:
    """Free space at `path`, for the download-budget checks in the roadmap."""
    import shutil
    return shutil.disk_usage(path).free


if __name__ == "__main__":
    import sys

    strict = "--require" in sys.argv
    if strict:
        try:
            root = require_data_root()
        except RuntimeError as exc:
            print(f"[data_root] {exc}")
            sys.exit(1)
    else:
        root = get_data_root()
    print(describe())
    for sub in SUBDIRS:
        print(f"  {root / sub}")
    print(f"  free: {free_bytes(root) / 1e9:.1f} GB")
