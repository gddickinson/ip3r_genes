"""s4_build_manifest.py — build results/genome_manifest.tsv, the S4 denominator.

Scope rule (roadmap S4, confirmed with the user 2026-09-04): one best
reference assembly per vertebrate order, UNION every margin species the S3
census leaves undecided — derived from the committed tables by
`s4_manifest_lib.derive_margins`, never hand-listed.

Every `datasets` call is cached under <data_root>/raw_api/ncbi_datasets/, so
a rerun is offline unless --refresh-assemblies is passed. The binary is not
on the bare PATH (S1's toolchain manifest records it as env-resident), so
its location is resolved rather than assumed.

Usage:
    python3 scripts/s4_build_manifest.py [--refresh-assemblies]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.data_root import get_data_root, free_bytes  # noqa: E402
from src.utils.family import MIN_LENGTH_AA  # noqa: E402
from scripts.s4_manifest_lib import (  # noqa: E402
    ANCHORS, EXPECTED_PARALOGS, ZIP_FACTOR, class_label, datasets_bin,
    derive_margins, load_jsonl, load_taxonomy, merge_rows, parse_assembly,
    rank_key, select_order_reps, write_manifest,
)
from scripts.s4_notes import write_notes  # noqa: E402

VERTEBRATA_TAXID = "7742"
TAXONOMY_CHUNK = 400


def run_datasets(args: list[str]) -> str:
    proc = subprocess.run([datasets_bin()] + args,
                          capture_output=True, text=True)
    if proc.returncode != 0 and not proc.stdout:
        raise RuntimeError(
            f"datasets {' '.join(args[:4])}… failed: {proc.stderr[:300]}")
    return proc.stdout


def ensure_assembly_dump(cache_dir: Path, refresh: bool) -> Path:
    dump = cache_dir / "vertebrata_reference_genomes.jsonl"
    if dump.exists() and dump.stat().st_size > 0 and not refresh:
        return dump
    print("Fetching all Vertebrata reference-genome summaries…")
    dump.write_text(run_datasets(
        ["summary", "genome", "taxon", VERTEBRATA_TAXID, "--reference",
         "--as-json-lines"]))
    return dump


def ensure_taxonomy(cache_dir: Path, taxids: set[int]) -> Path:
    """Fetch classification for any taxid not already in the cached dump."""
    dump = cache_dir / "vertebrata_taxonomy.jsonl"
    have = set(load_taxonomy(dump))
    missing = sorted(t for t in taxids if t and t not in have)
    if not missing:
        return dump
    print(f"Fetching taxonomy for {len(missing)} taxids…")
    with open(dump, "a") as fh:
        for i in range(0, len(missing), TAXONOMY_CHUNK):
            chunk = [str(t) for t in missing[i:i + TAXONOMY_CHUNK]]
            fh.write(run_datasets(
                ["summary", "taxonomy", "taxon"] + chunk + ["--as-json-lines"]))
            print(f"  {min(i + TAXONOMY_CHUNK, len(missing))}/{len(missing)}")
    return dump


def query_species_assemblies(cache_dir: Path, name: str,
                             taxid: int | None) -> list[dict]:
    """Every assembly for one margin species, cached per species.

    The reference dump only carries NCBI's chosen reference per species, so a
    margin species with no reference genome is invisible there. This asks for
    all of them — which is the difference between "no genome" and "no
    *reference* genome", and only the first is grounds to drop a species from
    the denominator.
    """
    mdir = cache_dir / "margins"
    mdir.mkdir(exist_ok=True)
    cache = mdir / f"{taxid or name.lower().replace(' ', '_')}.jsonl"
    if not cache.exists():
        try:
            out = run_datasets(["summary", "genome", "taxon",
                                str(taxid) if taxid else name,
                                "--as-json-lines"])
        except RuntimeError as exc:
            print(f"  ! no assemblies for {name}: {exc}")
            out = ""
        cache.write_text(out)
    if not cache.stat().st_size:
        return []
    return [parse_assembly(r) for r in load_jsonl(cache)]


def resolve_margins(cache_dir: Path, margins: list[dict],
                    assemblies: list[dict]):
    """Match each margin species to its best assembly."""
    # Best assembly per taxid, by the *same* rank function the order reps
    # use. A plain last-wins dict lets a species that has two reference
    # assemblies resolve to a different one here than it did there, which
    # puts the same species in the manifest twice under two accessions.
    by_taxid: dict[int, dict] = {}
    for a in assemblies:
        cur = by_taxid.get(a["taxid"])
        if cur is None or rank_key(a) > rank_key(cur):
            by_taxid[a["taxid"]] = a
    resolved, unresolved = [], []
    for m in margins:
        asm = by_taxid.get(m["taxid"])
        if asm is None:
            cands = [a for a in assemblies
                     if a["organism"] == m["organism"]
                     or a["organism"].startswith(m["organism"] + " ")]
            asm = max(cands, key=rank_key) if cands else None
        if asm is None:
            cands = query_species_assemblies(cache_dir, m["organism"],
                                             m["taxid"])
            asm = max(cands, key=rank_key) if cands else None
        if asm is None:
            unresolved.append(m)
            print(f"  ! margin species without any assembly: {m['organism']}")
        else:
            resolved.append(dict(asm, reasons=m["reasons"], notes=m["notes"]))
    return resolved, unresolved


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh-assemblies", action="store_true")
    args = ap.parse_args()

    cache_dir = get_data_root() / "raw_api" / "ncbi_datasets"
    cache_dir.mkdir(parents=True, exist_ok=True)
    results_dir = PROJECT_ROOT / "results"

    dump = ensure_assembly_dump(cache_dir, args.refresh_assemblies)
    assemblies = [parse_assembly(r) for r in load_jsonl(dump)]
    print(f"Reference assemblies in the dump: {len(assemblies):,}")

    margins = derive_margins(PROJECT_ROOT, MIN_LENGTH_AA)
    print(f"Margin species derived from the census: {len(margins)}")
    margin_rows, unresolved = resolve_margins(cache_dir, margins, assemblies)

    taxids = ({a["taxid"] for a in assemblies}
              | {m["taxid"] for m in margin_rows})
    taxmap = load_taxonomy(ensure_taxonomy(cache_dir, taxids))

    order_reps, no_order = select_order_reps(assemblies, taxmap)
    print(f"Vertebrate orders with a reference assembly: {len(order_reps)}")

    for m in margin_rows:  # margins need the class/order columns too
        tax = taxmap.get(m["taxid"], {})
        m.setdefault("vorder", tax.get("order", ""))
        m.setdefault("vclass", class_label(tax, m["vorder"]))

    rows = merge_rows(order_reps, margin_rows)
    out_tsv = results_dir / "genome_manifest.tsv"
    totals = write_manifest(rows, out_tsv)
    free_gb = free_bytes(get_data_root()) / 1e9
    write_notes(results_dir / "genome_manifest_notes.md", totals, rows,
                order_reps, no_order, unresolved, margin_rows, margins,
                free_gb, {"expected_paralogs": EXPECTED_PARALOGS,
                          "min_length_aa": MIN_LENGTH_AA,
                          "zip_factor": ZIP_FACTOR, "anchors": ANCHORS})
    # Dashboard live panel (session protocol): S4's product is the manifest
    # plus a fetch path proven end to end, so both are reported.
    fetched = len(list((get_data_root() / "genomes").glob("*/.done")))
    (results_dir / "session_live.json").write_text(json.dumps({
        "task": "S4",
        "steps": [
            {"label": f"assembly dump ({len(assemblies):,} vertebrate "
                      "reference genomes)", "done": True},
            {"label": f"margin species derived from the census "
                      f"({len(margins)}, {len(unresolved)} unresolved)",
             "done": True},
            {"label": f"manifest written ({totals['n']} genomes, "
                      f"{totals['total_fasta_gb']:.0f} GB FASTA)",
             "done": True},
            {"label": f"genomes fetched + md5-verified ({fetched}/"
                      f"{totals['n']})", "done": fetched >= totals["n"]},
        ],
    }, indent=1))
    print(f"Wrote {out_tsv} — {totals['n']} genomes; "
          f"FASTA ≈ {totals['total_fasta_gb']:.0f} GB, "
          f"zip ≈ {totals['total_zip_gb']:.0f} GB, "
          f"free on the data root {free_gb:.0f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
