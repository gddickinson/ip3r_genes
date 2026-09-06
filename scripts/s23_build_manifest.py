"""s23_build_manifest.py — build results/s23_scope/genome_manifest_s23.tsv.

The declared denominator for the non-vertebrate sweep: G1-G5 of
`s23_scope.py` applied to NCBI's eukaryote reference assemblies, with the size
guard from `s23_calibration.py`.

Every `datasets` call is cached under <data_root>/raw_api/ncbi_datasets/, so a
rerun is offline unless --refresh is passed. The binary is env-resident (S1's
toolchain manifest), so its location is resolved rather than assumed.

Usage:
    python3 scripts/s23_build_manifest.py
    python3 scripts/s23_build_manifest.py --refresh
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from src.utils.data_root import free_bytes, get_data_root      # noqa: E402
import s23_calibration as cal                                  # noqa: E402
import s23_manifest_lib as mlib                                # noqa: E402
import s23_scope as scope                                      # noqa: E402
from s23_scope import read_tsv                                 # noqa: F401,E402
from s4_manifest_lib import datasets_bin, load_jsonl, parse_assembly  # noqa: E402

EUKARYOTA_TAXID = "2759"
TAXONOMY_CHUNK = 400
OUT_DIR = PROJECT_ROOT / "results" / "s23_scope"


def cache_dir() -> Path:
    d = get_data_root() / "raw_api" / "ncbi_datasets"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_datasets(args: list[str]) -> str:
    proc = subprocess.run([datasets_bin()] + args, capture_output=True,
                          text=True)
    if proc.returncode != 0 and not proc.stdout:
        raise RuntimeError(
            f"datasets {' '.join(args[:4])}… failed: {proc.stderr[:300]}")
    return proc.stdout


def ensure_assembly_dump(refresh: bool) -> Path:
    dump = cache_dir() / f"s23_taxon_{EUKARYOTA_TAXID}_reference.jsonl"
    if dump.exists() and dump.stat().st_size and not refresh:
        return dump
    print("Fetching every Eukaryota reference-genome summary…")
    dump.write_text(run_datasets(
        ["summary", "genome", "taxon", EUKARYOTA_TAXID, "--reference",
         "--as-json-lines"]))
    return dump


def load_ncbi_taxonomy(path: Path) -> dict[int, dict]:
    """taxid -> {kingdom, phylum, class, order, parents} from a taxonomy dump.

    Richer than `s4_manifest_lib.load_taxonomy`, which keeps only order and
    class: this scope selects at phylum and kingdom rank and decides "is this a
    vertebrate" from the parent taxid list, so all of it is kept.
    """
    out: dict[int, dict] = {}
    if not path.exists():
        return out
    for rec in load_jsonl(path):
        w = rec.get("taxonomy", rec)
        cls = w.get("classification", {})
        entry = {rank: (cls.get(rank) or {}).get("name", "")
                 for rank in ("kingdom", "phylum", "class", "order", "genus",
                              "species")}
        entry["parents"] = set(w.get("parents") or [])
        entry["name"] = (w.get("current_scientific_name") or {}).get("name", "")
        for q in rec.get("query") or []:
            try:
                out[int(q)] = entry
            except (TypeError, ValueError):
                pass
        if w.get("tax_id"):
            out[int(w["tax_id"])] = entry
    return out


def ensure_taxonomy(taxids: set[int]) -> Path:
    dump = cache_dir() / "s23_eukaryota_taxonomy.jsonl"
    have = set(load_ncbi_taxonomy(dump))
    missing = sorted(t for t in taxids if t and t not in have)
    if not missing:
        return dump
    print(f"Fetching taxonomy for {len(missing)} taxids…")
    with open(dump, "a") as fh:
        for i in range(0, len(missing), TAXONOMY_CHUNK):
            chunk = [str(t) for t in missing[i:i + TAXONOMY_CHUNK]]
            fh.write(run_datasets(
                ["summary", "taxonomy", "taxon"] + chunk + ["--as-json-lines"]))
    return dump


def build_rules(by_clade: dict) -> dict:
    """G1-G5 as data, from the committed S20 tables plus the assembly pool."""
    presence = scope.load_presence()
    s20_tax = scope.load_taxonomy()
    big = scope.big_phyla(presence, s20_tax)
    absences = scope.absence_clades(presence, s20_tax)
    copies = scope.copy_number_questions(presence, s20_tax)

    # G1: every phylum the assembly pool actually holds. Derived from the pool
    # rather than from a phylum list, so a phylum NCBI has no eukaryote
    # reference genome for is absent from the denominator by fact, not by
    # omission.
    phyla = sorted({name for (rank, name) in by_clade if rank == "phylum"})

    # G2: the classes of the big phyla, taken from the pool for the same reason.
    class_phylum: dict[str, str] = {}
    for (rank, name), rows in by_clade.items():
        if rank == "class" and rows:
            class_phylum[name] = rows[0]["phylum"]
    big_phylum_classes = sorted(
        (ph, kl) for kl, ph in class_phylum.items() if ph in big)

    return {"phyla": phyla, "big_phyla": big,
            "big_phylum_classes": big_phylum_classes,
            "absence_clades": absences, "anchors": scope.ANCHORS,
            "copy_number": copies}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(cal.summary(), "\n")

    dump = ensure_assembly_dump(args.refresh)
    assemblies = [parse_assembly(json.loads(l)) for l in open(dump)
                  if l.strip()]
    print(f"Eukaryote reference assemblies: {len(assemblies):,}")

    taxmap = load_ncbi_taxonomy(
        ensure_taxonomy({a["taxid"] for a in assemblies}))
    by_clade, by_taxid, by_anc, skipped = mlib.index_assemblies(
        assemblies, taxmap)
    n_vert = sum(1 for s in skipped if s["why"].startswith("vertebrate"))
    print(f"  {n_vert:,} vertebrate (S4's scope), "
          f"{len(skipped) - n_vert:,} without a taxonomy record, "
          f"{len(assemblies) - len(skipped):,} in scope")

    rules = build_rules(by_clade)
    print(f"G1 {len(rules['phyla'])} phyla | "
          f"G2 {len(rules['big_phylum_classes'])} classes in "
          f"{len(rules['big_phyla'])} big phyla | "
          f"G3 {len(rules['absence_clades'])} absence clades | "
          f"G4 {len(rules['anchors'])} anchors | "
          f"G5 {len(rules['copy_number'])} copy-number species")

    picked, unresolved = mlib.select(by_clade, by_taxid, by_anc, rules)
    rows = list(picked.values())
    out_tsv = OUT_DIR / "genome_manifest_s23.tsv"
    totals = mlib.write_manifest(rows, out_tsv)

    reason_counts = Counter(r for row in rows for r in row["reasons"])
    multi = sum(1 for row in rows if len(set(row["reasons"])) > 1)
    free_gb = free_bytes(get_data_root()) / 1e9
    stats = {
        "assemblies_in_pool": len(assemblies),
        "vertebrate_excluded": n_vert,
        "no_taxonomy": len(skipped) - n_vert,
        "genomes": totals["n"],
        "total_gbp": round(totals["total_bp"] / 1e9, 1),
        "fasta_gb": round(totals["total_fasta_gb"], 1),
        "zip_gb": round(totals["total_zip_gb"], 1),
        "free_gb": round(free_gb, 1),
        "below_contiguity_bar": totals["below_contiguity_bar"],
        "by_group": totals["by_group"],
        "by_reason": dict(reason_counts),
        "genomes_with_multiple_reasons": multi,
        "unresolved": len(unresolved),
        "rules": {"big_phyla": rules["big_phyla"],
                  "n_phyla": len(rules["phyla"]),
                  "n_big_phylum_classes": len(rules["big_phylum_classes"]),
                  "n_absence_clades": len(rules["absence_clades"]),
                  "n_anchors": len(rules["anchors"]),
                  "n_copy_number": len(rules["copy_number"])},
        "calibration": cal.calibration(),
    }
    (OUT_DIR / "manifest_stats.json").write_text(json.dumps(stats, indent=1))

    with open(OUT_DIR / "unresolved_targets.tsv", "w") as fh:
        fh.write("rule\ttarget\tnote\n")
        for u in unresolved:
            fh.write(f"{u['rule']}\t{u['target']}\t{u['note']}\n")

    # absence-clade coverage: how many genomes each negative claim rests on
    per_clade: dict[str, int] = defaultdict(int)
    for row in rows:
        for reason, note in zip(row["reasons"], row["notes"]):
            if reason == "absence_clade":
                per_clade[note.split("proteomes of ")[-1].split(" and")[0]] += 1
    with open(OUT_DIR / "absence_coverage.tsv", "w") as fh:
        fh.write("rank\tclade\tswept_proteomes\tgenomes_in_manifest\n")
        for c in rules["absence_clades"]:
            fh.write(f"{c['rank']}\t{c['clade']}\t{c['swept']}\t"
                     f"{per_clade.get(c['clade'], 0)}\n")

    print(f"\nWrote {out_tsv} — {totals['n']} genomes, "
          f"{stats['total_gbp']:,} Gbp (FASTA ≈ {stats['fasta_gb']:.0f} GB, "
          f"zip ≈ {stats['zip_gb']:.0f} GB, free {free_gb:.0f} GB)")
    print("  by rule: " + ", ".join(f"{k} {v}"
                                    for k, v in reason_counts.most_common()))
    print("  by group: " + ", ".join(f"{k} {v}" for k, v in
                                     sorted(totals["by_group"].items())))
    print(f"  {totals['below_contiguity_bar']} genome(s) below their group's "
          "contiguity bar (kept and flagged — D4)")
    if unresolved:
        print(f"  {len(unresolved)} target(s) with no assembly at all "
              "→ unresolved_targets.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
