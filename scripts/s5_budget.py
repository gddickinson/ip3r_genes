"""s5_budget.py — project the full sweep's cost from the pilot's measured rates.

S5a's job is to hand S5b a budget that came off a stopwatch rather than a
guess. Rates are measured on clean recomputes (a cached `miniprot.gff` makes a
rerun look free, which is exactly how a budget gets written from a number that
means nothing), and the projection is applied to the genomes the manifest
still holds.

Outputs -> results/genome_ledger/s5b_budget.{tsv,json}

Usage:
  python scripts/s5_budget.py
  python scripts/s5_budget.py --workers 3
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from src.utils.data_root import get_data_root, free_bytes   # noqa: E402
import s5_sweep_lib as lib                                  # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "genome_ledger"
MANIFEST = PROJECT_ROOT / "results" / "genome_manifest.tsv"

#: Measured on clean recomputes in the S5a pilot: miniprot with 8 threads
#: over the 121,294-residue panel ran at 36 s/Gbp on both a 0.85 Gbp bird
#: (Todus, 30.6 s) and a 1.86 Gbp mammal (Ornithorhynchus, 67.5 s). Stated
#: here rather than re-derived so a rerun that drifts is visible.
MINIPROT_S_PER_GBP = 36.0

#: tblastn rescue only fires on a cell miniprot left empty, so its cost is
#: per *genome that needs it*, not per Gbp: 22.3 s for two empty cells on a
#: 0.85 Gbp assembly, including the makeblastdb build.
RESCUE_S_PER_GENOME = 30.0

#: Pilot fraction of genomes needing any rescue at all (2 of 6). The margin
#: set is dominated by fragmented bird assemblies, so the full sweep should
#: expect this to be higher, not lower.
RESCUE_FRACTION = 0.5


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=2,
                    help="parallel shards assumed for the wall-clock estimate")
    ap.add_argument("--max-fasta-gb", type=float, default=12.0,
                    help="the driver's default cut; giants are budgeted apart")
    args = ap.parse_args()

    manifest = read_tsv(MANIFEST)
    sweep_root = get_data_root() / "genome_sweep"
    done = {d.name for d in sweep_root.glob("*/") if (d / ".sweep.done").exists()}

    groups = {
        "swept (S5a pilot)": [r for r in manifest if r["accession"] in done],
        "remaining, standard": [r for r in manifest if r["accession"] not in done
                                and float(r["fasta_gb"]) <= args.max_fasta_gb],
        "remaining, giants": [r for r in manifest if r["accession"] not in done
                              and float(r["fasta_gb"]) > args.max_fasta_gb],
    }

    rows, totals = [], {}
    for label, rs in groups.items():
        gbp = sum(int(r["total_length_bp"]) for r in rs) / 1e9
        zip_gb = sum(float(r["est_zip_gb"]) for r in rs)
        fasta_gb = sum(float(r["fasta_gb"]) for r in rs)
        miniprot_h = gbp * MINIPROT_S_PER_GBP / 3600
        rescue_h = len(rs) * RESCUE_FRACTION * RESCUE_S_PER_GENOME / 3600
        rows.append({
            "group": label, "genomes": len(rs), "gbp": round(gbp, 1),
            "download_gb": round(zip_gb, 1), "fasta_gb": round(fasta_gb, 1),
            "miniprot_h": round(miniprot_h, 1),
            "rescue_h": round(rescue_h, 1),
            "compute_h": round(miniprot_h + rescue_h, 1),
            "wall_h_at_workers": round((miniprot_h + rescue_h) / args.workers, 1),
        })
        if label != "swept (S5a pilot)":
            for k in ("genomes", "gbp", "download_gb", "fasta_gb", "compute_h"):
                totals[k] = round(totals.get(k, 0) + rows[-1][k], 1)

    span = lib.itpr_span_stats()

    def below(r: dict) -> bool:
        return int(r["contig_n50"] or 0) < span["median"]

    below_bar = [r for r in manifest if below(r)]

    # D4's bar broken out, because *which* genomes fail it is the result.
    # The margin species were chosen for what their proteomes lack; if they
    # are also the assemblies that cannot carry the gene, the two signals are
    # confounded and the ledger has to say so before S15 reads any absence.
    contiguity: list[dict] = []
    seen: dict[str, list] = {}
    for r in manifest:
        seen.setdefault(r["vclass"], []).append(r)
    for vclass, rs in sorted(seen.items(), key=lambda kv: -len(kv[1])):
        n = sum(1 for r in rs if below(r))
        contiguity.append({"group": "class", "name": vclass, "below": n,
                           "total": len(rs),
                           "pct": round(100 * n / len(rs), 1)})
    for name, rs in (("margin species",
                      [r for r in manifest
                       if not r["reasons"].startswith("order_rep")]),
                     ("order reps",
                      [r for r in manifest
                       if r["reasons"].startswith("order_rep")])):
        n = sum(1 for r in rs if below(r))
        contiguity.append({"group": "scope", "name": name, "below": n,
                           "total": len(rs),
                           "pct": round(100 * n / len(rs), 1)})
    contiguity.append({"group": "total", "name": "all manifest genomes",
                       "below": len(below_bar), "total": len(manifest),
                       "pct": round(100 * len(below_bar) / len(manifest), 1)})
    write_tsv(OUT_DIR / "contiguity_bar.tsv",
              ["group", "name", "below", "total", "pct"], contiguity)
    capped = [r for r in manifest
              if lib.intron_rule_capped(int(r["total_length_bp"] or 0))]
    free_gb = free_bytes(get_data_root()) / 1e9

    stats = {
        "rates": {"miniprot_s_per_gbp": MINIPROT_S_PER_GBP,
                  "rescue_s_per_genome": RESCUE_S_PER_GENOME,
                  "rescue_fraction": RESCUE_FRACTION,
                  "measured_on": "S5a pilot clean recomputes, 8 threads, "
                                 "121,294-residue panel"},
        "groups": rows, "remaining_totals": totals,
        "workers_assumed": args.workers,
        "wall_h_remaining": round(totals.get("compute_h", 0) / args.workers, 1),
        "storage": {
            "free_gb_now": round(free_gb, 1),
            "fasta_gb_if_kept": totals.get("fasta_gb", 0),
            "fits_without_delete_after": free_gb > totals.get("fasta_gb", 0) * 1.2,
            "note": "--delete-after keeps the peak footprint at ~2 genomes; "
                    "keeping the FASTAs costs the fasta_gb total but saves "
                    "re-downloading for S8 synteny and S10 validation",
        },
        "flags": {
            "genomes_below_contiguity_bar": len(below_bar),
            "contiguity_breakdown": contiguity,
            "contiguity_bar_bp": span["median"],
            "genomes_with_capped_max_intron": len(capped),
            "capped_accessions": [r["accession"] for r in capped],
            "giants_needing_chunked_miniprot": [
                {"accession": r["accession"], "organism": r["organism"],
                 "fasta_gb": r["fasta_gb"]} for r in groups["remaining, giants"]],
        },
    }

    write_tsv(OUT_DIR / "s5b_budget.tsv",
              ["group", "genomes", "gbp", "download_gb", "fasta_gb",
               "miniprot_h", "rescue_h", "compute_h", "wall_h_at_workers"],
              rows)
    (OUT_DIR / "s5b_budget.json").write_text(json.dumps(stats, indent=1) + "\n")

    print(f"S5b budget (rates from the S5a pilot, {args.workers} workers)")
    for r in rows:
        print(f"  {r['group']:<22} {r['genomes']:>4} genomes  "
              f"{r['gbp']:>6.1f} Gbp  {r['download_gb']:>6.1f} GB download  "
              f"{r['compute_h']:>5.1f} h compute")
    print(f"  remaining total: {totals.get('genomes')} genomes, "
          f"{totals.get('compute_h')} h compute -> "
          f"{stats['wall_h_remaining']} h wall at {args.workers} workers")
    print(f"  storage: {totals.get('fasta_gb')} GB of FASTA if kept, "
          f"{free_gb:.0f} GB free "
          f"({'fits' if stats['storage']['fits_without_delete_after'] else 'use --delete-after'})")
    print(f"  {len(below_bar)}/{len(manifest)} manifest genomes are below the "
          f"contiguity bar ({span['median']:,} bp) — their absences are not "
          "claimable (D4)")
    for c in contiguity:
        if c["group"] != "class" or c["below"] == 0:
            continue
        print(f"      {c['name']:<18} {c['below']:>3}/{c['total']:<4} "
              f"{c['pct']:>5.0f}%")
    for c in contiguity:
        if c["group"] == "scope":
            print(f"      {c['name']:<18} {c['below']:>3}/{c['total']:<4} "
                  f"{c['pct']:>5.0f}%")
    print(f"  {len(capped)} genome(s) hit the -G cap: "
          + ", ".join(r["organism"] for r in capped))
    print(f"wrote {OUT_DIR / 's5b_budget.tsv'}")


if __name__ == "__main__":
    main()
