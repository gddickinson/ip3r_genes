"""S3 — calibrate the two profiles before the sweep uses them.

The census-v3 sweep assigns proteins to ITPR or RYR by which profile scores
them higher. That instrument has to be measured before it is trusted, and
S2 left the perfect ruler: 15,421 records whose call came from a *different*
kind of evidence — the presence and absence of Pfam signatures (D14b) —
audited there against 6,191 gene symbols with 0 disagreements.

So: run both profiles over S2's archived seeded-space FASTA, assign by
margin, and score the assignment against S2's architecture call. Two
independent instruments agreeing on 15,000 proteins is what lets the sweep's
novel hits — which have no architecture call, because they were never in
InterPro's result set — be assigned by profile alone.

The seeds are inside that set by construction, so the agreement is reported
twice: over everything, and with the 56 seed accessions removed. The second
number is the one that means anything.

Outputs (results/hmm_sweep/):
  calibration_assignments.tsv    every record: both scores, margin, both calls
  calibration_disagreements.tsv  every row where profile ≠ architecture
  calibration_summary.json       confusion counts, seed-excluded and not

Run:  python3 scripts/s3_calibrate.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_assign import (  # noqa: E402
    ASSIGN_FIELDS, MIN_SCORE, REL_MARGIN, assign, audit_against,
    load_profile_hits,
)
from scripts.s3_hmm_lib import (  # noqa: E402
    HMM_SWEEP_DIR, census_index, load_census_v2, read_tsv, write_tsv,
)
from src.utils.data_root import require_data_root  # noqa: E402

SEEDED_FASTA = "raw_api/uniprot/s2_seed_sweep.fasta"
PROFILE_NAMES = ("itpr", "ryr")


def log(msg: str) -> None:
    print(f"[s3_calib] {msg}", flush=True)


def search(name: str, db: Path, raw_dir: Path, cpu: int, evalue: str) -> Path:
    domtbl = raw_dir / f"calibration_{name}_vs_s2seeded.domtblout"
    if domtbl.exists():
        log(f"reusing {domtbl.name}")
        return domtbl
    logf = raw_dir / f"calibration_{name}_vs_s2seeded.log"
    t0 = time.time()
    proc = subprocess.run(
        ["hmmsearch", "--cpu", str(cpu), "-E", evalue, "--noali",
         "--domtblout", str(domtbl), "-o", str(logf),
         str(HMM_SWEEP_DIR / f"{name}.hmm"), str(db)],
        capture_output=True, text=True)
    log(f"hmmsearch {name}.hmm rc={proc.returncode} in {time.time() - t0:.0f}s")
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-2000:])
        raise SystemExit("hmmsearch failed")
    return domtbl


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--evalue", default="1e-5")
    args = ap.parse_args()

    db = require_data_root() / SEEDED_FASTA
    if not db.exists():
        raise SystemExit(f"S2 seeded-space FASTA missing: {db}")
    raw_dir = require_data_root() / "hmmer"
    raw_dir.mkdir(parents=True, exist_ok=True)

    domtbls = {n: search(n, db, raw_dir, args.cpu, args.evalue)
               for n in PROFILE_NAMES}
    rows = assign(load_profile_hits(domtbls["itpr"]),
                  load_profile_hits(domtbls["ryr"]))
    log(f"{len(rows)} records scored by both profiles")

    census = census_index(load_census_v2())
    # Two rulers, not one. `arch_call` is the pure architecture verdict of
    # D14b — the one that reads only which signatures a record carries.
    # `call` is that verdict with S2's low-confidence gene-symbol fallback
    # folded in for records the architecture could not call. The profiles
    # are scored against both, because agreeing with the architecture rule
    # and agreeing with a gene symbol are different claims.
    truth = {acc: r["arch_call"] for acc, r in census.items()}
    truth_call = {acc: r["call"] for acc, r in census.items()}
    seed_accs = {r["accession"] for r in
                 read_tsv(HMM_SWEEP_DIR / "seed_manifest.tsv")}

    # Carry S2's call and the record's own labels onto every row so the
    # committed table can be re-audited without re-joining anything.
    for row in rows:
        c = census.get(row["accession"], {})
        row["s2_call"] = c.get("call", "")
        row["s2_arch_call"] = c.get("arch_call", "")
        row["s2_reason"] = c.get("reason", "")
        row["s2_arch"] = c.get("n_itpr_arch", "")
        row["s2_group"] = c.get("group", "")
        row["is_seed"] = int(row["accession"] in seed_accs)

    non_seed_rows = [r for r in rows if not r["is_seed"]]
    full = audit_against(rows, truth, "census v2 architecture call")
    non_seed = audit_against(non_seed_rows, truth,
                             "census v2 architecture call, seeds excluded")
    vs_call = audit_against(non_seed_rows, truth_call,
                            "census v2 call incl. symbol fallback, "
                            "seeds excluded")

    fields = ASSIGN_FIELDS + ["s2_call", "s2_arch_call", "s2_reason",
                              "s2_arch", "s2_group", "is_seed"]
    write_tsv(HMM_SWEEP_DIR / "calibration_assignments.tsv", fields, rows)
    # One disagreements table, carrying both rulers' verdicts, so a row can
    # be read without re-joining anything.
    disagreements = {r["accession"]: r for r in non_seed["disagreements"]}
    for r in vs_call["disagreements"]:
        disagreements.setdefault(r["accession"], r)
    write_tsv(HMM_SWEEP_DIR / "calibration_disagreements.tsv",
              fields + ["truth"],
              sorted(disagreements.values(), key=lambda r: r["accession"]))

    # Where do the profiles decline to call, and where does S2?
    prof_by_s2 = Counter((r["s2_arch_call"] or "not in v2", r["assignment"])
                         for r in rows)
    unassigned_recovered = sum(
        1 for r in rows if r["s2_arch_call"] == "unassigned"
        and r["assignment"] != "unassigned")
    # Records S2 could only call from a gene symbol, which the profiles now
    # decide on sequence — and how often the profile contradicts the symbol.
    symbol_only = [r for r in rows if r["s2_arch_call"] == "unassigned"
                   and r["s2_call"] in ("ITPR", "RYR")]
    symbol_overturned = [r for r in symbol_only
                         if r["assignment"] in ("ITPR", "RYR")
                         and r["assignment"] != r["s2_call"]]
    summary = {
        "rel_margin": REL_MARGIN, "min_score": MIN_SCORE,
        "evalue": args.evalue,
        "records_scored": len(rows),
        "census_v2_records": len(census),
        "seeds_in_set": sum(r["is_seed"] for r in rows),
        "all_records": full["counts"],
        "seeds_excluded": non_seed["counts"],
        "vs_call_with_symbol_fallback": vs_call["counts"],
        "confusion_s2_vs_profile": {f"{a} → {b}": n
                                    for (a, b), n in sorted(prof_by_s2.items())},
        "s2_arch_unassigned_called_by_profile": unassigned_recovered,
        "symbol_fallback_records": len(symbol_only),
        "symbol_fallback_overturned": len(symbol_overturned),
        "profile_calls": dict(Counter(r["assignment"] for r in rows)),
    }
    (HMM_SWEEP_DIR / "calibration_summary.json").write_text(
        json.dumps(summary, indent=1))
    log(f"agreement (seeds excluded): {non_seed['counts']['agreement']} "
        f"({non_seed['counts']['agree']} agree, "
        f"{non_seed['counts']['disagree']} disagree, "
        f"{non_seed['counts']['profile_unassigned']} profile-unassigned)")
    log(f"agreement vs the symbol-fallback call: "
        f"{vs_call['counts']['agreement']} "
        f"({vs_call['counts']['disagree']} disagree)")
    log(f"architecture-unassigned records the profiles call: "
        f"{unassigned_recovered}")
    log(f"symbol-fallback calls the profiles overturn: "
        f"{len(symbol_overturned)} of {len(symbol_only)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
