"""S20 steps 1–2 & 5 — the two-profile sweep over the non-vertebrate groups.

Per group: `hmmsearch` both family profiles over that group's concatenated
reference-proteome DB, then S3's margin assignment (`s3_assign`) on the
result. Identical instrument, identical thresholds, different part of the
tree — that is what makes S3's vertebrate numbers and S20's comparable.

**One search, two sensitivities.** Every search runs at `-E 10` and the
primary (E ≤ 1e-5) call is taken by filtering the same domtblout, because
`-E` is a reporting threshold and does not touch hmmsearch's acceleration
filters. So the strict set is exactly what a strict run would have produced,
the relaxed set is a superset from the same search, and the negative claims
get their stated sensitivity for free rather than from a second experiment
that might have differed in some other way.

**The relaxed panel** (brief step 5) adds the four family Pfam domain models
over the groups where the negative claim lives — a 213-state IP₃-core model
can find something a 2,684-state channel model cannot. See `s20_profiles.py`.

Raw domtblout + logs → `<data_root>/hmmer/s20/`. Committed: per-group
assignments, the per-proteome presence table, and the relaxed-panel hits.

Run:  python3 scripts/s20_run_sweep.py --group fungi
      python3 scripts/s20_run_sweep.py --all --cpu 8
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_groups import ALL_GROUPS, EVALUE_PRIMARY, EVALUE_RELAXED, GROUPS  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    S20_DIR, accession_to_upid, group_paths, live, load_group_manifest, log,
    read_json, read_tsv, s20_dirs, write_json, write_tsv,
)
from scripts.s20_profiles import relaxed_panel  # noqa: E402
from scripts.s3_assign import ASSIGN_FIELDS, assign  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    acc_key, best_hits_by_target, parse_domtblout, parse_uniprot_header,
)

PROFILE_NAMES = ("itpr", "ryr")
SEARCH_EVALUE = EVALUE_RELAXED   # what hmmsearch is told; see module docstring

PRESENCE_FIELDS = [
    "upid", "organism", "taxid", "group", "protein_count",
    "n_itpr", "n_ryr", "n_unassigned", "best_itpr_score", "best_itpr_evalue",
    "best_itpr_acc", "itpr_status",
]


def run(cmd: list[str]) -> float:
    log("s20_sweep", "$ " + " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log("s20_sweep", f"  rc={proc.returncode} in {dt:.0f}s")
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-3000:])
        raise SystemExit(f"hmmsearch failed: {' '.join(cmd[:3])}")
    return dt


def search(hmm: Path, db: Path, out_stem: Path, cpu: int) -> tuple[Path, float]:
    """One profile over one DB at the search threshold; cached if present."""
    domtbl = out_stem.with_suffix(".domtblout")
    if domtbl.exists() and domtbl.stat().st_size > 0:
        log("s20_sweep", f"  reusing {domtbl.name}")
        return domtbl, 0.0
    tmp = domtbl.with_suffix(".domtblout.part")
    dt = run(["hmmsearch", "--cpu", str(cpu), "-E", SEARCH_EVALUE, "--noali",
              "--domtblout", str(tmp), "-o", str(out_stem.with_suffix(".log")),
              str(hmm), str(db)])
    tmp.rename(domtbl)          # atomic: an interrupted search leaves no
    return domtbl, dt           # file a later run would accept as complete


def filtered_hits(domtbl: Path, evalue: float) -> dict[str, dict]:
    """Best hit per target at a sequence E-value cut, from one domtblout."""
    rows = [r for r in parse_domtblout(domtbl) if r["full_evalue"] <= evalue]
    return {t: h for t, h in best_hits_by_target(rows).items()}


def sweep_group(group: str, cpu: int, raw_dir: Path) -> dict:
    """hmmsearch both profiles over one group, then the margin assignment."""
    paths = group_paths(group)
    if not paths["db"].exists():
        raise SystemExit(f"{group}: sweep DB missing ({paths['db']}) — "
                         "run scripts/s20_fetch.py first")
    stats: dict = {"group": group, "db": str(paths["db"]),
                   "search_evalue": SEARCH_EVALUE,
                   "primary_evalue": EVALUE_PRIMARY, "profiles": {}}
    domtbls: dict[str, Path] = {}
    for name in PROFILE_NAMES:
        hmm = S20_DIR.parent / "hmm_sweep" / f"{name}.hmm"
        stem = raw_dir / f"hmmsearch_{name}_vs_{group}"
        domtbl, dt = search(hmm, paths["db"], stem, cpu)
        domtbls[name] = domtbl
        n_relaxed = len(filtered_hits(domtbl, float(EVALUE_RELAXED)))
        n_primary = len(filtered_hits(domtbl, float(EVALUE_PRIMARY)))
        stats["profiles"][name] = {
            "targets_relaxed": n_relaxed, "targets_primary": n_primary,
            "elapsed_s": round(dt, 1), "domtblout": str(domtbl)}
        log("s20_sweep", f"{group}/{name}.hmm: {n_primary} targets at "
                         f"E<={EVALUE_PRIMARY} ({n_relaxed} at E<={EVALUE_RELAXED})")

    ecut = float(EVALUE_PRIMARY)
    rows = assign(_primary(domtbls["itpr"], ecut), _primary(domtbls["ryr"], ecut))
    for r in rows:
        r["group"] = group
    write_tsv(S20_DIR / f"assignments_{group}.tsv",
              ASSIGN_FIELDS + ["group"], rows)
    calls = {c: sum(1 for r in rows if r["assignment"] == c)
             for c in ("ITPR", "RYR", "unassigned")}
    stats["assignment"] = {"targets": len(rows), "calls": calls}
    log("s20_sweep", f"{group}: {len(rows)} targets → {calls}")
    return stats


def _primary(domtbl: Path, evalue: float) -> dict[str, dict]:
    """`s3_assign.load_profile_hits`, at an explicit E-value cut.

    Line for line the same collapse — accession key, best full-sequence
    score wins a duplicate — because the assignment `assign()` makes has to
    be the assignment S3 made. The only difference is that the domain rows
    are filtered to the primary threshold first, so one search can serve
    both sensitivities.
    """
    out: dict[str, dict] = {}
    for target, hit in filtered_hits(domtbl, evalue).items():
        meta = parse_uniprot_header(f"{target} {hit['description']}")
        acc = acc_key(meta["accession"])
        prev = out.get(acc)
        if prev is None or hit["full_score"] > prev["full_score"]:
            out[acc] = {**hit, **meta}
    return out


# ------------------------------------------------------------ presence table
def presence(groups: list[str]) -> list[dict]:
    """Per-proteome ITPR presence/absence — the range result's raw form.

    Hits are attributed to the proteome file they actually came from, not to
    their taxid: several dozen taxids across these groups carry two reference
    proteomes each, and a taxid key would credit both with either one's hits.
    """
    rows: list[dict] = []
    for group in groups:
        assigns = S20_DIR / f"assignments_{group}.tsv"
        if not assigns.exists():
            continue
        scored = [r for r in read_tsv(assigns) if r["accession"]]
        acc2upid = accession_to_upid(group, {r["accession"] for r in scored})
        by_upid: dict[str, list[dict]] = {}
        for r in scored:
            upid = acc2upid.get(r["accession"])
            if upid:
                by_upid.setdefault(upid, []).append(r)
        for m in load_group_manifest(group):
            if m["status"] != "swept":
                continue
            hits = by_upid.get(m["upid"], [])
            itpr = [h for h in hits if h["assignment"] == "ITPR"]
            best = max(itpr, key=lambda h: float(h["itpr_score"] or 0),
                       default=None)
            rows.append({
                "upid": m["upid"], "organism": m["organism"],
                "taxid": m["taxid"], "group": group,
                "protein_count": m["protein_count"],
                "n_itpr": len(itpr),
                "n_ryr": sum(1 for h in hits if h["assignment"] == "RYR"),
                "n_unassigned": sum(1 for h in hits
                                    if h["assignment"] == "unassigned"),
                "best_itpr_score": best["itpr_score"] if best else "",
                "best_itpr_evalue": best["itpr_evalue"] if best else "",
                "best_itpr_acc": best["accession"] if best else "",
                "itpr_status": "present" if itpr else "absent",
            })
    return rows


# ------------------------------------------------------------ relaxed panel
RELAXED_FIELDS = ["group", "profile", "accession", "target", "species",
                  "taxon_id", "gene", "length", "full_evalue", "full_score",
                  "hmm_coverage", "qlen", "positions", "protein_name"]


def relaxed(groups: list[str], cpu: int, raw_dir: Path) -> list[dict]:
    """Every profile in the relaxed panel over the negative-claim groups."""
    out: list[dict] = []
    for group in groups:
        db = group_paths(group)["db"]
        for label, hmm in relaxed_panel():
            # The two full-length profiles were already searched at the
            # relaxed threshold by the primary stage — `-E` is a reporting
            # threshold, so that output *is* the relaxed one. Re-running them
            # here would cost two full searches per group to reproduce a file
            # already on disk.
            stem = (raw_dir / f"hmmsearch_{label}_vs_{group}"
                    if label in PROFILE_NAMES
                    else raw_dir / f"relaxed_{label}_vs_{group}")
            domtbl, _ = search(hmm, db, stem, cpu)
            hits = filtered_hits(domtbl, float(EVALUE_RELAXED))
            log("s20_sweep", f"relaxed {group}/{label}: {len(hits)} targets "
                             f"at E<={EVALUE_RELAXED}")
            for target, h in hits.items():
                meta = parse_uniprot_header(
                    f"{target} {h['description']}")
                out.append({
                    "group": group, "profile": label,
                    "accession": meta["accession"], "target": target,
                    "species": meta["species"], "taxon_id": meta["taxon_id"],
                    "gene": meta["gene"], "length": h["tlen"],
                    "full_evalue": h["full_evalue"],
                    "full_score": h["full_score"],
                    "hmm_coverage": round(h["hmm_coverage"], 4),
                    "qlen": h["qlen"],
                    "positions": round(h["hmm_coverage"] * h["qlen"]),
                    "protein_name": meta["protein_name"],
                })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", action="append", choices=list(GROUPS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--stage", choices=["hmmsearch", "relaxed", "both"],
                    default="both")
    args = ap.parse_args()

    groups = list(ALL_GROUPS) if args.all else (args.group or [])
    if not groups:
        ap.error("give --group (repeatable) or --all")
    _, raw_dir = s20_dirs()

    stats = read_json(S20_DIR / "sweep_stats.json") or {"groups": {}}
    if args.stage in ("hmmsearch", "both"):
        for i, g in enumerate(groups):
            live([("fetch proteomes", True)]
                 + [(f"sweep {n}", n in groups[:i]) for n in groups]
                 + [("relaxed panel", False), ("verdicts", False),
                    ("jackhmmer", False), ("census v5", False)])
            stats["groups"][g] = sweep_group(g, args.cpu, raw_dir)
            write_json(S20_DIR / "sweep_stats.json", stats)
        # Always over every group with an assignment table, never just
        # the ones this invocation swept: the sweeps run one group at a
        # time, and a per-invocation presence table would overwrite the
        # file with a fraction of the denominator each run.
        rows = presence(list(ALL_GROUPS))
        write_tsv(S20_DIR / "proteome_presence.tsv", PRESENCE_FIELDS, rows)
        n_present = sum(1 for r in rows if r["itpr_status"] == "present")
        log("s20_sweep", f"presence: {n_present}/{len(rows)} swept proteomes "
                         f"carry an ITPR call")

    if args.stage in ("relaxed", "both"):
        targets = [g for g in groups if GROUPS[g]["relaxed"]]
        live([("fetch proteomes", True), ("sweep", True),
              ("relaxed panel", False), ("verdicts", False),
              ("jackhmmer", False), ("census v5", False)])
        rows = relaxed(targets, args.cpu, raw_dir)
        write_tsv(S20_DIR / "relaxed_hits.tsv", RELAXED_FIELDS, rows)
        log("s20_sweep", f"relaxed panel: {len(rows)} hits over "
                         f"{len(targets)} groups → relaxed_hits.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
