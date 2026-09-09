"""The S18 driver — ordered stages, each resumable, each writing tables only.

`loci → protein → zero → corrections → tables → figures → report`, with
`--only` / `--from` / `--list`. The self-test runs **before anything is
written** and the run refuses to continue if it fails. A failed stage does not
stop the ones after it that do not depend on it, and the report marks an
unfinished section *not run yet*, because waking up to a partial S18 that says
which parts are partial is worth more than waking up to nothing.

Everything but `protein`'s first run is offline: the locus half reads archived
GFFs and archived miniprot output, and `protein` caches its blastp table under
the data root, so a re-run is deterministic (D24).
"""

from __future__ import annotations

import argparse
import pickle
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_lib as L                                               # noqa: E402
import s18_locus_audit as A                                       # noqa: E402
import s18_locus_rules as R                                       # noqa: E402
import s18_corrections as C                                       # noqa: E402
import s18_protein_audit as P                                     # noqa: E402
import s18_tables as T                                            # noqa: E402
import s18_zero_hits as Z                                         # noqa: E402

STAGES = ["loci", "protein", "zero", "corrections", "tables", "figures",
          "report"]
#: The bar sweep the sensitivity grid is computed over, and the piece floors.
BARS = [0.30, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
PIECES = [0.02, 0.05, 0.10]
CACHE = L.OUT / ".measurements.pkl"       # gitignored; the measurement pass


def _rules_fingerprint() -> str:
    """What the cache is only valid for: the parsers, not the rules.

    The cache holds what the two GFF readers found, and nothing a rule
    decides — states and verdicts are recomputed on every run. But the readers
    themselves can change, so their source is hashed into the cache and a
    mismatch re-reads the annotations rather than silently serving a parse
    made by different code.
    """
    import hashlib
    h = hashlib.sha256()
    for name in ("s10_gff.py", "s18_locus_audit.py"):
        h.update(Path(__file__).with_name(name).read_bytes())
    return h.hexdigest()[:16]


def _cached_loci(refresh: bool = False) -> list[dict]:
    """The measurement pass, cached: 274 gzipped GFFs is four minutes.

    The cache holds only what the parsers found. Every state and every verdict
    is recomputed from it, so a rule change cannot survive in a committed
    table (it did once — see `s18_locus_audit.apply_verdicts`).
    """
    if CACHE.exists() and not refresh:
        blob = pickle.loads(CACHE.read_bytes())
        if isinstance(blob, dict) and blob.get("fp") == _rules_fingerprint():
            return blob["rows"]
        print("  [loci] cache is from different reader code — re-reading")
    genomes = L.genome_rows()
    summaries = L.load_summaries()
    t0 = time.time()
    rows = A.measure(summaries, genomes, on_progress=lambda i, n: (
        print(f"  [loci] {i}/{n} genomes  {time.time() - t0:.0f}s", flush=True),
        L.live("loci", [("loci", False), ("protein", False), ("zero", False),
                        ("corrections", False), ("tables", False),
                        ("figures", False), ("report", False)])))
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_bytes(pickle.dumps({"fp": _rules_fingerprint(), "rows": rows}))
    return rows


def stage_loci(args) -> dict:
    rows = _cached_loci(args.refresh)
    A.apply_verdicts(rows)
    A.apply_bar(rows, R.COMPLETE_FRAC)
    A.contiguity_join(rows, L.contiguity_index(), L.integrity_index())
    cal = R.calibrate_complete(rows, L.quantile)
    A.apply_bar(rows, R.COMPLETE_FRAC)          # the bar the tables use
    L.write_tsv(L.OUT / "locus_audit.tsv", rows, A.AUDIT_COLS)
    scorable = [r for r in rows if r["gene_set"] == "present"]
    L.write_tsv(L.OUT / "locus_state_counts.tsv",
                T.state_counts(rows, ["cell"])
                + T.state_counts(scorable, ["cell", "source"]),
                None)
    L.write_tsv(L.OUT / "state_by_cell.tsv",
                T.state_counts(scorable, ["cell"]), None)
    L.write_tsv(L.OUT / "state_by_class.tsv",
                T.state_counts([r for r in scorable if not r["is_control"]],
                               ["vclass"]), None)
    L.write_tsv(L.OUT / "by_source.tsv", T.by_source(rows), None)
    L.write_tsv(L.OUT / "family_vs_control.tsv", T.family_vs_control(rows),
                None)
    L.write_tsv(L.OUT / "complete_calibration.tsv", [cal], None)
    L.write_tsv(L.OUT / "state_sensitivity.tsv",
                A.sensitivity(rows, BARS, PIECES), None)
    L.write_tsv(L.OUT / "locus_name_verdicts.tsv",
                T.verdict_counts(scorable, "symbol_verdict", ["cell", "source"])
                + T.verdict_counts(scorable, "product_verdict",
                                   ["cell", "source"]), None)
    return {"rows": rows, "calibration": cal}


def stage_protein(args) -> dict:
    census = L.read_tsv(L.CENSUS_V6)
    scope, excluded = [], []
    for r in census:
        ok, why = P.in_scope(r)
        (scope if ok else excluded).append(
            r if ok else {"accession": r["accession"], "call": r.get("call", ""),
                          "source": r.get("source", ""),
                          "length": r.get("length", ""),
                          "excluded_by": why})
    work = L.data_root() / "s18"
    hits_path, _seqs, missing = P.resolve_and_blast(scope, work,
                                                    threads=args.threads)
    hits = P.parse_hits(hits_path)
    reach = L.panel_paralogs()
    labels = L.bait_labels()
    rows = [P.audit_record(r, hits.get(r["accession"], []), labels,
                           callable_by_panel=reach) for r in scope]
    for r in rows:
        r["sequence_resolved"] = int(r["accession"] not in set(missing))
    L.write_tsv(L.OUT / "protein_audit.tsv", rows,
                P.AUDIT_COLS + ["sequence_resolved"])
    L.write_tsv(L.OUT / "protein_scope_excluded.tsv", excluded,
                ["accession", "call", "source", "length", "excluded_by"])
    L.write_tsv(L.OUT / "protein_name_verdicts.tsv",
                T.verdict_counts(rows, "symbol_verdict", ["census_call", "group"])
                + T.verdict_counts(rows, "name_verdict",
                                   ["census_call", "group"]), None)
    L.write_tsv(L.OUT / "pfam_recall.tsv", P.pfam_recall(scope), None)
    return {"rows": rows, "n_missing": len(missing)}


def stage_zero(args) -> dict:
    rows = Z.resolve(L.read_tsv(L.ZERO_HIT), L.genome_rows(),
                     L.load_summaries())
    L.write_tsv(L.OUT / "zero_hit_proteomes.tsv", rows, Z.COLS)
    L.write_tsv(L.OUT / "zero_hit_summary.tsv", Z.summarise(rows), None)
    return {"rows": rows}


def stage_corrections(args, loci: list[dict], proteins: list[dict]) -> dict:
    rows = C.from_loci(loci) + C.from_proteins(proteins)
    L.write_tsv(L.OUT / "corrections.tsv", rows, C.COLS)
    L.write_tsv(L.OUT / "correction_summary.tsv", C.summarise(rows), None)
    return {"rows": rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--from", dest="from_stage")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--refresh", action="store_true",
                    help="re-read every annotation instead of the cache")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    if args.list:
        print("\n".join(STAGES))
        return 0

    todo = STAGES
    if args.from_stage:
        todo = STAGES[STAGES.index(args.from_stage):]
    if args.only:
        todo = [s for s in STAGES if s in args.only]

    if subprocess.run([sys.executable,
                       str(Path(__file__).with_name("s18_test_audit.py"))]
                      ).returncode:
        print("[s18] self-test failed — nothing written")
        return 1

    L.OUT.mkdir(parents=True, exist_ok=True)
    state: dict = {}
    failed = []
    done = {s: False for s in STAGES}

    def tick(stage):
        done[stage] = True
        L.live(stage, [(s, done[s]) for s in STAGES])

    for stage in todo:
        print(f"[s18] {stage}")
        try:
            if stage == "loci":
                state["loci"] = stage_loci(args)
            elif stage == "protein":
                state["protein"] = stage_protein(args)
            elif stage == "zero":
                state["zero"] = stage_zero(args)
            elif stage == "corrections":
                loci = state.get("loci", {}).get("rows") or _joined_loci()
                prot = state.get("protein", {}).get("rows") or \
                    L.read_tsv(L.OUT / "protein_audit.tsv")
                state["corrections"] = stage_corrections(args, loci, prot)
            elif stage == "tables":
                _write_stats(state, args)
            elif stage == "figures":
                import s18_figures
                s18_figures.main()
            elif stage == "report":
                import s18_report
                s18_report.main()
            tick(stage)
        except Exception as exc:                     # noqa: BLE001
            failed.append(f"{stage}: {exc}")
            print(f"[s18] stage {stage} FAILED: {exc}")
            import traceback
            traceback.print_exc()
    if failed:
        print("[s18] failed stages: " + "; ".join(failed))
        return 1
    print("[s18] done")
    return 0


def _joined_loci() -> list[dict]:
    """The committed locus table, re-typed, for a stage run on its own."""
    rows = L.read_tsv(L.OUT / "locus_audit.tsv")
    for r in rows:
        for k in ("coverage", "best_model_frac", "union_coding_frac",
                  "namer_frac_cds"):
            r[k] = float(r.get(k) or 0)
        for k in ("locus_idx", "is_control", "n_coding_models",
                  "n_noncoding_models", "locus_cds_bp"):
            r[k] = int(float(r.get(k) or 0))
    return rows


def _write_stats(state: dict, args) -> None:
    loci = state.get("loci", {}).get("rows") or _joined_loci()
    cal = state.get("loci", {}).get("calibration") or \
        (L.read_tsv(L.OUT / "complete_calibration.tsv") or [{}])[0]
    prot = state.get("protein", {}).get("rows") or \
        L.read_tsv(L.OUT / "protein_audit.tsv")
    for r in prot:
        r["paralog_askable"] = int(float(r.get("paralog_askable") or 0))
    corr = state.get("corrections", {}).get("rows") or \
        L.read_tsv(L.OUT / "corrections.tsv")
    for r in corr:
        r["vetoed"] = int(float(r.get("vetoed") or 0))
    zero = state.get("zero", {}).get("rows") or \
        L.read_tsv(L.OUT / "zero_hit_proteomes.tsv")
    for k in ("bar", "median", "frac_below_bar", "n"):
        if k in cal:
            cal[k] = float(cal[k])
    head = T.headline(loci, prot, corr, zero, cal)
    paths = {p.stem: p for p in sorted(L.OUT.glob("*.tsv"))}
    T.write_stats(paths, {
        "complete_frac": R.COMPLETE_FRAC,
        "complete_frac_source": "s5_classify.ANNOT_CDS_FRAC",
        "overcredit_frac": R.OVERCREDIT_FRAC,
        "min_piece_frac": R.MIN_PIECE_FRAC,
        "min_touch_frac": R.MIN_TOUCH_FRAC,
        "calib_min_coverage": R.CALIB_MIN_COVERAGE,
        "protein_min_length_aa": P.MIN_FULL_LENGTH,
        "protein_rel_margin": P.REL_MARGIN,
        "min_rename_bits": C.MIN_RENAME_BITS,
        "correction_min_recovery": C.MIN_RECOVERY,
        "correction_strong_recovery": C.STRONG_RECOVERY,
        "genome_found_coverage": Z.GENOME_FOUND_COVERAGE,
        "window_pad_bp": A.WINDOW_PAD,
        "bars": BARS, "pieces": PIECES,
    }, "all negative controls pass (s18_test_audit.py)", head)


if __name__ == "__main__":
    raise SystemExit(main())
