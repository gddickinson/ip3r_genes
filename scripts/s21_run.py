"""The S21 driver — ordered stages, each resumable, each writing tables only.

`measure → calibrate → frame → architecture → introns → annot → tandem →
tables → figures → report`, with `--only` / `--from` / `--list`. The negative
controls run **before anything is written** and the run refuses to continue if
they fail. A failed stage does not stop the ones after it that do not depend on
it, and the report marks an unfinished section *not run yet*, because waking up
to a partial S21 that says which parts are partial is worth more than waking up
to nothing.

Two measurement passes are expensive and cached under the data root: the block
pass (70 s over 309 retained GFFs and their genomes) and the annotation pass
(4 min over 262 gzipped GFF3s). Both caches hold only what the readers found —
every rule, every state and every verdict is recomputed on each run (S18's D55:
a rule change once survived into a committed table because the verdicts were
cached beside the measurement). Everything else is committed or cached, so a
re-run is offline and deterministic (D24).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_annot as AN                                            # noqa: E402
import s21_architecture as A                                      # noqa: E402
import s21_blocks as B                                            # noqa: E402
import s21_calibrate_intron as CAL                                # noqa: E402
import s21_frame as FR                                            # noqa: E402
import s21_introns as I                                           # noqa: E402
import s21_lib as L                                               # noqa: E402
import s21_measure as M                                           # noqa: E402
import s21_tables as T                                            # noqa: E402
import s21_tandem as TD                                           # noqa: E402

STAGES = ["measure", "calibrate", "frame", "architecture", "introns", "annot",
          "tandem", "tables", "figures", "report"]


def _annot_cache() -> tuple[Path, Path, Path]:
    d = M.cache_dir()
    return (d / "annot_concordance.tsv", d / "annot_termini.tsv",
            d / "annot_fragments.tsv")


def stage_measure(args) -> None:
    if not args.refresh and (M.cache_dir() / "block_pairs.tsv").exists():
        L.log("block pass already cached — use --refresh to re-read")
        return
    t0 = time.time()
    M.run(lambda i, n, s: L.log(f"[measure] {i}/{n} genomes {s:.0f}s"))
    L.log(f"block pass in {time.time() - t0:.0f}s")


def stage_calibrate(args) -> dict:
    return CAL.write(M.load_pairs())


def stage_frame(args) -> dict:
    anchors = FR.require_anchors()
    frames = FR.build(refresh=args.refresh)
    L.write_tsv(L.OUT / "frame_map.tsv", FR.frame_rows(frames), None)
    L.write_tsv(L.OUT / "frame_anchors.tsv", anchors, None)
    return frames


def stage_architecture(args, frames: dict) -> tuple[list[dict], list[dict]]:
    models = M.load_models()
    pairs = M.load_pairs()
    mi = B.min_intron_bp()
    arows, irows = A.build(models, pairs, frames, mi)
    L.write_tsv(L.OUT / "architecture_loci.tsv", arows, list(A.ARCH_COLS))
    L.write_tsv(L.OUT / "architecture_by_paralog.tsv",
                A.summarise(arows, ("cell",)), None)
    L.write_tsv(L.OUT / "architecture_by_class.tsv",
                A.summarise(arows, ("cell", "vclass")), None)
    L.write_tsv(L.OUT / "architecture_sensitivity.tsv",
                A.sensitivity(models, pairs, frames, mi), None)
    paired = A.paired(arows)
    for vclass in sorted({r["vclass"] for r in arows if r["in_scope"]}):
        rows = [r for r in arows if r["vclass"] == vclass]
        if len({r["accession"] for r in rows if r["in_scope"]}) >= 10:
            paired += A.paired(arows, vclass)
    L.write_tsv(L.OUT / "paired_comparisons.tsv", paired, None)
    L.write_tsv(L.OUT / "junction_quality.tsv",
                T.junction_quality(arows, pairs), None)
    L.write_tsv(L.OUT / "intron_positions.tsv", irows, list(A.INTRON_COLS))
    return arows, irows


def stage_introns(args, arows: list[dict], irows: list[dict],
                  frames: dict) -> None:
    pos = I.positions(irows)
    cons = I.conservation(pos, arows)
    L.write_tsv(L.OUT / "intron_conservation.tsv", cons, None)
    L.write_tsv(L.OUT / "intron_conservation_summary.tsv",
                I.conservation_summary(cons), None)
    shared = I.shared_table(pos, arows, frames)
    L.write_tsv(L.OUT / "shared_introns.tsv", shared, None)
    L.write_tsv(L.OUT / "shared_intron_summary.tsv",
                I.shared_summary(shared) + I.shared_summary(shared, True), None)


def stage_annot(args, arows: list[dict]) -> None:
    cp, tp, fp = _annot_cache()
    audit = L.audit_index()
    mi = B.min_intron_bp()
    if args.refresh or not tp.exists():
        t0 = time.time()
        crows, erows, _ = AN.run(
            arows, audit, mi,
            lambda i, n: L.log(f"[annot] {i}/{n} genomes "
                               f"{time.time() - t0:.0f}s") if i % 50 == 0
            else None)
        L.write_tsv(cp, crows, list(AN.CONC_COLS))
        L.write_tsv(tp, erows, list(AN.EDGE_COLS))
    else:
        L.log("annotation pass already cached — use --refresh to re-read")
        crows = L.read_tsv(cp)
        for r in crows:
            for k in ("n_annot_edges", "n_exact", "n_within_tol", "n_no_match",
                      "n_model_exons", "n_annot_blocks"):
                r[k] = int(float(r.get(k) or 0))
            for k in ("frac_exact", "frac_within_tol", "median_offset"):
                r[k] = float(r.get(k) or 0)
        erows = L.read_tsv(tp)
        for r in erows:
            for k in ("locus_idx", "model_index", "n_models", "position",
                      "distance_to_exon_edge", "splice_ok", "is_gene_terminus",
                      "both_axes_agree"):
                r[k] = int(float(r.get(k) or 0))
    for r in arows:
        key = (r["accession"], r["cell"], r["locus_idx"])
        r["state"] = (audit.get(key) or {}).get("state", "")
    # Verdicts are always recomputed from the cached termini (D55).
    lrows = AN.verdicts_from_termini(erows, arows)
    L.write_tsv(fp, lrows, list(AN.LOCUS_COLS))
    L.write_tsv(L.OUT / "boundary_concordance.tsv", crows, list(AN.CONC_COLS))
    L.write_tsv(L.OUT / "boundary_concordance_summary.tsv",
                AN.concordance_summary(crows), None)
    L.write_tsv(L.OUT / "fragment_termini.tsv", erows, list(AN.EDGE_COLS))
    L.write_tsv(L.OUT / "fragment_verdicts.tsv", lrows, list(AN.LOCUS_COLS))
    L.write_tsv(L.OUT / "fragment_summary.tsv", AN.fragment_summary(lrows), None)


def stage_tandem(args) -> None:
    prows, crows = TD.run(lambda i, n: L.log(f"[tandem] {i}/{n} genomes")
                          if i % 100 == 0 else None)
    L.write_tsv(L.OUT / "tandem_pairs.tsv", prows, list(TD.PAIR_COLS))
    L.write_tsv(L.OUT / "tandem_by_class.tsv", TD.class_summary(prows), None)
    L.write_tsv(L.OUT / "tandem_control.tsv", TD.control_summary(crows), None)
    L.write_tsv(L.OUT / "tandem_teleost_check.tsv", TD.teleost_check(crows),
                None)
    L.write_tsv(L.OUT / "tandem_cells.tsv", crows, list(TD.CELL_COLS))


def _loaded_arch() -> tuple[list[dict], list[dict]]:
    """The committed tables re-typed, for a stage run on its own."""
    arows = L.read_tsv(L.OUT / "architecture_loci.tsv")
    for r in arows:
        for k in ("locus_idx", "start", "end", "q_start", "q_end", "residues",
                  "n_blocks", "n_exons", "n_merges", "n_introns",
                  "frame_steps", "cds_bp", "span_bp", "total_intron_bp",
                  "max_intron_bp", "min_intron_bp_seen",
                  "canonical_junctions", "minor_junctions",
                  "non_canonical_junctions", "in_scope", "frame_via_cell"):
            r[k] = int(float(r.get(k) or 0))
        for k in ("coverage", "identity", "median_intron_bp", "mean_exon_bp",
                  "frac_canonical"):
            r[k] = float(r.get(k) or 0)
    irows = L.read_tsv(L.OUT / "intron_positions.tsv")
    for r in irows:
        for k in ("locus_idx", "intron_index", "length", "frame_ok",
                  "col_left", "col_right", "intron_phase", "q_left", "q_right"):
            r[k] = int(float(r.get(k) or 0))
    return arows, irows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--from", dest="from_stage")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--refresh", action="store_true",
                    help="re-read every GFF instead of the cached measurements")
    args = ap.parse_args()
    if args.list:
        print("\n".join(STAGES))
        return 0

    todo = STAGES
    if args.from_stage:
        todo = STAGES[STAGES.index(args.from_stage):]
    if args.only:
        todo = [s for s in STAGES if s in args.only]

    test_ok = subprocess.run(
        [sys.executable,
         str(Path(__file__).with_name("s21_test_arch.py"))]).returncode == 0
    if not test_ok:
        print("[s21] negative controls failed — nothing written")
        return 1

    L.OUT.mkdir(parents=True, exist_ok=True)
    L.FIGS.mkdir(parents=True, exist_ok=True)
    state: dict = {}
    failed: list[str] = []
    done = {s: False for s in STAGES}

    def tick(stage: str) -> None:
        done[stage] = True
        L.live(stage, [(s, done[s]) for s in STAGES])

    for stage in todo:
        print(f"[s21] {stage}")
        try:
            if stage == "measure":
                stage_measure(args)
            elif stage == "calibrate":
                state["calibration"] = stage_calibrate(args)
            elif stage == "frame":
                state["frames"] = stage_frame(args)
            elif stage == "architecture":
                frames = state.get("frames") or FR.build()
                state["frames"] = frames
                state["arch"], state["introns"] = stage_architecture(args,
                                                                     frames)
            elif stage == "introns":
                frames = state.get("frames") or FR.build()
                state["frames"] = frames
                arows, irows = (state.get("arch"), state.get("introns"))
                if arows is None:
                    arows, irows = _loaded_arch()
                stage_introns(args, arows, irows, frames)
            elif stage == "annot":
                arows = state.get("arch") or _loaded_arch()[0]
                stage_annot(args, arows)
            elif stage == "tandem":
                stage_tandem(args)
            elif stage == "tables":
                T.write_stats(test_ok, {"headline": T.headline()})
            elif stage == "figures":
                import s21_figures
                s21_figures.main()
            elif stage == "report":
                import s21_report
                s21_report.main()
            tick(stage)
        except Exception as exc:                                  # noqa: BLE001
            failed.append(f"{stage}: {exc}")
            print(f"[s21] stage {stage} FAILED: {exc}")
            import traceback
            traceback.print_exc()
    if failed:
        print("[s21] failed stages: " + "; ".join(failed))
        return 1
    print("[s21] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
