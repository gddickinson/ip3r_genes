"""s15b_run.py — the S15b driver.

Ordered stages `dollo -> sensitivity -> mk -> fossils -> tables -> figures
-> report`, with `--only` / `--from` / `--list`.  Runs
`s15b_test_counts.py` **before writing anything** and refuses to continue
if it fails.  A failed stage does not stop the ones after it that do not
depend on it, and the report marks an unfinished section *not run yet* —
the S9/S10/S11/S12/S15a rule: waking up to a partial S15b that says which
parts are partial is worth more than waking up to nothing.

Everything S15b reads is committed, so the whole task runs offline and a
rerun is deterministic.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback

import s15b_coding as coding
import s15b_lib as lib

STAGES = ["dollo", "sensitivity", "mk", "fossils", "tables", "figures",
          "report"]


def _load_inputs() -> dict:
    ctx = dict(
        character_matrix=lib.read_tsv(lib.MATRIX),
        loss_candidates=lib.read_tsv(lib.CANDIDATES),
        integrity_loci=lib.read_tsv(lib.INTEGRITY_LOCI),
        integrity_pairs=lib.read_tsv(lib.INTEGRITY_PAIRS),
        integrity_tests=lib.read_tsv(lib.INTEGRITY_TESTS),
        calibrations=lib.load_calibrations(),
        s15a_stats=json.loads(lib.S15A_STATS.read_text(encoding="utf-8")),
    )
    ctx["tree"] = lib.load_tree()
    return ctx


def stage_dollo(ctx: dict) -> None:
    import s15b_dollo as dollo
    import s15b_sensitivity as sens
    import s15b_tables as tables
    ctx["dollo_counts"] = sens.dollo_counts(ctx["character_matrix"],
                                            ctx["tree"])
    ctx["loss_edges"] = sens.loss_edges(ctx["character_matrix"], ctx["tree"])
    ctx["polytomy_profile"] = dollo.polytomy_profile(ctx["tree"])
    ctx["resolution_bound"] = tables.resolution_bound_rows(
        ctx["polytomy_profile"], ctx["dollo_counts"])
    base = coding.recode_all(ctx["character_matrix"], **coding.BASE_SETTING)
    mism = sum(1 for r, d in zip(ctx["character_matrix"], base)
               if r["state"] != d["state"])
    ctx["base_reproduction"] = dict(
        n_rows=len(base), n_mismatches=mism,
        setting=coding.BASE_SETTING,
        note="the recoder is S15a's own rule chain; a mismatch here would "
             "make every comparison in this task a comparison against "
             "something that never ran")


def stage_sensitivity(ctx: dict) -> None:
    import s15b_sensitivity as sens
    ctx["manufactured_losses"] = sens.manufactured_losses(
        ctx["character_matrix"])


def stage_mk(ctx: dict) -> None:
    import s15b_mk as mk
    import s15b_sensitivity as sens
    grid, fits = sens.matrix(ctx["character_matrix"], ctx["tree"],
                             ctx["calibrations"])
    ctx["sensitivity_matrix"] = grid
    ctx["mk_fits"] = fits
    ctx["mk_profile"] = sens.profile_rows(ctx["tree"], ctx["character_matrix"],
                                          ctx["calibrations"])
    ctx["mk_models"] = mk.MODELS


def stage_fossils(ctx: dict) -> None:
    import s15b_fossils as fos
    bar = ctx["s15a_stats"]["integrity_bar"]["bar"]
    rows, stats = fos.denominator(ctx["integrity_loci"],
                                  ctx["character_matrix"], bar)
    ctx["fossil_loci"] = rows
    ctx["fossil_stats"] = stats
    ctx["fossil_by_cell"] = fos.by_cell(rows)
    pair_rows, tests = fos.by_class(ctx["integrity_pairs"],
                                    ctx["character_matrix"], matched=True)
    ctx["lesion_by_class"] = tests
    ctx["lesion_class_controls"] = fos.class_controls(
        pair_rows, tests, ctx["character_matrix"])
    ctx["class_min_n"] = fos.MIN_CLASS_N


def stage_tables(ctx: dict) -> None:
    import s15b_tables as tables
    ctx["headline"] = tables.headline(ctx)
    written = tables.write_all(ctx)
    tables.write_stats(ctx, written, ctx["self_test"])
    ctx["written"] = written


def stage_figures(ctx: dict) -> None:
    import s15b_figures
    s15b_figures.main()


def stage_report(ctx: dict) -> None:
    import s15b_report
    s15b_report.main()


RUNNERS = {"dollo": stage_dollo, "sensitivity": stage_sensitivity,
           "mk": stage_mk, "fossils": stage_fossils, "tables": stage_tables,
           "figures": stage_figures, "report": stage_report}

DEPENDS = {"tables": ("dollo", "sensitivity", "mk", "fossils"),
           "figures": ("tables",), "report": ("tables",)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="S15b driver")
    ap.add_argument("--only", action="append", choices=STAGES)
    ap.add_argument("--from", dest="start", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--skip-self-test", action="store_true",
                    help="for debugging only; the driver refuses to write "
                         "tables when the controls have not passed")
    args = ap.parse_args(argv)
    if args.list:
        print("\n".join(STAGES))
        return 0

    import s15b_test_counts
    st = (dict(status="skipped", n_tests=0, n_passed=0, n_failed=0,
               failures=[], mutations_caught=[])
          if args.skip_self_test else
          s15b_test_counts.self_test(verbose=True))
    if st["status"] == "fail":
        print(f"\nself-test failed ({st['n_failed']} of {st['n_tests']}); "
              f"nothing written")
        return 1
    if not args.skip_self_test:
        st["mutations_caught"] = s15b_test_counts.mutate(verbose=True)

    todo = (args.only if args.only else
            STAGES[STAGES.index(args.start):] if args.start else STAGES)
    ctx = _load_inputs()
    ctx["self_test"] = st
    print(f"\ninputs: {len(ctx['character_matrix'])} matrix cells, "
          f"{len(lib.tips(ctx['tree']))} tree tips, "
          f"{len(ctx['integrity_loci'])} scored-or-not loci")

    done, failed = [], []
    for stage in STAGES:
        if stage not in todo:
            continue
        missing = [d for d in DEPENDS.get(stage, ()) if d in failed
                   or (d not in done and d in todo)]
        if missing:
            print(f"[{stage}] skipped — depends on {missing}")
            failed.append(stage)
            continue
        lib.live(stage, [(s, s in done) for s in STAGES])
        t0 = time.time()
        try:
            RUNNERS[stage](ctx)
            done.append(stage)
            print(f"[{stage}] ok ({time.time() - t0:.1f}s)")
        except Exception as exc:                        # noqa: BLE001
            failed.append(stage)
            print(f"[{stage}] FAILED: {exc}")
            traceback.print_exc(limit=3)
    lib.live("done" if not failed else "partial",
             [(s, s in done) for s in STAGES])
    if "headline" in ctx:
        print("\nheadline: " + json.dumps(ctx["headline"], default=str))
    print(f"\nstages ok: {done}")
    if failed:
        print(f"stages failed: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
