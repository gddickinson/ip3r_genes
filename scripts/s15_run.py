"""s15_run.py — the S15a driver.

Ordered stages `recon -> synteny -> integrity -> tree -> matrix -> tables
-> figures -> report`, with `--only` / `--from` / `--list`.  Runs
`s15_test_loss.py` **before writing anything** and refuses to continue if
it fails.  A failed stage does not stop the ones after it that do not
depend on it, and the report marks an unfinished section *not run yet* —
the S9/S10/S11/S12 rule: waking up to a partial S15 that says which parts
are partial is worth more than waking up to nothing.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback

import s15_lib as lib

STAGES = ["recon", "synteny", "integrity", "tree", "matrix", "tables",
          "figures", "report"]


def _load_inputs() -> dict:
    return dict(
        ledger=lib.read_tsv(lib.LEDGER),
        rescue=lib.read_tsv(lib.RESCUE),
        manifest=lib.read_tsv(lib.MANIFEST),
        summaries=lib.load_summaries(),
    )


def stage_recon(ctx: dict) -> None:
    import s15_reconstruct as recon
    import s15_calibrate_recon as cal
    ctx["reconstruction"] = recon.run(ctx["ledger"], ctx["rescue"],
                                      ctx["summaries"])
    table, stats = cal.calibrate(ctx["reconstruction"])
    ctx["recon_calibration"] = table
    ctx["recon_calibration_stats"] = stats
    ctx["bar"] = cal.bar(stats)


def stage_synteny(ctx: dict) -> None:
    import s15_synteny_reach as syn
    ctx["synteny_reach"] = syn.reach(ctx["rescue"])
    ctx["synteny_reach_summary"] = [
        syn.summarise(ctx["synteny_reach"], w)
        for w in ("informative10", "fixed10")]
    ctx["synteny_calls"] = syn.apply_caller(ctx["rescue"],
                                            ctx["synteny_reach"])
    ctx["caller_accuracy_by_keys"] = syn.accuracy_by_keys(
        lib.read_tsv(lib.RESULTS / "synteny" / "caller_calibration.tsv"))


def stage_integrity(ctx: dict) -> None:
    import s15_integrity as integ
    man = {r["accession"]: r for r in ctx["manifest"]}
    rows = integ.loci_table(ctx["summaries"], man)
    bar = integ.calibrate_bar(rows)
    ctx["integrity_bar"] = bar
    ctx["integrity_loci"] = integ.verdicts(rows, bar)
    ctx["integrity_covariates"] = integ.covariates(rows)
    pairs, tests = [], []
    for matched in (False, True):
        p, t = integ.paired_within_genome(rows, match_identity=matched)
        pairs.extend(p)
        tests.extend(t)
    ctx["integrity_pairs"] = pairs
    ctx["integrity_tests"] = tests


def stage_tree(ctx: dict) -> None:
    import s15_species_tree as tree
    root, tips, audit = tree.build(ctx["manifest"])
    ctx["tree_root"] = root
    ctx["tree_placement"] = audit
    ctx["tree_polytomies"] = tree.polytomy_report(root)
    ctx["tree_vs_s13"] = tree.compare_with_s13(root, tips, ctx["manifest"])
    ctx["newick"] = tree.newick(root)


def stage_matrix(ctx: dict) -> None:
    import s15_matrix as mx
    m = mx.build(ctx["ledger"], ctx["summaries"], ctx["reconstruction"],
                 ctx["synteny_calls"], ctx["synteny_reach"], ctx["bar"])
    ctx["character_matrix"] = m
    ctx["character_matrix_wide"] = mx.wide(m)
    ctx["state_counts"] = mx.state_counts(m)
    ctx["implied_copies"] = mx.implied_copies(m)
    ctx["loss_candidates"] = mx.loss_candidates(m)


def stage_tables(ctx: dict) -> None:
    import s15_tables as tables
    ctx["headline"] = tables.headline(ctx["character_matrix"],
                                      ctx["implied_copies"],
                                      ctx["synteny_reach_summary"])
    written = tables.write_all(ctx)
    tables.write_stats(ctx, written, ctx["self_test"])
    ctx["written"] = written


def stage_figures(ctx: dict) -> None:
    import s15_figures
    s15_figures.main()


def stage_report(ctx: dict) -> None:
    import s15_report
    s15_report.main()


RUNNERS = {"recon": stage_recon, "synteny": stage_synteny,
           "integrity": stage_integrity, "tree": stage_tree,
           "matrix": stage_matrix, "tables": stage_tables,
           "figures": stage_figures, "report": stage_report}

#: a stage cannot run unless these have produced their output
DEPENDS = {"matrix": ("recon", "synteny"), "tables": ("matrix",),
           "figures": ("tables",), "report": ("tables",)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="S15a driver")
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

    import s15_test_loss
    st = (dict(status="skipped", n_tests=0, n_passed=0, n_failed=0,
               failures=[], mutations_caught=[])
          if args.skip_self_test else s15_test_loss.self_test(verbose=True))
    if st["status"] == "fail":
        print(f"\nself-test failed ({st['n_failed']} of {st['n_tests']}); "
              f"nothing written")
        return 1

    todo = (args.only if args.only else
            STAGES[STAGES.index(args.start):] if args.start else STAGES)
    ctx = _load_inputs()
    ctx["self_test"] = st
    print(f"\ninputs: {len(ctx['ledger'])} ledger rows, "
          f"{len(ctx['rescue'])} rescue regions, "
          f"{len(ctx['summaries'])} genome summaries")

    done, failed = [], []
    for stage in STAGES:
        if stage not in todo:
            continue
        missing = [d for d in DEPENDS.get(stage, ()) if d in failed]
        if missing:
            print(f"[{stage}] skipped — depends on failed {missing}")
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
