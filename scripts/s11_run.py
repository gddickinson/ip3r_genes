"""The S11 driver — ordered stages, each resumable, each writing tables only.

    afdb -> panel -> tmalign -> plddt -> foldseek -> tables -> figures -> report

`--only` / `--from` / `--list` as in `s9_run.py` and `s10_run.py`, and for
the same reason: a stage that fails does not stop the ones after it that do
not depend on it, and the report renders *not run yet* for a section whose
table is absent, so a partial S11 says which parts are partial.

`foldseek` is off by default (`--with-foldseek` turns it on): the brief
marks it optional and it is the one stage that needs a multi-gigabyte index
downloaded before it can say anything.

The self-test (`s11_test_structures.py`) runs **before** anything is
written, and a failure stops the run.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_afdb_probe as afdb                                   # noqa: E402
import s11_foldseek as fold                                     # noqa: E402
import s11_panel as panel                                       # noqa: E402
import s11_plddt as plddt                                       # noqa: E402
import s11_tables as tb                                         # noqa: E402
import s11_test_structures as tests                             # noqa: E402
import s11_tmalign_run as tmr                                   # noqa: E402
from s11_lib import load_tsv                                    # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "structures"
STAGES = ["afdb", "panel", "tmalign", "plddt", "foldseek", "tables",
          "figures", "report"]


def _live(stage: str, completed: set[str], wanted: set[str]) -> None:
    """Feed the dashboard's live panel.

    A stage is `done` only if it actually ran and finished in *this*
    invocation. Marking every stage done at the end of a `--only panel`
    run — which the first version did — puts a full green bar on the
    dashboard for a run that did one eighth of the work.
    """
    try:
        (ROOT / "results" / "session_live.json").write_text(json.dumps({
            "task": "S11", "workers": 1,
            "steps": [{"label": s, "done": s in completed,
                       "skipped": s not in wanted} for s in STAGES],
            "current": stage,
            "progress": f"{len(completed)}/{len(wanted)}"}) + "\n")
    except OSError:
        pass


def _load(name: str) -> list[dict]:
    path = OUT / name
    return load_tsv(path) if path.exists() else []


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", choices=STAGES)
    ap.add_argument("--from", dest="start", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--with-foldseek", action="store_true")
    ap.add_argument("--foldseek-subset", default="afdb_swissprot")
    ap.add_argument("--foldseek-threads", type=int, default=2)
    args = ap.parse_args(argv)

    if args.list:
        for s in STAGES:
            print(s)
        return 0

    wanted = set(args.only) if args.only else set(STAGES)
    if args.start:
        wanted &= set(STAGES[STAGES.index(args.start):])
    if not args.with_foldseek and not args.only:
        wanted.discard("foldseek")

    problems = tests.run()
    if problems:
        for p in problems:
            print("FAIL", p)
        print("[s11] self-test failed — nothing written")
        return 2
    print(f"[s11] self-test: 0 failures")

    failures: list[str] = []
    written: list[Path] = []
    manifest: list[dict] = []
    tm_rows: list[dict] = []
    self_tests = {"s11_test_structures": "pass"}

    completed: set[str] = set()
    for stage in STAGES:
        if stage not in wanted:
            continue
        _live(stage, completed, wanted)
        t0 = time.time()
        try:
            if stage == "afdb":
                afdb.main(["--set", "both", "--workers",
                           str(max(4, args.workers))])

            elif stage == "panel":
                (manifest, cands, ctrl_cands, unfilled, slots,
                 selection) = panel.build_panel()
                written += [
                    tb.write_candidates(cands),
                    tb.write_control_candidates(ctrl_cands),
                    tb.write_unfilled(unfilled),
                    tb.write_slots(slots),
                    tb.write_manifest(manifest, panel.MANIFEST_COLUMNS),
                ]
                written.append(tb.write_selection(selection))
                print(f"[s11] panel: {len(manifest)} structures, "
                      f"{sum(1 for m in manifest if m.get('status') == 'ok')} usable")

            elif stage == "tmalign":
                manifest = manifest or _load("structure_manifest.tsv")
                tm_rows = tmr.run(
                    manifest, workers=args.workers,
                    on_progress=lambda d, n: print(f"[s11]   tmalign {d}/{n}",
                                                   flush=True))
                probs = tmr.self_test(tm_rows)
                self_tests["tmalign_self_test"] = (
                    "pass" if not probs else f"{len(probs)} failures")
                for p in probs[:10]:
                    print("[s11] tmalign self-test:", p)
                written.append(tb.write_tm_scores(tm_rows, tmr.PAIR_COLUMNS))
                written.append(tb.write_vs_reference(
                    tb.vs_reference(manifest, tm_rows)))
                print(f"[s11] tmalign: {len(tm_rows)} pairs")

            elif stage == "plddt":
                manifest = manifest or _load("structure_manifest.tsv")
                rows = plddt.run(manifest)
                written += tb.write_plddt(rows, plddt.COLUMNS,
                                          plddt.domain_summary(rows),
                                          plddt.SUMMARY_COLUMNS)
                print(f"[s11] plddt: {len(rows)} domain rows")

            elif stage == "foldseek":
                manifest = manifest or _load("structure_manifest.tsv")
                rows = fold.run(manifest, subset=args.foldseek_subset,
                                threads=args.foldseek_threads)
                written.append(fold.write(rows))
                print(f"[s11] foldseek: {len(rows)} hits")

            elif stage == "tables":
                # Every committed table in the directory, not only the ones
                # this invocation happened to write. A `--from tmalign` run
                # skips the panel stage, and the first version then recorded
                # a SHA-256 for 8 of the 14 tables — a manifest that is
                # silently partial is worse than none, because the missing
                # rows look like tables that do not exist.
                stats = tb.write_stats(
                    sorted(OUT.glob("*.tsv")),
                    {"panel_size": len(manifest or _load("structure_manifest.tsv")),
                     "full_coverage_bar": panel.FULL_COVERAGE,
                     "min_chain_residues": panel.MIN_CHAIN_RESIDUES,
                     "min_domain_cover": plddt.MIN_DOMAIN_COVER,
                     "foldseek_subset": (args.foldseek_subset
                                         if args.with_foldseek else "not run")},
                    self_tests)
                print(f"[s11] stats -> {stats}")

            elif stage == "figures":
                import s11_figures
                s11_figures.main([])

            elif stage == "report":
                import s11_report
                s11_report.main([])

            completed.add(stage)
            print(f"[s11] {stage} done in {time.time() - t0:.0f}s")
        except Exception as exc:                                # noqa: BLE001
            failures.append(stage)
            print(f"[s11] {stage} FAILED: {exc}")
            traceback.print_exc()

    _live("done", completed, wanted)
    if failures:
        print("[s11] failed stages:", ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
