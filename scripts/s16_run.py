"""The S16 driver — ordered stages, resumable, self-tested before it writes.

    python3 scripts/s16_run.py                # every stage
    python3 scripts/s16_run.py --only copies  # one stage
    python3 scripts/s16_run.py --from map     # this stage onwards
    python3 scripts/s16_run.py --list

`s16_test_dup.py` runs **before anything is written** and the driver refuses
to continue if it fails. A failed stage does not stop the ones after it that
do not depend on it, and the report marks an unfinished section *not run yet*
— waking up to a partial S16 that says which parts are partial is worth more
than waking up to nothing (S9's rule).

Everything except `map` is offline. `map` makes 25 BioMart requests on a cold
cache (~3 min) and none on a warm one, so a rerun is reproducible without
network.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402


def _copies(log):
    import s16_copy_number as M
    return M.run(log=log)


def _map(log):
    import s16_ensembl_map as M
    import s16_paralogon as P
    return M.run(extra_symbols=P.human_window_symbols(), log=log)


def _paralogon(log):
    import s16_paralogon as M
    return M.run(log=log)


def _quartet(log):
    import s16_quartet as M
    return M.run(log=log)


def _teleost(log):
    import s16_teleost as M
    return M.run(log=log)


def _tables(log):
    import s16_tables as M
    p = M.stats(self_test=_SELF_TEST)
    log(f"[tables] wrote {p}")
    return p


def _figures(log):
    import s16_figures as M
    return M.run(log=log)


def _report(log):
    import s16_report as M
    p = M.render()
    log(f"[report] wrote {p}")
    return p


STAGES = [
    ("copies", "per-locus copy table, merges, copy-number landscape", _copies),
    ("map", "human paralogy map from Ensembl BioMart (pinned archive)", _map),
    ("paralogon", "2R cross-window paralogy + null + block scan", _paralogon),
    ("quartet", "quartet, pooled test, 309-genome replication + its null",
     _quartet),
    ("teleost", "3R: groups, dispersion, DCS, cross-anchor blocks", _teleost),
    ("tables", "stats json with parameters and SHA-256 per table", _tables),
    ("figures", "the four figures", _figures),
    ("report", "render report.md from the committed tables", _report),
]

_SELF_TEST = "not run"


def main() -> int:
    global _SELF_TEST
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--from", dest="start", default="")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--skip-self-test", action="store_true",
                    help="for debugging a single stage only; the full run "
                         "must never use it")
    args = ap.parse_args()

    if args.list:
        for name, desc, _ in STAGES:
            print(f"  {name:10s} {desc}")
        return 0

    names = [s[0] for s in STAGES]
    todo = names
    if args.only:
        todo = [x for x in names if x in args.only]
    elif args.start:
        if args.start not in names:
            print(f"unknown stage {args.start!r}")
            return 2
        todo = names[names.index(args.start):]

    if not args.skip_self_test:
        import s16_test_dup
        print("=" * 62)
        rc = s16_test_dup.main()
        _SELF_TEST = "pass" if rc == 0 else "FAIL"
        if rc != 0:
            print("\nself-test failed — refusing to write anything")
            return 1
    else:
        _SELF_TEST = "skipped"

    failed: list[str] = []
    for i, (name, desc, fn) in enumerate(STAGES):
        if name not in todo:
            continue
        L.live(name, [(s[0], s[0] in todo[:todo.index(name)])
                      for s in STAGES])
        print("=" * 62)
        print(f"[{name}] {desc}")
        t0 = time.time()
        try:
            fn(print)
            print(f"[{name}] done in {time.time() - t0:.1f}s")
        except Exception:                                       # noqa: BLE001
            traceback.print_exc()
            failed.append(name)
            print(f"[{name}] FAILED after {time.time() - t0:.1f}s — "
                  f"continuing with the stages that do not depend on it")
    L.live("done", [(s[0], s[0] in todo and s[0] not in failed)
                    for s in STAGES])
    print("=" * 62)
    if failed:
        print(f"S16: {len(failed)} stage(s) failed: {', '.join(failed)}")
        return 1
    print("S16: all stages completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
