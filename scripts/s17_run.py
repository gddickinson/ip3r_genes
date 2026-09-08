"""S17 driver — constraint & function, in the one order that works.

The stages are not independent. `domains` puts the architecture into each
reference's numbering and is what every later stage's `element` column comes
from; `orthologs` builds the deep per-paralog alignments the scores are read
off; `conservation` needs both; `variants`, `fel` and `paint` all read the
per-site tables; `report` and `figures` render only from committed TSVs (D13).

`s17_test_constraint.py` runs **before anything is written** and the driver
refuses to continue if it fails. A failed stage does not stop the ones after it
that do not depend on it, and the report marks an unfinished section *not run
yet* — waking up to a partial S17 that says which parts are partial is worth
more than waking up to nothing.

    python scripts/s17_run.py                 # everything
    python scripts/s17_run.py --from variants
    python scripts/s17_run.py --only fel
    python scripts/s17_run.py --list

Measured stage costs on this machine: `domains` ~15 s (three pairwise MAFFT
transfers), `orthologs` ~4.5 min cold (two MAFFT passes per paralog — the
shape screen needs the first to measure on) and seconds warm, `conservation`
~25 s, `variants` ~2 min cold / ~10 s warm (ClinVar + UniProt, cached under
`<data_root>/raw_api/s17/`), **`fel` ~8 min** (three HyPhy runs; cached — pass
`--force` or delete `results/constraint/fel/*.json` to redo), `paint` ~90 s,
`tables`, `figures` and `report` seconds.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402


def _stages():
    import s17_conservation, s17_domains, s17_fel, s17_figures
    import s17_orthologs, s17_paint, s17_report, s17_tables, s17_variants

    return [
        ("domains", "the architecture in each reference's own numbering",
         lambda out: s17_domains.run(out)),
        ("orthologs", "deep per-paralog orthologue sets from the S5 models",
         lambda out: s17_orthologs.run(out)),
        ("conservation", "per-site constraint, four layers, + the element tests",
         lambda out: s17_conservation.run(out)),
        ("variants", "ClinVar/UniProt harvest, the AUC and the paralog audit",
         lambda out: s17_variants.run(out)),
        ("fel", "per-site dN/dS (HyPhy FEL) on the same coordinates",
         lambda out: s17_fel.run(out)),
        ("paint", "constraint into the B-factor column of the S11 structures",
         lambda out: s17_paint.run(out)),
        ("tables", "constraint_stats.json — rules, self-test, SHA-256 per table",
         lambda out: s17_tables.run(out)),
        ("figures", "the four S17 figures",
         lambda out: s17_figures.run(out)),
        ("report", "render report.md from the committed tables",
         lambda out: s17_report.build(out)),
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    ap.add_argument("--from", dest="start", default=None, help="resume here")
    ap.add_argument("--only", default=None, help="run just this stage")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--skip-self-test", action="store_true",
                    help="for debugging a stage only; never for a real build")
    args = ap.parse_args()

    stages = _stages()
    names = [n for n, _d, _f in stages]
    if args.list:
        for name, desc, _ in stages:
            print(f"  {name:<13} {desc}")
        return

    if not args.skip_self_test:
        import s17_test_constraint as T
        print(f"[s17] self-test: {len(T.TESTS)} constructed negative controls")
        if not T.self_test():
            raise SystemExit("[s17] self-test failed — refusing to write")

    if args.only:
        if args.only not in names:
            ap.error(f"unknown stage {args.only!r}; choose from {names}")
        stages = [s for s in stages if s[0] == args.only]
    elif args.start:
        if args.start not in names:
            ap.error(f"unknown stage {args.start!r}; choose from {names}")
        stages = stages[names.index(args.start):]

    args.out.mkdir(parents=True, exist_ok=True)
    t_all = time.time()
    failed = []
    for name, desc, fn in stages:
        t0 = time.time()
        print(f"\n=== {name} — {desc}")
        try:
            fn(args.out)
        except Exception as e:                                # noqa: BLE001
            failed.append(name)
            print(f"!!! {name} FAILED: {type(e).__name__}: {e}")
            continue
        print(f"--- {name} done in {time.time() - t0:.0f} s")
    print(f"\n[s17] {len(stages)} stage(s) in {time.time() - t_all:.0f} s "
          f"-> {args.out}")
    if failed:
        raise SystemExit(f"[s17] failed stages: {', '.join(failed)}")


if __name__ == "__main__":
    main()
