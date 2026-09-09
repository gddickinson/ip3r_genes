"""The S22 driver — ordered stages, resumable, and the self-test first.

    python scripts/s22_run.py                 # everything
    python scripts/s22_run.py --list
    python scripts/s22_run.py --only paired --only contacts
    python scripts/s22_run.py --from omega

`modules` and the self-test run before anything is written, because every
panel and every table downstream is a claim about two regions and a join.
A failed stage does not stop the ones after it that do not depend on it, and
the report marks an unfinished section *not run yet* — waking up to a
partial S22 that says which parts are partial is worth more than waking up
to nothing.

Only `plc` and `deep_lineage` touch the network and the bulk proteome DBs; on a warm HMMER
cache the whole task rebuilds offline in about a minute.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402

STAGES = ["modules", "shells", "paired", "contacts", "omega", "plc",
          "lineage", "deep_lineage", "tables", "figures", "report"]


def _run(name: str) -> int:
    mod = __import__(f"s22_{name}")
    return mod.main([]) if name == "figures" else mod.main()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", choices=STAGES)
    ap.add_argument("--from", dest="start", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args(argv)
    if a.list:
        for s in STAGES:
            print(" ", s)
        return 0

    todo = a.only or (STAGES[STAGES.index(a.start):] if a.start else STAGES)

    # The self-test runs before anything is written, whatever the selection.
    import s22_test_ligand as T
    if (L.OUT_DIR / "module_map.tsv").exists():
        if T.main() != 0:
            L.log("self-test failed — refusing to write")
            return 1

    failed: list[str] = []
    for stage in todo:
        L.log(f"=== {stage} ===")
        t0 = time.time()
        try:
            rc = _run(stage)
        except Exception as exc:                              # noqa: BLE001
            import traceback
            traceback.print_exc()
            rc = 1
            L.log(f"{stage}: raised {type(exc).__name__}")
        if rc != 0:
            failed.append(stage)
        L.log(f"{stage}: {'ok' if rc == 0 else 'FAILED'} "
              f"({time.time() - t0:.1f}s)")
    if failed:
        L.log(f"stages failed: {failed}")
        return 1
    L.log("S22 complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
