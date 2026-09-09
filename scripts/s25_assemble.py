"""S25 driver — build the thesis.

Stages, in order:

  guards   break every build guard on purpose and check that it fires
  assign   commit the chapter grouping and guard it (rule T6)
  controls derive the negative-control inventory from the test modules
  production measure how the project itself was made
  refs     audit every reference added for the thesis against a live record
  figures  copy every placed figure from the results tree it was drawn in
  claims   re-verify every load-bearing number, and that its chapter says it
  stitch   concatenate the chapters into thesis/thesis.md
  pdf      typeset it

  python scripts/s25_assemble.py                    # all stages
  python scripts/s25_assemble.py --only claims
  python scripts/s25_assemble.py --from stitch
  python scripts/s25_assemble.py --refs-offline     # re-audit from the archive

`guards` runs first and on every build, following the project's rule that a
self-test runs before anything is written: it breaks each of the build's own
checks in a sandboxed copy of `thesis/` and requires each to fire with the
message it is supposed to, and it verifies that it altered no committed file.

The rest is dependency order: the grouping is fixed before the figures
are collected, the figures exist before the text that places them is checked,
and the numbers are verified before a document quoting them is written. A
failing stage does not stop the ones after it that do not depend on it — a
thesis that builds and says which of its guards failed is worth more than no
thesis — but it does set the exit status, so `s25_assemble.py` exits non-zero
on a missing figure, a missing chapter, a cited key with no reference row, an
unaudited new reference, a dropped glyph or a failed claim.
"""

from __future__ import annotations

import argparse

import s25_assign
import s25_claims
import s25_controls
import s25_figures
import s25_production
import s25_pdf
import s25_refs
import s25_stitch
import s25_test_guards

STAGES = ["guards", "assign", "controls", "production", "refs",
          "figures", "claims", "stitch", "pdf"]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=STAGES)
    ap.add_argument("--from", dest="from_stage", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--refs-offline", action="store_true",
                    help="re-check the reference audit from the archived "
                         "responses instead of fetching")
    args = ap.parse_args()

    if args.list:
        for s in STAGES:
            print(s)
        return 0

    todo = ([args.only] if args.only
            else STAGES[STAGES.index(args.from_stage):] if args.from_stage
            else list(STAGES))

    status = 0
    for stage in todo:
        if stage == "guards":
            status |= s25_test_guards.main()
        elif stage == "assign":
            status |= s25_assign.run(verbose=args.verbose)
        elif stage == "controls":
            status |= s25_controls.run()
        elif stage == "production":
            status |= s25_production.run()
        elif stage == "refs":
            status |= s25_refs.audit(offline=args.refs_offline,
                                     verbose=args.verbose)
        elif stage == "figures":
            status |= s25_figures.run()
        elif stage == "claims":
            status |= s25_claims.run(verbose=args.verbose)
        elif stage == "stitch":
            status |= s25_stitch.run(verbose=args.verbose)
        elif stage == "pdf":
            status |= s25_pdf.run()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
