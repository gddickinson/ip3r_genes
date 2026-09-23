"""S26 driver: build the paper series.

Stages, in order:

  guards   every build guard broken on purpose, in a sandbox (D74)
  assign   P5: every results entry primary in one paper or excluded, and
           the cross-check against the thesis's assignment
  rules    P1, P2, P3, P6 and the citation graph
  figures  per paper: number by first mention, check P4/P5/R3, copy
  prose    per paper: the thesis's mechanical prose rules
  stitch   per paper: resolve figures, companions and citations
  claims   one ledger across the series, with a paper column
  pdf      per paper: typeset, with the thesis's page checks
  readme   papers/README.md, rendered from the series' committed tables
  deposit  per paper: the deposit manifest

  python scripts/s26_assemble.py                    # everything
  python scripts/s26_assemble.py --paper origin     # one paper's stages
  python scripts/s26_assemble.py --only claims
  python scripts/s26_assemble.py --from stitch --no-pdf

Series-level stages (`assign`, `rules`, `claims`) always run over the whole
series, because what they check is how the papers relate to one another; with
`--paper`, `claims` also writes that paper's own check table. Exit status is
non-zero if any stage failed. A failed stage does not stop the stages after
it, so one run reports every problem.
"""

from __future__ import annotations

import argparse
import sys

import s26_assign
import s26_claims
import s26_deposit
import s26_figures
import s26_lib as L
import s26_pdf
import s26_prose
import s26_readme
import s26_rules
import s26_stitch

STAGES = ["guards", "assign", "rules", "figures", "prose", "stitch",
          "claims", "pdf", "readme", "deposit"]
PER_PAPER = {"figures": s26_figures.run, "prose": s26_prose.run,
             "stitch": s26_stitch.run, "pdf": s26_pdf.run,
             "deposit": s26_deposit.run}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=STAGES)
    ap.add_argument("--from", dest="from_stage", choices=STAGES)
    ap.add_argument("--paper", choices=L.SERIES)
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    if a.list:
        print("\n".join(STAGES))
        return 0
    todo = ([a.only] if a.only else
            STAGES[STAGES.index(a.from_stage):] if a.from_stage else
            list(STAGES))
    if a.no_pdf and "pdf" in todo:
        todo.remove("pdf")
    if a.paper and "guards" in todo and not a.only:
        todo.remove("guards")          # the sandbox rebuilds the series
    papers = [a.paper] if a.paper else list(L.SERIES)
    status = 0
    for stage in todo:
        if stage == "guards":
            import s26_test_guards
            status |= s26_test_guards.run()
        elif stage == "assign":
            status |= s26_assign.run()
        elif stage == "rules":
            status |= s26_rules.run()
        elif stage == "readme":
            status |= s26_readme.run()
        elif stage == "claims":
            status |= s26_claims.run(verbose=a.verbose)
            if a.paper:
                status |= s26_claims.run(verbose=a.verbose, only=a.paper)
        else:
            for pid in papers:
                status |= PER_PAPER[stage](pid)
    print(f"[s26] {'ok' if status == 0 else 'FAILED'}")
    return status


if __name__ == "__main__":
    sys.exit(main())
