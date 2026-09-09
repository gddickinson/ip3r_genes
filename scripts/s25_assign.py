"""S25 stage `assign` — commit the chapter grouping, and guard it.

Brief step 1: *fix the chapter grouping and commit it as a table with the rule
that assigned each results directory to a chapter, before writing any prose.*

This module writes `thesis/chapter_assignment.tsv` from the declaration in
`s25_lib.ASSIGNMENT` and fails on four things, each of which is a way the
grouping could quietly stop being a grouping:

1. a results-tree entry that is neither assigned nor explicitly excluded —
   rule **T6**, enforced, so a chapter cannot become a bag of leftovers by
   accident and a new task's results cannot be silently left out of the
   thesis;
2. an assignment naming an entry that no longer exists (a renamed results
   directory);
3. an assignment naming a chapter with no file in `CHAPTERS`;
4. a chapter that is primary for nothing — which under **T1** means it is
   exposition rather than a chapter, and has to say so.

The fourth is a warning rather than an error: Chapters 1 and 15 are
introduction and discussion and are primary for nothing by design.
"""

from __future__ import annotations

import sys

import s25_lib as L

FIELDS = ["results_entry", "chapter", "chapter_title", "rule", "rule_text",
          "what_it_is", "why_not_elsewhere"]

#: Chapters that are primary for no results directory by design.
EXPOSITION_CHAPTERS = {1, 15}


def build() -> tuple[list[dict], list[str]]:
    rows, problems = [], []
    present = set(L.results_entries())
    assigned = set(L.ASSIGNMENT) | set(L.ASSIGNMENT_EXCLUDED)

    for entry in sorted(present - assigned):
        problems.append(
            f"results/{entry} is neither assigned to a chapter nor listed in "
            f"ASSIGNMENT_EXCLUDED (rule T6)")
    for entry in sorted(assigned - present):
        problems.append(f"results/{entry} is assigned but does not exist")

    for entry, (chapter, rule, what, why) in sorted(L.ASSIGNMENT.items()):
        if chapter not in L.FILE_OF_CHAPTER:
            problems.append(f"results/{entry} is assigned to chapter "
                            f"{chapter}, which has no file")
            continue
        if rule not in L.RULES:
            problems.append(f"results/{entry} cites rule {rule}, which is "
                            f"not in the rule vocabulary")
            continue
        rows.append({
            "results_entry": f"results/{entry}",
            "chapter": str(chapter),
            "chapter_title": L.chapter_title(chapter),
            "rule": rule,
            "rule_text": L.RULES[rule],
            "what_it_is": what,
            "why_not_elsewhere": why,
        })

    for entry, reason in sorted(L.ASSIGNMENT_EXCLUDED.items()):
        rows.append({
            "results_entry": f"results/{entry}",
            "chapter": "-",
            "chapter_title": "not in the thesis",
            "rule": "-",
            "rule_text": "-",
            "what_it_is": reason,
            "why_not_elsewhere": "-",
        })

    return rows, problems


def run(verbose: bool = False) -> int:
    rows, problems = build()
    L.write_tsv(L.TH / "chapter_assignment.tsv", rows, FIELDS)

    covered = {int(r["chapter"]) for r in rows if r["chapter"] != "-"}
    barren = sorted(set(L.FILE_OF_CHAPTER) - covered - EXPOSITION_CHAPTERS)
    for n in barren:
        print(f"  NOTE chapter {n} ({L.chapter_title(n)}) is primary for no "
              f"results directory", file=sys.stderr)

    for p in problems:
        print(f"  [FAIL] {p}", file=sys.stderr)
    if verbose:
        for r in rows:
            print(f"  {r['results_entry']:<46} ch {r['chapter']:<3} "
                  f"{r['rule']}")
    n_ass = len([r for r in rows if r["chapter"] != "-"])
    print(f"[s25 assign] {n_ass} results directories assigned across "
          f"{len(covered)} chapters, {len(L.ASSIGNMENT_EXCLUDED)} excluded "
          f"-> thesis/chapter_assignment.tsv")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(run(verbose="--verbose" in sys.argv))
