"""S26 stage `assign`: enforce P5 and write `papers/paper_assignment.tsv`.

Checks, each a build failure:

* an entry under `results/` that no row assigns or excludes (P5, strong
  direction: a future task's results cannot be left out silently);
* a row naming an entry that does not exist, an unknown paper or rule;
* a split directory (P5s) where a file matches no group or more than one;
* a departure from `thesis/chapter_assignment.tsv` that is not declared in
  `s26_assignment.DEPARTURES`, or a declared departure that has gone stale.
  This is the cross-check between the two documents' tables that the
  2026-09-09 emergent row asked for.

`primary_paper(path)` is the lookup the other stages use to enforce P2, P5
and P6: it answers which paper a committed file is primary in.
"""

from __future__ import annotations

import argparse
import sys
from fnmatch import fnmatch

import s26_assignment as A
import s26_lib as L

RULES = {"P1", "P2", "P3", "P5s", "excluded"}
FIELDS = ["results_entry", "files", "paper", "rule", "what", "nearest_paper",
          "why_not_there", "thesis_chapter", "thesis_paper", "departure"]


def _files(entry: str) -> list[str]:
    base = L.ROOT / entry
    if base.is_file():
        return [""]
    return sorted(str(p.relative_to(base)) for p in base.rglob("*")
                  if p.is_file() and p.name != ".DS_Store")


def _rows_by_entry() -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = {}
    for row in A.ASSIGNMENT:
        out.setdefault(row[0], []).append(row)
    return out


def primary_paper(path: str) -> str | None:
    """Paper a repo-relative results path is primary in; '-' if excluded,
    None if the path is outside results/ or unassigned."""
    parts = path.split("/")
    if len(parts) < 2 or parts[0] != "results":
        return None
    entry = "/".join(parts[:2])
    inner = "/".join(parts[2:])
    for e, glob, pid, *_ in A.ASSIGNMENT:
        if e == entry and (glob == "*" or fnmatch(inner, glob)):
            return pid
    return None


def _thesis_chapters() -> dict[str, int | None]:
    path = L.ROOT / "thesis" / "chapter_assignment.tsv"
    out = {}
    for r in L.read_tsv(path):
        ch = r["chapter"]
        out[r["results_entry"]] = int(ch) if ch.isdigit() else None
    return out


def build() -> tuple[list[dict], list[str]]:
    fails: list[str] = []
    by_entry = _rows_by_entry()
    entries = L.results_entries()
    thesis = _thesis_chapters()

    for entry in by_entry:
        if entry not in entries:
            fails.append(f"{entry}: assigned but does not exist")
    for e, glob, pid, rule, *_ in A.ASSIGNMENT:
        if pid != "-" and pid not in L.SERIES:
            fails.append(f"{e} [{glob}]: unknown paper {pid!r}")
        if rule not in RULES:
            fails.append(f"{e} [{glob}]: unknown rule {rule!r}")
        if (pid == "-") != (rule == "excluded"):
            fails.append(f"{e} [{glob}]: an exclusion must use rule "
                         f"'excluded' and nothing else may")

    rows: list[dict] = []
    for entry in entries:
        groups = by_entry.get(entry)
        if not groups:
            fails.append(f"{entry}: neither assigned to a paper nor excluded "
                         f"(P5)")
            continue
        split = len(groups) > 1
        if split:
            for g in groups:
                if g[1] == "*" or g[3] not in ("P5s", "excluded"):
                    fails.append(f"{entry}: a split directory needs file "
                                 f"globs under P5s, got {g[1]!r}/{g[3]}")
            for f in _files(entry):
                hits = [g for g in groups if fnmatch(f, g[1])]
                if len(hits) != 1:
                    fails.append(f"{entry}/{f}: matches {len(hits)} file "
                                 f"groups, must match exactly one (P5s)")
        ours = groups[0][2] if not split else "split"
        ch = thesis.get(entry)
        thesis_paper = (A.THESIS_CHAPTER_PAPER.get(ch, "?") if ch is not None
                        else None)
        thesis_paper = thesis_paper or "-"
        differs = ours != thesis_paper and not (
            not split and ours == "-" and thesis_paper == "-")
        declared = entry in A.DEPARTURES
        if differs and not declared:
            fails.append(f"{entry}: thesis places it in chapter {ch} "
                         f"(-> {thesis_paper}), the series in {ours}, and "
                         f"the departure is not declared")
        if declared and not differs:
            fails.append(f"{entry}: declared as a departure from the thesis "
                         f"but the two assignments agree")
        for e, glob, pid, rule, what, near, why in groups:
            rows.append({
                "results_entry": entry, "files": glob, "paper": pid,
                "rule": rule, "what": what, "nearest_paper": near,
                "why_not_there": why,
                "thesis_chapter": "" if ch is None else ch,
                "thesis_paper": thesis_paper,
                "departure": A.DEPARTURES.get(entry, ""),
            })
    for entry in A.DEPARTURES:
        if entry not in by_entry:
            fails.append(f"{entry}: declared departure for an entry the "
                         f"series does not assign")
    return rows, fails


def run(verbose: bool = False) -> int:
    rows, fails = build()
    L.write_tsv(L.PAPERS / "paper_assignment.tsv", rows, FIELDS)
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    by_paper: dict[str, int] = {}
    for r in rows:
        by_paper[r["paper"]] = by_paper.get(r["paper"], 0) + 1
    entries = {r["results_entry"] for r in rows}
    n_split = len({r["results_entry"] for r in rows if r["files"] != "*"})
    print(f"[s26 assign] {len(entries)} results entries, {len(rows)} rows "
          f"({n_split} split by P5s), "
          + ", ".join(f"{p}:{by_paper.get(p, 0)}" for p in L.SERIES + ["-"])
          + f", {len(A.DEPARTURES)} declared departures from the thesis: "
          f"{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    sys.exit(run(verbose=ap.parse_args().verbose))
