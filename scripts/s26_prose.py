"""S26 stage `prose`: the thesis's mechanical prose rules, applied per paper.

D78: a style rule is only applied when it is checked. The rules are the
thesis's, and so is the code: every section file of a paper goes through
`s25_prose._check_file` unchanged (no em-dash; a legend opens with a
statement and every sentence in it is one; a legend never says the same
thing twice or restates its neighbour; a bold lead-in is a statement), with
`prose_lex` as the one finite-verb detector. R3, that every figure is
described in the body, is enforced in `s26_figures`.

One check is added for papers. **A results subheading is a statement**
(D76): a `###` heading in the results section must carry a finite verb, since
a reader of a paper sees its results headings in a table of contents first.
"""

from __future__ import annotations

import re
import sys

import prose_lex as P
import s25_prose
import s26_lib as L

HEADING = re.compile(r"^###\s+(.+)$", re.M)


def check(pid: str, verbose: bool = False) -> tuple[list[str], int]:
    fails: list[str] = []
    n_legends = 0
    for name in L.SECTIONS:
        path = L.paper_dir(pid) / name
        if not path.exists():
            continue
        f, _refs = s25_prose._check_file(path, verbose)
        fails += [f"{pid}/{x}" for x in f]
        text = path.read_text(encoding="utf-8")
        n_legends += sum(1 for p in text.split("\n\n")
                         if L.LEGEND_HEAD.match(p.strip()))
        if name == "02_results.md":
            for head in HEADING.findall(text):
                if not P.has_finite_verb(head.strip()):
                    fails.append(f"{pid}/{name}: results heading has no "
                                 f"finite verb: '{head[:70]}'")
    return fails, n_legends


def run(pid: str, verbose: bool = False) -> int:
    fails, n = check(pid, verbose)
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    print(f"[s26 prose {pid}] {n} legends; no em-dash, no repetition, "
          f"every legend, lead-in and results heading a statement: "
          f"{len(fails)} failure(s)")
    return 1 if fails else 0
