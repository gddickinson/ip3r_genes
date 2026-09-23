"""S26 stage `rules`: P1, P2, P3 and P6, and the citation graph.

Each paper's configuration module declares its question, its controls, its
scope claims and its standalone answer. This stage holds each declaration
to its rule:

* **P1**: the question is one sentence ending in a question mark, with no
  "and" in it.
* **P2**: every declared control table exists and is primary in this paper.
* **P3**: every declared scope claim is in the ledger with this paper in its
  paper list (the claims stage then requires the paper to state it).
* **P6**: the standalone answer is written (at least 60 words), and every
  claim it rests on has its source in a results entry primary in this paper.
  A paper whose standalone claims rest on a sibling's table does not stand
  alone, and that is the one form of "slice" a build can detect.
* **The graph**: edges are read from `{paper:<id>}` citations in the
  sources, not declared. The graph must have no cycles, and a paper may cite
  only papers earlier in `s26_lib.SERIES`.

Writes `papers/dependency_graph.tsv` and `papers/standalone.tsv`.
"""

from __future__ import annotations

import re
import sys

import s26_assign
import s26_claims
import s26_lib as L

MIN_STANDALONE_WORDS = 60


def edges() -> list[dict]:
    rows = []
    for pid in L.SERIES:
        text = "\n".join(L.section_text(pid, n) for n in L.SECTIONS)
        counts: dict[str, int] = {}
        for other in L.PAPER_REF.findall(text):
            counts[other] = counts.get(other, 0) + 1
        for other, n in sorted(counts.items()):
            rows.append({
                "citing": pid, "citing_n": L.number(pid), "cited": other,
                "cited_n": L.number(other) if other in L.SERIES else "",
                "mentions": n,
                "order_ok": int(other in L.SERIES
                                and L.number(other) < L.number(pid)),
            })
    return rows


def _cycles(rows: list[dict]) -> list[str]:
    graph: dict[str, set[str]] = {p: set() for p in L.SERIES}
    for r in rows:
        if r["cited"] in graph:
            graph[r["citing"]].add(r["cited"])
    fails, state = [], {p: 0 for p in graph}

    def visit(p: str, path: list[str]) -> None:
        state[p] = 1
        for q in sorted(graph[p]):
            if state[q] == 1:
                fails.append("citation cycle: " + " -> ".join(
                    path[path.index(q):] + [q]) if q in path else q)
            elif state[q] == 0:
                visit(q, path + [q])
        state[p] = 2
    for p in L.SERIES:
        if state[p] == 0:
            visit(p, [p])
    return fails


def check() -> tuple[list[str], list[dict], list[dict]]:
    fails: list[str] = []
    merged, _ = s26_claims.ledger()
    graph = edges()
    for r in graph:
        if not r["order_ok"]:
            fails.append(f"{r['citing']} cites {r['cited']}, which is not "
                         f"earlier in the submission order")
    fails += _cycles(graph)
    table = []
    for pid in L.SERIES:
        cfg = L.paper(pid)
        q = cfg.QUESTION.strip()
        if not q.endswith("?") or q.count("?") != 1:
            fails.append(f"{pid}: the question is not one sentence ending "
                         f"in '?' (P1)")
        if re.search(r"\band\b", q, re.I):
            fails.append(f"{pid}: the question contains 'and' (P1): {q}")
        for what, path in cfg.CONTROLS:
            if not (L.ROOT / path).exists():
                fails.append(f"{pid}: control {what!r} names {path}, which "
                             f"does not exist")
            owner = s26_assign.primary_paper(path)
            if owner != pid:
                fails.append(f"{pid}: control {what!r} lives in {path}, "
                             f"primary in {owner!r} (P2)")
        for cid in cfg.SCOPE_CLAIMS:
            if cid not in merged or pid not in merged[cid]["papers"]:
                fails.append(f"{pid}: scope claim {cid} is not in the "
                             f"ledger against this paper (P3)")
        words = len(cfg.STANDALONE.split())
        if words < MIN_STANDALONE_WORDS:
            fails.append(f"{pid}: the standalone answer is {words} words; "
                         f"P6 needs one written out")
        if not cfg.STANDALONE_CLAIMS:
            fails.append(f"{pid}: the standalone answer rests on no claim")
        for cid in cfg.STANDALONE_CLAIMS:
            c = merged.get(cid)
            if c is None or pid not in c["papers"]:
                fails.append(f"{pid}: standalone claim {cid} is not in the "
                             f"ledger against this paper (P6)")
            elif c["primary"] != pid:
                fails.append(f"{pid}: standalone claim {cid} rests on "
                             f"{c['source']}, primary in {c['primary']!r}: "
                             f"the paper does not stand alone (P6)")
        cites = [r["cited"] for r in graph if r["citing"] == pid]
        cited_by = [r["citing"] for r in graph if r["cited"] == pid]
        table.append({
            "number": L.number(pid), "paper": pid, "title": cfg.TITLE,
            "question": q, "controls": len(cfg.CONTROLS),
            "standalone_claims": ",".join(cfg.STANDALONE_CLAIMS),
            "cites": ",".join(cites), "cited_by": ",".join(cited_by),
            "standalone": " ".join(cfg.STANDALONE.split()),
        })
    return fails, graph, table


def run(verbose: bool = False) -> int:
    fails, graph, table = check()
    L.write_tsv(L.PAPERS / "dependency_graph.tsv", graph,
                ["citing", "citing_n", "cited", "cited_n", "mentions",
                 "order_ok"])
    L.write_tsv(L.PAPERS / "standalone.tsv", table,
                ["number", "paper", "title", "question", "controls",
                 "standalone_claims", "cites", "cited_by", "standalone"])
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    print(f"[s26 rules] {len(L.SERIES)} papers, {len(graph)} citation "
          f"edge(s), P1/P2/P3/P6 and the graph: {len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(run())
