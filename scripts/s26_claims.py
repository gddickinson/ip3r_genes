"""S26 stage `claims`: one ledger for the whole series, with a paper column.

The brief's rule is that a number quoted in two papers is checked once and
cannot disagree between them. So there is one ledger, not six. Each paper's
data module (`s26_claims_<id>.py`) declares two things:

* `CARRY`: identifiers of claims already declared elsewhere: in the
  manuscript's ledger (`C…`), the thesis's (`T…`), or another paper's
  module. A carried claim keeps its identifier, source, op and expected
  value. Only its paper list grows.
* `CLAIMS`: numbers no earlier ledger declares, with an identifier carrying
  the paper's prefix (`PREFIX`).

The merged ledger runs through `s14_claims.check`, the engine the manuscript
and the thesis already share (D72), with two additions:

* **Declared once.** Two identifiers with the same source, op, selector and
  expected value are the same claim declared twice, which is a build error.
  One of them has to be carried instead.
* **Stated where declared (D73).** Every paper in a claim's paper list must
  state its value, rendered the way `s25_claims._renderings` allows. The
  paper column cannot pick up a paper that never quotes the number.

Writes `papers/claims_check.tsv` (the whole series, one row per claim) and
`papers/<id>/claims_check.tsv` (the rows each paper states).
"""

from __future__ import annotations

import importlib
import json
import sys

import s14_claims
import s25_claims
import s26_assign
import s26_lib as L

PREFIX = {"range": "RA", "origin": "OR", "retention": "RE",
          "archive": "AR", "constraint": "CO", "ligand": "LI"}


def _module(pid: str):
    return importlib.import_module(f"s26_claims_{pid}")


def _signature(c: dict) -> str:
    return json.dumps([c["source"], c["op"], c.get("column"), c.get("key"),
                       sorted(c.get("where", {}).items()), c["expect"]])


def ledger() -> tuple[dict[str, dict], list[str]]:
    """id -> claim (with a `papers` list), and the ledger's own failures."""
    fails: list[str] = []
    upstream = {c["id"]: c for c in s14_claims.CLAIMS}
    for c in s25_claims.CLAIMS:
        upstream.setdefault(c["id"], c)
    merged: dict[str, dict] = {}
    carries: list[tuple[str, str]] = []
    for pid in L.SERIES:
        try:
            mod = _module(pid)
        except ModuleNotFoundError:
            fails.append(f"{pid}: no claims module s26_claims_{pid}.py")
            continue
        for c in mod.CLAIMS:
            if not c["id"].startswith(PREFIX[pid]):
                fails.append(f"{pid}: new claim {c['id']} lacks the paper's "
                             f"prefix {PREFIX[pid]}")
            if c["id"] in merged or c["id"] in upstream:
                fails.append(f"{pid}: claim id {c['id']} declared twice")
                continue
            merged[c["id"]] = dict(c, papers=[pid], origin=pid)
        carries += [(pid, cid) for cid in mod.CARRY]
    for pid, cid in carries:
        if cid in merged:
            if pid not in merged[cid]["papers"]:
                merged[cid]["papers"].append(pid)
        elif cid in upstream:
            src = upstream[cid]
            merged[cid] = {k: v for k, v in src.items()
                           if k not in ("chapter", "section")}
            merged[cid].update(papers=[pid], origin="C" if
                               cid.startswith("C") else "T")
        else:
            fails.append(f"{pid}: carries {cid}, which no ledger declares")
    seen: dict[str, str] = {}
    for cid, c in merged.items():
        sig = _signature(c)
        if sig in seen:
            fails.append(f"{cid} and {seen[sig]} are the same claim declared "
                         f"twice; carry one of them instead")
        seen[sig] = cid
        c["papers"].sort(key=L.number)
        c["primary"] = s26_assign.primary_paper(c["source"]) or ""
    return merged, fails


_TEXT: dict[str, str] = {}


def paper_text(pid: str) -> str:
    if pid not in _TEXT:
        _TEXT[pid] = "\n".join(L.section_text(pid, n) for n in L.SECTIONS)
    return _TEXT[pid]


def _stated(claim: dict) -> str:
    if claim["op"] == "grep":
        return ""               # the needle is a phrase in a report
    missing = [p for p in claim["papers"]
               if not any(r in paper_text(p)
                          for r in s25_claims._renderings(claim["expect"]))]
    if missing:
        return (f"{', '.join(missing)} does not state {claim['expect']!r}: "
                f"the number is missing from the text or the claim is "
                f"padding")
    return ""


def run(verbose: bool = False, only: str | None = None) -> int:
    _TEXT.clear()
    merged, fails = ledger()
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    claims = []
    for cid, c in merged.items():
        if only and only not in c["papers"]:
            continue
        claims.append(dict(c, section=",".join(c["papers"])))
    out = (L.PAPERS / "claims_check.tsv" if only is None
           else L.paper_dir(only) / "claims_check.tsv")
    n_fail, rows = s14_claims.check(claims, out, "s26 claims"
                                    + (f" {only}" if only else ""),
                                    verbose=verbose, extra=_stated)
    if only is None:
        for pid in L.SERIES:
            mine = [r for r in rows if pid in r["section"].split(",")]
            L.write_tsv(L.paper_dir(pid) / "claims_check.tsv", mine,
                        s14_claims.CHECK_FIELDS)
        shared = sum(1 for c in merged.values() if len(c["papers"]) > 1)
        print(f"[s26 claims] {len(merged)} claims across {len(L.SERIES)} "
              f"papers, {shared} stated in more than one paper and checked "
              f"once")
    return 1 if (n_fail or fails) else 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--paper")
    a = ap.parse_args()
    sys.exit(run(verbose=a.verbose, only=a.paper))
