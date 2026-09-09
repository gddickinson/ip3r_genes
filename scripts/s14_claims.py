"""S14 stage `claims` — re-verify the manuscript's load-bearing numbers.

Every headline number in the manuscript is declared here with the committed
table it comes from and the operation that recovers it. Running this script
re-reads those tables and reports agreement, so the manuscript cannot drift
away from its own data — the same guarantee the per-analysis reports already
have, extended to the paper.

  python scripts/s14_claims.py            # check, print a summary
  python scripts/s14_claims.py --verbose  # print every claim

Writes `manuscript/claims_check.tsv`. Exit status is non-zero if any claim
fails or its source is missing.
"""

from __future__ import annotations

import argparse
import sys

import s14_claims_function
import s14_claims_history
import s14_claims_scope
import s14_lib as lib

# op:
#   count    - number of rows matching `where`
#   cell     - value of `column` in the single row matching `where`
#   nunique  - distinct values of `column` among rows matching `where`
#   sum      - total of `column` over rows matching `where`
#   json     - dotted key path `key` in a JSON file
#   grep     - `expect` appears literally in a text file (needle = expect)
#
# `where` keys support `col`, `col__startswith`, `col__in` (comma-separated).
#: The ledger itself lives in three data modules, one per part of the paper,
#: so no file exceeds the project's 500-line budget. Concatenating them here
#: keeps a single ordered list and a single `claims_check.tsv`.
CLAIMS: list[dict] = (
    s14_claims_scope.CLAIMS
    + s14_claims_history.CLAIMS
    + s14_claims_function.CLAIMS
)

def _matches(row: dict, where: dict) -> bool:
    for key, want in where.items():
        if key.endswith("__startswith"):
            if not (row.get(key[:-12]) or "").startswith(want):
                return False
        elif key.endswith("__in"):
            if (row.get(key[:-4]) or "") not in want.split(","):
                return False
        elif (row.get(key) or "") != want:
            return False
    return True


def _dig(obj, dotted: str):
    for part in dotted.split("."):
        obj = obj[part]
    return obj


def _evaluate(claim: dict) -> tuple[str, str]:
    """Return (got, error). `got` is '' when error is set."""
    path = lib.ROOT / claim["source"]
    if not path.exists():
        return "", "source missing"
    op = claim["op"]
    try:
        if op == "json":
            return str(_dig(lib.read_json(path), claim["key"])), ""
        if op == "grep":
            text = path.read_text(encoding="utf-8", errors="replace")
            return (claim["expect"] if claim["expect"] in text
                    else "<not found>"), ""
        rows = [r for r in lib.read_tsv(path)
                if _matches(r, claim.get("where", {}))]
        if op == "count":
            return str(len(rows)), ""
        if op == "nunique":
            return str(len({r[claim["column"]] for r in rows})), ""
        if op == "sum":
            total = sum(float(r[claim["column"]] or 0) for r in rows)
            return str(int(total) if total == int(total) else total), ""
        if op == "cell":
            if len(rows) != 1:
                return "", f"expected 1 matching row, got {len(rows)}"
            return (rows[0].get(claim["column"]) or ""), ""
    except KeyError as exc:
        return "", f"missing column/key {exc}"
    except Exception as exc:                      # noqa: BLE001 - reported
        return "", f"{type(exc).__name__}: {exc}"
    return "", f"unknown op {op}"


def _agrees(expect: str, got: str, tol: float) -> bool:
    if expect == got:
        return True
    try:
        return abs(float(expect) - float(got)) <= tol
    except ValueError:
        return False


def run(verbose: bool = False) -> int:
    rows, n_fail = [], 0
    for claim in CLAIMS:
        got, err = _evaluate(claim)
        ok = (not err) and _agrees(claim["expect"], got, claim.get("tol", 0.0))
        n_fail += 0 if ok else 1
        rows.append({
            "id": claim["id"],
            "section": claim["section"],
            "claim": claim["claim"],
            "expected": claim["expect"],
            "found": got,
            "verdict": "ok" if ok else "MISMATCH",
            "note": err,
            "source": claim["source"],
            "op": claim["op"],
        })
        if verbose or not ok:
            mark = "ok  " if ok else "FAIL"
            print(f"  [{mark}] {claim['id']} {claim['claim']}: "
                  f"expected {claim['expect']!r}, found {got!r}"
                  f"{' (' + err + ')' if err else ''}")

    lib.write_tsv(lib.MS / "claims_check.tsv", rows,
                  ["id", "section", "claim", "expected", "found", "verdict",
                   "note", "source", "op"])
    print(f"[s14 claims] {len(rows) - n_fail}/{len(rows)} load-bearing numbers "
          f"re-verified against the committed tables "
          f"-> manuscript/claims_check.tsv")
    n_fail += _checklist_total_matches(len(rows))
    return 1 if n_fail else 0


def _checklist_total_matches(total: int) -> int:
    """The reviewer checklist is hand-written; hold its headline count to
    the ledger's.

    Every other number in this package is generated. This one is typed, and
    at S14c it was four editions out of date -- it claimed 175 declared
    claims against a ledger of 180, which is exactly the kind of drift the
    ledger exists to prevent everywhere else. So the one number that says how
    many claims there are is itself checked.
    """
    path = lib.MS / "reviewer_checklist.md"
    if not path.exists():
        return 0
    text = path.read_text(encoding="utf-8")
    if f"{total} declared claims" in text and f"{total}/{total} pass" in text:
        return 0
    print(f"  [FAIL] reviewer_checklist.md does not say "
          f"'{total} declared claims' and '{total}/{total} pass' — the "
          f"hand-written self-audit has drifted from the ledger",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    sys.exit(run(verbose=args.verbose))
