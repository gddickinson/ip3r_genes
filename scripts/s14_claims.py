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
CLAIMS: list[dict] = [
    # Every load-bearing number in the manuscript gets a row here, naming the
    # committed table it comes from and the operation that recovers it. A
    # re-run that changes a table then fails loudly instead of leaving the
    # text stale. Add rows as the manuscript is written — an example of each
    # op, for reference:
    #
    # dict(id="C01", section="R1",
    #      claim="the sweep produced N genome x paralog cells",
    #      source="results/genome_ledger.tsv", op="count", expect="582"),
    # dict(id="C02", section="scope",
    #      claim="the declared genome scope is N assemblies",
    #      source="results/genome_manifest.tsv", op="nunique",
    #      column="accession", expect="194"),
    # dict(id="C03", section="R1",
    #      claim="ITPR2 is absent from N genomes",
    #      source="results/genome_ledger.tsv", op="count",
    #      where={"status": "absent", "paralog": "ITPR2"}, expect="0"),
    # dict(id="C04", section="R2",
    #      claim="the census holds N proteins",
    #      source="results/census_v5/census_v5_stats.json", op="json",
    #      key="v5_rows", expect="8329"),
    # dict(id="C05", section="R3",
    #      claim="mean omega on the ITPR3 stem",
    #      source="results/selection/omega_table.tsv", op="cell",
    #      where={"branch": "ITPR3_stem"}, column="omega",
    #      expect="0.0412", tol=1e-4),
]

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
    return 1 if n_fail else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    sys.exit(run(verbose=args.verbose))
