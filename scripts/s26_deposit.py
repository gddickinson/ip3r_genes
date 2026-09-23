"""S26 stage `deposit`: one paper's deposit manifest.

A paper deposits the results it is primary in and nothing it only cites, so
the results each paper deposits are the P5 assignment read from the other
side. For a directory split under P5s, only this paper's file groups go in.
Each deposit also carries the shared inputs every paper reads (the reference
table, the toolchain manifest), the scripts that regenerate its tables, and
the paper's own directory.

Rows, sizes and SHA-256 are recorded as in `s14_deposit`, and the notes
(including the bulk classes excluded, with the command that regenerates
each) are rendered by `s14_deposit._notes`, the manuscript's own renderer.
"""

from __future__ import annotations

import sys

import s14_deposit
import s14_lib
import s26_assign
import s26_lib as L

SHARED = ["results/s0_baseline/references.tsv",
          "results/toolchain_manifest.txt"]


def files(pid: str):
    for entry in L.results_entries():
        base = L.ROOT / entry
        paths = [base] if base.is_file() else sorted(base.rglob("*"))
        for path in paths:
            if not path.is_file() or path.name in s14_lib.DEPOSIT_SKIP_NAMES:
                continue
            rel = str(path.relative_to(L.ROOT))
            if s26_assign.primary_paper(rel) == pid:
                yield path
    for rel in SHARED:
        yield L.ROOT / rel
    for path in sorted((L.ROOT / "scripts").glob("*.py")):
        yield path
    for path in sorted(L.paper_dir(pid).rglob("*")):
        if (path.is_file() and path.name != "deposit_manifest.tsv"
                and not path.name.startswith(".")):
            yield path


def run(pid: str) -> int:
    rows, seen = [], set()
    for path in files(pid):
        rel = path.relative_to(L.ROOT)
        if rel in seen:
            continue
        seen.add(rel)
        rows.append({"path": str(rel),
                     "group": s14_lib.deposit_group(rel.parts),
                     "bytes": path.stat().st_size,
                     "sha256": L.sha256(path)})
    total = sum(r["bytes"] for r in rows)
    d = L.paper_dir(pid)
    L.write_tsv(d / "deposit_manifest.tsv", rows,
                ["path", "group", "bytes", "sha256"])
    (d / "deposit_notes.md").write_text(
        s14_deposit._notes(rows, total, generator="scripts/s26_deposit.py"),
        encoding="utf-8")
    n_results = sum(1 for r in rows if r["path"].startswith("results/"))
    print(f"[s26 deposit {pid}] {len(rows)} files ({n_results} results), "
          f"{s14_lib.size_str(total)} -> papers/{pid}/deposit_manifest.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(run(sys.argv[1]))
