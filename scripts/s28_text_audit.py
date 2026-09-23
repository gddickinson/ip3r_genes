"""S28 — `thesis/figure_text_audit.tsv`: every legend or label changed.

One row per legend, bold lead-in or in-figure text that S28 changed, with
what it said before and what it says now, so the revision is a table rather
than a claim. Three sources, each diffed the same way:

* **legends** — every `**{fig:slug}.**` paragraph in the chapter sources at
  the commit S28 started from (`git show <ref>:thesis/<file>`) against the
  working tree, keyed by slug;
* **lead-ins** — every bold paragraph opener, paired by file and by closest
  match (difflib) so a rewritten lead-in meets its original rather than its
  neighbour;
* **in-figure text** — the text inventory `figcheck` logged before any figure
  module was edited (committed as `thesis/figure_text_inventory_before.tsv`)
  against the inventory of the final render, per figure stem, paired by
  role and closest match.

  python scripts/s28_text_audit.py --after <inventory.tsv> [--ref <commit>]
"""

from __future__ import annotations

import argparse
import difflib
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s25_lib as L                                          # noqa: E402

BEFORE_INVENTORY = L.TH / "figure_text_inventory_before.tsv"
OUT = L.TH / "figure_text_audit.tsv"
FIELDS = ["kind", "where", "before", "after"]

LEGEND_RE = re.compile(r"\*\*\{fig:([a-z0-9_]+)\}\.\*\*(.*)", re.S)
LEADIN_RE = re.compile(r"^\*\*(.+?)\*\*", re.S)


def _git_show(ref: str, path: str) -> str:
    out = subprocess.run(["git", "show", f"{ref}:{path}"], cwd=L.ROOT,
                         capture_output=True, text=True)
    return out.stdout if out.returncode == 0 else ""


def _legends(text: str) -> dict[str, str]:
    out = {}
    for para in text.split("\n\n"):
        m = LEGEND_RE.match(para.lstrip())
        if m:
            out[m.group(1)] = " ".join(m.group(2).split())
    return out


def _leadins(text: str) -> list[str]:
    out = []
    for para in text.split("\n\n"):
        p = para.strip()
        if p.startswith("**") and not p.startswith("**{fig:"):
            m = LEADIN_RE.match(p)
            if m:
                out.append(" ".join(m.group(1).split()))
    return out


def _pair(before: list[str], after: list[str], floor: float = 0.3):
    """Pair each removed string with its closest added one, greedily."""
    rows, used = [], set()
    for b in before:
        best, score = None, floor
        for j, a in enumerate(after):
            if j in used:
                continue
            r = difflib.SequenceMatcher(None, b, a).ratio()
            if r > score:
                best, score = j, r
        if best is None:
            rows.append((b, ""))
        else:
            used.add(best)
            rows.append((b, after[best]))
    rows += [("", a) for j, a in enumerate(after) if j not in used]
    return rows


def _inventory(path: Path) -> dict[str, list[tuple[str, str]]]:
    """figure stem -> [(role, text)], problem rows excluded."""
    out: dict[str, list[tuple[str, str]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 6 or parts[1] == "PROBLEM":
            continue
        stem, role, text = parts[0], parts[1], parts[5]
        if role in ("tick", "axis", "legend", "off"):
            continue
        out.setdefault(stem, []).append((role, text))
    return out


def rows(after_inventory: Path, ref: str) -> list[dict]:
    out: list[dict] = []
    for name in L.CHAPTER_FILES:
        old = _git_show(ref, f"thesis/{name}")
        new = (L.TH / name).read_text(encoding="utf-8")
        lo, ln = _legends(old), _legends(new)
        for slug in sorted(set(lo) | set(ln)):
            if lo.get(slug, "") != ln.get(slug, ""):
                out.append({"kind": "legend", "where": f"{name}:{slug}",
                            "before": lo.get(slug, ""),
                            "after": ln.get(slug, "")})
        bo, bn = _leadins(old), _leadins(new)
        removed = [x for x in bo if x not in bn]
        added = [x for x in bn if x not in bo]
        for b, a in _pair(removed, added):
            out.append({"kind": "lead-in", "where": name, "before": b,
                        "after": a})
    before, after = _inventory(BEFORE_INVENTORY), _inventory(after_inventory)
    for stem in sorted(set(before) | set(after)):
        b_rows, a_rows = before.get(stem, []), after.get(stem, [])
        for role in sorted({r for r, _t in b_rows} | {r for r, _t in a_rows}):
            bt = [t for r, t in b_rows if r == role]
            at = [t for r, t in a_rows if r == role]
            removed = [t for t in bt if t not in at]
            added = [t for t in at if t not in bt]
            for b, a in _pair(removed, added):
                out.append({"kind": f"figure {role}", "where": stem,
                            "before": b, "after": a})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--after", type=Path, required=True,
                    help="figcheck text inventory of the final render")
    ap.add_argument("--ref", default="HEAD",
                    help="the commit S28 started from")
    args = ap.parse_args()
    table = rows(args.after, args.ref)
    L.write_tsv(OUT, table, FIELDS)
    kinds = {}
    for r in table:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
    print(f"[s28 audit] {len(table)} changes -> {OUT.relative_to(L.ROOT)}: "
          + ", ".join(f"{k} {n}" for k, n in sorted(kinds.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
