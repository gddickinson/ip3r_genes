"""S14 — resolve the manuscript's stable citation keys into numbered ones.

The section files cite with the project's stable keys (`[R22]`, `[R22, R25]`)
against `results/s0_baseline/references.tsv`, the same table the literature
review is built from. The *stitch* stage renumbers them into order of first
appearance and renders the bibliography, so the text and the bibliography
cannot drift and no author, title or year is ever retyped into a section file
(D13 applied to the reference list).

A cited key with no row in the table is a build error, as it is in
`scripts/s0_review_build.py`.
"""

from __future__ import annotations

import re

import s14_lib as lib

REFS_TSV = lib.RESULTS / "s0_baseline" / "references.tsv"

#: `[R22]` or `[R22, R25]` — keys only, so an ordinary bracket is untouched.
CITE = re.compile(r"\[(R\d{2,3}(?:\s*,\s*R\d{2,3})*)\]")

#: Where the bibliography is spliced into the stitched manuscript.
MARKER = "{references}"


def load_refs() -> dict[str, dict[str, str]]:
    return {r["ref_id"]: r for r in lib.read_tsv(REFS_TSV)}


def _format(num: int, row: dict[str, str]) -> str:
    bits = [f"{num}.", row["authors"] + ".", row["title"].rstrip(".") + "."]
    bits.append(f"*{row['journal']}* {row['year']}.")
    if row.get("doi"):
        bits.append(f"doi:{row['doi']}.")
    if row.get("pmid"):
        bits.append(f"PMID {row['pmid']}.")
    bits.append(f"[{row['ref_id']}]")
    return " ".join(bits)


def resolve(text: str) -> tuple[str, str, dict]:
    """Renumber every citation and render the bibliography.

    Returns (text with numbered citations, bibliography markdown, stats).
    Raises KeyError naming every cited key with no reference row.
    """
    refs = load_refs()
    order: list[str] = []

    def number(key: str) -> int:
        if key not in order:
            order.append(key)
        return order.index(key) + 1

    missing: list[str] = []

    def repl(m: re.Match) -> str:
        keys = [k.strip() for k in m.group(1).split(",")]
        for k in keys:
            if k not in refs:
                missing.append(k)
        return "[" + ", ".join(str(number(k)) for k in keys) + "]"

    out = CITE.sub(repl, text)
    if missing:
        raise KeyError("cited keys with no row in references.tsv: "
                       + ", ".join(sorted(set(missing))))

    lines = [_format(i + 1, refs[k]) for i, k in enumerate(order)]
    bib = "\n\n".join(lines)
    return out, bib, {"n_cited": len(order), "n_available": len(refs),
                      "keys": list(order)}
