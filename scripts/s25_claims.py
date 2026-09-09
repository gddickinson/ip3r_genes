"""S25 stage `claims` — the thesis's load-bearing numbers, checked twice.

The manuscript's ledger re-verifies every headline number against the
committed table it comes from. A thesis is a larger surface for drift, not a
licence for less checking, so this ledger does the same — through the *same*
engine (`s14_claims.check`), not a fork, so a number quoted in both documents
is recovered once and cannot disagree between them.

It adds one check the manuscript's ledger does not have. A claim declares the
chapter it belongs to, and the build requires the expected value to **appear
in that chapter's text**. That closes the ledger in both directions:

* a number in the thesis that no committed table produces fails, as in the
  paper;
* and a claim declared against a table that the thesis never actually states
  fails too — which is what stops a ledger being padded to look thorough.

The text check is deliberately generous about *rendering* and strict about
*place*: `1236` matches `1,236` and `0.50` matches `0.5`, because a document
writes numbers the way a reader reads them, but the match must be in the
chapter the claim names.

  python scripts/s25_claims.py --verbose
"""

from __future__ import annotations

import argparse
import re
import sys

import s14_claims
import s25_claims_carried
import s25_claims_thesis
import s25_lib as L

CLAIMS: list[dict] = s25_claims_carried.CLAIMS + s25_claims_thesis.CLAIMS

_TEXT_CACHE: dict[str, str] = {}


def _chapter_text(chapter: int) -> str | None:
    if chapter not in L.FILES_OF_CHAPTER:
        return None
    key = str(chapter)
    if key not in _TEXT_CACHE:
        _TEXT_CACHE[key] = L.chapter_text(chapter)
    return _TEXT_CACHE[key]


#: Prose spells small numbers out, so the appearance check has to read them.
#: Only up to twenty: past that a document writes digits, and admitting more
#: word forms would start matching ordinary English rather than a number.
WORDS = {
    0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
    6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
    12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen",
    16: "sixteen", 17: "seventeen", 18: "eighteen", 19: "nineteen",
    20: "twenty",
}


def _renderings(value: str) -> list[str]:
    """Every way the thesis might legitimately write this value."""
    out = {value}
    m = re.fullmatch(r"-?\d+", value)
    if m:
        n = int(value)
        out.add(f"{n:,}")
        if n in WORDS:
            out.add(WORDS[n])
            out.add(WORDS[n].capitalize())
    m = re.fullmatch(r"(-?\d+)\.(\d+)", value)
    if m:
        whole, frac = m.groups()
        out.add(f"{int(whole):,}.{frac}")
        out.add(value.rstrip("0").rstrip("."))          # 0.50 -> 0.5
        try:
            out.add(f"{float(value):g}")
        except ValueError:
            pass
    return sorted(out, key=len, reverse=True)


def _in_chapter(claim: dict) -> str:
    """'' when the claim's value appears in its chapter, else the error."""
    chapter = claim.get("chapter")
    if chapter is None:
        return "claim declares no chapter"
    text = _chapter_text(chapter)
    if text is None:
        return f"chapter {chapter} is not a chapter of the thesis"
    if not text:
        return f"chapter {chapter} has not been written yet"
    if claim["op"] == "grep":
        return ""            # the needle is a phrase in a report, not a value
    for rendering in _renderings(claim["expect"]):
        if rendering in text:
            return ""
    return (f"chapter {chapter} does not state {claim['expect']!r} — either "
            f"the number is missing from the text or the claim is padding")


def run(verbose: bool = False) -> int:
    _TEXT_CACHE.clear()
    n_fail, rows = s14_claims.check(
        CLAIMS, L.TH / "claims_check.tsv", "s25 claims",
        verbose=verbose, extra=_in_chapter)
    by_chapter: dict[str, int] = {}
    for c in CLAIMS:
        key = str(c.get("chapter"))
        by_chapter[key] = by_chapter.get(key, 0) + 1
    if verbose:
        print("  claims per chapter: " + ", ".join(
            f"{k}:{v}" for k, v in sorted(by_chapter.items(),
                                          key=lambda kv: int(kv[0]))))
    return 1 if n_fail else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    sys.exit(run(verbose=ap.parse_args().verbose))
