"""The priors S16's results are judged against, in one module imported by
both halves of the report so a prior cannot be stated twice and drift.

The verdict vocabulary is S8's five-valued one, imported from `s15_priors`
unchanged rather than restated. Two of its values do real work here:

* **`not corroborated`** — a related property measured cleanly that does not
  support the prior. S7's sister pair against the retained flank ohnologs is
  exactly this: a tree estimates the *order* of duplication, while a
  retained ohnolog records which copies survived deletion beside each gene,
  and 2R quartets lose flank copies independently of the duplication order.
* **`orthogonal`** — a different property of the same object. S13 places the
  duplications on the tree; S16 measures what the genome kept beside them.
  A tree node and a retained neighbour are not the same quantity, and
  calling a mismatch between them a contradiction would claim a
  disagreement between two things that never measured the same thing.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s15_priors import VERDICTS                                # noqa: E402,F401

#: key -> (what an earlier task concluded, where it said it)
PRIOR = {
    "paralogon_through_itpr1": (
        "the surviving 2R paralogon links run through ITPR1 — ITPR1 with "
        "ITPR2 and ITPR1 with ITPR3 each retain one shared root flank "
        "family (BHLHE, GRM), ITPR2 with ITPR3 retains none at any bar, and "
        "no ITPR x RyR pair retains one",
        "S8 §4"),
    "sister_pair": (
        "ITPR2 and ITPR3 are sisters — the AU test rejects ITPR1+ITPR2 "
        "(p-AU 1.8e-05) and ITPR1+ITPR3 (p-AU 1.65e-05) and does not reject "
        "ITPR2+ITPR3 (p-AU 0.476); the ML tree groups them at 100/100",
        "S7 §5.2"),
    "duplication_placement": (
        "the two duplications are not on the same branch — ITPR1 splits from "
        "the ITPR2/ITPR3 stem on the vertebrate stem (563 Ma, unbounded "
        "above), and ITPR2 from ITPR3 on the gnathostome stem (462-563 Ma)",
        "S13 §4"),
    "extra_copies_same_paralog": (
        "a cell holding two loci holds two copies of one gene rather than a "
        "misfiled second gene — 36 of 36 extra copies in multi-copy cells "
        "are placed in the same paralog as the cell they were filed under",
        "S8 §5.2"),
    "no_losses": (
        "no vertebrate lineage in this scope has lost an IP3 receptor — "
        "Dollo places 0 losses on the family character and 0 on every "
        "paralog character across 927 genome x paralog cells",
        "S15b §10"),
    "teleost_two_loci": (
        "the teleost 3R co-orthologs are in the representative set by rule "
        "and S5's ledger files both copies into one paralog cell, so a cell "
        "can hold two genes",
        "S6 §1 / S7 §5.5"),
    "ryr_is_the_hazard": (
        "the ryanodine receptors carry every ITPR-diagnostic Pfam domain and "
        "are inside every search this project runs; separating them is a "
        "positive test at every stage",
        "D14"),
}


def prior(key: str) -> tuple[str, str]:
    return PRIOR[key]


def line(key: str, verdict: str, detail: str) -> str:
    """One rendered prior line: what was said, where, what S16 measured.

    The verdict is checked against S8's committed vocabulary rather than
    formatted freely, so a typo cannot invent a sixth verdict.
    """
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not one of {VERDICTS}")
    what, where = PRIOR[key]
    return f"Prior: **{what}** — {where}. **{verdict}** — {detail}"
