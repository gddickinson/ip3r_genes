"""S8 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift. Each entry records **what an earlier task concluded
and where it said it**; the report computes S8's own answer beside it and
renders the verdict from the comparison, printing both numbers either way.

The verdict vocabulary is deliberately four-valued. `underpowered` is not a
polite word for a negative: a Jaccard measured on a genome whose flanking
genes are 73 % unnamed, or a consensus overlap the random-window null itself
reaches, is a measurement that did not happen, and reporting it as a
negative result would be reporting the annotation as biology.
"""

from __future__ import annotations

PRIOR = {
    "sister_pair": {
        "value": "ITPR2 + ITPR3",
        "where": "S7 §5.2 — the AU test rejects ITPR1+ITPR2 (p-AU 1.8e-05) "
                 "and ITPR1+ITPR3 (p-AU 1.65e-05) and does not reject "
                 "ITPR2+ITPR3 (p-AU 0.476); the unconstrained ML tree "
                 "groups ITPR2+ITPR3 at SH-aLRT 100 / UFBoot 100",
    },
    "cyclostome_loci": {
        "value": 6,
        "where": "S7 §5.4 — all 6 cyclostome loci sit in cyclostome-only "
                 "clades, so they are neither one expansion nor three 1:1 "
                 "ohnologs; S7 said explicitly that **which side of the "
                 "vertebrate duplication each lineage attaches to** is "
                 "what S8's synteny would settle",
    },
    "itpr_ryr_separate": {
        "value": 0.0,
        "where": "D14, measured at every stage since S1 — the ryanodine "
                 "receptors carry every ITPR-diagnostic Pfam domain, and "
                 "separating the two families is a positive test, never an "
                 "assumption. S1 measured ITPR-to-RyR covered identity at "
                 "0.249 against 0.828 within the family",
    },
    "paralog_cells": {
        "value": 812,
        "where": "S5 — 812 of the 1,236 genome × class cells are "
                 "`found_annotated`, and every downstream task treats a "
                 "cell's paralog label as orthology. Nothing before S8 "
                 "tested that with evidence outside the gene itself",
    },
    "itpr3_recovery": {
        "value": 232,
        "where": "S5 — ITPR3 is the most consistently recovered paralog "
                 "(232 `found_annotated` cells against 194 for ITPR1 and "
                 "170 for ITPR2), and S5's annotation-quality section "
                 "found the three within 6 points of each other once "
                 "contiguity was held constant",
    },
    "teleost_3r": {
        "value": "two loci per cell",
        "where": "S6 §1 and S7 §5.5 — the teleost 3R co-orthologs are in "
                 "the representative set by rule; S5's ledger files both "
                 "copies into the same paralog cell, so a cell can hold "
                 "two genes that are not the same gene",
    },
}


#: The five verdicts, and what each one means about the measurement.
#:
#:   confirmed        the same property, measured again, agrees
#:   contradicted     the same property, measured again, disagrees
#:   not corroborated a *related* property was measured cleanly and does not
#:                    support the prior — which is not the same as refuting
#:                    it, and saying "contradicted" here would claim a
#:                    disagreement between two things that were never
#:                    measuring the same quantity
#:   orthogonal       the prior stands; S8 measured a different property of
#:                    the same object and the two are not comparable
#:   underpowered     the measurement did not reach the question at all
VERDICTS = ("confirmed", "contradicted", "not corroborated", "orthogonal",
            "underpowered")


def verdict(name: str, holds, note: str = "") -> str:
    """Render one verdict. `holds` is True / False / a VERDICTS string /
    None (underpowered)."""
    if holds is True:
        tag = "confirmed"
    elif holds is False:
        tag = "contradicted"
    elif holds is None:
        tag = "underpowered"
    else:
        tag = str(holds)
        if tag not in VERDICTS:
            raise ValueError(f"{name}: unknown verdict {tag!r}")
    return f"**{tag}** — {note}" if note else f"**{tag}**"
