"""The priors S13's results are judged against, in one module.

Imported by both halves of the report so a prior cannot be stated twice and
drift. Each entry records what an earlier task concluded **and where it said
it**; the report computes this task's answer beside it and renders the verdict
from the comparison.

The verdict vocabulary is S8's five-valued one, unchanged. Two of its values
do real work here. `underpowered` is the honest verdict for a question this
task's *instrument* cannot reach — a reconciliation on a representative
alignment cannot count losses, and saying so is not a hedge. `orthogonal` is
for a prior measuring a different property of the same object: S8's retained
flanking ohnologs are a **deletion** record and a duplication placement is a
**branching** statement, so a mismatch between them is not a disagreement.
"""

from __future__ import annotations

VERDICTS = ("confirmed", "contradicted", "not corroborated", "orthogonal",
            "underpowered")

PRIOR = {
    "sister": dict(
        where="S7 §5.2",
        text="**ITPR2 + ITPR3 are sisters, ITPR1 outside** — the AU test "
             "rejects ITPR1+ITPR2 (p-AU 1.8e-05) and ITPR1+ITPR3 "
             "(p-AU 1.65e-05) and does not reject ITPR2+ITPR3 (p-AU 0.476); "
             "the unconstrained ML tree groups ITPR2+ITPR3 at 100/100."),
    "cyclostome": dict(
        where="S7 §5.4 and S8 §6",
        text="**the 6 cyclostome loci sit in cyclostome-only clades** — "
             "neither one lineage-specific expansion nor three 1:1 ohnologs. "
             "S7 asked which side of the vertebrate duplication each lineage "
             "attaches to; S8 measured it with flank synteny, returned "
             "**underpowered**, and named the instrument that could answer: "
             "\"the 2R paralogon reconstructed from a cyclostome-anchored "
             "gene tree\"."),
    "paralogon": dict(
        where="S8 §4",
        text="**the surviving paralogon links run through ITPR1** — ITPR1 "
             "shares one root-level flanking family with ITPR2 and one with "
             "ITPR3, while ITPR2 with ITPR3 retains none at any prevalence "
             "bar."),
    "ohnologs": dict(
        where="the review, §7.4",
        text="**not demonstrated** — whether ITPR1/2/3 are ohnologues from "
             "the two rounds of vertebrate whole-genome duplication has not "
             "been shown with synteny-backed, phylogeny-tested evidence, and "
             "no published support-annotated ML analysis with an RyR outgroup "
             "fixes the rooted topology."),
    "clock": dict(
        where="S9 §4.2",
        text="**no clock is available** — median pairwise dS within a single "
             "paralog set runs 4.6-13.5 and 84.6-94.2 % of within-paralog "
             "pairs exceed the dS = 1.5 saturation bar, so synonymous sites "
             "cannot date anything at this depth, let alone between paralogs."),
    "losses": dict(
        where="S5b",
        text="**the family is almost never lost in the vertebrate scope** — "
             "across 309 genomes the sweep records `absent` in 4 of 1,236 "
             "cells (ITPR2 and ITPR3 in the two cyclostomes) with 43 further "
             "cells holding a tblastn remnant and 7 a fragment."),
    "support": dict(
        where="S7 §5.6",
        text="**the deep arrangement is the weak point** — 69.5 % of internal "
             "nodes clear both support thresholds, but the node separating "
             "ITPR1 from ITPR2+ITPR3 within the vertebrate clade sits at "
             "SH-aLRT 17.4 / UFBoot 54."),
}


def verdict_line(key: str, computed: str, verdict: str) -> str:
    """`Prior: ... **verdict** — computed`, the shape every section uses."""
    p = PRIOR[key]
    assert verdict in VERDICTS, verdict
    return (f"Prior: {p['text']} ({p['where']})\n\n"
            f"**{verdict}** — {computed}")
