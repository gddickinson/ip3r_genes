"""S15 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift.  Each entry records **what an earlier task
concluded and where it said it**; the report computes S15a's own answer
beside it and renders the verdict from the comparison, printing both
numbers either way.

The verdict vocabulary is S8's five-valued one, reused unchanged.  Two of
its values do real work here.  `underpowered` is the honest answer for the
brief's own first step: a synteny comparison run where 273 of 432 regions
sit on a contig carrying no annotated gene at all is a measurement that
did not happen, and reporting it as "synteny does not corroborate the
attribution" would report an assembly's contig lengths as biology.  And
`orthogonal` is what stops the ORF screen from being read as a
disagreement with S9: lesion density is a count of indels in an assembly
and omega is a ratio of substitution rates, so the two can point the same
way without either corroborating the other.
"""

from __future__ import annotations

VERDICTS = ("confirmed", "contradicted", "not corroborated", "orthogonal",
            "underpowered")

PRIOR = {
    "corroborated_losses": {
        "value": 0,
        "where": "S13 §4 — the reconciliation implies 51-53 losses across "
                 "the vertebrate subtree and **0** survive being asked of "
                 "the S5 genome ledger: 26 of the 46 candidate cells are "
                 "`sampling_artefact` (the paralog is in the genome but "
                 "not in S6's representative set) and 4 are "
                 "`paralog_unassignable`. S13 handed S15 the explicit "
                 "instruction to build the matrix from the sweep and not "
                 "from tips of the tree",
    },
    "absent_cells": {
        "value": 4,
        "where": "S5b — 4 `absent` cells in 1,236, all four cyclostome "
                 "(ITPR2 and ITPR3 in *Petromyzon marinus* and *Myxine "
                 "glutinosa*), beside 43 `tblastn_trace` and 7 `fragment`",
    },
    "cyclostome_unassignable": {
        "value": 2,
        "where": "S13 D45 — both cyclostome genomes carry **three** ITPR "
                 "loci apiece, all filed in the ITPR1 cell because the S5 "
                 "bait panel has no cyclostome-labelled bait, so their "
                 "`absent` cells are a panel limit and not an absence",
    },
    "synteny_caller": {
        "value": 1.0,
        "where": "S8 §4 — the consensus paralog caller is 405/405 correct "
                 "on the loci it calls (leave-one-genome-out, on loci "
                 "whose paralog their own annotation establishes) at a "
                 "0.83 % false-call rate on matched random windows. S13 "
                 "and the S15 brief both name it as the instrument that "
                 "would disambiguate the traces",
    },
    "aves_contiguity": {
        "value": 0.66,
        "where": "S5b — 66 % of Aves assemblies fail D4's contiguity bar "
                 "against 11 % of Actinopteri, and 68 % of margin species "
                 "against 12 % of order representatives",
    },
    "annotation_right": {
        "value": 0.982,
        "where": "S10 — the annotation gets the gene right at 375 of 382 "
                 "loci (98.2 %) where it demonstrably could, and the 7 "
                 "failures are in 3 genomes",
    },
    "itpr1_constraint": {
        "value": "ITPR1 held roughly twice as tightly",
        "where": "S9b — three independent framings (one-ratio omega per "
                 "paralog, the two-ratio contrast, and RELAX) agree that "
                 "ITPR1 is under roughly twice the purifying constraint of "
                 "ITPR2 and ITPR3",
    },
    "ryr_control": {
        "value": 309,
        "where": "S5b — the RyR positive control fired in every one of the "
                 "309 genomes, so no genome is excluded on control grounds "
                 "and every empty ITPR cell is an empty cell in an "
                 "assembly the search demonstrably reached",
    },
}


def cite(key: str) -> str:
    p = PRIOR[key]
    return f"{p['where']}"


def verdict(prior, observed, *, tolerance: float = 0.0,
            orthogonal: bool = False, underpowered: str = "") -> str:
    """The comparison, rendered by the same rule everywhere in the report."""
    if underpowered:
        return "underpowered"
    if orthogonal:
        return "orthogonal"
    if isinstance(prior, (int, float)) and isinstance(observed, (int, float)):
        if abs(float(prior) - float(observed)) <= tolerance:
            return "confirmed"
        return "contradicted"
    return "confirmed" if str(prior) == str(observed) else "contradicted"
