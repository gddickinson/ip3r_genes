"""S18 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be stated
twice and drift. Each entry records what an earlier task concluded **and where
it said it**; the report computes S18's own answer beside it and renders the
verdict from the comparison, printing both numbers either way.

The verdict vocabulary is S8's five-valued one, imported from `s15_priors`
unchanged rather than restated. Two of its values do real work here.

`orthogonal` is what stops S10's 98.2 % from reading as a contradiction of
S18's 72.4 %. S10 measured the annotation over the loci where it *demonstrably
could* have delivered the gene — five eligibility rules, of which E5 asks
whether the annotation builds genes this long anywhere else in the same genome
— and S18 measures it over every locus in the scope. The two are different
denominators of the same numerator and the difference between them is itself a
result, not a disagreement.

`underpowered` is the honest verdict on the family-versus-sister comparison
once its numbers are in: above D4's contiguity bar the ITPR cells and the RyR
control differ by 1.6 points on 568 and 632 loci, which is a measurement that
cannot resolve a difference of the size it would take to matter.
"""

from __future__ import annotations

from s15_priors import VERDICTS                                   # noqa: F401

PRIOR = {
    "s10_annotation_correct": {
        "value": 0.982,
        "where": "S10 §3 — the annotation gets the gene right at 375 of 382 "
                 "loci (98.2 %) **where it demonstrably could**: 880 "
                 "recovered loci, 382 passing E1-E5, median annotation loss "
                 "0, 360 at exactly 0",
    },
    "s5b_dna_only_models": {
        "value": 318,
        "where": "S5b — census v4 adds 1,058 gene models from 224 genomes, of "
                 "which **318 ITPR models exist only as DNA** and 167 more "
                 "sit inside an annotated gene carrying no family name, "
                 "unreachable by any name-based search",
    },
    "s5b_paralog_naming_gap": {
        "value": 0.23,
        "where": "S5b — holding contiguity constant across 487 loci, ITPR3 is "
                 "88 % correctly named against ITPR1's 65 %, a 23-point gap "
                 "between genes of near-identical protein length in the same "
                 "genomes. Uncontrolled the gap looked like 54 % vs 14 %, "
                 "which was assembly quality",
    },
    "s5b_paralog_conflicts": {
        "value": 3,
        "where": "S5b — 3 loci are claimed by a cell other than the one their "
                 "annotation names (1 with a coherent sibling locus), "
                 "**supplied to S18, not adjudicated**",
    },
    "s3_zero_hit_proteomes": {
        "value": 15,
        "where": "S3 — 15 of 764 vertebrate reference proteomes returned no "
                 "family hit from either profile. S3 could not say whether "
                 "that was a gene or a gene caller",
    },
    "s2_naming_pfam_recall": {
        "value": 0.811,
        "where": "S2 — 2,911 of 15,417 seeded-space proteins do not carry "
                 "PF08709, the signature that names the family, so a "
                 "PF08709 query recovers 81.1 % of the space its own sister "
                 "signatures enumerate",
    },
    "s0_zebrafish_ryr_share": {
        "value": 0.49,
        "where": "S0 — 49 % (53/109) of zebrafish PF08709 records are "
                 "ryanodine receptors. D14's hazard, measured at the "
                 "signature level, is what makes wrong-family naming the "
                 "audit's sharpest question",
    },
    "s5b_contiguity_bar": {
        "value": 120,
        "where": "S5b/D4 — 120 of 309 manifest genomes cannot hold a median "
                 "ITPR gene on one contig (66 % of birds against 11 % of "
                 "fish), so any annotation failure rate has to be reported "
                 "with that control",
    },
}


def line(key: str, observed: str, verdict: str) -> str:
    """One prior, its source, this task's answer and the verdict between."""
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not in the committed vocabulary "
                         f"{VERDICTS}")
    p = PRIOR[key]
    return (f"- **{verdict}.** {p['where']}. S18 measures **{observed}**.")
