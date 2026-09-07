"""S7 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift. Each entry records **what an earlier task
concluded and where it said it**, and the report computes the tree's
own answer beside it and renders `confirmed` / `contradicted` /
`unresolved` from the comparison — printing both numbers either way. A
generator written to narrate the expected answers would print them
whatever the tree said, which is how a pipeline launders an assumption
into a result.
"""

from __future__ import annotations

PRIOR = {
    "sister_pair": {
        "value": "ITPR1 + ITPR2",
        "margin": 0.042,
        "where": "S6 §4.2 — mean covered-only identity over the trimmed "
                 "MSA ranks ITPR1×ITPR2 at 0.788 against 0.746 and 0.736, "
                 "a lead of 0.042 with non-overlapping interquartile "
                 "ranges. S6 reported it as a preview and said the AU "
                 "test was the answer",
    },
    "published_sister": {
        "value": None,
        "where": "the review §7.4 — which two of the three vertebrate "
                 "paralogues are sisters is **not fixed by any published, "
                 "support-annotated ML analysis with an RyR outgroup**",
    },
    "paralog_monophyly": {
        "value": 3,
        "where": "assumed by every stage since S2: the architecture call "
                 "assigns a paralog per record, S5's ledger has one cell "
                 "per paralog per genome, and S6 built the representative "
                 "grid on clade band × paralog",
    },
    "cyclostome_nearest": {
        "value": "ITPR1",
        "where": "S6 §4.3 — all 6 cyclostome loci fall nearest ITPR1 on "
                 "covered identity, which is equally what a "
                 "lineage-specific expansion and three fast-evolving 1:1 "
                 "orthologs look like. S6 said the tree distinguishes them",
    },
    "duplicate_pairs": {
        "value": "sisters",
        "where": "S6 §1 — the teleost 3R co-orthologs are in the "
                 "representative set by rule and are a stress test: if "
                 "the tree does not recover a species' two same-paralog "
                 "copies as sisters, the naming is wrong or the alignment "
                 "is",
    },
}
