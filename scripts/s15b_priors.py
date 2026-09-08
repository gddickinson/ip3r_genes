"""S15b — the priors this half's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift.  Each entry records **what an earlier task
concluded and where it said it**; the report computes S15b's answer beside
it and renders the verdict from the comparison, printing both numbers
either way.

The verdict vocabulary is S8's five-valued one, imported from
`s15_priors` unchanged rather than restated — S15a and S15b are two halves
of one task and a second copy of the vocabulary is a second thing to drift.

Two values earn their keep here.  `orthogonal` is what stops the
duplication-placement echo being read as corroboration of S13: Dollo's
gain node is the MRCA of the tips that *have* the character in this
matrix, which is a statement about sampling and about which tips the bait
panel could label, while S13's placement is a reconciliation of a gene
tree against a dated species tree.  They can agree without either being
evidence for the other, and here they do agree.  And `underpowered` is the
honest verdict on the bird ITPR3 result once its contiguity control is
read: 21 of its 27 pairs are in assemblies below D4's bar, and the 6 above
it all point the same way but cannot carry a test.
"""

from __future__ import annotations

from s15_priors import VERDICTS, verdict          # noqa: F401  (re-exported)

PRIOR = {
    "absent_cells": {
        "value": 0,
        "where": "S15a — 0 of 927 genome x paralog cells reaches the "
                 "`absent` state, the only state S15b may count, and that "
                 "state is reachable: `s15_test_loss.py` T8 constructs a "
                 "contiguous, controlled, empty, spare-free cell and "
                 "requires it to come back `absent`",
    },
    "min_gene_equivalents": {
        "value": 3.0,
        "where": "S15a — in all 189 assemblies above D4's contiguity bar "
                 "the minimum implied copy number is 3.00 gene-equivalents",
    },
    "loss_candidates": {
        "value": 4,
        "where": "S15a `loss_candidates.tsv` — four near-misses, all of "
                 "them cyclostome ITPR2/ITPR3 cells stopped by R5 (D45), "
                 "each row naming the rule that stopped it",
    },
    "family_coding_primary": {
        "value": "family-level presence per genome",
        "where": "S15a D46 — the `co_trace` population (20 regions, median "
                 "reassembly 0.631, sitting between the candidates and the "
                 "decoy) is the measured size of the paralog-attribution "
                 "problem in a shattered assembly, so S15b's primary "
                 "coding must be family-level and the paralog-resolved "
                 "matrix is the sensitivity axis",
    },
    "itpr1_gain_node": {
        "value": "Vertebrata",
        "where": "S13 — the split separating ITPR1 from ITPR2+ITPR3 is on "
                 "the **vertebrate stem**, older than crown Vertebrata, "
                 "and the split separating ITPR2 from ITPR3 is on the "
                 "**gnathostome stem**",
    },
    "itpr23_gain_node": {
        "value": "Gnathostomata",
        "where": "S13 — ITPR2 and ITPR3 separate on the gnathostome stem, "
                 "so no cyclostome carries either as such; S7 says the "
                 "same from the other side, with every cyclostome tip "
                 "outside all three paralog clades",
    },
    "itpr3_indel_excess": {
        "value": 0.0032,
        "where": "S15a D47 — ITPR3 carries an indel excess against its own "
                 "genome's identity-matched sibling loci, 39 genomes to "
                 "14, q = 0.0032, and nothing in this project explains it. "
                 "ITPR2's apparent excess disappears under identity "
                 "matching (q = 0.902) and the RyR control shows none "
                 "(q = 0.090)",
    },
    "no_fossils": {
        "value": 0,
        "where": "S15a §11.4 — every locus above the measured lesion bar "
                 "is at full coverage with an intact model, so the "
                 "shared-lesion Poisson test has no dead loci to run on",
    },
    "zero_stops_falsifies": {
        "value": 0,
        "where": "S10 — all 8 family loci across the two case genomes "
                 "carry 0 internal stops against 4-25 expected under "
                 "neutral drift at the observed divergence; the rule is "
                 "one-sided, zero stops falsifies a pseudogene call and a "
                 "handful does not establish one",
    },
    "aves_contiguity": {
        "value": 0.66,
        "where": "S5b — 66 % of Aves assemblies fail D4's contiguity bar "
                 "against 11 % of Actinopteri",
    },
    "polytomy_degree": {
        "value": 23,
        "where": "S15a — the 309-genome NCBI taxonomy tree carries 468 "
                 "internal nodes and 49 polytomies, the largest of degree "
                 "23 (Passeriformes)",
    },
}


def cite(key: str) -> str:
    return PRIOR[key]["where"]
