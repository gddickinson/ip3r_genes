"""S9 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift. Each entry records **what an earlier task concluded
and where it said it**; the report computes S9's own answer beside it and
renders the verdict from the comparison, printing both numbers either way.

The verdict vocabulary is S8's five-valued one, imported rather than
re-declared, and `orthogonal` earns its keep here more than anywhere. Two
of the priors S9 meets are about *the genome* — how well a paralog's
neighbourhood travels, how often its locus is recovered — and ω is about
the coding sequence. Calling a mismatch there "contradicted" would claim a
disagreement between two things that were never measuring the same
quantity.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s8_priors import VERDICTS, verdict  # noqa: F401,E402  (re-exported)

PRIOR = {
    "purifying": {
        "value": 0.910,
        "where": "S6 §4.1 — mean covered-only identity *within* a "
                 "vertebrate paralog group is 0.910, and the background "
                 "(`docs/ip3r_background.md` §1) describes a 2,700-residue "
                 "channel whose every domain is load-bearing. Nothing "
                 "before S9 measured selection: identity is a distance, "
                 "not a rate, and a conserved protein and a slowly "
                 "evolving one are not the same statement",
    },
    "paralog_divergence": {
        "value": 0.762,
        "where": "S6 §4.2 — mean covered-only identity *between* paralog "
                 "groups is 0.791 (ITPR1×ITPR2), 0.753 (ITPR2×ITPR3) and "
                 "0.741 (ITPR1×ITPR3), against 0.910 within. The three "
                 "paralogs are 2R products older than 500 Myr, so S9 "
                 "expects synonymous sites between them to be saturated "
                 "and every cross-paralog ω to be qualified by that",
    },
    "sister_pair": {
        "value": "ITPR2 + ITPR3",
        "where": "S7 §5.2 — the AU test rejects ITPR1+ITPR2 (p-AU 1.8e-05) "
                 "and ITPR1+ITPR3 (1.65e-05) and does not reject "
                 "ITPR2+ITPR3 (0.476). S9's branch models mark each "
                 "paralog's stem on that topology, so the tree is an "
                 "input here and not a result (D15's discipline applied "
                 "to the gene tree)",
    },
    "itpr3_neighbourhood": {
        "value": 0.062,
        "where": "S8 §4.4 — ITPR3's flanking neighbourhood is the one that "
                 "does not travel: cross-class Jaccard 0.062 against "
                 "0.160 (ITPR1) and 0.158 (ITPR2), a 2.5× gap, and in "
                 "human that ground is the MHC at 6p21. Whether the "
                 "*coding sequence* shows the same asymmetry is a "
                 "different measurement on the same object",
    },
    "itpr1_clinical": {
        "value": "ITPR1",
        "where": "`docs/ip3r_background.md` §5 and the review §9 — the "
                 "family's characterised disease alleles are "
                 "overwhelmingly ITPR1 (SCA15/SCA16 deletions, SCA29 and "
                 "Gillespie syndrome missense), with ITPR2 anhidrosis and "
                 "ITPR3 neuropathy each resting on far fewer families. A "
                 "dominant missense burden is what a gene under tight "
                 "constraint looks like clinically, so ITPR1 is the "
                 "paralog S9 would expect to carry the lowest ω",
    },
    "itpr3_recovery": {
        "value": 232,
        "where": "S5 — ITPR3 is the most consistently recovered paralog "
                 "(232 `found_annotated` cells against 194 ITPR1 and 170 "
                 "ITPR2). Recovery is a property of the assembly and the "
                 "bait panel, not of the gene's evolutionary rate",
    },
}
