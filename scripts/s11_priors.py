"""S11 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift. Each entry records **what an earlier task (or the
brief) concluded and where it said it**; the report computes S11's own
answer beside it and renders the verdict from the comparison, printing both
numbers either way.

The verdict vocabulary is S8's five-valued one, imported rather than
re-declared. `orthogonal` matters here for one reason in particular: most
of what S11 measures is a *fold* and most of what came before it measured a
*sequence*. Two proteins 25 % identical that share a fold are not a
contradiction of the 25 %; they are a different quantity, and calling that
disagreement would be a category error dressed as a result.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s8_priors import VERDICTS, verdict  # noqa: F401,E402  (re-exported)

PRIOR = {
    "afdb_poor": {
        "value": "poor",
        "where": "the S11 brief — \"AFDB coverage of the family: … Expect "
                 "it to be poor for ~2,700-residue proteins; that absence "
                 "is itself a result.\" This is the one prior in the list "
                 "that is a *guess written in advance*, which is exactly "
                 "why it is worth writing down before measuring",
    },
    "d14_separable": {
        "value": 1.0,
        "where": "D14, and every stage since S2 — the ITPR/RyR separation "
                 "has been a positive test at architecture (S2), profile "
                 "(S3), alignment score (S5), tree (S7) and neighbourhood "
                 "(S8) level, agreeing with the census call every time. "
                 "S11 asks the same question of a sixth instrument that "
                 "reads none of the same evidence: shape",
    },
    "sequence_separation": {
        "value": 0.249,
        "where": "S1 — ITPR-to-RyR covered-only sequence identity is "
                 "0.249, against 0.828 within the family; S6 confirmed "
                 "both on the representative alignment. Structure is a "
                 "different quantity from identity and the two families "
                 "are known to share the channel fold, so what S11 can "
                 "add is the *size* of the structural gap, not a repeat "
                 "of the sequence one",
    },
    "nonvertebrate_real": {
        "value": "real receptors",
        "where": "S20 and S23 — the family's range reaches Viridiplantae, "
                 "Fungi, SAR, Discoba and Amoebozoa, and S20's per-record "
                 "verdicts chased each plant and fungal record for "
                 "contamination rather than accepting it. Those records "
                 "are the ones a fold test is worth most on: they are "
                 "held on sequence evidence alone",
    },
    "state_spread": {
        "value": "gating, not fold",
        "where": "S0's structural measurement and the review §3 — the "
                 "family's cryo-EM series resolves apo, resting, "
                 "preactivated, activated and inhibited states of the "
                 "*same* protein. If those states moved TM-score as much "
                 "as paralog identity does, no fold claim in this task "
                 "would be readable, so the state panel is a control on "
                 "S11's own instrument",
    },
    "ip3_core_diagnostic": {
        "value": "PF08709",
        "where": "the review §2 and S22's scope — PF08709 is annotated on "
                 "*both* families (its Pfam name is \"Inositol "
                 "1,4,5-trisphosphate/ryanodine receptor\", and S0 "
                 "measured it on RYR1/2/3 as well as on ITPR1/2/3), but "
                 "only the IP3 receptors bind IP$_3$ there — the "
                 "β-trefoil is shared as a *domain* and unshared as a "
                 "*function*. S22 is scoped around it, so whether "
                 "AlphaFold models it well decides whether S17 and S22 "
                 "can stand on predicted structures at all",
    },
}


def where(name: str) -> str:
    return PRIOR[name]["where"]


def value(name: str):
    return PRIOR[name]["value"]
