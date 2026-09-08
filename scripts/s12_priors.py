"""S12 — the priors this task's results are judged against.

One module, imported by both halves of the report, so a prior cannot be
stated twice and drift.  Each entry records **what an earlier task
concluded and where it said it**; the report computes S12's own answer
beside it and renders the verdict from the comparison, printing both
numbers either way.

The verdict vocabulary is S8's five-valued one, imported rather than
re-declared.  Two of its values do real work here.

`underpowered` is the honest verdict for a deposit search whose denominator
is 10 records, and S10 already used it in prose: it wrote that a search of
43 records returning nothing "has not shown the gene is untranscribed — it
has shown the species has almost no transcript deposits".  S12 inherits
that sentence as a prior and has to be able to reproduce it as a verdict.

`orthogonal` matters because S12 measures **transcription** and almost
everything before it measured **sequence or annotation**.  A gene with an
intact reading frame that turns out to be silent in twenty libraries has
not contradicted the reading frame; it has answered a different question,
and the report must be able to say so rather than manufacture a conflict.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s8_priors import VERDICTS, verdict  # noqa: F401,E402  (re-exported)

PRIOR = {
    "handoff": {
        "value": "unanswerable from the deposits",
        "where": "S10 §3 and §4, both cases, verbatim: \"0 hits, 0 of them "
                 "contiguous across a junction. The denominator is the "
                 "result: this species has 43 transcript records in total "
                 "[10 for *D. eleginoides*]. … RNA-seq for it does exist "
                 "and reaching it needs the streaming aligner S12 builds; "
                 "these junctions are handed there rather than "
                 "half-answered here.\" This is the prior S12 was created "
                 "to settle, and it is the only one in the list that names "
                 "this task",
    },
    "annotation_failures": {
        "value": 7,
        "where": "S10 — of 382 loci in assemblies whose annotation could "
                 "reasonably have delivered the gene, 375 (98.2 %) are "
                 "delivered whole and 7 are not: 3 omissions, 2 "
                 "truncations, 2 fragmentations, across 3 genomes "
                 "(results/annotation_bugs/report.md §1)",
    },
    "not_pseudogenes": {
        "value": 0,
        "where": "S10 §3 — every family locus in the case genomes splices "
                 "into an uninterrupted reading frame: *Nibea albiflora* "
                 "ITPR2 is 2,673 codons with **0 internal stops** against "
                 "14.2 expected under neutral drift to its own divergence, "
                 "and *D. eleginoides* ITPR3 is 2,687 codons with 0 "
                 "against 11.4. An intact ORF is a necessary condition for "
                 "a transcribed gene and not a sufficient one, so S12's "
                 "answer can corroborate it but a silent locus would not "
                 "contradict it",
    },
    "boundaries_corroborated": {
        "value": 0.945,
        "where": "S10 §3 — 52 of 55 of *Nibea albiflora* ITPR2's internal "
                 "exon boundaries (94.5 %), and 59 of 59 of "
                 "*D. eleginoides* ITPR3's (100 %), are placed identically "
                 "by a majority of the 35 and 40 swept genomes whose own "
                 "annotation independently models the same paralog from "
                 "the same bait. That validates the boundaries as an "
                 "instrument; reads crossing them validate these loci",
    },
    "flank_position": {
        "value": 0.916,
        "where": "S10 §3 reading S8 — *Nibea albiflora*'s unannotated "
                 "ITPR2 lies between `SSPN` and `BHLHE41`, and `SSPN` is "
                 "this paralog's most conserved neighbour, found beside "
                 "ITPR2 in 197 of 215 swept vertebrates (91.6 %)",
    },
    "cell_status": {
        "value": "found_unannotated",
        "where": "S5b — the ledger calls all three ITPR cells of both "
                 "*Dissostichus* genomes and *Nibea albiflora*'s ITPR2 "
                 "`found_unannotated`, and S5 deliberately declined to "
                 "adjudicate: 'one instrument does not overturn a public "
                 "annotation'",
    },
    "decoy_floor": {
        "value": 0,
        "where": "the PIEZO project's S12, whose reversed-CDS decoys "
                 "collected 0 reads across all 78 runs of its four-species "
                 "panel. Method ported, result not: the decoy floor is "
                 "re-measured here on different species, different "
                 "libraries and much larger references",
    },
}


def where(key: str) -> str:
    return PRIOR[key]["where"]


def value(key: str):
    return PRIOR[key]["value"]
