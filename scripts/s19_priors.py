"""The priors S19's results are judged against.

One module imported by every half of the report, so a prior cannot be stated
twice and drift. The verdict vocabulary is S8's five-valued one, imported
from `s15_priors` unchanged rather than restated.

Two of those values do real work here. **`orthogonal`** stops a methods
number being read as a disagreement with the biology it was measured
alongside: S15b's zero losses and S19's 15 % cell-level false-negative rate
are not in conflict, they are the two halves of one design — the losses are
zero *because* the misses were chased, and the misses are what a sweep that
reported only its ledger would have called absences. **`underpowered`** is
the honest verdict wherever a channel's contribution is measured on a
denominator too small to carry it, which here is every non-vertebrate
database except the protists.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s15_priors import VERDICTS                                # noqa: E402

#: what an earlier task concluded, and where it said it.
PRIOR = {
    "false_negatives": (
        "S5b measured recovery at 98-99 % above D4's contiguity bar and "
        "57-70 % below it, over a scope in which S15b later reconstructed "
        "no losses at all",
        "roadmap S5b Results; results/loss_counts/report.md"),
    "contiguity_bar": (
        "D4's bar is the median measured ITPR genomic span, 142,212 bp of "
        "contig N50, chosen a priori from gene geometry and never "
        "calibrated against a measured error rate",
        "roadmap D4; scripts/s5_calibration.py:itpr_span_stats"),
    "span_bias": (
        "S5b confirmed a 13-point recovery gap below the bar in the "
        "direction of gene span — ITPR3 70 %, ITPR1 61 %, ITPR2 57 %",
        "roadmap S5b Results"),
    "hmm_over_pfam": (
        "S3 found the profile sweep adds 618 proteins the InterPro census "
        "never returned, all 209-942 aa fragmentary gene models, and that "
        "0 of 2,787 v2 ITPR records in a swept proteome were missed",
        "results/census_v3/report.md"),
    "genome_only": (
        "S5b's census v4 adds 1,058 gene models from 224 genomes, of which "
        "318 ITPR models exist only as DNA and 167 more sit inside an "
        "annotated gene carrying no family name",
        "roadmap S5b Results; results/census_v4/genome_models.tsv"),
    "d10b": (
        "S3 and S20b found off-family accretion *dilutes* the sister-family "
        "share, so K1 moves the wrong way while a run drifts and only the "
        "round ceiling catches it",
        "roadmap D10b; results/s20_sweep/report.md"),
    "bait_breadth": (
        "S5a found the unlabelled `vertebrate_basal` baits were competing "
        "as a fourth paralog until rescue attribution was fixed, and S5's "
        "panel leaves six slots unfilled for want of a labelled record",
        "roadmap S5a Results; results/s5_baits/unfilled_slots.tsv"),
    "d14_genomic": (
        "S5b found that across all 2,144 loci only one family's baits "
        "aligned at all — the ITPR and RyR panels never once contested a "
        "locus",
        "roadmap S5b Results"),
    "recon_losses": (
        "S13 audited its own reconciliation against the genome ledger and "
        "found the implied losses are overwhelmingly sampling, with four "
        "corroborated absences before D45 removed them",
        "results/reconciliation/report.md"),
    "synteny_reach": (
        "S15a found S8's caller is accurate where it acts and almost never "
        "reaches — 273 of 432 trace regions sit on a contig carrying no "
        "annotated gene at all, and 8 reach the four-key floor",
        "results/loss_dynamics/report.md"),
    "fel_power": (
        "S17 reported per-site selection on the same coordinates as its "
        "constraint layers, taking omega only over sites where the "
        "synonymous rate is identifiable",
        "results/constraint/report.md"),
}


def line(key: str, verdict: str, detail: str) -> str:
    """One rendered prior line: what was said, where, what S19 measured.

    The verdict is checked against the committed vocabulary rather than
    formatted freely, so a typo cannot invent a sixth verdict.
    """
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not one of {VERDICTS}")
    what, where = PRIOR[key]
    return f"Prior: **{what}** — {where}. **{verdict}** — {detail}"
