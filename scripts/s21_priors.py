"""The priors S21's results are judged against.

One module imported by every half of the report, so a prior cannot be stated
twice and drift. The verdict vocabulary is S8's five-valued one, imported from
`s15_priors` unchanged rather than restated.

Two of the five do real work here. **`orthogonal`** stops a gene-structure
number reading as a disagreement with a sequence one: S6 measured the three
paralogues at 0.83 mean protein identity and S21 measures their intron
positions, and two genes can share every intron and differ in sequence, or the
reverse — the exon structure and the coding sequence are different objects with
the same name in prose. **`underpowered`** is the honest verdict wherever a
comparison lands on a stratum too small to carry it, which here is every
vertebrate class outside the four the sweep samples deeply.

One prior is the *literature's*, and it is the reason S21 exists as a ledger
row: S0 audited "~58-60 exons, hundreds of kilobases" and struck the span half
of it outright — human ITPR3 spans 76 kb — leaving the exon count untested on
anything but human.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s15_priors import VERDICTS                                   # noqa: E402

#: what an earlier task (or the literature audit) concluded, and where.
PRIOR = {
    "exon_count": (
        "S0's audit of the literature's \"~58-60 exons\" left the count "
        "verified on human alone, with the accompanying span claim "
        "(\"hundreds of kb\") struck: human ITPR3 spans 76 kb",
        "roadmap S0 Results; docs/ip3r_background.md"),
    "span_varies": (
        "S0 measured the genomic span varying 6.5x across the three human "
        "paralogues while the protein length varies 3 %",
        "roadmap S0 Results; results/s0_baseline/gene_structure.tsv"),
    "boundary_concordance": (
        "S10 corroborated the sweep's exon boundaries for two loci against "
        "35 and 40 other genomes' independent annotations, at 94.5 % and "
        "100 % of boundaries shared by a majority",
        "results/annotation_bugs/report.md; boundary_concordance.tsv"),
    "annotation_states": (
        "S18 called 264 of the sweep's loci `fragmentary` and 27 `split` "
        "out of 1,874 scored, and found the family is not recorded worse "
        "than its sister",
        "results/annotation_audit/report.md"),
    "d9_source": (
        "S18 ran D9's RefSeq-versus-GenBank contrast twice, raw and above "
        "D4's contiguity bar, and found the gap survives the control",
        "results/annotation_audit/by_source.tsv"),
    "false_negatives": (
        "S19 measured 15.2 % of the sweep's cells as false negatives of the "
        "method and found every one is an assembly, with a 0.9 % residual "
        "above D4's bar over 189 genomes",
        "results/methods/report.md; contiguity_cells.tsv"),
    "bait_identity": (
        "S19 found a single bait at any identity above 0.5 recovers the "
        "locus, so which bait won a locus does not limit the exon "
        "boundaries the aligner places there",
        "roadmap ledger note before S21; results/methods/panel_min_identity.tsv"),
    "copy_number": (
        "S16 counted 174 gene copies over 167 clusters and found its "
        "split-model merge fired on 0 of 2,146 loci, with the 3R teleosts "
        "carrying two copies of a paralogue where the pre-3R ray-finned "
        "outgroup carries one",
        "results/duplication/report.md; loci.tsv, merges.tsv"),
    "paralog_identity": (
        "S6 measured the three paralogues at 0.83 mean within-family "
        "covered identity and ITPR-to-RyR at 0.249",
        "results/msa_v2/report.md"),
    "sister_pair": (
        "S7's AU test placed ITPR2 and ITPR3 as sisters, and S16 found the "
        "retained flanking ohnologs do not corroborate that order",
        "results/phylogeny/report.md; results/duplication/report.md"),
    "d14_genomic": (
        "S5b found that across all 2,144 loci only one family's baits "
        "aligned at all — the ITPR and RyR panels never once contested a "
        "locus",
        "roadmap S5b Results"),
}


def line(key: str, verdict: str, detail: str) -> str:
    """One rendered prior line: what was said, where, what S21 measured.

    The verdict is checked against the committed vocabulary rather than
    formatted freely, so a typo cannot invent a sixth verdict.
    """
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not one of {VERDICTS}")
    what, where = PRIOR[key]
    return f"Prior: **{what}** — {where}. **{verdict}** — {detail}"
