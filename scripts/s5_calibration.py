"""s5_calibration.py — read the committed calibrations, never retype them.

Three thresholds in the S5 sweep are measurements rather than settings, and
each one changes a call the ledger makes:

  `-G`          too small splits a gene, and a split ITPR reads out of the
                ledger as `fragment`;
  the cap on it a tractability choice, not biology, so the genomes it binds
                on have to be flagged;
  D4's bar      an assembly whose contigs are shorter than the gene cannot
                carry it, so its empty cells are evidence about the assembly.

All three are read back out of `results/s5_baits/intron_calibration.{json,tsv}`
rather than being written into the code, so a threshold cannot drift from the
measurement that justified it — the discipline D13 applies to reports, applied
to constants.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
S5_BAITS = PROJECT_ROOT / "results" / "s5_baits"

DEFAULT_MAX_INTRON = 200_000     # miniprot's own default (-G)


_INTRON_RULE: dict | None = None


def intron_rule() -> dict:
    """The -G rule, read from the committed calibration rather than retyped.

    `s5_intron_calibration.py` measured it (D13 applied to a threshold): the
    widest intron in either family across a vertebrate panel, doubled, scaled
    by genome size from the genome it was measured in, floored at miniprot's
    own default and capped for tractability.
    """
    global _INTRON_RULE
    if _INTRON_RULE is None:
        path = S5_BAITS / "intron_calibration.json"
        if not path.exists():
            raise SystemExit(f"intron calibration missing: {path}\n"
                             "  run: python scripts/s5_intron_calibration.py")
        _INTRON_RULE = json.loads(path.read_text())["max_intron_rule"]
    return _INTRON_RULE


def max_intron_for(genome_bp: int) -> int:
    """miniprot -G for a genome of this size, from the measured rule.

    A too-small -G does not lose a gene, it *splits* one — and a split ITPR
    reads out of the ledger as `fragment`, which is the deliverable's own
    failure mode. A too-large one joins neighbouring loci. Both directions
    are why this comes off a measurement.
    """
    rule = intron_rule()
    if genome_bp <= 0:
        return int(rule["floor_bp"])
    scaled = int(rule["scale_bp"] * genome_bp / rule["anchor_genome_bp"])
    return min(int(rule["cap_bp"]), max(int(rule["floor_bp"]), scaled))


_SPAN_STATS: dict | None = None


def itpr_span_stats() -> dict:
    """Measured ITPR genomic spans, for D4's contiguity bar.

    D4 says an absence claim needs a genome-wide contiguity floor as well as
    a local check, and leaves the floor to S19. This is not that floor — it
    is the weaker statement that can be made now without inventing a number:
    an assembly whose contig N50 is below the span of the gene cannot carry
    that gene on one contig, so `fragment` / `assembly_gap` / `no_locus`
    there says nothing about whether the gene is present. The spans come
    from this project's own `intron_calibration.tsv`, so the bar is measured
    rather than asserted.
    """
    global _SPAN_STATS
    if _SPAN_STATS is None:
        path = S5_BAITS / "intron_calibration.tsv"
        if not path.exists():
            raise SystemExit(f"intron calibration missing: {path}\n"
                             "  run: python scripts/s5_intron_calibration.py")
        spans = []
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["family"] == "ITPR" and row["span_bp"]:
                    spans.append(int(row["span_bp"]))
        spans.sort()
        _SPAN_STATS = {"n": len(spans), "min": spans[0], "max": spans[-1],
                       "median": spans[len(spans) // 2]}
    return _SPAN_STATS


def spans_a_gene(contig_n50: int) -> bool:
    """Can a contig of this N50 hold a whole ITPR gene? (D4's first bar.)"""
    return bool(contig_n50) and contig_n50 >= itpr_span_stats()["median"]


def intron_rule_capped(genome_bp: int) -> bool:
    """True when the cap — a tractability choice, not biology — is binding.

    The ledger flags these genomes: a `fragment` there could be a real intron
    the sweep declined to span.
    """
    rule = intron_rule()
    return int(rule["scale_bp"] * max(0, genome_bp)
               / rule["anchor_genome_bp"]) > int(rule["cap_bp"])
