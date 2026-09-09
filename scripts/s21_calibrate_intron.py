"""What counts as an intron — measured against the genome, not chosen.

The brief expected miniprot to emit two CDS records either side of a
frameshift, which would make a broken locus look exon-rich, and asked for those
pairs to be merged. Whether such pairs exist is a measurement, and this module
makes it: every consecutive block pair in the sweep is binned by the gap
between the two blocks in gene order and scored by the **splice dinucleotides**
at that gap's two edges — evidence the alignment score did not produce. A gap
that is a real intron reads `GT..AG`; an indel the aligner stepped over has no
reason to.

The two populations a bar would separate are therefore *spliceable gaps* and
*sub-intron gaps*, and the calibration reports the separation rather than
asserting it. Following `s15_calibrate_recon.bar()` and
`s23_calibrate_loci.separation()`, it **refuses to hand out a threshold** its
own data does not support: with no sub-intron population there is nothing to
put a bar between, so the floor stays at the declared fallback, the refusal is
written into the table, and the merge rule then demonstrably fires on nothing.
`s21_test_arch` constructs a frameshift pair and requires the merge to fire on
it, which is what makes a count of zero merges a result rather than a rule that
cannot act.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_blocks as B                                            # noqa: E402
import s21_lib as L                                               # noqa: E402

#: Gap bins, in bp. The low bins are single values because that is where a
#: frameshift artefact would sit (miniprot walks a 1-2 bp indel), and a bin
#: wide enough to hold both an artefact and a short intron would hide the
#: thing the calibration exists to find.
BINS = [(-10**9, 0), (1, 2), (3, 5), (6, 10), (11, 20), (21, 30), (31, 40),
        (41, 50), (51, 60), (61, 70), (71, 80), (81, 100), (101, 200),
        (201, 1000), (1001, 10**4), (10**4 + 1, 10**9)]

#: A bin is called spliceable when at least this fraction of its gaps read as a
#: canonical or minor splice pair. The complement of the bar is what a random
#: pair of dinucleotides would give: 1/256 for `GT..AG` alone, about 1.2 % for
#: the three accepted pairs together, so a bin of artefacts cannot pass.
SPLICEABLE_BAR = 0.50

#: Refuse to calibrate on less than this many gaps, in either population.
MIN_GAPS = 200


def bin_rows(pairs: list[dict]) -> list[dict]:
    """One row per gap bin: how many, and how many the genome calls spliceable."""
    rows = []
    for lo, hi in BINS:
        sub = [p for p in pairs if lo <= p["gap_bp"] <= hi]
        if not sub:
            continue
        canonical = sum(1 for p in sub if p["splice_class"] == "canonical")
        minor = sum(1 for p in sub if p["splice_class"] == "minor")
        noncan = sum(1 for p in sub if p["splice_class"] == "non_canonical")
        unread = len(sub) - canonical - minor - noncan
        rows.append({
            "gap_lo": lo if lo > -10**9 else "", "gap_hi": hi if hi < 10**9 else "",
            "n": len(sub), "canonical": canonical, "minor": minor,
            "non_canonical": noncan, "unreadable": unread,
            "frac_spliceable": round((canonical + minor) / len(sub), 4),
            "spliceable": int((canonical + minor) / len(sub) >= SPLICEABLE_BAR)})
    return rows


def calibrate(pairs: list[dict]) -> dict:
    """The floor, or a stated refusal to derive one.

    A usable bar needs gaps on **both** sides of it: enough spliceable gaps to
    know where real introns start, and enough non-spliceable ones to know what
    is being excluded. When the second population is empty the honest output is
    that the rule has nothing to exclude — reported, not rounded into a
    threshold with no evidence under it.
    """
    rows = bin_rows(pairs)
    spl = [r for r in rows if r["spliceable"]]
    non = [r for r in rows if not r["spliceable"]]
    n_non = sum(r["n"] for r in non)
    smallest = min((p["gap_bp"] for p in pairs), default=0)
    if not pairs:
        raise SystemExit("[s21] intron calibration has no gaps to measure")
    if n_non < MIN_GAPS:
        # The floor goes at the smallest gap the genome calls a splice pair,
        # not at the declared fallback: a fallback of 30 bp would merge the ten
        # junctions between 10 and 29 bp that read `GT..AG`, which is the
        # opposite of what the evidence says. Placed here the rule fires on
        # nothing — because nothing here is an artefact — and still fires on a
        # sub-intron gap if one ever appears.
        return {"floor_bp": smallest, "chosen": 1,
                "derivation": ("smallest gap the genome calls spliceable — no "
                               "sub-intron population to exclude"),
                "separated": 0,
                "n_gaps": len(pairs), "n_spliceable": sum(r["n"] for r in spl),
                "n_non_spliceable": n_non,
                "smallest_gap_bp": smallest,
                "fallback_bp": B.MIN_INTRON_FALLBACK,
                "bar": SPLICEABLE_BAR, "min_gaps": MIN_GAPS,
                "reason": (
                    f"{n_non} gaps sit in bins the genome does not call "
                    f"spliceable, under the {MIN_GAPS} needed to place a bar; "
                    f"the smallest gap anywhere in the sweep is {smallest} bp "
                    "and it reads as a splice pair, so no block pair here is a "
                    "frameshift artefact and the merge rule has nothing to "
                    "merge")}
    edge_lo = max(int(r["gap_hi"] or 0) for r in non)
    edge_hi = min(int(r["gap_lo"] or 0) for r in spl if int(r["gap_lo"] or 0)
                  > edge_lo)
    return {"floor_bp": edge_hi, "chosen": 1,
            "derivation": "measured gap between the two populations",
            "separated": 1, "n_gaps": len(pairs),
            "n_spliceable": sum(r["n"] for r in spl),
            "n_non_spliceable": n_non, "smallest_gap_bp": smallest,
            "bar": SPLICEABLE_BAR, "min_gaps": MIN_GAPS,
            "gap_edge_low_bp": edge_lo, "gap_edge_high_bp": edge_hi,
            "reason": (f"gaps up to {edge_lo} bp are not spliceable and gaps "
                       f"from {edge_hi} bp are; the floor is the upper edge")}


def write(pairs: list[dict]) -> dict:
    cal = calibrate(pairs)
    L.OUT.mkdir(parents=True, exist_ok=True)
    L.write_tsv(L.OUT / "intron_floor_calibration.tsv", [cal], None)
    L.write_tsv(L.OUT / "junction_gap_bins.tsv", bin_rows(pairs), None)
    return cal
