"""s5_bait_spec.py — the rules that choose the S5 genomic-sweep bait panel.

The brief asks for "one correctly labelled full-length bait per paralog per
major clade, screened for chimeras (D5), plus one RyR bait per genome as an
internal positive control". This module is that sentence written as rules the
builder *enforces*, in the shape S4 used for the margin species: the panel is
**derived from the committed census**, not hand-listed, so it rebuilds when
the census does and a typo cannot invent a bait.

Six rules, all checked in `s5_build_baits.py`:

 R1 **Census-derived label.** A bait's paralog comes from a gene symbol or
    protein name that names it *and* from the census call agreeing on the
    family (`call == ITPR` for an ITPR bait, `RYR` for a RyR bait). D5's
    "screened by label, not padded for breadth": a bait wearing the wrong
    clade name does real damage, so a record whose label and whose call
    disagree is never a bait.

 R2 **Full length, complete architecture.** An ITPR bait carries all five
    family signatures (`n_itpr_arch == 5`) and sits in `BAIT_BAND_AA`. S2
    found seven records with the complete architecture in 1,528-1,993 aa —
    truncated gene models UniProt does not flag — so the architecture test
    alone does not establish full length and the band is required with it.

 R3 **One bait per paralog per clade band**, two for the bands that dominate
    the genome scope (`DEEP_BANDS`), and the second from a different NCBI
    order than the first. S4's manifest is 152 Aves and 88 Actinopteri out of
    309 genomes; a single Galliform and a single cypriniform bait would put
    most of the sweep at the far end of its nearest bait.

 R4 **The chimera screen is this project's own instrument.** Every candidate
    is scored against `itpr.hmm` and `ryr.hmm` (S3) and must (a) be assigned
    to its own family under S3's relative margin — D14 as a positive test,
    the same rule the census used — and (b) have the winning profile cover at
    least `MIN_ENVELOPE_FRAC` of the bait with no unaligned internal segment
    over `MAX_INTERNAL_GAP_AA`. A fusion or a mis-joined model shows up as
    exactly that gap.

 R5 **The RyR baits are a control, not a census.** One per band, and the
    ledger reports RyR *presence* per genome rather than which RyR. A genome
    where the RyR control finds nothing has an assembly or pipeline problem,
    not a biological one. They are in the same miniprot panel as the ITPR
    baits for a second reason: without a RyR bait present, a genome's own RyR
    locus is claimed by an ITPR bait at partial coverage and reads out as an
    ITPR `fragment`.

 R6 **The S3 seeds are reused, not re-derived.** They already passed the
    stricter S3 selection rules, so they enter the panel as-is; S5 only adds
    the bands and paralogs S3 had no seed for. `vertebrate_basal` seeds
    (Callorhinchus, Petromyzon, Lepisosteus) carry no paralog label and stay
    unlabelled — in a lamprey or chimaera genome the paralog assignment is
    not established, and labelling a bait there would be R1's own failure.

 R7 **The panel is scoped to the sweep.** S4's manifest is 309 *vertebrate*
    genomes, so the S3 seeds from the `invertebrate` and `non_metazoan`
    grades are left out: they cost query residues in every miniprot run and
    cannot win a locus a vertebrate bait does not. S23 sweeps the
    invertebrate, protist and plant/fungal genomes and builds its own panel
    from the same rules with the bands that task needs.

The bands are the classes of `results/genome_manifest.tsv`, collapsed to the
resolution the bait panel can actually distinguish.
"""

from __future__ import annotations

import re

#: manifest `vclass` -> bait clade band.
BAND_OF_CLASS = {
    "Mammalia": "mammalia",
    "Aves": "aves",
    "Lepidosauria": "reptilia",
    "Testudines": "reptilia",
    "Crocodylia": "reptilia",
    "Amphibia": "amphibia",
    "Actinopteri": "actinopteri",
    "Cladistia": "actinopteri",
    "Chondrichthyes": "chondrichthyes",
    "Coelacanthimorpha": "sarcopterygian_fish",
    "Dipnoi": "sarcopterygian_fish",
    "Myxini": "cyclostomata",
    "Hyperoartia": "cyclostomata",
}

#: Fallback for records whose census `class` is empty. UniProt's ranked
#: lineage has no class rank for several deep vertebrate lineages — the
#: coelacanth is the one that matters here, and it silently dropped the
#: Latimeria ITPR1 seed out of the first build of this panel. Resolution goes
#: class -> order, and an unresolved *vertebrate* record is a build error
#: rather than a quiet omission.
BAND_OF_ORDER = {
    "Coelacanthiformes": "sarcopterygian_fish",
    "Ceratodontiformes": "sarcopterygian_fish",
    "Lepidosireniformes": "sarcopterygian_fish",
    "Petromyzontiformes": "cyclostomata",
    "Myxiniformes": "cyclostomata",
    "Squamata": "reptilia",
    "Rhynchocephalia": "reptilia",
    "Testudines": "reptilia",
    "Crocodylia": "reptilia",
}

#: Bands in the order reports and figures use them (crown -> base).
BANDS = ["mammalia", "aves", "reptilia", "amphibia", "actinopteri",
         "chondrichthyes", "sarcopterygian_fish", "cyclostomata"]


def band_of(row: dict) -> str:
    """The bait clade band of a census row, or "" if it is not in scope."""
    band = BAND_OF_CLASS.get((row.get("class") or "").strip())
    if band:
        return band
    return BAND_OF_ORDER.get((row.get("order") or "").strip(), "")

#: R3 — bands taking two baits per paralog, because they dominate the scope.
#: (S4: 152 Aves + 88 Actinopteri of 309 genomes.)
DEEP_BANDS = {"aves", "actinopteri"}

#: R2 — the bait length band. Tighter than `family.MIN_LENGTH_AA` ..
#: `MAX_LENGTH_AA` (2,000-3,600) on purpose: that band admits a real but
#: short family member, and a bait has a stricter job than a census row. The
#: floor is above S2's seven complete-architecture truncations (max 1,993 aa)
#: and the ceiling below the shortest RyR in the S3 seed set (4,715 aa).
BAIT_BAND_AA = (2_400, 3_100)

#: R5 — the RyR control baits' band, from the S3 RyR seeds (4,715-5,317 aa).
RYR_BAND_AA = (4_600, 5_400)

#: R4 — chimera screen thresholds, in profile match states / bait residues.
MIN_ENVELOPE_FRAC = 0.55      # winning profile must cover this much of the bait
MAX_INTERNAL_GAP_AA = 400     # no unaligned run this long inside the envelope

#: The paralogs a bait may be labelled with, per family.
ITPR_PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
RYR_PARALOGS = ("RYR1", "RYR2", "RYR3")

_GENE_RE = re.compile(r"^(itpr|ryr)([123])[ab]?$", re.IGNORECASE)
_NAME_TYPE_RE = re.compile(r"\btype[- ]?([123])\b", re.IGNORECASE)


def paralog_of(gene: str, protein_name: str) -> str:
    """R1 — the paralog a record's *label* names, or "" if it names none.

    Two routes, in order: an unambiguous gene symbol (`Itpr1`, `itpr1b`,
    `RYR2`), else a protein name that spells the family out and gives a type
    number. Nothing else counts — a locus tag, a "receptor-like" name, or a
    type number without the family name leaves the record unlabelled, and an
    unlabelled record cannot be a labelled bait.
    """
    m = _GENE_RE.match((gene or "").strip())
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    low = (protein_name or "").lower()
    m = _NAME_TYPE_RE.search(low)
    if not m:
        return ""
    if "trisphosphate receptor" in low:
        return f"ITPR{m.group(1)}"
    if "ryanodine receptor" in low:
        return f"RYR{m.group(1)}"
    return ""


def family_of(paralog: str) -> str:
    return "RYR" if paralog.startswith("RYR") else "ITPR"


def passes_shape(paralog: str, length: int, n_itpr_arch: str) -> tuple[bool, str]:
    """R2 — length band, and the complete architecture for an ITPR bait."""
    if family_of(paralog) == "RYR":
        lo, hi = RYR_BAND_AA
        if not lo <= length <= hi:
            return False, f"length {length} outside RyR bait band {lo}-{hi}"
        return True, ""
    lo, hi = BAIT_BAND_AA
    if not lo <= length <= hi:
        return False, f"length {length} outside bait band {lo}-{hi}"
    if str(n_itpr_arch) != "5":
        return False, f"architecture {n_itpr_arch}/5, not complete"
    return True, ""


#: R6/R7 — S3 seed clades that enter the S5 panel unlabelled. The deep
#: vertebrate grades where no paralog assignment is established.
UNLABELLED_CLADES = ("vertebrate_basal",)

#: R7 — S3 seed clades left out of this (vertebrate) sweep entirely.
OUT_OF_SCOPE_CLADES = ("invertebrate", "non_metazoan")


def quota(band: str, family: str = "ITPR") -> int:
    """R3/R5 — baits per slot. ITPR: per paralog per band, doubled in the
    bands that dominate the scope. RYR: one per band, whichever paralog,
    because the RyR cell is a presence control rather than a census."""
    if family == "RYR":
        return 1
    return 2 if band in DEEP_BANDS else 1


def rank_candidate(row: dict, reference_len: int) -> tuple:
    """R3 tie-break, best first.

    Reviewed (Swiss-Prot) records first, then the record whose length is
    closest to the paralog's reference length. Closest-to-reference rather
    than longest: the longest record in a band is the one most likely to be a
    fusion or a read-through, which is the failure R4 then has to catch.
    """
    reviewed = 1 if str(row.get("reviewed", "")) == "1" else 0
    length = int(row.get("length") or 0)
    return (-reviewed, abs(length - reference_len), row.get("accession", ""))
