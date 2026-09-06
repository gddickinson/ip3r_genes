"""The eight rules that choose the S6 representative set, and why.

D8 says representatives are chosen per clade *and* per kingdom by explicit
coded rules, never "the longest sequence per species" — a longest-first pick
reliably selects chimeric gene models. This module is those rules, written
out and enforced, in the style of `s3_seed_spec.py` / `s5_bait_spec.py` /
`s23_bait_spec.py`. The driver (`s6_select_reps.py`) applies them and writes
the audit; nothing here fetches, aligns or plots.

**What S6 inherits, and what it changes.**

*From S23 (D27/D29): identity to a bait is not a usable statistic outside
the vertebrates.* The inherited 0.40 floor discarded 87 confirmed loci, and
neither identity nor coverage separates the confirmed and contradicted
populations (Youden J 0.70/0.71). So **no rule here ranks or filters on
sequence identity.** Quality is judged on architecture, fragment status,
profile evidence and the record's own paperwork — evidence that does not
come from the same measurement the selection is trying to be independent of.

*From S23: copy number outside the vertebrates ranges 0 to 18.* A rule that
gives every non-vertebrate species one slot cannot represent a genome
carrying eighteen genes, and a tree built from one tip per expanded genome
cannot say whether the expansion is lineage-specific. **R6** spends slots on
the expansions deliberately.

*From S5/S23: the clade band, not the class rank.* UniProt's ranked lineage
has no class for Testudines, Crocodylia or Coelacanthiformes, and a
class-only lookup silently dropped the *Latimeria* seed from the first S5
panel. `s5_bait_spec.band_of()` is imported rather than re-derived.

**The one thing this set exists for.** S7 roots the family tree on the
ryanodine receptors and asks which two vertebrate paralogs are sisters. A
representative set without an outgroup has no root (**R8**), and one without
the pre-2R grade (**R3**) has nothing between the outgroup and the crown
trio to place the duplications against. Both are rules, not garnish.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from s5_bait_spec import BANDS, band_of, paralog_of  # noqa: E402,F401
from src.utils import family  # noqa: E402

# --------------------------------------------------------------- R1 grid

#: R1 — the vertebrate paralog grid: clade band x ITPR1/2/3.
PARALOGS = tuple(family.PARALOGS)

#: R1 — bands taking two records per paralog. The three that dominate the
#: vertebrate census (Actinopteri 2,292 / Aves 1,808 / Mammalia 1,349 ITPR
#: records of 6,112): a single tip for a clade holding a third of the
#: evidence is a sampling choice the tree would inherit as a topology.
DEEP_BANDS = {"mammalia", "aves", "actinopteri"}

#: R1 — forced in: the three human paralogs are the reference panel every
#: earlier task measured against (`family.REFERENCE_PANEL`), and the clinical
#: variant work (S11/S22) is stated on human numbering. A representative set
#: that resolved them away by score would leave every later task translating
#: between coordinate systems.
FORCED_ITPR = {acc: label for label, acc in family.REFERENCE_PANEL}

#: R8 — forced in: the human sister trio, which roots the tree.
FORCED_RYR = {acc: label for label, acc in family.SISTER_PANEL}

# ------------------------------------------------------------ R2 teleost

#: R2 — how many teleost species contribute both 3R co-orthologs. The
#: teleost genome duplication is the family's only other documented
#: duplication, and `itpr1a`/`itpr1b` pairs are in the census for 52
#: species x paralog cells. Two species is enough to show the 3R node
#: without the crown teleosts outvoting the rest of the vertebrates.
TELEOST_3R_SPECIES = 2

# ---------------------------------------------------- R3 unlabelled grade

#: R3 — the bands whose records enter unlabelled, because no paralog
#: assignment is established there. Cyclostomes are the pre-2R grade; the
#: chondrichthyan and sarcopterygian-fish grades sit either side of the
#: paralog radiation and S5 could fill no labelled ITPR1/ITPR2 bait for
#: them at all. The tree assigns these tips; this script does not.
UNLABELLED_BANDS = ("cyclostomata", "chondrichthyes", "sarcopterygian_fish")

#: R1/R3 — **when a bait attribution counts as a paralog label, and when
#: it does not.** The sweeps record which bait won each locus (`cell`), and
#: for most vertebrate clades that is a usable paralog label: the three
#: baits win different loci, so the winner discriminates. In the cyclostome
#: band it does not. The sea lamprey and the inshore hagfish carry three
#: ITPR loci each and **the ITPR1 bait wins all six** — the attribution is
#: constant, so it carries no information about which locus is which, and
#: filling an `ITPR1@cyclostomata` grid cell from it would assert an
#: orthology the data does not contain. `informative_attribution()` decides
#: this per band by counting distinct attributed cells, so it is measured
#: from the census rather than a band being named here as the exception.
def informative_attribution(cells: set[str]) -> bool:
    return len({c for c in cells if c in PARALOGS}) > 1


#: R3 — **one tip per distinct locus, not per species.** Measured in the
#: census: the sea lamprey carries three ITPR loci (LOC116947491,
#: LOC116951169, LOC116952798) and the inshore hagfish three, and S5's
#: sweep assigned *all six* to the ITPR1 cell — the ITPR1 bait wins every
#: one of them, so best-bait attribution does not separate them. Whether
#: those three are 1:1 orthologs of ITPR1/2/3 or a lineage-specific
#: expansion is the 2R question, and a set carrying one lamprey tip cannot
#: answer it: two thirds of the evidence would not be in the alignment.
#: The cap is the vertebrate paralog count rather than a typed 3, and a
#: locus is a distinct gene symbol (or, for a record with none, its own
#: accession).
UNLABELLED_LOCI_PER_SPECIES = len(PARALOGS)
UNLABELLED_SPECIES_PER_BAND = 2

# --------------------------------------------------------- R4/R5 phylum

#: R4 — extra slots in the non-vertebrate metazoan phyla that dominate that
#: part of the census (Arthropoda 831, Mollusca 331, Nematoda 292,
#: Platyhelminthes 152 of 1,990 records). Spread across classes.
METAZOAN_DEEP_PHYLA = {"Arthropoda": 4, "Mollusca": 2, "Nematoda": 2,
                       "Platyhelminthes": 2, "Chordata": 2}

#: R4/R5 — the default: one representative per phylum.
PER_PHYLUM = 1

#: R5 — the non-vertebrate, non-metazoan groups, in report order. Every one
#: of these is a range claim S20/S23 made, so each needs a tip: a tree that
#: omits the plant grade cannot be cited for where the family reaches.
EUKARYOTE_GROUPS = ["SAR", "Discoba", "Amoebozoa", "Viridiplantae",
                    "Fungi", "Eukaryota (other)"]

#: R5 — extra slots for the groups whose ITPR records are numerous enough
#: that one tip would be a coin toss (SAR 609 + 17 records, spread over
#: Ciliophora / Apicomplexa / Bacillariophyta / Oomycota).
EUKARYOTE_DEEP_GROUPS = {"SAR": 5, "Discoba": 3, "Viridiplantae": 3,
                         "Amoebozoa": 2, "Fungi": 2, "Eukaryota (other)": 3}

#: R5 — the S20 verdict a plant or fungal record must carry to be eligible.
#: S20 chased every plant and fungal ITPR record individually and wrote a
#: verdict on each; `contaminant_suspect` records are 95-100 % identical to
#: one particular animal, and putting one in the tree would place an
#: assembly contamination artifact on the plant branch and then read the
#: branch as a result. `fragment` and `no_genome_backing` are excluded for
#: the same reason a fragment is excluded everywhere else: it contributes
#: gaps to every column the tree is inferred from. A row with no verdict
#: (the genomic models S23 added, which S20 never saw) is eligible — the
#: gate rejects a *bad* verdict, it does not require a verdict.
PLANT_FUNGAL_REJECT = {"contaminant_suspect", "fragment", "no_genome_backing",
                       "module_only", "unresolved"}


def plant_fungal_ok(verdict: str) -> bool:
    return (verdict or "").strip() not in PLANT_FUNGAL_REJECT


# ------------------------------------------------------- R6 copy number

#: R6 — genomes whose measured copy number exceeds the vertebrate paralog
#: count get more than one tip. Threshold is `len(PARALOGS)` rather than a
#: typed 3, so it moves with the family definition.
COPY_NUMBER_THRESHOLD = len(PARALOGS)

#: R6 — how many genomes are sampled this way, and how many copies each
#: contributes. The cap exists because *Macrostomum lignano* alone carries
#: 18 full copies: uncapped, one flatworm would be 12 % of the alignment.
COPY_GENOMES = 4
COPY_PER_GENOME = 4

# ------------------------------------------------- R7 novel gene models

#: R7 — the project's own novel evidence. S5 and S23 built 1,058 + 183
#: genomic gene models, and census v4/v6 record *how* a database holds each
#: locus; `genome_only` means no database annotates it at all. 241 of those
#: are full-length. A representative set that omits them omits the only
#: sequences in this project that no public database contains, and S7 could
#: not then ask the question they exist for: does a locus nobody has named
#: fall where its bait says it should?
NOVEL_DB_STATUS = "genome_only"
NOVEL_VERTEBRATE = 8          # at most one per (clade band, paralog cell)
NOVEL_NONVERTEBRATE = 3       # at most one per census group

# ----------------------------------------------------------- R8 outgroup

#: R8 — non-vertebrate RyR outgroup tips, beyond the forced human trio.
#: A vertebrate-only outgroup cannot root the non-vertebrate ITPR grade:
#: the branch it defines would be inside the vertebrates.
RYR_NONVERT = 3

#: R8 — the RyR length band, from `s5_bait_spec.RYR_BAND_AA`'s reasoning:
#: a 1,000 aa "RyR" is a fragment, and a fragment cannot root anything.
RYR_MIN_AA = 4_000

# ------------------------------------------------------- shape and rank

#: The length band a representative must fall in. `family.MIN_LENGTH_AA` ..
#: `MAX_LENGTH_AA` (2,000-3,600) is the family's own band and is used
#: unchanged for the ITPR tips: unlike a bait, a representative is not
#: required to be exemplary, only to be a full-length gene. Records below
#: the floor are fragments, and a fragment contributes gaps to every column
#: the tree is inferred from.
ITPR_MIN_AA = family.MIN_LENGTH_AA
ITPR_MAX_AA = family.MAX_LENGTH_AA

#: The one deliberate exception, and it is the same one S3's seed manifest
#: makes: a clade whose *only* records are short still needs a tip, or the
#: range claim S20/S23 made about it has no representative in the tree.
#: A rule that admits such a record records it as `short_exception` and the
#: driver reports the count; it is never silent.
SHORT_EXCEPTION_MIN_AA = 1_200


def passes_shape(family_call: str, length: int, allow_short: bool
                 ) -> tuple[bool, str]:
    """Length gate. Returns (ok, reason-if-not | note-if-exception)."""
    if family_call == "RYR":
        if length < RYR_MIN_AA:
            return False, f"RyR {length} aa below {RYR_MIN_AA} outgroup floor"
        return True, ""
    if ITPR_MIN_AA <= length <= ITPR_MAX_AA:
        return True, ""
    if length > ITPR_MAX_AA:
        return False, f"length {length} above family band {ITPR_MAX_AA}"
    if allow_short and length >= SHORT_EXCEPTION_MIN_AA:
        return True, "short_exception"
    return False, f"length {length} below family band {ITPR_MIN_AA}"


#: Protein-existence ranks, best first. UniProt's own confidence in the
#: record existing at all; a `Predicted` entry is a gene-caller's opinion.
PE_RANK = {
    "Evidence at protein level": 4,
    "Evidence at transcript level": 3,
    "Inferred from homology": 2,
    "Predicted": 1,
    "": 0,
}


def rank_key(row: dict, length_target: float, bait_ids: set[str]) -> tuple:
    """The ordered quality key a cell's winner is chosen by.

    **Deliberately lexicographic, not a weighted sum.** A weighted sum hides
    which criterion decided, and the audit table has to be able to say. The
    components, in order:

      1. not a fragment          — a fragment aligns as gaps
      2. complete architecture   — all five ITPR signatures (D14b's own test)
      3. inside the length band  — support, never the call (D14)
      4. profile confidence high — both instruments agree (D23)
      5. Swiss-Prot reviewed     — a curator has read the record
      6. already a bait or seed  — R7, continuity with the searches
      7. protein-existence rank  — UniProt's confidence the gene is real
      8. length fit              — |len - the group's own median|, inverted;
                                   *not* longest-first (D8), and the target
                                   is the group's median so the plant grade
                                   is not scored against a vertebrate ruler
      9. profile score           — the last tie-break, and only that

    Identity appears nowhere: S23 (D29) retired it as a call gate for this
    family outside the vertebrates, and a selection rule stated in a
    statistic that does not separate the populations is not a rule.
    """
    from s6_lib import as_float, as_int
    length = as_int(row.get("length"))
    return (
        0 if (row.get("fragment") or "").strip() else 1,
        1 if str(row.get("n_itpr_arch")).strip() == "5" else 0,
        1 if str(row.get("in_band")).strip() == "1" else 0,
        1 if (row.get("profile_confidence") or "") == "high" else 0,
        1 if str(row.get("reviewed")).strip() == "1" else 0,
        1 if row.get("accession") in bait_ids else 0,
        PE_RANK.get((row.get("protein_existence") or "").strip(), 0),
        -abs(length - length_target),
        as_float(row.get("itpr_score")),
    )


#: The names of `rank_key`'s components, so the audit table can say which
#: one separated the winner from the runner-up.
RANK_COMPONENTS = ["not_fragment", "complete_architecture", "in_length_band",
                   "profile_high", "reviewed", "is_bait", "protein_existence",
                   "length_fit", "profile_score"]


def decisive_component(win: tuple, runner: tuple | None) -> str:
    """Which component of `rank_key` actually chose the winner."""
    if runner is None:
        return "sole_candidate"
    for name, a, b in zip(RANK_COMPONENTS, win, runner):
        if a != b:
            return name
    return "tie"


# --------------------------------------------------------------- groups

#: Figure/reporting group for a chosen row. Paralog identity where it is
#: established, the taxonomic grade otherwise — which is exactly what
#: `figstyle.GROUP` is built to colour.
def group_of(family_call: str, paralog: str, census_group: str,
             is_vertebrate: bool) -> str:
    """**A paralog label is only meaningful inside the vertebrates.**

    ITPR1/2/3 are a 2R product. A *Tetrabaena* or *Achlya* record whose
    UniProt protein name reads "receptor type 2" carries that number by
    annotation transfer from a vertebrate, not by descent from the
    vertebrate paralog, and colouring such a tip orange in the tree figure
    would assert the very thing S7 is being run to test. Outside the
    vertebrates the group is the taxonomic grade, whatever the name says;
    the name-derived label is still written to the audit table's `paralog`
    column with `paralog_source` naming where it came from, so the transfer
    is visible rather than deleted.
    """
    if family_call == "RYR":
        return "RYR"
    if is_vertebrate:
        return paralog if paralog in PARALOGS else "vertebrate_basal"
    if census_group == "Metazoa (non-vertebrate)":
        return "invert_metazoa"
    if census_group == "Viridiplantae":
        return "plant"
    if census_group == "Fungi":
        return "fungi"
    return "protist"


GROUP_ORDER = ["ITPR1", "ITPR2", "ITPR3", "vertebrate_basal",
               "invert_metazoa", "plant", "protist", "fungi", "RYR"]
