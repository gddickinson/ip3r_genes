"""s23_bait_spec.py — the rules that choose the S23 genomic-sweep bait panel.

S5's panel answers "which of the three vertebrate paralogs is at this locus".
Outside the vertebrates that question does not exist, so three things change,
and none of them is cosmetic.

 B1 **A bait's label is its clade band, never a paralog.** ITPR1/2/3 are a
    2R product; S5b's four `absent` cells are all cyclostome and the trio is
    not established below them. The census question here is **copy number** —
    how many ITPRs a genome has — so a bait carries the band it came from and
    the ledger counts loci rather than filling three named cells.

 B2 **The RyR positive control does not transfer, and pretending it does
    would silently void every negative claim in this task.** In S5 a genome
    with no RyR locus is a broken assembly, because every vertebrate has
    three. S20 measured what happens outside: at architecture level RyR
    appears in exactly two non-metazoan proteomes of 6,928 (*Salpingoeca*,
    *Capsaspora*). A land plant or a yeast with no RyR locus is a *correct*
    result, so RyR cannot be that genome's proof-of-search. This panel
    therefore carries a third role: the **MIR-domain sharer** (protein
    O-mannosyltransferase / dolichyl-phosphate-mannose mannosyltransferase),
    which S1's decoy panel was built around and which S20 used as its
    in-search positive control at proteome level — PF02815 returns 633 hits
    in land plants and 4,376 in Dikarya, the genomes where PF08709 returns
    zero. Taking that same control one level down to the assembly is what
    makes "no ITPR in this genome" a statement about the genome rather than
    about the search.

    It is a *second* control, not a replacement: RyR baits stay in the panel
    for the metazoan genomes, where they remain both the positive control and
    (D14) the decoy that must never win an ITPR locus.

 B3 **The length band is measured here, not inherited.** S5's 2,400-3,100 aa
    band is the vertebrate paralogs' band. The complete-architecture
    non-vertebrate records in census v5 are a different distribution, so the
    band is computed from them at build time by `measure_band()` and written
    into the manifest. `passes_shape()` **raises** until that measurement has
    been made, so a build that skipped it fails instead of quietly falling
    back to the vertebrate numbers.

 B4 **The chimera screen is S5's, unchanged** (D5) — every candidate scored
    against `itpr.hmm` and `ryr.hmm`, assigned to its own family under S3's
    relative margin, with the winning profile covering `MIN_ENVELOPE_FRAC` of
    the bait and no unaligned interior over `MAX_INTERNAL_GAP_AA`. Reusing
    the module rather than copying it is the point: one screen, one set of
    thresholds, and S5's own negative controls still run on every build.

 B5 **The deep bands get a stated exception, not a silent pass.** S0 found
    that *Dictyostelium* iplA — a characterised IP3 receptor — carries none
    of PF08709, PF02815 or PF00520, and S3 put it in the ITPR seed set
    anyway, "precisely so the profile can find its relatives". A band with no
    complete-architecture record may therefore take its best record **by
    profile score**, which is the instrument D14's call is actually made on;
    the row records `arch_exception` and the architecture count it was
    accepted at. A band with no record at all stays an unfilled slot and is
    reported as one.

 B6 **S3's seeds are reused where they fall in a band** (S5's R6): they have
    already passed the stricter S3 selection rules, and re-deriving them
    would let this panel drift from the profiles it is screened by. **A seed
    does not consume its band's quota**, it is added to it. The first build
    let it: Chlorophyta's single slot went to S3's *Cymbomonas* seed and the
    census's own *Chlamydomonas reinhardtii* record — 3,210 aa, complete
    5/5 architecture — never entered the panel. In the pilot the nearest bait
    to the *Chlamydomonas* genome was then an **amoeba**, its alignments came
    back at 25-28 % identity, and the genome read `no_locus`. "Reused, not
    re-derived" means the seeds are kept; it never meant they replace the
    derivation.

 B7 **The panel is scoped to this sweep.** S5's 30 vertebrate ITPR baits are
    left out: they cost query residues in every miniprot run and cannot win a
    locus a nearer bait does not. The reverse of S5's own R7, which left the
    invertebrate and non-metazoan seeds out for the same reason.
"""

from __future__ import annotations

import re

# ------------------------------------------------------------------ bands

#: The clade bands, crown-ward within each kingdom-level group. These are the
#: phyla census v5 actually holds ITPR calls in (`itpr_by_phylum.tsv`), which
#: is why the list is not a textbook phylum list: a band exists because there
#: is something to bait it with.
BANDS: dict[str, str] = {
    # metazoan bands -> the group used for control selection
    "Arthropoda": "metazoa",
    "Nematoda": "metazoa",
    "Mollusca": "metazoa",
    "Annelida": "metazoa",
    "Platyhelminthes": "metazoa",
    "Rotifera": "metazoa",
    "Brachiopoda": "metazoa",
    "Bryozoa": "metazoa",
    "Echinodermata": "metazoa",
    "Hemichordata": "metazoa",
    "Chordata": "metazoa",          # tunicates + cephalochordates only (B7)
    "Cnidaria": "metazoa",
    "Ctenophora": "metazoa",
    "Porifera": "metazoa",
    "Placozoa": "metazoa",
    "Priapulida": "metazoa",
    "Tardigrada": "metazoa",
    "Orthonectida": "metazoa",
    # non-metazoan eukaryote bands
    "Chlorophyta": "viridiplantae",
    "Chytridiomycota": "fungi",
    "Mucoromycota": "fungi",
    "Basidiobolomycota": "fungi",
    "Entomophthoromycota": "fungi",
    "Zoopagomycota": "fungi",
    "Ciliophora": "sar",
    "Oomycota": "sar",
    "Perkinsozoa": "sar",
    "Evosea": "amoebozoa",
    "Discosea": "amoebozoa",
    "Euglenozoa": "discoba",
    "Heterolobosea": "discoba",
    "Haptophyta": "other",
    "Preaxostyla": "other",
}

#: B2 — the groups where a RyR locus is *expected*, so a RyR bait is a
#: positive control there. Everywhere else RyR is a decoy and its absence is
#: a result (S20: 2 architecture-level RyR proteomes outside Metazoa in
#: 6,928). The MIR-sharer control (`CONTROL_MIR`) runs in every genome.
RYR_CONTROL_GROUPS = {"metazoa"}

#: B3 — bands taking two ITPR baits, because they dominate the scope. Set
#: from census v5's own complete-architecture supply: Arthropoda 565,
#: Mollusca 158 and Nematoda 140 records against 37 or fewer everywhere else,
#: and they are likewise the phyla with the most sequenced genomes.
DEEP_BANDS = {"Arthropoda", "Mollusca", "Nematoda"}

#: Band ordering for reports and figures.
BAND_ORDER = list(BANDS)

#: The three bait roles. `CONTROL_MIR` baits are not family members: they are
#: the proof-of-search (B2) and simultaneously the sharpest available decoy,
#: since MIR is the one family domain every eukaryote carries on something
#: else.
FAMILY_ITPR = "ITPR"
FAMILY_RYR = "RYR"
CONTROL_MIR = "MIR"

#: B4 — the chimera-screen thresholds. Imported from S5's spec rather than
#: retyped, so the two panels cannot be screened to different standards.
from s5_bait_spec import (MAX_INTERNAL_GAP_AA,  # noqa: E402,F401
                          MIN_ENVELOPE_FRAC)

#: B2 — the MIR control bait band, in aa. The floor is the PF02815 model's own
#: length (185 match states, rounded up): a bait shorter than the domain it is
#: supposed to carry is a fragment, and above that a short bait is an easier
#: miniprot target rather than a harder one.
#:
#: **It was 600 aa on the first build and that was wrong.** 600 selects the
#: protein O-mannosyltransferases, which is what the control is *usually*
#: made of — and it left five clades with no control at all, among them both
#: apicomplexan classes, i.e. exactly the negative claim S20a asked this task
#: to take to genome level. Those clades are not free of MIR proteins; they
#: carry short MIR-only ones (178-455 aa). A threshold that quietly deletes
#: the hardest negative claims is a threshold chosen for the wrong thing.
MIR_BAND_AA = (200, 2_400)

#: B2 — what makes a control *strong*. A named mannosyltransferase of at
#: least this length is the full PMT/DPM protein the control was designed
#: around; anything else is a short MIR-only protein, which still proves the
#: search runs but proves less of it. The distinction is recorded per clade
#: rather than resolved, because it is the reader's to weigh: an absence
#: claim standing on a 233 aa control in a clade where one proteome of 36
#: carries any MIR protein is a weaker claim than one standing on a 950 aa
#: mannosyltransferase found across a whole class.
CONTROL_STRONG_MIN_AA = 600

#: B2 — and how much of the clade the control has to be present in. Below
#: this fraction of the clade's swept proteomes carrying any PF02815 protein,
#: the control speaks for the clade only weakly.
CONTROL_STRONG_MIN_CLADE_FRAC = 0.25

#: B5 — the percentile of the *complete-architecture* records' own profile
#: scores that an architecture exception has to clear. Stated as a rule rather
#: than as a number: a record admitted without the architecture to prove it
#: must score at least as well against `itpr.hmm` as the weakest 5 % of the
#: records that do have it. Measured at build time by
#: `measure_arch_exception_floor()` — on census v5 it lands near 1,100 bits,
#: where a typed floor would have been set an order of magnitude too low.
ARCH_EXCEPTION_PERCENTILE = 0.05

_ARCH_FLOOR: float | None = None

# --------------------------------------------------------- measured band

_MEASURED_BAND: tuple[int, int] | None = None

#: B3 — the central fraction of the measured length distribution the band
#: spans, and the rounding applied to its edges.
BAND_CENTRAL_FRAC = 0.90
BAND_ROUND_AA = 50

#: B3/B5 — how many records **per band** enter the measurement of the band and
#: the floor. Both were first measured over the whole in-scope population, and
#: that population is 53 % Arthropoda: the band came out 2,550-3,100 aa and the
#: floor 1,095 bits, which are the arthropod family's numbers wearing the whole
#: tree's name. They excluded *Bodo saltans* — a 2,356-4,222 aa euglenozoan
#: ITPR at 4/5 architecture and 780 bits, which is exactly the deep-branching
#: record a non-vertebrate panel exists to carry.
#:
#: Capping the contribution per band is the same correction S20 made when it
#: stratified its bacterial sample by genus, and the same one G3 and G5 make
#: here. A threshold measured on the best-sampled clade and applied to the
#: whole tree does not describe the family; it describes the sampling.
MEASURE_CAP_PER_BAND = 20


def stratify(by_band: dict[str, list]) -> list:
    """The measurement population, capped at `MEASURE_CAP_PER_BAND` per band.

    Takes the *shortest* records within each band's cap for the length band and
    the *lowest-scoring* for the floor — whichever list is passed — because
    both thresholds are lower bounds and the question they answer is "how
    little may a real family member be", not "how much".
    """
    out = []
    for band, xs in sorted(by_band.items()):
        out.extend(sorted(xs)[:MEASURE_CAP_PER_BAND])
    return out


def measure_band(lengths: list[int]) -> tuple[int, int]:
    """B3 — the ITPR bait length band, from the census's own lengths.

    The central `BAND_CENTRAL_FRAC` of the complete-architecture
    non-vertebrate ITPR lengths, rounded *outward* to `BAND_ROUND_AA`. A
    percentile rather than min/max because the tails of that distribution are
    the very things a bait must not be: S2 found seven records carrying the
    complete architecture in 1,528-1,993 aa (truncated gene models UniProt
    does not flag), and the top of the range runs into read-through fusions.
    """
    if len(lengths) < 100:
        raise ValueError(f"band measured on only {len(lengths)} lengths; "
                         "the census should supply hundreds")
    # sanity: the band must admit the deepest bands' records, which is the
    # failure the stratification exists to prevent.
    xs = sorted(lengths)
    tail = (1.0 - BAND_CENTRAL_FRAC) / 2.0
    lo = xs[int(len(xs) * tail)]
    hi = xs[min(len(xs) - 1, int(len(xs) * (1.0 - tail)))]
    lo = (lo // BAND_ROUND_AA) * BAND_ROUND_AA
    hi = -((-hi) // BAND_ROUND_AA) * BAND_ROUND_AA
    return int(lo), int(hi)


def measure_arch_exception_floor(scores: list[float]) -> float:
    """B5 — the score floor, from the complete-architecture records' scores.

    `scores` is every in-band, complete-architecture, non-vertebrate ITPR
    record's `itpr.hmm` full-sequence bit score. The floor is the
    `ARCH_EXCEPTION_PERCENTILE` quantile of them.
    """
    xs = sorted(s for s in scores if s > 0)
    if len(xs) < 100:
        raise ValueError(f"floor measured on only {len(xs)} scores; "
                         "the census should supply hundreds")
    return float(xs[int(len(xs) * ARCH_EXCEPTION_PERCENTILE)])


def set_measured_band(band: tuple[int, int], arch_floor: float) -> None:
    global _MEASURED_BAND, _ARCH_FLOOR
    _MEASURED_BAND = band
    _ARCH_FLOOR = float(arch_floor)


def arch_exception_floor() -> float:
    if _ARCH_FLOOR is None:
        raise RuntimeError("the B5 architecture-exception floor has not been "
                           "measured; call set_measured_band(...) first")
    return _ARCH_FLOOR


def measured_band() -> tuple[int, int]:
    if _MEASURED_BAND is None:
        raise RuntimeError(
            "the S23 bait length band has not been measured. Call "
            "set_measured_band(measure_band(lengths)) with the census's "
            "complete-architecture non-vertebrate ITPR lengths first — B3 "
            "exists so this panel cannot silently inherit S5's vertebrate "
            "band (2,400-3,100 aa).")
    return _MEASURED_BAND


# ------------------------------------------------------------- band lookup

#: Census phylum names that are not the band name. `unclassified` rows carry
#: no phylum at all and are resolved by group instead.
_PHYLUM_ALIAS = {
    "Bacillariophyta": "",          # diatoms: no ITPR call to bait with
}


def band_of(row: dict) -> str:
    """The bait clade band of a census row, or "" if it is not in scope.

    Resolution is NCBI's lineage first (`tax_phylum`, written by S20 from the
    taxonomy API) and UniProt's ranked lineage second (`phylum`), because the
    two disagree on several protist groups and S20's is the one the sweep's
    own manifest will be built from. A **vertebrate** row is out of scope by
    B7 and returns "" rather than raising: S5 owns that panel.
    """
    if (row.get("tax_group") or row.get("group") or "") == "Vertebrata":
        return ""
    for key in ("tax_phylum", "phylum"):
        name = (row.get(key) or "").strip()
        name = _PHYLUM_ALIAS.get(name, name)
        if name in BANDS:
            # Chordata spans both sweeps; only the non-vertebrates are ours.
            return name
    return ""


def group_of_band(band: str) -> str:
    return BANDS.get(band, "")


def expects_ryr(band: str) -> bool:
    """B2 — is a RyR locus expected in this band's genomes?"""
    return group_of_band(band) in RYR_CONTROL_GROUPS


# ------------------------------------------------------------------ shape

_MIR_NAME_RE = re.compile(
    r"mannosyltransferase|\bpmt\d?\b|dolichyl-phosphate-mannose",
    re.IGNORECASE)


def control_strength(named: int, length: int, clade_frac: float
                     ) -> tuple[str, str]:
    """B2 — is this control strong, and if not, what is missing from it?"""
    problems = []
    if not named:
        problems.append("not named a mannosyltransferase")
    if length < CONTROL_STRONG_MIN_AA:
        problems.append(f"{length} aa — a short MIR-only protein, not the "
                        f"full PMT (< {CONTROL_STRONG_MIN_AA})")
    if clade_frac < CONTROL_STRONG_MIN_CLADE_FRAC:
        problems.append(f"PF02815 found in only {clade_frac:.0%} of the "
                        "clade's swept proteomes")
    if not problems:
        return "strong", ("a named mannosyltransferase of full length, "
                          f"present across {clade_frac:.0%} of the clade")
    return "weak", "; ".join(problems)


def is_mir_sharer_name(protein_name: str) -> bool:
    """B2 — does this protein's own name say it is the control, not the family?

    A positive test on the record's paperwork, used only to *choose* control
    baits. Their family status is then settled by the screen like every other
    bait: a MIR bait must be assigned to neither family, which is exactly what
    makes it a decoy.
    """
    return bool(_MIR_NAME_RE.search(protein_name or ""))


def passes_shape(family: str, length: int, n_itpr_arch: str) -> tuple[bool, str]:
    """B3/B5 — the length band, and the architecture an ITPR bait must carry.

    Signature matches `s5_bait_spec.passes_shape` so `s5_bait_screen` can run
    its own negative controls against this spec unchanged.
    """
    length = int(length or 0)
    if family == CONTROL_MIR:
        lo, hi = MIR_BAND_AA
        if not lo <= length <= hi:
            return False, f"length {length} outside MIR control band {lo}-{hi}"
        return True, ""
    if family == FAMILY_RYR:
        from s5_bait_spec import RYR_BAND_AA
        lo, hi = RYR_BAND_AA
        if not lo <= length <= hi:
            return False, f"length {length} outside RyR bait band {lo}-{hi}"
        return True, ""
    lo, hi = measured_band()
    if not lo <= length <= hi:
        return False, f"length {length} outside measured bait band {lo}-{hi}"
    if str(n_itpr_arch) != "5":
        return False, f"architecture {n_itpr_arch}/5, not complete"
    return True, ""


def passes_shape_with_exception(family: str, length: int, n_itpr_arch: str,
                                itpr_score: float) -> tuple[bool, str, bool]:
    """B5 — `passes_shape`, plus the deep-band architecture exception.

    Returns (ok, reason, is_exception). The exception relaxes **only** the
    architecture test, never the length band: a record that is the wrong size
    is not a deep homolog, it is a broken model, and that is the failure S2
    measured seven times.
    """
    ok, why = passes_shape(family, length, n_itpr_arch)
    if ok or family != FAMILY_ITPR:
        return ok, why, False
    if "architecture" not in why:
        return False, why, False
    # An exception needs something to be an exception *to*. A record whose
    # architecture was never measured — the 287 census rows that exist only
    # because this project's own profile sweep found them, with no InterPro
    # annotation behind them — is not a deep homolog the signatures missed,
    # it is a record the signatures were never run on, and the two must not
    # be spent from the same budget.
    if not str(n_itpr_arch or "").strip():
        return False, ("architecture never measured (no InterPro annotation); "
                       "B5's exception applies only to a *measured* "
                       "incomplete architecture"), False
    floor = arch_exception_floor()
    if float(itpr_score or 0.0) < floor:
        return False, (f"{why}, and itpr.hmm scores it only "
                       f"{float(itpr_score or 0):.0f} bits (< {floor:.0f}, "
                       "the 5th percentile of the complete-architecture "
                       "records' own scores)"), False
    return True, (f"architecture {n_itpr_arch}/5 accepted as a stated "
                  f"exception on itpr.hmm = {float(itpr_score):.0f} bits "
                  f"(>= {floor:.0f}) — the S0/S3 iplA case, where the "
                  "family's own defining signature does not reach a "
                  "characterised receptor"), True


#: B3 — the rank a band's baits must differ at. S5's R3 required its second
#: bait to come from a different NCBI **order** than the first, and that rule
#: was not ported here. It cost the first panel its *Chlamydomonas* bait:
#: Chlorophyta took two *Cymbomonas* records — same genus, a prasinophyte —
#: while *Chlamydomonas reinhardtii*, a chlorophycean with a documented
#: complete 5/5 receptor, got none. A band here is a whole **phylum**, far
#: wider than S5's class-level bands, so the spread matters more, not less.
SPREAD_RANK = "order"


def spread_key(row: dict) -> str:
    """The value two baits in one band must differ at (B3's spread rule)."""
    for key in ("tax_order", "order", "tax_class", "class"):
        v = (row.get(key) or "").strip()
        if v:
            return v
    # No ranked order (common for protists): fall back to the genus, taken as
    # the first token of the species name. A stratification key, not a claim.
    return (row.get("species") or "").split()[0] if row.get("species") else ""


def quota(band: str, family: str = FAMILY_ITPR) -> int:
    """B3/B2 — baits per slot."""
    if family == FAMILY_ITPR:
        return 2 if band in DEEP_BANDS else 1
    return 1


def rank_candidate(row: dict, reference_len: int) -> tuple:
    """B3 tie-break, best first.

    Reviewed (Swiss-Prot) records first, then the *highest-scoring* record
    against `itpr.hmm`, then the length closest to the band's centre.

    Score before length, which is the opposite of S5's rule and deliberate.
    There, every candidate was a labelled full-length paralog and the risk was
    picking a fusion, so closest-to-reference was right. Here most records
    carry no label at all and the risk is picking a well-sized fragment of
    something else, so the instrument that makes the call gets to rank the
    candidates it will later be asked to defend.
    """
    reviewed = 1 if str(row.get("reviewed", "")) == "1" else 0
    try:
        score = float(row.get("itpr_score") or 0.0)
    except ValueError:
        score = 0.0
    length = int(row.get("length") or 0)
    return (-reviewed, -score, abs(length - reference_len),
            row.get("accession", ""))
