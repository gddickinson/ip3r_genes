"""s23_scope.py — S23's genome scope, and why it is that scope.

S4 built the vertebrate denominator: one best assembly per order, union every
margin species the census left undecided, with every margin species *derived
from a committed table* rather than hand-listed. This is the same discipline
pointed at the rest of the eukaryotes, where the question is different.

S20 swept 6,928 reference **proteomes** and found the family absent from
Streptophyta (0/384), Dikarya (0/1,353), Apicomplexa (0/60) and a scatter of
fungal phyla. Every one of those is an *annotation* fact. A proteome is what
a gene-caller found; a genome is what is there. Taking those absences to
assembly level is the whole point of this task, and it sets the scope rule:

    **sample most finely where the negative claim is.**

Five rules. Each writes the reason it fired into the manifest row, and a
genome that several rules reach is one row carrying all of them (S4's
`merge_rows` pattern) — so the manifest's row count is the number of genomes
actually downloaded.

 G1 `phylum_rep`     One best assembly per non-vertebrate eukaryotic phylum.
                     Guarantees no phylum is unrepresented, including the ones
                     with a single sequenced species.

 G2 `class_rep`      One per class inside any phylum S20 swept at least
                     `BIG_PHYLUM_SWEPT` proteomes of. Resolution follows how
                     much of the tree the phylum actually is, measured on this
                     project's own sweep rather than on an opinion about which
                     phyla are big.

 G3 `absence_clade`  One per class inside every clade S20 swept at least
                     `MIN_SWEPT_FOR_CLAIM` proteomes of and returned **zero**
                     ITPR calls from. This is the rule that exists for this
                     task: it puts several independent assemblies under every
                     absence, so "no ITPR in the land plants" is not one
                     genome's word. Derived from
                     `results/s20_sweep/proteome_presence.tsv` joined to
                     `taxonomy.tsv`, never typed — which is how it also picks
                     up the absences S20a's summary did not name
                     (Bacillariophyta, Rhodophyta, Cestoda).

 G4 `anchor`         The reference organisms this project's claims are
                     calibrated on, each with its reason in the row. The one
                     hand-written list in the scope, and deliberately so: an
                     anchor is a genome whose answer is known from the
                     literature, which is exactly what no table can derive.

 G5 `copy_number`    The species with the most ITPR records in each class,
                     wherever that count exceeds `COPY_NUMBER_QUESTION` — more
                     copies than any *vertebrate* has paralogs. Outside the
                     vertebrates the expectation is a single Itpr (S1's
                     control panel is the "single-Itpr grade"), so an
                     expansion is a question, and S20a left one open by name
                     (*Cymbomonas*, 36 records). Stratified by class for the
                     same reason G3 is: unstratified it returns five
                     *Paramecium* species and calls that five findings.

**The size guard** is not in this module's rules because it is not a rule
about taxonomy. `s23_build_manifest.py` applies it: an assembly whose contig
N50 is below the measured non-vertebrate ITPR gene span cannot carry the gene
on one contig, so it is ranked below any assembly that can, and a clade whose
*best* assembly is still below the bar is admitted and flagged rather than
dropped — D4's discipline, with the bar measured for this scope by
`s23_span_calibration.py` instead of inherited from the vertebrate one.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
S20_DIR = PROJECT_ROOT / "results" / "s20_sweep"

#: The S20 groups that partition the non-vertebrate eukaryotes. The prokaryote
#: groups are out of scope: S20 swept 4,171 archaeal and bacterial proteomes
#: for 0 ITPR, and no one disputes that result hard enough to spend assemblies
#: on it.
EUKARYOTE_GROUPS = ("metazoa_nonvert", "fungi", "viridiplantae",
                    "protista_other")

#: G2 — a phylum this many swept proteomes deep is sampled per class.
BIG_PHYLUM_SWEPT = 100

#: G3 — the smallest swept sample a zero-ITPR clade needs before the absence
#: is worth taking to assembly level. Below it the proteome-level negative is
#: not yet a claim, so putting genomes under it would be spending compute on
#: noise.
MIN_SWEPT_FOR_CLAIM = 10

#: G5 — the copy number that makes a species a question. The vertebrate
#: paralog count: a non-vertebrate genome carrying **more ITPR records than
#: any vertebrate has paralogs** is an expansion, not a grade. Two would be
#: the literal contrast with S1's "single-Itpr grade", but it selects 264
#: species and most of them are one alternative gene model; above three the
#: signal is unambiguous and the set is small enough to sequence.
COPY_NUMBER_QUESTION = 4

#: G5 — the rank the copy-number rule is stratified at, matching G2/G3. One
#: genome per class, the species with the most records in it: without a
#: stratification the rule returns five *Paramecium* species at 36-39 copies
#: and calls that five findings.
COPY_NUMBER_RANK = ("class", "phylum", "kingdom")

#: G4 — the anchors, with why each one is here. Keyed by NCBI taxid so the
#: manifest resolves them the same way every other row is resolved.
#:
#: The brief names eight of these ("*Arabidopsis*, rice, a moss, a
#: chlorophyte, *S. cerevisiae*, *Neurospora*, a chytrid and a
#: microsporidian"); the positive controls are added because a panel of
#: nothing but expected-negatives cannot show that the pipeline works.
ANCHORS: dict[int, str] = {
    # positive controls — a characterised ITPR whose gene is known
    7227: "Drosophila melanogaster — Itpr, the S3 seed and S1 recall case",
    6239: "Caenorhabditis elegans — itr-1, the second S3 invertebrate seed",
    45351: "Nematostella vectensis — the cnidarian outgroup to bilaterian Itpr",
    7668: "Strongylocentrotus purpuratus — the echinoderm IP3R of S3's seeds",
    44689: "Dictyostelium discoideum — iplA, the receptor the family's own "
           "defining Pfam does not reach (S0)",
    3055: "Chlamydomonas reinhardtii — the chlorophyte reference, the green "
          "lineage where S20 does find the family",
    # the brief's named negatives — where the absence claim has to hold
    3702: "Arabidopsis thaliana — the land-plant reference; Streptophyta 0/384",
    39947: "Oryza sativa Japonica — the second land-plant reference, a "
           "monocot against Arabidopsis' eudicot",
    3218: "Physcomitrium patens — a moss, the deepest-branching land plant "
          "with a reference genome",
    559292: "Saccharomyces cerevisiae S288C — the Dikarya reference; "
            "Ascomycota 0/1,034",
    367110: "Neurospora crassa OR74A — the second Dikarya reference, a "
            "filamentous ascomycete against a yeast",
    109871: "Batrachochytrium dendrobatidis — a chytrid. Chytridiomycota is "
            "a phylum S20 *does* find the family in (9 records, 6 taxa), but "
            "this species is not one of them, so it tests a species-level "
            "absence inside a positive phylum — the patchy-loss pattern S20a "
            "found across the fungi",
    6035: "Encephalitozoon cuniculi — a microsporidian; Microsporidia 0/29, "
          "and the most reduced eukaryotic genome in the scope",
    5811: "Toxoplasma gondii — an apicomplexan; Apicomplexa 0/60",
}


# ------------------------------------------------------------------- tables

def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_presence() -> list[dict]:
    """S20's per-proteome presence table, eukaryotes only."""
    rows = read_tsv(S20_DIR / "proteome_presence.tsv")
    return [r for r in rows if r["group"] in EUKARYOTE_GROUPS]


def load_taxonomy() -> dict[str, dict]:
    """taxid -> the ranked lineage S20 wrote from NCBI's taxonomy API."""
    return {r["taxid"]: r for r in read_tsv(S20_DIR / "taxonomy.tsv")}


# --------------------------------------------------------------- derivation

def clade_counts(presence: list[dict], taxonomy: dict[str, dict],
                 rank: str) -> dict[str, tuple[int, int]]:
    """clade name -> (proteomes swept, proteomes with an ITPR call)."""
    agg: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for r in presence:
        name = (taxonomy.get(r["taxid"], {}).get(rank) or "").strip()
        if not name:
            continue
        agg[name][0] += 1
        if r["itpr_status"] == "present":
            agg[name][1] += 1
    return {k: (v[0], v[1]) for k, v in agg.items()}


def big_phyla(presence: list[dict], taxonomy: dict[str, dict]) -> dict[str, int]:
    """G2 — phyla swept at least `BIG_PHYLUM_SWEPT` proteomes deep."""
    counts = clade_counts(presence, taxonomy, "phylum")
    return {k: n for k, (n, _) in counts.items() if n >= BIG_PHYLUM_SWEPT}


def absence_clades(presence: list[dict],
                   taxonomy: dict[str, dict]) -> list[dict]:
    """G3 — every clade swept deep enough and returning zero ITPR calls.

    Reported at both phylum and class rank, because the two carry different
    claims: "Ascomycota 0/1,034" is the headline and "Saccharomycetes 0/61" is
    what makes it a statement about independent lineages rather than about one
    badly-annotated group. A class inside a zero-ITPR phylum appears in both,
    and the manifest merges the reasons onto one genome.
    """
    out = []
    for rank in ("phylum", "class"):
        for name, (swept, present) in sorted(
                clade_counts(presence, taxonomy, rank).items(),
                key=lambda kv: -kv[1][0]):
            if present == 0 and swept >= MIN_SWEPT_FOR_CLAIM:
                out.append({"rank": rank, "clade": name, "swept": swept,
                            "with_itpr": present})
    return out


def copy_number_questions(presence: list[dict],
                          taxonomy: dict[str, dict] | None = None,
                          threshold: int = COPY_NUMBER_QUESTION) -> list[dict]:
    """G5 — the highest-copy-number species per class, above `threshold`.

    Without `taxonomy` the rule is unstratified and returns every species over
    the threshold; the manifest always passes it, and the stratified form is
    what the scope means.
    """
    over = []
    for r in presence:
        try:
            n = int(r["n_itpr"])
        except (TypeError, ValueError):
            continue
        if n >= threshold:
            over.append({"taxid": r["taxid"], "organism": r["organism"],
                         "n_itpr": n, "upid": r["upid"], "group": r["group"],
                         "protein_count": r["protein_count"], "clade": ""})
    over.sort(key=lambda r: -r["n_itpr"])
    if taxonomy is None:
        return over
    best: dict[str, dict] = {}
    for r in over:
        t = taxonomy.get(r["taxid"], {})
        clade = next((t.get(rank) for rank in COPY_NUMBER_RANK if t.get(rank)),
                     "") or "unclassified"
        r["clade"] = clade
        if clade not in best:          # `over` is already sorted by count
            best[clade] = r
    return sorted(best.values(), key=lambda r: -r["n_itpr"])


def control_clades(presence: list[dict], taxonomy: dict[str, dict]
                   ) -> list[dict]:
    """The clades a positive-control bait has to be drawn from (bait rule B2).

    Every G3 absence clade at **class** rank, plus each zero-ITPR phylum that
    has no qualifying class. A control bait proves the search works in the
    genomes the negative claim is about, so it has to come from those genomes'
    own lineage — a chlorophyte mannosyltransferase is not a control for
    *Arabidopsis*.
    """
    clades = absence_clades(presence, taxonomy)
    by_rank = defaultdict(list)
    for c in clades:
        by_rank[c["rank"]].append(c)
    have_class = set()
    for c in by_rank["class"]:
        have_class.add(c["clade"])
    # a phylum keeps its own row only if none of its classes qualified
    phylum_of_class: dict[str, str] = {}
    for r in presence:
        t = taxonomy.get(r["taxid"], {})
        if t.get("class") and t.get("phylum"):
            phylum_of_class[t["class"]] = t["phylum"]
    covered_phyla = {phylum_of_class.get(c) for c in have_class}
    out = list(by_rank["class"])
    out += [c for c in by_rank["phylum"] if c["clade"] not in covered_phyla]
    return out
