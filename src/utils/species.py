"""Species name / taxon-id helpers.

The lookup table covers (a) the common model organisms, (b) the lineages
that decide the ITPR family's open questions: the teleosts that carry the
3R duplicates, and the non-metazoan eukaryotes where the family is reported
present (amoebozoa, ciliates, kinetoplastids) or absent (land plants,
fungi). It is intentionally small — the GUI lets the user type any species
name as a free-text override, and the census tasks work from taxon-id
manifests rather than from this table.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SpeciesInfo:
    scientific: str   # e.g. "Homo sapiens"
    common: str       # e.g. "Human"
    taxon_id: int     # NCBI taxonomy id
    ensembl_slug: str # e.g. "homo_sapiens" (path component for Ensembl REST)


_SPECIES: list[SpeciesInfo] = [
    SpeciesInfo("Homo sapiens",            "Human",          9606,  "homo_sapiens"),
    SpeciesInfo("Mus musculus",            "Mouse",          10090, "mus_musculus"),
    SpeciesInfo("Rattus norvegicus",       "Rat",            10116, "rattus_norvegicus"),
    SpeciesInfo("Danio rerio",             "Zebrafish",      7955,  "danio_rerio"),
    SpeciesInfo("Gallus gallus",           "Chicken",        9031,  "gallus_gallus"),
    SpeciesInfo("Xenopus tropicalis",      "Frog",           8364,  "xenopus_tropicalis"),
    SpeciesInfo("Drosophila melanogaster", "Fruit fly",      7227,  "drosophila_melanogaster"),
    SpeciesInfo("Caenorhabditis elegans",  "C. elegans",     6239,  "caenorhabditis_elegans"),
    SpeciesInfo("Saccharomyces cerevisiae","Yeast",          4932,  "saccharomyces_cerevisiae"),
    SpeciesInfo("Ornithorhynchus anatinus","Platypus",       9258,  "ornithorhynchus_anatinus"),
    SpeciesInfo("Anolis carolinensis",     "Anole lizard",   28377, "anolis_carolinensis"),
    SpeciesInfo("Takifugu rubripes",       "Fugu",           31033, "takifugu_rubripes"),
    SpeciesInfo("Callorhinchus milii",     "Elephant shark", 7868,  "callorhinchus_milii"),
    SpeciesInfo("Petromyzon marinus",      "Sea lamprey",    7757,  "petromyzon_marinus"),
    SpeciesInfo("Branchiostoma floridae",  "Amphioxus",      7739,  "branchiostoma_floridae"),
    # The margins that decide the family's range (S3/S20-class tasks):
    SpeciesInfo("Dictyostelium discoideum","Slime mould",    44689, "dictyostelium_discoideum"),
    SpeciesInfo("Trypanosoma brucei",      "Trypanosome",    5691,  "trypanosoma_brucei"),
    SpeciesInfo("Paramecium tetraurelia",  "Paramecium",     5888,  "paramecium_tetraurelia"),
    SpeciesInfo("Tetrahymena thermophila", "Tetrahymena",    5911,  "tetrahymena_thermophila"),
    SpeciesInfo("Capsaspora owczarzaki",   "Capsaspora",     595528,"capsaspora_owczarzaki"),
    SpeciesInfo("Monosiga brevicollis",    "Choanoflagellate",81824,"monosiga_brevicollis"),
    SpeciesInfo("Arabidopsis thaliana",    "Thale cress",    3702,  "arabidopsis_thaliana"),
    SpeciesInfo("Chlamydomonas reinhardtii","Chlamydomonas", 3055,  "chlamydomonas_reinhardtii"),
]


SPECIES_LOOKUP: dict[str, SpeciesInfo] = {}
for _sp in _SPECIES:
    SPECIES_LOOKUP[_sp.scientific.lower()] = _sp
    SPECIES_LOOKUP[_sp.common.lower()] = _sp
    SPECIES_LOOKUP[_sp.ensembl_slug] = _sp


DEFAULT_SPECIES_PANEL: list[str] = [
    "Homo sapiens",
    "Mus musculus",
    "Danio rerio",
    "Gallus gallus",
    "Xenopus tropicalis",
    "Ornithorhynchus anatinus",
    "Anolis carolinensis",
    "Callorhinchus milii",
]


def resolve_species(name: str) -> SpeciesInfo | None:
    """Resolve a free-text species name (common or scientific) to SpeciesInfo.

    Returns None if not found — caller decides whether to fall back to
    free-text query or skip.
    """
    if not name:
        return None
    return SPECIES_LOOKUP.get(name.strip().lower())


def ensembl_species_slug(name: str) -> str:
    """Convert any input species string to Ensembl's URL slug form."""
    info = resolve_species(name)
    if info:
        return info.ensembl_slug
    # fallback: lowercase + underscores; Ensembl uses this convention
    return name.strip().lower().replace(" ", "_")


def all_species() -> list[SpeciesInfo]:
    return list(_SPECIES)
