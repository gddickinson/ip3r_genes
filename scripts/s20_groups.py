"""S20 — the non-vertebrate sweep's scope: which proteomes, and why.

S3 swept `vertebrata` (763 reference proteomes) and found where the three
paralogs live. S20 asks the opposite question — the family's *range* — so
its scope has to cover everything vertebrates are not, and the negative
claims have to be made against a stated sample rather than against silence.

Six groups. The four eukaryote ones partition Eukaryota with S3's
`vertebrata`, so no proteome is ever swept twice and the two sweeps'
denominators add:

    vertebrata          (S3)            7742
    metazoa_nonvert     33208 − 7742
    fungi               4751
    viridiplantae       33090
    protista_other      2759 − (33208 ∪ 33090 ∪ 4751)

and two prokaryote ones carry the negative claim:

    archaea             2157            complete — 641 proteomes is tractable
    bacteria_genus      2               genus-stratified, see SAMPLE_RULE

**Why bacteria are sampled and archaea are not.** UniProt lists 17,981
bacterial reference proteomes (73.8 M proteins); sweeping them all would
cost more than the rest of the project's compute put together to answer a
question no one disputes. Archaea's 641 (1.8 M proteins) cost nothing, so
they are taken whole and no sampling caveat attaches to them.

**The bacterial sampling rule is deliberately generous to the hypothesis it
tests.** One proteome per genus, and within a genus the one with the *most*
proteins. A larger proteome is a more sensitive place to find a homolog, so
picking the largest makes a negative result harder to obtain — the opposite
of the choice that would flatter it. Stratifying by genus rather than
sampling at random keeps the sample from being 400 *Streptococcus* strains,
which is what UniProt's reference set looks like unstratified.

**Genus is taken as the first whitespace token of the organism name.** For
binomial names that is the genus; for the strain-heavy names in this set it
is still the genus, because UniProt writes them "Genus species (strain …)".
It is a stratification key, not a taxonomic claim, and the sample table
records the key beside the organism so the choice is auditable.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# --------------------------------------------------------------- the groups
# `query`  : the UniProt proteomes query (reference proteomes only)
# `domain` : the reference_proteomes/ FTP directory
# `sample` : None = take every proteome; otherwise the stratification rule
# `relaxed`: run the E<=10 relaxed panel over this group (brief step 5)
GROUPS: dict[str, dict] = {
    "metazoa_nonvert": {
        "query": "taxonomy_id:33208+NOT+taxonomy_id:7742+AND+reference:true",
        "domain": "Eukaryota",
        "sample": None,
        "relaxed": False,
        "note": "invertebrates — the single-Itpr grade",
    },
    "fungi": {
        "query": "taxonomy_id:4751+AND+reference:true",
        "domain": "Eukaryota",
        "sample": None,
        "relaxed": True,
        "note": "S2 called only early-diverging phyla; Dikarya contributed "
                "no records at all — the claim S20 has to test",
    },
    "viridiplantae": {
        "query": "taxonomy_id:33090+AND+reference:true",
        "domain": "Eukaryota",
        "sample": None,
        "relaxed": True,
        "note": "S2: every call was Chlorophyta, Streptophyta 0/15 records",
    },
    "protista_other": {
        "query": ("taxonomy_id:2759+NOT+taxonomy_id:33208"
                  "+NOT+taxonomy_id:33090+NOT+taxonomy_id:4751"
                  "+AND+reference:true"),
        "domain": "Eukaryota",
        "sample": None,
        "relaxed": False,
        "note": "amoebozoa, SAR, discoba — where iplA and the deep grade sit",
    },
    "archaea": {
        "query": "taxonomy_id:2157+AND+reference:true",
        "domain": "Archaea",
        "sample": None,
        "relaxed": True,
        "note": "complete, not sampled — 641 proteomes cost nothing",
    },
    "bacteria_genus": {
        "query": "taxonomy_id:2+AND+reference:true",
        "domain": "Bacteria",
        "sample": "largest_per_genus",
        "relaxed": True,
        "note": "genus-stratified sample of 17,981 reference proteomes",
    },
}

EUKARYOTE_GROUPS = ("metazoa_nonvert", "fungi", "viridiplantae",
                    "protista_other")
PROKARYOTE_GROUPS = ("archaea", "bacteria_genus")
ALL_GROUPS = EUKARYOTE_GROUPS + PROKARYOTE_GROUPS

# The sensitivity at which every claim in this task is made.
EVALUE_PRIMARY = "1e-5"      # S3's protocol, so the two sweeps compare
EVALUE_RELAXED = "10"        # brief step 5 — "not found at E<=10"

SAMPLE_RULE = {
    "largest_per_genus": (
        "one proteome per genus (first token of the organism name), the one "
        "with the most proteins — the most sensitive member of each genus, "
        "so the negative claim is made in the places most likely to break it"
    ),
}


def genus_key(organism: str) -> str:
    """The stratification key: first whitespace token of the organism name."""
    return organism.split()[0] if organism.split() else organism


def largest_per_genus(rows: list[dict]) -> tuple[list[dict], dict]:
    """Apply SAMPLE_RULE['largest_per_genus'].

    Returns (selected rows, stats). Every selected row gains `genus` and
    `genus_n_available`, so the sample table records what each pick stood
    for rather than only what was picked.
    """
    by_genus: dict[str, list[dict]] = {}
    for row in rows:
        by_genus.setdefault(genus_key(row["organism"]), []).append(row)
    selected = []
    for genus, members in sorted(by_genus.items()):
        best = max(members, key=lambda r: (r["protein_count"], r["upid"]))
        best = dict(best, genus=genus, genus_n_available=len(members))
        selected.append(best)
    stats = {
        "rule": "largest_per_genus",
        "rule_text": SAMPLE_RULE["largest_per_genus"],
        "n_available": len(rows),
        "n_genera": len(by_genus),
        "n_selected": len(selected),
        "proteins_available": sum(r["protein_count"] for r in rows),
        "proteins_selected": sum(r["protein_count"] for r in selected),
        "largest_genus": max(
            ((g, len(m)) for g, m in by_genus.items()),
            key=lambda gm: gm[1], default=("", 0)),
    }
    return selected, stats


SAMPLERS = {"largest_per_genus": largest_per_genus}


def apply_sample(group: str, rows: list[dict]) -> tuple[list[dict], dict]:
    """Sample a group's proteome list by its declared rule (or take it whole)."""
    rule = GROUPS[group].get("sample")
    if rule is None:
        return rows, {"rule": "none", "rule_text": "every reference proteome",
                      "n_available": len(rows), "n_selected": len(rows),
                      "proteins_available": sum(r["protein_count"] for r in rows),
                      "proteins_selected": sum(r["protein_count"] for r in rows)}
    return SAMPLERS[rule](rows)
