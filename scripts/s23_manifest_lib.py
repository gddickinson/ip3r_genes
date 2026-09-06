"""s23_manifest_lib.py — selection and writing for the S23 genome manifest.

Pure functions over NCBI `datasets` dumps plus this project's committed tables.
No network calls: the driver (`s23_build_manifest.py`) owns every `datasets`
invocation and its caching, exactly as S4 does.

The parsing and ranking primitives are **imported from `s4_manifest_lib`**
rather than copied — `parse_assembly`, `load_jsonl`, `load_taxonomy`,
`datasets_bin` and `LEVEL_RANK` carry no vertebrate assumptions, and two
copies of an assembly ranker is how two manifests come to disagree about which
assembly is best.

What is new here is the **size guard**, which S4 did not need. S4's vertebrate
genomes are all large and all recent; this scope reaches assemblies of 2.5 Mbp
(a microsporidian) alongside 20 Gbp plants, and a fragmentary assembly cannot
answer an absence question. The guard is D4 applied at selection time rather
than only at reporting time:

  * `rank_key_s23` prefers, in order: an assembly whose contigs can hold the
    group's gene (`s23_calibration.spans_a_gene`), then an annotated one, then
    RefSeq, then assembly level, then scaffold N50. The contiguity term goes
    **first** because it is the only one that decides whether the genome can
    answer the question at all — an annotated assembly that cannot carry the
    gene is worse evidence for an absence than an unannotated one that can.
  * a clade whose best assembly is still below the bar is **kept and flagged**,
    never dropped. Dropping it would silently shrink the denominator, which is
    the failure mode every absence claim in this project is built to avoid.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s23_calibration as cal                                  # noqa: E402
from s4_manifest_lib import (LEVEL_RANK, class_label,          # noqa: F401,E402
                             datasets_bin, load_jsonl,
                             load_taxonomy, parse_assembly)

#: An NCBI datasets zip of a genome FASTA is ~30 % of the sequence length.
#: S4's measured factor for vertebrate DNA; reused, and the manifest reports
#: both so a wrong factor is visible rather than load-bearing.
ZIP_FACTOR = 0.30

#: Kingdom-level group for the calibration lookup, from the taxonomy dump's
#: own ranks. These are the groups `s23_span_calibration.py` measured, which
#: are `s23_bait_spec.BANDS`' values — one vocabulary for baits, calibration
#: and manifest.
GROUP_OF_KINGDOM = {
    "Metazoa": "metazoa",
    "Fungi": "fungi",
    "Viridiplantae": "viridiplantae",
}

GROUP_OF_PHYLUM = {
    "Ciliophora": "sar", "Oomycota": "sar", "Perkinsozoa": "sar",
    "Apicomplexa": "sar", "Bacillariophyta": "sar",
    "Evosea": "amoebozoa", "Discosea": "amoebozoa", "Tubulinea": "amoebozoa",
    "Euglenozoa": "discoba", "Heterolobosea": "discoba",
}


def group_of(tax: dict) -> str:
    """The calibration group of one taxonomy record."""
    king = (tax.get("kingdom") or "").strip()
    if king in GROUP_OF_KINGDOM:
        return GROUP_OF_KINGDOM[king]
    phy = (tax.get("phylum") or "").strip()
    if phy in GROUP_OF_PHYLUM:
        return GROUP_OF_PHYLUM[phy]
    return "other"


#: NCBI taxid for Vertebrata — S4's scope, and this one's complement.
VERTEBRATA_TAXID = 7742


def is_vertebrate(tax: dict) -> bool:
    """S4 owns the vertebrates; this scope is everything else (bait rule B7).

    Decided on the **taxid lineage**, not on a class-name list: NCBI leaves
    several deep chordate lineages without a ranked class (S4 needed a
    `CLASS_FALLBACK` table for exactly this), and a name list would have let
    them through into a sweep whose bait panel has no vertebrate baits.
    """
    return VERTEBRATA_TAXID in (tax.get("parents") or ())


# ---------------------------------------------------------------- ranking

def rank_key_s23(asm: dict, group: str) -> tuple:
    """Best-assembly key for this scope. Contiguity first — see module doc."""
    holds, _bar, _why = cal.spans_a_gene(asm.get("contig_n50", 0), group)
    return (
        1 if holds else 0,
        1 if asm["annotated"] else 0,
        1 if asm["accession"].startswith("GCF_") else 0,
        LEVEL_RANK.get(asm["level"], 0),
        asm["scaffold_n50"],
    )


def best_of(assemblies: list[dict], group: str) -> dict | None:
    return max(assemblies, key=lambda a: rank_key_s23(a, group)) if assemblies \
        else None


# --------------------------------------------------------------- selection

def index_assemblies(assemblies: list[dict], taxmap: dict[int, dict]
                     ) -> tuple[dict, dict, dict, list]:
    """(by rank-name, by taxid, by ancestor taxid, skipped).

    One pass that builds every index the selection rules need, so the taxonomy
    join happens once. Vertebrates are dropped here and counted, because "how
    many assemblies did S4 already own" is a number the notes should state
    rather than leave to inference.

    **`by_ancestor` is not a convenience.** NCBI files a reference assembly
    under the taxid of the *strain* it was sequenced from, not the species:
    *Toxoplasma gondii* is taxid 5811 and its reference genome is *T. gondii*
    ME49, 508771. An anchor looked up by species taxid alone therefore misses,
    and on the first build it missed four of fourteen — *Dictyostelium*,
    *Toxoplasma*, *Encephalitozoon* and *Batrachochytrium*, i.e. the
    apicomplexan, the microsporidian and the amoebozoan the negative claims
    are named after. Indexing every ancestor taxid resolves them exactly,
    without name matching.
    """
    by_clade: dict[tuple, list[dict]] = defaultdict(list)
    by_taxid: dict[int, list[dict]] = defaultdict(list)
    by_ancestor: dict[int, list[dict]] = defaultdict(list)
    skipped = []
    for a in assemblies:
        tax = taxmap.get(a["taxid"])
        if not tax:
            skipped.append(dict(a, why="no taxonomy record"))
            continue
        if is_vertebrate(tax):
            skipped.append(dict(a, why="vertebrate — S4's scope"))
            continue
        grp = group_of(tax)
        row = dict(a, group=grp, kingdom=tax.get("kingdom", ""),
                   phylum=tax.get("phylum", ""), tclass=tax.get("class", ""),
                   torder=tax.get("order", ""))
        for rank, key in (("phylum", "phylum"), ("class", "tclass")):
            if row[key]:
                by_clade[(rank, row[key])].append(row)
        by_taxid[a["taxid"]].append(row)
        for anc in (tax.get("parents") or ()):
            by_ancestor[anc].append(row)
    return by_clade, by_taxid, by_ancestor, skipped


def select(by_clade: dict, by_taxid: dict, by_ancestor: dict,
           rules: dict) -> tuple[dict, list]:
    """Apply G1-G5. Returns (accession -> row with merged reasons, unresolved).

    A genome several rules reach is one row carrying all of them (S4's
    `merge_rows` pattern), so the row count is the number of genomes to
    download, not the number of reasons to want one.
    """
    picked: dict[str, dict] = {}
    unresolved: list[dict] = []

    def take(asm: dict | None, reason: str, note: str, where: str) -> None:
        if asm is None:
            unresolved.append({"rule": reason, "target": where, "note": note})
            return
        row = picked.setdefault(asm["accession"],
                                dict(asm, reasons=[], notes=[]))
        row["reasons"].append(reason)
        row["notes"].append(note)

    for phylum in sorted(rules["phyla"]):                       # G1
        pool = by_clade.get(("phylum", phylum), [])
        grp = pool[0]["group"] if pool else "other"
        take(best_of(pool, grp), "phylum_rep", phylum,
             f"best of {len(pool)} reference assemblies in {phylum}")

    for phylum, klass in sorted(rules["big_phylum_classes"]):   # G2
        pool = by_clade.get(("class", klass), [])
        grp = pool[0]["group"] if pool else "other"
        take(best_of(pool, grp), "class_rep", klass,
             f"class of {phylum}, a phylum S20 swept "
             f"{rules['big_phyla'][phylum]} proteomes of")

    for c in rules["absence_clades"]:                           # G3
        pool = by_clade.get((c["rank"], c["clade"]), [])
        grp = pool[0]["group"] if pool else "other"
        take(best_of(pool, grp), "absence_clade", c["clade"],
             f"S20 swept {c['swept']} proteomes of {c['clade']} and found "
             f"0 ITPR — the absence this genome tests")

    for taxid, why in rules["anchors"].items():                 # G4
        pool = (by_taxid.get(int(taxid))
                or by_ancestor.get(int(taxid), []))
        grp = pool[0]["group"] if pool else "other"
        take(best_of(pool, grp), "anchor", str(taxid), why)

    for r in rules["copy_number"]:                              # G5
        pool = (by_taxid.get(int(r["taxid"]))
                or by_ancestor.get(int(r["taxid"]), []))
        grp = pool[0]["group"] if pool else "other"
        take(best_of(pool, grp), "copy_number", r["organism"],
             f"{r['n_itpr']} ITPR records in its reference proteome "
             f"({r['clade']}); more copies than any vertebrate has paralogs")

    return picked, unresolved


# ------------------------------------------------------------------ output

MANIFEST_COLUMNS = [
    "accession", "organism", "taxid", "group", "kingdom", "phylum", "class",
    "order", "assembly_name", "level", "refseq_category", "annotated",
    "annotation_source", "release_date", "total_length_bp", "scaffold_n50",
    "contig_n50", "spans_gene", "contiguity_bar_bp", "bar_source",
    "fasta_gb", "est_zip_gb", "reasons", "notes",
]


def write_manifest(rows: list[dict], out_tsv: Path) -> dict:
    tot_bp = tot_fasta = tot_zip = 0.0
    below_bar = 0
    by_group: dict[str, int] = defaultdict(int)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(MANIFEST_COLUMNS)
        for r in sorted(rows, key=lambda r: (r.get("group", ""),
                                             r.get("phylum", ""),
                                             r["organism"])):
            holds, bar, why = cal.spans_a_gene(r.get("contig_n50", 0),
                                               r.get("group", "other"))
            below_bar += 0 if holds else 1
            by_group[r.get("group", "other")] += 1
            fasta_gb = r["total_length"] / 1e9
            zip_gb = fasta_gb * ZIP_FACTOR
            tot_bp += r["total_length"]
            tot_fasta += fasta_gb
            tot_zip += zip_gb
            w.writerow([
                r["accession"], r["organism"], r["taxid"], r.get("group", ""),
                r.get("kingdom", ""), r.get("phylum", ""), r.get("tclass", ""),
                r.get("torder", ""), r["assembly_name"], r["level"],
                r["refseq_category"], "Y" if r["annotated"] else "N",
                "RefSeq" if r["accession"].startswith("GCF_") else "GenBank",
                r["release_date"], r["total_length"], r["scaffold_n50"],
                r["contig_n50"], "Y" if holds else "N", bar, why,
                f"{fasta_gb:.3f}", f"{zip_gb:.3f}",
                ";".join(r["reasons"]), " | ".join(r.get("notes", [])),
            ])
    return {"n": len(rows), "total_bp": tot_bp,
            "total_fasta_gb": tot_fasta, "total_zip_gb": tot_zip,
            "below_contiguity_bar": below_bar, "by_group": dict(by_group)}
