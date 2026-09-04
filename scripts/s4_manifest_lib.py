"""s4_manifest_lib.py — parsing + selection for the S4 genome manifest.

Pure functions over NCBI `datasets` JSONL dumps plus this project's own
committed census tables. No network calls here — the driver
(`s4_build_manifest.py`) owns every `datasets` invocation and its caching.

**What is different from the PIEZO port.** There the margin species were a
hand-typed list of piezo3 cases. Here every margin species is *derived from
a committed table* by a rule with a stated threshold (D13 applied to scope
selection), so the denominator can be rebuilt from the census and cannot
drift from it:

  M1 `zero_hit_proteome`  a swept reference proteome with no ITPR record at
                          all (`results/census_v3/proteomes_without_hits.tsv`).
  M2 `missing_paralog`    fewer than 3 ITPR records — the vertebrate family
                          has three paralogs, so <3 is a candidate loss.
  M3 `fragment_only`      no ITPR record reaches the family's own length
                          floor (`family.MIN_LENGTH_AA`, 2,000 aa): the
                          species is represented only by broken gene models.
  M4 `anchor`             the reference species the whole project is
                          calibrated on.

M2 and M3 are *questions*, not findings. S3 showed copy number tracks
annotation depth (bird proteomes hold a third the proteins of mammalian
ones), so a low count is exactly what the genome sweep exists to resolve —
which is why these species are in the denominator rather than in a results
table.
"""

from __future__ import annotations

import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path

#: S1's toolchain manifest records `datasets` as env-resident, not on the
#: bare PATH, so every caller resolves it instead of assuming it.
ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")

# Assembly-level rank for ranking (higher = better).
LEVEL_RANK = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}

# Informal clade labels for orders NCBI leaves without a ranked class
# (turtles/crocodylians sit under unranked Archelosauria/Archosauria;
# lungfish and coelacanths under unranked Dipnomorpha/Coelacanthimorpha).
CLASS_FALLBACK = {
    "Ceratodontiformes": "Dipnoi",
    "Coelacanthiformes": "Coelacanthimorpha",
    "Crocodylia": "Crocodylia",
    "Testudines": "Testudines",
}

# Compression factor: an NCBI datasets zip of a genome FASTA is ~30 % of the
# uncompressed sequence length (empirical for vertebrate DNA).
ZIP_FACTOR = 0.30

# The vertebrate paralog count. Not a guess: `family.REFERENCE_ACCESSIONS`
# is human ITPR1/2/3, so three is what a completely annotated vertebrate
# proteome should carry.
EXPECTED_PARALOGS = 3

#: Reference species the project is calibrated on — the S0 background, the
#: S1 control panel and the S3 seeds all run through these.
ANCHORS = {
    9606: "human — the reference panel (ITPR1/2/3) and every structure claim",
    10090: "mouse — the genetic models the review's §9 rests on",
    9031: "chicken — the bird reference, and Aves is where the margin sits",
    7955: "zebrafish — the teleost 3R case and the S0 PF08709 query",
    8364: "Xenopus tropicalis — the tetrapod outgroup to amniotes",
}


def datasets_bin() -> str:
    """Absolute path to the NCBI `datasets` CLI, PATH first then the env."""
    found = shutil.which("datasets")
    if found:
        return found
    env = ENV_BIN / "datasets"
    if env.exists():
        return str(env)
    raise SystemExit("datasets CLI not found on PATH or in the piezo1 env")


def class_label(tax: dict, order: str) -> str:
    return tax.get("class") or CLASS_FALLBACK.get(order, "")


def load_jsonl(path: Path) -> list[dict]:
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_assembly(rec: dict) -> dict:
    """Flatten one `datasets summary genome` record to the fields we keep."""
    info = rec.get("assembly_info", {})
    stats = rec.get("assembly_stats", {})
    org = rec.get("organism", {})
    return {
        "accession": rec.get("accession", ""),
        "organism": org.get("organism_name", ""),
        "taxid": org.get("tax_id"),
        "common_name": org.get("common_name", ""),
        "assembly_name": info.get("assembly_name", ""),
        "level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "release_date": info.get("release_date", ""),
        "total_length": int(stats.get("total_sequence_length", 0) or 0),
        "scaffold_n50": int(stats.get("scaffold_n50", 0) or 0),
        "contig_n50": int(stats.get("contig_n50", 0) or 0),
        "annotated": "annotation_info" in rec,
        "source_db": rec.get("source_database", ""),
    }


def load_taxonomy(path: Path) -> dict[int, dict]:
    """taxid -> {'order': name, 'class': name} from a taxonomy JSONL dump."""
    taxmap: dict[int, dict] = {}
    if not path.exists():
        return taxmap
    for rec in load_jsonl(path):
        wrapped = rec.get("taxonomy", rec)
        cls = wrapped.get("classification", {})
        entry = {
            "order": (cls.get("order") or {}).get("name", ""),
            "class": (cls.get("class") or {}).get("name", ""),
        }
        for q in rec.get("query") or []:
            try:
                taxmap[int(q)] = entry
            except (TypeError, ValueError):
                pass
        tid = wrapped.get("tax_id")
        if tid:
            taxmap[int(tid)] = entry
    return taxmap


# ------------------------------------------------------- margin derivation
def derive_margins(project_root: Path, min_length_aa: int) -> list[dict]:
    """The margin species, computed from the committed S3 tables.

    Returns one row per species with every rule that fired on it, so a
    species that is both fragment-only and short of paralogs is one manifest
    row carrying both reasons rather than two rows.
    """
    v3 = project_root / "results" / "census_v3"
    census = read_tsv(v3 / "census_v3.tsv")
    swept = {int(m["Organism Id"]): m
             for m in read_tsv(project_root / "results" / "hmm_sweep"
                               / "proteome_manifest.tsv")}

    by_taxon: dict[int, list[dict]] = defaultdict(list)
    for r in census:
        if r["call"] != "ITPR" or not r["taxon_id"]:
            continue
        tid = int(r["taxon_id"])
        if tid in swept:
            by_taxon[tid].append(r)

    def longest(recs: list[dict]) -> int:
        out = 0
        for r in recs:
            try:
                out = max(out, int(r["length"]))
            except (TypeError, ValueError):
                pass
        return out

    hits: dict[int, dict] = {}

    def add(tid: int, name: str, reason: str, note: str) -> None:
        row = hits.setdefault(tid, {"taxid": tid, "organism": name,
                                    "reasons": [], "notes": []})
        row["reasons"].append(reason)
        row["notes"].append(note)

    for row in read_tsv(v3 / "proteomes_without_hits.tsv"):
        tid = int(row["Organism Id"])
        add(tid, row["Organism"].split("(")[0].strip(), "zero_hit_proteome",
            f"0 ITPR records in {row['Proteome Id']} "
            f"({int(row['Protein count']):,} proteins); absence in a proteome "
            "is not gene loss")

    for tid, recs in by_taxon.items():
        name = recs[0]["species"] or swept[tid]["Organism"]
        n_prot = int(swept[tid].get("Protein count") or 0)
        if len(recs) < EXPECTED_PARALOGS:
            add(tid, name, "missing_paralog",
                f"{len(recs)} of {EXPECTED_PARALOGS} expected paralogs "
                f"({n_prot:,} proteins in the reference proteome)")
        if longest(recs) < min_length_aa:
            add(tid, name, "fragment_only",
                f"longest ITPR record {longest(recs)} aa, below the family "
                f"floor of {min_length_aa:,} aa ({len(recs)} records)")

    for tid, why in ANCHORS.items():
        name = (by_taxon.get(tid, [{}])[0].get("species")
                or (swept.get(tid) or {}).get("Organism", "")
                or why.split("—")[0].strip())
        add(tid, name, "anchor", why)

    return sorted(hits.values(), key=lambda r: r["organism"])


# ---------------------------------------------------------------- ranking
def rank_key(asm: dict) -> tuple:
    """Sort key: annotated > RefSeq > assembly level > scaffold N50.

    D9: a RefSeq (`GCF_`) gene set and a submitter GenBank (`GCA_`) gene set
    are not comparable evidence, so annotation source is part of the choice
    and is carried into the manifest rather than discarded here.
    """
    return (
        1 if asm["annotated"] else 0,
        1 if asm["accession"].startswith("GCF_") else 0,
        LEVEL_RANK.get(asm["level"], 0),
        asm["scaffold_n50"],
    )


def select_order_reps(assemblies: list[dict],
                      taxmap: dict[int, dict]) -> tuple[dict, list]:
    """One best assembly per order. Returns (order -> asm, no_order list)."""
    by_order: dict[str, list[dict]] = defaultdict(list)
    no_order = []
    for asm in assemblies:
        tax = taxmap.get(asm["taxid"], {})
        order = tax.get("order", "")
        if not order:
            no_order.append(dict(asm, vclass=tax.get("class", "")))
            continue
        by_order[order].append(
            dict(asm, vorder=order, vclass=class_label(tax, order)))
    reps = {}
    for order, group in sorted(by_order.items()):
        best = max(group, key=rank_key)
        best["order_pool_size"] = len(group)
        reps[order] = best
    return reps, no_order


def merge_rows(order_reps: dict, margin_rows: list[dict]) -> list[dict]:
    """Union of order reps + margin assemblies, merging reasons on collision.

    A species can be its order's representative *and* a margin case; that is
    one genome, one row, two reasons — the union is over assemblies, so the
    manifest's row count is the number of genomes actually downloaded.
    """
    by_acc: dict[str, dict] = {}
    for order, asm in order_reps.items():
        row = dict(asm)
        row["reasons"] = [f"order_rep:{order}"]
        row["notes"] = [f"best of {asm.get('order_pool_size', 1)} reference "
                        f"assemblies in {order}"]
        by_acc[row["accession"]] = row
    for m in margin_rows:
        acc = m["accession"]
        if acc in by_acc:
            by_acc[acc]["reasons"].extend(m["reasons"])
            by_acc[acc]["notes"].extend(m["notes"])
        else:
            by_acc[acc] = dict(m)
    return sorted(by_acc.values(),
                  key=lambda r: (r.get("vclass", ""), r.get("vorder", ""),
                                 r["organism"]))


MANIFEST_COLUMNS = [
    "accession", "organism", "taxid", "vclass", "vorder", "assembly_name",
    "level", "refseq_category", "annotated", "annotation_source",
    "release_date", "total_length_bp", "scaffold_n50", "contig_n50",
    "fasta_gb", "est_zip_gb", "reasons", "notes",
]


def write_manifest(rows: list[dict], out_tsv: Path) -> dict:
    """Write the TSV; return the size totals the storage check needs."""
    tot_bp = tot_fasta = tot_zip = 0.0
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(MANIFEST_COLUMNS)
        for r in rows:
            fasta_gb = r["total_length"] / 1e9
            zip_gb = fasta_gb * ZIP_FACTOR
            tot_bp += r["total_length"]
            tot_fasta += fasta_gb
            tot_zip += zip_gb
            w.writerow([
                r["accession"], r["organism"], r["taxid"],
                r.get("vclass", ""), r.get("vorder", ""), r["assembly_name"],
                r["level"], r["refseq_category"],
                "Y" if r["annotated"] else "N",
                "RefSeq" if r["accession"].startswith("GCF_") else "GenBank",
                r["release_date"], r["total_length"], r["scaffold_n50"],
                r["contig_n50"], f"{fasta_gb:.2f}", f"{zip_gb:.2f}",
                ";".join(r["reasons"]), " | ".join(r.get("notes", [])),
            ])
    return {"n": len(rows), "total_bp": tot_bp,
            "total_fasta_gb": tot_fasta, "total_zip_gb": tot_zip}
