"""S10 step 5 — the assembly-version audit, and what the databases actually hold.

The brief's first question is the one that would invalidate everything else:
*do the two annotations even sit on the same build?* A gene model that looks
missing because it was called on a superseded assembly is a bookkeeping
artefact, not an annotation bug, and it is the first thing a reviewer will ask.

Three independent records are pulled and compared:

* the **assembly report** `datasets` wrote beside the genome — name, level,
  release date, submitter, and whether NCBI still calls this assembly
  `current`;
* the **GFF3 header**, which NCBI's annotwriter stamps with the build the
  annotation was written against (`genome-build-accession`). If that accession
  is not the assembly the sweep searched, the case is void;
* the **assembly list for the species**, fetched live from `datasets` and
  archived, so "this is the only assembly" or "a newer one exists" is a
  statement about NCBI's holdings rather than an assumption. A newer assembly
  does not void a case — the audited annotation is still what the databases
  serve — but it changes what the paper should say about it, so it is measured
  rather than hoped.

Then two inventories that turn a coordinate-level finding into a claim about
databases: every **family-named model in the whole annotation** with the
biotype that decides whether it emits a protein, and every **protein record
the census holds for the species**. The second is the consequence of the
first, and it is the number a reader cares about — how many IP3-receptor
proteins this species has in the protein databases, and why.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s10_gff import annotation_meta                            # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

#: Which annotated gene names count as family names. Sourced from S5's own
#: list, which is itself sourced from `src/utils/family.py` (CLAUDE.md's rule:
#: the family is defined in one place), so S10 cannot recognise a different set
#: of names from the sweep it is auditing.
try:                                                  # pragma: no cover
    from s5_classify import ITPR_NAME_HINTS, RYR_NAME_HINTS
except Exception:                                     # pragma: no cover
    ITPR_NAME_HINTS = ("itpr", "ip3r", "insp3r", "itr-1")
    RYR_NAME_HINTS = ("ryr", "ryanodine receptor")


def _as_int(v) -> int | str:
    """`datasets` types its gene counts inconsistently across releases — some
    come back as JSON numbers and some as strings. Arithmetic on the mixture
    raises; returning "" for a genuinely absent count keeps the distinction
    between zero and unknown."""
    if v in (None, ""):
        return ""
    try:
        return int(v)
    except (TypeError, ValueError):
        return ""


def _pick(d: dict, *keys, default=""):
    for k in keys:
        if d.get(k) not in (None, ""):
            return d[k]
    return default


def _datasets_bin() -> str:
    from shutil import which
    hit = which("datasets")
    if hit:
        return hit
    env = Path("/opt/anaconda3/envs/piezo1/bin/datasets")
    if env.exists():
        return str(env)
    return ""


def assemblies_for_taxid(taxid: int, archive: Path) -> list[dict]:
    """Every assembly NCBI lists for the species, cached on first fetch.

    Cached by taxid under the data root, in the same shape S4 archives its
    `datasets` calls, so a re-run of S10 is offline. A failure to reach NCBI
    returns an empty list and the audit records that it could not be done,
    rather than silently reporting "one assembly".
    """
    archive.mkdir(parents=True, exist_ok=True)
    cache = archive / f"taxon_{taxid}_all.jsonl"
    if not cache.exists() or cache.stat().st_size == 0:
        exe = _datasets_bin()
        if not exe:
            return []
        try:
            res = subprocess.run(
                [exe, "summary", "genome", "taxon", str(taxid),
                 "--as-json-lines"], check=True, capture_output=True,
                text=True, timeout=180)
        except Exception:                             # noqa: BLE001
            return []
        cache.write_text(res.stdout)
    out = []
    for line in cache.read_text().splitlines():
        if not line.strip():
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        # `datasets summary genome taxon` emits snake_case while the
        # per-genome `assembly_data_report.jsonl` emits camelCase for the same
        # fields. Reading only one spelling silently returns blanks for every
        # assembly, which reads as "no metadata" rather than as a parse bug.
        info = d.get("assembly_info") or d.get("assemblyInfo") or {}
        ann = d.get("annotation_info") or {}
        stats = d.get("assembly_stats") or d.get("assemblyStats") or {}
        counts = ((ann.get("stats") or {}).get("gene_counts") or {})
        out.append({
            "accession": d.get("accession", ""),
            "assembly_name": _pick(info, "assembly_name", "assemblyName"),
            "level": _pick(info, "assembly_level", "assemblyLevel"),
            "status": _pick(info, "assembly_status", "assemblyStatus"),
            "refseq_category": _pick(info, "refseq_category", "refseqCategory"),
            "release_date": _pick(info, "release_date", "releaseDate"),
            "submitter": info.get("submitter", ""),
            "annotated": "Y" if ann else "N",
            "annotation_name": ann.get("name", ""),
            "annotation_provider": ann.get("provider", ""),
            "annotation_pipeline": ann.get("pipeline", ""),
            "annotation_release_date": ann.get("release_date", ""),
            "n_protein_coding": _as_int(counts.get("protein_coding")),
            "n_pseudogene": _as_int(counts.get("pseudogene")),
            "n_genes_total": _as_int(counts.get("total")),
            "contig_n50": _pick(stats, "contig_n50", "contigN50"),
        })
    out.sort(key=lambda r: (r["release_date"], r["accession"]))
    return out


def assembly_audit(case: dict, genome_dir: Path, archive: Path) -> dict:
    """One row: does the annotation sit on the assembly the sweep searched?"""
    report = genome_dir / "assembly_data_report.jsonl"
    info: dict = {}
    if report.exists():
        d = json.loads(report.read_text())
        ai = d.get("assemblyInfo") or {}
        st = d.get("assemblyStats") or {}
        info = {
            "assembly_name": ai.get("assemblyName", ""),
            "assembly_level": ai.get("assemblyLevel", ""),
            "assembly_status": ai.get("assemblyStatus", ""),
            "refseq_category": ai.get("refseqCategory", ""),
            "release_date": ai.get("releaseDate", ""),
            "submitter": ai.get("submitter", ""),
            "sequencing_tech": ai.get("sequencingTech", ""),
            "assembly_method": ai.get("assemblyMethod", ""),
            "annotation_comment": ai.get("comments", ""),
            "bioproject": ai.get("bioprojectAccession", ""),
            "taxid": ((ai.get("biosample") or {}).get("description", {})
                      .get("organism", {}).get("taxId", "")),
            "contig_n50": st.get("contigN50", ""),
            "total_bp": st.get("totalSequenceLength", ""),
        }
    meta = annotation_meta(genome_dir / "genomic.gff.gz")
    build_acc = meta.get("build_accession", "").replace("NCBI_Assembly:", "")
    same = (build_acc == case["accession"]) if build_acc else None
    others = assemblies_for_taxid(int(info.get("taxid") or 0), archive) \
        if info.get("taxid") else []
    mine = next((a for a in others if a["accession"] == case["accession"]), {})
    newer = [a for a in others
             if a["release_date"] > (info.get("release_date") or "")
             and a["accession"] != case["accession"]]
    return {
        "case_id": case.get("case_id", ""), "accession": case["accession"],
        "organism": case["organism"], **info,
        "gff_build_name": meta.get("build_name", ""),
        "gff_build_accession": build_acc,
        "gff_processor": meta.get("processor", ""),
        "annotation_on_this_assembly": "" if same is None else str(bool(same)),
        "n_assemblies_for_species": len(others),
        "n_newer_assemblies": len(newer),
        "newer_assemblies": ";".join(a["accession"] for a in newer),
        "newer_annotated": ";".join(a["accession"] for a in newer
                                    if a["annotated"] == "Y"),
        "assembly_list_available": str(bool(others)),
        # NCBI's own gene counts for this annotation. The pseudogene fraction
        # is the systemic number behind a case: a teleost annotation filing a
        # third of its genes as pseudogenes has a policy, not an accident.
        "annot_protein_coding": mine.get("n_protein_coding", ""),
        "annot_pseudogene": mine.get("n_pseudogene", ""),
        "annot_genes_total": mine.get("n_genes_total", ""),
        "annot_pseudogene_frac": (
            round(mine["n_pseudogene"] / mine["n_genes_total"], 4)
            if isinstance(mine.get("n_genes_total"), int)
            and isinstance(mine.get("n_pseudogene"), int)
            and mine["n_genes_total"] else ""),
        "annotation_release_date": mine.get("annotation_release_date", ""),
        "annotation_provider": mine.get("annotation_provider", ""),
    }


def family_named_models(genes_slim: Path) -> list[dict]:
    """Every family-named gene in the annotation, with its biotype.

    Genome-wide rather than locus-local, because the question this answers is a
    database question: how many of this species' IP3-receptor gene models are
    in a form that produces a protein at all. A model filed as `pseudogene`
    carries the name, occupies the locus, and emits nothing.
    """
    rows = []
    with open(genes_slim) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            name = f[5] or ""
            low = name.lower()
            fam = ""
            if any(h in low for h in ITPR_NAME_HINTS):
                fam = "ITPR"
            elif any(h in low for h in RYR_NAME_HINTS):
                fam = "RYR"
            if not fam:
                continue
            rows.append({
                "contig": f[0], "start": int(f[1]), "end": int(f[2]),
                "strand": f[3], "gene_id": f[4], "name": name,
                "biotype": f[6], "family_from_name": fam,
                "span": int(f[2]) - int(f[1]) + 1,
                "emits_protein": int(f[6] == "protein_coding"),
            })
    rows.sort(key=lambda r: (r["contig"], r["start"]))
    return rows


def database_records(species: str, census: Path) -> list[dict]:
    """Every protein record the census holds for this species.

    Read from the committed census rather than re-queried, so the count in the
    report is the same one every other task in the project is using (D13).
    """
    out = []
    with open(census) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("species") != species:
                continue
            out.append({
                "accession": r.get("accession", ""), "gene": r.get("gene", ""),
                "protein_name": (r.get("protein_name") or "")[:120],
                "length": r.get("length", ""), "call": r.get("call", ""),
                "confidence": r.get("confidence", ""),
                "fragment": r.get("fragment", ""),
                "reviewed": r.get("reviewed", ""),
                "source": r.get("source", ""),
            })
    out.sort(key=lambda r: (r["call"], r["gene"]))
    return out


def flank_consensus_check(case: dict, neighbours: list[dict],
                          consensus: Path,
                          key_kind: str = "relaxed",
                          window: str = "fixed10") -> list[dict]:
    """Are this locus's flanking genes the ones this paralog usually has?

    A fifth line of evidence, and an entirely independent one: S8 measured, per
    paralog, which flanking gene symbols recur across the swept vertebrates,
    each against its own matched random-window background. If the genes
    bracketing an *unannotated* locus are the paralog's own consensus
    neighbours, the gene is not merely present and intact — it is in the
    position that paralog occupies in the rest of the vertebrates, which no
    property of this assembly's annotation could produce.

    Read from S8's committed table, never recomputed (D13). A flanking gene
    the consensus does not list is reported with fraction 0 rather than
    dropped, and a genome whose annotation carries no gene symbols returns
    rows that all read 0 — which is the honest outcome there, not a negative
    result: the check cannot run on locus tags.
    """
    if not consensus.exists():
        return []
    rows = []
    with open(consensus) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("cell") == case["cell"] and r.get("key_kind") == key_kind
                    and r.get("window") == window):
                rows.append(r)
    by_symbol = {r["symbol"].upper(): r for r in rows}
    out = []
    for n in sorted(neighbours, key=lambda x: int(x["distance_bp"]))[:6]:
        hit = by_symbol.get((n["name"] or "").upper())
        out.append({
            "case_id": case["case_id"], "accession": case["accession"],
            "cell": case["cell"], "side": n["side"], "symbol": n["name"],
            "distance_bp": n["distance_bp"],
            "in_paralog_consensus": int(bool(hit)),
            "consensus_fraction": hit["fraction"] if hit else "0",
            "consensus_species": hit["n_species"] if hit else "0",
            "consensus_species_total": hit["n_species_total"] if hit else "0",
            "consensus_rank": (1 + sorted(
                (float(r["fraction"]) for r in rows), reverse=True).index(
                    float(hit["fraction"]))) if hit else "",
        })
    return out


FLANK_CHECK_COLS = ["case_id", "accession", "cell", "side", "symbol",
                    "distance_bp", "in_paralog_consensus",
                    "consensus_fraction", "consensus_species",
                    "consensus_species_total", "consensus_rank"]

AUDIT_COLS = ["case_id", "accession", "organism", "assembly_name",
              "assembly_level", "assembly_status", "refseq_category",
              "release_date", "submitter", "sequencing_tech",
              "assembly_method", "annotation_comment", "bioproject", "taxid",
              "contig_n50", "total_bp", "gff_build_name",
              "gff_build_accession", "gff_processor",
              "annotation_on_this_assembly", "n_assemblies_for_species",
              "n_newer_assemblies", "newer_assemblies", "newer_annotated",
              "assembly_list_available", "annotation_provider",
              "annotation_release_date", "annot_protein_coding",
              "annot_pseudogene", "annot_genes_total",
              "annot_pseudogene_frac"]

FAMILY_MODEL_COLS = ["case_id", "accession", "contig", "start", "end", "strand",
                     "gene_id", "name", "biotype", "family_from_name", "span",
                     "emits_protein"]

DB_RECORD_COLS = ["case_id", "species", "accession", "gene", "protein_name",
                  "length", "call", "confidence", "fragment", "reviewed",
                  "source"]
