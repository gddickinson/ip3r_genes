"""S20 step 3 — an individual verdict on every plant and fungal record.

S2 found the family's plant and fungal records are not scattered but
phylogenetically clean: in Viridiplantae every ITPR call was Chlorophyta and
Streptophyta had none; in Fungi every call sat in an early-diverging phylum
and Dikarya contributed nothing. That is a striking pattern *if the records
are real*, and an artefact if they are contaminants or mis-annotations. The
brief says to chase them one at a time, and this is that chase.

**The candidate set** is every record called ITPR in a plant or fungal
lineage by *either* instrument — census v3 (InterPro enumeration + the
profiles) and S20's own proteome sweep — pooled and de-duplicated. Taking
only the sweep's would quietly drop the UniProtKB entries that are not in a
reference proteome, which is most of the interesting ones.

**The rules**, in order, first match wins. Each is a positive test on
recorded evidence and each writes the number it fired on into the row:

  R0 `unresolved`           no sequence and/or no UniProt record — said so
                            rather than guessed
  R1 `module_only`          the profile match spans < 200 match states
                            (D22's own floor): the hit rests on a module the
                            family shares, not on a family domain
  R2 `contaminant_suspect`  ≥ 95 % identical, over ≥ 50 % of its length, to
                            a protein from a different kingdom. Two genuine
                            IP₃ receptors a billion years apart are 20–40 %
                            identical; 95 % across kingdoms is a sequence in
                            the wrong assembly
  R3 `cross_kingdom_outlier` 80–95 % to another kingdom — flagged, not
                            called, because that band has no clean reading
  R4 `no_genome_backing`    no EMBL and no proteome cross-reference: a
                            prediction with no assembly behind it
  R5 `fragment`             UniProt flags it a fragment, or it is under the
                            family length floor
  R6 `real_gene`            survives all of the above

**Why the contamination test is the load-bearing one.** It is the only
evidence here that can separate "a green alga has an IP₃ receptor" from "an
animal sequence is sitting in a green alga's assembly", and it is a positive
test on a measured identity rather than an appeal to plausibility.

Every UniProt response is archived under the data root, so a re-run needs
no network — there is deliberately no "offline mode" that fills the
paperwork in with placeholders, because a verdict rendered from placeholder
evidence is worse than no verdict.

Run:  python3 scripts/s20_verdicts.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_evidence import (  # noqa: E402
    best_by_partition, blastp, gather_sequences, make_db, uniprot_detail,
    work_dir, write_fasta,
)
from scripts.s20_groups import ALL_GROUPS  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    S20_DIR, group_paths, live, log, read_tsv, require_data_root, write_json,
    write_tsv,
)
from scripts.s3_assign import MIN_PROFILE_POSITIONS  # noqa: E402
from src.utils.family import MIN_LENGTH_AA  # noqa: E402

CENSUS_V3 = PROJECT_ROOT / "results" / "census_v3" / "census_v3.tsv"
#: S3 scored the *whole* S2 seeded space against both profiles when it
#: calibrated the instrument, so every census candidate already has a
#: measured `itpr_positions` there. Without it R1 (the shared-module rule)
#: could only fire on the records this task's own sweep reached, and the
#: same record would be judged by different evidence depending on which
#: instrument happened to find it.
CALIBRATION = PROJECT_ROOT / "results" / "hmm_sweep" / "calibration_assignments.tsv"

#: The two kingdoms the chase is about, and the split inside each that makes
#: a verdict interesting. Phyla are UniProt's, as recorded in census v3.
TARGET_KINGDOMS = {"Viridiplantae", "Fungi"}
LINEAGE_CLASS = {
    "Streptophyta": "land_plant_lineage",
    "Chlorophyta": "green_alga",
    "Ascomycota": "dikarya",
    "Basidiomycota": "dikarya",
}

CONTAMINANT_PIDENT = 95.0     # R2
CONTAMINANT_QCOV = 0.50       # R2
OUTLIER_PIDENT = 80.0         # R3

VERDICT_FIELDS = [
    "accession", "verdict", "rule", "evidence",
    "kingdom", "phylum", "lineage_class", "species", "taxon_id",
    "gene", "protein_name", "length", "source", "call_by",
    "itpr_score", "ryr_score", "rel_margin", "itpr_positions",
    "fragment", "protein_existence", "reviewed", "pfams",
    "n_embl", "n_proteomes",
    "out_kingdom_best", "out_kingdom_pident", "out_kingdom_qcov",
    "out_kingdom_group", "within_best", "within_pident", "within_group",
]


# ------------------------------------------------------------- candidates
def census_candidates() -> dict[str, dict]:
    """Plant/fungal ITPR calls from census v3, with S3's measured positions."""
    positions = ({r["accession"]: r["itpr_positions"]
                  for r in read_tsv(CALIBRATION)}
                 if CALIBRATION.exists() else {})
    out: dict[str, dict] = {}
    for r in read_tsv(CENSUS_V3):
        if r["call"] != "ITPR" or r["kingdom"] not in TARGET_KINGDOMS:
            continue
        out[r["accession"]] = {
            "accession": r["accession"], "kingdom": r["kingdom"],
            "phylum": r["phylum"], "species": r["species"],
            "taxon_id": r["taxon_id"], "gene": r["gene"],
            "protein_name": r["protein_name"], "length": r["length"],
            "itpr_score": r["itpr_score"], "ryr_score": r["ryr_score"],
            "rel_margin": r["rel_margin"], "fragment": r["fragment"],
            "protein_existence": r["protein_existence"],
            "reviewed": r["reviewed"], "pfams": r["pfams"],
            "itpr_positions": positions.get(r["accession"], ""),
            "source": "census_v3", "call_by": r["instruments"],
        }
    return out


def sweep_candidates(taxa: dict[int, dict]) -> dict[str, dict]:
    """Plant/fungal ITPR calls from S20's own proteome sweep."""
    out: dict[str, dict] = {}
    for group in ("viridiplantae", "fungi"):
        path = S20_DIR / f"assignments_{group}.tsv"
        if not path.exists():
            continue
        for r in read_tsv(path):
            if r["assignment"] != "ITPR":
                continue
            tax = taxa.get(int(r["taxon_id"] or 0), {})
            out[r["accession"]] = {
                "accession": r["accession"],
                "kingdom": tax.get("kingdom", ""),
                "phylum": tax.get("phylum", ""),
                "species": r["species"], "taxon_id": r["taxon_id"],
                "gene": r["gene"], "protein_name": r["protein_name"],
                "length": r["length"], "itpr_score": r["itpr_score"],
                "ryr_score": r["ryr_score"], "rel_margin": r["rel_margin"],
                "itpr_positions": r["itpr_positions"],
                "fragment": "", "protein_existence": "",
                "reviewed": r["reviewed"], "pfams": "",
                "source": f"s20_sweep:{group}", "call_by": "profile",
            }
    return out


def merge_candidates(a: dict[str, dict], b: dict[str, dict]) -> dict[str, dict]:
    """Union, recording that a record both instruments found was found twice."""
    out = dict(a)
    for acc, row in b.items():
        if acc in out:
            prev = out[acc]
            prev["source"] = f"{prev['source']}+{row['source']}"
            for k in ("itpr_positions", "kingdom", "phylum"):
                if not prev.get(k) and row.get(k):
                    prev[k] = row[k]
        else:
            out[acc] = row
    return out


# ------------------------------------------------ the comparison background
def reference_set(taxa: dict[int, dict], per_species: int = 1) -> dict[str, dict]:
    """One ITPR representative per species, every group — the blast background.

    One per species rather than everything: the contamination test needs
    taxonomic *breadth* (is this record's nearest relative in another
    kingdom?), and 8,000 near-identical vertebrate paralogs add none of it
    while costing the whole runtime.
    """
    best: dict[tuple[str, str], dict] = {}
    for r in read_tsv(CENSUS_V3):
        if r["call"] != "ITPR":
            continue
        key = (r["group"], r["taxon_id"])
        cur = best.get(key)
        if cur is None or int(r["length"] or 0) > int(cur["length"] or 0):
            best[key] = r
    for group in ALL_GROUPS:
        path = S20_DIR / f"assignments_{group}.tsv"
        if not path.exists():
            continue
        for r in read_tsv(path):
            if r["assignment"] != "ITPR":
                continue
            tax = taxa.get(int(r["taxon_id"] or 0), {})
            key = (tax.get("group", group), r["taxon_id"])
            cur = best.get(key)
            if cur is None or int(r["length"] or 0) > int(cur["length"] or 0):
                best[key] = {**r, "group": tax.get("group", group),
                             "kingdom": tax.get("kingdom", ""),
                             "phylum": tax.get("phylum", "")}
    return {r["accession"]: r for r in best.values()}


def fasta_sources() -> list[Path]:
    """Every archived FASTA a sequence might be in, cheapest first."""
    root = require_data_root()
    srcs = [root / "raw_api" / "uniprot" / "s2_seed_sweep.fasta"]
    srcs += [group_paths(g)["db"] for g in ALL_GROUPS]
    srcs += [root / "proteomes" / "vertebrata_refprot.fasta"]
    return [p for p in srcs if p.exists()]


# -------------------------------------------------------------- the rules
def classify(cand: dict, detail: dict | None, nn: dict) -> tuple[str, str, str]:
    """(verdict, rule, evidence) for one candidate. First match wins."""
    out = nn.get("outside")
    within = nn.get("within")
    pos = int(cand.get("itpr_positions") or 0)

    if not cand.get("_has_sequence"):
        return ("unresolved", "R0",
                "no sequence in any archived FASTA — cannot be tested")
    if detail is None:
        return ("unresolved", "R0",
                "no UniProt record returned for the accession")
    if pos and pos < MIN_PROFILE_POSITIONS:
        return ("module_only", "R1",
                f"profile match spans {pos} match states, under the "
                f"{MIN_PROFILE_POSITIONS}-state family-domain floor")
    if out and out["pident"] >= CONTAMINANT_PIDENT and out["qcov"] >= CONTAMINANT_QCOV:
        return ("contaminant_suspect", "R2",
                f"{out['pident']:.1f} % identical over {out['qcov']:.0%} of "
                f"its length to {out['sseqid']} ({out['_group']}) — a "
                f"cross-kingdom identity no genuine deep homolog reaches")
    if out and out["pident"] >= OUTLIER_PIDENT:
        return ("cross_kingdom_outlier", "R3",
                f"{out['pident']:.1f} % to {out['sseqid']} ({out['_group']}) "
                f"— above the homology range, below the contamination call")
    n_embl = len([x for x in (detail.get("EMBL") or "").split(";") if x.strip()])
    n_prot = len([x for x in (detail.get("Proteomes") or "").split(";") if x.strip()])
    if n_embl == 0 and n_prot == 0:
        return ("no_genome_backing", "R4",
                "no EMBL and no proteome cross-reference — a prediction with "
                "no assembly behind it")
    frag = (detail.get("Fragment") or "").strip()
    length = int(cand.get("length") or 0)
    if frag or (length and length < MIN_LENGTH_AA):
        return ("fragment", "R5",
                f"{'UniProt flags it a fragment; ' if frag else ''}"
                f"{length} aa against the {MIN_LENGTH_AA} aa family floor")
    near = (f"nearest outside-kingdom relative {out['pident']:.1f} % "
            f"({out['sseqid']}, {out['_group']})" if out
            else "no outside-kingdom blastp hit at E<=1e-3")
    same = (f"; nearest within-kingdom {within['pident']:.1f} %"
            if within else "; unique in its kingdom at E<=1e-3")
    return ("real_gene", "R6",
            f"{length} aa, {n_embl} EMBL / {n_prot} proteome cross-refs; "
            f"{near}{same}")


def run() -> int:
    from scripts.s20_taxa import load_table
    taxa = load_table()
    cands = merge_candidates(census_candidates(), sweep_candidates(taxa))
    log("s20_verdicts", f"{len(cands)} plant/fungal ITPR records to chase")
    if not cands:
        return 1

    wd = work_dir()
    refs = reference_set(taxa)
    log("s20_verdicts", f"background: {len(refs)} one-per-species ITPR "
                        f"representatives")

    wanted = set(cands) | set(refs)
    seqs = gather_sequences(wanted, fasta_sources())
    for acc, row in cands.items():
        row["_has_sequence"] = acc in seqs

    write_fasta(wd / "candidates.faa",
                {a: seqs[a][1] for a in cands if a in seqs})
    write_fasta(wd / "background.faa",
                {a: seqs[a][1] for a in wanted if a in seqs})
    db = make_db(wd / "background.faa")
    hits = blastp(wd / "candidates.faa", db, wd / "candidates_vs_background.tsv")
    log("s20_verdicts", f"blastp: {len(hits)} HSP rows")

    part: dict[str, str] = {}
    for acc, r in refs.items():
        part[acc] = r.get("group") or taxa.get(int(r.get("taxon_id") or 0),
                                               {}).get("group", "?")
    for acc, r in cands.items():
        part[acc] = ("Viridiplantae" if r["kingdom"] == "Viridiplantae"
                     else "Fungi" if r["kingdom"] == "Fungi" else "?")
    nn = best_by_partition(hits, part, part)
    for rec in nn.values():
        for h in rec.values():
            h["_group"] = part.get(h["sseqid"], "?")

    detail = uniprot_detail(sorted(cands), wd)

    rows = []
    for acc in sorted(cands):
        cand = cands[acc]
        d = detail.get(acc)
        verdict, rule, evidence = classify(cand, d, nn.get(acc, {}))
        out = nn.get(acc, {}).get("outside")
        win = nn.get(acc, {}).get("within")
        rows.append({
            **{k: cand.get(k, "") for k in VERDICT_FIELDS if k in cand},
            "accession": acc, "verdict": verdict, "rule": rule,
            "evidence": evidence,
            "lineage_class": LINEAGE_CLASS.get(cand.get("phylum", ""),
                                               "other_early_diverging"),
            "fragment": (d or {}).get("Fragment", cand.get("fragment", "")),
            "protein_existence": (d or {}).get("Protein existence",
                                               cand.get("protein_existence", "")),
            "pfams": (d or {}).get("Pfam", cand.get("pfams", "")),
            "n_embl": len([x for x in ((d or {}).get("EMBL") or "").split(";")
                           if x.strip()]),
            "n_proteomes": len([x for x in ((d or {}).get("Proteomes") or "").split(";")
                                if x.strip()]),
            "out_kingdom_best": out["sseqid"] if out else "",
            "out_kingdom_pident": round(out["pident"], 1) if out else "",
            "out_kingdom_qcov": out["qcov"] if out else "",
            "out_kingdom_group": out["_group"] if out else "",
            "within_best": win["sseqid"] if win else "",
            "within_pident": round(win["pident"], 1) if win else "",
            "within_group": win["_group"] if win else "",
        })

    write_tsv(S20_DIR / "plant_fungal_verdicts.tsv", VERDICT_FIELDS, rows)
    tally: dict[str, int] = {}
    for r in rows:
        key = f"{r['kingdom']}/{r['lineage_class']}/{r['verdict']}"
        tally[key] = tally.get(key, 0) + 1
    write_json(S20_DIR / "verdict_summary.json",
               {"n_candidates": len(rows), "tally": tally,
                "thresholds": {"contaminant_pident": CONTAMINANT_PIDENT,
                               "contaminant_qcov": CONTAMINANT_QCOV,
                               "outlier_pident": OUTLIER_PIDENT,
                               "min_profile_positions": MIN_PROFILE_POSITIONS,
                               "min_length_aa": MIN_LENGTH_AA},
                "n_background": len(refs)})
    for k, n in sorted(tally.items()):
        log("s20_verdicts", f"  {k:60s} {n}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    live([("fetch proteomes", True), ("sweep", True), ("relaxed panel", True),
          ("plant/fungal verdicts", False), ("jackhmmer", False),
          ("census v5", False)])
    return run()


if __name__ == "__main__":
    sys.exit(main())
