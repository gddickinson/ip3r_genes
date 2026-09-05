"""S20 step 6 — census v5: the census with the whole eukaryotic tree in it.

Census v4 is the family as three instruments see it — InterPro's
architecture call (v2), the profile pair (v3) and the genomic sweep (v4) —
but its search space is the vertebrate-shaped one those tasks built. This
merge adds the fourth: `itpr.hmm`/`ryr.hmm` over every non-vertebrate
eukaryotic reference proteome plus the prokaryote sample.

**The merge rule is S3's, unchanged (D23).** A record both a previous
instrument and this sweep call keeps the agreed call and gains an
instrument; a record only this sweep reaches is enrolled on its profile
verdict alone, at the confidence a single instrument earns; a disagreement
becomes `conflict`, kept and reported rather than resolved by fiat. Writing
a new rule here would make v5's calls incomparable with v4's, which is the
whole reason for a versioned census.

**Every row gains lineage columns (D8)** from `s20_sweep/taxonomy.tsv`,
resolved from one source — UniProt taxonomy — for old rows and new alike.
v4's rows already carried `group`/`kingdom`/`phylum`; those are kept as
`*_v4` and the resolved columns sit beside them, so a taxonomy that has
moved since S2 is visible rather than overwritten.

**And every plant or fungal row carries its verdict** from step 3, because
"census v5 holds N Viridiplantae records" means nothing without saying how
many of them survived being chased.

Run:  python3 scripts/s20_census_v5.py
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_groups import ALL_GROUPS  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    CENSUS_V5_DIR, S20_DIR, live, log, read_tsv, write_json, write_tsv,
)
from scripts.s20_taxa import RANKS, load_table  # noqa: E402
from scripts.s3_hmm_lib import acc_key  # noqa: E402

CENSUS_V4 = PROJECT_ROOT / "results" / "census_v4" / "census_v4.tsv"

V4_FIELDS = [
    "accession", "call", "confidence", "reason", "instruments",
    "arch_call", "profile_call", "profile_confidence",
    "itpr_score", "ryr_score", "rel_margin", "searches", "source",
    "gene", "protein_name", "species", "taxon_id", "group", "kingdom",
    "phylum", "class", "order", "family", "length", "in_band", "fragment",
    "reviewed", "protein_existence", "in_alphafold", "seed_pfams", "pfams",
    "n_itpr_arch", "ryr_pfams", "source_detail", "genome_accession",
    "contig", "start", "end", "strand", "cell", "cell_status", "db_status",
    "coverage", "contig_spans_gene",
]
V5_EXTRA = (["s20_call", "s20_confidence", "s20_group", "s20_itpr_score",
             "s20_ryr_score", "s20_rel_margin", "s20_evidence"]
            + [f"tax_{r}" for r in RANKS]
            + ["tax_group", "plant_fungal_verdict", "plant_fungal_rule"])
V5_FIELDS = V4_FIELDS + V5_EXTRA


def merge_call(prev_call: str, prev_conf: str, s20_call: str,
               s20_conf: str) -> tuple[str, str, str, str]:
    """(call, confidence, reason, instruments) — S3's D23 rule, verbatim."""
    prev_spoke = prev_call in ("ITPR", "RYR")
    s20_spoke = s20_call in ("ITPR", "RYR")
    if prev_spoke and s20_spoke:
        if prev_call == s20_call:
            return (prev_call, "high",
                    f"census v4 and the S20 sweep both call {prev_call}",
                    "v4+s20")
        return ("conflict", "none",
                f"census v4 calls {prev_call}, the S20 sweep calls "
                f"{s20_call} — kept and reported", "v4+s20")
    if prev_spoke:
        return (prev_call, prev_conf,
                "census v4's call; the S20 sweep did not reach it", "v4")
    if s20_spoke:
        return (s20_call, s20_conf,
                "the S20 non-vertebrate sweep's profile call only", "s20")
    return ("unassigned", "none", "neither instrument called it", "none")


def s20_rows() -> dict[str, dict]:
    """Every S20 assignment, accession-keyed, best score wins a duplicate."""
    out: dict[str, dict] = {}
    for group in ALL_GROUPS:
        path = S20_DIR / f"assignments_{group}.tsv"
        if not path.exists():
            continue
        for r in read_tsv(path):
            acc = acc_key(r["accession"])
            prev = out.get(acc)
            if prev is None or float(r["itpr_score"] or 0) > float(
                    prev["itpr_score"] or 0):
                out[acc] = {**r, "group": group}
    return out


def verdict_index() -> dict[str, dict]:
    path = S20_DIR / "plant_fungal_verdicts.tsv"
    return ({acc_key(r["accession"]): r for r in read_tsv(path)}
            if path.exists() else {})


def build() -> dict:
    taxa = load_table()
    sweep = s20_rows()
    verdicts = verdict_index()
    log("census_v5", f"S20 sweep: {len(sweep)} scored targets; "
                     f"{len(verdicts)} plant/fungal verdicts; "
                     f"{len(taxa)} taxa resolved")

    # S5b's genomic gene models are keyed by genome accession and carry no
    # taxid, so their lineage is recovered from S4's declared manifest — the
    # same table that defined the genome scope in the first place.
    manifest = PROJECT_ROOT / "results" / "genome_manifest.tsv"
    genome_taxid = ({r["accession"]: r["taxid"] for r in read_tsv(manifest)}
                    if manifest.exists() else {})

    def lineage(taxid: str, genome_accession: str = "") -> dict:
        if not taxid and genome_accession:
            taxid = genome_taxid.get(genome_accession, "")
        t = taxa.get(int(taxid)) if taxid else None
        row = {f"tax_{r}": (t or {}).get(r, "") for r in RANKS}
        row["tax_group"] = (t or {}).get("group", "")
        return row

    rows: list[dict] = []
    seen: set[str] = set()
    for r in read_tsv(CENSUS_V4):
        acc = acc_key(r["accession"])
        seen.add(acc)
        s = sweep.get(acc)
        call, conf, reason, instruments = (
            merge_call(r["call"], r["confidence"], s["assignment"],
                       s["confidence"]) if s else
            (r["call"], r["confidence"], r["reason"], r["instruments"]))
        v = verdicts.get(acc, {})
        rows.append({
            **{k: r.get(k, "") for k in V4_FIELDS},
            "call": call, "confidence": conf, "reason": reason,
            "instruments": instruments,
            "s20_call": s["assignment"] if s else "",
            "s20_confidence": s["confidence"] if s else "",
            "s20_group": s["group"] if s else "",
            "s20_itpr_score": s["itpr_score"] if s else "",
            "s20_ryr_score": s["ryr_score"] if s else "",
            "s20_rel_margin": s["rel_margin"] if s else "",
            "s20_evidence": s["evidence"] if s else "",
            **lineage(r.get("taxon_id", ""),
                      r.get("genome_accession", "")),
            "plant_fungal_verdict": v.get("verdict", ""),
            "plant_fungal_rule": v.get("rule", ""),
        })

    novel = 0
    for acc, s in sorted(sweep.items()):
        if acc in seen:
            continue
        if s["assignment"] not in ("ITPR", "RYR"):
            continue          # an unassigned novel target enrols nothing
        novel += 1
        v = verdicts.get(acc, {})
        rows.append({
            **{k: "" for k in V4_FIELDS},
            "accession": acc, "call": s["assignment"],
            "confidence": s["confidence"],
            "reason": "the S20 non-vertebrate sweep's profile call only",
            "instruments": "s20", "profile_call": s["assignment"],
            "profile_confidence": s["confidence"],
            "itpr_score": s["itpr_score"], "ryr_score": s["ryr_score"],
            "rel_margin": s["rel_margin"],
            "searches": f"s20:{s['group']}", "source": "S20 sweep",
            "gene": s["gene"], "protein_name": s["protein_name"],
            "species": s["species"], "taxon_id": s["taxon_id"],
            "length": s["length"], "reviewed": s["reviewed"],
            "s20_call": s["assignment"], "s20_confidence": s["confidence"],
            "s20_group": s["group"], "s20_itpr_score": s["itpr_score"],
            "s20_ryr_score": s["ryr_score"], "s20_rel_margin": s["rel_margin"],
            "s20_evidence": s["evidence"],
            **lineage(s.get("taxon_id", "")),
            "plant_fungal_verdict": v.get("verdict", ""),
            "plant_fungal_rule": v.get("rule", ""),
        })

    CENSUS_V5_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(CENSUS_V5_DIR / "census_v5.tsv", V5_FIELDS, rows)
    log("census_v5", f"census v5: {len(rows)} records "
                     f"({novel} new from the S20 sweep)")

    # -------- delta and the range tables the report is rendered from -------
    v4_call = {acc_key(r["accession"]): r["call"] for r in read_tsv(CENSUS_V4)}
    changed = [r for r in rows
               if r["accession"] in v4_call
               and r["call"] != v4_call[r["accession"]]]
    write_tsv(CENSUS_V5_DIR / "call_changes.tsv",
              V5_FIELDS + ["v4_call"],
              [{**r, "v4_call": v4_call[r["accession"]]} for r in changed])
    write_tsv(CENSUS_V5_DIR / "conflicts.tsv", V5_FIELDS,
              [r for r in rows if r["call"] == "conflict"])

    by_group = Counter((r["tax_group"] or r["group"] or "unclassified",
                        r["call"]) for r in rows)
    groups = sorted({g for g, _ in by_group})
    write_tsv(CENSUS_V5_DIR / "calls_by_group.tsv",
              ["group", "ITPR", "RYR", "conflict", "unassigned", "total"],
              [{"group": g,
                "ITPR": by_group[(g, "ITPR")], "RYR": by_group[(g, "RYR")],
                "conflict": by_group[(g, "conflict")],
                "unassigned": by_group[(g, "unassigned")],
                "total": sum(n for (gg, _), n in by_group.items() if gg == g)}
               for g in groups])

    itpr = [r for r in rows if r["call"] == "ITPR"]
    # Fall back to census v4's own lineage columns where this task's
    # taxonomy resolution has no row: a record whose taxid UniProt taxonomy
    # does not return still has the placement S2 recorded for it, and
    # dropping it into "unclassified" would understate every group.
    def _key(r: dict) -> tuple[str, str]:
        return (r["tax_group"] or r["group"] or "unclassified",
                r["tax_phylum"] or r["phylum"] or "unclassified")

    by_phylum = Counter(_key(r) for r in itpr)
    taxa_by_phylum: dict[tuple[str, str], set] = {}
    for r in itpr:
        taxa_by_phylum.setdefault(_key(r), set()).add(r["taxon_id"])
    write_tsv(CENSUS_V5_DIR / "itpr_by_phylum.tsv",
              ["group", "phylum", "records", "taxa"],
              [{"group": g, "phylum": p, "records": n,
                "taxa": len(taxa_by_phylum[(g, p)])}
               for (g, p), n in sorted(by_phylum.items(),
                                       key=lambda kv: (kv[0][0], -kv[1]))])

    stats = {
        "n_records": len(rows), "n_novel_from_s20": novel,
        "calls": dict(Counter(r["call"] for r in rows)),
        "n_call_changes": len(changed),
        "n_conflicts": sum(1 for r in rows if r["call"] == "conflict"),
        "n_itpr": len(itpr),
        "n_itpr_taxa": len({r["taxon_id"] for r in itpr if r["taxon_id"]}),
        "itpr_groups": dict(Counter(r["tax_group"] or r["group"]
                                    or "unclassified" for r in itpr)),
    }
    write_json(CENSUS_V5_DIR / "census_v5_stats.json", stats)
    for k, v in stats.items():
        log("census_v5", f"  {k}: {v}")
    return stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    live([("fetch proteomes", True), ("sweep", True), ("relaxed panel", True),
          ("plant/fungal verdicts", True), ("jackhmmer", True),
          ("census v5", False)])
    build()
    live([("fetch proteomes", True), ("sweep", True), ("relaxed panel", True),
          ("plant/fungal verdicts", True), ("jackhmmer", True),
          ("census v5", True)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
