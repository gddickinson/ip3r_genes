"""s23_census_v6.py — merge the non-vertebrate sweep's gene models into census v6.

Census v5 is 17,882 records carrying, per row, an architecture call, a profile
call, S20's non-vertebrate sweep verdict and a lineage (D8). Like every census
before it, it can only hold a gene some database has already turned into a
protein entry. This adds the genes **S23 found in DNA** outside the
vertebrates, each with its genomic address.

**Same merge rule as v4 and v5, deliberately.** The profiles score every model
(`itpr.hmm` / `ryr.hmm`, D23) so a v6 row rests on the instrument that called
v3, not on whichever bait happened to find it; a model the profiles decline is
reported and not enrolled (D22's floor); and `db_status` records *how* a
database holds the locus rather than only whether:

  `genome_only`            no gene model at the locus, or the assembly carries
                           no annotation — this sweep is the only evidence
  `annotated_unnamed`      a gene model exists but names no family. This is
                           the big class here and it is not a defect of the
                           annotation: outside the vertebrates most gene
                           models carry locus tags, and S23a's pilot found
                           *Chlamydomonas*'s receptor filed as
                           `CHLRE_16g665450v5` and *Strongylocentrotus*'s as
                           `LOC594527`
  `annotated_other_family` a gene model exists and names the **ryanodine
                           receptor**. D14's sharpest possible disagreement,
                           so it is a class of its own rather than part of
                           `annotated_unnamed`
  `annotated_named`        already named for the family (counted, not
                           enrolled — the denominator is the result)

**What is different from v4.** There are no paralog cells. ITPR1/2/3 are a 2R
product and S5b found the trio absent below the cyclostomes, so a v6 row
carries a **copy index** within its genome instead of a cell, and the genome's
copy number beside it. That is the deliverable of S23, and putting it on the
census row is what lets a later task ask "how many ITPRs does this lineage
have" without going back to the ledger.

Outputs -> results/census_v6/
  census_v6.tsv          v5 rows + the non-vertebrate genome-derived rows
  genome_models_s23.tsv  every model: address, copy index, db_status, profile
  db_status_counts.tsv   how the databases hold each locus, by group
  unassigned_models.tsv  models the profiles declined to call
  census_v6_stats.json   totals and the v5 -> v6 delta
  (sequences -> <data_root>/census_v6/genome_models_s23.faa — bulk)

Usage:
  python3 scripts/s23_census_v6.py [--min-aa 200]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

for _extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(_extra).is_dir() and _extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = _extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root, require_data_root  # noqa: E402
import s5_sweep_lib as lib                                    # noqa: E402
from s5_census_v4 import read_tsv, write_tsv, score_models     # noqa: E402

RESULTS = PROJECT_ROOT / "results"
V5 = RESULTS / "census_v5" / "census_v5.tsv"
V6 = RESULTS / "census_v6"

#: Same floor S3 and v4 use for what counts as family evidence.
DEFAULT_MIN_AA = 200

#: Statuses at which no protein-database gene model exists at the locus.
GENOME_ONLY_STATUSES = {"found_no_annotation", "tblastn_trace", "no_locus"}


def db_status_of(genome_status: str, locus: dict) -> str:
    """How a protein database holds this locus — not merely whether."""
    if genome_status in GENOME_ONLY_STATUSES or not locus.get("annot_gene"):
        return "genome_only"
    if locus.get("annot_names_family"):
        return "annotated_named"
    if (locus.get("annot_family") or "") == "RYR":
        return "annotated_other_family"
    return "annotated_unnamed"


def collect_models(min_aa: int) -> tuple[list[dict], dict[str, str], Counter]:
    """Every S23 locus and its translation, from the completed sweeps."""
    root = get_data_root() / "s23_sweep"
    models: list[dict] = []
    seqs: dict[str, str] = {}
    tally: Counter = Counter()
    for d in sorted(root.glob("*/")):
        if not (d / ".sweep.done").exists():
            continue
        try:
            s = json.loads((d / "summary.json").read_text())
        except (json.JSONDecodeError, FileNotFoundError):
            continue
        translations = (lib.read_fasta(d / "novel_models.faa")
                        if (d / "novel_models.faa").exists() else {})
        by_addr = {k.split("|", 2)[2].split("|")[0]: (k, v)
                   for k, v in translations.items()}
        group = s.get("group", "")
        for rank, loc in enumerate(s.get("loci", [])):
            status = db_status_of(s.get("status", ""), loc)
            tally[(group, status)] += 1
            addr = f"{loc['contig']}:{loc['start']}-{loc['end']}"
            hit = by_addr.get(addr)
            if status == "annotated_named" or hit is None:
                continue
            _name, seq = hit
            if len(seq) < min_aa:
                tally[("_skipped_short", group)] += 1
                continue
            model_id = f"{s['accession']}|copy{rank + 1}|{addr}"
            seqs[model_id] = seq
            models.append({
                "model_id": model_id, "accession": s["accession"],
                "organism": s["organism"], "group": group,
                "kingdom": s.get("kingdom", ""), "phylum": s.get("phylum", ""),
                "class": s.get("class", ""), "order": s.get("order", ""),
                "genome_status": s.get("status", ""), "db_status": status,
                "copy_index": rank + 1,
                "genome_copy_number": s.get("n_full", 0),
                "contig": loc["contig"], "start": loc["start"],
                "end": loc["end"], "strand": loc["strand"],
                "length": len(seq), "grade": loc.get("grade", ""),
                "coverage": loc.get("coverage", ""),
                "identity": loc.get("identity", ""),
                "family_margin": loc.get("family_margin", ""),
                "span_bp": loc.get("span_bp", ""),
                "cds_footprint_bp": loc.get("cds_footprint_bp", ""),
                "frameshifts": loc.get("frameshifts", 0),
                "stop_codons": loc.get("stop_codons", 0),
                "bait": loc.get("bait", ""), "band": loc.get("band", ""),
                "annot_gene": (loc.get("annot_gene") or {}).get("name", ""),
                "annot_family": loc.get("annot_family", ""),
                "control_verdict": s.get("control_verdict", ""),
                "contig_spans_gene": int(bool(s.get("spans_gene"))),
                "annotated_assembly": s.get("annotated", ""),
            })
    return models, seqs, tally


def to_census_row(m: dict, v: dict, cols: list[str]) -> dict:
    """One model as a census row, in v5's own schema."""
    row = {c: "" for c in cols}
    row.update({
        "accession": m["model_id"],
        "call": v["profile_call"], "confidence": v["profile_confidence"],
        "reason": v["profile_reason"], "instruments": "s23_genome_sweep",
        "arch_call": "", "profile_call": v["profile_call"],
        "profile_confidence": v["profile_confidence"],
        "itpr_score": v["itpr_score"], "ryr_score": v["ryr_score"],
        "rel_margin": v["profile_rel_margin"],
        "source": "S23 genome sweep", "source_detail": m["bait"],
        "gene": m["annot_gene"], "protein_name": "",
        "species": m["organism"], "group": m["group"],
        "kingdom": m["kingdom"], "phylum": m["phylum"],
        "class": m["class"], "order": m["order"],
        "length": m["length"], "genome_accession": m["accession"],
        "contig": m["contig"], "start": m["start"], "end": m["end"],
        "strand": m["strand"], "cell": f"copy{m['copy_index']}",
        "cell_status": m["genome_status"], "db_status": m["db_status"],
        "coverage": m["coverage"],
        "contig_spans_gene": m["contig_spans_gene"],
        "tax_kingdom": m["kingdom"], "tax_phylum": m["phylum"],
        "tax_class": m["class"], "tax_order": m["order"],
        "tax_group": m["group"],
    })
    return {k: row.get(k, "") for k in cols}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-aa", type=int, default=DEFAULT_MIN_AA)
    args = ap.parse_args()
    require_data_root()

    if not V5.exists():
        raise SystemExit(f"census v5 missing: {V5}")
    v5 = read_tsv(V5)
    cols = list(v5[0].keys())

    models, seqs, tally = collect_models(args.min_aa)
    if not models:
        raise SystemExit("no S23 gene models — run scripts/s23_run_sweep.py")
    print(f"{len(models)} model(s) from the S23 sweep; scoring against "
          "itpr.hmm / ryr.hmm")
    verdicts = score_models(seqs)

    enrolled, declined = [], []
    for m in models:
        v = verdicts.get(m["model_id"])
        if not v or v["profile_call"] != "ITPR":
            declined.append(dict(m, why=(
                "no profile hit" if not v else
                f"profile called {v['profile_call']} "
                f"({v['profile_confidence']}) — {v['profile_reason']}")))
            continue
        m.update({k: v[k] for k in ("profile_call", "profile_confidence",
                                    "itpr_score", "ryr_score",
                                    "profile_rel_margin", "profile_envelope")})
        enrolled.append((m, v))

    V6.mkdir(parents=True, exist_ok=True)
    new_rows = [to_census_row(m, v, cols) for m, v in enrolled]
    write_tsv(V6 / "census_v6.tsv", cols, v5 + new_rows)
    write_tsv(V6 / "genome_models_s23.tsv",
              ["model_id", "accession", "organism", "group", "kingdom",
               "phylum", "class", "order", "genome_status", "db_status",
               "copy_index", "genome_copy_number", "contig", "start", "end",
               "strand", "length", "grade", "coverage", "identity",
               "family_margin", "span_bp", "cds_footprint_bp", "frameshifts",
               "stop_codons", "bait", "band", "annot_gene", "annot_family",
               "control_verdict", "contig_spans_gene", "annotated_assembly",
               "profile_call", "profile_confidence", "itpr_score",
               "ryr_score", "profile_rel_margin"],
              [m for m, _ in enrolled])
    write_tsv(V6 / "unassigned_models.tsv",
              ["model_id", "organism", "group", "phylum", "grade", "coverage",
               "identity", "length", "bait", "why"], declined)

    by_group: dict[str, Counter] = defaultdict(Counter)
    for (group, status), n in tally.items():
        if not group.startswith("_"):
            by_group[group][status] += n
    write_tsv(V6 / "db_status_counts.tsv",
              ["group", "status", "loci"],
              [{"group": g, "status": st, "loci": n}
               for g, c in sorted(by_group.items())
               for st, n in sorted(c.items())])

    bulk = get_data_root() / "census_v6"
    bulk.mkdir(parents=True, exist_ok=True)
    lib.write_fasta({m["model_id"]: seqs[m["model_id"]]
                     for m, _ in enrolled}, bulk / "genome_models_s23.faa")

    stats = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "census_v5_rows": len(v5), "census_v6_rows": len(v5) + len(new_rows),
        "models_found": len(models), "models_enrolled": len(enrolled),
        "models_declined": len(declined),
        "min_aa": args.min_aa,
        "db_status": {g: dict(c) for g, c in sorted(by_group.items())},
        "copy_number_rows": dict(Counter(m["genome_copy_number"]
                                         for m, _ in enrolled)),
    }
    (V6 / "census_v6_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"census v6: {len(v5)} + {len(new_rows)} = "
          f"{len(v5) + len(new_rows)} rows -> {V6}")
    print(f"  {len(declined)} model(s) the profiles declined")
    for g, c in sorted(by_group.items()):
        print(f"  {g:15s} " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))


if __name__ == "__main__":
    main()
