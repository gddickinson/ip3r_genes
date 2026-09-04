"""s5_census_v4.py — merge the sweep's genomic gene models into census v4.

Census v3 is 16,039 protein-database records carrying two independent
verdicts. It can only hold a gene some database already turned into a protein
entry. This adds the genes the sweep found **in DNA**, each keeping its
genomic address so the claim is auditable.

**The merge records how a database holds a gene, rather than only whether it
does.** The PIEZO port kept models by cell status alone — genome-only
evidence in, everything annotated out. That drops the most useful population
this project has: a locus that *is* an annotated gene but carries no family
name. S5a's very first pilot genome had one — *Takifugu* ITPR1's two 3R
duplicates side by side, one `itpr1b` and one `LOC101074739`. Both are in a
protein database; only one can be found by name, and only one is in census v3
under a family symbol. So `db_status` is a column:

  `genome_only`             no gene model at the locus at all, or the
                            assembly carries no annotation — the sweep is the
                            only evidence
  `annotated_unnamed`       a gene model exists but names no family
  `annotated_other_paralog` a gene model exists and names a *different*
                            paralog (or family) than the cell that won the
                            locus
  `annotated_named`         already named correctly (kept out of the census
                            rows; counted, because the denominator is the
                            result)

`annotated_other_paralog` is deliberately named for what was observed rather
than for a verdict. It says two labels disagree — the assembly's annotation
and this sweep's best-scoring bait — not that the annotation is wrong. One
instrument does not overturn a public annotation, and the paralog margins
this family affords are narrow (S5a measured 0.225-0.347 at complete loci).
`sibling_locus_for_annot_paralog` is what makes a case readable: when the
same genome *also* carries a separate locus for the paralog the annotation
names, the two cannot both be that paralog and the discrepancy is coherent.
S18 adjudicates these; S5 only supplies them.

Every model is scored against `itpr.hmm` and `ryr.hmm` — the same instrument
that called v3, so a v4 row rests on the same evidence as a v3 one rather
than on the bait that happened to find it (D23). A model the profiles cannot
call is reported and not enrolled, exactly as S3 did.

Outputs -> results/census_v4/
  census_v4.tsv         v3 rows + genome-derived rows, one schema
  genome_models.tsv     every model: address, cell, coverage, db_status, profile
  db_status_counts.tsv  how the databases hold each locus, by paralog and class
  unassigned_models.tsv models the profiles declined to call
  census_v4_stats.json  totals and the v3 -> v4 delta
  (sequences -> <data_root>/census_v4/genome_models.faa — bulk)

Usage:
  python scripts/s5_census_v4.py [--min-aa 200]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

for extra in ("/opt/anaconda3/envs/piezo1/bin", "/opt/homebrew/bin"):
    if Path(extra).is_dir() and extra not in os.environ.get("PATH", ""):
        os.environ["PATH"] = extra + os.pathsep + os.environ.get("PATH", "")

from src.utils.data_root import get_data_root, require_data_root  # noqa: E402
import s5_sweep_lib as lib                                  # noqa: E402
from s5_bait_screen import envelope_of, hmmsearch, profile_hits  # noqa: E402
from s3_assign import REL_MARGIN, MIN_SCORE, assign          # noqa: E402

RESULTS = PROJECT_ROOT / "results"
V3 = RESULTS / "census_v3" / "census_v3.tsv"
V4 = RESULTS / "census_v4"
HMM_DIR = RESULTS / "hmm_sweep"
ALL_CELLS = (*lib.CLASSES, lib.CONTROL_CLASS)

#: A model shorter than this is a scrap of alignment, not a gene model worth
#: enrolling. Matches S3's own floor for what counts as family evidence.
DEFAULT_MIN_AA = 200

#: Cell statuses where no protein-database gene model exists at the locus.
GENOME_ONLY_STATUSES = {"found_no_annotation", "tblastn_trace",
                        "tblastn_trace_ambiguous", "absent", "no_locus"}


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def db_status_of(cell_status: str, locus: dict, cell_class: str) -> str:
    """How a protein database holds this locus — not merely whether."""
    if cell_status in GENOME_ONLY_STATUSES or not locus.get("annot_gene"):
        return "genome_only"
    if locus.get("annot_paralog_matches"):
        return "annotated_named"
    if not locus.get("annot_family"):
        return "annotated_unnamed"
    return "annotated_other_paralog"


def collect_models(min_aa: int) -> tuple[list[dict], dict[str, str], Counter]:
    """Walk every completed sweep for its loci and their translations."""
    root = get_data_root() / "genome_sweep"
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
        translations = lib.read_fasta(d / "novel_models.faa") \
            if (d / "novel_models.faa").exists() else {}
        by_addr = {k.split("|", 2)[2].split("|")[0]: (k, v)
                   for k, v in translations.items()}
        # Which cells this genome actually filled, so a paralog-label
        # disagreement can be read against the rest of the genome.
        filled = {c for c in ALL_CELLS
                  if s["cells"][c]["status"].startswith("found")}
        for cell_class in ALL_CELLS:
            cell = s["cells"][cell_class]
            for rank, loc in enumerate(cell.get("loci", [])):
                status = db_status_of(cell["status"], loc, cell_class)
                tally[(cell_class, status)] += 1
                addr = (f"{loc['contig']}:{loc['start']}-{loc['end']}"
                        f"{loc['strand']}")
                hit = by_addr.get(addr)
                if status == "annotated_named" or hit is None:
                    continue
                name, seq = hit
                if len(seq) < min_aa:
                    tally[("_skipped_short", cell_class)] += 1
                    continue
                model_id = (f"{s['accession']}|{cell_class}|{addr}")
                seqs[model_id] = seq
                models.append({
                    "model_id": model_id, "accession": s["accession"],
                    "organism": s["organism"], "vclass": s.get("vclass", ""),
                    "vorder": s.get("vorder", ""),
                    "cell": cell_class, "cell_status": cell["status"],
                    "db_status": status, "locus_rank": rank,
                    "contig": loc["contig"], "start": loc["start"],
                    "end": loc["end"], "strand": loc["strand"],
                    "length": len(seq),
                    "coverage": loc.get("coverage", ""),
                    "identity": loc.get("identity", ""),
                    "family_margin": loc.get("family_margin", ""),
                    "paralog_margin": loc.get("paralog_margin", ""),
                    "frameshifts": loc.get("frameshifts", 0),
                    "stop_codons": loc.get("stop_codons", 0),
                    "bait": loc.get("bait", ""),
                    "assigned_via": loc.get("assigned_via", ""),
                    "annot_gene": (loc.get("annot_gene") or {}).get("name", ""),
                    "annot_family": loc.get("annot_family", ""),
                    "annot_paralog": loc.get("annot_paralog", ""),
                    "sibling_locus_for_annot_paralog": int(
                        bool(loc.get("annot_paralog"))
                        and loc.get("annot_paralog") in filled
                        and loc.get("annot_paralog") != cell_class),
                    "contig_spans_gene": int(lib.spans_a_gene(
                        int(s.get("contig_n50") or 0))),
                    "annotated_assembly": s.get("annotated", ""),
                })
    return models, seqs, tally


def score_models(seqs: dict[str, str]) -> dict[str, dict]:
    """The S3 profiles over every model — the instrument that called v3."""
    if not seqs:
        return {}
    verdicts: dict[str, dict] = {}
    with tempfile.TemporaryDirectory(prefix="s5_v4_") as tmp:
        work = Path(tmp)
        faa = work / "models.faa"
        lib.write_fasta(seqs, faa)
        hits, envs = {}, {}
        for fam, profile in (("itpr", HMM_DIR / "itpr.hmm"),
                             ("ryr", HMM_DIR / "ryr.hmm")):
            domtbl = work / f"{fam}.domtblout"
            hmmsearch(profile, faa, domtbl)
            hits[fam] = profile_hits(domtbl)
            envs[fam] = envelope_of(domtbl)
        rows = assign(hits["itpr"], hits["ryr"],
                      rel_margin=REL_MARGIN, min_score=MIN_SCORE)
        for r in rows:
            win = "itpr" if r["itpr_score"] >= r["ryr_score"] else "ryr"
            env = envs[win].get(r["accession"]) or {}
            verdicts[r["accession"]] = {
                "profile_call": r["assignment"],
                "profile_confidence": r["confidence"],
                "profile_evidence": r["evidence"],
                "profile_reason": r["reason"],
                "itpr_score": r["itpr_score"], "ryr_score": r["ryr_score"],
                "profile_rel_margin": r["rel_margin"],
                "profile_envelope": env.get("envelope_frac", 0.0),
            }
    return verdicts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-aa", type=int, default=DEFAULT_MIN_AA)
    args = ap.parse_args()

    v3 = read_tsv(V3)
    models, seqs, tally = collect_models(args.min_aa)
    if not models:
        raise SystemExit("no genome models found — run scripts/s5_run_sweep.py")
    verdicts = score_models(seqs)
    for m in models:
        m.update(verdicts.get(m["model_id"], {"profile_call": "no_hit",
                                              "profile_reason": "no profile hit"}))

    enrolled = [m for m in models if m["profile_call"] in ("ITPR", "RYR")]
    unassigned = [m for m in models if m["profile_call"] not in ("ITPR", "RYR")]

    # v4 rows: v3's schema, with the genome rows carrying their address.
    v3_cols = list(v3[0].keys())
    extra = ["source_detail", "genome_accession", "contig", "start", "end",
             "strand", "cell", "cell_status", "db_status", "coverage",
             "contig_spans_gene"]
    rows: list[dict] = []
    for r in v3:
        rows.append({**r, "source_detail": r.get("source", "census_v3")})
    for m in enrolled:
        rows.append({
            "accession": m["model_id"], "call": m["profile_call"],
            "confidence": m.get("profile_confidence", ""),
            "reason": m.get("profile_reason", ""),
            "instruments": "s5_genome_sweep",
            "profile_call": m["profile_call"],
            "profile_confidence": m.get("profile_confidence", ""),
            "itpr_score": m.get("itpr_score", ""),
            "ryr_score": m.get("ryr_score", ""),
            "rel_margin": m.get("profile_rel_margin", ""),
            "source": "S5_genome", "source_detail": m["db_status"],
            "gene": m["annot_gene"], "protein_name": "",
            "species": m["organism"], "class": m["vclass"],
            "order": m["vorder"], "length": m["length"],
            "genome_accession": m["accession"], "contig": m["contig"],
            "start": m["start"], "end": m["end"], "strand": m["strand"],
            "cell": m["cell"], "cell_status": m["cell_status"],
            "db_status": m["db_status"], "coverage": m["coverage"],
            "contig_spans_gene": m["contig_spans_gene"],
        })
    write_tsv(V4 / "census_v4.tsv", v3_cols + extra, rows)

    model_cols = ["model_id", "accession", "organism", "vclass", "vorder",
                  "cell", "cell_status", "db_status", "locus_rank", "contig",
                  "start", "end", "strand", "length", "coverage", "identity",
                  "family_margin", "paralog_margin", "frameshifts",
                  "stop_codons", "bait", "assigned_via", "annot_gene",
                  "annot_family", "annot_paralog",
                  "sibling_locus_for_annot_paralog", "annotated_assembly",
                  "contig_spans_gene", "profile_call", "profile_confidence",
                  "profile_evidence", "itpr_score", "ryr_score",
                  "profile_rel_margin", "profile_envelope", "profile_reason"]
    write_tsv(V4 / "genome_models.tsv", model_cols, models)
    write_tsv(V4 / "unassigned_models.tsv", model_cols, unassigned)

    status_rows = []
    for (cell, status), n in sorted(tally.items()):
        if str(cell).startswith("_"):
            continue
        status_rows.append({"cell": cell, "db_status": status, "n_loci": n})
    write_tsv(V4 / "db_status_counts.tsv", ["cell", "db_status", "n_loci"],
              status_rows)

    bulk = require_data_root() / "census_v4"
    bulk.mkdir(parents=True, exist_ok=True)
    lib.write_fasta({m["model_id"]: seqs[m["model_id"]] for m in enrolled},
                    bulk / "genome_models.faa")

    itpr_new = [m for m in enrolled if m["profile_call"] == "ITPR"]
    by_db = Counter(m["db_status"] for m in itpr_new)
    by_cell = Counter(m["cell"] for m in itpr_new)
    genomes = {m["accession"] for m in models}
    stats = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "min_aa": args.min_aa,
        "genomes_contributing": len(genomes),
        "v3_records": len(v3), "v4_records": len(rows),
        "models_collected": len(models), "models_enrolled": len(enrolled),
        "models_unassigned": len(unassigned),
        "itpr_models": len(itpr_new),
        "ryr_models": sum(1 for m in enrolled if m["profile_call"] == "RYR"),
        "itpr_by_db_status": dict(by_db),
        "itpr_by_cell": dict(by_cell),
        "loci_by_db_status": {f"{c}:{s}": n for (c, s), n in tally.items()
                              if not str(c).startswith("_")},
        "skipped_short": sum(v for k, v in tally.items()
                             if str(k[0]).startswith("_skipped")),
    }
    (V4 / "census_v4_stats.json").write_text(json.dumps(stats, indent=1) + "\n")

    print(f"census v4: {len(v3):,} v3 records + {len(enrolled):,} genome "
          f"models = {len(rows):,} ({len(genomes)} genomes)")
    print(f"  ITPR models {len(itpr_new):,}: "
          + ", ".join(f"{k} {v}" for k, v in by_db.most_common()))
    print(f"  by cell: " + ", ".join(f"{k} {v}" for k, v in by_cell.most_common()))
    if unassigned:
        print(f"  {len(unassigned)} model(s) the profiles declined to call — "
              "reported, not enrolled")
    print(f"wrote {V4}")


if __name__ == "__main__":
    main()
