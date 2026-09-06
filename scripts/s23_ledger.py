"""s23_ledger.py — assemble every completed summary.json into S23's tables.

D13: nothing is recomputed here. The per-genome `summary.json` files carry the
calls `s23_classify.py` made; this arranges them. A ledger that re-derived a
status would be a second classifier, and the two would eventually disagree.

Outputs -> results/s23_scope/
    copy_number_ledger.tsv     one row per genome: status, n_full, controls
    loci.tsv                   one row per ITPR locus, with its grade
    split_merges.tsv           every locus pair the conservative count folded
    control_ledger.tsv         the control verdict per genome, and why
    absence_at_genome.tsv      per absence clade: genomes, controls, findings
    ledger_stats.json          the counts the report is rendered from
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from src.utils.data_root import get_data_root                  # noqa: E402
import s23_scope as scope                                      # noqa: E402
from s5_build_baits import read_tsv, write_tsv                 # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "s23_scope"
MANIFEST = OUT_DIR / "genome_manifest_s23.tsv"

#: `n_full` is the **copy number** — distinct non-overlapping complete gene
#: models. `n_full_loci` is the older count of full-graded clusters, kept
#: beside it because the two differ exactly where a cluster chained two genes
#: or where several baits hit one, and that difference is a result.
LEDGER_COLS = ["accession", "organism", "group", "kingdom", "phylum", "class",
               "reasons", "status", "n_full", "n_full_loci", "n_fragment",
               "n_scrap", "n_below_floor", "best_below_floor_identity",
               "n_loci_raw", "n_merges", "best_coverage", "best_identity",
               "annotated_full", "family_named_full", "control_loci",
               "ryr_loci",
               "expects_ryr", "control_verdict", "control_why", "spans_gene",
               "contiguity_bar_bp", "genome_bp", "contig_n50", "annotated",
               "novel_model_count", "max_intron", "call_min_identity",
               "total_s"]

LOCI_COLS = ["accession", "organism", "group", "phylum", "contig", "start",
             "end", "strand", "grade", "bait", "band", "identity", "coverage",
             "score", "aligned_aa", "bait_len", "frameshifts", "stop_codons",
             "family_margin", "contig_edge", "n_gap", "annot_names_family",
             "span_bp", "cds_footprint_bp", "span_inflation"]

COPY_COLS = ["accession", "organism", "group", "phylum", "contig", "start",
             "end", "strand", "bait", "identity", "coverage", "score",
             "cds_footprint_bp"]


def load_summaries() -> list[dict]:
    root = get_data_root() / "s23_sweep"
    out = []
    for d in sorted(root.glob("*/summary.json")):
        if not (d.parent / ".sweep.done").exists():
            continue
        try:
            out.append(json.loads(d.read_text()))
        except json.JSONDecodeError:
            print(f"  ! unreadable summary: {d}")
    return out


def build(summaries: list[dict]) -> dict:
    ledger, loci, merges, controls, copies = [], [], [], [], []
    for s in summaries:
        ledger.append({c: s.get(c, "") for c in LEDGER_COLS})
        controls.append({
            "accession": s["accession"], "organism": s["organism"],
            "group": s.get("group", ""), "expects_ryr": s.get("expects_ryr"),
            "control_loci": s.get("control_loci", 0),
            "ryr_loci": s.get("ryr_loci", 0),
            "cross_kingdom_baits": ",".join(s.get("control_cross_groups")
                                            or []),
            "verdict": s.get("control_verdict", ""),
            "why": s.get("control_why", "")})
        for c in s.get("copies", []):
            copies.append({
                "accession": s["accession"], "organism": s["organism"],
                "group": s.get("group", ""), "phylum": s.get("phylum", ""),
                **{k: c.get(k, "") for k in COPY_COLS[4:]}})
        for d in s.get("loci", []):
            loci.append({
                "accession": s["accession"], "organism": s["organism"],
                "group": s.get("group", ""), "phylum": s.get("phylum", ""),
                **{k: d.get(k, "") for k in LOCI_COLS[4:]}})
        for m in s.get("merges", []):
            merges.append(dict(m, accession=s["accession"],
                               organism=s["organism"]))
    return {"ledger": ledger, "loci": loci, "merges": merges,
            "controls": controls, "copies": copies}


def absence_table(summaries: list[dict]) -> list[dict]:
    """Per absence clade: how many genomes, how many controlled, what was found.

    The table the negative claims are read off. A clade's absence counts only
    the genomes whose control fired — an uncontrolled genome contributes to
    neither column, which is the whole point of carrying a control at all.
    """
    presence = scope.load_presence()
    tax = scope.load_taxonomy()
    clades = scope.absence_clades(presence, tax)
    by_acc = {s["accession"]: s for s in summaries}
    manifest = {r["accession"]: r for r in read_tsv(MANIFEST)} \
        if MANIFEST.exists() else {}

    rows = []
    for c in clades:
        hits = []
        for acc, s in by_acc.items():
            row = manifest.get(acc, {})
            if (row.get("phylum") == c["clade"] and c["rank"] == "phylum") or \
               (row.get("class") == c["clade"] and c["rank"] == "class"):
                hits.append(s)
        controlled = [s for s in hits
                      if s.get("control_verdict", "").startswith("controlled")]
        with_itpr = [s for s in controlled if s.get("n_full", 0) > 0]
        traces = [s for s in controlled
                  if s.get("n_full", 0) == 0
                  and (s.get("n_fragment", 0) or s.get("n_scrap", 0)
                       or s.get("status") == "tblastn_trace")]
        rows.append({
            "rank": c["rank"], "clade": c["clade"],
            "swept_proteomes": c["swept"],
            "genomes_swept": len(hits),
            "genomes_controlled": len(controlled),
            "genomes_with_full_itpr": len(with_itpr),
            "genomes_with_trace_only": len(traces),
            "verdict": ("not yet swept" if not hits else
                        "uncontrolled" if not controlled else
                        "ITPR FOUND — proteome absence was an annotation fact"
                        if with_itpr else
                        "absence holds at assembly level"),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summaries = load_summaries()
    if not summaries:
        raise SystemExit("no completed genome summaries under "
                         f"{get_data_root() / 's23_sweep'}\n"
                         "  run: python3 scripts/s23_run_sweep.py --reasons anchor")
    t = build(summaries)
    write_tsv(OUT_DIR / "copy_number_ledger.tsv", LEDGER_COLS, t["ledger"])
    write_tsv(OUT_DIR / "loci.tsv", LOCI_COLS, t["loci"])
    write_tsv(OUT_DIR / "split_merges.tsv",
              ["accession", "organism", "contig", "strand", "kept_start",
               "kept_end", "merged_start", "merged_end", "gap_bp",
               "query_overlap", "kept_bait", "merged_bait", "why"],
              t["merges"])
    write_tsv(OUT_DIR / "copies.tsv", COPY_COLS, t["copies"])
    write_tsv(OUT_DIR / "control_ledger.tsv",
              ["accession", "organism", "group", "expects_ryr",
               "control_loci",
               "ryr_loci", "cross_kingdom_baits", "verdict", "why"],
              t["controls"])
    absences = absence_table(summaries)
    write_tsv(OUT_DIR / "absence_at_genome.tsv",
              ["rank", "clade", "swept_proteomes", "genomes_swept",
               "genomes_controlled", "genomes_with_full_itpr",
               "genomes_with_trace_only", "verdict"], absences)

    status = Counter(s["status"] for s in summaries)
    verdicts = Counter(s.get("control_verdict", "") for s in summaries)
    by_group = defaultdict(Counter)
    for s in summaries:
        by_group[s.get("group", "")][s["status"]] += 1
    stats = {
        "genomes": len(summaries),
        "manifest_genomes": len(read_tsv(MANIFEST)) if MANIFEST.exists() else 0,
        "status": dict(status),
        "control_verdicts": dict(verdicts),
        "uncontrolled": verdicts.get("uncontrolled", 0),
        "by_group": {g: dict(c) for g, c in sorted(by_group.items())},
        "total_full_loci": sum(s.get("n_full", 0) for s in summaries),
        "total_loci_raw": sum(s.get("n_loci_raw", 0) for s in summaries),
        "total_merges": sum(s.get("n_merges", 0) for s in summaries),
        "novel_models": sum(s.get("novel_model_count", 0) for s in summaries),
        "below_contiguity_bar": sum(1 for s in summaries
                                    if not s.get("spans_gene", True)),
        "copy_number": dict(Counter(s.get("n_full", 0) for s in summaries)),
        "total_full_loci_unmerged": sum(s.get("n_full_loci", 0)
                                        for s in summaries),
        "genomes_copy_differs_from_loci": sum(
            1 for s in summaries
            if s.get("n_full", 0) != s.get("n_full_loci", 0)),
        "loci_below_floor": sum(s.get("n_below_floor", 0) for s in summaries),
        "call_min_identity": (summaries[0].get("call_min_identity")
                              if summaries else None),
        "absence_clades_covered": sum(1 for a in absences
                                      if a["genomes_controlled"] > 0),
        "absence_clades": len(absences),
        "wall_clock_s": round(sum(s.get("total_s", 0) for s in summaries), 1),
    }
    (OUT_DIR / "ledger_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"{len(summaries)} genomes → {OUT_DIR}")
    print("  status: " + ", ".join(f"{k} {v}" for k, v in status.most_common()))
    print("  control: " + ", ".join(f"{k} {v}" for k, v in verdicts.most_common()))
    print(f"  {stats['total_full_loci']} full ITPR loci, "
          f"{stats['total_merges']} split-merges, "
          f"{stats['novel_models']} novel models")
    return 0


if __name__ == "__main__":
    sys.exit(main())
