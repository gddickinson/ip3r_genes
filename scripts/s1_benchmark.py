"""S1 control benchmark — is the discovery scorer a validated instrument?

Mirrors the production wiring exactly (`src/cli.py:run_headless`): live
InterPro domain lookup for every panel accession, `analyse(use_mafft=True)`,
the MSA-derived family-signature fallback, then `discover_novel_paralogs`.
One deliberate difference: `max_variants` is raised so ALL panel sequences
enter the MSA — this benchmarks the scorer, not the production truncation.

Four discovery runs: a baseline where ITPR1/2/3 are all known, plus one
hold-out per paralog. A positive control is scored in the run that hides its
own name; the invertebrate grade is scored in the baseline, because no
hold-out can protect a member that was never name-protected.

Writes tables only — `s1_report.py` renders the report from them (D13).

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s1_benchmark.py
Outputs → results/benchmark_controls/
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s1_lib import (MafftTracer, bait_margins, components,  # noqa: E402
                            nearest_known_identity, score_lookup,
                            sister_margin, write_tsv)
from scripts.s1_panels import DECOY_SPEC, POSITIVE_SPEC, fetch_panel  # noqa: E402
from src.analysis.distance import identity_matrix  # noqa: E402
from src.analysis.pipeline import _label_for, analyse  # noqa: E402
from src.databases.interpro import batch_fetch_family_signatures  # noqa: E402
from src.discovery.candidates import (DiscoveryConfig,  # noqa: E402
                                      build_signature_set,
                                      discover_novel_paralogs)
from src.utils.family import (FAMILY_PFAM_IDS, MAX_LENGTH_AA,  # noqa: E402
                              MIN_LENGTH_AA, PARALOGS)

OUT_DIR = PROJECT_ROOT / "results" / "benchmark_controls"
PROMOTION_THRESHOLD = 40
REFERENCE_ACC = "Q14643"          # human ITPR1
DOMAIN_CACHE = OUT_DIR / "domain_hits.tsv"


def load_domain_hits(accessions: list[str], log=print) -> dict[str, list]:
    """Live InterPro family-signature lookup, cached to TSV so the benchmark
    reruns offline. This is the production route (src/cli.py) — and it is
    what hands the RyR decoys their +20 domain point, by design."""
    from src.databases.interpro import DomainHit
    if DOMAIN_CACHE.exists():
        out: dict[str, list] = {a: [] for a in accessions}
        for line in DOMAIN_CACHE.read_text().splitlines()[1:]:
            parts = (line.split("\t") + ["", "", "", ""])[:5]
            acc, pfam, pfam_name, n_hits, coverage = parts
            if pfam:
                out.setdefault(acc, []).append(DomainHit(
                    accession=acc, pfam_id=pfam, pfam_name=pfam_name,
                    n_hits=int(n_hits or 1), coverage=float(coverage or 0.0)))
        log(f"[domains] loaded {sum(len(v) for v in out.values())} hits "
            f"from cache {DOMAIN_CACHE.name}")
        return out
    log(f"[domains] InterPro lookup for {len(accessions)} accessions "
        f"(~{len(accessions) * 1.2:.0f}s)…")
    hits = batch_fetch_family_signatures(accessions)
    rows = [{"accession": acc, "pfam_id": h.pfam_id,
             "pfam_name": h.pfam_name, "n_hits": h.n_hits,
             "coverage": f"{h.coverage:.4f}"}
            for acc in accessions for h in hits.get(acc, [])]
    write_tsv(DOMAIN_CACHE, rows,
              ["accession", "pfam_id", "pfam_name", "n_hits", "coverage"])
    n_with = sum(1 for a in accessions if hits.get(a))
    log(f"[domains] {n_with}/{len(accessions)} accessions carry >=1 family "
        f"Pfam ({','.join(FAMILY_PFAM_IDS)})")
    return hits


def run_discovery(variants, known, label_map, analysis, domain_hits, sset):
    cfg = DiscoveryConfig(known_paralogs=known)
    return discover_novel_paralogs(
        variants, cfg, analysis_label_for=label_map,
        distances=analysis.distances, domain_hits=domain_hits,
        signature_set=sset)


def main() -> int:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log = print

    positives = fetch_panel(POSITIVE_SPEC, "positive", OUT_DIR, log=log)
    decoys = fetch_panel(DECOY_SPEC, "decoy", OUT_DIR, log=log)
    variants = positives + decoys
    accs = [v.accession for v in variants]
    domain_hits = load_domain_hits(accs, log=log)

    log(f"[bench] {len(positives)} positives + {len(decoys)} decoys — "
        f"aligning with MAFFT…")
    with MafftTracer() as tracer:
        analysis = analyse(variants, identity_threshold=0.4, use_mafft=True,
                           max_variants=len(variants),
                           reference_accession=REFERENCE_ACC)
    if not tracer.calls or tracer.calls[0]["returncode"] != 0:
        log("[bench] FATAL: MAFFT was not invoked (or failed) — "
            "use_mafft=True fell back to the star alignment; benchmark void.")
        return 1
    width = len(analysis.msa[0].aligned) if analysis.msa else 0
    log(f"[bench] MAFFT ran: rc={tracer.calls[0]['returncode']}, "
        f"{analysis.n_analyzed} rows, alignment width {width}")
    (OUT_DIR / "mafft_trace.json").write_text(json.dumps({
        "calls": tracer.calls, "n_rows": analysis.n_analyzed,
        "alignment_width": width}, indent=1))

    label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                 for v in variants}
    sset = build_signature_set(variants, analysis.msa, label_map, PARALOGS)

    runs = {}
    for name, known in [("baseline", list(PARALOGS))] + [
            (f"holdout_{p}", [q for q in PARALOGS if q != p]) for p in PARALOGS]:
        log(f"[bench] discovery run '{name}' (known={known})…")
        runs[name] = run_discovery(variants, known, label_map, analysis,
                                   domain_hits, sset)

    # ---------------- D14: labelled-bait margin ----------------
    label_by_acc = {v.accession: label_map[f"{v.source}|{v.accession}|"
                                           f"{v.gene_symbol}"]
                    for v in variants}
    margins = bait_margins(analysis.distances, label_by_acc)
    # Same margin under the fragment-aware metric the census will use in S6
    # (identity over mutually covered columns only). RyR is ~1.8x the length
    # of an ITPR, so the full-alignment metric dilutes every RyR-vs-ITPR
    # identity toward zero; this is the sensitivity check on that.
    covered = identity_matrix(analysis.msa, covered_only=True)
    margins_cov = bait_margins(covered, label_by_acc)
    margin_rows = []
    for v in variants:
        m = margins.get(label_by_acc.get(v.accession, ""), {})
        mc = margins_cov.get(label_by_acc.get(v.accession, ""), {})
        if not m:
            continue
        truth = "RYR" if v.raw.get("panel_note", "").startswith("RyR") else (
            "ITPR" if v.raw.get("panel_kind") == "positive" else "other")
        margin_rows.append({
            "accession": v.accession, "gene_symbol": v.gene_symbol,
            "species": v.species, "length_aa": v.length_aa or "",
            "truth": truth, "panel_note": v.raw.get("panel_note", ""),
            "id_to_itpr_bait": f"{m['id_to_itpr']:.3f}",
            "id_to_ryr_bait": f"{m['id_to_ryr']:.3f}",
            "margin": f"{m['margin']:+.3f}", "call": m["call"],
            "correct": "" if truth == "other" else str(m["call"] == truth),
            "id_to_itpr_bait_cov": f"{mc.get('id_to_itpr', 0):.3f}",
            "id_to_ryr_bait_cov": f"{mc.get('id_to_ryr', 0):.3f}",
            "margin_cov": f"{mc.get('margin', 0):+.3f}",
            "call_cov": mc.get("call", ""),
            "correct_cov": "" if truth == "other" else str(
                mc.get("call") == truth),
        })
    write_tsv(OUT_DIR / "bait_margin.tsv", margin_rows)

    # ---------------- positive controls ----------------
    pos_rows = []
    for v in positives:
        grp = v.raw.get("panel_note", "")
        run_name = f"holdout_{grp}" if grp in PARALOGS else "baseline"
        cand = score_lookup(runs[run_name])(v)
        score = cand.score if cand else None
        pos_rows.append({
            "group": grp, "run": run_name, "accession": v.accession,
            "gene_symbol": v.gene_symbol, "species": v.species,
            "length_aa": v.length_aa or "",
            "in_band": str(bool(v.length_aa and
                                MIN_LENGTH_AA <= v.length_aa <= MAX_LENGTH_AA)),
            "n_family_pfam": len(domain_hits.get(v.accession, [])),
            "id_to_nearest_known_pct": nearest_known_identity(cand),
            "sister_margin": sister_margin(cand),
            "score": score if score is not None else "n/a",
            "pass": str((score or 0) >= PROMOTION_THRESHOLD),
            "components": components(cand),
        })
    write_tsv(OUT_DIR / "positive_controls.tsv", pos_rows)

    # ---------------- why each recall failure failed ----------------
    # A missed positive is only actionable if the missing component is
    # named, so record for each failure its nearest non-known sibling under
    # both identity metrics — the breadth component needs one at
    # >= breadth_identity_min, and the full-alignment metric is the one that
    # dilutes it.
    fail_rows = []
    breadth_min = DiscoveryConfig(known_paralogs=list(PARALOGS)).breadth_identity_min
    known_lbls = {label_by_acc[v.accession] for v in positives
                  if v.raw.get("panel_note", "") in PARALOGS}
    for r in pos_rows:
        if r["pass"] == "True":
            continue
        lbl = label_by_acc.get(r["accession"], "")
        best = {"full": ("", 0.0), "cov": ("", 0.0)}
        for tag, dt in (("full", analysis.distances), ("cov", covered)):
            if lbl not in dt.labels:
                continue
            i = dt.labels.index(lbl)
            for j, other in enumerate(dt.labels):
                if j == i or other in known_lbls:
                    continue
                if dt.identity[i][j] > best[tag][1]:
                    best[tag] = (other, dt.identity[i][j])
        fail_rows.append({
            "accession": r["accession"], "gene_symbol": r["gene_symbol"],
            "species": r["species"], "score": r["score"],
            "components_fired": r["components"],
            "breadth_identity_min": f"{breadth_min:.2f}",
            "nearest_sibling_full": best["full"][0],
            "id_full": f"{best['full'][1]:.3f}",
            "breadth_would_fire_full": str(best["full"][1] >= breadth_min),
            "nearest_sibling_covered": best["cov"][0],
            "id_covered": f"{best['cov'][1]:.3f}",
            "breadth_would_fire_covered": str(best["cov"][1] >= breadth_min),
        })
    write_tsv(OUT_DIR / "recall_failures.tsv", fail_rows)

    # ---------------- negative controls ----------------
    neg_rows = []
    for v in decoys:
        per_run = {n: (score_lookup(r)(v).score if score_lookup(r)(v) else 0)
                   for n, r in runs.items()}
        worst_run = max(per_run, key=per_run.get)
        worst = per_run[worst_run]
        cand = score_lookup(runs[worst_run])(v)
        neg_rows.append({
            "category": v.raw.get("panel_note", ""), "accession": v.accession,
            "gene_symbol": v.gene_symbol, "species": v.species,
            "length_aa": v.length_aa or "",
            "in_band": str(bool(v.length_aa and
                                MIN_LENGTH_AA <= v.length_aa <= MAX_LENGTH_AA)),
            "n_family_pfam": len(domain_hits.get(v.accession, [])),
            "id_to_nearest_known_pct": nearest_known_identity(cand),
            "sister_margin": sister_margin(cand),
            "score_baseline": per_run["baseline"], "score_max": worst,
            "worst_run": worst_run,
            "pass": str(worst < PROMOTION_THRESHOLD),
            "components": components(cand),
        })
    neg_rows.sort(key=lambda r: -r["score_max"])
    write_tsv(OUT_DIR / "negative_controls.tsv", neg_rows)

    # ---------------- machine-readable summary ----------------
    scored = [r for r in pos_rows if r["score"] != "n/a"]
    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M"),
        "elapsed_s": round(time.time() - t0, 1),
        "promotion_threshold": PROMOTION_THRESHOLD,
        "size_band_aa": [MIN_LENGTH_AA, MAX_LENGTH_AA],
        "family_pfam_ids": list(FAMILY_PFAM_IDS),
        "reference_accession": REFERENCE_ACC,
        "n_positives": len(positives), "n_decoys": len(decoys),
        "n_msa_rows": analysis.n_analyzed, "alignment_width": width,
        "mafft_argv": tracer.calls[0]["argv"],
        "mafft_returncode": tracer.calls[0]["returncode"],
        "mafft_n_calls": len(tracer.calls),
        "signature_blocks": (len(sset.signatures) if sset else 0),
        "runs": list(runs),
        "recall_n": len(scored),
        "recall_hit": sum(1 for r in scored if r["pass"] == "True"),
        "specificity_n": len(neg_rows),
        "specificity_hit": sum(1 for r in neg_rows if r["pass"] == "True"),
        "ryr_n": sum(1 for r in neg_rows if r["category"].startswith("RyR")),
        "ryr_pass": sum(1 for r in neg_rows
                        if r["category"].startswith("RyR") and r["pass"] == "True"),
        "bait_margin_n": sum(1 for r in margin_rows if r["correct"]),
        "bait_margin_correct": sum(1 for r in margin_rows
                                   if r["correct"] == "True"),
        "bait_margin_correct_covered": sum(1 for r in margin_rows
                                           if r["correct_cov"] == "True"),
        "sister_margin_setting": 0.10,
        "n_recall_failures": len(fail_rows),
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1))
    log(f"[bench] recall {summary['recall_hit']}/{summary['recall_n']}, "
        f"specificity {summary['specificity_hit']}/{summary['specificity_n']} "
        f"(RyR {summary['ryr_pass']}/{summary['ryr_n']}), "
        f"bait margin {summary['bait_margin_correct']}/{summary['bait_margin_n']}"
        f" — {summary['elapsed_s']}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
