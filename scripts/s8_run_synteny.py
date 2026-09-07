#!/usr/bin/env python3
"""S8 driver -- flanks, Jaccard matrices, the matched null, the paralogon
test and the consensus caller.

Usage:
  python scripts/s8_run_synteny.py [--flank-n 10] [--min-keys 4]
                                   [--limit N] [--skip-flanks-tsv]

Outputs (results/synteny/): see INTERFACE.md. Nothing here writes prose; the
report is rendered from these tables by s8_report.py (D13).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

import s8_control as ctl        # noqa: E402
import s8_flank_lib as lib      # noqa: E402
import s8_paralogon as par      # noqa: E402
import s8_tables as tab         # noqa: E402
import s8_test_flanks as tests  # noqa: E402

OUT = PROJECT_ROOT / "results" / "synteny"
LIVE = PROJECT_ROOT / "results" / "session_live.json"

WINDOWS = ("fixed10", "informative10")
PRIMARY_WINDOW = "fixed10"        # the brief's rule
CALL_WINDOW = "informative10"     # equal key-set size across genomes
PRIMARY_KEY = "relaxed"


def live(stage: str, **kw):
    try:
        LIVE.write_text(json.dumps(
            {"task": "S8", "stage": stage, "at": time.strftime("%H:%M:%S"),
             **kw}))
    except OSError:
        pass


# ------------------------------------------------------------ extraction

def build(loci, flank_n, limit=None):
    """One pass over the genomes: flank sets for every locus and window,
    plus the matched control windows. The gene table is loaded once."""
    by_acc = defaultdict(list)
    for L in loci:
        by_acc[L.accession].append(L)
    accs = sorted(by_acc)
    if limit:
        accs = accs[:limit]
    flank_sets: dict[str, list] = {w: [] for w in WINDOWS}
    controls: dict[str, dict] = {w: {} for w in WINDOWS}
    with_table, skipped = set(), []
    for i, acc in enumerate(accs, 1):
        genes = lib.load_genes(acc)
        if genes is None:
            skipped.extend(by_acc[acc])
            continue
        with_table.add(acc)
        index = lib.index_by_contig(genes)
        first = by_acc[acc][0]
        for w in WINDOWS:
            for L in by_acc[acc]:
                flank_sets[w].append(
                    lib.extract_flanks(index, L, n=flank_n, window=w))
            cw = ctl.sample_windows(index, acc, first.organism, first.vclass,
                                    n=flank_n, window=w)
            if cw:
                controls[w][acc] = cw
        if i % 25 == 0:
            live("flanks", done=i, total=len(accs))
    return flank_sets, controls, with_table, skipped


def one_per_species(sets_: list[lib.FlankSet]) -> list[lib.FlankSet]:
    """The best-covered locus per species per cell.

    Pair statistics are computed on this subset so that a within-ITPR1 pair
    compares like with like: a teleost ITPR1 cell holds itpr1a and itpr1b,
    and scoring itpr1a in one species against itpr1b in the next measures
    the 3R duplication, not orthology.
    """
    best: dict[tuple, lib.FlankSet] = {}
    for fs in sets_:
        k = (fs.locus.organism, fs.locus.cell)
        cur = best.get(k)
        if cur is None or ((fs.locus.coverage, fs.locus.identity, -fs.locus.idx)
                           > (cur.locus.coverage, cur.locus.identity, -cur.locus.idx)):
            best[k] = fs
    return [best[k] for k in sorted(best, key=lambda k: (k[1], k[0]))]


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flank-n", type=int, default=lib.DEFAULT_N)
    ap.add_argument("--min-keys", type=int, default=par.MIN_KEYS)
    ap.add_argument("--limit", type=int, default=None,
                    help="debug: only the first N genomes")
    ap.add_argument("--skip-flanks-tsv", action="store_true")
    args = ap.parse_args()
    par.MIN_KEYS = args.min_keys

    # the rules that decide what a flank is, what two symbols count as the
    # same gene, and when a consensus call is evidence are all negative-
    # controlled before anything is measured with them
    tests.self_test()

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    t0 = time.time()

    loci, excluded, missing_summary = lib.load_loci()
    live("load", n_loci=len(loci))
    flank_sets, controls, with_table, skipped = build(loci, args.flank_n,
                                                      args.limit)

    tab.write_loci(OUT / "loci.tsv", loci, with_table)
    if not args.skip_flanks_tsv:
        tab.write_flanks(OUT / "flanks.tsv",
                         flank_sets[PRIMARY_WINDOW] + flank_sets[CALL_WINDOW])

    # ---- the matrix / pair-statistics subset
    subset = {w: {c: one_per_species([fs for fs in flank_sets[w]
                                      if fs.locus.cell == c])
                  for c in lib.ALL_CLASSES} for w in WINDOWS}
    in_matrix = {(fs.locus.label, w)
                 for w in WINDOWS for c in lib.ALL_CLASSES
                 for fs in subset[w][c]}
    tab.write_locus_sets(OUT / "locus_sets.tsv",
                         flank_sets[PRIMARY_WINDOW] + flank_sets[CALL_WINDOW],
                         in_matrix)
    live("jaccard", **{c: len(subset[PRIMARY_WINDOW][c])
                       for c in lib.ALL_CLASSES})

    # ---- matrices (primary window, relaxed keys) + every pair class
    vclass_of = {fs.locus.label: fs.locus.vclass
                 for w in WINDOWS for c in lib.ALL_CLASSES
                 for fs in subset[w][c]}
    stats_rows = []
    within_pairs = {}
    for c in lib.ALL_CLASSES:
        within_pairs[c] = tab.write_matrix(
            OUT / f"jaccard_{c}.tsv", subset[PRIMARY_WINDOW][c], PRIMARY_KEY)

    def add(name, pairs, w, kind):
        stats_rows.append({**ctl.pair_class_row(name, pairs, controls[w], kind),
                           "window": w, "stratum": "all"})
        for stratum, sub in ctl.stratify_pairs(pairs, vclass_of).items():
            stats_rows.append({**ctl.pair_class_row(name, sub, controls[w], kind),
                               "window": w, "stratum": stratum})

    for w in WINDOWS:
        for kind in ("relaxed", "root"):
            for c in lib.ALL_CLASSES:
                pairs = (within_pairs[c] if (w == PRIMARY_WINDOW
                                             and kind == PRIMARY_KEY)
                         else lib.pairwise_jaccard(subset[w][c], None, kind=kind))
                add(f"within_{c}", pairs, w, kind)
            for i, a in enumerate(lib.ALL_CLASSES):
                for b in lib.ALL_CLASSES[i + 1:]:
                    pairs = lib.pairwise_jaccard(subset[w][a], subset[w][b],
                                                 kind=kind)
                    add(f"cross_{a}_vs_{b}", pairs, w, kind)
            # The RYR control cell is one cell holding three genes: the best
            # locus per genome is RYR2 in most species and RYR1 in others, so
            # the pooled within-RYR number is a mixture and understates the
            # method. Split by which bait won, and the sister family becomes
            # a second family the method is demonstrated on.
            for rp in ("RYR1", "RYR2"):
                sub = [fs for fs in subset[w]["RYR"]
                       if fs.locus.bait_paralog == rp]
                if len(sub) >= 10:
                    add(f"within_RYR_{rp}",
                        lib.pairwise_jaccard(sub, None, kind=kind), w, kind)
    tab.write_pair_stats(OUT / "pair_stats.tsv", stats_rows)

    # ---- consensus + paralogon
    bg = {w: {k: ctl.background_prevalence(controls[w], k)
              for k in ("relaxed", "root")} for w in WINDOWS}
    cons_rows, par_rows, par_summary = [], [], []
    for w in WINDOWS:
        for kind in ("relaxed", "root"):
            prev_bg, bg_n = bg[w][kind]
            for c in lib.ALL_CLASSES:
                prev, n_sp, _ = par.side_prevalence(subset[w][c], kind)
                # symbol as the final tiebreak (see s8_paralogon.shared_roots)
                ranked = sorted(prev.items(), key=lambda kv: (-kv[1], kv[0]))
                for sym, frac in ranked:
                    if frac < 0.10:
                        break
                    b = prev_bg.get(sym, 0.0)
                    cons_rows.append(dict(
                        cell=c, key_kind=kind, window=w, symbol=sym,
                        n_species=round(frac * n_sp), n_species_total=n_sp,
                        fraction=frac, background=b,
                        enrichment=(frac / b) if b else float("inf")))
        prev_bg, bg_n = bg[w]["root"]
        for i, a in enumerate(lib.ALL_CLASSES):
            for b_ in lib.ALL_CLASSES[i + 1:]:
                rows = par.shared_roots(subset[w][a], subset[w][b_], prev_bg)
                name = f"{a}_vs_{b_}"
                for r in rows:
                    par_rows.append({**r, "pair": f"{name}@{w}"})
                par_summary.append({**par.paralogon_summary(name, rows, bg_n),
                                    "pair": f"{name}@{w}"})
    tab.write_consensus(OUT / "flank_consensus.tsv", cons_rows)
    tab.write_paralogon(OUT / "paralogon_shared.tsv", par_rows)
    tab.write_paralogon_summary(OUT / "paralogon_summary.tsv", par_summary)

    # ---- the consensus caller: calibrate, null-test, then call the rest
    (calib_rows, calib_summary, unplaced_rows, null_rows, null_summary,
     sweep_rows) = call_loci(flank_sets, subset, controls, args)
    tab.write_calibration(OUT / "caller_calibration.tsv", calib_rows,
                          lib.ITPR_CLASSES)
    tab.write_unplaced(OUT / "unplaced_loci.tsv", unplaced_rows,
                       lib.ITPR_CLASSES)
    tab.write_caller_null(OUT / "caller_null.tsv", null_rows, lib.ITPR_CLASSES)
    tab.write_frac_sweep(OUT / "caller_frac_sweep.tsv", sweep_rows,
                         lib.ITPR_CLASSES)

    stats = dict(
        generated=time.strftime("%Y-%m-%d %H:%M"),
        runtime_s=round(time.time() - t0, 1),
        flank_n=args.flank_n, min_keys=args.min_keys,
        primary_window=PRIMARY_WINDOW, call_window=CALL_WINDOW,
        primary_key=PRIMARY_KEY,
        n_control_replicates=ctl.N_REPLICATES, control_seed=ctl.CONTROL_SEED,
        n_loci=len(loci), n_excluded_cells=len(excluded),
        n_missing_summary=len(missing_summary),
        n_genomes_with_gene_table=len(with_table),
        n_loci_without_gene_table=len(skipped),
        n_control_genomes=len(controls[PRIMARY_WINDOW]),
        self_tests=[n for n, _ in tests.self_test(verbose=False)],
        excluded=excluded,
        subset_sizes={c: len(subset[PRIMARY_WINDOW][c]) for c in lib.ALL_CLASSES},
        caller=calib_summary, caller_null=null_summary,
        unplaced_supported=sum(1 for r in unplaced_rows
                               if r["null_verdict"] == "supported"),
        unplaced_within_null=sum(1 for r in unplaced_rows
                                 if r["null_verdict"] == "within_null"),
        paralogon_min_frac=par.PARALOGON_MIN_FRAC,
        background_windows={w: bg[w]["root"][1] for w in WINDOWS},
    )
    # the SHA-256 of every committed table, so a rerun that drifts is
    # visible in the data (S3/S6's align_stats.json discipline). The one
    # source of drift this catches is real: ranked tables built by walking
    # a set are hash-seeded per process unless the sort has a final
    # tiebreak, and `flank_consensus.tsv` did not have one.
    stats["table_sha256"] = {
        f.name: hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(OUT.glob("*.tsv"))}
    (OUT / "synteny_stats.json").write_text(json.dumps(stats, indent=1))
    live("done")
    print(f"[s8] {len(loci)} loci, {len(with_table)} gene tables, "
          f"caller accuracy {calib_summary['accuracy']:.3f} on "
          f"{calib_summary['n_called']} calibration loci "
          f"(false-call rate {null_summary['false_call_rate']:.3f} on "
          f"{null_summary['n']} random windows); "
          f"tables → {OUT}")


def call_loci(flank_sets, subset, controls, args):
    """Calibrate the consensus caller on annotation-confirmed loci, then run
    it on every locus the sweep could not label."""
    w = CALL_WINDOW
    sets_by_class = {c: subset[w][c] for c in lib.ITPR_CLASSES}
    caller = par.ConsensusCaller(sets_by_class, "relaxed")

    labelled = [(fs, fs.locus.cell) for fs in
                [x for c in lib.ITPR_CLASSES for x in subset[w][c]]
                if fs.locus.annot_paralog in lib.ITPR_CLASSES
                and fs.locus.annot_paralog == fs.locus.cell]

    # the consensus threshold is measured, not typed: sweep first, pick by a
    # stated rule, then rebuild the caller at the chosen value
    sweep_rows = par.frac_sweep(sets_by_class, labelled,
                                lib.ITPR_CLASSES, controls[w])
    par.CONSENSUS_FRAC, frac_rule = par.select_frac(sweep_rows)
    caller = par.ConsensusCaller(sets_by_class, "relaxed")

    calib_rows, summary = par.calibrate(caller, labelled, lib.ITPR_CLASSES)
    null_rows, null_summary = par.caller_null(caller, controls[w],
                                              lib.ITPR_CLASSES)
    summary["consensus_frac"] = par.CONSENSUS_FRAC
    summary["consensus_frac_rule"] = frac_rule
    summary["consensus_sizes"] = {c: len(caller.consensus(c))
                                  for c in lib.ITPR_CLASSES}

    # every locus whose paralog the sweep did not establish
    targets = []
    for fs in flank_sets[w]:
        L = fs.locus
        reason = ""
        if L.cell == "RYR":
            continue
        if L.bait_paralog == "vertebrate_basal":
            reason = "won by an unlabelled basal-vertebrate bait"
        elif L.vclass in ("Myxini", "Hyperoartia"):
            reason = "cyclostome locus (no 2R paralog label available)"
        elif L.status in ("fragment", "assembly_gap"):
            reason = f"S5 status {L.status}"
        elif L.idx > 0:
            reason = "extra locus in a multi-copy cell"
        elif not L.annot_paralog:
            reason = "no paralog-named annotation at the locus"
        if reason:
            targets.append((fs, reason))

    tail, max_null = par.null_tail(null_summary)
    rows = []
    for fs, reason in targets:
        s = caller.score(fs, lib.ITPR_CLASSES, leave_out=True)
        L = fs.locus
        rows.append(dict(
            label=L.label, organism=L.organism, vclass=L.vclass, cell=L.cell,
            locus_idx=L.idx, status=L.status, contig=L.contig, start=L.start,
            end=L.end, bait_paralog=L.bait_paralog, annot_gene=L.annot_gene,
            n_keys=s["n_keys"], call=s["call"], best=s["best"],
            best_score=s["best_score"], margin=s["margin"],
            p_null=tail.get(s["best_score"], 0.0),
            null_verdict=par.null_verdict(s["best_score"], s["call"], max_null),
            **{f"score_{c}": s["scores"][c] for c in lib.ITPR_CLASSES},
            reason=reason,
            keys=";".join(sorted(lib.keys_of(fs, "relaxed")))))
    rows.sort(key=lambda r: (r["vclass"], r["organism"], r["cell"],
                             r["locus_idx"]))
    return calib_rows, summary, rows, null_rows, null_summary, sweep_rows


if __name__ == "__main__":
    main()
