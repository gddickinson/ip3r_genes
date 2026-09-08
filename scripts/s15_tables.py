"""s15_tables.py — every committed S15a table, written once, plus
`loss_dynamics_stats.json` (rules, parameters, self-test status and the
SHA-256 of every table).

The discipline is S6/S8/S9/S10's: nothing downstream re-derives a state, a
bar or a test statistic — the report and the figures render from these
files (D13).
"""

from __future__ import annotations

import json
import time

import s15_lib as lib
import s15_calibrate_recon as cal
import s15_integrity as integ
import s15_matrix as mx
import s15_reconstruct as recon
import s15_species_tree as tree
import s15_states as states
import s15_synteny_reach as syn

FILES = {
    "reconstruction": "reconstruction.tsv",
    "recon_calibration": "recon_calibration.tsv",
    "synteny_reach": "synteny_reach.tsv",
    "synteny_reach_summary": "synteny_reach_summary.tsv",
    "synteny_calls": "synteny_calls.tsv",
    "caller_accuracy_by_keys": "caller_accuracy_by_keys.tsv",
    "integrity_loci": "integrity_loci.tsv",
    "integrity_covariates": "integrity_covariates.tsv",
    "integrity_pairs": "integrity_pairs.tsv",
    "integrity_tests": "integrity_tests.tsv",
    "character_matrix": "character_matrix.tsv",
    "character_matrix_wide": "character_matrix_wide.tsv",
    "state_counts": "state_counts.tsv",
    "implied_copies": "implied_copies.tsv",
    "loss_candidates": "loss_candidates.tsv",
    "tree_placement": "tree_placement.tsv",
    "tree_polytomies": "tree_polytomies.tsv",
    "tree_vs_s13": "tree_vs_s13.tsv",
}


def write_all(payload: dict) -> list:
    out = lib.OUT
    out.mkdir(parents=True, exist_ok=True)
    written = []

    def w(key, rows, cols):
        written.append(lib.write_tsv(out / FILES[key], rows, cols))

    w("reconstruction", payload["reconstruction"], recon.COLS)
    w("recon_calibration", payload["recon_calibration"], cal.COLS)
    w("synteny_reach", payload["synteny_reach"], syn.REACH_COLS)
    w("synteny_reach_summary", payload["synteny_reach_summary"],
      list(payload["synteny_reach_summary"][0].keys())
      if payload["synteny_reach_summary"] else ["window"])
    w("synteny_calls", payload["synteny_calls"], syn.CALL_COLS)
    w("caller_accuracy_by_keys", payload["caller_accuracy_by_keys"],
      syn.ACC_COLS)
    w("integrity_loci", payload["integrity_loci"], integ.LOCUS_COLS)
    w("integrity_covariates", payload["integrity_covariates"], integ.CORR_COLS)
    w("integrity_pairs", payload["integrity_pairs"], integ.PAIR_COLS)
    w("integrity_tests", payload["integrity_tests"], integ.TEST_COLS)
    w("character_matrix", payload["character_matrix"], states.COLS)
    w("character_matrix_wide", payload["character_matrix_wide"], mx.WIDE_COLS)
    w("state_counts", payload["state_counts"], mx.COUNT_COLS)
    w("implied_copies", payload["implied_copies"], mx.COPY_COLS)
    w("loss_candidates", payload["loss_candidates"], mx.CAND_COLS)
    w("tree_placement", payload["tree_placement"], tree.AUDIT_COLS)
    w("tree_polytomies", payload["tree_polytomies"], tree.POLY_COLS)
    w("tree_vs_s13", payload["tree_vs_s13"], tree.CMP_COLS)
    (out / "species_tree_309.nwk").write_text(payload["newick"] + ";\n")
    written.append(out / "species_tree_309.nwk")
    return written


def write_stats(payload: dict, written: list, self_test: dict) -> None:
    stats = dict(
        task="S15a",
        written_at=time.strftime("%Y-%m-%d %H:%M:%S"),
        parameters=dict(
            rescue_evalue=lib.RESCUE_E,
            known_locus_pad_bp=lib.KNOWN_PAD,
            cov_found=states.COV_FOUND,
            recon_min_positives=cal.MIN_POS,
            recon_min_negatives=cal.MIN_NEG,
            integrity_bar_quantile=integ.BAR_QUANTILE,
            integrity_min_aligned_aa=integ.MIN_ALIGNED_AA,
            integrity_identity_window=integ.IDENTITY_WINDOW,
            synteny_window="informative10",
            synteny_min_keys=payload["synteny_reach_summary"][0]
            .get("min_keys_required") if payload["synteny_reach_summary"]
            else None,
        ),
        recon_calibration=payload["recon_calibration_stats"],
        integrity_bar=payload["integrity_bar"],
        state_vocabulary=list(states.STATES),
        headline=payload["headline"],
        self_test=self_test,
        tables={p.name: {"sha256": lib.sha256(p), "bytes": p.stat().st_size}
                for p in written},
    )
    (lib.OUT / "loss_dynamics_stats.json").write_text(
        json.dumps(stats, indent=1, default=str))


def headline(matrix: list[dict], copies: list[dict],
             reach_summary: list[dict]) -> dict:
    """The numbers the report is not allowed to compute for itself."""
    import collections
    st = collections.Counter(r["state"] for r in matrix)
    above = [c for c in copies if c["contig_spans_gene"]]
    below = [c for c in copies if not c["contig_spans_gene"]]
    rs = next((r for r in reach_summary
               if r.get("window") == "informative10"), {})
    return dict(
        n_cells=len(matrix), n_genomes=len(copies),
        n_absent=st.get("absent", 0),
        n_present=sum(st.get(s, 0) for s in states.PRESENT),
        n_undecided=sum(st.get(s, 0) for s in states.UNDECIDED),
        n_fragmented=st.get("present_fragmented", 0),
        n_paralog_unassignable=st.get("paralog_unassignable", 0),
        n_genomes_above_bar=len(above),
        n_genomes_below_bar=len(below),
        min_implied_above_bar=min((c["implied_copies"] for c in above),
                                  default=float("nan")),
        median_implied_above_bar=lib.median(
            [c["implied_copies"] for c in above]),
        n_below_2_5_above_bar=sum(1 for c in above
                                  if c["implied_copies"] < 2.5),
        n_below_2_5_below_bar=sum(1 for c in below
                                  if c["implied_copies"] < 2.5),
        synteny_regions=rs.get("n_regions", 0),
        synteny_reached=rs.get("n_reached", 0),
        synteny_no_neighbourhood=rs.get("n_no_neighbourhood", 0),
    )
