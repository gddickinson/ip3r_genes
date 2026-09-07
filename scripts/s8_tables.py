#!/usr/bin/env python3
"""S8's table writers. Every number the report renders comes from one of
these files (D13); nothing downstream recomputes from the gene tables."""
from __future__ import annotations

from pathlib import Path

import s8_flank_lib as lib


def _w(fh, cols, row):
    fh.write("\t".join(_fmt(row[c]) for c in cols) + "\n")


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.4f}"
    if isinstance(v, bool):
        return str(int(v))
    return "" if v is None else str(v)


def write_rows(path: Path, cols: list[str], rows: list[dict]):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            _w(fh, cols, r)


# --------------------------------------------------------------- loci

LOCI_COLS = ["label", "accession", "organism", "vclass", "vorder", "cell",
             "locus_idx", "status", "contig", "start", "end", "span_bp",
             "identity", "coverage", "bait", "bait_paralog", "annot_gene",
             "annot_paralog", "has_gene_table"]


def write_loci(path: Path, loci: list[lib.Locus], with_table: set[str]):
    rows = [dict(
        label=L.label, accession=L.accession, organism=L.organism,
        vclass=L.vclass, vorder=L.vorder, cell=L.cell, locus_idx=L.idx,
        status=L.status, contig=L.contig, start=L.start, end=L.end,
        span_bp=L.end - L.start, identity=L.identity, coverage=L.coverage,
        bait=L.bait, bait_paralog=L.bait_paralog, annot_gene=L.annot_gene,
        annot_paralog=L.annot_paralog,
        has_gene_table=L.accession in with_table) for L in loci]
    write_rows(path, LOCI_COLS, rows)


# -------------------------------------------------------------- flanks

FLANK_COLS = ["label", "accession", "organism", "vclass", "cell", "window",
              "side", "rank", "symbol", "informative", "relaxed_key",
              "root_key", "gene_id", "biotype", "contig", "start", "end",
              "strand", "distance_bp"]


def write_flanks(path: Path, flank_sets: list[lib.FlankSet]):
    with open(path, "w") as fh:
        fh.write("\t".join(FLANK_COLS) + "\n")
        for fs in flank_sets:
            L = fs.locus
            for side, genes in (("up", fs.upstream), ("down", fs.downstream)):
                for rank, g in enumerate(genes, 1):
                    inf = lib.informative(g[5])
                    dist = L.start - g[2] if side == "up" else g[1] - L.end
                    _w(fh, FLANK_COLS, dict(
                        label=L.label, accession=L.accession,
                        organism=L.organism, vclass=L.vclass, cell=L.cell,
                        window=fs.window, side=side, rank=rank, symbol=g[5],
                        informative=inf,
                        relaxed_key=lib.relaxed_key(g[5]) if inf else "",
                        root_key=lib.root_key(g[5]) if inf else "",
                        gene_id=g[4], biotype=g[6], contig=g[0], start=g[1],
                        end=g[2], strand=g[3], distance_bp=dist))


# ---------------------------------------------------------- locus sets

SET_COLS = ["label", "accession", "organism", "vclass", "vorder", "cell",
            "locus_idx", "status", "window", "n_flanks", "scanned_up",
            "scanned_down", "n_informative", "n_root_keys", "in_matrix",
            "keys_relaxed"]


def write_locus_sets(path: Path, flank_sets: list[lib.FlankSet],
                     in_matrix: set[tuple[str, str]]):
    rows = []
    for fs in flank_sets:
        L = fs.locus
        rows.append(dict(
            label=L.label, accession=L.accession, organism=L.organism,
            vclass=L.vclass, vorder=L.vorder, cell=L.cell, locus_idx=L.idx,
            status=L.status, window=fs.window, n_flanks=fs.n_flanks,
            scanned_up=fs.scanned_up, scanned_down=fs.scanned_down,
            n_informative=fs.n_informative, n_root_keys=len(fs.keys_root),
            in_matrix=(L.label, fs.window) in in_matrix,
            keys_relaxed=";".join(sorted(fs.keys_relaxed))))
    write_rows(path, SET_COLS, rows)


# ------------------------------------------------------------- matrix

def write_matrix(path: Path, sets_: list[lib.FlankSet], kind: str):
    labels = [fs.locus.label for fs in sets_]
    pairs = lib.pairwise_jaccard(sets_, None, kind=kind)
    m = lib.matrix_from_pairs(labels, pairs)
    with open(path, "w") as fh:
        fh.write("label\t" + "\t".join(labels) + "\n")
        for a in labels:
            fh.write(a + "\t" + "\t".join(f"{m.get((a, b), 0.0):.4f}"
                                          for b in labels) + "\n")
    return pairs


# ---------------------------------------------------------- pair stats

PAIR_COLS = ["pair_class", "stratum", "key", "window", "n_pairs", "mean_j",
             "median_j",
             "max_j", "frac_positive", "control_mean_j", "control_median_j",
             "control_frac_positive", "excess_mean", "ratio", "n_used",
             "n_ties", "frac_beats_control", "z"]


def write_pair_stats(path: Path, rows: list[dict]):
    write_rows(path, PAIR_COLS, rows)


# ---------------------------------------------------------- consensus

CONS_COLS = ["cell", "key_kind", "window", "symbol", "n_species",
             "n_species_total", "fraction", "background", "enrichment"]


def write_consensus(path: Path, rows: list[dict]):
    write_rows(path, CONS_COLS, rows)


# ---------------------------------------------------------- paralogon

PARALOGON_COLS = ["pair", "key", "frac_a", "frac_b", "n_species_a",
                  "n_species_b", "n_vclass_a", "n_vclass_b", "background",
                  "enrichment", "vclass_a", "vclass_b"]
PARALOGON_SUMMARY_COLS = ["pair", "n_shared", "n_at_10", "n_at_25", "n_at_50",
                          "n_shared_bg_lt_05", "n_shared_bg_lt_01",
                          "n_control_windows", "top"]


def write_paralogon(path: Path, rows: list[dict]):
    write_rows(path, PARALOGON_COLS, rows)


def write_paralogon_summary(path: Path, rows: list[dict]):
    write_rows(path, PARALOGON_SUMMARY_COLS, rows)


# -------------------------------------------------------- calibration

def calib_cols(classes):
    return (["label", "organism", "vclass", "cell", "truth", "call", "best",
             "best_score", "margin", "n_keys"]
            + [f"score_{c}" for c in classes])


def write_calibration(path: Path, rows: list[dict], classes):
    write_rows(path, calib_cols(classes), rows)


# ------------------------------------------------------- unplaced loci

def unplaced_cols(classes):
    return (["label", "organism", "vclass", "cell", "locus_idx", "status",
             "contig", "start", "end", "bait_paralog", "annot_gene",
             "n_keys", "call", "best", "best_score", "margin", "p_null",
             "null_verdict"]
            + [f"score_{c}" for c in classes] + ["reason", "keys"])


def write_unplaced(path: Path, rows: list[dict], classes):
    write_rows(path, unplaced_cols(classes), rows)


# ------------------------------------------------------- caller controls

def null_cols(classes):
    return (["label", "organism", "vclass", "accession", "anchor_gene",
             "n_keys", "call", "best", "best_score", "margin"]
            + [f"score_{c}" for c in classes])


def write_caller_null(path: Path, rows: list[dict], classes):
    write_rows(path, null_cols(classes), rows)


def sweep_cols(classes):
    return (["consensus_frac", "n", "n_called", "call_rate", "accuracy",
             "n_control", "n_control_called", "false_call_rate", "margin"]
            + [f"consensus_size_{c}" for c in classes]
            + [f"call_rate_{c}" for c in classes]
            + [f"accuracy_{c}" for c in classes])


def write_frac_sweep(path: Path, rows: list[dict], classes):
    write_rows(path, sweep_cols(classes), rows)
