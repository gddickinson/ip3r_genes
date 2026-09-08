"""s15_synteny_reach.py — the brief's first step, and the number that
answers it.

S15 was asked to disambiguate the sweep's undecided `tblastn_trace` cells
by synteny, with measured accuracy against known loci.  S8 already built
that instrument — a per-paralog flank consensus, calibrated leave-one-
genome-out on the loci whose paralog their own assembly's annotation
establishes, with a matched random-window null — so the work here is not
to build a second one but to point it at the trace regions and measure
whether it arrives.

**It does not, and that is the result.**  A flank comparison needs a
neighbourhood, and a trace region does not have one: the reason a gene is
a trace is that its assembly is shattered, and the contigs the pieces sit
on are shorter than the gene.  `reach()` measures exactly that, per
region — how many coding genes lie on the region's own contig at all, how
many flank it, and how many informative symbols that yields — and
`accuracy_by_keys()` puts S8's own calibration on the same axis, so the
accuracy quoted is the accuracy available *at the number of keys the
trace regions actually have* rather than the accuracy S8 measured on
chromosome-level assemblies.

Two rules worth stating.  The window rules and symbol vocabularies are
S8's, imported unchanged: a second normalisation here would let the two
tasks disagree about what a flank is.  And a region on a contig carrying
no annotated gene is reported as `no_neighbourhood`, never as a
disagreement with the caller — the caller was never asked.
"""

from __future__ import annotations

import collections

import s15_lib as lib

import s8_flank_lib as flib
import s8_paralogon as pg


def _locus_from_region(r: dict) -> flib.Locus:
    return flib.Locus(
        accession=r["accession"], organism=r["organism"], vclass=r["vclass"],
        vorder="", cell=r["cell"], status=r["cell_status"],
        contig=r["contig"], start=int(r["start"]), end=int(r["end"]),
        idx=0, bait="", bait_paralog="", identity=0.0, coverage=0.0,
        annot_gene=r.get("overlapping_genes", ""), annot_paralog="")


def reach(rescue_rows: list[dict], windows=("fixed10", "informative10")
          ) -> list[dict]:
    """Per trace region: is there a neighbourhood to compare at all?"""
    by_acc: dict[str, list[dict]] = collections.defaultdict(list)
    for r in rescue_rows:
        if r["cell_status"].startswith("tblastn_trace"):
            by_acc[r["accession"]].append(r)
    out: list[dict] = []
    for acc in sorted(by_acc):
        genes = flib.load_genes(acc)
        if genes is None:
            for r in by_acc[acc]:
                out.append(dict(_base(r), window="", n_genes_on_contig=0,
                                contig_extent_bp=0, n_flanks=0, n_keys=0,
                                reached=0, why="no_gene_table"))
            continue
        idx = flib.index_by_contig(genes)
        extent: dict[str, int] = {}
        ngene: dict[str, int] = collections.Counter()
        for g in genes:
            extent[g[0]] = max(extent.get(g[0], 0), int(g[2]))
            ngene[g[0]] += 1
        for r in by_acc[acc]:
            loc = _locus_from_region(r)
            for win in windows:
                fs = flib.extract_flanks(idx, loc, window=win)
                nk = len(fs.keys_relaxed)
                why = ("ok" if nk >= pg.MIN_KEYS else
                       "no_neighbourhood" if ngene[r["contig"]] == 0 else
                       "too_few_keys")
                out.append(dict(_base(r), window=win,
                                n_genes_on_contig=ngene[r["contig"]],
                                contig_extent_bp=extent.get(r["contig"], 0),
                                n_flanks=fs.n_flanks, n_keys=nk,
                                reached=int(nk >= pg.MIN_KEYS), why=why))
    return out


def _base(r: dict) -> dict:
    return dict(accession=r["accession"], organism=r["organism"],
                vclass=r["vclass"], cell=r["cell"],
                cell_status=r["cell_status"], contig=r["contig"],
                start=int(r["start"]), end=int(r["end"]),
                region_kind=r.get("region_kind", ""),
                assigned_clade=r.get("assigned_clade", ""))


def summarise(rows: list[dict], window: str = "informative10") -> dict:
    g = [r for r in rows if r["window"] == window]
    n = len(g)
    if not n:
        return dict(window=window, n_regions=0)
    return dict(
        window=window, n_regions=n,
        n_reached=sum(r["reached"] for r in g),
        n_no_neighbourhood=sum(1 for r in g if r["why"] == "no_neighbourhood"),
        n_too_few_keys=sum(1 for r in g if r["why"] == "too_few_keys"),
        n_no_gene_table=sum(1 for r in g if r["why"] == "no_gene_table"),
        median_keys=lib.median([float(r["n_keys"]) for r in g]),
        median_genes_on_contig=lib.median(
            [float(r["n_genes_on_contig"]) for r in g]),
        median_contig_extent_kb=lib.median(
            [r["contig_extent_bp"] / 1000.0 for r in g
             if r["contig_extent_bp"]]),
        min_keys_required=pg.MIN_KEYS,
    )


def accuracy_by_keys(calibration_rows: list[dict],
                     bins=((0, 0), (1, 3), (4, 7), (8, 12), (13, 99))
                     ) -> list[dict]:
    """S8's own leave-one-genome-out calibration, binned by keys available.

    This is what makes the accuracy statement usable: S8 reported one
    number over its whole labelled set, and the question here is what the
    caller achieves with the handful of keys a shattered assembly offers.
    """
    out = []
    for lo, hi in bins:
        g = [r for r in calibration_rows
             if lo <= int(r["n_keys"]) <= hi]
        called = [r for r in g if r["call"] != "no_call"]
        correct = [r for r in called if r["call"] == r["truth"]]
        out.append(dict(
            keys_lo=lo, keys_hi=hi, n=len(g), n_called=len(called),
            n_correct=len(correct),
            call_rate=(len(called) / len(g)) if g else float("nan"),
            accuracy=(len(correct) / len(called)) if called else float("nan"),
        ))
    return out


def apply_caller(rescue_rows: list[dict], reach_rows: list[dict],
                 window: str = "informative10") -> list[dict]:
    """Run S8's caller on the trace regions it can reach.

    Built from S8's committed loci so the consensus is the same object S8
    calibrated; a consensus rebuilt from a different locus set would be a
    different instrument wearing S8's calibration.
    """
    reachable = {(r["accession"], r["cell"], r["contig"], int(r["start"]))
                 for r in reach_rows
                 if r["window"] == window and r["reached"]}
    if not reachable:
        return []
    loci, _, _ = flib.load_loci()
    sets_by_class: dict[str, list] = collections.defaultdict(list)
    cache: dict[str, dict] = {}
    for L in loci:
        if L.cell not in flib.ALL_CLASSES:
            continue
        if L.accession not in cache:
            genes = flib.load_genes(L.accession)
            cache[L.accession] = (flib.index_by_contig(genes)
                                  if genes is not None else {})
        idx = cache[L.accession]
        if not idx:
            continue
        sets_by_class[L.cell].append(
            flib.extract_flanks(idx, L, window=window))
    caller = pg.ConsensusCaller(sets_by_class, "relaxed")
    out = []
    for r in rescue_rows:
        key = (r["accession"], r["cell"], r["contig"], int(r["start"]))
        if key not in reachable:
            continue
        idx = cache.get(r["accession"])
        if idx is None:
            genes = flib.load_genes(r["accession"])
            idx = cache[r["accession"]] = (flib.index_by_contig(genes)
                                           if genes is not None else {})
        if not idx:
            continue
        fs = flib.extract_flanks(idx, _locus_from_region(r), window=window)
        sc = caller.score(fs, flib.ITPR_CLASSES, leave_out=True)
        out.append(dict(_base(r), window=window, n_keys=sc["n_keys"],
                        call=sc["call"], best=sc["best"],
                        best_score=sc["best_score"], margin=sc["margin"],
                        agrees_with_bitscore=int(
                            sc["call"] == r.get("assigned_clade"))))
    return out


REACH_COLS = ["accession", "organism", "vclass", "cell", "cell_status",
              "contig", "start", "end", "region_kind", "assigned_clade",
              "window", "n_genes_on_contig", "contig_extent_bp", "n_flanks",
              "n_keys", "reached", "why"]
CALL_COLS = ["accession", "organism", "vclass", "cell", "cell_status",
             "contig", "start", "end", "region_kind", "assigned_clade",
             "window", "n_keys", "call", "best", "best_score", "margin",
             "agrees_with_bitscore"]
ACC_COLS = ["keys_lo", "keys_hi", "n", "n_called", "n_correct",
            "call_rate", "accuracy"]
