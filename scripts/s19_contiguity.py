"""S19 — assembly contiguity as a confounder, measured on a family with no losses.

How much apparent gene loss is an assembly artefact? Usually that question
has to be answered with a control paralog and an argument. Here it does not:
**S15b reconstructs no losses anywhere in the 309-genome scope**, so every
one of the 923 assignable ITPR cells holds a gene that is there, and every
ledger cell that is not `found_*` is a false negative of the search. The
false-negative rate is therefore a direct measurement rather than an
estimate, and it can be read as a function of assembly contiguity.

The second series is the sister family. The RyR cell was swept in the same
pass by the same aligner, is present in every vertebrate as three genes, and
was never used to call an ITPR — so it is an independent replicate of the
same measurement, and the two agreeing is what stops the ITPR rate being an
artefact of how S15a assigned its states.

What that buys, and what the brief asks for: **D4's contiguity bar was chosen
a priori** as the median measured ITPR genomic span (142,212 bp of contig
N50) from gene geometry alone. `floor_scan` calibrates it — the residual
false-negative rate at every floor, against the genomes each floor costs —
so the bar is either validated or replaced by a measurement.

Two things this module refuses to do. It does not treat the 4
`paralog_unassignable` cyclostome cells as either present or absent, because
S15a declined to and re-deciding here would be a fork. And it does not read
the residual misses as biology without asking whether the *neighbourhood* is
there: contig N50 is a genome-wide statistic and a regional assembly defect
is invisible to it, so `neighbourhood_check` asks S8's own consensus flanks
whether the region survived at all.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_floor as FL                                          # noqa: E402
import s19_lib as S                                             # noqa: E402

# Re-exported so `s19_tables` and the self-test can read the calibration's
# parameters from the module that owns the measurement they belong to.
FLOOR_GRID = FL.FLOOR_GRID
FN_TARGETS = FL.FN_TARGETS
NEIGHBOURHOOD_PRESENT = FL.NEIGHBOURHOOD_PRESENT

NEGATIVE_GRADES = {
    "absent": "absent", "tblastn_trace": "trace",
    "tblastn_trace_ambiguous": "trace", "assembly_gap": "assembly_gap",
    "fragment": "fragment", "no_locus": "no_locus",
}


def grade(status: str) -> str:
    return "found" if S.is_found(status) else NEGATIVE_GRADES.get(status, status)


# ------------------------------------------------------------------- 1. cells

def cell_table(log=S.log) -> list[dict]:
    """Every cell, with whether it is a control and whether it was missed."""
    bar = S.contiguity_bar()
    rows = []
    for c in S.ledger_cells():
        if c["cell"] == S.CONTROL_CELL:
            control, why = "ryr_sister", "present in every vertebrate"
        elif c["s15_state"] in S.PRESENT_STATES:
            control, why = "itpr_present", f"S15a {c['s15_state']} ({c['s15_rule']})"
        elif c["s15_state"] in S.UNDECIDABLE_STATES:
            control, why = "", f"S15a {c['s15_state']} — not counted either way"
        else:
            control, why = "", f"S15a {c['s15_state'] or 'no state'}"
        rows.append({
            "accession": c["accession"], "organism": c["organism"],
            "vclass": c["vclass"], "cell": c["cell"], "status": c["status"],
            "grade": grade(c["status"]), "found": c["found"],
            "control": control, "control_reason": why,
            "false_negative": int(bool(control) and not c["found"]),
            "s15_state": c["s15_state"], "s15_rule": c["s15_rule"],
            "level": c["level"], "annotated_assembly": c["annotated_assembly"],
            "annotation_source": c["annotation_source"],
            "contig_n50": c["contig_n50"], "scaffold_n50": c["scaffold_n50"],
            "total_length_bp": c["total_length_bp"],
            "spans_gene": int(c["contig_n50"] >= bar),
            "best_coverage": S.fnum(c.get("best_coverage"), float, 0.0),
            "recon_coverage": c["recon_coverage"],
        })
    S.write_tsv(S.out_dir() / "contiguity_cells.tsv", list(rows[0].keys()),
                [list(r.values()) for r in rows])
    for name in ("itpr_present", "ryr_sister"):
        sub = [r for r in rows if r["control"] == name]
        fn = sum(r["false_negative"] for r in sub)
        log(f"contiguity_cells: {name} {fn}/{len(sub)} false negatives "
            f"({fn / max(1, len(sub)):.1%})")
    return rows


# -------------------------------------------------------------------- 2. bins

BINS = [(0, 50_000), (50_000, 142_212), (142_212, 500_000),
        (500_000, 2_000_000), (2_000_000, 10_000_000),
        (10_000_000, 10 ** 12)]


def _si(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1e6:.0f}M" if n % 1_000_000 == 0 else f"{n / 1e6:.2f}M"
    return f"{n / 1e3:.0f}k" if n % 1000 == 0 else f"{n / 1e3:.1f}k"


def _bin_label(lo: int, hi: int) -> str:
    if hi >= 10 ** 12:
        return f">={_si(lo)}"
    return f"<{_si(hi)}" if lo == 0 else f"{_si(lo)}-{_si(hi)}"


def bins(cells: list[dict], log=S.log) -> list[list]:
    """False-negative rate against contig N50, per series and per cell.

    Binned rather than smoothed, and the bin edges are placed **on** D4's bar
    rather than across it: a window straddling the threshold reports a rate no
    genome in it has, which is the one thing this table exists to show
    (`s5_figures.py`'s rule).
    """
    rows = []
    for lo, hi in BINS:
        label = _bin_label(lo, hi)
        sub = [c for c in cells if lo <= c["contig_n50"] < hi]
        series = [("itpr_present", [c for c in sub if c["control"] == "itpr_present"]),
                  ("ryr_sister", [c for c in sub if c["control"] == "ryr_sister"])]
        series += [(p, [c for c in sub if c["cell"] == p and c["control"]])
                   for p in S.PARALOGS]
        for name, s in series:
            if not s:
                continue
            fn = [c for c in s if c["false_negative"]]
            wlo, whi = S.wilson(len(fn), len(s))
            rows.append([label, lo, hi if hi < 10 ** 12 else "", name,
                         len({c["accession"] for c in s}), len(s), len(fn),
                         round(len(fn) / len(s), 4), round(wlo, 4),
                         round(whi, 4)])
    S.write_tsv(S.out_dir() / "contiguity_bins.tsv",
                ["bin", "contig_n50_lo", "contig_n50_hi", "series",
                 "n_genomes", "n_cells", "n_false_negative", "rate",
                 "wilson_lo", "wilson_hi"], rows)
    log(f"contiguity_bins: {len(rows)} rows")
    return rows


def by_stratum(cells: list[dict], log=S.log) -> list[list]:
    """The same rate broken out by assembly level, source (D9) and class."""
    rows = []
    ctrl = [c for c in cells if c["control"]]
    for field, label in (("level", "assembly_level"),
                         ("annotation_source", "annotation_source"),
                         ("vclass", "vertebrate_class")):
        for value in sorted({c[field] or "(none)" for c in ctrl}):
            sub = [c for c in ctrl if (c[field] or "(none)") == value]
            for name in ("itpr_present", "ryr_sister"):
                s = [c for c in sub if c["control"] == name]
                if not s:
                    continue
                fn = sum(c["false_negative"] for c in s)
                above = [c for c in s if c["spans_gene"]]
                fn_above = sum(c["false_negative"] for c in above)
                rows.append([label, value, name,
                             len({c["accession"] for c in s}), len(s), fn,
                             round(fn / len(s), 4), len(above),
                             round(fn_above / len(above), 4) if above else ""])
    S.write_tsv(S.out_dir() / "contiguity_by_stratum.tsv",
                ["stratum", "value", "series", "n_genomes", "n_cells",
                 "n_false_negative", "rate", "n_cells_above_bar",
                 "rate_above_bar"], rows)
    log(f"contiguity_by_stratum: {len(rows)} rows")
    return rows


# ------------------------------------------------------------------- 3. tests

def tests(cells: list[dict], log=S.log) -> list[list]:
    """Is the miss a contiguity effect, and is it the same in both series?"""
    from scipy import stats

    rows: list[list] = []
    series = {"itpr_present": [c for c in cells if c["control"] == "itpr_present"],
              "ryr_sister": [c for c in cells if c["control"] == "ryr_sister"]}
    for name, ctrl in series.items():
        found = [c["contig_n50"] for c in ctrl if not c["false_negative"]]
        miss = [c["contig_n50"] for c in ctrl if c["false_negative"]]
        if not miss:
            continue
        u = stats.mannwhitneyu(found, miss, alternative="greater")
        rows.append([f"mannwhitney_contig_n50_found_vs_missed", name,
                     f"n_found={len(found)}, n_missed={len(miss)}",
                     f"median {S.median(found):,.0f} vs {S.median(miss):,.0f}",
                     f"U={u.statistic:.0f}", f"{u.pvalue:.3g}"])
        chrom = [c for c in ctrl if c["level"] == "Chromosome"]
        other = [c for c in ctrl if c["level"] != "Chromosome"]
        table = [[sum(1 for c in chrom if c["false_negative"]),
                  sum(1 for c in chrom if not c["false_negative"])],
                 [sum(1 for c in other if c["false_negative"]),
                  sum(1 for c in other if not c["false_negative"])]]
        odds, p = stats.fisher_exact(table, alternative="less")
        rows.append(["fisher_chromosome_vs_lower", name,
                     f"chromosome {table[0][0]}/{sum(table[0])} missed; "
                     f"other {table[1][0]}/{sum(table[1])} missed",
                     f"odds={odds:.3g}", "", f"{p:.3g}"])
        # A logistic fit in stdlib: Newton-Raphson on one covariate. scipy is
        # already a dependency but statsmodels is not, and the odds ratio per
        # 10x of contig N50 is the number the report quotes.
        beta, se = _logit_log10(ctrl)
        if beta is not None:
            z = beta / se
            p = _two_sided_z(z)
            rows.append(["logit_found~log10(contig_n50)", name,
                         f"n={len(ctrl)}", f"beta={beta:.3f}",
                         f"OR per 10x={math.exp(beta):.2f}, z={z:.1f}",
                         f"{p:.3g}" if p > 0 else "<1e-300"])
        fn = sum(1 for c in ctrl if c["false_negative"])
        lo, hi = S.wilson(fn, len(ctrl))
        rows.append(["false_negative_rate", name, f"n={len(ctrl)}",
                     f"{fn / len(ctrl):.4f}", f"95% CI {lo:.4f}-{hi:.4f}", ""])

    a, b = series["itpr_present"], series["ryr_sister"]
    tab = [[sum(1 for c in a if c["false_negative"]),
            sum(1 for c in a if not c["false_negative"])],
           [sum(1 for c in b if c["false_negative"]),
            sum(1 for c in b if not c["false_negative"])]]
    odds, p = stats.fisher_exact(tab)
    rows.append(["fisher_itpr_vs_ryr_false_negative", "the two control series",
                 f"ITPR {tab[0][0]}/{sum(tab[0])}; RyR {tab[1][0]}/{sum(tab[1])}",
                 f"odds={odds:.3g}", "two-sided", f"{p:.3g}"])

    # The span-bias question S5b answered on recovery, asked of the misses.
    ctrl = [c for c in cells if c["control"] == "itpr_present"]
    below = [c for c in ctrl if not c["spans_gene"]]
    for p_ in S.PARALOGS:
        s = [c for c in below if c["cell"] == p_]
        if not s:
            continue
        fn = sum(c["false_negative"] for c in s)
        lo, hi = S.wilson(fn, len(s))
        rows.append(["false_negative_rate_below_bar", p_, f"n={len(s)}",
                     f"{fn / len(s):.4f}", f"95% CI {lo:.4f}-{hi:.4f}", ""])
    S.write_tsv(S.out_dir() / "contiguity_tests.tsv",
                ["test", "series", "n", "effect", "detail", "p"], rows)
    log(f"contiguity_tests: {len(rows)} tests")
    return rows


def _two_sided_z(z: float) -> float:
    """P(|Z| > |z|) via erfc, not 1 - Phi(z).

    `1 - 0.5*(1 + erf(z/sqrt2))` cancels to exactly 0.0 above |z| ~= 6 in
    double precision, and the first version of this table printed `0` for a
    z of 9.8 — a p-value of zero is not a strong result, it is a lost one.
    """
    return math.erfc(abs(z) / math.sqrt(2))


def _logit_log10(ctrl: list[dict], iters: int = 60):
    """found ~ b0 + b1*log10(contig_n50), Newton-Raphson, stdlib only.

    Returns (b1, se_b1) or (None, None) if the fit does not move — a
    separated or degenerate design should not silently report a coefficient.
    """
    xs = [[1.0, math.log10(max(1, c["contig_n50"]))] for c in ctrl]
    ys = [1.0 - c["false_negative"] for c in ctrl]
    b = [0.0, 0.0]
    for _ in range(iters):
        g = [0.0, 0.0]
        h = [[0.0, 0.0], [0.0, 0.0]]
        for x, y in zip(xs, ys):
            eta = b[0] * x[0] + b[1] * x[1]
            mu = 1 / (1 + math.exp(-max(-30.0, min(30.0, eta))))
            w = mu * (1 - mu)
            for i in range(2):
                g[i] += (y - mu) * x[i]
                for j in range(2):
                    h[i][j] += w * x[i] * x[j]
        det = h[0][0] * h[1][1] - h[0][1] * h[1][0]
        if abs(det) < 1e-12:
            return (None, None)
        inv = [[h[1][1] / det, -h[0][1] / det], [-h[1][0] / det, h[0][0] / det]]
        step = [inv[0][0] * g[0] + inv[0][1] * g[1],
                inv[1][0] * g[0] + inv[1][1] * g[1]]
        b = [b[0] + step[0], b[1] + step[1]]
        if max(abs(s) for s in step) < 1e-9:
            break
    return (b[1], math.sqrt(abs(inv[1][1])))


def run(log=S.log) -> dict:
    cells = cell_table(log)
    bins(cells, log)
    by_stratum(cells, log)
    tests(cells, log)
    _, rec = FL.floor_scan(cells, log)
    FL.clade_composition(cells, log)
    resid = FL.neighbourhood_check(cells, log)
    summary = {"n_cells": len(cells), "d4_bar": S.contiguity_bar(),
               "series": {}, "floor": rec}
    for name in ("itpr_present", "ryr_sister"):
        s = [c for c in cells if c["control"] == name]
        fn = sum(c["false_negative"] for c in s)
        lo, hi = S.wilson(fn, len(s))
        summary["series"][name] = {
            "n_cells": len(s), "n_false_negative": fn,
            "rate": round(fn / max(1, len(s)), 4),
            "wilson": [round(lo, 4), round(hi, 4)],
            "above_bar": _rate([c for c in s if c["spans_gene"]]),
            "below_bar": _rate([c for c in s if not c["spans_gene"]]),
        }
    summary["not_counted"] = sum(1 for c in cells if not c["control"])
    summary["residual"] = {
        "tested": len(resid),
        "neighbourhood_missing": sum(1 for r in resid
                                     if r[9] == "neighbourhood_missing"),
        "neighbourhood_present": sum(1 for r in resid
                                     if r[9] == "neighbourhood_present"),
        "no_gene_table": sum(1 for r in resid if r[9] == "no_gene_table"),
    }
    S.write_json(S.out_dir() / "contiguity_summary.json", summary)
    return summary


def _rate(sub: list[dict]) -> dict:
    fn = sum(c["false_negative"] for c in sub)
    return {"n": len(sub), "n_false_negative": fn,
            "rate": round(fn / max(1, len(sub)), 4)}


if __name__ == "__main__":
    run()
