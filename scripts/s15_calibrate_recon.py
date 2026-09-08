"""s15_calibrate_recon.py — the bar a reconstruction has to clear before a
shattered gene counts as present, measured rather than chosen.

`s15_reconstruct.py` gives every candidate cell two numbers: how much of
its reference the assembly holds outside every locus the aligner found
(`coverage`), and how much family exonic sequence that is in
gene-equivalents (`gene_equiv`).  Neither is interpretable until something
says what an *absent* gene scores.

The negative is the one the sweep's own output supplies without any new
search: regions the full 38-bait panel attributes to a paralog whose gene
the aligner **already placed at a locus in the same genome**.  That
paralog is accounted for, so the fragment set cannot be that gene, and
what it recovers is what cross-paralog similarity delivers on its own.

Three things this module refuses to do.

It does not quote a threshold without its separation.  `separation()`
reports Youden's J beside the bar (`s23_calibrate_loci.separation()`'s
rule), so a pair of populations that does not separate is visible as such
— and the first pass here *did not* separate, at J = 0.52, until the
decoy stopped counting fragments of co-shattered paralogs as negatives.

It does not put the operating point at Youden's own threshold.  Youden
returns the smallest cut with maximal J, which on a perfectly separated
pair lands exactly on the lowest positive — the most permissive bar
consistent with the data.  `operating_point()` takes the midpoint of the
gap instead and commits the gap's two edges, so a later run whose
populations have drifted into contact is visible in the table.

And it does not pass vacuously.  Fewer than `MIN_POS` positives or
`MIN_NEG` negatives is a failure with its own message
(`s20_test_sensitivity.py`'s rule): a bar read off two negatives is not a
calibration, and the `co_trace` population is reported beside both
precisely because its size is the measured extent of the problem the bar
cannot solve.
"""

from __future__ import annotations

import s15_lib as lib
import s15_reconstruct as recon

#: refuse to write a calibration read off fewer than this many rows
MIN_POS = 20
MIN_NEG = 5

#: the metrics the bar is measured on, both committed
METRICS = ("coverage", "gene_equiv")


def populations(rows: list[dict], metric: str
                ) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {s: [] for s in recon.SCOPES}
    for r in rows:
        if not r.get("bait"):
            continue
        out[r["scope"]].append(float(r[metric]))
    return out


def operating_point(pos: list[float], neg: list[float]) -> dict:
    """The midpoint of the gap between the populations, with both edges."""
    if not pos or not neg:
        return dict(operating_point=float("nan"), gap_lo=float("nan"),
                    gap_hi=float("nan"), gap=float("nan"))
    lo, hi = max(neg), min(pos)
    return dict(operating_point=(lo + hi) / 2.0, gap_lo=lo, gap_hi=hi,
                gap=hi - lo)


def separation(rows: list[dict], metric: str) -> dict:
    pops = populations(rows, metric)
    pos, neg, mid = (pops["own_clade"], pops["decoy_accounted"],
                     pops["co_trace"])
    y = lib.youden(pos, neg)
    op = operating_point(pos, neg)
    row = dict(metric=metric, n_positive=len(pos), n_negative=len(neg),
               n_co_trace=len(mid),
               positive_median=lib.median(pos), positive_min=min(pos) if pos else float("nan"),
               negative_median=lib.median(neg), negative_max=max(neg) if neg else float("nan"),
               co_trace_median=lib.median(mid),
               youden_threshold=y["threshold"], youden_j=y["j"],
               sensitivity=y["sens"], specificity=y["spec"], **op)
    row["separated"] = 1 if (y["j"] >= 0.99 and op["gap"] > 0) else 0
    row["verdict"] = _verdict(row)
    return row


def _verdict(row: dict) -> str:
    if row["n_positive"] < MIN_POS or row["n_negative"] < MIN_NEG:
        return (f"underpowered: {row['n_positive']} positives and "
                f"{row['n_negative']} negatives, floors are "
                f"{MIN_POS}/{MIN_NEG}")
    if row["separated"]:
        return (f"separated: J = {row['youden_j']:.2f}, gap "
                f"{row['gap_lo']:.3f}-{row['gap_hi']:.3f}")
    return (f"not separated: J = {row['youden_j']:.2f}; the bar cannot be "
            f"read off these populations")


def calibrate(rows: list[dict]) -> tuple[list[dict], dict]:
    """One row per metric, plus the chosen bar and whether it may be used."""
    table = [separation(rows, m) for m in METRICS]
    chosen = next((r for r in table if r["metric"] == "coverage"), None)
    usable = bool(chosen and chosen["separated"]
                  and chosen["n_positive"] >= MIN_POS
                  and chosen["n_negative"] >= MIN_NEG)
    stats = dict(
        metric=chosen["metric"] if chosen else "",
        bar=chosen["operating_point"] if chosen else float("nan"),
        gap_lo=chosen["gap_lo"] if chosen else float("nan"),
        gap_hi=chosen["gap_hi"] if chosen else float("nan"),
        youden_j=chosen["youden_j"] if chosen else float("nan"),
        n_positive=chosen["n_positive"] if chosen else 0,
        n_negative=chosen["n_negative"] if chosen else 0,
        n_co_trace=chosen["n_co_trace"] if chosen else 0,
        usable=int(usable), min_pos=MIN_POS, min_neg=MIN_NEG,
        verdict=chosen["verdict"] if chosen else "no rows",
    )
    return table, stats


def bar(stats: dict) -> float:
    """The committed bar, read back rather than retyped (D13 for constants).

    Refuses to hand out a bar the calibration marked unusable — a
    threshold nobody could measure must not silently become a default.
    """
    if not stats.get("usable"):
        raise RuntimeError(
            "the reconstruction bar was not established: " + str(
                stats.get("verdict", "no calibration")))
    return float(stats["bar"])


COLS = ["metric", "n_positive", "n_negative", "n_co_trace",
        "positive_median", "positive_min", "negative_median", "negative_max",
        "co_trace_median", "youden_threshold", "youden_j", "sensitivity",
        "specificity", "gap_lo", "gap_hi", "gap", "operating_point",
        "separated", "verdict"]
