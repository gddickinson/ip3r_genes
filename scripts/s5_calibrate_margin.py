"""s5_calibrate_margin.py — is the inherited attribution margin usable here?

`s5_rescue.ATTRIBUTION_REL_MARGIN` decides whether a rescue trace is credited
to the paralog whose absence is being tested or reported as ambiguous. It
arrived as the PIEZO port's 1.5x bit-score ratio restated in D7's units
(0.333), and that project's paralogs are 40-50 % identical where this
family's are 61-68 % (S1). A threshold tuned on better-separated paralogs is
not evidence about these ones, so this measures what separation is actually
available.

**What is measured, and what it does and does not establish.** The sweep
records, at every locus it finds, how far the winning paralog's bait beats the
next paralog's (`paralog_margin`), and at loci the assembly's own annotation
names, that identity is established independently of the alignment. So these
are paralog margins at loci of *known* identity, with *complete* evidence —
a full-length gene, every exon available to the aligner.

That is an **upper bound** on what a rescue region can achieve, not a direct
calibration of it: a rescue region is a fragment of a gene, scored by tblastn
rather than by spliced alignment, so it separates paralogs less well than a
complete locus does. The bound is what makes the measurement decisive in one
direction only — if complete loci separate by less than the threshold, then
fragments certainly do, and every rescue trace would be reported ambiguous.
It cannot, on its own, say what the threshold *should* be. Rescue regions of
known identity are what would settle that, and they arrive with the full
sweep: a genome carrying two paralogs and genuinely missing the third gives
traces whose identity the other two cells establish.

Outputs -> results/genome_ledger/margin_calibration.tsv + .json

Usage:
  python scripts/s5_calibrate_margin.py
"""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import s5_sweep_lib as lib                                # noqa: E402
from s5_rescue import ATTRIBUTION_REL_MARGIN              # noqa: E402

LEDGER_DIR = PROJECT_ROOT / "results" / "genome_ledger"


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


def describe(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    v = sorted(values)
    return {"n": len(v), "min": round(v[0], 4), "max": round(v[-1], 4),
            "median": round(statistics.median(v), 4),
            "p05": round(v[max(0, int(0.05 * len(v)) - 1)], 4),
            "p25": round(v[max(0, int(0.25 * len(v)) - 1)], 4),
            "mean": round(statistics.fmean(v), 4)}


def main() -> None:
    path = LEDGER_DIR / "locus_margins.tsv"
    if not path.exists():
        raise SystemExit(f"{path} missing — run scripts/s5_ledger.py first")
    rows = read_tsv(path)

    def fval(r: dict, key: str) -> float | None:
        try:
            return float(r[key])
        except (KeyError, ValueError):
            return None

    # Loci whose paralog identity the annotation settles independently, and
    # where a second paralog's bait actually competed (margin 1.0 means no
    # other paralog aligned at all — real separation, but it carries no
    # information about how close a *contested* locus gets).
    itpr = [r for r in rows if r["class"] in lib.CLASSES]
    known = [r for r in itpr if r["annot_paralog"]]
    agree = [r for r in known if r["annot_agrees"] == "1"]
    disagree = [r for r in known if r["annot_agrees"] != "1"]
    contested = [r for r in agree if (fval(r, "paralog_margin") or 1.0) < 1.0]

    groups = {
        "all_itpr_loci": [fval(r, "paralog_margin") for r in itpr],
        "annotation_known": [fval(r, "paralog_margin") for r in known],
        "annotation_agrees": [fval(r, "paralog_margin") for r in agree],
        "annotation_agrees_contested": [fval(r, "paralog_margin")
                                        for r in contested],
        "annotation_disagrees": [fval(r, "paralog_margin") for r in disagree],
        "family_margin_all": [fval(r, "family_margin") for r in rows],
    }
    summary = {k: describe([v for v in vals if v is not None])
               for k, vals in groups.items()}

    contested_vals = [v for v in groups["annotation_agrees_contested"]
                      if v is not None]
    n_below = sum(1 for v in contested_vals if v < ATTRIBUTION_REL_MARGIN)
    # A threshold every correct, contested, complete locus clears — rounded
    # down to two places so it is a stated number rather than a data point.
    suggested = (int(min(contested_vals) * 100) / 100 if contested_vals
                 else None)

    fam_vals = [v for v in groups["family_margin_all"] if v is not None]
    fam_contested = [v for v in fam_vals if v < 1.0]

    verdict = []
    if contested_vals:
        verdict.append(
            f"{n_below}/{len(contested_vals)} contested loci of known, "
            f"correctly-called identity sit below the inherited threshold "
            f"{ATTRIBUTION_REL_MARGIN:.3f}.")
        if n_below:
            verdict.append(
                "Complete loci are an upper bound on rescue fragments, so a "
                "threshold that complete evidence fails cannot be met by a "
                "fragment: under it every rescue trace would be reported "
                "`tblastn_trace_ambiguous`, and no absence claim would ever "
                "be attributed. The inherited constant does not transfer.")
        verdict.append(
            f"Observed separation at contested loci: median "
            f"{summary['annotation_agrees_contested']['median']}, min "
            f"{summary['annotation_agrees_contested']['min']}. A threshold of "
            f"{suggested} is cleared by every one of them.")
    else:
        verdict.append("No contested locus in this run — every ITPR locus was "
                       "won by its own paralog's bait with no other paralog "
                       "aligning there at all, so the pilot puts no lower "
                       "bound on the margin.")
    if fam_contested:
        verdict.append(
            f"D14: {len(fam_contested)}/{len(fam_vals)} loci had both families' "
            f"baits align, the closest at a margin of {min(fam_contested):.3f}.")
    else:
        verdict.append(
            f"D14: at all {len(fam_vals)} loci only one family's baits aligned "
            "— the ITPR and RyR panels never contested a locus, so the family "
            "separation was decided before any margin had to be applied.")

    out_rows = [{"group": k, **v} for k, v in summary.items()]
    write_tsv(LEDGER_DIR / "margin_calibration.tsv",
              ["group", "n", "min", "p05", "p25", "median", "mean", "max"],
              out_rows)
    (LEDGER_DIR / "margin_calibration.json").write_text(json.dumps({
        "inherited_threshold": ATTRIBUTION_REL_MARGIN,
        "inherited_from": "PIEZO port, 1.5x bit-score ratio",
        "n_contested_below_threshold": n_below,
        "n_contested": len(contested_vals),
        "suggested_floor": suggested,
        "measures": "paralog margin at complete loci — an upper bound on "
                    "rescue fragments, not a direct calibration of them",
        "groups": summary, "verdict": verdict}, indent=2) + "\n")

    print(f"margin calibration over {len(rows)} loci "
          f"({len(known)} with an annotation-established paralog, "
          f"{len(contested_vals)} of those contested)")
    for k, v in summary.items():
        if v["n"]:
            print(f"  {k:<30} n={v['n']:<4} min={v['min']:<8} "
                  f"median={v['median']:<8} max={v['max']}")
    print()
    for line in verdict:
        print(f"  {line}")
    print(f"wrote {LEDGER_DIR / 'margin_calibration.tsv'}")


if __name__ == "__main__":
    main()
