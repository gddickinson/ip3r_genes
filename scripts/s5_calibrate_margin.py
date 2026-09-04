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
It cannot, on its own, say what the threshold *should* be.

**S5b closes that gap directly.** The full sweep produces rescue regions that
overlap an annotated gene the assembly itself names for a paralog, so their
identity is established independently of the bait scores — the same trick as
above, one evidence level down, on exactly the fragments the threshold acts
on. Those are reported here as `fragment_*` groups, and they are what the
threshold should actually be set from.

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

#: The value this project inherited from the PIEZO port — that project's 1.5x
#: bit-score ratio restated in D7's units. Kept as a literal so the report can
#: say what was inherited even after `ATTRIBUTION_REL_MARGIN` has been moved
#: off it; comparing the data against the live constant alone would quietly
#: rename whatever is currently in force as "inherited".
INHERITED_REL_MARGIN = 0.333


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


def fragment_calibration(fval) -> tuple[dict, list[str], list[float]]:
    """The measurement S5a could not make: margins on rescue *fragments*
    whose paralog the assembly's own annotation names."""
    path = LEDGER_DIR / "rescue_regions.tsv"
    if not path.exists():
        return {}, ["No rescue-region table yet — run scripts/s5_ledger.py "
                    "after a sweep that produced rescues."], []
    regions = read_tsv(path)
    known = [r for r in regions if r["annot_paralog"]]
    agree = [r for r in known if r["annot_agrees"] == "1"]
    wrong = [r for r in known if r["annot_agrees"] != "1"]
    vals = [v for v in (fval(r, "attribution_rel_margin") for r in agree)
            if v is not None]
    wrong_vals = [v for v in (fval(r, "attribution_rel_margin") for r in wrong)
                  if v is not None]
    groups = {
        "fragment_all": [fval(r, "attribution_rel_margin") for r in regions],
        "fragment_annotation_known": [fval(r, "attribution_rel_margin")
                                      for r in known],
        "fragment_annotation_agrees": vals,
        "fragment_annotation_disagrees": wrong_vals,
    }
    summary = {k: describe([v for v in vv if v is not None])
               for k, vv in groups.items()}
    lines = []
    if vals:
        below = sum(1 for v in vals if v < ATTRIBUTION_REL_MARGIN)
        would_lose = sum(1 for v in vals if v < INHERITED_REL_MARGIN)
        lines.append(
            f"**Fragment-level, the measurement S5a deferred**: {len(known)} "
            f"rescue region(s) overlap a gene the assembly names for a "
            f"paralog. The attribution agrees with that name on "
            f"{len(agree)}/{len(known)}, over margins {min(vals):.3f}-"
            f"{max(vals):.3f}.")
        lines.append(
            f"{below}/{len(vals)} of the correct fragment attributions fall "
            f"below the threshold in force ({ATTRIBUTION_REL_MARGIN}); "
            f"{would_lose}/{len(vals)} would have fallen below the inherited "
            f"{INHERITED_REL_MARGIN} and been reported ambiguous.")
        if wrong_vals:
            lines.append(
                f"{len(wrong)} region(s) were attributed to a paralog the "
                f"annotation names differently, at margins up to "
                f"{max(wrong_vals):.3f} — the population any threshold has to "
                "exclude, and the reason the floor is not simply the minimum "
                "correct margin.")
        else:
            lines.append(
                "No region was attributed against its annotation, so this "
                "evidence bounds the threshold from below only: it says how "
                "low the floor must be to keep correct calls, not how high it "
                "may go before wrong ones enter.")
    return summary, lines, vals


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
    n_below_inherited = sum(1 for v in contested_vals
                            if v < INHERITED_REL_MARGIN)
    n_below_live = sum(1 for v in contested_vals if v < ATTRIBUTION_REL_MARGIN)
    # A threshold every correct, contested, complete locus clears — rounded
    # down to two places so it is a stated number rather than a data point.
    suggested = (int(min(contested_vals) * 100) / 100 if contested_vals
                 else None)

    fam_vals = [v for v in groups["family_margin_all"] if v is not None]
    fam_contested = [v for v in fam_vals if v < 1.0]

    verdict = []
    if contested_vals:
        verdict.append(
            f"{n_below_inherited}/{len(contested_vals)} contested loci of "
            f"known, correctly-called identity sit below the **inherited** "
            f"threshold {INHERITED_REL_MARGIN:.3f}; "
            f"{n_below_live}/{len(contested_vals)} sit below the threshold "
            f"**in force** ({ATTRIBUTION_REL_MARGIN:.3f}).")
        if n_below_inherited:
            verdict.append(
                "Complete loci are an upper bound on rescue fragments, so a "
                "threshold that complete evidence fails cannot be met by a "
                "fragment: under the inherited value every rescue trace would "
                "be reported `tblastn_trace_ambiguous` and no absence claim "
                "would ever be attributed. That constant does not transfer.")
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

    frag_summary, frag_lines, frag_vals = fragment_calibration(fval)
    summary.update(frag_summary)
    verdict.extend(frag_lines)

    out_rows = [{"group": k, **v} for k, v in summary.items()]
    write_tsv(LEDGER_DIR / "margin_calibration.tsv",
              ["group", "n", "min", "p05", "p25", "median", "mean", "max"],
              out_rows)
    (LEDGER_DIR / "margin_calibration.json").write_text(json.dumps({
        "inherited_threshold": INHERITED_REL_MARGIN,
        "inherited_from": "PIEZO port, 1.5x bit-score ratio",
        "threshold_in_force": ATTRIBUTION_REL_MARGIN,
        "n_contested_below_inherited": n_below_inherited,
        "n_contested_below_in_force": n_below_live,
        "n_contested": len(contested_vals),
        "suggested_floor": suggested,
        "fragment_n_correct": len(frag_vals),
        "fragment_min_margin": (min(frag_vals) if frag_vals else None),
        "fragment_below_threshold": sum(
            1 for v in frag_vals if v < ATTRIBUTION_REL_MARGIN),
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
