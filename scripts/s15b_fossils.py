"""s15b_fossils.py — the pseudogene-fossil analysis, and what to report
when the thing it is a test of does not exist.

The brief asks for a shared-lesion Poisson test: if a paralog died once on
a deep branch, its descendants' dead copies should share lesions, and the
number they share should exceed what independent decay would give.  That
test needs dead loci.  S15a found none — every locus above its measured
lesion bar is at full coverage — so this module's first job is to report
the **denominator**: how many loci were scored, how many cleared the bar,
and how many of those are dead by any of three independent readings.  An
omitted section would be indistinguishable from a section nobody ran.

Its second job is the lead the brief points at instead.  S15a's paired
within-genome test found that ITPR3 carries an indel excess against its
own genome's identity-matched sibling loci (39 genomes to 14,
q = 0.0032) and nothing in this project explains it.  `by_class()`
stratifies the same committed pairs by vertebrate class — no new
computation, D47's own instrument — and asks whether the excess is one
lineage or the whole family.  Every test is the identity-matched form,
because D47 says a statement about disabling lesions is made on that form
or not at all, and the whole stratified set is BH-corrected together
(S9's rule: one family, one correction).

`dead_locus()` is deliberately generous to the hypothesis it fails to
find.  Three readings are offered and a locus counts as a fossil if *any*
of them fires — a lesion density above the bar together with an incomplete
locus, a premature stop, or a cell whose paralog state is not a live gene.
A test made hard to pass would make the zero uninformative.
"""

from __future__ import annotations

import s15b_coding as coding
import s15b_lib as lib

#: the sweep's coverage bar for a complete locus (s15_states.COV_FOUND)
COV_FOUND = 0.70
#: S10's rule: zero internal stops falsifies a pseudogene call, and a
#: handful does not establish one, so a stop is a reading and not a verdict
MIN_STOPS_FOR_FOSSIL = 1


def _f(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def dead_locus(locus: dict, cell_state: str) -> dict:
    """Three independent readings of 'this locus is a dead gene'.

    The third reading is only *available* for a cell the character matrix
    holds.  The RyR control is not an ITPR cell and has no state there, so
    asking it whether its cell is coded present would answer `no` for
    every RyR locus in the sweep and flag the entire control as fossils.
    """
    elevated = locus.get("verdict") == "elevated_lesions"
    incomplete = _f(locus.get("coverage")) < COV_FOUND
    stops = int(_f(locus.get("stop_codons")))
    in_matrix = cell_state != "not_in_matrix"
    not_live = in_matrix and cell_state not in coding.PRESENT_STATES
    reasons = []
    if elevated and incomplete:
        reasons.append("above the lesion bar and below the coverage bar")
    if elevated and stops >= MIN_STOPS_FOR_FOSSIL:
        reasons.append(f"above the lesion bar with {stops} internal stops")
    if elevated and not_live:
        reasons.append(f"above the lesion bar in a cell coded {cell_state}")
    return dict(is_fossil=int(bool(reasons)),
                reading=";".join(reasons) or "not a fossil under any reading",
                elevated=int(elevated), incomplete=int(incomplete),
                stop_codons=stops, not_live=int(not_live),
                state_reading_available=int(in_matrix))


def denominator(loci: list[dict], matrix_rows: list[dict],
                bar: float) -> tuple[list[dict], dict]:
    """Every scored locus, its fossil reading, and the counts behind them."""
    state = {(r["accession"], r["cell"]): r["state"] for r in matrix_rows}
    out, n_scored, n_elev, n_fossil = [], 0, 0, 0
    for lo in loci:
        if not int(_f(lo.get("scored"))):
            continue
        n_scored += 1
        cs = state.get((lo["accession"], lo["cell"]), "not_in_matrix")
        d = dead_locus(lo, cs)
        n_elev += d["elevated"]
        n_fossil += d["is_fossil"]
        out.append(dict(accession=lo["accession"], organism=lo["organism"],
                        vclass=lo["vclass"], cell=lo["cell"],
                        contig=lo.get("contig", ""),
                        coverage=lo.get("coverage"),
                        identity=lo.get("identity"),
                        aligned_aa=lo.get("aligned_aa"),
                        frameshifts=lo.get("frameshifts"),
                        lesion_density=lo.get("lesion_density"),
                        verdict=lo.get("verdict"), cell_state=cs, **d))
    elevated = [r for r in out if r["elevated"]]
    stats = dict(
        lesion_bar=bar, n_loci=len(loci), n_scored=n_scored,
        n_above_bar=n_elev, n_fossils=n_fossil,
        n_above_bar_full_coverage=sum(1 for r in elevated
                                      if not r["incomplete"]),
        n_above_bar_zero_stops=sum(1 for r in elevated
                                   if r["stop_codons"] == 0),
        n_above_bar_live_cell=sum(1 for r in elevated
                                  if r["state_reading_available"]
                                  and not r["not_live"]),
        n_above_bar_state_unavailable=sum(
            1 for r in elevated if not r["state_reading_available"]),
        n_fossils_itpr=sum(1 for r in elevated
                           if r["is_fossil"] and r["cell"] != "RYR"),
        median_identity_above_bar=lib.median([_f(r["identity"])
                                              for r in elevated]),
        poisson_test="not run: it needs two dead loci sharing a lesion, "
                     "and there are no dead loci")
    return out, stats


def by_cell(rows: list[dict]) -> list[dict]:
    """Fossil readings tallied per cell, so the RyR control is visible."""
    out = []
    for cell in ("ITPR1", "ITPR2", "ITPR3", "RYR"):
        g = [r for r in rows if r["cell"] == cell]
        el = [r for r in g if r["elevated"]]
        out.append(dict(cell=cell, n_scored=len(g), n_above_bar=len(el),
                        n_fossils=sum(r["is_fossil"] for r in el),
                        frac_above_bar=(round(len(el) / len(g), 4)
                                        if g else ""),
                        median_identity=round(lib.median(
                            [_f(r["identity"]) for r in g]), 4) if g else "",
                        median_identity_above_bar=(round(lib.median(
                            [_f(r["identity"]) for r in el]), 4)
                            if el else "")))
    return out


# ---------------------------------------------------------------- D47's lead

MIN_CLASS_N = 8


def by_class(pairs: list[dict], matrix_rows: list[dict],
             matched: bool = True) -> tuple[list[dict], list[dict]]:
    """S15a's paired within-genome test, stratified by vertebrate class.

    Same pairs, same sign test, same identity matching — the only new
    thing is the stratum.  A class with fewer than `MIN_CLASS_N` untied
    pairs is reported with its n and no p-value rather than a p-value
    nobody should read; it stays out of the BH family for the same reason.
    """
    vclass = {r["accession"]: r["vclass"] for r in matrix_rows}
    use = [p for p in pairs if int(_f(p.get("matched"))) == int(matched)]
    rows = []
    for p in use:
        rows.append(dict(p, vclass=vclass.get(p["accession"], "unknown")))
    tests = []
    cells = sorted({r["cell"] for r in rows})
    classes = sorted({r["vclass"] for r in rows})
    for cell in cells:
        for cl in classes:
            g = [r for r in rows if r["cell"] == cell and r["vclass"] == cl]
            if not g:
                continue
            st = lib.sign_test([_f(r["diff"]) for r in g])
            tests.append(dict(
                cell=cell, vclass=cl, matched=int(matched),
                n_genomes=len(g), n=st["n"], n_pos=st["n_pos"],
                n_neg=st["n_neg"], n_ties=st["n_ties"],
                median_diff=round(lib.median([_f(r["diff"]) for r in g]), 5),
                direction=("excess" if st["n_pos"] > st["n_neg"] else
                           "deficit" if st["n_neg"] > st["n_pos"] else "tied"),
                p=(st["p"] if st["n"] >= MIN_CLASS_N else ""),
                underpowered=int(st["n"] < MIN_CLASS_N)))
    testable = [t for t in tests if t["p"] != ""]
    qs = lib.bh([t["p"] for t in testable])
    for t, q in zip(testable, qs):
        t["q_bh"] = q
    for t in tests:
        t.setdefault("q_bh", "")
    tests.sort(key=lambda t: (t["cell"], t["vclass"]))
    return rows, tests


def class_controls(pair_rows: list[dict], tests: list[dict],
                   matrix_rows: list[dict], q_cut: float = 0.05
                   ) -> list[dict]:
    """The two things that could produce a stratified result without a
    lineage effect, measured on the same committed pairs.

    **Contiguity.** A class whose assemblies are mostly below D4's bar
    could be reporting fragmentation, so each significant stratum is
    re-tested inside and outside the bar and both halves are printed.

    **The paired design itself.** `others` for a cell is the same genome's
    *other* family loci, so within one genome an excess at one paralog and
    a deficit at its siblings are one observation seen from two sides.
    `sibling_of` names the strata that are the other side of this one, so
    a report cannot present them as independent corroboration.
    """
    contig = {r["accession"]: int(_f(r["contig_spans_gene"]))
              for r in matrix_rows}
    out = []
    for t in tests:
        if t["p"] == "" or _f(t.get("q_bh"), 1.0) > q_cut:
            continue
        g = [r for r in pair_rows
             if r["cell"] == t["cell"] and r["vclass"] == t["vclass"]]
        above = [r for r in g if contig.get(r["accession"], 0)]
        below = [r for r in g if not contig.get(r["accession"], 0)]
        sa = lib.sign_test([_f(r["diff"]) for r in above])
        sb = lib.sign_test([_f(r["diff"]) for r in below])
        siblings = sorted({s["cell"] for s in tests
                           if s["vclass"] == t["vclass"] and s["p"] != ""
                           and _f(s.get("q_bh"), 1.0) <= q_cut
                           and s["cell"] != t["cell"]})
        out.append(dict(
            cell=t["cell"], vclass=t["vclass"], direction=t["direction"],
            n=t["n"], q_bh=t["q_bh"],
            n_above_bar=sa["n"], pos_above=sa["n_pos"], neg_above=sa["n_neg"],
            p_above=(sa["p"] if sa["n"] >= MIN_CLASS_N else ""),
            n_below_bar=sb["n"], pos_below=sb["n_pos"], neg_below=sb["n_neg"],
            p_below=(sb["p"] if sb["n"] >= MIN_CLASS_N else ""),
            median_identity=round(lib.median([_f(r["identity"])
                                              for r in g]), 4),
            median_others_identity=round(lib.median(
                [_f(r["others_identity_median"]) for r in g]), 4),
            sibling_of=";".join(siblings),
            note=("this stratum and its siblings are the same within-genome "
                  "comparison read from opposite sides"
                  if siblings else "")))
    return out
