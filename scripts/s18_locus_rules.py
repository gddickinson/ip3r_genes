"""How an annotation holds a gene — the five states, and the two verdicts.

The brief asks every gene-scale locus to be scored *complete / split /
fragmentary / noncoding / unannotated, same strand only*, with a name verdict
and a symbol verdict applied identically here and on the protein side. Three
things make that a measurement rather than a labelling exercise:

**Same strand only.** A gene on the other strand overlapping an ITPR locus is
not a model of this gene, however much sequence it shares with it. The rule
is a filter, not a tiebreak.

**Scored against CDS blocks, never gene spans.** A gene span covers its own
introns, so scoring against spans credits an annotation with every base it
never called — and with 152 kb introns in this family (S5's intron
calibration) a small passenger gene sitting *inside* an ITPR intron would
score as covering it. S10 established this at two loci; here it runs at 2,110.

**The completeness bar is inherited, and then measured.** This project already
has one definition of "this annotated gene is the model of that locus" —
`s5_classify.ANNOT_CDS_FRAC`, half the alignment's coding footprint — and S5,
S10 and S23 all read a locus through it. A second bar chosen here would be a
fork in what the project means by an annotated gene, so S18 uses that one and
spends its calibration on the question that is actually open: *is 0.50 the
right place for it?* `calibrate_complete()` answers that by measuring the
distribution the bar sits in — best-single-model coverage over loci whose own
annotation names the correct paralog with a coding model, whose alignment
recovered ≥ `CALIB_MIN_COVERAGE` of its bait and whose contig clears D4's bar.
The answer is a result rather than a formality (`s5_calibrate_margin`'s shape
applied to a second inherited constant), and `s18_locus_audit.sensitivity()`
re-counts every state across the whole bar range so what the choice buys is in
the data and not in a sentence.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_gff                                                    # noqa: E402
from s5_classify import (ANNOT_CDS_FRAC, ITPR_NAME_HINTS,     # noqa: E402
                         RYR_NAME_HINTS, name_family, name_paralog)
from s23_locus_evidence import is_placeholder                     # noqa: E402

# A gene has to reach this share of the locus's coding footprint to count as a
# *piece* of it at all; below it the overlap is a passenger, not a model.
MIN_PIECE_FRAC = 0.05
# Anything at all on the same strand: separates `unannotated` from `noncoding`.
MIN_TOUCH_FRAC = 0.01
# The calibration population's coverage admission.
CALIB_MIN_COVERAGE = 0.90
#: The bar, inherited: the project's one definition of "this annotated gene is
#: the model of that locus" (S5, S10, S23). S18 validates it, never replaces it.
COMPLETE_FRAC = ANNOT_CDS_FRAC
#: A locus above the bar but below this still leaves a fifth of the gene with
#: no model — reported as the price of the bar, not silently absorbed by it.
OVERCREDIT_FRAC = 0.80

STATES = ("cds_unavailable", "gff_unavailable", "no_gene_set", "unannotated",
          "noncoding", "complete", "split", "fragmentary")
NAME_VERDICTS = ("correct_paralog", "paralog_unspecified",
                 "correct_family_wrong_paralog", "wrong_family",
                 "family_ambiguous", "family_unnamed", "placeholder", "absent")


# --------------------------------------------------------------------------
# what the annotation puts over a locus
# --------------------------------------------------------------------------
def _product_of(gene: dict) -> str:
    """The descriptive name: the gene's own, else a transcript's product.

    RefSeq writes `description=` on the gene; a submitter-deposited GenBank
    gene set usually writes nothing there and puts the product on the CDS, so
    reading only the gene attribute would report most of the GenBank half of
    the scope as unnamed — which is a D9 difference manufactured by the
    parser rather than found in the data.
    """
    if gene.get("description"):
        return gene["description"]
    for tx in (gene.get("mrnas") or {}).values():
        if tx.get("product"):
            return tx["product"]
    return ""


def models_over(window: dict, locus_cds: list[tuple[int, int]], strand: str
                ) -> list[dict]:
    """Every same-strand annotated gene touching the locus's coding blocks.

    Coding and non-coding overlappers are both returned, tagged: a locus held
    only by a pseudogene is a different annotation failure from a locus held by
    nothing, and the deliverable is a correction list that has to say which.
    """
    total = sum(e - s + 1 for s, e in locus_cds)
    out = []
    for g in window.get("genes", []):
        if g["strand"] != strand:
            continue
        # A gene with real CDS features is a coding model whatever its
        # biotype string says, and a gene without them is not one however it
        # is labelled. *Podiceps*' ITPR3 is annotated `gene_biotype=other`,
        # `Note=contains frameshift`, named "Inositol 1,4,5-trisphosphate
        # receptor type 3", with no CDS anywhere — correctly identified and
        # serving no protein. Reading the biotype string would have called
        # that coding; reading the CDS calls it what a user of the database
        # would find.
        coding = bool(g["cds"]) and not g["pseudo"]
        blocks = g["cds"] if g["cds"] else g["exons"]
        if not blocks:
            blocks = [(g["start"], g["end"])]
        ov = s10_gff.overlap_bp(s10_gff.merge(blocks), locus_cds)
        if ov <= 0:
            continue
        out.append({
            "gene_id": g["gene_id"], "symbol": g["name"],
            "product": _product_of(g), "biotype": g["biotype"],
            "pseudo": bool(g["pseudo"]), "coding": coding,
            "start": g["start"], "end": g["end"], "strand": g["strand"],
            "overlap_bp": ov, "frac_cds": round(ov / total, 4) if total else 0.0,
            # The reciprocal direction: how much of the annotated model's own
            # coding sequence is inside this locus. A full-length gene next
            # door overlapping by 300 bp and a truncated model wholly inside
            # the gene have the same `frac_cds` and are different failures.
            "model_bp": sum(e - s2 + 1 for s2, e in s10_gff.merge(blocks)),
            "own_frac": round(ov / max(1, sum(e - s2 + 1 for s2, e
                                              in s10_gff.merge(blocks))), 4),
            "blocks": s10_gff.merge(blocks),
        })
    out.sort(key=lambda m: -m["overlap_bp"])
    return out


def naming_model(models: list[dict]) -> dict:
    """Which model the name verdicts are read off.

    Deliberately generous to the annotation, because the deliverable is a
    correction list and accusing a database of failing to name a gene it did
    name is the expensive error. So: the largest same-strand overlapper whose
    symbol *or* product claims either family; failing that the largest coding
    overlapper; failing that the largest overlapper of any kind.

    Reading the largest *coding* model first would have scored *Podiceps
    cristatus* ITPR3 `absent`, when its annotation names the gene exactly
    right and merely files it as non-coding — two different failures that a
    correction list has to keep apart.
    """
    named = [m for m in models
             if name_family(m.get("symbol", "")) or
             name_family(m.get("product", ""))]
    for pool in (named, [m for m in models if m["coding"]], models):
        if pool:
            return max(pool, key=lambda m: m["overlap_bp"])
    return {}


# --------------------------------------------------------------------------
# the state
# --------------------------------------------------------------------------
def locus_state(models: list[dict], locus_cds: list[tuple[int, int]],
                gene_set: str, complete_frac: float,
                min_piece: float = MIN_PIECE_FRAC) -> dict:
    """The five states as ordered positive tests; first to fire wins.

    Each state writes the number it fired on into its row, so a state is never
    a bare label: `fragmentary` beside `best_model_frac = 0.31` says what the
    annotation delivered, and a later reader can move the bar and see what
    moves with it (`s18_locus_audit.sensitivity`).
    """
    total = sum(e - s + 1 for s, e in locus_cds)
    coding = [m for m in models if m["coding"]]
    noncoding = [m for m in models if not m["coding"]]
    pieces = [m for m in coding if m["frac_cds"] >= min_piece]
    union = s10_gff.merge([b for m in coding for b in m["blocks"]])
    union_frac = round(s10_gff.overlap_bp(union, locus_cds) / total, 4) \
        if total else 0.0
    best = max((m["frac_cds"] for m in coding), default=0.0)

    d = {"n_models": len(models), "n_coding_models": len(coding),
         "n_coding_pieces": len(pieces), "n_noncoding_models": len(noncoding),
         "best_model_frac": round(best, 4), "union_coding_frac": union_frac,
         "locus_cds_bp": total,
         "noncoding_biotype": ";".join(sorted({m["biotype"] for m in noncoding}))
         if noncoding else ""}

    if gene_set == "cds_unavailable" or not locus_cds:
        return {**d, "state": "cds_unavailable",
                "state_reason": "the sweep's alignment for this locus could "
                                "not be recovered from the archived GFF — "
                                "unscorable, never scored against its span"}
    if gene_set == "unavailable":
        return {**d, "state": "gff_unavailable",
                "state_reason": "the manifest records a gene set this project "
                                "does not hold — not scorable either way"}
    if gene_set != "present":
        return {**d, "state": "no_gene_set",
                "state_reason": "the assembly ships no gene set (D9)"}
    if union_frac < MIN_TOUCH_FRAC and not noncoding:
        # Not "no same-strand feature": on 9 of 2,144 loci there *is* one and
        # it reaches under 1 % of the footprint — a passenger inside an
        # intron, not a model of this gene. A correction list is read by the
        # curator who owns the record, so the reason has to say which of the
        # two it is rather than assert the stronger one.
        near = (f"; the nearest same-strand coding model reaches "
                f"{best:.3f} of it" if coding else "")
        return {**d, "state": "unannotated",
                "state_reason": f"no same-strand coding model reaches "
                                f"{MIN_TOUCH_FRAC:g} of the locus's {total} "
                                f"coding bp{near}"}
    if union_frac < MIN_TOUCH_FRAC:
        return {**d, "state": "noncoding",
                "state_reason": f"{len(noncoding)} same-strand non-coding "
                                f"feature(s) only ({d['noncoding_biotype']}); "
                                f"coding models reach {union_frac:.3f} of the "
                                f"footprint"}
    if best >= complete_frac:
        return {**d, "state": "complete",
                "state_reason": f"one model covers {best:.2f} of the coding "
                                f"footprint (bar {complete_frac:.2f})"}
    if len(pieces) >= 2 and union_frac >= complete_frac:
        return {**d, "state": "split",
                "state_reason": f"{len(pieces)} models cover {union_frac:.2f} "
                                f"together, best alone {best:.2f}"}
    return {**d, "state": "fragmentary",
            "state_reason": f"coding models cover {union_frac:.2f} of the "
                            f"footprint, best alone {best:.2f}"}


# --------------------------------------------------------------------------
# the two verdicts — one rule, both sides of the audit
# --------------------------------------------------------------------------
def _names_both(low: str) -> bool:
    """A name that puts the two families on either side of a slash.

    `ITPR_NAME_HINTS` are substrings of the *IP3 receptor's* own names, and
    the superfamily names UniProt actually uses for non-vertebrate records
    break every one of them: "Inositol 1,4,5-trisphosphate/ryanodine receptor
    domain-containing protein" has no "…trisphosphate receptor" in it, because
    the word "receptor" is shared. Matching on the inositol half separately is
    what makes 15 records read as the superfamily name they are rather than as
    the database calling an IP3 receptor a ryanodine receptor.
    """
    if not any(h in low for h in RYR_NAME_HINTS):
        return False
    return "inositol" in low and ("trisphosphate" in low
                                  or "triphosphate" in low
                                  or "ip3r" in low or "insp3" in low)


def name_verdict(text: str, expected_cell: str, family: str) -> tuple[str, str]:
    """What a name or a symbol claims, against what the evidence calls it.

    Applied unchanged to a genome gene's symbol, a genome gene's product name,
    a UniProt `gene` field and a UniProt protein name — the brief's "applied
    identically", which is what makes the genome and protein halves of the
    audit comparable at all. `expected_cell` is the paralog where one is
    defined (the vertebrate ITPR1/2/3 and RYR1/2/3 cells) and `family` the
    family call; outside the vertebrates the paralog question is not asked and
    `expected_cell` is empty.

    Two verdicts exist because the first version manufactured errors without
    them, and both are annotation states in their own right:

    **`family_ambiguous`.** UniProt's commonest name for a non-vertebrate
    family record is "RyR/IP3R Homology associated domain-containing protein",
    and "Inositol 1,4,5-trisphosphate/ryanodine receptor, putative" is close
    behind. Those name *both* families. `name_family` resolves them to RYR
    because it tests the RyR hints first, which is right for its own purpose
    and would score 66 records here as the database calling an IP3 receptor a
    ryanodine receptor — a wrong-family error rate invented by a precedence
    rule.

    **`paralog_unspecified`.** "Inositol 1,4,5-trisphosphate receptor" with no
    type number is not a wrong paralog; it is a correct family name that
    declines the paralog. Folding it into `correct_family_wrong_paralog` put
    6,660 records in the wrong-paralog cell, which is most of the vertebrate
    scope and would have been the audit's headline.
    """
    if not (text or "").strip():
        return "absent", ""
    low = text.lower()
    hits_itpr = any(h in low for h in ITPR_NAME_HINTS)
    hits_ryr = any(h in low for h in RYR_NAME_HINTS)
    if (hits_itpr and hits_ryr) or _names_both(low):
        return "family_ambiguous", "ITPR/RYR"
    fam = name_family(text)
    par = name_paralog(text)
    if not fam:
        if is_placeholder(text):
            return "placeholder", ""
        return "family_unnamed", ""
    if fam != family:
        return "wrong_family", par or fam
    if not expected_cell or expected_cell == family:
        return "correct_paralog", par or fam
    if not par:
        return "paralog_unspecified", fam
    return (("correct_paralog" if par == expected_cell
             else "correct_family_wrong_paralog"), par)


# --------------------------------------------------------------------------
# the bar
# --------------------------------------------------------------------------
def calibration_population(rows: list[dict]) -> list[dict]:
    """Loci the annotation demonstrably holds, in assemblies that carry them.

    Four admissions, each a positive test the coverage number did not produce:
    a same-strand model whose *name or symbol* names the locus's own paralog,
    that model being a **coding** one, an alignment recovering ≥ 90 % of its
    bait, and D4's contiguity bar. Everything else — including every locus
    whose annotation is silent — stays out.

    The coding admission is not bookkeeping. Without it the population's lower
    tail runs to 0.0, because a gene the annotation names perfectly and files
    as `gene_biotype=other` with no CDS anywhere covers none of itself (the
    *Podiceps cristatus* ITPR3 case) — and those are the failures the bar
    exists to detect, not evidence about where it belongs.

    What survives is still not a clean population: 5 % of it sits below 0.91
    and a thin tail runs to 0.002, which is what a correctly-*named* locus
    whose model is a neighbouring gene's looks like. That is reported rather
    than filtered away, because it is the reason the bar is validated here and
    not derived here.
    """
    return [r for r in rows
            if r.get("state") not in ("no_gene_set", "gff_unavailable")
            and float(r.get("coverage") or 0) >= CALIB_MIN_COVERAGE
            and str(r.get("contig_spans_gene")) in ("1", "True", "true")
            and str(r.get("namer_coding")) in ("1", "True", "true")
            and (r.get("symbol_verdict") == "correct_paralog"
                 or r.get("product_verdict") == "correct_paralog")]


def calibrate_complete(rows: list[dict], quantile_fn,
                       bar: float = COMPLETE_FRAC) -> dict:
    """Where the inherited bar sits in the distribution it is applied to.

    Reports the distribution, the bar's percentile in it, and the two error
    rates the choice buys on that population: how many demonstrably-annotated
    loci the bar calls *incomplete* (its false-fragmentary rate), and how many
    it calls complete while a full 20 % of the gene reaches no model
    (`OVERCREDIT_FRAC` — the ones a correction list would still want to see).

    Refuses to report from fewer than `MIN_CALIB` loci: a threshold check
    measured on a handful is the handful's number, and this project has been
    caught by exactly that before — `s23_calibrate_loci` wrote 0.25 from two
    loci and the next sweep read it back.
    """
    MIN_CALIB = 50
    pop = calibration_population(rows)
    vals = sorted(float(r["best_model_frac"]) for r in pop)
    if len(vals) < MIN_CALIB:
        return {"usable": 0, "n": len(vals), "bar": bar,
                "reason": f"only {len(vals)} calibration loci (< {MIN_CALIB})"}
    below = sum(1 for v in vals if v < bar)
    overcredited = sum(1 for v in vals if bar <= v < OVERCREDIT_FRAC)
    pct = sum(1 for v in vals if v <= bar) / len(vals)
    return {
        "usable": 1, "n": len(vals), "bar": bar,
        "bar_source": "s5_classify.ANNOT_CDS_FRAC (S5, S10, S23)",
        "bar_percentile": round(pct, 4),
        "n_below_bar": below, "frac_below_bar": round(below / len(vals), 4),
        "n_overcredited": overcredited,
        "frac_overcredited": round(overcredited / len(vals), 4),
        "overcredit_frac": OVERCREDIT_FRAC,
        "p01": round(quantile_fn(vals, 0.01), 4),
        "p05": round(quantile_fn(vals, 0.05), 4),
        "p10": round(quantile_fn(vals, 0.10), 4),
        "p25": round(quantile_fn(vals, 0.25), 4),
        "median": round(quantile_fn(vals, 0.50), 4),
        "min": round(min(vals), 4), "max": round(max(vals), 4),
        "reason": (f"best-single-model coverage over {len(vals)} "
                   f"correctly-named, fully-recovered loci; the inherited bar "
                   f"{bar:g} is that distribution's {pct:.1%} point"),
    }
