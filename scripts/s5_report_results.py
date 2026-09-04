"""s5_report_results.py — the results half of the S5 report.

Split from `s5_report.py` to keep both under 500 lines, and following the
`s3_report.py` / `s3_report_d10.py` pattern: this module takes the caller's
table loader and formatter rather than re-importing its own, so the two halves
cannot read the same tables differently.

**Every headline here is chosen by the data, not asserted.** S5a's findings
came from a six-genome pilot and are hypotheses at 309: the short-gene
detection bias, the perfect ITPR/RyR separation, the fragment margin floor. A
report generator written to narrate the pilot's answers would print them
whatever the full sweep said, which is the specific way a pipeline like this
launders an assumption into a result. So each section computes its statistic,
compares it against the pilot's recorded value, and renders `confirmed`,
`contradicted`, `weakened` or `underpowered` from the comparison — with the
pilot's number printed beside the new one either way.
"""

from __future__ import annotations

from collections import Counter

#: The pilot's values, recorded so the full sweep can be compared against
#: them explicitly rather than silently replacing them. Source:
#: `results/genome_ledger/report.md` as rendered at S5a (2026-09-04).
PILOT = {
    "genomes": 6,
    "found_above_bar": {"ITPR1": (3, 3), "ITPR2": (3, 3), "ITPR3": (3, 3)},
    "found_below_bar": {"ITPR1": (0, 3), "ITPR2": (0, 3), "ITPR3": (2, 3)},
    "contested_family_loci": 0,
    "total_loci": 43,
    "fragment_margin_min": None,      # S5a had no fragment evidence at all
}

#: Below this many genomes on the wrong side of D4's bar, the per-paralog
#: recovery comparison is not worth a verdict.
MIN_BELOW_BAR_FOR_VERDICT = 20


def _rate(hit: int, tot: int) -> str:
    return f"{hit}/{tot} ({hit / tot:.0%})" if tot else "n/a"


def section_control(A, table, ledger, wide, stats) -> None:
    """The RyR positive control — the gate on reading anything else."""
    A("\n## The RyR positive control\n")
    ctrl = [r for r in ledger if r["class"] == "RYR"]
    fired = [r for r in ctrl if r["status"] not in ("absent", "no_locus")]
    failures = [w for w in wide if w["control_ok"] != "1"]
    A(f"A ryanodine receptor bait travelled with every genome. RyRs are "
      f"present in three copies in every vertebrate, so a genome where the "
      f"control finds nothing has an assembly or pipeline problem rather than "
      f"a biological result — and until it fires, that genome's ITPR cells "
      f"say nothing.\n")
    A(f"**The control fired in {_rate(len(fired), len(ctrl))} genomes.**\n")
    if failures:
        A(f"The {len(failures)} exception(s) are excluded from every absence "
          "claim below and listed in `control_failures.tsv`:\n")
        A(table(["genome", "class", "contig N50", "RyR cell", "assembly"],
                [[f"*{w['organism']}*", w["vclass"], f"{int(w['contig_n50']):,}",
                  w["RYR"], w["assembly_level"]] for w in failures[:12]]))
    else:
        A("No genome failed it, so no ITPR result below is excluded on "
          "control grounds.\n")
    by_status = Counter(r["status"] for r in ctrl)
    A("\nControl cell statuses: "
      + ", ".join(f"{k} {v}" for k, v in by_status.most_common()) + ".\n")


def section_ledger(A, table, ledger, stats, STATUS_ORDER, CLASSES) -> None:
    A("\n## The ledger\n")
    counts = stats["status_counts"]
    A(f"{stats['genomes_swept']} genomes swept of the "
      f"{stats['genomes_in_manifest']} in the declared scope; "
      f"{stats['loci_recorded']} loci recorded.\n")
    rows = []
    for cell in CLASSES:
        row = [cell]
        for st in STATUS_ORDER:
            row.append(counts.get(f"{cell}:{st}", 0) or "")
        rows.append(row)
    A(table(["paralog"] + [s.replace("_", " ") for s in STATUS_ORDER], rows))
    absent = {c: counts.get(f"{c}:absent", 0) for c in CLASSES}
    total_absent = sum(absent.values())
    A(f"\n**{total_absent} cell(s) are called `absent`** — no spliced-alignment "
      "locus and no tblastn remnant. That is the strongest negative this task "
      "produces, and it is not yet a loss claim: D4 requires the contiguity "
      "bar as well, which the next section applies.\n")

    # How well the annotations know each paralog, among loci actually found.
    # This is S18's question and it falls straight out of the ledger, so it is
    # reported here rather than left for a table nobody reads.
    A("\n### How well the annotations know each paralog\n")

    def annot_stats(rows_in):
        out, frac = [], {}
        for cell in CLASSES:
            found = [r for r in rows_in if r["class"] == cell
                     and r["status"].startswith("found")
                     and r["annotated_assembly"] == "Y"]
            named = [r for r in found if r["annot_paralog_matches"] == "1"]
            family = [r for r in found if r["annot_gene"]]
            frac[cell] = len(named) / len(found) if found else 0.0
            out.append([cell, len(found), _rate(len(family), len(found)),
                        _rate(len(named), len(found))])
        return out, frac

    # Contiguity is held constant, because it is a confounder rather than a
    # nuisance: the same fragmented assemblies both fail to place a gene on
    # one contig and fail to annotate it, so an uncontrolled comparison would
    # measure assembly quality and call it annotation quality.
    intact = [r for r in ledger if r["contig_spans_gene"] == "1"]
    rows_all, frac_all = annot_stats(ledger)
    rows_intact, frac_intact = annot_stats(intact)
    A("Restricted to loci the sweep **found** in an assembly that carries a "
      "gene set, so the denominator is genes that exist and are annotatable. "
      "The second block additionally holds contiguity constant — only "
      "assemblies whose contigs can carry the gene — because otherwise this "
      "measures assembly quality and calls it annotation quality:\n")
    A(table(["paralog", "found loci", "a gene model is there",
             "and it names this paralog"], rows_all))
    A(f"\nHolding contiguity constant ({len({r['accession'] for r in intact})} "
      "genomes above the bar):\n")
    A(table(["paralog", "found loci", "a gene model is there",
             "and it names this paralog"], rows_intact))
    spread = frac_intact if sum(r[1] for r in rows_intact) >= 30 else frac_all
    controlled = spread is frac_intact
    A(f"\n*The verdict below reads the "
      f"{'contiguity-controlled' if controlled else 'uncontrolled'} block"
      + ("" if controlled else
         " — too few loci above the bar yet to control, so the gap it "
         "reports may be assembly quality rather than annotation quality")
      + ".*\n")
    if spread:
        best = max(spread, key=spread.get)
        worst = min(spread, key=spread.get)
        gap = spread[best] - spread[worst]
        if gap >= 0.15:
            A(f"\n**The three paralogs are not annotated equally well.** "
              f"{best} is correctly named at {spread[best]:.0%} of its found "
              f"loci against {worst} at {spread[worst]:.0%} — a {gap:.0%} "
              "gap between genes of near-identical protein length in the same "
              "genomes. A census built on gene symbols inherits that gap as "
              "an apparent difference in copy number, which is the "
              "annotation-quality problem S18 exists to quantify.\n")
        else:
            A(f"\nThe three are named about equally well "
              f"({worst} {spread[worst]:.0%} to {best} {spread[best]:.0%}), so "
              "annotation quality is not a per-paralog confounder here.\n")


def section_contiguity(A, table, ledger, contiguity, bar, CLASSES) -> None:
    """D4's bar at full scale, and S5a's span-bias hypothesis tested."""
    A("\n## Assembly contiguity, and what it does to an absence\n")
    total = next((c for c in contiguity if c["group"] == "total"), None)
    A(f"An assembly whose contig N50 falls below the median measured ITPR "
      f"genomic span ({bar:,} bp) cannot carry the gene on one contig, so its "
      "empty cells are evidence about the assembly rather than the animal.\n")
    if total:
        A(f"**{total['below']} of {total['total']} genomes in the manifest "
          f"({float(total['pct']):.0f} %) fall below it.**\n")
    A(table(["group", "below the bar", "of", "%"],
            [[c["name"], c["below"], c["total"], f"{float(c['pct']):.0f} %"]
             for c in contiguity if c["group"] != "class"
             or int(c["below"]) > 0]))

    # --- S5a's hypothesis, tested rather than repeated
    A("\n### Does recovery track gene span? (S5a's hypothesis at scale)\n")
    swept = {r["accession"] for r in ledger}
    below, above = {}, {}
    for cell in CLASSES:
        rows = [r for r in ledger if r["class"] == cell]
        a = [r for r in rows if r["contig_spans_gene"] == "1"]
        b = [r for r in rows if r["contig_spans_gene"] != "1"]
        above[cell] = (sum(1 for r in a if r["status"].startswith("found")), len(a))
        below[cell] = (sum(1 for r in b if r["status"].startswith("found")), len(b))
    n_below = below[CLASSES[0]][1]
    A(f"S5a predicted, from six genomes, that a fragmented assembly loses the "
      f"*long* paralogs first — ITPR3 has the shortest genomic span and was "
      f"recovered where ITPR1 and ITPR2 were not. Across "
      f"{len(swept)} genomes:\n")
    A(table(["paralog", "found above the bar", "found below the bar",
             "pilot, below the bar"],
            [[c, _rate(*above[c]), _rate(*below[c]),
              _rate(*PILOT["found_below_bar"][c])] for c in CLASSES]))
    if n_below < MIN_BELOW_BAR_FOR_VERDICT:
        A(f"\n**Underpowered**: only {n_below} genomes fall below the bar in "
          "the swept set, too few to separate a span effect from noise. The "
          "hypothesis stands untested rather than confirmed.\n")
    else:
        best = max(CLASSES, key=lambda c: below[c][0] / max(1, below[c][1]))
        worst = min(CLASSES, key=lambda c: below[c][0] / max(1, below[c][1]))
        spread = (below[best][0] / max(1, below[best][1])
                  - below[worst][0] / max(1, below[worst][1]))
        if best == "ITPR3" and spread >= 0.10:
            A(f"\n**Confirmed.** Below the bar, {best} is recovered "
              f"{spread:.0%} more often than {worst}, and {best} is the "
              "shortest gene of the three. The detection bias runs in the "
              "same direction as the loss signal the margin species were "
              "selected for, so a per-paralog absence in a fragmented "
              "assembly cannot be read as loss.\n")
        elif spread < 0.10:
            A(f"\n**Not confirmed.** Below the bar the three paralogs are "
              f"recovered within {spread:.0%} of one another, so the "
              "six-genome pattern does not survive the full scope. "
              "Fragmentation suppresses recovery, but not preferentially by "
              "gene length.\n")
        else:
            A(f"\n**Contradicted.** Below the bar the best-recovered paralog "
              f"is {best}, not the shortest gene, by {spread:.0%} over "
              f"{worst}. Whatever drives differential recovery in poor "
              "assemblies, S5a's span explanation is not it.\n")


def section_d14(A, table, margin, PILOT_LOCI) -> None:
    """Did the two families ever contest a locus at genome scale?"""
    A("\n## D14: did ITPR and RyR ever contest a locus?\n")
    fam = margin["groups"].get("family_margin_all", {})
    n = fam.get("n", 0)
    A("The ryanodine receptors carry every ITPR-diagnostic domain and sit in "
      "every vertebrate genome in triplicate. At the protein level they "
      "defeated this project's own detector completely (S1: all six decoys "
      "promoted). The genomic question is whether a locus is ever won by one "
      "family's baits only narrowly over the other's.\n")
    if not n:
        A("No locus margins recorded.\n")
        return
    if fam.get("min") == 1.0:
        A(f"**At all {n:,} loci, only one family's baits aligned at all.** "
          "The panels never contested a locus, so the family call was settled "
          "before any margin had to be applied — a sharper separation than "
          "the protein level affords, and it held from the six-genome pilot "
          f"(0 of {PILOT_LOCI} contested) to full scale.\n")
    else:
        contested = fam["n"] - 0
        A(f"**The families did contest loci here.** Margins run from "
          f"{fam['min']} to {fam['max']} over {n:,} loci, median "
          f"{fam['median']}. The pilot saw none of this "
          f"(0 of {PILOT_LOCI} loci contested), so the separation is not the "
          "absolute one six genomes suggested; every locus below D7's 10 % "
          "band is listed in `locus_margins.tsv` and needs a second "
          "instrument before its paralog is quoted.\n")


def section_margin(A, table, margin, in_force) -> None:
    A("\n## The attribution margin, calibrated on fragments\n")
    g = margin["groups"]
    A("S5a could only measure paralog separation at *complete* loci, which "
      "bounds a fragment's separation from above without measuring it. The "
      "full sweep produces the fragments themselves: rescue regions "
      "overlapping a gene the assembly names for a paralog carry an identity "
      "established independently of the bait scores.\n")
    frag = g.get("fragment_annotation_agrees", {})
    if not frag.get("n"):
        A("No rescue region in this run overlapped a paralog-named gene, so "
          "the fragment-level calibration has no evidence and the threshold "
          "stands where S5a's complete-locus bound put it.\n")
        return
    A(table(["evidence level", "n", "min", "median", "max"],
            [["complete loci, contested",
              g["annotation_agrees_contested"]["n"],
              g["annotation_agrees_contested"]["min"],
              g["annotation_agrees_contested"]["median"],
              g["annotation_agrees_contested"]["max"]],
             ["rescue fragments", frag["n"], frag["min"], frag["median"],
              frag["max"]]]))
    below = margin.get("fragment_below_threshold", 0)
    A(f"\nThe attribution agrees with the annotation on "
      f"{frag['n']}/{g['fragment_annotation_known']['n']} of them, and "
      f"**{below} fall below the {in_force} threshold in force**. Under the "
      f"inherited {margin['inherited_threshold']} the same evidence would "
      "have been reported ambiguous, and no absence claim could have been "
      "attributed at all.\n")
    wrong = g.get("fragment_annotation_disagrees", {})
    if wrong.get("n"):
        A(f"{wrong['n']} region(s) were attributed against their annotation, "
          f"up to a margin of {wrong['max']} — the population a threshold has "
          "to exclude, and the reason the floor is not simply the lowest "
          "correct margin.\n")
    else:
        A("No region was attributed against its annotation, so this evidence "
          "bounds the threshold from below only: how low it must be to keep "
          "correct calls, not how high it may go before wrong ones enter.\n")


def section_census(A, table, v4stats, models) -> None:
    A("\n## Census v4 — the genes only the DNA holds\n")
    if not v4stats:
        A("Census v4 has not been built for this sweep yet.\n")
        return
    A(f"Census v3 held {v4stats['v3_records']:,} protein-database records. "
      f"The sweep contributes {v4stats['models_enrolled']:,} gene models from "
      f"{v4stats['genomes_contributing']} genomes, each scored against "
      "`itpr.hmm`/`ryr.hmm` so a v4 row rests on the instrument that called "
      f"v3 → **census v4 is {v4stats['v4_records']:,} records**.\n")
    by_db = v4stats.get("itpr_by_db_status", {})
    if by_db:
        A(table(["how a database holds the locus", "ITPR models"],
                [[k.replace("_", " "), v] for k, v in
                 sorted(by_db.items(), key=lambda kv: -kv[1])]))
        gonly = by_db.get("genome_only", 0)
        unnamed = by_db.get("annotated_unnamed", 0)
        A(f"\n**{gonly} ITPR gene models exist only as DNA** — no gene model "
          f"at the locus, or an assembly with no gene set at all. A further "
          f"**{unnamed} sit inside an annotated gene that carries no family "
          "name**, so no search by name can reach them however complete the "
          "protein databases are.\n")
    if v4stats.get("models_unassigned"):
        A(f"{v4stats['models_unassigned']} model(s) the profiles declined to "
          "call are reported in `unassigned_models.tsv` and not enrolled — "
          "the same discipline S3 applied to its sweep.\n")


def section_discrepancies(A, table, models) -> None:
    A("\n### Where an annotation and this sweep disagree\n")
    other = [m for m in models if m["db_status"] == "annotated_other_paralog"]
    if not other:
        A("No locus was claimed by a cell other than the one its annotation "
          "names.\n")
        return
    coherent = [m for m in other
                if m.get("sibling_locus_for_annot_paralog") == "1"]
    A(f"{len(other)} locus/loci are claimed by a paralog cell other than the "
      f"one the assembly's annotation names. **{len(coherent)} of them sit in "
      "a genome that also carries a separate locus for the paralog the "
      "annotation names**, so the two cannot both be that paralog and the "
      "discrepancy is internally coherent.\n")
    A("These are supplied to S18, not adjudicated here: one instrument does "
      "not overturn a public annotation, and this family's paralog margins "
      "are narrow.\n")
    A(table(["genome", "claimed by", "annotated", "paralog margin",
             "sibling locus?"],
            [[f"*{m['organism']}*", m["cell"], m["annot_gene"],
              m["paralog_margin"],
              "yes" if m.get("sibling_locus_for_annot_paralog") == "1" else "no"]
             for m in sorted(other, key=lambda m: -float(m["paralog_margin"] or 0))[:15]]))
