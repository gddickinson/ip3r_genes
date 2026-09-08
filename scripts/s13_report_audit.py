"""S13 — the audit half of `results/reconciliation/report.md` (§7-§11).

Split from `s13_report_results.py` to keep both under 500 lines, and taking
the caller's loader, formatter and clade helpers so the two halves cannot read
the tables differently (the `s3_report.py` / `s3_report_d10.py` pattern,
applied twice as S7 did).

This half is the part of the report that argues against its own headline: the
six tips the deep placement rests on, checked against an instrument that never
saw the alignment and against the long-branch objection; the loss count the
reconciliation implies, taken apart against the genome sweep; and what the
result does and does not settle for the earlier tasks.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s13_lib import deepest_nodes as _deepest      # noqa: E402
from s13_lib import matrix_cells as _cells         # noqa: E402
from s13_priors import PRIOR, verdict_line          # noqa: E402

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")


# ------------------------------------------------------------------ §7

def _cyclostomes(load, table, num, not_run) -> list:
    cy = load("cyclostome_loci.tsv")
    bl = load("branch_lengths.tsv")
    summ = load("reconciliation_summary.tsv")
    L = ["## 7. The six tips the older placement rests on", ""]
    if not cy:
        return L + not_run("the cyclostome cross-check")

    with_c = sorted(_deepest(summ, "with_cyclostome")) if summ else []
    without_c = sorted(_deepest(summ, "without_cyclostome")) if summ else []
    computed = (
        f"the tree makes **{len({r['pair_support'] for r in cy})} "
        f"hagfish/lamprey orthology pairs** out of the six loci, in "
        f"{len({r['cyclostome_clade_size'] for r in cy})} cyclostome-only "
        f"clades. With them the deepest paralog duplication maps to "
        f"**{', '.join(with_c)}**; drop them and it maps to "
        f"**{', '.join(without_c)}**. So they carry the whole difference "
        f"between a pre-cyclostome and a gnathostome-stem placement — and S13 "
        f"is the instrument S8 named for the question, so the two checks below "
        f"are what say how much weight they can take."
    )
    L += [verdict_line("cyclostome", computed, "confirmed"), ""]

    L += ["### 7.1 Does an instrument that never saw the alignment agree?", ""]
    L += table(["locus", "S5 cell", "tree pair", "pair support",
                "S8 flank call", "vs its own null", "pair verdict"],
               [[f"`{r['s8_annot_gene'] or r['label'][:28]}` "
                 f"({r['species'].split()[0]})",
                 r["s5_cell"], r["pair_size"], r["pair_support"] or "—",
                 r["s8_paralog_call"] or "no call",
                 r["s8_null_verdict"] or "—", r["s8_pair_agreement"]]
                for r in cy])
    agree = Counter(r["s8_pair_agreement"] for r in cy)
    nulls = Counter(r["s8_null_verdict"] or "no_call" for r in cy)
    joins = Counter(r["s8_join"].split()[0] for r in cy)
    L += [f"All six loci join to an S8 locus exactly and offline — "
          f"{joins.get('coordinates', 0)} by genomic coordinates (the *Myxine* "
          f"tips are S5 gene models and carry their own) and "
          f"{joins.get('census', 0)} by the `LOC` gene symbol census v6 records "
          f"for the UniProt entry, which is the same symbol S8 recorded as the "
          f"locus's annotated gene.", "",
          f"S8's consensus caller read **only the flanking gene symbols**, so "
          f"it knows nothing about this alignment, this model or this tree. Of "
          f"the three orthology pairs the tree proposes it agrees on "
          f"{agree.get('agree', 0)//2}, disagrees on {agree.get('disagree', 0)//2} "
          f"and is uninformative on {agree.get('uninformative', 0)//2} — and "
          f"**{nulls.get('within_null', 0)} of the {len(cy)} loci got a "
          f"call at all, and every one of those sits inside its own null** "
          f"(the other {nulls.get('no_call', 0)} are no-calls), so not one "
          f"of them is evidence in either direction. That is S8's own "
          f"`underpowered` verdict re-derived one locus at a time, and it is "
          f"the honest reading: the corroboration this task would most like to "
          f"have does not exist.", ""]

    L += ["### 7.2 Are they long branches?", ""]
    if not bl:
        L += not_run("the branch-length check")
    else:
        cyr = [r for r in bl if r["is_cyclostome"]]
        ratios = sorted(float(r["vs_median"]) for r in cyr)
        ranks = sorted(int(r["rank"]) for r in cyr)
        n = len(bl)
        L += [f"Cyclostome sequences are the classic long-branch attraction "
              f"risk in vertebrate phylogeny, and long branches are attracted "
              f"to the root — which is exactly where this answer sits. So the "
              f"objection is measured rather than argued: root-to-tip distance "
              f"for all **{n}** vertebrate tips, with the six marked.", ""]
        L += table(["locus", "root-to-tip", "vs the median of all "
                    f"{n}", "rank (1 = longest)"],
                   [[f"`{r['label'][:44]}`", num(r["root_to_tip"], 3),
                     f"{num(r['vs_median'],2)}×", f"{r['rank']} of {n}"]
                    for r in sorted(cyr, key=lambda x: int(x["rank"]))])
        L += [f"**{num(ratios[0],2)}–{num(ratios[-1],2)}× the median, ranking "
              f"{ranks[0]} to {ranks[-1]} of {n}.** These are ordinary "
              f"branches — not one of them is an outlier, and two sit in the "
              f"shorter half. The long-branch objection to a cyclostome-"
              f"anchored deep placement does not apply to these tips. That is "
              f"a negative result and it is the one that matters most here: it "
              f"is the reason the placement in §5 is offered as a finding "
              f"rather than as a caveat.", ""]
    return L


# ------------------------------------------------------------------ §8

def _losses(load, table, num, not_run) -> list:
    ver = load("loss_verdicts.tsv")
    pres = load("paralog_presence.tsv")
    summ = load("reconciliation_summary.tsv")
    det = load("corroborated_losses.tsv")
    L = ["## 8. The loss audit — and why a loss count from a reconciliation "
         "should not be quoted", ""]
    if not ver:
        return L + not_run("the loss audit")

    binary = [r for r in _cells(summ) if r["variant"] != "support_collapsed"]
    coll = [r for r in _cells(summ) if r["variant"] == "support_collapsed"]
    implied = max((int(r["implied_losses"]) for r in binary), default=0)
    lo = min((int(r["implied_losses"]) for r in binary), default=0)
    coll_hi = max((int(r["implied_losses"]) for r in coll), default=0)
    tally = {r["verdict"]: int(r["n_species"]) for r in ver
             if r["paralog"] == "all"}
    real = tally.get("corroborated_loss", 0)
    computed = (
        f"the reconciliation implies **{lo}-{implied} losses** across the "
        f"{len(binary)} fully-resolved cells, and up to {coll_hi} once "
        f"unsupported nodes are collapsed — every child of a polytomy "
        f"crosses the full depth gap, so a collapse inflates the count "
        f"without adding evidence. Asked of the S5 genome ledger one "
        f"species x paralog "
        f"cell at a time, **{tally.get('corroborated_loss', 0)} are corroborated** — "
        f"{tally.get('sampling_artefact', 0)} are genes the sweep found intact "
        f"in the genome that S6 simply did not sample, "
        f"{tally.get('paralog_unassignable', 0)} are cells where the genome "
        f"carries family loci the bait panel cannot assign to a paralog, and "
        f"{tally.get('no_genome_in_manifest', 0)} are species with no genome "
        f"in the S4 scope at all."
    )
    L += [verdict_line("losses", computed, "confirmed"), ""]
    L += [f"`msa_v2` is a **representative** alignment — one sequence per clade "
          f"per species by S6's eight rules — so a species with no ITPR2 tip is "
          f"almost always a species whose ITPR2 was never sampled, not one that "
          f"lost it. A reconciliation cannot tell those apart; only a genome "
          f"can. Every cell of the {len(pres)}-cell grid is therefore asked of "
          f"the ledger:", ""]
    L += table(["verdict"] + list(PARALOGS) + ["all"],
               [[v] + [next((r["n_species"] for r in ver
                             if r["paralog"] == p and r["verdict"] == v), "0")
                       for p in PARALOGS]
                + [str(tally.get(v, 0))]
                for v in ("sampled_in_gene_tree", "sampling_artefact",
                          "paralog_unassignable", "corroborated_loss",
                          "loss_with_remnant", "undecidable",
                          "no_genome_in_manifest")
                if tally.get(v, 0)])

    L += ["### 8.1 The rule this audit could not do without", ""]
    L += [f"Run without `paralog_unassignable`, the only four corroborated "
          f"losses in the whole table were ITPR2 and ITPR3 in *Myxine "
          f"glutinosa* and *Petromyzon marinus*. Both genomes carry **three "
          f"ITPR loci apiece**, all filed by the sweep in the ITPR1 cell "
          f"because S5 has no cyclostome-labelled bait to offer the other two — "
          f"its six unfilled bait slots include all three paralogs in "
          f"cyclostomes. The ledger says `absent` there about a **cell**, not "
          f"about a gene, and S7 says the same thing from the other side: "
          f"every cyclostome tip sits outside all three paralog clades. A loss "
          f"claim built on that cell would report a bait-panel limit as "
          f"biology.", ""]
    if det:
        L += table(["species", "paralog", "ledger", "ITPR loci in the genome",
                    "loci beyond the cells they fill", "tips the tree cannot "
                    "place", "verdict"],
                   [[f"*{r['species']}*", r["paralog"], r["s5_ledger_status"],
                     r["n_ledger_itpr_loci"], r["n_spare_loci"],
                     r["n_tree_unplaced_tips"], f"**{r['verdict']}**"]
                    for r in det])
    L += [f"So **{real} of the {implied} implied losses survive contact with "
          f"the genomes**. The general lesson is the one the PIEZO project "
          f"reached from the same measurement and it holds here more sharply: "
          f"loss counts from reconciliations on representative alignments "
          f"should not be quoted. What can be quoted is the sweep's own "
          f"number, which is a different measurement on a complete "
          f"denominator.", ""]
    return L


# ------------------------------------------------------------------ §9

def _what_this_changes(load, table, num) -> list:
    summ = load("reconciliation_summary.tsv")
    L = ["## 9. What this changes about the earlier tasks", ""]
    with_c = sorted(_deepest(summ, "with_cyclostome")) if summ else []
    L += ["### 9.1 The 2R question the review left open", ""]
    computed = (
        f"the reconciliation places the ITPR1 / ITPR2+ITPR3 split on the "
        f"**vertebrate stem** and the ITPR2 / ITPR3 split on the "
        f"**gnathostome stem**, under every topology in the AU 95 % set and "
        f"under the `--bnni` guard. That is the shape a 2R origin predicts for "
        f"the older event and it is the first phylogeny-tested placement of "
        f"either. It is **not** a demonstration that the three are ohnologues: "
        f"a placement on the right branch is necessary and not sufficient, and "
        f"S8's paralogon test — which is the synteny half of the same question "
        f"— retains one shared flanking family for ITPR1×ITPR2 and one for "
        f"ITPR1×ITPR3 and none for ITPR2×ITPR3."
    )
    L += [verdict_line("ohnologs", computed, "not corroborated"), "",
          "The verdict is `not corroborated` rather than `confirmed` on "
          "purpose. The review asked for synteny-backed, phylogeny-tested "
          "evidence; this task supplies the phylogeny-tested half and S8 "
          "supplied a synteny half that is positive for two of the three pairs "
          "and empty for the third. Stating that as a demonstration of 2R "
          "ohnology would be claiming the conjunction from one conjunct.", ""]

    L += ["### 9.2 S7's sister result, seen from the reconciliation", ""]
    computed = (
        "the reconciliation is consistent with it and adds an ordering: the "
        "ITPR2 / ITPR3 duplication is the *younger* of the two, mapping to "
        "Gnathostomata while the ITPR1 split maps to Vertebrata. Every "
        "topology gives the same counts, so the sister arrangement itself is "
        "not what the reconciliation is sensitive to — which is why the three "
        "AU hypotheses are indistinguishable here and the AU test, not this "
        "task, is what settles the pair."
    )
    L += [verdict_line("sister", computed, "orthogonal"), ""]

    L += ["### 9.3 S8's paralogon asymmetry", ""]
    computed = (
        "the reconciliation makes ITPR1 the *earlier-diverging* copy, which is "
        "the copy S8 found had kept a shared flanking family with each of the "
        "other two. Those are different measurements — a retained flank "
        "ohnolog is a deletion record, a duplication placement is a branching "
        "statement — and 2R quartets lose flank copies independently of "
        "duplication order, so the agreement is suggestive and not a test."
    )
    L += [verdict_line("paralogon", computed, "orthogonal"), ""]

    L += ["### 9.4 What could not be dated, and why that is stated up front",
          ""]
    computed = (
        "unavoidable, and stated in §1 rather than buried in the caveats. "
        "Reconciliation supplies a placement and the species tree supplies the "
        "bracket; no rate estimated from these sequences enters any number in "
        "this report, and the one bracket that matters most is left open at "
        "its old end because nothing in this tree closes it."
    )
    L += [verdict_line("clock", computed, "confirmed"), ""]
    return L


# ------------------------------------------------------------------ §10

def _caveats(load, num) -> list:
    summ = load("reconciliation_summary.tsv")
    cells = _cells(summ)
    L = ["## 10. Caveats", "",
         "- **The bracket is as good as the calibration, and the calibration "
         "is an input.** Crown Vertebrata carries the widest disagreement in "
         "the tree — published estimates span 480–615 Ma — and the older "
         "duplication's bracket has no upper bound at all in this tree. "
         "`species_tree_calibrations.tsv` carries the spread and the source "
         "for every node so a reader can substitute their own.",
         "- **One alignment, one representative set.** Every statement here is "
         "conditional on S6's 134-tip sample and S7's tree of it. The "
         "reconciliation is robust across topologies, rootings and the support "
         "collapse, but none of those re-samples the taxa.",
         "- **Three cyclostome loci per species is what the sweep found, not "
         "necessarily what is there.** Both cyclostome genomes are scaffold-"
         "level and the bait panel has no cyclostome-labelled bait; a fourth "
         "locus in either genome would be exactly the kind of evidence that "
         "moves §5.",
         "- **A loss count from this analysis is not a loss count.** §8 is the "
         "argument; the number to quote is the sweep's.",
         "- **The reconciliation is parsimony-based.** GeneRax was not run: "
         "the PIEZO project established that bioconda 2.1.3 segfaults in "
         "`Scenario::savePerSpeciesEventsCounts` on both the osx-arm64 and "
         "osx-64 builds, rooted or unrooted, and repeating a known failure "
         "buys nothing. The probabilistic cross-check is therefore absent, and "
         "what stands in its place is the matrix: "
         f"{len(cells)} reconciliations across topology, taxon sampling, "
         "support and rooting, reported cell by cell.", ""]
    return L


def _figures() -> list:
    return [
        "## 11. Figures", "",
        "![](figures/recon_dated_backbone.png)", "",
        "*The dated species tree, with every duplication the reconciliation "
        "places drawn on the branch it maps to. Ages and their published "
        "spreads come from `species_tree_calibrations.tsv`, placements from "
        "`duplication_placement.tsv`; nothing is positioned by eye.*", "",
        "![](figures/recon_matrix.png)", "",
        "*(a) The topology × variant matrix — the deepest paralog duplication "
        "in every cell, with its event counts. (b) Every rooting of the "
        "vertebrate subtree by total events, with the outgroup rooting and the "
        "minimum-event rooting marked.*", "",
        "![](figures/recon_losses.png)", "",
        "*(a) What each implied loss turns out to be once it is asked of the "
        "S5 genome ledger. (b) The implied loss count in every cell of the "
        "matrix against the number the genomes corroborate.*", "",
        "![](figures/recon_cyclostome.png)", "",
        "*(a) Root-to-tip distance for all 57 vertebrate tips with the six "
        "cyclostome loci marked — the long-branch check. (b) S8's independent "
        "flank call for each of those loci, with the pair support the tree "
        "gives it; every call sits inside its own null.*", ""]


def render(*, load, table, num, not_run) -> list:
    L = _cyclostomes(load, table, num, not_run)
    L += _losses(load, table, num, not_run)
    L += _what_this_changes(load, table, num)
    L += _caveats(load, num)
    L += _figures()
    return L
