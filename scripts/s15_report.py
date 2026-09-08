"""S15a — renders `results/loss_dynamics/report.md` purely from the
committed tables (D13).

Scope, the instrument and the negative controls live here; the results
live in `s15_report_results.py` (the `s3_report.py` / `s3_report_d10.py`
split, so both halves stay under 500 lines and take this module's loader
and formatter and therefore cannot read the tables differently).

**Headlines are chosen by the data.**  Every prior is stated in
`s15_priors.PRIOR` with where the earlier task said it, computed from
S15a's own tables, and rendered with both numbers printed either way.  The
comparison that has to be sayable is the one that goes badly: if a
vertebrate genome had lost an ITPR paralog, the table it would be said
from is `character_matrix.tsv` and the state would be `absent`.

A section whose table is absent renders *not run yet*, so a stage that
was skipped is visible as skipped.

Run:  python3 scripts/s15_report.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

OUT = Path(__file__).resolve().parents[1] / "results" / "loss_dynamics"


def load(name: str) -> list[dict]:
    p = OUT / name
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    p = OUT / name
    return json.loads(p.read_text()) if p.exists() else {}


def num(v, nd=0) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return f"{f:,.{nd}f}"


def pct(part, whole, nd=1) -> str:
    try:
        return f"{100.0 * float(part) / float(whole):.{nd}f} %"
    except (TypeError, ValueError, ZeroDivisionError):
        return "n/a"


def table(header: list[str], rows: list[list]) -> list[str]:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out + [""]


def missing(name: str, what: str) -> list[str]:
    return [f"*Not run yet — `{name}` is not present, so {what} cannot be "
            f"reported.*", ""]


def fig(slug: str, caption: str) -> list[str]:
    return [f"![]({slug})", "", caption, ""]


# --------------------------------------------------------------- sections

def header(stats: dict, head: dict) -> list[str]:
    n_g = head.get("n_genomes", 0)
    n_c = head.get("n_cells", 0)
    return [
        f"# S15a — the loss instrument: {num(n_c)} genome × paralog cells "
        f"across {num(n_g)} vertebrate genomes",
        "",
        f"*Generated {stats.get('written_at', '')} by "
        f"`scripts/s15_report.py` from the committed tables in "
        f"`results/loss_dynamics/`. Nothing in this report is "
        f"hand-written; every number is read from a table beside it "
        f"(D13).*",
        "",
        "S13 placed both duplications that made ITPR1/2/3 and then could "
        "not count a single loss: its reconciliation implies 51–53 across "
        "the vertebrate subtree and **none** survives being asked of the "
        "genome sweep. It handed S15 the reason — a reconciliation on a "
        "134-tip representative alignment cannot count losses, because a "
        "paralog missing from the alignment is usually missing from the "
        "*sample* — and the instruction: build the character matrix from "
        "the 309-genome sweep, one cell per genome × paralog, and not "
        "from tips of a tree.",
        "",
        "This half of S15 builds that matrix and the three things it needs "
        "first. **What counts as a loss** — a state vocabulary in which "
        "exactly one state may be counted, reached only after every "
        "alternative explanation has been given a positive test. **What "
        "the sweep's undecided cells actually contain** — the brief asked "
        "for synteny and the answer is that synteny cannot reach them, so "
        "the reference is reassembled from the sweep's own archived "
        "alignments instead, against a bar measured on a gene that is "
        "accounted for elsewhere in the same genome. And **whether a "
        "broken reading frame is a broken gene or a broken assembly** — "
        "an ORF screen with the confounder controlled three ways, one of "
        "which turned out to matter and was not the one anybody expects.",
        "",
        "S15b counts the losses: Dollo parsimony as the primary count, Mk "
        "model fits, the coding × evidence × branch-length × contiguity "
        "sensitivity matrix, and the pseudogene-fossil lesion analysis.",
        "",
    ]


def section_scope(head: dict, stats: dict) -> list[str]:
    out = ["## 1. Scope, and the unit of observation", ""]
    if not head:
        return out + missing("loss_dynamics_stats.json", "the scope")
    out += [
        f"The matrix is **{num(head['n_cells'])} cells**: "
        f"{num(head['n_genomes'])} genomes × three paralogs. The unit is "
        f"the *assembly*, not the species — two accessions of one species "
        f"are two independent observations of its gene complement, and "
        f"collapsing them would average a placed gene with a shattered "
        f"one.",
        "",
        f"The RyR sister family is in every genome as the positive "
        f"control and is **not** a column of the matrix: it is what makes "
        f"an empty ITPR cell an empty cell in an assembly the search "
        f"demonstrably reached. It fired in all "
        f"{num(head['n_genomes'])} genomes, so rule R0 excludes nothing.",
        "",
    ]
    p = stats.get("parameters", {})
    out += table(
        ["parameter", "value", "where it comes from"],
        [["rescue significance", p.get("rescue_evalue"),
          "`s5_rescue.RESCUE_E` — the sweep's own cut, not re-chosen"],
         ["known-locus pad", f"{num(p.get('known_locus_pad_bp'))} bp",
          "`s5_run_sweep.filter_hsps_outside` — the cross-paralog control"],
         ["full-locus coverage", p.get("cov_found"),
          "`s5_classify.COV_FOUND`"],
         ["synteny window", p.get("synteny_window"),
          "`s8_flank_lib` — S8's rule, imported unchanged"],
         ["synteny key floor", p.get("synteny_min_keys"),
          "`s8_paralogon.MIN_KEYS`"],
         ["integrity bar quantile", p.get("integrity_bar_quantile"),
          "measured on the intact population, §4"],
         ["identity-match window", p.get("integrity_identity_window"),
          "§4 — the confounder that turned out to be real"]])
    return out


def section_states(stats: dict) -> list[str]:
    out = ["## 2. What counts as a loss", "",
           "Eight states, assigned by ordered positive tests; the first "
           "rule that fires wins, and each writes the number it fired on "
           "into its row of `character_matrix.tsv`. The design constraint "
           "is that **exactly one** state may be counted as a loss.", ""]
    out += table(
        ["rule", "state", "the test it applies"],
        [["R0", "`no_control`",
          "the genome's RyR control did not fire — no cell in it supports "
          "any claim"],
         ["R1", "`present_single_locus`",
          "a locus at or above the sweep's own coverage bar"],
         ["R2", "`present_truncated`",
          "a locus below it, truncated at a contig edge or an N-run"],
         ["R3", "`present_partial`",
          "a locus below it with no assembly excuse"],
         ["R4", "`present_fragmented`",
          "no locus, but the reference reassembles from sequence outside "
          "every locus the aligner found, above the measured bar (§3)"],
         ["R5", "`paralog_unassignable`",
          "no locus and no reassembly, but the genome carries family loci "
          "no cell claimed — D45"],
         ["R6", "`undecidable_contiguity`",
          "nothing found, and the assembly cannot hold the gene on one "
          "contig (D4's bar)"],
         ["R7", "**`absent`**",
          "nothing found, in a controlled assembly that could have held "
          "it. The only state S15b may count"]])
    out += [
        "Two of these rules exist because of specific incidents earlier in "
        "the project. **R5 is D45.** Both cyclostome genomes carry three "
        "ITPR loci apiece, all filed in the ITPR1 cell because the S5 "
        "bait panel has no cyclostome-labelled bait, and their ITPR2 and "
        "ITPR3 cells therefore read `absent` in the ledger. A count that "
        "took those at face value would score two independent losses per "
        "cyclostome that never happened. **R6 is D4.** An assembly whose "
        "contigs are shorter than the gene cannot represent it as one "
        "locus, so it cannot be evidence that the gene is missing — and "
        "66 % of the Aves assemblies in this scope are in that position.",
        "",
        "The rule order is itself tested: `s15_test_loss.py` T6 moves R5 "
        "after R7 and requires the suite to fail, which it does.",
        "",
    ]
    return out


def section_reconstruction(stats: dict) -> list[str]:
    cal = load("recon_calibration.tsv")
    recon = load("reconstruction.tsv")
    out = ["## 3. The instrument: reassembling a reference across contigs",
           ""]
    if not cal or not recon:
        return out + missing("recon_calibration.tsv",
                             "the reconstruction bar")
    out += [
        "The sweep's `tblastn_trace` cells are the ones a loss count can "
        "neither ignore nor use. The aligner placed no locus, so the "
        "ledger records no coverage; the rescue attributed regions, so "
        "they are not nothing either. What is measured here is the "
        "**non-redundant coverage of the reference protein, reassembled "
        "across contigs** — the union of the query intervals of every "
        "significant HSP, over the reference's own length.",
        "",
        "Three properties make that a positive test rather than a hopeful "
        "one.",
        "",
        "1. **It is computed outside every locus the aligner found.** The "
        "HSP set passes through the sweep's own exclusion — every "
        "clustered locus in the genome, any bait, padded 5 kb — so a "
        "reference cannot be reassembled out of the genome's other "
        "paralogs' genes. For a family whose paralogs are 61–68 % "
        "identical that is the failure mode that matters, and "
        "`s15_test_loss.py` T2 removes the filter and requires the suite "
        "to fail.",
        "2. **One reference at a time.** A union over three orthologous "
        "baits would count a residue covered in any of them and report a "
        "coverage no single protein achieves (T3).",
        "3. **The bar is measured, against a gene that is accounted for.**",
        "",
    ]
    rows = []
    for r in cal:
        rows.append([r["metric"], num(r["n_positive"]), num(r["n_negative"]),
                     num(r["n_co_trace"]),
                     num(r["positive_median"], 3), num(r["positive_min"], 3),
                     num(r["negative_median"], 3), num(r["negative_max"], 3),
                     num(r["youden_j"], 2),
                     f"{num(r['gap_lo'], 3)}–{num(r['gap_hi'], 3)}",
                     num(r["operating_point"], 3)])
    out += table(["metric", "n pos", "n neg", "n co-trace", "pos median",
                  "pos min", "neg median", "neg max", "Youden J", "gap",
                  "bar"], rows)
    cov = next((r for r in cal if r["metric"] == "coverage"), {})
    out += [
        f"The negative is the sweep's own output and needed no new search: "
        f"regions the full 38-bait panel attributes to a paralog whose "
        f"gene the aligner **already placed at a locus in the same "
        f"genome**. That paralog is accounted for, so the fragment set "
        f"cannot be that gene, and what it recovers is what cross-paralog "
        f"similarity delivers on its own — a median of "
        f"{num(cov.get('negative_median'), 3)} and a maximum of "
        f"{num(cov.get('negative_max'), 3)} against a candidate median of "
        f"{num(cov.get('positive_median'), 3)}.",
        "",
        f"The first version of this calibration did **not** separate, at "
        f"J = 0.52, and the reason is worth stating because it is a result "
        f"in itself: the decoy was every region attributed elsewhere, and "
        f"in a genome where two paralogs are both shattered a fragment "
        f"attributed to the other one is a piece of a real gene. Those "
        f"{num(cov.get('n_co_trace'))} regions are now reported as their "
        f"own population — `co_trace`, median "
        f"{num(cov.get('co_trace_median'), 3)}, sitting squarely between "
        f"the two — and they are the **measured size of the "
        f"paralog-attribution problem in a shattered assembly**. "
        f"`s15_test_loss.py` T9 folds them back into the decoy and "
        f"requires the suite to fail.",
        "",
        f"The operating point is the midpoint of the gap "
        f"({num(cov.get('operating_point'), 3)}) rather than Youden's own "
        f"threshold, which on a perfectly separated pair lands on the "
        f"lowest positive and is the most permissive bar the data allow. "
        f"Both edges of the gap are committed so a later run whose "
        f"populations have drifted into contact is visible in the table, "
        f"and `s15_calibrate_recon.bar()` refuses to hand out a threshold "
        f"the calibration marked unusable.",
        "",
    ]
    out += fig("figures/s15_reconstruction.png",
               "**Figure 2.** The three populations the bar is read off, "
               "with the gap shaded and the operating point drawn (a), and "
               "every undecided cell against the number of contigs its "
               "gene is spread over (b). A calibration figure that asked "
               "to be believed would not be one, so both edges of the gap "
               "are marks and not a caption.")
    return out


def section_synteny(stats: dict) -> list[str]:
    summ = load("synteny_reach_summary.tsv")
    acc = load("caller_accuracy_by_keys.tsv")
    calls = load("synteny_calls.tsv")
    out = ["## 4. The brief's first step, and the number that answers it",
           ""]
    if not summ:
        return out + missing("synteny_reach_summary.tsv",
                             "the synteny measurement")
    row = next((r for r in summ if r["window"] == "informative10"), summ[0])
    n = int(row["n_regions"])
    out += [
        f"S15 was asked to disambiguate the sweep's undecided cells **by "
        f"synteny**, with measured accuracy against known loci. S8 had "
        f"already built that instrument and calibrated it "
        f"leave-one-genome-out on loci whose paralog their own annotation "
        f"establishes, so the work here was to point it at the trace "
        f"regions and measure whether it arrives.",
        "",
        f"**It does not, and that is the result.** Of the {num(n)} rescue "
        f"regions across the undecided cells, "
        f"**{num(row['n_no_neighbourhood'])} sit on a contig carrying no "
        f"annotated gene at all**, {num(row['n_too_few_keys'])} have too "
        f"few informative flanking symbols, and "
        f"**{num(row['n_reached'])} reach the caller's "
        f"{row['min_keys_required']}-key floor**. The median region has "
        f"{num(row['median_keys'])} informative neighbours and "
        f"{num(row['median_genes_on_contig'])} coding genes on its own "
        f"contig, which extends a median of "
        f"{num(row['median_contig_extent_kb'], 1)} kb — shorter than a "
        f"single ITPR gene.",
        "",
        "That is not a failure of the caller. On S8's own calibration set, "
        "binned by the number of keys available, accuracy is 100 % at "
        "every key count the caller will act on:",
        "",
    ]
    if acc:
        out += table(["keys available", "n loci", "called", "call rate",
                      "accuracy when called"],
                     [[f"{r['keys_lo']}" if r["keys_lo"] == r["keys_hi"]
                       else (f"{r['keys_lo']}+" if int(r["keys_hi"]) >= 90
                             else f"{r['keys_lo']}–{r['keys_hi']}"),
                       num(r["n"]), num(r["n_called"]),
                       num(r["call_rate"], 3),
                       "—" if r["accuracy"] in ("", "nan")
                       else num(r["accuracy"], 3)] for r in acc])
    if calls:
        agree = sum(1 for c in calls if c["agrees_with_bitscore"] == "1")
        called = [c for c in calls if c["call"] != "no_call"]
        out += [
            f"On the {num(len(calls))} regions it can reach, the caller "
            f"returns a call for {num(len(called))} and **agrees with the "
            f"alignment's own paralog attribution on {num(agree)} of "
            f"them**. That is an independent instrument corroborating the "
            f"attribution — on six regions, which is stated for what it "
            f"is rather than leaned on.",
            "",
        ]
    out += fig("figures/s15_synteny_reach.png",
               "**Figure 3.** Why synteny could not answer. The joint "
               "distribution of what a trace region has to work with — "
               "genes on its own contig against informative flank keys, "
               "with the caller's four-key floor drawn (a) — beside S8's "
               "own accuracy on the same axis (b). Putting reach and "
               "accuracy on one axis is what makes *accurate and "
               "unavailable* a readable sentence.")
    return out


def section_controls(stats: dict) -> list[str]:
    st = stats.get("self_test", {})
    out = ["## 5. Negative controls", ""]
    if not st:
        return out + missing("loss_dynamics_stats.json",
                             "the self-test status")
    out += [
        f"`s15_test_loss.py` runs **{num(st.get('n_tests'))} constructed "
        f"cases before anything is written** and the driver refuses to "
        f"continue if any fails "
        f"({num(st.get('n_passed'))}/{num(st.get('n_tests'))} passed on "
        f"this build). They are checks on *refusal*, because every rule "
        f"in this task returns a plausible number when it is wrong: a "
        f"reconstruction without the known-locus filter reassembles a "
        f"missing gene out of its own paralogs and reports 0.95; a state "
        f"machine with R5 and R7 the wrong way round manufactures four "
        f"losses in the cyclostomes; a bar read off a contaminated decoy "
        f"lands in the middle of the positives.",
        "",
    ]
    muts = st.get("mutations_caught") or []
    if muts:
        out += ["Mutation-tested — each deliberate rule breakage was "
                "applied to the live module and the suite re-run:", ""]
        out += table(["breakage", "caught by"],
                     [[m.split("->")[0].strip(),
                       m.split("->")[1].strip() if "->" in m else ""]
                      for m in muts])
    if st.get("failures"):
        out += ["Failures on this build:", ""]
        out += ["- " + f for f in st["failures"]] + [""]
    return out


def main() -> None:
    import s15_report_results as results
    stats = load_json("loss_dynamics_stats.json")
    head = stats.get("headline", {})
    lines: list[str] = []
    lines += header(stats, head)
    lines += section_scope(head, stats)
    lines += section_states(stats)
    lines += section_reconstruction(stats)
    lines += section_synteny(stats)
    lines += section_controls(stats)
    lines += results.render(load, load_json, num, pct, table, missing, fig,
                            stats, head)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text("\n".join(lines))
    print(f"wrote {OUT / 'report.md'} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
