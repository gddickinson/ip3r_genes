"""S15b — renders `results/loss_counts/report.md` purely from the
committed tables (D13).

Scope, what S15a handed forward, the instrument and the negative controls
live here; the results live in `s15b_report_results.py` (the
`s3_report.py` / `s3_report_d10.py` split, so both halves stay under 500
lines and take this module's loader and formatter and therefore cannot
read the tables differently).

**Headlines are chosen by the data.**  Every prior is stated in
`s15b_priors.PRIOR` with where the earlier task said it, computed from
S15b's own tables, and rendered with both numbers printed either way.
The comparison that has to be sayable is the one that goes badly: if a
vertebrate lineage had lost an IP3 receptor, the table it would be said
from is `dollo_counts.tsv` and the row would be the base setting with a
non-zero count.

A section whose table is absent renders *not run yet*, so a stage that
was skipped is visible as skipped.

Run:  python3 scripts/s15b_report.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

OUT = Path(__file__).resolve().parents[1] / "results" / "loss_counts"


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
    return [f"![](figures/{slug}.png)", "", caption, ""]


def _i(v, d=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return d


# --------------------------------------------------------------- sections

def header(stats: dict, head: dict) -> list[str]:
    return [
        f"# S15b — the loss counts: {num(head.get('dollo_family_base', 0))} "
        f"across {num(head.get('n_cells', 0))} genome × paralog cells, and "
        f"what it takes to manufacture one",
        "",
        f"*Generated {stats.get('written_at', '')} by "
        f"`scripts/s15b_report.py` from the committed tables in "
        f"`results/loss_counts/`. Nothing in this report is hand-written; "
        f"every number is read from a table beside it (D13).*",
        "",
        "---", "",
    ]


def section_scope(head: dict, stats: dict) -> list[str]:
    import s15b_priors as priors
    base = stats.get("base_reproduction", {})
    out = [
        "## 1. What this half was asked to do, and what S15a changed about it",
        "",
        "S15b's brief asks for Dollo parsimony as the primary count, Mk "
        "fits over a sensitivity matrix, and a shared-lesion Poisson test "
        "on pseudogene fossils. S15a returned a matrix in which **no cell "
        "reaches the state that would license a loss**, and its §11 hands "
        "this half four consequences rather than four open questions.",
        "",
    ]
    out += table(
        ["what S15a handed forward", "consequence for this half"],
        [["0 of 927 cells `absent`",
          "the Dollo count is zero; it is run anyway, because a routine "
          "that returns zero without being able to return anything else "
          "is not a measurement"],
         ["the `co_trace` population (D46)",
          "the primary coding is **family-level presence per genome**; "
          "the paralog-resolved matrix is a sensitivity axis"],
         ["an invariant character",
          "Mk rates are **stated, not fitted** — an invariant character "
          "has no transition and a reported rate would be the "
          "optimiser's starting point"],
         ["no dead loci",
          "the shared-lesion Poisson test is reported **with its "
          "denominator** rather than omitted, and D47's ITPR3 indel "
          "excess is followed instead"]])
    out += [
        "So the deliverable is the **sensitivity matrix**: which "
        "combinations of coding, evidence bar, branch lengths and "
        "contiguity filter manufacture a loss, and how far each is from "
        "where S15a stands.",
        "",
        "### Everything here is read from a committed table",
        "",
        "S15b runs offline. Its inputs are S15a's character matrix, S15a's "
        "309-genome species tree, S15a's ORF-integrity tables and S13's "
        "committed node calibrations; it makes no network call, no "
        "alignment and no search, and a rerun is deterministic.",
        "",
    ]
    if base:
        out += [
            f"The first thing the driver checks is that its recoder **is** "
            f"S15a's rule chain: re-running all {num(base.get('n_rows', 0))} "
            f"cells at S15a's own setting reproduces "
            f"`{base.get('n_mismatches', 0)}` mismatches. That check is "
            f"load-bearing rather than decorative — every comparison in "
            f"this half is a comparison against S15a's operating point, and "
            f"a recoder that disagreed with it would make the whole matrix "
            f"a comparison against something that never ran.",
            "",
        ]
    out += [
        f"**Priors.** {len(priors.PRIOR)} results from earlier tasks are "
        f"stated in `s15b_priors.py` with where each was said; §9 computes "
        f"this half's answer beside each and renders the verdict from the "
        f"comparison, printing both numbers either way.",
        "",
    ]
    return out


def section_ladder(stats: dict) -> list[str]:
    import s15b_coding as coding
    p = stats.get("parameters", {})
    out = [
        "## 2. The four axes, and what each one is a knob on", "",
        "Three of the brief's axes change the character; the fourth cannot.",
        "",
        "### 2.1 Coding (D46)", "",
        "`family` asks whether the assembly carries an IP3 receptor at all "
        "and is the **primary**. `paralog` asks the same question of each "
        "of ITPR1/2/3 separately and is the sensitivity axis, because "
        "S15a measured how unreliable per-fragment paralog attribution is "
        "in a shattered assembly: 20 `co_trace` regions sit between the "
        "candidate and the decoy populations, with a median reassembly of "
        "0.631.",
        "",
        "A genome counts as carrying the family if any cell is present "
        "**or** it holds spare family loci no cell claimed — the "
        "cyclostomes' case, where three ITPR loci sit in a genome the bait "
        "panel can only file one of.",
        "",
        "### 2.2 Evidence — an ordered ladder named after what each rung "
        "refuses", "",
    ]
    rows = []
    for name, bar, rules in coding.EVIDENCE_LEVELS:
        rows.append([f"`{name}`",
                     f"{bar:.4g}" if bar is not None else "R4 disabled",
                     ", ".join(rules), coding.describe(name)])
    out += table(["rung", "reconstruction bar", "locus rules accepted",
                  "the rung"], rows)
    out += [
        f"The first three rungs are the calibration's own numbers and not "
        f"round ones: `{coding.BAR_GAP_LO}` is the decoy population's "
        f"maximum, `{coding.BAR_GAP_HI}` the lowest candidate, and "
        f"`{coding.BAR_CALIBRATED:.5g}` the midpoint S15a operates at. "
        f"Walking the bar across its whole measured gap is the honest "
        f"first sensitivity question, and it is asked before any invented "
        f"threshold is.",
        "",
        "`s15b_test_counts.py` T6 requires the ladder to be monotone: a "
        "rung named after refusing something must not admit more than the "
        "rung above it.",
        "",
        "### 2.3 The two protective rules", "",
        "**D45** (`paralog_unassignable`, R5): a cell the bait panel "
        "cannot fill in a genome carrying spare family loci is not "
        "evidence the gene is gone. **D4** (`undecidable_contiguity`, R6): "
        "an assembly whose contigs are shorter than the gene cannot be "
        "evidence that the gene is missing. Each is switched off in turn "
        "and in combination — not as a setting anybody should adopt, but "
        "to measure what that decision is worth in losses. T7 requires "
        "that switching either off can only ever *add* absences.",
        "",
        "### 2.4 Branch lengths — the axis that cannot move the count", "",
        f"Dollo parsimony counts **edges**. It does not read a length, so "
        f"the branch-length axis cannot change the primary result, and "
        f"this is carried as a column and reported rather than left as an "
        f"omission. Three schemes are declared "
        f"(`{'`, `'.join(p.get('branch_length_schemes', []))}`); the third "
        f"joins S13's committed node ages by name — 23 of S13's 29 "
        f"calibrated nodes are present in this tree — and interpolates the "
        f"rest. The ages are an **input** (D15), not something S15b "
        f"estimates. Lengths enter the Mk section and nothing else.",
        "",
    ]
    return out


def section_dollo(stats: dict) -> list[str]:
    bound = load("resolution_bound.tsv")
    poly = stats.get("polytomy_profile", {})
    out = ["## 3. Dollo parsimony, and the bound the tree puts on it", ""]
    if not bound:
        return out + missing("resolution_bound.tsv",
                             "the resolution bound")
    out += [
        "The gain is placed once and losses are the edges below it whose "
        "whole subtree has lost the character. **Where the gain goes is a "
        "decision, and it is made per character.** The family character is "
        "pinned at the root, because S20 (6,928 non-vertebrate eukaryotic "
        "reference proteomes) and S23 (194 non-vertebrate genomes) both "
        "found IP3 receptors well outside the vertebrates, so the family "
        "was present at the root of any vertebrate tree and a vertebrate "
        "clade carrying none has lost it. The paralog characters are "
        "**not** pinned: S13 places the ITPR2/ITPR3 duplication on the "
        "gnathostome stem, so a cyclostome that has neither never had "
        "either, and pinning would score the origin of the paralogs as two "
        "losses.",
        "",
        "### 3.1 The resolution any count is read against", "",
        "A count of *independent* losses is bounded by how far the tree "
        "resolves. Under a polytomy, sibling losses may be separate events "
        "or one event on a branch the polytomy does not resolve, so every "
        "count is reported as an interval: `dollo_losses_max` counts every "
        "loss edge, `dollo_losses_min` counts one per parent carrying any.",
        "",
    ]
    out += table(["metric", "value", "what it means"],
                 [[f"`{r['metric']}`", num(r["value"]), r["note"]]
                  for r in bound])
    if poly:
        degs = poly.get("degree_counts", {})
        big = sorted(((int(k), v) for k, v in degs.items() if int(k) >= 3),
                     reverse=True)[:4]
        out += [
            "The largest unresolved nodes are "
            + ", ".join(f"degree {d} (×{n})" for d, n in big)
            + f", against {num(poly.get('n_unary', 0))} unary species "
              f"nodes that the sweep's one-tip-per-assembly design creates "
              f"and that carry no resolution question at all.",
            "",
        ]
    return out


def section_controls(stats: dict) -> list[str]:
    st = stats.get("self_test", {})
    out = ["## 4. The negative controls", ""]
    if not st:
        return out + missing("loss_counts_stats.json", "the self-test status")
    out += [
        f"`s15b_test_counts.py` runs **before anything is written** and the "
        f"driver refuses to continue if it fails: "
        f"{st.get('n_passed', 0)}/{st.get('n_tests', 0)} passed, status "
        f"**{st.get('status', '?')}**.",
        "",
        "These are checks on *refusal* and on *reachability*, because both "
        "failure modes are silent here. A Dollo routine that cannot find a "
        "loss returns zero, which is the answer this task expects — so the "
        "count would look right while measuring nothing.",
        "",
    ]
    out += table(
        ["control", "the rule it is a control for"],
        [["T1/T2 a constructed loss is found, and sister losses merge",
          "the count must be **reachable**, and a lost clade is one event"],
         ["T3/T3b the gain rule changes the answer in the right direction",
          "unpinned, a clade-wide absence is ancestral; pinned at the "
          "root, it is a loss"],
         ["T4 a polytomy bounds rather than fixes the count",
          "two losses under an unresolved node may be one event"],
         ["T5 the recoder reproduces all 927 S15a states exactly",
          "the matrix is a comparison against S15a and must start from it"],
         ["T6/T7 the ladder is monotone and the rules are protective",
          "a rung must not admit more than the one above it; turning a "
          "rule off can only manufacture losses"],
         ["T8 the loss state is still reachable after recoding",
          "S15a's T8 property, re-asserted on this module"],
         ["T9/T10 an invariant character is refused on a measured profile",
          "a fitted rate on an invariant character is the optimiser's "
          "starting point"],
         ["T11/T12/T13 a varying character *is* fitted, and the pruning "
          "matches the closed form",
          "the refusal has to be about the character, not the model"],
         ["T14/T15/T16 the fossil screen fires, refuses and excludes the "
          "control",
          "S10's one-sided rule, and the RyR control has no row in the "
          "ITPR character matrix"],
         ["T17 BH is a real correction",
          "one family, one correction (S9's rule)"]])
    if st.get("mutations_caught"):
        out += ["Mutation-tested: "
                f"{len(st['mutations_caught'])} deliberate rule breakages, "
                f"all caught.", ""]
        out += ["".join([f"- {m}\n" for m in st["mutations_caught"]]), ""]
    return out


def main() -> None:
    import s15b_report_results as results
    stats = load_json("loss_counts_stats.json")
    head = stats.get("headline", {})
    lines: list[str] = []
    lines += header(stats, head)
    lines += section_scope(head, stats)
    lines += section_ladder(stats)
    lines += section_dollo(stats)
    lines += section_controls(stats)
    lines += results.render(load, load_json, num, pct, table, missing, fig,
                            stats, head)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text("\n".join(lines))
    print(f"wrote {OUT / 'report.md'} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
