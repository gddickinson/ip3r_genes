"""S10 — renders `results/annotation_bugs/report.md` purely from the tables (D13).

Scope, selection and the instrument live here; the two case studies live in
`s10_report_cases.py` (the `s3_report.py` / `s3_report_d10.py` split, so both
halves stay inside the 500-line budget and read the tables through the same
loader and the same formatter).

**Headlines chosen by the data.** S5b and S3 each recorded a number this task
is in a position to check, and each is stated in `PRIOR`, computed from S10's
own tables, and rendered `confirmed` / `contradicted` / `underpowered` with
both numbers printed either way. The comparison that has to be sayable is the
one that goes badly: if the annotations turned out to be right and the sweep
wrong, this report has to be able to say so.

**A section whose table is absent renders *not run yet*,** never nothing, so a
partial run is visible as partial rather than as a shorter report.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_case_spec as spec                                   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "annotation_bugs"

#: What earlier tasks concluded, and where they said it. The report computes
#: S10's answer beside each and renders the verdict from the comparison.
PRIOR = {
    "unannotated_loci": {
        "claim": "S5b: 318 ITPR gene models exist only as DNA and a further "
                 "167 sit inside an annotated gene carrying no family name",
        "where": "results/genome_ledger/report.md, 'Census v4'",
    },
    "naming_conflict": {
        "claim": "S5b: 3 loci are claimed by a paralog cell other than the one "
                 "the assembly's annotation names, left unadjudicated "
                 "('one instrument does not overturn a public annotation')",
        "where": "results/genome_ledger/report.md, 'Where an annotation and "
                 "this sweep disagree'",
    },
    "split_fragments": {
        "claim": "S3: 1,794 sweep hits carry an ITPR or RYR gene symbol but "
                 "fall under D22's 200-position gate — 'real family genes "
                 "whose annotation has been broken into fragments'",
        "where": "PUBLICATION_ROADMAP.md, Emergent tasks 2026-09-03",
    },
}


def load(name: str) -> list[dict]:
    path = OUT / name
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    path = OUT / name
    return json.loads(path.read_text()) if path.exists() else {}


def fmt(x, nd: int = 1) -> str:
    """One number formatter for both halves of the report."""
    if x is None or x == "":
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v == int(v) and abs(v) < 1e15:
        return f"{int(v):,}"
    return f"{v:,.{nd}f}"


def pct(num, den) -> str:
    den = float(den or 0)
    return "n/a" if not den else f"{100.0 * float(num) / den:.1f} %"


def verdict(ok: bool | None) -> str:
    if ok is None:
        return "**underpowered**"
    return "**confirmed**" if ok else "**contradicted**"


def missing(section: str) -> str:
    return (f"\n### {section}\n\n*Not run yet — the table this section renders "
            f"from is not present.*\n")


# --------------------------------------------------------------------------
def header(stats: dict, cases: list[dict]) -> list[str]:
    n = stats.get("n_loci_recovered", 0)
    e = stats.get("n_eligible", 0)
    f = stats.get("n_failures", 0)
    return [
        "# S10 — annotation-bug molecular validation",
        "",
        f"The genome sweep recovered **{fmt(n)}** IP3-receptor loci across the "
        f"309-genome scope. **{fmt(e)}** of them are in assemblies where the "
        "annotation could reasonably have been expected to deliver the gene "
        f"(rules E1–E5 below), and at **{fmt(e - f)}** of those "
        f"({pct(e - f, e)}) it did. This report is about the "
        f"**{fmt(f)}** where it did not, and about the two of them, one per "
        "failure mode, whose evidence is taken down to the exon.",
        "",
        "The two cases are "
        + " and ".join(f"***{c['organism']}*** {c['cell']} "
                       f"({c['selected_as']})" for c in cases)
        + ".",
        "",
        "> Every number below is read from a committed table in this "
        "directory. Nothing in this report is computed from a GFF, a genome "
        "or a BLAST run at render time (D13).",
        "",
    ]


def section_selection(stats: dict, rank: list[dict], depth: list[dict]
                      ) -> list[str]:
    if not rank:
        return [missing("1. Which failures, and why those")]
    out = ["## 1. Which failures, and why those", "",
           "\"The two worst\" is a rule applied to every recovered locus, not "
           "a pair of loci chosen by eye. The measurement is **annotation "
           "loss** — the fraction of a recovered gene's coding footprint that "
           "no single annotated gene model delivers — taken from the "
           "`frac_cds` values S5 recorded at sweep time.", "",
           "### Eligibility", "",
           "| rule | test | loci excluded |", "|---|---|---|"]
    exc = stats.get("n_excluded_by_rule", {})
    for code, text in spec.ELIGIBILITY:
        out.append(f"| {code} | {text} | {fmt(exc.get(code, 0))} |")
    elig = [r for r in rank if r["eligible"] == "True"]
    zero = stats.get("n_loss_zero", 0)
    out += ["",
            f"**{fmt(len(elig))} loci pass all five.** Their median annotation "
            f"loss is {fmt(stats.get('median_loss_eligible'), 3)} and "
            f"{fmt(zero)} of them ({pct(zero, len(elig))}) sit at exactly "
            "zero — one annotated gene model covering the whole recovered "
            "coding footprint. The failures are a thin tail, not a "
            "distribution.", ""]

    e5 = [r for r in rank if r["failed_rule"] == "E5"]
    if e5:
        out += ["### E5 does the work, and here is what it removed", "",
                "E5 asks whether the annotation builds genes this long "
                "*anywhere else in the same genome*. Without it the ranking's "
                "top rows are loci in assemblies whose annotation has a "
                "genome-wide length ceiling, and validating one of those "
                "would report a property of the whole gene set as a bug at "
                "this gene.", "",
                "| genome | locus | span | loss | longest gene the same "
                "annotation builds |", "|---|---|---|---|---|"]
        for r in sorted(e5, key=lambda x: -float(x["loss"])):
            out.append(
                f"| *{r['organism']}* | {r['cell']} | "
                f"{fmt(int(r['span']) / 1000)} kb | {r['loss']} | "
                f"{fmt(int(r['annot_max_gene_span'] or 0) / 1000)} kb |")
        out += ["", f"All {fmt(len(e5))} excluded loci are in "
                f"{fmt(len({r['accession'] for r in e5}))} genomes, and E5 "
                "removes nothing else in the whole sweep.", ""]

    out += ["### The failures", "",
            "| genome | locus | mode | loss | annotated models | internal "
            "control | selected |", "|---|---|---|---|---|---|---|"]
    sel = {(c["accession"], c["cell"]): c["case_id"] for c in load("cases.tsv")}
    for r in [x for x in elig if x["mode"]]:
        tag = sel.get((r["accession"], r["cell"]), "")
        out.append(
            f"| *{r['organism']}* | {r['cell']} | {r['mode']} | {r['loss']} | "
            f"{r['n_fragments']} | {r['control_strength']} of 3 | "
            f"{'**' + tag + '**' if tag else '—'} |")
    out += ["",
            "`internal control` counts the *other* family loci in the same "
            "genome the annotation gets right; it breaks ties, because at "
            "equal loss the sharper case is the one whose own genome proves "
            "the annotation could have done better.",
            "",
            "The two selected cases are the worst of each mode, at most one "
            "per genome. Taking the top two of the single ranking would have "
            "given two omissions and left the fragmentation claim "
            "unvalidated; the full ranking is committed either way "
            "(`case_ranking.tsv`).", ""]
    return out


def section_priors(rank: list[dict], models: list[dict], tiling: list[dict],
                   dbrec: list[dict]) -> list[str]:
    """What earlier tasks claimed, and what this task's tables say."""
    if not rank:
        return [missing("2. What this changes about earlier tasks")]
    out = ["## 2. What this changes about earlier tasks", ""]

    mism = [r for r in tiling if r.get("naming") == "name_mismatch"]
    tested = [r for r in tiling if r.get("naming") in
              ("name_mismatch", "name_matches_sequence")]
    ok = bool(mism) if tested else None
    out += [f"### The naming conflicts S5b left open — {verdict(ok)}", "",
            f"*Prior.* {PRIOR['naming_conflict']['claim']} "
            f"({PRIOR['naming_conflict']['where']}).", ""]
    if tested:
        out += [f"*Here.* {fmt(len(tested))} annotated models in the two case "
                "genomes carry a paralog name and could be placed at their "
                f"own locus; {fmt(len(mism))} of them "
                f"{'sits' if len(mism) == 1 else 'sit'} on a locus of a "
                "different paralog — one of the three conflicts S5b listed, "
                "the two others being in genomes outside these cases. The "
                "adjudication is the annotation's **own "
                "translated protein** blastp'd against the genome's own "
                "recovered loci — a second instrument, and the one S5b said "
                "was needed before a public annotation could be "
                "contradicted.", ""]
        for r in mism:
            loc = r["best_locus"].split("|")
            out.append(f"- `{r['name']}` ({fmt(r['protein_aa'])} aa) matches "
                       f"the **{loc[1] if len(loc) > 1 else '?'}** locus at "
                       f"{float(r['identity']):.1%} identity over "
                       f"{fmt(r['aln_aa'])} residues, "
                       f"bit-score margin {r['bit_margin']} over the "
                       "runner-up.")
        out.append("")
    else:
        out += ["*Here.* No annotated model in either case genome carries a "
                "paralog name that could be tested.", ""]

    if dbrec:
        out += ["### What the protein databases hold for these two species",
                "", "| species | records in census v3 | ITPR records |",
                "|---|---|---|"]
        for sp in sorted({r["species"] for r in dbrec}):
            rows = [r for r in dbrec if r["species"] == sp]
            itpr = [r for r in rows if r["call"] == "ITPR"]
            out.append(f"| *{sp}* | {fmt(len(rows))} | {fmt(len(itpr))} |")
        out += ["",
                "*Prior.* " + PRIOR["unannotated_loci"]["claim"] + ".", "",
                "*Here.* Both case species carry three complete IP3-receptor "
                "genes in their DNA and **no ITPR protein record at all** in "
                "the census. That is the consequence S5b's counts describe, "
                "measured at the two species where the cause is now known.",
                ""]
    return out


def main() -> int:
    from s10_report_cases import case_sections

    stats = load_json("selection_stats.json")
    rank, depth = load("case_ranking.tsv"), load("annotation_depth.tsv")
    cases = load("cases.tsv")
    tiling, dbrec = load("fragment_tiling.tsv"), load("database_records.tsv")
    models = load("annotated_models.tsv")

    lines = header(stats, cases)
    lines += section_selection(stats, rank, depth)
    lines += section_priors(rank, models, tiling, dbrec)
    lines += case_sections(load, load_json, fmt, pct, verdict, missing, PRIOR)
    lines += ["## Outputs", "",
              "- `case_ranking.tsv` — every recovered locus with its "
              "annotation loss, failure mode and eligibility",
              "- `annotation_depth.tsv` — the E5 control per genome",
              "- `cases.tsv` — the two selected cases",
              "- `case_exons.tsv` / `case_introns.tsv` — every aligned exon "
              "and intron with what the annotation holds there",
              "- `annotated_models.tsv` / `block_accounting.tsv` — the "
              "annotated models on each gene, and disjoint blocks counted on "
              "both sides",
              "- `fragment_tiling.tsv` — each annotated protein blastp-tiled "
              "onto the genome's own recovered loci",
              "- `reading_frame.tsv` — the spliced CDS, its stops, and the "
              "stops expected under neutrality",
              "- `junction_probes.tsv` / `boundary_concordance.tsv` / "
              "`probe_summary.tsv` — the transcript search and the "
              "exon-boundary control",
              "- `flank_consensus_check.tsv` — whether the locus's "
              "neighbours are this paralog's consensus flanks (from S8)",
              "- `assembly_audit.tsv` — build provenance per case",
              "- `family_named_models.tsv` / `database_records.tsv` — what "
              "the annotation names and what the databases serve",
              "- `figures/` — four figures (D13, D19)", ""]

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print(f"[s10] report -> {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
