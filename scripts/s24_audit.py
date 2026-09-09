"""The figure audit (D11): open every figure and read it against its legend.

Two halves, and only one of them can be code.

**The mechanical half is checked here.** Every figure the manifest lists must
have a legend that names it; every legend must name a figure that exists; and
for each Extended Data figure the panel letters its legend uses must match the
panel files the manifest holds — the check that catches "the legend says four,
the figure draws three" for a multi-file figure. These run on every build and
a failure is reported, not warned about.

**The other half is an inspection, and it is recorded rather than computed.**
Whether a legend's *description* matches what the panel draws cannot be
derived from the files: someone has to look. `FINDINGS` is that record — one
row per thing the S24 inspection found wrong, naming the figure, what the
legend said, what the figure shows, how it was verified, and what was changed.
Writing it as data rather than prose means the audit is a table in `results/`
with the same standing as any other, and the next session can see what was
looked at as well as what was found.

Every numeric finding here was re-derived from the committed table named in
its `verified_against` field before it was written down.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s24_lib as L                                            # noqa: E402

MANUSCRIPT = L.PROJECT_ROOT / "manuscript"
MANIFEST = MANUSCRIPT / "figure_manifest.tsv"
LEGENDS_MAIN = MANUSCRIPT / "11_figure_legends.md"
LEGENDS_ED = MANUSCRIPT / "12_extended_data.md"
MANUSCRIPT_MD = MANUSCRIPT / "manuscript.md"

#: What looking at the figures found. `status` is `fixed_figure` when the
#: figure was redrawn, `fixed_legend` when the legend was the thing that was
#: wrong, and `noted` when neither needed changing but the reader does.
FINDINGS: list[dict] = [
    dict(figure="Fig. 2", panel="key",
         legend_said="the two lighter blues are a locus found in an assembly "
                     "whose annotation misses it or that has no gene set; "
                     "grey is a locus the assembly is too fragmented to place",
         figure_shows="four blues, not two — found+annotated, found "
                      "unannotated, found with no gene set, and partial "
                      "locus — and the two greys are assembly gap and "
                      "ambiguous trace; the fragmentary class is the palest "
                      "blue, not the grey",
         verified_against="scripts/figstyle.py:STATUS_ORDER, "
                          "results/genome_ledger/ledger_status_counts.tsv",
         status="fixed_legend"),
    dict(figure="Fig. 2", panel="whole",
         legend_said="evidence class ... as a fraction of the genomes of "
                     "each vertebrate class",
         figure_shows="six of the scope's thirteen classes, 302 of 309 "
                      "genomes; the figure keeps only classes with at least "
                      "two genomes",
         verified_against="results/genome_ledger/genome_ledger.tsv, "
                          "scripts/s5_figures.py:fig_ledger_by_class",
         status="fixed_legend"),
    dict(figure="Fig. 2", panel="key",
         legend_said="the absent colour appears nowhere at this scale: four "
                     "cells (both cyclostomes' ITPR2 and ITPR3) are treated "
                     "as paralogue-unassignable rather than absent",
         figure_shows="the ledger does hold those four cells as absent; they "
                      "are not drawn because Myxini and Hyperoartia have one "
                      "genome each and the panel drops classes below two. "
                      "Paralogue-unassignable is S15's re-statement of them, "
                      "not this figure's",
         verified_against="results/genome_ledger/genome_ledger.tsv "
                          "(4 rows with status=absent, both cyclostomes)",
         status="fixed_legend"),
    dict(figure="Fig. 2", panel="whole",
         legend_said="Aves carry the most non-blue area",
         figure_shows="Lepidosauria does, on the fraction the panel plots "
                      "(0.27 mean against Aves' 0.23); Aves carries the most "
                      "non-blue cells because it is the largest class",
         verified_against="results/genome_ledger/genome_ledger.tsv, "
                          "results/genome_ledger/contiguity_bar.tsv",
         status="fixed_legend"),
    dict(figure="Fig. 3", panel="backbone",
         legend_said="the four 100/100 labels on the backbone are the "
                     "extended paralogue clades",
         figure_shows="five labels: the two ancestors of each boxed clade, "
                      "which are the three extended paralogue clades, an "
                      "inner ITPR1 node of 15 tips, and the 32-tip node "
                      "uniting ITPR2 and ITPR3",
         verified_against="results/phylogeny/rooted.nwk re-walked with "
                          "scripts/s7_figure.py's own annotate rule",
         status="fixed_legend"),
    dict(figure="Fig. 6", panel="a",
         legend_said="(no panel is referenced)",
         figure_shows="a bold panel letter 'a' with no panel b anywhere in "
                      "the figure",
         verified_against="scripts/s17_figures.py:fig_channel_profile",
         status="fixed_figure"),
    dict(figure="Fig. 1", panel="grouping",
         legend_said="lineages grouped by kingdom",
         figure_shows="the sweep's taxonomic groups, three of which — SAR, "
                      "Discoba, Amoebozoa — are not kingdoms",
         verified_against="results/s20_sweep/proteome_presence.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 1", panel="a",
         legend_said="the long right tail is three genomes",
         figure_shows="nine genomes carry more than three complete gene "
                      "models, and four carry six or more",
         verified_against="results/s23_scope/copy_number_ledger.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 1", panel="b",
         legend_said="the grey bar reaching the light bar in every row is "
                     "the result",
         figure_shows="the light bar is hidden behind the grey one in every "
                      "row, so the coincidence the panel exists to show "
                      "cannot be seen",
         verified_against="results/s23_scope/absence_at_genome.tsv",
         status="fixed_figure"),
    dict(figure="Extended Data Fig. 1", panel="c",
         legend_said="every record sits between 20 % and 46 %",
         figure_shows="19.9 % to 45.8 %, and a second dashed line at 80 % "
                      "that the legend does not explain",
         verified_against="results/s20_sweep/plant_fungal_verdicts.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 2", panel="c",
         legend_said="against three to six copies for the control",
         figure_shows="the ryanodine control's most common count is two "
                      "(33 % of genomes), then three (30 %) and six (21 %)",
         verified_against="results/genome_ledger/genome_ledger.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 3", panel="a",
         legend_said="per-column conservation",
         figure_shows="a rolling mean of 25 columns, and an unexplained "
                      "dashed mean line",
         verified_against="scripts/s6_figures.py",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 7", panel="a",
         legend_said="each calibrated node's published age spread drawn as "
                     "a band",
         figure_shows="a band on the two nodes a duplication is placed at, "
                      "and on no other",
         verified_against="scripts/s13_figures.py "
                          "(bands are drawn for `used` nodes only)",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 7", panel="c",
         legend_said="each of the 51–53 implied losses",
         figure_shows="51 to 102: the two support-collapsed variants imply "
                      "83 and 102",
         verified_against="results/reconciliation/reconciliation_summary.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 8", panel="a",
         legend_said="(the second panel is not described)",
         figure_shows="two panels: the calibration, and how far each "
                      "undecided cell's gene is spread across contigs",
         verified_against="results/loss_dynamics/figures/"
                          "s15_reconstruction.png",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 8", panel="c",
         legend_said="(not a legend fault)",
         figure_shows="the two 'below the caller's floor' annotations "
                      "overlap into unreadable text",
         verified_against="results/loss_dynamics/figures/"
                          "s15_synteny_reach.png",
         status="fixed_figure"),
    dict(figure="Extended Data Fig. 8", panel="d",
         legend_said="(not a legend fault)",
         figure_shows="the two heat maps carry independent colour scales, so "
                      "the single manufactured loss in panel a is drawn as "
                      "dark as the 45 in panel b",
         verified_against="results/loss_counts/sensitivity_matrix.tsv",
         status="fixed_figure"),
    dict(figure="Extended Data Fig. 10", panel="a",
         legend_said="the 29-structure panel",
         figure_shows="30 bars: 29 usable structures plus the ITPR2 record "
                      "AlphaFold DB serves as a 181-residue isoform",
         verified_against="results/structures/structure_manifest.tsv "
                          "(30 rows, 29 status=ok, 1 too_short)",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 10", panel="b",
         legend_said="(marker shape is not mentioned)",
         figure_shows="triangles, circles and squares carry meaning that "
                      "appears in no key",
         verified_against="scripts/s11_figures.py",
         status="fixed_figure"),
    dict(figure="Extended Data Fig. 11", panel="a",
         legend_said="each protein's own linker mean drawn",
         figure_shows="one dashed line per panel — the mean over all three "
                      "proteins' linker controls, not three lines",
         verified_against="scripts/s17_figures.py:fig_elements",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 11", panel="c",
         legend_said="(not a legend fault)",
         figure_shows="the four AUCs are printed at mixed precision "
                      "(0.758 beside 0.6842)",
         verified_against="results/constraint/variant_constraint_test.tsv",
         status="fixed_figure"),
    dict(figure="Extended Data Fig. 14", panel="c",
         legend_said="almost entirely below 1,000 residues, except in the "
                     "protists",
         figure_shows="the fungal 1,000–1,999 aa bar is 0.58 and the plant "
                      "≥ 2,000 aa bar 0.16; both rest on a dozen to two "
                      "dozen records",
         verified_against="results/methods/head_to_head.tsv",
         status="fixed_legend"),
    dict(figure="Extended Data Fig. 7", panel="d",
         legend_said="(not a legend fault)",
         figure_shows="the panel title read 'every one inside its null' when "
                      "two of the six loci made no call at all, and the "
                      "right-hand support numbers were unlabelled",
         verified_against="results/reconciliation/cyclostome_loci.tsv "
                          "(4 within_null, 2 no_call)",
         status="fixed_figure"),
    dict(figure="Supplementary Fig. 1", panel="c",
         legend_said="(not a legend fault)",
         figure_shows="the palest grey in the quality scale vanished at "
                      "printed size in the assembled PDF, so one of the two "
                      "histograms was invisible on the page",
         verified_against="manuscript/itpr_family_manuscript.pdf p. 47",
         status="fixed_figure"),
    dict(figure="every figure in the project", panel="pdf output",
         legend_said="(not a legend fault)",
         figure_shows="matplotlib stamps the wall clock into a PDF's "
                      "/CreationDate, so a figure rebuilt from the same code "
                      "on the same data differed from itself by two bytes "
                      "and no SHA-256 recorded against a pdf meant anything",
         verified_against="two runs of scripts/s24_run.py --only figures; "
                          "12 of 12 files now byte-identical",
         status="fixed_figure"),
    dict(figure="Fig. 7, Extended Data Figs 9, 12, 13", panel="letters",
         legend_said="(**A**) / (**B**)",
         figure_shows="uppercase panel letters in four figures where every "
                      "other figure in the set uses lowercase",
         verified_against="scripts/s9_figures.py, scripts/s10_figures.py, "
                          "scripts/s18_figures.py",
         status="fixed_figure"),
]


def _panel_letters(text: str) -> list[str]:
    """The panel letters a legend paragraph references, in order."""
    return re.findall(r"\(\*\*([a-zA-Z])\*\*\)", text)


def split_legends(path: Path, pattern: str) -> dict[int, str]:
    """{figure number: legend text} from a legends markdown file."""
    body = path.read_text()
    out: dict[int, str] = {}
    hits = list(re.finditer(pattern, body))
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(body)
        out[int(m.group(1))] = body[m.start():end]
    return out


def citation_order() -> tuple[list[dict], list[str]]:
    """Every Extended Data figure cited, and cited in ascending order.

    A figure nobody points at is a figure the reader never opens, and a
    numbering that does not follow first mention is one a copy-editor will
    renumber for you. Both were true of this package before S14c: two figures
    were never cited by any sentence, and the methods figure was numbered last
    and first cited in the third Results section. Checked mechanically here so
    it cannot come back — the citations are read from the stitched body, which
    is where a reader meets them, and the legend section is excluded so a
    legend heading cannot count as a citation of itself.
    """
    problems: list[dict] = []
    if not MANUSCRIPT_MD.exists():
        return [], ["manuscript.md not built — citation order not checked"]
    flat = re.sub(r"\s+", " ", MANUSCRIPT_MD.read_text(encoding="utf-8"))
    body = flat.split("Extended Data figure legends")[0]
    legends = {int(n) for n in
               re.findall(r"\*\*Extended Data Fig\. (\d+) \|", flat)}
    order: list[int] = []
    for n in (int(m) for m in re.findall(r"Extended Data Fig\. (\d+)", body)):
        if n not in order:
            order.append(n)
    msgs = []
    for num in sorted(legends - set(order)):
        msgs.append(f"Extended Data Fig. {num}: has a legend but no sentence "
                    f"in the paper cites it")
    for num in sorted(set(order) - legends):
        msgs.append(f"Extended Data Fig. {num}: cited but has no legend")
    if order != sorted(order):
        msgs.append(f"Extended Data figures are not numbered in order of "
                    f"first mention: {order}")
    rows = [{"figure": f"Extended Data Fig. {n}", "first_mention_rank": i + 1,
             "in_order": int(order == sorted(order))}
            for i, n in enumerate(order)]
    return rows, msgs


def audit() -> tuple[list[dict], list[str]]:
    """The mechanical half. Returns (rows, problems)."""
    manifest = [r for r in L.read_tsv(MANIFEST) if r["format"] == "png"]
    main = split_legends(LEGENDS_MAIN, r"\*\*Fig\. (\d+) \|")
    ed = split_legends(LEGENDS_ED, r"\*\*Extended Data Fig\. (\d+) \|")

    rows, problems = [], []
    for kind, legends, label in (("main", main, "Fig."),
                                 ("extended_data", ed, "Extended Data Fig.")):
        panels: dict[int, list[dict]] = {}
        for r in manifest:
            if r["kind"] == kind:
                panels.setdefault(int(r["number"]), []).append(r)
        for num in sorted(set(panels) | set(legends)):
            files = panels.get(num, [])
            text = legends.get(num, "")
            letters = _panel_letters(text)
            row = {"kind": kind, "number": num,
                   "panel_files": len(files),
                   "legend_panel_letters": "".join(letters),
                   "legend_present": bool(text),
                   "figures_present": all(
                       (MANUSCRIPT / "figures" /
                        f"{f['figure']}.png").exists() for f in files),
                   "note": ""}
            if not text:
                problems.append(f"{label} {num}: figure files but no legend")
                row["note"] = "no legend"
            elif not files:
                problems.append(f"{label} {num}: legend but no figure file")
                row["note"] = "no figure"
            elif not row["figures_present"]:
                problems.append(f"{label} {num}: a manifest file is missing")
                row["note"] = "missing file"
            elif kind == "extended_data" and letters and \
                    len(letters) != len(files):
                problems.append(
                    f"{label} {num}: legend references {len(letters)} panels "
                    f"({''.join(letters)}) but the manifest holds "
                    f"{len(files)} panel files")
                row["note"] = "panel count mismatch"
            elif letters and letters != [c for c in
                                         "abcdefgh"[:len(letters)]]:
                problems.append(
                    f"{label} {num}: panel letters {''.join(letters)} are not "
                    f"lowercase and in order")
                row["note"] = "panel letters"
            rows.append(row)
    return rows, problems


def write(out_dir: Path) -> dict:
    rows, problems = audit()
    order_rows, order_problems = citation_order()
    problems = problems + order_problems
    L.write_tsv(out_dir / "figure_audit.tsv", rows,
                ["kind", "number", "panel_files", "legend_panel_letters",
                 "legend_present", "figures_present", "note"])
    L.write_tsv(out_dir / "figure_citation_order.tsv", order_rows,
                ["figure", "first_mention_rank", "in_order"])
    findings = [dict(f) for f in FINDINGS]
    L.write_tsv(out_dir / "figure_findings.tsv", findings,
                ["figure", "panel", "legend_said", "figure_shows",
                 "verified_against", "status"])
    return {"figures_audited": len(rows), "problems": problems,
            "citations_checked": len(order_rows),
            "findings": len(findings),
            "findings_by_status": {
                s: sum(1 for f in findings if f["status"] == s)
                for s in sorted({f["status"] for f in findings})}}
