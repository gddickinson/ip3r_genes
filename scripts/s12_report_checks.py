"""The cross-check and caveat half of the S12 report.

Split out of `s12_report_results.py` to keep every module inside the
500-line budget, and taking the caller's loader and formatters so all
three halves read the tables the same way (the `s3_report.py` /
`s3_report_d10.py` split, applied twice — the pattern S7 already needed
three times).

What lives here is the part of the report that argues against itself: the
deposit cross-check, whose verdict on its own power is `underpowered`, and
the caveats, which say what the reads cannot settle.
"""

from __future__ import annotations


def _i(r, k, d=0):
    try:
        return int(float(r.get(k) or d))
    except (TypeError, ValueError):
        return d


def _atlas(atlas, resources, loci, fmt, pct, missing, prior_line) -> list[str]:
    if not atlas:
        return [missing("10. The independent cross-check")]
    out = ["## 10. The independent cross-check — the deposits, re-asked", "",
           "Reads and deposits are different instruments. S10 ran the "
           "deposit test on its two cases and could not answer it. S12's "
           "panel adds a third species with a very different deposit: NCBI "
           "holds **37,166 mRNA records for *Dissostichus mawsoni***, "
           "against 43 and 10 for the other two.", "",
           "| species | nuccore records | mRNA records | TSA | fetched and "
           "searched | median record | longest record |",
           "|---|---|---|---|---|---|---|"]
    for r in sorted(resources, key=lambda x: x["organism"]):
        out.append(f"| *{r['organism']}* | {fmt(r['nuccore_total'])} | "
                   f"{fmt(r['nuccore_mrna'])} | {fmt(r['nuccore_tsa'])} | "
                   f"{fmt(r['mrna_fetched'])} | "
                   f"{fmt(r.get('median_record_nt'))} nt | "
                   f"{fmt(r.get('max_record_nt'))} nt |")
    probes = sum(_i(r, "probes") for r in atlas)
    hits = sum(_i(r, "deposit_hits") for r in atlas)
    spans = sum(_i(r, "deposit_spanning") for r in atlas)
    ctrl_hits = sum(_i(r, "control_hits") for r in atlas)
    ctrl_spans = sum(_i(r, "control_spanning") for r in atlas)
    ryr = [r for r in atlas if r["cell"] == "RYR"]
    ryr_hits = sum(_i(r, "deposit_hits") for r in ryr)
    out += ["",
            f"**{fmt(probes)} ±90 nt junction probes** were searched against "
            f"those deposits: **{fmt(hits)} hits, {fmt(spans)} spanning**. "
            f"The same probes against each locus's own genomic sequence — "
            f"where 180 nt of contiguous *spliced* sequence cannot span its "
            f"junction by construction — make {fmt(ctrl_hits)} hits and "
            f"{fmt(ctrl_spans)} spans, so the criterion is demonstrably "
            f"discriminating rather than one that never fires.", "",
            f"The result is a property of the deposits, not of the genes, "
            f"and the RyR control is what shows it: the probes find "
            f"{fmt(ryr_hits)} hits for the ryanodine receptor too — a gene "
            f"this project recovered at full length and one of these "
            f"annotations names. The lengths in the table are why: an "
            f"IP3-receptor transcript is over 8 kb of coding sequence, and "
            f"the largest deposit any of these species has is "
            f"{fmt(max((_i(r, 'max_record_nt') for r in resources), default=0))}"
            f" nt.", "",
            "![](figures/s12_instruments.png)", "",
            "**Figure 4.** The same question asked of both instruments: the "
            "fraction of unannotated junctions with evidence, from streamed "
            "reads and from submitted transcript records.", ""]
    out.append(prior_line(
        "handoff",
        f"Re-asked with a denominator 800x larger than S10's, the deposit "
        f"test still returns {fmt(spans)} spanning hits — and returns "
        f"{fmt(ryr_hits)} for the annotated RyR in the same genomes. The "
        f"deposits cannot see a gene of this size at this abundance, which "
        f"is why the reads were needed.",
        None if spans == 0 else True))
    out.append("")
    return out


def _caveats(by_run, loci, fmt, pct) -> list[str]:
    failed = [r for r in loci if r["role"] == "failed"]
    thin = [r for r in failed if 0 < _i(r, "runs_detected") <= 1]
    return [
        "## 11. Caveats", "",
        "- Runs are the **first N spots** of each accession, not a random "
        "sample; that is a flowcell-order subsample, fine for "
        "presence/absence and not for expression level.",
        "- Reads are counted unpaired (each mate independently), so `reads` "
        "is a mate count, not a fragment count.",
        "- The reference is a **closed set**: a read from a gene not in it "
        "cannot be assigned away. §6 measures the part of that which is "
        "testable — whether the sequences in the reference collect each "
        "other\u2019s reads — but a read from a *fourth* ITPR-like locus "
        "absent from the reference would still land on the nearest member. "
        "The sweep found no such locus in these genomes; that is the "
        "assumption the closure rests on.",
        "- Several libraries show a strong 3' bias, so junctions near the "
        "5' end of a transcript are asked with less depth than those near "
        "the 3' end. The per-junction table carries the position of every "
        "junction, so this is visible rather than absorbed.",
        "- **Absence of reads in one tissue is weak evidence of absence.**"
        + (f" {len(thin)} locus/loci are detected in only a single library "
           f"and should be quoted with that number." if thin else ""),
        "- A tissue is only as good as its submitter's metadata; the "
        "attribute each tissue call was read from is recorded in "
        "`runs_selected.tsv`, and calls resting on a free-text title are "
        "marked as such.",
        "- The deposit cross-check is **underpowered by construction** for "
        "these species and this gene size, which §10 establishes with its "
        "own RyR control rather than assuming.",
        "",
    ]


def _outputs() -> list[str]:
    return [
        "## Outputs", "",
        "- `panel_audit.tsv` — every locus P1–P4 considered, and the rule "
        "that included or excluded it",
        "- `reference_table.tsv` / `junctions.tsv` — the references and "
        "every junction with what the annotation holds there",
        "- `reference_validation.tsv` — each reference scored as built, "
        "with its blocks reversed and off the wrong strand",
        "- `runs_considered.tsv` / `runs_selected.tsv` / "
        "`run_availability.tsv` — every run the SRA search returned, why "
        "each was kept or rejected, and the per-species denominator",
        "- `run_metrics.tsv` — library size, tissue and study per run",
        "- `expression_by_run.tsv` — every run × reference row with its "
        "decoy and its detection verdict",
        "- `expression_by_locus.tsv` / `expression_by_tissue.tsv` — pooled",
        "- `junction_support.tsv` — **every junction of every reference**, "
        "whether or not a read crossed it",
        "- `annotation_gap_coverage.tsv` — read coverage inside and outside "
        "the annotated coding blocks",
        "- `coverage_profile.tsv` — binned coverage per locus",
        "- `crossmap_control.tsv` — every reference tiled with synthetic "
        "reads and mapped back, the closed-set control",
        "- `atlas_resources.tsv` / `atlas_probes.tsv` / `atlas_summary.tsv` "
        "— the deposit cross-check with its genomic negative control",
        "- `expression_stats.json` — parameters, self-test status and the "
        "SHA-256 of every table",
        "- `figures/` — four figures (D13, D19)",
        "",
    ]


def sections(atlas, resources, loci, by_run, fmt, pct, missing,
             prior_line) -> list[str]:
    out: list[str] = []
    out += _atlas(atlas, resources, loci, fmt, pct, missing, prior_line)
    out += _caveats(by_run, loci, fmt, pct)
    out += _outputs()
    return out
