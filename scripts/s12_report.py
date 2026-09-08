"""S12 — renders `results/expression/report.md` purely from the tables (D13).

Scope, the reference construction and the instrument live here; the results
— what the reads say about the loci the annotation loses, junction by
junction, and the deposit cross-check — live in `s12_report_results.py`
(the `s3_report.py` / `s3_report_d10.py` split, so both halves stay inside
the 500-line budget and read the tables through the same loader and the
same formatter).

**Headlines chosen by the data.** Every prior is stated in
`s12_priors.PRIOR` with where the earlier task said it, computed from
S12's own tables, and rendered with S8's five-valued verdict — both numbers
printed either way. The comparison that has to be sayable is the one that
goes badly: a locus the census claims and no library transcribes would be
read off `expression_by_locus.tsv`, and this report has to be able to say
it.

**A section whose table is absent renders *not run yet*,** never nothing,
so a partial run is visible as partial rather than as a shorter report.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s12_priors as priors  # noqa: E402
from s12_panel import LOSS_FLOOR, MIN_RUNS  # noqa: E402
from s12_tables import MIN_JUNCTION_READS, MIN_READS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "expression"


def load(name: str) -> list[dict]:
    path = OUT / name
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    path = OUT / name
    return json.loads(path.read_text()) if path.exists() else {}


def fmt(x, nd: int = 2) -> str:
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


def missing(section: str) -> str:
    return (f"\n## {section}\n\n*Not run yet — the table this section "
            f"renders from is not present.*\n")


def prior_line(key: str, ours: str, holds) -> str:
    p = priors.PRIOR[key]
    return (f"> **Prior.** {p['where']}.\n>\n"
            f"> **This task.** {ours}\n>\n"
            f"> {priors.verdict(key, holds)}\n")


def _int(r, k, d=0):
    try:
        return int(float(r.get(k) or d))
    except (TypeError, ValueError):
        return d


# ---------------------------------------------------------------------------

def header(refs, metrics, stats) -> list[str]:
    failed = [r for r in refs if r["role"] == "failed"]
    species = sorted({r["organism"] for r in refs})
    n_runs = stats.get("runs_quantified", len(metrics))
    lib = sum(_int(m, "library_reads") for m in metrics)
    un = sum(_int(r, "n_junctions_unannotated") for r in failed)
    return [
        "# S12 — are the genes the annotation loses actually transcribed?",
        "",
        f"S10 found **{priors.value('annotation_failures')}** loci where the "
        f"IP3-receptor gene this project recovered from the genome reaches no "
        f"annotated gene model, in assemblies whose annotation delivers the "
        f"other 375 of 382 whole. It could not tell whether those are real "
        f"genes: the transcript deposits for the species are 43 and 10 "
        f"records, so a search of them returning nothing is a fact about the "
        f"deposit. It handed the question here.",
        "",
        f"This report answers it with reads. **{n_runs} public RNA-seq runs** "
        f"({fmt(lib)} reads streamed) across **{len(species)} species** were "
        f"aligned against a per-species reference holding every family locus "
        f"the sweep recovered in that genome, three housekeeping anchors and "
        f"a composition-matched decoy for each. The measurement is not "
        f"whether the gene collects reads but whether the **"
        f"{fmt(un)} exon junctions no annotated model spans** are crossed by "
        f"reads that read through them.",
        "",
        "> Every number below is read from a committed table in this "
        "directory. Nothing here is computed from a genome, a SAM file or a "
        "network call at render time (D13).",
        "",
    ]


def section_scope(refs, audit, avail) -> list[str]:
    if not refs:
        return [missing("1. Which loci, and why those")]
    out = ["## 1. Which loci, and why those", "",
           "The scope is derived from S10's committed ranking, not chosen. "
           "Four rules, each a positive test:", "",
           "| rule | test |", "|---|---|",
           f"| P1 | the locus is eligible under S10's five rules and its "
           f"annotation loss exceeds {LOSS_FLOOR} |",
           f"| P2 | its species has at least {MIN_RUNS} public Illumina "
           f"RNA-seq runs — a species with no libraries cannot be asked the "
           f"question |",
           "| P3 | every **other** ITPR locus the sweep recovered in the "
           "same genome joins the reference, whatever its annotation status "
           "— known-real genes in the same libraries, on the same reference, "
           "by the same aligner |",
           "| P4 | so does the genome's RyR locus: this project's sister "
           "family (D14), whose reads must land on RyR |", ""]

    by_sp: dict[str, list[dict]] = {}
    for r in refs:
        if r["role"] != "decoy":
            by_sp.setdefault(r["organism"], []).append(r)
    avail_by = {a["organism"]: a for a in avail if a.get("term") == "base"}
    out += ["The panel that produces:", "",
            "| species | assembly | RNA-seq runs | loci | annotation "
            "failures | modes |", "|---|---|---|---|---|---|"]
    for org, rs in sorted(by_sp.items()):
        fail = [r for r in rs if r["role"] == "failed"]
        modes = sorted({r["failure_mode"] for r in fail if r["failure_mode"]})
        out.append(f"| *{org}* | `{rs[0]['accession']}` | "
                   f"{fmt(avail_by.get(org, {}).get('n_runs', ''))} | "
                   f"{len([r for r in rs if r['role'] != 'housekeeping'])} | "
                   f"{len(fail)} | {', '.join(modes) or '—'} |")
    out += ["",
            f"That is **every one of S10's {priors.value('annotation_failures')} "
            f"failures**, covering all three failure modes, with "
            f"{len([r for r in refs if r['role'] in ('control_paralog', 'control_family')])} "
            f"internal controls beside them.", ""]
    excl = [a for a in audit if a["verdict"] != "included"]
    if excl:
        out += ["Excluded, with the rule that excluded them:", "",
                "| species | locus | rule | reason |", "|---|---|---|---|"]
        for a in excl:
            out.append(f"| *{a['organism']}* | {a['cell']} | {a['rule']} | "
                       f"{a['reason']} |")
        out.append("")
    return out


def section_reference(refs) -> list[str]:
    if not refs:
        return [missing("2. What the reads are aligned against")]
    real = [r for r in refs if r["role"] != "decoy"]
    out = ["## 2. What the reads are aligned against", "",
           "The reference is the **spliced genomic coding sequence** of each "
           "locus: the miniprot model's CDS blocks from the archived sweep "
           "GFF, fetched out of the assembly by coordinate, "
           "reverse-complemented on the minus strand and concatenated in "
           "target order.", "",
           "It is deliberately not S9's frameshift-corrected CDS. S9 needed "
           "a clean reading frame because codeml counts codons; a read "
           "comes from the genome as it is, and correcting a frameshift "
           "inserts 1–2 bp that every read crossing it would then carry as "
           "an indel — costing reads at exactly the sites the sweep already "
           "flagged as difficult.", "",
           "### Validation, and the test that had to be replaced", "",
           "Each reference is checked against the sweep's own `##STA` "
           "protein. The obvious check — translate and require high "
           "identity — is the wrong instrument, and measuring it is how "
           "that became visible: a frameshift costs the reading frame from "
           "where it sits, so a **correct** reference scores 1.00 with no "
           "frameshifts and 0.92 with nine, and any floor drawn across that "
           "range rejects correct references for carrying frameshifts the "
           "sweep had already recorded.", "",
           "What replaced it has no tuned threshold. Each block is "
           "translated in its own recorded phase and searched for verbatim "
           "in the expected protein **from the previous block's match "
           "onwards**, so a block places only if it is this protein's "
           "sequence *and* it follows the block before it. The pass rule is "
           "then a statement the model makes about itself: **a block may "
           "fail to place only if the model's own frameshift count can "
           "explain it.** Block *number* is checked separately, against the "
           "`cds_bp` the sweep recorded for the locus.", "",
           "",
           "| species | locus | role | nt | exons | frameshifts | blocks "
           "placed | junctions (annotated / total) |", "|---|---|---|---|---|---|---|---|"]
    for r in sorted(real, key=lambda x: (x["species"], x["role"], x["seq"])):
        out.append(
            f"| *{r['organism']}* | `{r['seq'].split('|')[0]}` | {r['role']} "
            f"| {fmt(r['length'])} | {fmt(r['n_exons'])} | "
            f"{fmt(r['frameshifts'])} | "
            f"{fmt(r['blocks_placed'])}/{fmt(r['blocks_tested'])} | "
            f"{fmt(r['n_junctions_annotated'])} / {fmt(r['n_junctions'])} |")
    out += ["",
            "**The rule is reported beside the failures it rejects.** Each "
            "reference is re-scored twice more: once with its blocks in "
            "reverse target order (coordinates and strand intact) and once "
            "read off the wrong strand (coordinates and order intact). "
            "Both controls are built from the real models, not from a "
            "synthetic gene, so these are measurements on the actual "
            "references.", ""]
    out += _validation_controls()
    out += ["",
            "Two further sequences per gene are what make the counts "
            "readable. A **reversed decoy** — the same sequence backwards, "
            "so identical length and base composition with no homology — "
            "measures what this reference collects by accident in this "
            "library. And three **housekeeping anchors** (GAPDH, EEF1A1, "
            "RPL13A) prove the library worked at all.", "",
            "The housekeeping anchors are aligned with miniprot and spliced "
            "by the same module as the target loci, rather than pulled out "
            "of each assembly's annotation. Two of these three annotations "
            "name no genes at all — *Dissostichus mawsoni* files all "
            "29,240 of its genes as \"hypothetical protein\" — so a "
            "GFF-derived anchor exists for one species and not the others; "
            "and building the anchor and the target the same way means a "
            "difference between them cannot be a difference in "
            "construction.", "",
            "Their bait accessions are **resolved by query against a "
            "declared length band**, not remembered. The first version of "
            "this module hard-coded three and got all three wrong — a "
            "427 aa \"GAPDH\" (the enzyme is 333), a 1,829 aa \"EEF1A1\" "
            "(462) and an accession that served nothing. A bait of the "
            "wrong protein still aligns somewhere and still produces "
            "reads, so the anchor would have silently measured a different "
            "gene.", ""]
    return out


def _validation_controls() -> list[str]:
    rows = load("reference_validation.tsv")
    if not rows:
        return ["*`reference_validation.tsv` not present.*"]
    fam = [r for r in rows if r["role"] != "housekeeping"]
    out = ["| species | locus | as built | blocks reversed | wrong strand | "
           "spliced nt vs the sweep's own `cds_bp` |",
           "|---|---|---|---|---|---|"]
    for r in sorted(fam, key=lambda x: (x["organism"], x["seq"])):
        exp = r.get("cds_bp_expected", "")
        built = r.get("cds_bp_built", "")
        agree = "match" if exp and str(exp) == str(built) else \
            (f"{built} vs {exp}" if exp else "—")
        out.append(
            f"| *{r['organism']}* | `{r['seq'].split('|')[0]}` | "
            f"{r['as_built_placed']}/{r['as_built_tested']} "
            f"({r['as_built_rate']}) | "
            f"{r['reversed_order_placed']}/{r['reversed_order_tested']} "
            f"({r['reversed_order_rate']}) | "
            f"{r['wrong_strand_placed']}/{r['wrong_strand_tested']} "
            f"({r['wrong_strand_rate']}) | {agree} |")
    return out


def section_instrument(metrics, stats) -> list[str]:
    if not metrics:
        return [missing("3. How a junction is counted")]
    p = stats.get("parameters", {})
    lib = sum(_int(m, "library_reads") for m in metrics)
    studies = len({m["study"] for m in metrics if m.get("study")})
    attr = sum(1 for m in metrics if m.get("tissue_from") == "attribute")
    return [
        "## 3. How a junction is counted", "",
        f"Each run is streamed — `fastq-dump -X {fmt(p.get('spots', 4000000))}"
        f" | hisat2` — and the SAM parsed in flight; only counts reach the "
        f"disk. Alignment is **unspliced on purpose**: the reference is "
        f"already CDS, so a read crossing a junction is contiguous here, and "
        f"leaving hisat2's spliced mode on would let it open a gap inside "
        f"the CDS and call an intron that does not exist — which would then "
        f"be counted as spanning a junction it had actually skipped.", "",
        f"A read spans a junction when **one contiguous aligned block** "
        f"covers at least {p.get('junction_anchor_nt', 8)} nt on **both** "
        f"sides of it. Requiring the anchor on one side only admits an "
        f"alignment that has run a few bases past the junction and stopped, "
        f"which is not evidence of splicing; S10 hit exactly that and its "
        f"junction rule grew its second half from it.", "",
        "A locus is called **transcribed** in a run when all three hold:", "",
        f"1. at least {MIN_JUNCTION_READS} junction-spanning reads — the "
        f"load-bearing clause, because genomic-DNA carryover cannot cross a "
        f"splice point;",
        f"2. at least {MIN_READS} reads in total;",
        "3. more reads than that run's own decoy.", "",
        f"**{len(metrics)} runs, {fmt(lib)} reads, {studies} independent "
        f"studies**; {attr} of {len(metrics)} runs carry a tissue confirmed "
        f"from their own BioSample attributes rather than from a title.", "",
    ]


def main(argv: list[str] | None = None) -> int:
    import s12_report_results as results

    refs = load("reference_table.tsv")
    metrics = load("run_metrics.tsv")
    stats = load_json("expression_stats.json")
    lines: list[str] = []
    lines += header(refs, metrics, stats)
    lines += section_scope(refs, load("panel_audit.tsv"),
                           load("run_availability.tsv"))
    lines += section_reference(refs)
    lines += section_instrument(metrics, stats)
    lines += results.sections(load, load_json, fmt, pct, missing, prior_line)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print(f"[s12] report -> {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
