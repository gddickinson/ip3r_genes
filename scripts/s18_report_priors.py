"""S18's priors and hand-off — the third piece of the report split.

Kept apart from `s18_report_results.py` for the 500-line budget (the
`s3_report.py` / `s3_report_d10.py` split, applied twice as S7, S9 and S12 do)
and taking the caller's loader and formatters, so no half of the report can
read the tables differently from another.

Two verdicts here are the ones worth reading first. **`orthogonal`** on S10's
98.2 %: that number is measured over the loci where the annotation
demonstrably *could* have delivered the gene, and S18's is measured over every
locus in the scope, so the distance between them is the size of the population
S10 excluded rather than a disagreement. And **`contradicted`** on S5b's
23-point paralog naming gap, which does not survive reading the name off
whichever model the annotation actually places at the locus.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_priors as PR                                           # noqa: E402


def _row(rows, **kw):
    for r in rows:
        if all(str(r.get(k)) == str(v) for k, v in kw.items()):
            return r
    return {}


def priors(load, num, pct, H) -> list[str]:
    fvc = load("family_vs_control.tsv")
    pf = load("pfam_recall.tsv")
    audit = load("locus_audit.tsv")
    sc = [r for r in audit if r["gene_set"] == "present"
          and r["is_control"] == "0"]
    above = [r for r in sc if r["contig_spans_gene"] == "1"]

    def named(pool, cell):
        c = [r for r in pool if r["cell"] == cell]
        ok = sum(1 for r in c
                 if r["symbol_verdict"] == "correct_paralog"
                 or r["product_verdict"] == "correct_paralog")
        return ok, len(c)

    gaps = {}
    for cell in ("ITPR1", "ITPR2", "ITPR3"):
        ok, n = named(above, cell)
        gaps[cell] = ok / n if n else 0.0
    spread = max(gaps.values()) - min(gaps.values())
    both = _row(pf, axis="which instruments called it", bucket="both")
    d = [r for r in fvc if r["scope"].startswith("loci on a contig")
         and r["measure"] == "any annotation failure"]
    dd = d[0] if d else {}

    L = ["## 11. The priors this task was judged against", "",
         "Each is stated with where the earlier task said it, computed on "
         "S18's own tables, and rendered with both numbers printed either "
         "way.", ""]
    L += [PR.line("s10_annotation_correct",
                  f"{pct(H.get('frac_itpr_complete'))} complete over every "
                  f"scorable locus, and "
                  f"{pct(1 - float(H.get('frac_itpr_failure_above_d4', 0)))} "
                  f"over the loci above D4's bar",
                  "orthogonal")]
    L += ["  S10's denominator is the loci where the annotation *demonstrably "
          "could* have delivered the gene (five eligibility rules, of which "
          "E5 asks whether that annotation builds genes this long anywhere "
          "else in the genome); S18's is every locus in the scope. The gap "
          "between 98.2 % and "
          f"{pct(H.get('frac_itpr_complete'))} is the size of the population "
          "S10 excluded, and is a result rather than a disagreement."]
    L += [PR.line("s5b_dna_only_models",
                  f"{H.get('n_itpr_unannotated')} ITPR loci with no "
                  f"same-strand annotated feature and "
                  f"{H.get('n_itpr_noncoding')} held only by a non-coding "
                  f"one, in the {num(H.get('n_scorable_itpr_loci'))} loci "
                  f"whose assembly ships a gene set",
                  "confirmed")]
    L += [PR.line("s5b_paralog_naming_gap",
                  f"a {spread * 100:.1f}-point spread across ITPR1/2/3 "
                  f"({', '.join(f'{k} {v:.1%}' for k, v in gaps.items())}) "
                  f"once contiguity is held constant",
                  "contradicted")]
    L += ["  S5b measured naming against its own `annot_paralog_matches`, "
          "which requires a model covering half the locus; S18 reads the name "
          "off whichever model the annotation places there, coding or not. "
          "On that reading the paralogs are named equally well and S5b's "
          "23-point gap does not survive — what differs between the paralogs "
          "is whether a *coding* model exists, not whether the gene is named."]
    conflicts = [r for r in sc
                 if r["symbol_verdict"] == "correct_family_wrong_paralog"
                 or r["product_verdict"] == "correct_family_wrong_paralog"]
    raised = [r for r in load("corrections.tsv")
              if r["cls"] == "C5_wrong_paralog_name"]
    L += [PR.line("s5b_paralog_conflicts",
                  f"{len(conflicts)} loci whose covering model names a "
                  f"different paralog, of which {len(raised)} clears the "
                  f"recovery bar a correction needs",
                  "confirmed")]
    L += ["  The other "
          f"{len(conflicts) - len(raised)} are recovered at under "
          f"{min(float(r['coverage']) for r in conflicts if float(r['coverage']) < 0.5):.2f} "
          "of the bait, where an alignment covering a tenth of the gene is "
          "not evidence for renaming the model that covers the rest."]
    L += [PR.line("s3_zero_hit_proteomes",
                  f"{H.get('n_zero_hit_gene_caller')} of "
                  f"{H.get('n_zero_hit_proteomes')} resolved as gene-caller "
                  f"failures, 0 as absences",
                  "confirmed")]
    L += [PR.line("s2_naming_pfam_recall",
                  f"{pct(both.get('frac_with_naming_pfam'))} of the records "
                  f"both instruments call family carry PF08709",
                  "confirmed")]
    L += [PR.line("s0_zebrafish_ryr_share",
                  f"{H.get('n_name_wrong_family')} wrong-family names in "
                  f"{num(H.get('n_protein_records'))} full-length records, "
                  f"and 0 at any of the {num(H.get('n_loci'))} genomic loci",
                  "orthogonal")]
    L += ["  S0 measured what a *signature* returns; S18 measures what a "
          "*name* claims. A Pfam shared by both families says nothing about "
          "whether either is named correctly, and on this evidence both are."]
    L += [PR.line("s5b_contiguity_bar",
                  f"the ITPR failure rate falls from "
                  f"{pct(H.get('frac_itpr_any_failure'))} to "
                  f"{pct(dd.get('frac_itpr_failing'))} across that bar",
                  "confirmed")]
    L += [""]
    return L


def handoff(pct, H) -> list[str]:
    return [
        "## 12. What this settles, and what it does not",
        "",
        "**Settled.**",
        "",
        "- How the family is recorded across the declared genome scope, with "
        "the sister family measured beside it in the same assemblies: "
        f"{pct(H.get('frac_itpr_complete'))} of ITPR loci delivered as one "
        f"complete model against {pct(H.get('frac_control_complete'))} of RyR "
        f"loci, no difference surviving correction.",
        "- That the archive matters far more than the family: RefSeq and "
        "submitter-deposited GenBank gene sets differ on every state, and "
        "only part of that gap is assembly quality.",
        "- That an ITPR-shaped hole in a vertebrate reference proteome is, in "
        "every case in this scope, a gene caller and not a gene.",
        "- That the family's naming in the protein databases is essentially "
        "correct where a name exists, and that a name often does not: over "
        "half the full-length records have no usable gene symbol.",
        "",
        "**Not settled.**",
        "",
        "- **Whether a proposed correction is right.** The audit proposes; it "
        "does not validate. S10 validated two cases to the exon and found "
        f"both; the remaining {H.get('n_corrections')} items here carry "
        "alignment evidence and nothing more.",
        "- **The RYR3 cell.** The bait panel has no RYR3 bait, so 74 records "
        "whose symbol says RYR3 are outside the instrument's reach rather "
        "than adjudicated. That is a limit of S5's panel, recorded as one.",
        "- **Whether the non-coding demotions are wrong.** An annotation that "
        "files a gene as `gene_biotype=other` with a note about a frameshift "
        "may be right about the frameshift. S15's ORF screen vetoes the ones "
        "it can; the rest are proposed with that stated.",
        "- **Anything about the 33 assemblies with no gene set**, which are "
        "reported as unscorable rather than as failures.",
        "",
        "**What a reader should hold against it.** The completeness bar is "
        "inherited from S5 and sits in the far lower tail of the distribution "
        "it is applied to, so `complete` is a generous verdict — 15 loci read "
        "complete while a fifth of the gene reaches no model. And the whole "
        "genome half rests on the S5 alignments being right about where the "
        "exons are; S10's exon-boundary control validated that instrument at "
        "two loci in two genomes, not at 2,144.",
        "",
    ]
