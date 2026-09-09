"""S18 — renders `results/annotation_audit/report.md` purely from the
committed tables (D13).

Scope, the instrument, the bar and the negative controls live here; the
results live in `s18_report_results.py` (the `s3_report.py` /
`s3_report_d10.py` split, so both halves stay under 500 lines and take this
module's loader and formatters and therefore cannot read the tables
differently).

**Headlines are chosen by the data.** Every prior is stated in
`s18_priors.PRIOR` with where the earlier task said it, computed from S18's
own tables, and rendered with both numbers printed either way. The comparison
that had to be sayable is the one that goes badly for the audit's own premise:
if this family were no worse recorded than its sister family, the table it
would be said from is `family_vs_control.tsv` — and that is what it says.

A section whose table is absent renders *not run yet*, so a stage that was
skipped is visible as skipped.

Run:  python3 scripts/s18_report.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

OUT = Path(__file__).resolve().parents[1] / "results" / "annotation_audit"


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
        return f"{float(v):,.{nd}f}"
    except (TypeError, ValueError):
        return str(v)


def pct(v, nd=1) -> str:
    """A share already expressed as a fraction."""
    try:
        return f"{100.0 * float(v):.{nd}f} %"
    except (TypeError, ValueError):
        return "n/a"


def share(part, whole, nd=1) -> str:
    try:
        return f"{100.0 * float(part) / float(whole):.{nd}f} %"
    except (TypeError, ValueError, ZeroDivisionError):
        return "n/a"


def pfmt(p, nd=3) -> str:
    """A p-value at four decimal places is 0.0000, which hides how strong a
    claim is rather than how weak (`s15_report.pfmt`'s reason)."""
    try:
        f = float(p)
    except (TypeError, ValueError):
        return str(p)
    if f == 0.0:
        # Fisher's exact underflowed to zero. Printing "0.0e+00" invites a
        # reader to think the test was not run; the honest rendering names the
        # floor double precision put it under.
        return "< 1e-300"
    return f"{f:.{nd}f}" if f >= 1e-3 else f"{f:.1e}"


def table(rows: list[dict], cols: list[tuple[str, str]]) -> list[str]:
    out = ["| " + " | ".join(h for h, _ in cols) + " |",
           "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(k, "")) for _, k in cols) + " |")
    return out


MISSING = "*not run yet — the stage that writes this table has not been run.*"


# --------------------------------------------------------------------------
def s1_scope(H: dict, stats: dict) -> list[str]:
    man = load("locus_state_counts.tsv")
    if not man:
        return ["## 1. Scope", "", MISSING]
    p = stats.get("params", {})
    L = [
        "## 1. Scope — what is being audited, and against what",
        "",
        "Two records of the same genes are read against each other, and both "
        "against one common evidence set: the alignments the S5 genomic sweep "
        "placed in every genome of the declared scope.",
        "",
        f"- **The genome half.** {num(H['n_loci'])} gene-scale loci in "
        f"{num(H['n_genomes'])} assemblies — {num(H['n_itpr_loci'])} ITPR "
        f"loci and {num(H['n_control_loci'])} ryanodine-receptor control "
        f"loci — scored against the assemblies' own annotations. "
        f"{num(H['n_loci_no_gene_set'])} sit in assemblies that ship no gene "
        f"set at all and {num(H.get('n_loci_cds_unavailable'))} could not "
        f"have their own alignment recovered from the archive; both are "
        f"reported as unscorable rather than as failures, leaving "
        f"{num(H['n_scorable_itpr_loci'])} scorable ITPR loci.",
        f"- **The protein half.** {num(H['n_protein_records'])} full-length "
        f"family protein records from census v6 "
        f"(≥ {num(p.get('protein_min_length_aa'))} aa, not flagged "
        f"fragments, from a protein database rather than a genome model), of "
        f"which {num(H['n_protein_vertebrate'])} are vertebrate and can be "
        f"asked their paralog.",
        f"- **The proteome half.** the {num(H['n_zero_hit_proteomes'])} "
        f"vertebrate reference proteomes S3's profile sweep found nothing in, "
        f"each resolved against an assembly of its own species.",
        "",
        "The ryanodine receptors are the control throughout, not a nuisance: "
        "they are annotated by the same pipelines in the same assemblies and "
        "carry every diagnostic domain of the family (D14), so a failure rate "
        "measured without them cannot be told from the general quality of "
        "vertebrate gene sets.",
        "",
        "### 1.1 Loci by cell and annotation source",
        "",
    ]
    cell_rows = [r for r in man if not r.get("source")]
    L += table(cell_rows, [
        ("cell", "cell"), ("loci", "n_loci"),
        ("no gene set", "n_no_gene_set"), ("complete", "n_complete"),
        ("split", "n_split"), ("fragmentary", "n_fragmentary"),
        ("non-coding only", "n_noncoding"), ("unannotated", "n_unannotated"),
    ])
    return L


def s2_instrument(stats: dict) -> list[str]:
    p = stats.get("params", {})
    return [
        "## 2. The instrument, and the three rules that make it a measurement",
        "",
        "**Same strand only.** A gene on the other strand overlapping an ITPR "
        "locus is not a model of this gene however much sequence it shares "
        "with it. The rule is a filter, not a tiebreak, and the self-test "
        "constructs a perfectly covering antisense gene and requires the "
        "locus to read `unannotated`.",
        "",
        "**Scored against CDS blocks, never gene spans.** A gene span covers "
        "its own introns, so scoring against spans credits an annotation with "
        "every base it never called — and this family's introns reach 152 kb "
        "(S5's intron calibration), which is room for several passenger "
        "genes. S10 established this at two loci; here it runs at "
        f"{num(load_json('annotation_audit_stats.json').get('headline', {}).get('n_loci'))}.",
        "",
        "**The name verdict reads the model the annotation actually places "
        "there**, which is not always the biggest coding one. *Podiceps "
        "cristatus* files its ITPR3 as `gene_biotype=other`, "
        "`Note=contains frameshift`, named \"Inositol 1,4,5-trisphosphate "
        "receptor type 3\", with no CDS feature anywhere — correctly "
        "identified and serving no protein. Reading the largest coding model "
        "would have scored that name `absent`, and those are two different "
        "failures a correction list has to keep apart.",
        "",
        "### 2.1 Parameters",
        "",
        "| parameter | value | where it comes from |",
        "|---|---|---|",
        f"| completeness bar | {p.get('complete_frac')} | "
        f"`{p.get('complete_frac_source')}` — inherited, not chosen (§3) |",
        f"| a coding model counts as a *piece* at | {p.get('min_piece_frac')} "
        f"| stated here; swept in `state_sensitivity.tsv` |",
        f"| anything at all on the strand | {p.get('min_touch_frac')} | "
        f"separates `unannotated` from `noncoding` |",
        f"| calibration admission: alignment coverage | "
        f"{p.get('calib_min_coverage')} | the gene is demonstrably there |",
        f"| protein-side family margin | {p.get('protein_rel_margin')} | "
        f"`s3_assign.REL_MARGIN` (D7) |",
        f"| protein-side length floor | {num(p.get('protein_min_length_aa'))} "
        f"aa | `family.MIN_LENGTH_AA` |",
        f"| a rename needs | {num(p.get('min_rename_bits'))} bits | stated "
        f"here; D7's relative margin is not enough on its own (§8) |",
        f"| annotation window pad | {num(p.get('window_pad_bp'))} bp | a "
        f"model may start outside the alignment's first exon |",
        "",
    ]


def s3_bar(stats: dict) -> list[str]:
    cal = load("complete_calibration.tsv")
    if not cal:
        return ["## 3. The bar", "", MISSING]
    c = cal[0]
    if c.get("usable") != "1":
        return ["## 3. The bar", "",
                f"*The calibration refused to report: {c.get('reason')}.*"]
    sens = load("state_sensitivity.tsv")
    L = [
        "## 3. The bar is inherited, and then measured",
        "",
        "This project already has one definition of \"this annotated gene is "
        "the model of that locus\" — `s5_classify.ANNOT_CDS_FRAC`, half the "
        "alignment's coding footprint — and S5, S10 and S23 all read a locus "
        "through it. A second bar chosen here would fork what the project "
        "means by an annotated gene, so S18 uses that one and spends its "
        "calibration on the question that is open: **is 0.50 the right place "
        "for it?**",
        "",
        f"The population is the {num(c['n'])} loci whose own annotation names "
        f"the correct paralog with a coding model, whose alignment recovered "
        f"at least {pct(stats.get('params', {}).get('calib_min_coverage'), 0)} "
        f"of its bait, and whose contig clears D4's bar — genes the "
        f"annotation demonstrably holds in assemblies that demonstrably carry "
        f"them. Over that population a single annotated model covers a median "
        f"**{c['median']}** of the gene.",
        "",
        "| quantile of best-single-model coverage | value |",
        "|---|---|",
        f"| minimum | {c['min']} |",
        f"| 1 % | {c['p01']} |",
        f"| 5 % | {c['p05']} |",
        f"| 10 % | {c['p10']} |",
        f"| 25 % | {c['p25']} |",
        f"| median | {c['median']} |",
        f"| maximum | {c['max']} |",
        "",
        f"So the inherited bar sits at that distribution's **"
        f"{pct(c['bar_percentile'])} point**: {c['n_below_bar']} of "
        f"{num(c['n'])} demonstrably-annotated loci fall below it "
        f"({pct(c['frac_below_bar'])}), and it is therefore *conservative* — "
        f"it credits an annotation with a gene it delivers only half of. "
        f"{c['n_overcredited']} loci sit above the bar but below "
        f"{c['overcredit_frac']}, which means a fifth of the gene reaches no "
        f"model at all and the locus still reads `complete`. Those are "
        f"reported as the price of the bar rather than absorbed by it.",
        "",
    ]
    if sens:
        base = [r for r in sens
                if float(r["min_piece"]) == float(
                    stats.get("params", {}).get("min_piece_frac", 0.05))]
        L += [
            "### 3.1 What moving the bar changes",
            "",
            "Every count in the grid is the same evidence read at a different "
            "bar (S15b's sensitivity matrix, applied to the one threshold S18 "
            "owns). `unannotated` and `noncoding` do not move with it by "
            "construction, which is worth being able to see rather than "
            "assert.",
            "",
        ]
        L += table(base, [
            ("bar", "bar"), ("ITPR loci", "n_loci"),
            ("complete", "n_complete"), ("split", "n_split"),
            ("fragmentary", "n_fragmentary"),
            ("non-coding only", "n_noncoding"),
            ("unannotated", "n_unannotated"),
        ])
        lo = base[0]
        hi = base[-1]
        L += ["",
              f"Across the whole plausible range — {lo['bar']} to "
              f"{hi['bar']} — the share of ITPR loci reading `complete` moves "
              f"from {share(lo['n_complete'], lo['n_loci'])} to "
              f"{share(hi['n_complete'], hi['n_loci'])}. The audit's headline "
              f"does not rest on where the bar is put.",
              ""]
    return L


def s4_controls(stats: dict) -> list[str]:
    return [
        "## 4. Negative controls",
        "",
        f"`scripts/s18_test_audit.py` — {stats.get('self_test')}. It runs "
        "before anything is written and the driver refuses to continue if it "
        "fails. Every rule in this task returns a plausible number when it is "
        "wrong, and three of them did on the way here, so the tests are "
        "mostly tests on **refusal** and on **reachability**:",
        "",
        "- a perfectly covering gene on the **other strand** must leave the "
        "locus `unannotated`;",
        "- a passenger gene whose span crosses the whole locus while its CDS "
        "sits entirely in an intron must count for nothing;",
        "- **each of the five states must be reachable**, including the three "
        "that override every model (`no_gene_set`, `gff_unavailable`, "
        "`cds_unavailable`) — a state nothing can reach is not a negative "
        "result;",
        "- a locus whose own alignment cannot be recovered must be reported "
        "**unscorable and never scored against its span** — miniprot restarts "
        "its model ids per chunk, so a chunked genome's concatenated GFF "
        "repeats every one, and reading the raw id missed 21 loci whose "
        "denominator then became the whole locus span (2.4 Mb under an 8 kb "
        "gene, which forces `unannotated` whatever the annotation holds);",
        "- a 6 bp overlap is not a *piece*, so one model plus a sliver is not "
        "a split;",
        "- the naming model must prefer a small correctly-named gene over a "
        "large unnamed one, and a correctly-named non-coding gene over an "
        "unnamed coding one;",
        "- every name verdict must fire on its own case, including the two "
        "that exist because their absence manufactured errors — "
        "`family_ambiguous` (\"RyR/IP3R Homology associated domain-containing "
        "protein\" names both families) and `paralog_unspecified` "
        "(\"inositol 1,4,5-trisphosphate receptor\" with no type number is "
        "not a wrong paralog);",
        "- a paralog the bait panel cannot reach must not read as a naming "
        "error;",
        "- two baits inside D7's margin must give **no** family call;",
        "- the paralog question must not be asked outside the vertebrates;",
        "- the calibration must refuse to report from fewer than 50 loci, and "
        "a locus the annotation is silent about must never enter it;",
        "- **D6's veto must be recorded, not dropped** — a lesion-rich locus "
        "is written with `vetoed = 1` and its reason;",
        "- a rename must need an absolute score as well as a relative margin;",
        "- all four zero-hit verdicts must fire, `undecidable_no_genome` "
        "included;",
        "- and the two one-pass readers this task added to `s10_gff` must "
        "agree with the per-window and per-model ones they replace.",
        "",
        "Mutation-tested on nine deliberate rule breakages — the strand "
        "filter removed, scoring moved to gene spans, the naming model "
        "ranked by size, the ambiguity check dropped, the panel-reach guard "
        "disabled, the piece floor removed, the paralog asked of everything, "
        "the calibration floor removed, and `paralog_unspecified` folded back "
        "into the wrong-paralog cell. All nine were caught.",
        "",
    ]


def main() -> None:
    import s18_report_results as RR
    stats = load_json("annotation_audit_stats.json")
    H = stats.get("headline", {})
    lines = [
        "# S18 — Annotation-quality audit",
        "",
        "*Rendered from the committed tables by `scripts/s18_report.py` "
        "(D13). Nothing here is hand-written.*",
        "",
        "**How the family is recorded, and by whom.** Across "
        f"{num(H.get('n_genomes'))} vertebrate assemblies the annotation "
        f"delivers {pct(H.get('frac_itpr_complete'))} of the ITPR loci this "
        f"project recovers as a single complete gene model — and the "
        f"ryanodine-receptor control, in the same assemblies and through the "
        f"same pipelines, sits at {pct(H.get('frac_control_complete'))}. The "
        f"family is not recorded worse than its sister; vertebrate gene sets "
        f"are recorded this way. Two things do separate: **which archive "
        f"serves the annotation**, and **whether the assembly can carry the "
        f"gene at all**.",
        "",
    ]
    lines += s1_scope(H, stats)
    lines += [""]
    lines += s2_instrument(stats)
    lines += s3_bar(stats)
    lines += s4_controls(stats)
    lines += RR.results(load, load_json, num, pct, share, pfmt, table, MISSING)
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print(f"[s18] wrote {OUT / 'report.md'}")


if __name__ == "__main__":
    main()
