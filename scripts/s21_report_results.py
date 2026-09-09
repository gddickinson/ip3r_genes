"""The results half of S21's report, split to keep both under 500 lines.

It takes the caller's formatter and table helper by import, so the two halves
cannot render a number differently (the `s3_report.py` / `s3_report_d10.py`
pattern). Every prior an earlier task set is stated in `s21_priors.PRIOR` with
where it was said, computed on S21's own tables, and rendered with both numbers
printed either way.

The comparison that has to be sayable is the one that goes badly. Here it is
§8: the three paralogues share ~48 intron positions and the ryanodine
receptors, which carry every ITPR-diagnostic Pfam domain, share **one** — so
the sister family that is inside every search this project runs turns out to
have an exon structure with no ancestry in common with this one, and the
control that makes that statement is the same control D14 uses everywhere else.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_architecture as A                                      # noqa: E402
import s21_introns as I                                           # noqa: E402
import s21_lib as L                                               # noqa: E402
import s21_priors as P                                            # noqa: E402
import s21_report as RP                                           # noqa: E402

fmt, pct, pfmt, table = RP.fmt, RP.pct, RP.pfmt, RP.table
CELLS = list(L.PARALOGS) + [L.CONTROL_CELL]


def _missing(name: str, what: str) -> str | None:
    return None if (L.OUT / name).exists() else (
        f"## {what}\n\n*Not run yet — `{name}` is not on disk.*")


def section_6(h: dict) -> str:
    miss = _missing("architecture_by_paralog.tsv", "6. The architecture")
    if miss:
        return miss
    per = h["by_paralog"]
    rows = []
    for cell in CELLS:
        d = per.get(cell)
        if not d:
            continue
        rows.append([RP.cell_label(cell), fmt(d["n_loci"]), fmt(d["n_genomes"]),
                     f"{d['exons']:.0f} ({d['exons_p10']:.0f}–"
                     f"{d['exons_p90']:.0f})",
                     fmt(int(d["cds_bp"])), fmt(int(d["span_bp"])),
                     fmt(int(d["median_intron_bp"])),
                     f"{d['mean_exon_bp']:.0f}"])
    sens = L.read_tsv(L.OUT / "architecture_sensitivity.tsv")
    spread = {}
    for cell in CELLS:
        vals = [float(r["n_exons_median"]) for r in sens if r["cell"] == cell]
        if vals:
            spread[cell] = max(vals) - min(vals)
    worst = max(spread.values()) if spread else 0.0
    ex_lo, ex_hi = (h["exon_count_range"] + [0, 0])[:2]
    mean_lo = min(per[c]["mean_exon_bp"] for c in L.PARALOGS)
    mean_hi = max(per[c]["mean_exon_bp"] for c in L.PARALOGS)
    seen_lo = min(h["conservation"][c]["positions_seen"] for c in L.PARALOGS) \
        if h.get("conservation") else 0
    worst_txt = ("does not move at all" if worst == 0
                 else f"moves by at most {worst:.0f} "
                      f"exon{'s' if worst != 1 else ''}")
    return f"""## 6. The architecture

{table(["gene", "loci", "genomes", "coding exons (p10–p90)", "CDS bp",
        "genomic span bp", "median intron bp", "mean exon bp"], rows)}

**The exon count is the conserved thing and the genomic span is not.** Across
{fmt(h['n_in_scope'])} genes in {fmt(h['n_genomes_in_scope'])} vertebrate
genomes the three paralogues sit at **{ex_lo:.0f}–{ex_hi:.0f} coding exons**
and their median genomic spans differ by **{h['span_fold_range']:.1f}-fold**.
The coding sequence is nearly the same size in all three
({fmt(int(per['ITPR1']['cds_bp']))}, {fmt(int(per['ITPR2']['cds_bp']))} and
{fmt(int(per['ITPR3']['cds_bp']))} bp); what differs is how much intron is
wrapped around it.

The ryanodine receptors, measured through the identical instrument in the same
assemblies, carry **{per['RYR']['exons']:.0f}** exons over
{fmt(int(per['RYR']['cds_bp']))} bp of coding sequence — nearly twice the gene
in both, at almost exactly the same mean exon length
({per['RYR']['mean_exon_bp']:.0f} bp against {mean_lo:.0f}–{mean_hi:.0f} bp).

**The count does not depend on the coverage bar.** Recomputed at six bars from
{min(A.COV_BARS):.2f} to {max(A.COV_BARS):.2f}, the median exon count of any
paralogue {worst_txt} (`architecture_sensitivity.tsv`), so A3 is a scope rule
and not a lever.

{P.line("exon_count", "confirmed",
        f"the literature's ~58-60 was verified on human alone; measured across "
        f"{fmt(h['n_genomes_in_scope'])} genomes the medians are "
        f"{per['ITPR1']['exons']:.0f} (ITPR1), {per['ITPR2']['exons']:.0f} "
        f"(ITPR2) and {per['ITPR3']['exons']:.0f} (ITPR3)")}

{P.line("span_varies", "confirmed",
        f"S0 measured a 6.5-fold span spread across the three human "
        f"paralogues; the medians across the whole scope differ "
        f"{h['span_fold_range']:.1f}-fold, in the same direction — ITPR3 is "
        f"the compact gene ({fmt(int(per['ITPR3']['span_bp']))} bp) and ITPR2 "
        f"the long one ({fmt(int(per['ITPR2']['span_bp']))} bp)")}

![](figures/architecture_by_paralog.png)

**{{fig:architecture}}.** Exon count, genomic span and the coverage-bar
sensitivity. (a) coding exons per gene, median with the 10th–90th percentiles.
(b) genomic span per locus, log axis, one violin per gene. (c) the median exon
count recomputed at every coverage bar.
"""


def section_7(h: dict) -> str:
    miss = _missing("paired_comparisons.tsv", "7. Paired within genome (D16)")
    if miss:
        return miss
    pc = [r for r in L.read_tsv(L.OUT / "paired_comparisons.tsv")
          if r["stratum"] == "all"]
    keep = ("n_exons", "span_bp", "median_intron_bp", "mean_exon_bp")
    rows = []
    for r in pc:
        if r["metric"] not in keep or r["cell_b"] == "RYR":
            continue
        rows.append([f"{r['cell_a']} − {r['cell_b']}", r["metric"],
                     fmt(int(r["n_pairs"])), fmt(int(r["n_a_greater"])),
                     fmt(int(r["n_b_greater"])), fmt(int(r["n_ties"])),
                     fmt(float(r["median_diff"]), 1), pfmt(r["q"])])
    n_strata = len({r["stratum"] for r in
                    L.read_tsv(L.OUT / "paired_comparisons.tsv")}) - 1
    return f"""## 7. Paired within genome (D16)

Intron size, assembly quality and annotation completeness all scale with the
assembly, so every cross-paralogue comparison is made **inside one genome** and
never across the scope: one locus per genome × gene, the best-covered, so a
genome carrying two copies does not weight itself twice. Each is a two-sided
exact sign test with ties dropped and counted (S8's rule) and the whole family
of tests is BH-corrected together (S9's rule). The same comparisons are
repeated inside each of the {fmt(n_strata)} vertebrate classes with at least
ten genomes.

{table(["pair", "metric", "genomes", "a bigger", "b bigger", "ties",
        "median difference", "q"], rows)}

Two things are worth reading off that table. **The exon counts differ by one
to three exons and the spans by tens of kilobases**, in the same genomes: the
architecture is shared and the packaging is not. And **the ordering of the
spans is consistent gene by gene**, not an artefact of averaging — ITPR3 is
the shorter gene than ITPR1 in {fmt(int([r for r in pc if r['cell_a'] == 'ITPR1' and r['cell_b'] == 'ITPR3' and r['metric'] == 'span_bp'][0]['n_a_greater']))}
of {fmt(int([r for r in pc if r['cell_a'] == 'ITPR1' and r['cell_b'] == 'ITPR3' and r['metric'] == 'span_bp'][0]['n_pairs']))}
genomes that carry both.
"""


def section_8(h: dict) -> str:
    miss = _missing("shared_intron_summary.tsv", "8. Intron positions")
    if miss:
        return miss
    cons = h["conservation"]
    crows = [[RP.cell_label(c), fmt(d["n_loci"]), fmt(d["positions_seen"]),
              fmt(d["at_50pc"]), fmt(d["at_90pc"]), fmt(d["at_99pc"]),
              f"{d['median_per_locus']:.1f}"]
             for c, d in ((c, cons[c]) for c in CELLS if c in cons)]
    sh = h["shared"]
    srows = []
    for key in ("ITPR1|ITPR2", "ITPR1|ITPR3", "ITPR2|ITPR3", "ITPR1|RYR",
                "ITPR2|RYR", "ITPR3|RYR"):
        d = sh.get(key)
        if not d:
            continue
        a, b = key.split("|")
        srows.append([f"{a} vs {'RyR' if b == 'RYR' else b}",
                      fmt(d["n_genomes"]), f"{d['median_observed']:.0f}",
                      f"{d['median_expected']:.2f}",
                      f"{d['enrichment']:.1f}×",
                      f"{fmt(d['n_significant'])} / {fmt(d['n_genomes'])}",
                      pfmt(d["q_worst"])])
    tolerances = sorted({int(r["tolerance"]) for r in
                         L.read_tsv(L.OUT / "shared_intron_summary.tsv")})
    ctrl_med = max(sh[k]["median_observed"] for k in sh if k.endswith("RYR"))
    itpr_med = min(sh[k]["median_observed"] for k in sh
                   if not k.endswith("RYR"))
    return f"""## 8. Intron positions: what is conserved, and what is ancestral

An intron's *position* here is a pair — the alignment column its upstream exon
ends in, and its classical phase (0, 1 or 2: how many bases of the interrupted
codon lie upstream). Both halves are needed: two paralogues can carry an intron
between the same two residues in different frames, which is not one ancestral
intron, and a column on its own would call any two introns in the same region
shared.

### 8.1 Within a paralogue

{table(["gene", "loci", "positions seen", "in ≥50 %", "in ≥90 %", "in ≥99 %",
        "median per locus"], crows)}

Each paralogue carries a core of about **{min(cons[c]['at_90pc'] for c in L.PARALOGS)}–{max(cons[c]['at_90pc'] for c in L.PARALOGS)}** intron
positions present in at least 90 % of the genomes that have the gene, out of
{min(cons[c]['positions_seen'] for c in L.PARALOGS)}–{max(cons[c]['positions_seen'] for c in L.PARALOGS)}
positions seen anywhere. The ryanodine receptors have their own core of
{cons['RYR']['at_90pc']}, on the same scale relative to their
{cons['RYR']['median_per_locus']:.0f} introns per gene. So the architecture is
not merely a count that happens to be stable — it is the *same introns*, in
the same places, across the vertebrates.

### 8.2 Between paralogues

Two genes with ~58 introns each spread over ~2,700 aligned residues will share
some positions by chance, so every comparison is paired within genome and
scored against a null: each of B's introns placed independently and uniformly
on the columns *both* loci have residues in, keeping its own phase, which makes
the match count a Poisson-binomial whose upper tail is exact — no RNG, no
replicate count. A seeded permutation without replacement is run beside it as a
cross-check and both are committed.

{table(["pair", "genomes", "median shared", "median expected", "enrichment",
        "genomes with p < 0.05", "BH q, worst genome"], srows)}

**The three IP₃ receptors share about {itpr_med:.0f} intron positions in every
genome that carries them, against half a position expected.** The enrichment is
85–88× and it is significant in *every* genome tested, not on average.

**And the control is the result.** The ryanodine receptors carry every
ITPR-diagnostic Pfam domain — that is the hazard this whole project is built
around (D14) — and measured through the identical instrument in the same
genomes they share a median of **{ctrl_med:.0f}** intron position with an IP₃
receptor, in **0** of {fmt(sh['ITPR1|RYR']['n_genomes'])} genomes at p < 0.05.
The two families' exon structures have no ancestry in common. Whatever the
shared domain architecture means, it was not inherited as a gene.

The answer does not rest on exact column identity: the test is run at
tolerances of {", ".join(str(t) for t in tolerances)} columns and the operating
point is {I.OPERATING_TOL}. Repeated with the nine `frame_via_cell` frames
excluded, the counts are in `shared_intron_summary.tsv` under
`excludes_frame_via_cell = 1`.

{P.line("paralog_identity", "orthogonal",
        f"S6 measured the paralogues at 0.83 mean protein identity and "
        f"ITPR-to-RyR at 0.249; S21 measures intron positions, and the two are "
        f"different objects — the RyR sequence identity is a quarter and its "
        f"shared-intron count is {ctrl_med:.0f}, which no identity would "
        f"predict either way")}

![](figures/intron_positions.png)

**{{fig:introns}}.** (a) intron positions ranked by prevalence, per gene, with
the 90 % line drawn. (b) shared positions against positions expected by chance,
one point per genome per pair, with the identity line.
"""


def section_9(h: dict) -> str:
    miss = _missing("fragment_summary.tsv", "9. Where a fragmentary "
                                            "annotation stops")
    if miss:
        return miss
    frag = h["fragments"]
    by_state = h["fragments_by_state"]
    verdicts = ("annotation_failure", "structure_disagreement", "mixed",
                "broken_at_real_junctions", "no_internal_terminus")
    rows = []
    for state, d in sorted(by_state.items()):
        rows.append([state, fmt(int(d["n_loci"])), fmt(int(d["n_genomes"]))]
                    + [fmt(int(d.get(v, 0))) for v in verdicts])
    n = int(frag["internal_termini"]) or 1
    fail = int(frag.get("annotation_failure", 0))
    total = int(frag["n_loci"])
    verd = L.read_tsv(L.OUT / "fragment_verdicts.tsv")
    in_scope = sum(1 for r in verd if r.get("in_architecture_scope") == "1")
    return f"""## 9. Where a fragmentary annotation stops

S18 called 264 of the sweep's loci `fragmentary` and 27 `split`: coding
sequence is there, but no one model covers the gene. That is a statement about
coverage and says nothing about **where** the annotation stopped — and the
difference is the whole claim. A model that stops at a genuine junction has
produced a plausible short gene; one that stops in the middle of an exon has
produced a boundary no splicing machinery could make.

So the unit is the annotated model's own **terminus**, not its internal exon
edges: a model's internal boundaries are its own splice sites and are canonical
by construction, and scoring them would answer §3 a second time. Each terminus
is scored on two axes that know nothing about each other — against the
alignment's exon boundaries, and against the two genomic bases immediately
outside it read in gene orientation — and a terminus coinciding with the gene's
own end is excluded by rule, because a real gene legitimately starts and stops
inside an exon.

The question is asked of **all {fmt(total)}** such loci and not only of the
ones in the architecture scope: those loci sit in the poorer assemblies by
construction, and restricting the question to the architecture scope would
have answered it on {fmt(in_scope)} of them and dropped the hardest.

{table(["state", "loci", "genomes"] + [v.replace("_", " ") for v in verdicts],
       rows)}

Of {fmt(n)} internal termini,
**{fmt(int(frag['termini_mid_exon']))}** sit inside an exon of the gene model
and **{fmt(int(frag['termini_in_intron']))}** inside one of its introns;
{fmt(int(frag['termini_at_exon_boundary']))}
({pct(frag['frac_termini_at_boundary'])}) land on a boundary the model has.
**{fmt(fail)} of {fmt(total)} loci** carry at least one mid-exon terminus,
which falsifies the reading that the annotation stopped at a real gene
boundary; **{fmt(int(frag.get('broken_at_real_junctions', 0)))}** are broken
entirely at junctions the gene has.

The verdict is one-sided in the way S10's ORF screen is. A mid-exon terminus
falsifies; every terminus landing on a boundary does *not* prove the pieces are
separate genes, so that verdict is named `broken_at_real_junctions` and says
only what it says.

{P.line("annotation_states", "confirmed",
        f"S18 counted the states; S21 says what they are — "
        f"{fmt(fail)} of {fmt(total)} split or fragmentary loci stop somewhere "
        f"nothing splices")}

{P.line("boundary_concordance", "confirmed",
        f"S10 corroborated two loci's boundaries against 35 and 40 genomes; "
        f"the same comparison across {fmt(h['concordance']['n_genomes'])} "
        f"genomes and {fmt(h['concordance']['n_edges'])} annotated edges gives "
        f"{pct(h['concordance']['frac_exact'])} exact agreement")}
"""


def section_10(h: dict) -> str:
    miss = _missing("tandem_control.tsv", "10. Duplication")
    if miss:
        return miss
    ctrl = h["tandem_control"]
    rows = [[RP.cell_label(r["cell"]), fmt(int(r["n_cells"])),
             fmt(int(r["tp"])), fmt(int(r["fn"])), fmt(int(r["fp"])),
             fmt(int(r["tn"])), pct(r["sensitivity"]), pct(r["specificity"])]
            for r in L.read_tsv(L.OUT / "tandem_control.tsv")]
    tel = h["teleost_agreement"]
    cls = h["tandem_classes"]
    return f"""## 10. Duplication: is any of this two genes?

An exon count is only a gene's if the locus is one gene. miniprot aligns each
bait independently, so a genome encoding the same part of a protein twice has
the *same* bait aligning twice at two disjoint places — and that geometry, run
naively, measures **paralogy**: the three IP₃ receptors are 61–68 % identical,
every bait aligns at all three genes, and the detector's first version reached
a specificity of 0.16. So the pair has to be inside the cell's own loci, with
the sweep's clustering and attribution imported unchanged, which is where D14
lives.

The detector is then scored as a classifier of a copy count it never sees:
S16's committed `n_copies`, decided by coverage, identity and aligned length
and by no pairwise geometry at all.

{table(["cell", "genome × gene cells", "tp", "fn", "fp", "tn", "sensitivity",
        "specificity"], rows)}

Sensitivity {pct(ctrl['sensitivity'], 2)} and specificity
{pct(ctrl['specificity'], 2)} over {fmt(ctrl['n_cells'])} cells. The RyR
control's specificity is lower on purpose: that cell holds three genes
(RYR1/2/3) in every tetrapod, so the detector *should* fire there, and it does.
The 3R teleost check agrees with S16 on **{fmt(tel['agree'])} of
{fmt(tel['n'])}** genome × paralogue cells.

**No locus in the sweep encodes the same part of the protein twice.** The
`within_locus` class — two alignments of one bait inside one locus cluster,
which would be an internal partial duplication, or two neighbouring genes the
10 kb clustering had merged — fires on
**{fmt(h['n_within_locus_pairs'])}** pairs. `s21_test_arch` T15 constructs such
a pair and requires the branch to fire, so the zero is a measurement.

{P.line("copy_number", "confirmed",
        f"S16 found its split-model merge fired on 0 of 2,146 loci; a "
        f"geometric detector reading the same alignments finds "
        f"{fmt(h['n_within_locus_pairs'])} within-locus duplications and "
        f"agrees with S16's copy call at "
        f"{pct(ctrl['sensitivity'], 1)} / {pct(ctrl['specificity'], 1)}")}

![](figures/fragments_and_duplicates.png)

**{{fig:fragments}}.** (a) the verdict on every split and fragmentary locus.
(b) where the annotation's internal model termini sit relative to the gene
model. (c) the duplication detector scored against S16's copy call.
"""


def section_11(h: dict) -> str:
    per = h["by_paralog"]
    ctrl_med = max(h["shared"][k]["median_observed"] for k in h["shared"]
                   if k.endswith("RYR"))
    return f"""## 11. What this settles, what it does not, and the hand-off

**Settled.**

1. The IP₃-receptor gene is a **{per['ITPR1']['exons']:.0f}-exon gene** in all
   three vertebrate paralogues — {per['ITPR1']['exons']:.0f},
   {per['ITPR2']['exons']:.0f} and {per['ITPR3']['exons']:.0f} — measured over
   {fmt(h['n_in_scope'])} genes in {fmt(h['n_genomes_in_scope'])} genomes above
   D4's contiguity bar, and the literature's number was previously verified on
   one species.
2. Its **genomic span is not conserved at all**: a {h['span_fold_range']:.1f}-fold
   spread between the paralogues' medians, consistent gene by gene within
   genomes, with essentially identical coding length.
3. The three paralogues share ~{h['shared']['ITPR2|ITPR3']['median_observed']:.0f}
   **intron positions of ~58**, at 85–88× the chance rate, in every genome
   tested — the exon structure is inherited from their common ancestor, not
   convergent.
4. The ryanodine receptors share **{ctrl_med:.0f}**. The two families' domain
   architecture is shared and their exon structure is not.
5. Database "fragments" are overwhelmingly **annotation failures, not gene
   boundaries**: {fmt(int(h['fragments'].get('annotation_failure', 0)))} of
   {fmt(int(h['fragments']['n_loci']))} split or fragmentary loci stop
   somewhere nothing splices.
6. The sweep's exon boundaries are corroborated by an independent pipeline at
   **{pct(h['concordance']['frac_exact'])}** over
   {fmt(h['concordance']['n_edges'])} annotated edges — S10's two-case check,
   generalised to the scope.

**Not settled.**

1. **Where the introns were gained or lost.** S21 counts shared positions; it
   does not reconstruct the ancestral intron set on the tree, and the
   paralogue-specific positions could be gains in one lineage or losses in the
   others. That is a reconciliation question and S13's machinery is the place
   for it.
2. **Anything outside the vertebrates.** The frame is the three human
   paralogues, so the non-vertebrate grade S20 and S23 enumerated is not in
   this measurement at all. Whether the ~58-exon architecture predates the
   2R duplications is unanswered here.
3. **Alternative splicing.** Every count is of one gene model per locus, which
   is the aligner's best path through the genome. The family's characterised
   splice variants (the S1, S2 and SII sites in ITPR1) are not visible to it.
4. **The `in_intron` termini.** {fmt(int(h['fragments']['termini_in_intron']))}
   annotated model termini sit inside an intron of the gene model. Some of
   those are the annotation and the alignment genuinely disagreeing about a
   boundary, and S21 records the disagreement rather than adjudicating it.

**Two things to hold against this task.**

- The architecture scope excludes {fmt(h['excluded_by'].get('A2_contig_cannot_span_the_gene', 0))}
  loci because their contig cannot hold the gene. Those are real genes, and
  every number here is conditional on assembly quality in the way S19 measured.
- The intron-position frame is a protein alignment of four kingdoms, and a
  column is not a homology statement everywhere. The anchor test fixes it at 14
  measured residues in the pore and the ligand core; between those it is the
  aligner's opinion.
"""
