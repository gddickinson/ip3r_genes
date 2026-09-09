"""S19's report, part three: the limits (§7-§10).

The half that argues against the task's own instrument — where iterative
search drifts and the rule written to catch it does not fire, where the
inference methods ran out, what the whole thing changes for anyone doing the
same work, and what it does not settle. Split from `s19_report_results` to
keep both inside the 500-line budget (the `s3_report.py` /
`s3_report_d10.py` pattern, applied twice as S7, S9, S12 and S15b do).

Takes the caller's headline dict, so the three halves cannot read the tables
differently.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_priors as P                                          # noqa: E402
from s19_report import fmt, pct, table                          # noqa: E402
from s19_report_results import _missing, _t                     # noqa: E402


# --------------------------------------------------------------- 7. the drift

def section_7(h: dict) -> str:
    val = _t("kill_criterion_validation.tsv")
    if not val:
        return "## 7. Iterative-search drift\n\n" + _missing(
            "kill_criterion_validation.tsv")
    trace = _t("kill_criterion_trace.tsv")
    outcome = _t("drift_outcome.tsv")
    runs = _t("jackhmmer_runs.tsv")
    yield_rows = _t("iteration_yield.tsv")
    seed_rows = _t("seed_effect.tsv")
    vert = next(r for r in yield_rows if r["database"] == "vertebrata")
    prot = next(r for r in yield_rows if r["database"] == "protista_other")
    drifted = [r["run"] for r in outcome if r["drifted"] == "1"]
    return f"""## 7. The rule written to catch iterative drift never fires

Seven jackhmmer runs, {fmt(h['drift']['drifted'])} of which walked out of the
family. The outcome is measured on the **finished model** rather than taken
from any session's prose: the share of a run's final included set that
neither profile scores at all. The seven separate cleanly —
{", ".join(f"{r['run']} {float(r['offfamily_frac_final']):.2f}" for r in outcome
           if r["drifted"] == "1")}
against
{", ".join(f"{r['run']} {float(r['offfamily_frac_final']):.2f}" for r in outcome
           if r["drifted"] == "0")}
— with nothing between {max(float(r["offfamily_frac_final"]) for r in outcome
                             if r["drifted"] == "0"):.2f} and
{min(float(r["offfamily_frac_final"]) for r in outcome
     if r["drifted"] == "1"):.2f}, so the labels are not a judgement call.

Scoring each of D10's rules as a classifier of that outcome:

{table(["rule", "fires on a drifted run", "fires on one that did not",
        "sensitivity", "specificity"],
       [[r["rule"], f"{r['true_positive']}/{h['drift']['drifted']}",
         f"{r['false_positive']}/{int(r['n_runs']) - h['drift']['drifted']}",
         pct(float(r["sensitivity"])), pct(float(r["specificity"]))]
        for r in val])}

**K1 — the rule written for exactly this hazard — fires on none of the seven
runs, including all three that drifted.** The reason is D10b, and it is now
measured rather than described: off-family accretion *dilutes* the
sister-family share, so K1's own statistic moves the wrong way while a run
drifts. The sister share **fell** over the run in
{fmt(h['drift']['sister_fell'])} of {fmt(h['drift']['n_runs'])} runs.

K3, the round ceiling, catches all three but also flags three runs that did
not drift — it counts rounds rather than content, which is right for a budget
and useless as a diagnosis. That is why S20b had to separate `protista_other`
from `metazoa_nonvert` in prose: both hit K3, one had accreted 22,913
off-family targets and the other had simply not finished.

### 7.1 One change fixes it, and it is offered as a proposal

Moving K1's own threshold — the same 0.10 — from the sister-family share to
the **off-family** share separates the seven runs completely, firing at round
2 or 3 in each of {", ".join(drifted)} and never in the other four. It is
reported here and **not applied**: it is validated after the fact on seven
runs with an inherited threshold, and applying it would overturn committed
verdicts, which is a decision for a session that owns those tables.

{P.line("d10b", "confirmed",
        f"K1 fires on 0 of 7 runs and the sister share falls in "
        f"{h['drift']['sister_fell']} of them; the off-family axis "
        f"separates the same seven runs at sensitivity 1.00 and "
        f"specificity 1.00")}

### 7.2 What iteration bought

{table(["run", "database", "one-pass family calls", "accepted targets",
        "final targets", "verdict", "rule"],
       [[r["run"], r["database"], fmt(int(r["hmmsearch_called_family"])),
         fmt(int(r["accepted_targets"])), fmt(int(r["final_targets"])),
         r["verdict"], r["kill_rule"]] for r in runs])}

What those targets *are* is the question a completeness argument turns on:

{table(["database", "one-pass family calls", "accepted targets across runs",
        "of which called family", "of which not"],
       [[r["database"], fmt(int(r["hmmsearch_called_family"])),
         fmt(int(r["accepted_targets_union"])) if r["accepted_targets_union"]
         else "—",
         fmt(int(r["of_which_called_family"]))
         if r["of_which_called_family"] else "—",
         fmt(int(r["of_which_not_called_family"]))
         if r["of_which_not_called_family"] else "—"]
        for r in yield_rows if r["accepted_targets_union"]])}

In the vertebrate database, three iterated runs between them hold
{fmt(int(vert['of_which_called_family']))} of the
{fmt(int(vert['hmmsearch_called_family']))} records one hmmsearch pass calls
family, and {fmt(int(vert['of_which_not_called_family']))} that it does not.
So iteration returned **no record the profile pair calls family that one
pass had not already returned**, and twenty thousand that it does not call
family at all. The same shape holds in the protists
({fmt(int(prot['of_which_called_family']))} family against
{fmt(int(prot['of_which_not_called_family']))} not). Only the run that
converged cleanly stays close to its one-pass set. Whether any of those
non-family targets is a real family member the profiles failed to score is
not answerable from this table, and S3 and S20 both examined them: they are
module-only matches to the shared domains.

### 7.3 Seed choice does not change the family, only the debris

The three vertebrate runs were seeded from a human paralog, a fly gene and an
amoebozoan gene — as unlike each other as this family allows.

{table(["run", "seed", "accepted targets", "of which family",
        "family unique to this seed", "family calls this run lacks"],
       [[r["run"], r["seed"], fmt(int(r["accepted_targets"]))
         if r["accepted_targets"] else "—",
         fmt(int(r["of_which_family"])),
         fmt(int(r["family_unique_to_this_seed"]))
         if r["family_unique_to_this_seed"] != "" else "—",
         fmt(int(r["hmmsearch_family_not_in_this_run"]))]
        for r in seed_rows])}

The family content is the same set three times over: the three runs
intersect on {fmt(int(seed_rows[-1]["of_which_family"]))} family records and
differ by at most {fmt(max(int(r["family_unique_to_this_seed"]) for r in
                           seed_rows[:-1]
                           if r["family_unique_to_this_seed"] != ""))}.
What the seed changes is everything that is *not* the family: the amoebozoan
run carries {fmt(int(seed_rows[2]["of_which_not_family"]))} non-family
targets against
{fmt(int(seed_rows[0]["of_which_not_family"]))} and
{fmt(int(seed_rows[1]["of_which_not_family"]))} for the other two. A distant
seed does not reach further into the family; it reaches further out of it.
"""


# ----------------------------------------------------------- 8. inference

def section_8(h: dict) -> str:
    if not _t("codon_model_power.tsv"):
        return "## 8. Where the inference ran out\n\n" + _missing(
            "codon_model_power.tsv")
    inf = h["inference"]
    syn = _t("synteny_power.tsv")
    reach = next((r for r in syn if r["measure"] == "reach"), {})
    return f"""## 8. Where the inference methods ran out, not the search

Four limits this project met while doing something else. Each is a general
result about the method.

**A reconciliation's loss count is mostly sampling.** Losses implied by a
gene tree over a 134-tip representative alignment, checked cell by cell
against the 309-genome ledger: {fmt(inf['reconciliation']['implied'])} implied
cells, {fmt(inf['reconciliation']['artefact'])} of them the paralog present in
the genome and absent only from the sample, and
{fmt(inf['reconciliation']['corroborated'])} corroborated.

{P.line("recon_losses", "confirmed",
        f"{inf['reconciliation']['artefact']} of "
        f"{inf['reconciliation']['implied']} implied-loss cells are "
        f"sampling; the recomputation agrees with S13's own aggregate row")}

**Synteny disambiguation is accurate and unavailable.** S8's caller is right
wherever it acts and almost never acts on the cells that need it: of
{fmt(int(reach.get('n', 0)))} trace regions it reaches
{fmt(int(reach.get('n_positive', 0)))}
({pct(float(reach.get('frac', 0)))}), because a fragment's contig carries no
neighbours to read. Accuracy and reach are different numbers and quoting only
the first would report contig lengths as method performance.

{P.line("synteny_reach", "confirmed",
        f"{int(reach.get('n_positive', 0))} of {int(reach.get('n', 0))} "
        f"regions reach the caller's floor, at 100 % accuracy in every key "
        f"bin where it acts")}

**Per-site codon models are past their power at these divergences.** The
synonymous rate is unidentifiable at
{fmt(inf['fel']['n_unidentifiable'])} of {fmt(inf['fel']['n_sites'])} sites
({pct(inf['fel']['frac'])}), and
{fmt(inf['saturation']['n_saturated'])} of
{fmt(inf['saturation']['n_pairs'])} pairwise comparisons
({pct(inf['saturation']['frac'])}) are flagged saturated. Any per-site or
per-element omega formed as a ratio of sums is dominated by sites whose
denominator is not estimable.

{P.line("fel_power", "orthogonal",
        f"S17 already took omega only over identifiable sites; this "
        f"recomputation states the size of what that excludes "
        f"({pct(inf['fel']['frac'])} of scored sites) rather than "
        f"disagreeing with it")}

**A likelihood tree need not resolve the question asked of it.** The AU test
leaves {fmt(inf['tree_confidence_set'])} topologies in the 95 % confidence
set, so the sister arrangement is reported as what the data supports rather
than as the ML tree's answer.
"""


# ------------------------------------------------------------ 9. what it means

def section_9(h: dict) -> str:
    return f"""## 9. What this changes for anyone doing the same thing

Five results, in the order they would change a protocol.

1. **A genome-scale ortholog sweep misses about one cell in seven, and every
   miss is an assembly.** {pct(h['itpr_present']['rate'])} of control cells
   over the whole scope, {pct(h['itpr_present']['rate_at_d4_bar'])} above a
   contiguity bar set at the gene's own median span. A survey that reports
   absences without a contiguity floor is reporting assembly quality.
2. **A contiguity bar can be set from gene geometry before any error is
   measured.** D4's was, and it lands where the calibration would have put
   it. The bar is not free: it costs
   {pct(1 - h['itpr_present']['genomes_at_d4_bar'] / h['n_genomes'])} of the
   genomes, disproportionately in the clades a loss survey is most interested
   in.
3. **Phylogenetic breadth in a protein bait panel is nearly worthless within
   the vertebrates; paralog coverage is everything.** Four human baits match
   thirty-eight. One bait recovers its gene at any identity above 0.5. Spend
   the panel budget on paralogs and on clades where no labelled record
   exists, not on sampling depth within a paralog.
4. **A family profile HMM over reference proteomes is not a discovery method
   for whole genes.** At gene scale it returns what domain annotation already
   returns. Its gain is fragments — and, in the least-annotated clades,
   gene-scale records the annotation never carried.
5. **The rule everyone writes to catch iterative-search drift measures the
   wrong thing.** Sister-family contamination is diluted by the drift it is
   meant to detect. The off-family share is what moves.

And the one that is about the archive rather than the method:
**{pct(h['gene_recovery']['invisible'] / h['gene_recovery']['n_genes'])} of
the genes this project demonstrated are not reachable from any protein
database.** Not because they are hard to find — because no protein record of
them exists.
"""


# ------------------------------------------------------------ 10. caveats

def section_10(h: dict) -> str:
    return """## 10. What this does not settle

- **The control rests on S15b's zero.** If a loss were later established
  anywhere in the scope, that cell would move from the numerator of the
  false-negative rate to a real absence. The RyR series is carried precisely
  so the reader can see how much of the number depends on that: it is
  measured independently of S15a's states and agrees to within two points.
- **`present_fragmented` and `present_partial` states were reached using
  alignment evidence.** They are outside the loci the sweep placed and
  cleared a calibrated bar, so they are not the sweep's own answer restated
  — but they are not an independent instrument either. The stricter reading
  is the RyR series, which does not use them at all.
- **The panel ablation is exact for locus placement and silent about
  rescue.** Removing baits changes what tblastn rescue attributes, and that
  path needs the genome FASTA, which the sweep deleted after searching. So
  the ablation measures what the panel was worth to miniprot, not to the
  whole pipeline.
- **The single-bait design rule is conditional on the gene being there.** It
  is measured at loci the full panel placed. It bounds what one bait can do
  at a known locus; it does not license a one-bait panel in an unsearched
  clade.
- **The drift criterion is validated after the fact on seven runs**, with a
  threshold inherited from K1 rather than fitted. Seven runs cannot establish
  a specificity of 1.00. It is a proposal with its evidence attached, and it
  is deliberately not applied to any committed verdict.
- **The head-to-head is a comparison of two channels' *calls*, not of two
  algorithms in the abstract.** The profile pair's gate (D22) is part of the
  method being scored, and a different gate would move the tail columns.
"""
