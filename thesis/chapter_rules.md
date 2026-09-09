# Chapter rules

*Written before any chapter prose, and committed with the assignment table it
produced (`chapter_assignment.tsv`). S26 carves the same body of results into
papers; whichever task ran first was to commit its assignment and the other to
adopt it or record every departure. This one ran first.*

The manuscript is 17,807 words and cites 29 references. It is the compressed
form of 109,243 words of committed task reports, 75 recorded methodological
decisions (D0-D71, four of them lettered sub-decisions) and a 137-reference literature review. Compression is not the
thesis's problem — the paper does that well. The thesis's problem is the
opposite one: with the page budget removed, what stops the document becoming
the same paper with padding between the paragraphs?

These seven rules are the answer. They were written to be *testable*, and
four of them are enforced by `scripts/s25_assign.py` and the build rather than
asserted here.

**T1 — one thread of argument.** A chapter follows one line of reasoning from
a question to an answer. It may take several questions where they are steps in
one argument; it may not take two arguments because they happened in the same
session. This is deliberately weaker than S26's P1 ("one question, statable in
a sentence without an *and*"), because a chapter is allowed to be a stage in a
longer construction and a paper is not.

**T2 — the instrument is explained where it is first used.** A results
directory is primary in the chapter that first depends on its measurement, and
that chapter carries the reasoning behind it: what the threshold is, what it
was measured against, and what it would have cost to guess. This is the rule
that separates this document from the paper. A paper cites an instrument; a
chapter builds it.

**T3 — one sitting.** Four to twelve figures, and a chapter a reader can
finish. A grouping that yields two figures is a section of its neighbour; one
that needs twenty is two chapters. This is what split the channel in two: the
constraint map and the ligand-site comparison are one subject and sixteen
thousand words of report, so Chapter 11 takes the receptor as a whole and
Chapter 12 takes the one part the ryanodine receptors do not share.

The band is deliberately wider than S26's P4, which allows a paper four to
seven main figures. That is the difference between the formats stated as a
number: a paper's figure budget is set by a journal and a chapter's by a
reader's attention, and a chapter is allowed to be longer. The introduction,
the methods chapter and the general discussion are exempt — they are
exposition and carry whatever they need, which for the methods chapter is
one figure.

**T4 — a results directory is primary in exactly one chapter.** Enforced.
Every other chapter cites it. The rule exists so that the thesis and the paper
series can be compared: S26's P5 is the same rule, so a directory that is
primary in Chapter 7 here and in two papers there is a visible disagreement
rather than an invisible one.

**T5 — a measured dead end belongs with the instrument it was measured on.**
There is no chapter of failures. ModelFinder's exhaustive scan, abandoned at
11 of up to 1,232 models and projected at 21 hours, is in the chapter that
chooses the substitution model. The synteny caller threshold that optimising
call rate alone would have selected is in the chapter that calibrates the
caller. The duplication detector that reached a specificity of 0.16 before it
was scoped to the cell's own loci is in the chapter that measures gene
structure. A collected chapter of dead ends would read as an apology; each one
in its place reads as the reason a number is trustworthy.

**T6 — a grouping that only holds leftovers is an appendix.** Enforced, and
in the strong direction: every entry under `results/` is either assigned to a
chapter with a rule and a reason, or named in an exclusion list with why it is
not a result. A results directory that exists and is not assigned fails the
build. One entry is excluded — `results/session_live.json`, the dashboard's
progress panel, which is session state and not a measurement.

**T7 — every figure is copied from a committed results directory, never
re-plotted.** D13 and D19 applied to the thesis. The thesis may carry figures
the paper had no room for; it may not carry a figure that is not already
committed beside the data it was drawn from. The build fails on a figure whose
png or pdf is missing, on a placed figure that is not declared, and on a
declared figure that is never placed.

## What the rules produced

Fifteen chapters, five appendices, and the departure below.

| # | Chapter | Primary for |
|---|---------|-------------|
| 1 | The receptor, and the question | — (exposition) |
| 2 | The baseline, and what it could be trusted to say | `s0_baseline`, the four S0 smoke bundles |
| 3 | Two families, one architecture | `benchmark_controls`, `census_v2`, `census_v3`, `hmm_sweep` |
| 4 | The vertebrate genomic sweep, and what the search was worth | `genome_manifest.tsv`, `s5_baits`, `genome_ledger`, `census_v4`, `methods` |
| 5 | How far the family reaches | `s20_sweep`, `s23_baits`, `s23_scope`, `census_v5`, `census_v6` |
| 6 | The alignment and the tree | `msa_v2`, `phylogeny` |
| 7 | Where ITPR1, ITPR2 and ITPR3 came from | `synteny`, `reconciliation`, `duplication` |
| 8 | The gene: exons, introns, and what a fragment is | `gene_architecture` |
| 9 | Counting a loss that never happened | `loss_dynamics`, `loss_counts` |
| 10 | Selection across the vertebrate family | `selection` |
| 11 | The channel: shape, constraint and the clinic | `structures`, `constraint` |
| 12 | The part the ryanodine receptors do not share | `ligand_site` |
| 13 | An archive that cannot find the gene | `annotation_bugs`, `expression`, `annotation_audit` |
| 14 | Methods, and the reasoning behind them | `supplementary`, `toolchain_manifest.txt` |
| 15 | General discussion | — (exposition) |

## The one departure to expect from S26

S26's starting proposal splits S19 (`results/methods`) across two papers: the
retention paper takes its false-negative rate, the archive paper its
contribution half. T4 forbids that here, so the whole of S19 sits in Chapter 4
with the search it measures, and Chapters 9 and 13 cite it.

The disagreement is real and worth stating in both directions. S19's
false-negative rate is a property of *the search* — 15.0 % of cells whose gene
is independently known present were not found by the ledger — so a thesis that
explained the search in Chapter 4 and then withheld the number that says what
that search was worth until Chapter 9 would be hiding the instrument's error
bar from the chapter that builds the instrument. A paper is under the opposite
pressure: a retention claim without its own false-negative rate is a paper
whose central control lives somewhere else, which is exactly what S26's P2
forbids. Both are right about their own format. Appendix E states it as a
disagreement rather than resolving it, and S26 is free to keep its split as
long as it records this table as the reason it had to choose.
