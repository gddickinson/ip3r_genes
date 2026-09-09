# 10. Measuring selection across the vertebrate family

## 10.1 An identity is a distance, and this chapter measures a rate

Chapter 6 measured identity within a paralogue at 0.910 and concluded the
family is conserved. That is a distance. A conserved protein and a slowly
evolving one are not the same statement: identity says how far two sequences
have moved apart, and a selection test says whether the moves that did happen
were the ones selection tolerates.

This chapter measures the rate. It is the shortest analytical chapter in the
thesis and the one with the most traps in it, because a codon model returns a
perfectly plausible number when its input is wrong, when its optimiser has
found a local peak, or when its denominator has stopped existing.

## 10.2 Every coding sequence had to be proved to encode the aligned protein

A ratio of non-synonymous to synonymous substitution is a statement about
codons, so every number here rests on a nucleotide sequence that provably
encodes the exact protein the alignment holds. Nothing is a translated
database record taken on trust.

Three routes supply candidates, comprising a genome database, cross-references
from the protein record, and the genomic gene models the sweep produced, and
each is a **generator of candidates** rather than an answer. A protein entry
lists every transcript of its gene and the entry's own sequence is one
particular isoform. Human ITPR1 lists five and ITPR2 two, and in **both** the
first is not the isoform the alignment holds. Returning the first coding
sequence that downloads would silently substitute a different isoform for the
protein the tree was built on, which is a wrong codon alignment rather than a
missing one.

Every candidate is translated against the aligned protein and **every**
disagreement is masked, not only internal stops. A codon that encodes
something else is not that protein's codon at that site whatever the cause,
and the model must see missing data rather than a residue that is not there.

Fifty-seven of 57 vertebrate tips have a validated coding sequence, from two
routes, with 24 codons masked across the whole set.

**The genomic route has its own trap.** The sweep ran the aligner in a mode
whose emitted protein skips frameshift residues, so the alignment's coding
blocks run out of frame against it and splicing them is not offered at all.
The locus is realigned in a mode that keeps those residues, and the
reconstruction is accepted only if it reproduces the sweep's protein. The one
concession is measured rather than assumed: a locus re-run indexes a fraction
of the genome the sweep indexed, so a terminal exon can be placed differently,
by nine residues of 2,666 in one frog. Same length and within 1 % is the same
gene model and its disagreeing codons are masked, and anything beyond that is
a different model and is refused.

The codon alignment is then produced by a protein-guided aligner [R148] and
**cross-checked nucleotide by nucleotide against an independent in-house
mapping**, with a disagreement aborting the build. A residue masked in step
one is written as an unknown in the protein rows too, so the aligner is never
asked to reconcile a pair that disagrees, which it does by dropping the
sequence silently. Trimming is chosen on the protein and applied codon-aware,
whole triplets only.

![](figures/supp_codon_alignment.png)

**{fig:supp_codon_alignment}.** The trimmed codon alignment behind every
estimate in this chapter, plotted in codons, because an axis in nucleotides
would make a three-fold difference look like a property of the data.

**Twelve constructed negative controls run on every build**, and their
character is the point. A codon alignment is the one artefact in this project
where a silent error is invisible downstream, because the model will fit a
frame-shifted alignment and return a perfectly plausible rate. So the controls
are checks on refusal: a frame-shifted sequence, one encoding a different
protein, an internal stop that must be masked and not refused, whole-triplet
gap stripping, and a non-monophyletic foreground being refused.

## 10.3 The tree decides which tips are in which selection set

A branch model marks a node. A foreground that is not a clade does not fail.
It silently marks a larger one and returns a well-formed rate for a hypothesis
nobody asked.

So the three selection sets are Chapter 6's extended paralogue clades,
re-derived here from the rooted tree with that chapter's own rule and
cross-checked against its committed table, with a disagreement in size or
membership treated as a hard failure.

Two consequences follow. A tip the tree recorded as unplaced is placed by the
tree or by nothing, never by a database symbol. And the seven unlabelled tips
the tree nests inside a paralogue clade join that set, because outside the
vertebrates a paralogue label means nothing, but inside a maximally supported
vertebrate paralogue clade the tree has just supplied one.

Six tips sit in no paralogue clade. They stay in the whole-tree analyses,
because dropping them would change the branch lengths every other estimate is
made on, and they are in no foreground, because a locus the tree could not
place is not evidence about a paralogue.

They are the six cyclostome loci. Chapter 6 placed them in two cyclostome-only
clades and handed the question on, Chapter 7 reported the neighbourhood
underpowered because after 550 million years there is no shared flank
vocabulary left, and a third instrument now declines the same question for a
third reason: a codon model can only ask about a paralogue whose clade the
tree defines, and for these six it defines none.

## 10.4 All three paralogues are under strong purifying selection

One-ratio estimates [R145] give ITPR1 at 0.0238, ITPR2 at 0.0430 and ITPR3 at
0.0415. The highest of the three is **23 times below neutrality**.

![](figures/omega_by_paralog.png)

**{fig:omega_by_paralog}.** Per-paralogue rate with the curated-sequence
sensitivity estimate beside it, and neutrality drawn rather than described.
The axis is logarithmic, because at a rate of 0.03 a linear axis puts every
bar on the floor and hides the one thing a reader wants, which is how far
below one it sits.

The identity Chapter 6 measured says the same thing far less sharply. This is
a 2,700-residue channel accumulating one non-synonymous change per 23
synonymous ones.

The sensitivity subsets, meaning the same paralogues with every genomic gene
model removed so that no masked codon contributes, move the estimate by at
most 0.0112. The estimates are not an artefact of the reconstructed models.

## 10.5 Synonymous saturation qualifies every other number in this chapter

Synonymous sites are **saturated across the vertebrate span, inside a single
paralogue and not only between them.** In the worst set, 94.2 % of
within-paralogue pairs exceed the conventional saturation bar, and the medians
run from 4.6 to 13.5.

![](figures/dnds_saturation.png)

**{fig:dnds_saturation}.** Pairwise rates within each paralogue, drawn log-log
with the neutral diagonal and the saturation bar. On a logarithmic-x, linear-y
plot the neutral diagonal is not a line at all: the first draft's neutrality
ran off the panel within the first pixel and left the saturation bar as the
only line on the figure, which reads as neutrality.

The expectation going in was saturation between the paralogues, which are
older than 500 million years. It is already reached within them, because a
single paralogue set spans shark to teleost to mammal.

That has one immediate consequence and one that runs through the rest of the
chapter. **The pairwise rate matrix is a diagnostic here rather than an
estimate**, so every rate quoted comes from a tree-based model, which
distributes substitutions over branches instead of asking one pair to carry
450 million years. And every rate should be read as a lower bound on precision
rather than a point estimate with a small error, because a ratio whose
denominator is soft is soft.

## 10.6 ITPR1 is the most constrained of the three

Ranked from most constrained to least, the order is ITPR1, then ITPR3, then
ITPR2.

Each paralogue clade tested against the rest of the family as a two-ratio
branch model confirms it, with every test corrected across the family. ITPR1's
foreground rate is 0.0241 against a background of 0.0432, ITPR2's is 0.0435
against 0.0317, and ITPR3's is 0.0455 against 0.0308.

![](figures/branch_contrast.png)

**{fig:branch_contrast}.** Background against foreground rate per paralogue
clade, with the relaxation coefficient beside it.

ITPR1 is the paralogue that carries the family's dominant missense disease
burden, and it is the most constrained. It is also the paralogue whose
neighbourhood Chapter 7 found best preserved and the one that alone kept its
teleost duplicate.

**ITPR3's neighbourhood is the one that does not travel, while its coding
sequence is not the least constrained.** Those are different quantities,
because neighbourhood conservation is rearrangement history and this is coding
sequence rate, so the mismatch is not a disagreement, and the project's
verdict vocabulary has a word for that rather than forcing a choice.

## 10.7 The stem branches, and three optimiser failures that a single run would have hidden

A duplicate's fate is decided on its stem, meaning the interval between the
duplication and the first surviving split of the new copy. If a paralogue was
ever free to change, that is when. Branch-site model A [R146] asks whether a
class of sites on that one branch has a rate above one while the rest of the
tree does not.

**Model A is restarted from four initial rates by construction here** rather
than repaired afterwards. A nested alternative cannot have a lower optimum
than its own null, and yet **3 of the twelve restarts converged below their
own null, one on every stem.** Each of those is a local-optimum failure that a
single-start run would have reported as its answer.

A different starting value fails on each stem, so no single initial value
would have been safe. That is the case for running several rather than for
choosing a better one.

![](figures/bs_restarts.png)

**{fig:bs_restarts}.** Every restart against its own null. A point below the
line is a local optimum rather than a result.

**All three stems are significant after correction, and one of the three
carries a rate the data actually determine.** That distinction is the
substance of this section, and it comes from a parameter printed outside the
test rather than from the test.

On the ITPR2 and ITPR3 stems the foreground rate is pinned at the optimiser's
upper bound of 999. That is not an estimate of 999. It is the optimiser
reporting that the foreground has no synonymous signal left to normalise a
rate against, which is exactly what §10.5 predicts for a branch this old.
Restarts reaching the same likelihood put one of them at 162 and at 999, which
is a six-fold spread at an unchanged likelihood and is the definition of an
unidentified parameter. On one of them, no site clears the posterior bar at
all.

On the ITPR1 stem the rate is **5.53 on 11 % of sites**, well inside the
estimable range, stable across restarts, with 8 sites clearing the posterior
bar [R147] and 5 clearing the stricter one.

**So the reportable branch-site finding is ITPR1 alone**: a class of sites on
its stem evolving several times faster than neutrally while the rest of the
tree sits near 0.03. The other two stems' tests are significant and their
effect size is unmeasurable, and those are different sentences. A pipeline
that printed the three p-values would have reported its strongest signal on
the stem whose parameter is least determined.

## 10.8 A significant site-model test can mean something quieter than its p-value

Site models ask whether any site in a paralogue has a rate above one across
the whole clade, which is a different question from the stem, and the one that
would find recurrent positive selection anywhere in the vertebrate history of
that copy.

The pattern that recurs is a significant likelihood-ratio test whose fitted
parameter says something quieter than the p-value suggests: an extra rate
class at exactly 1.000 on 0.3 % of sites, or a class above 1 with a share of
zero, and no site clearing the posterior bar. A significant test with no site
above the posterior bar is the signature of a fit driven by saturation rather
than by identifiable substitutions, so every significant test in this chapter
is printed beside its posterior column.

That discipline is the reason three of this chapter's significant tests are
reported as meaning something other than what their p-values suggest, and in
every case the tell was a parameter printed outside the test.

## 10.9 Selection on ITPR1 is intensified, and on its sisters relaxed

A branch model asks whether a foreground's rate differs. The relaxation test
[R150], run locally in the same framework as the per-site model [R149], asks
whether the whole distribution on the test branches is pulled towards
neutrality or away from it, which in a family where every rate is far below
one is the sharper question.

ITPR1 is intensified at k = 9.36, ITPR2 relaxed at k = 0.908, and ITPR3
relaxed at k = 0.836.

That is the same ordering the one-ratio estimates give, arrived at by a
different statistic on a different model, in that one compares point estimates
and the other compares whole distributions, so the two are a check on each
other rather than one number told twice.

**The three copies have not been held to the same standard since they were
made.** ITPR1 is under roughly twice the purifying selection of its sisters
and is intensifying relative to them, and it is the copy the teleost
duplication kept, the copy whose neighbourhood is best preserved, and the copy
that carries the dominant disease burden.

The unlabelled tips the tree places in no paralogue clade are left unlabelled
in these runs rather than swept into the reference set, because a branch whose
paralogue identity is unresolved is not evidence about either side of a
contrast.

One implementation note is worth keeping because the failure it causes is
silent in the wrong direction. The descriptive model estimates a per-branch
rate, and a branch with no synonymous change gets an infinite one, which the
program writes as a bare token that is not legal JSON. It is normal output,
the analysis succeeds and the parse throws, and on the first run that turned
one paralogue's entire result into a blank row.

## 10.10 What this chapter does not establish

**Saturation is not repealed by using a tree.** Nearly nine tenths of all
within-paralogue pairs exceed the saturation bar. Tree-based models handle it
far better than pairwise estimation, but a rate estimated where the synonymous
rate is poorly determined is a ratio with a soft denominator.

**Fourteen of the coding sequences come from genomic reconstructions** with
masked codons. The curated subsets show the estimates barely move without
them, but those models are also the only evidence for several lineages, so the
subset is a control rather than a replacement.

**The tree is conditioned on.** Every branch test runs on Chapter 6's
topology. If the sister arrangement were wrong, the branches marked here would
be the wrong branches, which is why that chapter ran a topology test over all
three arrangements before this one started, and why the dependency is recorded
rather than quietly relied on. A branch test cannot corroborate the topology
it is conditioned on.

**This is a vertebrate result.** The non-vertebrate grade is not in the codon
alignment at all, so nothing here bears on the constraint acting on the
single-copy receptors Chapter 5 found across the eukaryotes.

**Two of the three branch-site effect sizes are not identified**, which is not
the same as being large.

The whole suite ran unattended in 22.7 hours, and every stage renders whatever
has landed and marks an unfinished section as unfinished, because waking up to
a partial analysis that says which parts are partial is worth more than waking
up to nothing.
