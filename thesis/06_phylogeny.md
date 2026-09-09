# 6. Building the alignment and the tree that later chapters stand on

## 6.1 What this chapter has to produce, and why the care is warranted

Chapters 7, 10, 11 and 12 all rest on a single multiple sequence alignment. It
is worth the care that goes into it, and it is worth saying which decisions in
building it are measurements and which are choices.

The requirement is a set of representatives spanning the whole family, meaning
every clade and every kingdom in which Chapter 5 found the receptor, with the
ryanodine receptors included to root the tree.

## 6.2 Eight rules choose the representatives, and none of them uses length or identity

Eight rules choose the set, and the first thing to say about them is what none
of them does. **No rule ranks or filters on sequence identity, and no rule
takes the longest record per species.**

Longest-per-species is the obvious rule and it reliably selects chimeric gene
models, because a mis-joined model is longer than the real gene, so a rule
that prefers length prefers the error. Identity was retired as a call gate
outside the vertebrates in Chapter 5, after it was measured and found to
separate neither confirmed from contradicted loci, and a selection rule stated
in a statistic that does not separate the populations is not a rule.

What the eight rules do select is the following. They take the vertebrate
paralogue grid of clade band by paralogue, with human forced in. They take the
teleost co-orthologue pairs as a stress test. They take the deep vertebrate
grade, meaning cyclostomes, cartilaginous fish and the coelacanth grade,
entered **unlabelled**, because a paralogue label is a vertebrate concept and
assigning one here would assert what the tree is being run to test. They take
non-vertebrate metazoa by phylum. They take non-metazoan eukaryotes, gated on
Chapter 5's per-record plant and fungal verdicts, so that a suspected
contaminant cannot become the plant branch's tip. They take the copy-number
expansions, because Chapter 5 found up to 18 copies in one genome and a single
tip would misrepresent that. They take the novel gene models no database
annotates. And they take the ryanodine outgroup.

The quality key is **lexicographic rather than a weighted sum**, precisely so
the audit can name the component that decided each pick. Every chosen row
carries its rule, the cell it filled, the runner-up it beat and which
component separated them.

**Length targets are measured rather than borrowed.** Each group is scored
against its own median over its non-fragment records, because the family's
size varies far more outside the vertebrates than within, and one target would
score the plant and ciliate grades against a vertebrate ruler. The bar for a
record to count toward that median is four of five signatures rather than
five, and the record floor is three rather than five. Both were loosened for a
measured reason. At the strict setting two groups had too few qualifying
records and fell back to the global median: the green plants had two and the
amoebozoa one, so both were scored at 2,694 residues, a metazoan number, and
in the plant grade that decided which tip was chosen. At the looser setting
the well-populated groups barely move and the starved ones gain a real ruler,
in that the green plants go from two qualifying records to thirteen and their
target from a metazoan 2,694 to 3,182. The tier and record count behind every
target is committed, so a number computed from four records is visible as
such.

**One deviation from the brief is stated rather than absorbed.** The brief asks
for every novel gene model from both sweeps. There are 424 full-length models
with no database record, and including all of them would triple the alignment,
put it past the size where the chosen aligner is affordable single-threaded,
and, the real objection, fill the tree with vertebrate gene models, since 206
of the 241 vertebrate ones are birds and ray-finned fish. The rule instead
takes one per clade band and paralogue cell and one per non-vertebrate group,
which is what the requirement is for: that sequences no public database
contains are in the tree, and that each clade where they are the only evidence
has one.

The result is 134 representatives, 375,137 residues, the longest of them a
5,317-residue ryanodine receptor.

Six of the selection rules have constructed negative controls that run on
every build, covering a fragment, a ryanodine-length record offered as a
family member, a suspected contaminant, a protist whose UniProt name reads
"type 2" by annotation transfer, one species under two spellings, and a
score-ranked single-clade pool. Two further checks cover properties the step
before the aligner must have for reproducibility to mean anything downstream:
determinism, and Newick-legal unique tip labels.

## 6.3 The alignment is pinned to one thread because otherwise it is not reproducible

The aligner is MAFFT in its most accurate iterative mode [R138]: 61 minutes,
producing 11,777 columns at 76.23 % gaps.

**It is run single-threaded, by decision.** At automatic thread count it is not
reproducible on this machine, because the same seeds aligned twice in Chapter
3 gave 8,510 and 8,468 columns, since the iterative refinement stage combines
partial results in whatever order the threads finish. Everything downstream is
built from this file, so it is pinned and the SHA-256 of input and output are
recorded. A ragged alignment is a hard failure rather than a warning, because
a silent aligner failure degrades to a star alignment with no other symptom.

An automated trimming heuristic [R139] then keeps **1,797 of 11,777 columns,
15.3 %**. The choice of heuristic is recorded rather than tuned, because a
trimming threshold picked by looking at the resulting tree is a threshold
fitted to the answer.

![](figures/supp_representative_alignment.png)

**{fig:supp_representative_alignment}.** The alignment, with the columns the
tree actually saw marked in the input's own coordinates. The trimmed
alignment is the input with 9,980 columns deleted and the two share no
x-axis, so only one raster can honestly be drawn and the cuts marked beneath
it. The raster bins columns and plots occupancy rather than residue
identity, because at 11,777 columns one printed pixel is nine columns and a
residue palette would draw whichever residue happened to land on it. The
importance of this figure is that it lets a reader see what the tree was
actually computed on. Trimming removed 85 % of the columns, and the ones it
removed were the sparse ones, at a median occupancy of 0.03 against 0.98 for
those it kept. Every tree, selection estimate and constraint layer
downstream rests on the kept set, and a reader who wants to know how much of
the alignment that is can read it here rather than take it on trust.

![](figures/msa_coverage.png)

**{fig:msa_coverage}.** Per-sequence coverage of the trimmed alignment. The
median is 0.96 and one tip of 134 covers less than half. A tip below half is
not wrong, but it contributes gaps to every column the tree is inferred
from, so it is named rather than left inside a median. Coverage matters
because a tip that covers half the alignment contributes gaps to every
column the tree is inferred from, and gaps are not neutral for a likelihood
model. Naming the one tip below half, rather than reporting only the median,
is what allows a reader to judge whether any deep node depends on it.

## 6.4 What the alignment says before any tree is built

**The family separation, re-measured.** The ryanodine receptors are in this
alignment on purpose, and their separation from the family is a positive test
at every stage rather than an assumption. Measured here, mean identity within
a vertebrate paralogue group is 0.910 and from each paralogue group to the
outgroup 0.253, giving a separation of 0.656 identity units. Chapter 3's
benchmark, using a different estimator on a different panel, measured 0.828
and 0.249. The comparison is made on the separation rather than on either
absolute value, and it is confirmed.

That separation is what every stage of this project works inside. A
5,000-residue ryanodine receptor and a 2,700-residue IP₃ receptor are 25 %
identical over the columns they share, which is not far enough apart to trust
a heuristic with.

![](figures/msa_identity_heatmap.png)

**{fig:msa_identity_heatmap}.** All-pairs identity across the representative
set. Two identity matrices are committed, one scoring only mutually covered
columns and one counting gaps, because they answer different questions and
disagree systematically where a fragment is involved. The importance of
committing two identity matrices rather than one is that they answer
different questions and disagree systematically wherever a fragment is
involved. A metric that counts gaps as mismatches makes every partial
sequence look divergent, which would confound sequence divergence with
assembly quality, and that confound runs through this entire project.

![](figures/msa_group_identity.png)

**{fig:msa_group_identity}.** Identity within and between groups. The three
vertebrate paralogues sit at 0.74 to 0.79 to each other and 0.91 within
themselves. This figure fixes the scale every later comparison is read
against. The three paralogues sit at 0.74 to 0.79 to one another and at 0.91
within themselves, which is close enough that a bait from one aligns at all
three genes. That single fact is why the duplication detector of Chapter 8
needed the family call before its geometry meant anything, and why the
attribution margin inherited from a less similar family had to be
overturned.

**The alignment's preview of the sister question turns out to be wrong.** The
alignment can rank the three between-paralogue identities: ITPR1 with ITPR2
leads at 0.791, against 0.753 and 0.741, and the leading pair's interquartile
range does not overlap either of the others. That is a real ranking and it is
not a phylogenetic estimate, because it ignores the outgroup, the rate
variation and the branch lengths, so it was reported as a preview with the
tree named as the answer. Section 6.7 reports that the tree contradicts it.

**The alignment's preview of the cyclostome loci is wrong in a different way.**
Every cyclostome in this set carries three loci, all won by the same bait, so
nothing before the alignment could say which was which. On identity all six
fall nearest ITPR1, at margins of 0.014 to 0.066.

That table needs a control, because leaning towards one paralogue could be a
fact about the cyclostomes or a fact about the metric. If ITPR1 is simply the
slowest-evolving of the three, everything deep is nearest it. The same
statistic over groups that are certainly not vertebrate paralogues gives the
no-signal baseline: invertebrate metazoa lean toward ITPR1 at a median margin
of 0.008, protists at 0.002, and the ryanodine receptors at 0.002. The
cyclostome margin is five times that, so the lean is real. What it means is a
different question, and §6.8 gives an answer neither reading predicted.

![](figures/msa_conservation.png)

**{fig:msa_conservation}.** Per-column conservation with human ITPR1's
domain architecture mapped onto it through the alignment, by walking the
human row and counting ungapped positions, rather than by scaling residue
coordinates onto column coordinates. The latter is the mistake that puts a
domain boundary in the wrong place by exactly the gap content of the
sequence. The importance of mapping the architecture through the alignment
rather than by scaling coordinates is that the second method places every
domain boundary wrong by exactly the gap content of the sequence. Since
Chapters 11 and 12 both make claims about which element a residue belongs
to, a systematic offset here would have propagated into every one of them.

## 6.5 The exhaustive model scan was measured and abandoned

The obvious command is an exhaustive model scan [R141]. It was started, timed
and abandoned, and the measurement is the reason: **11 models of up to 1,232
in 672 seconds, projecting about 20.9 hours.**

The cost is the free-rate models, and they are not optional on this alignment,
because the best free-rate model beats the best gamma model by 682.3 BIC
units. Dropping them to buy the time back would have returned a measurably
worse model.

The replacement is greedy in two stages, with both tables committed: every
exchangeability matrix at fixed rate heterogeneity, then every rate model on
the winning matrix. Stage A chose its matrix by 614.6 BIC units over the
runner-up, and stage B bought a further 1,109.4 by changing only the rate
model.

**What that costs is stated rather than buried.** The matrix that wins under a
gamma model need not be the matrix that would win under the final rate model.
That is the bet a greedy selection makes, and the margin it was made on is in
a committed table rather than in a sentence, so how close the second-best
matrix came is checkable.

## 6.6 The tree search is pinned and guarded, and its rooting doubles as a control

The analysis used one maximum-likelihood search [R140], taking 71 minutes,
with 1,000 ultrafast bootstrap replicates [R142] and 1,000 approximate
likelihood-ratio replicates [R143]. **Threads are pinned rather than
automatic**, because automatic thread count reads the machine's current load
and the search is reproducible only at a fixed seed and a fixed thread count.
The number pinned is what the program's own benchmark returned for this
alignment, so pinning costs nothing.

Support is reported as both branch tests together, and a node counts as well
supported only when it clears both thresholds. One alone is not enough,
because the bootstrap is optimistic under model violation, which is exactly
the condition a 134-tip alignment spanning four kingdoms is in.

Two guards on the search itself came out of one incident. Every stage resumes
on the report file rather than the tree file, because the program writes a
tree at checkpoints and a search killed part-way leaves a tree of an
unfinished topology that the next run would silently reuse. And a run refuses
a prefix that a live process already owns. Both exist because a previous
session's constrained searches were interrupted, kept running orphaned for
fifty minutes, and a fresh run started three more on the same prefixes, after
which the topology test was handed whichever checkpoint tree happened to be on
disk.

**Rooting doubles as a control on the whole alignment.** The tree is rooted on
the ryanodine receptors, and if they did not come back as a clade the
alignment underneath every downstream result would be the thing to doubt. They
do, at maximal support.

![](figures/tree_ml_rooted.png)

**{fig:tree_ml_rooted}.** The rooted maximum-likelihood phylogram. Branches
are in neutral ink, the three paralogue clades and the outgroup are boxed,
and a filled dot marks every node clearing both support thresholds. No tip
is ringed, because the relabelling rule described in §6.7 fires on none of
them, and the legend entry for a ring appears only when a ring does, since a
key naming a marker the figure does not carry asserts a correction that was
never made. The importance of this tree is that it answers a question the
literature did not have an answer to. The review's own audit found no
published, support-annotated maximum-likelihood analysis with a ryanodine
outgroup that fixes which two of the three paralogues are sisters. The
rooting doubles as a control on the whole alignment, since a superfamily
that did not come back as a clade would put every downstream result in
doubt.

## 6.7 The three paralogues are clades, but only after the tree's own corrections

Every stage since Chapter 3 has assumed they are clades: the architecture call
assigns a paralogue per record, the genome ledger has one cell per paralogue
per genome, and the representative grid is built on clade band by paralogue.

Asked on the census's own labels with nothing corrected, **none of the three
is monophyletic.** ITPR1's largest pure clade holds 13 of 15 tips, ITPR2's 10
of 11, and ITPR3's 16 of 18.

A label is not evidence, and a single mis-annotated tip breaks a clade of
forty. So the same question is asked after the tree's own corrections, derived
by coded rule and never by hand. A tip inside its group's largest pure clade
keeps its label. A tip outside it whose smallest containing clade holds at
least three members of exactly one other paralogue, at a well-supported node,
is moved. Anything else is left out of all three sets rather than assumed
either way.

Under tree-corrected membership **all three are monophyletic.** The difference
between the two tables is the size of the family's naming problem rather than
a change in the tree.

**The correction rule fires on nothing.** 39 tips agree with their label, 0
are reassigned, and 5 are left unplaced. That is worth saying plainly: the
tree does not overturn a single census name in this set.

The rule that produces that result was itself wrong on its first version, and
the way it was wrong is instructive. It climbed the tree from a mislabelled
tip until some paralogue dominated and then assigned the tip there, which
hands a tip to whichever clade is largest three nodes up, and that is a fact
about clade sizes rather than about the tip. It now stops at the tip's own
nearest neighbourhood. Two constructed controls, a positive and a negative
half, caught it: a weakly supported placement must not overrule a census
label, and a strongly supported one must.

**Every naming call the tree disputes is checked by a different instrument.**
Neither the alignment nor the model can tell a mis-annotated record from a
tree error, because both produced the placement in question. So each of the
five unplaced tips is put through reciprocal best hits against the human
reference proteome and back, using an instrument that knows nothing about this
alignment, this model or this tree. **All five uphold the census name.** For
each of them the tree offered no alternative, only a refusal to place, so what
is upheld is the name: the record is correctly named, and the conflict is this
tree's own uncertainty about where to hang it.

That verification is set up so that it can fire. An earlier version filtered
to the reassigned tips only, which made it report that the tree contradicted
no census label while the table it had just read held five tips the tree
declined to place. A verification that cannot fire is not a verification.

![](figures/paralog_placement.png)

**{fig:paralog_placement}.** Where the tree places each tip against its
census label, with the five unplaced tips named. The importance of this
figure is that the tree overturns no census name. Thirty-nine tips agree
with their label, none is reassigned, and five are left unplaced, and
reciprocal best hits against an independent database uphold the census name
for all five. A gene tree that had disagreed with the annotation would have
cast doubt on the paralogue cells that every later chapter is built on.

## 6.8 ITPR2 and ITPR3 are sisters, and the other two arrangements are rejected

Three rooted hypotheses were tested, one per way of pairing the three
paralogues, each realised as a constrained search under the same model, and
each making all three paralogue clades monophyletic, so that the test compares
the sister arrangement and nothing else.

Two traps had to be avoided and the second was walked into first.

Membership must be **tree-corrected**, because a constraint built from census
labels asks the test about a topology the data rejects for reasons unrelated
to the sister question, and all three hypotheses then fail together. That trap
is documented, and the correction here is derived from the tree rather than
hand-listed.

The second trap is that the constraint mechanism places freely exactly the
taxa a constraint omits, so a tip that is listed, even in a top-level
polytomy, is pinned outside every group the constraint declares. The first
constraint files named all 134 tips and thereby forced the unlabelled
vertebrate tips out of the paralogue clades this tree nests them in,
identically in all three hypotheses. The test then rejected every arrangement
at a log-likelihood difference of about 1,500, **including the one the tree
itself holds at maximal support.** Each constraint now names the three
paralogue cores and the outgroup and nothing else, and a constructed control
fails any constraint that names a free tip.

![](figures/sister_au.png)

**{fig:sister_au}.** The three hypotheses under the approximately unbiased
test [R144], with log-likelihood difference and p-value on two panels
sharing one row of categories rather than on two y-axes. They are different
measures on different scales, and a twin axis invites a comparison that has
no meaning. This figure is the answer to Chapter 6's central question, and
its importance is that it excludes rather than merely prefers. Two of the
three possible sister arrangements fall outside the 95 % confidence set of
topologies for this alignment, at p below 2 × 10⁻⁵, so the result is not a
ranking that a slightly different alignment could reverse. Chapter 7 then
shows the genomic neighbourhood pointing the other way, and the two are
reported as a disagreement rather than reconciled.

**The unconstrained tree groups ITPR2 and ITPR3 at maximal support**, and the
test agrees. The arrangement pairing ITPR1 with ITPR2 is rejected at
p = 1.8 × 10⁻⁵, ITPR1 with ITPR3 at 1.65 × 10⁻⁵, and ITPR2 with ITPR3 is not
rejected at p = 0.476.

The alignment's identity preview predicted ITPR1 with ITPR2 and the test
returns ITPR2 with ITPR3, so **the preview is contradicted**, which is the
reason it was reported as a preview.

This is an answer the literature did not have. The review's own audit found no
published, support-annotated maximum-likelihood analysis with a ryanodine
outgroup that fixes the pair.

## 6.9 The cyclostome loci are two ancient lineages, which neither earlier reading predicted

The alignment left two hypotheses. Three one-to-one orthologues from the
vertebrate duplications would put each cyclostome locus inside a different
paralogue clade, while a lineage-specific expansion would put all of them
together, outside all three.

The tree gives a third answer. **Every one of the six loci sits in a
cyclostome-only clade.** There are two such clades, each well supported, and
both hold hagfish and lamprey together. So the loci are neither one expansion
nor three one-to-one orthologues. They are two anciently separate cyclostome
lineages, and a clade spanning the two species is a duplication older than the
hagfish and lamprey split rather than a copy either genome made on its own.

The identity preview is contradicted as a reading rather than as a
measurement. Identity ranked all six against the same paralogue because it can
only measure distance to the labelled clades, and these loci are outside all
of them. Which side of the vertebrate duplication each lineage attaches to is
Chapter 7's question, and the neighbourhood is what would settle it.

## 6.10 A stress test on the naming failed, for the right reason

The teleost co-orthologue pairs are in the representative set as a stress test
on the naming: if the tree does not recover a species' two same-paralogue
copies as sisters, either the naming is wrong or the alignment is.

**Neither set comes back as a clade.** Both are broken only by other tips of
the same paralogue, in that the copies sit in a small clade made entirely of
that paralogue, with each copy nearer another species' copy than its own
genome's.

That is exactly what a duplication older than the species looks like, and the
teleost genome duplication is precisely that. The naming passes. What fails is
the expectation that co-orthologues of a shared duplication should be sisters,
which was the wrong prior rather than the wrong tree.

## 6.11 How much of the tree is resolved, and what the model-violation guard says

A sister question answered on a tree whose deep nodes are unsupported is not
answered. Of 131 internal nodes, **91 clear both thresholds, which is 69.5 %**,
with a median bootstrap of 100 and a median branch-test value of 99.

![](figures/support_profile.png)

**{fig:support_profile}.** Node support drawn as a scatter rather than as
two histograms, because the claim is about the joint condition and a pair of
marginals cannot show it. Drawing support as a joint condition matters
because a node is only trustworthy when both branch tests agree, and two
marginal histograms cannot show that. The bootstrap alone is known to be
optimistic under model violation, which a 134-tip alignment spanning four
kingdoms guarantees, so every well-supported claim in this thesis requires
both.

Every node a claim rests on is listed with its own support. The ITPR2 with
ITPR3 clade the sister result depends on is at maximal support. The one weakly
supported claim node is ITPR1's core clade, meaning the census's labelled set
taken bare, at a branch-test value of 47.8, which the tree prefers to group
slightly differently. The same clade including the unlabelled tips is at
maximal support.

**A model-violation guard was run over the whole search.** Because the
bootstrap is optimistic under model violation, the whole search was re-run with an extra optimisation round on
every bootstrap tree, and the same clade questions were asked of that tree,
with membership read from a committed table rather than re-derived so that a
bookkeeping difference cannot surface as a topology change.

**No claim is weakened or lost.** Nine of ten clear both thresholds in the
guarded tree as well. The tenth was not well supported in the main tree
either, so the guard took nothing away. It is listed rather than folded into
the count, and the guard does push it further down, which is the honest
reading.

The tree reported is also the best tree found, because no constrained search
reached a higher likelihood than the unconstrained one. That check exists
because a constrained search can reach a better optimum than the free one, and
if it did, the tree every table is built from would not be the global optimum.

## 6.12 What the tree settles, and what it hands to Chapter 7

The tree settles that the three vertebrate paralogues are clades, that ITPR2
and ITPR3 are sisters with the other two arrangements outside the 95 %
confidence set, and that the cyclostome loci are an ancient cyclostome-specific
pair of lineages rather than orthologues of the three.

It does not settle when the duplications happened, where on the vertebrate
tree they sit, or whether they are the two rounds of whole-genome duplication
at all. A gene tree gives branching order, it does not give dates, and it
cannot distinguish a genome duplication from a tandem one. That is Chapter 7,
and it needs a different kind of evidence entirely: the genomic
neighbourhood.
