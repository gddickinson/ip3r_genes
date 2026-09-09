# 7. Where ITPR1, ITPR2 and ITPR3 came from

## 7.1 A gene tree cannot answer the questions this chapter asks

Chapter 6 established that the three vertebrate paralogues are clades, that
ITPR2 and ITPR3 are sisters, and that the cyclostome loci are something else
again. It could not say when the duplications happened, where on the
vertebrate tree they sit, or whether they are whole-genome duplications at
all.

Those are three different questions and they need three different kinds of
evidence. Where a gene sits relative to its neighbours is a statement about
the genome that no alignment contains. When a duplication happened requires a
species tree with dates on it, which is an input rather than a result. And
whether a duplication was part of a whole-genome event requires showing that
the neighbourhood was duplicated too, which needs a paralogy map rather than a
synteny map.

This chapter does all three, and the answers do not entirely agree. Where they
disagree, the disagreement is reported rather than resolved by preference, and
§7.13 says what it is.

## 7.2 Reading the genomic neighbourhood needs three decisions and a matched control

Every IP₃ and ryanodine receptor locus of every swept genome was taken with
its flanking genes, and every measurement below is a comparison between two
such neighbourhoods.

Three decisions make that comparison mean anything.

**Loci come from the per-genome sweep output rather than from the ledger.**
The ledger holds one row per genome by cell and therefore only each cell's
best locus. In a teleost a cell holds two genes, so a ledger-driven comparison
would score one species' first copy against the next species' second copy and
report the mismatch as a synteny result.

**Two window rules are committed rather than one**, because naming density is
not constant: human annotation names 97 % of coding genes and sea lamprey
27 %. One rule takes a fixed number of coding genes each side, and the other
takes the nearest informative symbols each side, so that a difference in
neighbourhood similarity cannot be a difference in annotation depth.

**Three key vocabularies are used.** The first is a strict symbol, the second
a relaxed one with teleost duplicate suffixes stripped, and the third a
**root** key with the trailing digit run removed. The last is the only level
at which a whole-genome-duplication ohnologue pair is visible at all, because
after 500 million years the two copies almost never carry the same symbol.

**Every real comparison is scored against a matched random control.** A
neighbourhood similarity of 0.21 says nothing on its own, so every pair of
loci is scored against control pairs of random coding genes in the same two
genomes, under the same window and key rules, which holds the two genomes,
their annotation depth, their naming conventions and the normalisation
constant fixed. Sampling is seeded per accession, and a contig too short to
hold a window is refused rather than allowed to degrade the null.

Eleven constructed negative controls run on every build. The most instructive
requires that **two empty neighbourhoods score zero rather than one**,
because otherwise two unannotated genomes are a perfect synteny match, and a
survey across 309 genomes of very unequal annotation depth would find its
strongest signal in its worst data. Another caught a real hash-seeded
nondeterminism in a ranked output table.

## 7.3 The neighbourhood confirms the paralogue cells at two orders of magnitude

Every stage since Chapter 3 has treated a cell's paralogue label as orthology.
Nothing before this had tested that with evidence outside the gene itself.

**It holds, by two orders of magnitude.** Within-paralogue neighbourhood
similarity runs at 216 to 413 times the matched null, with 98 to 99.8 % of
individual pairs beating their own control.

![](figures/synteny_pair_classes.png)

**{fig:synteny_pair_classes}.** Every pair class against its matched
random-window control. This is the first test in the thesis of an assumption
every earlier chapter had made, since the paralogue cells were assigned from
sequence and nothing before this had checked them against evidence outside
the gene itself. Neighbourhood similarity within a paralogue runs at two to
four hundred times a matched null while every cross-paralogue class sits at
or below it, so the cells are orthology groups rather than similarity bins.

**Across the family boundary there is nothing.** Over 168,241 pairs of IP₃ and
ryanodine loci the highest mean similarity of any class is 0.0002, which is
below the random-window control itself. An instrument that shares nothing with
the sequence evidence returns the same family separation.

**The three paralogues are not equal, and pooling hides it.** A pooled
within-paralogue mean answers two questions at once: whether two mammals share
a neighbourhood, which they nearly trivially do, and whether a mammal shares
one with a teleost. Only the second is about the locus.

![](figures/synteny_clade_decay.png)

**{fig:synteny_clade_decay}.** How far a neighbourhood travels, split into
within-class and cross-class comparisons. Within a vertebrate class the
three paralogues span 0.254 to 0.347. Across classes, ITPR3 retains 0.062
against 0.160 and 0.158, which is a 2.5-fold gap and the one asymmetry in
this panel that a pooled mean would have hidden.

**ITPR3's neighbourhood is the one that does not travel.** It is still 157
times its own null and 96.9 % of its cross-class pairs still beat their
control, so this is decay rather than absence. It is also worth noting what it
is not. ITPR3 is the best-recovered paralogue in the genome sweep and the
best-annotated. Ease of finding a gene and stability of the ground it sits on
are different properties, and in human that ground is the major
histocompatibility region at 6p21.31, one of the most rearranged and most
polymorphic neighbourhoods in the genome.

## 7.4 The surviving paralogon links run through ITPR1

If the three paralogues are the product of whole-genome duplication [R179,
R180], their neighbourhoods should be paralogous rather than identical,
meaning the flanking genes should be surviving copies of the same ancestral
families under different names. That is exactly what a same-symbol test
cannot see, and it is why cross-paralogue symbol similarity above is a flat
zero.

At the root-key level two links survive, and each was scored against the rate
at which random windows in the same genomes share a root family.

**ITPR1 with ITPR3 share a metabotropic glutamate receptor family**, present
in 42 % of one side and 53 % of the other against a background of 0.45 %,
which is an enrichment of 93-fold. **ITPR1 with ITPR2 share a basic
helix-loop-helix family**, at 62 % and 85 % against 0.74 %, which is 84-fold.

**ITPR2 with ITPR3 share nothing, at any threshold**, and nor does any pair of
IP₃ and ryanodine neighbourhoods.

![](figures/synteny_paralogon.png)

**{fig:synteny_paralogon}.** The three human neighbourhoods drawn as gene
tracks with the shared ohnologue families linked. The tracks are read out of
a committed table rather than named in the figure code, so the figure cannot
show a gene the data does not have. This figure carries the neighbourhood
half of the duplication argument. If the three paralogues arose in
whole-genome duplications, their neighbourhoods should be paralogous rather
than identical, and the surviving shared families are what remains of that.
The importance of the asymmetry is that both surviving links run through
ITPR1 and none connects ITPR2 to ITPR3, which is the first of five
independent measurements that put ITPR1 on its own side.

That is a clean measurement, since both families are under 1 % of random
windows at every threshold from 10 % to 50 %, and it is the first result in
this thesis that does not corroborate an earlier one. Chapter 6's tree makes
ITPR2 and ITPR3 sisters, and ITPR2 and ITPR3 are the one pair whose
neighbourhoods share nothing.

That is not a contradiction, and §7.13 says why at length. In short, a tree
estimates the order of duplication, while a retained flanking ohnologue
records which copies survived deletion beside each gene, and the quartets
left by whole-genome duplication lose flank copies lineage by lineage in a
way that carries no memory of the duplication order [R205, R183]. What can
be said is that the synteny does not corroborate the tree's pair, and that
whichever pair is sister, the ITPR1 neighbourhood is the one that kept its
ohnologues.

## 7.5 A caller built only from the neighbours places 131 unlabelled loci

The neighbourhood can also be used forwards. If each paralogue has a
characteristic set of neighbours, an unlabelled locus can be assigned by its
neighbours alone.

A consensus is built per paralogue from the keys carried by a stated fraction
of the species holding it, **with the query's own species dropped from every
consensus first**, so a locus cannot be scored against evidence it supplied.
The caller is calibrated on loci whose paralogue their own assembly's
annotation establishes, which is evidence the caller never sees, since it
reads the symbols of the neighbours and the truth is the symbol of the gene
itself.

**The threshold is measured, and the obvious way to measure it is wrong.** The
sweep reports what each setting buys, meaning the call rate on the confirmed
loci, and what it costs, meaning the rate at which random neighbourhoods in
the same genomes are called a paralogue. Optimising call rate alone selects
the loosest setting on offer, which is how an instrument gets tuned into
agreeing with itself. The rule is to maximise the difference, with ties broken
toward the stricter setting.

![](figures/synteny_caller.png)

**{fig:synteny_caller}.** The caller calibrated and swept across its
threshold range, with the random-window false-call rate drawn beside the
call rate. At the chosen setting it is correct on 405 of 405 calls over 503
annotation-confirmed loci, calls 80.5 % of them, and calls 6 of 726 random
control windows. Drawing both curves is what makes the operating point
defensible, since the two move together and only their difference says what
a setting is worth. The caller then places 131 loci that the alignment could
not label, each with a null tail probability attached.

Accuracy is 1.000 across the whole range, so it separates nothing and is not
what is optimised. It is reported rather than used.

The reporting bar is the **null's own maximum** rather than a probability cut.
The highest consensus overlap any random window reached is 2, so a call at or
below that score is reported as inside the null however clean it looks. With
726 control windows, a one-in-726 tail is not a rate to build a claim on.

**131 loci gain a paralogue assignment that no random neighbourhood could have
produced**, comprising 73 ITPR1, 42 ITPR2 and 16 ITPR3, each carrying its
score, its margin and its null tail probability. Thirty-six of them are second
and later copies inside one paralogue cell, and 36 of 36 are placed in the
same paralogue as the cell they were filed under. A cell holding two loci
holds two copies of one gene rather than a misfiled second gene.

## 7.6 The neighbourhood cannot reach the cyclostome question

Chapter 6 handed this chapter a specific question: the six cyclostome loci sit
in two ancient cyclostome-only lineages, and which side of the vertebrate
duplication each attaches to is undecided.

The neighbourhood cannot answer it. Of the six loci, four get a call at all,
and **every one of those sits inside its own null.** The other two are
no-calls.

That is not a negative result, it is an absence of measurement, and the
distinction matters enough to have its own verdict in this project's
vocabulary. A neighbourhood comparison taken where 27 to 29 % of the flanking
genes carry a symbol is a measurement that did not happen. The corroboration
this chapter would most like to have does not exist, and saying so is better
than reporting the four inside-null calls as evidence.
