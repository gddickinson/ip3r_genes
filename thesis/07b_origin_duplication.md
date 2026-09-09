## 7.7 Dating the duplications requires a species tree, which is an input

A reconciliation maps a gene tree onto a species tree and reads duplications
off the mapping. The species tree is therefore an input, and it is declared as
one rather than estimated here.

It is a hand-curated topology for the 31 vertebrate species the alignment
samples, with 29 named internal nodes, every one carrying a literature age,
the spread of published estimates, and its source [R188, R189, R190, R191,
R192]. A validator fails the build on a gene-tree species missing from the
tree, a taxon with no gene-tree tip, an uncalibrated node, a point age outside
its own spread, or any age inversion.

Two things it does that a cladogram with ages written on it does not.

**Where the literature does not resolve an arrangement, the node is left a
polytomy** rather than given a shape the sources do not carry.

**Every node carries a stem age beside its crown age.** The sampled tree omits
every unsequenced lineage, so a node's parent in this tree is usually older
than the split that created it. Without the stem age, a duplication on the
placental mammal branch would be bracketed at 96 to 319 Ma instead of 96 to
99.

The reconciliation itself uses the **non-binary** duplication rule [R174]
rather than the familiar binary one, and that is load-bearing rather than
fastidious. The binary rule is valid only on a binary tree, and on the
support-collapsed root, which is a four-way polytomy, it reads a speciation,
which would have reported the collapse as the duplications disappearing.
Losses are counted by the standard method [R173] and the mapping is the
standard least-common-ancestor one [R172].

**Fourteen groups of constructed negative controls run before anything is
written**, and two of them caught real errors: the two halves of the
non-binary duplication rule, and the support-collapse routine refusing a tree
with no support labels. That last exists because a constrained search writes
no support values, and the first run dissolved all three constrained
topologies into one 134-tip polytomy.

## 7.8 The two duplications sit on different branches

Twelve reconciliations were run. Five topologies, comprising the
maximum-likelihood tree, the three constrained sister hypotheses and the
model-violation re-search, were crossed with three variants: all tips, the
cyclostome loci pruned, and every node below the support bar collapsed. Three
combinations are refused with their reason.

![](figures/recon_matrix.png)

**{fig:recon_matrix}.** The topology by variant matrix, with the deepest
paralogue duplication in every cell, and beside it every rooting of the
vertebrate subtree scored by total events, with the outgroup rooting and the
minimum-event rooting marked. The matrix matters because it separates the
placement from the topology it was read off. Twelve reconciliations across
five topologies and three taxon treatments put the younger duplication on
the gnathostome stem in every cell, so that placement does not depend on
which sister arrangement is correct. The rooting panel beside it shows the
same for the choice of root.

Read plainly, on the tree with every tip in, two results follow.

**The split that separated ITPR1 from ITPR2 and ITPR3 is on the vertebrate
stem**, meaning older than crown Vertebrata [R186, R187], whose published
estimates run 480 to 615 Ma. Nothing in this tree bounds it from above, and
the bracket is reported open rather than closed silently.

**The split that separated ITPR2 from ITPR3 is on the gnathostome stem**,
meaning between crown Gnathostomata at 462 Ma and crown Vertebrata at 563 Ma.
It maps there in every cell of the matrix, cyclostomes in or out.

![](figures/recon_dated_backbone.png)

**{fig:recon_dated_backbone}.** The dated backbone on a linear time axis
with each calibration's spread drawn as a band, rather than a cladogram with
ages written on it, because the result is an interval and a cladogram cannot
show one. Placements are deduplicated to one marker per distinct
arrangement, since the matrix repeats the same placement across cells and
drawing each would make agreement look like weight. The deepest bracket is
left open at its old end because nothing in this tree closes it, and the
widest disagreement in the whole calibration set sits at exactly the node
the older duplication maps to.

**The two duplications are not on the same branch.** That is the finding, and
it is what makes ITPR1 the earlier-diverging copy.

## 7.9 Two objections to the deep placement were measured rather than argued

Two objections a reader should raise are measured rather than argued.

**Does the placement rest on the weakest node in the tree?** The gene-tree
node separating ITPR1 from ITPR2 and ITPR3 sits at a branch-test value of 17.4
and a bootstrap of 54, which is the weakest node in the whole vertebrate
subtree and apparently the one carrying the headline. It is not. Collapsing
every node below the support bar dissolves that arrangement and the placement
survives, because the root becomes a polytomy in which more than one child
still maps into the gnathostome subtree, and the non-binary rule reads a
duplication at Vertebrata from that alone. What the collapse costs is
resolution, in that two separate vertebrate-stem duplications merge into one,
rather than the placement.

**What forces the older placement is one hagfish and lamprey pair.** Of the six
cyclostome loci, two sit in a cyclostome-only clade that is itself nested
among the gnathostome paralogues, sister to ITPR2 and ITPR3 at maximal
support. A lineage holding both cyclostome and gnathostome copies forces every
duplication above it to predate the split between cyclostomes and
gnathostomes. The same tip is what keeps the younger event on the gnathostome
stem, because that pair is sister to ITPR2 and ITPR3 rather than inside
either, so the node uniting ITPR2 with ITPR3 holds gnathostomes only. Had it
fallen inside either, that split too would have been pushed onto the
vertebrate stem. **The two answers are decided by where one two-tip clade
attaches**, and saying so is the honest form of the result.

**Are those tips long branches?** Cyclostomes are the classic long-branch
attraction risk in vertebrate phylogeny, and long branches are attracted to
the root, which is exactly where this answer sits. So the objection is
measured. Root-to-tip distance for all 57 vertebrate tips puts the six
cyclostome loci at **0.96 to 1.04 times the median, ranking 18th to 55th of
57.** Not one is an outlier and two sit in the shorter half.

![](figures/recon_cyclostome.png)

**{fig:recon_cyclostome}.** The six cyclostome loci marked on the
distribution of root-to-tip distance across the tree, and beside it the
independent neighbourhood call for each locus with the pair support the tree
gives it. The first panel removes the standard objection to this placement
by measurement rather than by argument. The second matters because the
neighbourhood knows nothing about the alignment, so its agreement is
corroboration from an instrument that could have disagreed.

That negative result is the one that matters most here. It is the reason the
placement is offered as a finding rather than as a caveat.

**A loss count from this analysis is not a loss count.** The reconciliation
implies 47 lost cells. Checked one at a time against the genome ledger, 26 of
them are the paralogue present in the genome and absent only from the 134-tip
sample, and none is corroborated.

![](figures/recon_losses.png)

**{fig:recon_losses}.** What each implied loss turns out to be once it is
asked of the genome ledger, and the implied count in every cell of the
matrix against the number the genomes corroborate. The figure is the reason
a loss count is not read off a reconciliation: a reconciliation over a
representative sample counts sampling, not biology. That is why Chapter 9
counts losses from 309 genomes instead, and why the two chapters' numbers
are not alternative estimates of the same quantity.

## 7.10 Testing whole-genome duplication needs dated paralogy, not synteny

Neighbourhood similarity can show that two blocks are related. It cannot show
when they became related, and that distinction is what separates an ohnologue
pair from the two rounds of duplication from an older duplication whose two
copies happen to sit beside these genes.

So the human paralogy map was taken from a comparative-genomics database
[R161] through a **pinned dated archive** rather than the rolling one [R162].
That is the same reproducibility discipline applied to a database release,
since the rolling host follows the release cycle and a re-run would score the
same windows against a different tree. The archive's own registry is committed
beside the map, so the release is evidence in the results directory rather
than a sentence in a report.

Three rules make it a test rather than a description.

**Every IP₃ and ryanodine receptor gene leaves every window.** Counting the
ITPR1 and ITPR2 paralogy as evidence that their blocks are paralogous is
circular, and leaving a ryanodine receptor in an IP₃ window would import the
control's answer into the test.

**The null is drawn from real genomic windows** at the matched gene count, so
the clustering of gene families that a shuffled gene set throws away is kept,
and it inflates the real window and the null alike.

**The links are dated.** The database's duplication node is what separates an
ohnologue pair from the two rounds from an older duplication, and it is read
in two nested vocabularies.

**Dating the links changed the answer.** Two paralogue links survive between
the three human neighbourhoods at every window size, and they are exactly the
two families the neighbourhood analysis found by a completely different
instrument. But one of them, the basic helix-loop-helix pair linking ITPR1 and
ITPR2, is dated to Opisthokonta. It is a pair far older than the vertebrates
whose two copies happen to sit beside these genes. It is real paralogy and it
is not an ohnologue pair from the two rounds. The metabotropic glutamate
receptor pair linking ITPR1 and ITPR3 is dated to Vertebrata, which is what
the two rounds mean.

**The single human genome is underpowered, and the report says so twice.** At
the primary window each family retains exactly one vertebrate-dated ohnologue
pair between two of its three neighbourhoods, each clears its own permutation
null when asked once per family rather than three times per pair, and neither
survives correction across all 135 tests the stage ran. The ryanodine trio's
two-round origin is not in question, so what the human test measures is the
instrument's ceiling on one genome rather than a difference between the
families.

**The replication is where the answer is.** The same human map, asked of the
flank sets in every swept genome, against matched random neighbourhoods in the
same genomes drawn by the same sampler, gives the following. The ITPR1
neighbourhood carries a paralogue of the ITPR2 neighbourhood in 141 of 175
genomes, which is 80.6 %, and of the ITPR3 neighbourhood in 89 of 152, against
2.6 % of matched random windows. **ITPR2 and ITPR3 sit at exactly the
background rate, at 2.7 % against 2.6 %, p = 0.554.**

![](figures/s16_paralogon.png)

**{fig:s16_paralogon}.** The two-round test with its null drawn across the
bars rather than quoted in a caption, and the ryanodine trio run through the
identical instrument in the same genomes beside it. What carries the claim
is not the human window but the replication across 309 genomes against
matched random neighbourhoods, and what makes it readable is the sister
family, whose two-round origin is not in question and which behaves the same
way through the same instrument.

The dated column splits the two links cleanly: the ITPR1 to ITPR3 link is
vertebrate-dated in 84 genomes and the ITPR1 to ITPR2 link in 2.

**The fourth slot cannot be found**, and current reconstructions of what the
vertebrate duplications left behind [R182, R183] are what a positive answer
would have to be read against. A genome-wide block scan returns no block that
both carries no family gene and looks like these blocks' missing sibling. With
one dated ohnologue pair surviving between the blocks that do carry a gene, a
block that lost the gene too has nothing left to be recognised by. That is a
limit of the evidence rather than a claim that no fourth slot existed. The
quartet test additionally flags its circular pairs in the data, because a
window tested against its own top-scoring block tests the selection rather
than the quartet, and all six such pairs are significant by construction.

## 7.11 The teleost duplication doubled ITPR1 and only ITPR1

Ray-finned fish underwent a further whole-genome duplication [R181, R184], and
the copy-number landscape says that is where all this family's variation is.

![](figures/s16_copy_number.png)

**{fig:s16_copy_number}.** Copy number grouped by whole-genome-duplication
status rather than by taxonomy, because the result is a contrast between
lineages defined by which duplications they have been through, and a
per-class bar buries it inside the ray-finned fish. Grouping by duplication
history rather than by taxonomy is what makes this figure a test. The
prediction is specific: lineages that diverged before the teleost
duplication should carry one copy of each paralogue, teleosts two of one,
and lineages with a further duplication more again. All three hold, and a
second copy in the outgroup or an undoubled sister family would each have
falsified the reading.

**Above the contiguity bar, every non-teleost gnathostome genome in the sweep
carries exactly one of each paralogue and three ryanodine receptors.**

Five predictions were checked.

**Pre-duplication lineages carry one copy.** The bichir, gar and bowfin carry
1.00 of each and 3.00 ryanodine receptors, which is the gnathostome state.

**Teleosts carry two of one.** They average 1.97 ITPR1, 1.04 ITPR2, 1.04 ITPR3
and 5.82 ryanodine receptors, with 97.3 % of 73 genomes carrying more than one
ITPR1.

**Lineages with a further duplication carry more again** [R185]. Sturgeon and
salmon carry 3.00, 2.00, 2.00 and 8.00.

A second copy appearing in the outgroup, or a teleost ryanodine count that had
not doubled, would each have killed the reading, and neither does. **The
ryanodine control is what makes the singletons interpretable**: the
duplication duplicated ITPR2 and ITPR3 as surely as ITPR1, and the sister
family in the same genomes kept all six of its copies. So the teleost ITPR2
and ITPR3 singletons are a statement about retention rather than about the
sweep's ability to find a duplicate.

**The copies are not tandem, and this is measured rather than asserted.** Of
71 teleost genomes with two copies above the bar, 20 put them on different
contigs and the rest sit a median 8,870,578 bp apart, with the closest pair
anywhere at 535,374 bp.

**Both copies keep part of the ancestral neighbourhood and between them
account for it.** Of 49 two-copy genomes with both copies flanked, 46 keep
ancestral symbols on both copies against a tetrapod reference and 47 against a
pre-duplication ray-finned one, and in 45 and 46 respectively **the two sets
are disjoint**. Disjoint is the load-bearing word: the copies do not merely
each resemble the ancestor, they partition it, which is what reciprocal gene
loss after one duplication produces and what a pair of independent later
duplications would not.

![](figures/s16_dcs.png)

**{fig:s16_dcs}.** Double-conserved synteny, with the two copies plotted
against each other so that disjointness is a geometric fact on the figure
rather than a number in a table. The load-bearing word in this figure is
disjoint. Two copies that each resemble the ancestral neighbourhood could be
two independent later duplications, but two copies that partition it between
them are what reciprocal gene loss after a single duplication produces. That
distinction is what turns a copy count into a claim about one shared event.

Two references were used rather than one, because teleost gene symbols diverge
from tetrapod ones even after normalisation, so a tetrapod consensus
under-counts. Both give the same answer.

## 7.12 The two teleost copies are one ancestral duplication, not many

The five predictions are all consistent with a single duplication in the
teleost ancestor. They are also consistent with a series of independent
lineage-specific duplications, and separating those is the sixth check.

The obvious statistic has no power, because each genome's copy labels are
arbitrary, so a matched-versus-crossed comparison measures nothing. Instead,
anchors from different orders assign every genome independently, and the
question is whether the anchors agree.

**705 of 705 assignments agree**, across six anchors from six orders. Under
independent lineage-specific duplications the anchors carry no shared
information and the statistic sits at 0.5, and a constructed control builds
exactly that case and requires the statistic to land there, so the perfect
agreement is a measurement rather than a property of the routine.

![](figures/s16_blocks.png)

**{fig:s16_blocks}.** Cross-anchor block identity, and the sensitivity of
every count to the coverage bar. Seven bars are scanned, because a
duplication claim that survives only one bar is a claim about the bar. The
anchor panel is what separates one ancestral duplication from a series of
lineage-specific ones, which the copy counts alone cannot do, because
anchors from six different orders assign every genome independently and
would carry no shared information if each lineage had duplicated on its own.

It is corroborated by evidence of a different kind. Which bait won each copy
is a sequence call made with no neighbourhood input at all, and in the six
genomes whose two copies won different baits the sequence call and the
neighbourhood block agree six times out of six. The reference for that check
is chosen as the first anchor whose own two copies won different baits,
because anchors are ranked on flank richness and taking the first one blindly
makes the check unrunnable.

## 7.13 What this chapter settles, and the disagreement it leaves standing

**Six things are settled.** The three paralogue neighbourhoods are shared
within a paralogue
at 216 to 413 times a matched null and share nothing across the family
boundary. The blocks are paralogous and the paralogy runs through ITPR1. One
of the two surviving links is vertebrate-dated and the other is far older. The
ITPR1 against ITPR2-plus-ITPR3 split sits on the vertebrate stem and the ITPR2
against ITPR3 split on the gnathostome stem, and the placement survives
collapsing every weakly supported node. The teleost duplication doubled ITPR1
and only ITPR1, in 97 % of teleost genomes, while doubling all three ryanodine
receptors in the same genomes, and the two copies are one ancestral
duplication.

**One disagreement is left standing.** The tree makes ITPR2 and ITPR3
sisters. The
neighbourhood makes ITPR1 the block that kept its ohnologues with both of the
others, and ITPR2 with ITPR3 the one pair that retains nothing above
background. Two instruments give two answers, and the project's verdict is
*not corroborated* rather than *contradicted*.

The reason is that they are not the same quantity. A tree estimates the order
in which copies diverged. A retained flanking ohnologue records which copies
survived deletion beside each gene, hundreds of millions of years later, and
quartets from whole-genome duplication lose flank copies independently of the
duplication order. The reconciliation adds a third view that is consistent
with both, because it makes ITPR1 the earlier-diverging copy, which is the
copy whose neighbourhood kept a shared family with each of the other two. That
is suggestive rather than a test.

Stating that as a disagreement rather than resolving it is a decision. The
alternative, which is quietly preferring whichever instrument agreed with the
headline, is how a project of this size talks itself into a wrong answer.

**Three things are not settled.** The chapter does not settle which of the two
duplication
rounds made which split, because the dated link is dated to Vertebrata, which
is both rounds, and separating them needs the cyclostome side of the quartet,
which §7.6 found underpowered. It does not settle whether the quartet had a
fourth slot. And it does not settle **why ITPR1 alone kept its teleost
duplicate**, because the observation is clean and the cause is not in this
chapter's evidence. Chapters 10 and 11 are where a dosage or
subfunctionalisation argument would have to be made.
