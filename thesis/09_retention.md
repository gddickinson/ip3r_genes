# 9. Counting a loss that never happened

## 9.1 A zero is the hardest result in this thesis to state

Gene families of this age normally lose copies [R193, R194]. The prior for a
three-paralogue family present since the origin of the vertebrates is that
somewhere in 309 genomes at least one lineage has dropped one.

This chapter reports that **not one of the 927 genome by paralogue cells in
the vertebrate scope reaches the state this project defines as a loss**, that
Dollo parsimony [R175] places no loss anywhere on the tree, and that no rate
can be fitted to a character that never changes.

A zero is the easiest result to obtain by accident and the hardest to defend.
It is what a broken pipeline produces, what an over-conservative rule
produces, and what a search too insensitive to find anything produces. So most
of this chapter is not the count. It is the machinery that makes the count
mean something: a state vocabulary in which exactly one state can be counted
as a loss and every other explanation must be positively excluded first; an
instrument for the cells the sweep could not decide; a false-negative rate,
measured in Chapter 4, that says how often the search misses a gene that is
there; and a sensitivity matrix that names, cell by cell, the analytical
settings that would manufacture a loss out of this data.

## 9.2 Eight states, of which only one can be counted as a loss

The genome ledger's own statuses were not built to license a loss count. In
that ledger "absent" means the rescue search attributed no region, which in a
shattered assembly is a statement about contig lengths. So every cell is
restated under a chain of ordered positive tests, where the first rule that
fires wins and each writes the number it fired on into its own row.

The design constraint is that **exactly one state may be counted as a loss**,
reachable only by passing every test that could explain the absence otherwise.

The chain runs as follows. The genome's positive control did not fire, so no
cell in it supports any claim. Or there is a locus at or above the coverage
bar. Or a locus below it, truncated at a contig edge or a run of ambiguous
bases. Or a locus below it with no assembly excuse. Or no locus, but the
reference reassembles from sequence outside every locus the aligner found. Or
no locus and no reassembly, but the genome carries family loci that no cell
claimed. Or nothing found, and the assembly cannot hold the gene on one
contig. And finally, nothing found in a controlled assembly that could have
held it.

**Two of those rules exist because of specific incidents.**

The sixth is the cyclostome rule. Both cyclostome genomes carry three family
loci apiece, all filed in the ITPR1 cell because the bait panel has no
cyclostome-labelled bait, so their ITPR2 and ITPR3 cells read as absent in the
ledger. A count taking that at face value would score two independent losses
per cyclostome that never happened.

The seventh is the contiguity rule. An assembly whose contigs are shorter than
the gene cannot represent it as one locus, so it cannot be evidence that the
gene is missing, and two thirds of the bird assemblies in this scope are in
that position.

The rule order is itself tested. A constructed control moves the cyclostome
rule after the absence rule and requires the suite to fail.

## 9.3 Reassembling a reference across contigs decides the undecided cells

The sweep's undecided cells are the ones a loss count can neither ignore nor
use. The aligner placed no locus, so the ledger records no coverage, and the
rescue search attributed regions, so they are not nothing either.

What is measured is the **non-redundant coverage of the reference protein,
reassembled across contigs**. Three properties make that a positive test
rather than a hopeful one.

**It is computed outside every locus the aligner found**, passing through the
sweep's own exclusion set of every clustered locus in the genome, any bait,
padded. For a family whose paralogues are 61 to 68 % identical that is the
failure mode that matters, because without the filter a missing gene
reassembles out of its own paralogues and reports 0.95 coverage. A constructed
control removes the filter and requires the suite to fail.

**It uses one reference at a time**, never a union over baits, because a union
over three orthologous baits counts a residue covered in any of them and
reports a coverage no single protein achieves.

**The bar is measured against a gene that is accounted for**, with the
separation statistic [R177] reported beside the threshold. The negative
population needed no new search: it is the regions the full bait panel
attributes to a paralogue whose gene the aligner already placed at a locus in
the same genome. That paralogue is accounted for, so the fragment set cannot
be that gene, and what it recovers is what cross-paralogue similarity delivers
on its own, at a median of 0.030 against a candidate median of 0.795.

![](figures/reconstruction.png)

**{fig:reconstruction}.** The three populations the bar is read off, with the
gap shaded and the operating point drawn, and every undecided cell against the
number of contigs its gene is spread over. A calibration figure that asked to
be believed would not be one, so both edges of the gap are marks rather than a
caption.

**The first version of this calibration did not separate**, at a Youden index
of 0.52, and why is a result in itself. The decoy was every region attributed
elsewhere, but in a genome where two paralogues are both shattered, a fragment
attributed to the other one is a piece of a real gene. Those 20 regions are
now reported as their own population, with a median of 0.631 sitting squarely
between the two, and they are the **measured size of the paralogue-attribution
problem in a shattered assembly**. A constructed control folds them back into
the decoy and requires the suite to fail.

The operating point is the **midpoint of the gap** rather than the Youden
threshold, which on a perfectly separated pair lands on the lowest positive
and is the most permissive bar the data allow. Both edges are committed, so a
later run whose populations have drifted into contact is visible in the table,
and the calibration refuses to hand out a threshold it has marked unusable.

## 9.4 Synteny could not reach the undecided cells, and that is the result

The task was asked to disambiguate the undecided cells **by synteny**, with
measured accuracy. Chapter 7's caller already existed and was already
calibrated, so the work was to point it at these regions and measure whether
it arrives.

**It does not, and that is the result.** Of 432 rescue regions, **273 sit on a
contig carrying no annotated gene at all**, 151 have too few informative
flanking symbols, and 8 reach the caller's floor. The median region has zero
informative neighbours and zero coding genes on its own contig, which extends
a median of 39.9 kb, shorter than a single IP₃ receptor gene.

That is not a failure of the caller. Binned by the number of keys available,
its accuracy on its own calibration set is 100 % at every key count at which
it acts.

![](figures/synteny_reach.png)

**{fig:synteny_reach}.** Why synteny could not answer, showing the joint
distribution of what a trace region has to work with, with the caller's floor
drawn, beside the caller's own accuracy on the same axis. Putting reach and
accuracy on one axis is what makes "accurate and unavailable" a readable
sentence rather than an excuse.

On the 8 regions it can reach the caller returns a call for 6 and agrees with
the alignment's own attribution on all 6. That is an independent instrument
corroborating the attribution, on six regions, stated for what it is.

## 9.5 No cell in 927 reaches the absence state

Of the 927 cells, 783 are a single locus at full coverage, 90 a locus
truncated by the assembly, 7 a partial locus with no assembly excuse, 43
reassembled across contigs, and 4 family loci the panel cannot file. Zero fall
in each of the three states that would license a loss claim or invalidate a
genome.

![](figures/character_matrix.png)

**{fig:character_matrix}.** Every genome by paralogue cell, ordered by
assembly contiguity with the bar drawn, and the gene-equivalents each assembly
holds. The panel exists so a reader can see that the red the genome ledger
showed is gone, and see where it went.

The 43 reassembled cells are the substantive change this chapter makes to the
ledger, because those were the cells a naive count would have read as
candidate absences. Every one of them holds a gene.

The four unfilable cells are ITPR2 and ITPR3 in the two cyclostomes, which are
exactly the four the reconciliation flagged, recovered here by a rule that
reads the sweep's own per-genome locus counts rather than that chapter's
table.

**Reference coverage cannot count copies, and genomic sequence can.** The
paralogues are similar enough that a genome holding only ITPR1 recovers most
of the ITPR2 reference, which is why the calibration needed a decoy at all. A
placed locus and a reassembly occupy different places, so the per-cell
contributions add.

**In every one of the 189 assemblies contiguous enough to carry this gene, all
three paralogues are there**, at a median of 3.00 gene-equivalents and a
minimum of 3.00. Below the bar, 11 of 120 fall short, and the shortfall tracks
contig N50 rather than taxonomy.

## 9.6 Lesion counts are mostly an alignment statistic, and the confounder is not the obvious one

The sweep records how many frameshifts and in-frame stops each locus required.
Read naively that is a pseudogene screen. Read honestly it is mostly a
sequencing and alignment statistic, and the family's own control makes the
point: the ryanodine receptors, which are a 5,000-residue gene nobody claims
is dead, have the lowest share of lesion-free loci of the four cells.

Lesions are counted per kilo-aligned-residue rather than per gene, or the
1.8-times-longer ryanodine reference leads every ranking by construction. The
bar is the upper quantile of a population the screen never scores, comprising
loci at full coverage, in contiguous assemblies, whose own annotation names
them as family, which are genes a second pipeline independently calls
functional. 72 % of them carry no lesion at all.

**The confounder that matters is not the one anybody expects.** Measured over
every scored locus, assembly contiguity barely moves the lesion count
(ρ = −0.077) and the locus's identity to its bait moves it a great deal
(ρ = −0.397). A lesion count is substantially a measure of how far the
reference is from the gene, because a poorly matched bait buys alignment with
frameshifts.

![](figures/integrity.png)

**{fig:integrity}.** Lesion density against the two confounders that could
produce it without a gene being dead, and the paired within-genome test that
removes both. The identity panel is drawn first because contiguity is the
confounder everyone expects and identity is the one that turned out to be
real. Density is logarithmic with an explicit zero band, since 72 % of intact
loci carry no lesion and a linear axis puts the whole calibration population
on one pixel.

So the paired within-genome test, in which each cell is compared against the
same genome's other family loci and which removes the assembly entirely, is
run twice, once identity-matched, with ties dropped and counted and the four
tests corrected together by the false-discovery-rate procedure this project
uses throughout [R170].

**One result survives, and it is paralogue-specific.** ITPR3 carries more
disabling lesions than its own genome's identity-matched sibling loci, in 39
genomes to 14. ITPR2's excess in the unmatched test disappears once identity
is matched, so it was an alignment artefact. ITPR1's deficit does not survive
correction, and the ryanodine control shows no excess, so the ITPR3 signal is
not a property of the family's gene structure as a whole.

What that is **not** is a pseudogene finding. All 44 loci above the bar in
contiguous assemblies are at full coverage with an intact gene model, and the
one sound direction here was already established: zero stops falsifies a
pseudogene call, and a handful does not establish one. The verdict column
stays one-sided.

## 9.7 The tree the count is placed on is a taxonomy, declared as an input

A reconciliation needs dates and a hand-curated tree, and a loss count needs
coverage. Chapter 7's 31-species tree cannot serve here, because 278 of the
species whose cells this matrix holds are not in it.

So the topology is NCBI taxonomy over all 309 assemblies, taken from the same
archived dumps the manifest was built from: 309 genomes placed, 468 internal
nodes, **49 of them polytomies**. There is one tip per accession rather than
per species, because the sweep's unit of observation is an assembly.

It is an input rather than a result. Its polytomies are real and are not
resolved, which for parsimony makes a placement less confident and never
wrong, but a count of independent losses under a polytomy is bounded by the
resolution, so the degree distribution is committed. The tree's compatibility
with the curated one is checked rather than assumed, and the check is asked so
that the two trees' different sampling cannot register as a disagreement. The
first version counted assemblies of species the curated tree never sampled as
intruders and reported 21 of 29 clades as unrecovered while every matched node
carried the right name.

## 9.8 Dollo parsimony places no loss, and the gain nodes agree with Chapter 7

Under the operating point, on both codings, the family character and all three
paralogue characters are present in every genome that can be scored, absent in
none, and Dollo places **zero losses**, at both the maximum and the minimum.

The count is zero because no cell reaches the state that would license it,
rather than because the routine cannot return anything else. A constructed
control puts a loss on a known edge and requires it found, another requires
two sister losses merged into their parent, and a third constructs a
contiguous, controlled, spare-free empty cell and requires it to come back as
an absence. **A count that can only ever return zero is not a measurement**,
and these are what stop this one being that.

**The gain nodes are a consistency check nobody asked for.** Dollo places the
single gain at the most recent common ancestor of the tips carrying the
character, and for the paralogue characters that node is not pinned. The three
answers are Vertebrata, Gnathostomata and Gnathostomata, which is exactly
where Chapter 7 placed the two duplications, recovered by a method that reads
no gene tree, no alignment and no reconciliation.

It is worth printing and it is **not independent evidence**, because the
reason ITPR2 and ITPR3 have no cyclostome tip is that those cells are
unfilable for want of a cyclostome-labelled bait, and the reason the
reconciliation placed the duplication there is a reconciliation.

## 9.9 What it would take to manufacture a loss out of this data

With no loss to place, what is worth reporting is the shape of the zero: how
far the rules must move before a loss appears, which rule has to move, and
what each move buys.

Thirty-two settings were walked, comprising an ordered eight-rung evidence
ladder, each rung named after what it refuses, crossed with the two protective
rules on and off, on both codings, under three branch-length schemes. The
first three rungs of the ladder are the calibration's own measured gap edges
and midpoint rather than round numbers, so the first sensitivity question
asked is whether moving the bar across its own uncertainty changes anything.

![](figures/sensitivity_matrix.png)

**{fig:sensitivity_matrix}.** The loss count in every cell of the grid, for
the family-level coding and the paralogue-resolved one, with the operating
point marked. A zero is drawn as an explicit zero and never as an empty cell,
because an empty cell reads as "not measured" and the zero is the result.

**Moving the reconstruction bar across the whole gap the calibration measured
manufactures no loss on either coding.** The bar's position inside its own
uncertainty is not what any result here rests on.

Read across the two codings, the result is an asymmetry rather than a number.
The **family-level coding manufactures a loss in 2 of 32 settings**, and its
worst case is one genome in 309, reached only by refusing everything except a
complete locus and ignoring the contiguity bar at the same time. The
**paralogue-resolved coding manufactures one in 18 of 32**, up to 45 loss
edges. A per-paralogue absence is fragile to every one of these knobs and a
family-level absence is not.

![](figures/reconstruction_bar.png)

**{fig:reconstruction_bar}.** The calibrated bar with the gap's two edges
drawn. The within-row offset of each point is its rank in its own row, using
no hash and no random number generator, so the figure is reproducible.

**The branch-length axis changes nothing, and that is reported rather than
omitted.** Across all 32 settings and all three schemes, the number of
settings at which a branch-length scheme changes a parsimony count is zero.
Parsimony counts edges and does not read a length, so this could not have come
out any other way, but the brief names branch lengths as an axis and a matrix
that quietly dropped one axis would be indistinguishable from one that had
tested it.

## 9.10 The rate models are refused, and the refusal is measured

The brief asks for equal-rates, all-rates-different and irreversible Markov
models [R176]. On the primary character all three are **refused**, and the
refusal is measured rather than asserted.

An invariant character contains no transition to estimate. Profiling each
likelihood along a rate grid over eight orders of magnitude gives, at every one
of the twelve model by branch-length by axis combinations, a monotone curve
with its maximum on the grid's boundary. The loss axis falls to the floor and
the all-rates-different model's gain axis rises to the ceiling, which is what
unidentifiability looks like when it is drawn rather than argued.

![](figures/mk_profile.png)

**{fig:mk_profile}.** The likelihood along a rate grid for every model, axis
and branch-length scheme. The all-rates-different model is profiled on its
gain axis rather than on its diagonal, because the diagonal is the equal-rates
model by construction and would put the same curve on the figure twice under
two names.

So no rate is reported for the primary character, because a fitter run on it
would return its own starting point.

Where the sensitivity matrix does produce variation the fits are real, and
they are reported for what they are: a property of the filter rather than of
the family. In the most extreme cell of the matrix the three branch-length
schemes give rates spanning a factor of 495 while describing the same
character, which is the branch-length axis doing the only thing it can do
here, namely setting the units a rate is quoted in.

## 9.11 The fossil test has no dead loci to run on, and it hands on one lead

If a paralogue died in an ancestor, its descendants should share lesions more
often than independent decay would give. That test needs dead loci.

Of 1,760 scored loci, 44 clear the measured lesion bar, and **all 44 are at
full coverage, delivering a complete gene model.** The screen was made
deliberately generous, counting a locus as a candidate if any of three
readings fires, because a test made hard to pass would make the zero
uninformative. Even so, seven family loci in 874 fire any reading at all, and
all seven do so on one or two internal stops in a complete model.

One or two stops in a 2,700-residue model at full coverage is not the
signature of decay. **The shared-lesion test has no dead loci to run on, and
that is reported with its denominator** rather than omitted, because an
omitted section is indistinguishable from a section nobody ran.

**The one lead this chapter hands on is a bird result.** The ITPR3 lesion
excess of §9.6, stratified by vertebrate class and corrected across the
strata, is 25 genomes to 2 in birds and 7 to 6 in ray-finned fish, which is no
signal at all on a sample that is smaller but not small.

![](figures/lesion_strata.png)

**{fig:lesion_strata}.** The identity-matched sign test stratified by class,
the strongest stratum split by the contiguity bar, and the fossil denominator.

The design makes a cell and its siblings the same observation twice, because
the comparison is within one genome against that genome's other family loci,
so a bird ITPR3 that is elevated makes its own ITPR1 and ITPR2 look deficient
by construction. The three bird rows are one result rather than three.

**The control has to be printed.** Twenty-one of the 27 bird pairs are in
assemblies below the contiguity bar, and that is where the test has its power.
Above the bar all six pairs point the same way and none points against, but
six pairs cannot carry a test. Two thirds of bird assemblies are below the
bar, the worst of any class, so this is precisely the class where an indel
signal is hardest to separate from an assembly signal. The identity control is
clean, so the confounder that mattered in §9.6 is not what this is. **The
lineage is now named and the mechanism is not, and the contiguous half of the
evidence is too small to settle it.**

## 9.12 What the zero means, and the two ways to escape it

Every vertebrate genome in a declared scope of 309, searched with an
instrument whose false-negative rate is measured at 15.2 % overall and 0.9 %
in contiguous assemblies, carries all three IP₃ receptors. No lineage has lost
one. The count is robust to moving every threshold across its measured
uncertainty, and the settings that would manufacture a loss are named and
counted rather than avoided.

That is a strong claim, and the strongest form in which it can be made is the
one that names its own escape routes. There are two. A per-paralogue absence
is fragile to the analytical settings in a way a family-level absence is not,
and a reader who wants a loss can have one by refusing all reconstruction
evidence and ignoring assembly contiguity, which is a description of the
setting rather than a defence of it. And the scope is 309 vertebrate genomes,
so nothing here extends to a lineage with no assembly.
