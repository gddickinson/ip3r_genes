# 8. The gene: exons, introns, and what a fragment is

## 8.1 The question a failed background claim left

Chapter 2 struck a claim: that each IP₃ receptor is a 58-to-60-exon gene
spanning hundreds of kilobases. The exon counts turned out to be 62, 57 and
58, and the spans 354, 498 and **76** kb. Three paralogues encoding proteins
that differ in length by 3 % occupy genomic spans that differ by a factor of
six and a half.

That leaves a question with two possible answers. Either the exon structure
differs between the paralogues, or the introns do — and only one of those is a
statement about how the gene has evolved. Answering it across 309 genomes also
answers a second question the annotation audit needs: when a database calls a
locus *fragmentary*, is that a real gene boundary or an annotation that
stopped somewhere nothing splices?

## 8.2 Reading exons out of a genomic alignment

The measurement is made on the sweep's own alignments rather than on
annotations, because Chapter 4 found that most of these genes have no usable
annotation. Four things a naive reader of that output gets wrong, each with a
constructed control that fires on it.

**The model keys.** The aligner numbers its alignments from one per run, so
the two assemblies large enough to be processed in chunks repeat every
identifier. A reader that keys models the way the aligner names them merges
the blocks of two different genes, and the merged gene has more exons and
longer introns than either.

**Gene order.** Blocks come in target order, and an intron computed as the
next block's start minus the previous block's end reverses every minus-strand
gene.

**Splice dinucleotides.** Reading an intron's two ends without reverse-
complementing both turns every canonical minus-strand intron into a
non-canonical one — a systematic signal that would discredit exactly the
alignments this chapter is defending.

**A block is not an exon.** The aligner emits blocks; some adjacent pairs are
separated by a gap too small to be an intron.

**Frameshifts turn out to be inside blocks, not between them**, which is a
measurement rather than an assumption, and it means the merge this chapter
needs is a sub-intron-gap merge rather than a frameshift merge.

## 8.3 What an intron is, measured against the genome

Rather than choose a minimum intron length, every block pair was binned by gap
size and scored by its splice dinucleotides — evidence the alignment score did
not produce.

The calibration then **refuses to hand out the threshold it was written to
produce.** The sub-intron population is empty: there is no population of small
gaps that fail the splice test. So the floor goes at the smallest gap the
genome calls spliceable, 10 bp, rather than at a declared 30 bp fallback,
which would have merged the ten junctions between 10 and 29 bp that read as
canonical splice pairs.

## 8.4 The instrument's own error rate

Before measuring exon structure with an aligner, it is worth knowing how often
the aligner puts a boundary where no splice site is. **99.89 % of 112,254
junctions read as a canonical splice pair.**

![](figures/junction_quality.png)

**{fig:junction_quality}.** The aligner's junction error rate, drawn on the
**complement**: at 99.9 % agreement a bar of the agreement is four full bars
and shows nothing. Each gap bin carries the count it rests on.

That is corroborated by a second, independent pipeline: the sweep's boundaries
against **188,146 annotated coding-sequence edges from 164 genomes**, at
94.5 % exact agreement. An earlier chapter made the same comparison for two
loci against 35 and 40 genomes; this is that check generalised.

## 8.5 One coordinate frame for three genes

An exon boundary is a position in the *bait's* numbering, and the sweep used
38 baits, so comparing an ITPR1 gene's boundaries with an ITPR2 gene's needs a
frame both can reach. That frame is Chapter 6's alignment: all three human
paralogues and a human ryanodine receptor are tips of it, so a boundary
travels from the bait to its cell's human reference and then to an alignment
column.

Ninety-eight frames were built and none refused.

**The frame is checked, not assumed.** All fourteen residues measured on the
structure — the ten IP₃ contacts, the two selectivity-filter residues and the
two gate residues, each in its own paralogue's numbering — land in the **same
alignment column in all three paralogues**. A frame that had slipped anywhere
in the pore or the ligand core would fail there, and a slipped frame is exactly
the error that makes a shared-intron count look like a result. A constructed
control shifts one paralogue's aligned row by a single column and requires the
test to refuse.

Twenty-one constructed controls run before anything is written, and the suite
was mutation-tested on six deliberate rule breakages, all six caught. Four of
this chapter's rules return a *good-looking* number when they are wrong, which
is why the suite is as large as it is: the model-key collision, the
strand branch, an inverted intron-phase convention that swaps two phases
everywhere without changing a single count, and the duplication detector
discussed in §8.9.

## 8.6 The architecture: the count is conserved, the packaging is not

Across 1,378 genes in 189 vertebrate genomes:

The three paralogues sit at **57 to 58 coding exons** with median genomic
spans differing by **4.2-fold**, and coding sequence lengths of 8,250, 8,100
and 8,008 bp. What differs between them is how much intron is wrapped around
nearly the same protein.

![](figures/architecture_by_paralog.png)

**{fig:architecture_by_paralog}.** Exon count and genomic span on **separate
axes**, span logarithmic and count not. The result is that one is conserved
and the other is not, and a shared scale would hide it.

The ryanodine receptors, measured through the identical instrument in the same
assemblies, carry **104 exons over 14,910 bp** of coding sequence — nearly
twice the gene in both dimensions, at almost exactly the same mean exon length,
144 bp against 138 to 142. The control tells us the instrument is measuring
gene structure rather than a property of this family's alignments.

Every count is recomputed at six coverage bars, because a count that survives
only one bar is a claim about the bar.

**The comparison is paired inside one genome.** Intron size, assembly quality
and annotation completeness all scale with the assembly, so every
cross-paralogue comparison uses one locus per genome per gene, sign-tested with
ties dropped and counted, and the whole family of tests corrected together.

Two things fall out. **The exon counts differ by one to three exons and the
spans by tens of kilobases, in the same genomes.** And **the ordering of the
spans is consistent gene by gene** rather than an artefact of averaging: ITPR3
is the shorter gene than ITPR1 in 161 of the 182 genomes carrying both.

## 8.7 The same introns, in the same places

An intron's position here is a pair: the alignment column its upstream exon
ends in, and its phase. Both halves are needed. Two paralogues can carry an
intron between the same two residues in different frames, which is not one
ancestral intron, and a column alone would call any two introns in the same
region shared.

**Within a paralogue**, each carries a core of about 48 to 49 intron positions
present in at least 90 % of the genomes that have the gene, out of 175 to 185
positions seen anywhere. The ryanodine receptors have their own core of 91, on
the same scale relative to their 101 introns per gene. So the architecture is
not merely a count that happens to be stable. It is the *same introns*, in the
same places, across the vertebrates.

![](figures/intron_positions.png)

**{fig:intron_positions}.** Intron positions by alignment column and phase,
per paralogue. The columns line up.

**Between paralogues** is the test this chapter exists for. Two genes with
about 58 introns each spread over about 2,700 aligned residues will share some
positions by chance, so every comparison is paired within a genome and scored
against a null in which each of one gene's introns is placed independently and
uniformly on the columns *both* loci have residues in, keeping its own phase.
That makes the match count a Poisson-binomial whose upper tail is exact — no
random number generator and no replicate count — and a seeded permutation
without replacement is run beside it as a cross-check.

**The three paralogues share 46 to 49 of their roughly 58 intron positions in
every genome that carries them**, an enrichment of 85 to 88-fold, significant
in every one of 181 to 183 genomes.

**The ryanodine receptors share one.** An enrichment of 2.1 to 2.3-fold,
significant in none of 188 genomes.

That is the sharpest single result in this chapter and it needs the control to
be readable. The ryanodine receptors carry every diagnostic domain of this
family; they are the same superfamily; they have been in every search. And
their intron positions have nothing to do with the IP₃ receptors', while the
three IP₃ receptors' have almost everything to do with each other. The three
paralogues did not converge on a shared architecture — they inherited one,
from a common ancestor much younger than the split from the ryanodine
receptors.

## 8.8 Where a fragmentary annotation stops

Chapter 13's audit calls 264 loci *fragmentary* and 27 *split*: coding
sequence is present, and no single model covers the gene. That is a statement
about coverage and says nothing about **where** the annotation stopped, and
the difference is the whole claim. A model that stops at a genuine junction has
produced a plausible short gene. One that stops in the middle of an exon has
produced a boundary no splicing machinery could make.

So the unit is the annotated model's own **terminus**, not its internal exon
edges — a model's internal boundaries are its own splice sites and are
canonical by construction, so scoring them would answer §8.4 a second time.
Each terminus is scored on two axes that know nothing about each other:
against the alignment's exon boundaries, and against the two genomic bases
immediately outside it read in gene orientation. A terminus coinciding with
the gene's own end is excluded by rule, because a real gene legitimately
starts and stops inside an exon.

The question is asked of all 291 such loci rather than only those in the
architecture scope, because those loci sit in the poorer assemblies by
construction and restricting the question would have answered it on 50 of them
and dropped the hardest.

Of 1,448 internal termini, 406 sit inside an exon of the gene model and 353
inside one of its introns; 689 land on a boundary the model has. **228 of the
291 loci carry at least one mid-exon terminus**, which falsifies the reading
that the annotation stopped at a real gene boundary. Three are broken entirely
at junctions the gene has.

The verdict is one-sided by design. A mid-exon terminus falsifies. A terminus
landing on a boundary does *not* prove the pieces are separate genes, so that
verdict is named for what it is and claims nothing more.

## 8.9 Is any of this two genes?

An exon count is only a gene's if the locus is one gene. The aligner aligns
each bait independently, so a genome encoding the same part of a protein twice
has the *same* bait aligning twice at two disjoint places.

Run naively, that geometry measures **paralogy** rather than duplication: the
three receptors are 61 to 68 % identical, every bait aligns at all three
genes, and the detector's first version reached a specificity of **0.16**. The
pair has to be inside the cell's own loci, with the sweep's clustering and
family attribution imported unchanged.

Scored as a classifier of a copy count it never sees — decided by coverage,
identity and aligned length, and by no pairwise geometry at all — the detector
reaches **99.7 % sensitivity and 97.7 % specificity over 1,236 cells.**

The ryanodine control's specificity is lower on purpose: that cell holds three
genes in every tetrapod, so the detector *should* fire there, and it does.

**No locus in the sweep encodes the same part of the protein twice.** The
within-locus class — two alignments of one bait inside one cluster, which
would be an internal partial duplication — fires on zero pairs, and a
constructed control builds such a pair and requires the branch to fire, so the
zero is a measurement rather than an unreachable code path.

![](figures/fragments_and_duplicates.png)

**{fig:fragments_and_duplicates}.** The verdict on every split and fragmentary
locus; where the annotation's internal termini sit relative to the gene model;
and the duplication detector scored against an independent copy call.

## 8.10 What this chapter settles

The three vertebrate paralogues inherited one gene architecture and have kept
it: 57 to 58 coding exons, of which 46 to 49 positions are shared between any
two of them in every genome that carries both, at nearly two orders of
magnitude above a null drawn from the alignment itself. The ryanodine
receptors, the same superfamily, share one.

What varies is intron length, and it varies by a factor of four in the median
and consistently gene by gene. The 6.5-fold span difference that Chapter 2's
failed background claim exposed is entirely intronic.

And a database's *fragmentary* is usually not a gene boundary: 228 of 291 such
loci stop somewhere nothing splices.
