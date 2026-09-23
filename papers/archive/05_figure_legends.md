## Figure legends

**{fig:contribution}.** Most demonstrated genes have no protein record,
and a better protein search would not change that. Panel (a) compares two
protein-level searches of the same reference proteomes, an exhaustive
enumeration of the family's Pfam signatures and a pair of profile hidden
Markov models. Each bar gives the share of the records either search calls
family that only one of them found, split by record length, for each swept
group. The profile adds nothing at gene scale in the vertebrates or the other
animals, so its gain lies in the fragment tail. Panel (b) asks, for
each paralogue cell and for the ryanodine receptor control cell, what
fraction of the genes the genome sweep demonstrates no protein record
resolves to, with 95 % Wilson intervals. The whiskers do not reach below
0.68 for any IP₃ receptor paralogue.

**{fig:by_source}.** Loci in RefSeq annotations are delivered whole far
more often than loci in submitter annotations. Panel (a) stacks
the five locus states for every scorable locus, split by whether the
assembly's gene set is a curated RefSeq annotation or a GenBank annotation
deposited by the submitter. Panel (b) repeats the contrast on the loci
whose contig is long enough to hold the gene. The gap narrows in (b)
without closing, which separates the archive effect from the assembly
effect: GenBank assemblies are less contiguous, yet the GenBank gene sets
still fail on good contigs.

**{fig:family_control}.** The IP₃ receptors are recorded at the rate their
sister family is recorded in the same assemblies. Panel (a) stacks the five
locus states for each of *ITPR1*, *ITPR2* and *ITPR3* beside the ryanodine
receptor control. Panel (b) compares the share of loci with any annotation
failure between the family and the control over all scorable loci and over
loci on contigs able to carry the gene, with the Benjamini–Hochberg
corrected p-value above each pair. Neither pair differs after correction,
so a quarter of loci missing a complete model is the rate for vertebrate
genes of this size, not a defect peculiar to this family.

**{fig:protein_side}.** Protein records that carry a name are almost
always named correctly, and most carry none. In panel (a), each of the
11,402 full-length family records is scored twice, on its gene symbol and
on its protein name, against the paralogue its sequence is assigned to.
The categories run from a correct paralogue, through the family named
without a paralogue and names outside the panel's reach, to a wrong
paralogue, the sister family, both families at once, a locus tag, or
nothing. The pale bands of placeholder and absent symbols fill more than
half the symbol bar. Panel (b) sorts the 15 vertebrate reference proteomes
with no family hit by what an assembly of the same species holds, and
every one falls where the genome holds the gene and the proteome does
not.

**{fig:exon_tracks}.** Two annotation failures are drawn at true genomic
width, each exon at its real size. Panel (a) is the *Nibea albiflora*
*ITPR2* locus, where all 56 exons are coloured as touched by no annotated
gene. Panel (b) is the *Dissostichus eleginoides* *ITPR3* locus, with the
three annotated models drawn beneath the exons they cover and 23 of the 60
exons outside all of them. Flanking named genes mark each end. Exons are
never widened for visibility, because a 56-exon gene over 52 kb averages
143 bp an exon and a widened drawing would misstate how little of the
locus the missed sequence occupies.

**{fig:junctions}.** Reads cross the junctions that no annotated model
spans. Each row is one of the seven IP₃ receptor genes the annotation
loses, and each junction of its spliced coding sequence is a bar at its
position along the gene, with height equal to the number of reads that
read through it. Colour marks whether an annotated model contains the
junction, and a cross on the baseline marks a junction no read crossed, so
the denominator is visible. Transcription alone would say little about a
partly annotated gene; spliced molecules across the missing introns show
that the unmodelled exon structure exists.
