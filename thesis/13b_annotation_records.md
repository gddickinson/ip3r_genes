## 13.8 Two annotation failures were validated down to the exon

An audit that counts states is a survey. Proving that a particular annotation
is wrong about a particular gene needs molecular evidence, and this section
takes two cases down to the exon.

**Which two cases, and why those.** The measurement is annotation loss,
meaning the share of a recovered gene's coding footprint that no single
annotated model delivers, read straight off the sweep's own coverage rather
than re-derived. Five eligibility rules apply, each a positive test, and the
fifth is the one that does the work: at least ten annotated protein-coding
genes elsewhere in the same genome must be longer than the locus.

Without it, the top of the ranking is loci in assemblies with a genome-wide
gene-length ceiling, and validating one would report a property of the whole
gene set as a bug at this gene. Measured over the sweep it removes six loci in
two genomes, one whose longest annotated gene anywhere is 42 kb and one whose
longest is 138 kb against a 726 kb locus, and nothing else.

Three failure modes are labelled and **two are selectable, one case each**,
because taking the top two of a single ranking gives two omissions and leaves
the fragmentation claim unvalidated. Ties are broken on **control strength**,
meaning how many other family loci in the same genome the annotation gets
right.

![](figures/annotation_loss.png)

**{fig:annotation_loss}.** Annotation loss across every recovered locus: 880
loci, 382 eligible, and 360 of them at zero loss. The full distribution is
drawn because it is the denominator the two validated cases are chosen from,
and a case study without its denominator is an anecdote. The failures are a
tail rather than a norm, so the two cases taken forward are the extreme of a
measured distribution rather than the two examples that happened to be
noticed.

**Case A is an omission.** It is a 52.5 kb locus in a chromosome-level fish
assembly whose contigs are 84 times the length of the gene. The gene has 56
coding exons totalling 8,021 bp, and **the annotation places zero gene models
over them.** Not a partial model and not a mis-named one: the annotation has
nothing there, while working on both sides of the gap, with named genes 4.2 kb
downstream and 5.5 kb upstream.

**Case B is a fragmentation.** It is a 68 kb locus in another chromosome-level
fish assembly, with 60 coding exons totalling 8,064 bp, over which the
annotation places **three gene models** contributing 38 disjoint coding blocks
and covering 58 % of the coding footprint. 3,400 bp in 28 blocks has no coding
model at all and reaches no protein.

![](figures/exon_tracks.png)

**{fig:exon_tracks}.** The two loci at true genomic width, with exons never
widened to be visible. A 56-exon gene over 52 kb averages 143 bp an exon, so
fattening the exons would draw a gene whose coding fraction looks like 40 %
when it is 15 %. That geometry is the point: what the annotation missed is a
few per cent of the locus distributed over dozens of small exons across tens
of kilobases, which is why the gene is hard to annotate and why the failure
is systematic rather than careless.

## 13.9 What the DNA says at those two loci, and five checks on it

**Case A's locus splices into 2,673 codons with zero internal stops** and a
terminal stop present. Under neutral drift to its observed divergence, 14.2
stops would be expected, computed from the locus's own codon usage. **Case B
gives 2,687 codons, zero internal stops, and 11.4 expected.**

That comparison is the right shape for the claim being made. The annotation's
implicit assertion, that this is not a protein-coding gene, is falsifiable,
and zero nonsense codons over 2,700 falsifies it. The rule is deliberately
one-sided: zero stops falsifies a pseudogene call, and a handful would not
establish one.

Five independent checks follow, and they are stated in the order a reader
needs them rather than in the order that flatters the conclusion.

**Splice sites.** Every intron in both loci carries a canonical splice pair.
An alignment is not a gene, and an exon structure whose introns are spliceable
is evidence that this one is.

**Exon boundaries against other genomes' annotations.** Thirty-five and 40
swept genomes carry the same paralogue, aligned from the same bait, with their
own annotation independently placing one model over the whole alignment.
94.5 % and 100 % of each locus's internal boundaries are shared by a majority
of them. This validates the instrument at this gene rather than the individual
locus, which is why it is reported alongside the reading frame and the splice
sites rather than instead of them.

**Transcript evidence at the junction was searched for, and the denominator
is the result.** Probes across every junction the annotation does not model
were searched against every transcript record for the species, and returned
nothing. The species has 10 transcript records in total. A search of that
many returning nothing has not shown the gene is untranscribed. It has shown
the species has almost no deposits.

**That criterion has its own negative control**, and it needed a second half.
Run against the locus's own genomic DNA, where a probe of contiguous spliced
sequence cannot span its junction by construction, the probes make 50 hits and
none spans. So the zero above is a discriminating zero rather than a test that
never fires. On its first run that control returned six false spans, which is
how the rule acquired its second half: an anchor requirement alone admits an
alignment that has run a dozen bases past the junction into the intron.

**The rest of the family in the same genome is intact.** Every other family
locus in both assemblies also splices into an uninterrupted reading frame.
Whatever the annotation is doing here, it is not responding to a damaged
gene.

**The neighbourhood check is available for one case and not the other.**
Case A's missing gene sits between the two neighbours its paralogue carries
in 197 and 183 of 215 swept vertebrates. For Case B the
check is not available at all, because that annotation names its genes by
locus tag and there is nothing to match, and that is reported as a limitation
of the comparison rather than as a negative result.

![](figures/fragment_tiling.png)

**{fig:fragment_tiling}.** The annotated proteins tiled back onto the
genome's own recovered loci, which is what says whether the annotation
missed a gene or broke it into pieces. The proteins are translated from the
assembly's own annotation and genome rather than downloaded, because a locus
filed as a pseudogene emits no protein record and those loci are exactly the
ones under dispute. The guard that matters is that the subject set is the
genome's own loci: a query whose own locus is missing from it lands on its
nearest paralogue instead, and the row then reads as a confident naming
disagreement, which is exactly what happened before the check existed.

![](figures/case_validation.png)

**{fig:case_validation}.** The five checks on each case, together. The
importance of showing the five checks together is that no single one of them
is decisive. Splice dinucleotides, an intact reading frame, junction probes,
exon-boundary concordance and the neighbourhood consensus each fail in
different circumstances, and they do not fail together. A reader asking
whether a recovered gene is an alignment artefact is asking whether every
one of these could be wrong at once, which is what this panel is arranged to
answer.

**One number frames both cases.** One of the two assemblies files 7,380 of its
23,345 genes as pseudogenes, which is 31.6 % of the gene set. Its 14
family-named gene models include 8 filed as pseudogenes and 6 emitting a
protein.

## 13.10 Reads show that all seven lost genes are transcribed and spliced

The junction question the case studies could not answer needs reads rather
than deposits, and this is where it is answered.

**The scope is derived from the ranking rather than chosen.** Four rules, each
a positive test, apply. Every locus the eligibility rules admit whose
annotation loss clears a floor is included, read straight off the committed
ranking rather than re-derived. The species must have enough public runs,
which is a floor on askability rather than on power. And every other family
locus in the same genome joins as an internal control, including the ryanodine
locus, whose reads must land on the ryanodine gene.

A control locus carries the loss measured for it or none at all, in that the
field is empty rather than zero, because a locus never ranked and a locus
measured at zero loss must not share a cell.

Reads come from the public sequence archive [R163] and are mapped with a
splice-aware aligner [R165], run in unspliced mode for the reason given below.

**The reference every read is mapped against is built rather than
downloaded.** It is the spliced genomic coding sequence, deliberately not the
frameshift-corrected version, because correcting a frameshift inserts one or
two bases that every read crossing it would carry as an indel, costing reads
at exactly the difficult sites.

Validation is **colinear block placement**, which has no tuned threshold. Each
block is translated in its own phase and searched verbatim in the aligner's
own protein from the previous block's match onwards, and a block may fail to
place only if the model's own frameshift count explains it. The identity test
it replaced was the wrong instrument, and measuring it showed why: a correct
reference scores 1.00 with no frameshifts and 0.92 with nine.

Three housekeeping anchors are aligned **by the same aligner and spliced by
the same module as the targets** rather than pulled from each annotation,
because two of the three annotations name no genes at all, one filing all
29,240 of its models as hypothetical proteins. Building anchor and target
alike means a difference between them cannot be a difference in construction.
The anchor accessions are **resolved by query against a declared length band**
rather than remembered, because the first version hard-coded three and got all
three wrong, including a 427-residue protein named for a 335-residue enzyme.

**Seven of seven loci that the annotation loses are transcribed and spliced**,
meeting all three detection criteria in at least one library, against their
own reversed decoys at zero.

![](figures/s12_detection.png)

**{fig:s12_detection}.** Detection per locus against its own reversed decoy,
which is the spurious-mapping floor of this particular reference rather than
a threshold. The decoy replaces a chosen cut with a measured one: a sequence
of identical length and composition and no homology collects whatever this
reference collects by accident, so a locus above its own decoy is detected
on evidence. Counts are on logarithmic axes because these libraries differ
roughly forty-fold in depth, and the comparison that means anything is
between a locus and its decoy in the same run.

**298 of 314 junctions that no annotated model spans are crossed by reads**,
which is 94.9 %, against 95.2 % of the annotated junctions in the same genes.

![](figures/s12_junctions.png)

**{fig:s12_junctions}.** Every junction of every reference, crossed or not,
scored per junction rather than per gene, with the decoy floor drawn rather
than stated. That unit is what makes the read evidence answer the annotation
question rather than a weaker one: that a gene is transcribed says little
when part of it is already annotated, whereas the particular junctions no
model contains being crossed by reads is direct evidence for the exon
structure the annotation is missing.

**The confound is printed beside the ratio.** Several of these libraries are
strongly 3′-biased, and in a truncated or fragmented gene the annotated
junctions are the ones at the 5′ end, where coverage is thinnest. Comparing
the two recovery rates without saying so would report a property of library
preparation as a property of the annotation. **The claim this makes is the
absolute one**, that these particular junctions are crossed, rather than the
ratio.

![](figures/s12_gap_coverage.png)

**{fig:s12_gap_coverage}.** Read coverage over the sequence the annotation
loses, which is the read-level form of the loss measurement. Two independent
instruments, one aligning a protein to a genome and one mapping reads to a
spliced reference, are asked about the same missing sequence, and their
agreement is what rules out the loss being an artefact of either.

## 13.11 Four controls make the read evidence readable

**A spurious-mapping floor is measured per sequence.** Every reference carries
a reversed decoy of identical length and composition with no homology, so what
it collects is this reference's own floor rather than a number from a paper.

**Cross-mapping between paralogues is measured rather than argued away.** The
reference is a closed set, so if the paralogues were similar enough at
nucleotide level for a read to cross between them, the whole exercise would be
measuring the wrong gene, and these paralogues are more similar to each other
than the family this method was ported from, whose version dismissed the risk
in a caveat. Every reference was tiled exhaustively with synthetic reads at
fixed steps and mapped back under the same settings. **Zero of 12,500 reads
are assigned elsewhere**, and the result is committed whether it is zero or
not. It is tiled rather than sampled, so it is exhaustive over positions and
deterministic.

**Spliced alignment is switched off deliberately.** The reference is already
coding sequence, so a junction read is contiguous here, and leaving spliced
mode on would let the aligner open a gap inside the coding sequence and call
an intron that does not exist.

**Junctions with no reads are written as zeros rather than omitted**, because
the denominator is half the result.

![](figures/s12_instruments.png)

**{fig:s12_instruments}.** The deposit cross-check with its genomic negative
control. The same junction-probe test the case studies ran is re-asked on a
panel that now includes a species with 37,166 transcript records against the
other two species' 43 and 10, and it still cannot answer. A search returning
nothing in a species with tens of thousands of deposits and nothing in a
species with ten is two very different situations, and the negative control
is what shows the probes themselves work. Reporting the result as
underpowered is what keeps a database's emptiness from being read as a
biological absence.

Eleven constructed negative controls run before anything is written, and were
mutation-tested on five deliberate breakages, all five caught. The most
instructive requires junction classes to be matched **by genomic coordinate
rather than by index**, on a two-intron minus-strand gene, which is the
smallest case that can tell a coordinate lookup from an index lookup, because
a coordinate-sorted and a target-sorted list have the same length and opposite
order on a minus-strand gene. A one-intron test cannot distinguish them, and a
mutation test found exactly that mapping slipping through.

## 13.12 The correction list is 297 items addressed to somebody else

The deliverable of this chapter is addressed to somebody else: **297 correction
items, 52 of them high priority, 18 withheld.**

Each row carries the assembly, the coordinates, the current state, the current
name, the proposal, an archived evidence file a curator can open, and the
numbers the class fired on.

**Priority is a rule rather than an impression.** High priority needs the gene
demonstrably present, the assembly demonstrably able to carry it, and the
reading frame intact. A partial recovery or an unscored reading frame is
medium. Anything below the contiguity bar is low, because there the
annotation's silence may be the assembly's fault rather than the annotator's.

**The integrity veto is applied as a column rather than a filter.** Where the
lesion screen scored a locus as lesion-rich, the audit does not propose
resurrecting it, but the row is still written, flagged, with its reason,
because a locus this audit declined to correct is evidence about the audit.

A protein-record rename additionally needs an absolute score as well as a
relative margin, since the candidates score 90 to 191 bits over roughly 2,700
residues.

## 13.13 What this chapter settles about the record

The family's public record is substantially worse than the family. Three
quarters of the genes this project demonstrated are unreachable from any
protein database, more than half of the full-length protein records that do
exist carry no usable gene symbol, and a quarter of the loci in annotated
assemblies are not delivered as one complete model.

**And it is not this family's fault.** The sister family, in the same
assemblies through the same pipelines, fails at the same rate. What this
chapter measures is how vertebrate genes of this size and structure are
recorded, for which this family is an unusually well-controlled probe, because
it comes with a size-matched, domain-sharing sister to compare against.

Two of the failures are proved at the exon and confirmed by reads. Both are
genes with intact reading frames, canonical splice sites, boundaries
corroborated by dozens of other genomes' independent annotations, and
transcript evidence crossing the junctions no model spans, in chromosome-level
assemblies whose contigs are dozens of times the length of the gene.

Two thirds of the apparent annotation problem, however, is assembly quality,
and every claim here is made with the contiguity control beside it.
