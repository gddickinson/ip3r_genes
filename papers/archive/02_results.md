## Results

### Four in five demonstrated IP₃ receptor genes have no protein record that reaches them

We first asked the reader's question directly. For every genome and every
paralogue cell in which the genomic sweep demonstrates a gene, is there a
protein record, in a reference proteome or elsewhere in the protein
databases, that resolves to that gene? The genome sweep demonstrates 923
IP₃ receptor genes across the three paralogue cells of the 309 genomes.
**744 of them (80.6 %) are reachable by no protein-database search**
({fig:contribution}b). The rate is similar in all three paralogues: 257 of
309 *ITPR1* genes, 260 of 307 *ITPR2* genes and 227 of 307 *ITPR3* genes
have no record that resolves to them. The ryanodine receptor control cell,
asked the same question in the same genomes, is unreachable in 196 of 309.
The control's rate is lower but still well above half.

The 744 are missing for four different reasons, and they have different
remedies. 289 sit in species with no reference proteome at all. 186 sit in
species whose only family records are fragments. In 254 the species has
full-length family records, but none resolves to that paralogue, and in 15
the species has no family record of any kind. Only 179 IP₃ receptor genes
are held both as DNA and as a protein record assignable to their cell. The
deficit does not come from the margin species the scope was extended to
include. Genomes chosen as order representatives lose 413 of 556 genes and
margin species 465 of 592, counting both families.

Changing the search method would not recover them. Inside one database,
the vertebrate reference proteomes, an exhaustive enumeration of the
family's Pfam signatures returns 3,135 gene-scale records and a pair of
profile hidden Markov models returns 3,136, of which 3,135 are shared
({fig:contribution}a). At gene scale the two protein-level searches find
the same set. What they cannot find is a gene that no record holds. The
genome sweep itself contributed 1,058 gene models from 224 genomes to the
census, and of the 488 that are IP₃ receptors, 318 exist only as DNA, 167
lie inside an annotated gene that carries no family name, and 3 lie inside
a gene named for a different paralogue.

### A quarter of annotated loci are not delivered as one complete model

The sweep places 2,144 gene-scale loci in the 309 assemblies: 1,059 IP₃
receptor loci and the ryanodine receptor control loci beside them. 255 loci
sit in assemblies that ship no gene set and are reported as unscorable
rather than as failures. Of the 932 scorable IP₃ receptor loci, **689
(73.9 %) are delivered as one complete annotated model** and 243 are not
({fig:family_control}a). 149 are fragmentary, 14 are split across models
that each deliver a piece, 42 are held only by a feature that emits no
protein, and 38 have no same-strand annotated feature at all.

The completeness bar was inherited from the sweep, which calls a locus
annotated when one model covers half its coding footprint, and was
validated rather than re-derived. Over 1,077 loci whose own annotation
names the correct paralogue on an assembly able to carry the gene, the best
single model covers a median 0.993 of the gene, and only 14 of the 1,077
fall below the bar ({fig:calibration}a). Moving the bar from 0.3 to 0.95
changes the complete count from 712 to 635 of 932 ({fig:calibration}b), so
no conclusion here depends on where it sits.

Assembly contiguity is the largest single determinant. On contigs long
enough to hold the gene, the failure rate falls to 38 of 568 loci (6.7 %).
Two thirds of what looks like an annotation problem is therefore an
assembly too fragmented to carry a 2,700-codon gene. By vertebrate class,
cartilaginous fish and mammals are nearly complete, and birds are the worst
recorded: 243 of 433 bird loci are complete and 124 are fragmentary.

### The family is recorded no worse than its sister family

A 58-exon gene with a sister family sharing its diagnostic domains might be
expected to be annotated badly. The control shows that it is not. In the
same assemblies, through the same pipelines, the ryanodine receptors are
complete at 734 of 942 loci and fail at 208 ({fig:family_control}b). The
difference from the IP₃ receptors does not survive correction for multiple
testing (q = 0.13), and on contigs able to carry the gene the two fail at
38 of 568 and 35 of 632 loci (q = 0.82). How this family is recorded is how
vertebrate genes of this size are recorded.

**One state does separate the families, and it is the one no search can
see.** An IP₃ receptor locus is held only by a non-coding feature in 42
cases against 16 for the ryanodine receptors (q = 0.006). Such a gene is
identified, often named correctly, and filed under a biotype that emits no
protein, so no name- or sequence-based protein search reaches it. The
excess is confined to broken assemblies: on contigs able to carry the gene
the comparison is 4 against 5 (q = 1). The measurement cannot say whether
the excess is a fact about the family or about the assemblies its loci
happen to sit in.

### RefSeq annotations deliver the gene whole far more often than submitter annotations do

Curated and submitter-deposited gene sets are not comparable evidence
[R159]. Across all scorable loci of both families, **1,161 of 1,175 loci in
RefSeq annotations are complete, against 262 of 699 in GenBank annotations
deposited by the submitter** ({fig:by_source}a). Every failure state is
concentrated in the GenBank sets: 263 of their loci are fragmentary and 92
are unannotated, against 1 and 10 in RefSeq.

Submitter assemblies are also less contiguous, so the contrast was run
again on loci whose contig can hold the gene ({fig:by_source}b). It
narrows and does not close: 1,009 of 1,015 RefSeq loci are complete against
118 of 185 GenBank loci. Roughly two fifths of the gap between the archives
is assembly quality, and the rest is the gene set.

### Protein records are named correctly where they exist but half carry no usable symbol

We scored each of the 11,402 full-length family protein records the
databases hold against a labelled panel of 38 bait proteins by blastp. Each
record was assigned to a family only on a stated relative margin, and to a
paralogue only within the winning family. Its own gene symbol and protein
name were then read through the same verdict rule applied to the genome
annotations.

**The family call is not in dispute.** Sequence disagrees with the census
family call on 4 of 11,402 records, and 5 records carry a name from the
sister family. All of these are non-vertebrate records whose best family
score is under 200 bits, where the sequence barely separates the families
either. **Nor is the paralogue call.** 52 of 8,306 vertebrate records carry
a symbol naming a paralogue the panel assigns elsewhere, and a further 74
name a paralogue the panel carries no bait for. Those 74 are reported as
outside the instrument's reach, not as errors.

**What is in dispute is whether a record can be found by name.** 3,872
records carry a placeholder gene symbol and 2,395 carry none, so 55.0 % of
the family's full-length records have no usable gene symbol
({fig:protein_side}a). A further 66 protein names claim the superfamily in
a form that names both families and separates neither. A search by the
family's defining signature does better: PF08709 [R155] is carried by
95.6 % of the 9,841 records both of this project's protein instruments call
family, and by 81.2 % of the 1,554 that only one instrument found.

### Every empty reference proteome is a gene-caller failure

Of 763 vertebrate reference proteomes swept with the family profiles, 15
returned no family protein at all. On its own an empty proteome cannot
distinguish a species that lacks the gene from a gene caller that missed
it. Each was therefore resolved against an assembly of its own species,
searched with the genomic sweep. **All 15 are gene-caller failures**
({fig:protein_side}b). Every one of the 15 species has an assembly in
scope, and in every one the sweep recovers at least one IP₃ receptor locus
at more than half the bait's length while the proteome holds none. The
verdict vocabulary includes an *undecidable* value for a species with no
assembly in scope, and it fired on none of the 15.

### Two annotation failures are intact genes at the exon

A survey of states shows how often the record fails; proving that a
particular annotation is wrong about a particular gene takes molecular
evidence. The cases were chosen by rule. The sweep recovered 880 IP₃
receptor loci, and five eligibility rules admit the 382 that an annotation
could reasonably have been expected to deliver
({fig:cases}a). The decisive rule requires the same annotation to build at
least ten longer genes elsewhere in the same genome. Without it the top of
the ranking consists of loci in gene sets with a genome-wide length
ceiling, and the rule removes 6 loci in 2 genomes and nothing else. At 375
of the 382 eligible loci the annotation delivers the gene, and 360 have an
annotation loss of exactly zero. We took the worst omission and the worst
fragmentation, one per genome, to the exon ({fig:exon_tracks}).

**The *Nibea albiflora* *ITPR2* gene has no annotated model at all.** It is
a 52.5 kb locus with 56 coding exons totalling 8,021 bp, in a
chromosome-level assembly whose contig N50 is 84.24 times the locus span.
The annotation places no gene model over any exon, while annotating named
genes 4.2 kb downstream and 5.5 kb upstream. The locus splices into 2,673
codons with no internal stop, against 14.2 stops expected under neutral
drift to its observed divergence from the bait. All 55 introns carry a
spliceable dinucleotide (54 GT–AG, 1 minor), and 52 of 55 exon boundaries
are shared by a majority of 35 independently annotated genomes carrying
the same paralogue ({fig:case_checks}). The gene sits beside *SSPN*, found next
to this paralogue in 197 of 215 swept vertebrates. The same annotation
files 7,380 of its 23,345 genes as pseudogenes, and it names its *ITPR3*
locus *ITPR2* ({fig:cases}b).

**The *Dissostichus eleginoides* *ITPR3* gene is annotated as three
pieces.** Its 60 coding exons totalling 8,064 bp are covered by 3 gene
models that deliver 4,664 bp between them, leaving 3,400 bp in 28 blocks
with no coding model. The locus gives 2,687 codons with no internal stop,
against 11.4 expected, and all 59 introns are GT–AG. 59 of 59 exon
boundaries are shared by a majority of 40 independently annotated genomes.
Neither species holds any IP₃ receptor protein record, although each
genome carries three complete IP₃ receptor genes.

Transcript deposits could not settle whether these genes are expressed.
Probes across every junction the annotation does not model found nothing,
but the two species hold only 43 and 10 transcript records. A negative
control shows that the probes can fire: against the locus's own genomic
DNA they return 112 and 50 hits, none spanning a junction.

### Reads cross the junctions that no annotated model spans

Reads can answer the question that transcript records cannot
({fig:expression}). The panel is every locus that passes the eligibility
rules with more than half its coding footprint unannotated, in species
with public RNA-seq: all seven of the annotation failures, in three
species, with every other family locus of the same genomes as internal
controls. We aligned 67 public RNA-seq runs, 536,000,000 reads from 16
studies, to a per-species reference of spliced genomic coding sequence.
**All 7 lost loci are transcribed**: each meets all three detection
criteria, including at least two junction-spanning reads and more reads
than its own reversed decoy, in 13 to 31 runs. Read coverage outside
the annotated coding blocks tracks the exon-level loss locus by locus
({fig:gap_coverage}).

The claim that answers the annotation is per junction
({fig:junctions}). Of the 314 junctions in these genes that no annotated
model spans, 298 (94.9 %) are crossed by reads that read through them. The
same genes' annotated junctions are crossed at a similar rate. We do not
rest anything on that ratio, because several libraries are strongly 3′
biased and the annotated junctions of a truncated gene lie at its 5′ end.
The claim is the absolute one: these particular junctions are spliced.

### The audit ends in 297 corrections

The audit's deliverable is addressed to the archives. **297 corrections**
each give the assembly, contig coordinates, current state, current name, a
proposal and an archived evidence file. The corrections comprise 20
unannotated genes, 54 genes held only by a non-coding feature, 217
incomplete models, 1 wrong paralogue name, and 5 protein records named for
the wrong family. Priority is assigned by rule. High priority, 52 items,
needs the gene demonstrably present, the assembly able to carry it and an
intact reading frame. Anything on a contig too short to hold the gene is
low priority, because there the annotation's silence may be the assembly's
fault. 18 items whose loci the reading-frame screen scores as lesion-rich
are kept in the list with the reason, marked as withheld rather than
removed.
