### The family is not badly recorded; the archive is

Every one of the 2,144 gene-scale loci was scored against the assembly's own
annotation — same strand only, on coding blocks and never on gene spans —
as complete, split, fragmentary, held only by a non-coding feature, or
unannotated, with the ryanodine receptors scored the same way in the same
assemblies through the same pipelines as the control.

**The audit's own premise is contradicted by its control.** IP₃ receptor
loci are 73.9 % complete and 26.1 % failing; the ryanodine receptor control
is 77.9 % and 22.1 %, and no overall difference survives multiple-testing
correction (q = 0.13 over all scorable loci, q = 0.82 above the contiguity
bar). One state does separate, and it is the family-specific one: an IP₃
receptor locus is 2.7 times more likely than a ryanodine receptor locus to
be held *only* by a non-coding feature (42 against 16, q = 0.006) — a gene
the annotation identifies and names correctly, files as a non-coding
biotype, and therefore serves no protein for. Above the contiguity bar that
excess is 4 against 5 (q = 1.0), so it is confined to assemblies too broken
to carry the gene.

**What decides whether the gene is recorded is the archive, not the gene**
(Fig. 7). RefSeq annotations are 98.8 % complete against submitter-deposited
GenBank annotations at 37.5 %, with every state separating at q < 10⁻³⁰⁰.
Controlling for contiguity narrows the gap without closing it: 99.4 %
against 63.8 %. Assembly contiguity is the second-largest effect — the IP₃
receptor failure rate falls from 26.1 % to 6.7 % across the bar. The
completeness bar itself was inherited and validated rather than re-derived:
over 1,077 correctly named, fully recovered loci a single annotated model
covers a median 0.993 of the alignment, so the 0.50 bar sits at that
distribution's 1.3 % point (Extended Data Fig. 15).

**A "fragment" in the archive is an annotation that stopped where nothing
splices.** The 264 fragmentary and 27 split loci are a statement about
coverage and say nothing about *where* the annotation stopped, and the
difference is the whole claim: a model ending at a genuine junction has
produced a plausible short gene, one ending inside an exon has produced a
boundary no splicing machinery could make. Scoring each model's own internal
termini against the alignment's exon boundaries and against the two genomic
bases immediately outside them, **228 of the 291 loci carry at least one
terminus inside an exon** and **3 are broken entirely at junctions the gene
has** (Extended Data Fig. 9). Of 1,448 internal termini, 406 sit inside an
exon and 353 inside an intron, against 689 (47.6 %) on a boundary the gene
model carries. The verdict is one-sided by design — a mid-exon terminus
falsifies the reading that the annotation found a real gene end, while a
terminus on a boundary does not prove the pieces are separate genes — and on
that one-sided reading almost every fragmentary record in this family is an
annotation failure rather than a gene boundary.

**On the protein side the family call is not in dispute; findability is.**
Of 11,402 full-length family records searched against the labelled bait
panel, 4 disagree with the census on family by sequence and 5 carry a name
from the wrong family, all non-vertebrate and all under 200 bits; 52 of
8,306 vertebrate symbols name the wrong paralogue. But **3,872 records carry
a placeholder symbol and 2,395 carry none at all — 55.0 % of full-length
family records have no usable gene symbol** — and 66 more are named for the
superfamily rather than the family. All 15 reference proteomes that returned
no family protein at all are gene-caller failures: every one of those
species has a genome in the scope and every one of those genomes carries the
gene.

**The consequence, counted.** Of 1,232 genes demonstrated to exist in these
assemblies, **940 (76.3 %) cannot be reached by any protein-database
search** — 386 in species with no reference proteome, 248 where only
fragmentary records exist, and 286 where records exist but none resolves to
that paralogue. This is not an artefact of the margin species the scope was
extended for: 74.3 % for order representatives against 78.5 % for margin
species.

**Two failures validated to the exon, and both genes are transcribed.**
Where the annotation demonstrably could have delivered the gene — 382 loci
passing five eligibility rules, of which the decisive one asks whether the
same annotation builds genes that long anywhere else in the same genome —
it gets it right at **375 of 382 (98.2 %)**, with a median annotation loss
of zero and 360 loci at exactly zero. The seven failures fall in three
genomes and cover all three failure modes. We validated one omission and one
fragmentation at nucleotide resolution (Extended Data Fig. 16).

*Nibea albiflora ITPR2* has 56 coding exons and **not one annotated gene
model over any of them**, in a chromosome-level assembly with 84-fold
headroom whose annotation gets *ITPR1* and *ITPR3* right. The locus splices
into 2,673 codons with **zero internal stops**, against 14.2 expected under
neutral drift at the observed divergence; all 55 introns carry a spliceable dinucleotide (54 GT-AG, 1 minor);
52 of 55 exon boundaries are shared by a majority of 35 independently
annotated genomes; and it sits between `SSPN` and `BHLHE41`, this
paralogue's two most conserved neighbours across 197 and 183 of 215 swept
vertebrates. The same genome files every complete family gene as a
pseudogene and names its *ITPR3* locus *ITPR2*. *Dissostichus eleginoides
ITPR3* is one 68 kb gene called as three protein-coding models tiling
residues 1–67, 52–467 and 468–1593, with the 3′ 41 % unmodelled; all 59 introns are
GT-AG and 59 of 59 boundaries are shared by a majority of 40
genomes. Both species have **zero IP₃ receptor protein records in any
database and three complete IP₃ receptor genes in their DNA.**

Reads settle it. 67 public RNA-seq runs (536 million reads) across three
species were aligned against a per-species reference of spliced genomic
coding sequence with a reversed, composition-matched decoy for every
sequence. Each failed locus is detected in 13–31 of its runs and in 6–9
tissues, and **298 of the 314 junctions (94.9 %) that no annotated model
spans are crossed by reads that read through them**, against 79 of 83
(95.2 %) of the annotated junctions in the same genes. The closed-set risk
was measured rather than argued: every reference tiled exhaustively with
synthetic reads gives 0 of 12,500 cross-mapped.

The audit ends in a list addressed to somebody else: **297 corrections** —
assembly, coordinates, current state, current name, proposal and archived
evidence — of which 52 are high priority and 18 are recorded as withheld,
with the reason, rather than dropped.

