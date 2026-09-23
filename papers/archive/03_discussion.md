## Discussion

The public record of the IP₃ receptor genes is much worse than the genes.
Four in five of the genes the vertebrate genomes demonstrably carry cannot
be reached from a protein database. More than half of the full-length
protein records that do exist carry no usable gene symbol. A quarter of the
loci in annotated assemblies are not delivered as one complete model. None
of this is a statement about the biology of the family, and the paper's
central control is what makes that sayable: the ryanodine receptors, in
the same assemblies and through the same pipelines, fail at a rate
indistinguishable from the IP₃ receptors'. This family is therefore an
unusually well-controlled probe of how vertebrate genes of this size and
structure are recorded, and the answer is general.

Three things decide whether a gene reaches the record, and none of them is
the gene. The first is whether a reference proteome exists for the
species: 289 of the unreachable IP₃ receptor genes are in species without
one. The second is which archive annotated the genome. Curated RefSeq
annotations deliver nearly every locus whole, while submitter-deposited
annotations deliver about a third, and assembly contiguity explains less
than half of that difference. The third is contiguity itself, which removes
two thirds of the apparent annotation failures once loci on contigs too
short to hold the gene are set aside. A study that samples a gene family
from protein databases therefore samples annotation effort. A survey of
this family built that way would have lost most of the birds, whose
assemblies are fragmented and whose gene sets are largely submitter
deposits.

The failures that remain on good assemblies are real failures of
annotation, not of assembly or biology. The two cases taken to the exon are
in chromosome-level assemblies with contigs 60 to 85 times the length of
the gene. Their reading frames are open, their introns canonical, and their
exon boundaries are the ones dozens of other genomes' annotations draw for
the same paralogue. Reads cross almost every junction no model spans. One
of the two annotations files nearly a third of its genes as pseudogenes, so
its silence at a locus that splices into an uninterrupted 2,673-codon
reading frame looks like a pipeline setting, not a judgement about this
gene. The companion analysis of gene architecture reaches the same
conclusion from the other side: most fragmentary models of this family
stop inside an exon rather than at a junction, which is not where a real
gene boundary can lie ({paper:origin}).

The protein-record findings point to a cheap remedy. The databases are not
confusing the families or the paralogues; the sequence and the name
disagree on a handful of records out of 11,402. What they lack is symbols.
A record that carries a placeholder or no symbol is invisible to a
name-based query and legible only to a sequence search. Assigning symbols
to records whose sequence already places them would recover more
findability than reannotating any genome.

**The audit has limits that bound what it can claim.** Its ground truth is
one instrument, a protein-to-genome sweep, and a locus that instrument did
not find cannot appear here as an annotation failure. The sweep misses
genes in fragmented assemblies ({paper:retention}), so the failure rates
reported are, if anything, conservative on poor assemblies. The
completeness bar was inherited rather than chosen, but it was validated
against the population it acts on, and every state is re-counted across
the whole plausible range. The excess of IP₃ receptor loci held only by
non-coding features is real in the record set, yet it disappears on
assemblies able to carry the gene. We therefore do not claim it as a
property of the family. The expression evidence covers the seven failures
in three fish species that have public RNA-seq. It shows that those
particular unannotated junctions are spliced, not that every unannotated
locus is transcribed. The 297 corrections are proposals from one
instrument about genes that instrument recovered; two are validated to the
exon and by reads, and the rest rest on the same instrument applied at
scale.

The deliverable is addressed to somebody else. Each correction carries
coordinates, the current state, a proposal and an evidence file a curator
can open, and the priority rule says which ones a curator can act on
without re-examining the assembly.
