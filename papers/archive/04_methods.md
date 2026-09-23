## Methods

### Scope

The genome scope is 309 vertebrate assemblies from NCBI Datasets [R159],
declared before any search as the best assembly of every vertebrate order,
ranked by annotation status, RefSeq category, assembly level and scaffold
N50, together with margin species chosen by stated rules from the protein
census. 168 of the assemblies carry a RefSeq annotation and 141 a GenBank
annotation or none; every locus is scored against the gene set that ships
with its own assembly and against no other. The genes audited are the
sweep's gene-scale loci, 2,144 in all, of which 1,059 are IP₃ receptor
loci. Loci in assemblies that ship no gene set are counted and reported as
unscorable, leaving 932 scorable IP₃ receptor loci and 942 scorable
ryanodine receptor loci. The protein scope is the 11,402 full-length family
records (at least 2,000 residues, not flagged as fragments, from a protein
database rather than a genome model), of which 8,306 are vertebrate. The
reference-proteome scope is the 763 vertebrate UniProt reference proteomes
[R157].

**Contiguity is carried as a stratum everywhere.** A locus is on a contig
able to carry the gene when that contig is at least as long as the median
measured IP₃ receptor gene span, 142,212 bp. Every comparison is reported
over all scorable loci and again over loci that clear this bar, because a
contig shorter than the gene cannot be annotated completely by anybody.

### The ground truth

Each assembly was searched with miniprot [R164] using a panel of 38 bait
proteins, 30 IP₃ receptors across the vertebrate clades and 8 ryanodine
receptors as a presence control. Loci were clustered, called to family on
alignment score before any paralogue was considered, and assigned to a
paralogue cell inside the winning family. Empty cells were re-searched by
whole-genome tblastn [R158]. The sweep reads no annotation when it places
a gene, which is what allows it to serve as an independent reference for
the annotation. Each locus carries its coding blocks in gene order, its
bait coverage and identity, and its frameshift and stop counts.

### Locus states

Every scorable locus was assigned one of five states against the
assembly's own GFF3 annotation: *complete* (one coding model covers at
least half the locus's aligned coding footprint), *split* (two or more
coding models each contribute a piece and none reaches the bar),
*fragmentary* (the coding models present deliver less than the bar),
*non-coding only* (a same-strand feature covers the locus but emits no
protein, most often a pseudogene) and *unannotated* (no same-strand
feature). Three rules make the state a measurement. Only same-strand
features count, because an antisense gene is not a model of this one.
Loci are scored against coding blocks and never gene spans, because these
genes' introns reach 152 kb and a passenger gene's span inside one would
otherwise cover the locus. The name is read from whichever model the
annotation places over the locus, preferring a small correctly named
model to a large unnamed one, because accusing a database of failing to
name a gene it did name is the costlier error.

The name verdict has two values that exist because their absence
manufactured errors. A name that claims the family but no paralogue is
recorded as *paralogue unspecified*, not as wrong. A name that claims the
superfamily in a form naming both families is recorded as *family
ambiguous*, not as a wrong-family call.

### Controls

**The sister family is the central control.** Every ryanodine receptor locus the sweep placed was
scored through the same states, names and bars in the same assemblies.
Family and control were compared state by state with Fisher's exact test
and Benjamini–Hochberg correction [R170] across all states within each
scope. This is the comparison that separates a statement about this family
from a statement about vertebrate gene sets. Its outcome is reported in the
Results: no overall difference in either scope.

**The completeness bar was validated, not chosen.** It was inherited from the sweep, where it
decides whether a found gene counts as annotated, and was validated on
1,077 loci whose own annotation names the correct paralogue, whose
alignment recovers at least 90 % of the bait, and whose contig clears the
contiguity bar. The bar sits at that distribution's 1.3 % point. Every
state was re-counted at bars from 0.3 to 0.95.

**The source contrast was run twice**, over all scorable loci and again over loci
that clear the contiguity bar, with the same test and correction, so that
the archive effect and the assembly effect are reported side by side.

### Protein records and proteomes

Each full-length record was compared by blastp against the bait panel,
assigned to a family only when its best score in one family beats the
other by more than a relative margin of 0.1, and to a paralogue within
the winning family by the same rule. The family question was asked of
every record and the paralogue question of vertebrate records only. A
record whose symbol names a paralogue the panel carries no bait for is
reported as outside the instrument's reach. A protein-record rename was
proposed only when the record's best score also cleared 200 bits. Each
reference proteome with no family hit in the profile sweep [R153] was
resolved against an assembly of its own species under four verdicts:
gene caller missed it, genome also empty, assembly cannot carry the gene,
and undecidable because no assembly is in scope.

### Reachability

For each genome and each paralogue cell in which the sweep demonstrates a
gene (1,232 genome-by-cell genes, 923 of them IP₃ receptors), a gene counts
as reachable when at least one full-length protein record of the same
species resolves to that cell. Unreachable genes are sorted by reason: no
reference proteome for the species, only fragmentary records, records
present but none resolving to the cell, or no family record at all. The
same classification was applied to the ryanodine receptor cell. The
protein-level search channels were compared inside the vertebrate
reference proteomes, counting only records each channel calls family and
only accessions the swept files hold.

### Case selection and validation

Annotation loss is the fraction of a recovered gene's coding footprint that
no single annotated model delivers. It was ranked over every recovered IP₃
receptor locus that passes five eligibility rules: the assembly carries a
gene set; the sweep recovered the gene at coverage of at least 0.7 with no
contig edge or gap; the contig holds the locus; contig N50 is at least ten
times the locus span; and at least ten annotated protein-coding genes
elsewhere in the same genome are longer than the locus. The fifth rule is
the per-genome annotation-depth control; it removes 6 loci in 2 genomes.
Failures were labelled omission, truncation or fragmentation, and the
worst of two modes was taken, one case per genome, with ties broken by how
many other family loci in the same genome the annotation gets right.

For each case the spliced coding sequence was reconstructed from the
aligned blocks and its internal stops compared with the number expected
under neutral drift to the observed divergence, computed from the locus's
own codon usage. The comparison is one-sided: zero stops falsifies a
pseudogene call, and a handful would not establish one. Splice
dinucleotides were read from the genome with the strand corrected.
**Exon-boundary concordance** is the control on the instrument: every
swept genome whose same paralogue was aligned from the same bait, and
whose own annotation independently places one model over the whole
alignment, contributes its boundaries, and each case boundary is scored by
the fraction of those genomes sharing it. Annotated proteins were
translated from each assembly's own annotation and tiled onto the genome's
own recovered loci by blastp. Junction probes of 90 nucleotides either
side were searched against every transcript record of the species, counting
only hits contiguous across the junction with at least 8 nucleotides on
each side. The same probes were searched against the locus's own genomic
DNA as a negative control that cannot span a junction.

### Expression

The expression panel took every eligible locus with annotation loss above
0.5 in species with at least four public Illumina RNA-seq runs, plus every
other family locus of the same genomes as internal controls. Runs were
chosen from the Sequence Read Archive [R163] by tissue keyword, verified
against each run's own sample attributes, spread across studies, and
streamed through HISAT2 [R165] with spliced alignment disabled, because the
reference is already coding sequence. Each reference sequence is the
spliced genomic coding sequence of its locus, deliberately not corrected
for frameshifts, validated by placing every block colinearly in the
sweep's own protein. Three housekeeping anchors were built the same way.

**Two controls make the counts readable.** Every reference carries a
reversed decoy of identical length and composition, and a locus is
detected in a run only when it has at least two junction-spanning reads,
at least five reads, and more reads than its own decoy. Across all runs
the decoys collected 42 reads in 2 of 469 run-by-locus comparisons, and
every affected locus still cleared its decoy. Cross-mapping between
paralogues was measured rather than assumed: every reference was tiled
with 100-nucleotide synthetic reads at 10-nucleotide steps and mapped back
under the same settings, and 0 of 12,500 reads were assigned to a sequence
other than their source. A read spans a junction when one contiguous
aligned block covers at least 8 nucleotides on both sides. Junctions with
no reads are written as zeros, not omitted.

### Corrections

A correction is written for every scorable locus not delivered complete
and for every protein record whose name and sequence disagree on family.
Priority is high when the gene is demonstrably present, its contig clears
the contiguity bar and its reading frame is intact; medium for a partial
recovery or an unscored reading frame; and low below the contiguity bar.
Where the reading-frame screen scores a locus as lesion-rich, the row is
written, marked as withheld and given its reason, so the list records what
the audit declined to propose.

### Negative controls and reproducibility

Every stage runs a suite of constructed negative controls before writing
anything: 44 for the locus and protein audit, 47 for the case validation
and 11 for the expression analysis. Each suite was mutation-tested by
breaking rules deliberately. Every table and figure is generated by
script from committed tables, and every number in this paper is declared
in a claims ledger with the table it comes from and re-verified on every
build.
