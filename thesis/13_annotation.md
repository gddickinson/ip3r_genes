# 13. An archive that cannot find the gene

## 13.1 Three quarters of these genes are unreachable from any protein database

Chapter 4 measured that **940 of 1,232 demonstrated genes, or 76.3 %, are not
reachable by any protein-database search.** They are not hard to find, because
no protein record of them exists.

That is a statement about databases rather than about biology, and it is
easily the most immediately useful result in the thesis. It is also a claim
that has to be made carefully, because there are three quite different ways a
gene can be missing from a public record, and they have different remedies.

The gene may not be annotated at all. It may be annotated and split across
several models, none of which delivers the protein. Or it may be annotated
perfectly, named perfectly, and filed as something that emits no protein.

All three happen. This chapter measures how often, audits whether it is worse
for this family than for its sister, validates two cases down to the exon,
confirms with reads that the lost genes are transcribed, and produces a
correction list addressed to somebody else.

## 13.2 Three rules decide whether the audit is measuring anything

Every gene-scale locus in the 309-genome scope is scored into one of five
states, comprising complete, split, fragmentary, non-coding only and
unannotated, under three rules that decide whether the audit is measuring
anything.

**Only same-strand features count**, because an antisense gene overlapping the
locus perfectly is not an annotation of it.

**Loci are scored against coding blocks rather than gene spans**, because a
gene span covers its own introns and these genes' introns reach 152 kb, so a
passenger gene inside one would score as covering the locus.

**The name is read off whichever model the annotation actually places there**,
deliberately generously, because accusing a database of failing to name a gene
it did name is the expensive error.

The name verdict itself has two values that exist because their absence
manufactured errors. A name that claims neither family, and a name that claims
the family but no paralogue, are recorded as such rather than as wrong.
Without those two values the audit produced 66 wrong-family and 6,660
wrong-paralogue errors out of names that claim neither.

**The bar for "complete" is validated rather than re-derived.** The audit
inherits the sweep's coverage threshold and scores it against the distribution
it sits in, and every state is re-counted across the whole bar range so the
reader can see what the choice costs.

Forty-four constructed negative controls run before anything is written, and
the suite was mutation-tested on nine deliberate rule breakages, all nine
caught. The measurement is cached and the verdicts are not, because the first
version cached both and a later rule change survived into a committed table.

## 13.3 A quarter of annotated loci are not delivered as one complete model

Of the 932 loci in assemblies that ship a gene set, **73.9 % are delivered as
one complete annotated model and 26.1 % are not.** Thirty-eight have no
same-strand annotated feature at all, and 42 are held only by a non-coding
one, so the database serves no protein for either.

**The single largest determinant is not the paralogue but the assembly.**
Above the contiguity bar the failure rate falls from 26.1 % to **6.7 %**, so
two thirds of what looks like an annotation problem is a contig too short to
hold a 2,700-residue gene.

By class the spread is wide. Cartilaginous fish and mammals are
near-perfect, ray-finned fish good, and **birds are complete at 56 % with 29 %
fragmentary**, which is the same class the contiguity bar removes two thirds
of.

![](figures/s18_calibration.png)

**{fig:s18_calibration}.** The completeness bar drawn inside the
distribution it sits in, because a number in a legend cannot show a tail.
The importance of drawing the bar inside its own distribution is that the
whole audit turns on it. Every locus in this chapter is called complete or
not against one coverage threshold, and a threshold quoted as a number in a
sentence cannot show whether it sits in a gap or in the middle of a
population. The bar was validated rather than replaced here, and this panel
is the evidence for that decision.

## 13.4 Curated and submitter-deposited gene sets differ, and only partly because of assembly quality

A curated reference gene set and a submitter-deposited one are not comparable
evidence, and they are not close.

**98.8 % of loci in curated annotations are complete, against 37.5 % in
submitter-deposited ones.**

![](figures/s18_by_source.png)

**{fig:s18_by_source}.** Locus state by annotation source, raw and above the
contiguity bar, with the control drawn beside the raw contrast rather than
instead of it. The importance of this figure is that it separates a database
effect from an assembly effect, which no raw comparison of annotation
sources can do. Curated and submitter-deposited gene sets are not applied to
the same assemblies, so a difference between them partly measures which
genomes each was run on. Drawing the contiguity-controlled contrast beside
the raw one shows how much of the gap survives that control and how much
does not.

The confounder is obvious and is controlled rather than argued. Submitter
assemblies are less contiguous, and a locus on a contig too short to hold the
gene cannot be annotated completely by anybody. Re-run over the loci that
clear the bar, the contrast **narrows and does not close**, at 99.4 % against
63.8 %. About 42 % of the gap between the two archives is assembly quality and
the rest is the gene set.

## 13.5 The family is not recorded worse than its sister family

The audit's premise was that a 2,700-residue, 58-exon gene with a sister
family sharing every diagnostic domain should be badly recorded.

**The control says otherwise, and that is the result.** The IP₃ receptor cells
fail on 26.1 % of loci, and the ryanodine receptors, in the same assemblies
through the same pipelines, on 22.1 %. Above the contiguity bar the two are
6.7 % and 5.5 %. **No overall difference survives correction. How this family
is recorded is how vertebrate genes of this size are recorded.**

![](figures/s18_family_vs_control.png)

**{fig:s18_family_vs_control}.** This family against its sister in the same
assemblies. Without this comparison, a failure rate is not a statement about
this family at all. This figure carries the result that overturns the
chapter's own premise. The audit was built expecting a large, many-exon gene
with a confusable sister family to be recorded badly, and the sister family
in the same assemblies through the same pipelines fails at a comparable
rate. Without a control of this kind, a quarter of loci not delivered as one
model reads as an indictment of this family; with it, it reads as the
failure rate for vertebrate genes of this size.

**One state does separate, and it is the family-specific one.** An IP₃
receptor locus is 2.7 times more likely than a ryanodine locus to be held only
by a non-coding feature, in 42 loci against 16. That is a gene the annotation
identifies correctly, names correctly, and files as non-coding, so no protein
record is ever created and no name-based search can reach it.

**It does not survive the contiguity control.** Above the bar the same
comparison is 4 against 5, so the effect is confined to assemblies too broken
to carry the gene, which is where a submitter has most reason to file a model
as non-coding in the first place. The excess is real in the record set, and
this measurement cannot say it is a fact about the family rather than about
the assemblies its loci happen to sit in.

## 13.6 Half the protein records carry no usable gene symbol

Each of the 11,402 full-length family protein records was scored against the
committed bait panel, assigned to a family only on a stated margin, then to a
paralogue inside the winning family, with the record's own gene symbol and
protein name read through the same verdict rule the genome half uses.

**The family call is not in dispute.** The sequence disagrees with the census
call on 4 of 11,402 records, and the databases name the sister family at only
5, every one a non-vertebrate record where the sequence barely separates the
families either. The hazard measured at 49 % of records in one query in
Chapter 2 does not appear as a naming error in the full-length record set.

**The paralogue call is not in dispute either.** Fifty-two of 8,306 vertebrate
records carry a symbol naming a paralogue the panel assigns elsewhere, and a
further 74 claim a paralogue **this panel cannot call**, because the bait
table records no bait for it, so those are reported as outside the
instrument's reach rather than as errors. Without that guard they would have
been 74 database naming errors manufactured by a missing bait.

**What is in dispute is whether the records are findable.** 3,872 records
carry a placeholder gene symbol and 2,395 carry none at all, so **55.0 % of
the family's full-length protein records have no usable gene symbol.** On the
protein-name side, 66 are named for the superfamily in a way that names the
family and its sister together and therefore separates neither.

![](figures/s18_protein_side.png)

**{fig:s18_protein_side}.** The protein records, showing what the name
claims against what the sequence is. The importance of separating what a
name claims from what a sequence is is that the two failure modes have
different remedies. The family and paralogue calls turn out not to be in
dispute, so this is not a record full of misidentified proteins; more than
half of it simply carries no usable gene symbol. A record that is correct
and unfindable needs a name, not a reannotation, and that distinction is
what the panel makes visible.

**Would a signature query have found them?** The signature that names this
family reaches 95.6 % of the records both instruments agree on and 81.2 % of
those only one instrument found. It is a good but not complete index of its
own family, and the deficit is concentrated outside the vertebrates, falling
to 64 % in the amoebozoa and the stramenopiles.

## 13.7 All fifteen empty proteomes are gene-caller failures

Chapter 3 swept 764 vertebrate reference proteomes and 15 returned no family
hit at all. On its own that is uninterpretable, because a receptor-shaped hole
in a proteome is either a gene the species lacks or a gene its gene caller did
not find.

Each was resolved against an assembly of its own species, searched with an
instrument that owes the gene caller nothing.

**All 15 are gene-caller failures.** Every one of the 15 species has a genome
in scope, and in every one the genomic sweep recovers at least one locus at
over half the bait's length while the proteome holds none.

The verdict vocabulary carries a fourth value that fires on none of them, and
that is the point of having it: a species with no assembly in scope would be
recorded as undecidable rather than as an absence. Being able to say it fires
on none is worth more than assuming it would.
