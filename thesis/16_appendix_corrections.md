# Appendix A — the correction list

This appendix is the one deliverable of the thesis addressed to somebody else.
It is a list of specific, coordinate-level proposals to the public archives,
each with the evidence a curator would need to act on it, and it is committed
in full as `results/annotation_audit/corrections.tsv`.

**297 items. 52 high priority, 19 medium, 226 low. 18 vetoed.**

## What a row contains

Every row carries the assembly accession, the annotation source, the organism
and its vertebrate class, the paralogue cell, the contig with start, end and
strand, the current annotation state, the current name and biotype, the
proposal, the numbers the class fired on, the priority with the rule that
assigned it, and an archived evidence file a curator can open.

## The six classes

**Unannotated (20 items).** No same-strand annotated feature overlaps the
recovered gene at all. The proposal is a coding model over the aligned exons.

**Held only by a non-coding feature (54 items).** The annotation places a
feature there and files it as something that emits no protein — most often a
pseudogene. This is the class where the family differs from its sister, and
where the failure is most invisible to a user: the gene is identified,
frequently named correctly, and simply produces no protein record, so no
name-based or sequence-based protein search can ever reach it.

**Incomplete (217 items).** One or more coding models are present and between
them deliver a fraction of the gene. The proposal names how many models are
there, what fraction they deliver together, what the best single one delivers,
and asks for extension or merging to the aligned exons. The typical row
delivers a tenth of the gene across two models.

**A wrong paralogue name (1 item).** A locus whose annotation names a
paralogue the alignment assigns elsewhere, in a genome that also carries a
separate locus for the paralogue the annotation names — so the two cannot both
be that gene.

**A protein record naming the wrong family (5 items).** Full-length protein
records whose sequence the bait panel assigns to one family and whose name
claims the other. All five are non-vertebrate records where the sequence
barely separates the families either, and each additionally required an
absolute score as well as a relative margin before it was written, because
these candidates score 90 to 191 bits over roughly 2,700 residues.

Two hundred and ninety-two of the items are about a genome annotation and five
about a protein record.

## Priority is a rule

**High** requires three things simultaneously: the gene demonstrably present,
the assembly demonstrably able to carry it, and the reading frame intact. Only
52 items meet all three, and they are concentrated in ray-finned fish, with a
handful in birds, mammals and amphibians.

**Medium** is a partial recovery or an unscored reading frame.

**Low** is anything below the contiguity bar — because there the annotation's
silence may be the assembly's fault rather than the annotator's, and a
correction proposed on that basis would be asking a curator to annotate a gene
the assembly cannot represent.

## The veto, applied as a column

Eighteen items are vetoed. Where the lesion screen scored a locus as
lesion-rich, this audit does not propose resurrecting it.

**The row is still written**, flagged, with its reason. A locus the audit
declined to correct is evidence about the audit, and filtering those rows out
would leave a list that looks cleaner than the evidence behind it.

## What the list is not

It is not a claim that these annotations are wrong in a sense a curator must
accept. It is a set of specific, checkable proposals with coordinates and
evidence, produced by one instrument, about genes that instrument recovered.
Two of them — one omission and one fragmentation — are validated in Chapter 13
down to the exon and confirmed by junction-spanning reads. The other 295 rest
on the same instrument that produced those two, applied at scale.
