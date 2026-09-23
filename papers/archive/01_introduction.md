## Introduction

A gene that no protein database holds is a gene that gene-family studies,
variant interpreters and structure predictors never see. Most comparative
work starts from protein records, whether by name, by domain signature or by
sequence similarity against a proteome, and every one of those routes
depends on a genome annotation having delivered the gene as a coding model
and a database having served its product. Automated annotation is known to
fail [R194], and new assemblies arrive faster than curation can follow them
[R195]. What is rarely available is a measurement, for one gene family
across a declared set of genomes, of how often the record fails, in which
way, and whether the failure belongs to the family or to the archive.

The inositol 1,4,5-trisphosphate (IP₃) receptors are a demanding test of
that record. They are the endoplasmic reticulum's ligand-gated calcium
release channels [R07], encoded in vertebrates by three paralogues, *ITPR1*,
*ITPR2* and *ITPR3* [R32]. Each gene is about 2,700 codons spread over some
58 exons and tens to hundreds of kilobases of genomic sequence. Mutations in
them cause ataxias and other neurological disease [R122], so a missing or
mis-assembled record has practical cost. The family also has a built-in
control. The ryanodine receptors are the IP₃ receptors' sister family
[R04, R37]: similar in size, sharing every diagnostic domain [R58], and
annotated by the same pipelines in the same assemblies. A failure rate
measured on the IP₃ receptors alone cannot be told from the general quality
of vertebrate gene sets. Measured beside the ryanodine receptors, it can.

Here we ask how completely the public databases record the IP₃ receptor
genes that vertebrate genomes carry. The ground truth is a protein-to-genome
alignment sweep over 309 vertebrate assemblies [R164], which recovers each
gene from the DNA without consulting any annotation, so it can be used to
audit the annotations. Every locus it places is scored against the
assembly's own gene set as complete, split, fragmentary, held only by a
non-coding feature, or unannotated, and the ryanodine receptors are scored
the same way in the same genomes. We then ask the question from the reader's
side: which demonstrated genes can be reached from a protein database at
all ({fig:contribution}). We compare curated RefSeq annotations with
submitter-deposited GenBank ones ({fig:by_source}) and the family with its
sister ({fig:family_control}), and score every full-length protein record
the databases hold for its family, its paralogue and its name
({fig:protein_side}). Finally we take two failures down to the exon
({fig:exon_tracks}) and ask public RNA-seq whether the junctions no
annotated model spans are transcribed ({fig:junctions}).

The answer divides the blame differently from the way the question was
posed. The family is recorded no worse than its sister, and the protein
records that exist are almost always named correctly. The failure lies in
which archive annotated the genome, in assembly contiguity, and above all
in records that do not exist.
