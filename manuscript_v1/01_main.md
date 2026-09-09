## Introduction

Inositol 1,4,5-trisphosphate releases calcium from a non-mitochondrial
intracellular store [R01], and the channel that does it is a
~2,700-residue, four-fold symmetric protein of the endoplasmic reticulum
membrane [R02, R03]. From its first sequence it was recognised as the
smaller relative of the ryanodine receptor [R04], and the two families
have been read against each other ever since: they share the MIR, RIH and
RIH-associated domains, the six-transmembrane pore module, and — as this
work repeatedly had to accommodate — every Pfam signature used to
diagnose either one. In vertebrates the receptor is encoded by three
genes, *ITPR1*, *ITPR2* and *ITPR3* [R32], with distinct tissue
distributions [R29, R33], distinct single-channel behaviour [R67] and
distinct disease associations: *ITPR1* in spinocerebellar ataxia 15 and 29
and in Gillespie syndrome [R38, R42, R43, R44], *ITPR2* in isolated
anhidrosis [R45], and *ITPR3* in Charcot-Marie-Tooth neuropathy and a
multisystem disorder [R46, R47, R48].

Almost everything known about the family's distribution outside the
vertebrates rests on searches that were never scoped. Surveys of calcium
signalling machinery report IP₃ receptors in animals, amoebae and various
protists and their absence from land plants and from the yeasts
[R35, R36, R107, R108], but an absence is a claim about a search space,
and the search spaces behind those statements are not declared, not
controlled, and not separable from the annotation quality of the databases
they were run in. The same is true of the vertebrate paralogues. Their
number is known; whether any lineage has lost one is not, because no study
has asked the question of genomes rather than of gene sets, and a gene
missing from a gene set is more often an annotation than a biological
fact.

Three properties of this family make it a hard case, and each of them
shapes what follows. It is long, so its locus is longer than the contigs
of a large fraction of published assemblies. It is paralogous with a
sister family that carries all of its diagnostic domains, so any search
sensitive enough to find divergent IP₃ receptors also finds ryanodine
receptors. And it is poorly named: much of its record sits under
placeholder symbols or under no symbol at all, so name-based retrieval
systematically under-counts it.

We therefore built the census as a genomic one, with the sister family
searched alongside as an internal control at every stage, and with every
absence required to carry a positive control demonstrating that the search
reached the assembly. The scope is declared and finite throughout: 763
vertebrate and 6,928 other reference proteomes, 309 vertebrate genomes
selected as one per order plus every species the protein record left in
doubt, and 194 non-vertebrate genomes sampled most finely where the
negative claims are. Every claim below is a claim about that space.

