## 6. Three paralogues, and the diversity within each

### 6.1 The genes

Humans carry three *ITPR* genes. Their protein products are near-identical in
length — IP<sub>3</sub>R1 2,758 aa (UniProt Q14643), IP<sub>3</sub>R2 2,701 aa
(Q14571), IP<sub>3</sub>R3 2,671 aa (Q14573), a spread of 3% — but the genes
that encode them are not built alike. Re-derived from Ensembl release 15.12 for
this review:

| Gene | Locus | Ensembl gene | Canonical exons | Genomic span |
|------|-------|--------------|-----------------|--------------|
| *ITPR1* | 3p26.1 | ENSG00000150995 | 62 | 354,174 bp |
| *ITPR2* | 12p11.23 | ENSG00000123104 | 57 | 497,888 bp |
| *ITPR3* | 6p21.31 | ENSG00000096433 | 58 | 76,245 bp |

The genomic span varies **6.5-fold** across three genes of essentially
identical exon count and protein length. *ITPR3* achieves the same architecture
in roughly one-fifth of the DNA that *ITPR2* uses. Whether this reflects
lineage-specific intron gain, differential loss, or selection related to the
paralogues' distinct expression programmes is not established, and it is a
question the genomic literature has largely not asked.

![](figures/gene_architecture.png)

**{fig:gene_architecture}.** Three genes that make near-identical proteins out of very different amounts of DNA. (**a**) Genomic span; the exon count and protein length that stay constant are printed with each gene. (**b**) The asymmetry as one number: *ITPR2* spends 184 bp of locus per residue of protein and *ITPR3* spends 29. Re-derived from Ensembl release 15.12 for this review — the claim this replaced, that all three span "hundreds of kb", is one of the four that failed audit (§12.2).


### 6.2 Functional differences between paralogues

The paralogues are not interchangeable. Single-channel analysis showed
isoform-specific gating [R67]; the type 2 receptor was identified and
reconstituted separately and shown to differ functionally [R105]; ATP
regulation differs between types [R70]; the receptors differ in
IP<sub>3</sub> sensitivity, in Ca<sup>2+</sup> dependence and in their
regulation by the modulators of §4. The cleanest whole-cell demonstration
remains the finding that Ca<sup>2+</sup> signal *encoding* — the shape and
frequency of the response — is set by which subtypes a cell expresses [R101],
which recasts subtype composition as a tuning parameter rather than redundancy.
They also differ in the ER–mitochondrial contacts they support [R12].

![](figures/conservation_profile.png)

**{fig:conservation_profile}.** Constraint across the family, and how far apart its members actually are. (**a**) Per-column conservation from an alignment of the 25 control positives — IP<sub>3</sub>R1/2/3 across a vertebrate panel plus the invertebrate and non-metazoan single-*itpr* grade — mapped onto human IP<sub>3</sub>R1 numbering, over that protein's domain architecture. The ligand-contact and pore positions measured in {fig:channel_structure} are marked; both sit on conservation maxima. The gap trace shows where the alignment is carrying indels, and the deepest troughs are the linkers between domains rather than the domains themselves. (**b**, **c**) All-pairs identity over mutually covered columns. The three human paralogues are 66–72% identical to each other and about 25% identical to a ryanodine receptor — and the fly receptor is closer to human IP<sub>3</sub>R1 (62%) than any IP<sub>3</sub>R is to any RyR. Paralogue identity is high enough that the differences of §6.2 are differences of tuning, not of kind.


### 6.3 Heterotetramers

Because the paralogues co-assemble, most cells contain a mixed population of
tetramers rather than three separate channel species. This has two
consequences that matter throughout this review. It expands the functional
repertoire combinatorially. And it converts a heterozygous missense variant in
one gene into a defect distributed across a large fraction of a cell's
channels, which — combined with the all-four-sites requirement of §3.3 — is why
dominant-negative disease mechanisms are so prominent in this family (§9).

### 6.4 Alternative splicing

*ITPR1* carries three alternatively spliced segments. SI and SII were defined
in mouse, with tissue- and development-specific usage [R29]; the SII region
distinguishes neuronal from non-neuronal forms and the two differ in
phosphorylation [R30]; SII splicing was characterised in rat brain versus
peripheral tissues [R100]; and a third spliced segment, SIII, was identified on
cloning the human cDNA [R31]. Because SII lies in the regulatory region between
the ligand-binding and channel domains, splicing here changes how the receptor
responds rather than whether it conducts.

### 6.5 Tissue distribution

The classical picture — reviewed comprehensively for the three subtypes [R33] —
is of broad expression with strong paralogue bias:

- **IP<sub>3</sub>R1** dominates the cerebellum, where it was discovered as the
  Purkinje-cell P400 antigen [R55, R21], and is the principal neuronal isoform
  [R103].
- **IP<sub>3</sub>R2** is broadly expressed, prominent in secretory epithelia
  [R45], cardiac myocytes [R116] and glia — including oligodendrocytes, where it
  is required for myelination [R117].
- **IP<sub>3</sub>R3** is the epithelial and secretory isoform, concentrated
  apically in pancreatic and salivary gland cells [R102] and central to
  hepatobiliary function [R34, R125].

Two cautions. First, this distribution is not fixed: subtype expression is
remodelled in disease, as shown for the two intracellular release-channel types
in end-stage heart failure [R104], and quantitative modern atlases have not been
applied systematically across the family. Second, and more importantly, §5.3 shows that
the functionally relevant quantity is not how much receptor a cell contains but
how much of it is licensed and correctly placed — so expression level is a weak
predictor of signalling capacity.
