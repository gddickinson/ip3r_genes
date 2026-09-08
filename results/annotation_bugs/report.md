# S10 — annotation-bug molecular validation

The genome sweep recovered **880** IP3-receptor loci across the 309-genome scope. **382** of them are in assemblies where the annotation could reasonably have been expected to deliver the gene (rules E1–E5 below), and at **375** of those (98.2 %) it did. This report is about the **7** where it did not, and about the two of them, one per failure mode, whose evidence is taken down to the exon.

The two cases are ***Nibea albiflora*** ITPR2 (omission) and ***Dissostichus eleginoides*** ITPR3 (fragmentation).

> Every number below is read from a committed table in this directory. Nothing in this report is computed from a GFF, a genome or a BLAST run at render time (D13).

## 1. Which failures, and why those

"The two worst" is a rule applied to every recovered locus, not a pair of loci chosen by eye. The measurement is **annotation loss** — the fraction of a recovered gene's coding footprint that no single annotated gene model delivers — taken from the `frac_cds` values S5 recorded at sweep time.

### Eligibility

| rule | test | loci excluded |
|---|---|---|
| E1 | assembly carries a gene set and the sweep indexed it | 103 |
| E2 | sweep recovered the gene: coverage >= 0.7, no contig edge, no N-gap | 314 |
| E3 | the contig holds the whole locus (D4) | 23 |
| E4 | contig N50 >= 10x the locus span | 52 |
| E5 | >= 10 annotated protein-coding genes elsewhere in the same genome are longer than the locus | 6 |

**382 loci pass all five.** Their median annotation loss is 0 and 360 of them (94.2 %) sit at exactly zero — one annotated gene model covering the whole recovered coding footprint. The failures are a thin tail, not a distribution.

### E5 does the work, and here is what it removed

E5 asks whether the annotation builds genes this long *anywhere else in the same genome*. Without it the ranking's top rows are loci in assemblies whose annotation has a genome-wide length ceiling, and validating one of those would report a property of the whole gene set as a bug at this gene.

| genome | locus | span | loss | longest gene the same annotation builds |
|---|---|---|---|---|
| *Cirrhinus mrigala* | ITPR3 | 41.5 kb | 0.971 | 42.0 kb |
| *Cirrhinus mrigala* | ITPR1 | 38.9 kb | 0.942 | 42.0 kb |
| *Cirrhinus mrigala* | ITPR2 | 116.3 kb | 0.942 | 42.0 kb |
| *Saguinus oedipus* | ITPR1 | 332.6 kb | 0.909 | 137.8 kb |
| *Saguinus oedipus* | ITPR2 | 726.2 kb | 0.835 | 137.8 kb |
| *Saguinus oedipus* | ITPR3 | 310.0 kb | 0.541 | 137.8 kb |

All 6 excluded loci are in 2 genomes, and E5 removes nothing else in the whole sweep.

### The failures

| genome | locus | mode | loss | annotated models | internal control | selected |
|---|---|---|---|---|---|---|
| *Nibea albiflora* | ITPR2 | omission | 1.0 | 0 | 2 of 3 | **case_a** |
| *Dissostichus mawsoni* | ITPR1 | omission | 1.0 | 0 | 1 of 3 | — |
| *Dissostichus eleginoides* | ITPR2 | omission | 1.0 | 0 | 1 of 3 | — |
| *Dissostichus mawsoni* | ITPR3 | truncation | 0.989 | 1 | 1 of 3 | — |
| *Dissostichus mawsoni* | ITPR2 | truncation | 0.876 | 1 | 1 of 3 | — |
| *Dissostichus eleginoides* | ITPR3 | fragmentation | 0.581 | 3 | 1 of 3 | **case_b** |
| *Dissostichus eleginoides* | ITPR1 | fragmentation | 0.537 | 2 | 1 of 3 | — |

`internal control` counts the *other* family loci in the same genome the annotation gets right; it breaks ties, because at equal loss the sharper case is the one whose own genome proves the annotation could have done better.

The two selected cases are the worst of each mode, at most one per genome. Taking the top two of the single ranking would have given two omissions and left the fragmentation claim unvalidated; the full ranking is committed either way (`case_ranking.tsv`).

## 2. What this changes about earlier tasks

### The naming conflicts S5b left open — **confirmed**

*Prior.* S5b: 3 loci are claimed by a paralog cell other than the one the assembly's annotation names, left unadjudicated ('one instrument does not overturn a public annotation') (results/genome_ledger/report.md, 'Where an annotation and this sweep disagree').

*Here.* 11 annotated models in the two case genomes carry a paralog name and could be placed at their own locus; 1 of them sits on a locus of a different paralog — one of the three conflicts S5b listed, the two others being in genomes outside these cases. The adjudication is the annotation's **own translated protein** blastp'd against the genome's own recovered loci — a second instrument, and the one S5b said was needed before a public annotation could be contradicted.

- `ITPR2` (2,521 aa) matches the **ITPR3** locus at 89.3% identity over 2,683 residues, bit-score margin 0.1363 over the runner-up.

### What the protein databases hold for these two species

| species | records in census v3 | ITPR records |
|---|---|---|
| *Dissostichus eleginoides* | 4 | 0 |
| *Nibea albiflora* | 1 | 0 |

*Prior.* S5b: 318 ITPR gene models exist only as DNA and a further 167 sit inside an annotated gene carrying no family name.

*Here.* Both case species carry three complete IP3-receptor genes in their DNA and **no ITPR protein record at all** in the census. That is the consequence S5b's counts describe, measured at the two species where the cause is now known.

## 3. Case A — *Nibea albiflora* ITPR2 (omission)

`GCA_014281875.1` · CM024800.1:12,447,838–12,500,306- · 52.5 kb · recovered from bait `A0AAX7UHL1` at 86.6% identity over 100% of the bait.

### Is this even the right build?

The annotation's GFF3 header names build `ASM1428187v1` (`GCA_014281875.1`), written by `NCBI annotwriter`. The sweep searched `GCA_014281875.1`, so the annotation and the alignment sit on the **same assembly**: `True`. NCBI lists 3 assemblies for this species; 1 is newer, of which none carries an annotation. The assembly is chromosome-level, released 2020-08-25 by Marine Fishery Institute of Zhejiang Province, contig N50 4,420.1 kb — 84.2× this locus.

NCBI's own counts for this annotation: 15,965 protein-coding genes and 7,380 pseudogenes out of 23,345 — 31.6% of the gene set is filed as pseudogene.

### What the annotation put on the gene

The recovered gene has **56 coding exons** totalling 8,021 bp. The annotation places **0 gene models** over them, contributing 0 disjoint coding blocks and covering 0 bp (0.0%) of the coding footprint. 8,021 bp in 56 blocks has no coding model at all, and 8,021 bp reaches no protein.

Blocks rather than genes, because that is what a reader can recover: two gene models a few hundred bp apart are two genes but their merged coding evidence is still two pieces.

| exon class | exons | share |
|---|---|---|
| intergenic | 56 | 100.0 % |

**No annotated gene model overlaps any exon of this gene.** Not a partial model, not a mis-named one: the annotation has nothing here.

The annotation is working on both sides of the gap. The nearest annotated genes are `PLEKHG7` (4.2 kb downstream), `SSPN` (5.5 kb upstream), `BHLHE41` (12.2 kb upstream), `EEA1` (17.8 kb downstream).

### What the DNA says

The locus splices into **8,019 nt of coding sequence, 2,673 codons**, translating with **0 internal stop codons** and a terminal stop present. Under neutral drift to the observed 86.6% protein identity, 14.2 stops would be expected (computed from this locus's own codon usage). miniprot reports 2 frameshifts, which at this divergence is not on its own evidence of pseudogeny; the absence of nonsense over 2,673 codons is.

Genome-wide, the annotation carries 14 family-named gene models, of which 8 are filed as pseudogenes and 6 emit a protein.

| model | biotype | span | emits a protein |
|---|---|---|---|
| `RYR3.5` | pseudogene | 81.5 kb | no |
| `ITPR1` | pseudogene | 53.2 kb | no |
| `RYR3.6` | pseudogene | 44.0 kb | no |
| `RYR2.3` | pseudogene | 42.5 kb | no |
| `RYR2` | pseudogene | 42.0 kb | no |
| `RYR3.3` | pseudogene | 40.2 kb | no |
| `RYR2.2` | protein_coding | 22.0 kb | yes |
| `ITPR2` | pseudogene | 20.3 kb | no |

### The name and the sequence disagree

`ITPR2` (2,521 aa, pseudogene) is named for one paralog and its translated protein matches the **ITPR3** locus at 89.3% identity over 2,683 residues — a bit-score margin of 0.1363 over the next locus. The comparison is the annotation's own sequence against this genome's own loci, and the model sits at the coordinates of the locus it matches.

### Five checks on this evidence

1. **Splice sites.** 55 of 55 introns (100.0 %) carry canonical or minor splice dinucleotides (54 GT–AG, 1 minor). An alignment is not a gene; an exon structure whose introns are spliceable is evidence that this one is.

2. **Exon boundaries against independently annotated genomes.** 35 swept genomes carry this paralog, aligned from the same bait, with their own annotation independently placing one gene model over the whole alignment. 52 of 55 of this locus's internal exon boundaries (94.5%) are shared by a majority of them; the median boundary is shared by 35. This validates the instrument at this gene, not the individual locus — which is why it is reported alongside the reading frame and the splice sites rather than instead of them.

3. **Transcript evidence at the junction.** 55 ±90 nt probes were built across the junctions the annotation does not model and searched against every transcript record NCBI holds for this species. 0 hits, 0 of them contiguous across a junction.

    The denominator is the result: this species has 43 transcript records in total. A search of that many returning nothing has not shown the gene is untranscribed — it has shown the species has almost no transcript deposits. RNA-seq for it does exist and reaching it needs the streaming aligner S12 builds; these junctions are handed there rather than half-answered here.

    The criterion has its own negative control. Run against the locus's own genomic DNA — where a probe of contiguous *spliced* sequence cannot span its junction by construction — the probes make 112 hits and 0 of them span. So the zero above is a discriminating zero, not a test that never fires. (On its first run that control returned six false spans, which is how the rule acquired its second half: an 8 nt anchor alone admits an alignment that has run a dozen bases past the junction into the intron.)

4. **The rest of the family in the same genome.** 3 of 3 other family loci in this assembly also splice into an uninterrupted reading frame. Whatever the annotation is doing here, it is not responding to a damaged gene.

5. **The neighbourhood.** 5 of 6 of the nearest flanking genes are in S8's consensus flank set for ITPR2 — including `SSPN`, this paralog's most conserved neighbour, found beside it in 197 of 215 swept vertebrates (91.6%). The gene is not merely present and intact; it is in the position this paralog occupies across the vertebrates, which nothing about this assembly's annotation could produce.

## 4. Case B — *Dissostichus eleginoides* ITPR3 (fragmentation)

`GCA_031216635.1` · CM062267.1:8,174,094–8,242,096- · 68.0 kb · recovered from bait `A0A8C6KWV8` at 88.8% identity over 100% of the bait.

### Is this even the right build?

The annotation's GFF3 header names build `KU_De_1.0` (`GCA_031216635.1`), written by `NCBI annotwriter`. The sweep searched `GCA_031216635.1`, so the annotation and the alignment sit on the **same assembly**: `True`. NCBI lists 2 assemblies for this species; none is newer. The assembly is chromosome-level, released 2023-09-08 by Korea University, contig N50 4,229.0 kb — 62.2× this locus.

### What the annotation put on the gene

The recovered gene has **60 coding exons** totalling 8,064 bp. The annotation places **3 gene models** over them, contributing 38 disjoint coding blocks and covering 4,664 bp (57.8%) of the coding footprint. 3,400 bp in 28 blocks has no coding model at all, and 3,400 bp reaches no protein.

Blocks rather than genes, because that is what a reader can recover: two gene models a few hundred bp apart are two genes but their merged coding evidence is still two pieces.

| exon class | exons | share |
|---|---|---|
| in annotated cds | 37 | 61.7 % |
| intergenic | 23 | 38.3 % |

| annotated model | biotype | span | share of the coding footprint | residues delivered |
|---|---|---|---|---|
| `KUDE01_018466` | protein_coding | 26.2 kb | 40.8% | 468–1593 |
| `KUDE01_017523` | protein_coding | 11.3 kb | 15.0% | 52–467 |
| `KUDE01_018062` | protein_coding | 13.0 kb | 2.0% | 1–67 |

The annotation is working on both sides of the gap. The nearest annotated genes are `KUDE01_017692` (12.7 kb downstream), `KUDE01_018563` (13.9 kb upstream), `KUDE01_017602` (30.9 kb downstream), `KUDE01_018254` (42.8 kb downstream).

### What the DNA says

The locus splices into **8,061 nt of coding sequence, 2,687 codons**, translating with **0 internal stop codons** and a terminal stop present. Under neutral drift to the observed 88.8% protein identity, 11.4 stops would be expected (computed from this locus's own codon usage). miniprot reports 2 frameshifts, which at this divergence is not on its own evidence of pseudogeny; the absence of nonsense over 2,687 codons is.

### Five checks on this evidence

1. **Splice sites.** 59 of 59 introns (100.0 %) carry canonical or minor splice dinucleotides (59 GT–AG, 0 minor). An alignment is not a gene; an exon structure whose introns are spliceable is evidence that this one is.

2. **Exon boundaries against independently annotated genomes.** 40 swept genomes carry this paralog, aligned from the same bait, with their own annotation independently placing one gene model over the whole alignment. 59 of 59 of this locus's internal exon boundaries (100.0%) are shared by a majority of them; the median boundary is shared by 39. This validates the instrument at this gene, not the individual locus — which is why it is reported alongside the reading frame and the splice sites rather than instead of them.

3. **Transcript evidence at the junction.** 25 ±90 nt probes were built across the junctions the annotation does not model and searched against every transcript record NCBI holds for this species. 0 hits, 0 of them contiguous across a junction.

    The denominator is the result: this species has 10 transcript records in total. A search of that many returning nothing has not shown the gene is untranscribed — it has shown the species has almost no transcript deposits. RNA-seq for it does exist and reaching it needs the streaming aligner S12 builds; these junctions are handed there rather than half-answered here.

    The criterion has its own negative control. Run against the locus's own genomic DNA — where a probe of contiguous *spliced* sequence cannot span its junction by construction — the probes make 50 hits and 0 of them span. So the zero above is a discriminating zero, not a test that never fires. (On its first run that control returned six false spans, which is how the rule acquired its second half: an 8 nt anchor alone admits an alignment that has run a dozen bases past the junction into the intron.)

4. **The rest of the family in the same genome.** 3 of 3 other family loci in this assembly also splice into an uninterrupted reading frame. Whatever the annotation is doing here, it is not responding to a damaged gene.

5. **The neighbourhood — not testable here.** None of the flanking genes carries a gene symbol (this annotation names its genes by locus tag), so S8's consensus flank set has nothing to match. That is a limitation of the comparison, not a negative result.

## Outputs

- `case_ranking.tsv` — every recovered locus with its annotation loss, failure mode and eligibility
- `annotation_depth.tsv` — the E5 control per genome
- `cases.tsv` — the two selected cases
- `case_exons.tsv` / `case_introns.tsv` — every aligned exon and intron with what the annotation holds there
- `annotated_models.tsv` / `block_accounting.tsv` — the annotated models on each gene, and disjoint blocks counted on both sides
- `fragment_tiling.tsv` — each annotated protein blastp-tiled onto the genome's own recovered loci
- `reading_frame.tsv` — the spliced CDS, its stops, and the stops expected under neutrality
- `junction_probes.tsv` / `boundary_concordance.tsv` / `probe_summary.tsv` — the transcript search and the exon-boundary control
- `flank_consensus_check.tsv` — whether the locus's neighbours are this paralog's consensus flanks (from S8)
- `assembly_audit.tsv` — build provenance per case
- `family_named_models.tsv` / `database_records.tsv` — what the annotation names and what the databases serve
- `figures/` — four figures (D13, D19)

