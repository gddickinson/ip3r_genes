# S24 — supplementary figures, and the figure audit

*Rendered 2026-09-09 from the committed tables (D13). Nothing in this report is hand-written.*

## 1. What this task is

Two deliverables. Six supplementary figures showing the alignments and structures the main figures rest on — none of which re-aligns or re-renders anything, each drawn from the same committed file as the main figure it supports. And a figure-by-figure audit (**D11**): every main and Extended Data figure opened and read against its own legend, because errors of the kind *the legend says four and the figure draws three* are invisible to every table check in the project.

## 2. The two guards, run before any figure was drawn

A supplementary figure exists so a reader can check a join. Both joins were therefore checked first, in code, as hard failures.

| guard | what it requires | scale |
|---|---|---|
| the trimAl column map | every trimmed column must be the input column the map names, for every sequence | 1,797 columns x 134 sequences |
| the residue at the column | every variant's reference residue must be the residue its own paralogue's table holds there, and every aligned partner the residue the *other* paralogue's table holds | 1,780 variants, 2,699 partners |

Both passed. trimAl writes no column map of its own; S6 recovered one with `-colnumbering` and S24 does not trust it — the walk is exhaustive and an off-by-one would have no other symptom, because every column would still map to a column.

17 constructed negative controls run before anything is written (`scripts/s24_test_supp.py`), and 5 deliberate rule breakages were mutation-tested against them; all 5 were caught. Two of them were not, at first, and both are recorded in the test. The duplicate-column case has to be built where the two input columns hold the *same* residues, or the content walk catches it first and the one-to-one guard is never exercised. And a byte-identity check on a saved figure passes vacuously whenever both saves land in the same second, so the property is tested directly instead: the pdf must carry no creation timestamp.

## 3. The six supplementary figures

| # | figure | what it shows | drawn from |
|---|---|---|---|
| S1 | `SuppFig1_representative_alignment` | The representative alignment, and the columns the tree actually saw | `aln.fasta`, `trimmed.fasta`, `column_map.tsv`, `representatives.tsv` |
| S2 | `SuppFig2_labelled_positions` | The ligand core and the pore module at residue resolution across the three paralogues | `variants.tsv`, `paralog_variant_positions.tsv`, `constraint_ITPR1_Q14643.tsv`, `constraint_ITPR2_Q14571.tsv`, `constraint_ITPR3_Q14573.tsv` |
| S3 | `SuppFig3_paralog_alignments` | The within-paralogue alignments the constraint map is computed on | `ortholog_shape.tsv`, `layer_sizes.tsv`, `orthologs_manifest.tsv` |
| S4 | `SuppFig4_codon_alignment` | The trimmed codon alignment behind every omega estimate | `codon_aln.fasta`, `codon_trimmed.fasta`, `tip_codes.tsv`, `cds_status.tsv` |
| S5 | `SuppFig5_constraint_on_channel` | The constraint map painted onto the channel | `constraint_reference_ITPR1_7LHF.pdb`, `constraint_reference_ITPR2_9YKK.pdb`, `constraint_reference_ITPR3_8TKG.pdb`, `selection_reference_ITPR3_8TKG.pdb` |
| S6 | `SuppFig6_variants_on_structure` | Every labelled variant, and the per-element enrichment test | `constraint_reference_ITPR2_9YKK.pdb`, `constraint_reference_ITPR3_8TKG.pdb`, `variant_by_element.tsv`, `variants.tsv` |

**Supplementary Fig. S1** — the L-INS-i alignment is 11,777 columns and the tree, the selection tests and the constraint map were all computed on the 1,797 trimAl kept. What trimAl removed was the sparse columns: median occupancy 0.9776 kept against 0.0299 cut.

**Supplementary Fig. S2** — 22 pathogenic positions fall in the ligand core or the pore module, and 20 of them carry the same residue in all three paralogues. Every letter in the panel was checked against the paralogue's own per-residue table.

**Supplementary Fig. S3** — the deep layers run 249–265 orthologues per paralogue against msa_v2's 13–19. 779 sequences went through the shape screen and 1 was dropped.

**Supplementary Fig. S4** — 2,459 of 3,253 codons survive trimming, and that is the alignment every ω in the paper was estimated on.

**Supplementary Fig. S5** — the constraint map as S17 painted it, on all three cryo-EM references, with the same structure painted with the selection layer beside it. Unscored residues are grey, never the low end of the scale: S17 wrote them −1 precisely so that *we could not score this* and *this is the least conserved part of the receptor* cannot look the same, and the luminal loop is both.

**Supplementary Fig. S6** — a structure may carry a human variant position only if *every* residue it shares with the human per-residue table carries the same amino acid. ITPR2, ITPR3 pass; `constraint_reference_ITPR1_7LHF.pdb`, `constraint_model_Q14643.pdb` do not, so *ITPR1*'s pathogenic positions are not drawn on coordinates that are not theirs. The per-element enrichment test runs over 15 element × gene cells, of which 4 clear q < 0.05 after Benjamini–Hochberg.

## 4. The figure audit (D11)

23 manuscript figures — 7 main and 14 Extended Data — were opened and read against their own legends. The mechanical half is checked in code on every build: every figure must have a legend, every legend a figure, and for each Extended Data figure the panel letters its legend uses must match the panel files the manifest holds.

26 findings, 10 fixed figure, 16 fixed legend.

| figure | the legend said | the figure shows | fixed |
|---|---|---|---|
| Fig. 2 (key) | the two lighter blues are a locus found in an assembly whose annotation misses it or that has no gene set; grey is a locus the assembly is too fragmented to place | four blues, not two — found+annotated, found unannotated, found with no gene set, and partial locus — and the two greys are assembly gap and ambiguous trace; the fragmentary class is the palest blue, not the grey | fixed legend |
| Fig. 2 (whole) | evidence class ... as a fraction of the genomes of each vertebrate class | six of the scope's thirteen classes, 302 of 309 genomes; the figure keeps only classes with at least two genomes | fixed legend |
| Fig. 2 (key) | the absent colour appears nowhere at this scale: four cells (both cyclostomes' ITPR2 and ITPR3) are treated as paralogue-unassignable rather than absent | the ledger does hold those four cells as absent; they are not drawn because Myxini and Hyperoartia have one genome each and the panel drops classes below two. Paralogue-unassignable is S15's re-statement of them, not this figure's | fixed legend |
| Fig. 2 (whole) | Aves carry the most non-blue area | Lepidosauria does, on the fraction the panel plots (0.27 mean against Aves' 0.23); Aves carries the most non-blue cells because it is the largest class | fixed legend |
| Fig. 3 (backbone) | the four 100/100 labels on the backbone are the extended paralogue clades | five labels: the two ancestors of each boxed clade, which are the three extended paralogue clades, an inner ITPR1 node of 15 tips, and the 32-tip node uniting ITPR2 and ITPR3 | fixed legend |
| Fig. 6 (a) | (no panel is referenced) | a bold panel letter 'a' with no panel b anywhere in the figure | fixed figure |
| Fig. 1 (grouping) | lineages grouped by kingdom | the sweep's taxonomic groups, three of which — SAR, Discoba, Amoebozoa — are not kingdoms | fixed legend |
| Extended Data Fig. 1 (a) | the long right tail is three genomes | nine genomes carry more than three complete gene models, and four carry six or more | fixed legend |
| Extended Data Fig. 1 (b) | the grey bar reaching the light bar in every row is the result | the light bar is hidden behind the grey one in every row, so the coincidence the panel exists to show cannot be seen | fixed figure |
| Extended Data Fig. 1 (c) | every record sits between 20 % and 46 % | 19.9 % to 45.8 %, and a second dashed line at 80 % that the legend does not explain | fixed legend |
| Extended Data Fig. 2 (c) | against three to six copies for the control | the ryanodine control's most common count is two (33 % of genomes), then three (30 %) and six (21 %) | fixed legend |
| Extended Data Fig. 3 (a) | per-column conservation | a rolling mean of 25 columns, and an unexplained dashed mean line | fixed legend |
| Extended Data Fig. 7 (a) | each calibrated node's published age spread drawn as a band | a band on the two nodes a duplication is placed at, and on no other | fixed legend |
| Extended Data Fig. 7 (c) | each of the 51–53 implied losses | 51 to 102: the two support-collapsed variants imply 83 and 102 | fixed legend |
| Extended Data Fig. 8 (a) | (the second panel is not described) | two panels: the calibration, and how far each undecided cell's gene is spread across contigs | fixed legend |
| Extended Data Fig. 8 (c) | (not a legend fault) | the two 'below the caller's floor' annotations overlap into unreadable text | fixed figure |
| Extended Data Fig. 8 (d) | (not a legend fault) | the two heat maps carry independent colour scales, so the single manufactured loss in panel a is drawn as dark as the 45 in panel b | fixed figure |
| Extended Data Fig. 10 (a) | the 29-structure panel | 30 bars: 29 usable structures plus the ITPR2 record AlphaFold DB serves as a 181-residue isoform | fixed legend |
| Extended Data Fig. 10 (b) | (marker shape is not mentioned) | triangles, circles and squares carry meaning that appears in no key | fixed figure |
| Extended Data Fig. 11 (a) | each protein's own linker mean drawn | one dashed line per panel — the mean over all three proteins' linker controls, not three lines | fixed legend |
| Extended Data Fig. 11 (c) | (not a legend fault) | the four AUCs are printed at mixed precision (0.758 beside 0.6842) | fixed figure |
| Extended Data Fig. 14 (c) | almost entirely below 1,000 residues, except in the protists | the fungal 1,000–1,999 aa bar is 0.58 and the plant ≥ 2,000 aa bar 0.16; both rest on a dozen to two dozen records | fixed legend |
| Extended Data Fig. 7 (d) | (not a legend fault) | the panel title read 'every one inside its null' when two of the six loci made no call at all, and the right-hand support numbers were unlabelled | fixed figure |
| Supplementary Fig. 1 (c) | (not a legend fault) | the palest grey in the quality scale vanished at printed size in the assembled PDF, so one of the two histograms was invisible on the page | fixed figure |
| every figure in the project (pdf output) | (not a legend fault) | matplotlib stamps the wall clock into a PDF's /CreationDate, so a figure rebuilt from the same code on the same data differed from itself by two bytes and no SHA-256 recorded against a pdf meant anything | fixed figure |
| Fig. 7, Extended Data Figs 9, 12, 13 (letters) | (**A**) / (**B**) | uppercase panel letters in four figures where every other figure in the set uses lowercase | fixed figure |

Every numeric finding was re-derived from the committed table named in its `verified_against` column before it was written down; the table is `figure_findings.tsv`.

## 5. What this changes, and what it does not

**It changes the manuscript text.** Every legend finding above is corrected in `manuscript/11_figure_legends.md` or `manuscript/12_extended_data.md`, and the corrections are numbers, not wording: a legend that said four backbone labels where the tree draws five, a legend that named the wrong colour for the fragmentary class, a legend that claimed a range of 51–53 for a panel that plots up to 102.

**It changes 10 things about the figures themselves, against 16 legend corrections.** Panel letters were uppercase in four figures and lowercase everywhere else; a lone panel letter `a` sat on a figure with no panel b; two annotations overlapped into unreadable text; two heat maps side by side carried independent colour scales; a bar the legend asked the reader to compare against was hidden behind the bar it was being compared with; a threshold was drawn and never named; marker shape carried meaning that appeared in no key; a panel title asserted an agreement two of its six rows did not reach; four AUCs were printed at mixed precision; and one colour that read clearly on screen vanished at printed size.

**And it makes every figure in the project byte-reproducible.** matplotlib stamps the wall clock into a PDF's creation date, so a figure rebuilt from the same code on the same data differed from itself and no checksum recorded against a pdf meant anything. `figstyle.save` now drops the field; twelve of twelve S24 files rebuild byte-identically, and the property is a self-test rather than a claim.

**It does not change a result.** No table was recomputed and no analysis was re-run. Every fix is to a legend or to how a committed number is drawn.

**What it leaves open.** The audit is an inspection, and an inspection finds what it looks for. The mechanical checks that now run on every build cover existence, pairing and panel counts; whether a legend's *description* matches what a panel draws is not checkable in code and will need looking at again whenever a figure or a legend changes — which is what D11 says.

