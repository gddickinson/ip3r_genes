# S18 — Annotation-quality audit

*Rendered from the committed tables by `scripts/s18_report.py` (D13). Nothing here is hand-written.*

**How the family is recorded, and by whom.** Across 309 vertebrate assemblies the annotation delivers 73.9 % of the ITPR loci this project recovers as a single complete gene model — and the ryanodine-receptor control, in the same assemblies and through the same pipelines, sits at 77.9 %. The family is not recorded worse than its sister; vertebrate gene sets are recorded this way. Two things do separate: **which archive serves the annotation**, and **whether the assembly can carry the gene at all**.

## 1. Scope — what is being audited, and against what

Two records of the same genes are read against each other, and both against one common evidence set: the alignments the S5 genomic sweep placed in every genome of the declared scope.

- **The genome half.** 2,144 gene-scale loci in 309 assemblies — 1,059 ITPR loci and 942 ryanodine-receptor control loci — scored against the assemblies' own annotations. 255 sit in assemblies that ship no gene set at all and 0 could not have their own alignment recovered from the archive; both are reported as unscorable rather than as failures, leaving 932 scorable ITPR loci.
- **The protein half.** 11,402 full-length family protein records from census v6 (≥ 2,000 aa, not flagged fragments, from a protein database rather than a genome model), of which 8,306 are vertebrate and can be asked their paralog.
- **The proteome half.** the 15 vertebrate reference proteomes S3's profile sweep found nothing in, each resolved against an assembly of its own species.

The ryanodine receptors are the control throughout, not a nuisance: they are annotated by the same pipelines in the same assemblies and carry every diagnostic domain of the family (D14), so a failure rate measured without them cannot be told from the general quality of vertebrate gene sets.

### 1.1 Loci by cell and annotation source

| cell | loci | no gene set | complete | split | fragmentary | non-coding only | unannotated |
|---|---|---|---|---|---|---|---|
| ITPR1 | 401 | 52 | 260 | 6 | 53 | 16 | 11 |
| ITPR2 | 323 | 33 | 196 | 2 | 75 | 2 | 13 |
| ITPR3 | 335 | 35 | 233 | 6 | 21 | 24 | 14 |
| RYR | 1085 | 135 | 734 | 13 | 115 | 16 | 64 |

## 2. The instrument, and the three rules that make it a measurement

**Same strand only.** A gene on the other strand overlapping an ITPR locus is not a model of this gene however much sequence it shares with it. The rule is a filter, not a tiebreak, and the self-test constructs a perfectly covering antisense gene and requires the locus to read `unannotated`.

**Scored against CDS blocks, never gene spans.** A gene span covers its own introns, so scoring against spans credits an annotation with every base it never called — and this family's introns reach 152 kb (S5's intron calibration), which is room for several passenger genes. S10 established this at two loci; here it runs at 2,144.

**The name verdict reads the model the annotation actually places there**, which is not always the biggest coding one. *Podiceps cristatus* files its ITPR3 as `gene_biotype=other`, `Note=contains frameshift`, named "Inositol 1,4,5-trisphosphate receptor type 3", with no CDS feature anywhere — correctly identified and serving no protein. Reading the largest coding model would have scored that name `absent`, and those are two different failures a correction list has to keep apart.

### 2.1 Parameters

| parameter | value | where it comes from |
|---|---|---|
| completeness bar | 0.5 | `s5_classify.ANNOT_CDS_FRAC` — inherited, not chosen (§3) |
| a coding model counts as a *piece* at | 0.05 | stated here; swept in `state_sensitivity.tsv` |
| anything at all on the strand | 0.01 | separates `unannotated` from `noncoding` |
| calibration admission: alignment coverage | 0.9 | the gene is demonstrably there |
| protein-side family margin | 0.1 | `s3_assign.REL_MARGIN` (D7) |
| protein-side length floor | 2,000 aa | `family.MIN_LENGTH_AA` |
| a rename needs | 200 bits | stated here; D7's relative margin is not enough on its own (§8) |
| annotation window pad | 20,000 bp | a model may start outside the alignment's first exon |

## 3. The bar is inherited, and then measured

This project already has one definition of "this annotated gene is the model of that locus" — `s5_classify.ANNOT_CDS_FRAC`, half the alignment's coding footprint — and S5, S10 and S23 all read a locus through it. A second bar chosen here would fork what the project means by an annotated gene, so S18 uses that one and spends its calibration on the question that is open: **is 0.50 the right place for it?**

The population is the 1,077 loci whose own annotation names the correct paralog with a coding model, whose alignment recovered at least 90 % of its bait, and whose contig clears D4's bar — genes the annotation demonstrably holds in assemblies that demonstrably carry them. Over that population a single annotated model covers a median **0.993** of the gene.

| quantile of best-single-model coverage | value |
|---|---|
| minimum | 0.0032 |
| 1 % | 0.2657 |
| 5 % | 0.9542 |
| 10 % | 0.9725 |
| 25 % | 0.9826 |
| median | 0.993 |
| maximum | 1 |

So the inherited bar sits at that distribution's **1.3 % point**: 14 of 1,077 demonstrably-annotated loci fall below it (1.3 %), and it is therefore *conservative* — it credits an annotation with a gene it delivers only half of. 5 loci sit above the bar but below 0.8, which means a fifth of the gene reaches no model at all and the locus still reads `complete`. Those are reported as the price of the bar rather than absorbed by it.

### 3.1 What moving the bar changes

Every count in the grid is the same evidence read at a different bar (S15b's sensitivity matrix, applied to the one threshold S18 owns). `unannotated` and `noncoding` do not move with it by construction, which is worth being able to see rather than assert.

| bar | ITPR loci | complete | split | fragmentary | non-coding only | unannotated |
|---|---|---|---|---|---|---|
| 0.3 | 932 | 712 | 8 | 132 | 42 | 38 |
| 0.5 | 932 | 689 | 14 | 149 | 42 | 38 |
| 0.6 | 932 | 680 | 6 | 166 | 42 | 38 |
| 0.7 | 932 | 674 | 8 | 170 | 42 | 38 |
| 0.8 | 932 | 668 | 5 | 179 | 42 | 38 |
| 0.9 | 932 | 652 | 3 | 197 | 42 | 38 |
| 0.95 | 932 | 635 | 2 | 215 | 42 | 38 |

Across the whole plausible range — 0.3 to 0.95 — the share of ITPR loci reading `complete` moves from 76.4 % to 68.1 %. The audit's headline does not rest on where the bar is put.

## 4. Negative controls

`scripts/s18_test_audit.py` — all negative controls pass (s18_test_audit.py). It runs before anything is written and the driver refuses to continue if it fails. Every rule in this task returns a plausible number when it is wrong, and three of them did on the way here, so the tests are mostly tests on **refusal** and on **reachability**:

- a perfectly covering gene on the **other strand** must leave the locus `unannotated`;
- a passenger gene whose span crosses the whole locus while its CDS sits entirely in an intron must count for nothing;
- **each of the five states must be reachable**, including the three that override every model (`no_gene_set`, `gff_unavailable`, `cds_unavailable`) — a state nothing can reach is not a negative result;
- a locus whose own alignment cannot be recovered must be reported **unscorable and never scored against its span** — miniprot restarts its model ids per chunk, so a chunked genome's concatenated GFF repeats every one, and reading the raw id missed 21 loci whose denominator then became the whole locus span (2.4 Mb under an 8 kb gene, which forces `unannotated` whatever the annotation holds);
- a 6 bp overlap is not a *piece*, so one model plus a sliver is not a split;
- the naming model must prefer a small correctly-named gene over a large unnamed one, and a correctly-named non-coding gene over an unnamed coding one;
- every name verdict must fire on its own case, including the two that exist because their absence manufactured errors — `family_ambiguous` ("RyR/IP3R Homology associated domain-containing protein" names both families) and `paralog_unspecified` ("inositol 1,4,5-trisphosphate receptor" with no type number is not a wrong paralog);
- a paralog the bait panel cannot reach must not read as a naming error;
- two baits inside D7's margin must give **no** family call;
- the paralog question must not be asked outside the vertebrates;
- the calibration must refuse to report from fewer than 50 loci, and a locus the annotation is silent about must never enter it;
- **D6's veto must be recorded, not dropped** — a lesion-rich locus is written with `vetoed = 1` and its reason;
- a rename must need an absolute score as well as a relative margin;
- all four zero-hit verdicts must fire, `undecidable_no_genome` included;
- and the two one-pass readers this task added to `s10_gff` must agree with the per-window and per-model ones they replace.

Mutation-tested on nine deliberate rule breakages — the strand filter removed, scoring moved to gene spans, the naming model ranked by size, the ambiguity check dropped, the panel-reach guard disabled, the piece floor removed, the paralog asked of everything, the calibration floor removed, and `paralog_unspecified` folded back into the wrong-paralog cell. All nine were caught.

## 5. What the annotation does with the gene

Of the 932 ITPR loci in assemblies that ship a gene set, **73.9 %** are delivered as one complete annotated model and **26.1 %** are not. The failures are not evenly shaped: 38 loci have no same-strand annotated feature at all and 42 are held only by a non-coding one, so the database serves no protein for either.

*(Figure `s18_fig3_family_vs_control` panel A.)*

| cell | loci | complete | split | fragmentary | non-coding only | unannotated | complete (share) |
|---|---|---|---|---|---|---|---|
| ITPR1 | 346 | 260 | 6 | 53 | 16 | 11 | 0.7514 |
| ITPR2 | 288 | 196 | 2 | 75 | 2 | 13 | 0.6806 |
| ITPR3 | 298 | 233 | 6 | 21 | 24 | 14 | 0.7819 |
| RYR | 942 | 734 | 13 | 115 | 16 | 64 | 0.7792 |

The single largest determinant is not the paralog but the assembly: above D4's contiguity bar the ITPR failure rate falls from 26.1 % to **6.7 %**. Two thirds of what looks like an annotation problem is a contig too short to hold a 2,700-residue gene.

### 5.1 By vertebrate class

The best and worst three of the 5 classes with at least 20 scorable loci:

| class | loci | complete | fragmentary | unannotated |
|---|---|---|---|---|
| Chondrichthyes | 33 | 1 | 0 | 0 |
| Mammalia | 95 | 0.9474 | 0.0316 | 0.0105 |
| Actinopteri | 317 | 0.8675 | 0.0599 | 0.0379 |
| Amphibia | 27 | 0.7778 | 0.1111 | 0.037 |
| Aves | 433 | 0.5612 | 0.2864 | 0.0554 |

## 6. Who produced the annotation (D9)

D9 says a RefSeq gene set and a submitter-deposited GenBank one are not comparable evidence. They are not close.

Over every scorable locus, **98.8 %** of the 1,175 loci in RefSeq (`GCF_`) annotations are complete against **37.5 %** of the 699 in submitter GenBank (`GCA_`) ones. Every state differs at q ≤ 1.000.

*(Figure `s18_fig1_by_source`.)*

| state | RefSeq n | RefSeq | GenBank n | GenBank | q (BH) |
|---|---|---|---|---|---|
| cds_unavailable | 0 | 0 | 0 | 0 | 1 |
| unannotated | 10 | 0.0085 | 92 | 0.1316 | 0 |
| noncoding | 2 | 0.0017 | 56 | 0.0801 | 0 |
| complete | 1161 | 0.9881 | 262 | 0.3748 | 0 |
| split | 1 | 0.0009 | 26 | 0.0372 | 0 |
| fragmentary | 1 | 0.0009 | 263 | 0.3763 | 0 |

### 6.1 The confounder, controlled

GenBank assemblies are less contiguous, and a locus on a contig too short to hold the gene cannot be annotated completely by anybody — so the contrast is run again over the loci that clear D4's bar. It narrows and does not close: **99.4 %** against **63.8 %** on 1,015 and 185 loci (q = < 1e-300). About 42.1 % of the gap between the two archives is assembly quality; the rest is the gene set.

| state | RefSeq n | RefSeq | GenBank n | GenBank | q (BH) |
|---|---|---|---|---|---|
| cds_unavailable | 0 | 0 | 0 | 0 | 1 |
| unannotated | 3 | 0.003 | 9 | 0.0486 | 7e-06 |
| noncoding | 1 | 0.001 | 8 | 0.0432 | 3e-06 |
| complete | 1009 | 0.9941 | 118 | 0.6378 | 0 |
| split | 1 | 0.001 | 18 | 0.0973 | 0 |
| fragmentary | 1 | 0.001 | 32 | 0.173 | 0 |

## 7. Is it this family, or is it vertebrate annotation?

The audit's premise is that a ~2,700-residue, ~58-exon gene with a sister family sharing every diagnostic domain should be badly recorded. The control says otherwise, and that is the result.

The ITPR cells fail on **26.1 %** of 932 loci; the ryanodine receptors, in the same assemblies through the same pipelines, on **22.1 %** of 942 (q = 0.128). Above D4's bar the two are 6.7 % and 5.5 % (q = 0.819). **No overall difference survives correction.** How this family is recorded is how vertebrate genes of this size are recorded.

*(Figure `s18_fig3_family_vs_control` panel B.)*

| measure | ITPR n | ITPR | RyR n | RyR | q (BH) |
|---|---|---|---|---|---|
| any annotation failure | 243 | 0.2607 | 208 | 0.2208 | 0.127691 |
| cds_unavailable | 0 | 0 | 0 | 0 | 1 |
| unannotated | 38 | 0.0408 | 64 | 0.0679 | 0.074942 |
| noncoding | 42 | 0.0451 | 16 | 0.017 | 0.006342 |
| complete | 689 | 0.7393 | 734 | 0.7792 | 0.127691 |
| split | 14 | 0.015 | 13 | 0.0138 | 1 |
| fragmentary | 149 | 0.1599 | 115 | 0.1221 | 0.093763 |

One state does separate, and it is the family-specific one: an ITPR locus is **2.7×** more likely than a RyR locus to be held only by a non-coding feature — 42 loci against 16 (q = 0.006). That is the *Podiceps* failure at scale: a gene the annotation identifies correctly, names correctly, and files as non-coding, so no protein record is ever created and no name-based search can reach it.

**And it does not survive the contiguity control.** Above D4's bar the same comparison is 4 against 5 on 568 and 632 loci (q = 1.000) — the effect is confined to assemblies too broken to carry the gene, which is where a submitter has most reason to file a model as non-coding in the first place. The excess is real in the record set and this measurement cannot say it is a fact about the family rather than about the assemblies the family's loci happen to sit in.

## 8. What the protein databases call these records

Each of the 11,402 full-length family protein records was scored against the committed 38-bait panel by blastp: best bit score in each family, assigned to the winner only when it beats the loser by more than D7's relative margin, then the paralog inside the winning family — and the record's own gene symbol and protein name read through the *same* verdict rule the genome half uses.

**The family call is not in dispute.** The sequence disagrees with the census call on 4 of 11,402 records, and the databases name the sister family at only 5 — every one of them a non-vertebrate record whose best family score is under 200 bits, i.e. where the sequence barely separates the families either. D14's hazard, which S0 measured at 49 % of zebrafish PF08709 records, does not appear as a *naming* error in the full-length record set.

**The paralog call is not in dispute either.** 52 of 8,306 vertebrate records carry a symbol naming a paralog the panel assigns elsewhere, and a further 74 claim RYR3 — which this panel cannot call, because S5's slot table records no RYR3 bait, so those are reported as outside the instrument's reach rather than as errors.

**What is in dispute is whether the records are findable.** 3,872 records carry a placeholder gene symbol and 2,395 carry none at all: **55.0 % of the family's full-length protein records have no usable gene symbol.** On the protein-name side 66 are named for the superfamily — "RyR/IP3R Homology associated domain-containing protein" — which names the family and its sister together and therefore separates neither.

*(Figure `s18_fig4_protein_side` panel A.)*

### 8.1 Pfam recall — would a signature query have found them?

The signature that *names* this family is PF08709. A reader looking for the IP3 receptors would query it, so the recall question is what that query returns out of the records the project's own two instruments call family.

| axis | bucket | records | carry PF08709 | complete architecture |
|---|---|---|---|---|
| census call | ITPR | 5833 | 0.9523 | 0.903 |
| census call | RYR | 5569 | 0.9196 | 0.8655 |
| which instruments called it | both | 9841 | 0.9564 | 0.9291 |
| which instruments called it | profile | 7 | 0.2857 | 0 |
| which instruments called it | v4+s20 | 1554 | 0.8121 | 0.6075 |

A PF08709 query reaches 95.6 % of the records both instruments agree on, and 81.2 % of those only one instrument found. The signature is a good but not complete index of its own family, and the deficit is concentrated outside the vertebrates.

| group | records | carry PF08709 | complete architecture |
|---|---|---|---|
| Amoebozoa | 11 | 0.6364 | 0.0909 |
| SAR | 407 | 0.6462 | 0.2457 |
| Viridiplantae | 24 | 0.7917 | 0.25 |
| Eukaryota (other) | 34 | 0.8235 | 0.4118 |
| Metazoa (non-vertebrate) | 2552 | 0.9244 | 0.8672 |
| Discoba | 42 | 0.9524 | 0.3571 |
| Vertebrata | 8306 | 0.9553 | 0.9302 |
| Fungi | 26 | 0.9615 | 0.4615 |

## 9. The proteomes that returned nothing

S3 swept 764 vertebrate reference proteomes with both profiles and 15 came back with no family hit at all. On its own that is uninterpretable — an ITPR-shaped hole in a proteome is either a gene the species lacks or a gene its gene caller did not find — so each was resolved against an assembly of its own species, which the S5 sweep searched with an instrument that owes the gene caller nothing.

**All 15 are gene-caller failures.** Every one of the 15 species has a genome in the S4 scope, and in every one the genomic sweep recovers at least one ITPR cell at over half the bait's length while the proteome holds none. Not one is `genome_also_empty`, and not one is `undecidable_no_genome` — the verdict that exists so a species with no assembly in scope could not be reported as an absence.

*(Figure `s18_fig4_protein_side` panel B.)*

| proteome | species | class | proteins | assembly | ITPR loci | cells recovered | verdict |
|---|---|---|---|---|---|---|---|
| UP001529510 | Cirrhinus mrigala | Actinopteri | 57290 | GCA_036247105.2 | 4 | 3 | gene_caller_missed_it |
| UP000053858 | Charadrius vociferus | Aves | 14257 | GCF_000708025.1 | 3 | 3 | gene_caller_missed_it |
| UP000053283 | Nipponia nippon | Aves | 14132 | GCF_047370975.1 | 3 | 3 | gene_caller_missed_it |
| UP000053286 | Aptenodytes forsteri | Aves | 13704 | GCF_000699145.1 | 3 | 3 | gene_caller_missed_it |
| UP000053760 | Cuculus canorus | Aves | 13582 | GCF_017976375.1 | 3 | 3 | gene_caller_missed_it |
| UP000536381 | Semnornis frantzii | Aves | 10801 | GCA_013399775.1 | 1 | 1 | gene_caller_missed_it |
| UP000537747 | Mystacornis crossleyi | Aves | 10563 | GCA_013400655.1 | 2 | 1 | gene_caller_missed_it |
| UP000660247 | Todus mexicanus | Aves | 10042 | GCA_013389965.1 | 1 | 1 | gene_caller_missed_it |
| UP000053537 | Acanthisitta chloris | Aves | 9650 | GCF_000695815.1 | 3 | 3 | gene_caller_missed_it |
| UP000053238 | Phalacrocorax carbo | Aves | 8990 | GCF_963921805.1 | 3 | 3 | gene_caller_missed_it |
| UP000295264 | Sousa chinensis | Mammalia | 8382 | GCA_007760645.1 | 3 | 3 | gene_caller_missed_it |
| UP000054232 | Eurypyga helias | Aves | 7979 | GCF_000690775.1 | 2 | 2 | gene_caller_missed_it |
| UP000031515 | Chaetura pelagica | Aves | 3162 | GCF_000747805.2 | 3 | 3 | gene_caller_missed_it |
| UP000631465 | Cervus hanglu | Mammalia | 713 | GCA_010411085.1 | 3 | 3 | gene_caller_missed_it |
| UP000681084 | Bothrops jararaca | Lepidosauria | 84 | GCA_018340635.1 | 2 | 2 | gene_caller_missed_it |

## 10. The correction list

297 correction items, 52 of them `high` priority, 18 withheld under **D6**. Each row carries the assembly, the coordinates, the current state, the current name, the proposal, the archived evidence file a curator can open, and the numbers the class fired on. → `results/annotation_audit/corrections.tsv`

**Priority is a rule, not an impression.** `high` needs the gene demonstrably present *and* the assembly demonstrably able to carry it *and* the reading frame intact; a partial recovery or an unscored ORF is `medium`; anything below D4's contiguity bar is `low`, because there the annotation's silence may be the assembly's fault and not the annotator's.

**D6 is applied as a column, not a filter.** If S15's ORF screen scored a locus lesion-rich, the audit does not tell RefSeq to resurrect it — but the row is still written, flagged, with its reason, because a locus this audit declined to correct is evidence about the audit.

| class | priority | items | vetoed (D6) | RefSeq | GenBank |
|---|---|---|---|---|---|
| C1_unannotated | high | 10 | 0 | 3 | 7 |
| C1_unannotated | medium | 2 | 1 | 0 | 2 |
| C1_unannotated | low | 8 | 0 | 0 | 8 |
| C2_noncoding_only | high | 7 | 0 | 0 | 7 |
| C2_noncoding_only | medium | 1 | 1 | 1 | 0 |
| C2_noncoding_only | low | 46 | 0 | 0 | 46 |
| C3_incomplete | high | 34 | 0 | 1 | 33 |
| C3_incomplete | medium | 16 | 16 | 1 | 15 |
| C3_incomplete | low | 167 | 0 | 0 | 167 |
| C5_wrong_paralog_name | high | 1 | 0 | 0 | 1 |
| C6_protein_wrong_family | low | 5 | 0 | 0 | 0 |

### 10.1 The first ten `high`-priority items

Sorted as the table is written; the full list is the file.

| assembly | species | cell | locus | state | current name | proposal |
|---|---|---|---|---|---|---|
| GCA_011823955.1 | Dissostichus mawsoni | ITPR1 | JAAKFY010000006.1 | unannotated | (none) | add a protein-coding gene model over the 8228 bp of aligned coding sequence on JAAKFY010000006.1:25372485-25435846 (-) |
| GCA_011823955.1 | Dissostichus mawsoni | ITPR2 | JAAKFY010000007.1 | unannotated | (none) | add a protein-coding gene model over the 8065 bp of aligned coding sequence on JAAKFY010000007.1:17781036-17839072 (-) |
| GCA_013398455.1 | Loxia curvirostra | ITPR1 | VZSM01004910.1 | fragmentary | Itpr1 | 1 coding model(s) deliver 0.03 of the gene, best alone 0.03; extend or merge to the aligned exons |
| GCA_013398455.1 | Loxia curvirostra | ITPR2 | VZSM01004385.1 | fragmentary | Itpr2_1 | 2 coding model(s) deliver 0.14 of the gene, best alone 0.11; extend or merge to the aligned exons |
| GCA_013398455.1 | Loxia curvirostra | ITPR3 | VZSM01006023.1 | noncoding | LOXCUR_R12629 | the locus is held only by 1 non-coding feature(s) (pseudogene); the alignment gives a coding model over the same exons |
| GCA_013398455.1 | Loxia curvirostra | RYR | VZSM01002731.1 | fragmentary | Ryr3_1 | 4 coding model(s) deliver 0.08 of the gene, best alone 0.03; extend or merge to the aligned exons |
| GCA_014281875.1 | Nibea albiflora | ITPR1 | CM024790.1 | fragmentary | ITPR1 | 1 coding model(s) deliver 0.02 of the gene, best alone 0.02; extend or merge to the aligned exons |
| GCA_014281875.1 | Nibea albiflora | ITPR1 | CM024790.1 | unannotated | (none) | add a protein-coding gene model over the 8250 bp of aligned coding sequence on CM024790.1:3836243-3913687 (+) |
| GCA_014281875.1 | Nibea albiflora | ITPR2 | CM024800.1 | unannotated | (none) | add a protein-coding gene model over the 8021 bp of aligned coding sequence on CM024800.1:12447838-12500306 (-) |
| GCA_014281875.1 | Nibea albiflora | ITPR3 | CM024793.1 | noncoding | ITPR2 | the covering model claims ITPR2; the alignment assigns ITPR3 |

## 11. The priors this task was judged against

Each is stated with where the earlier task said it, computed on S18's own tables, and rendered with both numbers printed either way.

- **orthogonal.** S10 §3 — the annotation gets the gene right at 375 of 382 loci (98.2 %) **where it demonstrably could**: 880 recovered loci, 382 passing E1-E5, median annotation loss 0, 360 at exactly 0. S18 measures **73.9 % complete over every scorable locus, and 93.3 % over the loci above D4's bar**.
  S10's denominator is the loci where the annotation *demonstrably could* have delivered the gene (five eligibility rules, of which E5 asks whether that annotation builds genes this long anywhere else in the genome); S18's is every locus in the scope. The gap between 98.2 % and 73.9 % is the size of the population S10 excluded, and is a result rather than a disagreement.
- **confirmed.** S5b — census v4 adds 1,058 gene models from 224 genomes, of which **318 ITPR models exist only as DNA** and 167 more sit inside an annotated gene carrying no family name, unreachable by any name-based search. S18 measures **38 ITPR loci with no same-strand annotated feature and 42 held only by a non-coding one, in the 932 loci whose assembly ships a gene set**.
- **contradicted.** S5b — holding contiguity constant across 487 loci, ITPR3 is 88 % correctly named against ITPR1's 65 %, a 23-point gap between genes of near-identical protein length in the same genomes. Uncontrolled the gap looked like 54 % vs 14 %, which was assembly quality. S18 measures **a 1.8-point spread across ITPR1/2/3 (ITPR1 90.9%, ITPR2 91.7%, ITPR3 89.9%) once contiguity is held constant**.
  S5b measured naming against its own `annot_paralog_matches`, which requires a model covering half the locus; S18 reads the name off whichever model the annotation places there, coding or not. On that reading the paralogs are named equally well and S5b's 23-point gap does not survive — what differs between the paralogs is whether a *coding* model exists, not whether the gene is named.
- **confirmed.** S5b — 3 loci are claimed by a cell other than the one their annotation names (1 with a coherent sibling locus), **supplied to S18, not adjudicated**. S18 measures **3 loci whose covering model names a different paralog, of which 1 clears the recovery bar a correction needs**.
  The other 2 are recovered at under 0.11 of the bait, where an alignment covering a tenth of the gene is not evidence for renaming the model that covers the rest.
- **confirmed.** S3 — 15 of 764 vertebrate reference proteomes returned no family hit from either profile. S3 could not say whether that was a gene or a gene caller. S18 measures **15 of 15 resolved as gene-caller failures, 0 as absences**.
- **confirmed.** S2 — 2,911 of 15,417 seeded-space proteins do not carry PF08709, the signature that names the family, so a PF08709 query recovers 81.1 % of the space its own sister signatures enumerate. S18 measures **95.6 % of the records both instruments call family carry PF08709**.
- **orthogonal.** S0 — 49 % (53/109) of zebrafish PF08709 records are ryanodine receptors. D14's hazard, measured at the signature level, is what makes wrong-family naming the audit's sharpest question. S18 measures **5 wrong-family names in 11,402 full-length records, and 0 at any of the 2,144 genomic loci**.
  S0 measured what a *signature* returns; S18 measures what a *name* claims. A Pfam shared by both families says nothing about whether either is named correctly, and on this evidence both are.
- **confirmed.** S5b/D4 — 120 of 309 manifest genomes cannot hold a median ITPR gene on one contig (66 % of birds against 11 % of fish), so any annotation failure rate has to be reported with that control. S18 measures **the ITPR failure rate falls from 26.1 % to 6.7 % across that bar**.

## 12. What this settles, and what it does not

**Settled.**

- How the family is recorded across the declared genome scope, with the sister family measured beside it in the same assemblies: 73.9 % of ITPR loci delivered as one complete model against 77.9 % of RyR loci, no difference surviving correction.
- That the archive matters far more than the family: RefSeq and submitter-deposited GenBank gene sets differ on every state, and only part of that gap is assembly quality.
- That an ITPR-shaped hole in a vertebrate reference proteome is, in every case in this scope, a gene caller and not a gene.
- That the family's naming in the protein databases is essentially correct where a name exists, and that a name often does not: over half the full-length records have no usable gene symbol.

**Not settled.**

- **Whether a proposed correction is right.** The audit proposes; it does not validate. S10 validated two cases to the exon and found both; the remaining 297 items here carry alignment evidence and nothing more.
- **The RYR3 cell.** The bait panel has no RYR3 bait, so 74 records whose symbol says RYR3 are outside the instrument's reach rather than adjudicated. That is a limit of S5's panel, recorded as one.
- **Whether the non-coding demotions are wrong.** An annotation that files a gene as `gene_biotype=other` with a note about a frameshift may be right about the frameshift. S15's ORF screen vetoes the ones it can; the rest are proposed with that stated.
- **Anything about the 33 assemblies with no gene set**, which are reported as unscorable rather than as failures.

**What a reader should hold against it.** The completeness bar is inherited from S5 and sits in the far lower tail of the distribution it is applied to, so `complete` is a generous verdict — 15 loci read complete while a fifth of the gene reaches no model. And the whole genome half rests on the S5 alignments being right about where the exons are; S10's exon-boundary control validated that instrument at two loci in two genomes, not at 2,144.

