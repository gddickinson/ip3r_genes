# S8 — Synteny of the ITPR loci across 309 genomes

*Generated 2026-09-07 06:35; 10 flanking genes per side, `fixed10` as the primary window and `informative10` for the caller, 3 control windows per genome (seed 20260907), runtime 56.9 s.*

**What this task asks.** S5 assigned a paralog to each locus from sequence: which bait won the alignment. Synteny asks the same question with evidence the alignment never sees — what the *neighbouring* genes are called — so that an orthology claim does not rest twice on the same measurement.

## 1. What is in the analysis

Loci come from the per-genome `summary.json` files rather than from the ledger, because the ledger carries one row per genome × cell and therefore only the best locus. A teleost ITPR1 cell holds *itpr1a* and *itpr1b*; the sea lamprey's holds three. Reading the best one would compare *itpr1a* in one species against *itpr1b* in the next and report the mismatch as a synteny result.

| cell | loci | of which extra copies in the same cell |
|---|---|---|
| ITPR1 | 401 | 109 |
| ITPR2 | 323 | 32 |
| ITPR3 | 335 | 38 |
| RYR | 1,085 | 777 |

- Loci with coordinates across 309 swept genomes: **2,144**
- Genomes carrying an annotation gene table: **274** (the rest are unannotated assemblies; **270** loci have no flanks and are absent from every statistic below)
- Genomes contributing control windows: **242** — a genome whose longest contig carries fewer than 21 coding genes cannot host a fair control window and is excluded from the null rather than allowed to degrade it
- Cells excluded by name: **0** (the mechanism exists so that an exclusion cannot be made silently; S7's relabelling rule fired on no tip, so none is due)

## 2. The instrument

### 2.1 Two windows, because annotation naming density is not constant

Matching flanks across species means matching gene *symbols*, and a symbol is only usable if the annotation gave the gene one. Across this genome set that varies about fourfold — 97 % of human coding genes carry an informative symbol against 27 % of the sea lamprey's — so a fixed ±10-gene window hands a well-named genome 20 usable keys and a poorly-named one 5. A Jaccard difference would then be an annotation difference.

So both windows are committed and every table carries a `window` column:

| window | rule | mean informative keys per locus |
|---|---|---|
| `fixed10` | the brief's rule — the 10 nearest coding genes each side | 9.5 |
| `informative10` | the 10 nearest *informative* symbols each side, scanning at most 60 genes | 13.1 |

`fixed10` is primary for the pair statistics because it is the brief's rule and it makes no assumption about what a symbol is; `informative10` is what the consensus caller uses, because a caller whose evidence depends on the assembly's naming conventions cannot be calibrated across them.

### 2.2 Three key vocabularies

A symbol is normalised before it is compared, and the three levels answer different questions:

| key | rule | what it can see |
|---|---|---|
| `strict` | the symbol, uppercased | the same gene, same name |
| `relaxed` | lowercase teleost/amphibian duplicate suffixes stripped (`gnasa` → `GNAS`, trailing `.2` dropped); uppercase symbols untouched, so `GNB1` and `BAK1` survive | the same gene across clades that name it differently |
| `root` | the relaxed key with its trailing digit run removed (`BHLHE40`, `BHLHE41` → `BHLHE`), floored at 3 characters so `TP53` and `C3` keep their key | the same gene **family** — which is the only level at which a 2R ohnolog pair is visible, because the two copies almost never carry the same symbol |

The root key over-merges (every zinc finger reaches `ZNF`). That is tolerable only because it is scored against a null built with the *same* rule on random windows, so over-merging inflates the signal and its background together.

### 2.3 The control: matched random neighbourhoods

A Jaccard of 0.21 means nothing on its own. Two neighbourhoods in two well-annotated mammals share vocabulary for reasons that have nothing to do with orthology, and two in a hagfish and a lamprey share almost none whatever their history.

So every real pair *(locus in genome A, locus in genome B)* is scored against matched control pairs *(random coding gene in A, random coding gene in B)* — 3 replicates per genome, drawn with the same window rule and the same key rule, seeded per accession so the draw is reproducible. The control therefore holds constant the two genomes, their annotation depth, their naming conventions, the window and the normalisation. What is left is the locus.

The paired comparison is reported as a **sign test** rather than a difference of means: Jaccard is bounded, zero-inflated and nowhere near normal, and pairs where both the locus and its control score zero are counted as ties and dropped rather than scored as failures, because counting them against the locus would make an unannotated genome look like a negative result.

### 2.4 The negative controls, run on every build

The rules that decide what a flank is, when two symbols are the same gene, and when a consensus call is evidence are all constructed-case tested before anything is measured with them (`scripts/s8_test_flanks.py`, the pattern of `s5_bait_screen.self_test()` and `s7_test_tree.py`). A failure aborts the run.

| test | what must fail |
|---|---|
| **T1** | the placeholder vocabularies — `LOC116952798`, `FN964_004414`, `ENSG…`, `si:ch211-…` — are rejected and real symbols kept |
| **T2** | `gnasa` reaches `GNAS`; `GNB1`, `BAK1`, `ctsa`, `vapb` are left alone |
| **T3** | `BHLHE40`/`BHLHE41` and `GRM7`/`GRM4` merge; `TP53` and `C3` are not reduced below the character floor |
| **T4** | the locus itself and its overlappers stay out of its own flank set, non-coding biotypes are filtered, and the two windows behave differently on the same table |
| **T5** | Jaccard is symmetric and **empty-vs-empty is 0, not 1** — otherwise two unannotated genomes would score as a perfect synteny match |
| **T6** | the pair subset takes one locus per species, the highest-covered, so a teleost *a*/*b* pair cannot be compared across species |
| **T7** | the control draw is reproducible and refuses a contig too short to hold a window |
| **T8** | the caller's two halves — a genuine neighbourhood is called, an unrelated one and a key-poor one are refused, and a locus cannot be scored against a consensus it voted into |
| **T9** | the threshold rule refuses a setting that calls everything, including 90 % of random windows |
| **T10** | a call is only `supported` above the highest score any random window reached |
| **T11** | a ranked table's row order does not depend on which order the input arrived in — Python's set iteration is hash-seeded per process, and a rank without a final tiebreak differs between two runs of identical data |

`synteny_stats.json` also records the SHA-256 of every committed table, so a rerun that drifts is visible in the data rather than only in a diff.

## 3. Does the neighbourhood carry the paralog?

### 3.1 Within a paralog against a matched random neighbourhood

One locus per species per cell (the highest-covered), `fixed10` window, `relaxed` keys. `ratio` is the mean Jaccard divided by the mean of that pair's own matched control; `beats control` is the fraction of pairs that individually exceed their own matched control pair, ties dropped.

| pair class | pairs | mean J | median J | J > 0 | control mean J | ratio | beats control | z |
|---|---|---|---|---|---|---|---|---|
| within ITPR1 | 27,261 | 0.211 | 0.059 | 0.71 | 0.0005 | 411× | 0.997 | 139 |
| within ITPR2 | 27,028 | 0.219 | 0.129 | 0.79 | 0.0005 | 413× | 0.998 | 146 |
| within ITPR3 | 27,261 | 0.125 | 0.000 | 0.44 | 0.0006 | 216× | 0.981 | 106 |
| within RYR | 28,920 | 0.096 | 0.000 | 0.33 | 0.0006 | 173× | 0.977 | 94 |
| ITPR1 vs ITPR2 | 54,293 | 0.000 | 0.000 | 0.00 | 0.0005 | 0.01× | 0.022 | -24 |
| ITPR1 vs ITPR3 | 54,529 | 0.000 | 0.000 | 0.00 | 0.0005 | 0.30× | 0.176 | -18 |
| ITPR1 vs RYR | 56,160 | 0.000 | 0.000 | 0.00 | 0.0005 | 0.29× | 0.268 | -14 |
| ITPR2 vs ITPR3 | 54,295 | 0.000 | 0.000 | 0.00 | 0.0006 | 0.00× | 0.000 | -25 |
| ITPR2 vs RYR | 55,920 | 0.000 | 0.000 | 0.00 | 0.0005 | 0.06× | 0.052 | -23 |
| ITPR3 vs RYR | 56,161 | 0.000 | 0.000 | 0.00 | 0.0006 | 0.00× | 0.000 | -26 |

Every ITPR paralog's neighbourhood is shared far beyond what two random neighbourhoods in the same two genomes share: **216×–413× the matched null**, with 98.1%–99.8% of individual pairs beating their own control. The spread across the trio is real — ITPR2 at 0.219 against ITPR3 at 0.125 — and §3.2 is where it comes from.

Prior: **812 `found_annotated` cells** — S5 — 812 of the 1,236 genome × class cells are `found_annotated`, and every downstream task treats a cell's paralog label as orthology. Nothing before S8 tested that with evidence outside the gene itself. **confirmed** — the paralog cells are backed by genomic neighbourhood, which is evidence the bait alignment never saw

### 3.2 Within a class, and across classes

A pooled within-paralog mean answers two questions at once: do two mammals share the neighbourhood (they do, nearly trivially), and does a mammal share it with a teleost. Only the second is about the locus, and a clade-restricted signal reported pooled reads as a vertebrate-wide one.

| cell | same class, mean J | different classes, mean J | ratio | cross-class J > 0 | beats control |
|---|---|---|---|---|---|
| ITPR1 | 0.323 | 0.160 | 0.50 | 0.71 | 0.997 |
| ITPR2 | 0.347 | 0.158 | 0.46 | 0.79 | 0.999 |
| ITPR3 | 0.254 | 0.062 | 0.24 | 0.32 | 0.969 |
| RYR | 0.185 | 0.053 | 0.29 | 0.30 | 0.973 |

**ITPR3's neighbourhood is the one that does not travel.** Across vertebrate classes it retains 0.062 against 0.160 and 0.158 for the other two — a 2.5-fold gap — while within a class the three span only 0.254–0.347. It is still 157× its own null and 96.9% of its cross-class pairs still beat their control, so this is decay, not absence.

Prior: **232 `found_annotated` cells for ITPR3** — S5 — ITPR3 is the most consistently recovered paralog (232 `found_annotated` cells against 194 for ITPR1 and 170 for ITPR2), and S5's annotation-quality section found the three within 6 points of each other once contiguity was held constant. **orthogonal** — ITPR3 is the *best*-recovered paralog and the one whose neighbourhood is *least* conserved. The prior stands; the two measure different properties of the same gene — how easy it is to find, and how stable the ground it sits on. In human that ground is the MHC region at 6p21

### 3.3 The sister family, as both controls at once

The ryanodine receptors ride the same genomes through the same code. They are a **positive** control — a second family the method has to work on — and a **negative** one: ITPR and RyR neighbourhoods must not be shared (D14).

| RyR pair class | pairs | mean J | ratio to null |
|---|---|---|---|
| within RYR1 | 1,770 | 0.166 | 329× |
| within RYR2 | 14,028 | 0.164 | 343× |

The pooled `within_RYR` number is lower than either of these because the RYR cell is one cell holding three genes: the best locus per genome is RYR2 in most species and RYR1 in others, so the pooled figure is a mixture of two neighbourhoods, not a measurement of one.

Across the family boundary, over **168,241** ITPR × RyR pairs, the highest mean Jaccard of any class is **0.0002** — below the random-window control itself.

Prior: **ITPR and RyR are separate families** — D14, measured at every stage since S1 — the ryanodine receptors carry every ITPR-diagnostic Pfam domain, and separating the two families is a positive test, never an assumption. S1 measured ITPR-to-RyR covered identity at 0.249 against 0.828 within the family. **confirmed** — an instrument that shares nothing with the sequence evidence returns the same separation

## 4. The 2R paralogon

ITPR1/2/3 are a 2R product, so their neighbourhoods should be *paralogous* rather than identical: the flanking genes should be the surviving copies of the same ancestral families under different names. §3.1 shows symbol Jaccard between paralogs at zero, which is what that prediction looks like to a test that is looking for the same word twice. The root key is what makes an ohnolog pair visible.

A root is reported when at least 10% of the species on **both** sides carry it. The summary counts at three bars so the answer does not rest on where one line is drawn.

| pair | shared root families | at ≥10 % | at ≥25 % | at ≥50 % | background < 1 % |
|---|---|---|---|---|---|
| ITPR1 vs ITPR2 | 1 | 1 | 1 | 1 | 1 |
| ITPR1 vs ITPR3 | 1 | 1 | 1 | 0 | 1 |
| ITPR1 vs RYR | 0 | 0 | 0 | 0 | 0 |
| ITPR2 vs ITPR3 | 0 | 0 | 0 | 0 | 0 |
| ITPR2 vs RYR | 0 | 0 | 0 | 0 | 0 |
| ITPR3 vs RYR | 0 | 0 | 0 | 0 | 0 |

Every root that passes, with the random-window background it had to beat:

| pair | root family | in side A | in side B | classes A / B | random-window background | enrichment |
|---|---|---|---|---|---|---|
| ITPR1 vs ITPR3 | `GRM` | 0.42 | 0.53 | 6 / 8 | 0.0045 | 93× |
| ITPR1 vs ITPR2 | `BHLHE` | 0.62 | 0.85 | 11 / 10 | 0.0074 | 84× |

**The surviving paralogon links run through ITPR1.** ITPR1 with ITPR2 and ITPR1 with ITPR3 each retain one shared flanking family; **ITPR2 with ITPR3 retains none at any bar**, and no ITPR × RyR pair retains one (0 of 3 cross-family pairs with any shared root).

Prior: **ITPR2 + ITPR3 are sisters** — S7 §5.2 — the AU test rejects ITPR1+ITPR2 (p-AU 1.8e-05) and ITPR1+ITPR3 (p-AU 1.65e-05) and does not reject ITPR2+ITPR3 (p-AU 0.476); the unconstrained ML tree groups ITPR2+ITPR3 at SH-aLRT 100 / UFBoot 100. **not corroborated** — the tree's sister pair is the one pair whose neighbourhoods share nothing. These are not the same measurement and one does not overturn the other: a tree estimates the order of duplication, while a retained flanking ohnolog records which copies *survived deletion* beside each gene, and 2R quartets are known to lose flank copies independently of the duplication order. What can be said is that the synteny does not corroborate it, and that whichever pair is sister, the ITPR1 neighbourhood is the one that kept its ohnologs. The measurement itself is clean — two families, both under 1 % of random windows, at every bar from 10 % to 50 %

## 5. Placing the loci the sweep could not label

### 5.1 A paralog caller built only from the neighbours

Each paralog's flank consensus is the set of keys carried by at least a stated fraction of the species holding that paralog. A locus is scored by how many consensus keys its own flanks carry, **with its own species dropped from every consensus first**, so a locus cannot be scored against evidence it supplied.

It is calibrated on the loci whose paralog identity their own assembly's annotation already establishes — evidence the caller never sees, since it reads the symbols of the *neighbours* and the truth is the symbol of the gene itself.

The consensus threshold is measured, not typed. The sweep reports what each setting buys (call rate on the confirmed loci) and what it costs (the rate at which random neighbourhoods in the same genomes are called a paralog), and the rule is to maximise their difference — one call gained is worth one false call avoided — with ties broken toward the stricter setting.

| consensus frac | consensus sizes 1/2/3 | call rate | accuracy | false-call rate | difference |
|---|---|---|---|---|---|
| 0.20 | 26/28/20 | 0.811 | 1.000 | 0.0152 | 0.7960 |
| 0.30 | 16/15/14 | 0.805 | 1.000 | 0.0083 | 0.7969 |
| 0.40 | 15/11/9 | 0.805 | 1.000 | 0.0083 | 0.7969 ← |
| 0.50 | 8/7/2 | 0.767 | 1.000 | 0.0041 | 0.7633 |
| 0.60 | 5/4/0 | 0.547 | 1.000 | 0.0014 | 0.5453 |
| 0.70 | 1/2/0 | 0.539 | 1.000 | 0.0014 | 0.5374 |

Chosen: **0.4** — maximises call rate - random-window false-call rate (0.805 - 0.008 = 0.797); accuracy on labelled loci 1.000. Accuracy is 1.000 across the whole sweep, so it separates nothing and is not what is optimised; it is reported.

At the chosen setting the caller is **405 correct of 405 calls on 503 annotation-confirmed loci (1.000), call rate 0.805**, and it calls **6 of 726 random control windows** (0.008).

| paralog | confirmed loci | called | correct | accuracy |
|---|---|---|---|---|
| ITPR1 | 136 | 121 | 121 | 1.000 |
| ITPR2 | 163 | 155 | 155 | 1.000 |
| ITPR3 | 204 | 129 | 129 | 1.000 |

The highest consensus overlap any random window reached is **2**, so a call at or below that score is reported `within_null` however clean it looks. That bar is the null's own maximum rather than a probability cut, because with 726 control windows a 1-in-726 tail is not a rate to build a claim on.

### 5.2 What the caller adds

Every locus whose paralog the sweep did not establish is scored: loci in cells with no paralog-named annotation, the extra copies in multi-copy cells, the fragment and assembly-gap cells, and the cyclostome loci.

| why it was unplaced | loci | supported call | within the null | no call |
|---|---|---|---|---|
| S5 status assembly_gap | 147 | 8 | 2 | 137 |
| S5 status fragment | 7 | 0 | 0 | 7 |
| cyclostome locus (no 2R paralog label available) | 6 | 0 | 4 | 2 |
| extra locus in a multi-copy cell | 81 | 36 | 19 | 26 |
| no paralog-named annotation at the locus | 211 | 87 | 32 | 92 |

**131 loci gain a paralog assignment that no random neighbourhood could have produced** — ITPR1 73, ITPR2 42, ITPR3 16. Each carries its score, its margin over the runner-up and its null tail probability in `unplaced_loci.tsv`.

Of these, **36** are the second and later copies inside one paralog cell — the teleost 3R co-orthologs and their relatives — and **36** of them are placed in the same paralog as the cell they were filed under.

Prior: **two loci per cell** — S6 §1 and S7 §5.5 — the teleost 3R co-orthologs are in the representative set by rule; S5's ledger files both copies into the same paralog cell, so a cell can hold two genes that are not the same gene. **confirmed** — the extra copies are the same paralog as their cell, so a cell holding two loci holds two copies of one gene rather than a misfiled second gene

### 5.3 The cyclostome loci — the question S7 handed to S8

| species | cell.copy | informative keys | ITPR1 | ITPR2 | ITPR3 | call | null tail | verdict |
|---|---|---|---|---|---|---|---|---|
| *Petromyzon marinus* | ITPR1.0 | 20 | 0 | 0 | 0 | no_call | 1.0000 | `no_call` |
| *Petromyzon marinus* | ITPR1.1 | 20 | 1 | 0 | 0 | ITPR1 | 0.0083 | `within_null` |
| *Petromyzon marinus* | ITPR1.2 | 20 | 1 | 0 | 0 | ITPR1 | 0.0083 | `within_null` |
| *Myxine glutinosa* | ITPR1.0 | 20 | 0 | 0 | 1 | ITPR3 | 0.0083 | `within_null` |
| *Myxine glutinosa* | ITPR1.1 | 20 | 0 | 0 | 0 | no_call | 1.0000 | `no_call` |
| *Myxine glutinosa* | ITPR1.2 | 20 | 2 | 1 | 0 | ITPR1 | 0.0055 | `within_null` |

**Synteny does not answer it.** All 6 loci carry at least 20 informative flank symbols, so the window is not the limitation — but the highest overlap any of them reaches with a gnathostome paralog consensus is **2**, and random neighbourhoods reach 2. 0 of 6 clear the null. The two species also disagree: the calls that do fire point at different paralogs.

Prior: **6 cyclostome loci in cyclostome-only clades** — S7 §5.4 — all 6 cyclostome loci sit in cyclostome-only clades, so they are neither one expansion nor three 1:1 ohnologs; S7 said explicitly that **which side of the vertebrate duplication each lineage attaches to** is what S8's synteny would settle. **underpowered** — roughly 550 My of independent rearrangement, and cyclostome annotations that name 27–29 % of their coding genes, leave no shared vocabulary to measure. S7 asked S8 to settle which side of the vertebrate duplication each cyclostome lineage attaches to; the answer is that flanking-gene synteny cannot, and the question needs an instrument that does not depend on orthologous gene *names* — Compara-style orthology calls on the flanks, or the 2R paralogon reconstructed from a cyclostome-anchored gene tree

## 6. Caveats

- **Flank orthology is by gene symbol.** RefSeq nomenclature is ortholog-derived, so this is not circular, but an unnamed gene cannot match anything: every Jaccard here is a **floor**, not a point estimate. That is exactly why nothing is read off a raw Jaccard — every claim is a comparison against the matched null.
- **270 loci have no gene table at all** (unannotated assemblies) and are in no statistic above. They are in `loci.tsv` with `has_gene_table=0`.
- **The root key over-merges.** `ZNF3` and `ZNF800` reach the same key. The paralogon result survives it because the background is built with the same rule, and because both surviving roots are carried by fewer than 1 % of random windows.
- **Leave-one-species-out is not leave-one-clade-out.** A mammalian locus is still scored against a consensus its close relatives voted into, so the caller's accuracy is an upper bound for a locus with no close relative in the set — which is exactly the cyclostome case. That is why §5.3 reports the null tail and not the accuracy.
- **A retained flanking ohnolog is a deletion record, not a phylogeny.** §4 says which ITPR neighbourhoods kept a shared family, not which paralogs are sisters; §4's verdict on S7's sister pair is stated in exactly those terms.
- Loci on scaffold-level assemblies can have fewer than 20 flanks (contig ends); `locus_sets.tsv` records the count per locus.

## 7. Figures

| figure | what it shows |
|---|---|
| `figures/synteny_pair_classes` | mean Jaccard per pair class beside its own matched random-window control |
| `figures/synteny_clade_decay` | the same, split into same-class and cross-class pairs |
| `figures/synteny_paralogon` | the three human neighbourhoods as gene tracks, with the shared ohnolog families linked, and their prevalence against background |
| `figures/synteny_caller` | the caller's two score distributions and the sweep that chose its operating point |

