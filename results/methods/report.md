# S19 — methods results: what the search was worth

*Generated 2026-09-08 by `scripts/s19_run.py` from the
committed tables in `results/methods/`. Nothing here is hand-written.*

## 1. What this task can measure that most cannot

Every methods section that reports a search's sensitivity has to estimate it,
because the genes the search missed are the ones nobody can count. This
project does not have that problem. S15b reconstructed **no losses anywhere
in the 309-genome scope**: all
923 assignable genome × paralog cells hold a
gene that is there, established by evidence the sweep's locus placement did
not produce. So every ledger cell that is not `found` is a **false negative
of the method**, and the sensitivity of a genome-scale ortholog sweep is a
direct measurement.

There are two such series and they were measured by the same instrument in
the same assemblies:

| series | cells | missed | false-negative rate | 95 % CI |
|---|---|---|---|---|
| ITPR cells S15 states present | 923 | 140 | 15.2 % | 13.0 % – 17.6 % |
| RyR sister cell, present in every vertebrate | 309 | 42 | 13.6 % | 10.2 % – 17.9 % |

The two rates are not distinguishable (Fisher exact
*p* = 0.578), which is the point of carrying the second
one: the ITPR figure rests on S15a's state assignments and the RyR figure
does not, so agreement between them is evidence that the number is a
property of the search rather than of how the states were called.

Exactly 4 cells are in **neither** series,
and they are the cyclostome ITPR2/ITPR3 cells S15a filed
`paralog_unassignable`: those genomes carry spare ITPR loci the bait panel
cannot label, so calling them present or absent would decide a question S15a
declined.

## 2. Scope and inputs

Nothing in this task downloads or re-aligns. Every input is a committed
artefact of an earlier session, with one derived asset built here:

- the **swept accession universe** — one `grep '^>'` pass over each of the
  seven reference-proteome FASTAs, cached under the data root. Without it,
  "the InterPro enumeration holds a record the profile HMM did not return"
  cannot be told apart from "that record was never in the database the
  profile HMM searched", and the head-to-head in §5 would be an accounting
  artefact rather than a comparison.
- the **retained per-genome miniprot alignments** — 309 GFFs. miniprot
  aligns each bait independently, so dropping baits from the file and
  re-clustering reproduces exactly what the sweep would have reported had
  those baits never been in the panel. The ablation in §6 is therefore
  exact, not a model of one.

Total: 1,232 control cells, 309 genomes,
7 jackhmmer runs, and
15.2 recorded compute-hours over eight search
channels.

## 3. Negative controls

`s19_test_methods.py` runs before anything is written and refuses the build
if it fails. The checks are on *refusal* and on *reachability*, because every
rule in this task returns a plausible number when it is wrong.

| check | what it refuses |
|---|---|
| T1–T2 | a FASTA header must yield a bare accession in all three shapes, and a genome model id must survive intact — stripping one at its first `.` would collapse every model of an assembly onto its accession prefix |
| T3 | a zero-row parse must refuse to cache a universe; an empty universe makes every record read as absent from the database it was found in |
| T4 | a `paralog_unassignable` cell must enter neither control series |
| T5 | a false negative must be *constructible* — a control that can only ever return zero misses is not a measurement |
| T6 | raising the contiguity floor must never retain more genomes, and D4's own bar must be on the scanned grid |
| T7 | the full-panel simulation must reproduce the committed ledger cell for cell |
| T8 | a locus with no labelled paralog bait must fill no paralog cell, and must be offered once one aligns there |
| T9 | D14 — a locus won by the RyR baits must be offered to no ITPR cell |
| T10 | a killed jackhmmer run must contribute only its accepted rounds |
| T11 | the proposed drift rule must fire on a drifting run and decline on a stable one, both constructed |
| T12 | the drift label must come from the finished model, and the two populations must separate |
| T13 | the head-to-head must compare family *calls*, not raw profile targets |
| T14 | every manifest genome must land in exactly one scope class |
| T15 | the suite itself must alter no committed table — the first build's T11 called the real trace routine and overwrote `kill_criterion_trace.tsv` with its two constructed rows, and the report then read 2 jackhmmer runs where there are 7 |

All 32 checks pass. Five were
mutation-tested by breaking the rule they guard — counting
`paralog_unassignable` as present, comparing raw targets instead of calls,
pointing the proposed drift rule back at K1's axis, letting an empty universe
cache, and letting the suite write to the results directory — and every
mutation was caught by its own check.


## 4. Contiguity is the whole of the false-negative rate

The misses are not distributed across the scope. A missed control cell's
assembly has a median contig N50 of
23,460 bp
against
3,396,515 bp
for a found one, and the odds of finding the gene rise
8.10×
per tenfold of contig N50 in the ITPR series and
19.99×
in the RyR series. On chromosome-level assemblies the ITPR series misses
3/512
cells and the RyR series misses
0/172.

### 4.1 D4's bar was chosen a priori, and it survives calibration

D4's contiguity bar is **142,212 bp** of contig N50 — the median measured
ITPR genomic span, taken from gene geometry alone, with no error rate in its
derivation. This is the first time it has been scored against one:

| series | genomes kept | cells | missed | residual rate | upper 95 % bound |
|---|---|---|---|---|---|
| ITPR cells | 189 | 563 | 5 | 0.9 % | 2.1 % |
| RyR sister cell | 189 | 189 | 0 | 0.0 % | 2.0 % |

The bar the project has been using since S5a lands at under one per cent
residual error in the series that can have any, and at zero in the other. It
is *conservative* against the conventional five-per-cent target, which the
scan reaches at a floor of
100,000 bp,
and about right against a one-per-cent target.

Prior: **D4's bar is the median measured ITPR genomic span, 142,212 bp of contig N50, chosen a priori from gene geometry and never calibrated against a measured error rate** — roadmap D4; scripts/s5_calibration.py:itpr_span_stats. **confirmed** — scored against a measured false-negative rate for the first time, the a-priori bar gives 0.9 % residual error on the ITPR series and 0.0 % on the RyR sister, retaining 189 of the 309 genomes

Prior: **S5b measured recovery at 98-99 % above D4's contiguity bar and 57-70 % below it, over a scope in which S15b later reconstructed no losses at all** — roadmap S5b Results; results/loss_counts/report.md. **confirmed** — S5b's 98-99 % / 57-70 % split is reproduced as a false-negative rate: 15.2 % over the whole scope, 0.9 % above the bar

### 4.2 What the floor costs, and which clades it removes

A floor that reaches one per cent by keeping a third of the scope has not
made the survey more reliable, it has made it smaller. At D4's own bar
189 of 309
genomes are retained; `floor_clade_composition.tsv` records which classes
each floor removes, and the answer is the one S5a already warned about — the
margin species S4 added for their *proteome* gaps are very largely the same
genomes whose assemblies cannot hold the gene, so a contiguity filter removes
the part of the scope the scope was extended for.

### 4.3 Below the bar, the misses run with gene span

| cell | cells below the bar | false-negative rate | 95 % CI |
|---|---|---|---|
| ITPR1 | 120 | 39.2 % | 30.9 % – 48.1 % |
| ITPR2 | 120 | 43.3 % | 34.8 % – 52.3 % |
| ITPR3 | 120 | 30.0 % | 22.5 % – 38.7 % |

Prior: **S5b confirmed a 13-point recovery gap below the bar in the direction of gene span — ITPR3 70 %, ITPR1 61 %, ITPR2 57 %** — roadmap S5b Results. **confirmed** — below the bar the ordering is the same one S5b measured and points the same way as gene span: the shortest gene is the one the fragmented assembly still yields

### 4.4 The residual, diagnosed locally

Contig N50 is a genome-wide statistic and cannot see a regional assembly
defect. Asking S8's own consensus flanks whether the neighbourhood survived
at all: **57 of 182** false negatives sit in an
assembly that carries fewer than half of that paralog's usual neighbours, so
what is missing there is the region and the cell says nothing about the gene.


## 5. What each search channel contributed

### 5.1 How the census accumulated

| census | channel that produced it | records | added |
|---|---|---|---|
| v1 | targeted database search (the app: NCBI/Ensembl/UniProt/Compara) | 1,571 | 1,571 |
| v2 | InterPro/Pfam exhaustive enumeration | 15,421 | 15,421 |
| v3 | profile HMM + jackhmmer, vertebrate reference proteomes | 16,039 | 618 |
| v4 | genome sweep (miniprot, 309 vertebrate assemblies) | 17,097 | 1,058 |
| v5 | profile HMM + jackhmmer, 6,928 other reference proteomes | 17,882 | 785 |
| v6 | genome sweep (miniprot, 194 non-vertebrate assemblies) | 18,065 | 183 |

### 5.2 Profile HMM against domain annotation, inside one database

Both channels are counted by what they **call family** — the enumeration by
carrying a family signature, the sweep by clearing D22's gate — and both are
restricted to accessions the swept FASTAs actually hold. Comparing a curated
set against a raw hit list would credit the vertebrate sweep with the
13,371 targets its own gate declined.

At **gene scale** they return nearly the same set:

| database | Pfam enumeration | profile HMM | shared | HMM only | Pfam only |
|---|---|---|---|---|---|
| vertebrata | 3,135 | 3,136 | 3,135 | 1 | 0 |
| fungi | 26 | 26 | 26 | 0 | 0 |
| viridiplantae | 25 | 21 | 21 | 0 | 4 |
| metazoa_nonvert | 1,021 | 1,023 | 1,021 | 2 | 0 |
| protista_other | 488 | 577 | 488 | 89 | 0 |

In the **fragment tail** they do not:

| database | Pfam enumeration | profile HMM | HMM only | Pfam only |
|---|---|---|---|---|
| vertebrata | 1,681 | 1,509 | 594 | 766 |
| fungi | 15 | 17 | 16 | 14 |
| viridiplantae | 22 | 19 | 11 | 14 |
| metazoa_nonvert | 652 | 916 | 469 | 205 |
| protista_other | 82 | 152 | 102 | 32 |

So a family profile HMM over reference proteomes is not a discovery method
for whole genes — in the vertebrates it adds one gene-scale record to
3,135, and in the non-vertebrate metazoa two to
1,021. Its entire gain is in fragments. The
one exception is worth naming: in the **protists** it adds
89 gene-scale records the enumeration
never returned, which is where the family's architecture is least well
annotated.

The Pfam-only column is the same statement from the other side:
766 vertebrate records carry a family
signature and are declined by the profile pair, all of them under 1,000 aa
— D22's 200-match-state gate doing what it was added for.

Prior: **S3 found the profile sweep adds 618 proteins the InterPro census never returned, all 209-942 aa fragmentary gene models, and that 0 of 2,787 v2 ITPR records in a swept proteome were missed** — results/census_v3/report.md. **confirmed** — S3's fragment characterisation is now a measurement rather than a description: of everything the vertebrate profile pair returns that the enumeration did not, exactly 1 record is gene-scale and 594 are under 1,000 aa

### 5.3 Per gene: what a protein-database search would have missed

The question asked per genome × cell rather than per record, because the
genome sweep is the ground truth for where the genes are:

| how the gene is reachable | cells |
|---|---|
| genome only — no reference proteome | 386 |
| protein database and genome | 292 |
| genome only — no record resolves to this paralog | 286 |
| genome only — records exist but none full length | 248 |
| genome only — no family record for species | 20 |

**940 of
1,232** demonstrated genes
(76.3 %) are
not reachable by any protein-database search. That is not an artefact of the
margin species: split by why S4 put each genome in scope, the rate is
74.3 % for order
representatives against
78.5 % for margin species.

Prior: **S5b's census v4 adds 1,058 gene models from 224 genomes, of which 318 ITPR models exist only as DNA and 167 more sit inside an annotated gene carrying no family name** — roadmap S5b Results; results/census_v4/genome_models.tsv. **confirmed** — S5b's census v4 counts are the per-locus form of this; per cell the figure is 76.3 %, and the sister family is the best-served cell at 63.4 %

### 5.4 Cost

15.2 recorded compute-hours over eight channels.
4 channels carry no wall clock at all and
are reported blank rather than as zero: `method_cost.tsv` names them and
says why.


## 6. The bait panel: breadth is nearly free, paralog coverage is not

### 6.1 The ablation is exact, and validated first

The simulation reproduces the committed ledger **1,236
cells out of 1,236** under the full panel
(0 disagreements). Every delta below is
therefore measured against the sweep itself, not against a model of it.

### 6.2 What each ablation costs

| panel | baits | ITPR cells found | change | description |
|---|---|---|---|---|
| full38 | 38 | 783 | +0 | the panel as run (S5b) |
| labelled_only | 35 | 783 | +0 | no unlabelled `vertebrate_basal` bait — the gar, chimaera and lamprey seeds that carry no paralog assignment |
| no_ryr_control | 30 | 783 | +0 | the ITPR baits alone, no sister control |
| drop_bird | 31 | 785 | +2 | no bird bait |
| drop_mammal | 31 | 783 | +0 | no mammal bait |
| drop_ray_finned_fish | 30 | 782 | -1 | no ray-finned fish bait |
| amniote_only | 18 | 781 | -2 | mammal + bird + reptile baits only — a tetrapod-centric panel |
| s3_seeds_only | 25 | 782 | -1 | the S3 seed set alone, before S5a's additions |
| human_only | 4 | 782 | -1 | the human baits alone |
| one_per_cell_nonhuman | 4 | 782 | -1 | one non-human bait per cell — the longest of each |
| drop_ITPR1_baits | 29 | 545 | -238 | no bait labelled ITPR1 — can the other paralogs' baits and the unlabelled ones still place its gene? |
| drop_ITPR2_baits | 30 | 543 | -240 | no bait labelled ITPR2 — can the other paralogs' baits and the unlabelled ones still place its gene? |
| drop_ITPR3_baits | 28 | 542 | -241 | no bait labelled ITPR3 — can the other paralogs' baits and the unlabelled ones still place its gene? |
| basal_only | 3 | 0 | -783 | only the three unlabelled baits |

Two readings, and they point opposite ways.

**Phylogenetic breadth buys almost nothing.** Four human baits — one per
cell — recover 782 of the
783 cells the 38-bait panel recovers. Dropping any
single clade band costs at most two cells. The three unlabelled
`vertebrate_basal` baits cost **nothing** when removed, and on their own
recover **zero** cells: an unlabelled bait cannot fill a paralog cell by
itself, because the sweep offers an unresolved-clade locus to whichever
*labelled* paralog scores highest there, and with no labelled bait present
there is nobody to offer it to. Their function is to keep a real shark or
lamprey gene inside the family when it outranks every labelled bait, not to
place it in a cell.

**Paralog coverage buys everything.** Removing one paralog's own baits costs
241, 240
and 238 cells — about a quarter of the
recovery each — and no amount of breadth in the other paralogs replaces it.

**Removing the sister-family control changes no ITPR call.** D14 at genome
scale costs nothing in recall, which is the cheapest possible price for a
positive family test.

### 6.3 An ablation can *gain* a cell, and that is a property of the rule

68 of the 2,179 call changes across all panels
are gains. The mechanism is the same every time: a cell's coverage is read
off its **top-scoring** bait, not its best-covering one, so removing a
competitor can promote a bait whose coverage is fractionally higher and push
a marginal cell across the 0.7 coverage bar.
The two gains under `drop_bird` are both a chicken ITPR1 bait at coverage
0.6956 giving way to a lizard bait at 0.7022 — seven thousandths either side
of the threshold. They are reported here rather than folded away, because an
ablation table with unexplained positive deltas invites the reading that a
smaller panel searches better.

### 6.4 The six unfilled slots

| cell | band | classes | cells | found by the full panel | found without the unlabelled baits |
|---|---|---|---|---|---|
| ITPR1 | chondrichthyes | Chondrichthyes | 12 | 12 | 12 |
| ITPR1 | cyclostomata | Hyperoartia;Myxini | 2 | 2 | 2 |
| ITPR2 | chondrichthyes | Chondrichthyes | 12 | 12 | 12 |
| ITPR2 | cyclostomata | Hyperoartia;Myxini | 2 | 0 | 0 |
| ITPR2 | sarcopterygian_fish | Coelacanthimorpha;Dipnoi | 2 | 2 | 2 |
| ITPR3 | cyclostomata | Hyperoartia;Myxini | 2 | 0 | 0 |

S5's panel has no labelled bait for these six slots because no labelled
full-length record exists there. The chondrichthyan and coelacanth cells are
recovered anyway, by labelled baits from other clades — those genes are
within reach of a bony-fish or tetrapod bait. The four cyclostome cells are
not recovered by any panel, which is S5b's own result and the reason S15a
filed them `paralog_unassignable`.

Prior: **S5a found the unlabelled `vertebrate_basal` baits were competing as a fourth paralog until rescue attribution was fixed, and S5's panel leaves six slots unfilled for want of a labelled record** — roadmap S5a Results; results/s5_baits/unfilled_slots.tsv. **not corroborated** — the unlabelled baits change no miniprot cell call in either direction (0 cells), so S5a's fix mattered to rescue attribution — a different code path — and not to locus placement; and the six unfilled slots cost only the four cyclostome cells

### 6.5 The design rule

Every (genome, cell, bait) triple the sweep produced is a measurement of what
one bait alone recovers at a known sequence distance:

| bait-to-target identity | measurements | recovered | recall | lower 95 % bound |
|---|---|---|---|---|
| no alignment | 38 | 0 | 0.0 % | 0.0 % |
| 0.50-0.60 | 14 | 14 | 100.0 % | 78.5 % |
| 0.60-0.70 | 1,145 | 1,143 | 99.8 % | 99.4 % |
| 0.70-0.80 | 2,396 | 2,396 | 100.0 % | 99.8 % |
| 0.80-0.90 | 3,047 | 3,047 | 100.0 % | 99.9 % |
| 0.90-0.95 | 1,303 | 1,303 | 100.0 % | 99.7 % |
| 0.95-1.01 | 1,256 | 1,256 | 100.0 % | 99.7 % |

A single ITPR bait recovers the gene essentially always at any identity above
0.5. The scope of that claim matters: it is measured **at loci the full panel
already placed**, so it says what one bait can do where a gene is known to
be, not what a one-bait panel would find in an unsearched genome. Read that
way it is still the useful number — it is why a four-bait panel matches a
38-bait one here.


## 7. The rule written to catch iterative drift never fires

Seven jackhmmer runs, 3 of which walked out of the
family. The outcome is measured on the **finished model** rather than taken
from any session's prose: the share of a run's final included set that
neither profile scores at all. The seven separate cleanly —
s20_fungi 0.99, s20_protista_other 0.97, itpr_acanthamoeba 0.81
against
s20_metazoa_nonvert 0.25, itpr1_human 0.34, itpr_fly 0.34, s20_viridiplantae 0.23
— with nothing between 0.34 and
0.81, so the labels are not a judgement call.

Scoring each of D10's rules as a classifier of that outcome:

| rule | fires on a drifted run | fires on one that did not | sensitivity | specificity |
|---|---|---|---|---|
| K1 (sister-family share rises > 0.10 in a round) | 0/3 | 0/4 | 0.0 % | 100.0 % |
| K2 (a round includes > 10x the previous) | 1/3 | 0/4 | 33.3 % | 100.0 % |
| K3 (the 10-round ceiling) | 3/3 | 3/4 | 100.0 % | 25.0 % |
| proposed: off-family share rises > 0.10 in a round | 3/3 | 0/4 | 100.0 % | 100.0 % |

**K1 — the rule written for exactly this hazard — fires on none of the seven
runs, including all three that drifted.** The reason is D10b, and it is now
measured rather than described: off-family accretion *dilutes* the
sister-family share, so K1's own statistic moves the wrong way while a run
drifts. The sister share **fell** over the run in
4 of 7 runs.

K3, the round ceiling, catches all three but also flags three runs that did
not drift — it counts rounds rather than content, which is right for a budget
and useless as a diagnosis. That is why S20b had to separate `protista_other`
from `metazoa_nonvert` in prose: both hit K3, one had accreted 22,913
off-family targets and the other had simply not finished.

### 7.1 One change fixes it, and it is offered as a proposal

Moving K1's own threshold — the same 0.10 — from the sister-family share to
the **off-family** share separates the seven runs completely, firing at round
2 or 3 in each of s20_fungi, s20_protista_other, itpr_acanthamoeba and never in the other four. It is
reported here and **not applied**: it is validated after the fact on seven
runs with an inherited threshold, and applying it would overturn committed
verdicts, which is a decision for a session that owns those tables.

Prior: **S3 and S20b found off-family accretion *dilutes* the sister-family share, so K1 moves the wrong way while a run drifts and only the round ceiling catches it** — roadmap D10b; results/s20_sweep/report.md. **confirmed** — K1 fires on 0 of 7 runs and the sister share falls in 4 of them; the off-family axis separates the same seven runs at sensitivity 1.00 and specificity 1.00

### 7.2 What iteration bought

| run | database | one-pass family calls | accepted targets | final targets | verdict | rule |
|---|---|---|---|---|---|---|
| itpr1_human | vertebrata | 5,130 | 7,226 | 7,222 | killed | K3 |
| itpr_fly | vertebrata | 5,130 | 7,235 | 7,228 | clean | - |
| itpr_acanthamoeba | vertebrata | 5,130 | 24,902 | 26,266 | killed | K3 |
| s20_viridiplantae | viridiplantae | 55 | 70 | 70 | clean | - |
| s20_fungi | fungi | 55 | 42 | 6,302 | killed | K2 |
| s20_protista_other | protista_other | 850 | 22,913 | 26,148 | killed | K3 |
| s20_metazoa_nonvert | metazoa_nonvert | 2,350 | 2,908 | 2,921 | killed | K3 |

What those targets *are* is the question a completeness argument turns on:

| database | one-pass family calls | accepted targets across runs | of which called family | of which not |
|---|---|---|---|---|
| vertebrata | 5,130 | 24,929 | 4,960 | 19,969 |
| fungi | 55 | 42 | 33 | 9 |
| viridiplantae | 55 | 70 | 54 | 16 |
| metazoa_nonvert | 2,350 | 2,908 | 2,198 | 710 |
| protista_other | 850 | 22,913 | 772 | 22,141 |

In the vertebrate database, three iterated runs between them hold
4,960 of the
5,130 records one hmmsearch pass calls
family, and 19,969 that it does not.
So iteration returned **no record the profile pair calls family that one
pass had not already returned**, and twenty thousand that it does not call
family at all. The same shape holds in the protists
(772 family against
22,141 not). Only the run that
converged cleanly stays close to its one-pass set. Whether any of those
non-family targets is a real family member the profiles failed to score is
not answerable from this table, and S3 and S20 both examined them: they are
module-only matches to the shared domains.

### 7.3 Seed choice does not change the family, only the debris

The three vertebrate runs were seeded from a human paralog, a fly gene and an
amoebozoan gene — as unlike each other as this family allows.

| run | seed | accepted targets | of which family | family unique to this seed | family calls this run lacks |
|---|---|---|---|---|---|
| itpr1_human | Q14643 | 7,226 | 4,787 | 0 | 343 |
| itpr_fly | P29993 | 7,235 | 4,797 | 6 | 333 |
| itpr_acanthamoeba | L8GF85 | 24,902 | 4,953 | 162 | 177 |
| ALL THREE (intersection) |  | — | 4,785 | — | 170 |

The family content is the same set three times over: the three runs
intersect on 4,785 family records and
differ by at most 162.
What the seed changes is everything that is *not* the family: the amoebozoan
run carries 19,949 non-family
targets against
2,439 and
2,438 for the other two. A distant
seed does not reach further into the family; it reaches further out of it.


## 8. Where the inference methods ran out, not the search

Four limits this project met while doing something else. Each is a general
result about the method.

**A reconciliation's loss count is mostly sampling.** Losses implied by a
gene tree over a 134-tip representative alignment, checked cell by cell
against the 309-genome ledger: 47 implied
cells, 26 of them the paralog present in
the genome and absent only from the sample, and
0 corroborated.

Prior: **S13 audited its own reconciliation against the genome ledger and found the implied losses are overwhelmingly sampling, with four corroborated absences before D45 removed them** — results/reconciliation/report.md. **confirmed** — 26 of 47 implied-loss cells are sampling; the recomputation agrees with S13's own aggregate row

**Synteny disambiguation is accurate and unavailable.** S8's caller is right
wherever it acts and almost never acts on the cells that need it: of
432 trace regions it reaches
8
(1.8 %), because a fragment's contig carries no
neighbours to read. Accuracy and reach are different numbers and quoting only
the first would report contig lengths as method performance.

Prior: **S15a found S8's caller is accurate where it acts and almost never reaches — 273 of 432 trace regions sit on a contig carrying no annotated gene at all, and 8 reach the four-key floor** — results/loss_dynamics/report.md. **confirmed** — 8 of 432 regions reach the caller's floor, at 100 % accuracy in every key bin where it acts

**Per-site codon models are past their power at these divergences.** The
synonymous rate is unidentifiable at
671 of 7,368 sites
(9.1 %), and
373 of
420 pairwise comparisons
(88.8 %) are flagged saturated. Any per-site or
per-element omega formed as a ratio of sums is dominated by sites whose
denominator is not estimable.

Prior: **S17 reported per-site selection on the same coordinates as its constraint layers, taking omega only over sites where the synonymous rate is identifiable** — results/constraint/report.md. **orthogonal** — S17 already took omega only over identifiable sites; this recomputation states the size of what that excludes (9.1 % of scored sites) rather than disagreeing with it

**A likelihood tree need not resolve the question asked of it.** The AU test
leaves 2 topologies in the 95 % confidence
set, so the sister arrangement is reported as what the data supports rather
than as the ML tree's answer.


## 9. What this changes for anyone doing the same thing

Five results, in the order they would change a protocol.

1. **A genome-scale ortholog sweep misses about one cell in seven, and every
   miss is an assembly.** 15.2 % of control cells
   over the whole scope, 0.9 % above a
   contiguity bar set at the gene's own median span. A survey that reports
   absences without a contiguity floor is reporting assembly quality.
2. **A contiguity bar can be set from gene geometry before any error is
   measured.** D4's was, and it lands where the calibration would have put
   it. The bar is not free: it costs
   38.8 % of the
   genomes, disproportionately in the clades a loss survey is most interested
   in.
3. **Phylogenetic breadth in a protein bait panel is nearly worthless within
   the vertebrates; paralog coverage is everything.** Four human baits match
   thirty-eight. One bait recovers its gene at any identity above 0.5. Spend
   the panel budget on paralogs and on clades where no labelled record
   exists, not on sampling depth within a paralog.
4. **A family profile HMM over reference proteomes is not a discovery method
   for whole genes.** At gene scale it returns what domain annotation already
   returns. Its gain is fragments — and, in the least-annotated clades,
   gene-scale records the annotation never carried.
5. **The rule everyone writes to catch iterative-search drift measures the
   wrong thing.** Sister-family contamination is diluted by the drift it is
   meant to detect. The off-family share is what moves.

And the one that is about the archive rather than the method:
**76.3 % of
the genes this project demonstrated are not reachable from any protein
database.** Not because they are hard to find — because no protein record of
them exists.


## 10. What this does not settle

- **The control rests on S15b's zero.** If a loss were later established
  anywhere in the scope, that cell would move from the numerator of the
  false-negative rate to a real absence. The RyR series is carried precisely
  so the reader can see how much of the number depends on that: it is
  measured independently of S15a's states and agrees to within two points.
- **`present_fragmented` and `present_partial` states were reached using
  alignment evidence.** They are outside the loci the sweep placed and
  cleared a calibrated bar, so they are not the sweep's own answer restated
  — but they are not an independent instrument either. The stricter reading
  is the RyR series, which does not use them at all.
- **The panel ablation is exact for locus placement and silent about
  rescue.** Removing baits changes what tblastn rescue attributes, and that
  path needs the genome FASTA, which the sweep deleted after searching. So
  the ablation measures what the panel was worth to miniprot, not to the
  whole pipeline.
- **The single-bait design rule is conditional on the gene being there.** It
  is measured at loci the full panel placed. It bounds what one bait can do
  at a known locus; it does not license a one-bait panel in an unsearched
  clade.
- **The drift criterion is validated after the fact on seven runs**, with a
  threshold inherited from K1 rather than fitted. Seven runs cannot establish
  a specificity of 1.00. It is a proposal with its evidence attached, and it
  is deliberately not applied to any committed verdict.
- **The head-to-head is a comparison of two channels' *calls*, not of two
  algorithms in the abstract.** The profile pair's gate (D22) is part of the
  method being scored, and a different gate would move the tail columns.

