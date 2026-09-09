## 4.10 How the search can be used as its own control series

Everything above is a search. What follows is a measurement of that search,
and it is placed here rather than in a methods appendix because the numbers it
produces are the error bars on Chapters 5, 9 and 13.

The design is what makes it possible, and it depends on a result that has not
been reported yet. Chapter 9 reconstructs no losses anywhere in this scope. If
no gene in the scope is absent, then every cell whose gene is independently
known to be present and which the ledger did not find is a **false negative of
the method**. It is not a biological absence and not an ambiguous case, but a
miss with a known right answer. That turns the whole sweep into its own
control series.

Two independent series are carried, because one of them would be circular on
its own. The first is every cell that Chapter 9's state rules call present.
The second is the ryanodine receptor sister cell, which uses none of those
rules at all: every vertebrate has three ryanodine receptors, so every
genome's control cell has a known answer that no part of this project's state
logic produced.

Before any of it is written, 32 constructed negative controls run and refuse
the build if they fail. Their character is worth noting, because it is the
character of every test suite in this project: they are checks on refusal and
on reachability. One requires that a false negative be **constructible** at
all, because a control that can only ever return zero misses is not a
measurement. Another requires the full-panel simulation to reproduce the
committed ledger cell for cell. A third requires the suite itself to alter no
committed table, and it exists because an earlier build's test called the real
routine and overwrote a committed table with its two constructed rows, after
which the report read two runs where there are seven.

## 4.11 The search misses about one cell in seven, and every miss is an assembly

Over the whole 309-genome scope the search misses **140 of 923 control cells,
15.2 %**.

The misses are not distributed evenly. A missed cell's assembly has a median
contig N50 of 23,460 bp against 3,396,515 bp for a found one, and the odds of
finding the gene rise 8.10-fold per tenfold of contig N50 in the IP₃ receptor
series and 19.99-fold in the ryanodine series. On chromosome-level assemblies
the IP₃ series misses 3 of 512 cells and the ryanodine series misses none of
172.

![](figures/s19_contiguity.png)

**{fig:s19_contiguity}.** The measured false-negative rate against assembly
contiguity, with the bar drawn. This is the figure that licenses the
retention result, because a survey reporting zero losses is worth nothing
unless somebody has measured how often the same search fails to find a gene
that is demonstrably present. The two series are independent, since the
ryanodine sister series uses none of the state rules the family series
depends on, and they agree, which is what stops the number being read as the
search agreeing with itself.

**The contiguity bar was chosen a priori and survives calibration.** The
142,212 bp floor came from gene geometry, meaning the median measured genomic
span of the gene, with no error rate anywhere in its derivation, and it had
never been scored against one. Scored now, the IP₃ series has a residual miss
rate of 0.9 % over 563 cells above the bar and the ryanodine series 0.0 % over
189. That is conservative against a conventional five-per-cent target, which
the scan reaches at a floor of 100,000 bp, and about right against a
one-per-cent target.

The bar is not free. It removes 38.8 % of the genomes, and it removes them
disproportionately from exactly the clades a loss survey is most interested
in. That cost is quantified rather than absorbed: the clade composition of
what each candidate floor retains is committed, so a reader can see which
question each floor makes unanswerable.

The residual is also diagnosed locally rather than dismissed. Contig N50 is a
genome-wide statistic and a regional assembly defect is invisible to it, so
each remaining miss is checked against the flanking-gene consensus of its own
paralogue to ask whether the region survived at all.

## 4.12 The bait panel was ablated exactly rather than modelled

How much of the panel was doing work? The question is normally answered by
argument. Here it is answered exactly, because the aligner aligns each bait
independently: dropping baits from a retained alignment and re-clustering
reproduces exactly what the sweep would have reported with a smaller panel.

The full chain is reproduced rather than approximated. It runs cluster, then
the identity floor, then the cell assignment where the family call lives,
because scoring the best coverage over every alignment instead would let an
ITPR3 bait's hit at the ITPR1 gene stand in as ITPR3 evidence and make every
panel look alike. The simulation reproduces the committed ledger 1,236 cells
out of 1,236, so every delta below is measured against the sweep itself rather
than against a model of it.

![](figures/s19_panel.png)

**{fig:s19_panel}.** Nineteen panels, drawn as change from the full panel on
a symmetric-log axis. Every panel scores between 0.58 and 0.85 of the cells,
so on an absolute linear axis the entire breadth result would be one pixel,
and the one ablation that removes every labelled bait sits 783 cells away
from the rest. The importance of this figure is that it overturns the usual
intuition about how to build a bait panel. Phylogenetic breadth, which is
what a panel is normally padded with, buys almost nothing inside the
vertebrates, while paralogue coverage buys everything. Anyone designing a
comparable survey should spend their panel budget on paralogues and on
clades where no labelled record exists, and this figure is the measurement
that says so.

Two readings fall out, and they point in opposite directions.

**Phylogenetic breadth buys almost nothing.** Four human baits, one per cell,
recover 782 of the 783 cells that the 38-bait panel recovers. Dropping any
single clade band costs at most two cells. The three unlabelled baits cost
nothing when removed and recover zero cells on their own, because the sweep
offers an unresolved-clade locus to whichever labelled paralogue scores
highest there, and with no labelled bait present there is nobody to offer it
to. Their function is to keep a real shark or lamprey gene inside the family
when it outranks every labelled bait, not to place it in a cell.

**Paralogue coverage buys everything.** Removing one paralogue's own baits
costs 238, 240 and 241 cells, which is about a quarter of the recovery each,
and no amount of breadth in the other two replaces it.

**Removing the sister-family control changes no call.** The positive family
test at genome scale costs nothing in recall, which is the cheapest possible
price for it.

**An ablation can gain a cell, and that is a property of the rule.**
Sixty-eight of the 2,179 call changes across all panels are gains, by the same
mechanism every time. A cell's coverage is read off its top-scoring bait
rather than its best-covering one, so removing a competitor can promote a bait
whose coverage is fractionally higher and push a marginal cell across the
coverage bar. The two gains under one ablation are a chicken bait at 0.6956
coverage giving way to a lizard bait at 0.7022, which is seven thousandths
either side of the threshold. They are reported rather than folded away,
because an ablation table with unexplained positive deltas invites the reading
that a smaller panel searches better.

## 4.13 What each search channel actually contributed, on three separate denominators

Three denominators are kept apart, because conflating them is the standard way
this question gets answered wrongly.

**How the census accumulated.** The targeted database search that started the
project returned 1,571 records. Exhaustive signature enumeration took it to
15,421. The vertebrate profile sweep added 618 and the genome sweep added
1,058. The non-vertebrate proteome sweep added 785 and the non-vertebrate
genome sweep 183, giving 18,065.

**Profile models against domain annotation, inside one database.** Both
channels are counted by what they call family rather than by raw hits, and
both are restricted to accessions the swept files actually hold. Comparing a
curated set against a raw hit list would credit the vertebrate sweep with the
13,371 targets its own gate declined.

At gene scale the two return nearly the same set. In the vertebrates the
enumeration returns 3,135 gene-scale records and the profile sweep 3,136, of
which 3,135 are shared, so the profile adds one. In the non-vertebrate metazoa
it adds two to 1,021. Its entire gain is in the fragment tail. The one
exception is worth naming: in the protists it adds 89 gene-scale records the
enumeration never returned, which is where the family's architecture is least
well annotated.

The other side of that statement is the same measurement seen from the other
direction: 766 vertebrate records carry a family signature and are declined by
the profile pair, all of them under 1,000 residues, which is the span gate
doing exactly what it was added for.

![](figures/s19_contribution.png)

**{fig:s19_contribution}.** What each channel contributed, drawn as
disagreement rates rather than as stacked counts. The claim is a proportion,
and stacking proportions on a logarithmic axis misreads them by
construction. The figure matters because it prices the standard approach. A
family profile searched over reference proteomes, which is how most
gene-family surveys are done, returns at gene scale almost exactly what
domain annotation already returns. Its entire gain is fragments, and the
only place it finds whole genes that annotation missed is the clades where
annotation is worst.

**Per gene, what a protein-database search would have missed.** Asked per
genome by cell rather than per record, with the genome sweep as ground truth
for where the genes are, **940 of 1,232 demonstrated genes, or 76.3 %, are not
reachable by any protein-database search.** They are not hard to find. No
protein record of them exists.

That is not an artefact of the margin species. Split by why each genome is in
scope, the rate is 74.3 % for order representatives against 78.5 % for margin
species. Chapter 13 is what that number turns into.

## 4.14 The rule written to catch iterative-search drift never fires

Chapter 3 left a question open: two of three iterative searches were killed by
a round ceiling rather than by the rule written for the hazard. With seven such
runs now available across the whole project, the rules can be scored as
classifiers.

The outcome is measured on the finished model rather than taken from any
session's prose, using the share of a run's final included set that neither
profile scores at all. The seven runs separate cleanly, with three at 0.81,
0.97 and 0.99 against four at 0.23 to 0.34 and nothing between 0.34 and 0.81,
so the labels are not a judgement call.

![](figures/s19_drift.png)

**{fig:s19_drift}.** Each kill rule scored as a classifier of a drift
outcome measured on the finished model. The figure records an instrument
failing at the one job it was built for: the rule written to catch
sister-family drift has a sensitivity of zero across seven runs, including
all three that drifted, because the drift dilutes the very quantity the rule
measures. The correction is a one-line change to a different axis, and it is
reported rather than applied, because it was validated after the fact on the
runs it would reclassify.

**The sister-family rule fires on none of the seven runs, including all three
that drifted.** The reason is now measured rather than described: off-family
accretion dilutes the sister-family share, so the rule's own statistic moves
the wrong way while a run drifts. The sister share fell over the run in four
of the seven.

The round ceiling catches all three drifted runs and also flags three that did
not. It counts rounds rather than content, which is right for a budget and
useless as a diagnosis.

**One change fixes it.** Moving the same threshold from the sister-family
share to the off-family share separates all seven runs completely, firing at
round 2 or 3 in each of the three that drifted and never in the other four, at
sensitivity 1.00 and specificity 1.00. It is reported and **not applied**,
because it is validated after the fact on seven runs with an inherited
threshold, and applying it would overturn committed verdicts. That is a
decision for a session that owns those tables rather than for the one that
noticed.

**Seed choice does not change the family, only the debris.** The three
vertebrate runs were seeded from a human paralogue, a fly gene and an
amoebozoan gene. Their family content is the same set three times over: they
intersect on 4,785 family records and differ by at most 162. What the seed
changes is everything that is not the family, in that the amoebozoan run
carries 19,949 non-family targets against roughly 2,439 for each of the other
two. A distant seed does not reach further into the family. It reaches further
out of it.

## 4.15 The inference ran out in four places, rather than the search

Four limits were met while doing something else, and each is a general result
about the method rather than about this family.

**A reconciliation's loss count is mostly sampling.** Losses implied by a gene
tree over a 134-tip alignment, checked cell by cell against the 309-genome
ledger, give 47 implied cells, 26 of them the paralogue present in the genome
and absent only from the sample, and none corroborated.

**Synteny disambiguation is accurate and unavailable.** The neighbourhood
caller built in Chapter 7 is right wherever it acts and almost never acts on
the cells that need it: of 432 trace regions it reaches 8, because a
fragment's contig carries no neighbours to read. Accuracy and reach are
different numbers, and quoting only the first would report contig lengths as
method performance.

**Per-site codon models are past their power at these divergences.** The
synonymous rate is unidentifiable at 671 of 7,368 sites, and 373 of 420
pairwise comparisons are flagged saturated. Any per-site rate formed as a
ratio of sums is dominated by sites whose denominator is not estimable.

**A likelihood tree need not resolve the question asked of it.** The test
that compares the three sister arrangements leaves two topologies inside its
95 % confidence set, so the sister question is answered by what the data
excludes rather than by reading the best tree. Chapter 6 reports both.

## 4.16 Five results here apply to anyone doing the same thing

1. **A genome-scale orthologue sweep misses about one cell in seven, and every
   miss is an assembly.** A survey that reports absences without a contiguity
   floor is reporting assembly quality.
2. **A contiguity bar can be set from gene geometry before any error is
   measured**, and it lands where a calibration would have put it, at a cost
   of 38.8 % of the genomes, disproportionately in the clades a loss survey
   most cares about.
3. **Phylogenetic breadth in a protein bait panel is nearly worthless within
   the vertebrates, and paralogue coverage is everything.** Spend the budget
   on paralogues and on clades where no labelled record exists.
4. **A family profile over reference proteomes is not a discovery method for
   whole genes.** At gene scale it returns what domain annotation already
   returns, its gain is fragments, and it adds gene-scale records only where
   the annotation is worst.
5. **The obvious rule for catching iterative-search drift measures the wrong
   thing.** Sister-family contamination is the quantity this project and its
   predecessor both chose to monitor, and it is diluted by the very drift it
   is meant to detect.
