# 14. Methods, and the reasoning behind them

## 14.1 What a methods chapter is for here

A paper's methods section says what was done. It has to, because a reader
needs to reproduce it. It cannot say *why* each threshold is the number it is,
and for this project that is where most of the work went.

Seventy-five methodological decisions were recorded across twenty-six working
sessions, numbered as they were made. Several of them changed an answer.
Several were made, tested and overturned. None of them has ever been written
out as prose, and this chapter is that: not a list, but the arguments in the
order they matter, with the incident that produced each one.

The chapter is organised by what a decision is *about* — thresholds, controls,
refusals, reproducibility, and the rules that keep a document from drifting
away from its data.

## 14.2 Thresholds: measured, not chosen, and not ported

**A threshold ported from another project must be re-measured before it is
used, not merely restated in this project's units.** This project inherited a
methodology from a sister project on a different gene family, and the
inheritance is deliberate: the methods were paid for over two dozen sessions
and re-deriving them buys nothing.

The boundary of that rule is a number. The attribution margin arrived as a
bit-score ratio restated in this project's units, which made it *look* derived
while leaving it tuned to a family whose paralogues are 40 to 50 % identical.
These are 61 to 68 % identical. Measured against loci whose identity an
independent annotation establishes, the inherited value would have reported
every rescue fragment ambiguous and made every absence claim unattributable.

**A method ports. A number that encodes how far apart that family's paralogues
sit does not** — and neither does one that encodes its gene sizes. The same
session found the aligner's maximum-intron parameter needed measuring for the
same reason. Where a ported constant cannot be re-measured yet, it is used
with its limitation stated in the code and the per-row value recorded, so the
cut can be revisited without re-running anything.

**A threshold cannot be calibrated from the population it has already
filtered.** The sweeps deliberately record every candidate down to a floor far
below the one they apply, precisely so that the floor is not measured on data
it has already removed.

**A threshold the data has no population for is not lowered to a round
number.** One brief asked for a merge rule that would join the pair of coding
blocks an aligner emits either side of a frameshift. Measured across all 309
genomes, **no such pair exists**: every one of 149,148 consecutive block pairs
has contiguous query spans and the smallest genomic gap anywhere is 10 bp. So
the calibration refuses to derive a bar — there is nothing on the other side
of it — and the floor goes at the smallest gap the genome calls spliceable
rather than at the declared fallback, which would have merged ten real
junctions. The merge then fires on zero junctions, and a constructed control
builds such a pair and requires the branch to fire, which is what makes the
zero a measurement rather than a rule that cannot act.

**And a calibration must be able to refuse.** More than one in this project
does: it reports its separation statistic beside its threshold and raises
rather than handing out a bar its own data does not support. One refused with
a Youden index of 0.52 and sent a chapter back to find the population that was
contaminating its decoy.

**A calibration must not be able to pass vacuously.** A smoke-test run over
one genome once produced a threshold from two loci, wrote it, and the next
sweep read it back and moved three genomes between statuses. Calibrations now
carry floors on the number of observations they may be computed from, and one
equivalence test refuses to report success on a comparison of fewer than 25
targets after its first run went green on a group that had no hits at either
setting.

**An operating point is not the same as a separation statistic.** Youden's
index on a perfectly separated pair of populations lands on the lowest
positive, which is the most permissive bar the data allow. Where a gap exists,
the operating point goes at its midpoint and **both edges are committed**, so
a later run whose populations have drifted into contact is visible in the
table rather than only in a verdict.

## 14.3 Controls: positive, negative, and matched

**The sister family is this project's positive control, not its nuisance.**
The ryanodine receptors carry every diagnostic domain of the IP₃ receptors,
are twice their length, and are inside every search this project ran. Rather
than filter them out, they are carried through every stage and measured
through the identical instrument: as a presence control in the genome sweep,
as a decoy in the discovery benchmark, as the outgroup that roots the tree, as
a second family the synteny method has to work on, as the family whose
triplication makes the teleost singletons interpretable, as the annotation
audit's comparison, and as the ligand-site chapter's positive control for a
receptor with no IP₃ site.

**A positive control does not port across the scope it controls.** The
vertebrate sweep can write "the ryanodine control fired in all 309 genomes"
because every vertebrate has three. Outside the metazoa that sentence is not
available. The replacement was a domain-sharing protein family — and *that*
failed too, in the one clade where it matters, because a single profile fixed
in advance for every clade is the mistake rather than the choice of profile.
The fix is to choose the control per clade by measurement, running six
candidates over each clade's own proteomes and taking the one that is
actually there.

**A within-genome paired test makes a cell and its siblings one observation.**
Where a comparison is between paralogues inside one genome, the strata from it
are not independent corroboration: an elevated cell makes its own siblings
look deficient by construction, and three significant rows can be one result.

**A confounder that can be matched is matched, not argued away.** Two chapters
here found their headline effect disappear under identity matching, and both
report the matched result as the answer. A third measured a risk that a sister
project had dismissed in a caveat — cross-mapping between paralogues — by
tiling every reference exhaustively with synthetic reads and mapping them
back. **An argument in a caveat and a measurement in a table are not the same
evidence, and the second costs minutes.**

**A ported control's *result* is re-measured, never inherited.** The sister
project's reversed decoys collected zero reads; here they collected 42, in two
of 469 comparisons, confined to 64 bp of one decoy. It changes no detection
call. Writing "zero" because zero was the ported answer is exactly the failure
the rule exists to prevent.

## 14.4 Refusals

**A refusal is not a disagreement.** This is the rule that appears in the most
instruments. A structural test declines six models; their agreement column is
left **blank, not zero**, because scoring a refusal as a disagreement would
turn *this model is too partial to place* into *shape contradicts the census*.
The tree's unplaced tips, the sweep's uncontrolled genomes and the census's
unassigned records are the same rule in earlier instruments.

The piece that makes a refusal checkable was added late: the report prints the
**arithmetic ceiling** each declined structure could have reached, so "too
short to pass" becomes a claim that can be falsified. Here it was: every one
of the six had the headroom and fell short on similarity instead.

**A verdict vocabulary needs a word for a measurement that did not happen.**
This project's is five-valued — confirmed, contradicted, not corroborated,
orthogonal, underpowered — and the last three do the work that stops a
comparison being overclaimed. *Not corroborated* is a related property
measured cleanly that does not support the prior. *Orthogonal* is a different
property of the same object: a tree's branching order against a retained
flanking ohnologue, a column's dispersion against a rate on a tree, a search's
recovery rate against a gene's evolutionary rate. And *underpowered* is not a
polite word for a negative — a neighbourhood comparison taken where 73 % of
the flanking genes are unnamed is a measurement that did not happen.

**A ratio is estimated on a tree, never from a pair, once the pairs are
saturated.** The counting method a pairwise estimate rests on [R152] is used
in the application's own exploratory tool and nowhere in this thesis's
results, for the reason Chapter 10 measures: nearly nine tenths of
within-paralogue pairs are past the saturation bar, so a pairwise matrix is a
diagnostic here and not an estimate.

**A significant test is not a claim; the fitted parameter is.** Three of the
selection chapter's significant tests mean something other than what their
p-values suggest, and in every case the tell was a parameter printed *outside*
the test: an effect size pinned at the optimiser's bound, a rate class at
exactly one with a share of zero, a likelihood-ratio test with no site
clearing its own posterior threshold. Every significant test in this thesis is
printed beside the parameter it is about.

**A likelihood-ratio test whose null sits on a boundary needs the right
null distribution**, and both are written out rather than one chosen.

## 14.5 Reproducibility, and three things that were not reproducible

**Anything a committed artefact is built from is aligned single-threaded.**
The alignment program's [R138] iterative refinement combines partial results in
whatever order the threads finish, so at automatic thread count it is not
reproducible: the same seeds aligned twice gave 8,510 and 8,468 columns, and
the profiles built from them 4,933 and 4,908 match states. This was discovered
by rebuilding profiles after an unrelated edit and noticing the match-state
count had moved. **A profile that changes when it is rebuilt cannot be the
profile a committed result was produced with.** Single-threaded costs 87
seconds instead of 15, once.

**Thread counts and seeds are pinned for the tree search [R140] too**, because
automatic thread selection reads the machine's current load and the search is
reproducible only at a fixed seed *and* a fixed thread count. The number
pinned is what the program's own benchmark returned, so pinning costs nothing.

**A database release is pinned like a parameter.** The paralogy map is taken
from a dated archive rather than the rolling host, because the rolling host
follows the release cycle and a re-run would score the same windows against a
different tree. The archive's own registry is committed beside the map.

**A figure is not reproducible until its file is.** Figures were being written
with a creation timestamp in the vector output, so the same code on the same
data produced different bytes and no checksum recorded against a figure meant
anything. The timestamp is now dropped at save time. A mutation test written
to catch this initially *passed on the broken code*, because both saves landed
in the same second.

**A failed run is never cached, and a running stage owns a lock.** An
interrupted structural run left its parent process orphaned while its children
were killed; it went on spawning work beside the driver that replaced it, and
the killed pairs were written to the cache as zero-score results with empty
error text — ten permanent false negatives that only the stage's own
self-test caught. A result that is not marked successful is no longer written,
and a stage claims a lock checked against a *live* process, so a crash does
not block the next run for ever.

**And a cache may hold what a parser found, never what a rule decided.** One
audit cached both, and a later rule change survived into a committed table.

## 14.6 Making a report unable to lie

**Reports are rendered from the committed tables, never written by hand
alongside them.** Every `report.md` in this project's results tree is
generated, and the numbers in it come from the tables in the same directory.
That is why the claims ledger can use a text search against a report as a
check that is not weaker than a table lookup: the report carries no
hand-written numbers.

**Headlines are chosen by the data.** Where a chapter tests a prediction an
earlier chapter made, the prior is stated in code with where it was said, the
new statistic is computed, and the verdict is rendered from the comparison —
with both numbers printed either way. A generator written to narrate the
expected answer would print it whatever the data said, which is how a pipeline
launders an assumption into a result.

**A comparison that has to be sayable is the one that goes badly.** Several
chapters here report a result against their own design: the constraint chapter
finds the deep layer it was built around is the second best of four; the
ligand-site chapter finds the module that names the family is the less
constrained of two; the annotation chapter finds the family is not recorded
worse than its sister. Each of those had to be as easy to print as its
opposite, and the report generators are written so that it is.

**Every load-bearing number is declared with the table it comes from and the
operation that recovers it.** The manuscript carries 276 such declarations and
this thesis carries its own, through the same engine rather than a fork, so a
number quoted in both documents is recovered once and cannot disagree between
them.

![](figures/census_growth.png)

**{fig:census_growth}.** How the census grew across its six editions, and on
what evidence each addition rests. It is in the methods chapter rather than in
a results chapter because it is a statement about the search rather than about
the family: each edition is a channel, and the figure is the shape of what
each channel was worth. The targeted database search that started the project
returned about a tenth of what the exhaustive enumeration did, and the two
genome sweeps together added the records no protein database holds at all.

**Citations are keys resolved at build time and the bibliography is rendered,
never typed.** No author, title, year or journal is written into a chapter
file. A cited key with no reference row is a build error. This is a decision
rather than a convenience because a bibliography typed beside the text is the
one part of a document no table check can reach.

**A typeset build is not finished until its own log has been read for dropped
glyphs.** A missing character is a warning to the typesetter and a hole in the
page to a reader. One early build printed the receptor's own name without its
subscript because the document font has no subscript glyphs, and nothing in
the build said so.

## 14.7 The rules that are specific to this family

Four of the seventy-five are about this family and could not have been ported
from anywhere.

**Separating IP₃ from ryanodine receptors is a positive test at every stage.**
Never assume a domain hit, a similarity hit or a gene model is a family
member. Assign by best profile or by a labelled-bait margin, record the
margin, and treat the length band as a filter that supports the call rather
than as the call itself.

**A paralogue label is only meaningful inside the vertebrates.** ITPR1, ITPR2
and ITPR3 are a whole-genome-duplication product, so a protist record whose
protein name reads "receptor type 2" carries that number by annotation
transfer, not by descent. Four representatives were affected.

**And a bait attribution is only a label where it varies.** The sweeps record
which bait won each locus, and for most clades that discriminates because the
three baits win different loci. In the cyclostomes it does not: each genome
carries three full-length loci and the same bait wins all of them. A constant
attribution is no information, and filling a paralogue grid cell from it would
assert an orthology the data does not contain. The rule counts distinct
attributed cells per band, so the exception is *found* rather than named.

**In a family whose members are 61 to 68 % identical, "the same query aligns
twice" measures paralogy, not duplication.** Run on raw alignments, a
duplication detector fires in essentially every vertebrate genome because
every bait aligns at all three genes: specificity 0.16. The fix is the family
call applied to a pairwise test — the pair must sit inside the cell's own loci
— and specificity goes to 0.977 against a copy call the detector never sees.
**A geometric signature is only evidence once the family assignment has been
made.**

## 14.8 Two rules about scope

**An absence claim must pass two bars, not one:** a genome-wide contiguity
floor *and* a local check that the region itself is present. Residual false
negatives are regional rather than random, and a genome missing a whole
syntenic block is not evidence of gene loss.

That bar was chosen a priori from gene geometry — the median measured genomic
span of the gene, with no error rate anywhere in its derivation — and has now
been scored against a measured false-negative rate for the first time. It
gives 0.9 % residual error on the family series and 0.0 % on the sister
series, and it **stands**. The cost is reported beside it: it removes 38.8 %
of the scope, disproportionately from the clades a loss survey most wants.

**A search's sensitivity is measurable, not estimable, when the family has no
losses.** That is the design that makes Chapter 4's second half possible, and
it needs three rules to be a measurement rather than a restatement: the cells
the state rules decline to call are in neither the numerator nor the
denominator; the sister cell is carried as an independent replicate using none
of those state rules; and a false negative must be *constructible*, because a
control that can only ever return zero misses is not a measurement.
