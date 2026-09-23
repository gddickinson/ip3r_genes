## Methods

### Scope

The denominator is 309 vertebrate genome assemblies, declared before the
search and fixed. It is the union of two rules applied to NCBI assembly
metadata [R159]. The taxonomic rule takes one best assembly per vertebrate
order, 161 orders, ranked by a stated function that prefers an annotated
assembly, then RefSeq over GenBank, then assembly level, then scaffold N50.
The margin rule adds 169 species computed from the protein census by four
positive tests: the species' reference proteome returned no family record,
carried fewer than three paralogues, held only records below the family's
minimum length, or the species was an anchor needed for another reason. The
same rank function resolves each margin species to an assembly, and a species
selected by both rules is one row. The manifest records accession, taxonomy,
assembly level, annotation source, contig and scaffold N50 and the rules that
admitted each genome. The scope totals 552,535,418,825 bp. The deliverable
unit is the genome × cell: three IP₃R paralogues plus the ryanodine receptor
control cell in each genome, 1,236 cells, of which the 927 paralogue cells
enter the loss analysis.

### Bait panel and paralogue labels

The panel holds 38 baits, 30 IP₃ receptors and 8 ryanodine receptors,
121,294 residues, chosen by seven rules: a bait's paralogue label comes from
the protein census rather than its own gene symbol; it is full-length with
the complete domain architecture; there is one bait per paralogue per clade
band, two in the bands that dominate the scope; every candidate passes a
chimera screen against the IP₃R and ryanodine receptor profiles; the
ryanodine baits serve as a presence control; existing profile seeds are
reused; and the panel is scoped to this sweep. Six slots have no labelled,
full-length record and stay empty (*ITPR1* and *ITPR2* in cartilaginous fish,
*ITPR2* in the coelacanth grade, all three in cyclostomes). Three unlabelled
baits from those grades keep their genes inside the family without assigning
a paralogue. The chimera screen is itself tested on every build with three
constructed failures (an IP₃R–ryanodine fusion, a mis-joined model and a
truncation), each of which must be rejected by the rule responsible.

Paralogue identity in this paper therefore rests on the labelled baits, and a
cell is called for the paralogue whose bait scores highest at its locus,
after the family has been settled by which family's baits align at all. No
cell call uses flanking-gene synteny or a gene tree. The family-level coding
of the loss analysis uses no paralogue label at all.

### Spliced alignment and rescue

Baits were aligned to each genome with miniprot [R164]. Its maximum intron
parameter decides whether a long gene is split, and a split IP₃R reads as a
fragment, so the parameter was measured rather than taken as the default of
200 kb. Across an eleven-species Ensembl panel [R160], the widest intron of
any IP₃R gene is 152,216 bp (human *ITPR1*) and of any ryanodine receptor
gene 227,927 bp (human *RYR2*). The sweep uses twice the widest measured
intron scaled by genome size, floored at the aligner's default and capped for
tractability; the two genomes on which the cap binds are flagged. Cells left
empty were re-asked with a whole-genome tblastn search [R158], and each
region was attributed first to a family and then to a paralogue against the
full panel. The paralogue-attribution margin arrived from a sister project as
0.333; measured at 542 loci whose paralogue the assembly's own annotation
establishes, 368 of them fell below it, so the margin in force is 0.22, the
highest value that rejects none of those calls. All 53 rescue regions that
overlap an annotated gene of known paralogue were attributed in agreement
with it.

### The positive control

The ryanodine receptors are present in every vertebrate and are searched in
the same run by the same aligner. A genome whose control cell is not found
has an assembly or pipeline problem, and its IP₃R cells support no claim.
The control failed in 0 of 309 genomes.

### States

Each paralogue cell was restated by eight ordered rules, first rule to fire
wins, each recording the number it fired on (Results). Coverage is aligned
bait length over bait length, with the bar inherited from the sweep. The
cyclostome rule reads the sweep's own per-genome locus counts: a cell with no
locus in a genome carrying family loci no cell claimed is
paralogue-unassignable. The contiguity rule is the bar described below. The
rule order is tested: moving the cyclostome rule after the absence rule must
fail the suite.

### Reconstruction across contigs

For each cell with no locus, each reference bait was taken in turn, never as
a union, and its non-redundant coverage by translated-search hits outside
every locus the aligner clustered anywhere in the genome (any bait, padded)
was computed, with the fraction of reference residues covered once
(tiling) reported beside the total (pile-up), and the genomic exonic length
on the subject side as gene-equivalents. The bar was calibrated against 44
positive regions (cells with no other explanation), 9 negative regions
(attributed to a paralogue already placed at a locus in the same genome) and
the 20 co-shattered regions kept as their own population. Youden's index
[R177] is reported with the gap's two edges, and the operating point is the
gap's midpoint. The calibration refuses to return a threshold on too few
regions or on unseparated populations.

### Synteny reach

Each rescue region was offered to a flanking-gene paralogue caller with a
four-key floor, whose accuracy was measured on its own leave-one-genome-out
calibration and binned by the number of informative keys, so that reach and
accuracy are read on one axis.

### Contiguity bar and false-negative rate

The contiguity bar is 142,212 bp of contig N50, the median genomic span of
28 IP₃R genes measured across the Ensembl panel. It was fixed from gene
geometry before any error rate was computed. Because no loss is
reconstructed, the IP₃R control series is every cell the states call
present, 923 cells (the four unassignable cells are excluded), and a
false negative is a cell in that series the ledger did not grade as found.
The ryanodine receptor series is the 309 control cells, which use none of the
state rules. For each series we report the rate with a Wilson interval, a
logistic fit of found against log₁₀ contig N50, a Mann–Whitney comparison of
contig N50 between found and missed cells, a Fisher test of chromosome-level
against lower assemblies and of one series against the other, and the
residual rate at every candidate floor with the genomes each floor retains
by clade. A constructed case requires a false negative to be reachable.

### Panel ablation

Eighteen reduced panels were simulated beside the full one by dropping baits from the retained
miniprot output and re-running clustering, the identity floor and the cell
assignment in full. Scoring best coverage over all alignments instead would
let one paralogue's bait stand in for another at the same gene, so the full
chain is reproduced. The full-panel simulation must reproduce the committed
ledger cell for cell, and it does in 1,236 of 1,236 cells. Every call change,
gains included, is recorded with both baits and both coverages.

### Tree, parsimony and rate models

The tree is NCBI taxonomy over all 309 assemblies, one tip per accession,
unary nodes kept so that the edge a loss is placed on does not move. It is an
input; its 49 polytomies are reported, not resolved. Losses were counted by
Dollo parsimony [R175] as the maximum (every loss edge) and the minimum (one
per parent carrying any), which differ only under a polytomy. The family
character's gain is pinned at the root and each paralogue's gain is left
free. The sensitivity grid is the eight-rung evidence ladder crossed with the
cyclostome and contiguity rules on and off, 32 settings, on the family and
paralogue codings and under unit, ultrametric and calibrated branch lengths.
Equal-rates, all-rates-different and irreversible Mk models [R176] were
profiled along a rate grid for every combination, and a fit is refused, with
its reason, when the character does not vary.

### Reading-frame integrity

Frameshifts and internal stops were counted per thousand aligned residues at
every placed locus. The bar is the upper quantile of a population the screen
never scores: full-coverage loci above the contiguity bar whose own
annotation names them as family. Density was correlated with contig N50,
bait identity and length. Each cell was compared with the same genome's other
family loci in a paired sign test, ties dropped and counted, once unmatched
and once identity-matched, and the tests were corrected together by the
Benjamini–Hochberg procedure [R170], then stratified by vertebrate class with
the contiguity bar applied within each stratum. The verdict is one-sided: a
locus without stops is not a pseudogene, and a handful of stops does not make
one.

### Negative controls and reproducibility

Every stage runs a suite of constructed negative controls before writing any
table, each of which must be rejected by the rule responsible (for example,
a reconstruction without the locus exclusion must reassemble a missing gene
out of its paralogues and fail). Everything downstream of the alignments is
offline and deterministic, and every committed table carries a SHA-256 in its
stage's statistics file.
