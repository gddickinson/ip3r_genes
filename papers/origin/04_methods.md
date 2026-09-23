## Methods

**Scope:** Four denominators are used, and each is declared here. The tree is
built on 134 representative proteins chosen from a census of 18,065 IP₃ and
ryanodine receptor records assembled from reference proteomes, InterPro and
genome sweeps (the census and its family call are described in
{paper:range}). The neighbourhood analyses use every IP₃ and ryanodine
receptor locus recovered by a protein-to-genome sweep of 309 vertebrate
genome assemblies, one per vertebrate order plus every species the protein
record left in doubt: 2,144 loci, of which those in the 274 assemblies that
carry a gene annotation have a readable neighbourhood. The reconciliation
uses the 57 vertebrate tips of the tree and a species tree of 29 named
internal nodes. The gene architecture uses the 1,378 genes in the 189
assemblies whose contigs are long enough to carry the whole gene.

**Representatives:** Eight coded rules choose the set: a vertebrate grid of
clade band by paralogue with human forced in, the teleost co-orthologue pairs
as a stress test, the deep vertebrate grade (cyclostomes, cartilaginous fish,
coelacanth) entered unlabelled because a paralogue label is a vertebrate
concept, non-vertebrate metazoa by phylum, non-metazoan eukaryotes gated on
per-record contamination verdicts, copy-number expansions, gene models no
database holds, and the ryanodine outgroup. The quality key is lexicographic,
so the audit names the component that decided each pick, and each group's
length target is its own measured median. Six selection rules have
constructed negative controls that run on every build.

**Alignment:** MAFFT L-INS-i [R138] was run single-threaded, because at an
automatic thread count its refinement is not reproducible on this machine,
and the SHA-256 of input and output are recorded. A ragged alignment is a hard
failure, since a silent MAFFT failure degrades to a star alignment. trimAl
`-automated1` [R139] chose the kept columns; the heuristic was fixed in
advance rather than tuned against the resulting tree.

**Tree:** An exhaustive ModelFinder scan [R141] was started and abandoned
after it projected about 21 hours, and replaced by a two-stage greedy scan
(every matrix at fixed rate heterogeneity, then every rate model on the
winning matrix) with both tables committed. Free-rate models were kept
because they beat gamma models by 682 BIC units on this alignment.
IQ-TREE 2 [R140] ran with 1,000 ultrafast bootstrap [R142] and 1,000
SH-aLRT [R143] replicates, threads and seed pinned. A node is well supported
only when it clears both SH-aLRT ≥ 80 and UFBoot ≥ 95. Paralogue membership
was corrected by a coded rule that keeps a tip inside its group's largest
pure clade, moves a tip only into a well-supported clade of at least three
members of one other paralogue, and otherwise leaves it unplaced; every tip
the tree disputes was checked by reciprocal best hits against the human
proteome.

**Sister test:** Each of the three rooted pairings
was a constrained maximum-likelihood search under the same model, with the
constraint naming only the three paralogue cores and the outgroup. A
constraint that names a tip pins it outside every group it declares, and a
first version that named all 134 tips rejected every pairing, including the
one the tree holds at maximal support; a constructed control now fails any
constraint that names a free tip. The approximately unbiased test [R144] used
10,000 RELL replicates, and the 95 % confidence set of topologies was read off
it. As a model-violation guard the search was repeated with `--bnni`, and
every clade claim was re-asked of that tree from committed tip sets.

**Neighbourhoods:** Flanking genes were read from each assembly's own
annotation under two window rules (±10 coding genes, and the 10 nearest
informative symbols each side, since human annotation names about 97 % of
coding genes and sea lamprey about 27 %) and three symbol vocabularies, the
third stripping trailing digits so that an ohnologue pair is visible at all.
Loci were taken from the per-genome sweep output rather than a per-cell
table, so that a teleost cell holding two genes contributes both. Every real
locus pair was scored against matched random-window pairs drawn in the same
two genomes, seeded per accession, and compared by a sign test with ties
dropped and counted. Two empty neighbourhoods score zero, not one, so that two
unannotated genomes cannot form a perfect match; this is one of 11
constructed controls run on every build. The paralogue caller builds a
consensus per paralogue with the query's own species removed, and its
threshold maximises call rate minus the random-window false-call rate; a call
at or below the highest score any random window reached is reported as inside
the null.

**Dated paralogy:** Human paralogy came from a pinned Ensembl BioMart
archive (Ensembl Genes 116) [R161, R162], with the archive's registry
committed, so a re-run scores the same windows against the same gene trees.
Every IP₃ and ryanodine receptor gene was removed from every window. The null
for the human test is a permutation over real genomic windows at the matched
gene count, 5,000 permutations at each of three window sizes, which keeps the
gene-family clustering a shuffled gene set discards. Each link carries its
Compara duplication node. The replication across 309 genomes applies the
human map to the flank sets of every swept genome and compares against matched
random neighbourhoods drawn by the neighbourhood analysis's own sampler,
unchanged. The ryanodine receptors pass through every step as the positive
control, since their origin in the early duplications is not in question.

**Reconciliation:** The species tree covers the 31 sampled vertebrate
species. Every internal node carries a crown age, the spread of published
estimates, a stem age and a source, and a validator fails on an uncalibrated
node, a point age outside its own spread, an age inversion or a species
mismatch. Where the
literature does not resolve an arrangement the node is left a polytomy. The
stem age matters because the sampled tree omits unsequenced lineages, so a
node's parent here is usually older than the split that created it.
Duplications were called by the non-binary least-common-ancestor rule [R174]
and losses counted after Zmasek and Eddy [R173]. Five topologies were
crossed with three tip treatments (all tips, cyclostome tips pruned, nodes
below the support bar collapsed), and the collapse refuses a tree with no
support labels, since a constrained search writes none. As a rooting control,
the vertebrate subtree was re-rooted at each of its edges and scored by total
events beside the outgroup rooting. Long-branch attraction was assessed by
root-to-tip distance over all vertebrate tips. Every implied loss was checked
against the genome sweep. Fourteen groups of constructed controls run before
anything is written.

**Teleost copies:** A copy is a locus whose model covers at least half of
its bait over at least 500 aligned residues at 0.4 identity or more, and
counts were recomputed at seven coverage bars. Copy dispersion was measured on
the genome. Double-conserved synteny was scored against a tetrapod and a
pre-duplication ray-finned reference, because teleost symbols diverge from
tetrapod ones. Block consistency was tested by letting anchor genomes from
different orders assign every genome's two copies independently, against a
two-sided binomial at one half; a constructed control builds independent
lineage-specific duplications and requires the statistic to land at one half.

**Gene architecture:** Exon blocks were read from the sweep's spliced
alignments in target order, with model keys disambiguated across the two
assemblies processed in chunks and splice dinucleotides read
strand-corrected. The intron floor was calibrated by binning block gaps and
scoring them against the genome's splice dinucleotides. Intron positions were
carried into the representative alignment's columns through a map from each
bait to its cell's human reference; the map is refused unless all 14
structurally measured anchor residues land in one column in all three
paralogues. The shared-intron null is exact: each intron of one locus is
placed independently and uniformly on the columns both loci resolve, keeping
its phase, so the shared count is Poisson-binomial and its tail needs no
simulation. A seeded permutation without replacement is run beside it. Every
cross-paralogue comparison is paired inside one genome, sign-tested with ties
dropped, and Benjamini-Hochberg corrected across the family of tests. The
ryanodine receptors pass through the same steps as the control. Twenty-one
constructed controls run before anything is written, and six deliberate
breakages of the rules were all caught.

**Reproducibility:** Every table is written by a script and every report is
rendered from committed tables. Every load-bearing number in this paper is
declared in a claims ledger with its source table and the operation that
recovers it, and is re-verified on every build.
