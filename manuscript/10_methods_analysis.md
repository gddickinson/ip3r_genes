## Methods — alignment, phylogeny and downstream analysis

**Representatives and alignment.** 134 representatives were chosen from the
18,065-record census by eight coded rules, per clade and per kingdom, never
longest-per-species: a longest-first pick reliably selects chimeric gene
models. No rule ranks or filters on sequence identity, because identity was
retired as a call gate outside the vertebrates; a paralogue label is
vertebrate-only, because four non-vertebrate records carry a type number by
annotation transfer; and each group's length target is its own measured
median rather than a global one, since the family's size varies far more
outside the vertebrates than within. MAFFT L-INS-i (`--thread 1`) gave
11,777 columns and trimAl `-automated1` kept 1,797 (15.3 %, 7.7 % gaps),
96.8 % of them parsimony-informative. A ragged alignment is a hard failure:
a silent MAFFT failure degrades to a star alignment with no other symptom.

**Phylogeny.** ModelFinder's exhaustive scan was started, measured at 11 of
up to 1,232 models in 11 min, projected at ~21 h and abandoned for a
two-stage greedy scan with both stages committed; `Q.insect+R7` was chosen
(the free-rate models are not optional here — `LG+R5` beats `LG+I+G4` by 682
BIC units). IQ-TREE 2.3.6 with 1,000 ultrafast bootstrap and 1,000 SH-aLRT
replicates, threads and seed pinned, gave log-likelihood −215,452.0. Sister
hypotheses were tested with the approximately unbiased test over 10,000 RELL
replicates on three constrained maximum-likelihood trees, each constraining
only the taxa its hypothesis is about — a constraint naming every tip pins
tips outside the groups it declares and makes all three hypotheses fail
together. The search was repeated with an extra
nearest-neighbour-interchange round per bootstrap tree as a model-violation
guard, and every clade claim re-asked of that tree from the same committed
tip sets.

**Synteny.** Flanking genes were read from each assembly's own annotation
under two window rules (±10 coding genes; the 10 nearest *informative*
symbols each side) and three symbol vocabularies, the third stripping
trailing digits so that a two-round ohnologue pair is visible at all. Every
real locus pair is scored against matched random-window control pairs drawn
in the same two genomes, seeded per accession, and compared by a sign test
with ties dropped and counted. The paralogue caller's operating point was
chosen by maximising call rate minus random-window false-call rate, not call
rate alone.

**Duplication.** Human paralogy came from a pinned Ensembl BioMart archive
(Ensembl Genes 116) with the archive's registry committed beside the map, so
a re-run scores the same windows against the same gene tree. Every IP₃ *and*
ryanodine receptor gene is removed from every window before testing, and the
null is a permutation over *real* genomic windows at the matched gene count.
Links carry their Compara duplication node, which is what separates a
two-round ohnologue pair from an older duplication whose copies happen to
lie in these blocks. The 309-genome replication uses the synteny task's own
control sampler unchanged.

**Reconciliation.** The species tree is an input, not a result: 31 species,
29 named internal nodes, each with a published crown age, the spread of
published estimates, a stem age and its source, validated by a checker that
fails on an uncalibrated node, an age inversion or a species mismatch.
Duplications were called by the non-binary LCA rule and losses per Zmasek &
Eddy; the familiar binary rule reads a support-collapsed four-way node as a
speciation. Five topologies × three tip variants were reconciled, and the
vertebrate subtree re-rooted at all 112 of its edges.

**Gene architecture.** Exon blocks come from the same spliced alignments the
census was built on, read in target order so a minus-strand gene's first
block is the one carrying residue 1. What counts as an intron was calibrated
against the genome rather than chosen: block gaps were binned and scored by
their splice dinucleotides, and the floor was put at the smallest gap the
genome calls spliceable. Intron positions are a pair — the alignment column
the upstream exon ends in, and the classical phase — carried into one frame
through a bait-to-human-reference-to-column map whose anchor test requires
all 14 residues measured on the reference structure to land in one column in
all three paralogues. The shared-intron null is exact: each intron of one
gene placed independently and uniformly on the columns both loci resolve,
keeping its own phase, which makes the match count a Poisson-binomial whose
upper tail needs no simulation; a seeded permutation is run beside it and
both are committed. Every comparison is paired inside one genome.

**Selection.** Every vertebrate tip carries a coding sequence that provably
encodes the exact protein aligned: three independent routes each generate
*candidates*, and the first that validates by translation wins, because a
UniProt entry lists every transcript of its gene and in both human *ITPR1*
and *ITPR2* the first is not the aligned isoform. Every translation
disagreement is masked to `NNN`. PAL2NAL output was cross-checked
nucleotide-by-nucleotide against an independent in-house codon mapping, and
trimAl columns chosen on the protein were applied codon-aware, whole
triplets only. Selection sets are the tree's own extended paralogue clades,
cross-checked against the committed clade table, because a foreground that
is not a clade does not fail — it silently marks a larger one. Branch-site
model A was restarted from four initial ω by construction; the branch-site
likelihood-ratio test is reported against a 50:50 mixture as well as χ²₁;
and Benjamini-Hochberg correction runs across the whole 12-test family.

**Constraint and variants.** Conservation was computed per residue in each
human paralogue's own numbering on four layers, sequence-weighted before any
column statistic, with gaps excluded from the distribution and the
occupancy each score is conditional on reported beside it. Pfam element
coordinates were measured per accession and not transferred; the
structural elements (selectivity filter, gate, IP₃ contacts) were
transferred from PDB 6DQN with an anchor test that aborts rather than
writing a coordinate. The luminal loop has no annotation and was located
geometrically from the structure's own membrane span. ClinVar records were
placed only after each cited transcript's own translated coding sequence was
aligned to the canonical and the reference amino acid required to match.
Classifier performance is reported on three contrasts, the decisive one
restricted to the positions every layer scores.

**The ligand site.** The two functional modules are defined here rather than
taken from Pfam, twice each: a primary definition derived from measurement
(the ligand core as the minimal span holding every measured IP₃ contact; the
pore as PF00520 less the geometrically located luminal loop) and a
sensitivity definition taken from how the field draws the same region (the
published binding core; PF00520 as InterPro draws it). Every test is run
under all four combinations. Each definition is checked against something it
does not contain — the core must hold all ten contacts and no filter or gate
residue, the pore both filter and both gate residues and no contact — and a
definition that fails raises rather than being written. The pocket was
measured all-atom, not on Cα, within 15 Å of the ligand in six independent
IP₃-bound depositions, with recovery of the ten published contacts required
as a positive control. The module comparison is paired per orthologue, one
core and one pore number per sequence, with a tip required to resolve half
of *both* modules to enter. Phospholipase C presence requires both halves of
the catalytic barrel in one protein; the lineage comparison is matched on
pore identity within 0.03, and its power is reported against the shift the
ryanodine receptor control produces on the same instrument.

**Structures.** References were resolved by an RCSB query over the union of
the family's Pfam signatures — a PF08709-only query misses this project's
own IP₃ receptor reference — and assigned to a family by the census, never
by entry title. Every structure is reduced to its largest chain before
scoring, because both families are homotetramers and an assembly is ~11,000
residues. AlphaFold DB coverage is reported three ways, and only the third
— resolved residues over the census length — detects that the database
answers a canonical accession with an isoform (181 residues of a
2,701-residue *ITPR2*).

**Expression.** Per-species references are spliced *genomic* coding
sequence, deliberately not frameshift-corrected, with three housekeeping
anchors built by the same pipeline and a reversed, composition-matched decoy
for every sequence. Runs were selected by keyword and verified against each
run's own sample attributes. Alignment ran with spliced alignment disabled,
because the reference is already coding sequence and a junction read is
contiguous in it.

**Negative controls.** Every stage from the bait screen onward runs a suite
of constructed negative controls before it writes anything, and refuses to
write if one fails: 11 for synteny, 12 for the codon alignment, 14 for
constraint, 15 for the loss instrument, 21 for the loss counts and 21 for
gene architecture, 31 for duplication, 32 for the methods results, 44 for
the annotation audit and 44 for the ligand site, and 47 for the
annotation-bug validation. Each suite was itself mutation-tested by
breaking a rule deliberately and requiring the suite to catch it. Several
caught real defects, including a within-protein control that had swept in
the 225-residue domain the ligand question is about, a hash-seeded
non-determinism in a ranked table, and a self-test that was overwriting the
committed table it was checking.

**Reproducibility.** Reports are rendered from committed tables by script
and contain no hand-written numbers; figures are drawn only from committed
tables and are copied, never re-plotted, into the manuscript; and every
load-bearing number in this paper is declared in a claims ledger with its
source table and the operation that recovers it, re-verified on every build
(`manuscript/claims_check.tsv`).

