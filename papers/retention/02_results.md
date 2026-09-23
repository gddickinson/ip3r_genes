## Results

### The sweep searched every genome in a declared scope with a control that never failed

The scope is 309 vertebrate assemblies, 552,535,418,825 bp in total: one
best assembly for each of 161 vertebrate orders, and 169 margin species that
the protein record left in doubt (a species that is both is one row). The
margin species were derived from the protein census by stated rules rather
than listed by hand: their reference proteome returned no family record,
carried fewer than three paralogues, held only fragments, or was an anchor
species. They are therefore enriched, by construction, for genomes in which
a paralogue looks missing, which is the population a loss survey most needs
and the one in which a false absence is most likely.

A panel of 38 baits (30 IP₃ receptors and 8 ryanodine receptors, 121,294
residues) was aligned to every genome, and cells the alignment left empty
were re-asked with a whole-genome translated search [R158]. All 309 genomes
completed. The ryanodine receptor control, which should fire in every
vertebrate, failed in 0 of them, so no genome is excluded on control
grounds. The sweep recorded 1,236 genome × cell results over 2,144 loci
({fig:census}; {fig:status}). At every one of the 2,144 loci only one
family's baits aligned, so the family question is settled by the alignment
before any threshold is applied.

Read as a raw ledger, the result is already almost complete. Four cells in
1,236 were called absent, meaning no aligned locus and no remnant, and all 4
are in the two cyclostomes, *Petromyzon marinus* and *Myxine glutinosa*,
which lack *ITPR2* and *ITPR3* in the ledger. The rest of the missing
evidence is not absence but weaker presence: truncated loci, fragments and
translated-search traces, concentrated in particular vertebrate classes
({fig:census}). The remainder of this paper asks whether any of that is a
loss.

### The search misses a gene only where the assembly is fragmented

The central control of any absence claim is how often the same search fails
to find a gene that is present. We measured it on the search's own output.
Anticipating the result below, that no loss is reconstructed anywhere in the
scope, every cell whose gene is independently known to be present and which
the ledger did not grade as found is a false negative of the method, a miss
with a known answer. Over the whole scope that is 140 of 923 cells, a rate
of 0.1517 ({fig:false_negatives}a).

That series depends on the state rules of the next sections, so we carried a
second one that does not. Every vertebrate has ryanodine receptors, so every
genome's control cell has a known answer that no part of the state logic
produced. It misses 42 of 309, a rate of 0.1359, indistinguishable from the
IP₃R series (Fisher p = 0.578). The miss rate is a property of the method
applied to these assemblies, not of the paralogue calls.

Every miss tracks assembly contiguity. A missed cell's assembly has a median
contig N50 of 23,460 bp against 3,396,515 bp for a found one; the odds of
finding the gene rise 8.10-fold per tenfold increase in contig N50 in the
IP₃R series and 19.99-fold in the ryanodine series; and chromosome-level
assemblies miss 3 of 512 IP₃R cells and 0 of 172 ryanodine cells
({fig:false_negatives}a; {fig:confound}). Below the contiguity bar the miss
rate is 0.3917 for *ITPR1*, 0.4333 for *ITPR2* and 0.3 for *ITPR3*.

The contiguity bar was fixed before any error rate was known: 142,212 bp of
contig N50, the median measured genomic span of the gene, on the grounds
that an assembly whose contigs are shorter than the gene cannot represent it
on one contig. Scored against the miss rate afterwards, it holds. Above it,
189 genomes carry 563 IP₃R control cells, of which 5 are missed (a rate of
0.0089), and the ryanodine series misses none of 189
({fig:false_negatives}b). The bar has a price, which is stated rather than
absorbed: it sets aside 120 of 309 genomes, disproportionately birds and
margin species ({fig:false_negatives}c). Every count below is therefore
given over all 309 genomes and over the 189 above the bar.

### The bait panel's recall depends on paralogue coverage, not on phylogenetic breadth

The spliced aligner treats each bait independently, so dropping baits from
the retained alignments and re-running the cell-assignment chain reproduces
exactly what the sweep would have reported with a smaller panel. The
simulation reproduces the committed ledger in 1,236 of 1,236 cells, with 0
mismatches, so each of 19 panels, the full one and 18 ablations, is measured rather
than modelled
({fig:panel}).

The full panel recovers 783 cells. Four human baits, one per cell, recover
782 of them, and removing any single vertebrate clade band changes the count
by at most two cells in either direction. Removing one paralogue's own baits
is what costs recall: the deltas are -238, -240 and -241 cells for *ITPR1*,
*ITPR2* and *ITPR3* ({fig:panel}a). A single bait recovers its gene at any
identity above 0.5 ({fig:panel}b). Removing the ryanodine receptor baits
changes no IP₃R call, so the positive family test costs nothing in recall.
Of the 2,179 call changes across all panels, 68 are gains, each produced by a
marginal cell whose coverage is read from its top-scoring bait and which a
competitor's removal pushed across the coverage bar. They are properties of
the assignment rule, not evidence that a smaller panel searches better.

### No genome × paralogue cell reaches the state that licenses a loss

The ledger's own statuses were not designed to license a loss count. There,
"absent" means that the rescue search attributed no region, which in a
shattered assembly is a statement about contig lengths. We therefore
restated each of the 927 genome × paralogue cells under eight ordered rules,
in which the first rule to fire assigns the state and exactly one state is
countable as a loss. A cell is uninformative if its genome's control failed.
It is present if it carries a locus at or above the coverage bar; truncated
if its locus ends at a contig edge or a run of ambiguous bases; partial if
its locus is short with no assembly excuse; fragmented if the reference
reassembles from sequence outside every locus the aligner placed; and
paralogue-unassignable if the genome carries family loci that no cell
claimed. Only then can a cell be undecidable because its assembly cannot
hold the gene on one contig, or absent. The absent state is reachable by
construction: a self-test builds a contiguous, controlled, empty, spare-free
cell and requires it to come back absent.

Of the 927 cells, 783 are a single full-coverage locus, 90 a locus truncated
by the assembly, 7 a partial locus, 43 a gene reassembled across contigs and
4 paralogue-unassignable. None is absent, and 923 hold a gene that is present
({fig:matrix}a).

The four unassignable cells are the ledger's four absences. Both cyclostome
genomes carry three full-length family loci, all of which fall to the
*ITPR1* baits because the panel holds no cyclostome-labelled *ITPR2* or
*ITPR3* bait. Filing those genomes' *ITPR2* and *ITPR3* cells as absent would
score two losses per cyclostome that the genomes themselves contradict. They
are counted neither as present nor as lost.

Genomic sequence, unlike reference coverage, can count copies, because a
placed locus and a reassembly occupy different places in the genome. Summed
per genome, every one of the 189 assemblies above the contiguity bar holds a
minimum of 3.0 gene-equivalents, all three paralogues. Below the bar, 11 of
120 fall short of 2.5, and the shortfall follows contig N50 rather than
taxonomy ({fig:matrix}b).

### Every undecided cell reassembles into one gene above a calibrated bar

The 43 fragmented cells are the ones a naive count would read as candidate
losses: the aligner placed no locus, but the translated search found
something. For each we measured the non-redundant coverage of one reference
protein reassembled across contigs, computed only from sequence outside every
locus the aligner clustered in that genome. Both restrictions are necessary
in a family whose paralogues are 61 to 68 % identical: without the exclusion
set, a missing gene reassembles out of its own paralogues, and a union over
several baits reports a coverage that no single protein achieves.

The bar a reassembly must clear was measured against a gene known to be
elsewhere. The negative population is the set of regions the full panel
attributes to a paralogue whose gene the aligner already placed at a locus in
the same genome, so the fragments cannot be that gene and measure only what
cross-paralogue similarity delivers on its own. The 44 positive regions
reach a median coverage of 0.7948 and the 9 negative ones 0.0303
({fig:reconstruction}a). The first version of this calibration did not
separate, and the reason is itself a result: in a genome where two
paralogues are both shattered, a fragment attributed to the other one is a
piece of a real gene. Those 20 co-shattered regions, at a median of 0.6308,
are now reported as their own population, the measured size of the
paralogue-attribution problem in a fragmented assembly. With them separated,
Youden's index [R177] is 1 and the gap runs from 0.1 to 0.1737. The
operating point is the gap's midpoint, 0.13685, rather than Youden's own
threshold, which on a perfectly separated pair is the most permissive bar the
data allow. Every undecided cell clears the bar, whether its gene is spread
over two contigs or twelve ({fig:reconstruction}b).

We also asked whether flanking-gene synteny could independently place these
fragments, using a neighbourhood caller that is correct at every key count at
which it acts. It could not reach them. Of 432 rescue regions, 273 sit on a
contig that carries no annotated gene at all and 8 reach the caller's
four-key floor; the median region's contig extends 39.936 kb, shorter than
the gene ({fig:synteny_reach}). Synteny is accurate here and unavailable, a
property of the assemblies rather than of the method.

### Parsimony places no loss, and the settings that would manufacture one are named

Under the operating point, on the family-level coding and on each paralogue,
the character is present in every genome that can be scored, and Dollo
parsimony [R175] over the NCBI-taxonomy tree of all 309 assemblies places 0
losses. The count is zero because no cell reaches the loss state, not because
the routine cannot return anything else: constructed controls put a loss on a
known edge and require it found, and require two sister losses merged into
their parent. The tree is an input, not a result, and its 468 internal nodes
include 49 polytomies, which bound any count of independent losses from
above but cannot create one where there is none.

With no loss to place, the informative quantity is how far the rules must
move before a loss appears. We walked 32 settings: an ordered eight-rung
evidence ladder, each rung named after what it refuses, crossed with the
cyclostome rule and the contiguity rule each switched on and off, on both
codings and under three branch-length schemes ({fig:sensitivity}). The
first three rungs are the calibration's own gap edges and midpoint, so the
first question is whether moving the reconstruction bar across its measured
uncertainty changes anything; it changes no cell in 927
({fig:reconstruction_bar}).

The two codings then behave differently. The family-level coding
manufactures a loss in 2 of 32 settings, at worst 1 genome in 309, and only
when everything except a complete locus is refused and the contiguity bar is
ignored at the same time ({fig:sensitivity}a). The paralogue-resolved coding
manufactures one in 18 of 32 settings, up to 45 loss edges
({fig:sensitivity}b). A per-paralogue absence is fragile to every one of
these settings; a family-level absence is not. Branch lengths change the
count in 0 of the 32 settings, as parsimony requires.

Continuous-time Markov models of the character [R176] are stated and not
fitted. An invariant character contains no transition to estimate, and at
all twelve model, axis and branch-length combinations the likelihood profile
is monotone to the edge of an eight-order-of-magnitude rate grid
({fig:mk_profile}). No rate is reported for it.

### No reading frame marks a dead gene

Frameshifts and premature stops counted along each placed locus would, read
naively, be a pseudogene screen. We counted them per kilo-aligned-residue
and found that the confounder is not the obvious one: assembly contiguity
barely moves the lesion count, whereas a locus's identity to its bait moves
it a great deal ({fig:integrity}), because a distant reference buys its
alignment with frameshifts. Compared within each genome against the same
genome's other family loci and matched on identity, one paralogue shows an
excess: *ITPR3* carries more lesions than its identity-matched siblings in 39
genomes to 14, while the unmatched *ITPR2* excess disappears once identity is
matched.

That excess is not decay. Of 1,760 scored loci, 44 clear the lesion bar, and
all 44 are full-coverage loci with an intact gene model; 7 family loci fire
any of three deliberately generous fossil readings, each on one or two stops
in a complete model. There are therefore no dead loci on which a shared-lesion
test of an ancestral pseudogene could be run. Stratified by class, the *ITPR3*
excess is a bird result, 25 genomes to 2, but 21 of those 27 pairs lie in
assemblies below the contiguity bar, and above it 6 pairs cannot carry a test
({fig:lesion_strata}). We report it as unresolved rather than as a lineage
in which the gene is decaying.
