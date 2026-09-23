## Results

### A tree rooted on the ryanodine receptors makes *ITPR2* and *ITPR3* sisters

**The alignment is built to root the tree, not to find the family.** We chose
134 representatives from the census by eight coded rules, one per clade and
paralogue cell, with the non-vertebrate grades and six ryanodine receptors
included. No rule ranks on length or identity, because a longest-first pick
favours chimeric gene models. MAFFT L-INS-i gave 11,777 columns at 76.23 %
gaps, and trimAl kept 1,797 of them (15.3 %), 96.8 % of which are
parsimony-informative ({fig:alignment}; {fig:representative_alignment}). The
family separation is measured again on this alignment, as a check that the
outgroup is where it should be: mean identity is 0.910 within a vertebrate
paralogue group and 0.253 from each group to the ryanodine receptors
({fig:identity}).

**All three paralogues and the outgroup are clades.** Maximum likelihood under
Q.insect+R7 recovers the ryanodine receptors as a clade at the root and
*ITPR1*, *ITPR2* and *ITPR3* each as a clade inside the vertebrates
({fig:tree}). Of the 131 internal nodes, 91 (69.5 %) clear both SH-aLRT ≥ 80
and UFBoot ≥ 95, the joint bar every claim below is held to
({fig:support}a). The tree overturns no census name: 39 vertebrate tips agree
with their label, 0 are reassigned, and 5 are left unplaced. For each of the
five, a reciprocal best-hit search against the human proteome upholds the
database name, so the refusal to place reflects the tree's uncertainty rather
than a mislabelled record ({fig:support}b).

**The sister pair is *ITPR2* with *ITPR3*, and the other two pairings are
rejected.** Each of the three rooted pairings was realised as a constrained
maximum-likelihood search in which all three paralogue clades are
monophyletic, so the test compares the sister arrangement and nothing else.
Under the approximately unbiased test over 10,000 RELL replicates,
*ITPR1*+*ITPR2* costs 114.02 log-likelihood units and is rejected
(p = 1.8 × 10⁻⁵), *ITPR1*+*ITPR3* costs 113.98 and is rejected
(p = 1.65 × 10⁻⁵), and *ITPR2*+*ITPR3* is not rejected (p = 0.476), nor is
the unconstrained tree (p = 0.525) ({fig:sister}). The 95 % confidence set of
topologies therefore holds 2 trees, and both carry the same pair. Mean
identity had pointed the other way: *ITPR1* and *ITPR2* are the most similar
pair on the alignment ({fig:group_identity}). That ranking ignores the outgroup and rate variation,
which is why it was treated as a preview and the tree as the answer.

**The answer does not depend on the model.** Ultrafast bootstrap is
optimistic under model violation, and an alignment that spans four kingdoms
violates any single model. We re-ran the search with an extra
nearest-neighbour-interchange round on every bootstrap tree and asked the same
clade questions of the result, with clade membership read from a committed
table so that bookkeeping cannot masquerade as topology. Nine of the claim
clades are held and none is weakened or lost. The one clade that is poorly
supported in the main tree, *ITPR1*'s census-labelled core, is poorly
supported in both, while the same clade with the unlabelled tips included is
at maximal support.

**The cyclostome loci form two ancient lineages of their own.** Hagfish and
lamprey each carry three IP₃ receptor loci, and all six sit in
cyclostome-only clades rather than inside any gnathostome paralogue clade.
There are two such clades, each holding both species, so the cyclostome copies
are two lineages that separated before hagfish and lamprey did, not an
expansion either genome made on its own.

### The three neighbourhoods are paralogous, and the paralogy runs through *ITPR1*

**Each paralogue has a neighbourhood of its own.** We read the flanking genes
around all 2,144 IP₃ and ryanodine receptor loci in the 274 genomes that
carry an annotation, and scored every pair of loci against matched random
windows in the same two genomes, which holds annotation depth, naming
convention and window rule constant. Within-paralogue similarity runs 216- to
411-fold above its own null, and every cross-paralogue and cross-family class
sits at or below the null; over 168,241 pairs of IP₃ and ryanodine receptor
loci the largest class mean is 0.0002 ({fig:synteny}a). A paralogue caller
built on these neighbourhoods alone, with its operating point set by
maximising call rate minus false-call rate, is correct on 405 of 405 calls
over 503 loci whose paralogue the assembly's own annotation establishes, and
calls 6 of 726 random windows ({fig:synteny}b). It then places 131 loci that
the sequence sweep could not label. The same instrument shows that the
*ITPR3* neighbourhood, which lies in the major histocompatibility region in
human, is the one that decays fastest across vertebrate classes
({fig:synteny}c).

**Two families of flanking genes survive between the three neighbourhoods,
and both touch *ITPR1*.** Shared gene names are no test of whole-genome
duplication, because after about 500 million years the flanking copies are
paralogues of one another under different names. At the level of the gene
family, with the trailing digits of each symbol removed, two links survive
across the vertebrates: a metabotropic glutamate receptor family between
*ITPR1* and *ITPR3*, and a basic helix-loop-helix family between *ITPR1* and
*ITPR2* ({fig:neighbourhood}). *ITPR2* and *ITPR3* share nothing at any
threshold, and neither does any IP₃ and ryanodine receptor pair.

**The paralogy replicates across 309 genomes against a real-window null.**
To test paralogy rather than shared names, we took the human paralogy map from
a pinned Ensembl Compara release, removed every IP₃ and ryanodine receptor
gene from every window so that the family's own paralogy could not count as
evidence, and dated each link by its Compara duplication node. Asked of the
flank sets in every swept genome, the *ITPR1* neighbourhood carries a
paralogue of the *ITPR2* neighbourhood in 141 of 175 genomes (80.6 %, ten
vertebrate classes) and of the *ITPR3* neighbourhood in 89 of 152 (58.6 %,
nine classes). A null of 932 matched random neighbourhoods drawn by the same
sampler in the same genomes carries a link 2.6 % of the time. *ITPR2* against
*ITPR3* is 4 of 149 (2.7 %), which is the background (p = 0.554)
({fig:paralogon}a).

**Dating the links removes one of them from the two-round account.** The
glutamate receptor link between *ITPR1* and *ITPR3* is dated to Vertebrata,
which is what an ohnologue pair from the early duplications looks like, and it
is carried in 84 genomes. The helix-loop-helix link between *ITPR1* and
*ITPR2* is dated to Opisthokonta: it is real paralogy, far older than the
vertebrates, whose copies happen to lie beside these genes, and it is dated
to Vertebrata in only 2 genomes. The ryanodine receptors pass through the
identical instrument in the same genomes as the positive control, and in the
single human genome they fare no better: at the primary window each family
keeps one vertebrate-dated link between two of its three neighbourhoods, and
neither survives correction across every test the stage ran
({fig:paralogon}b). What the human genome alone measures is the instrument's
ceiling, which is why the replication carries the claim.

**The neighbourhood does not corroborate the tree's sister pair.** The tree
makes *ITPR2* and *ITPR3* sisters, and they are the one pair whose
neighbourhoods keep nothing above background. The two instruments do not
measure the same quantity. A tree estimates the order in which the copies
diverged. A retained flanking ohnologue records which copies survived
deletion beside each gene, and the quartets left by whole-genome duplication
lose their flanking copies lineage by lineage with no memory of the
duplication order [R205, R183]. We therefore report the disagreement as a
lack of corroboration and not as a contradiction.

### Reconciliation places the two duplications on different branches

**The species tree is an input with its uncertainty attached.** We mapped the
gene tree onto a hand-curated species tree for the 31 sampled vertebrate
species, with 29 named internal nodes, each carrying a published crown age, the
spread of published estimates, a stem age and its source
[R188, R189, R190]. Duplications were called by the non-binary rule [R174],
which is required because the familiar binary rule reads a
support-collapsed polytomy as a speciation, and losses were counted by the
standard method [R172, R173].

**The earlier split is on the vertebrate stem, the later one on the
gnathostome stem.** On the reported tree with every tip included, the split
separating *ITPR1* from *ITPR2* and *ITPR3* maps to the vertebrate stem, older
than crown Vertebrata. Published estimates for that crown run from 480 to
615 Ma around a point age of 563 Ma [R186, R187], and nothing in the sampled
tree bounds the duplication from above. The split separating *ITPR2* from
*ITPR3* maps to the gnathostome stem, between crown Gnathostomata at 462 Ma
and crown Vertebrata at 563 Ma ({fig:dated}). The two duplications sit on
different branches, and that makes *ITPR1* the earlier-diverging copy.

**The placement survives every perturbation we could make to the gene tree.**
We reconciled five topologies (the reported tree, the three constrained sister
hypotheses and the model-violation re-search) under three tip treatments,
twelve reconciliations in all. With the cyclostome loci kept, the deepest
paralogue duplication maps to Vertebrata in all 5 topologies, including the
two sister arrangements the AU test rejects, so the dating does not depend on
which pair is sister ({fig:reconciliation}a). Collapsing every node below the
support bar keeps the placement, even though the gene-tree node that unites
*ITPR2* and *ITPR3* with *ITPR1* is the weakest in the vertebrate subtree
(17.4/54). Re-rooting the vertebrate subtree at every one of its 112 edges
leaves the placement where it is: the minimum-event rooting needs 64 events
against 67 for the outgroup rooting, and both put the deepest duplication at
Vertebrata ({fig:reconciliation}a).

**What forces the older placement is where one cyclostome lineage attaches.**
One of the two cyclostome clades sits among the gnathostome paralogues, sister
to *ITPR2* and *ITPR3* at maximal support. A lineage holding both cyclostome
and gnathostome copies forces every duplication above it to predate the split
between the two groups. The same attachment keeps the later duplication on the
gnathostome stem, because the node uniting *ITPR2* and *ITPR3* then holds
gnathostomes only. Dropping the six cyclostome tips moves the deepest
duplication to Gnathostomata, so the result rests on six tips and we measured
the standard objection to them. Cyclostomes are the classic long-branch risk,
and long branches are drawn to the root, where this answer sits. Their
root-to-tip distances lie at 0.957 to 1.045 times the median of all 57
vertebrate tips, and none is an outlier ({fig:reconciliation}b).

**The implied losses count sampling, not biology.** The
reconciliation implies 47 lost cells. Checked one at a time against the genome
sweep, 26 are the paralogue present in the genome and absent only from the
134-tip sample, and 0 are corroborated ({fig:reconciliation}c). We report no
losses from it; a loss count needs the genomes themselves.

### The teleost duplication doubled *ITPR1* alone, in one ancestral event

**Five predictions of the third genome duplication hold.** Above the
contiguity bar, every non-teleost gnathostome genome carries one copy of each
paralogue and three ryanodine receptors. The four pre-duplication ray-finned
genomes (bichir, gar and bowfin) carry one of each and three ryanodine
receptors. The 73 teleost genomes carry a mean of 1.973 copies of *ITPR1*,
1.041 of *ITPR2* and 1.041 of *ITPR3*, with 97.3 % carrying two copies of
*ITPR1*, while the same genomes carry 5.822 ryanodine receptors
({fig:duplication}a). The lineages with a further duplication, sturgeon and
salmon, carry three, two and two copies and eight ryanodine receptors
[R185]. The ryanodine count is what makes the singletons interpretable: the
duplication copied *ITPR2* and *ITPR3* as surely as *ITPR1*, and the sister
family kept its copies in the same genomes, so the missing second copies are
a fact about retention, not about the search.

**The two copies are dispersed, not tandem.** Of 71 teleost genomes carrying
two copies, 20 place them on different contigs, and the rest a median
8,870,578 bp apart, the closest pair anywhere at 535,374 bp.

**The two copies partition one ancestral block.** In 49 two-copy genomes with
both copies flanked, the two copies' ancestral flanking symbols are disjoint
in 45 against a tetrapod reference and in 46 against a pre-duplication
ray-finned reference ({fig:teleost}). Disjoint is the load-bearing word. Two
copies that each resemble the ancestor could be two later duplications; two
copies that divide the ancestral neighbourhood between them are what
reciprocal loss after a single duplication produces.

**One duplication made both copies.** Each genome's copy
labels are arbitrary, so we let anchor genomes from 6 orders assign every
genome's two copies independently and asked whether the anchors agree. Under
independent lineage-specific duplications they would agree half the time, and
a constructed control requires the statistic to land there. All 705 of 705
assignments agree (p = 1.2 × 10⁻²¹²) ({fig:duplication}b). A sequence call
made with no neighbourhood input agrees too: in the six genomes whose two
copies won different baits, the bait and the block agree six times out of six.

### The three paralogues share one intron set, and the ryanodine receptors share none of it

**The exon structure is read from the genome and checked against it.** We
measured the exon structure of every locus in the 189 assemblies contiguous
enough to carry the whole gene, 1,378 genes and 112,254 junctions, from the
same spliced alignments the sweep produced. What counts as an intron was set
by the genome rather than chosen: every block gap was scored by its splice
dinucleotides, and since no population of short unspliceable gaps exists, the
floor went at the smallest gap the genome calls spliceable. 99.89 % of
junctions read as a spliceable pair, with 1,035 at minor splice sites, and
177,710 of 188,146 coding edges annotated by independent pipelines in 164
genomes (94.5 %) fall exactly on a boundary the alignment placed
({fig:architecture}b).

**The exon count is conserved and the genomic span is not.** The three
paralogues carry a median 58, 57 and 58 coding exons around 8,250, 8,100 and
8,008 bp of coding sequence, at mean exon lengths of 141.4, 142.2 and
137.8 bp. Their median genomic spans are 147,445, 243,700 and 58,175 bp, a
4.19-fold range at a 3 % difference in coding length ({fig:architecture}a).
The ordering holds inside single genomes: *ITPR3* is the shorter gene than
*ITPR1* in 161 of 182 genomes and than *ITPR2* in 146 of 181. The ryanodine
receptors, measured through the same instrument, are 104-exon genes around
14,910 bp of coding sequence, at a mean exon length of 144.5 bp.

**The intron positions are shared.** We placed every intron in one
coordinate frame, the representative alignment's columns, reached through
each bait's human reference. The frame was checked before use: all 14
residues measured on the IP₃-bound structure (the ten ligand contacts and the
filter and gate residues) land in one column in all three paralogues. An
intron's position is the column its upstream exon ends in together with its
phase, since two introns in the same place in different frames are not one
ancestral intron. Scored inside each genome against an exact null that places
each intron independently and uniformly on the columns both loci resolve,
*ITPR1* and *ITPR2* share a median 48 positions, *ITPR1* and *ITPR3* 46, and
*ITPR2* and *ITPR3* 49, where 0.55 are expected. These are 87.6-, 85.6- and
88.4-fold enrichments, significant in every one of the 183, 182 and 181
genomes tested ({fig:introns}). Each paralogue carries a core of 48 to 49
positions present in over 90 % of its genomes.

**The sister family shares a fold and none of the introns.** Measured through
the same instrument in the same genomes, an IP₃ receptor and a ryanodine
receptor share a median of one intron position, significant in 0 of 188, 0 of
184 and 0 of 183 genomes ({fig:introns}b). The ryanodine receptors have a core
of their own, 91 positions present in over 90 % of their genomes, so the
instrument reads their structure well. The two families share four Pfam
signatures and a pore, and whatever their shared domain architecture means, it
was not inherited as a gene structure.

**The same measurement shows where fragmentary annotations stop.** The
assemblies' own annotations hold 291 of these loci as split or fragmentary
gene models. Of 1,448 internal termini of those models, 689 fall on an exon
boundary the gene carries, 406 fall inside an exon and 353 inside an intron,
and 228 of the 291 loci carry a terminus that no splice site explains
({fig:architecture}c). Only 3 are broken entirely at real junctions. A
duplication detector built on the same blocks, scored against the committed copy
count it never sees, finds 342 of 343 duplicated cells and wrongly calls 21 of
893 single-copy ones.
