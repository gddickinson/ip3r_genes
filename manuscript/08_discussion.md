## Discussion

**An ancient family that vertebrates cannot do without and other lineages
have dropped repeatedly.** The IP₃ receptor is present across the eukaryotic
tree, in animals, amoebae, discobids, ciliates, oomycetes, dinoflagellates,
haptophytes and green algae, and its duplication from the ryanodine
receptors predates the animals: two unicellular holozoans carry both
families at full length. Set against that background, the absences are the
informative part, and they are not one event. Land plants have none, and
their sister green algae do. The Dikarya have none, and neither do four
other fungal phyla that branch apart from them, so at least five independent
fungal losses are implied. Apicomplexa, Microsporidia, Rhodophyta, diatoms
and one tapeworm class complete a list of lineages that lost a channel
otherwise kept for a billion years. These are eukaryotes with reduced or
re-routed calcium signalling — obligate parasites, walled autotrophs, yeasts
— a pattern consistent with wider surveys of the eukaryotic calcium toolkit
[R107, R108, R109], which makes the correlation worth naming and not worth
over-reading from a presence table. What can be said, and what previous surveys could not say,
is that each absence is a claim about a declared space in which a comparably
long, deeply conserved gene *was* found by the same instrument in the same
assembly.

**Complete retention within the vertebrates is the sharper result.** In 309
genomes spanning every vertebrate order, no paralogue is lost: not once in
927 genome × paralogue cells, and in every assembly contiguous enough to
carry the gene the genome holds at least three gene-equivalents. Because a
count of zero is easy to obtain by not looking hard enough, the design was
inverted — the search's own sensitivity was measured against an independent
replicate (the ryanodine receptor cells, 13.6 % missed against 15.2 %), the
loss state was shown to be reachable by construction, and instead of a
robustness check we report the settings that manufacture a loss: the
family-level count breaks in 2 of 32 and only by discarding two decisions
made for reasons unrelated to the count. A gene family that a whole
vertebrate lineage has never dispensed with, against a eukaryotic background
of repeated loss, is a statement about vertebrate physiology rather than
about the protein.

**Retention is not symmetric, and the asymmetry has one direction.**
*ITPR1* is held about twice as tightly as its sisters by three independent
framings (one-ratio ω, whole-tree two-ratio contrasts, and a
distribution-level relaxation test), is the only paralogue whose stem
carries an interpretable episodic signal, is the only one doubled and kept
after the teleost genome duplication, and is the paralogue through which
both surviving two-round paralogon links run. The one decay signal in the
data points the other way, at *ITPR3*, and does not survive its own
contiguity control. The obvious reading — that *ITPR1* carries the
ancestral, least substitutable role and its sisters have been freer to
diverge — is consistent with what is known of their expression and channel
behaviour [R29, R33, R67, R101] and with the phenotypes of the single-gene
knockouts [R114, R115], and this work supplies the evolutionary half of it
rather than proving the functional half.

**The sister question, and a tension worth keeping.** *ITPR2* and *ITPR3*
are sisters on sequence: the alternatives are rejected at p < 2 × 10⁻⁵ and
the grouping survives a model-violation guard. Yet those are the two
neighbourhoods that share no dated ohnologue link, while both surviving
links attach *ITPR1* to one of them. The two instruments are measuring
different things — a tree estimates the order of duplication, a retained
ohnologue records which copies escaped deletion afterwards — so this is a
tension rather than a contradiction, and the honest form of the result is to
report both.

**The two families share a fold and not a gene.** The ryanodine receptors
carry every Pfam signature used to diagnose this family, and that is the
single fact that shaped every search in this paper. Intron positions are a
character no protein alignment produces, and measured on them the two
families have nothing in common: an IP₃ receptor and a ryanodine receptor in
the same genome share a median of one intron position, in none of 183 to 188
genomes at p < 0.05, while the three IP₃ receptor paralogues share 46 to 49
in every genome tested. The shared domain architecture is real and it is
old, but it was not inherited as a gene, and the practical consequence for
anyone searching this superfamily is that domain content is the wrong
evidence to separate them on — which is the conclusion this project reached
by a different route at every earlier stage.

**The part that names the family is not the part selection holds hardest.**
The IP₃-binding core is what distinguishes these receptors from the
ryanodine receptors functionally and what the family's diagnostic Pfam
signature is named for, so it is the natural place to expect the tightest
constraint. Paired inside single orthologues against about 250 sequences per
paralogue it is the *less* constrained of the two functional modules, and
the finding that matters more than the direction is that the direction
reverses on where the pore is taken to stop: fifty residues of luminal loop,
which InterPro includes in PF00520, are enough to flip all three paralogues.
A comparison between two protein modules is only as good as the boundaries
drawn around them, and those boundaries are usually inherited from a
database rather than declared. Within the site, what selection holds is a
pocket about fifteen ångström across rather than the ten residues that touch
the ligand — visible only because the site was measured as a distance in six
independent structures rather than taken as a contact list from one. And
losing the enzyme that makes IP₃ does not relax that pocket: matched on
divergence the effect is a bounded null, against a positive control the same
test detects almost always. Whatever holds the binding site in these
lineages, the canonical route to its ligand is not required for it.

**A gene the archive cannot find.** Three quarters of the genes demonstrated
here are unreachable by any protein-database search, more than half of the
full-length records that do exist carry no usable gene symbol, and 318 IP₃
receptor genes have no protein record of any kind. The audit that was
supposed to show this family is poorly recorded instead showed, through its
own sister-family control, that it is recorded about as well as its
neighbour, and that what predicts the record is the archive the assembly
went to and how broken that assembly is. Two consequences follow. For this
family specifically, any survey based on protein records — including every
previous survey of its distribution — is working from about a quarter of
what exists. For gene-family surveys generally — including the ones this family's own
distribution has been read from [R35, R36] — the result argues that the
denominator has to be assemblies, that absence claims need a positive
control in the same assembly, and that a sister family searched alongside is
worth more than any threshold.

**Limits.** The vertebrate scope is one genome per order plus the species
the protein record left in doubt, not every vertebrate genome; a paralogue
lost in a single species inside a sampled order would not be seen. 120 of
309 assemblies cannot hold the gene on one contig, and every result is
therefore reported twice, above and below that bar. The absence claims
outside the vertebrates are claims about 6,928 reference proteomes and 194
genomes, sampled most finely where the claims are, and a clade represented
by one genome is a clade whose absence rests on one genome. Synonymous
saturation means every ω here is tree-based and no pairwise distance is
interpretable. Predicted structures cover 0.2 % of the full-length records,
so nothing structural in this paper rests on a model of a whole subunit. And
the deepest constraint layers built here are not the best variant
classifier: the family-wide layer is, which is a result about what
conservation-based variant scoring should use, and a caution about assuming
that more orthologues of the same paralogue is the same thing as more
signal. The module and pocket comparisons run on those same deep
alignments, which are vertebrate by construction, so the ligand core's
position relative to the pore is a vertebrate statement; and everything
within 15 Å of the ligand is also inside the fold that holds it, so no
measurement here separates constraint on ligand binding from constraint on
domain packing. That would take a mutational or binding dataset; sequence
will not do it.

