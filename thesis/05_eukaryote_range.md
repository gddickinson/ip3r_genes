# 5. How far the family reaches across the eukaryotes

## 5.1 Two questions that look like one, and have to be answered separately

"Where does this gene family occur?" is normally answered from a database and
reported as a range. That answer conflates two questions that have to be kept
apart, and this chapter answers them separately.

The first question is about the record: in how many reference proteomes does a
family member appear? That is a statement about what gene callers have found,
and it is the one a database can support on its own.

The second question is about the DNA: in genomes where no gene caller found
one, is there one? Answering it requires searching assemblies and, crucially,
a **positive control that demonstrates the search could have found the gene
had it been there**. Without that control, "no IP₃ receptor in *Arabidopsis*"
is indistinguishable from "the search did not run properly on *Arabidopsis*".

Both questions are answered here, over 6,928 reference proteomes and 194
genomes.

## 5.2 The proteome sweep covers 6,928 proteomes and partitions the eukaryotes

The sweep covers six groups: non-vertebrate metazoa, fungi, green plants,
other protists, archaea and a genus-stratified bacterial sample. Together they
comprise **6,928 reference proteomes, 63,144,898 proteins and 24.93 gigabases
of residue.** The four eukaryotic groups partition Eukaryota with Chapter 3's
vertebrate sweep, so no proteome is swept twice and the two denominators add.

Two sampling decisions are stated rather than assumed. Bacteria are sampled at
one proteome per genus, taking the one with the most proteins, which is the
most sensitive member of each genus, so a negative claim is made in the places
most likely to break it. Archaea's reference set is small enough to take
whole, so no sampling caveat attaches to it.

Twenty-three proteomes that UniProt lists are not published in the release
tree and return a permanent 404. They are excluded from the denominator and
recorded, together with what they took out of it, rather than retried. A
listed but unpublished proteome is not a transient failure and retrying cannot
fix it.

The instrument is Chapter 3's, unchanged: the same two profiles, the same
30-bit floor, the same 200-match-state span gate and the same relative margin.
Reusing it rather than building a new one is what makes the vertebrate and
non-vertebrate numbers comparable at all.

## 5.3 The family is present in 662 of 6,928 proteomes, and in no prokaryote

**662 of 6,928 reference proteomes carry an IP₃ receptor call**, and 45 of the
135 clades swept do.

![](figures/range_by_phylum.png)

**{fig:range_by_phylum}.** The family across the eukaryotes, drawn as a
fraction of swept proteomes rather than of records. That denominator is the
whole point: counted by records, a single well-sequenced alga outvotes a
sparsely sampled phylum, and the figure would report sequencing effort.

The shape is a family that is ancestrally eukaryotic and has been lost
repeatedly. It is in 94 % of arthropod proteomes, 96 % of nematode, and 100 %
of molluscan, cnidarian and sponge. It is in 65 % of oomycetes, 67 % of
euglenozoans, 94 % of ciliates and 69 % of amoebozoans. It is in none of 634
archaeal and none of 3,537 bacterial proteomes.

![](figures/profile_separation_euk.png)

**{fig:profile_separation_euk}.** The same two-profile separation outside the
vertebrates, with the no-call band drawn. 28,137 targets were scored by at
least one profile and 2,769 by both above the floor, which are the only ones
where the two families can be said to compete at all. Of those, one falls
inside the band where the instrument declines to choose. The family separation
built for the vertebrates transfers to the eukaryotes without modification.

## 5.4 Land plants and Dikarya have no IP₃ receptor, and their early-diverging relatives do

Chapter 2 recorded an anomaly: forty plant and forty-one fungal proteins carry
the IP₃-binding-core signature while *Arabidopsis* and budding yeast have
none. Chapter 3 sharpened it. In the enumerated space every green-plant call
was in Chlorophyta and **Streptophyta, the land-plant lineage, had none**,
while every fungal call was in an early-diverging phylum and **Dikarya
contributed no records at all**.

Those were statements about a space a protein enters only by already carrying
a family domain annotation. Sweeping the proteomes themselves asks the
question without that filter, and both statements hold.

For land plants, **0 of 384 Streptophyta reference proteomes carry a call**,
against 15 of 48 Chlorophyta. For the fungi, **0 of 1,034 Ascomycota and 0 of
319 Basidiomycota** carry one, against 18 of 34 Mucoromycota, 6 of 16
Chytridiomycota and both Basidiobolomycota.

That is the shape of a loss rather than of an absence: the family is present
in the early-diverging lineage of both kingdoms and gone from the derived one.

**Every plant and fungal record was then chased individually**, 64 and 35 of
them, through four independent lines of evidence, because a handful of records
in a kingdom whose model organisms have none is exactly the shape of a
contamination artefact.

![](figures/plant_fungal_chase.png)

**{fig:plant_fungal_chase}.** Every plant and fungal record chased, with the
cross-kingdom identities the contamination test rests on. That is the
load-bearing axis: a genuine deep homologue is 20 to 40 % identical to its
metazoan relatives, while an assembly contaminant is 95 to 100 % identical to
one particular animal. Nothing here is above 40 %.

Twenty-two plant records and 25 fungal records survive every test, and the
rest are fragments of real genes. The surviving plant records run from 19.9 %
to 39.9 % identity to anything outside their kingdom, and the fungal ones from
20.5 % to 33.5 %. None is a contaminant.

They are also not evenly spread. *Cymbomonas tetramitiformis*, a green alga,
contributes twelve of the twenty-two plant records against a median of one per
species. Whether that is a real copy-number expansion or a duplicated assembly
is not a question a proteome can answer, and §5.8 answers it.

## 5.5 Every negative claim was made twice, at two stated sensitivities

Every negative claim in this chapter was made twice: once with the two
full-length family profiles at a strict threshold, and once at a much looser
one with those profiles plus the family's four individual Pfam domain models.
A 213-state domain model can reach something a 2,684-state channel model
cannot, and a negative claim should be made with the most sensitive instrument
available rather than the most specific.

The design behind that is worth stating. The same search is run once and
filtered twice: every search runs at the loose threshold and the strict call
is taken by filtering the same output, so the strict set is exactly what a
strict run would have produced and the loose set is a superset from the same
search rather than a second experiment.

That equivalence was tested rather than assumed, and **the assumption behind
it turned out to be wrong.** The reporting threshold is applied to the
conditional expectation value, which is normalised by how many sequences
passed, so a looser setting inflates every conditional value by a constant
factor and pushes marginal domain rows out of the report. Measured across 798
protist targets, the sequence sets and all 798 assignments agree and no gate
decision moves, but 11 domain rows differ. A second effect appears only at
scale: the sequence expectation value is printed to two significant figures,
so a target whose true value sits just above the cut prints as if it were at
the cut, and 42 of 13,770 plant targets on the ryanodine profile are admitted
by a less-than-or-equal filter that a strict run would never have reported.
All 42 are declined by the span gate anyway.

The pass criterion is therefore about the call rather than about set equality:
no assignment moves, no gate decision moves, and no one-sided target is
called. The equivalence test also refuses to pass vacuously. Its first run
went green on the archaeal group, which has no hit at either threshold, so a
comparison of fewer than 25 targets is now a failure with its own message.

**Each group was also searched again iteratively**, from a seed native to that
group, under the same coded kill criterion Chapter 3 used. The completeness
question this answers is precise. A model seeded in one lineage and iterated
to convergence either reaches a neighbouring lineage or it does not, and the
targets that only iteration found are the ones that matter. If they sit in the
same lineage as the seed, the absence next door is not a sensitivity artefact.

![](figures/jackhmmer_s20.png)

**{fig:jackhmmer_s20}.** Per-group convergence with the sister-family trace
beside it. One of four groups converges and the other three hit the kill
criterion, and Chapter 4 explains why the rule that catches them is the wrong
one.

Nine targets entered an accepted model that the single pass never reported.
They are named rather than counted, and two of them sit in a lineage this
chapter calls empty, so the absence there is not quite absolute and what those
two are decides whether that matters. **Both are mannosyltransferases**,
meaning the MIR-domain sharer that Chapter 3's decoy panel was built around
rather than a receptor. The iterated search reaches these lineages exactly far
enough to pick up the known false positive and no further, which is the most
informative result the iteration could have produced.

## 5.6 Taking the absences from the proteome to the genome

A proteome absence is a fact about a gene caller. The 194-genome sweep turns
the ones that matter into facts about genomes.

**The scope is chosen to sample most finely where the negative claim is.**
Five rules were applied to 24,596 NCBI eukaryotic reference assemblies, minus
the 6,216 vertebrate ones: one assembly per non-vertebrate eukaryotic phylum;
one per class in the five phyla with more than 100 swept proteomes; one per
clade with at least ten swept proteomes and zero calls; the anchor organisms
whose answer is known from the literature; and the highest-copy species per
class. That gives 194 genomes and 100.4 Gbp.

The third rule re-derives the target list from the presence table rather than
repeating a summary, and in doing so it **found three absences that summary
never named**: diatoms at 0 of 16, red algae at 0 of 12, and **Cestoda at 0 of
11, a metazoan clade with no record at all.**

**Both genomic thresholds had to be re-measured, and neither transferred.** A
vertebrate IP₃ receptor gene spans 76 to 498 kb and a *Drosophila* one spans
22 kb. Measured across sixteen genes in sixteen species from NCBI's own gene
annotations, and deliberately not from the aligner, whose own intron setting
shapes the loci it reports so that a measurement taken that way cannot falsify
the setting it calibrates, spans run from 3,739 to 324,840 bp with a median of
19,435. That is a hundred-fold range against the 6.5-fold measured inside the
vertebrates, so **the contiguity bar is set per group**. Metazoan genes are
about nine times longer than protist and fungal ones, at 83 kb against 7 to 9
kb, and the genomes carrying this chapter's negative claims are judged against
the smaller bar, which almost any modern assembly clears.

One gap is stated in the output rather than hidden. No green-plant gene span
could be measured, because the chlorophyte census records carry locus tags
with no gene record, so every plant genome is judged against the global
median. The calibration function returns that fallback in its own return value
so that a caller cannot fail to notice.

## 5.7 The positive control had to be replaced twice

Chapter 4 could write "the ryanodine control fired in all 309 genomes" because
every vertebrate has three ryanodine receptors. Outside the metazoa that
sentence is not available, because the proteome sweep found architecture-level
ryanodine receptors in **2 of 6,928** non-metazoan proteomes. A land plant
with no ryanodine locus is the correct answer, so it cannot be that genome's
proof of search.

The first replacement was the MIR-domain sharer, meaning the
O-mannosyltransferase family, which is in the family's own signature set,
present in the clade each claim is about, and simultaneously the sharpest
available decoy. Where a query on the IP₃-core signature returns nothing in
land plants, a query on MIR returns 633.

**It was still one profile fixed in advance for every clade, and that is what
failed.** Both apicomplexan classes carry a single MIR protein each across 60
swept proteomes, and the pilot's *Toxoplasma gondii* came back with neither a
receptor nor a control, making it the one genome in the pilot whose absence
claim therefore rested on nothing.

The fix is not a different profile. It is not fixing the profile in advance. A
control's job is to fire in the clade whose absence is the claim, so which
family makes the best control is a property of the clade and is measurable.
Six candidate profiles, all large, deeply conserved, multi-exon eukaryotic
families, were run over each clade's own swept proteomes, and each clade takes
the one that is actually there, with every candidate's coverage committed
rather than only the winner's. The MIR bait stays in the panel whatever the
measurement says, because dropping it would buy a proof of search and sell the
family's negative control.

**193 of 194 genomes are controlled**, and 115 of them across a kingdom
boundary, meaning that a control bait from a different kingdom-level group
aligns across the same locus. That is the tier an absence claim needs, because
it shows the search reaches across the divergence any receptor here would have
to be found across rather than merely that the assembly is readable. In the
pilot, 1 of 14 genomes was uncontrolled. Here, none.

![](figures/absence_at_genome.png)

**{fig:absence_at_genome}.** Each absence claim with the number of controlled
genomes drawn beside the number searched. A claim whose control bar is short
is standing on nothing, and drawing the two together is the only honest way to
present a table of zeros.

**Every one of the 35 clade-level absences holds at assembly level, and none
fails.** Ascomycota has none in 31 controlled assemblies, Streptophyta none in
25, and Basidiomycota none in 17. Apicomplexa, Microsporidia, Glomeromycota,
diatoms, red algae and Cestoda likewise.

## 5.8 Copy number outside the vertebrates reaches eighteen

Outside the vertebrates the three paralogue cells do not exist, because ITPR1,
ITPR2 and ITPR3 are a vertebrate whole-genome-duplication product and Chapter
4 found the trio absent below the cyclostomes. The question is simply how many
receptors a genome has.

Of 193 controlled genomes, 116 carry none, 43 carry one, and 34 carry more
than one. The top of the distribution is a flatworm, *Macrostomum lignano*,
with 18 complete gene models, then a ciliate, *Stentor coeruleus*, with 13,
then two sponges at 8 and 6.

![](figures/nonvert_copy_number.png)

**{fig:nonvert_copy_number}.** Copy number outside the vertebrates, with the
vertebrate paralogue count drawn as a reference line rather than as a
category. The family's copy number is not a vertebrate story: several
invertebrate and protist lineages carry more receptors than any vertebrate.

*Cymbomonas*, the green alga that contributed twelve of the twenty-two
surviving plant protein records, carries **three** complete gene models across
two clusters. Twelve records and three genes is what a duplicated assembly and
a multi-isoform gene set look like from the protein side, and it is the answer
the proteome could not give.

## 5.9 Two further thresholds had to be measured for this scope

Two thresholds in this sweep are measurements rather than settings, and both
had to be made here because neither transfers from the vertebrates.

**A cluster of alignments is not a gene.** With the intron parameter set high,
as §5.6 requires, a cluster is a large object, and two rules pull in opposite
directions. Neighbouring same-strand clusters whose bait spans are
complementary are one gene the aligner broke, and one cluster containing two
distinct complete models is two genes. Both were implemented, both record the
numbers they fired on, and their measured effect over the sweep is seven genes
recovered in seven genomes.

**The identity floor below which a locus is not a locus** was measured against
loci whose identity the assembly's own annotation establishes, which is
evidence the alignment score did not produce. The sweep deliberately records
every cluster down to a much lower floor precisely so that the floor is not
measured from the population it has already filtered.

![](figures/identity_floor.png)

**{fig:identity_floor}.** The measured identity floor and the two populations
it separates. The negative result here is what changed the instrument: outside
the vertebrates the annotation axis is nearly empty, with 21 of 917 clusters
carrying an informative gene name, so a second axis was added in which every
recorded cluster is scored against both family profiles. Neither identity nor
coverage separates the confirmed and contradicted populations cleanly.

That calibration also refuses to write below a floor of genomes and confirmed
loci, and the reason is an incident. A smoke-test run over one genome produced
a threshold of 0.25 from two loci, wrote it, and the next sweep read it back
and moved three genomes from one status to another. A calibration that can
pass vacuously is worse than no calibration.

![](figures/span_inflation.png)

**{fig:span_inflation}.** Locus span against coding footprint, by group. A
locus is much larger than the gene inside it, and by a factor that varies by
group, which is the reason the intron parameter is set per group and the
reason a locus is not a copy.

## 5.10 The family's range, stated with its scope attached

The family is **ancestrally eukaryotic**. It is present across the metazoa,
across several protist lineages including ciliates, oomycetes, euglenozoans
and amoebozoans, in the green algae, and in the early-diverging fungi. It is
absent from land plants, from Dikarya, and from every archaeal and bacterial
proteome swept.

Every one of those absences is a claim about a declared space. The
proteome-level claims are about 6,928 reference proteomes, and the
genome-level claims are about 194 assemblies in which a measured positive
control recovered a comparably long, deeply conserved gene. Where the two
disagree, which they do not in any of the 35 clades tested, the genome would
win.

The merged census is 18,065 records, of which 8,990 are called IP₃ receptor
across 1,402 taxa. Chapter 12 returns to the absences with a different
question: whether the enzyme that makes the ligand went with the receptor.
