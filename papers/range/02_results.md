## Results

### Three independent instruments separate the IP₃ receptors from the ryanodine receptors

Before any search ran, the discovery scorer was benchmarked on 25 positive
controls, the three vertebrate paralogues across a species panel plus the
invertebrate and non-metazoan single-receptor grade, and 31 decoys chosen to
be hard: six ryanodine receptors, which carry all four family signatures, the
MIR-sharing O-mannosyltransferases POMT1 and POMT2, fifteen in-band channels
and non-channels, and out-of-band giants. The first run called all six
ryanodine receptors IP₃ receptors. The fix, used unchanged for the rest of the
work, is a labelled-bait identity margin: a candidate is assigned to whichever
family's labelled baits it is closer to, and only when it wins by more than
0.10. The margin never reads a gene name, and length never decides a call.
With it, recall is 24 of 25 and specificity 31 of 31, and all 6 ryanodine
receptors are rejected. Applied to every record in the census, the margin
leaves an empty gap between the two families ({fig:margin}). The one missed positive, *Drosophila* Itpr, fell 0.008
short of the identity its novelty component needs, which measures the
instrument's worst case: a lone family member with no close sibling in the
set.

The search space was then enumerated to exhaustion rather than sampled. Three
family signatures (PF08709, PF01365 and PF08454) were walked through InterPro
[R155, R156] to cursor exhaustion, because the advertised counts are wrong in
both directions: PF08709 advertises 12,339 records and serves 12,507, while
PF08454 serves 11,890, fewer than it advertises. The union is 15,421 proteins
across 1,488 taxa, and only 10,256 carry all three seeds. A census built on
PF08709 alone, the signature that names the family, would have missed 2,914
proteins, 758 of them IP₃ receptors across 385 taxa ({fig:enumeration}a,
b). A search of the public databases by name, the obvious alternative, had
returned almost none of the records outside the vertebrates
({fig:census_growth}a).

Every enumerated record was called by domain architecture. A record is a
ryanodine receptor if it carries any of four RyR-specific signatures and an
IP₃ receptor if it carries the complete five-signature IP₃R architecture and
none of them. That rule calls 6,433 IP₃ receptors, 6,807 ryanodine receptors
and leaves 2,181 records uncalled. Because the rule never sees a gene symbol,
symbols are an independent label to audit it against: it agrees with the
symbol on 2,479 of 2,479 decided IP₃ receptor records and 2,492 of 2,492
decided ryanodine receptor records. Length separates the two called families
cleanly, but the uncalled records are short fragments where length says
nothing, which is why length never makes the call ({fig:census_growth}b).

A rule that reads annotation cannot call a partial record, so a second
instrument reads residues. Two profile hidden Markov models [R153], one built
from 34 IP₃ receptor seeds (2,684 match states) and one from 22 ryanodine
receptor seeds, assign each target to the profile that beats the other by more
than 10 % of its own score, provided the winner clears 30 bits over at least
200 match states. Calibrated against the architecture call before use, with
the seeds removed, the profiles agree on 11,875 records and disagree on 1
({fig:separation}). They also call 2,074 records the architecture rule left
unassigned as IP₃ receptors and 240 as ryanodine receptors. Merged by a
stated rule, census v3 holds 16,039 records: 8,000 IP₃ receptors, 7,432
ryanodine receptors and 605 unassigned ({fig:instruments}).

### The profile search misses no receptor inside the proteomes it searched

Swept over 763 vertebrate reference proteomes, 14,414,821 proteins, the
profiles declined 12,609 hits as module-only matches, proteins that share one
small domain with a channel model, and found 618 targets the InterPro walk had
never returned, 188 of them IP₃ receptors. The completeness check was run in
the expensive direction. Of the enumerated IP₃ receptor records in swept taxa
that the profiles did not report, all 2,787 were looked up in the search
database itself: none was in the database and missed. Every one is a UniProtKB
entry that its taxon's reference proteome does not contain. Fifteen
vertebrate proteomes carry no family record at all, and copy number per
proteome rises with the size of the gene set ({fig:instruments}b), which is
why every absence claim below is taken to an assembly.

Iterative search gives a third line of evidence. Three jackhmmer runs [R154]
seeded from a human, a fly and an amoebozoan receptor intersect on 4,785
family records and differ by at most 162; what the seed changes is only the
off-family debris ({fig:completeness}). Iteration returned no record the
profile pair calls family that one pass had not already returned. The coded
kill criterion written for this family's hazard, a rise in the ryanodine
receptor share of a round, fired on none of the three runs that drifted,
because accreting unrelated proteins dilutes that share rather than raising
it. Measured on the finished models, the seven runs separate cleanly, with the
drifted ones at an off-family share of 0.81 or more and the rest at 0.34 or
less. The round ceiling caught all three drifted runs, so no drifted round
entered any count, and moving the same rule onto the off-family share
classifies all seven runs correctly (3 of 3 drifted, 4 of 4 clean).

### The family is present in 662 of 6,928 non-vertebrate proteomes and in no prokaryote

The same instrument, unchanged, was swept over 6,928 further reference
proteomes: 63,144,898 proteins in four eukaryotic groups that partition
Eukaryota with the vertebrate set, all 634 archaeal proteomes, and 3,537
bacterial proteomes sampled as the largest proteome per genus. The family
separation transfers intact: of 28,137 targets scored by at least one profile,
2,769 were scored by both above the floor, and 1 of those falls inside the
no-call band ({fig:separation_outside}).

**Presence spans the animals and several protist and algal lineages.** 662 of
the 6,928 proteomes carry an IP₃ receptor call, and 45 of the 135 clades swept
do ({fig:range}). It is in 292 of 311 arthropod, 106 of 110 nematode and all
molluscan, cnidarian and sponge proteomes, and outside the animals in
Amoebozoa, Discoba, ciliates, oomycetes, dinoflagellates, haptophytes and
green algae. It is in no archaeal or bacterial proteome. The choanoflagellate
*Salpingoeca rosetta* carries a full-length ryanodine receptor that is not a
profile seed, so the split between the two families predates Metazoa, as
comparative surveys suggested [R35, R37]. With the gene models from the genome
sweep added, the census holds 18,065 records, of which 8,990 are IP₃ receptors
across 1,402 taxa.

### Land plants and Dikarya have lost the receptor that their early-diverging relatives retain

**The two absences sit beside presences in the same kingdom.** None of 384
Streptophyta (land plant) reference proteomes carries a call, against 15 of 48
Chlorophyta. None of 1,034 Ascomycota or 319 Basidiomycota proteomes carries
one, against 18 of 34 Mucoromycota and 6 of 16 Chytridiomycota. The family is
also absent from all 60 apicomplexan, 29 microsporidian, 27 glomeromycote, 18
mortierellomycote and 35 kickxellomycote proteomes, and from 16 diatom and 12
red algal ones.

**Each negative claim has a positive control inside the same search.** Every
negative was made twice, at E ≤ 1e-5 with the two full-length profiles and at
E ≤ 10 with those profiles plus the family's four Pfam domain models, because
a short domain model can reach a sequence a 2,684-state channel model cannot.
At the relaxed sensitivity, PF08709 returns no substantial match in land
plants against 26 in Chlorophyta, and none in Dikarya against 16 in
Mucoromycota. PF02815, the MIR domain that every eukaryote carries on
unrelated proteins, returns 633 substantial matches in land plants and 4,376
in Dikarya. The instrument finds thousands of the promiscuous domain in the
proteomes where it finds none of the family-defining one. Iterating a model
seeded inside each group added 9 targets the single pass never reported, and
the two that sit in an empty lineage are O-mannosyltransferases, the MIR
decoy, rather than receptors.

### Every plant and fungal record is a real gene or a fragment of one

A handful of records in kingdoms whose model organisms have none is the shape
of a contamination artefact, so all 99 plant and fungal IP₃ receptor records,
64 plant and 35 fungal, were chased individually through four lines of
evidence: sequence, the record's own provenance, the nearest relative outside
its kingdom and the nearest within it. 47 are real genes and 52 are fragments
of real genes; none is a contamination suspect and none lacks genome backing
({fig:chase}). The contamination test had full power. A genuine deep homologue
is 20–40 % identical to its metazoan relatives and an assembly contaminant is
95–100 % identical to one particular animal. The 22 surviving plant records
sit at 19.9–39.9 % identity to anything outside their kingdom and the 25
fungal ones at 20.5–33.5 %, and no record of the 99 exceeds 45.8 %.

### Thirty-five clade-level absences hold in controlled genome assemblies

A proteome absence is a fact about a gene caller, so we re-asked each one of a
genome. We selected 194 non-vertebrate assemblies, 100.4 Gbp, from the 24,596
NCBI eukaryotic reference assemblies [R159] less the 6,216 vertebrate ones,
by five rules that sample most finely where the negative claims are: one
assembly per phylum, one per class in the most-sampled phyla, one per clade
that the proteome sweep called empty, anchor organisms with a known answer,
and the highest-copy species per class. The third rule re-derives the target
list from the presence table and found three absences no summary had named:
diatoms (0 of 16), red algae (0 of 12) and Cestoda (0 of 11), a metazoan
clade.

**The positive control was chosen per clade by measurement.** The vertebrate
convention of using the ryanodine receptors as proof of search does not
transfer, because architecture-level ryanodine receptors occur in 2 of the
6,928 non-vertebrate proteomes and a silent RyR in a land plant is the correct
answer. A single fixed replacement also failed: both apicomplexan classes carry
a single MIR protein across their 60 proteomes. Six candidate families, each
large, deeply conserved and multi-exon, were therefore run over every
clade's own proteomes, and each clade takes the one that is actually present.
Apicomplexa take the myosin head domain, present in 36 of 36 Aconoidasida and
23 of 23 Conoidasida proteomes; red algae take the SMC N-terminal domain,
present in 12 of 12. The panel carries 108 baits, of which 37 are IP₃
receptors, 16 are ryanodine receptors and 55 are clade controls, and every
one of the 56 control clades is filled.

**No genome in the sweep is uncontrolled.** In 115 genomes the control is
recovered and baits from another kingdom-level group align across the same
locus, showing the search crosses the divergence a receptor would have to be
found across; in 77 the receptor itself was found; in two ciliate genomes the
control is partial or absent, and no absence claim rests on either. Against
that, **all 35 clade-level absences hold at assembly level** and none fails
({fig:absence}). Ascomycota has no receptor in 31 controlled assemblies,
Streptophyta none in 25 and Basidiomycota none in 17. Two land-plant
assemblies carry only scattered tblastn traces, short of a gene. The genomic
thresholds were measured for this scope rather than transferred: sixteen
annotated gene spans run from 3,739 to 324,840 bp with a median of 19,435, so
the contiguity bar is set per group, and a cluster of alignments is not
counted as a gene until it passes rules for split and chained models
({fig:genome_instrument}).

### Copy number outside the vertebrates runs from zero to eighteen

Outside the vertebrates there are no paralogue cells to fill, so the question
is how many receptor genes a genome carries. Of the 194 genomes, 117 carry
none, 43 carry one, and the rest carry two or more ({fig:copies}). The flatworm
*Macrostomum lignano* carries 18 complete gene models, the ciliate *Stentor
coeruleus* 13 and the sponge *Dysidea avara* 8. The vertebrate count of three
is not the family's maximum. The genome also corrects the protein side: the
green alga *Cymbomonas tetramitiformis*, which contributes twelve of the 22
surviving plant records, carries three complete gene models.

### The fold separates the two families where a structure exists

A structural check is independent of every instrument above, because it reads
no gene symbol, domain annotation or alignment score. It is also limited by
what exists. AlphaFold DB [R167, R168] holds a usable model for 20.9 % of the
census's protein records, for 9 of the 134 representatives the downstream
analyses use, and for 13 of the 5,861 full-length records (0.2 %). It answers
an accession with whatever record it holds, which for human ITPR2 is a
181-residue isoform ({fig:afdb}). The structural test therefore rests
mainly on cryo-electron microscopy depositions, chosen by seven rules with the
family call taken from this census rather than from entry titles.

The comparison scale was calibrated on this panel before anything was called
({fig:fold}). None of 81 pairs involving a negative control reaches TM-align's
same-fold bar of 0.5 [R166]. Two structures of the same paralogue in
different conformations score a median of 0.78 against 0.43 for two different
IP₃ receptors, so conformation alone moves the score by more than the
difference any fold claim here relies on. Scored against both reference
families with a 10 % relative margin, all 20 structures the test could call
are called as the census calls them, and none of the 3 negative controls is
called at all. The test does not reach the deep lineages: 2 of 7
non-vertebrate models clear the same-fold bar, and the rest prefer the IP₃
receptors by about two to one while falling short on similarity, not on
length. Model confidence is highest in the IP₃-binding core (median pLDDT
83.9) and lowest in the pore among named domains (71.0), against 69.5 outside
annotated domains ({fig:structures}b).
