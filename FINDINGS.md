# FINDINGS.md — what this project has discovered

The running biological story, one dated entry per completed roadmap task.
Findings before methods, plain language, no file paths. Claims that are not
yet confirmed are marked *(pending: which task confirms it)*.

`SESSION_LOG.md` is the chronological record of what was *run*; this file is
the record of what is *known*.

---

## 2026-08-18 — Project set up; nothing measured yet

No analysis has run. What exists is the question and the instrument.

**The subject.** The IP3 receptor is the endoplasmic reticulum's
ligand-gated calcium-release channel — the source of most agonist-evoked
calcium signals in non-muscle cells. Vertebrates have three of them
(ITPR1, ITPR2, ITPR3), each about 2,700 residues, assembled as tetramers,
and each associated with a different human disease: cerebellar ataxia and
Gillespie syndrome for ITPR1, an inability to sweat for ITPR2, and a
peripheral neuropathy for ITPR3.

**What is already visible in the databases** (checked 2026-08-18, and to be
re-derived properly in S2 before any of it is quoted):

- The family looks overwhelmingly like an animal family: 12,149 of the
  12,338 proteins carrying its core domain are metazoan.
- Yet there are 40 plant and 41 fungal records in a family that textbooks
  say plants and fungi lack — while *Arabidopsis* and baker's yeast have
  none at all. Something there is either a real and under-appreciated
  branch of the family, or a set of database errors. Both would be worth
  knowing. *(pending: S2, S20, S23)*
- Zebrafish carries four: itpr1a, itpr1b, itpr2 and itpr3. The first two
  look like a teleost duplicate pair. *(pending: S16)*

**The complication that shapes the whole project.** Ryanodine receptors —
the other big calcium-release channel, twice the size — carry every domain
that marks an IP3 receptor. They will appear in every search this project
runs. That is a nuisance and an opportunity: it means the two families can
be counted by the same instrument and compared directly, and RyR gives the
IP3R tree a proper outgroup instead of a guess.

**What the project intends to find out.** Whether the family really is
absent from plants and fungi; where the three vertebrate paralogs came from
and whether they are the same age as the three ryanodine receptors; whether
any of them has ever been lost; whether the conserved core explains the
three diseases; and how much of what the databases say about this family is
wrong.

---

## 2026-08-18 — S0: the baseline, checked

**Nearly all of what we thought we knew about this family survived
checking. Four things did not.**

**The textbook picture of the IP3 receptor holds up.** It is the ER's
calcium-release channel; IP3 and calcium together open it, neither alone;
the channel is a four-subunit assembly with the ligand-binding sites at the
top and the pore at the bottom, connected by a tail that runs back up
through the molecule. Every one of those statements now has a paper behind
it rather than a memory. One detail matters more than it looks: **all four
subunits must bind IP3 before the channel opens**. That is why a single
faulty subunit can disable a whole channel — the mechanism behind the
dominant disease variants.

**The three genes are not built the same way, and that was a surprise.**
They make near-identical proteins — 2,671 to 2,758 residues, within 3 % of
each other — but the genes themselves differ enormously. ITPR2 is spread
across 498,000 letters of chromosome 12. ITPR3 does the same job in 76,000
letters of chromosome 6. That is a **6.5-fold difference in size for genes
that build the same machine with the same number of parts** (57–62 pieces
each). Something has been adding or removing the non-coding filler between
those parts, in one paralog and not another, and we do not yet know what or
when. *(pending: S21)*

This also corrects the project's own starting document, which said all three
genes span "hundreds of thousands of letters". ITPR3 does not.

**Half of what a standard search returns is the wrong family.** This was the
sharpest result of the session. Asking the databases for every zebrafish
protein carrying the IP3 receptor's own signature domain — the one that
binds IP3, the thing that *defines* the family — returns 109 proteins. **53
of them are ryanodine receptors**, the sister family that shares the domain
but is nearly twice the size. The search that is supposed to find IP3
receptors is 49 % wrong before anyone looks at it. Every count this project
publishes has to survive that, and the rule that separates the two families
has to be tested rather than assumed.

One of those 53 has no name at all — a 4,900-residue protein filed only as
`LOC101884734`. It is the first concrete example of how much of this family
is sitting in the databases unlabelled. *(pending: S18)*

**We nearly assumed our own answer.** The starting document stated that the
IP3 receptors and the ryanodine receptors *independently* grew from one gene
to three in vertebrates. The shared ancestry is real and well supported. The
independence is not established anywhere — and it happens to be one of the
questions this project set out to answer. It has been moved out of the
"known" column and back into the "to find out" column, where it belongs.
*(pending: S7, S13, S16)*

**The disease picture is richer than we recorded.** Gillespie syndrome —
ataxia with a partial iris — turns out to be caused two different ways: by
inheriting two broken copies of ITPR1, *or* by a single new mutation that
poisons the four-subunit channel. Recording only the second would have
skewed the analysis of which parts of the protein matter. And ITPR3, which
we had down as a nerve disease gene, causes something considerably worse in
some people: a multi-system disorder including immune failure. The clinical
panel this family offers is broader, and more informative, than we thought.
*(pending: S17)*

**Nothing changed about the strangest thing in the data.** Forty plant
proteins and forty-one fungal proteins still carry the IP3-binding domain,
in a family that plants and fungi are supposed to lack — while thale cress
and brewer's yeast have none at all. Re-checking reproduced those numbers
exactly. Whether they are real genes, mistakes, or contamination is still
the sharpest question on the list. *(pending: S2, S20)*

---

## 2026-08-18 — S1: the instrument was tested, and it failed the first time

**Our search would have called every ryanodine receptor an IP3 receptor.**
This is the finding of the session, and it is a finding about ourselves.
We knew the two families were confusable — the ryanodine receptors carry all
four of the protein signatures we use to recognise an IP3 receptor. What we
did not know was how the search behaved when actually handed one. We built a
panel of known answers: 25 real IP3 receptors from across the animals, and 31
impostors, six of them ryanodine receptors. Every single ryanodine receptor
came back flagged as a promising new IP3 receptor. Not one was rejected.

The reason is almost funny. The strongest piece of evidence our search uses
is "does this protein carry the family's signature domains?" — and the
ryanodine receptors carry all of them. The one thing that most confidently
identifies a member of this family is the one thing that cannot distinguish
it from its sister.

**What fixed it was a comparison, not a rule.** The obvious patch — throw out
anything called "RYR" — is worthless, because the problem is precisely the
proteins nobody has named yet. (We already have one: an unnamed zebrafish
gene of ryanodine-receptor size that our earlier search returned as an IP3
receptor candidate.) So the test we built ignores names entirely. It asks
each candidate a comparative question: *are you more similar to a known IP3
receptor, or to a known ryanodine receptor?* The answer separates the two
families completely and with room to spare — genuine IP3 receptors lean
toward IP3 receptors by 7 to 87 percentage points, ryanodine receptors lean
the other way by 54 to 87. There is no overlap anywhere in the middle.

**But the test has an edge, and we found where.** *Dictyostelium*, a soil
amoeba, has a true IP3 receptor — and it leans toward the IP3 receptors by
only 6.5 points. It is called correctly, but barely; by our own standard for
what counts as a confident call, it does not qualify. The lesson is specific
and useful: this comparison works throughout the animals and fails at the
base of the tree. Any claim about whether some deep-branching organism has an
IP3 receptor will need a different and more sensitive instrument, not this
one. *(pending: S3, S20)*

**A quieter result about what recognition rests on.** The three human IP3
receptors are 61–68 % identical to one another. That is close enough that
none of our "this looks like an interesting outlier" tests ever fire for a
real family member — they are all too similar to each other to look like
anything new. In practice, every true IP3 receptor in the test was recognised
by its signature domains and by nothing else. So the census's reach is set by
one thing: whether a database record can be matched to a domain signature at
all. A real gene with no such record attached will be missed. That is now a
measured property of the search, not a worry about it. *(pending: S2)*

**One real IP3 receptor was missed.** The fruit fly's, by the narrowest
margin imaginable. The search wants corroboration — a second, similar protein
in another species — and the fly's closest relative in the panel was the
roundworm's, at 34.2 % similarity against a 35 % requirement. It failed by
eight-tenths of a percentage point. Measured a fairer way (comparing only the
regions the two proteins actually share), the same pair scores 42 % and the
fly would have been recognised. This is partly an artefact of a small test
panel — in the real census the invertebrates will corroborate each other —
but the worst case is now known: **a lone, moderately-diverged family member
with no close relative in the dataset can be missed.** *(pending: S6)*

---

## 2026-08-19 — What illustrating the review taught us

Drawing the figures meant measuring things the review had only described, and
three of those measurements changed what it says.

**The gene that defines the family does not find all of the family.** The
IP₃-binding core (PF08709) is the domain this whole project uses to recognise
an IP₃ receptor. *Dictyostelium* iplA — a slime-mould IP₃ receptor that has
been studied for decades, and which sits in our own list of known-good test
cases — does not carry it. Nor does it carry the MIR domain or the pore
domain: three of the family's five signatures are simply absent from its
database record. A census built on that signature would not return it at all.

We already knew the signature returns too much: in zebrafish, half of what it
returns are ryanodine receptors. Now we know it also returns too little. The
family's official badge is unreliable in both directions, so the census
cannot be one domain query — it needs profile models built to recognise the
deep branches, and until it has them, any claim that a lineage *lacks* an
IP₃ receptor is a claim about a database. *(pending: S2, S3)*

**The pore was found without being looked for.** We took the deposited
coordinates of a human IP₃ receptor and asked a purely geometric question:
going along the channel's axis, how close does any atom come to it? Two
narrow points appeared. The calculation had no access to the protein's
sequence — and yet the narrow points fall exactly on the two features the
literature names: the GGGVGD selectivity filter, and a gate made by a
phenylalanine and an isoleucine one helical turn apart. That the geometry and
the biochemistry agree, having been derived independently, is a good sign for
every structural measurement built on top of this.

The filter turns out to be **wide** — about 5 Å to the nearest atom, against
2.5 Å at the gate. That fits what the channel is for: it dumps calcium down a
steep gradient, so it needs throughput, not discrimination. This is not a
precision filter and the structure says so.

**The famous "100 Å" is only part of the distance.** The receptor's defining
puzzle is that IP₃ binds far from the gate it opens, and the number quoted is
about 100 Å. Measured on the ligand-bound structure, the binding site is
103 Å *above* the gate — but also 62 Å *sideways* from the channel's axis, so
the actual through-space separation is **120 Å**. The signal does not travel
straight down; it travels diagonally, through the armadillo solenoid and into
a neighbouring subunit. The larger number is the one the protein has to solve.

**The disease evidence is thinner than it reads.** Laying the clinical
variants onto the three genes forced us to record how precisely each is
actually located. Of nine catalogued entries, **two** name a specific amino
acid; the rest are placed only within a domain, or not placed at all. Both
precise ones are in ITPR3 — and one of them, the variant causing the severe
multisystem disease, sits **seven residues past the gate we measured**, on
the stretch of protein that mechanically couples ligand to pore. That is a
testable structural explanation for why the variant is so damaging, rather
than a coincidence. *(pending: S17)*

We also corrected a claim we were about to make: SCA29 is dominant and
missense, but its *direction* of effect is not established — one variant is a
gain of function — so it is now marked unresolved rather than lumped in with
the dominant-negative mechanisms.


---

## 2026-09-03 — S2: how many IP3 receptors are there, and where do they live?

**We counted the family.** Asking the domain databases for every protein
carrying one of the IP3 receptor's three diagnostic domains returns
**15,421 proteins from 1,488 organisms**. Just under half of them
(**6,433**) are IP3 receptors; slightly more (**6,807**) are ryanodine
receptors, the sister family that shares every one of those domains; and
**2,181** are fragments too partial to call either way. The near-even split
is worth pausing on: a search built to find IP3 receptors returns as many
ryanodine receptors as targets, which is exactly the hazard this project
was designed around.

**Where the family lives is more interesting than how big it is.** Outside
the animals the records are not scattered noise — they fall into a clean
pattern:

- In **green plants**, every single IP3-receptor call is in the **green
  algae** (*Chlamydomonas*, *Volvox*, and relatives — 11 species, one with
  the complete five-domain architecture). The **land-plant lineage has
  none**, despite contributing 15 proteins to the search.
- In **fungi**, every call sits in an **early-diverging phylum** —
  Mucoromycota, Chytridiomycota, Basidiobolomycota, Entomophthoromycota.
  The familiar fungi — yeasts, moulds, mushrooms — contribute **no
  proteins to the search at all**.

So the receptor is not simply "absent from plants and fungi", as the
textbook framing goes. It is present in the early branches of both
kingdoms and gone from the derived ones — the shape you would expect from
two independent losses. That is still a statement about what the protein
databases hold rather than about genomes themselves. *(pending: S20, S23)*

**We nearly built the census on the wrong domain.** The obvious way to
enumerate this family is to ask for the domain that defines it — the
IP3-binding core. Doing that would have missed **2,914 proteins, 758 of
them real IP3 receptors, across 385 species**, because plenty of family
members are annotated with only one or two of the three domains. The
sharpest case is *Dictyostelium* iplA, a characterised IP3 receptor that
carries neither the binding core nor two other family domains; last
session we predicted the census would miss it, and it does not, because the
census asks for any of three domains rather than the one. It comes back —
and then sits in the uncallable pile, because two domains out of five is
not enough evidence to call it anything. Resolving records like that is
what next session's profile search is for. *(pending: S3)*

**How we know the calls are right.** Each record is assigned by what
domains it carries, never by its name — so the names are free to serve as
an independent check. Across **6,191 proteins with an unambiguous gene
symbol, the domain-based call and the name agree every single time**. We
then re-ran the separation a second, completely different way, on raw
sequence similarity to six labelled human reference proteins, over one
representative from each of 32 phyla: **agreement wherever that test has
the power to decide (27 of 27)**, and the handful of coin-flips all fall in
the band where the sequence test explicitly declines to choose — which is
precisely where the deepest, oldest branches sit.

**A byproduct: seven broken database entries.** Because protein length is
recorded but never used to make a call, it was free to catch something
else. Seven records carry the *complete* five-domain architecture of a
full-length IP3 receptor packed into 1,528–1,993 residues — roughly 700
short of the shortest real one. An intact architecture in two-thirds of the
length means the gene model is truncated, and none of the seven is flagged
as a fragment by the database, because a truncated model submitted as a
whole protein is not marked as one. The biology is right and the record is
wrong. *(pending: S18)*

---

## S3 — 2026-09-04 — a second opinion on every protein, and it agrees

**Two completely different ways of asking "is this an IP3 receptor?" now
agree on 11,875 of 11,876 proteins.** The first, built last session, reads
what a database says about a protein: which domain signatures it carries,
and — crucially — which ryanodine-receptor-specific ones it does not. The
second, built this session, ignores annotation entirely and reads the
protein's own sequence against two statistical models, one distilled from
34 known IP3 receptors spanning humans to amoebae, the other from 22 known
ryanodine receptors. Every protein is scored against both; whichever model
fits better wins, but only if it wins by a clear margin.

They should agree, and they do — but that is worth something precisely
because they could have disagreed. One instrument reads a label, the other
reads the molecule. The two families they are separating share every domain
that defines either of them, which is the hazard this whole project is
built around.

**The single disagreement is a genuine correction.** A slime mould protein
from *Tieghemostelium lacteum*, 2,845 residues — IP3-receptor-sized, not
ryanodine-receptor-sized — was filed as a ryanodine receptor last session
because it carries a domain literally named "Ryanodine Receptor TM 4-6".
That domain is the channel's *pore*, which is the part the two families
share; naming it after one of them was a historical accident of which
protein got sequenced first. The sequence models are unambiguous: it scores
303 against the IP3-receptor model and 133 against the ryanodine one. Last
session's own write-up flagged that one of its domain rules was only
conditionally trustworthy. This is that caution turning out to be right,
once, in 11,876 chances.

**The uncallable pile shrinks by 70 %.** Last session left 3,361 records
that the domain rule could not call, because you cannot conclude anything
from a protein having *some* of the domains — an absence might mean the
protein lacks it or that nobody annotated it. A sequence model has no such
problem. **2,314 of those 3,361 now have a call** (1,379 IP3 receptors,
197 ryanodine receptors among the records already in the census). What is
left is 605 records that neither method can call — genuinely too fragmentary
to be anything to anyone.

**The search found 618 proteins the domain search never returned.** These
are proteins sitting in curated vertebrate reference proteomes — hagfish,
plaice, a bat, an Antarctic icefish — that the domain-based census missed
entirely. They are all short, 209–942 residues, which is the explanation:
a fragmentary gene model carries too few domains to be enumerated by a
domain query, but a sequence model recognises what is there. This is a
measure of what an annotation-derived census costs. *(pending: S18)*

**And in the other direction, nothing was missed.** The stronger claim is
the one that could have failed: of 2,787 IP3-receptor records from the
searched species that the sweep did *not* return, every single one was
checked against the search database itself — and **not one of them was in
the database and overlooked**. All 2,787 are entries the curated
per-species protein sets simply do not contain, so they were never searched
at all. Within its declared search space the sweep has no known blind spot.

**The trap we walked into, and out of.** The first run of the search called
**14,981 vertebrate proteins ryanodine receptors** — an absurd number, about
twenty per species. The cause: ryanodine receptors contain a small module
called SPRY that also sits in thousands of entirely unrelated proteins, so
troponins, calcium-binding proteins and ubiquitin ligases were all scoring
against the ryanodine model and winning by default, because the IP3 model
has no SPRY to compete with. The fix is a rule with a measured basis rather
than a chosen one: a match only counts as family evidence if it spans at
least 200 positions of the model, because 200 is the shortest we have ever
measured the family's own defining domain to be, and 137 is the longest we
have measured SPRY to be. There is a clean gap between those two numbers and
the threshold sits in it.

**A side-effect of that rule is the project's biggest annotation lead yet.**
1,794 of the proteins the gate turns away carry an IP3- or ryanodine-receptor
gene name of their own. They are not impostors — they are real family genes
whose database entries have been broken into pieces too small to recognise.
That is a concrete, per-species list of broken gene models to check against
the actual DNA. *(pending: S18, S10)*

**Three receptors per species is the norm.** Across 758 vertebrate species,
the commonest number of IP3-receptor records is exactly three — the three
paralogues — and species with many more turn out to be the ones with the
deepest annotation, not the ones with extra genes. Fifteen species have
none at all, and those fifteen have unusually thin protein sets (a median of
10,042 proteins, against a well-annotated ~20,000+). On this evidence they
are gaps in annotation, not absences of the gene, and the project's own
rules forbid calling them losses until the genomes themselves are searched.
*(pending: S4, S5)*

**One methodological embarrassment, caught and fixed.** Rebuilding the two
models after an unrelated edit produced *different models* from identical
input — the alignment program had been left to use as many processor cores
as it liked, and combines its results in whatever order they finish. A model
that changes when you rebuild it cannot be the model any published result
was produced with. Both are now built single-threaded and byte-identical
across repeated runs, and every file records its own checksum.

---

## 2026-09-04 — S3 addendum: the third seed finished

**The family looks the same from anywhere on the tree of life.** The
completeness argument was built by starting an iterative search from a
single receptor sequence and letting it pull in everything it could reach.
Three such searches were run, from a human, a fruit fly and *Acanthamoeba*
— an amoeba whose lineage split from ours long before animals existed. The
third took 7.6 hours and only finished after the last session ended. Its
result: all three searches recover **exactly the same family**. The same
1,656 full-architecture IP3 receptors, the same 952 partial ones, in every
run. Which sequence you start from does not change what the family turns
out to contain, and that is a far stronger statement about completeness
than any one search converging.

**What the deep-branching seed did change is how much rubbish it dragged
in.** The amoeba search ended up resting on 26,266 sequences against the
human search's 7,222, and 18,670 of those extras are proteins that neither
of the project's two receptor models recognises at all. Starting from a
sequence 1.5 billion years from the reference does not find more receptors
— it finds the same receptors plus a great deal of unrelated protein.

**And that failure had a shape nobody had anticipated.** The safeguard
built into this project watches for a search sliding into the ryanodine
receptors, the sister family that shares the same architecture. This search
did the opposite: it filled up with unrelated protein so fast that the
ryanodine fraction *fell* by 27 percentage points, making a badly drifting
search look, by that measure, cleaner than the ones that behaved. It was
caught only by the blunt rule that a search which has not settled after ten
rounds cannot support a completeness claim. The lesson is recorded as a
decision: keep the blunt rule, and always report what a search is actually
built from rather than only whether it stopped growing.

---

## 2026-09-04 — S4: the denominator, and the discovery that this is a bird problem

**Every claim this project makes about a missing gene will be a claim about
309 genomes, and about nothing else.** That list is now fixed and published:
one best assembly for each of 161 vertebrate orders, plus every species the
protein-level census could not settle. It comes to 552 billion bases of DNA.

**What the census could not settle turned out to be, overwhelmingly, birds.**
Of the 169 species pulled in because their receptor complement looks wrong —
too few copies, only broken fragments, or nothing at all — **130 are birds**,
against 21 ray-finned fishes and 10 mammals. This is the same finding S3
reported from the other side: bird genomes are described by roughly a third
as many proteins as mammal genomes, so the family's three receptors are
routinely missing from bird protein sets. Whether that means a bird has lost
a receptor gene or merely that nobody has annotated it is precisely what
searching the DNA will decide, and every one of those 130 birds is now in
the list to be searched. *(pending: S5)*

**Thirty-three of the 309 genomes have no gene list at all.** They are raw
sequence: assembled DNA nobody has yet marked up with where the genes are.
For those species the project may say what the DNA contains but must never
quote an annotation, and the manifest records which ones they are so the
distinction cannot be lost later.

**Three genomes are enormous and will dominate the work.** The African
lungfish assembly alone is 40 billion bases — thirteen times the human
genome — with a newt at 23 billion and a fire-bellied toad at 10 billion
close behind. Amphibian and lungfish genomes are famously bloated with
repetitive DNA, and searching them for a receptor gene costs accordingly.

**Nothing was dropped for want of a genome.** All 169 uncertain species have
a sequenced assembly available. Every absence this project eventually reports
will therefore be an absence measured in DNA, not an absence of data.

---

## 2026-09-04 — S5a: what a fragmented genome does to a missing gene

**The genomes chosen because their proteins looked incomplete are largely the
genomes that were never assembled well enough to show the gene.** This
project's list of 309 genomes was built in two halves: one representative per
vertebrate order, plus 169 "margin" species flagged because their protein
sets were missing an IP3 receptor, or held only fragments of one. Almost all
of those margin species are birds.

Before searching any of them, we asked a simple question of the genomes
themselves: are the assembled pieces even long enough to contain an IP3
receptor gene? An IP3 receptor spans roughly 80,000 to 500,000 bases of DNA.
If a genome has been assembled into pieces shorter than that, the gene cannot
sit whole on any one of them, and its apparent absence says nothing about the
animal.

**120 of the 309 genomes fail that test — and not at random.** Two thirds of
the birds fail it, against one in nine of the ray-finned fishes. Sixty-eight
per cent of the margin species fail it, against twelve per cent of the order
representatives. The suspicion that made these species interesting and the
technical shortcoming that would manufacture that same suspicion are very
largely the same set of genomes. Nothing about bird IP3 receptor loss can be
read from this evidence until the two are told apart.

**Worse, the bias points the same way as the biology would.** The three human
paralogs encode near-identical proteins but occupy wildly different amounts
of DNA — ITPR3 is compact, ITPR1 and ITPR2 sprawl over two to three times as
much. In the six-genome pilot, every paralog was found in every well-assembled
genome. In the poorly-assembled ones, the compact ITPR3 was still found in
two of three, while the two sprawling paralogs were found in none. A broken
assembly does not lose receptor genes at random: **it loses the big ones
first** — which is exactly the pattern that would be mistaken for birds having
lost ITPR1 and ITPR2. *(pending: S5b, which tests this across all 309.)*

**A ryanodine receptor was carried into every search as a tripwire, and it
worked.** Ryanodine receptors are the IP3 receptor's sister family, present in
every vertebrate genome in three copies. If a search finds none of them, the
problem is the search or the assembly, not the animal. In *Todus mexicanus* —
the Puerto Rican tody, the bird whose protein set contained no IP3 receptor at
all — the tripwire found a single ryanodine receptor where there should be
three, and even that one unnamed. That genome is not missing calcium channels.
It is missing continuous DNA.

**Where the sequences are good, the two families do not confuse each other at
all.** The ryanodine receptors have shadowed this project from the start:
they carry every protein signature that identifies an IP3 receptor, and at the
protein level they fooled an early version of our own detector completely.
Searched against genomic DNA instead, the confusion vanishes — across all 43
gene loci examined, each was matched by one family's sequences and simply not
by the other's. The ambiguity that dogs protein databases is not intrinsic to
the two families; it is a property of searching fragments rather than genes.

**One thing the databases are quietly getting wrong.** In both pufferfish and
pink cusk-eel, the two copies of ITPR1 that fish carry sit side by side in the
genome, one properly named and the other left with a placeholder identifier.
The same is true of two of their ryanodine receptors. These are real, complete,
correctly-placed genes that no search by name will ever find — which is why
every locus whose annotation does not already name it correctly is being
carried forward rather than discarded. *(pending: S18, the annotation audit.)*

**And a limit worth stating plainly.** Below the well-studied vertebrates, the
paralog names simply run out. Across the entire protein census there is no
complete, full-length, correctly-named ITPR1 or ITPR2 from any shark, ray or
chimaera, none from lampreys or hagfish, and no ITPR2 from the coelacanth
lineage. The search can still find those genes — it just cannot say which of
the three they are. That question belongs to the family tree, not to a name.

---

## 2026-09-05 — S5b: what 309 genomes say, read straight from the DNA

**Every one of the 309 genomes was searched, and the tripwire held in all of
them.** A ryanodine receptor — the IP3 receptor's sister, present in every
vertebrate — travelled with every search as a check that the search worked. It
was found in all 309. Nothing below is a case of a search failing quietly.

**Searched as DNA, the two families never once got confused.** Across 2,144
gene locations, every single one was matched by one family's sequences and
simply not by the other's. Not one was a close call. This is the family that
defeated our own detector completely at the protein level, where ryanodine
receptors carry every signature that identifies an IP3 receptor. The ambiguity
that has dogged these two families in protein databases is not a property of
the proteins — it is a property of searching fragments of them.

**Lampreys and hagfish have one IP3 receptor where other vertebrates have
three.** These are the only two genomes out of 309 where the search found
nothing at all for a paralog, and it found nothing twice in each: both carry
the equivalent of ITPR1 and neither carries ITPR2 or ITPR3. Both genomes are
well assembled and both passed the tripwire, so this is not a technical
failure. Jawless fishes split from the rest of the vertebrates before the
genome duplications that are thought to have produced most vertebrate gene
trios, so a single receptor there is what that history predicts. The honest
statement is *one findable receptor*, not *two specific absences*: no properly
named lamprey or hagfish ITPR2 or ITPR3 exists anywhere in the protein
databases for the search to use as a template. *(pending: S7, the family tree,
which is what turns this into a statement about genome duplication.)*

**Everywhere else, the receptors are there.** In genomes assembled well enough
to hold a gene of this size on one piece of DNA, all three paralogs were found
in 98–99 % of cases — sharks, rays, a coelacanth, a lungfish, a reedfish, and
every bird, fish, mammal, amphibian and reptile in the scope. **There is no
evidence of IP3 receptor loss anywhere in the jawed vertebrates.**

**And the apparent exceptions are a property of the assemblies.** In the 120
genomes too fragmented to hold the gene on one contig, recovery drops to
57–70 %. That is not biology; it is the same shortcoming reappearing, and it
falls hardest on exactly the species that were flagged as suspicious in the
first place. The bias even has a direction: the compact ITPR3 is recovered 13
points more often than the sprawling ITPR1 and ITPR2, because a short gene
fits on a short contig. A broken assembly loses the big receptors first, which
is precisely the pattern that would be misread as birds having lost two of
their three.

**485 IP3 receptor genes are effectively invisible to anyone searching by
name.** 318 exist only as DNA — no gene model at the locus, or an assembly
with no gene list at all — and another 167 sit inside a gene the databases do
record but never named, so no search for "ITPR" will ever return them. They
are real, complete, correctly placed genes.

**The databases are not equally good at the three receptors.** Comparing only
well-assembled genomes, so that assembly quality cannot masquerade as
annotation quality, a gene model names the right paralog 88 % of the time for
ITPR3 and 87 % for ITPR2 — but only **65 % for ITPR1**. These are genes of
near-identical protein length sitting in the same genomes, so a count built
from gene names will report a difference in copy number that does not exist.
*(pending: S18, the annotation audit.)*

**A caution about our own method, worth recording.** In the two largest
genomes — a newt and a lungfish, twenty and forty times the size of a typical
bird genome — the search initially reported up to ten copies of a receptor
where there is one. Each phantom was a fragment of the pore module, the part
these channels share with unrelated proteins, stitched across millions of
bases of empty DNA. The real gene was never in doubt, but the copy count was,
and copy count is a result. It is now filtered on sequence identity, using a
threshold with a wide empty gap on either side of it.

---

## 2026-09-05 — S20: where the IP₃ receptor actually lives

**The family is not a metazoan invention, and it is not universal either.**
Searching 6,928 reference proteomes — 63 million proteins, every non-vertebrate
eukaryote UniProt calls a reference, plus all the archaea and one bacterium
per genus — the receptor turns up in 662 of them, and the pattern is not
random. It is in every animal phylum we looked at, in the amoebae, in the
ciliates and oomycetes and dinoflagellates, in the haptophytes, in the
euglenids, in the green algae, and in several early-diverging fungal phyla.
It is in the two closest single-celled relatives of animals — a filasterean
and a choanoflagellate — which already carry **both** this receptor and its
sister the ryanodine receptor, full length. The two channels had already
split before animals existed.

**The land plants really do not have it.** This was the question the task was
set to answer, and it had an unglamorous possible answer: that the handful of
plant records in the databases were there because someone had annotated them,
and *Arabidopsis* had none because nobody had. That is not what happened.
Across **384 land-plant reference proteomes and 16 million proteins — every
flowering plant, moss, fern and conifer in the set — there is not one.** Its
green-algal cousins have it: 15 of 48 chlorophyte proteomes, including
*Chlamydomonas*. So the gene was there in the green lineage and the plants
lost it, somewhere on the way onto land.

**The yeasts and moulds have lost it too, and they are not alone.** Zero in
1,353 Ascomycota and Basidiomycota proteomes. But the fungal picture is
patchier than "early fungi kept it": it is present in the Mucoromycota (18 of
34), the chytrids (6 of 16), *Basidiobolus*, the Zoopagomycota and
*Entomophthora* — and absent from the Glomeromycota, the Mortierellomycota,
the Kickxellomycota and the microsporidia. The receptor has been lost
repeatedly and independently across the fungi, not once at their base.

**And it is absent from bacteria and archaea entirely** — 0 of 3,537 bacterial
and 0 of 634 archaeal proteomes.

**Why these absences can be believed.** An absence is only as good as the
search that failed to find anything, so each one was re-made with a
deliberately blunter instrument — the family's own Pfam domain models at a
sensitivity a thousand-fold looser — and checked against a positive control
inside the same search. In land plants the IP₃-binding domain that names the
family returns **nothing** substantial, while in the green algae next door it
returns 26. In the Dikarya it returns nothing, against 16 in the Mucoromycota.
And the shared MIR domain, which every eukaryote carries on other proteins,
comes back **633 times in land plants and 4,376 times in the Dikarya** — so
the search is demonstrably working in those very genomes. It is finding
everything except the receptor.

**The plant and fungal records are real genes, not database mistakes.** All 99
of them were chased one at a time. The obvious worry was contamination: a
sequencing project picks up an animal, its DNA lands in the assembly, and a
"green algal IP₃ receptor" turns out to be 99 % identical to a mouse. Every
record was therefore compared to its nearest relative outside its own kingdom.
The identities run from **20 % to 46 %, median 24 %** — the ordinary range for
genes that parted a billion years ago, and nowhere near the 95 % that would
signal a sequence in the wrong assembly. Not one record was a contaminant, and
not one lacked a genome to sit in. Forty-seven are complete genes; the other
fifty-two are real but fragmentary gene models, which says something about how
well algal and fungal genomes are annotated rather than about the receptor.

**One oddity worth chasing.** A single green alga, *Cymbomonas
tetramitiformis*, contributes twelve of the twenty-two complete plant genes.
Either that lineage has genuinely expanded the family, or its assembly is
duplicated. A proteome cannot tell those apart *(pending: S23, the
assembly-level search.)*

**Searching harder finds nothing new anywhere.** In all four groups — plants,
fungi, protists and invertebrates — the deeper iterative search settles on
**exactly** the number of receptors the ordinary search had already found: 54,
35, 729 and 1,194. It then kept going, in the protists' case swallowing another
twenty-four thousand proteins that are not receptors at all, and never found a
fifty-fifth green algal receptor or a seven-hundred-and-thirtieth protist one.
The ordinary search had them all.

**Searching harder does not find them either.** For each group the search was
run again the hard way — starting from a single receptor native to that group
and letting the model rebuild itself from whatever it found, round after round,
until it stopped changing. In the green algae it settled after five rounds on
seventy proteins, **every one of them a green alga**, including the six that
only this deeper search could reach. In the fungi it reached the yeasts and
moulds exactly twice in ten rounds, and both times what it reached was a
mannosyltransferase — a protein that shares one module with the receptor and
was already in this project's list of known look-alikes. So the absences are
not an artefact of not looking hard enough. They are absences.

**In the protists the deeper search found nothing new either, and failed in a
different way again.** It settled on 729 receptors by its third round — exactly
what the ordinary search had already found — and then spent seven more rounds
growing from two thousand proteins to twenty-six thousand without adding a
single further receptor. Most of what it swallowed were **apicomplexan
proteins**: the malaria parasites and their relatives, a group in which this
same sweep found no receptor in any of sixty genomes.

**And the deeper search failed in an instructive way.** The fungal run did not
settle at all: at round three the model grew thirty-four-fold, from 42 proteins
to 1,442, and by round ten it rested on 6,302 — while the number of actual
receptors in it went from 33 to 35. It had stopped being a model of the family
and become a model of fungal proteins in general. The rule that watches for
this project's expected failure — drifting into the sister family — saw
nothing at all, because the flood was of proteins belonging to *neither*
family, which dilutes the sister share rather than raising it. A second rule,
watching for explosive growth, caught it.

The protist run went wrong the same way but *slowly* — never more than
threefold in any single round — so the growth rule missed it too, and its
sister-family share actually **fell** while it drifted. Only the hard limit on
the number of rounds stopped it. The invertebrate run failed in yet another way — or rather, did not fail at
all. It also ran out of rounds, but nothing had gone wrong with it: three
quarters of what it held were still receptors of one family or the other, and
it was simply still finding a few more each round when the limit stopped it.
Four groups, four outcomes from the same rules: one settled, one caught by
growth, one caught by the ceiling after drifting badly, and one caught by the
ceiling while perfectly healthy. The rule cannot tell the last two apart,
because it counts rounds rather than looking at what the model holds — which is
the right way for a safety limit to behave, and a reason to look at the model
as well as the verdict. It is worth stating plainly: a model that has gone wrong does not
always go wrong in the way you are watching for, and the value of writing the
kill rules down in advance is that they still catch it when it does not.

**A caution about our own counting.** The sweep reports 142 "ryanodine
receptor" hits in fungi, green algae and protists, where no ryanodine receptor
is thought to exist. They are proteins of a few hundred to two thousand
residues matching **4–9 %** of a five-thousand-residue channel model — the
module these families share with each other and with unrelated proteins, not
the gene. Counted as genes they would invent a receptor family across half the
eukaryotic tree. Every call in this task is therefore reported with how much of
the model it actually covers: 63 % of the IP₃-receptor calls span at least half
their model, against 1 % of the ryanodine ones.

---

## 2026-09-05 — S23a: taking the absences to the genome, and what it takes to be allowed to

Everything this project has said about where the IP3 receptor is *missing* has
so far been said about protein databases. A protein database holds what a
gene-finder found. This task builds the machinery to ask the genomes
themselves, and the first thing it found is how easy it would have been to ask
badly.

**A negative needs a witness, and the witness we had does not live in plants.**
When we searched vertebrate genomes, every one of them could be checked: the
ryanodine receptor, the IP3 receptor's big cousin, is in every vertebrate
genome, so finding it proved the search had worked and a genome where it was
missing was a broken search rather than a discovery. Outside the animals that
cousin is essentially absent — we found it at full length in two of nearly
seven thousand non-animal species — so in a plant or a yeast its absence is
the right answer and it can witness nothing. Had we carried it over anyway,
"no IP3 receptor in *Arabidopsis*" and "the search never ran on *Arabidopsis*"
would have looked identical, and every absence in this task would have been
unfalsifiable while appearing carefully controlled.

The replacement is a protein that shares one of the receptor's four signature
domains and is carried by essentially every eukaryote on something else
entirely — the enzyme that attaches sugars to proteins. It is the same
molecule this project used as its hardest decoy in the very first control
experiment. Finding it proves the search reached the assembly; it must never
be mistaken for the receptor, and it wasn't.

**The pilot immediately showed why witnesses have to be graded rather than
counted.** In yeast, in *Neurospora*, in a microsporidian and in a chytrid the
witness appeared and the receptor did not — those absences are now genome
facts, not annotation facts. But in *Toxoplasma*, one of the headline targets,
**neither the receptor nor the witness turned up**, so that genome tells us
nothing at all. The apicomplexans, it turns out, barely carry the witness
protein either: one copy between them across sixty species. The honest
conclusion is that the apicomplexan absence cannot yet be taken to the genome,
and saying so is the point of having a control.

**The receptor's gene is a very different shape outside the vertebrates.** In
us it sprawls: three-quarters of a megabase for one gene in the largest case,
and it was the single biggest caveat on the vertebrate survey, because a third
of those assemblies are too fragmentary to hold a gene that long. Elsewhere it
is compact — a fly's is twenty-two thousand bases, a *Perkinsus*'s under four
thousand — and the whole range across the tree spans a hundredfold where the
three human copies span sixfold. Animal genes are about nine times longer than
those of protists and fungi. The practical consequence is good news for the
absences: the assemblies carrying them need only a fraction of the contiguity
vertebrate assemblies did, and almost all of them have it.

**Some lineages carry the receptor and no one has looked.** Deriving the list
of empty lineages from our own measurements rather than restating the previous
summary turned up three nobody had named: diatoms, red algae, and — the
interesting one — **tapeworms**, eleven species with no record between them.
Tapeworms are animals. Every animal absence this project has seen so far has
been a fragmentary genome rather than a missing gene, so if this one survives
the genome search it would be the first real loss inside the animals.
*(pending: S23b)*

**And the copies keep multiplying in odd places.** A flatworm carries
sixty-two records of a gene most animals have one of; a ciliate fifty-one;
*Paramecium* thirty-nine; a sponge thirty-seven; a green alga thirty-six.
Whether those are real gene families or the same gene written down many times
is exactly what a genome search settles, and thirty-six of these species are
now in the queue. *(pending: S23b)*

A closing note on our own instrument, in the spirit of the last one. Four
numerical thresholds were carried over from the vertebrate survey; **all four
turned out to describe the vertebrates rather than the family.** Three were
caught before the pilot — they had been measured on a set that is half insects
and, applied to the whole tree, they quietly discarded a genuine receptor from
a deep-branching flagellate. The fourth was caught by the pilot itself: in
*Chlamydomonas*, an alga we know carries a complete receptor, the search found
nine pieces of one and threw all nine away, because they resembled the bait
only a quarter as closely as a vertebrate gene resembles a vertebrate bait —
and the nearest bait available to it was an amoeba's. A threshold measured on
the best-studied corner of a tree describes the studying, not the tree.

---

## 2026-09-06 — What it takes to be allowed to say "it isn't there"

Most of this project's remaining claims are absences. The family is missing
from the flowering plants, from the mushrooms and moulds of Dikarya, from the
malaria parasites and their relatives — or so the protein databases say. A
database absence is a statement about what a gene-finding program noticed. To
turn it into a statement about the organism you have to go to the DNA, and to
be allowed to report *that* absence you need one more thing: proof the search
would have found the gene if it had been there.

That proof is called a positive control, and this session it nearly failed us
in the one place it mattered most.

**The malaria parasites.** The control we had been using is a different
protein that shares one domain with the receptor — a sugar-transferring enzyme
almost every organism carries. Find it, and the search demonstrably works.
Except that the apicomplexans, the group containing *Toxoplasma* and the
malaria parasites, barely have it: one such protein across sixty sequenced
species. So when *Toxoplasma* came back with no receptor *and* no control, we
could not tell whether the receptor is absent or the search had failed. The
absence claim about a clade of major medical importance was stuck.

The fix turned out not to be a better control protein. It was to stop
choosing the control in advance. **Which protein makes a good control is a
fact about the group you are asking about, and it can simply be measured**:
take a handful of large, ancient, well-conserved proteins, count how many of
the group's species actually have each one, and use the winner. For the
apicomplexans that is myosin — the motor protein muscles are built from,
which these parasites use to glide into host cells. Every one of the fifty-nine
apicomplexan species we checked has it, and it is a big multi-part gene, so
finding it proves the search can find a big multi-part gene here. *Toxoplasma*
is now controlled, and "no IP₃ receptor in Apicomplexa" is a claim we are
entitled to test.

The same measurement caught something we would otherwise have got wrong in the
opposite direction. **Red algae have largely lost myosin** — only one in twelve
species has one — so myosin would have been a terrible control for them. They
get a chromosome-maintenance protein instead, present in all twelve. Neither
choice was ours; both fell out of counting.

**And the control now has two strengths.** Because every group's control
protein is searched in every genome, we can see not just whether a control
fired but *how far away* the bait that found it came from. A control found
using a bait from the same group shows the assembly is readable. A control
found using a bait from a different kingdom shows the search reaches across
the kind of evolutionary distance any hidden receptor would have to be found
across. Eighty-seven of the first hundred and twenty-seven genomes clear the
higher bar.

**Meanwhile, the search itself.** Two thirds of the way through, the picture
so far: **every one of the thirty-four absence groups we can currently speak
for holds up in the DNA.** Thirty-one ascomycete genomes, seventeen
basidiomycete ones, six land plants — no receptor in any of them, with the
control firing in every case. These are still small genomes, since we are
working from the smallest up, and the animals are mostly still to come.
*(pending: S23c)*

Two positives stand out. A marine worm, *Dimorphilus gyrociliatus*, carries
**four** complete copies of a gene most invertebrates have one of. And one
fungus does have the receptor — *Basidiobolus*, an early-branching fungus that
sits well outside the Dikarya, so it sharpens rather than undermines the
fungal story: the family looks lost on the branch leading to the familiar
fungi, not absent from fungi as a whole. *(pending: S23c)*

**One more note on our own instrument**, continuing last session's. A
threshold we needed to measure could not be measured, because the pipeline was
throwing away the evidence before anyone could look at it — the search
discarded weak matches at the very cut-off we were trying to justify. It now
keeps everything and decides later. And a first attempt at the measurement,
run on a single genome as a test, produced a number, wrote it down, and was
promptly believed by the next run, which moved three organisms into a
different category on the strength of two data points. The measurement now
refuses to answer until it has enough to answer with. It is a small thing, but
it is the same failure as a fixed control: a procedure that always returns
something will always be believed.

---

## 2026-09-06 (cont.) — The answer, and the ruler that turned out to be wrong

The search finished: 194 genomes, a hundred billion letters of DNA, no
failures. Two kinds of thing came out of it.

**First, the absences are real.** Every one of the thirty-five groups the
databases said had no IP₃ receptor still has none when you go to the DNA, and
every one of them has a working positive control — so in each case we know the
search would have found the gene had it been there. The ascomycete fungi
(thirty-one genomes), the mushrooms and their relatives (seventeen), the land
plants (twenty-five), the flowering plants, the malaria parasites and
*Toxoplasma*, the microsporidian parasites: none of them has this receptor.
**Not one genome was excluded for want of a control**, which is the thing that
was in doubt yesterday.

So the picture the databases gave was, in this respect, right — and now it is
a statement about genomes rather than about gene-finding software.

**Second, the copy numbers.** Most invertebrates have one. But the range is far
wider than anyone would guess from the vertebrates' tidy three:

- a flatworm, *Macrostomum lignano* — **eighteen** complete genes
- a giant single-celled ciliate, *Stentor coeruleus* — **thirteen**
- sponges — six and eight
- a green alga, *Cymbomonas* — three

The flatworm number settles a question we flagged last session. The databases
held **sixty-two** records for that species, and there was no way to tell
whether that meant sixty-two genes or one gene written down sixty-two times.
It is eighteen genes. Something in that lineage has been copying this receptor
repeatedly, and eighteen is six times what a human has.

**And then the ruler.** This is the part worth telling on ourselves.

To decide whether a stretch of DNA is a copy of the gene, the vertebrate survey
used a simple test: how closely does it resemble the reference gene we searched
with? Below 40 % similar, discard it. That works when your references come from
close relatives. Across the whole eukaryotic tree it does not, and we measured
how badly. **The genuine genes reach down to 19 % similarity; the junk reaches
up to 32 %. There is no line you can draw between them.** Applying the
inherited cut-off would have thrown away eighty-seven real genes — fifty-six of
them complete, intact ones. A third of everything we found.

The reason is almost embarrassing once stated. Similarity to the nearest
reference measures *how far away the nearest reference is*. In the vertebrates
that is a few tens of millions of years. Out here it can be an entire phylum.
The number was never measuring what we thought it measured; it was measuring
our own sampling.

We checked whether a different simple measure would do better — how much of the
reference gene is covered — and it is no better. So the test was retired
rather than retuned, and the decision handed to the tool that was already
making every other family call in this project: a statistical model of the
receptor built from thirty examples. Against the only independent evidence
available — what the genome's own curators named each gene — that model agrees
with ten out of ten real ones and rejects ten of eleven impostors.

The eleventh is a nice ending. It is a gene in the alga *Emiliania huxleyi*
that the curators called **IPR1**, which our name-matching did not recognise
and so filed as an impostor. The model scored it as a receptor, emphatically.
`IPR1` is almost certainly short for *inositol trisphosphate receptor 1*. Our
list of names was short; the model was right. We have left the disagreement on
the record rather than quietly fixing the list, because a name that matches
"IPR1" also matches every InterPro accession number in existence, and a fix
that creates a thousand new errors is not a fix.

---

## 2026-09-06 — Lampreys have three of them, and we cannot yet say which

*(from S6, the alignment every later result stands on)*

The three human IP₃ receptors are supposed to come from the two rounds of
whole-genome duplication that happened at the base of the vertebrates. If
that is right, an animal that branched off *before* those duplications
should have one receptor, and everything after them should have three.

Lampreys and hagfishes are the animals that sit closest to that boundary.
They have **three**.

That was not the surprise. The surprise was that our own search cannot tell
them apart. Every one of the six loci we found — three in the sea lamprey,
three in the hagfish — is matched best by the same reference gene, human
ITPR1. Not one of them looks like an ITPR2 or an ITPR3. Six for six.

That could mean two very different things. Either the lamprey's three copies
are the same three genes we have, and the three simply have not diverged as
far in that lineage; or the lamprey duplicated its single ancestral receptor
on its own, separately from us, and the resemblance to three-ness is a
coincidence of counting. The first would confirm the textbook story. The
second would mean the textbook story rests on a coincidence.

We checked the obvious objection first: perhaps ITPR1 is just the most
conserved of the three, so *everything* distantly related looks most like
it. That would make the result meaningless. So we ran the same measurement
over animals that certainly do not have vertebrate paralogues — insects,
molluscs, worms, amoebae, algae — and over the ryanodine receptors. They do
drift towards ITPR1 slightly more often than chance, but by a margin of
about seven parts in a thousand: no signal at all. The lamprey and hagfish
margin is six times larger. The lean is real.

So the question is genuinely open, and it is now sharp enough to answer. The
tree we build next, and the gene-neighbourhood evidence after it, are what
separate the two stories. *(pending: S7 for the tree, S8 for the synteny)*

Two smaller things came out of the same work.

**Which two of the three are closest relatives** — a question the published
literature has never settled with a proper rooted analysis — has a first
answer from the alignment: **ITPR1 and ITPR2**, and not narrowly. Across
every one of the 165 comparisons between them, they are more similar to each
other than either is to ITPR3, and the middle half of those comparisons does
not overlap the middle half of the other two pairings. This is a strong
hypothesis, not yet a result: similarity is not ancestry, and the formal
test comes with the tree. *(pending: S7)*

**Some databases have put a vertebrate gene number on an amoeba.** Four
records outside the vertebrates are labelled "receptor type 2" or ITPR1 —
in a pond alga, an amoeba, a lancelet and a sea squirt. None of those
animals has vertebrate paralogues to be a type 2 *of*; the number has been
copied across from a human entry by automatic annotation. We have kept the
labels visible in our own tables rather than deleting them, and grouped
those sequences by what they actually are. It is a small thing, but it is
the same failure mode as the unnamed genes and the mis-typed pseudogenes:
a database asserting more than it knows. *(pending: S18, the annotation
audit)*

And one methodological note worth recording, because it nearly became a
wrong figure. When you cut the poorly-aligned columns out of an alignment,
the positions renumber. Our first version of the domain diagram forgot
that, and drew the receptor's binding core, its pore and everything between
at roughly two thirds of their true positions — with the pore itself pushed
off the end of the plot entirely. Nothing in the numbers looked wrong; the
picture just was. After the fix, the pore lands exactly on the single most
conserved stretch of the whole alignment, across a billion years of
evolution — which is where it should be, and which nothing in the
correction knew to aim for.

## 2026-09-07 — Which two of the three are closest relatives, at last

**ITPR2 and ITPR3 are each other's nearest relatives. ITPR1 is the odd one
out.** This is the question the published literature has never settled with
a proper rooted analysis, and the tree settles it decisively: the two
alternatives — pairing ITPR1 with either of the others — are rejected with
p-values around two in a hundred thousand, while the ITPR2/ITPR3 pairing
survives comfortably. Two independent lines of evidence in the same run
agree, and the ryanodine receptors root the tree so the statement is about
ancestry and not merely similarity.

**This overturns what our own alignment predicted.** Last session, raw
sequence similarity ranked ITPR1 and ITPR2 as the closest pair, and by a
margin that looked convincing. The tree rejects exactly that pairing
hardest. The lesson is an old one worth restating: similarity is not
ancestry. Two proteins can resemble each other because they are close
relatives, or because one of them has simply changed less than its
siblings, and only a model of how sequences actually evolve can tell those
apart. We flagged the similarity result as a hypothesis rather than an
answer at the time, which is why finding it wrong costs nothing.

**The hagfish and lamprey receptors are older than they looked.** These
jawless fish sit just outside the group where the three-receptor system
arose, so what they carry is a direct clue to how the family got its three
members. Each of the two species has three receptors. They do not, as we
half-expected, correspond one-to-one with ITPR1, ITPR2 and ITPR3 — nor are
they a recent local expansion. Instead they fall into two ancient
groupings, each containing both a hagfish and a lamprey gene. That means
the duplications that produced them happened *before* hagfish and lampreys
went their separate ways, several hundred million years ago. Which of them
corresponds to which of our three receptors is still open. *(pending: S8,
the gene-neighbourhood analysis)*

**A test we set ourselves failed, and the test was wrong, not the data.**
We had asked whether a fish carrying two copies of the same receptor puts
them side by side on the tree, treating anything else as a sign our naming
was wrong. Neither fish passed. But looking at what actually separates the
copies, each one sits with the *matching* copy from another fish species —
which is precisely what you see when the duplication is older than the
species themselves. It is: these copies date from a whole-genome
duplication early in fish evolution. The naming was right; our expectation
was not.

**Every disputed name held up.** Five records sit somewhere the tree cannot
confidently place. Checked against an entirely separate method — matching
each protein to the human genome and back again — all five come back
correctly named. So these are cases of the tree being unsure where to hang
a branch, not of a database being wrong about what a gene is.

**And a green alga that was missing from the picture.** A reader asked why
*Volvox carteri*, a well-studied alga known to carry one of these
receptors, was absent from our figure. It was not absent from our data —
we had found it, and confirmed it as a genuine gene. It had been dropped
from the illustration by a rule that judged how long a "typical" family
member should be. For the algae, that yardstick had quietly been borrowed
from animals, which are built to a different scale, and the alga was
discarded for being the wrong size against the wrong ruler. Fixed, the
whole analysis re-run from scratch, and *Volvox* now sits where it belongs,
beside its close relative *Chlamydomonas*. Every conclusion above survived
the re-run unchanged — which is the reassuring half of the story. The
sobering half is that the error was invisible in every table we had, and
surfaced only because someone looked at the picture and asked where a
familiar name had gone.

## 2026-09-07 — The three receptors are where they should be, and the neighbours say who is related to whom

Until now every statement about which of the three human receptor genes is
which has rested on one kind of evidence: how similar the proteins are. That
is a good instrument, but it is a single one, and the whole family is built
from copies of the same gene — so it is fair to ask whether the labels would
survive being checked a completely different way.

This week we checked them by looking at the *neighbours*. In a genome, genes
sit in stretches that stay together over very long spans of time. If our
ITPR1 in a mouse and our ITPR1 in a lizard really are the same gene inherited
from a common ancestor, they should still be surrounded by many of the same
neighbouring genes — even though nothing about the neighbours was used to
identify them.

**They are.** Across 274 annotated genomes, two loci we call the same
receptor share their neighbourhood two to four *hundred* times more than two
randomly chosen stretches of the same two genomes do, and this holds for
98–99.8 % of individual comparisons. Two loci we call *different* receptors
share essentially nothing — and neither do our receptors and the ryanodine
receptors, the look-alike family that has shadowed this project from the
start. The labels hold up under an instrument that never saw the protein.

**The neighbourhood also remembers the ancient duplication.** The three
vertebrate receptors were made by two rounds of whole-genome duplication
more than 500 million years ago. When a genome is duplicated, the genes
*beside* the receptor are duplicated too — so a family should have left
cousins beside two different receptors. We found exactly two such surviving
pairs, and both of them link ITPR1 to one of the others: a circadian-clock
gene (BHLHE40 beside ITPR1, BHLHE41 beside ITPR2) and a glutamate receptor
(GRM7 beside ITPR1, GRM4 beside ITPR3). Between ITPR2 and ITPR3 there is
nothing left at all.

That is worth stating carefully, because last session's family tree
concluded that ITPR2 and ITPR3 are each other's closest relatives — and they
are the one pair whose neighbourhoods retain nothing in common. The two
findings are not in conflict: a tree records the *order* in which the copies
were made, while a shared neighbour records which copies happened to *escape
deletion* afterwards, and those are different histories. But it does mean
the ancestral block that ITPR1 sits in is the one that survived best, and
that reconstructing the full picture will need the whole quartet of
duplicated regions, not the two pairs that happen to have lasted. *(pending:
a dedicated paralogon reconstruction — logged as an emergent task)*

**One receptor's ground is much less stable than the others'.** ITPR3's
neighbourhood is conserved within a class of animals but breaks down between
classes about two and a half times faster than ITPR1's or ITPR2's. In humans
ITPR3 sits inside the major histocompatibility region on chromosome 6 —
famously the most rearranged, most variable neighbourhood in the vertebrate
genome. So the gene that is *easiest* to find across genomes is the one
sitting on the shakiest ground, and those two facts are measuring different
things.

**131 genes gained a name from their neighbours.** Where a genome's own
annotation left a receptor unlabelled — or filed a second copy without
saying which of the three it was — the neighbourhood can often say. After
demanding that each assignment beat the best any random stretch of genome
achieved, 131 loci now carry a receptor identity they did not have, 87 of
them where the database says nothing at all.

**And the lamprey question stays open.** Last session's tree found that
lampreys and hagfish carry three receptors of their own, in lineages older
than the split between the two animals, and left it to this analysis to say
which of the vertebrate three each one corresponds to. It cannot. Matching
neighbours across species means matching gene *names*, and after roughly 550
million years of rearrangement — with only about a quarter of lamprey and
hagfish genes carrying a name at all — there is no shared vocabulary left to
compare. Every one of the six scores no better than a random stretch of the
same genome. That is not a negative answer; it is the honest report that
this instrument cannot reach the question, and the next one will have to
compare the neighbours by what they *are* rather than by what they are
called. *(pending: a name-independent synteny instrument)*

## 2026-09-07 (cont.) — How hard evolution has been holding on to these genes

Everything so far has measured *similarity*: how alike two receptors are.
That is a good instrument but it cannot separate two very different
histories. A protein can look almost unchanged because it has been protected
by natural selection for half a billion years — or because not much time has
passed. To tell those apart you have to compare the changes that alter the
protein against the changes that do not, and the ratio between them is the
closest thing biology has to a direct read-out of how strongly a gene is
being held in place.

**These receptors are held about as tightly as a gene can be.** For every
change that alters the protein, roughly **twenty-five to forty** silent
changes have accumulated instead. Put another way: almost every mutation that
would have changed one of these 2,700 amino acids has been removed by
selection. That is not a mild constraint. It is the profile of a protein
where nearly every position matters — consistent with a channel that has to
fold into a four-part assembly, bind its messenger at one end and open a
gate 100 Å away at the other.

**And the three copies are not held equally.** ITPR1 — the receptor that
carries almost all of this family's known human disease mutations — is the
most tightly constrained of the ones measured so far, ITPR2 the least. If
that ordering holds when the third receptor lands, it says something worth
saying: the gene where a single amino-acid change most often causes disease
in people is the same gene evolution has been most reluctant to let change at
all, across 450 million years of vertebrates. Those two facts come from
completely different places — one from clinics, one from genomes — and they
agree. *(pending: S9b, which finishes the third receptor and the statistical
tests)*

**The silent positions, meanwhile, have been completely rewritten.** This is
the surprise of the session. We expected the silent sites to be saturated —
turned over so many times that they carry no more information — *between* the
three receptors, which were made by whole-genome duplications more than 500
million years ago. They turn out to be saturated **within a single receptor**,
comparing a shark to a fish to a mammal. Over 90 % of such comparisons are
past the point where the silent clock can still be read.

The picture that gives is striking. Take human ITPR1 and its shark
counterpart: the protein is nearly the same protein, while the DNA
underneath it has been rewritten at essentially every position that was free
to change. The gene has been running in place, very fast, for a very long
time — and the protein has not moved. It also means every number in this
analysis has to be read as a firm *direction* rather than a precise value,
and we have logged the follow-up that would sharpen it: repeat the
measurement inside a shallower group, such as the mammals alone, where the
clock has not yet run out.

**The lampreys and hagfish drop out again — for a third distinct reason.**
Six receptor genes in those two ancient lineages have now been asked three
times which of ITPR1/2/3 they correspond to. The family tree put them in
their own separate groups and could not attach them. The neighbouring genes
could not say, because after 550 million years there is no shared vocabulary
of gene names left to compare. And now the selection analysis has to decline
as well, for a reason of its own: it can only ask a question about a group
the tree actually defines, and for these six it defines none. Three
instruments, three different kinds of silence. That is not a failure to
answer — it is the same answer arriving from three directions, that these
genes sit outside the vertebrate three rather than inside them.

**A footnote on how nearly this went wrong.** Building this analysis needs a
DNA sequence for each protein, and the databases offer several for the same
gene — one per alternative version of the protein. Taking the first one
offered silently substituted the *wrong version* of human ITPR1 and human
ITPR2, the two best-studied genes in the whole family. Nothing would have
crashed; the numbers would simply have been about a slightly different
protein than the one every other result in this project is about. The fix is
to translate every candidate and keep the one that actually matches. Human
ITPR1 needed nine tries.

## 2026-09-08 — One of the three has been guarded twice as closely as the others

Yesterday's entry could only compare two of the receptors. All three are now
measured, three different ways, and they agree.

**ITPR1 is held about twice as tightly as ITPR2 and ITPR3.** For every change
that alters the ITPR1 protein, roughly forty silent changes have accumulated
instead; for the other two it is about twenty-three. Measured on each
receptor by itself, measured again by testing each one against the other two
on the same tree, and measured a third time by a method that compares the
whole spread of rates rather than an average — the answer is the same every
time, and the third method puts it most sharply: **selection on ITPR1 is
*intensified* relative to its siblings, while ITPR2 and ITPR3 are *relaxed*
relative to theirs.**

That is worth pausing on. ITPR1 is the receptor that carries almost all of
this family's known human disease — the ataxias, Gillespie syndrome — while
ITPR2 and ITPR3 have a handful of families each. The clinical record and the
evolutionary record are completely independent kinds of evidence, and they
point at the same gene. A gene where a single amino-acid change causes
disease in people is a gene where a single amino-acid change was removed by
selection in every other vertebrate for 450 million years.

**Something happened on the branch that made ITPR1, and only on that one.**
When the ancestral receptor was duplicated, each new copy had a period alone
before it began diversifying — the branch on which a duplicate's fate is
decided. Asking whether any part of the protein was changing *faster than
neutral* on those three branches gives an answer for one of them: on the
ITPR1 branch, about a ninth of the protein was evolving several times faster
than the neutral rate, with eight individual positions identifiable with
high confidence — while everything else in the family sat at a fortieth of
neutral. This is the one place in the whole analysis where the receptors
look like they were being *changed* rather than merely preserved.

We are deliberately not making the same claim for ITPR2 and ITPR3. Their
tests come out significant too, but the rate itself cannot be measured on
those branches — the estimate runs off the end of what the method can
express, which means the data cannot say *how* fast, only that the model
prefers "fast". Significant and unmeasurable are different things, and only
the ITPR1 answer is a number. *(pending: S17, which will place those eight
positions on the cryo-EM structure — "a few sites changed fast" means
something very different in the messenger-binding pocket than in a floppy
linker)*

**And nothing in this family is under positive selection today.** The tests
that ask whether any position in a living receptor is being actively driven
come back with a clear negative once you look at what they actually fitted.
Two of them appear significant, and both turn out to be describing a small
fraction of positions — under one percent — that are simply *unguarded*,
drifting freely rather than being pushed. In a 2,700-residue protein where
everything else is locked down, less than one percent of slack is itself a
statement about how little of this channel is expendable.

**A note on how nearly the biggest number went missing.** The analysis that
found ITPR1's intensified selection wrote its answer into a file that could
not be read back, because one branch had an unmeasurably large rate and the
program spelled that as a word its own output format does not allow. The
analysis had succeeded; the reading of it failed, and silently — the result
came back blank, and a blank sitting beside two real answers looks like a
negative finding rather than a broken pipe. It is the same lesson as the
plant genome that vanished from a figure last week: the failures that matter
are the ones that produce something plausible rather than an error.

---

## 2026-09-08 — S10: the databases are missing genes that are demonstrably there

**The reassuring result first: the annotations are mostly right.** Across the
309 vertebrate genomes, there are 382 IP3-receptor genes sitting in
assemblies good enough that the annotation had every chance of describing
them properly. It described 375 of them correctly — one gene model covering
the whole gene. The family's public record is in good shape almost
everywhere, and the failures are not a gradient of sloppiness but a short
list: seven genes, in three genomes.

**But where it fails, it fails completely, and the gene disappears from
biology.** Two of those failures were taken apart exon by exon.

**A croaker with an invisible receptor.** *Nibea albiflora* is a food fish
with a good chromosome-level genome. Its ITPR2 gene is there — 56 coding
exons, a complete uninterrupted protein of 2,673 amino acids, spliced at
every junction the way a real gene is. The annotation has **nothing** on it.
Not a partial model, not a mislabelled one: the gene-finding step walked
right past it, having correctly described the genes on either side. And the
gene is not in some anonymous stretch of chromosome — it sits exactly
between SSPN and BHLHE41, the two genes that flank ITPR2 in 197 and 183 of
the 215 vertebrates we have looked at. This is the most recognisable
address in the genome for this gene, and the annotation still missed it.

**The same fish files its whole receptor family as broken.** Every complete
IP3-receptor and ryanodine-receptor gene in that genome is recorded as a
*pseudogene* — a dead gene, one that makes no protein. All of them are
intact: each reads through from start to stop with not a single premature
stop codon, where a genuinely dead gene of that length would have collected
a dozen or more. Nearly a third of the entire gene set of this fish is filed
as pseudogene. And one further error compounds it: the gene the annotation
does label "ITPR2" is not ITPR2 at all — its own protein sequence matches the
ITPR3 gene. So the species' ITPR3 carries ITPR2's name, and ITPR2 carries no
name at all. The net effect is stark: **search any protein database for an
IP3 receptor in this species and you get nothing.** Three complete,
functional receptor genes, zero protein records.

**An Antarctic toothfish with one gene recorded as three.** *Dissostichus
eleginoides* has its ITPR3 gene chopped into three separate "genes" that
tile the protein end to end — one covering residues 1–67, the next 52–467,
the third 468–1,593 — and then the final 41 % of the protein, including the
entire channel that does the actual work, has no gene model at all. Anyone
downloading this species' proteins gets three short fragments and no
receptor. Here too: zero IP3-receptor protein records for the species.

**Why this matters beyond two fish.** Every count of where a gene family
exists, every claim that a lineage has lost a gene, every alignment built
from database proteins, rests on annotations like these. A gene that is
present, intact and in the right place can be entirely absent from the
protein databases — and nothing on the protein side reveals it. That is why
this project searched DNA rather than proteins, and these two cases are what
that decision was for.

**A note on catching your own mistakes.** The test for whether a gene's
splice junctions appear in transcript data was checked against a case where
the answer has to be zero — the gene's own raw DNA, which by definition
cannot contain a spliced junction. It came back six. The test was too
lenient: the alignment program runs a little way past the true junction, far
enough to fool the rule. The rule was tightened and the control now reads
zero, which is what makes the real zero meaningful. A test that never fires
and a test that finds nothing look identical in the output.

## 2026-09-08 (cont.) — S11: the shape is right, and almost nobody has looked

**The receptors fold like receptors, and nothing else does.** Every complete
IP3 receptor structure in this project's panel — three cryo-EM references,
six conformational states of one of them, and eleven predicted models across
mammals, fish and a starfish — is identified as an IP3 receptor and not as a
ryanodine receptor by shape alone, with no gene name, no domain annotation
and no alignment involved. Twenty out of twenty. Three deliberately chosen
decoys — a molecular motor, a calcium channel and a cytoskeletal protein,
each about the same size as a receptor — are correctly identified as none of
the above. This is the sixth completely different kind of evidence to
separate these two families in this project, and the first that looks only
at the physical arrangement of atoms in space.

**The part that does the receptor's defining job is the part that is
modelled best.** AlphaFold's confidence, measured domain by domain along the
protein, is highest at the IP3-binding pocket at the N-terminus and *lowest*
at the channel pore. That ordering is fortunate: the pocket is what makes an
IP3 receptor an IP3 receptor rather than a ryanodine receptor, and it is
where the project's later questions about constraint and disease variants
are aimed. Those questions can be asked of predicted structures. Questions
about the pore cannot, or not as confidently.

**Almost none of this family has ever been folded.** AlphaFold DB — the
public database of predicted protein structures, which covers most of the
protein universe — holds a usable model for **13 of the 5,861 full-length IP3
receptor records this project has catalogued**. Two tenths of one percent.
The reason is mundane and complete: the prediction pipeline that built the
database stops at about 2,700 residues, and a vertebrate IP3 receptor is
about 2,700 residues. The sister family, at 5,000, is not in the database at
all. So the structural coverage of this family is not poor because it is
obscure — a great deal is known about it — but because it is *long*. Where
the database does hold models of these proteins, they are overwhelmingly the
short fragments: the median modelled record is 392 residues and the median
unmodelled one is 2,674.

**And where it does hold something, it may not be what you asked for.** Ask
AlphaFold DB for the human ITPR2 structure and it returns a model — of a
181-residue splice variant, not the 2,701-residue protein. Ask for ITPR1 and
it returns a 2,695-residue isoform rather than the 2,758-residue canonical
form. Only ITPR3 comes back as the protein you asked for. Anyone counting
"does this protein have a structure?" as a yes/no question gets three yeses
and one badly wrong picture of the family.

**A structure-based search of the whole reviewed proteome finds this family
and one other thing.** Searching every predicted structure of every reviewed
protein for anything shaped like an IP3 receptor returns the receptors
themselves and, apart from them, exactly one class of protein: SDF2 and
SDF2L1, in humans, mice, cattle, a plant and a slime mould. These are small
proteins built almost entirely from a single domain — the MIR domain — that
the receptors also carry. They match about 9 % of a receptor, extremely well.
Nothing else in the reviewed proteome is shaped like an IP3 receptor. It is
a modest negative, because the database being searched is missing most of
this family, but within what was searched it is clean.

**The deep branches remain unsettled, and honestly so.** The plant, fungal
and protist receptors — the ones this project has argued hardest for, since
they are held on sequence evidence alone — cannot be confirmed structurally,
because the only models available for them are partial (about 1,100 of 2,700
residues) and low-confidence. They lean towards the IP3 receptors over the
ryanodine receptors by about two to one, which is the right direction, but
they do not reach the bar. Crucially, this is a statement about the models
and not a verdict on the genes: each one had enough length to have passed if
it had matched well enough, and each got roughly two-thirds of the way. The
test did not fail here; it did not arrive. *(pending: predicted models
generated for the project would settle it — an emergent task.)*

**A note on negative controls doing their job.** The rule for deciding
whether a structure is an IP3 receptor or a ryanodine receptor was inherited
from earlier stages of this project, where it works. Applied to structures,
it confidently identified a molecular motor, a calcium channel and a talin
as IP3 receptors — because it asked only whether one answer beat the other,
and a clear winner between two wrong answers is still a clear winner. The
fix was to require the winner to be a *good* match before comparing at all.
The controls were the only thing that revealed this; without them the
report would have said twenty-six out of twenty-six agreed — and would have
called a dynein motor, a calcium channel and a talin IP3 receptors on the
same page.

---

## 2026-09-08 — S12: the genes the databases lose are real, and they are switched on

Earlier in this project a genome-by-genome search found IP3-receptor genes
that the public annotations of those genomes do not contain. Seven of them,
in three Antarctic and East Asian fish, sat in assemblies whose annotation
handles the other 375 of 382 family genes perfectly well. The obvious worry
about a finding like that is that the computer has drawn a gene where there
isn't one. This session put the question to biology instead of to software:
**if these are real genes, cells should be making messenger RNA from them.**

Half a billion sequencing reads later, from 67 public experiments across 16
independent studies and 11 tissues, the answer is that **all seven are real
genes, and they are being transcribed**.

**The evidence is the joins, not the gene.** A gene is written in the genome
in pieces, and the cell splices those pieces together to make a message.
The joins only exist in the message, never in the DNA — so a sequencing read
that runs straight across a join can only have come from a real transcript.
Of the 314 joins that no annotated gene model accounts for, **298 (95 %) are
crossed by reads**. That is the same rate as the joins the annotations *do*
account for, in the very same genes. Whatever these annotations are doing,
they are not responding to anything different about the DNA.

The genes are also not obscure. Each is switched on in six to nine different
tissues — brain, gill, gut, kidney, liver, muscle, ovary, skin, spleen — in
the amounts you would expect from an ordinary working gene.

**One species has no IP3 receptor on record at all.** *Nibea albiflora*, the
yellow drum, is a farmed food fish with a chromosome-level genome. Its three
IP3-receptor genes are all present in its DNA and all three are transcribed
in essentially every library examined. But one of them is missing from the
annotation entirely, and the other two are filed as **pseudogenes** — the
label for a broken relic of a gene, which means the databases serve no
protein for them. So a biologist searching any protein database for this
species' IP3 receptors finds nothing, and would reasonably conclude the
family had been lost. All three genes are there, intact, and working.

**Why nobody noticed.** There is another way to check whether a gene is
transcribed: look at the sequences researchers have already deposited for
that species. We tried it, on the species that should have been able to
answer — 37,166 deposited messenger RNAs for the Antarctic toothfish
*Dissostichus mawsoni*. It found nothing. Not the missing genes, and not the
ryanodine receptor either, a closely related gene that this species' own
annotation does name. The reason is in the lengths: those deposits are short
fragments, half a kilobase on average, and an IP3 receptor message is more
than eight kilobases long. The deposits were never capable of showing this.
That is worth stating plainly, because it is the shape of the whole problem:
**these genes are absent from the records not because they are absent from
the animals, but because they are long, and long genes are where automated
annotation and short-read deposition both fail.**

**On being wrong about a borrowed result.** This method was ported from an
earlier project on a different gene family, where a control designed to
measure false signal returned exactly zero, every time. It did not do that
here. Two libraries produced a small false signal, tightly confined to one
64-base stretch of one control sequence — real, bounded, and not enough to
change any conclusion. It would have been easy to write "zero" because zero
was what the method returned last time. The number is 42, and the report
says 42. Separately, an assumption that the method's earlier authors had
argued away in a footnote — that these very similar sibling genes cannot be
confused for one another — was tested directly here rather than argued, and
holds: not one of 12,500 test reads landed on the wrong gene.

---

## 2026-09-08 — Two duplications, two different moments (S13)

Vertebrates have three IP3 receptors where every invertebrate has one. That
much has been in textbooks for decades. What has never been shown is *when*
the two duplications that made the trio happened — or even which of the two
came first. This session placed both.

**They did not happen together.** The first split separated the ancestor of
ITPR1 from the ancestor of ITPR2 and ITPR3, and it sits on the branch
leading to all vertebrates — before the lineage that became hagfish and
lampreys parted from the one that became everything with jaws. That is
older than about 563 million years, and nothing in this analysis puts an
upper limit on it. The second split, which separated ITPR2 from ITPR3,
came later, on the branch leading to jawed vertebrates: somewhere between
about 462 and 563 million years ago. So the trio was not made in one event.
There was a two-receptor stage that lasted long enough to be visible in the
tree.

**The whole answer rests on six genes in two animals**, and this session
spent most of its effort on whether they can carry it. Hagfish and lamprey
each have three IP3 receptor genes. The tree pairs them up — one hagfish
gene with one lamprey gene, three times over — which means all three
lineages are older than the hagfish/lamprey split itself, roughly 459
million years. None of the six is a copy of ITPR1, ITPR2 or ITPR3. One pair
sits immediately beside the ITPR2/ITPR3 group; the other four form a group
of their own, off to one side. Remove those six genes from the analysis and
the older duplication slides forward onto the jawed-vertebrate branch. They
are the evidence, and there is no second copy of it.

So the obvious objection had to be answered: fast-evolving sequences get
pulled to the base of trees by a well-known artefact, and hagfish and
lamprey are the textbook case. They are not fast here. Measured against all
57 vertebrate sequences in the tree, the six sit within four per cent of the
median — two of them in the slower half. The artefact that would explain
away this result is not present in these sequences.

**A third ancestral lineage that only cyclostomes still have.** The
reconstruction says the vertebrate ancestor carried *three* IP3 receptor
lineages, not two: the one that became ITPR1, the one that became ITPR2 and
ITPR3, and a third represented today by four hagfish and lamprey genes and
by nothing at all in any jawed vertebrate. Either a fourth ancient copy was
lost early on the jawed-vertebrate branch and survives in these two animals,
or that group is something cyclostomes made themselves. Distinguishing them
means going looking for it in jawed genomes, which nothing in this project
has done, because every search bait it owns is labelled ITPR1, ITPR2, ITPR3
or ryanodine receptor *(pending: S15/S16)*.

**What this does not show.** The three receptors are widely assumed to come
from the two rounds of whole-genome duplication early in vertebrate history.
The older split lands exactly where that idea predicts. But the younger one
does not: if both rounds of duplication happened before hagfish and lampreys
branched off, ITPR2 and ITPR3 should each have a hagfish and lamprey
counterpart, and neither does. That leaves two possibilities the data cannot
yet separate — the ITPR2/ITPR3 duplication really did happen later, or
cyclostomes lost both copies. Calling the trio a product of whole-genome
duplication remains a reasonable expectation and not a demonstrated fact
*(pending: S16)*.

**On counting losses from a tree.** The reconstruction implies that this
gene has been lost about fifty times across the species examined. It has
not. Every one of those apparent losses was checked against the genome
sweep, and all of them are bookkeeping: mostly species whose gene simply was
not chosen for the alignment, plus a handful of species with no sequenced
genome in the project at all. Not one is a gene that is actually gone. The
number to remember is the genome sweep's, not the tree's.

**And one that looked like biology and was not.** Four of the apparent
losses initially survived the check: ITPR2 and ITPR3 in hagfish and in
lamprey, where the project's own genome ledger says "absent". Both animals
have three IP3 receptor genes. They read as absent because the search had no
hagfish or lamprey sequence to search *with*, so all three of each animal's
genes were filed under the same name and the other two slots were left
empty. An absence in a filing cabinet is not an absence in a genome, and the
analysis now says so explicitly rather than counting four losses that never
happened.

---

## 2026-09-08 — S15a: nobody has lost this gene

**The headline is a negative, and it is a strong one.** Across 309
vertebrate genomes — from hagfish to hummingbirds — there is not a single
case of a species that has lost one of the three IP3 receptor genes. Every
one of the 927 species-by-gene slots examined comes back with the gene
present. Not once does the evidence support absence.

That is worth stating carefully, because the raw data looked otherwise. The
genome search had left roughly fifty slots unresolved: four reading "absent",
forty-four where the search found scattered pieces of the gene but could not
assemble them into one, and a hundred and twenty where the gene was cut off
part-way. On the face of it that is a lot of candidate gene loss. All of it
dissolves once you ask what the genome assembly was capable of showing.

**The forty-four scattered cases are whole genes, shattered by the assembly
rather than by evolution.** These are almost all bird genomes sequenced a
decade ago, and their problem is that the DNA fragments they are built from
are shorter than the gene itself. An IP3 receptor gene runs to a hundred
thousand base pairs or more; the typical fragment in these genomes stops
after forty thousand. So the gene is genuinely in there — it just is not in
one piece anywhere. Reassembling it across the fragments recovers, on
average, four-fifths of the protein, spread over about six separate pieces
of DNA. Not one of the forty-four is missing.

**And the check that this is not wishful thinking.** The three IP3 receptors
resemble each other closely enough — about two residues in three — that
fragments of one can be mistaken for another, which would make any gene look
present whether it was or not. So the test was calibrated against cases
where that mistake is impossible: fragments belonging to a gene the search
had *already* found sitting intact elsewhere in the same genome. Those score
3 % of the protein. The real candidates score 80 %. There is no overlap
between the two, and the calibration was rebuilt once, after the first
version put real shattered genes into the comparison group and blurred the
answer.

**How many copies each animal actually has.** In every one of the 189
genomes good enough to hold the gene on a single piece of DNA, all three
copies are there — three, exactly, with no exceptions. Where more turn up it
is either a known extra round of genome duplication (Atlantic salmon has
eight) or a fragmented assembly counting the same gene twice.

**One thing did not work, and the number matters.** The plan was to settle
the unresolved cases by looking at the genes on either side — a
neighbourhood, which is stable over hundreds of millions of years. That
instrument works: where it can be applied it is right every time, and it
agreed with the sequence evidence on every case it reached. It reached eight
of four hundred and thirty-two. In the rest there is no neighbourhood to
look at — for two hundred and seventy-three, the piece of DNA carrying the
gene fragment has no other gene on it at all. The reason the gene is
unresolved is the same reason its neighbourhood cannot be read.

**And one genuinely new lead.** The gene models recovered for ITPR3 carry
more small disruptions to their reading frame — single-base insertions and
deletions — than ITPR1 or ITPR2 do in the same animal. That comparison was
made within each genome, against sibling genes matched for how well the
search sequence fit, so it is not an artefact of one gene being harder to
find. ITPR2 showed the same pattern until that matching was applied, and
then it vanished; ITPR3's survives. Whether this means ITPR3's coding
sequence genuinely tolerates more disruption, or that its gene structure is
harder to read through, is not yet known *(pending: S21/S22)*.

**What has not been shown.** No gene loss found is not the same as no gene
loss possible — nothing here rules out losses in the many vertebrate species
with no sequenced genome, or in a lineage the project's search sequences
cannot reach, which is still the case for hagfish and lampreys. And these
are three genes whose loss in a mouse is lethal or nearly so, so a complete
absence of loss across the vertebrates is what one would have predicted.
What is new is that it has now been checked, one genome at a time, with the
assembly's own limitations counted as limitations rather than as absences.

---

## 2026-09-08 — S15b: how hard you have to try to find a loss

The previous session established that no vertebrate genome supplies evidence
that an IP3 receptor gene is missing. This session asked the obvious next
question — what would it take to make one appear? — and the answer turns out
to be more interesting than the count.

**The count.** Zero. Placed formally on a tree of all three hundred and nine
genomes, with the method that would have found a loss had one been there: it
was given constructed losses on known branches and found them, merged two
losses in sister lineages into the single event they would have been, and
refused to over-count losses under parts of the tree that are not fully
resolved. Nine hundred and twenty-seven gene-in-genome questions, no absence.

**How hard it is to manufacture one.** The rules that decide whether a gene
counts as present were loosened and tightened across thirty-two different
settings, from the most generous the evidence allows to the strictest — one
that accepts nothing except a complete, uninterrupted gene model and ignores
whether the assembly could even hold such a thing. Asking only "does this
animal have an IP3 receptor at all", a loss appears in **two of thirty-two**
settings, and the worst case is one animal in three hundred and nine. Asking
"does this animal have *this particular* one of the three", a loss appears in
**eighteen of thirty-two**, and at the extreme setting there are forty-five
of them. None of those forty-five is real. What the comparison shows is that
a statement about the family is robust and a statement about an individual
paralog in a broken assembly is not — which is why the family-level statement
is the one this project makes.

**Two things that change nothing.** The threshold for accepting a gene
reassembled from pieces scattered across an assembly was moved across the
entire range its own calibration leaves open, and not a single one of the
nine hundred and twenty-seven answers changed. And the assumed timescale of
the vertebrate tree — whether branches are counted as steps, spaced evenly,
or dated from the fossil record — cannot change a count of losses at all, for
the simple reason that the counting method looks at branching and not at
time. That was worth checking rather than assuming.

**One thing that changes a great deal.** Refusing to make one specific
allowance — that a genome carrying spare, unidentifiable copies of the gene
family is not evidence that a particular copy is gone — immediately produces
four losses, in the hagfish and the lamprey. Those two animals each carry
three IP3 receptor genes. The four "losses" would have been an artefact of
the search's own inability to tell which of their three is which.

**The one cell that a sceptic should be handed.** Of all nine hundred and
twenty-seven, exactly one becomes an absence when the evidence bar is
tightened while every protective allowance is left in place: the second IP3
receptor of the jararaca pit viper. Its gene is there — it reassembles across
six pieces of the assembly at four-fifths of full length, in a genome
contiguous enough to have held it whole — but it is the only gene in the whole
survey whose presence rests on that reassembly and on nothing else.

**No rate of loss can be quoted, and that is a result rather than a gap.**
Standard models for gene gain and loss estimate how often the gene is lost.
On a character where it has never been lost there is nothing to estimate, and
a model asked anyway will return whatever number it started from, dressed up
as an answer. This was measured rather than asserted: the fit was traced
across eight orders of magnitude of loss rate and it slides monotonically to
the edge every time. Where the loosened settings *do* manufacture losses, the
fitted rate simply tracks how many were manufactured — a property of the
filter, not of the animals.

**No dead genes.** A gene that died long ago should leave a corpse: a
recognisable but broken copy, its reading frame full of stops. Of one thousand
seven hundred and sixty gene models examined, forty-four carry enough small
disruptions to be worth a second look — and all forty-four are complete
genes. Deliberately using the most generous definition of a corpse available,
seven IP3 receptor loci qualify, and every one of them is a full-length gene
carrying one or two stray stops. There is no test to run for shared damage
between related dead copies, because there are no dead copies.

**The one new biological lead, now localised.** The last session found that
ITPR3's gene models carry more small reading-frame disruptions than its
sibling genes do in the same animal, and could not say why. Splitting that
comparison by animal group gives a sharp answer: **it is a bird result.**
In twenty-five of twenty-seven birds, ITPR3 carries more disruptions than its
siblings; in ray-finned fishes the same test is seven to six, which is
nothing. The catch is that birds have the most fragmented genome assemblies
of any group here, and twenty-one of those twenty-seven birds sit below the
quality bar this project uses. The six that sit above it all point the same
way and none points against — but six is not enough to settle it. So the
lineage is now named and the mechanism is not, and the next step is more
birds with better assemblies rather than a cleverer statistic
*(pending: S18/S19)*.

**An echo worth noting, not leaning on.** Asked where on the tree each of the
three genes first appears, this method — which reads no gene tree and no
alignment, only which animals have which gene — puts ITPR1 at the origin of
the vertebrates and ITPR2 and ITPR3 at the origin of the jawed vertebrates.
That is exactly where the reconciliation two sessions ago placed the two
duplications. It agrees because it is partly reading the same evidence from
another side, so it is a consistency check and not a second confirmation.

---

## 2026-09-08 — Where the three receptors came from, and what a fish did to one of them (S16)

**The three human IP3 receptors are not three unrelated genes that happen to
be similar. They are three pieces of one ancient neighbourhood, and you can
still see the seams.** When the vertebrates were young, the whole genome was
duplicated twice. If the three receptors came out of that, then the stretches
of chromosome around them should not merely look alike — they should be
*related copies* of one another, still carrying descendants of the same
neighbouring genes under different names. They do. In 141 of 175 vertebrate
genomes the genes flanking ITPR1 have relatives among the genes flanking
ITPR2, and in 89 of 152 they have relatives among the genes flanking ITPR3.
In matched random stretches of the same genomes, picked the same way and read
with the same dictionary, that happens 2.6 % of the time.

**But ITPR2 and ITPR3 are not each other's neighbours in this sense at all.**
Their two neighbourhoods share relatives at exactly the background rate — 4
genomes in 149, against a background of 2.6 %. Everything the genome
remembers about this family's origin runs through ITPR1. That is a strange
result to sit beside the family tree, which says ITPR2 and ITPR3 are the
closest pair. The two are not in conflict: a tree records the order in which
genes split, while a surviving neighbour records which copies happened to
escape deletion afterwards, and those are different histories. But it is
worth saying plainly that the genome's memory and the gene's memory point
different ways.

**One of those two links is old enough to count, and one is not.** The
previous session found two families of neighbouring genes shared between the
receptors and could not say when they were duplicated. This session could.
The pair beside ITPR1 and ITPR3 — the glutamate receptors GRM7 and GRM4 —
was duplicated at the origin of the vertebrates, which is exactly when the
whole-genome duplications happened, and the same link shows up in 84
genomes. The pair beside ITPR1 and ITPR2 — the clock genes BHLHE40 and
BHLHE41 — was duplicated hundreds of millions of years earlier, before
animals and fungi parted company. They are relatives sitting in the right
places for the wrong reason. Only the first is evidence about the vertebrate
duplications.

**The best control this project has is the family it spends most of its time
trying not to confuse itself with.** The ryanodine receptors carry every
domain the IP3 receptors do, they come in threes in vertebrates, and they are
just as old. Everything above was measured on them too, in the same genomes,
in the same run. They give the same answer at the same strength — which is
how we know that a weak signal here is a limit of what one human genome can
show, not a sign that the IP3 receptors are different.

**Then there is the fish.** About 320 million years ago the ancestor of the
teleost fishes duplicated its entire genome again. Every gene was copied,
including all three receptors. Most duplicates are quietly lost afterwards,
and that is what happened to two of them: in 73 well-assembled ray-finned
genomes, ITPR2 and ITPR3 are back to single copies. **ITPR1 is not. It kept
both copies, in 97 % of those genomes.** The ryanodine receptors in the very
same fish kept all six of theirs, so this is not the search failing to find
things — it is a choice the fish genome made about this one gene.

**The two fish copies really are the two halves of one event.** Each of them
sits in part of the neighbourhood their single ancestor occupied, and — the
part that matters — they *divide* it between them rather than both keeping
the same bits: 45 of 49 genomes show a clean partition against a
land-vertebrate reference and 46 of 49 against gar and bowfin, the fish that
split off before the duplication and still carry one copy. And they are the
*same* two halves everywhere. Take any well-mapped fish as a template and ask
which of its two copies each other fish's copies correspond to, and all 705
comparisons agree. If ITPR1 had simply been duplicated over and over in
different fish lineages, those answers would have been coin flips.

**What is left open.** Nobody knows why ITPR1 alone was worth keeping twice.
It is now the paralog that stands out in two independent ways — the one whose
neighbourhood kept its ancient relatives, and the one whose fish duplicate
survived — and whether those are one fact or two is a question about what the
protein is doing, not about where it sits *(pending: S17)*. Which of the two
ancient duplication rounds made which split is also unresolved, and it needs
the lampreys and hagfish, whose gene names are too sparse to read
neighbourhoods from.

---

## 2026-09-08 — The five residues that never changed, and the fifty that did (S17)

**The IP3 receptor is a gate at the end of a very long lever.** Calcium
leaves the ER through a pore about 2,500 residues downstream of the site
where IP3 binds. So the obvious question about the protein is which parts of
it are actually load-bearing — and now that there are 250 orthologues of each
of the three human copies, that question has an answer per residue rather
than per domain.

**The gate has not changed.** Five residues line the point where the channel
closes, and across the whole vertebrate sample they are the most invariant
sequence in the protein. More striking: **they are identical in all three
human receptors.** ITPR1, ITPR2 and ITPR3 have been separate genes since
before jawed vertebrates existed, they are only 64–70 % identical overall,
and they still share those five residues exactly. Whatever the three copies
were free to diverge in over half a billion years, it was not the gate. The
selectivity filter just below it is nearly as conserved (86–100 % identical
between copies), and the ten residues that actually touch IP3 are more
conserved than the rest of the domains they sit in.

**The domain named after IP3 does not contain the IP3 site.** The protein
family's signature domain — the one the databases call "Inositol
1,4,5-trisphosphate/ryanodine receptor" — is the N-terminal β-trefoil, and
**not one of the ten measured IP3 contacts is inside it.** They sit in the
two domains next door, MIR and RIH, which are also the two domains the
ryanodine receptors share. So if you had trusted the domain name and measured
constraint on it, you would have measured the suppressor domain and called it
the ligand site. This is worth stating plainly because the name is in every
database record for every one of the ~15,000 proteins in this project's
census.

**And there is one stretch that has changed enormously.** Inside the pore
domain, on the far side of the membrane from the cytoplasm, sits a
50-residue loop that hangs into the ER lumen. It is the **least** conserved
sequence in the whole receptor: where the protein overall is 64–70 %
identical between the three copies, this loop is **13–31 %**. Two completely
different instruments agree — residue conservation across 250 orthologues,
and the rate of protein-changing substitutions estimated from DNA — and it is
also the region cryo-EM maps resolve worst, and the region the standard
alignment-cleaning step deletes.

So the most conserved sequence in the receptor and the least conserved sit
about fifty residues apart, inside the same domain. Before that loop was
separated out, its variability was dragging the whole pore domain's average
below the receptor's floppy inter-domain linkers — which would have meant
reporting that the pore of an ion channel is less constrained than its
spacers. It is not; the average was hiding two opposite things. What the loop
is *for* is not something this project has looked at *(pending: an emergent
task)*; a paralog-specific luminal insert sitting next to an invariant gate
is the shape of a functional difference between the three copies, but that is
a hypothesis, not a result.

**Does any of this help read a patient's variant?** Partly, and the honest
answer has two halves.

The first half: **conservation does separate the variants that matter from
the ones that do not.** Taking ClinVar's own labels, positions called
pathogenic are markedly more conserved than positions called benign — well
enough to be a usable signal. But the layer that does it best is not the one
this task was built around. Conservation measured across *all* IP3 receptors,
from protists to humans, beats conservation measured deeply within a single
human gene. A position that has survived both 500 million years of the three
copies diverging *and* the whole eukaryotic history of the family is a better
predictor than a position that merely never varies among 260 vertebrate
copies of the same gene. That is worth knowing before anyone builds a
prediction tool for this family.

The second half is less comfortable. **88 % of the family's variant record is
uncertain.** Of 1,753 missense variants in ClinVar across the three genes,
1,546 have no interpretation at all. ITPR1 carries 49 pathogenic ones; ITPR3
five; **ITPR2 has one.** That is not because ITPR2 tolerates change — its
gate is identical to the others', its IP3 contacts are the most conserved of
the three, and 8 of those 10 residues are invariant across 249 species. It is
because nobody has been looking at ITPR2. It has no dominant childhood
syndrome pulling sequencing towards it, so its variants stay unclassified.
The prediction that follows is testable: variants at ITPR2's gate, filter and
IP3-contact positions should turn out to be pathogenic when someone looks
*(pending: nothing in this project can settle it — it needs patients)*.

**One small check worth recording.** The recurrent ITPR3 variant that causes
a multisystem disease with immunodeficiency sits seven residues past the
gate. It is not in a tolerant position: that residue is invariant across all
264 orthologues measured. And the two variants in the whole published
literature that name a specific residue were both found again in the database
harvest — which is the kind of check that is only worth anything because it
could have failed.

---

## 2026-09-08 — The gene is fine. The paperwork is not (S18)

We went looking for a badly-recorded gene family and found a badly-recorded
*archive*.

The expectation was reasonable. An IP3 receptor is about 2,700 amino acids
long, spread over dozens of exons across tens or hundreds of kilobases of
DNA, and it has a sister family — the ryanodine receptors — that carries
every one of its diagnostic domains. If any gene family were going to be
systematically mangled by automatic gene-finding software, it should be this
one.

So we checked, across 309 vertebrate genomes: for every IP3-receptor gene we
had found ourselves, what does that genome's own official gene list say is
there? And crucially, we asked the same question of the ryanodine receptors
sitting in the same genomes, annotated by the same software.

**They come out the same.** About 74 % of IP3-receptor genes are delivered as
one complete, correct gene model, and about 78 % of ryanodine-receptor genes
are too. The difference between them does not survive statistical
correction. Whatever is wrong here is not wrong about
this family — it is how genes of this size get recorded.

**What does differ, enormously, is who wrote the gene list.** Genomes whose
annotation comes from the curated RefSeq pipeline get this gene right about
99 % of the time. Genomes whose annotation was deposited by the sequencing
group get it right about 37 % of the time. Part of that is the assemblies
themselves — submitter genomes are more broken, and no annotation can build a
whole gene across a piece of DNA too short to hold it — but only about a
third of the gap. The rest is the gene list.

The assembly matters enormously on its own, too. If we look only at genomes
whose DNA is contiguous enough to actually carry the gene, the failure rate
falls from 26 % to 7 %. Three quarters of what looks like careless annotation
is really a fragmented genome.

**There is one failure that is this family's own**, and it is a strange one.
An IP3-receptor gene is two and a half times more likely than a ryanodine
receptor in the same genome to be recorded as *non-coding*. The annotation
finds the gene. It names the gene correctly — "Inositol 1,4,5-trisphosphate
receptor type 3", spelled out in full. And then it files it as a
non-protein-coding feature, sometimes with a note that it contains a
frameshift. The consequence is that no protein record is ever created, and
nobody searching a protein database by name will ever find it. Forty-two
genes in this survey are in that state, against sixteen ryanodine receptors.
The honest qualification is that almost all of them sit in badly fragmented
genomes — restrict the comparison to good assemblies and it is four against
five, which is nothing — so this may be a fact about which assemblies these
genes happen to be in rather than about the genes *(pending: the exon-level
machinery built in S10 has never been pointed at one, and it would also
settle whether the submitters are right about the frameshift)*.

**The proteins that do exist are named correctly and cannot be found.** We
took every full-length IP3-receptor and ryanodine-receptor protein record in
the databases — 11,402 of them — and asked whether the name matches the
sequence. It almost always does. Only five records in eleven thousand are
named for the wrong family, and all five are obscure single-celled organisms
where the sequence itself barely separates the two families. The specific
type — 1, 2 or 3 — is right too, in all but a few dozen cases.

But **more than half of these records have no usable gene name at all.** They
carry an internal database code, or nothing. A researcher searching for "the
IP3 receptors" by name reaches under half of what the databases actually
hold. That is the same problem we found on the DNA side, arriving from the
other direction.

**And a small, clean result about missing genes.** An earlier stage of this
project searched 764 vertebrate protein catalogues and found fifteen with no
IP3 receptor whatsoever — a gene every vertebrate is supposed to have. We
went back to the raw DNA of each of those fifteen species. **All fifteen have
the gene.** Not one is a genuine absence. Eleven of the fifteen are birds,
whose genomes are notoriously fragmented. When a protein catalogue says a
vertebrate is missing an IP3 receptor, it is telling you about the software,
not the animal.

The practical output is a list of 297 specific corrections — this gene, this
genome, these coordinates, this is what the record says, this is what it
should say — of which 52 are cases where the gene is demonstrably present,
the DNA is demonstrably intact, and the reading frame is unbroken. Eighteen
more were deliberately withheld: where our own earlier analysis says the gene
really is dead, we are not going to ask a database to bring it back.

---

## 2026-09-08 — S19: how much of what we found, we nearly missed

This entry is about the search rather than the receptor, but it changes how
much weight every earlier finding can carry, so it belongs in the story.

**We can say exactly how often our own search failed, and almost nobody can.**
Measuring a gene hunt's sensitivity normally means estimating it, because the
genes you missed are the ones you cannot count. This family gave us a way
round that. Across 309 vertebrate genomes we found that not one of the three
IP3 receptors has ever been lost — every genome has all three. So any cell in
our results where the search came up empty is not a missing gene. It is a
missing answer. **The search failed on one cell in seven: 140 of 923.**

**We checked that number against a second gene family and got the same
answer.** The ryanodine receptors were searched in the same genomes by the
same tool and are present in every vertebrate too. They failed at 42 of 309 —
13.6 % against 15.2 %, statistically the same. The two measurements share no
reasoning, so their agreement is real evidence rather than the method
agreeing with itself.

**Every single failure is the genome, not the gene.** Where the search missed
a gene, the assembly was on average a hundred and fifty times more
fragmented: a typical missed gene sits in an assembly whose contiguous pieces
average 23,000 letters, a typical found one in an assembly averaging 3.4
million. On the best-assembled genomes — the ones built chromosome by
chromosome — the search missed three IP3 receptor genes out of 512, and no
ryanodine receptor at all. **Above the quality bar this project set at the
outset, the failure rate is 0.9 %.** That bar was chosen from the size of the
gene itself, before anybody had measured an error rate. It turns out to land
almost exactly where the measurement would have put it.

**The gene we recover best from a broken genome is the shortest one.** Below
the quality bar, ITPR3 is found 70 % of the time, ITPR1 61 %, ITPR2 57 % —
and that is precisely the order of their lengths on the chromosome. This
matters because it points the same way as a real biological signal would. Any
future claim that one paralog is more often lost than another has to clear
this bias first.

**Three quarters of the IP3 receptor genes we found are not in any protein
database.** Of 1,232 genes we demonstrated in these genomes, **940 cannot be
reached by searching protein records at all.** For 386 the species has no
protein catalogue; for 248 the only records are fragments; for 286 records
exist but none can be assigned to that particular paralog. This is not an
artefact of having deliberately chosen poorly-served species — the rate is
74 % even among the genomes picked simply as representatives of their order.
The genes are there in the DNA. Nobody has written them down.

**A surprise about how to search for a gene family.** We built a 38-protein
bait panel spanning the vertebrates — mammal, bird, reptile, amphibian,
fish, shark, lamprey — on the assumption that a search needs a close relative
of whatever it is looking for. Re-running the entire search with only the
three human receptors plus one ryanodine receptor recovers **782 of the 783
genes the full panel recovers.** Removing every bird bait, or every mammal
bait, or every fish bait, costs at most two genes each. What does matter is
having a bait for each paralog: remove one paralog's baits and a quarter of
that gene's copies disappear. A single bait recovers its gene reliably even
when it is only 50 % identical to it. For anyone searching a gene family
across genomes, the lesson is that breadth is nearly free and specificity is
not — spend the effort on covering the paralogs, not on sampling the tree.

**And a caution about a standard safeguard.** Iterative profile searches can
wander out of the family they started in, and the standard guard is to watch
for contamination by the nearest sister family. We ran seven such searches
and three of them wandered badly. **The guard fired on none of them.** The
reason is almost funny: when a search wanders, it picks up so much unrelated
material that the sister family's *share* of what it has collected goes
*down*. The alarm is wired to a quantity that falls during the emergency it
is meant to detect. Watching the share of genuinely unrelated material
instead separates all seven searches perfectly. We report that as a proposal
rather than a change, because seven searches is not enough to prove a rule
*(pending: validation on an independent set of runs)*.
