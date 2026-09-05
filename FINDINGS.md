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

**A caution about our own counting.** The sweep reports 142 "ryanodine
receptor" hits in fungi, green algae and protists, where no ryanodine receptor
is thought to exist. They are proteins of a few hundred to two thousand
residues matching **4–9 %** of a five-thousand-residue channel model — the
module these families share with each other and with unrelated proteins, not
the gene. Counted as genes they would invent a receptor family across half the
eukaryotic tree. Every call in this task is therefore reported with how much of
the model it actually covers: 63 % of the IP₃-receptor calls span at least half
their model, against 1 % of the ryanodine ones.
