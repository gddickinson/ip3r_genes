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
