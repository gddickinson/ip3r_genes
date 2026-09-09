# 3. Two families, one architecture

## 3.1 The problem stated precisely

Chapter 2 measured the hazard: 49 % of the records returned by a query on the
IP₃-binding-core signature in one well-annotated genome are ryanodine
receptors. This chapter builds the instrument that separates them, and it does
so in three passes, because no single test reaches the whole family.

The requirement is stated as a rule rather than a preference. Every record in
this project is assigned to one family or the other by a **positive test** —
evidence that it *is* one thing, not evidence that it is not the other — and
where no positive test reaches the required margin, the record is left
uncalled and counted as uncalled. The gene's own name is never an input to
that test at any stage. A length filter is never the call.

Three instruments follow, in the order they were built: a labelled-bait
identity margin, benchmarked against 56 controls; a domain-architecture rule,
applied to an exhaustive enumeration of the search space and audited against
gene symbols it never sees; and a pair of profile hidden Markov models, one
per family, calibrated against the architecture rule before being used.

## 3.2 A benchmark before an instrument

Nothing was searched until the discovery scorer had been run against a panel
of known answers: 25 positive controls — the three vertebrate paralogues
across a species panel plus the invertebrate and non-metazoan single-receptor
grade — and 31 decoys.

The decoys were chosen to be hard rather than plausible. Six are the human and
mouse ryanodine receptors: the sharp decoy, carrying all four
family-diagnostic signatures. Two are the O-mannosyltransferases POMT1 and
POMT2, which share the MIR domain. Fifteen are in-band channels and
non-channels — proteins of the right size and the wrong identity, from
CACNA1A to talin. Four are out-of-band giants.

Recall was 24 of 25, and specificity 31 of 31.

**Specificity was not 31 of 31 on the first run.** It was 25 of 31, and every
one of the six failures was a ryanodine receptor. Each carried all four family
signatures, each therefore took the domain component, which also satisfies the
evidence gate that stops size and novelty alone promoting a candidate, and
each scored 45 against a threshold of 40. The instrument, run as designed,
called the sister family the family.

What fixed it is the labelled-bait margin, and its three properties are the
reason it is used unchanged for the rest of the project. It compares a
candidate's identity to the nearest *labelled* IP₃ receptor bait against its
identity to the nearest *labelled* ryanodine receptor bait, and assigns to
whichever wins by more than a stated margin. It is a positive test on
distances, so an unnamed 4,900-residue locus is called exactly as a named one
is. The length band contributes only to the size component and never to the
call. And the margin is recorded on every candidate, so every exclusion can be
audited afterwards.

![](figures/census_margin.png)

**{fig:census_margin}.** The labelled-bait margin measured across the census,
inside and outside the band in which the project declines to call. The two
families do not overlap: the narrowest true positive sits at +0.065 and the
narrowest ryanodine receptor at −0.536, a gap of 0.601 with nothing in it.

That gap is comfortable, and one true positive sits inside the no-call band
anyway. *Dictyostelium* iplA, a characterised IP₃-gated channel, has a margin
of +0.065 — the right sign, and smaller than the 0.10 the rule requires for a
call. The margin is not wrong about it; the margin simply cannot reach it.
That is a scope limit on this instrument, measured on the first day, and it is
why the profile route in §3.5 is not optional.

**Why the obvious identity metric understated the risk.** A ryanodine receptor
is about 1.8 times the length of an IP₃ receptor, so identity computed over
the full alignment dilutes every cross-family comparison toward zero: the
ryanodine decoys sit at 0.105 to 0.110 identity to the nearest IP₃ receptor
bait, *below* the floor at which the scorer's novelty component fires. Under a
fragment-aware metric that scores only mutually covered columns, the same
comparisons rise to 0.249 to 0.258 — inside the twilight zone, worth a further
20 points. A later switch to the better metric, which the one recall failure
argues for, would have promoted every ryanodine receptor to 65 without the
sister test. Both columns are committed, which is how that is known rather
than suspected.

**The single miss is worth its line.** *Drosophila* Itpr scored 35. It needs
one non-known sibling at 0.35 identity or better; its nearest is the worm
receptor at 0.342 under the metric the scorer uses — short by 0.008. Under the
fragment-aware metric the same pair scores 0.420 and the component fires. The
honest worst case is therefore measured rather than assumed: **a lone
true family member, 40 to 60 % diverged, with no sibling in the set, scores
35 and is missed.** In the census that follows the invertebrate grade is
represented by many lineages that corroborate one another, so the artefact
does not recur — but it is the shape of this instrument's failure and it is on
the record.

## 3.3 Enumerating the space, and not trusting the count

The census was built by walking three Pfam signatures [R155] to exhaustion,
through the domain database that serves them [R156]:
PF08709 (the IP₃-binding core), PF01365 (the RyR–IP₃R homology domain) and
PF08454 (RIH-associated). The MIR domain PF02815 is deliberately *not* a seed
— it is carried by the O-mannosyltransferases as well as by both receptor
families, so it widens the space without adding evidence — and is kept as an
annotation column instead.

Three things about that walk are worth a thesis's space.

**InterPro's own record count is not a completeness criterion, and it is wrong
in both directions.** PF08709 advertises 12,339 records and serves 12,507, an
excess of 168. PF01365 advertises 13,177 and serves 13,233. PF08454 advertises
12,066 and serves 11,890 — 176 fewer than advertised. The advertised values
were stable across re-queries, so this is not an edit during the walk; and in
every case what the interface *serves* matches UniProt's independent count for
the same signature to within two records. It is the advertised number that is
wrong.

Both directions bite. Stopping at an understated count drops records in
silence. Trusting an overstated one declares a complete walk incomplete —
which is what happened: this task's own completeness test failed PF08454 and
exited non-zero on a walk that had run its cursor chain to the end. The
recorded criterion is therefore cursor exhaustion, with both counts kept
beside it.

**The documented interface was down for the duration.** The database answers the
same queries on two hosts: the documented one and the one its own website
uses. During this work the documented host returned HTTP 500 to 11 of 12 probe
requests while the other served 12 of 12, with identical payloads and the same
counts. The client now tries both before it sleeps, which is why the walk
completed at all.

**One signature would not have been enough.** A census built on PF08709 alone
— the signature that names the family — would have enumerated 12,507 of the
15,421 proteins here and missed 2,914, including 758 that this census calls
IP₃ receptor across 385 taxa. Only 10,256 records, 66.5 %, carry all three
seeds. *Dictyostelium* iplA is among the records a single-signature census
would have missed: it carries PF01365 and PF08454 and none of PF08709,
PF02815 or PF00520 — a characterised receptor that the family's defining
domain does not annotate.

![](figures/census_space.png)

**{fig:census_space}.** The enumerated search space by signature, with the
overlaps drawn. The union is 15,421 proteins across 1,488 taxa; the
intersection of all three seeds is 10,256. The difference between those two
numbers is the argument for enumerating a union.

## 3.4 The architecture call, and an audit that could have failed

Every record is called by domain architecture. A record is a **ryanodine
receptor** if it carries any of four RyR-specific signatures: PF02026,
PF06459, PF21119 or PF00622. It is an **IP₃ receptor** if it carries the
complete five-signature architecture — PF08709, PF01365, PF08454, PF02815 and
PF00520 — and none of those four. Length is recorded on every row and enters
only as support for a medium-confidence call. The size band was chosen to
exclude ryanodine receptors, so letting it decide would beg the question the
call exists to answer.

The result is 6,433 IP₃ receptors, 6,807 ryanodine receptors and 2,181 records
that satisfy neither positive test.

![](figures/census_lineage.png)

**{fig:census_lineage}.** The census by lineage with both family calls, and
the sister family drawn beside the family everywhere. The vertebrate bar
carries 10,936 of the 15,421 records — which is a fact about sequencing effort
and is why Chapter 5 counts taxa and proteomes rather than records when it
asks about range.

![](figures/census_lengths.png)

**{fig:census_lengths}.** Length distributions of the two families as called.
The medians are 2,671 and 4,856 residues and the distributions barely touch.
This figure is the reason the length band is a good *filter* and the reason it
is not used as the call: it separates the two families beautifully in a set
whose annotations are already good, which is precisely the population in which
no separation is needed.

**The audit.** The architecture rule never sees a gene symbol, which makes
symbols an independent label to score it against — and every symbol-labelled
record in the census is used, not a sample. The rule agrees with the symbol on
2,479 of 2,479 decided IP₃ receptor cases and 2,492 of 2,492 decided ryanodine
receptor cases, disagreeing on none. Four subsidiary claims about individual
signatures are scored the same way and each is 100 % correct across two to
three thousand records.

One of those subsidiary claims deserves suspicion and gets it. The rule leans
on SPRY (PF00622) being diagnostic of ryanodine receptors, and SPRY sits in
roughly 114,000 UniProt proteins and is in no sense RyR-specific. The claim
being made is conditional — *among proteins that already carry a family seed
signature*, SPRY is diagnostic — and the audit tests exactly that conditional
claim, on 2,335 records, and it holds. It is stated here as the row of the
audit table a reader should check first, because it is the one that is a claim
about this search space rather than about a domain.

**A note on the 1,220 calls that rest on a gene symbol.** Some records have an
architecture too partial to decide and a name that is not. Those are called
from the name at explicitly low confidence, with the reason written into every
affected row and the pure architecture call kept in its own column so the
audit above can score the rule with no symbols anywhere in its input. Section
3.5 revisits all of them with an instrument that reads residues, and overturns
one.

**What the unassigned pile is.** The 2,181 records that satisfy neither test
are overwhelmingly fragments: median 678 residues against 2,671 for a called
receptor, with only 143 of them inside the size band at all. A rule that reads
the absence of a domain as evidence cannot call a partial annotation, and
refusing to is the intended behaviour rather than a shortfall.

**What the size band caught instead.** Because length never decides a call,
the band was free to catch something else: seven records carry the complete
five-signature architecture inside 2,000 residues, from 1,528 to 1,993, against
2,000 for the shortest real family member. An intact architecture in
two-thirds of the length is a truncated gene model. None of the seven is
flagged as a fragment by UniProt — a truncated model submitted as a whole
protein is not marked as one — and all seven are unnamed locus tags from three
species. They are the first evidence in this project for the annotation
failures Chapter 13 audits.

## 3.5 Two profiles, calibrated before use

The architecture rule reads annotation. A profile reads residues. Building
both and requiring them to agree is what makes the census defensible, because
they can fail in different ways.

Two profile hidden Markov models [R153] were built — one of 2,684 match
states from 34 seeds and one of 4,930 from 22 — each seed drawn from the census and carrying that census's
call, so the profiles are labelled by a rule rather than by a gene name. Five
selection rules are enforced in code rather than trusted: a seed whose census
call, architecture count or length has drifted since the manifest was written
aborts the build, and on the first run it did.

Five seeds carry an incomplete architecture and are in the set deliberately,
each with its reason written into the manifest row. The most important is
*Dictyostelium* iplA, at two signatures of five — put into the IP₃ receptor
seed set precisely so that the profile can find its relatives, since the
architecture rule cannot.

Two build decisions are measurements rather than settings. MAFFT is run
**single-threaded**, because at automatic thread count it is not reproducible:
the same ryanodine seeds aligned twice gave 8,510 and 8,468 columns and
profiles of 4,933 and 4,908 match states. And MAFFT's return code is checked
and a ragged alignment aborts the build, because the benchmark in §3.2
established that a silent MAFFT failure degrades to a star alignment with no
other symptom at all. The SHA-256 of every seed set, alignment and profile is
recorded, so a rebuild that drifts is visible in the data rather than only in
a count.

**The instrument is calibrated before it is used.** Both profiles were run
over the whole archived search space — 15,381 records whose call came from a
rule that read annotation rather than residues — and the two were scored
against each other. With the seed sequences removed, because they are in that
set by construction, the profile assignment agrees with the architecture call
on 11,875 records and disagrees on 1.

![](figures/profile_separation_vert.png)

**{fig:profile_separation_vert}.** The two profiles' scores against each
other, with the band inside which no call is made drawn rather than described.
A call requires the winning profile to clear 30 bits, to span at least 200
match states, and to beat the loser by more than 10 % of its own score. The
margin is **relative** and not absolute because bit scores scale with
alignable length: a fixed gap would call every full-length protein confidently
and no fragment at all.

The span floor is measured rather than tuned. In this project's own structural
domain coordinates the shortest observed PF08709 — the IP₃-binding core, which
names the family — is 200 residues, and the longest observed SPRY is 137, so
the floor sits in the gap between them. What it costs is on the record: 89
records that the architecture rule could call lose their profile verdict to
it.

And the second instrument earns its place exactly where the first fails.
**2,314 of the 3,360 records the architecture rule leaves unassigned get a
call from the profiles**, because a partial architecture defeats a rule that
reads absence as evidence and does not defeat a sequence profile. Of the 1,220 calls that had rested
on a gene symbol, the profiles scored 1,219 and overturned exactly one. The
symbol fallback was sound, and it is now checked rather than assumed.

![](figures/instrument_agreement.png)

**{fig:instrument_agreement}.** Which instrument calls each record. The census
keeps the two verdicts side by side and merges them by a stated rule: both
agreeing gives a high-confidence call, one speaking gives that one's, a
disagreement is kept and reported as a conflict rather than resolved by
preference, and neither speaking leaves the record unassigned. Across 16,039
records there are two conflicts.

## 3.6 The sweep, and the gate that had to exist

The profiles were then run over 763 vertebrate reference proteomes [R157]:
14,414,821 canonical proteins, one per gene, because isoform sets add isoforms
of genes already counted and this census counts genes.

The union of the two searches is 18,501 targets. Before that number means
anything, **12,609 hits are declined as module-only matches.** Both profiles
are full-length channel models, so a protein sharing one small domain with
either scores against it — and the ryanodine profile carries SPRY, already
flagged as sitting in roughly 114,000 proteins. Without the span gate this
sweep called 14,981 vertebrate proteins ryanodine receptors, headed by
troponin, DEAD-box helicases and calcium-binding proteins. That is what a
sensitive profile search returns if nothing stops it, and it is the single
clearest demonstration in this project that a threshold on score alone is not
an instrument.

Of the declined hits, 1,794 carry a gene symbol from this family or its
sister. They are not domain-sharing proteins that happen to score: they are
pieces of split or truncated gene models, sitting under the gate because there
is not enough of the protein left to span a family domain. They are kept as an
annotation lead rather than discarded, and Chapter 13 takes them up.

**618 sweep hits are absent from the enumeration entirely** — proteins in a
vertebrate reference proteome that the exhaustive InterPro walk never
returned, 188 of them called IP₃ receptor. A signature-based enumeration and a
profile sweep of the same organisms do not return the same set, and the
difference is not small.

## 3.7 Convergence, and a kill criterion that had to be rewritten

The completeness argument needs a third line: iterate a model from a single
sequence until it stops finding new things, and see whether the family it
converges on is the one the census holds.

Three iterative searches [R154] were run to convergence, from a human IP₃R1, a
*Drosophila* Itpr and an *Acanthamoeba* receptor — a vertebrate, an insect and
an amoebozoan. Only the fly run converged. The other two reached the ten-round
ceiling and were killed.

Killing them is the point. The failure mode for this family is specific: an
IP₃-seeded iterative search drifts across the shared architecture into the
ryanodine receptors and reports a converged, confident, wrong answer. So the
criterion is coded and evaluated per round on that round's own inclusion list,
rather than judged by eye afterwards.

**The first version of that criterion was wrong, and how it was wrong is a
result.** It was coded as a flat ceiling of 5 % sister-family content, and it
fired on every run at round one — because a single IP₃ receptor sequence
searched at a sensitive threshold already returns 32 % to 37 % ryanodine
receptors before any iteration has happened. That is not contamination. The
two families are genuine homologues sharing the entire pore, and a search
sensitive enough to reach *Acanthamoeba* is necessarily sensitive enough to
reach RYR1. The *level* is a fact about shared ancestry; only the *change* in
it can be attributed to iterating the model. The rule measures the rise.

![](figures/jackhmmer_convergence.png)

**{fig:jackhmmer_convergence}.** Convergence for each seed, with the
sister-family trace beside it. Two runs never converge and are excluded from
any completeness claim; the third does. The curves decay to an asymptote of a
few new targets a round rather than to zero, and §3.8 says what that tail is
made of.

## 3.8 What an iterated model is actually built from

A convergence curve says a search has stopped finding things. It does not say
what it found. Every target supporting each run's final model was therefore
classified by the sweep's own two-profile verdict.

The result is the strongest completeness statement in this chapter, and it
does not depend on convergence at all. **The family core is the same whichever
seed finds it.** Three starting points as far apart as a human, a fly and an
amoeba produce final models carrying identical family-called content: 1,656
architecture-complete IP₃ receptors in all three, 1,565 to 1,569
architecture-complete ryanodine receptors, 952 partial IP₃ receptors in all
three.

What differs between them is everything that is *not* the family. The
*Acanthamoeba* run's final model rests on 26,266 targets against the human
run's 7,222, and 18,670 of those are proteins that neither profile scores at
all. And the rise-based kill criterion **cannot see that drift**: off-family
accretion *dilutes* the sister-family share rather than raising it — across
that run it falls from 36.6 % to 8.9 % while the model grows five-fold — so
the rule that catches the failure it was designed for reads this one
backwards. The ten-round ceiling is what caught it, and the composition table
is why the run is reported rather than quietly dropped. Chapter 4 returns to
this and scores all three rules as classifiers of an outcome measured on the
finished model.

## 3.9 Completeness, asked in the expensive direction

A sweep is a completeness argument only if it is also checked the other way
round: how many records that the enumeration called IP₃ receptor, in a
proteome that was swept, did the profiles fail to find?

2,787. Each was then looked up in the 8.8-gigabyte search database itself, in
one streaming pass, rather than explained away by its gene name. Of those,
1,987 were absent from the database though another entry of the same gene was
hit, and 800 were absent and never searched.

**Not one was in the database and missed.** Every apparent miss is a UniProtKB
entry that its taxon's reference proteome does not contain — a reference
proteome holds one canonical protein per gene, and the enumeration walked all
of UniProtKB — so it was never searched at all. The profiles' sensitivity
failures over 763 vertebrate reference proteomes number zero.

Fifteen of those 763 proteomes carry no IP₃ receptor record. That is a list of
leads for the genome sweep, not a list of losses: a reference proteome is an
annotation, and Chapter 4 tests all fifteen against the genomes themselves.

![](figures/proteome_copy_number.png)

**{fig:proteome_copy_number}.** Copy number per proteome against annotation
depth. The relationship is the warning label on every proteome-level count in
this thesis: a gene set with fewer genes has fewer of these genes, and a
census of annotations measures annotation as much as biology. Chapter 4 is the
answer to it.

## 3.10 Census v3, and what it does not settle

The merged census is 16,039 records: 8,000 IP₃ receptors, 7,432 ryanodine
receptors, 605 unassigned and 2 conflicts. Both instruments speak for 79.0 %
of records; the profiles alone for 13.7 %; the architecture rule alone for
3.6 %.

Three limits carry forward, and each becomes a later chapter.

The sweep is vertebrates only. Every count here is a statement about 763
vertebrate reference proteomes, and nothing in it bears on the plant, fungal
or protist questions the enumeration opened. That is Chapter 5.

A reference proteome is a gene set, not a genome. A gene missing from one may
be missing from the annotation rather than from the DNA. That is Chapter 4.

And profile assignment inherits its seeds' breadth. The deep branches are
represented by seven non-metazoan seeds, and a lineage more diverged than any
of them is outside what this instrument has been *shown* to call. Chapter 5
builds a second panel by different rules for exactly that reason.
