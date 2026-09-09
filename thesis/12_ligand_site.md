# 12. The part the ryanodine receptors do not share

## 12.1 The one asymmetry in the superfamily

Every chapter so far has had to separate two families that share a fold, a
pore, four domains and half their titles. There is one thing they do not
share: the ryanodine receptors do not bind IP₃.

That makes the ligand site the single place where a functional difference
between the families is visible in sequence, and it makes three questions
askable that are not askable anywhere else. Is the IP₃-binding core under
different constraint from the pore? Do the ten residues that touch the ligand
behave differently from the rest of the site? And in lineages that have lost
the enzyme that *makes* IP₃, does the receptor's binding site relax?

The answers are, in order: yes and in the opposite direction to the
expectation; no, and the reason is informative; and no, with the null bounded.

## 12.2 Defining the modules, which Pfam cannot do

Neither module is a Pfam domain, and that is the whole reason this stage
exists.

Pfam calls the ligand core *MIR* and *RIH_N* and knows nothing about the
ligand — Chapter 11 established that none of the ten measured contacts lies in
the signature named after IP₃. And Pfam's pore domain contains fifty residues
of luminal loop that Chapter 11 measured as the least conserved sequence in
the receptor.

So the **primary ligand core** is the minimal contiguous span containing every
measured IP₃ contact, and the **primary pore module** is the Pfam pore domain
*less* the luminal loop. A sensitivity definition is committed beside each and
every test is run under all four combinations: the published binding core
transferred through the same routine, and the pore domain as the database
draws it.

Every definition is checked against something it does not contain — the core
must hold all ten contacts and no pore site, the pore both filter residues and
both gate residues and no contact — and a failure raises rather than being
written.

![](figures/supp_labelled_positions.png)

**{fig:supp_labelled_positions}.** The two modules at residue resolution
across all three paralogues, with every pathogenic position's residue
**printed**. A residue-level panel exists so that a reader can check the
claim, and every letter in it has been through the join guard described in
Chapter 14: each variant's reference amino acid must be the residue its own
paralogue's per-residue table holds there, and each aligned partner the
residue the other paralogue's table holds.

## 12.3 Distance to the ligand, not a contact label

A binary contact-or-not label discards the one thing a structure can say that
an alignment cannot. So every residue within 15 Å of the bound ligand was
measured, all-atom, in **six independent depositions**: the structure used
throughout this thesis and five further IP₃-bound entries from other groups
and other gating states.

The reader is deliberately small and all-atom. Chapter 11's structural reader
is backbone-only, and a backbone trace puts an arginine 8 Å from a phosphate
its side chain is hydrogen-bonded to.

**The ten contacts are the positive control and a hard failure**: the reader
must recover all ten in the reference structure or the stage raises. A residue
past the search radius is absent from the table rather than pooled into an
open bin, which would be a population defined by the geometry of the search.

## 12.4 The pairing, and what it controls

The two modules sit in the same protein, so every comparison is **one core
number and one pore number per orthologue**, and a difference between them
cannot be a difference in taxon sampling.

A tip must resolve half of **both** modules to enter, and the dropped ones name
the module that lost them — without which an N-terminally truncated gene model
enters as an extreme ligand-core divergence.

Forty-four constructed negative controls run before anything is written, split
across two suites with one entry point so the driver cannot run half of them,
and the suite was mutation-tested on ten deliberate rule breakages, all ten
caught. Four of this chapter's rules return a number that looks *better* when
they are wrong: a core drawn without checking it holds the contacts still
gives a clean p-value; a paired test admitting a tip that covers the pore and
not the core reports a truncated gene model as ligand-core divergence; a
permutation that resamples the whole module instead of the label answers a
different question significantly; and a panel filing "no reference proteome"
as "no enzyme" manufactures its own test set.

## 12.5 The answer, and it is not the one the question expects

**The pore module is more conserved than the ligand core, not less.**

Paired per orthologue inside each paralogue's own alignment, across roughly
250 orthologues each: ITPR1 puts the pore ahead by 0.024 identity units with
the tips running seven to one; ITPR3 by 0.020 with the tips running five to
one; ITPR2 shows no difference at all.

The module that gives this family its name and its defining signature is the
*less* constrained of the two across a vertebrate orthologue set.

## 12.6 And it reverses on a boundary nobody states

**All three paralogues flip sign** when the luminal loop is left inside the
pore module — which is how the domain database draws it, and therefore how the
comparison would be made by default.

Include the loop and the ligand core wins by about five identity points at
overwhelming significance in all three paralogues. Exclude it and the pore
wins in two and ties in the third.

![](figures/s22_modules.png)

**{fig:s22_modules}.** The core against the pore under **both** pore
definitions on one axis with zero marked, because the answer reverses between
them and a figure showing only the primary would assert the choice instead of
showing what it costs.

Neither answer is wrong about its own module. They are answers about different
modules, and the difference between them is one boundary that a comparison
drawn on database spans would never have to declare. **Any version of this
comparison that does not state the boundary is not interpretable.**

**A second disagreement is reported rather than dropped.** At column level the
divergence metric sees no difference between the modules while the
per-orthologue identity does. That is not a contradiction: the divergence
metric is measured against a background amino-acid table, so a transmembrane
module scores lower than a soluble one at equal conservation. The
composition-free column metric agrees with the paired test. The null result is
a property of the metric, and it is printed here because a reader coming from
Chapter 11's tables would otherwise find two of this project's own numbers
pointing different ways with no explanation.

## 12.7 The constrained unit is the pocket, not the contacts

**The ten residues that touch IP₃ are more constrained than the rest of the
binding core, and not more constrained than the rest of the pocket.**

Against the core the contacts win in all three paralogues. Against every
*other* residue the structure places within 15 Å of the ligand they win in
none.

The permutation resamples **which positions carry the contact label**, keeping
the constraint values where they are, because the hypothesis is about these
particular residues.

![](figures/s22_shells.png)

**{fig:s22_shells}.** Constraint against distance from the ligand, drawn as a
scatter rather than as a bar of shells: the claim is the *absence* of a step
at the contact radius, and four bars cannot show an absence.

Every shell out to 15 Å sits above the whole-protein mean, and the step a
contact-driven model predicts at the contact radius is not there. The gradient
across the neighbourhood is real and shallow, and reaches significance in one
paralogue.

**The result replicates on an independent axis.** Every contact site in every
paralogue is under detectable purifying selection by the per-site rate model,
and the share of significantly constrained sites falls with distance from the
ligand. It is not monotone: the outermost shell sits a little above the third,
so what the data show is a step down from the ligand's first two shells onto a
floor, rather than a gradient running all the way out — and it is reported as
such rather than as a gradient.

Two instruments, two alignments, two different statistics, the same pocket.

![](figures/s22_omega.png)

**{fig:s22_omega}.** Per-site rates by module and by shell.

**The module comparison does not replicate on the rate axis, and in one
paralogue it points the other way.** Both are computed correctly and they are
not the same measurement — 250 sweep orthologues of one paralogue against 57
vertebrate tips, a column's dispersion against a rate on a tree — but a reader
should have the disagreement in front of them rather than only the half that
agrees. With a rate that is exactly zero at most sites and a quarter of the
sequences, this is the underpowered version of the same comparison, and the
module result is read off the identity axis.

## 12.8 The upstream pathway, measured rather than assumed

The third question needs a list of lineages whose IP₃-producing pathway is
reduced or absent. A list taken from reading would be an assumption dressed as
a scope, so it is derived from the same reference proteomes the family was
swept over.

**What counts as a phosphoinositide-specific phospholipase C**: a protein
carrying *both* halves of the catalytic barrel in one sequence. Either half
alone is not one — one of them in particular turns up in unrelated proteins,
and counting it would report a pathway where there is none. Both halves are
counted separately, so a proteome scoring one and not the other is visible as
such.

**The vertebrate group is swept as the positive control for the search**: 760
of 763 vertebrate reference proteomes carry the enzyme. A vertebrate proteome
scoring zero would be an instrument failure, not a result.

**Sixty-four proteomes in the whole sweep carry an IP₃ receptor and no
phosphoinositide-specific phospholipase C.** They are not scattered: the
oomycetes, the early-diverging fungi, three *Perkinsus* species, two ciliates,
four prasinophyte algae and a group of flatworms. Two clades account for more
than half.

**The internal control for that list costs nothing and is strong.** Every one
of those proteomes carries a full-length IP₃ receptor — a 2,700-residue
multi-exon gene the sweep found in it. A gene set complete enough to hold this
receptor is complete enough to hold a phospholipase C.

A separate stratum records proteomes with no reference set at all, because
filing *we could not look* as *there is none* would manufacture the test set
out of missing data.

## 12.9 The lineage test, and why the pooled answer is wrong

Asked of the representative alignment, the test has two tips on one side. That
alignment was built to span four kingdoms with 134 sequences, and these taxa
are not the taxa a diversity rule picks.

**The positive control fires there.** The ryanodine receptors — same pore,
same diagnostic domains, no IP₃ site — move the paired difference by about a
tenth of an identity unit on six tips.

So the test is repeated on the **whole population the question is asked of**:
all 662 non-vertebrate reference proteomes carrying a receptor, each aligned
to the human reference through the same module residues, of which 486 resolve
both modules and enter.

**Pooled, the enzyme-absent records have a wider gap, and that number should
not be believed.** The covariate table says why: enzyme-absent proteomes sit
at a median pore identity of 0.358 against 0.589 for the rest. They are
oomycetes and early-diverging fungi, far further from the human reference than
the arthropods and nematodes that dominate the other cell — and the paired
difference is itself correlated with divergence. A pooled comparison of absent
against present is largely a comparison of distant against near.

Inside a phylum, only three strata have both cells populated at all. And the
one that cannot be tested is the most informative: **every early-diverging
fungal proteome in this set that carries a receptor lacks the enzyme**, so the
phylum has no internal control.

## 12.10 Matched on divergence, the effect is gone — and bounded

Each enzyme-absent record is matched to up to three enzyme-present records
within a stated distance of its own pore identity and scored against their
median. That is the identity-matched paired design of Chapter 9 applied to the
same kind of claim, and it removes the confound rather than arguing with it.

All 36 matched. **The median within-pair difference is −0.0064, with a
confidence interval from −0.0161 to 0.0152, the tips splitting 19 to 17, and
p = 0.87.**

![](figures/s22_lineage.png)

**{fig:s22_lineage}.** The lineage strata with the power curve beside them,
because a lineage test that finds nothing is only readable next to what it
could have found. The ligand core is drawn in a neutral dark rather than a
hue, because the palette reserves its accent colour for the ryanodine
receptors and this figure draws them.

**And the null is bounded rather than empty.** Measured through the *same*
pairwise instrument at the *same* divergence as the test group, the ryanodine
receptors' paired difference is a shift of −0.062. The matched test has
essentially full power at that shift and 59 % at half of it, and its
confidence interval excludes anything larger than about 0.016.

**Losing the enzyme that makes IP₃ does not relax the receptor's IP₃-binding
core.** Not *we could not tell*: the test can see a ligand-free core when
there is one, at the same evolutionary distance, and it does not see one here.

## 12.11 What this settles, and what it does not

**The ligand core is not the most constrained part of this receptor.** Paired
per orthologue against roughly 250 sequences per paralogue, the pore module is
ahead in two and level in the third — and that answer depends on one boundary
and reverses when fifty residues of luminal loop are counted as pore.

**The constrained unit at the ligand site is a pocket about 15 Å across, not a
contact set.** The ten measured contacts beat the rest of the core and do not
beat the rest of the pocket, and the same shape appears on the rate axis.

**Losing the upstream enzyme does not relax the site**, with the null bounded
by a positive control measured on the same instrument.

**What is not settled** is whether the pocket's constraint is *about* IP₃ at
all. Everything within 15 Å of the ligand is also within the fold that holds
it, and no measurement here separates ligand binding from domain packing; a
mutational or binding dataset would, and sequence will not. Nor what these
lineages' receptors are gated by — a search for one enzyme family says the
canonical route to IP₃ is missing, not that the receptor has no ligand, and
the bounded null is consistent with the site being held by something else.

And there is one prior this chapter does not corroborate, which is worth
stating plainly. The signature that names the family is the IP₃-binding core,
and it is the one diagnostic domain the ryanodine receptors do not carry
*functionally*. It is also the less constrained of the two modules measured
here. **Naming a family after its diagnostic domain is a statement about what
the domain identifies, not about what selection holds most tightly**, and this
is the first measurement in the project to ask the second question.
