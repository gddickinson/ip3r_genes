# 11. The channel: shape, constraint and the clinic

## 11.1 Two questions about the protein itself

Everything so far has been about genes: where they are, where they came from,
whether they are still there. This chapter is about the protein.

Two questions. Does what the census calls an IP₃ receptor *fold* like one, and
can it be told from a ryanodine receptor by shape alone? And which parts of
the receptor cannot change — where is the constraint, does it sit where the
structure says the machine works, and does it explain the clinical variants?

The first question has an answer that is mostly about a database, and it is
worth having in advance because it constrains what the second can be built on.

## 11.2 AlphaFold does not hold this family

A vertebrate IP₃ receptor subunit is about 2,700 residues, which is where the
monomer prediction pipeline [R167] stops; the sister family at about 5,000 is past it
outright. So structural coverage is measured before anything is compared, and
it is measured three ways that do not agree: the cross-reference in the
protein record, an actual probe of the prediction database [R168], and *model
coverage* against the length the census holds.

Only the third can see the trap.

**The database is keyed on an accession and answers with whatever record it
holds, and that record may be an isoform.** Asked for the three human
paralogues it returns a 2,695-residue model of a 2,758-residue protein, the
canonical 2,671-residue ITPR3 — and **181 residues of the 2,701-residue
ITPR2.** A probe that took the first record and called it a hit would report
the family's own reference paralogues as fully modelled.

**20.9 % of the census's protein-database records have a usable model**, and
of the 134 representatives every alignment, tree and selection result in this
thesis stands on, **nine** do.

![](figures/s11_afdb_coverage.png)

**{fig:s11_afdb_coverage}.** Coverage by group, and against record length. The
second panel is the result: the usable-model mass sits below about 1,300
residues while the peak at about 2,700 — a full-length subunit — is almost
entirely unmodelled.

The median modelled record is 392 residues and the median unmodelled one
2,674. Of the 5,861 census records at or above the family's own length floor —
the plausibly full-length ones — **13, or 0.2 %, have a usable model.** The
coverage this database offers the family is concentrated on its fragments,
which are the records a structural argument can do least with.

This was written down as a guess before it was measured, which is the only
reason it is worth reporting as a confirmed prediction rather than as an
observation.

## 11.3 Building a panel by rule

The references are resolved by query rather than from memory, and seven rules
choose what goes in, each a positive test.

**The family call comes from this project's census on the entity's accession,
never from an entry title** — the two families share every diagnostic domain
and half their titles. Structures must be cryo-electron microscopy with a
recorded resolution and full length in the family's declared band, so a
ryanodine receptor cannot enter through the IP₃ band because the band was
defined to exclude it. The conformational state must be readable from the
title against a stated vocabulary, ordered so that compound states are tested
before the words they contain. One reference per paralogue, and a paralogue
the census does not name gets **no** label rather than falling back to the
family name. A state panel on the paralogue offering the most conformations,
because a model scored against one conformation confounds fold with gating.
And the negative controls come from Chapter 3's committed decoy panel, one per
class, held to the same rules and size-matched.

One further instrument is available and is deliberately optional: a
structure-based homology search over the whole prediction database [R169],
which is the one stage in the project needing a multi-gigabyte index and is
off by default. What it returns is the shared domain rather than the family,
which is the same statement §11.5 makes with coordinates.

Enumerating candidates on the **union** of the four family signatures is
load-bearing. The structural database's own domain annotation of human ITPR3
— the project's own IP₃ receptor reference — carries four signatures and
**not** the one that names the family, the same absence Chapter 3 measured on
2,911 of 15,417 proteins. A query on that signature alone would have missed
the reference.

**Nothing enters unread.** Every file is parsed, reduced to **one chain** —
both families are homotetramers and an assembly is about 11,000 residues, so
scoring a monomer against a deposited assembly would answer a question about
quaternary structure — and measured. Twenty-nine of thirty structures are
usable, and the one rejection is the 181-residue isoform model, caught by a
rule rather than by inspection.

Where a taxonomic slot has no representative with a model, the best-covered
census record stands in, and the audit records that the slot was *filled*
rather than met. Without that, the panel is six proteins.

![](figures/s11_panel.png)

**{fig:s11_panel}.** Every structure in the panel: the pale bar is the length
the record claims, the filled bar what the structure delivers.

## 11.4 Calibrating the scale before using it

Before any structure is called, the scale is calibrated on this panel rather
than taken from the literature. Five classes of pair, each answering a
different question.

**The negative controls span 0.09 to 0.25 and none of 81 control pairs reaches
the same-fold bar** under either normalisation. That is the floor every
positive result is read against, measured on proteins chosen as decoys for a
*sequence* scorer and re-used here without reselection.

**Conformation is worth 0.22.** Two structures of the same paralogue in
different states score a median 0.78; two different IP₃ receptors score 0.43.
So a fold difference of that size or less is not readable on this panel and is
not claimed — which is what the state panel exists to establish.

![](figures/s11_tm_calibration.png)

**{fig:s11_tm_calibration}.** The family call structurally, and the
calibration behind it. Both of the alignment method's published bars are
**drawn** — the random-similarity floor and the same-fold bar [R166] — rather
than described. A calibration figure that asked to be believed would not be
one.

Both normalisations are kept, because on this panel they say different things:
a subunit resolves to about 2,200 residues and a ryanodine receptor to about
4,300, so **which chain a cross-family score is normalised by decides whether
it reads as a fold result or a size one.**

## 11.5 The family separation, asked of shape

This is the sixth instrument in the thesis to separate the two families, and
the first that reads no gene symbol, no domain architecture, no alignment
score and no tree — only coordinates.

The test is the best score against the IP₃ references, the best against the
ryanodine references, and a **relative** margin between them, with scores
normalised **by the reference**, which is the question being asked: does this
structure account for an IP₃ receptor, or for a ryanodine receptor?

A margin is only consulted once the winner clears the same-fold bar, and that
gate is not decoration: all three negative controls beat their own runner-up
by about 30 % of their score while scoring 0.17 to 0.27 against everything, so
a rule gated on the margin alone calls a dynein heavy chain an IP₃ receptor.

**Twenty of twenty structures the instrument could call are called the way the
census calls them, and none of three negative controls receives a family call
at all.**

The instrument declines six further census-named structures, and a refusal is
not a disagreement. Every one of them still prefers the IP₃ receptors over the
ryanodine receptors by roughly two to one; what they cannot do is account for
enough of a reference to earn a call.

## 11.6 Do the deep records fold like receptors?

This is where the structural work is worth most. The plant, fungal, protist
and non-vertebrate metazoan records were called IP₃ receptor on sequence
evidence alone — profile score, architecture, and a per-record contamination
chase. None of that is shape.

**Two of seven deep models reach the same-fold bar.** The rest fall short.

The interesting column is the **ceiling**. A reference-normalised score
divides by the reference's length, so a model of *n* residues scored against
an *L*-residue reference cannot exceed *n*/*L* before similarity is considered
at all. Most of these models are 1,000 to 1,300 residues — the only records
the prediction database holds for these groups — against references of about
2,050.

**None of them is excluded by arithmetic alone.** Every one had the headroom
to pass and did not, reaching 0.56 to 0.74 of its own ceiling: far above the
negative controls and far below a full-length receptor. That is what a
genuinely divergent homologue modelled at low confidence should look like.

The honest verdict is that the fold test neither confirms nor contradicts
these records. It does not reach them, and the project's vocabulary has a word
for a measurement that did not happen.

## 11.7 Where the models are confident, and where the receptor is

A mean confidence score over a 2,700-residue multi-domain channel averages a
well-predicted β-trefoil with hundreds of residues of linker, and the
resulting number is high enough to look reassuring while saying nothing about
the part any claim rests on. So confidence is reported **per domain**, with
the boundaries transferred by pairwise alignment, and a domain landing on too
little of its reference span is reported unplaced rather than averaged over
whatever aligned.

![](figures/s11_plddt_domains.png)

**{fig:s11_plddt_domains}.** Confidence per domain, ordered along the subunit,
with the prediction method's own confident and very-high bands drawn. Every
model gets an *outside annotated domains* contrast row, without which "the
pore is at 85" has nothing to be high against.

**The IP₃-binding core is the best-modelled domain of the receptor** at a
median of 83.9, against 69.5 outside the annotated domains — and **the pore,
the part the two families genuinely share as a working channel, is the worst
of the named domains** at 71.0. That matters for the next chapter, which is
scoped around the ligand site: whether the prediction is good there decides
whether a structural argument about it can stand at all.

## 11.8 What the constraint map is built on

The rest of this chapter is per-residue conservation, and it needs depth the
representative alignment does not have. Thirteen to nineteen tips per
paralogue is a phylogeny; it is not a per-site estimate.

So one locus per genome per paralogue was taken out of the 309 sweep outputs:
**249 to 265 orthologues of each individual paralogue.** Five rules select
them, and the fifth is one the sweep's own bars cannot do.

Bait coverage is aligned residues over bait length, so a model that covers the
bait *and* carries two thousand extra residues passes everything — which is
exactly what a large intron setting manufactures in the giant genomes. So the
fraction of each sequence's **own** residues landing in reference columns is
measured, and the bar is put in that distribution's own gap with both edges
committed, refusing to fire at all when nothing sits below the lowest curated
record.

**It drops one sequence, and that sequence was inflating one paralogue's
alignment from 3,380 to 5,676 columns.**

![](figures/supp_paralog_alignments.png)

**{fig:supp_paralog_alignments}.** The within-paralogue alignments the
constraint map is computed on, as per-residue occupancy of the human
reference rather than per alignment column — the three alignments have three
widths and no shared coordinate.

Three conventions run through every score. **Conservation is
sequence-weighted** [R171] before any column statistic, and the column metric
is an information-theoretic divergence [R178], because ray-finned
fish are about a third of the genome scope and unweighted, every
teleost-specific residue in a 260-sequence alignment reads as conserved.
**Gaps are missing data, not a 21st state**, and the occupancy each score is
conditional on is reported. And the paralogue of each alignment tip is read
from Chapter 6's committed table rather than re-derived from the tree, because
a third derivation that disagreed would be a silent fork in what this project
means by "ITPR2".

## 11.9 The Pfam name for the ligand site does not contain the ligand site

Joining the domain coordinates to the structural measurements says something
that changes how the rest of this chapter has to be read.

**None of the ten measured IP₃ contacts lies in the Pfam signature named
*Inositol 1,4,5-trisphosphate/ryanodine receptor*.** They sit in MIR and RIH —
the two domains the family shares with the ryanodine receptors.

So the N-terminal β-trefoil is not the ligand site, and the ligand question is
asked throughout of the **measured contacts** and never of the Pfam label. A
report that had taken the label at face value would have measured the
suppressor domain and called it the IP₃-binding core.

The structural elements — the filter, the gate, the ten contacts — are in one
paralogue's numbering and *are* transferred, each with an **anchor test that
can fail**: the filter must arrive on the same motif, the gate on the same
lining residues, and a failed anchor aborts rather than writing a coordinate.
The Pfam elements are measured per accession and are not transferred at all,
so the whole class of transfer error does not arise for them.

One element no annotation carries is located **by geometry**: the luminal
loop, defined as the contiguous run of channel residues whose backbone lies
beyond the membrane on the luminal side of the measured axial span. A boundary
drawn on the conservation profile would be the profile explaining itself.

And residues outside every named element are given **named linkers** rather
than left blank, because a single unassigned bucket would pool the
N-terminus, five inter-domain stretches and the whole C-terminal tail and then
use it as the within-protein control.
