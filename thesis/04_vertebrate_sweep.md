# 4. Searching 309 vertebrate genomes, and measuring what that search was worth

## 4.1 Why the reference proteomes were not enough

Chapter 3 ended with two limits. A reference proteome is a gene set rather
than a genome, so a gene missing from one may be missing from an annotation
rather than from the DNA, and fifteen of 763 vertebrate reference proteomes
carried no IP₃ receptor record at all. Neither of those is a biological
statement. Turning them into one requires searching the assemblies.

This chapter does that over a declared scope of 309 vertebrate genomes, and
then does something less usual: it measures what the search was worth. The
second half is not a methods appendix. The numbers it produces are the error
bars on Chapters 5, 9 and 13, and a census that reports zero losses, as
Chapter 9 does, is worth nothing unless somebody has measured how often the
same search fails to find a gene that is demonstrably there. That measurement
belongs with the instrument rather than with the result it qualifies.

## 4.2 The denominator was declared before the search ran

The scope is 309 assemblies, and it is the union of two rules.

The first rule is taxonomic. It takes one best assembly per vertebrate order
[R159], 161 of them, ranked by a stated function that prefers annotated over
unannotated, RefSeq over GenBank, then assembly level, then scaffold N50,
which is the length such that half the assembly sits in pieces at least that
long and is the standard summary of how contiguous an assembly is.

The second rule is a margin rule, and it is computed from the census rather
than hand-listed. Four positive tests put a species in scope: its reference
proteome returned no family record at all, or it carries fewer than three
paralogues, or every record it has is below the family's minimum length, or it
is an anchor species the project needs for another reason. That produces 169
margin species, and the union with the order representatives is 309 genomes
and 552.5 Gbp.

The point of computing the margin set from the census is that the denominator
rebuilds itself. If the census changes, the scope changes with it, and the
scope document is rendered from the manifest rather than written.

That design has a consequence discovered later and worth stating here, because
it shapes every absence claim in the thesis. The margin species were selected
for what their *proteomes* lack. They turn out to be very largely the same
genomes whose *assemblies* cannot hold the gene: 68 % of margin species fall
below the contiguity bar against 12 % of order representatives. The two
signals are confounded by construction, and separating them is what §4.7 and
Chapter 9 exist to do.

## 4.3 Seven rules chose the bait panel, and six slots stayed empty

The sweep aligns proteins to genomes, so it needs baits. There are 38 of them,
comprising 30 IP₃ receptors and 8 ryanodine receptors as an internal positive
control, 121,294 residues in total, derived from the census by seven stated
rules rather than assembled by hand.

The rules are worth listing because each is a positive test. A bait's
paralogue label must come from the census rather than from its own gene
symbol. It must be full length and carry the complete architecture. There is
one bait per paralogue per clade band, and two in the bands that dominate the
scope. Every candidate passes a chimera screen. The ryanodine baits are a
presence control. The existing profile seeds are reused rather than
re-derived. And the panel is scoped to the vertebrate genomes of this sweep
and no other.

**Six slots stay empty, and that is a finding rather than an oversight.** They
are ITPR1 and ITPR2 in cartilaginous fish, ITPR2 in the coelacanth grade, and
all three in cyclostomes. The census holds no labelled, full-length,
complete-architecture record for those combinations at all, because the
paralogue labels simply run out below the well-annotated clades. Rather than
promote an unlabelled record into a labelled slot, three unlabelled baits
cover those genomes without making a paralogue claim, and Chapter 6 settles
the assignment.

**The chimera screen has its own negative control, run on every build.** The
screen rejected none of the 57 candidates, and a screen that cannot fire looks
identical from outside to a shortlist that is clean. So every build constructs
three synthetic failures from real panel baits and requires each to be
rejected by the rule responsible: an IP₃ receptor N-terminus joined to a
ryanodine receptor C-terminus, caught by the family assignment; a mis-joined
gene model with 435 unaligned residues inside its envelope, caught by the
shape rule; and a truncated model, caught by the length band.

That test earned its place on its first run by failing. It originally asked
the envelope rule to reject the truncation, which the envelope cannot and
should not do, because a protein truncated to 40 % of its length aligns 100 %
of itself to the profile in one clean segment. Truncation is the length band's
job and the envelope's job is fusion.

## 4.4 The aligner's intron parameter is a threshold on the call, not a performance knob

The protein-to-genome aligner [R164] takes a maximum intron length. Its
default is 200 kb, and it is not a performance knob: a gene whose largest
intron exceeds it is split into pieces, and a split IP₃ receptor reads out of
this ledger as a fragment. Setting it wrongly manufactures exactly the
observation the thesis is trying to measure.

Chapter 2 supplied the wrong number for this. It measured genomic *spans*,
with human ITPR2 at 498 kb, and a span over 57 exons is not an intron. So the
largest single intron of every IP₃ and ryanodine receptor gene was measured
directly across an eleven-species panel.

**No IP₃ receptor gene in the panel has an intron over the default.** The
widest is 152,216 bp, in human ITPR1. The ryanodine control exceeds 200 kb
twice, the widest being human RYR2 at 227,927 bp. So the default is adequate
for the deliverable and marginal for the control, which is the opposite of
what the span figures suggested.

The rule the sweep uses takes the widest measured intron, doubles it, scales
it by genome size, floors it at the aligner's own default and caps it for
tractability. It is deliberately not the intron-per-gigabase ratio, because
that statistic is largest in the smallest genomes measured, and extrapolating
it linearly asks for a 7.4 Mbp setting on the lungfish, which is an artefact
of a small denominator rather than a measurement. The cap binds on two genomes
and the ledger flags both, because a fragment there could be a real intron the
sweep declined to span.

## 4.5 What the sweep found across 309 genomes

The sweep covered 309 genomes of the declared 309, recorded 2,144 loci, and
carried a ryanodine receptor bait with every one of them.

**The control fired in 309 of 309 genomes.** Ryanodine receptors are present
in three copies in every vertebrate, so a genome where the control finds
nothing has an assembly or a pipeline problem rather than a biological result,
and until it fires that genome's IP₃ receptor cells say nothing. None failed,
so no result below is excluded on control grounds.

![](figures/ledger_by_class.png)

**{fig:ledger_by_class}.** The three-paralogue ledger across 309 vertebrate
genomes, by class. This is the raw ledger before any of Chapter 9's
restatement, so it shows the problem that chapter had to solve. The family
is found nearly everywhere it is looked for, and what red there is falls in
particular vertebrate classes rather than scattered across them. Those
classes turn out to be the ones whose assemblies are worst rather than the
ones whose biology is different, so reading this figure directly as a map of
gene loss is the mistake the rest of the thesis is built to avoid.

![](figures/ledger_status.png)

**{fig:ledger_status}.** Cell status across the sweep, by paralogue. Four
cells of 927 are called absent, meaning no spliced-alignment locus and no
remnant from the rescue search. Those four are the strongest negative the
search itself can produce, and they are still not a loss claim: Chapter 9
shows that every one sits in an assembly that could not hold the gene or in
a genome whose paralogue labels the bait panel cannot resolve. The distance
between this figure and that conclusion is the distance between a search
result and a biological result.

**The two families never contested a locus.** At all 2,144 loci, only one
family's baits aligned at all. That is a far sharper separation than the
protein level afforded, where, before the sister test existed, all six
ryanodine decoys were promoted as candidate family members, and it holds from
the six-genome pilot through to full scale. Spliced alignment against a
labelled panel settles the family before any margin has to be applied.

## 4.6 A threshold inherited from another project was measured and overturned

The rescue step [R158] attributes a fragment to a paralogue when its best bait
beats the runner-up by a stated relative margin. That margin arrived in this
project as an inheritance: 0.333, a 1.5-fold bit-score ratio, from a sister
project on a gene family whose paralogues are 40 to 50 % identical. These are
61 to 68 % identical.

The question is whether the threshold transfers, and it can be answered
without any new search. Loci whose paralogue identity the assembly's own
annotation establishes carry a margin measured on evidence the bait scores did
not produce. Across 542 such loci the margin runs from 0.136 to 0.411 with a
median of 0.314, and **368 of the 542 fall below the inherited threshold.**

A complete locus bounds a fragment's separation from above, so a threshold
that complete evidence fails cannot be met by a fragment. Under the inherited
value every rescue trace in this project would have been reported ambiguous
and no absence claim could ever have been attributed to a paralogue at all.
The threshold is now 0.22, which is the highest value that rejects none of
those correct calls, and its limitation is stated rather than hidden: derived
from complete loci and applied to fragments, it is a ceiling on the right
answer rather than the right answer.

The full sweep then produced the fragments themselves. Fifty-three rescue
regions overlap a gene the assembly names for a paralogue, so their identity
is established independently of the bait scores. The attribution agrees with
the annotation on 53 of 53, and none falls below the threshold in force. That
bounds the threshold from below, meaning it establishes how low the threshold
must be to keep correct calls, and not from above. That limitation is said
plainly because no region was attributed against its annotation, so nothing
measures where wrong calls would start.

**A second inherited constant failed the same way.** Rescue attribution
originally ranked every bait clade against every other, including the
unlabelled baits. Those are the gar, chimaera and lamprey seeds, which are
evidence that a region is an IP₃ receptor rather than a rival hypothesis about
which one. Ranking them as competitors collapsed real margins. In one bird
genome an ITPR1 trace scoring 916.1 against ITPR2's 700.6, a margin of 0.235,
was reported at 0.027 because an unlabelled bait sat between them at 891.6.
Every ITPR1 and ITPR2 trace in that genome was deflated the same way.

## 4.7 Assembly contiguity is why four absences are not four losses

An assembly whose contigs are shorter than the gene cannot carry the gene on
one contig, so an empty cell in such an assembly is evidence about the
assembly.

The bar this project uses is the median measured IP₃ receptor genomic span,
142,212 bp of contig N50, taken from gene geometry with no error rate in its
derivation. **120 of 309 genomes fall below it, which is 39 %, and they are
not a random 39 %.** Assembly quality across the vertebrates is uneven by
clade and by sequencing era, which is the problem the reference genome
consortia were set up to address [R195], and this family's genes are long
enough to feel it. Two thirds of birds fail it against 11 % of
ray-finned fish, and 68 % of margin species against 12 % of order
representatives.

![](figures/contiguity_confound.png)

**{fig:contiguity_confound}.** Recovery against assembly contiguity, binned
on the bar rather than across it. A sliding window straddling the threshold
would report a recovery rate that no genome in the window has, smoothing
away the very discontinuity the panel exists to show. The relationship drawn
here is the confounder every absence claim in this thesis has to survive,
and it is why the contiguity floor is applied before any cell is read as
evidence about biology.

The pilot had predicted, from six genomes, that a fragmented assembly loses
the long paralogues first, because a 231 kb ITPR2 needs a contig that an 82 kb
ITPR3 does not. At 309 genomes that is confirmed: above the bar all three
paralogues are found in 98 to 99 % of genomes, while below it ITPR1 is
recovered in 61 %, ITPR2 in 57 % and ITPR3 in 70 %.

This matters more than it looks. The detection bias runs in the same direction
as the loss signal the margin species were selected for. A per-paralogue
absence in a fragmented assembly cannot be read as a loss, and the confound is
structural rather than incidental.

![](figures/ledger_copy_number.png)

**{fig:ledger_copy_number}.** Copy number per genome across the sweep. The
excess above three copies is not spread across the vertebrates but sits
almost entirely in the ray-finned fish, which is what a whole-genome
duplication confined to one lineage looks like from a copy count. Chapter 7
turns that observation into a tested claim about which paralogue was doubled
and retained.

## 4.8 The three paralogues are not annotated equally well

The sweep also produces something the protein databases cannot: for every gene
it finds, whether an annotation is there and whether that annotation names the
right paralogue. The measurement is restricted to loci the sweep found in an
assembly that carries a gene set, so the denominator is genes that exist and
are annotatable, and it then additionally holds contiguity constant, because
otherwise it measures assembly quality and calls it annotation quality.

Above the bar, a gene model is present at 96 % of found loci for all three
paralogues. It names the paralogue correctly at 65 % for ITPR1, 87 % for ITPR2
and 88 % for ITPR3.

A 23-point gap between genes of near-identical protein length in the same
genomes is not a biological difference. A census built on gene symbols
inherits that gap as an apparent difference in copy number, which is why
Chapter 13 exists.

The uncontrolled version of that table is instructive about the control.
Without holding contiguity constant, ITPR3 appears to be 21 points better
annotated than ITPR2. Above the bar the two are within one point, so the
apparent annotation gap between those two paralogues was assembly quality.

## 4.9 The sweep found 318 genes that exist only as DNA

The sweep contributes 1,058 gene models from 224 genomes, each scored against
both profiles so that a new census row rests on the instrument that called the
old ones. The census becomes 17,097 records.

The interesting column is not how many but **how a database holds each locus**.
318 IP₃ receptor gene models exist only as DNA, meaning there is no gene model
at the locus, or the assembly has no gene set at all. A further 167 sit inside
an annotated gene that carries no family name, so no search by name can reach
them however complete the protein databases are.

Three loci are claimed by a paralogue cell other than the one the assembly's
annotation names, and one of those sits in a genome that also carries a
separate locus for the paralogue the annotation names, so the two cannot both
be that paralogue and the disagreement is internally coherent. They are handed
to Chapter 13 rather than adjudicated here, because one instrument does not
overturn a public annotation and this family's paralogue margins are narrow.
