## 11.10 The gate and the selectivity filter are the least changeable elements

Every element is tested against the same protein's own linkers rather than
against its whole-protein mean, which would contain the element being tested.

**The gate and the selectivity filter are the most constrained elements of the
protein, on both metrics and in all three paralogues.** The RIH-associated
domain is the most constrained domain. Every named element except one sits
above its own protein's linker mean.

![](figures/s17_elements.png)

**{fig:s17_elements}.** Constraint by element, with the composition-free
metric drawn beside the divergence metric rather than instead of it, so a
reader can see the two agree rather than being asked to believe it. A
divergence from a background frequency table [R178] scores a transmembrane
element lower at equal conservation, which is exactly the artefact §11.11 had
to rule out.

## 11.11 One element inside the channel is the exception, and finding it changed the result

The luminal loop scores far below every other element and below the linker
mean, making it by a wide margin the least conserved sequence in the receptor,
on the composition-free metric as well.

**Before it was separated out, the pore domain read as the least constrained
named element in the receptor, below the linkers.** That is a strange thing to
report about the pore of an ion channel, and the sort of result that should
not be printed without asking what else could produce it. Three explanations
were checked and all three were measured rather than argued.

The first candidate explanation was a **metric artefact**, since the
divergence measure is computed against a background frequency table and a
transmembrane domain is built of common residues. It was ruled out, because
the composition-free metric gave the same picture.

The second was a **gene-model artefact**, since the deep alignment is about
95 % genomic translations and the exon-dense transmembrane region is where a
mis-placed boundary would do most damage. It was ruled out, because
recomputing on the curated subset alone still put the domain below the
linkers.

The third was **an unresolved element**, and it was confirmed. The domain's
mean was the average of the most conserved stretch in the protein and the
least, and a bar of that mean is a number no residue has.

With the loop separated, **the pore module is more constrained than the
linkers, and the fifty residues of luminal loop inside it are the most
variable sequence in the receptor.** Those two live about fifty residues apart
in the same Pfam domain.

This is what the geometric definition was for. Had the loop boundary been
drawn on the conservation profile, this section would be circular. Drawn on
the membrane's own axial span in the structure, it is an independent
prediction that landed on the dip.

![](figures/s17_channel_profile.png)

**{fig:s17_channel_profile}.** Constraint along the channel, binned, with the
bins never crossing an element boundary. A sliding window across the luminal
loop's edge would draw the curve straight through the boundary the panel
exists to show.

## 11.12 The ligand contacts and the gate are the two extremes of paralogue divergence

Each site class is tested against two controls: against the whole protein,
and, the sharper one, against the rest of its own element. A gate residue
beating the average residue of a 2,700-residue receptor is nearly guaranteed,
and beating the rest of the pore domain is not.

**The measured IP₃ contacts are more constrained than the rest of the domains
that carry them, in all three paralogues.** That is a statement about ten
residues and is reported as such. The filter and gate sets are two residues
each, and their means are given because the element-level test is where those
elements are actually powered, so a p-value on two residues is not offered.

![](figures/s17_functional_sites.png)

**{fig:s17_functional_sites}.** The measured functional residues against the
whole protein and against the rest of their own element.

A second and completely independent instrument on the same question is
identity between the paralogues, which needs no alignment depth and no
conservation metric at all.

**The gate is identical in all three pairs.** The luminal loop retains 13 to
31 %. The whole protein sits at 0.64 to 0.70.

So the two extremes of paralogue divergence in this receptor are **five
residues that have not changed since the duplications** and **fifty that
retain a tenth to a third of their identity**, and they sit inside the same
domain. Whatever the three copies were free to differ in, it was not the gate.

**One number looks like a disagreement with Chapter 6 and is not.** Chapter 6
reports between-paralogue identity at 0.741 to 0.791 and this chapter reports
0.64 to 0.70. Neither is wrong, because Chapter 6 measured on the trimmed
alignment, which kept 15 % of the columns. So every per-element row here
carries the same pair measured Chapter 6's way as well, and it reproduces that
chapter's committed matrix exactly.

This chapter cannot use the trimmed alignment, and the reason is the result.
Trimming deletes divergent and gappy columns, and the element this section's
headline rests on is the most divergent and gappiest in the receptor: 12 of
the luminal loop's 51 residues survive trimming, against all of the gate and
all of the filter. A per-element identity read off the trimmed alignment would
be a measurement of what survived trimming.

## 11.13 Conservation works as a variant classifier, and the best layer is not the deep one

The harvest returned 1,753 missense records across the three genes, and
**none was dropped.**

Numbering is checked rather than assumed. Every cited transcript's own
translated coding sequence is fetched, aligned to the canonical, positions
transferred through that alignment, and the reference amino acid then required
to match. All three genes file their records on a transcript whose translation
is the canonical, so the check passes with nothing dropped.

That is a result of the check rather than a reason it was unnecessary. The
same code on another channel family found 512 of 773 positions carrying a
different amino acid in the canonical sequence, and a resource built without
the check would have been wrong everywhere. Here it is right everywhere, and
the report can say which.

**The harvest carries a positive control it could fail.** Chapter 2's curated
table localises exactly two variants to a named residue, and both must come
back out of a query returning about 1,750 records. Both are recovered.

**The clinical record is mostly uncertain, and that is a finding.** 1,546 of
1,753 records, or 88 %, are of uncertain significance, and the labelled sets a
classifier can be scored on are 49, 1 and 5 pathogenic positions. One
paralogue's entire clinical record is a single pathogenic missense variant, so
no per-gene classifier score is reported for it.

**The four layers are compared on one fixed set of variants.** They do not
cover the same residues, because the representative-alignment layers lose
gappy columns the deep alignment fills and the reverse, so ranking four scores
measured on four slightly different variant sets would compare the sets as
much as the layers. The comparison is restricted to the 44 pathogenic and 34
benign positions where every layer has a reliable score.

![](figures/s17_variant_classifier.png)

**{fig:s17_variant_classifier}.** Four conservation layers as classifiers,
drawn as full curves rather than as a bar of summary scores, because the
layers cross.

Two things fall out, and one is against this chapter's own design.

The shallow control is the worst layer, so the depth was worth building,
because without the deep ortholog sets that is the number this thesis would
have had.

But **the best layer is the family-wide one rather than the deep
within-paralogue one.** Taxonomic breadth beats within-gene depth here,
because a position conserved across 500 million years of paralogue divergence
and across the eukaryotes discriminates pathogenic from benign better than a
position merely invariant across 260 vertebrate orthologues of the same gene.
A resource for this family should quote the family-wide layer, and the
instrument this chapter was built around is the second best of four.

**The whole-protein control is why the first column cannot be read alone.**
Every layer also separates pathogenic positions from the average residue, but
by less than it separates them from the benign set, so part of the signal is
benign calls being enriched in unconstrained positions rather than only
pathogenic calls being enriched in constrained ones. Both directions are real
and the control is what distinguishes them.

## 11.14 What a constraint resource can offer an uncertain variant

Eighty-eight per cent of the record is uncertain, and a per-site score is what
a resource of this kind can offer such a variant. Each is placed against the
labelled distributions of the same gene, on thresholds that are those
distributions' own medians, so no cut was chosen to make a count.

**This is a stratification rather than a call.** It says where a variant sits
on an axis the labelled variants separate on, which is a different claim from
pathogenicity, and the per-gene numbers inherit the problem above, in that one
paralogue's pathogenic median is a single position's score.

![](figures/supp_variants_on_structure.png)

**{fig:supp_variants_on_structure}.** Every labelled variant with the
per-element enrichment test beside it. An infinite odds ratio is drawn at the
ceiling with a marker rather than allowed off the axes, where a significant
result would simply vanish.

The deliverable is one table per human paralogue, giving every residue in its
own numbering, four conservation layers, the occupancy each is conditional on,
its element, and whether it is a measured contact, filter or gate residue.

## 11.15 Per-site selection rates land on the same coordinates

Chapter 10 answered the whole-gene question and not where the constraint sits.
A per-site rate model [R151] on the same codon alignments gives that, and
every site is carried onto the human reference through a **validated
transfer**, in which a site whose amino acid disagrees after realignment is
written with no residue rather than with a plausible wrong one.

The gate and the filter come back at a non-synonymous rate of zero with 100 %
of sites called significantly constrained, in every paralogue where the test
has power.

Three cautions are carried forward rather than restated. A per-site ratio is
taken only over sites where the synonymous rate is identifiable, with the
excluded count as a column. The constrained fraction is read at the method's
own significance threshold. And the obvious aggregate, meaning the sum of one
rate over the sum of the other, is not usable here, because Chapter 10
measured synonymous saturation on these alignments and the method duly pins
the synonymous rate at its bound on a share of sites, so the sum is dominated
by sites whose denominator is unidentifiable rather than large.

## 11.16 The constraint map is painted onto the experimental structures

Both layers are written into the B-factor column of every IP₃ receptor
structure in the panel, on one scale read the same way round, so they colour
with one command. Sixteen structures are painted at 98 to 100 % of their
resolved residues.

**Unscored residues are written as −1 rather than 0**, and that distinction is
load-bearing for exactly one region. The luminal loop is both the least
conserved element in the receptor and the one the maps resolve worst, and on a
coloured structure those two must not look the same. A chain that maps less
than half of itself to its paralogue is refused rather than painted with
somebody else's profile.

![](figures/supp_constraint_on_channel.png)

**{fig:supp_constraint_on_channel}.** The constraint map painted onto the
channel, drawn from the file the painting step wrote and coloured from its own
B-factor column, so a disagreement with the per-residue tables would be a bug
in the painting step, which is the point of drawing it. Unscored residues are
grey rather than the low end of the scale.

**The panel is carried by experimental structures rather than predicted ones.**
The human ITPR2 prediction is the 181-residue isoform of §11.2 and has no file
to paint, so that paralogue's constraint is shown on a cryo-electron
microscopy entry. A structural resource for this family cannot be built out of
the prediction database.

One further check has to be recorded because it fails on half its candidates.
A structure may carry a human variant position only if every residue it shares
with the human table carries the same amino acid, and two of four candidates
fail, comprising a rat reference and an isoform model. One paralogue's 55
pathogenic positions are therefore not placed on a structure at all, which is
a limit stated rather than worked around.

## 11.17 What this chapter settles about the channel

The gate and the selectivity filter are the least changeable elements of the
receptor, on two metrics, in all three paralogues, and the gate is identical
between all three. The measured IP₃ contacts are more constrained than the
rest of the domains carrying them, and those domains are MIR and RIH rather
than the Pfam signature named after the ligand. The pore module is more
constrained than the linkers once the fifty residues of luminal loop are
separated from it, and that loop is the most variable sequence in the
receptor.

Conservation is a usable variant classifier for this family, at a
discrimination the shallow control does not reach, and the layer that works
best is not the one this chapter was built to produce.

Fourteen constructed negative controls run before anything is written, and one
of them found a real bug on its first run: the within-protein control was
being selected by a prefix test that swept 225 residues of the N-terminal
β-trefoil, which is the element the ligand question is about, into the control
set. The suite was mutation-tested on four deliberate rule breakages and
caught all four.

What it does not settle is what the ligand site is for evolutionarily, which
is the next chapter.
