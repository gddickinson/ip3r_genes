# 1. The IP₃ receptor, and the questions this thesis answers

## 1.1 The receptor is a channel that reads a second messenger

Almost every eukaryotic cell keeps a store of calcium inside its endoplasmic
reticulum at a concentration three or four orders of magnitude above the
cytosol, and almost every eukaryotic cell has a way of opening that store on
demand. In animals the principal way is the inositol 1,4,5-trisphosphate
receptor, a ligand-gated channel in the endoplasmic reticulum membrane that
opens when it binds the soluble second messenger IP₃ and releases calcium into
the cytosol.

The pathway that produces the ligand is one of the best-characterised in cell
biology. A receptor at the plasma membrane activates a phospholipase C, which
cleaves phosphatidylinositol 4,5-bisphosphate into diacylglycerol and IP₃, and
the IP₃ then diffuses to the endoplasmic reticulum and opens the receptor
({fig:signal_hierarchy}). What follows, meaning the amplitude, the timing and
the spatial pattern of the calcium signal, is largely set by the channel
itself. That is why the receptor has been studied for forty years as a
signalling device rather than merely as a pore.

![](figures/signal_hierarchy.png)

**{fig:signal_hierarchy}.** The receptor's place in the phosphoinositide
pathway, drawn as a hierarchy of scales running from the messenger to the
cellular output. The layer this thesis is about is a single one of these
boxes. The figure is here to show what the rest of the pathway is doing while
the channel opens, because two of the results in later chapters concern the
pathway rather than the channel: the loss of the family in whole eukaryotic
kingdoms, and the enzyme repertoire that survives beside it.

The receptor was identified in the decade after IP₃ was shown to release
calcium from a non-mitochondrial store in pancreatic acinar cells [R01].
Purification and reconstitution demonstrated that a single protein was
sufficient for the flux [R03], and the primary structure, obtained from the
cerebellar P400 protein, revealed a very large polypeptide of about 2,700
residues [R02, R21]. The title of one of those first reports named the problem
this thesis had to solve before it could measure anything: *Putative receptor
for inositol 1,4,5-trisphosphate similar to ryanodine receptor* [R04].

![](figures/discovery_timeline.png)

**{fig:discovery_timeline}.** Four decades of work on the receptor as the
literature records it, with each milestone drawn at the year of its published
source. The figure is generated from a curated table in which no year is
typed: each is read from the reference row it cites, so a milestone cannot be
dated differently from the paper it rests on. The density on the right is the
cryo-electron microscopy era, which is what made the structural half of this
project possible.

## 1.2 The receptor shares its entire domain architecture with a much larger sister family

An IP₃ receptor is a homotetramer. Each subunit is around 2,700 residues, of
which roughly nine tenths are cytosolic. Reading from the N-terminus, a
subunit consists of a β-trefoil suppressor domain, then the IP₃-binding core
(a β-trefoil with an armadillo fold), then a long central solenoid of armadillo
repeats, and then, near the C-terminus, a six-transmembrane pore module of the
voltage-gated-channel superfamily, followed by a tail that folds back into the
cytosolic mass ({fig:channel_structure}). The ligand binds some 100 Å from the
gate, and everything between the two is the machine that couples them.

![](figures/channel_structure.png)

**{fig:channel_structure}.** The channel measured rather than drawn. The
figure shows one C4-symmetric subunit's Cα trace from PDB 6DQN, a human IP₃R3
in the IP₃-bound state at 3.33 Å [R24], together with the four-fold axis, the
axial extent of the membrane, the pore radius profile and the two
constrictions recovered from the coordinates. The narrowest luminal point
lands on the GGGVGD selectivity-filter motif and the cytosolic constriction on
the gate residues, neither of which the geometry was told about. That
agreement is what makes the measurement usable as a coordinate system in
Chapters 11 and 12.

In domain-annotation terms the diagnostic signatures of the family are the MIR
domains (PF02815) in the suppressor region, the IP₃-binding core Ins145_P3_rec
(PF08709), the RyR–IP₃R homology domain RIH (PF01365), the RIH-associated
domain (PF08454), and the generic Ion_trans pore (PF00520).

Every one of the first four is also carried by every human ryanodine receptor
({fig:domain_architecture}). This is not an artefact of database annotation.
The two families are one structural superfamily: their N-terminal regions are
conserved to the point of direct comparison [R58], and the near-atomic
structures published for both in 2015 [R22, R26, R27, R28] show the same
organisation, in which a vast cytosolic solenoid transduces ligand binding to
a C-terminal pore. They differ in scale, the ryanodine receptors being nearly
twice the size at around 4,900 to 5,000 residues.

![](figures/domain_architecture.png)

**{fig:domain_architecture}.** Every signature that defines an IP₃ receptor is
also carried by a ryanodine receptor, in the same copy number, down to the two
RIH domains. What separates the families is what the ryanodine receptors carry
in addition: four further domains and some 2,200 extra residues. A search
built on the shared signatures cannot tell the two apart, which is the
practical form of a statement about descent.

Three consequences of that shared architecture run through every chapter of
this thesis.

**Any similarity search for one family returns the other.** The scale of the
problem is easy to underestimate. A single query for the IP₃-binding-core
signature in zebrafish returns 109 protein records across 10 gene symbols, of
which 53, or 49 %, are ryanodine receptors. Seven of those records come from a
4,900-residue locus that carries no gene name at all. That number is measured
in Chapter 2, in the exact query this project's own planning document had
quoted as a clean result.

**The separation between the two families must therefore be a positive test.**
It cannot be a filter that removes what looks wrong. It has to be evidence
that a record is one thing rather than the other, in the form of best-profile
assignment with a stated margin, or a labelled-bait identity margin, applied
at every stage from the first enumeration to the last structural comparison.
Where the evidence does not reach the required bar, the honest output is a
refusal, and this project produces a great many of them.

**Size is a filter and not evidence.** The two families do separate cleanly by
length in a well-annotated genome. That is a fact about annotation quality in
that genome, and using it as the call would make every result downstream a
restatement of how good the annotation already was, which is one of the things
this thesis set out to measure.

## 1.3 What the channel is for, and why the three paralogues are not interchangeable

The channel's behaviour is not a simple function of ligand concentration. Its
response to cytosolic calcium is biphasic, activating at low concentrations
and inhibiting at high, which under IP₃ produces the bell-shaped
calcium-response curve that has organised thinking about the receptor since it
was measured [R13]. IP₃ binding and calcium binding are cooperative rather
than independent [R14], and the stoichiometry required to initiate release has
been resolved, in that more than one subunit must be occupied [R15]. ATP
modulates the channel at two sites [R16], protein kinase A phosphorylates it
[R17], and it is a hub for protein interactions, including with Bcl-2 [R18]
and IRBIT [R19].

![](figures/gating_logic.png)

**{fig:gating_logic}.** What opens the channel and what closes it, drawn as
the logic rather than as a mechanism. The point for this thesis is the last
row: the same channel is used to produce puffs, waves and oscillations, and
the properties that distinguish those outputs are properties of the individual
paralogue.

Those output differences are the reason the paralogues matter. The three
vertebrate receptors differ in IP₃ affinity, in calcium sensitivity, in
regulation and in tissue distribution [R20, R29, R30], and they are not
redundant. Calcium transfer from the endoplasmic reticulum to the
mitochondrion, and the apoptotic decisions downstream of it, depend on which
paralogue is present [R12]. A cell's calcium signalling repertoire is in part
a statement about which of ITPR1, ITPR2 and ITPR3 it expresses.

## 1.4 Four things about the family that were not established

For a family this well studied, a surprising amount was not established at the
level a genome-scale question requires.

**Its taxonomic range.** The receptor is described as animal machinery with
scattered occurrences elsewhere. Those scattered occurrences had never been
enumerated against a declared search space, and the conspicuous absences,
meaning no IP₃ receptor in *Arabidopsis* and none in budding yeast, had never
been tested at the level of an assembly rather than a gene set. An absence in
a proteome is a statement about what a gene caller found. Whether it is also a
statement about the genome is a different measurement, and nobody had made it.

**Where the three vertebrate paralogues came from.** That vertebrates carry
three is a database fact. That the three arose in the two rounds of
whole-genome duplication at the base of the vertebrates [R179, R180] is a
reasonable inference from their number and their age, and it had not been
tested against the genomic neighbourhoods that a duplication of that kind
leaves behind. Which two of the three are sisters had no published answer at
all. Worse, the proposition that the IP₃ and ryanodine receptor triplications
happened independently, which this project's own background document asserted,
turned out to have no primary source and had to be downgraded to an open
question before any of the work could begin. Chapters 6 and 7 answer it.

**What has happened to the paralogues since.** Gene families of this age
normally lose copies [R193]. Whether this one has, and where, is a question
that cannot be answered without a false-negative rate, because a gene that a
search fails to find looks exactly like a gene that is not there. That single
sentence determined the structure of a third of this thesis.

**How well the family is recorded.** Every one of the questions above is
answered using public databases, and the quality of those databases for this
family had never been audited. The audit turned out to be a result in its own
right.

## 1.5 How every claim in this thesis is scoped

Every result here is scoped to a declared search space. The claim is never
"all IP₃ receptors" but always a specific set: 7,691 UniProt reference
proteomes and 503 NCBI genome assemblies, each set enumerated by a stated rule
and committed as a manifest before any search ran. An absence is an absence
from that space, measured with a positive control that demonstrates the search
could have found the gene had it been there.

That discipline is why the negative results in this document can be read as
results. The family is absent from all 384 land-plant reference proteomes, and
it is also absent from those genomes, in assemblies where a measured control
recovers a comparably long, deeply conserved gene. Those are two different
claims and both are made separately, because the second is the one that means
something biological and the first is what a database can support on its own.

The same discipline runs in the other direction. Chapter 9 reports that not
one of the 927 genome by paralogue cells in the vertebrate scope reaches the
state this project defines as a loss. A zero of that kind is worth nothing
unless the instrument that produced it could have produced something else, so
Chapter 4 measures the search's own false-negative rate on cells whose gene is
independently known to be present, and Chapter 9 reports, cell by cell, the
analytical settings that would manufacture a loss out of this data. The
strongest form in which that result can be stated is the one that names its
own escape routes.

## 1.6 What each chapter of this thesis does

Chapter 2 audits the literature this project started from, claim by claim, and
reports the four statements that did not survive. Chapter 3 builds the
instrument that separates the two receptor families and enumerates the search
space. Chapter 4 takes that instrument to 309 vertebrate genomes and then
measures what the search was worth. Chapter 5 takes it outside the vertebrates
and reports the family's range and its absences.

Chapter 6 builds the alignment and the tree that everything downstream stands
on. Chapter 7 asks where the three paralogues came from, using the genomic
neighbourhood, the species tree and the duplication record. Chapter 8 measures
the gene itself, meaning its exons, its introns, and what a fragmentary
annotation actually is. Chapter 9 counts the losses.

Chapter 10 measures selection across the vertebrate family. Chapter 11 maps
constraint onto the channel and asks whether it explains the clinical
variants. Chapter 12 takes the one part of the receptor the ryanodine
receptors do not share. Chapter 13 audits the archive. Chapter 14 sets out the
methods as the argument behind each instrument rather than as a procedure, and
Chapter 15 discusses what the whole says.
