<!-- Assembled by scripts/s0_review_build.py — do not edit docs/ip3r_review_2026.md directly. -->

# The inositol 1,4,5-trisphosphate receptor family

## Abstract

The inositol 1,4,5-trisphosphate receptor (IP<sub>3</sub>R) is the endoplasmic
reticulum's ligand-gated calcium-release channel and the terminal element of
one of the most widely deployed signalling pathways in eukaryotic biology.
Four decades after IP<sub>3</sub> was shown to release Ca<sup>2+</sup> from an
intracellular store, the receptor is among the best-characterised large
channels in structural biology: a 1.3-MDa homotetramer, resolved by cryo-EM in
apo, ligand-bound and lipid-embedded states, whose ligand-binding sites sit
some 100 Å from the gate they control. It is also, in three of its vertebrate
paralogues, a disease gene — for cerebellar ataxia, for aniridia with ataxia,
for an inability to sweat, and for a demyelinating neuropathy that shades into
multisystem immunodeficiency.

This review consolidates that literature with a specific downstream purpose: it
is the verified baseline for a genome-scale census of the *ITPR* family, and it
is therefore written to distinguish what is established from what is repeated.
Every statement carries a citation to a source retrieved for this purpose;
statements that could not be verified were removed rather than softened, and
the four claims that failed audit are listed explicitly in §12. Two themes run
through the article and matter for anything built on top of it. First, the
IP<sub>3</sub> receptors and the ryanodine receptors are one structural
superfamily, not two families that resemble each other — a fact visible in the
very first sequence and confirmed domain by domain, and one that makes every
sequence-similarity search for one family return the other. Second, the three
vertebrate paralogues are not redundant: they differ in ligand sensitivity,
in regulation, in subcellular placement, in the organelle contacts they
support, and in the diseases they cause.

**Scope.** Molecular architecture and structural biology (§2), ligand binding
and gating (§3), regulation (§4), cellular physiology (§5), paralogue and
splice diversity (§6), evolution (§7), organismal physiology from genetic
models (§8), human disease genetics (§9), pharmacology and experimental tools
(§10), and the questions the field cannot currently answer from its records
(§11). Nomenclature: *ITPR1/2/3* for human genes, *Itpr1/2/3* for other
vertebrates, IP<sub>3</sub>R1/2/3 for proteins; InsP<sub>3</sub>R appears in
quoted titles. RyR denotes the ryanodine receptors throughout.

## 1. From a soluble messenger to a channel

### 1.1 The discovery sequence

The pathway was assembled backwards, from the message to the receiver. In 1983
IP<sub>3</sub> was shown to release Ca<sup>2+</sup> from a non-mitochondrial
intracellular store in permeabilised pancreatic acinar cells [1], and within
a year the same messenger was linked to agonist-evoked Ca<sup>2+</sup>
mobilisation in intact cells [2]. A receptor followed: high-affinity
IP<sub>3</sub> binding sites were characterised in brain membranes and shown to
be regulated by pH and Ca<sup>2+</sup> [3], and mapped across brain
regions, with by far the highest density in cerebellum [4].

The protein those sites belonged to had already been purified, from a different
direction entirely. P400 was purified from mouse cerebellum as a
glycoprotein characteristic of Purkinje cells [5]; it was then shown to be the
IP<sub>3</sub>-binding protein itself, isolated as a large oligomeric complex
[6]. Cloning in 1989 closed the loop from both ends in the same journal
issue: the primary structure of P400 established it as the IP<sub>3</sub>
receptor [7], and an independent cDNA was reported under a title that named
the field's central complication at the outset — a "putative receptor for
inositol 1,4,5-trisphosphate **similar to ryanodine receptor**" [8].

That the purified protein was the channel, and not an accessory to one, was
settled by reconstitution: purified receptor conferred IP<sub>3</sub>-dependent
Ca<sup>2+</sup> flux on reconstituted lipid vesicles [9]; an equivalent
receptor was isolated and characterised from smooth muscle [10].

### 1.2 The first two decades

Three lines of work then developed in parallel, and their conclusions still
frame the subject.

*Gating is not a simple ligand switch.* Single-channel recordings from
cerebellar membranes revealed a bell-shaped dependence on cytosolic
Ca<sup>2+</sup> — activation at low concentrations, inhibition at high [11] —
which made the receptor a candidate oscillator rather than a valve.

*The signal is spatially organised.* Ca<sup>2+</sup> release proved to be
built from discrete elementary events rather than a graded whole-cell response
[12], organised into local and global patterns within single cells [13] and
propagated between cells as intercellular waves [14].

*There is more than one receptor.* Cloning of additional subtypes [15] and the
demonstration that they are expressed tissue- and development-specifically
[16] converted a single channel into a family — one whose members would prove
to differ in ligand sensitivity, regulation and cellular role (§6).

The modern consolidations of this literature are the two *Physiological
Reviews* treatments [17,18], a structural-functional synthesis [19], and
the historical account by one of the field's founders [20].

## 2. Molecular architecture

### 2.1 Domain organisation

Each IP<sub>3</sub>R subunit is ~2,700 residues, of which roughly nine-tenths
are cytosolic. Reading from the N-terminus, the conserved architecture is: a
β-trefoil suppressor domain; the IP<sub>3</sub>-binding core, itself a β-trefoil
plus an armadillo-repeat fold; an extended central region built largely of
armadillo solenoid repeats; and, at the C-terminus, a six-transmembrane pore
module of the voltage-gated-channel superfamily fold, followed by a C-terminal
tail that returns into the cytosolic mass.

In domain-annotation terms the diagnostic signatures are the MIR domains
(PF02815) in the suppressor region, the IP<sub>3</sub>-binding core
(Ins145_P3_rec, PF08709), the RyR–IP<sub>3</sub>R homology domain (RIH,
PF01365), the RIH-associated domain (PF08454) and the Ion_trans pore domain
(PF00520). Three of those five are shared with the ryanodine receptors, and the
fourth — the pore — is shared with most of the cation-channel world. This is
not a technicality of database annotation but a statement about descent
(§7.1), and it is why domain content alone cannot assign a sequence to one
family or the other.

### 2.2 The cryo-EM structures

Low-resolution reconstructions established the overall shape — a large
cytosolic "mushroom cap" over a comparatively small membrane domain — well
before the resolution revolution [21]. The near-atomic era began in 2015 with
the apo structure of tetrameric rat IP<sub>3</sub>R1 at 4.7 Å, which traced
~85% of the backbone and identified the elements involved in gating [22]. It
was contemporaneous with, and interpretable alongside, the three independent
near-atomic RyR1 structures published the same year [23,24,25] — the
comparison that made the superfamily relationship structural rather than
sequence-based.

Subsequent work has filled in the functional states. Ligand-induced allosteric
rearrangements were resolved for IP<sub>3</sub>R1 [26]; human IP<sub>3</sub>R3
was solved in Ca<sup>2+</sup>- and IP<sub>3</sub>-bound states, giving the first
paralogue comparison [27]; the channel was imaged in a lipid bilayer rather
than detergent [28]; and activation and gating were reconstructed across a
ligand series [29], with a companion study resolving the conformational motions
that couple ligand binding to the pore [30]. Structure-guided mutagenesis then
identified which of the candidate Ca<sup>2+</sup> sites are functionally
required for activation [31]. A concise structural overview of the
pre-2018 state of the field is available [32].

### 2.3 The pore and permeation

The pore is formed between TM5 and TM6, with a short selectivity filter and a
gate at the cytosolic end of the TM6 bundle [22,27,29]. The channel is not
a precision Ca<sup>2+</sup> filter in the sense that voltage-gated
Ca<sup>2+</sup> channels are: single-channel recordings from cerebellar
preparations established a large conductance with modest discrimination among
divalent and monovalent cations [33]. Functionally this is appropriate — the
receptor releases Ca<sup>2+</sup> down a steep gradient from a store held at
high concentration, so throughput matters more than selectivity.

### 2.4 Long-range allosteric coupling

The defining structural problem of this receptor is distance. IP<sub>3</sub>
binds at the N-terminal core, roughly 100 Å from the gate. The cryo-EM series
resolves the coupling path: ligand binding closes the clam-shell of the
binding core, that motion is transmitted through the armadillo solenoid, and
the C-terminal tail — which runs from beyond TM6 back up into the cytosolic
domain, forming a left-handed helical bundle at the four-fold axis and
contacting the N-terminal domains of *adjacent* subunits — communicates the
change to the gate [22,26,29,30].

Two consequences follow, and both recur later in this review. Because the
coupling element contacts neighbouring subunits, the channel's allostery is
intrinsically inter-subunit, which is the structural basis for the requirement
that all four ligand sites be occupied (§3.3) and for the dominant-negative
behaviour of heterozygous variants in a tetramer (§9). And because the
C-terminal tail is a load-bearing mechanical element rather than a tail in the
dispensable sense, variants there are not peripheral — they sit on the path
that makes the channel a channel.

## 3. Ligand binding and gating

### 3.1 The IP<sub>3</sub>-binding core and the suppressor domain

The structure of the IP<sub>3</sub>-binding core with ligand bound [34] showed
a clam-shell of two domains closing around the inositol phosphate, with the
4,5-bisphosphate pair making the specificity contacts. Apo and bound structures
of the isolated ligand-binding domain confirmed the conformational cycle
directly [35].

The N-terminal suppressor domain, solved separately [36], reduces
IP<sub>3</sub> affinity — deleting it raises affinity but abolishes gating —
so it is not an inhibitor but part of the transduction apparatus. Its
importance is underlined from the clinical direction: a gain-of-function
variant in the suppressor domain causes cerebellar ataxia [37], and the
Gillespie-syndrome variants cluster in this same N-terminal region and in the
channel domain (§9.1).

### 3.2 Ca<sup>2+</sup> as obligatory co-agonist

The bell-shaped Ca<sup>2+</sup> dependence [11] is the receptor's signature
property, and its interpretation matured over two decades. Marchant and Taylor
showed that the two ligands bind *sequentially* — IP<sub>3</sub> first,
Ca<sup>2+</sup> second — and argued that this ordering is what protects cells
from spontaneous release [38]; spontaneous release does occur when the
safeguard is stressed [39], and even quiescent channels show low-probability
openings [40]. Systematic analysis of gating in permeabilised preparations
established the co-agonist model as an allosteric scheme rather than a
description [41]; the state of the activation problem before the structures
arrived is reviewed in [42].

The structural counterpart arrived with the human IP<sub>3</sub>R3 structures
solved with Ca<sup>2+</sup> and IP<sub>3</sub> [27] and was then tested
functionally: of the candidate Ca<sup>2+</sup>-binding sites visible in the
maps, mutagenesis identified which are actually required for activation [31].
The inhibitory limb of the bell — high Ca<sup>2+</sup> closing the channel —
remains the less well resolved half of the mechanism.

### 3.3 Stoichiometry: all four sites

Concatenated-subunit constructs, in which individual binding sites could be
disabled within a defined tetramer, showed that Ca<sup>2+</sup> release
requires IP<sub>3</sub> occupancy at **all four** subunits [43]. This single
result carries much of the interpretive weight in this review. It explains why
the receptor behaves as a high-order, switch-like detector of IP<sub>3</sub>
rather than a graded one; it explains why a mutant subunit that cannot bind or
transduce ligand poisons the whole tetramer; and it is the mechanistic reason
that heterozygous missense variants in *ITPR1* and *ITPR3* produce dominant
disease (§9), while loss of one gene copy in isolation is comparatively well
tolerated.

The question of whether the four subunits act in strict concert or with some
independence is still argued, including in recent work using pore-directed
peptides to disable subunits selectively.

### 3.4 ATP and the energy state of the cell

ATP modulates the channel at more than one site, shifting open probability
without being a permissive requirement [44], and the effect differs between
paralogues — the type 3 receptor's ATP regulation was characterised separately
and is not identical to the type 1 [45]. Because the concentrations involved
lie within the physiological range, this couples release probability to
cellular energy status, and it is one of several mechanisms by which the same
IP<sub>3</sub> signal produces different outputs in different metabolic
contexts.

## 4. Regulation: a channel as a signalling hub

The IP<sub>3</sub> receptor is regulated less like an ion channel and more like
a scaffold that happens to conduct. Its cytosolic mass presents a large surface
for protein and post-translational input, and the collected partner list runs
to dozens of proteins [46]. What follows is organised by mechanism rather
than by partner, and is deliberately selective.

### 4.1 Phosphorylation

PKA-dependent phosphorylation was the first covalent modification described,
and it *decreased* Ca<sup>2+</sup> release in brain preparations [47];
protein kinase C and Ca<sup>2+</sup>/calmodulin-dependent kinase were shown to
phosphorylate the receptor shortly afterwards [48], and cGMP-dependent kinase
adds a further site [49]. The functional sign of PKA phosphorylation is
paralogue- and cell-type-dependent, which is a recurring pattern: the same
modification on a different family member need not do the same thing.

Phosphorylation also acts upstream of the receptor to shape its input. Control
of Ca<sup>2+</sup> oscillation frequency by phosphorylation of a metabotropic
glutamate receptor [50] is an early and clean demonstration that oscillation
frequency is an encoded variable, not a byproduct.

### 4.2 The Bcl-2 family and the apoptotic set-point

The interaction between anti-apoptotic Bcl-2 proteins and the IP<sub>3</sub>
receptor is the best-developed example of regulation with a clear
pathophysiological reading. Bcl-2 binds the receptor and inhibits
IP<sub>3</sub>-evoked ER Ca<sup>2+</sup> release [51], acting through its BH4
domain on the central modulatory region [52]; Mcl-1 and Bcl-2 modulate the
channel in ways that protect against apoptosis [53]; and the BH4 domains of
Bcl-2 and Bcl-X<sub>L</sub> are not functionally equivalent at this target
[54]. Competition experiments indicate that Bcl-2 and IP<sub>3</sub> can
compete for the ligand-binding domain itself, modulating signalling output
[55], and dual BH3-like domains have been proposed to underlie biphasic
regulation of gating [56]. The area is reviewed comprehensively [57].

The therapeutic corollary is direct: peptides that disrupt the
Bcl-2–IP<sub>3</sub>R interaction reverse Bcl-2's inhibition of apoptotic
Ca<sup>2+</sup> signalling [58], and constitutive IP<sub>3</sub> signalling is
what makes B-cell malignancies sensitive to the BIRD-2 disruptor [59]. The
same axis is a proposed mechanism of Bcl-2's role in cancer more broadly
[60].

### 4.3 Luminal redox and ER environment

ERp44 binds the receptor in a redox- and pH-dependent manner from the ER lumen
and inhibits it, giving the channel a direct read-out of the oxidising state of
the compartment it drains — and it does so subtype-specifically [61]. This is
one of the clearest cases of regulation that is impossible to detect from the
cytosolic side alone.

### 4.4 IRBIT and competitive occupancy

IRBIT binds the IP<sub>3</sub>-binding region and is released when
IP<sub>3</sub> rises [62]; mechanistically it suppresses receptor activity by
competing with IP<sub>3</sub> for the same site [63]. It therefore sets a
threshold: the receptor does not respond to IP<sub>3</sub> until the messenger
has displaced a competitor, which sharpens an already high-order response
(§3.3).

### 4.5 Calcium sensors and calmodulin

A family of neuronal Ca<sup>2+</sup>-binding proteins (CaBPs) was identified as
protein ligands of the receptor, capable of activating it independently of
IP<sub>3</sub> in some conditions [64], and calmodulin modulates the cardiac
type 2 receptor directly [65]. Together with the intrinsic Ca<sup>2+</sup>
sites (§3.2), this means the channel senses Ca<sup>2+</sup> through several
routes at once.

### 4.6 Degradation and quality control

The receptor is not merely modulated but disposed of: sustained stimulation
triggers ubiquitination and down-regulation [66], and an early report
described Ca<sup>2+</sup>-induced proteolysis of the channel. The
ERAD machinery involved is medically relevant in its own right — a point
mutation in the ubiquitin ligase RNF170, which acts on the IP<sub>3</sub>
receptor, causes autosomal dominant sensory ataxia [67]. That a *disposal*
pathway for this receptor causes an ataxia, as loss of the receptor itself does
(§9.1), is a striking convergence.

### 4.7 Disease proteins as regulators

Presenilins regulate IP<sub>3</sub>R gating, and familial Alzheimer-disease
presenilin mutants produce a gain-of-function enhancement of channel activity
[68]. This positions the receptor downstream of a major neurodegenerative
locus, and is part of the wider argument that ER Ca<sup>2+</sup> handling is a
convergence point in neurodegeneration [69].

## 5. Cellular physiology

### 5.1 Elementary events and the hierarchy of signals

Ca<sup>2+</sup> release is quantal. Imaging in *Xenopus* oocytes resolved
"puffs" — discrete, localised release events from small groups of receptors —
as the elementary building blocks of the IP<sub>3</sub>-evoked signal [12].
Below the puff sits the single-channel "blip"; above it, recruitment and
regenerative coupling between release sites generate cell-wide oscillations and
propagating waves [70]. Exocrine cells display the local-to-global transition
particularly clearly, with agonist-evoked release beginning in a restricted
apical pole before spreading [13], and waves propagate between coupled cells
via IP<sub>3</sub> diffusion through gap junctions [14].

### 5.2 Oscillations: frequency as the encoded variable

That repetitive Ca<sup>2+</sup> spikes are not an epiphenomenon was established
early: pulsatile release does not require oscillations in IP<sub>3</sub> itself
[71], and oscillation frequency scales with agonist concentration
[72]. Two classes of model followed — one in which the oscillator is built
from Ca<sup>2+</sup>-induced Ca<sup>2+</sup> release between two pools [73],
and one in which the bell-shaped Ca<sup>2+</sup> dependence of a single pool of
IP<sub>3</sub> receptors suffices [74]. The second, built directly on the
gating properties in §3.2, became the standard framework. Receptor-level
control of frequency was then demonstrated experimentally by manipulating
upstream phosphorylation [50].

The functional significance is that the *same* second messenger can carry
different instructions at different frequencies, and that the receptor's own
gating kinetics are part of the encoder.

### 5.3 Placement: clusters, mobility and licensing

Where the receptors are matters as much as how they gate. Most IP<sub>3</sub>
receptors are mobile and do not participate in the initial response; the
Ca<sup>2+</sup> signal instead initiates at a small population of *immobile*
receptors positioned adjacent to ER–plasma-membrane junctions, effectively a
licensed subset [75]. This resolves a long-standing discrepancy between the
abundance of the protein and the sparseness of release sites, and it means
receptor number is a poor proxy for release capacity — a caveat that applies
directly to any inference from expression data.

### 5.4 Membrane contact sites: mitochondria

The receptor's best-characterised organellar partnership is with mitochondria.
Close ER–mitochondria contacts determine the amplitude of the mitochondrial
Ca<sup>2+</sup> response [76]; the linkage is a defined structural apposition
with a measurable gap and tethering elements [77]; and the resulting transfer
sets the threshold for apoptosis [78]. Tethering machinery has been identified
in yeast [79] and in mammalian neurons [80], and the interface is a signalling
platform in its own right, hosting mTORC2–Akt signalling that feeds back on the
receptor [81].

Critically for a review concerned with paralogues, **the three isoforms are not
equivalent at this interface**: they differ in the ER–mitochondrial contacts
they support and in the efficiency of local Ca<sup>2+</sup> transfer [82].
The physiological consequences are visible at the organ level — chronic
enrichment of hepatic ER–mitochondria contact drives mitochondrial
Ca<sup>2+</sup> overload and metabolic dysfunction in obesity [83] — and the
interface is remodelled in Alzheimer disease [68].

### 5.5 Apoptosis and immune-cell fate

The link between the receptor and cell death was established genetically before
it was structurally: T cells deficient in the IP<sub>3</sub> receptor resist
apoptosis [84], and increased type 3 receptor expression mediates lymphocyte
apoptosis [85]. The Bcl-2 axis (§4.2) supplies the regulatory logic, and the
transfer geometry of §5.4 supplies the mechanism. This is also where paralogue
identity becomes clinically consequential: IP<sub>3</sub>R3 is the isoform most
associated with efficient pro-apoptotic transfer, and it is the isoform whose
overexpression is repeatedly reported in cancers [86,87].

### 5.6 Other acceptors: lysosomes

Mitochondria are not the only acceptor compartment. IP<sub>3</sub> receptors
associate preferentially with ER–lysosome contact sites and deliver
Ca<sup>2+</sup> selectively to lysosomes [88], placing the receptor upstream of
lysosomal Ca<sup>2+</sup> handling and, by extension, of autophagic control.

## 6. Three paralogues, and the diversity within each

### 6.1 The genes

Humans carry three *ITPR* genes. Their protein products are near-identical in
length — IP<sub>3</sub>R1 2,758 aa (UniProt Q14643), IP<sub>3</sub>R2 2,701 aa
(Q14571), IP<sub>3</sub>R3 2,671 aa (Q14573), a spread of 3% — but the genes
that encode them are not built alike. Re-derived from Ensembl release 15.12 for
this review:

| Gene | Locus | Ensembl gene | Canonical exons | Genomic span |
|------|-------|--------------|-----------------|--------------|
| *ITPR1* | 3p26.1 | ENSG00000150995 | 62 | 354,174 bp |
| *ITPR2* | 12p11.23 | ENSG00000123104 | 57 | 497,888 bp |
| *ITPR3* | 6p21.31 | ENSG00000096433 | 58 | 76,245 bp |

The genomic span varies **6.5-fold** across three genes of essentially
identical exon count and protein length. *ITPR3* achieves the same architecture
in roughly one-fifth of the DNA that *ITPR2* uses. Whether this reflects
lineage-specific intron gain, differential loss, or selection related to the
paralogues' distinct expression programmes is not established, and it is a
question the genomic literature has largely not asked.

### 6.2 Functional differences between paralogues

The paralogues are not interchangeable. Single-channel analysis showed
isoform-specific gating [89]; the type 2 receptor was identified and
reconstituted separately and shown to differ functionally [90]; ATP
regulation differs between types [45]; the receptors differ in
IP<sub>3</sub> sensitivity, in Ca<sup>2+</sup> dependence and in their
regulation by the modulators of §4. The cleanest whole-cell demonstration
remains the finding that Ca<sup>2+</sup> signal *encoding* — the shape and
frequency of the response — is set by which subtypes a cell expresses [91],
which recasts subtype composition as a tuning parameter rather than redundancy.
They also differ in the ER–mitochondrial contacts they support [82].

### 6.3 Heterotetramers

Because the paralogues co-assemble, most cells contain a mixed population of
tetramers rather than three separate channel species. This has two
consequences that matter throughout this review. It expands the functional
repertoire combinatorially. And it converts a heterozygous missense variant in
one gene into a defect distributed across a large fraction of a cell's
channels, which — combined with the all-four-sites requirement of §3.3 — is why
dominant-negative disease mechanisms are so prominent in this family (§9).

### 6.4 Alternative splicing

*ITPR1* carries three alternatively spliced segments. SI and SII were defined
in mouse, with tissue- and development-specific usage [16]; the SII region
distinguishes neuronal from non-neuronal forms and the two differ in
phosphorylation [92]; SII splicing was characterised in rat brain versus
peripheral tissues [93]; and a third spliced segment, SIII, was identified on
cloning the human cDNA [94]. Because SII lies in the regulatory region between
the ligand-binding and channel domains, splicing here changes how the receptor
responds rather than whether it conducts.

### 6.5 Tissue distribution

The classical picture — reviewed comprehensively for the three subtypes [95] —
is of broad expression with strong paralogue bias:

- **IP<sub>3</sub>R1** dominates the cerebellum, where it was discovered as the
  Purkinje-cell P400 antigen [5,6], and is the principal neuronal isoform
  [96].
- **IP<sub>3</sub>R2** is broadly expressed, prominent in secretory epithelia
  [97], cardiac myocytes [98] and glia — including oligodendrocytes, where it
  is required for myelination [99].
- **IP<sub>3</sub>R3** is the epithelial and secretory isoform, concentrated
  apically in pancreatic and salivary gland cells [100] and central to
  hepatobiliary function [101,87].

Two cautions. First, this distribution is not fixed: subtype expression is
remodelled in disease, as shown for the two intracellular release-channel types
in end-stage heart failure [102], and quantitative modern atlases have not been
applied systematically across the family. Second, and more importantly, §5.3 shows that
the functionally relevant quantity is not how much receptor a cell contains but
how much of it is licensed and correctly placed — so expression level is a weak
predictor of signalling capacity.

## 7. Evolution

### 7.1 One superfamily, not two families

The IP<sub>3</sub> receptors and the ryanodine receptors are a single
structural superfamily. The relationship was visible in the first
IP<sub>3</sub>R sequence, whose reporting title stated it [8], and it has
been confirmed at every subsequent level of resolution: the receptors share the
MIR, RIH, RIH-associated and Ion_trans domains; their N-terminal regions are
structurally and functionally conserved to the point that domains can be
compared directly [103]; and the 2015 near-atomic structures of both families
[22,23,24,25] show the same overall organisation — a vast cytosolic
solenoid cap transducing ligand binding to a C-terminal pore module — at
different scales (RyR subunits are ~4,900–5,000 aa, nearly twice the size).

Three consequences run through the rest of this article, and through any
computational study of the family.

1. **Any similarity search for one family returns the other.** This is not a
   nuisance to be filtered away; it is a statement about descent, and the
   separation must be made by positive evidence — best-profile assignment or a
   labelled-bait margin — rather than by assumption. The scale of the problem is
   easy to underestimate: a single query for the IP<sub>3</sub>-binding-core
   domain in zebrafish returns 109 protein records across 10 gene symbols, of
   which 53 (49%) are ryanodine receptors, including one unnamed
   4,900-residue locus.
2. **RyR is the natural outgroup** for rooting an IP<sub>3</sub>R phylogeny —
   better conditioned than any invertebrate IP<sub>3</sub>R, because it is a
   genuine sister clade rather than a long branch within the ingroup.
3. **Size is a filter, not evidence.** The two families separate cleanly by
   length in well-annotated genomes, but that is a property of annotation
   quality in those genomes, not a phylogenetic argument.

### 7.2 Origins: the split predates animals

Comparative genomics places the machinery early. Analysis of genomes flanking
the animal–fungal divergence — including the apusozoan *Thecamonas trahens*,
from the putative unicellular sister group to Opisthokonta — finds many
components of animal and fungal Ca<sup>2+</sup> signalling already present in
the common ancestor, together with lineage-specific expansions of
Ca<sup>2+</sup> channels in the unicellular ancestors of animals and in basal
fungi [104]. Homologues of both intracellular release-channel families are
present in parasitic protists [105]. The broader comparative literature
places Ca<sup>2+</sup> signalling as an ancient and elaborate system rather
than a metazoan invention [106,107], with acidic-store channels following
their own evolutionary trajectory [108]; the RyR side is reviewed separately
[109].

### 7.3 The invertebrate single-gene state

Most invertebrates carry a single *itpr*. *Drosophila* is the best-developed
model: disruption of the gene affects larval metamorphosis and
ecdysone release [110], genetic dissection assigns a vital requirement to aminergic neurons
[111], and hypomorphs lose flight and the associated neuronal rhythmicity
[112]. The single *Drosophila* receptor has been characterised biophysically
and behaves recognisably like its vertebrate counterparts [113]. The
*Xenopus* receptor, cloned early, similarly established conservation of
structure and function across vertebrates [114].

A single-gene invertebrate state carrying out the functions distributed across
three vertebrate paralogues is the strongest available argument that
vertebrate paralogue specialisation is subfunctionalisation of an ancestral
repertoire rather than the acquisition of new capabilities.

### 7.4 The vertebrate expansion — and what is not established

Vertebrates carry three paralogues; so, independently or not, do the ryanodine
receptors. Teleosts add a further layer: zebrafish carries four
IP<sub>3</sub> receptor genes — *itpr1a*, *itpr1b*, *itpr2* and *itpr3* — with
the *itpr1a*/*itpr1b* pair bearing the signature of the teleost-specific
whole-genome duplication (longest isoform per gene 2,635–2,819 aa, re-derived
here).

Three questions about this history are **not settled by the existing
literature**, and they are stated as open rather than glossed:

- Whether *ITPR1/2/3* are ohnologues from the two rounds of vertebrate
  whole-genome duplication has not been demonstrated with synteny-backed,
  phylogeny-tested evidence.
- Whether the IP<sub>3</sub>R and RyR triplications were **independent** events
  is frequently asserted and, as far as this review's search could establish,
  nowhere demonstrated. It is an attractive parallel, not a result.
- Which two of the three vertebrate paralogues are sisters — the rooted
  topology of the family — is not fixed by any published, support-annotated
  maximum-likelihood analysis with an RyR outgroup.

### 7.5 The taxonomic range problem

Textbook accounts hold that land plants and fungi lack IP<sub>3</sub>
receptors, and the model organisms support this: *Arabidopsis thaliana* and
*Saccharomyces cerevisiae* carry no protein annotated with the
IP<sub>3</sub>-binding core. Yet the InterPro protein set for that same
signature contains 40 Viridiplantae and 41 Fungi entries (re-derived
2026-08-18, alongside 12,149 Metazoa of 12,338 total). The lineage-specific
expansion of Ca<sup>2+</sup> channels reported in basal fungi [104] means a
real-gene explanation cannot be dismissed a priori — but neither can
mis-annotation or contamination. **Whether the famous absences are facts about
genomes or facts about proteome databases is, at present, unresolved**, and it
cannot be settled from database counts alone.

## 8. Organismal physiology: what the genetic models show

### 8.1 *Itpr1* — the cerebellar phenotype

Mice lacking the type 1 receptor develop ataxia and epileptic seizures and die
young [115], a phenotype that matches both the receptor's Purkinje-cell
enrichment [6,4] and the human ataxias caused by *ITPR1* variants (§9.1).
The concordance between mouse null, spontaneous mouse mutants and human
heterozygous deletion [116] is unusually tight for a large channel gene, and it
is the main reason IP<sub>3</sub>R1's cerebellar role is considered settled.

### 8.2 *Itpr2* and *Itpr3* — secretion, energy balance and glia

Single knockouts of *Itpr2* or *Itpr3* are comparatively mild; the double
knockout is not. *Itpr2*<sup>−/−</sup>*Itpr3*<sup>−/−</sup> mice show severely
impaired exocrine secretion — saliva, pancreatic juice, gastric acid — leading
to difficulties in nutrient assimilation and an energy-metabolism phenotype
[117]. This is the clearest demonstration of paralogue redundancy *within* a
tissue coexisting with paralogue specialisation *between* tissues, and it
matches the apical concentration of these isoforms in secretory epithelia
[100].

Individually, the two genes have distinct assignments. *Itpr2* loss abolishes
endothelin-1-induced arrhythmogenic Ca<sup>2+</sup> signalling in atrial
myocytes [98] and, in oligodendrocytes, disturbs Ca<sup>2+</sup> homeostasis
and impairs myelination with behavioural consequences [99]. In humans, loss
of IP<sub>3</sub>R2 function abolishes sweating [97] — a phenotype reproduced
in *Itpr2*<sup>−/−</sup> mice, which is the strongest single cross-species
validation in the family.

### 8.3 Immune and haematopoietic function

T cells lacking the receptor resist apoptosis [84], and increased
IP<sub>3</sub>R3 expression mediates lymphocyte apoptosis [85]. The human
genetics has since caught up with the mouse work from an unexpected direction:
a recurrent dominant *ITPR3* variant produces severe immunodeficiency with
CD4<sup>+</sup> lymphopenia as part of a multisystem disorder [118] (§9.3).

### 8.4 Invertebrate physiology

The single-gene invertebrate systems (§7.3) make the pleiotropy of this
receptor unusually legible: one gene, and mutants show defects in moulting and
metamorphosis [110], in aminergic neuronal function [111] and in the
rhythmic motor output required for flight [112]. Whatever the three vertebrate
paralogues have specialised into, the ancestral gene was already doing a great
deal.

## 9. Human disease

The *ITPR* family offers an unusually informative clinical panel: three
paralogues, three largely distinct phenotypes, and dominant, dominant-negative
and recessive mechanisms all represented. This makes it a test case for the
proposition that pathogenic variants concentrate in the constrained parts of a
protein — a proposition that can be checked against the structures of §2.

### 9.1 *ITPR1*: cerebellar ataxia and Gillespie syndrome

**SCA15/16.** Heterozygous *deletions* of *ITPR1* cause spinocerebellar ataxia
type 15, identified simultaneously in human families and in the corresponding
ataxic mouse [116]. The finding replicated in Japanese [119] and large European
[120] cohorts, and the separately named SCA16 was shown to be the same entity —
heterozygous *ITPR1* deletion — rather than a distinct disorder [121,122].
The mechanism here is haploinsufficiency; the clinical context is reviewed
among the dominant ataxias [123].

**SCA29.** *Missense* variants cause a congenital, non-progressive, autosomal
dominant ataxia [124], consolidated as SCA29 in a multi-family case series
[125]. At least one such variant is a demonstrated *gain* of function in the
suppressor domain [37] — so the *ITPR1* ataxias are not a single mechanism
with a single direction of effect.

**Gillespie syndrome.** Partial aniridia, cerebellar ataxia and intellectual
disability. The genetics has two arms, and reporting only one of them is a
common error: Gillespie syndrome arises from **both** biallelic recessive
*ITPR1* variants **and** de novo heterozygous variants [126]. Only the latter
act dominant-negatively, through a restricted repertoire of mutations
concentrated in the channel domain [127] — exactly what §3.3 and §6.3 predict,
since one poisoned subunit disables a tetramer that requires all four.

### 9.2 *ITPR2*: autosomal recessive isolated anhidrosis

Loss of IP<sub>3</sub>R2 function causes generalised isolated anhidrosis —
an inability to sweat, with morphologically normal eccrine glands [97]. Three
qualifications belong with the claim. The human evidence is a *single*
consanguineous family with five affected members. The variant is a homozygous
*missense* in the pore-forming region that abolishes Ca<sup>2+</sup> release —
functionally null, not a null allele. And the strongest support is the mouse:
*Itpr2*<sup>−/−</sup> animals show markedly reduced sweat secretion with a
corresponding Ca<sup>2+</sup> defect in the gland. It is a well-supported claim
resting on a narrow base, and it is the only recessive, non-neurological
phenotype in the family.

### 9.3 *ITPR3*: neuropathy, and more than neuropathy

Dominant *ITPR3* variants cause demyelinating Charcot–Marie–Tooth disease
(CMT1J). The gene–disease relationship was established with segregation and
functional evidence in two families, including a variant with a
dominant-negative effect on Ca<sup>2+</sup> transients [128]; further families
followed [129]; and a recurrent p.Thr1424Met was found in 33 affected
individuals from nine unrelated families, with unusually variable age of onset
and severity even within families [130]. A spontaneous canine homozygous
nonsense variant reproduces the neuropathy and adds enamel defects, and reveals
that loss of IP<sub>3</sub>R3 also reduces IP<sub>3</sub>R1 and IP<sub>3</sub>R2
protein levels — evidence of co-regulation across the family [131].

The phenotype is **broader than neuropathy**. The recurrent de novo
p.Arg2524Cys causes a complex multisystem disorder with severe
immunodeficiency, CD4<sup>+</sup> lymphopenia, mitochondrial dysfunction,
ectodermal dysplasia and bone-marrow failure, of which CMT is one component,
acting through a dominant-negative mechanism [118]. Any analysis treating the
*ITPR3* variant set as purely neuropathic will mis-specify it.

### 9.4 Cancer

IP<sub>3</sub>R3 is repeatedly implicated in tumour biology, with
anti-apoptotic and proliferative roles reported in cancer cells [86] and
increased expression enhancing malignant properties in cholangiocarcinoma
[87]. The mechanistic frame is §5.4–5.5: ER-to-mitochondrial Ca<sup>2+</sup>
transfer sets the apoptotic threshold, and Bcl-2-family proteins tune the
receptor to raise it [60]. The therapeutic reading is the BIRD-2 peptide,
which kills B-cell malignancies dependent on constitutive IP<sub>3</sub>
signalling [59].

### 9.5 What the variant set is good for

Because the mechanisms are heterogeneous — haploinsufficiency, gain of
function, dominant-negative poisoning of a tetramer, and recessive loss — the
family supports a question most disease-gene families cannot answer: *do
pathogenic variants sit preferentially in the most constrained regions, and is
the constrained core the same in all three paralogues?* The structural
prerequisites (§2) and the clinical panel (§9.1–9.3) are both in place; what is
missing is a family-wide, alignment-based constraint analysis mapped onto the
cryo-EM coordinates.

## 10. Pharmacology and experimental tools

The IP<sub>3</sub> receptor has no clinically used drug and, for most of its
history, no genuinely selective pharmacological tool. This is worth stating
plainly, because a good deal of the literature rests on reagents whose
limitations are known.

**Heparin** is the classical competitive antagonist at the IP<sub>3</sub>
binding site — potent, reversible and well characterised [132] — and was used
to establish IP<sub>3</sub> dependence of Ca<sup>2+</sup> release in intact
preparations [133]. It is membrane-impermeant and promiscuous (it binds many
proteins), so it is a tool for permeabilised or microinjected preparations
only.

**2-APB** is the most widely used membrane-permeant blocker and the most widely
misinterpreted. It does inhibit IP<sub>3</sub>-induced Ca<sup>2+</sup> release
[134], but it is not selective: it blocks store-operated Ca<sup>2+</sup> entry
at comparable concentrations and affects several TRP channels, and the case
that it should not be used as a diagnostic for IP<sub>3</sub>R involvement was
made explicitly [135]. Results resting on 2-APB alone should be read as
consistent-with rather than demonstrating IP<sub>3</sub>R involvement.

**Xestospongins** were introduced as potent membrane-permeable IP<sub>3</sub>R
blockers [136] and have been valuable, though subsequent work has raised
questions about selectivity and off-target actions on the ER.

**Adenophostin A** and its disaccharide-polyphosphate analogues are the
high-potency agonists — substantially more potent than IP<sub>3</sub> itself —
and remain the main pharmacological probes of the binding site [137].

**Peptide and genetic tools** have proved more selective than small molecules.
BIRD-2 disrupts the Bcl-2–IP<sub>3</sub>R interaction and has a defined
mechanism and a therapeutic rationale [58,59]. Pore-directed peptides are
being used to disable subunits within a tetramer, addressing the
concerted-versus-independent question of §3.3. On the genetic side, the
paralogue knockouts of §8 and concatenated-tetramer constructs [43] have
answered questions no small molecule could.

**Caged IP<sub>3</sub> and imaging** underpin the entire elementary-event
literature: photolysis of caged IP<sub>3</sub> with confocal or TIRF imaging is
what resolved puffs [12] and what localised the licensed receptor population
to ER–plasma-membrane junctions [75].

The gap is conspicuous. For a channel with three paralogues, three diseases and
a documented role in tumour cell survival, there is no isoform-selective
inhibitor. Structures of all three paralogues now exist (§2.2), so
structure-guided design is finally tractable.

## 11. What the field cannot currently answer

This section is deliberately separated from the rest of the review. Everything
above is supported by cited work; what follows is a list of things that are
*asserted more often than they are demonstrated*, or simply unknown. They are
the questions a genome-scale census of the family is designed to answer.

**Q1 — Range.** What is the family's true taxonomic distribution? The
databases say the family is overwhelmingly metazoan (12,149 of 12,338
IP<sub>3</sub>-binding-core proteins), yet hold 40 Viridiplantae and 41 Fungi
records in a family textbooks say plants and fungi lack, while *Arabidopsis*
and *S. cerevisiae* have none. Are those real genes, mis-annotations or
contamination — and are the famous absences facts about genomes or about
proteome databases? Nothing in the current literature settles this, and the
finding of expanded Ca<sup>2+</sup>-channel repertoires in basal fungi [104]
means the answer is not obvious.

**Q2 — Origin.** Are *ITPR1/2/3* ohnologues of the two rounds of vertebrate
whole-genome duplication? Are teleost *itpr1a*/*itpr1b* products of the
teleost-specific third round? Where exactly does the IP<sub>3</sub>R/RyR split
fall relative to the origin of animals? And — the claim most often repeated
without support — were the IP<sub>3</sub>R and RyR expansions to three
vertebrate paralogues **independent** events? A rooted, support-annotated
maximum-likelihood phylogeny with an RyR outgroup, reconciled against a
calibrated species tree, would answer all four.

**Q3 — Fates.** Across a declared vertebrate genome scope, has any paralogue
been lost, how often, and independently or once? Is any surviving as a decaying
pseudogene? Loss claims in this family are currently anecdotal, and an absence
claim is only as good as the assembly and annotation behind it.

**Q4 — Records.** How often is a real *ITPR* locus missing, fragmentary, split
across gene models, unnamed, or filed under the wrong paralogue — or as a
ryanodine receptor? The last risk is not hypothetical: the domain that
*defines* this family returns ryanodine receptors at close to a 1:1 ratio in a
routine query (§7.1), and at least one of them is an unnamed locus.

**Q5 — Constraint and mechanism.** Do the pathogenic variants of §9 sit in the
most constrained parts of the channel, and is the constrained core the same in
all three paralogues? Is the IP<sub>3</sub>-binding core — the one module whose
function RyR does not share — under different constraint from the pore, which
RyR does? Does it change in lineages that have lost the upstream
PLC/IP<sub>3</sub> pathway?

**Q6 — Anything unnamed.** Is there a fourth vertebrate IP<sub>3</sub>
receptor, or an IP<sub>3</sub> receptor in a lineage reported to lack one?

**Q7 — Gene architecture.** Why does *ITPR3* achieve the same 58-exon
architecture in 76 kb that *ITPR2* spreads over 498 kb (§6.1)? A 6.5-fold span
difference among paralogues of identical protein length is a substantial,
unexplained observation.

**Q8 — Selective pharmacology.** No isoform-selective inhibitor exists (§10),
despite three paralogues with distinct disease associations and now-available
structures for all three.

## 12. How this review was built, and what failed verification

### 12.1 Method

Candidate literature was harvested from Europe PMC across ~60 topic queries
spanning the sections above, each run twice — ranked by citation count and by
recency — to avoid a purely canonical or purely recent bias. Candidates were
curated by hand; bibliographic metadata for every retained reference was then
fetched programmatically from Europe PMC by PMID, so no author list, title,
journal or year in the bibliography was transcribed manually. The final set is
**137 references — 117 primary research articles and 20 reviews** — spanning
1983 to 2025, with publication type taken from the Europe PMC record rather
than assigned by hand. The reference
table is committed as `results/s0_baseline/references.tsv`, and this document
is assembled from its section files by `scripts/s0_review_build.py`, which
renumbers the stable citation keys into order of first appearance and renders
the bibliography from that table. A cited key with no reference row is a build
error, so the text and the bibliography cannot drift apart.

Database figures quoted in §6.1, §7.1, §7.5 and §11 were re-derived for this
review rather than repeated: InterPro protein counts and taxonomic
distributions from the InterPro API, protein lengths and domain architectures
from UniProt, and gene structure from Ensembl release 15.12. The scripts are
`scripts/s0_db_snapshot.py` and `scripts/s0_gene_structure.py`; the outputs are
the committed tables in `results/s0_baseline/`.

**Limitations.** This is a narrative review, not a systematic one: there is no
pre-registered protocol, no PRISMA flow, and no formal inclusion criteria
beyond topical relevance and the preference for primary sources over reviews
for specific claims. Search was restricted to Europe PMC and to
English-language records. Coverage is deliberately deeper on structure,
evolution and genetics than on cell-type-specific physiology, reflecting the
downstream use.

### 12.2 The audit: four claims that did not survive

This review replaced a planning document whose statements were tagged by how
far they could be trusted. Nineteen atomic claims were checked against the
literature. Twelve were verified as written and three were kept with the
wording qualified; the remaining four are recorded here because a review that
silently drops what it could not confirm is less useful than one that says so.

| Claim as originally written | Verdict | Resolution |
|---|---|---|
| Loci *ITPR1* 3p26.1, *ITPR2* 12p11.23, *ITPR3* 6p21.31 | Confirmed by database | Re-derived from Ensembl; all three exact. Promoted from literature claim to database fact (§6.1). |
| "Each is a ~58–60 exon gene spanning hundreds of kb" | **False as written** | Measured: 57–62 canonical exons, and *ITPR3* spans 76 kb, not "hundreds of kb". Replaced by the measured table in §6.1; the 6.5-fold span asymmetry became open question Q7. |
| "Both families independently expanded to three vertebrate paralogues" | **Not established** | The three-paralogue state of each family is a database fact; the *independence* of the two triplications has no primary source. Moved to open question Q2 (§7.4, §11). |
| Gillespie syndrome "behaves as dominant-negative in the tetramer" | **Incomplete** | Both biallelic recessive and de novo dominant-negative mechanisms exist [126,127]. Corrected in §9.1. |

Three further claims are retained but qualified: "most agonist-evoked
Ca<sup>2+</sup> signals in non-muscle cells" is a review-level generalisation
with no measured fraction behind it (§1.2, §5.1); the paralogue tissue-emphasis
summary is qualitative and not a quantitative ranking (§6.5); and the *ITPR2*
anhidrosis claim rests on a single family plus a mouse model (§9.2). One
numerical correction was also applied: the tetramer is **1.3 MDa** as reported
in the primary structural work [22] — four copies of the 2,758-residue human
subunit give ≈1.25 MDa — and should not be quoted as "~1.2 MDa" unattributed.

### 12.3 Reading this review

Cite the sources, not this document. Where a statement here is a database
figure rather than a literature claim it says so and names the release; those
figures should be re-derived before being quoted, because they will move.

## References

1. Streb H, Irvine RF, Berridge MJ, Schulz I. Release of Ca2+ from a nonmitochondrial intracellular store in pancreatic acinar cells by inositol-1,4,5-trisphosphate. *Nature* **1983**. PMID [6605482](https://pubmed.ncbi.nlm.nih.gov/6605482/). doi:[10.1038/306067a0](https://doi.org/10.1038/306067a0).
2. Berridge MJ, Heslop JP, Irvine RF, Brown KD. Inositol trisphosphate formation and calcium mobilization in Swiss 3T3 cells in response to platelet-derived growth factor. *Biochem J* **1984**. PMID [6089758](https://pubmed.ncbi.nlm.nih.gov/6089758/). doi:[10.1042/bj2220195](https://doi.org/10.1042/bj2220195).
3. Worley PF, Baraban JM, Supattapone S, Wilson VS, Snyder SH. Characterization of inositol trisphosphate receptor binding in brain. Regulation by pH and calcium. *J Biol Chem* **1987**. PMID [3040730](https://pubmed.ncbi.nlm.nih.gov/3040730/). doi:[10.1016/s0021-9258(18)45326-4](https://doi.org/10.1016/s0021-9258(18)45326-4).
4. Worley PF, Baraban JM, Colvin JS, Snyder SH. Inositol trisphosphate receptor localization in brain: variable stoichiometry with protein kinase C. *Nature* **1987**. PMID [3027583](https://pubmed.ncbi.nlm.nih.gov/3027583/). doi:[10.1038/325159a0](https://doi.org/10.1038/325159a0).
5. Maeda N, Niinobe M, Nakahira K, Mikoshiba K. Purification and characterization of P400 protein, a glycoprotein characteristic of Purkinje cell, from mouse cerebellum. *J Neurochem* **1988**. PMID [3141586](https://pubmed.ncbi.nlm.nih.gov/3141586/). doi:[10.1111/j.1471-4159.1988.tb01151.x](https://doi.org/10.1111/j.1471-4159.1988.tb01151.x).
6. Maeda N, Niinobe M, Mikoshiba K. A cerebellar Purkinje cell marker P400 protein is an inositol 1,4,5-trisphosphate (InsP3) receptor protein. Purification and characterization of InsP3 receptor complex. *EMBO J* **1990**. PMID [2153079](https://pubmed.ncbi.nlm.nih.gov/2153079/). doi:[10.1002/j.1460-2075.1990.tb07386.x](https://doi.org/10.1002/j.1460-2075.1990.tb07386.x).
7. Furuichi T, Yoshikawa S, Miyawaki A, Wada K, Maeda N, Mikoshiba K. Primary structure and functional expression of the inositol 1,4,5-trisphosphate-binding protein P400. *Nature* **1989**. PMID [2554142](https://pubmed.ncbi.nlm.nih.gov/2554142/). doi:[10.1038/342032a0](https://doi.org/10.1038/342032a0).
8. Mignery GA, Sudhof TC, Takei K, De Camilli P. Putative receptor for inositol 1,4,5-trisphosphate similar to ryanodine receptor. *Nature* **1989**. PMID [2554146](https://pubmed.ncbi.nlm.nih.gov/2554146/). doi:[10.1038/342192a0](https://doi.org/10.1038/342192a0).
9. Ferris CD, Huganir RL, Supattapone S, Snyder SH. Purified inositol 1,4,5-trisphosphate receptor mediates calcium flux in reconstituted lipid vesicles. *Nature* **1989**. PMID [2554143](https://pubmed.ncbi.nlm.nih.gov/2554143/). doi:[10.1038/342087a0](https://doi.org/10.1038/342087a0).
10. Chadwick CC, Saito A, Fleischer S. Isolation and characterization of the inositol trisphosphate receptor from smooth muscle. *Proc Natl Acad Sci USA* **1990**. PMID [2156261](https://pubmed.ncbi.nlm.nih.gov/2156261/). doi:[10.1073/pnas.87.6.2132](https://doi.org/10.1073/pnas.87.6.2132).
11. Bezprozvanny I, Watras J, Ehrlich BE. Bell-shaped calcium-response curves of Ins(1,4,5)P3- and calcium-gated channels from endoplasmic reticulum of cerebellum. *Nature* **1991**. PMID [1648178](https://pubmed.ncbi.nlm.nih.gov/1648178/). doi:[10.1038/351751a0](https://doi.org/10.1038/351751a0).
12. Yao Y, Choi J, Parker I. Quantal puffs of intracellular Ca2+ evoked by inositol trisphosphate in Xenopus oocytes. *J Physiol* **1995**. PMID [7738847](https://pubmed.ncbi.nlm.nih.gov/7738847/). doi:[10.1113/jphysiol.1995.sp020703](https://doi.org/10.1113/jphysiol.1995.sp020703).
13. Thorn P, Lawrie AM, Smith PM, Gallacher DV, Petersen OH. Local and global cytosolic Ca2+ oscillations in exocrine cells evoked by agonists and inositol trisphosphate. *Cell* **1993**. PMID [8395347](https://pubmed.ncbi.nlm.nih.gov/8395347/). doi:[10.1016/0092-8674(93)90513-p](https://doi.org/10.1016/0092-8674(93)90513-p).
14. Boitano S, Dirksen ER, Sanderson MJ. Intercellular propagation of calcium waves mediated by inositol trisphosphate. *Science* **1992**. PMID [1411526](https://pubmed.ncbi.nlm.nih.gov/1411526/). doi:[10.1126/science.1411526](https://doi.org/10.1126/science.1411526).
15. Ross CA, Danoff SK, Schell MJ, Snyder SH, Ullrich A. Three additional inositol 1,4,5-trisphosphate receptors: molecular cloning and differential localization in brain and peripheral tissues. *Proc Natl Acad Sci USA* **1992**. PMID [1374893](https://pubmed.ncbi.nlm.nih.gov/1374893/). doi:[10.1073/pnas.89.10.4265](https://doi.org/10.1073/pnas.89.10.4265).
16. Nakagawa T, Okano H, Furuichi T, Aruga J, Mikoshiba K. The subtypes of the mouse inositol 1,4,5-trisphosphate receptor are expressed in a tissue-specific and developmentally specific manner. *Proc Natl Acad Sci USA* **1991**. PMID [1648733](https://pubmed.ncbi.nlm.nih.gov/1648733/). doi:[10.1073/pnas.88.14.6244](https://doi.org/10.1073/pnas.88.14.6244).
17. Foskett JK, White C, Cheung KH, Mak DO. Inositol trisphosphate receptor Ca2+ release channels. *Physiol Rev* **2007**. PMID [17429043](https://pubmed.ncbi.nlm.nih.gov/17429043/). doi:[10.1152/physrev.00035.2006](https://doi.org/10.1152/physrev.00035.2006).
18. Berridge MJ. The Inositol Trisphosphate/Calcium Signaling Pathway in Health and Disease. *Physiol Rev* **2016**. PMID [27512009](https://pubmed.ncbi.nlm.nih.gov/27512009/). doi:[10.1152/physrev.00006.2016](https://doi.org/10.1152/physrev.00006.2016).
19. Prole DL, Taylor CW. Structure and Function of IP3 Receptors. *Cold Spring Harb Perspect Biol* **2019**. PMID [30745293](https://pubmed.ncbi.nlm.nih.gov/30745293/). doi:[10.1101/cshperspect.a035063](https://doi.org/10.1101/cshperspect.a035063).
20. Mikoshiba K. IP3 receptor/Ca2+ channel: from discovery to new signaling concepts. *J Neurochem* **2007**. PMID [17697045](https://pubmed.ncbi.nlm.nih.gov/17697045/). doi:[10.1111/j.1471-4159.2007.04825.x](https://doi.org/10.1111/j.1471-4159.2007.04825.x).
21. Serysheva II, Bare DJ, Ludtke SJ, Kettlun CS, Chiu W, Mignery GA. Structure of the type 1 inositol 1,4,5-trisphosphate receptor revealed by electron cryomicroscopy. *J Biol Chem* **2003**. PMID [12714606](https://pubmed.ncbi.nlm.nih.gov/12714606/). doi:[10.1074/jbc.c300148200](https://doi.org/10.1074/jbc.c300148200).
22. Fan G, Baker ML, Wang Z, Baker MR, Sinyagovskiy PA, Chiu W *et al.* Gating machinery of InsP3R channels revealed by electron cryomicroscopy. *Nature* **2015**. PMID [26458101](https://pubmed.ncbi.nlm.nih.gov/26458101/). doi:[10.1038/nature15249](https://doi.org/10.1038/nature15249).
23. Zalk R, Clarke OB, des Georges A, Grassucci RA, Reiken S, Mancia F *et al.* Structure of a mammalian ryanodine receptor. *Nature* **2015**. PMID [25470061](https://pubmed.ncbi.nlm.nih.gov/25470061/). doi:[10.1038/nature13950](https://doi.org/10.1038/nature13950).
24. Yan Z, Bai X, Yan C, Wu J, Li Z, Xie T *et al.* Structure of the rabbit ryanodine receptor RyR1 at near-atomic resolution. *Nature* **2015**. PMID [25517095](https://pubmed.ncbi.nlm.nih.gov/25517095/). doi:[10.1038/nature14063](https://doi.org/10.1038/nature14063).
25. Efremov RG, Leitner A, Aebersold R, Raunser S. Architecture and conformational switch mechanism of the ryanodine receptor. *Nature* **2015**. PMID [25470059](https://pubmed.ncbi.nlm.nih.gov/25470059/). doi:[10.1038/nature13916](https://doi.org/10.1038/nature13916).
26. Fan G, Baker MR, Wang Z, Seryshev AB, Ludtke SJ, Baker ML *et al.* Cryo-EM reveals ligand induced allostery underlying InsP3R channel gating. *Cell Res* **2018**. PMID [30470765](https://pubmed.ncbi.nlm.nih.gov/30470765/). doi:[10.1038/s41422-018-0108-5](https://doi.org/10.1038/s41422-018-0108-5).
27. Paknejad N, Hite RK. Structural basis for the regulation of inositol trisphosphate receptors by Ca2+ and IP3. *Nat Struct Mol Biol* **2018**. PMID [30013099](https://pubmed.ncbi.nlm.nih.gov/30013099/). doi:[10.1038/s41594-018-0089-6](https://doi.org/10.1038/s41594-018-0089-6).
28. Baker MR, Fan G, Seryshev AB, Agosto MA, Baker ML, Serysheva II. Cryo-EM structure of type 1 IP3R channel in a lipid bilayer. *Commun Biol* **2021**. PMID [34035440](https://pubmed.ncbi.nlm.nih.gov/34035440/). doi:[10.1038/s42003-021-02156-4](https://doi.org/10.1038/s42003-021-02156-4).
29. Schmitz EA, Takahashi H, Karakas E. Structural basis for activation and gating of IP3 receptors. *Nat Commun* **2022**. PMID [35301323](https://pubmed.ncbi.nlm.nih.gov/35301323/). doi:[10.1038/s41467-022-29073-2](https://doi.org/10.1038/s41467-022-29073-2).
30. Fan G, Baker MR, Terry LE, Arige V, Chen M, Seryshev AB *et al.* Conformational motions and ligand-binding underlying gating and regulation in IP3R channel. *Nat Commun* **2022**. PMID [36376291](https://pubmed.ncbi.nlm.nih.gov/36376291/). doi:[10.1038/s41467-022-34574-1](https://doi.org/10.1038/s41467-022-34574-1).
31. Arige V, Terry LE, Wagner LE, Malik S, Baker MR, Fan G *et al.* Functional determination of calcium-binding sites required for the activation of inositol 1,4,5-trisphosphate receptors. *Proc Natl Acad Sci USA* **2022**. PMID [36122240](https://pubmed.ncbi.nlm.nih.gov/36122240/). doi:[10.1073/pnas.2209267119](https://doi.org/10.1073/pnas.2209267119).
32. Baker MR, Fan G, Serysheva II. Structure of IP3R channel: high-resolution insights from cryo-EM. *Curr Opin Struct Biol* **2017**. PMID [28618351](https://pubmed.ncbi.nlm.nih.gov/28618351/). doi:[10.1016/j.sbi.2017.05.014](https://doi.org/10.1016/j.sbi.2017.05.014).
33. Bezprozvanny I, Ehrlich BE. Inositol (1,4,5)-trisphosphate (InsP3)-gated Ca channels from cerebellum: conduction properties for divalent cations and regulation by intraluminal calcium. *J Gen Physiol* **1994**. PMID [7876825](https://pubmed.ncbi.nlm.nih.gov/7876825/). doi:[10.1085/jgp.104.5.821](https://doi.org/10.1085/jgp.104.5.821).
34. Bosanac I, Alattia JR, Mal TK, Chan J, Talarico S, Tong FK *et al.* Structure of the inositol 1,4,5-trisphosphate receptor binding core in complex with its ligand. *Nature* **2002**. PMID [12442173](https://pubmed.ncbi.nlm.nih.gov/12442173/). doi:[10.1038/nature01268](https://doi.org/10.1038/nature01268).
35. Lin CC, Baek K, Lu Z. Apo and InsP₃-bound crystal structures of the ligand-binding domain of an InsP₃ receptor. *Nat Struct Mol Biol* **2011**. PMID [21892169](https://pubmed.ncbi.nlm.nih.gov/21892169/). doi:[10.1038/nsmb.2112](https://doi.org/10.1038/nsmb.2112).
36. Bosanac I, Yamazaki H, Matsu-Ura T, Michikawa T, Mikoshiba K, Ikura M. Crystal structure of the ligand binding suppressor domain of type 1 inositol 1,4,5-trisphosphate receptor. *Mol Cell* **2005**. PMID [15664189](https://pubmed.ncbi.nlm.nih.gov/15664189/). doi:[10.1016/j.molcel.2004.11.046](https://doi.org/10.1016/j.molcel.2004.11.046).
37. Casey JP, Hirouchi T, Hisatsune C, Lynch B, Murphy R, Dunne AM *et al.* A novel gain-of-function mutation in the ITPR1 suppressor domain causes spinocerebellar ataxia with altered Ca2+ signal patterns. *J Neurol* **2017**. PMID [28620721](https://pubmed.ncbi.nlm.nih.gov/28620721/). doi:[10.1007/s00415-017-8545-5](https://doi.org/10.1007/s00415-017-8545-5).
38. Marchant JS, Taylor CW. Cooperative activation of IP3 receptors by sequential binding of IP3 and Ca2+ safeguards against spontaneous activity. *Curr Biol* **1997**. PMID [9210378](https://pubmed.ncbi.nlm.nih.gov/9210378/). doi:[10.1016/s0960-9822(06)00230-5](https://doi.org/10.1016/s0960-9822(06)00230-5).
39. Missiaen L, Taylor CW, Berridge MJ. Spontaneous calcium release from inositol trisphosphate-sensitive calcium stores. *Nature* **1991**. PMID [1857419](https://pubmed.ncbi.nlm.nih.gov/1857419/). doi:[10.1038/352241a0](https://doi.org/10.1038/352241a0).
40. Mak DO, McBride SM, Foskett JK. Spontaneous channel activity of the inositol 1,4,5-trisphosphate (InsP3) receptor (InsP3R). Application of allosteric modeling to calcium and InsP3 regulation of InsP3R single-channel gating. *J Gen Physiol* **2003**. PMID [14581584](https://pubmed.ncbi.nlm.nih.gov/14581584/). doi:[10.1085/jgp.200308809](https://doi.org/10.1085/jgp.200308809).
41. Hirose K, Kadowaki S, Iino M. Allosteric regulation by cytoplasmic Ca2+ and IP3 of the gating of IP3 receptors in permeabilized guinea-pig vascular smooth muscle cells. *J Physiol* **1998**. PMID [9490868](https://pubmed.ncbi.nlm.nih.gov/9490868/). doi:[10.1111/j.1469-7793.1998.407bw.x](https://doi.org/10.1111/j.1469-7793.1998.407bw.x).
42. Taylor CW, Tovey SC. IP(3) receptors: toward understanding their activation. *Cold Spring Harb Perspect Biol* **2010**. PMID [20980441](https://pubmed.ncbi.nlm.nih.gov/20980441/). doi:[10.1101/cshperspect.a004010](https://doi.org/10.1101/cshperspect.a004010).
43. Alzayady KJ, Wang L, Chandrasekhar R, Wagner LE, Van Petegem F, Yule DI. Defining the stoichiometry of inositol 1,4,5-trisphosphate binding required to initiate Ca2+ release. *Sci Signal* **2016**. PMID [27048566](https://pubmed.ncbi.nlm.nih.gov/27048566/). doi:[10.1126/scisignal.aad6281](https://doi.org/10.1126/scisignal.aad6281).
44. Bezprozvanny I, Ehrlich BE. ATP modulates the function of inositol 1,4,5-trisphosphate-gated channels at two sites. *Neuron* **1993**. PMID [7686381](https://pubmed.ncbi.nlm.nih.gov/7686381/). doi:[10.1016/0896-6273(93)90319-m](https://doi.org/10.1016/0896-6273(93)90319-m).
45. Mak DO, McBride S, Foskett JK. ATP regulation of recombinant type 3 inositol 1,4,5-trisphosphate receptor gating. *J Gen Physiol* **2001**. PMID [11331355](https://pubmed.ncbi.nlm.nih.gov/11331355/). doi:[10.1085/jgp.117.5.447](https://doi.org/10.1085/jgp.117.5.447).
46. Prole DL, Taylor CW. Inositol 1,4,5-trisphosphate receptors and their protein partners as signalling hubs. *J Physiol* **2016**. PMID [26830355](https://pubmed.ncbi.nlm.nih.gov/26830355/). doi:[10.1113/jp271139](https://doi.org/10.1113/jp271139).
47. Supattapone S, Danoff SK, Theibert A, Joseph SK, Steiner J, Snyder SH. Cyclic AMP-dependent phosphorylation of a brain inositol trisphosphate receptor decreases its release of calcium. *Proc Natl Acad Sci USA* **1988**. PMID [2847175](https://pubmed.ncbi.nlm.nih.gov/2847175/). doi:[10.1073/pnas.85.22.8747](https://doi.org/10.1073/pnas.85.22.8747).
48. Ferris CD, Huganir RL, Bredt DS, Cameron AM, Snyder SH. Inositol trisphosphate receptor: phosphorylation by protein kinase C and calcium calmodulin-dependent protein kinases in reconstituted lipid vesicles. *Proc Natl Acad Sci USA* **1991**. PMID [1848697](https://pubmed.ncbi.nlm.nih.gov/1848697/). doi:[10.1073/pnas.88.6.2232](https://doi.org/10.1073/pnas.88.6.2232).
49. Komalavilas P, Lincoln TM. Phosphorylation of the inositol 1,4,5-trisphosphate receptor by cyclic GMP-dependent protein kinase. *J Biol Chem* **1994**. PMID [8132598](https://pubmed.ncbi.nlm.nih.gov/8132598/). doi:[10.1016/s0021-9258(17)37024-2](https://doi.org/10.1016/s0021-9258(17)37024-2).
50. Kawabata S, Tsutsumi R, Kohara A, Yamaguchi T, Nakanishi S, Okada M. Control of calcium oscillations by phosphorylation of metabotropic glutamate receptors. *Nature* **1996**. PMID [8779726](https://pubmed.ncbi.nlm.nih.gov/8779726/). doi:[10.1038/383089a0](https://doi.org/10.1038/383089a0).
51. Chen R, Valencia I, Zhong F, McColl KS, Roderick HL, Bootman MD *et al.* Bcl-2 functionally interacts with inositol 1,4,5-trisphosphate receptors to regulate calcium release from the ER in response to inositol 1,4,5-trisphosphate. *J Cell Biol* **2004**. PMID [15263017](https://pubmed.ncbi.nlm.nih.gov/15263017/). doi:[10.1083/jcb.200402193](https://doi.org/10.1083/jcb.200402193).
52. Rong YP, Bultynck G, Aromolaran AS, Zhong F, Parys JB, De Smedt H *et al.* The BH4 domain of Bcl-2 inhibits ER calcium release and apoptosis by binding the regulatory and coupling domain of the IP3 receptor. *Proc Natl Acad Sci USA* **2009**. PMID [19706527](https://pubmed.ncbi.nlm.nih.gov/19706527/). doi:[10.1073/pnas.0907555106](https://doi.org/10.1073/pnas.0907555106).
53. Eckenrode EF, Yang J, Velmurugan GV, Foskett JK, White C. Apoptosis protection by Mcl-1 and Bcl-2 modulation of inositol 1,4,5-trisphosphate receptor-dependent Ca2+ signaling. *J Biol Chem* **2010**. PMID [20189983](https://pubmed.ncbi.nlm.nih.gov/20189983/). doi:[10.1074/jbc.m109.096040](https://doi.org/10.1074/jbc.m109.096040).
54. Monaco G, Decrock E, Akl H, Ponsaerts R, Vervliet T, Luyten T *et al.* Selective regulation of IP3-receptor-mediated Ca2+ signaling and apoptosis by the BH4 domain of Bcl-2 versus Bcl-Xl. *Cell Death Differ* **2012**. PMID [21818117](https://pubmed.ncbi.nlm.nih.gov/21818117/). doi:[10.1038/cdd.2011.97](https://doi.org/10.1038/cdd.2011.97).
55. Ivanova H, Wagner LE, Tanimura A, Vandermarliere E, Luyten T, Welkenhuyzen K *et al.* Bcl-2 and IP3 compete for the ligand-binding domain of IP3Rs modulating Ca2+ signaling output. *Cell Mol Life Sci* **2019**. PMID [30989245](https://pubmed.ncbi.nlm.nih.gov/30989245/). doi:[10.1007/s00018-019-03091-8](https://doi.org/10.1007/s00018-019-03091-8).
56. Yang J, Vais H, Gu W, Foskett JK. Biphasic regulation of InsP3 receptor gating by dual Ca2+ release channel BH3-like domains mediates Bcl-xL control of cell viability. *Proc Natl Acad Sci USA* **2016**. PMID [26976600](https://pubmed.ncbi.nlm.nih.gov/26976600/). doi:[10.1073/pnas.1517935113](https://doi.org/10.1073/pnas.1517935113).
57. Ivanova H, Vervliet T, Monaco G, Terry LE, Rosa N, Baker MR *et al.* Bcl-2-Protein Family as Modulators of IP3 Receptors and Other Organellar Ca2+ Channels. *Cold Spring Harb Perspect Biol* **2020**. PMID [31501195](https://pubmed.ncbi.nlm.nih.gov/31501195/). doi:[10.1101/cshperspect.a035089](https://doi.org/10.1101/cshperspect.a035089).
58. Rong YP, Aromolaran AS, Bultynck G, Zhong F, Li X, McColl K *et al.* Targeting Bcl-2-IP3 receptor interaction to reverse Bcl-2's inhibition of apoptotic calcium signals. *Mol Cell* **2008**. PMID [18657507](https://pubmed.ncbi.nlm.nih.gov/18657507/). doi:[10.1016/j.molcel.2008.06.014](https://doi.org/10.1016/j.molcel.2008.06.014).
59. Bittremieux M, La Rovere RM, Akl H, Martines C, Welkenhuyzen K, Dubron K *et al.* Constitutive IP3 signaling underlies the sensitivity of B-cell cancers to the Bcl-2/IP3 receptor disruptor BIRD-2. *Cell Death Differ* **2019**. PMID [29899382](https://pubmed.ncbi.nlm.nih.gov/29899382/). doi:[10.1038/s41418-018-0142-3](https://doi.org/10.1038/s41418-018-0142-3).
60. Akl H, Vervloessem T, Kiviluoto S, Bittremieux M, Parys JB, De Smedt H *et al.* A dual role for the anti-apoptotic Bcl-2 protein in cancer: mitochondria versus endoplasmic reticulum. *Biochim Biophys Acta* **2014**. PMID [24768714](https://pubmed.ncbi.nlm.nih.gov/24768714/). doi:[10.1016/j.bbamcr.2014.04.017](https://doi.org/10.1016/j.bbamcr.2014.04.017).
61. Higo T, Hattori M, Nakamura T, Natsume T, Michikawa T, Mikoshiba K. Subtype-specific and ER lumenal environment-dependent regulation of inositol 1,4,5-trisphosphate receptor type 1 by ERp44. *Cell* **2005**. PMID [15652484](https://pubmed.ncbi.nlm.nih.gov/15652484/). doi:[10.1016/j.cell.2004.11.048](https://doi.org/10.1016/j.cell.2004.11.048).
62. Ando H, Mizutani A, Matsu-ura T, Mikoshiba K. IRBIT, a novel inositol 1,4,5-trisphosphate (IP3) receptor-binding protein, is released from the IP3 receptor upon IP3 binding to the receptor. *J Biol Chem* **2003**. PMID [12525476](https://pubmed.ncbi.nlm.nih.gov/12525476/). doi:[10.1074/jbc.M210119200](https://doi.org/10.1074/jbc.M210119200).
63. Ando H, Mizutani A, Kiefer H, Tsuzurugi D, Michikawa T, Mikoshiba K. IRBIT suppresses IP3 receptor activity by competing with IP3 for the common binding site on the IP3 receptor. *Mol Cell* **2006**. PMID [16793548](https://pubmed.ncbi.nlm.nih.gov/16793548/). doi:[10.1016/j.molcel.2006.05.017](https://doi.org/10.1016/j.molcel.2006.05.017).
64. Yang J, McBride S, Mak DO, Vardi N, Palczewski K, Haeseleer F *et al.* Identification of a family of calcium sensors as protein ligands of inositol trisphosphate receptor Ca(2+) release channels. *Proc Natl Acad Sci USA* **2002**. PMID [12032348](https://pubmed.ncbi.nlm.nih.gov/12032348/). doi:[10.1073/pnas.102006299](https://doi.org/10.1073/pnas.102006299).
65. Bare DJ, Kettlun CS, Liang M, Bers DM, Mignery GA. Cardiac type 2 inositol 1,4,5-trisphosphate receptor: interaction and modulation by calcium/calmodulin-dependent protein kinase II. *J Biol Chem* **2005**. PMID [15710625](https://pubmed.ncbi.nlm.nih.gov/15710625/). doi:[10.1074/jbc.m414212200](https://doi.org/10.1074/jbc.m414212200).
66. Wojcikiewicz RJ, Ernst SA, Yule DI. Secretagogues cause ubiquitination and down-regulation of inositol 1, 4,5-trisphosphate receptors in rat pancreatic acinar cells. *Gastroenterology* **1999**. PMID [10220512](https://pubmed.ncbi.nlm.nih.gov/10220512/). doi:[10.1016/s0016-5085(99)70023-5](https://doi.org/10.1016/s0016-5085(99)70023-5).
67. Wright FA, Lu JP, Sliter DA, Dupré N, Rouleau GA, Wojcikiewicz RJ. A Point Mutation in the Ubiquitin Ligase RNF170 That Causes Autosomal Dominant Sensory Ataxia Destabilizes the Protein and Impairs Inositol 1,4,5-Trisphosphate Receptor-mediated Ca2+ Signaling. *J Biol Chem* **2015**. PMID [25882839](https://pubmed.ncbi.nlm.nih.gov/25882839/). doi:[10.1074/jbc.m115.655043](https://doi.org/10.1074/jbc.m115.655043).
68. Cheung KH, Shineman D, Müller M, Cárdenas C, Mei L, Yang J *et al.* Mechanism of Ca2+ disruption in Alzheimer's disease by presenilin regulation of InsP3 receptor channel gating. *Neuron* **2008**. PMID [18579078](https://pubmed.ncbi.nlm.nih.gov/18579078/). doi:[10.1016/j.neuron.2008.04.015](https://doi.org/10.1016/j.neuron.2008.04.015).
69. Foskett JK. Inositol trisphosphate receptor Ca2+ release channels in neurological diseases. *Pflugers Arch* **2010**. PMID [20383523](https://pubmed.ncbi.nlm.nih.gov/20383523/). doi:[10.1007/s00424-010-0826-0](https://doi.org/10.1007/s00424-010-0826-0).
70. Berridge MJ. Inositol trisphosphate and calcium signalling. *Nature* **1993**. PMID [8381210](https://pubmed.ncbi.nlm.nih.gov/8381210/). doi:[10.1038/361315a0](https://doi.org/10.1038/361315a0).
71. Wakui M, Potter BV, Petersen OH. Pulsatile intracellular calcium release does not depend on fluctuations in inositol trisphosphate concentration. *Nature* **1989**. PMID [2498663](https://pubmed.ncbi.nlm.nih.gov/2498663/). doi:[10.1038/339317a0](https://doi.org/10.1038/339317a0).
72. Harootunian AT, Kao JP, Paranjape S, Tsien RY. Generation of calcium oscillations in fibroblasts by positive feedback between calcium and IP3. *Science* **1991**. PMID [1986413](https://pubmed.ncbi.nlm.nih.gov/1986413/). doi:[10.1126/science.1986413](https://doi.org/10.1126/science.1986413).
73. Goldbeter A, Dupont G, Berridge MJ. Minimal model for signal-induced Ca2+ oscillations and for their frequency encoding through protein phosphorylation. *Proc Natl Acad Sci USA* **1990**. PMID [2304911](https://pubmed.ncbi.nlm.nih.gov/2304911/). doi:[10.1073/pnas.87.4.1461](https://doi.org/10.1073/pnas.87.4.1461).
74. De Young GW, Keizer J. A single-pool inositol 1,4,5-trisphosphate-receptor-based model for agonist-stimulated oscillations in Ca2+ concentration. *Proc Natl Acad Sci USA* **1992**. PMID [1329108](https://pubmed.ncbi.nlm.nih.gov/1329108/). doi:[10.1073/pnas.89.20.9895](https://doi.org/10.1073/pnas.89.20.9895).
75. Thillaiappan NB, Chavda AP, Tovey SC, Prole DL, Taylor CW. Ca2+ signals initiate at immobile IP3 receptors adjacent to ER-plasma membrane junctions. *Nat Commun* **2017**. PMID [29138405](https://pubmed.ncbi.nlm.nih.gov/29138405/). doi:[10.1038/s41467-017-01644-8](https://doi.org/10.1038/s41467-017-01644-8).
76. Rizzuto R, Pinton P, Carrington W, Fay FS, Fogarty KE, Lifshitz LM *et al.* Close contacts with the endoplasmic reticulum as determinants of mitochondrial Ca2+ responses. *Science* **1998**. PMID [9624056](https://pubmed.ncbi.nlm.nih.gov/9624056/). doi:[10.1126/science.280.5370.1763](https://doi.org/10.1126/science.280.5370.1763).
77. Csordas G, Renken C, Varnai P, Walter L, Weaver D, Buttle KF *et al.* Structural and functional features and significance of the physical linkage between ER and mitochondria. *J Cell Biol* **2006**. PMID [16982799](https://pubmed.ncbi.nlm.nih.gov/16982799/). doi:[10.1083/jcb.200604016](https://doi.org/10.1083/jcb.200604016).
78. Pinton P, Giorgi C, Siviero R, Zecchini E, Rizzuto R. Calcium and apoptosis: ER-mitochondria Ca2+ transfer in the control of apoptosis. *Oncogene* **2008**. PMID [18955969](https://pubmed.ncbi.nlm.nih.gov/18955969/). doi:[10.1038/onc.2008.308](https://doi.org/10.1038/onc.2008.308).
79. Kornmann B, Currie E, Collins SR, Schuldiner M, Nunnari J, Weissman JS *et al.* An ER-mitochondria tethering complex revealed by a synthetic biology screen. *Science* **2009**. PMID [19556461](https://pubmed.ncbi.nlm.nih.gov/19556461/). doi:[10.1126/science.1175088](https://doi.org/10.1126/science.1175088).
80. Hirabayashi Y, Kwon SK, Paek H, Pernice WM, Paul MA, Lee J *et al.* ER-mitochondria tethering by PDZD8 regulates Ca2+ dynamics in mammalian neurons. *Science* **2017**. PMID [29097544](https://pubmed.ncbi.nlm.nih.gov/29097544/). doi:[10.1126/science.aan6009](https://doi.org/10.1126/science.aan6009).
81. Betz C, Stracka D, Prescianotto-Baschong C, Frieden M, Demaurex N, Hall MN. Feature Article: mTOR complex 2-Akt signaling at mitochondria-associated endoplasmic reticulum membranes (MAM) regulates mitochondrial physiology. *Proc Natl Acad Sci USA* **2013**. PMID [23852728](https://pubmed.ncbi.nlm.nih.gov/23852728/). doi:[10.1073/pnas.1302455110](https://doi.org/10.1073/pnas.1302455110).
82. Bartok A, Weaver D, Golenar T, Nichtova Z, Katona M, Bansaghi S *et al.* IP3 receptor isoforms differently regulate ER-mitochondrial contacts and local calcium transfer. *Nat Commun* **2019**. PMID [31427578](https://pubmed.ncbi.nlm.nih.gov/31427578/). doi:[10.1038/s41467-019-11646-3](https://doi.org/10.1038/s41467-019-11646-3).
83. Arruda AP, Pers BM, Parlakgül G, Güney E, Inouye K, Hotamisligil GS. Chronic enrichment of hepatic endoplasmic reticulum-mitochondria contact leads to mitochondrial dysfunction in obesity. *Nat Med* **2014**. PMID [25419710](https://pubmed.ncbi.nlm.nih.gov/25419710/). doi:[10.1038/nm.3735](https://doi.org/10.1038/nm.3735).
84. Jayaraman T, Marks AR. T cells deficient in inositol 1,4,5-trisphosphate receptor are resistant to apoptosis. *Mol Cell Biol* **1997**. PMID [9154798](https://pubmed.ncbi.nlm.nih.gov/9154798/). doi:[10.1128/mcb.17.6.3005](https://doi.org/10.1128/mcb.17.6.3005).
85. Khan AA, Soloski MJ, Sharp AH, Schilling G, Sabatini DM, Li SH *et al.* Lymphocyte apoptosis: mediation by increased type 3 inositol 1,4,5-trisphosphate receptor. *Science* **1996**. PMID [8662540](https://pubmed.ncbi.nlm.nih.gov/8662540/). doi:[10.1126/science.273.5274.503](https://doi.org/10.1126/science.273.5274.503).
86. Rezuchova I, Hudecova S, Soltysova A, Matuskova M, Durinikova E, Chovancova B *et al.* Type 3 inositol 1,4,5-trisphosphate receptor has antiapoptotic and proliferative role in cancer cells. *Cell Death Dis* **2019**. PMID [30796197](https://pubmed.ncbi.nlm.nih.gov/30796197/). doi:[10.1038/s41419-019-1433-4](https://doi.org/10.1038/s41419-019-1433-4).
87. Ueasilamongkol P, Khamphaya T, Guerra MT, Rodrigues MA, Gomes DA, Kong Y *et al.* Type 3 Inositol 1,4,5-Trisphosphate Receptor Is Increased and Enhances Malignant Properties in Cholangiocarcinoma. *Hepatology* **2020**. PMID [31251815](https://pubmed.ncbi.nlm.nih.gov/31251815/). doi:[10.1002/hep.30839](https://doi.org/10.1002/hep.30839).
88. Atakpa P, Thillaiappan NB, Mataragka S, Prole DL, Taylor CW. IP3 Receptors Preferentially Associate with ER-Lysosome Contact Sites and Selectively Deliver Ca2+ to Lysosomes. *Cell Rep* **2018**. PMID [30540949](https://pubmed.ncbi.nlm.nih.gov/30540949/). doi:[10.1016/j.celrep.2018.11.064](https://doi.org/10.1016/j.celrep.2018.11.064).
89. Ramos-Franco J, Fill M, Mignery GA. Isoform-specific function of single inositol 1,4,5-trisphosphate receptor channels. *Biophys J* **1998**. PMID [9675184](https://pubmed.ncbi.nlm.nih.gov/9675184/). doi:[10.1016/s0006-3495(98)77572-1](https://doi.org/10.1016/s0006-3495(98)77572-1).
90. Perez PJ, Ramos-Franco J, Fill M, Mignery GA. Identification and functional reconstitution of the type 2 inositol 1,4,5-trisphosphate receptor from ventricular cardiac myocytes. *J Biol Chem* **1997**. PMID [9295347](https://pubmed.ncbi.nlm.nih.gov/9295347/). doi:[10.1074/jbc.272.38.23961](https://doi.org/10.1074/jbc.272.38.23961).
91. Miyakawa T, Maeda A, Yamazawa T, Hirose K, Kurosaki T, Iino M. Encoding of Ca2+ signals by differential expression of IP3 receptor subtypes. *EMBO J* **1999**. PMID [10064596](https://pubmed.ncbi.nlm.nih.gov/10064596/). doi:[10.1093/emboj/18.5.1303](https://doi.org/10.1093/emboj/18.5.1303).
92. Danoff SK, Ferris CD, Donath C, Fischer GA, Munemitsu S, Ullrich A *et al.* Inositol 1,4,5-trisphosphate receptors: distinct neuronal and nonneuronal forms derived by alternative splicing differ in phosphorylation. *Proc Natl Acad Sci USA* **1991**. PMID [1849282](https://pubmed.ncbi.nlm.nih.gov/1849282/). doi:[10.1073/pnas.88.7.2951](https://doi.org/10.1073/pnas.88.7.2951).
93. Schell MJ, Danoff SK, Ross CA. Inositol (1,4,5)-trisphosphate receptor: characterization of neuron-specific alternative splicing in rat brain and peripheral tissues. *Brain Res Mol Brain Res* **1993**. PMID [8389956](https://pubmed.ncbi.nlm.nih.gov/8389956/). doi:[10.1016/0169-328x(93)90004-9](https://doi.org/10.1016/0169-328x(93)90004-9).
94. Nucifora FC, Li SH, Danoff S, Ullrich A, Ross CA. Molecular cloning of a cDNA for the human inositol 1,4,5-trisphosphate receptor type 1, and the identification of a third alternatively spliced variant. *Brain Res Mol Brain Res* **1995**. PMID [7500840](https://pubmed.ncbi.nlm.nih.gov/7500840/). doi:[10.1016/0169-328x(95)00194-w](https://doi.org/10.1016/0169-328x(95)00194-w).
95. Taylor CW, Genazzani AA, Morris SA. Expression of inositol trisphosphate receptors. *Cell Calcium* **1999**. PMID [10668562](https://pubmed.ncbi.nlm.nih.gov/10668562/). doi:[10.1054/ceca.1999.0034](https://doi.org/10.1054/ceca.1999.0034).
96. Sharp AH, Nucifora FC, Blondel O, Sheppard CA, Zhang C, Snyder SH *et al.* Differential cellular expression of isoforms of inositol 1,4,5-triphosphate receptors in neurons and glia in brain. *J Comp Neurol* **1999**. PMID [10096607](https://pubmed.ncbi.nlm.nih.gov/10096607/). doi:[10.1002/(sici)1096-9861(19990405)406:23.0.co;2-7](https://doi.org/10.1002/(sici)1096-9861(19990405)406:23.0.co;2-7).
97. Klar J, Hisatsune C, Baig SM, Tariq M, Johansson AC, Rasool M *et al.* Abolished InsP3R2 function inhibits sweat secretion in both humans and mice. *J Clin Invest* **2014**. PMID [25329695](https://pubmed.ncbi.nlm.nih.gov/25329695/). doi:[10.1172/JCI78173](https://doi.org/10.1172/JCI78173).
98. Li X, Zima AV, Sheikh F, Blatter LA, Chen J. Endothelin-1-induced arrhythmogenic Ca2+ signaling is abolished in atrial myocytes of inositol-1,4,5-trisphosphate(IP3)-receptor type 2-deficient mice. *Circ Res* **2005**. PMID [15933266](https://pubmed.ncbi.nlm.nih.gov/15933266/). doi:[10.1161/01.res.0000172556.05576.4c](https://doi.org/10.1161/01.res.0000172556.05576.4c).
99. Zhang M, Zhi N, Feng J, Liu Y, Zhang M, Liu D *et al.* ITPR2 Mediated Calcium Homeostasis in Oligodendrocytes is Essential for Myelination and Involved in Depressive-Like Behavior in Adolescent Mice. *Adv Sci* **2024**. PMID [38476116](https://pubmed.ncbi.nlm.nih.gov/38476116/). doi:[10.1002/advs.202306498](https://doi.org/10.1002/advs.202306498).
100. Lee MG, Xu X, Zeng W, Diaz J, Wojcikiewicz RJ, Kuo TH *et al.* Polarized expression of Ca2+ channels in pancreatic and salivary gland cells. Correlation with initiation and propagation of [Ca2+]i waves. *J Biol Chem* **1997**. PMID [9188472](https://pubmed.ncbi.nlm.nih.gov/9188472/). doi:[10.1074/jbc.272.25.15765](https://doi.org/10.1074/jbc.272.25.15765).
101. Mangla A, Guerra MT, Nathanson MH. Type 3 inositol 1,4,5-trisphosphate receptor: A calcium channel for all seasons. *Cell Calcium* **2020**. PMID [31790953](https://pubmed.ncbi.nlm.nih.gov/31790953/). doi:[10.1016/j.ceca.2019.102132](https://doi.org/10.1016/j.ceca.2019.102132).
102. Go LO, Moschella MC, Watras J, Handa KK, Fyfe BS, Marks AR. Differential regulation of two types of intracellular calcium release channels during end-stage heart failure. *J Clin Invest* **1995**. PMID [7860772](https://pubmed.ncbi.nlm.nih.gov/7860772/). doi:[10.1172/jci117739](https://doi.org/10.1172/jci117739).
103. Seo MD, Velamakanni S, Ishiyama N, Stathopulos PB, Rossi AM, Khan SA *et al.* Structural and functional conservation of key domains in InsP3 and ryanodine receptors. *Nature* **2012**. PMID [22286060](https://pubmed.ncbi.nlm.nih.gov/22286060/). doi:[10.1038/nature10751](https://doi.org/10.1038/nature10751).
104. Cai X, Clapham DE. Ancestral Ca2+ signaling machinery in early animal and fungal evolution. *Mol Biol Evol* **2012**. PMID [21680871](https://pubmed.ncbi.nlm.nih.gov/21680871/). doi:[10.1093/molbev/msr149](https://doi.org/10.1093/molbev/msr149).
105. Prole DL, Taylor CW. Identification of intracellular and plasma membrane calcium channel homologues in pathogenic parasites. *PLoS One* **2011**. PMID [22022573](https://pubmed.ncbi.nlm.nih.gov/22022573/). doi:[10.1371/journal.pone.0026218](https://doi.org/10.1371/journal.pone.0026218).
106. Plattner H, Verkhratsky A. Ca2+ signalling early in evolution--all but primitive. *J Cell Sci* **2013**. PMID [23729741](https://pubmed.ncbi.nlm.nih.gov/23729741/). doi:[10.1242/jcs.127449](https://doi.org/10.1242/jcs.127449).
107. Plattner H. Molecular aspects of calcium signalling at the crossroads of unikont and bikont eukaryote evolution--the ciliated protozoan Paramecium in focus. *Cell Calcium* **2015**. PMID [25601027](https://pubmed.ncbi.nlm.nih.gov/25601027/). doi:[10.1016/j.ceca.2014.12.002](https://doi.org/10.1016/j.ceca.2014.12.002).
108. Patel S, Cai X. Evolution of acidic Ca²⁺ stores and their resident Ca²⁺-permeable channels. *Cell Calcium* **2015**. PMID [25591931](https://pubmed.ncbi.nlm.nih.gov/25591931/). doi:[10.1016/j.ceca.2014.12.005](https://doi.org/10.1016/j.ceca.2014.12.005).
109. Mackrill JJ. Ryanodine receptor calcium release channels: an evolutionary perspective. *Adv Exp Med Biol* **2012**. PMID [22453942](https://pubmed.ncbi.nlm.nih.gov/22453942/). doi:[10.1007/978-94-007-2888-2_7](https://doi.org/10.1007/978-94-007-2888-2_7).
110. Venkatesh K, Hasan G. Disruption of the IP3 receptor gene of Drosophila affects larval metamorphosis and ecdysone release. *Curr Biol* **1997**. PMID [9273145](https://pubmed.ncbi.nlm.nih.gov/9273145/). doi:[10.1016/s0960-9822(06)00221-1](https://doi.org/10.1016/s0960-9822(06)00221-1).
111. Joshi R, Venkatesh K, Srinivas R, Nair S, Hasan G. Genetic dissection of itpr gene function reveals a vital requirement in aminergic cells of Drosophila larvae. *Genetics* **2004**. PMID [15020420](https://pubmed.ncbi.nlm.nih.gov/15020420/). doi:[10.1534/genetics.166.1.225](https://doi.org/10.1534/genetics.166.1.225).
112. Banerjee S, Lee J, Venkatesh K, Wu CF, Hasan G. Loss of flight and associated neuronal rhythmicity in inositol 1,4,5-trisphosphate receptor mutants of Drosophila. *J Neurosci* **2004**. PMID [15356199](https://pubmed.ncbi.nlm.nih.gov/15356199/). doi:[10.1523/jneurosci.0656-04.2004](https://doi.org/10.1523/jneurosci.0656-04.2004).
113. Srikanth S, Wang Z, Tu H, Nair S, Mathew MK, Hasan G *et al.* Functional properties of the Drosophila melanogaster inositol 1,4,5-trisphosphate receptor mutants. *Biophys J* **2004**. PMID [15189860](https://pubmed.ncbi.nlm.nih.gov/15189860/). doi:[10.1529/biophysj.104.040121](https://doi.org/10.1529/biophysj.104.040121).
114. Kume S, Muto A, Aruga J, Nakagawa T, Michikawa T, Michikawa T *et al.* The Xenopus IP3 receptor: structure, function, and localization in oocytes and eggs. *Cell* **1993**. PMID [8387895](https://pubmed.ncbi.nlm.nih.gov/8387895/). doi:[10.1016/0092-8674(93)90142-d](https://doi.org/10.1016/0092-8674(93)90142-d).
115. Matsumoto M, Nakagawa T, Inoue T, Nagata E, Tanaka K, Takano H *et al.* Ataxia and epileptic seizures in mice lacking type 1 inositol 1,4,5-trisphosphate receptor. *Nature* **1996**. PMID [8538767](https://pubmed.ncbi.nlm.nih.gov/8538767/). doi:[10.1038/379168a0](https://doi.org/10.1038/379168a0).
116. van de Leemput J, Chandran J, Knight MA, Holtzclaw LA, Scholz S, Cookson MR *et al.* Deletion at ITPR1 underlies ataxia in mice and spinocerebellar ataxia 15 in humans. *PLoS Genet* **2007**. PMID [17590087](https://pubmed.ncbi.nlm.nih.gov/17590087/). doi:[10.1371/journal.pgen.0030108](https://doi.org/10.1371/journal.pgen.0030108).
117. Futatsugi A, Nakamura T, Yamada MK, Ebisui E, Nakamura K, Uchida K *et al.* IP3 receptor types 2 and 3 mediate exocrine secretion underlying energy metabolism. *Science* **2005**. PMID [16195467](https://pubmed.ncbi.nlm.nih.gov/16195467/). doi:[10.1126/science.1114110](https://doi.org/10.1126/science.1114110).
118. Molitor A, Lederle A, Radosavljevic M, Sapuru V, Zavorka Thomas ME, Yang J *et al.* A pleiotropic recurrent dominant ITPR3 variant causes a complex multisystemic disease. *J Exp Med* **2024**. PMID [39270020](https://pubmed.ncbi.nlm.nih.gov/39270020/). doi:[10.1084/jem.20232178](https://doi.org/10.1084/jem.20232178).
119. Hara K, Shiga A, Nozaki H, Mitsui J, Takahashi Y, Ishiguro H *et al.* Total deletion and a missense mutation of ITPR1 in Japanese SCA15 families. *Neurology* **2008**. PMID [18579805](https://pubmed.ncbi.nlm.nih.gov/18579805/). doi:[10.1212/01.wnl.0000311912.05593.1e](https://doi.org/10.1212/01.wnl.0000311912.05593.1e).
120. Marelli C, van de Leemput J, Johnson JO, Tison F, Thauvin-Robinet C, Picard F *et al.* SCA15 due to large ITPR1 deletions in a cohort of 333 white families with dominant ataxia. *Arch Neurol* **2011**. PMID [21555639](https://pubmed.ncbi.nlm.nih.gov/21555639/). doi:[10.1001/archneurol.2011.81](https://doi.org/10.1001/archneurol.2011.81).
121. Iwaki A, Kawano Y, Miura S, Shibata H, Matsuse D, Li W *et al.* Heterozygous deletion of ITPR1, but not SUMF1, in spinocerebellar ataxia type 16. *J Med Genet* **2008**. PMID [17932120](https://pubmed.ncbi.nlm.nih.gov/17932120/). doi:[10.1136/jmg.2007.053942](https://doi.org/10.1136/jmg.2007.053942).
122. Novak MJ, Sweeney MG, Li A, Treacy C, Chandrashekar HS, Giunti P *et al.* An ITPR1 gene deletion causes spinocerebellar ataxia 15/16: a genetic, clinical and radiological description. *Mov Disord* **2010**. PMID [20669319](https://pubmed.ncbi.nlm.nih.gov/20669319/). doi:[10.1002/mds.23223](https://doi.org/10.1002/mds.23223).
123. Durr A. Autosomal dominant cerebellar ataxias: polyglutamine expansions and beyond. *Lancet Neurol* **2010**. PMID [20723845](https://pubmed.ncbi.nlm.nih.gov/20723845/). doi:[10.1016/s1474-4422(10)70183-6](https://doi.org/10.1016/s1474-4422(10)70183-6).
124. Huang L, Chardon JW, Carter MT, Friend KL, Dudding TE, Schwartzentruber J *et al.* Missense mutations in ITPR1 cause autosomal dominant congenital nonprogressive spinocerebellar ataxia. *Orphanet J Rare Dis* **2012**. PMID [22986007](https://pubmed.ncbi.nlm.nih.gov/22986007/). doi:[10.1186/1750-1172-7-67](https://doi.org/10.1186/1750-1172-7-67).
125. Zambonin JL, Bellomo A, Ben-Pazi H, Everman DB, Frazer LM, Geraghty MT *et al.* Spinocerebellar ataxia type 29 due to mutations in ITPR1: a case series and review of this emerging congenital ataxia. *Orphanet J Rare Dis* **2017**. PMID [28659154](https://pubmed.ncbi.nlm.nih.gov/28659154/). doi:[10.1186/s13023-017-0672-7](https://doi.org/10.1186/s13023-017-0672-7).
126. Gerber S, Alzayady KJ, Burglen L, Bremond-Gignac D, Marchesin V, Roche O *et al.* Recessive and Dominant De Novo ITPR1 Mutations Cause Gillespie Syndrome. *Am J Hum Genet* **2016**. PMID [27108797](https://pubmed.ncbi.nlm.nih.gov/27108797/). doi:[10.1016/j.ajhg.2016.03.004](https://doi.org/10.1016/j.ajhg.2016.03.004).
127. McEntagart M, Williamson KA, Rainger JK, Wheeler A, Seawright A, De Baere E *et al.* A Restricted Repertoire of De Novo Mutations in ITPR1 Cause Gillespie Syndrome with Evidence for Dominant-Negative Effect. *Am J Hum Genet* **2016**. PMID [27108798](https://pubmed.ncbi.nlm.nih.gov/27108798/). doi:[10.1016/j.ajhg.2016.03.018](https://doi.org/10.1016/j.ajhg.2016.03.018).
128. Ronkko J, Molchanova S, Revah-Politi A, Pereira EM, Auranen M, Toppila J *et al.* Dominant mutations in ITPR3 cause Charcot-Marie-Tooth disease. *Ann Clin Transl Neurol* **2020**. PMID [32949214](https://pubmed.ncbi.nlm.nih.gov/32949214/). doi:[10.1002/acn3.51151](https://doi.org/10.1002/acn3.51151).
129. Cabello-Murgui J, Jimenez-Jimenez J, Vilchez JJ, Azorin I, Marti-Martinez P, Millet E *et al.* ITPR3-associated neuropathy: Report of a further family with adult onset intermediate Charcot-Marie-Tooth disease. *Eur J Neurol* **2024**. PMID [39287469](https://pubmed.ncbi.nlm.nih.gov/39287469/). doi:[10.1111/ene.16466](https://doi.org/10.1111/ene.16466).
130. Beijer D, Dohrn MF, Rebelo A, Danzi MC, Grosz BR, Ellis M *et al.* A recurrent missense variant in ITPR3 causes demyelinating Charcot-Marie-Tooth with variable severity. *Brain* **2025**. PMID [38938188](https://pubmed.ncbi.nlm.nih.gov/38938188/). doi:[10.1093/brain/awae178](https://doi.org/10.1093/brain/awae178).
131. Hytonen MK, Ronkko J, Hundi S, Jokinen TS, Suonto E, Teravainen E *et al.* IP3 receptor depletion in a spontaneous canine model of Charcot-Marie-Tooth disease 1J with amelogenesis imperfecta. *Dis Model Mech* **2025**. PMID [39804930](https://pubmed.ncbi.nlm.nih.gov/39804930/). doi:[10.1242/dmm.052078](https://doi.org/10.1242/dmm.052078).
132. Ghosh TK, Eis PS, Mullaney JM, Ebert CL, Gill DL. Competitive, reversible, and potent antagonism of inositol 1,4,5-trisphosphate-activated calcium release by heparin. *J Biol Chem* **1988**. PMID [3136153](https://pubmed.ncbi.nlm.nih.gov/3136153/). doi:[10.1016/s0021-9258(18)37923-7](https://doi.org/10.1016/s0021-9258(18)37923-7).
133. Kobayashi S, Kitazawa T, Somlyo AV, Somlyo AP. Cytosolic heparin inhibits muscarinic and alpha-adrenergic Ca2+ release in smooth muscle. Physiological role of inositol 1,4,5-trisphosphate in pharmacomechanical coupling. *J Biol Chem* **1989**. PMID [2509451](https://pubmed.ncbi.nlm.nih.gov/2509451/). doi:[10.1016/s0021-9258(19)84670-7](https://doi.org/10.1016/s0021-9258(19)84670-7).
134. Peppiatt CM, Collins TJ, Mackenzie L, Conway SJ, Holmes AB, Bootman MD *et al.* 2-Aminoethoxydiphenyl borate (2-APB) antagonises inositol 1,4,5-trisphosphate-induced calcium release, inhibits calcium pumps and has a use-dependent and slowly reversible action on store-operated calcium entry channels. *Cell Calcium* **2003**. PMID [12767897](https://pubmed.ncbi.nlm.nih.gov/12767897/). doi:[10.1016/s0143-4160(03)00026-5](https://doi.org/10.1016/s0143-4160(03)00026-5).
135. Bootman MD, Collins TJ, Mackenzie L, Roderick HL, Berridge MJ, Peppiatt CM. 2-aminoethoxydiphenyl borate (2-APB) is a reliable blocker of store-operated Ca2+ entry but an inconsistent inhibitor of InsP3-induced Ca2+ release. *FASEB J* **2002**. PMID [12153982](https://pubmed.ncbi.nlm.nih.gov/12153982/). doi:[10.1096/fj.02-0037rev](https://doi.org/10.1096/fj.02-0037rev).
136. Gafni J, Munsch JA, Lam TH, Catlin MC, Costa LG, Molinski TF *et al.* Xestospongins: potent membrane permeable blockers of the inositol 1,4,5-trisphosphate receptor. *Neuron* **1997**. PMID [9331361](https://pubmed.ncbi.nlm.nih.gov/9331361/). doi:[10.1016/s0896-6273(00)80384-0](https://doi.org/10.1016/s0896-6273(00)80384-0).
137. Marchant JS, Beecroft MD, Riley AM, Jenkins DJ, Marwood RD, Taylor CW *et al.* Disaccharide polyphosphates based upon adenophostin A activate hepatic D-myo-inositol 1,4,5-trisphosphate receptors. *Biochemistry* **1997**. PMID [9335535](https://pubmed.ncbi.nlm.nih.gov/9335535/). doi:[10.1021/bi971397v](https://doi.org/10.1021/bi971397v).
