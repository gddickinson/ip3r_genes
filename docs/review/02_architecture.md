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
family or the other ({fig:domain_architecture}).

![](figures/domain_architecture.png)

**{fig:domain_architecture}.** Every signature that defines an IP<sub>3</sub> receptor is also carried by a ryanodine receptor. (**a**) The six human proteins drawn to one scale from InterPro coordinates. (**b**) Copies per subunit; a blank cell means the signature is absent. The four shared signatures and the generic pore are identical in content between the families — including copy number, down to the two RIH domains. What separates them is what RyR carries *in addition*: four further domains and some 2,200 extra residues. A search built on the domains in the blue and amber columns cannot distinguish the two families, which is the practical form of §7.1's argument about descent.


### 2.2 The cryo-EM structures

Low-resolution reconstructions established the overall shape — a large
cytosolic "mushroom cap" over a comparatively small membrane domain — well
before the resolution revolution [R62]. The near-atomic era began in 2015 with
the apo structure of tetrameric rat IP<sub>3</sub>R1 at 4.7 Å, which traced
~85% of the backbone and identified the elements involved in gating [R22]. It
was contemporaneous with, and interpretable alongside, the three independent
near-atomic RyR1 structures published the same year [R26, R27, R28] — the
comparison that made the superfamily relationship structural rather than
sequence-based.

Subsequent work has filled in the functional states. Ligand-induced allosteric
rearrangements were resolved for IP<sub>3</sub>R1 [R23]; human IP<sub>3</sub>R3
was solved in Ca<sup>2+</sup>- and IP<sub>3</sub>-bound states, giving the first
paralogue comparison [R24]; the channel was imaged in a lipid bilayer rather
than detergent [R25]; and activation and gating were reconstructed across a
ligand series [R59], with a companion study resolving the conformational motions
that couple ligand binding to the pore [R60]. Structure-guided mutagenesis then
identified which of the candidate Ca<sup>2+</sup> sites are functionally
required for activation [R61]. A concise structural overview of the
pre-2018 state of the field is available [R63].

### 2.3 The pore and permeation

The pore is formed between TM5 and TM6, with a short selectivity filter and a
gate at the cytosolic end of the TM6 bundle [R22, R24, R59]. Both
constrictions are recovered directly from the ligand-bound type-3 structure
({fig:channel_structure}c): the luminal one sits on the GGGVGD filter motif,
the cytosolic one on Phe2513 and Ile2517. The channel is not
a precision Ca<sup>2+</sup> filter in the sense that voltage-gated
Ca<sup>2+</sup> channels are: single-channel recordings from cerebellar
preparations established a large conductance with modest discrimination among
divalent and monovalent cations [R64]. Functionally this is appropriate — the
receptor releases Ca<sup>2+</sup> down a steep gradient from a store held at
high concentration, so throughput matters more than selectivity.

### 2.4 Long-range allosteric coupling

The defining structural problem of this receptor is distance. IP<sub>3</sub>
binds at the N-terminal core, roughly 100 Å from the gate — measured on that
structure, 103 Å along the pore axis, and 120 Å through space, because the
site also sits 62 Å out from the axis ({fig:channel_structure}a). The cryo-EM series
resolves the coupling path: ligand binding closes the clam-shell of the
binding core, that motion is transmitted through the armadillo solenoid, and
the C-terminal tail — which runs from beyond TM6 back up into the cytosolic
domain, forming a left-handed helical bundle at the four-fold axis and
contacting the N-terminal domains of *adjacent* subunits — communicates the
change to the gate [R22, R23, R59, R60].

Two consequences follow, and both recur later in this review. Because the
coupling element contacts neighbouring subunits, the channel's allostery is
intrinsically inter-subunit, which is the structural basis for the requirement
that all four ligand sites be occupied (§3.3) and for the dominant-negative
behaviour of heterozygous variants in a tetramer (§9). And because the
C-terminal tail is a load-bearing mechanical element rather than a tail in the
dispensable sense, variants there are not peripheral — they sit on the path
that makes the channel a channel.

![](figures/channel_structure.png)

**{fig:channel_structure}.** The channel measured from its own coordinates — PDB 6DQN, the human type-3 receptor with IP<sub>3</sub> bound at 3.33 Å [R24]. (**a**) Cα trace of the tetramer, coloured by Pfam domain, with the four bound ligands in red; the mushroom is 178 Å tall and the pore domain occupies the shaded band. The distance IP<sub>3</sub> acts across resolves into two numbers rather than one: the site sits **103 Å above the gate along the pore axis** — the quantity §2.4 quotes as "roughly 100 Å" — but also 62 Å out from it, so the through-space separation is **120 Å**. (**b**) The same particle down the four-fold axis; the deposit is C4-symmetric to 0.06 Å RMSD, so one subunit is committed and the other three are drawn by rotation. (**c**) Minimum distance from any heavy atom to the axis, along it. Two constrictions appear, and neither was given to the calculation: the luminal one falls on the GGGVGD selectivity-filter motif, and the cytosolic one on Phe2513 and Ile2517 — one helical turn apart, the arrangement reported for the IP<sub>3</sub>R1 gate. The filter is wide (5.1 Å to the nearest atom), which is what a high-conductance, weakly selective channel should look like (§2.3).

