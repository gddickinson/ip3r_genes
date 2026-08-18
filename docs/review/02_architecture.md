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
gate at the cytosolic end of the TM6 bundle [R22, R24, R59]. The channel is not
a precision Ca<sup>2+</sup> filter in the sense that voltage-gated
Ca<sup>2+</sup> channels are: single-channel recordings from cerebellar
preparations established a large conductance with modest discrimination among
divalent and monovalent cations [R64]. Functionally this is appropriate — the
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
change to the gate [R22, R23, R59, R60].

Two consequences follow, and both recur later in this review. Because the
coupling element contacts neighbouring subunits, the channel's allostery is
intrinsically inter-subunit, which is the structural basis for the requirement
that all four ligand sites be occupied (§3.3) and for the dominant-negative
behaviour of heterozygous variants in a tetramer (§9). And because the
C-terminal tail is a load-bearing mechanical element rather than a tail in the
dispensable sense, variants there are not peripheral — they sit on the path
that makes the channel a channel.
