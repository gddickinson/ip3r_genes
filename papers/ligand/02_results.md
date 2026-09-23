## Results

### Both modules are defined by measurement and checked against what they must not contain

We defined the primary ligand core as the minimal contiguous span containing
every IP₃ contact measured on the 6DQN structure [R24], residues 265–569 in
*ITPR1*, and the primary pore module as the Pfam channel domain less the
geometrically located luminal loop, residues 2363–2608 in *ITPR1* with the
loop removed ({fig:modules}a). Each definition was checked against residues it
must not hold: the core holds all ten contacts and no filter or gate residue,
and the pore both filter and both gate residues and no contact. A failing
definition stops the analysis rather than being written. The sensitivity
definitions are the published IP₃-binding core [R05], transferred to each
paralogue, and the channel domain as InterPro draws it, loop included. Both
modules are laid out at residue resolution in all three paralogues in
{fig:labelled_positions}.

### The pore is more conserved than the ligand core in two paralogues and level in the third

Every comparison is paired inside one sequence: each orthologue contributes
one core identity and one pore identity to its own paralogue's reference, and
only their difference enters the test, so a generally divergent sequence
cannot favour either module. The tests ran on 260 *ITPR1*, 246 *ITPR2* and 262
*ITPR3* orthologues.

The pore module is ahead in *ITPR1* (mean difference -0.0239; 223 orthologues
favour the pore against 32 the core; q = 4 × 10⁻³⁴) and in *ITPR3* (-0.0199;
218 against 42; q = 2 × 10⁻²⁸). *ITPR2* shows no difference (q = 0.13)
({fig:modules}b). The module that binds the ligand, and whose function
separates this family from its sister, is the less conserved of the two
across vertebrate orthologue sets.

### The answer reverses in all three paralogues when the luminal loop counts as pore

Leaving the luminal loop inside the pore module, as the database draws the
domain and as a comparison would be made by default, reverses the sign in
every paralogue. The ligand core then leads by 0.055 in *ITPR1*, 0.0445 in
*ITPR2* and 0.0558 in *ITPR3*, each at q < 10⁻³⁷ ({fig:modules}b). The
loop's fifty residues are the least conserved stretch of the receptor, and
they sit inside the same domain as the gate. Neither answer is wrong about its
own module. They are answers about different modules, and a comparison that
does not state where the pore ends cannot be read.

Swapping in the published binding core moves one paralogue: with the loop
excluded, *ITPR2* tips from level to a slight lead for the core (0.0049, not
significant), while
*ITPR1* and *ITPR3* keep the pore ahead. The luminal-loop boundary decides the
sign everywhere, and the core boundary decides it only where the two modules
were already level.

One further disagreement is reported rather than dropped. Scored column by
column with the Jensen–Shannon divergence [R178], the two modules do not
differ, while the composition-free share of the modal residue agrees with the
paired test. The divergence is measured against a background amino-acid
table, so a transmembrane module scores lower than a soluble one at equal
conservation; the column-level null is a property of that metric.

### Six independent structures place the same ligand pocket

We measured every residue within 15 Å of IP₃, using all atoms rather than Cα,
in six IP₃-bound human IP₃R3 depositions from several groups and
gating states [R23, R24, R59, R60]. Recovery of the ten published contacts in the
6DQN structure was required before anything downstream ran, and all six
depositions recover all ten. They also agree on two more residues, Ala276 and
Arg411, each within 4.5 Å in four of the six depositions and outside it in
6DQN. The pocket so defined holds 125 residues, binned into a contact shell of
twelve residues (≤ 4.5 Å) and three outer shells out to 15 Å, the outermost
holding 59 ({fig:shells}a).

### Selection holds a pocket rather than a contact set

The ten contacts are more constrained than the rest of the ligand core in all
three paralogues (mean deep-layer divergence 0.813 against 0.744 in *ITPR1*;
q = 0.048, 0.017 and 0.041), with the permutation resampling which positions
carry the contact label while the constraint values stay where they are.
Against every other residue within 15 Å of the ligand they are more
constrained in none (q ≥ 0.12).

Every shell out to 15 Å sits above the whole-protein mean (0.7366 in *ITPR1*), and the
step a contact-driven model predicts at 4.5 Å is absent ({fig:shells}b).
Conservation falls with distance across the pocket, strongly in *ITPR2*
(ρ = -0.436) and weakly in the other two, where the trend does not reach
significance. What selection holds at the ligand site is a neighbourhood
about 15 Å across, which a contact label cannot see.

### Per-site substitution rates reproduce the pocket but not the module ranking

Column conservation and a substitution rate on a tree are both called
constraint and are different measurements, so we asked the same two questions
of per-site rates estimated by FEL [R149] on the vertebrate codon alignments
({fig:omega}). Every contact site in every paralogue is under detectable
purifying selection. That share falls to 0.875, 0.75 and 0.825 in the third
shell of *ITPR1*, *ITPR2* and *ITPR3*. The fall is not monotone, since the
outermost shell sits slightly above the third (0.8966 in *ITPR1*), so the data
show a step down from the first two shells onto a floor rather than a
gradient.

The module ranking does not replicate on this axis, and in *ITPR1* it points
the other way: the pore's mean nonsynonymous rate (0.0446) exceeds the core's
(0.0251). In *ITPR3* the core evolves faster (q = 0.058), and *ITPR2* shows
nothing. The rate axis uses about a quarter as many sequences as the identity
axis and a rate that is zero at most sites, so we read the module result off
the identity axis and report the disagreement beside it.

### The enzyme that makes IP₃ is missing from 64 proteomes that carry the receptor

We swept both halves of the phosphoinositide-specific phospholipase C
(PI-PLC) catalytic barrel over 3,527 eukaryotic reference proteomes and
counted a proteome as carrying the enzyme only when one protein held both
halves. The 763 vertebrate proteomes are the search's positive control: 760
carry a PI-PLC, and none carries a receptor without one.

Across the whole sweep, 64 proteomes carry an IP₃ receptor and no PI-PLC
({fig:lineage}a): 22 fungal, 28 from other protists, 10 non-vertebrate
animal and 4 green-plant proteomes. They concentrate in the oomycetes and the
early-diverging fungi. Every one of them carries a full-length receptor, and a
gene set complete enough to hold a 2,700-residue multi-exon gene is complete
enough to hold a phospholipase.

### Matched on divergence, losing the enzyme does not relax the ligand core

The representative alignment holds only two of these lineages, which is too
few to test, although the ryanodine receptor control already fires there
(shift -0.0996 on six tips). We therefore aligned all 662 non-vertebrate
receptors from reference proteomes to the human reference through the same
module residues; 486 resolve both modules and enter.

Pooled, the 36 enzyme-absent receptors show a wider core-minus-pore gap than
the 450 others (p = 9.7 × 10⁻⁶). The pooled comparison is confounded by clade:
enzyme-absent proteomes sit at a median pore identity of 0.3582 against 0.5886
for the rest, and the paired difference itself correlates with divergence
({fig:lineage}b). Inside a phylum, where the clade is held constant, the
comparison is possible in only three strata. All 15 early-diverging fungal
(Mucoromycota) receptors lack the enzyme, so that phylum has no internal
control.

We therefore matched each enzyme-absent receptor to up to three
enzyme-carrying receptors within 0.03 pore identity of its own. All 36
matched. The median within-pair difference is -0.0064 (95 % CI -0.0161 to
0.0152), with 19 pairs on one side and 17 on the other (sign test p = 0.87).

The null is bounded. Measured through the same pairwise instrument at the
same divergence, the ryanodine receptors shift the paired difference by
-0.0624. The matched test detects a shift of 0.05 with power 0.99 and one of
0.03 with power 0.59 ({fig:lineage}c), and its interval excludes shifts larger
than about 0.016. Losing the enzyme that makes IP₃ does not relax the
receptor's IP₃-binding core, and a ligand-free core at this divergence would
have been detected.
