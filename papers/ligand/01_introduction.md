## Introduction

Inositol 1,4,5-trisphosphate (IP₃) releases calcium from the endoplasmic
reticulum by opening the IP₃ receptor, a tetrameric channel of about 2,700
residues per subunit [R01, R50]. Vertebrates carry three paralogues, *ITPR1*,
*ITPR2* and *ITPR3*. The ligand binds at the N-terminus, in a core formed by a
β-trefoil and an armadillo-like domain [R05, R57], and the signal travels
some two thousand residues to the pore at the C-terminus [R22, R24, R59].

The ryanodine receptors are the IP₃ receptor's sister family. They share the
receptor's architecture, its pore and every domain used to diagnose it
[R04, R58], and they do not bind IP₃; their gating is organised around other
ligands [R37]. The ligand core is therefore the one module whose function
separates the two families, which makes it the natural place to expect the
tightest constraint in the receptor. That expectation has not been measured,
for three reasons this paper addresses.

First, neither module is a database domain. Pfam names the ligand region
after the families it identifies rather than after the ligand, and its pore
domain includes a luminal loop that is the least conserved sequence in the
receptor. A comparison drawn on database spans therefore compares something
other than the ligand core against something other than the pore. We define
both modules by measurement, and we run every test under a second,
conventional definition of each.

Second, the ligand site is usually described as a set of contact residues.
A contact label discards what a structure measures, the distance of every
residue from the ligand, and it cannot say whether selection holds the
contacts or the neighbourhood around them. We measure the site as a distance,
all-atom, in six independent IP₃-bound depositions.

Third, if the core is constrained because it binds IP₃, lineages that no
longer make IP₃ should relax it. Which lineages those are is usually taken
from reading. We derive them by sweeping every eukaryotic reference proteome
for the enzyme that makes IP₃, and we test the prediction against a positive
control, the ryanodine receptors, whose ligand-free core shows what relaxation
looks like on the same instrument.

The question this paper answers is whether the IP₃-binding core is under
different constraint from the pore it gates ({fig:modules}). The per-residue
conservation layers it stands on are built and described in
{paper:constraint}. The comparison, its controls and its scope are all
measured here.
