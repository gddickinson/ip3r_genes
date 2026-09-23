## Methods

### Scope

The module comparison runs on three per-paralogue orthologue alignments of
vertebrate genome-sweep gene models (the construction is described in
{paper:constraint}): 260 *ITPR1*, 246 *ITPR2* and 262 *ITPR3* orthologues
enter the paired test after the coverage rule below, which drops three, two
and two. The pocket is measured in six IP₃-bound human IP₃R3 depositions
(6DQN, 8TKG, 8TKF, 8TKH, 7T3P, 8TLA). The enzyme panel covers 3,527
eukaryotic reference proteomes, the 763 vertebrate sets and 2,764 others. The
lineage test covers all 662 non-vertebrate reference proteomes that carry an
IP₃ receptor, of which 486 enter. Prokaryotes are outside the question and
were not swept.

### Module definitions

Each module has a primary definition derived from measurement and a
sensitivity definition taken from how the field draws the region. The primary
ligand core is the minimal contiguous span containing every measured IP₃
contact, carried to each paralogue through the alignment. The primary pore
module is the Pfam channel domain (PF00520 [R155]) less the luminal loop,
which is located geometrically as the contiguous run of channel residues
whose Cα lies beyond the membrane on the luminal side of the measured axial
span of 6DQN. The sensitivity definitions are the published binding core,
*ITPR1* residues 224–604 [R05], and PF00520 as InterPro [R156] draws it. Every
test was run under all four combinations. Each definition must hold what it
is defined to hold and none of what the other holds (all ten contacts and no
pore site; both filter and both gate residues and no contact), and a failure
raises rather than being written.

### Ligand shells

Every residue within 15 Å of IP₃ was measured in each deposition with a
minimal all-atom mmCIF reader, taking the first model and altloc A and
measuring within one subunit. A Cα trace puts an arginine several ångström
further from a phosphate than the side chain that binds it, so all atoms are
used. The positive control is a hard failure: all ten contacts published for
6DQN must be recovered at ≤ 4.5 Å in 6DQN or the stage stops. Residues were
assigned to four shells by their median distance across the six structures
(≤ 4.5, 4.5–8, 8–11.5 and 11.5–15 Å), and residues beyond 15 Å are absent
rather than pooled into an open bin. Shells were carried from IP₃R3 to the
other two paralogues by alignment, and a residue that does not transfer is
dropped and counted rather than assigned to a neighbour.

### Paired module comparison

Each orthologue's identity to its paralogue's human reference was computed
over each module's reference residues. An orthologue enters only if it
resolves at least half of both modules; the excluded ones are listed with the
module that lost them, because an N-terminally truncated gene model would
otherwise enter as an extreme ligand-core divergence. Core-minus-pore
differences were tested by Wilcoxon signed-rank and sign tests, and all tests
in the study were corrected together by Benjamini–Hochberg [R170]. Column-level
contrasts used the sequence-weighted [R171] Jensen–Shannon divergence [R178]
beside the composition-free share of the modal residue.

### Contacts and pocket

Contact constraint was tested by permutation (200,000 iterations), resampling
which positions carry the contact label against two backgrounds: the rest of
the core and the rest of the pocket. The gradient was tested by Spearman
correlation of distance with conservation across the 125 pocket residues.
Per-site rates were taken from FEL [R149] runs on the vertebrate codon
alignments, reporting the share of sites called purifying at q ≤ 0.05 per
shell and comparing module β by Mann–Whitney test.

### Enzyme panel

Profiles for both halves of the PI-PLC catalytic barrel (PF00387 and PF00388)
were searched with hmmsearch [R153] over each eukaryotic group's reference
proteomes at E ≤ 10, and the call taken by filtering the same output at
E ≤ 10⁻⁵. Hits were attributed to proteomes by a measured accession-to-proteome
map rather than by taxon. A proteome carries a PI-PLC only when one protein
holds both halves; each half is also counted separately. The vertebrate group
is the positive control for the search. Species with no reference proteome
are a separate stratum and never enter the test as absences.

### Lineage test and its control

Each non-vertebrate receptor was aligned pairwise with MAFFT [R138] to human
IP₃R3 and scored through the same module residues as the vertebrate test.
Enzyme-absent receptors were compared with enzyme-carrying ones pooled, within
phylum where both cells held three records, and matched: each absent record
against the median of up to three present records within 0.03 pore identity.
The ryanodine receptors, carrying the same pore and no IP₃ site, were scored
through the same instrument as the positive control. Power was simulated at
the matched design's size across a grid of true shifts and read against the
shift the control produces.

### Controls on the analysis itself

Forty-four constructed negative controls ran before any table was written,
each built to break one rule (a core missing a contact, a pore holding the
loop, a Cα-only reader, a permutation scoring the whole module, an
unreferenced proteome filed as an absence, a matched test that could not
fire), and ten deliberate code mutations were each caught by the control
responsible.
