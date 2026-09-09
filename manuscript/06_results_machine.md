### What the receptor cannot change, and what that says about its variants

Per-residue constraint was computed in each human paralogue's own numbering
on four nested alignment layers: a **deep** layer of 249–265 orthologues of
that one paralogue drawn from the genome sweep, a **vertebrate** layer, a
**family-wide** layer of all 128 IP₃ receptor tips, and a **shallow** layer
of that paralogue's representative sequences kept as the control for what
the deep sets bought. Scores are sequence-weighted before any column
statistic and gaps are treated as missing data, not as a twenty-first state.

**The gate and the filter are the least changeable parts of the receptor,
and the gate is identical between the paralogues.** On both a divergence
metric and a composition-free one, the gate and the selectivity filter rank
at the top of every element in all three paralogues (gate mean JSD 0.799,
0.822, 0.822; filter 0.787, 0.797, 0.811, against a linker control of
0.721–0.726; Extended Data Fig. 12), and the gate's five lining residues are
**100 % identical between all three human copies**. The most constrained whole domain is
RIH-associated (0.798–0.804).

**A 50-residue luminal loop is the least conserved sequence in the
receptor, and it sits about fifty residues from the most conserved.**
Unresolved, the channel Pfam domain scored *below* the linker control
(p = 0.96), which survived both a metric control and a gene-model control.
The cause is a stretch of the channel domain that no annotation carries: we
located it geometrically, as the contiguous run of channel residues whose
Cα lies beyond the membrane on the luminal side of the measured axial span
(*ITPR3* 2401–2450). It is the least constrained element in the receptor
(mean JSD 0.452–0.485 against 0.721–0.726 for the linkers), the least
conserved between paralogues (identity 0.13–0.31 against 0.64–0.70
whole-protein), and the one region where site-wise selection is weak
(FEL median β 0.33–0.39, 27–33 % of sites purifying against 71–86 %
protein-wide). With it separated out, the pore module clears the linker
control (p = 0.005). The most and least constrained sequence in the receptor
lie about fifty residues apart in the same Pfam domain (Fig. 6), and the
next result turns on which of them counts as pore.

**The ligand site is not in the domain named after the ligand, and what
selection holds there is a pocket rather than a contact set.** None of the
ten IP₃ contacts measured on the 6DQN structure [R24] (≤ 4.5 Å) falls inside
PF08709, the signature called *Inositol 1,4,5-trisphosphate/ryanodine
receptor*; all ten sit in MIR and in the N-terminal RIH domain, which the
family shares with the ryanodine receptors — the same domains the crystal
structures of the isolated binding core identified [R05, R57] and the region
whose architecture is conserved between the two families [R58]. To ask what
those residues are special *against*, we re-measured the site as a distance
rather than a label: every residue within 15 Å of the ligand, all-atom, in
six independent IP₃-bound human IP₃R₃ depositions. All six recover all ten
of the published contacts — the positive control the measurement was
required to pass — and agree on **two more**, Ala276 and Arg411, which sit
inside 4.5 Å in a majority of structures and outside it in 6DQN alone. The
ten contacts are more constrained than the receptor as a whole in all three
paralogues (p = 0.013, 6 × 10⁻⁴, 4.2 × 10⁻³) and more constrained than the
rest of the binding core in all three (q = 0.048, 0.017, 0.041). Against
every *other* residue the structure places within 15 Å of the ligand they
are more constrained in **none** (q ≥ 0.12). Every shell out to 15 Å sits
above the whole-protein mean and there is no step at 4.5 Å, and site-wise
selection replicates it on an independent axis: **100 % of contact sites in
all three paralogues are under detectable purifying selection**, and that
share falls by 0.13, 0.25 and 0.18 to each paralogue's lowest shell
(Extended Data Fig. 13). The fall is not monotone — the outermost shell sits
slightly above the third in all three — so what the data show is a step down
from the ligand's first two shells onto a floor rather than a gradient
running all the way out. The constrained unit is a neighbourhood about
fifteen ångström across, which a contact/not label cannot see.

**The module that names the family is the less constrained of the two, and
the answer reverses on a boundary nobody states.** The IP₃-binding core is
the one part of this receptor the ryanodine receptors do not share
functionally, which makes it the obvious place to expect the tightest
constraint. It is not. Measured per orthologue and paired inside each
paralogue's own deep alignment — one core number and one pore number per
sequence, 246–262 sequences per paralogue, so a generally divergent tip
contributes a divergent core *and* a divergent pore and only the difference
enters — the **pore module is ahead by about two identity points in *ITPR1*
(−0.0239, 223 tips to 32, q = 4 × 10⁻³⁴) and *ITPR3* (−0.0199, 218 to 42,
q = 2 × 10⁻²⁸), and level in *ITPR2*** (−0.0005, q = 0.13). That result
depends entirely on one boundary. Leave the 50-residue luminal loop inside
PF00520, which is how InterPro draws the domain and therefore how the
comparison would be made by default, and **all three paralogues flip**: the
ligand core wins by five identity points at q < 10⁻³⁷ in every one. Neither
answer is wrong about its own module; they are answers about different
modules, and any version of this comparison that does not declare where the
pore stops is not interpretable.

**Losing the enzyme that makes the ligand does not relax the site that binds
it.** The lineages with a reduced upstream pathway were derived rather than
read: phospholipase C profiles for both halves of the catalytic barrel were
swept over all 3,527 eukaryotic reference proteomes, and **64 carry an IP₃
receptor and no phosphoinositide-specific phospholipase C** — concentrated
in the oomycetes and the early-diverging fungi, every one of them holding a
full-length receptor, in a gene set complete enough to have held the enzyme.
Pooled, their ligand core looks relaxed (p = 10⁻⁵), and that is a clade
artefact: they sit at a median pore identity of 0.358 against 0.589 for the
rest, so the comparison is largely distant against near. Matched on
divergence — each record scored against phospholipase-carrying records
within 0.03 pore identity of itself, all 36 matched — **the effect is gone**:
median within-pair difference −0.0064 (95 % CI −0.016 to +0.015), 19 tips
against 17, p = 0.87. This is a bounded null rather than an absent
measurement. Measured through the same instrument at the same divergence,
the ryanodine receptors — same pore, same diagnostic domains, no IP₃ site —
shift by −0.062, and the matched test detects a shift that size essentially
always and one half that size 59 % of the time (Extended Data Fig. 13).

**Constraint separates the clinical variants, and the deepest layer is not
the best classifier.** 1,753 ClinVar missense records — the family's disease record runs from
spinocerebellar ataxia and Gillespie syndrome at *ITPR1* [R38, R43, R44]
through anhidrosis at *ITPR2* [R45] to neuropathy at *ITPR3* [R46, R47] —
were placed on the canonical sequences with every cited transcript's own coding sequence
translated and checked residue by residue; none was dropped, and both of the
literature's residue-level citations were recovered by a positive control
the harvest could have failed. **1,546 (88 %) are variants of uncertain
significance**, and *ITPR2*'s entire pathogenic missense record is a single
variant. Scored as a classifier of pathogenic/likely-pathogenic against
benign/likely-benign on the one fixed set of 44 versus 34 positions that
every layer scores, the ranking is **family 0.872 > vertebrate 0.854 > deep
0.758 > shallow 0.684** (Extended Data Fig. 12) — against this task's own
design, the 249–265-orthologue sets are the second-best of four, and the
whole-family layer that ignores paralogue identity is the best. Site-wise
selection agrees with constraint and adds nothing new: FEL finds 5,766
purifying sites across the three paralogues and **one** diversifying site,
inside its own false-discovery budget.

**Structures confirm the family call by shape, and show why they cannot do
more.** AlphaFold DB holds a usable (≥ 95 % covered) model for 20.9 % of the
census's 8,319 UniProt-shaped IP₃ receptor records — but of the 5,861 at or
above the family's 2,000 aa floor, **13 (0.2 %)** do, and of the 134
representatives every alignment and tree here stands on, **9** do. The
median modelled record is 392 aa and the median unmodelled one 2,674 aa: the
coverage is concentrated on fragments. A panel of 29 structures — 5
cryo-EM references resolved by query and assigned by the census rather than
by entry title, a 5-state conformational panel, 3 size-matched negative
controls and 16 predicted models — gives 20 of 20 agreement between the
structural family call and the census, and no family call at all for any of
the 3 negative controls (Extended Data Fig. 14). Calibration: two
conformations of the same protein score a median TM of 0.78, two IP₃
receptors 0.43, an IP₃ receptor against a ryanodine receptor 0.39, and the
negative controls 0.19 with 0 of 81 pairs reaching the 0.50 same-fold bar.
Per-domain confidence explains what predicted structures are good for here:
the IP₃-binding core is the best-modelled domain (median pLDDT 83.9) and the
pore the worst (71.0), against 69.5 outside annotated domains.
