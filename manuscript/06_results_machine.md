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
0.721–0.726), and the gate's five lining residues are **100 % identical
between all three human copies**. The most constrained whole domain is
RIH-associated (0.798–0.804).

**The ligand site is not in the domain named after the ligand.** None of the
ten IP₃ contacts measured on the 6DQN structure [R24] (≤ 4.5 Å) falls inside
PF08709, the signature called *Inositol 1,4,5-trisphosphate/ryanodine
receptor*; all ten sit in MIR and in the N-terminal RIH domain, which the
family shares with the ryanodine receptors — the same domains the crystal
structures of the isolated binding core identified [R05, R57] and the region
whose architecture is conserved between the two families [R58]. Those ten residues are more
constrained than the receptor as a whole in all three paralogues
(p = 0.013, 6 × 10⁻⁴, 4.2 × 10⁻³) and more constrained than the rest of the
domains that carry them in two of three (p = 0.070, 2.1 × 10⁻³,
2.5 × 10⁻²).

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
lie about fifty residues apart in the same Pfam domain (Fig. 6).

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
0.758 > shallow 0.684** (Extended Data Fig. 11) — against this task's own
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
the 3 negative controls (Extended Data Fig. 10). Calibration: two
conformations of the same protein score a median TM of 0.78, two IP₃
receptors 0.43, an IP₃ receptor against a ryanodine receptor 0.39, and the
negative controls 0.19 with 0 of 81 pairs reaching the 0.50 same-fold bar.
Per-domain confidence explains what predicted structures are good for here:
the IP₃-binding core is the best-modelled domain (median pLDDT 83.9) and the
pore the worst (71.0), against 69.5 outside annotated domains.

