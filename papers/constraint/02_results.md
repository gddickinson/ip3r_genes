## Results

### All three paralogues are under strong purifying selection

We recovered a coding sequence for all 57 vertebrate IP₃ receptor tips of
the representative alignment, 43 from protein-database cross-references and
14 from genome gene models, and accepted each only if its translation
reproduced the aligned protein; 24 codons that disagreed were masked rather
than kept ({fig:codon_alignment}). The trimmed codon alignment holds 2,459
codons. One-ratio models give ω = 0.02375 for ITPR1, 0.04303 for ITPR2 and
0.04152 for ITPR3 ({fig:omega}; Table 1). The least constrained paralogue
therefore accumulates about one non-synonymous substitution for every 23
synonymous ones across a 2,700-residue channel. Removing every genome gene
model moves no estimate by more than 0.0112 (curated-only ω 0.02118, 0.0318
and 0.03714), so the rates are not an artefact of the reconstructed models.

**Synonymous sites are saturated inside a single paralogue.** Across
pairwise comparisons within a paralogue, 161 of 171 ITPR1 pairs, and 373 of
420 pairs over all three, exceed a synonymous distance of 1.5, with median dS
of 5.721, 4.53 and 13.483 for ITPR1, ITPR2 and ITPR3 ({fig:codon_checks}a).
A single paralogue set spans sharks, teleosts and mammals, so saturation is
reached within the copies and not only between them. Every rate we report is
therefore from a tree-based model, and the pairwise matrix is used only as
this diagnostic.

| | ITPR1 | ITPR2 | ITPR3 |
|---|---|---|---|
| one-ratio ω | 0.02375 | 0.04303 | 0.04152 |
| one-ratio ω, curated sequences only | 0.02118 | 0.0318 | 0.03714 |
| RELAX k (test clade against the other two) | 9.3587 | 0.9085 | 0.8356 |
| branch-site ω₂ on the stem | 5.52937 | 999 (at bound) | 999 (at bound) |
| share of sites in the ω₂ class | 0.1102 | 0.14149 | 0.278 |

*Table 1. Rates of the three paralogues, read from the committed codon-model
tables. Two stems carry an ω₂ at the optimiser's upper bound, which is not an
estimate.*

### ITPR1 is the most constrained paralogue and its selection is intensified

Ranked from most to least constrained, the order is ITPR1, ITPR3, ITPR2.
Two-ratio branch models, each paralogue clade against the rest of the
family, separate every clade from its background after correction across the
whole family of 12 likelihood-ratio tests: ITPR1's foreground ω of 0.0241
against a background of 0.0432 (q = 2.9 × 10⁻⁶³), ITPR2's 0.0435 against
0.0317 and ITPR3's 0.0455 against 0.0308 ({fig:branches}a). A relaxation test,
which compares whole ω distributions rather than point estimates [R150],
gives the same ordering by a different statistic: selection on ITPR1 is
intensified relative to the other two (k = 9.3587), and relaxed on ITPR2
(k = 0.9085) and ITPR3 (k = 0.8356) relative to theirs ({fig:branches}b).
ITPR1 is also the paralogue that carries most of the family's pathogenic
missense record.

**Only the ITPR1 stem carries a branch-site rate the data determine.**
Branch-site model A asks whether a class of sites evolved faster than
neutrally on the stem of a paralogue, the interval in which a new copy was
most free to change [R146]. We ran it from four initial ω on every stem, and
3 of the 12 restarts converged below their own nested null, one on every
stem and from a different starting value each time ({fig:codon_checks}b).
Each of those is a local optimum a single-start run would have reported. All
three stem tests are significant after correction, but only one is
interpretable. On ITPR1 the foreground class has ω₂ = 5.52937 on a share of
0.1102 of sites, stable across the restarts that reach the best likelihood,
with 8 sites above a posterior of 0.95 [R147]. On the ITPR2 and ITPR3 stems ω₂
sits at codeml's bound of 999, and on ITPR2 restarts at the same likelihood
place it anywhere from 162 to 999: the stem has no synonymous signal left to
normalise a rate against. Those two tests are significant and their effect
sizes are not measured.

**No site in any paralogue is positively selected.** Within each paralogue,
M2a never gives its positive-selection class a non-zero share, and M8 fits
better than M7 in all three only because it adds a class at exactly ω = 1 on
under 1 % of sites, with 0, 0 and 1 sites reaching a posterior of 0.95. A
class of sites free to drift is not a class being driven.

### The gate and the selectivity filter are the most constrained elements

For per-residue conservation we took one locus per genome per paralogue from
the vertebrate genome sweep, giving 264 ITPR1, 249 ITPR2 and 265 ITPR3
orthologues ({fig:paralog_alignments}), and scored each column with a
sequence-weighted Jensen–Shannon divergence [R171, R178]. Every element was
tested against the same protein's own linkers, since a whole-protein mean
contains the element being tested. The linker means are 0.7256, 0.7209 and
0.7214.

On the modal-residue fraction, which uses no background model, the gate
and the selectivity filter are the two most conserved elements in all three
paralogues (gate 0.9846, 1.0 and 1.0), with the ITPR1 filter tied with the
RIH-associated domain ({fig:elements}; Table 2). On the divergence metric
the gate ranks first in ITPR2 and ITPR3 and second in ITPR1, and the filter
ranks second to fourth. The two are five and seven residues long, so their
tests against the linkers reach p < 0.05 only in ITPR2 and ITPR3 for the gate
(0.0154, 0.0224) and in ITPR3 for the filter (0.0237). The most constrained
whole domain is the RIH-associated domain (0.8037, 0.8016 and 0.7976), and
every named element but one sits above its protein's linker mean.

| element | ITPR1 | ITPR2 | ITPR3 |
|---|---|---|---|
| gate | 0.7994 | 0.822 | 0.822 |
| selectivity filter | 0.7866 | 0.7966 | 0.8114 |
| RIH-associated | 0.8037 | 0.8016 | 0.7976 |
| pore domain, luminal loop removed | 0.7532 | 0.7497 | 0.7489 |
| luminal loop | 0.4883 | 0.4521 | 0.4846 |
| linkers (control) | 0.7256 | 0.7209 | 0.7214 |

*Table 2. Mean Jensen–Shannon divergence per element on the deep layer.
Higher is more conserved.*

### A luminal loop is the least conserved sequence in the receptor

The one exception sits inside the pore domain. We located it from the
structure rather than from any annotation, as the contiguous run of channel
residues whose Cα lies beyond the membrane on the luminal side of the axial
span measured on 6DQN [R24]: residues 2401 to 2450 of ITPR3, transferred by
alignment onto ITPR1 and ITPR2, where 49 and 40 of its residues carry a
score. It scores far below every other element and below the
linkers in all three paralogues (0.4883, 0.4521 and 0.4846), and its
modal-residue fraction (0.5468, 0.5632 and 0.5236) says the same without a
background model. On the constraint profile along the pore-forming half, the
curve collapses in the loop and recovers on both sides of it
({fig:profile}).

**Leaving the loop inside the pore domain inverts the pore's result.** With
the loop included, the pore domain read as the least constrained named
element in ITPR1, below the linkers (JSD 0.697, p = 0.96), a strange result
for the pore of an ion channel. We checked three explanations. A metric
artefact was ruled out by the composition-free metric, which gave the same
picture. A gene-model artefact was ruled out by recomputing on the curated
sequences alone, where the loop still scored lowest (0.6191 in ITPR1). The
third explanation held: the domain's mean averaged the most and least
conserved stretches of the receptor. With the loop separated, the pore domain
clears the linker control in all three paralogues (p = 0.0054, 0.0011 and
0.0083). Because the loop was drawn on the membrane's geometry, that result
is not the profile explaining itself.

### The measured functional residues are conserved beyond their own domains

The ten residues within 4.5 Å of IP₃ in 6DQN all fall in MIR and the
N-terminal RIH domain, and none in PF08709, the Pfam signature named after
the ligand. Tested against the whole protein, they are more constrained in
all three paralogues (p = 0.013, 0.0006 and 0.0042). Tested against the rest
of the domains that carry them, the harder control, they are more constrained
in ITPR2 (p = 0.0021) and ITPR3 (p = 0.0253) and higher but not significant
in ITPR1 (p = 0.0697) ({fig:sites}a). The filter-lining and gate-lining sets
are two residues each, so we report their means and leave their tests to the
element level, where those elements are powered.

**The gate has not changed in any paralogue since the duplications.**
Identity between the human paralogues, which needs no alignment depth and no
conservation metric, is 1.0 across the five gate residues in all three pairs
({fig:sites}b). The luminal loop retains 0.1304, 0.2041 and 0.3111 identity,
against 0.6405, 0.6628 and 0.7026 for the whole protein. The two extremes of
divergence between the copies therefore lie about fifty residues apart in
the same domain. Measured on the trimmed representative alignment instead,
the whole-protein identities are 0.7526, 0.7728 and 0.812, which is the
difference trimming makes: it keeps 12 of the luminal loop's 51 residues.

**Per-site rates agree on the same coordinates.** Carried onto the human
proteins, FEL [R151] calls 5,766 sites purifying across the three paralogues
and one site diversifying, which is inside the procedure's own false
discovery budget. Every gate, filter and contact site in ITPR1 and ITPR3 is
called purifying; in the luminal loop the share falls to 0.2667, 0.3333 and
0.2667 of the codons the trimmed alignment kept, against 0.8603, 0.712 and
0.7733 over the whole protein, with median non-synonymous rates of 0.3275,
0.3902 and 0.3946 where elsewhere they are zero. Painted into the B-factor
column of 16 experimental and predicted structures, both layers show the same
pattern on the channel ({fig:constraint_on_channel}).

### Conservation separates pathogenic from benign variants, and depth is not the best layer

We placed 1,753 ClinVar missense records on the three canonical sequences
after checking each cited transcript's own translation against the canonical
residue by residue; none was dropped (1,006 ITPR1, 299 ITPR2 and 448 ITPR3
records). A positive control the harvest could have failed passed: both
residue-level variants the literature localises, ITPR3 p.Thr1424Met and
p.Arg2524Cys [R46, R48], were recovered. The record is mostly uncertain:
1,546 records are of uncertain significance (875, 261 and 410), and ITPR2's
entire pathogenic record is one position.

The four layers score slightly different sets of residues, so ranking them
on their own coverage would compare the sets as much as the layers. We
therefore scored them on the one fixed set that every layer scores, 44
pathogenic and 34 benign positions ({fig:classifier}a). The family-wide layer
of 128 IP₃ receptor tips is the best classifier (AUC 0.8723), then the
vertebrate layer (0.8536), then the deep within-paralogue layer (0.758), and
the shallow control is worst (0.6842). Depth was worth building, since the
shallow layer is what the representative alignment alone would have given.
Breadth was worth more: a position conserved across the family's eukaryotic
range and 500 million years of paralogue divergence discriminates better than
one invariant across 260 orthologues of a single gene.

**The ranking holds inside the genes that can be scored.** Without the
restriction to shared positions, the deep layer separates 45 pathogenic from
42 benign positions at AUC 0.7926, higher than on the fixed set, which is why
the restriction matters: each layer's own coverage flatters it differently.
Within ITPR1, which holds most of the labelled positions, the family layer
reaches 0.9013 on the shared set, and within ITPR3 it reaches 0.8533 on five
pathogenic positions. We report no per-gene score for ITPR2, whose single
pathogenic position cannot define a curve.

**The whole-protein control separates two effects.** Each layer also
separates pathogenic positions from the average residue, but by less (0.7574
for the family layer, against 0.8723 for the benign contrast). Part of the
signal is therefore benign calls concentrating in unconstrained positions,
not only pathogenic calls concentrating in constrained ones. Pathogenic
positions are enriched in the gating elements: 2 of ITPR1's pathogenic
positions fall in its five gate residues (odds ratio 48.9, p = 0.0019)
({fig:variants_on_structure}). ITPR1's pathogenic positions cannot be placed
on a structure at all, because neither its rat cryo-electron microscopy
reference nor the AlphaFold DB isoform model carries the human residue at
every shared position.

**Uncertain variants can be stratified, not called.** Against the labelled
distributions of its own gene on the family layer, 58 of 570 scored ITPR1
variants of uncertain significance, 11 of 196 in ITPR2 and 34 of 299 in ITPR3
sit at or above that gene's pathogenic median ({fig:classifier}b). That says
where a variant sits on an axis the labelled variants separate on. It is not
a pathogenicity call, and the ITPR2 row rests on a single pathogenic
position.
