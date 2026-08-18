# The IP3 receptor family — verified literature baseline

*Built in S0 (2026-08-18) from `docs/ip3r_background.md`. Every statement
below carries a citation to a source that was retrieved and read for this
purpose. The claim-by-claim audit trail — what was checked, what passed, what
did not — is `results/s0_baseline/lit_claims.tsv`; the bibliography is
`results/s0_baseline/references.tsv`; the rendered summary is
`results/s0_baseline/report.md`.*

**Reading rules for this file.**

- A statement with a citation has been verified against that source. Cite the
  source, never this file.
- A statement marked **⚠ qualified** is supported but is a review-level
  generalisation, not a measured quantity. It may be used to motivate, never
  to conclude.
- Claims that failed verification are **not softened here — they were struck
  from `ip3r_background.md`**. What replaced them is in §7.
- References are `R01`–`R51`, matching `references.tsv`.

---

## 1. What the protein is

The inositol 1,4,5-trisphosphate receptor is the endoplasmic reticulum's
ligand-gated calcium-release channel. Phospholipase C cleaves PIP2 to IP3;
IP3 binds an N-terminal core; the channel opens and releases Ca²⁺ from the ER
into the cytosol. Each step of that sentence has its own primary source:
IP3 releases Ca²⁺ from a non-mitochondrial intracellular store [R01]; the
receptor was cloned as the cerebellar P400 protein and functionally expressed
[R02]; the purified receptor alone conducts Ca²⁺ in reconstituted vesicles,
so it is the channel and not merely a regulator of one [R03]; and the
N-terminal IP3-binding core was solved with its ligand bound [R05], sitting
below an N-terminal suppressor domain that tunes affinity [R06].

**⚠ qualified.** The receptor is described as the origin of most
agonist-evoked cytosolic Ca²⁺ signals in non-muscle cells. The *hierarchy* is
established — quantal "puffs" are the elementary release events [R08], and
they summate into the oscillations and waves that carry the signal [R07] —
but "most" is a review-level generalisation with no measured fraction behind
it. Do not quote it as a number.

The receptor also governs ER-to-mitochondrion Ca²⁺ transfer. Close ER–
mitochondria contacts determine the mitochondrial Ca²⁺ response [R09]; the
physical linkage has been measured structurally and functionally [R10]; that
transfer sets the threshold for apoptosis [R11]. Importantly for this project,
**the three isoforms are not interchangeable in this role** — they differ in
how they regulate ER–mitochondrial contacts and local transfer [R12]. That is
a functional asymmetry between paralogs of exactly the kind S17 and S22 test.

### Gating

Gating is co-agonist and the Ca²⁺ dependence is biphasic: low Ca²⁺ activates,
high Ca²⁺ inhibits, giving the bell-shaped response curve that makes the
receptor a wave generator rather than a valve [R13]. IP3 alone does not open
it — IP3 and Ca²⁺ bind sequentially, which is what protects the cell from
spontaneous release [R14] — and **all four** IP3-binding sites of the
tetramer must be occupied before the channel opens [R15]. The last point
matters here beyond physiology: a channel that requires all four subunits to
be competent is a channel in which one bad subunit can poison the tetramer,
which is the mechanism behind the dominant-negative disease variants in §5.

Verified modulators, one primary source each: ATP [R16], PKA-dependent
phosphorylation [R17], Bcl-2 [R18], IRBIT [R19] and ERp44 [R20]. *PKC was
listed in the planning document but was not separately verified and has been
dropped.*

### Architecture

The subunit assembles as a homotetramer — purified as a tetrameric complex
from cerebellum [R21] and resolved by cryo-EM as a large cytosolic
"mushroom cap" over a small membrane domain [R22, R25]. The apo rat IP3R1
structure is at 4.7 Å and the paper calls the channel **1.3 MDa** [R22]; four
copies of the 2,758-residue human ITPR1 subunit (Q14643) give ≈1.25 MDa. The
planning document's "~1.2 MDa" is the same number rounded down — quote 1.3 MDa
with [R22], or the arithmetic with the accession, but not "~1.2 MDa" on its
own.

The pore lies between TM5 and TM6 with a short selectivity filter, and the
C-terminal tail after TM6 runs back up into the cytosolic domain, coupling
ligand binding at the top of the molecule to the gate at the bottom
[R22, R23, R24, R25]. This is the structural reason a C-terminal variant is
dangerous: it sits on the allosteric path, in a tetramer that needs all four
subunits. Both established dominant-negative disease variants are
C-terminal/pore-proximal — ITPR1 in Gillespie syndrome [R44] and ITPR3
p.Arg2524Cys [R48].

## 2. The genes

Human ITPR1 (Q14643, 2,758 aa), ITPR2 (Q14571, 2,701 aa) and ITPR3 (Q14573,
2,671 aa) are the three paralogs. **Their loci are no longer a literature
claim**: 3p26.1, 12p11.23 and 6p21.31 were re-derived from Ensembl and are
recorded in `results/s0_baseline/gene_structure.tsv`.

**⚠ qualified — tissue emphasis.** The paralogs are conventionally
distinguished as ITPR1 (cerebellar Purkinje cells, broad), ITPR2 (broad,
secretory epithelia, glia) and ITPR3 (epithelial/secretory, pancreas, liver).
This is supported qualitatively and from several directions: ITPR1 is the
Purkinje-cell P400 antigen [R21]; the subtypes are expressed tissue- and
development-specifically [R29] and localise differentially in brain and
periphery [R32], as reviewed across tissues [R33]; ITPR2 carries secretory
function in sweat-gland epithelium [R45]; ITPR3 predominates in pancreas,
liver and biliary epithelium [R34]. **No quantitative ranking is asserted
here** — GTEx/HPA quantification is S12's job, and the qualitative version
must not be used to support a claim about relative abundance.

ITPR1 carries alternatively spliced segments that distinguish neuronal from
peripheral isoforms: SI and SII were defined in mouse [R29]; the
SII-containing neuronal and non-neuronal forms differ in phosphorylation
[R30]; a third spliced segment, SIII, was identified in human ITPR1 [R31].

**The exon/span claim failed and was struck** — see §7.

## 3. The sister family, and why it matters here

The ryanodine receptors are not a separate problem to be ignored; they are
inside every search this project runs. This was visible from the beginning:
the IP3 receptor was reported in 1989 as a "putative receptor for inositol
1,4,5-trisphosphate **similar to ryanodine receptor**" [R04], in the same
issue as the cloning paper [R02]. The two families descend from a single
ancestral Ca²⁺-release channel, and the split predates animals: the Ca²⁺-
signalling machinery including both channel families is reconstructed to
before the animal/fungal divergence [R35], and both IP3R and RyR homologues
are present in parasitic protists [R36]; the RyR side of that history is
reviewed in [R37].

Human RYR1 (P21817, 5,038 aa), RYR2 (Q92736, 4,967 aa) and RYR3 (Q15413,
4,870 aa) carry **all four** of the ITPR-diagnostic Pfam signatures — PF02815,
PF08709, PF01365, PF08454 — verified per protein in
`results/s0_baseline/pfam_architecture.tsv`.

S0 measured what that costs. The zebrafish query the planning document quoted
as evidence for "four IP3 receptors" returns 109 records across 10 gene
symbols, and **53 of the 109 (49%) are ryanodine receptors**. That is a
measured floor on the contamination any Pfam-driven enumeration inherits, and
it is why D14 requires a positive ITPR/RYR call at every stage. In this
particular dataset the length band separates the two families cleanly
(≤2,819 aa vs ≥4,863 aa), but that is a fact about zebrafish annotation
quality, not a rule — hence "filter, never evidence".

**What is *not* established**, and has been retagged `[open]`: that the two
families' expansions to three vertebrate paralogs were *independent* events.
The three-paralog state of each family is a database fact [R32], but the
independence of the two triplications has no verified primary source and is
precisely open question **Q2**. Assuming it in the background document would
have made this project's own conclusion (S7, S13, S16) an input.

## 4. What the databases hold

Re-derived 2026-08-18 and unchanged from the planning snapshot; see
`results/s0_baseline/report.md` §2 for the tables. The headline for the
project is the Q1 anomaly, which survived re-derivation exactly: **40
Viridiplantae and 41 Fungi proteins carry PF08709**, the IP3-binding core,
in a family textbooks say plants and fungi lack — while *Arabidopsis
thaliana* and *Saccharomyces cerevisiae* have none. Whether those records are
real genes, mis-annotations or contamination is S2/S20's question. It is worth
noting that [R35] finds lineage-specific expansion of Ca²⁺ channels in basal
fungi, so a real-gene answer is not a priori absurd.

## 5. Disease — the clinical variant set

**ITPR1.** Heterozygous *deletions* cause spinocerebellar ataxia type 15: the
deletion was identified in humans and the corresponding ataxic mouse [R38],
then replicated in independent Japanese [R39] and European [R40] cohorts.
*Missense* variants cause SCA29, an autosomal dominant congenital
non-progressive ataxia [R41], with the phenotype consolidated in a case series
[R42].

ITPR1 also causes **Gillespie syndrome** (partial aniridia, cerebellar ataxia,
intellectual disability) — but **the planning document's version of this was
incomplete and has been corrected**. Gillespie syndrome arises from *both*
biallelic recessive ITPR1 variants *and* de novo heterozygous variants [R43];
only the latter act dominant-negatively, through a restricted repertoire of
channel-domain mutations [R44]. Both mechanisms must be carried into S17:
a constraint analysis that assumes only dominant-negative C-terminal variants
would mis-weight the recessive half of the variant set.

**ITPR2.** Homozygous loss of function causes autosomal recessive isolated
anhidrosis [R45] — with three qualifications that the one-line version loses.
The evidence is a *single* consanguineous family (five affected members); the
variant is a homozygous **missense** in the pore-forming region that abolishes
Ca²⁺ release, i.e. functionally null rather than a null allele; and the human
genetics is supported by *Itpr2*⁻/⁻ mice showing reduced sweat secretion. It
is a well-supported claim resting on a narrow base.

**ITPR3.** Dominant variants cause demyelinating Charcot–Marie–Tooth disease
(CMT1J). The gene–disease relationship was established with genetic and
functional evidence [R46] and is now strong: a recurrent p.Thr1424Met was
found in 33 individuals from nine unrelated families [R47], a further family
was reported independently [R51], and a spontaneous canine *ITPR3* nonsense
variant reproduces the neuropathy [R49].

The ITPR3 phenotype is **broader than the planning document stated**: the
recurrent de novo p.Arg2524Cys causes a complex multisystemic disease with
immunodeficiency, of which CMT is one component [R48]. This is an addition,
not a correction — but S17 should treat the ITPR3 variant set as multisystem,
not purely neuropathic.

Taken together this remains an unusually good clinical panel: three paralogs,
three phenotypes, and dominant, dominant-negative and recessive mechanisms all
represented. It makes S17 a test with an answer rather than a description.

## 6. Structures

Multiple cryo-EM structures of the tetramer exist. Rat IP3R1: the apo
structure at 4.7 Å [R22], ligand-bound states revealing the allostery [R23],
and the channel in a lipid bilayer [R25]. Human IP3R3: Ca²⁺- and IP3-bound
states [R24]. RyR1 was solved to near-atomic resolution three times in 2015
[R26, R27, R28] and remains among the best-resolved large channels available.
The structural picture as a whole, and how it maps onto function, is reviewed
in [R50].

**No PDB IDs are recorded in this file**, deliberately. S11's first step is an
RCSB query, so the reference panel is discovered rather than remembered; a
wrong four-character code that silently fetches a different protein is exactly
the error class that survives to publication.

## 7. What changed in `ip3r_background.md`

Four of the nineteen audited claims did not survive as written. All four edits
have been applied to the background document.

| Claim | Verdict | What changed |
|-------|---------|--------------|
| Loci 3p26.1 / 12p11.23 / 6p21.31 | `upgraded_db` | Re-derived from Ensembl; retagged `[lit]` → `[db]`. All three exact. |
| "~58–60 exon gene spanning hundreds of kb" | `corrected` | **False as written.** Measured: ITPR1 62 exons / 354 kb, ITPR2 57 / 498 kb, ITPR3 58 / **76 kb**. Range is 57–62 exons; ITPR3 does not span hundreds of kb. Struck and replaced with the measured values. |
| "both independently expanded to three vertebrate paralogs" | `downgraded_open` | The independence of the two triplications is unverified and is question Q2. Retagged `[lit]` → `[open]`. |
| Gillespie syndrome as dominant-negative | `corrected` | Incomplete. Both recessive biallelic and de novo dominant-negative mechanisms exist [R43, R44]. Amended. |

Three further claims are kept but **qualified** in this file and marked
`⚠ qualified` in the background document: "most agonist-evoked Ca²⁺ signals"
(§1), the tissue-emphasis summary (§2), and the ITPR2 anhidrosis evidence base
(§5).

One incidental finding from the exon/span correction is worth carrying
forward: **genomic span varies 6.5-fold across the three paralogs (76 kb to
498 kb) while protein length varies by 3%**. That is a target for S21, not a
footnote.
