# The IP3 receptor family — the project's literature and database baseline

*Written 2026-08-18 by the project-setup session. Every statement carries a
tag saying how far you may trust it:*

- **[db]** — checked against a live database on 2026-08-18, with the query
  recorded. Safe to build on; re-derive the number if you cite it.
- **[lit]** — standard textbook/review-level knowledge, not yet verified
  against a primary source in this repo. **S0 verifies these.** Do not put a
  `[lit]` statement in a manuscript until it has a checked citation.
- **[open]** — an open question. These are what the project is for.

---

## 1. What the protein is

**[lit]** The inositol 1,4,5-trisphosphate receptor (IP3R, gene symbol
`ITPR`) is the endoplasmic reticulum's ligand-gated calcium-release channel.
Phospholipase C cleaves PIP2 to IP3; IP3 binds the receptor's N-terminal
core; the channel opens and releases Ca(2+) from the ER into the cytosol.
It is the origin of most agonist-evoked cytosolic Ca(2+) signals in
non-muscle cells — puffs, oscillations and waves — and of the ER-to-
mitochondrion Ca(2+) transfer at MAM contact sites that sets the threshold
for apoptosis.

**[lit]** Gating is co-agonist: IP3 alone does not open the channel, and the
Ca(2+) dependence is biphasic (low Ca(2+) activates, high Ca(2+) inhibits),
which is what makes the receptor a wave generator rather than a valve. ATP,
PKA/PKC phosphorylation, Bcl-2-family proteins, IRBIT and ERp44 all modulate
it.

**[db]** Architecture, from InterPro for human ITPR1 (`Q14643`, 2,758 aa):

| Pfam | Name | Role |
|------|------|------|
| PF02815 | MIR | N-terminal, in the suppressor region |
| PF08709 | Ins145_P3_rec | the IP3-binding core (beta-trefoil) |
| PF01365 | RYDR_ITPR (RIH) | RyR and IP3R homology domain |
| PF08454 | RIH_assoc | RIH-associated |
| PF00520 | Ion_trans | the six-TM pore module, C-terminal |

**[lit]** The subunit assembles as a homotetramer (~1.2 MDa) with a large
cytosolic "mushroom cap" over a small membrane domain; the pore lies between
TM5 and TM6 with a short selectivity filter. The C-terminal tail after TM6
runs back up into the cytosolic domain and couples ligand binding to the
gate — which is why C-terminal variants can be dominant-negative in a
tetramer.

## 2. The genes

**[db]** The three human paralogs, with lengths confirmed from UniProt:

| Gene | UniProt | Length | **[lit]** locus | **[lit]** tissue emphasis |
|------|---------|--------|-----------------|---------------------------|
| ITPR1 | Q14643 | 2,758 aa | 3p26.1 | cerebellar Purkinje cells, broad |
| ITPR2 | Q14571 | 2,701 aa | 12p11.23 | broad; secretory epithelia, glia |
| ITPR3 | Q14573 | 2,671 aa | 6p21.31 | epithelial/secretory, pancreas, liver |

**[lit]** Each is a ~58–60 exon gene spanning hundreds of kb. ITPR1 carries
the alternatively spliced SI, SII and SIII segments that distinguish
neuronal from peripheral isoforms.

**[db]** Zebrafish carries four: `itpr1a`, `itpr1b`, `itpr2`, `itpr3`
(2,634–2,819 aa; UniProt query `taxonomy_id:7955 AND xref:pfam-PF08709`).
The `itpr1a`/`itpr1b` pair is the signature of a teleost-specific
duplication — the 3R test in S16 has a real target.

## 3. The sister family, and why it matters here

**[db]** The ryanodine receptors are not a separate problem to be ignored —
they are *inside every search this project will run*. Human RYR1 `P21817`
(5,038 aa), RYR2 `Q92736` (4,967 aa), RYR3 `Q15413` (4,870 aa) carry **all
four** of the Pfam domains listed above (checked for RYR1 on 2026-08-18),
plus RyR-specific ones: PF02026 (RyR domain), PF06459 (RyR TM 4-6), PF21119
(junctional solenoid repeat), PF00622 (SPRY).

Three consequences, and they are load-bearing:

1. A Pfam enumeration of `PF08709` returns IP3Rs **and** RyRs. Every census
   step needs a positive ITPR/RYR call, not an assumption.
2. Length separates them cleanly in the data seen so far: ~2,630–2,830 aa
   versus ~4,850–5,110 aa. `src/utils/family.py` sets the band at
   2,000–3,600 aa for exactly this reason. Length alone is not a
   phylogenetic argument, so it is a filter, never the evidence.
3. RyR is the natural **outgroup** for rooting the ITPR tree, which is a
   better-conditioned root than any invertebrate ITPR.

**[lit]** The IP3R/RyR split predates animals; both families descend from a
single ancestral Ca(2+)-release channel, and both independently expanded to
three vertebrate paralogs.

## 4. What the databases actually hold (the honest starting scale)

**[db]** InterPro protein counts per Pfam signature, 2026-08-18:

| Pfam | proteins |
|------|----------|
| PF08709 Ins145_P3_rec | 12,338 |
| PF01365 RIH | 13,177 |
| PF08454 RIH_assoc | 12,062 |
| PF02815 MIR | 23,453 (also in O-mannosyltransferases — not family-specific) |
| PF00520 Ion_trans | 206,115 (generic channel domain) |

**[db]** Taxonomic distribution of PF08709 (and PF01365 in brackets):

| Taxon | proteins |
|-------|----------|
| Metazoa | 12,149 (12,730) |
| SAR (stramenopiles/alveolates/rhizaria) | 309 (215) |
| Discoba (incl. kinetoplastids) | 48 (35) |
| Fungi | 41 (16) |
| Viridiplantae | 40 (56) |
| Amoebozoa | 7 (11) |
| Bacteria | 0 (1) |
| Archaea | 0 (0) |
| *Arabidopsis thaliana* | 0 (0) |
| *Saccharomyces cerevisiae* | 0 (0) |
| *Paramecium tetraurelia* | 21 (14) |
| *Trypanosoma brucei* | 1 (1) |
| *Chlamydomonas reinhardtii* | 1 (1) |
| *Dictyostelium discoideum* | 0 (1) |

These are API snapshots for planning, not results. **Re-derive them in S2
before any of them is quoted anywhere.** The queries are
`GET /api/protein/uniprot/entry/pfam/{PF}/taxonomy/uniprot/{taxid}` at
`https://www.ebi.ac.uk/interpro`.

## 5. Disease — the clinical variant set

**[lit]**, all to be verified against ClinVar/OMIM in S0 and S17:

- **ITPR1** — heterozygous deletions cause **spinocerebellar ataxia type 15
  (SCA15)**; missense variants cause **SCA29**, a congenital non-progressive
  ataxia; and specific, often C-terminal/pore-proximal variants cause
  **Gillespie syndrome** (partial aniridia, cerebellar ataxia, intellectual
  disability), which behaves as dominant-negative in the tetramer.
- **ITPR2** — homozygous loss of function causes **autosomal recessive
  isolated anhidrosis** (inability to sweat).
- **ITPR3** — dominant variants are reported in **Charcot–Marie–Tooth
  disease** / demyelinating neuropathy.

This is a better clinical panel than most families offer: three paralogs,
three different phenotypes, dominant and recessive mechanisms both present.
It makes the constraint task (S17) a test with an answer, not a description.

## 6. Structures

**[lit]** Multiple cryo-EM structures of the tetramer exist for rat IP3R1
(apo and ligand-bound states) and human IP3R3, and RyR1 is one of the
best-resolved large channels in the PDB. **Do not hard-code PDB IDs from
this document.** S11's first step is an RCSB query for the family, so the
reference panel is discovered rather than remembered; a wrong four-character
code that silently fetches a different protein is exactly the kind of error
that survives to publication.

## 7. The open questions this project can answer

- **[open] Q1 — Range.** What is the family's true taxonomic distribution?
  The database snapshot above says the family is overwhelmingly metazoan,
  yet holds 40 Viridiplantae and 41 Fungi records in a family that textbooks
  say plants and fungi lack, and zero in *Arabidopsis* and *S. cerevisiae*.
  Are those records real genes, mis-annotations, or contamination — and are
  the famous absences facts about **genomes** or about **proteome
  databases**? (S2, S3, S20)
- **[open] Q2 — Origin.** Are ITPR1/2/3 2R ohnologs? Are teleost
  itpr1a/itpr1b from 3R? Where does the IP3R/RyR split fall? (S7, S13, S16)
- **[open] Q3 — Fates.** Across a declared vertebrate genome scope, is any
  paralog lost, how often, and independently or once? Is any surviving as a
  decaying pseudogene? (S5, S15)
- **[open] Q4 — Records.** How often is a real ITPR locus missing,
  fragmentary, split across models, unnamed, or filed under the wrong
  paralog — with the wrong-paralog risk raised by a sister family sharing
  every diagnostic domain? (S10, S18, S21)
- **[open] Q5 — Machine.** Do the pathogenic variants sit in the most
  constrained parts of the channel, and is the constrained core the same in
  all three paralogs? (S9, S11, S17)
- **[open] Q6 — Anything unnamed.** Is there a fourth vertebrate ITPR, or an
  ITPR in a lineage reported to lack one? This is the app's discovery path
  applied to a family where the answer is genuinely unknown. (S2, S3, S5)

## 8. Reference material to add here

This project has no equivalent of a single preprint to test. The literature
baseline is a review-level synthesis, and S0 builds it: pull 15–25 primary
references covering sections 1–6, verify every `[lit]` statement above
against one of them, and record the verified set as
`docs/ip3r_review_2026.md` with a citation on every claim. Anything that
fails verification gets struck from this file, not softened.
