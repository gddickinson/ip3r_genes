# Analysis catalogue — what this data will be able to answer

Enumerating every ITPR gene is the project's primary aim (S0–S7). But the
harvest chain produces data assets that answer much broader questions about
the family — its **origin, range, fates, constraint and function** — plus a
few meta-scientific questions about how gene families are catalogued at all.

This file is the menu behind ledger rows S15–S22. Each entry states the
**question**, the **inputs** it needs, the **method**, the **deliverable**,
and — importantly — the **caveat that could invalidate it**.

Feasibility, as of project setup (2026-08-18), when no task has yet run:
✅ ready · 🟡 needs an upstream task · 🔴 needs data we do not have yet.
**Everything is 🟡 or 🔴 today.** Update the marks as tasks complete.

---

## The data assets (what each upstream task leaves behind)

| Asset | What it will hold | Path | From |
|---|---|---|---|
| Census v2 | every UniProt protein carrying a family Pfam, with an ITPR/RYR call and reason per record | `results/census_v2/` | S2 |
| Census v3 | v2 + profile-HMM hits over the vertebrate reference proteomes | `results/census_v3/` | S3 |
| Census v4 | v3 + genome-only miniprot gene models | `results/census_v4/` | S5 |
| Census v5 | v4 + the non-vertebrate sweep, lineage columns on every row | `results/census_v5/` | S20 |
| Genome ledger | per genome × paralog status, coordinates, coverage, identity, **frameshifts, stop codons**, contig-edge/N-run flags, annotation overlap, rescue regions with per-bait bitscores | `results/genome_ledger.tsv` | S5 |
| Gene models | **exon/CDS coordinates** for every locus in every swept genome | `<data_root>/genome_sweep/<acc>/miniprot.gff` | S5 |
| Neighbourhoods | every annotated gene interval per assembly | `<data_root>/genome_sweep/<acc>/genes_slim.tsv` | S5 |
| Assembly metadata | level, **scaffold/contig N50**, size, release date, annotation source (`GCF_`/`GCA_`) | `results/genome_manifest.tsv` | S4 |
| HMM sweep | per-protein bitscores + profile coverage, ITPR-vs-RYR margins, convergence series | `results/hmm_sweep/` | S3, S20 |
| Controls | measured recall/specificity of the scorer, including the RyR decoys | `results/benchmark_controls/` | S1 |
| Codon alignments | translation-validated CDS per tip, PAL2NAL codon alignment | `results/selection/` | S9 |
| Structures | AFDB + cryo-EM panel, TM-scores, per-domain pLDDT | `results/structures/` | S11 |

---

## S15 — Loss dynamics: what happened to each paralog

**A1. Ancestral-state reconstruction of presence/absence.** 🟡 (needs S5, S13)
- *Question.* Is any of ITPR1/2/3 lost in any vertebrate lineage, how often
  independently, and where? A family this central to calcium signalling may
  turn out never to lose a paralog — **that is a result**, and the analysis
  is the same either way.
- *Inputs.* The genome ledger + a dated species tree (D15).
- *Method.* Code each cell present / absent / **missing**, then Dollo
  parsimony plus ML binary-trait models (ER / ARD / irreversible).
  `assembly_gap` and ambiguous traces are **missing data, not absence** —
  that coding choice is the whole analysis, so run the full coding ×
  evidence × branch-length × contiguity sensitivity matrix.
- *Deliverable.* Loss-map figure and an event table with per-event support.
- *Caveat.* Absence inherits assembly quality. Re-run excluding genomes
  below the contiguity floor and report whether the count is stable (A15).

**A2. ORF-integrity screen.** 🟡 (needs S5)
- *Question.* Which surviving loci are intact ORFs and which are decaying?
- *Inputs.* Per-locus frameshift and stop counts from the sweep.
- *Method.* Score integrity per locus; **paired within-genome** comparison
  against that genome's other paralogs (D16); calibrate the disabling-density
  threshold on the controls; correlate density against contig N50.
- *Caveat.* miniprot frameshifts also arise from assembly error. No
  pseudogenisation claim rests on this alone — confirm against a second
  assembly or RNA-seq (S10, S12).

## S16 — Duplication history: where ITPR1, ITPR2 and ITPR3 came from

**A3. The 2R test.** 🟡 (needs S7, S8)
- *Question.* Are the three paralogs ohnologs from the vertebrate genome
  duplications, or older/younger tandem or segmental events? The ryanodine
  receptors also have exactly three vertebrate paralogs — if both families
  triplicated at the same time, that is strong 2R evidence; if not, the
  coincidence needs another explanation. **This comparison is available to
  this project and to almost nobody else**, because the searches return both
  families anyway.
- *Inputs.* Human paralogy map (pinned Ensembl BioMart archive), the ITPR and
  RYR locus coordinates, S8 flank sets.
- *Method.* Cross-window paralogy at ±10/20/30 genes against a permutation
  null of **real** genomic windows (D17), with the family genes removed from
  their own windows; a genome-wide block scan ranking blocks by shared gene
  *families*; a quartet test; cross-species replication.
- *Caveat.* Paralogon signal decays; a negative is weak evidence of absence.
  Report the block scan's rank, not just a p-value.

**A4. The 3R test.** 🟡 (needs S5, S8)
- *Question.* Are zebrafish `itpr1a` / `itpr1b` (confirmed to exist, 2,696–
  2,819 aa) a teleost-3R ohnolog pair, and did the other two paralogs keep
  duplicates anywhere?
- *Method.* Pre-3R outgroup / teleost / extra-WGD grouping with
  Acipenseriformes and Salmoniformes held out; double-conserved synteny
  against both a tetrapod and a pre-3R ray-finned reference; cross-anchor
  block identity across ≥ 6 orders.
- *Caveat.* Each genome's copy labels are arbitrary, so a matched-vs-crossed
  pairing statistic has no power. Ask whether independent anchors agree.

## S17 — Constraint and function: does the sequence explain the disease?

**A5. Constraint mapped onto the channel.** 🟡 (needs S6, S9, S11)
- *Question.* Which residues are most constrained, and is the constrained
  core the same in all three paralogs?
- *Method.* Four conservation layers on the same residues (deep
  within-paralog / vertebrate / family-wide / representative-alignment as the
  control), Henikoff-weighted Jensen-Shannon against a BLOSUM background;
  transfer the domain architecture into each reference's own numbering by
  pairwise alignment; test each element against the *same protein's* distal
  regions rather than against a global mean.

**A6. The clinical-variant test.** 🟡 (needs A5)
- *Question.* Do the pathogenic alleles — SCA15/SCA29 and Gillespie syndrome
  in ITPR1, recessive anhidrosis in ITPR2, dominant neuropathy in ITPR3 —
  sit in the most constrained positions, and are they enriched in particular
  structural elements (the pore module, the IP3-binding core, the C-terminal
  coupling tail)?
- *Method.* ClinVar plus UniProt features, with **transcripts mapped, not
  assumed**: fetch each cited transcript's own translated CDS, align to the
  canonical, require the reference amino acid to match. ROC-AUC of
  constraint separating pathogenic from benign, with a protein-wide control
  and a Fisher test for element enrichment.
- *Deliverable.* The figure that makes this project useful to a clinical
  geneticist: constraint, structure and variant class on one axis system.
- *Caveat.* ClinVar's benign set is ascertainment-biased toward
  well-sequenced genes. State it, and keep the protein-wide control.

**A7. Three phenotypes, three paralogs.** 🟡
- *Question.* The three paralogs cause three unrelated diseases. Is that
  because they differ in sequence where it matters, or only in where they
  are expressed? Compare cross-paralog identity per structural element: if
  the machine is identical and only the regulatory regions differ, the
  answer is expression, and that is a clean, quotable result.

## S18 — Annotation quality (a result in its own right)

**A8. How the family is recorded.** 🟡 (needs S5, S15)
- *Question.* How often is a real ITPR locus missing, fragmentary, split
  across models, unnamed — or filed as a ryanodine receptor?
- *Why this family.* ~2,700 residues, ~58 exons, and a sister family sharing
  every diagnostic domain. If any family is systematically mis-recorded,
  it is this one.
- *Method.* Score every gene-scale locus (complete / split / fragmentary /
  noncoding / unannotated, same strand only); name and symbol verdicts
  applied identically to genome gene names and UniProt gene fields; join
  S15's integrity call as the **veto** on any correction (D6); break results
  down by annotation source (D9); protein-side audit with a recorded
  bait margin (D7).
- *Deliverable.* A correction list with assembly, coordinates, current
  state, proposal, evidence file and priority.

## S19 — Methods: what the field should learn from how we searched

**A9. Per-method contribution.** 🟡 — what would proteome-only searching have
missed? Restrict the comparison to accessions actually present in each
searched database, or the answer conflates "the HMM missed it" with "it was
never there". Ask it **per gene**, not per record.

**A10. Assembly contiguity as a confounder.** 🟡 — use whichever paralog turns
out never to be lost as the control whose every non-`found` cell is a false
negative; scan the contiguity floor that D4 asserts a priori.

**A11. Bait-panel ablation.** 🟡 — re-cluster retained alignments from bait
subsets to reproduce the sweep under a different panel. **Validate the
reproduction cell-for-cell against the committed ledger before trusting any
ablation.**

## S20 / S23 — Range: the family beyond animals

**A12. The plant and fungal records.** 🟡 (needs S2, S20, S23)
- *Question.* The database snapshot holds 40 Viridiplantae and 41 Fungi
  PF08709 records in a family textbooks say plants and fungi lack, and zero
  in *Arabidopsis* and *S. cerevisiae*. What are those records — real genes
  in early-diverging lineages, mis-annotations, or contamination? And is the
  absence in land plants and dikarya a fact about **genomes** or about
  **proteome databases**?
- *Method.* Individual verdicts per record (S20), then genome-level searches
  in the specific lineages (S23) at a stated sensitivity.
- *Deliverable.* The project's headline range result either way.
- *Caveat.* A negative needs its sensitivity stated, and a contaminant needs
  positive evidence (assembly provenance, neighbouring genes), not just
  surprise.

**A13. Where the IP3R/RyR split falls.** 🟡 — with both families enumerated by
the same searches, the split can be placed on the eukaryote tree rather than
assumed. Ask which lineages have one, both, or neither.

## S21 — Gene architecture

**A14. The exon map.** 🟡 (needs S5) — a ~58-exon architecture across the
genome scope, from CDS blocks the sweep already computes. Are conserved
intron positions shared between paralogs at the same alignment column *and*
phase (test against a permutation null drawn from columns where both have
residues)? Are database "fragments" bounded by real splice sites?

## S22 — Ligand-site evolution

**A15. The IP3-binding core.** 🟡 (needs S9, S17) — the one part RyR does not
share functionally. Is it under different constraint from the pore? Do the
IP3-contacting residues stand out? And in lineages with a reduced upstream
PLC/IP3 pathway, is the core relaxed relative to the pore? Report the taxon
list and the power honestly; this test is easy to underpower.

---

## Report skeleton (S14) — where each analysis lands

| § | Section | Fed by |
|---|---------|--------|
| R1 | Enumerated census & completeness argument | S2, S3 |
| R2 | The family across the eukaryotes; the plant/fungal question | S20, S23 |
| R3 | Genome-wide ledger: present / lost / assembly-gap | S5 |
| R4 | Novel and unannotated ITPRs | S5, census v4 |
| R5 | Phylogeny, the RyR root, and which paralogs are sisters | S7 |
| R6 | Duplication history: 2R, 3R, and the RyR parallel | S16, S13 |
| R7 | Paralog fates: retention, loss, decay | S15 |
| R8 | Selection and constraint | S9, S17 |
| R9 | Structure-mapped conservation and the clinical variants | S11, S17 |
| R10 | Expression evidence | S12 |
| R11 | Annotation-quality audit and corrections | S18, S10, S21 |
| R12 | Methods: contribution, confounders, panel design | S19 |
| D | Discussion: range, ligand-site evolution, variant interpretation | S20, S22 |
