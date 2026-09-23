# The IP₃ receptor family across 503 genomes

**A census of the IP₃ receptor (*ITPR*) gene family across the eukaryotes. The
search is declared and controlled throughout, and every number in it can be
traced back to a committed table.** The repository also holds the
database-search app the census was built on.

> *Retained in every vertebrate, lost repeatedly elsewhere: a 503-genome
> census of the IP₃ receptor family* — manuscript in
> [`manuscript/`](manuscript/README.md), long-form thesis in
> [`thesis/`](thesis/README.md), literature review in
> [`docs/ip3r_review_2026.pdf`](docs/ip3r_review_2026.pdf).

**Status:** analysis complete (34 of 37 ledger tasks done). The paper (60 pp) and
the thesis (194 pp) both build end to end from the committed tables. Still to
do: a revision pass on figure text and typesetting (S28, in progress), a
paper series (S26), and deposit and release (S14b, which needs a human). See
[Project status](#project-status).

---

## Contents

- [Background](#background)
- [What was done](#what-was-done)
- [Key results](#key-results)
- [What the repository produces](#what-the-repository-produces)
- [Getting started](#getting-started)
- [Reproducing the analyses](#reproducing-the-analyses)
- [How the evidence is held to account](#how-the-evidence-is-held-to-account)
- [Repository layout](#repository-layout)
- [Project status](#project-status)
- [How this project was made](#how-this-project-was-made)
- [Citing, licence, contact](#citing-licence-contact)

---

## Background

When a cell-surface receptor activates phospholipase C, the cell makes
inositol 1,4,5-trisphosphate (IP₃). IP₃ opens the **IP₃ receptor**, a
~2,700-residue tetrameric channel in the endoplasmic reticulum membrane, and
calcium floods out of the ER store. Most agonist-evoked calcium signals in
non-muscle cells start at this channel, as do the calcium waves that pattern
development and the ER-to-mitochondrion calcium transfer that sets the
apoptotic threshold.

![The IP₃ receptor: the pathway, the channel and its sequence](docs/figures/receptor_overview.png)

*(a) The pathway that makes the ligand (schematic). (b) The channel's
dimensions measured on cryo-EM structure 6DQN, with the domain architecture
of one subunit. (c) The sequence at an IP₃ contact and at the selectivity
filter, in the three human paralogues and the fly receptor.*

Vertebrates have three paralogues, and each is linked to a different disease:

| Gene | UniProt | Length | Disease association |
|------|---------|-------:|---------------------|
| *ITPR1* | [Q14643](https://www.uniprot.org/uniprotkb/Q14643) | 2,758 aa | spinocerebellar ataxia 15/29, Gillespie syndrome |
| *ITPR2* | [Q14571](https://www.uniprot.org/uniprotkb/Q14571) | 2,701 aa | isolated anhidrosis |
| *ITPR3* | [Q14573](https://www.uniprot.org/uniprotkb/Q14573) | 2,671 aa | Charcot–Marie–Tooth neuropathy |

Almost everything known about where the family occurs comes from searches
whose scope was never declared. An absence is only a claim about a search
space, so it can only be tested once that space is stated. Three properties
make this family a hard case:

1. **Its sister family carries all of its diagnostic domains.** The
   ryanodine receptors (RYR1/2/3, ~5,000 aa) carry every Pfam signature used
   to recognise an IP₃ receptor. Any search sensitive enough to find a
   divergent IP₃ receptor also finds them: 6,807 of the first 15,421 records
   enumerated were ryanodine receptors. Every stage of this project therefore
   makes the ITPR/RyR call as a positive test, and the ryanodine receptors
   are searched alongside as an internal control.
2. **The gene is long.** It has ~58 exons over tens to hundreds of
   kilobases, which is longer than the contigs of many published assemblies.
   A "missing" gene is often just a broken assembly.
3. **It is poorly named.** Much of the record sits under placeholder
   symbols or no symbol at all, so searching by name undercounts the family.

The biology baseline, with every statement tagged by how far it can be
trusted, is [`docs/ip3r_background.md`](docs/ip3r_background.md).

---

## What was done

The census is **genomic**, not a keyword search. The scope is declared up
front, and every absence must come with a positive control showing that the
search actually reached that assembly.

| Search space | Size |
|---|---|
| Reference proteomes swept with family and sister-family profile HMMs | **7,691** (77.6 M proteins): 763 vertebrate, 6,928 other eukaryotes, plus 634 archaea and 3,537 bacteria |
| Vertebrate genomes searched directly (one per order, plus every species the protein record left in doubt) | **309** |
| Non-vertebrate genomes, sampled most finely where the absence claims are | **194** |
| ITPR genes demonstrated in genomes | 1,232, with the RyR control measured alongside |

On top of the census the project builds a 134-tip ML phylogeny rooted on the
ryanodine receptors, synteny and paralogon tests, codon-model selection
analyses (codeml, HyPhy RELAX/FEL), gene-tree/species-tree reconciliation
with dating, loss counting under Dollo and Mk models, exon–intron
architecture, residue-level constraint mapped onto cryo-EM structures,
a ClinVar variant test, an annotation audit of every locus, and a
measurement of the method's own false-negative rate.

---

## Key results

### 1. The family is ancestrally eukaryotic and has been lost many times, independently

The family is found in animals, amoebae, green algae, early-branching fungi
and several protist lineages. It is absent from **all 384 land-plant** and
**all 1,353 Dikarya** reference proteomes, and from Apicomplexa,
Microsporidia, red algae, diatoms and four further fungal phyla. It is also
absent from every archaeal and bacterial proteome. All 35 clade-level
absences hold at the genome level, in assemblies where a measured positive
control recovers a comparably long, deeply conserved gene.

![The ITPR family across the eukaryotes](results/s20_sweep/figures/range_by_phylum.png)

### 2. No vertebrate has ever lost one of the three paralogues

None of the 927 genome × paralogue cells across 309 vertebrate genomes
reaches the absence state. In the 189 assemblies contiguous enough to hold
the gene on one contig, every genome carries at least three gene-equivalents.
Dollo parsimony places no loss anywhere on the tree, so the project reports
instead which settings would *manufacture* one.

![Every genome × paralogue cell, ordered by assembly contiguity](results/loss_dynamics/figures/s15_character_matrix.png)

### 3. *ITPR2* and *ITPR3* are sisters, and the duplications are dated

The approximately unbiased (AU) test favours *ITPR2*+*ITPR3* (p = 0.476) and
rejects both alternatives at p < 2 × 10⁻⁵. The *ITPR1* / (*ITPR2*, *ITPR3*)
split sits on the vertebrate stem. The *ITPR2*/*ITPR3* split sits on the
gnathostome stem, at 462–563 Ma. The surviving two-round (2R) paralogon links
run through *ITPR1*.

![ML phylogeny rooted on the ryanodine receptors](results/phylogeny/figures/tree_ml_rooted.png)

![The 2R paralogon around the three receptors](results/duplication/figures/s16_paralogon.png)

### 4. *ITPR1* is the special one

*ITPR1* is the only paralogue that was doubled and kept after the teleost
genome duplication (97.3 % of 73 teleost genomes carry both copies). It is
also held under roughly twice the purifying selection of its sisters:
ω = 0.024 against 0.043 and 0.042, and RELAX finds intensified selection
(k = 9.36).

![Purifying selection on every paralogue](results/selection/figures/s9_omega_by_paralog.png)

### 5. One gene architecture, inherited intact

In every genome that carries them, the three paralogues share 46–49 of their
~58 intron positions, an 85–88-fold enrichment over chance. The ryanodine
receptors carry every diagnostic domain of the family but share just one.

![Shared intron positions](results/gene_architecture/figures/intron_positions.png)

### 6. What the channel cannot change

The gate and the selectivity filter are the least changeable elements of the
receptor, and the gate is identical across paralogues. The ten measured IP₃
contacts lie *outside* the Pfam domain named after the ligand. At the ligand
site, the constrained unit is a 15 Å pocket rather than the contacts
themselves. The IP₃-binding core comes out *less* constrained than the pore,
but that comparison reverses depending on whether a 50-residue luminal loop,
the least conserved sequence in the receptor, is counted as part of the pore.

![Constraint along the pore-forming half](results/constraint/figures/s17_channel_profile.png)

### 7. The public record is worse than the family, and the archive decides it

55.0 % of full-length ITPR protein records carry no usable gene symbol. Of
the 1,232 genes demonstrated in genomes, 940 (76.3 %) cannot be reached by
any protein-database search. Whether a gene is recorded correctly depends on
the archive that holds it, not on the gene. RefSeq annotations are 98.8 %
complete, submitter GenBank annotations 37.5 %. Above the contiguity bar the
figures are 99.4 % and 63.8 %. The ryanodine receptors, measured the same way
in the same assemblies, are annotated no better, so this is a property of the
archive rather than of this family. A correction list addressed to the
databases is in
[`results/annotation_audit/corrections.tsv`](results/annotation_audit/corrections.tsv).

![Annotation state by archive](results/annotation_audit/figures/s18_fig1_by_source.png)

### 8. The method's own error rate is measured

No IP₃ receptor gene has been lost anywhere in the vertebrate scope, so
every cell where the search came up empty is a **false negative of the
method**: 140 of 923, or 15.2 %. The ryanodine receptors, measured
independently, agree at 13.6 %. Every such failure traces to the assembly,
and above the contiguity bar the miss rate falls to 0.9 %. Two findings bear
on how a survey like this should be run. Bait-panel breadth is almost free:
four human baits recover 782 of the 783 genes that the 38-bait panel finds.
And the standard guard against drift in iterative searches fired on none of
the seven runs, including the three that did drift.

The plain-language story, task by task, is in [`FINDINGS.md`](FINDINGS.md).
Each analysis also has its own rendered `results/<task>/report.md`.

---

## What the repository produces

| Output | Where | Built by |
|---|---|---|
| **Manuscript**: 17 sections, 7 main + 16 Extended Data + 6 Supplementary figures, 276 load-bearing numbers re-verified on every build, 60-page PDF | [`manuscript/`](manuscript/README.md) · [PDF](manuscript/itpr_family_manuscript.pdf) | `scripts/s14_assemble.py` |
| **Thesis**: 16 chapters + 5 appendices, ~73,000 words, 104 figures, 116 audited references, 194-page PDF | [`thesis/`](thesis/README.md) · [PDF](thesis/itpr_family_thesis.pdf) | `scripts/s25_assemble.py` |
| **Literature review**: 32 pages, 137 references, 12 figures, with a claim-by-claim audit | [`docs/ip3r_review_2026.pdf`](docs/ip3r_review_2026.pdf) | `scripts/s0_review_build.py --pdf` |
| **Analysis results**: one directory per task, holding tables, figures and a rendered `report.md` | [`results/`](results/) | `scripts/s<n>_*.py` |
| **Deposit manifest**: 1,951 files with SHA-256, plus the commands that regenerate the excluded bulk data | [`manuscript/deposit_manifest.tsv`](manuscript/deposit_manifest.tsv) | `s14_assemble.py --only deposit` |
| **Protein Variant Finder**: GUI/CLI app | [`src/`](src/), [`run.py`](run.py) | `python run.py` |

### The Protein Variant Finder app

The app the project started from lists every protein isoform and ortholog of
a gene across **NCBI, Ensembl, UniProt, Ensembl Compara, BLAST, InterPro,
AlphaFold DB and Foldseek**, normalises the results into one record type, and
then analyses them: MSA → NJ tree → clusters → novel-paralogue scoring, with
dN/dS, presence/absence, ESMFold fold checks and per-accession
"investigation" case files. Network calls never run on the GUI thread. The
family is defined in a single file,
[`src/utils/family.py`](src/utils/family.py), so pointing the app at another
gene family means editing that file and the presets.

---

## Getting started

```bash
git clone git@github.com:gddickinson/ip3r_genes.git
cd ip3r_genes
pip install -r requirements.txt          # biopython, requests

# The GUI
python run.py --email you@example.com    # NCBI asks for a contact address

# Headless, against the live APIs
python run.py --headless --preset ip3r --save-results
python run.py --headless --preset ip3r_zebrafish --save-results
```

Bundled presets in [`presets/`](presets/): `ip3r`, `ip3r_zebrafish`,
`ip3r_discovery`, `ip3r_paralog_mine`, `ip3r_domain_scan` (enumerates every
protein carrying a family Pfam) and `ip3r_all` (exhaustive harvest plus every
analysis).

### Analysis toolchain

The publication pipeline shells out to pinned versions of: MAFFT 7.526,
HMMER 3.4, BLAST+ 2.16.0, miniprot 0.18, trimAl 1.5, IQ-TREE 2.3.6, PAML,
HyPhy, TM-align, HISAT2, sra-tools, NCBI `datasets` 18.35 and Foldseek.
Exact versions and paths are in
[`results/toolchain_manifest.txt`](results/toolchain_manifest.txt), which
`python scripts/s1_toolchain.py` regenerates. Documents are typeset with
pandoc + xelatex.

### Bulk data

Genomes, proteomes, BLAST databases and structures (far too large for git) never
enter the repository. They live under a data root named in
[`data_root.txt`](data_root.txt):

```bash
python -m src.utils.data_root --require   # exits non-zero if the data root is missing
```

---

## Reproducing the analyses

Each ledger task has a driver with ordered, resumable stages (`--list`,
`--only <stage>`, `--from <stage>`). Each driver runs its negative controls
before it writes anything.

```bash
python scripts/s7_run.py --list           # phylogeny
python scripts/s17_run.py                 # constraint & clinical variants
python scripts/s14_assemble.py            # figures → claims → stitch → pdf → deposit
python scripts/s25_assemble.py            # the thesis
python3 scripts/dashboard.py --open       # project dashboard (stdlib only)
```

Most analyses from S6 onwards read only committed tables and rebuild
offline. The sweeps (S2–S5, S20, S23) need the network and the bulk data
root. [`INTERFACE.md`](INTERFACE.md) maps every script and module, and
[`docs/session_briefs.md`](docs/session_briefs.md) gives each task's goal,
steps and completion criteria.

---

## How the evidence is held to account

The project's rules are recorded as numbered decisions in the
[roadmap's Decisions log](PUBLICATION_ROADMAP.md). The ones that shape every
result:

- **Reports render from tables, never from memory (D13).** Every
  `report.md`, every figure and the manuscript's numbers are generated from
  committed TSV/JSON. The manuscript's claims ledger re-checks 276 numbers
  against their source tables on every build and fails on any drift.
- **The ITPR/RyR call is a positive test (D14).** It uses best-profile
  assignment or a labelled-bait margin. Length is supporting evidence, never
  the call.
- **An absence needs a control.** A genome whose positive control fails
  supports no absence claim, and an assembly too fragmented to hold the gene
  is `undecidable`, not `absent` (D4).
- **Thresholds are measured, not chosen.** Intron limits, identity floors,
  reassembly bars and completeness bars are each calibrated against evidence
  the instrument did not produce, and the calibration is committed next to
  the result.
- **Negative controls run before anything is written.** 437 constructed
  controls across the test modules check that each rule *refuses* what it
  should. Most suites are also mutation-tested by deliberately breaking a
  rule and confirming the suite catches it.
- **Reproducible by construction (D24).** Aligners and tree searches run at
  a pinned thread count and seed, databases are pinned to dated releases,
  and figures are byte-identical on rebuild.
- **References are audited on entry.** A new citation is declared by PMID or
  DOI plus a phrase its title must contain, and is resolved against the live
  record. This caught 9 references that had been written from memory and
  pointed to entirely different papers.

---

## Repository layout

```
run.py                     app entry point (GUI / --headless)
src/                       the app: core/, databases/, analysis/, discovery/,
                           investigation/, gui/, utils/family.py
scripts/                   one module family per ledger task (s0_… s28_),
                           plus figstyle.py, dashboard.py and the build drivers
results/<task>/            committed tables, figures and report.md per task
manuscript/                the paper: hand-written sections + generated build
manuscript_v1/             the first draft, frozen
thesis/                    the long form: chapter sources + generated build
docs/                      background, literature review, session briefs,
                           analysis catalogue, review figures
presets/                   bundled app queries
PUBLICATION_ROADMAP.md     task ledger, Decisions log, session protocol
FINDINGS.md                the biological story, task by task
SESSION_LOG.md             what ran in each session
INTERFACE.md               module-by-module navigation map
```

---

## Project status

The work ran as a ledger of single-task sessions. Dependencies, results and
decisions for every task are in
[`PUBLICATION_ROADMAP.md`](PUBLICATION_ROADMAP.md).

| Phase | Tasks | Status |
|---|---|---|
| Baseline & controls | S0 literature audit · S1 toolchain and control benchmark | ✅ |
| Enumeration | S2 InterPro census · S3 profile-HMM sweep · S4 genome manifest · S5 309-genome vertebrate sweep · S20 non-vertebrate proteome sweep · S23 194-genome non-vertebrate sweep | ✅ |
| Evolution | S6 alignment · S7 phylogeny · S8 synteny · S9 selection · S13 reconciliation & dating · S15 loss dynamics · S16 duplication history · S21 gene architecture | ✅ |
| Function & records | S10 annotation-bug validation · S11 structures · S12 expression · S17 constraint & variants · S18 annotation audit · S19 methods · S22 ligand site | ✅ |
| Writing | S14a/S14c manuscript · S24 supplementary figures & figure audit · S25 thesis · S27 editorial pass | ✅ |
| Revision | S28 legends, in-figure text, typesetting | 🔄 in progress |
| Series | S26 the results regrouped as individual papers | ⏳ pending |
| Release | S14b Zenodo DOI, public repository, preprint | ⏳ pending, needs a human |

What still needs a human before submission is listed in
[`manuscript/reviewer_checklist.md`](manuscript/reviewer_checklist.md).

---

## How this project was made

The project was carried out with [Claude Code](https://claude.com/claude-code),
following the protocol in the roadmap: one ledger task per session, a written
brief with completion criteria, a literature baseline audited claim by claim,
a download scope decided by rule and committed as a manifest, analyses left
to run unattended, and every report, figure and document generated from
committed tables. The app, the tooling and the protocol were ported from a
companion project on the PIEZO family; methods were ported, results were
not. Chapter 16 of the thesis describes the process, and its numbers are
measured from this repository on every build.

---

## Citing, licence, contact

The manuscript is in preparation. A Zenodo DOI and preprint will be added
here at release (S14b). Until then, please cite the repository:

> Dickinson, G. (2026). *Retained in every vertebrate, lost repeatedly
> elsewhere: a 503-genome census of the IP₃ receptor family.*
> https://github.com/gddickinson/ip3r_genes

**Licence:** not yet chosen. It will be set at release, and until then all
rights are reserved.

**Contact:** George Dickinson — george.dickinson@gmail.com
