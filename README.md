# Protein Variant Finder — and the ITPR census

Two things live in this repository:

1. **A GUI/CLI tool** that enumerates every known protein isoform / ortholog
   of a gene across NCBI, Ensembl, UniProt, Compara, BLAST, AlphaFold DB and
   Foldseek, then analyses what it finds (MSA → tree → clusters →
   novel-paralog scoring).
2. **A publication project** using that tool: a systematic census of the
   **IP3 receptor (ITPR) family** — the endoplasmic reticulum's ligand-gated
   calcium-release channel — across a declared set of genomes and proteomes,
   to the evidence standard of an MBE / GBE / Genome Research paper. Plan and
   task ledger: [`PUBLICATION_ROADMAP.md`](PUBLICATION_ROADMAP.md);
   per-session narrative: [`SESSION_LOG.md`](SESSION_LOG.md); the biology:
   [`docs/ip3r_background.md`](docs/ip3r_background.md).

> **Reading order for a new session:** `PUBLICATION_ROADMAP.md` (protocol +
> ledger) → `docs/session_briefs.md` (the current task's brief) →
> `INTERFACE.md` (module map). This README is the state-of-the-project
> summary and is refreshed at the end of every session.

**Status: set up, nothing measured.** `results/` is empty by design. The next
session runs **S0**.

---

## Quickstart

```bash
# GUI
python run.py --email you@example.com

# Headless, against live APIs
python run.py --headless --preset ip3r --save-results
python run.py --headless --preset ip3r_zebrafish --save-results

# The project dashboard (stdlib only, no env needed)
python3 scripts/dashboard.py --open
python3 scripts/dashboard.py --watch &

# Is the bulk-storage drive attached? (exit 1 = no; that is a stop)
python -m src.utils.data_root --require
```

Dependencies: `pip install -r requirements.txt` (`biopython`, `requests`).
The analysis toolchain (MAFFT, HMMER, BLAST+, miniprot, trimAl, IQ-TREE) is
installed and version-pinned in **S1**.

---

## The science

The IP3 receptor releases calcium from the ER when IP3 binds it — the origin
of most agonist-evoked calcium signals in non-muscle cells, of the calcium
waves that pattern development, and of the ER-to-mitochondrion transfer that
sets the apoptotic threshold. Vertebrates have three paralogs of about 2,700
residues each, assembled as tetramers, and each is associated with a
different human disease:

| Gene | UniProt | Length | Disease association *(verified in S0)* |
|------|---------|--------|-----------------------------------------|
| ITPR1 | Q14643 | 2,758 aa | spinocerebellar ataxia (SCA15/SCA29), Gillespie syndrome |
| ITPR2 | Q14571 | 2,701 aa | autosomal recessive isolated anhidrosis |
| ITPR3 | Q14573 | 2,671 aa | dominant Charcot–Marie–Tooth neuropathy |

**The complication that shapes the project.** The ryanodine receptors
(RYR1/2/3, ~5,000 aa) carry *every* Pfam domain that marks an IP3 receptor —
checked directly: human RYR1 carries PF08709, PF01365, PF08454 and PF02815.
They are inside every search this project runs. So separating ITPR from RYR
is a positive test at every stage (roadmap Decisions **D14**) — and, turned
around, it is an opportunity: both families get counted by the same
instrument, and RyR gives the ITPR tree a proper outgroup.

### What is already visible, before any analysis

Checked against InterPro and UniProt on 2026-08-18, recorded in
[`docs/ip3r_background.md`](docs/ip3r_background.md), and to be re-derived
properly in S2 before any of it is quoted:

- The family looks overwhelmingly like an animal family — **12,149 of the
  12,338** proteins carrying PF08709 are metazoan.
- But there are **40 plant and 41 fungal records** in a family textbooks say
  plants and fungi lack, while *Arabidopsis* and *S. cerevisiae* have
  **none**. Real branch of the family, mis-annotation, or contamination?
- Zebrafish carries **four**: `itpr1a`, `itpr1b`, `itpr2`, `itpr3` — the
  first two look like a teleost-duplication pair.

### The questions the project answers

- **Range.** What is the family's true taxonomic distribution, and are the
  famous absences facts about genomes or about proteome databases?
- **Origin.** Are ITPR1/2/3 2R ohnologs? Are the teleost pairs from 3R? Did
  the ryanodine receptors triplicate at the same time?
- **Fates.** Has any paralog ever been lost, and how often independently?
- **Records.** How often is a real ITPR locus missing, split, unnamed, or
  filed as a ryanodine receptor?
- **Machine.** Do the pathogenic variants sit in the most constrained parts
  of the channel, and is the constrained core the same in all three paralogs?
- **Anything unnamed.** Is there a fourth vertebrate ITPR, or an ITPR in a
  lineage reported to lack one?

---

## Project status board

One task per session. Full ledger with dependencies and results in
[`PUBLICATION_ROADMAP.md`](PUBLICATION_ROADMAP.md).

| ID | Task | Status |
|----|------|--------|
| S0 | Literature baseline + scope confirmation | ⏳ pending — **next** |
| S1 | Toolchain + positive/negative controls (RyR is the sharp decoy) | ⏳ pending |
| S2 | Uncapped InterPro enumeration → census v2 | ⏳ pending |
| S3 | Profile-HMM sweep (itpr.hmm + ryr.hmm) + convergence argument | ⏳ pending |
| S4 | Genome scope manifest (the denominator) | ⏳ pending |
| S5 | Genomic sweep → per-genome ledger + novel gene models | ⏳ pending |
| S20 | Non-vertebrate sweep — the family's true range | ⏳ pending |
| S23 | Invertebrate / protist / plant / fungal genome sweep | ⏳ pending |
| S6 | Alignment upgrade (MAFFT L-INS-i + trimAl) | ⏳ pending |
| S7 | ML phylogeny, rooted on RyR; which paralogs are sisters | ⏳ pending |
| S8 | Synteny across the ITPR loci | ⏳ pending |
| S9 | ML selection (codeml, HyPhy RELAX) | ⏳ pending |
| S10 | Annotation-bug molecular validation | ⏳ pending |
| S11 | Structures + TM-align vs cryo-EM references | ⏳ pending |
| S12 | Expression evidence (SRA junction reads) | ⏳ pending |
| S13 | Reconciliation & dating | ⏳ pending |
| S15 | Loss dynamics | ⏳ pending |
| S16 | Duplication history (2R / 3R, and the RyR parallel) | ⏳ pending |
| S17 | Constraint & function — the clinical-variant test | ⏳ pending |
| S18 | Annotation-quality audit + correction list | ⏳ pending |
| S19 | Methods results | ⏳ pending |
| S21 | Gene architecture (~58 exons) | ⏳ pending |
| S22 | Ligand-site evolution | ⏳ pending |
| S14a | Manuscript assembly | ⏳ pending |
| S24 | Supplementary figures + figure audit | ⏳ pending |
| S14c | Manuscript rewrite pass | ⏳ pending |
| S14b | Deposit + release (Zenodo, public repo, preprint) | ⏳ pending — human-gated |

---

## Findings so far

Nothing yet — see [`FINDINGS.md`](FINDINGS.md), which gets a plain-language
entry after every completed task. Each session also adds its headline result
here with the figure that shows it.

---

## Storage

Bulk data (genome assemblies, reference proteomes, BLAST databases,
structures, SRA) never enters the repository. It lives under the path in
[`data_root.txt`](data_root.txt) — currently `/Volumes/FANTOM/IP3R_DATA`, an
**external drive that must be attached before any bulk session**. Moving it
is a one-line edit; the directory layout recreates itself.

`python -m src.utils.data_root --require` exits non-zero when the drive is
absent, and the session protocol runs it at the start of every session, so a
400 GB download cannot land on the laptop by accident.

Expected footprint, from the equivalent PIEZO-project measurements:
assemblies ~400 GB (fetch → search → delete keeps the live footprint near
20 GB), reference proteomes ~10 GB, sweep evidence ~30 GB, HMMER output
~12 GB, structures ~2 GB.

---

## Provenance

The application, the figure style, the dashboard, the manuscript-assembly
and claim-checking tooling and the session protocol are ported from the
PIEZO project (`../piezo_genes`), which took 24 sessions to build and
harden them. Its `INTERFACE.md` maps a reference implementation for most
tasks here, and its methodological rules are carried forward as **D0–D17**
in the roadmap's Decisions log rather than being re-derived.

The family is defined in exactly one place — `src/utils/family.py` — so
nothing else in `src/` hard-codes a gene name. **No result is ported.**
Every number in this repository must come from this project's own runs.
