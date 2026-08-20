# PUBLICATION_ROADMAP.md — publication-grade ITPR census & discovery

**Goal.** Enumerate and validate every IP3-receptor (ITPR) gene across a
declared genome scope, and answer what the family's records cannot currently
answer: its true taxonomic range, where the three vertebrate paralogs came
from, what happened to each of them, whether its constrained core explains
its clinical variants, and how badly it is recorded. To the evidence
standard of an MBE / GBE / Genome Research paper: an exhaustive enumerated
search space with a completeness argument, an ML phylogeny with support
values, synteny-backed orthology, ML selection tests, molecular validation
of every annotation-error claim, structural confirmation, and one-command
reproducibility.

**How the claim is scoped.** "A systematic census across N vertebrate
genomes and M eukaryotic reference proteomes" (N fixed in S4, M in S3/S20),
never an unbounded "all IP3 receptors". Every absence claim is a claim about
a *declared* search space.

**The biology this rests on** is `docs/ip3r_background.md`, which tags every
statement `[db]` (verified against a live database), `[lit]` (literature,
verified in S0) or `[open]` (a question this project answers). Read it before
S1.

**The one thing that will bite you.** The ryanodine receptors are inside
every search this project runs — RYR1 carries all four of the ITPR-diagnostic
Pfam domains. Separating ITPR from RYR is a positive test at every stage, not
an assumption (see Decisions **D14**).

---

## Session protocol

Claude: follow this protocol in every session that touches this project.

### Start of session
0. Open the project dashboard: `python3 scripts/dashboard.py --open` and
   start its watcher in the background (`python3 scripts/dashboard.py
   --watch` as a background task) so the page stays live through the
   session. Stdlib-only — no conda env needed.
1. Read this file top to bottom. Read `docs/session_briefs.md` for the brief
   of the task you are about to work.
2. `git pull`.
3. Run `python -m src.utils.data_root --require`. It exits non-zero if the
   bulk-storage drive is not attached. **If it fails, stop and tell the
   user** — do not download anything, and do not "temporarily" use the
   internal disk. Note free space (`df -h`) before any bulk task.
4. Pick the task: the single ledger row whose Status is `in_progress`, else
   the topmost `pending` row whose dependencies are all `completed`.
   Announce it to the user. Set its Status to `in_progress` (commit at
   session end).
5. Work **only that task** until its completion criteria pass. Do not start
   the next task even if time remains — spend surplus on tests, docs, or
   hardening of the current task.

### During the session
- Big/raw files → `get_data_root()` subdirs, never the repo. Committed
  artefacts: summaries, TSV/CSV ≤ ~20 MB, figures, manifests, code.
- Record every load-bearing number and file path in the task's Results
  column (link out to a file under `results/` if long).
- Write `results/session_live.json` from any long-running driver
  (`{"task": "S5", "workers": 2, "steps": [{"label": ..., "done": ...}]}`)
  so the dashboard shows real progress.
- New leads, surprises, or scope changes → add a row to **Emergent tasks**.
  Never silently expand the current task.
- A task that proves too large for one session: split it in the ledger
  (S5 → S5a/S5b), complete the first part properly, leave the rest `pending`
  with a note. Splitting is normal; half-done is not.
- Anything needing an interactive login or sudo: ask the user to run it via
  the `!` prefix.

### End of session (checklist — do all of it)
1. Completion criteria all pass → Status `completed` + date. Not all passing
   → stays `in_progress` with a NEXT note giving the exact resume point.
2. Update this file (ledger, Results, Emergent, Decisions).
3. Append a dated entry to `SESSION_LOG.md` (what ran, what resulted, what's
   next).
3a. **Append a biological-findings entry to `FINDINGS.md`** — what the task
   changed in the biological story, in plain language, findings before
   methods, unconfirmed claims marked *(pending: which task confirms it)*.
   Keep methods and paths out of it.
3b. **Refresh `README.md`** — it is the state-of-the-project summary for
   anyone arriving cold. Update the status board row, and if the session
   produced a result worth showing, add it to **Findings so far** with a link
   and an embedded figure. Delete what the session made obsolete.
4. `git add -A && git commit` (message: `S<n>: <one-line outcome>`) and
   `git push` if a remote is configured.
5. Tell the user: task status, headline results, what the next session does.

---

## Storage

- Active bulk-data root: see `data_root.txt`, read via
  `src/utils/data_root.py:get_data_root()`. Currently
  `/Volumes/FANTOM/IP3R_DATA` — **an external drive that must be attached
  before any bulk session**. Moving it later is a one-line edit; the
  directory layout recreates itself.
- `require_data_root()` raises rather than falling back. Anything that
  downloads or sweeps must call it, so an unplugged drive stops the session
  at the top instead of after 200 GB has landed on the laptop.
- Expected bulk footprint, from the PIEZO project's measured equivalents:
  vertebrate assemblies ~400 GB (fetch → search → delete keeps the live
  footprint near 20 GB), reference proteomes ~10 GB, per-genome sweep
  evidence ~30 GB, HMMER raw output ~12 GB, structures ~2 GB.
- Budget rule: any single download > 20 GB needs a note in Decisions with
  the running disk total.

---

## Task ledger

Statuses: `pending` / `in_progress` / `completed YYYY-MM-DD`.
Full step-by-step briefs: `docs/session_briefs.md`.

> **Next session: S2.** S0 and S1 are complete. Read
> `results/benchmark_controls/report.md` before S2 — it changes how the
> census must be run.
>
> **Three S1 results S2 depends on.** (1) The scorer promoted every RyR
> decoy at 45 before S1 added the sister-family test; the test is now the
> thing standing between the census and 49 % RyR contamination, so **every
> S2 analysis set must carry the labelled RyR bait panel**
> (`src/utils/family.py:SISTER_PANEL`) or the test has nothing to measure
> against. (2) **Recall rests entirely on the domain component** — ITPR1/2/3
> are 61–68 % identical, so the twilight-zone component never fires for a
> vertebrate paralog; a census record with no attachable Pfam hit sits at 30
> and is missed. (3) The bait margin **cannot call the non-metazoan grade**
> (*Dictyostelium* iplA at +0.065, inside the D7 no-call band) — S2 records
> such rows as `unassigned` and leaves the call to S3's profiles.
>
> **The external drive was not attached during S0** (S0 needs no bulk
> storage). S4/S5 cannot start until it is.
>
> **This project is the PIEZO project's machinery pointed at a different
> family.** The app (`src/`), the figure style, the dashboard, the
> manuscript-assembly and claim-checking tooling and the session protocol are
> ported from `../piezo_genes`, which took 24 sessions to build them. The
> **Decisions** section below carries that project's methodological scar
> tissue forward as pre-agreed rules: they were paid for in sessions and
> should not be re-derived. What is *not* ported is any result — every
> number in this repo must come from this project's own runs.
>
> **Where the analysis scripts come from.** Each task writes its own
> `scripts/s<n>_*.py`. Where the PIEZO project solved the same problem, its
> script is a reference implementation worth reading before writing a new
> one (`../piezo_genes/scripts/`, mapped in `../piezo_genes/INTERFACE.md`) —
> but the family differs in ways that matter (a sister family inside every
> search, a tetramer not a trimer, a ligand-binding domain with a clinical
> variant set), so port the method, not the constants.

| ID | Task (one session each) | Depends | Status | Results (headline) |
|----|-------------------------|---------|--------|--------------------|
| S0 | Literature baseline + scope confirmation: verify every `[lit]` claim in `docs/ip3r_background.md`, build `docs/ip3r_review_2026.md`, smoke-test the app against live APIs | — | completed 2026-08-18 | **19 `[lit]` claims audited vs 51 refs (45 primary): 12 verified, 3 qualified, 2 struck, 1 → `[open]`, 1 → `[db]`.** `docs/ip3r_review_2026.md` written. All `[db]` numbers re-derived exactly. **D14 measured: 49 % (53/109) of zebrafish PF08709 records are RyRs.** Exon/span claim false — ITPR3 spans 76 kb, not "hundreds of kb"; span varies 6.5× across paralogs, length 3 %. Ensembl `/xrefs/symbol/homo_sapiens/` stalls (species-specific, `BRCA2` too); client switched to `/lookup/symbol/` → human ITPR1 **0 → 24 variants**. Smoke 3/3 on human (UniProt 34 / NCBI 44 / Ensembl 41) and on zebrafish; the 8-species panel is 2/3 on Ensembl latency alone. → `results/s0_baseline/report.md` **Extended 2026-08-19 (user request): the review is now illustrated** — 12 generated figures, 24 → 32 pages, rendered by `scripts/s0_review_figures.py` from committed tables only (D13, D19). Measuring the structure for the figures produced three results the text did not have: the pore's two constrictions recovered blind and landing on the GGGVGD filter motif and on Phe2513/Ile2517; the IP₃→gate distance resolved into 103 Å axial vs 120 Å through space; and ***Dictyostelium* iplA carries none of PF08709, PF02815 or PF00520** — a characterised receptor the family's defining signature does not find. |
| S1 | Toolchain install + positive/negative control benchmark (RyR is the sharp decoy) | S0 | completed 2026-08-18 | **Recall 24/25 (96 %); specificity 31/31 (100 %) — but only after S1 had to implement D14.** All 12 binaries resolve with recorded versions (`results/toolchain_manifest.txt`); MAFFT proven invoked (rc=0, 56 seqs → 9,494 cols). **On the first run all six RyR decoys were promoted at 45** — they carry all four diagnostic Pfams, so `pfam+20` also satisfies the D3 gate; specificity was 25/31 and every failure was a RyR. Added the **labelled-bait sister-family test** to `discovery/candidates.py` (D14/D7, margin 0.10, positive test on distances, name never consulted, length band not the call) → all six now capped at 39 with the margin recorded. The margin is 31/31 correct under both identity metrics (true-ITPR +0.065…+0.874 vs RyR −0.868…−0.536). **Caveat: the deepest true member, *Dictyostelium* iplA, has a margin of +0.065 — inside the D7 no-call band**, so profile-based assignment is mandatory for the deep branches. One recall miss: fly `Itpr` at 35, its nearest sibling 0.342 vs the 0.35 breadth threshold (0.420 under `covered_only`). → `results/benchmark_controls/report.md` |
| S2 | Uncapped InterPro enumeration of the family Pfams → census v2, with a positive ITPR/RYR call on every record | S1 | pending | |
| S3 | Profile-HMM sweep (itpr.hmm + ryr.hmm) over vertebrate reference proteomes + jackhmmer-to-convergence completeness argument → census v3 | S1, S2 | pending | |
| S4 | Genome scope: declared assembly manifest (the denominator) + download tooling | S1 | pending | |
| S5 | Genomic tblastn + miniprot sweep → per-genome ledger (found / lost / assembly-gap) + novel gene models → census v4 | S4 | pending | |
| S20 | **Non-vertebrate sweep** — the family's true range across eukaryotic reference proteomes, and whether the land-plant / dikarya absence is a genome fact or a database fact | S3 | pending | |
| S23 | **Targeted invertebrate + protist + plant/fungal genome sweep** (triggered by S20): the negative claims at genome level, and copy number outside vertebrates | S20 | pending | |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl on census representatives per clade × kingdom, RyR outgroup included | S2, S3, S5, S20 | pending | |
| S7 | ML phylogeny: IQ-TREE 2 + ModelFinder + 1000 UFBoot, rooted on RyR; resolve which two of ITPR1/2/3 are sisters (AU test on the three rooted topologies) | S6, S23 | pending | |
| S8 | Synteny: shared flanking-gene analysis across the ITPR loci | S2, S5 | pending | |
| S9 | ML selection: PAL2NAL + codeml branch/site models, HyPhy RELAX | S6 | pending | |
| S10 | Annotation-bug molecular validation: assembly-version audit, RNA-seq spanning evidence, miniprot gene models, for the two worst cases S5/S18 surface | S5 | pending | |
| S11 | Structures: AFDB coverage, TM-align vs the cryo-EM IP3R and RyR references, Foldseek AFDB-wide sweep | S2, S6 | pending | |
| S12 | Expression evidence (SRA junction-spanning reads + atlases) for whichever paralog or lineage the census leaves in doubt | S5 | pending | |
| S13 | Gene-tree/species-tree reconciliation — date the duplications that made ITPR1/2/3 | S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3, S5, S7, S8, S9, S10, S11, S20, S23, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Loss dynamics** — ancestral-state reconstruction of per-paralog presence/absence, ORF-integrity screen from miniprot frameshift/stop counts, dating any pseudogene fossils | S5, S13 | high | pending | |
| S16 | **Duplication history** — are ITPR1/2/3 2R ohnologs; are teleost itpr1a/itpr1b from 3R; copy-number landscape and retention asymmetry | S7, S8, S13 | high | pending | |
| S17 | **Constraint & function** — per-site conservation mapped onto the cryo-EM channel; do the SCA15/SCA29/Gillespie, anhidrosis and neuropathy variants sit in the constrained core? Is the IP3-binding core more constrained than the pore? | S6, S9, S11 | high | pending | |
| S18 | **Annotation-quality audit** — how often a real ITPR locus is missing, fragmentary, split, unnamed or filed under the wrong paralog (or as a RyR) across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5, S15 | high | pending | |
| S19 | **Methods results** — per-method contribution ("what would proteome-only searching have missed?"), assembly contiguity as a confounder of loss claims, bait-panel design sensitivity | S5, S15, S18 | medium | pending | |
| S21 | **Gene architecture** — exon/intron structure across the genome scope from the miniprot CDS blocks; are database "fragments" real exon boundaries or annotation failures; is the ~58-exon architecture conserved | S5, S18 | medium | pending | |
| S22 | **Ligand-site evolution** — the IP3-binding core is the one part RyR does not share functionally. Is it under different constraint from the pore, does it differ between paralogs, and does it change in lineages that lost the upstream PLC/IP3 pathway | S9, S17 | medium | pending | |
| S14c | **Manuscript rewrite pass** — one full pass over the draft once every analysis has landed, with the figure audit's lessons applied | S14a, S24 | medium | pending | |

---

## Emergent tasks & new aims

Anything discovered mid-session that deserves its own work goes here rather
than expanding the task in progress.

| Added | From | Task | Status |
|-------|------|------|--------|
| 2026-08-18 | setup | Confirm the external drive is attached and `data_root.txt` points at it before S4/S5 | open |
| 2026-08-18 | setup | 40 Viridiplantae + 41 Fungi PF08709 records exist in a family textbooks say plants and fungi lack. Identify what they are (real gene / mis-annotation / contamination) — this is Q1's sharpest edge | open (S2 → S20) |
| 2026-08-18 | setup | ~~**Ensembl REST is unreliable right now**~~ — **diagnosed and fixed in S0.** Not general flakiness: `/xrefs/symbol/**homo_sapiens**/{symbol}` stalls indefinitely (no response, no error) for `BRCA2` as well as `ITPR1`, while the same endpoint answers in 0.6 s for `danio_rerio` and `/lookup/symbol/homo_sapiens/` answers normally. `src/databases/ensembl.py:_symbol_to_ids` now uses `lookup/symbol` with `xrefs` as fallback. Evidence: `results/s0_baseline/ensembl_endpoint_probe.tsv` | closed 2026-08-18 |
| 2026-08-18 | S0 | **Ensembl is slow enough to be a scheduling problem, separate from the stall.** Measured: 14 s for a 451-byte `lookup/symbol`, 11 s for a `lookup/id?expand=1`; one gene in one species costs ~95 s end to end. The default 8-species panel × 3 genes therefore needs ~38 min, so `run_headless`'s budget went 300 s → 900 s and a full panel sweep still needs `--species`. The real fix is to parallelise `EnsemblClient.search`'s per-species loop (the other clients already return in seconds) — a client change, not a session's worth of work, but out of S0's scope | open |
| 2026-08-18 | S0 | **Genomic span varies 6.5× across the three human paralogs (ITPR3 76 kb → ITPR2 498 kb) while protein length varies 3 %.** Found while correcting a false `[lit]` claim. Intron-content asymmetry between paralogs of identical architecture is a result, not a footnote — and D16 says the comparison must be paired within genome | open (S21) |
| 2026-08-18 | S0 | **Zebrafish `LOC101884734` (4,900 aa, 7 records) is an unnamed RyR-sized locus** returned by an ITPR-diagnostic Pfam query. First concrete instance of the unnamed-locus problem; keep it as a worked example for the correction list | open (S18) |
| 2026-08-18 | S0 | **ITPR3's clinical phenotype is broader than neuropathy** — the recurrent de novo p.Arg2524Cys causes a multisystemic disease with immunodeficiency. S17 must treat the ITPR3 variant set as multisystem, and Gillespie syndrome has **both** recessive and dominant-negative mechanisms, not only the latter | open (S17) |
| 2026-08-18 | S1 | **The discovery scorer's distance components run on full-alignment identity, which dilutes every comparison between proteins of unequal length.** Measured both ways on the S1 panel: RyR-vs-ITPR identity is 0.105–0.110 full-alignment (below the 0.20 novelty floor, so the twilight-zone component never fires) but 0.249–0.258 fragment-aware — *inside* the 0.20–0.40 zone, worth a further +20. The same dilution caused S1's only recall miss (fly `Itpr`: nearest sibling 0.342 vs the 0.35 breadth threshold, 0.420 fragment-aware). Switching `analyse()`/discovery to `identity_matrix(covered_only=True)` would fix the miss **and** raise every RyR to 65 — safe only because the sister test now caps them. A scorer re-weighting, so it is logged here rather than done silently | open (S6) |
| 2026-08-18 | S1 | **The bait margin cannot call the non-metazoan grade.** *Dictyostelium* `iplA`, a true family member, sits at +0.065 — inside the D7 10 % no-call band — while every metazoan positive is ≥ +0.256. The labelled-bait route is a vertebrate/invertebrate instrument; the deep branches need best-profile assignment (`itpr.hmm` vs `ryr.hmm`) before any presence/absence claim rests on them | open (S3 → S20) |
| 2026-08-18 | S1 | **The MSA-signature fallback, the fold component and the split-annotation component are all still untested.** Every S1 panel member had a real InterPro record, so the domain component always took the Pfam route; no Foldseek run and no split gene models were in a UniProt-only panel. Three of the eight scorer criteria therefore carry no validation, and the fallback is exactly what a Compara/BLAST-sourced candidate depends on | open (S2/S5) |
| 2026-08-18 | S1 | **`src/discovery/candidates.py` is 483 lines** — under the 500-line rule with 17 to spare. The next component added to the scorer must split the file (candidate scoring vs report rendering is the natural seam) | open |
| 2026-08-19 | review figures | ***Dictyostelium* iplA (Q9NA13) carries neither PF08709 nor PF02815 nor PF00520** — three of the five family signatures, including the one that *defines* the family. It is a characterised IP₃ receptor and it is in this project's own positive panel, so the InterPro enumeration of S2 would not return it. The census cannot be a single-domain query in either direction: PF08709 over-returns RyRs (49 % in zebrafish) **and** under-returns real members. S3's profile HMMs must be built to find it, and any absence claim resting on signature counts is unsafe until they are | open (S2 → S3 → S20) |
| 2026-08-19 | review figures | **The literature variant set is thinner than §9 reads.** Of the nine disease entries the review catalogues, only two carry a residue the cited source names (both *ITPR3*: p.Thr1424Met, p.Arg2524Cys); the rest are localised to a domain or not at all. S17's constraint analysis therefore cannot be built from the review — it needs a real variant table (ClinVar / gnomAD) with its own provenance | open (S17) |
| 2026-08-19 | review figures | **p.Arg2524Cys sits 7 residues past the measured gate** (Phe2513/Ile2517 in 6DQN). The recurrent multisystem variant is on the C-terminal load-bearing stretch of §2.4, which is a structural prediction S17/S12 can test rather than a coincidence to note | open (S17) |
| 2026-08-18 | setup | AlphaFold DB returned models for 8/8 human ITPR queries in the smoke test — better coverage than the PIEZO family had. Worth checking early whether AFDB covers full-length ITPRs or only fragments, since it changes S11's scope | open (S11) |

---

## Decisions log

**D0 — This project inherits the PIEZO project's methodological decisions.**
They were paid for over 24 sessions and are not to be re-derived. D3–D17
below are those rules, restated for this family. A session may overturn one,
but must say so explicitly here with its reason.

**D1 — Bulk storage lives on the external drive**, `data_root.txt`
(`/Volumes/FANTOM/IP3R_DATA`). Sessions call `require_data_root()`; an
unplugged drive stops the session rather than filling the internal disk.

**D2 — The repository stays private until S14b**, which flips it public
alongside the Zenodo DOI. Nothing in the repo may assume a public URL before
then. The remote was created after S1 (2026-08-19):
`git@github.com:gddickinson/ip3r_genes.git`, **visibility PRIVATE** —
sessions now end with `git push` as well as `git commit`, and the S14b step
that flips it public is `gh repo edit --visibility public`.

**D3 — The discovery scorer needs an evidence gate.** A score of ≥ 40 needs
at least one family-specific component (domain / outlier / fold / split
annotation); size and novelty alone cap at 39. Implemented in
`src/discovery/candidates.py`.

**D4 — An absence claim must pass two bars, not one.** A genome-wide
contiguity floor (calibrated in S19, ~50 kb contig N50 in the PIEZO project)
*and* a local check that the paralog's own genomic neighbourhood is present.
Residual false negatives are regional, not random: a genome missing a whole
syntenic block is not evidence of gene loss.

**D5 — Bait panels are screened by label, not padded for breadth.** One
correctly labelled bait per paralog above ~40 % identity recovers nearly
everything; a bait wearing the wrong clade name does real damage. Screen
every bait for chimeras (a domain-envelope check against the profile) before
it enters a panel.

**D6 — Never issue a database correction without the integrity veto.** If
the ORF-integrity screen (S15) says a locus is genuinely dead, the audit
does not tell RefSeq to resurrect it.

**D7 — Best-hit paralog assignment needs a margin.** Two baits scoring
within ~10 % of each other is not a call. This applies with extra force
here, where ITPR and RYR share every diagnostic domain.

**D8 — Representatives are chosen per clade × per kingdom**, with explicit
selection rules in a script, not "the longest sequence per species". A
longest-first pick reliably selects chimeric gene models.

**D9 — Quote an annotation claim with its annotation source.** RefSeq
(`GCF_`) gene sets and submitter-deposited GenBank (`GCA_`) gene sets are
not comparable evidence; some assemblies ship no gene set at all.

**D10 — Iterative searches need a coded kill criterion.** jackhmmer runs
that diverge must be killed by a rule recorded in code and reported, not by
eye.

**D11 — Look at the figure.** Any session that changes a figure looks at it;
any session that writes a legend looks at the figure it describes. Four
errors in the PIEZO manuscript were of a kind no table check could catch.

**D12 — Every load-bearing number in the manuscript needs a claim row** in
`scripts/s14_claims.py`, naming its source table and the operation that
recovers it. A re-run that changes a table then fails loudly.

**D13 — Reports are rendered from the committed tables**, never written by
hand alongside them, so a report and its data cannot drift.

**D14 — ITPR vs RYR is a positive test at every stage** (new, family-specific).
Never assume a Pfam hit, a BLAST hit or a gene model is an ITPR. Assign by
best profile (`itpr.hmm` vs `ryr.hmm`) or by a labelled-bait margin, record
the margin, and treat the length band (2,000–3,600 aa) as a filter that
supports the call rather than as the call itself.

**D14a — D14 is implemented as a labelled-bait margin in the scorer**
(S1, 2026-08-18). `DiscoveryConfig.sister_paralogs` / `sister_margin` (0.10,
D7's number) / `exclude_sister_family`: a candidate whose identity to the
nearest labelled sister bait exceeds its identity to the nearest known
paralog by more than the margin is assigned to the sister family and capped
at 39, with the margin written into its evidence. It is a positive test on
distances — the candidate's own gene symbol is never consulted, so unnamed
RyR-sized loci are called the same way — and the length band contributes
only the size component, never the call. It needs at least one *labelled*
sister bait in the analysis set to measure against; a search that omits the
RyR panel silently loses the test. Measured: without it, 6/6 RyR decoys
promoted at 45.

**D19 — The review's figures are generated, never drawn, and every one
declares its provenance.** `scripts/s0_review_figures.py` renders all 12
from committed tables under `results/s0_baseline/review_figures/`; no
figure module queries a database or reads a structure, so a figure cannot
disagree with the data behind it (D13 applied to pictures). Each carries a
corner tag — **measured** (a live database or structure), **computed**
(derived here from committed sequences), **schematic** (a drawing of a
cited mechanism, no data, not to scale) or **curated** (transcribed from
cited literature) — which is the `[db]`/`[lit]`/`[open]` discipline of
`ip3r_background.md` carried into the figures. Two further rules follow.
A curated row must carry the review's stable citation keys and
`s0_figdata_curated.py` fails if one has no reference row. And a variant
is drawn only at the resolution its source gives — `point` / `domain` /
`gene` — because a figure that promotes a domain-level claim to a residue
is inventing data, which is precisely what this review exists not to do.

**D18 — The project runs in the reused `piezo1` conda env**
(`/opt/anaconda3/envs/piezo1`, python 3.11.15), not a cloned `ip3r` env.
BLAST+ 2.16.0+, the NCBI `datasets` CLI and Foldseek live there rather than
on the bare PATH; MAFFT, HMMER, miniprot, trimAl and IQ-TREE 2 come from
Homebrew and work in any shell. Reused rather than cloned to avoid
duplicating several GB of an identical toolchain — the env holds tools, not
project data, so nothing about a result depends on its name. Exact versions:
`results/toolchain_manifest.txt`, regenerated by `scripts/s1_toolchain.py`.

**D15 — The species tree is an input, not a result.** Reconciliation uses a
hand-curated, literature-calibrated topology with a source on every
calibrated node, validated by a `--check` mode before use.

**D16 — Cross-paralog comparisons are paired within genome.** Intron size,
assembly quality and annotation completeness all scale with the assembly, so
an unpaired comparison measures the assemblies.

**D17 — Permutation nulls are drawn from real genomic windows**, not from a
uniform shuffle, and the genes under test are removed from their own windows
or the test is circular.
