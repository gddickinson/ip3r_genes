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

> **Next session: S5b.** S0–S4 and S5a are complete. Read
> `results/genome_ledger/report.md` before S5b — it carries the measured
> budget (303 genomes, 6.7 h compute, 164 GB download), the two constants
> S5a had to re-derive, and **D4's contiguity bar, which 120 of the 309
> manifest genomes fail**. The bar is the thing to hold on to: the margin
> set and the low-contiguity set are very largely the same genomes, so no
> bird absence is readable as loss until they are separated.
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
| S2 | Uncapped InterPro enumeration of the family Pfams → census v2, with a positive ITPR/RYR call on every record | S1 | completed 2026-09-03 | **15,421 proteins across 1,488 taxa; ITPR 6,433 / RYR 6,807 / unassigned 2,181 (14.1 %).** All three seeds walked their cursor chains to the end (63 + 67 + 60 pages, 0 restarts). **InterPro's advertised `count` is wrong in both directions** — PF08709 +168, PF01365 +56, PF08454 −176 — while what it *serves* matches UniProt's independent count to within 2 records every time (D20); the first version of this task's own completeness test failed PF08454 on a complete walk because of it. The documented API host answered 1 request in 12 during the run while `/interpro/wwwapi/` answered 12 of 12, so the client now fails over between them. **The call is a positive architecture test audited against a label it never sees (D14b): 6,191 symbol-labelled records, 0 disagreements**, and re-checked on sequence alone over the per-phylum core panel — 27/27 and 26/26 agreement wherever the bait margin decides, with all four sign flips inside D7's no-call band. **A PF08709-only census would have missed 2,914 records, 758 of them called ITPR across 385 taxa** (D21) — including *Dictyostelium* iplA, which the union recovers and correctly leaves `unassigned` at 2/5 signatures. **The plant and fungal records resolve phylogenetically**: every Viridiplantae call is Chlorophyta (Streptophyta 0/15 records, 13 taxa) and every fungal call is in an early-diverging phylum (Dikarya contributes no records at all). Delta vs census v1 ×9.8 with **no enumeration holes** — all 866 absences verdicted. 7 records carry the complete architecture in 1,528–1,993 aa: truncated gene models UniProt does not flag. → `results/census_v2/report.md` |
| S3 | Profile-HMM sweep (itpr.hmm + ryr.hmm) over vertebrate reference proteomes + jackhmmer-to-convergence completeness argument → census v3 | S1, S2 | completed 2026-09-04 | **Census v3: 16,039 records — ITPR 8,000, RYR 7,432, unassigned 605, conflict 2 — with two independent verdicts on every row (D23).** Two profiles built from the S2 census under enforced selection rules (34 ITPR / 22 RyR seeds, 2,684 / 4,930 match states) and swept over **763 vertebrate reference proteomes, 14,414,821 canonical proteins, 6.97 G residues**. **The instrument was calibrated before use**: profile assignment vs S2's architecture call over the whole seeded space agrees **11,875 / 11,876** (seeds excluded), and resolves **2,314 of the 3,361 records the architecture rule could not call**; S2's gene-symbol fallback on 1,219 records is overturned exactly once. **The one genuine disagreement is a correction**: *Tieghemostelium lacteum* A0A152A7I8 (2,845 aa), called RYR by S2 on PF06459 — "Ryanodine Receptor TM 4-6", which is the pore both families share — beaten 303 to 133 bits by `itpr.hmm`. That is D14b's own conditional-specificity caveat firing, once. **The sweep needed a gate the roadmap did not have (D22)**: the first run called **14,981 vertebrate proteins RYR**, 12,874 under a tenth of the profile and named *Tnnc2* / *Rspry1* / *Cabp1* / *Rnf123*, because `ryr.hmm` carries SPRY and `itpr.hmm` has nothing to match it. The floor is **200 match states**, measured from this project's own S0 coordinates (shortest PF08709 200 aa, longest SPRY 137 aa); it costs 86 of 12,020 records their profile verdict and keeps split gene models. **Completeness measured in both directions**: of 2,787 v2 ITPR records from a swept proteome the sweep did not return, **zero were in the database and missed** — all 2,787 are UniProtKB entries the reference proteomes do not contain; and the sweep adds **618 proteins the InterPro census never returned**, all 209–942 aa fragmentary gene models. **jackhmmer, all three seeds** (the *Acanthamoeba* run folded in 2026-09-04 after 7.6 h, against 2.4 h and 2.5 h): `itpr_fly` converged in 10 rounds (D10 clean); `itpr1_human` and `itpr_acanthamoeba` both reached the ceiling (K3). **The family core is seed-independent** — a vertebrate, an insect and an amoebozoan seed recover *identical* family content: ITPR/architecture **1,656 and ITPR/partial 952 in all three models**, RYR/architecture 1,565–1,569. What differs is everything that is not the family: the *Acanthamoeba* model rests on 26,266 targets against human's 7,222, and **18,670 of them are proteins neither profile scores at all** (7 in each of the other two runs). **That drift is invisible to K1** — off-family accretion *dilutes* the sister-family share instead of raising it, and it falls 27 points across this run — so K3 is the rule that caught it (**D10b**). Module-only matches are 10–34 % of these models, which is why they asymptote rather than reach zero. **Three method fixes**: MAFFT is not reproducible at `--thread -1` (**D24**), K1 measured sister-family *level* where only its *rise* is evidence (**D10a**), and K3 punished a run that converged on its last round. → `results/census_v3/report.md` |
| S4 | Genome scope: declared assembly manifest (the denominator) + download tooling | S1 | completed 2026-09-04 | **309 genomes: 161 vertebrate orders ∪ 169 margin species, 552.5 Gbp (≈553 GB FASTA, ≈166 GB as zip, against 1,434 GB free).** Scope confirmed with the user: the *full* margin union, so the bird question is settled per species rather than by sample. **The margin species are derived from the committed census tables, not hand-listed** (unlike the PIEZO port) — four rules with stated thresholds in `s4_manifest_lib.derive_margins`: `zero_hit_proteome` 15, `missing_paralog` (<3 paralogs) 101, `fragment_only` (longest record < `family.MIN_LENGTH_AA` = 2,000 aa) 100, `anchor` 5; 61 genomes carry more than one reason and 11 carry three. All 169 resolved to an assembly — **no margin species is excluded for want of a genome**. **The margin is a bird problem**: 130 of 169 margin species are Aves, against 21 Actinopteri and 10 Mammalia, which is S3's annotation-depth result reappearing as a scope requirement. **D9 is live in the denominator**: 168 RefSeq / 141 GenBank, and **33 genomes carry no gene set at all**, so no annotation claim may be quoted for them. Three assemblies exceed 10 Gbp (*Protopterus annectens* 40.1, *Lissotriton helveticus* 23.2, *Bombina bombina* 10.0) and dominate both download and S5 runtime. 67 reference species have no NCBI order rank (65 percomorphs + 2 sharks) and cannot represent one; the margin rules still reach them. `fetch_genomes.py` tested end to end on 3 small genomes — md5-verified against NCBI's own `md5sum.txt`, `.done` resume, fetch → search → purge, and a failing search correctly keeping its genome. **Two bugs found by testing**: margin species were matched by a last-wins taxid lookup while order reps used the rank function, putting one species in twice under two accessions; and `--verify` skipped any checksummed file the install does not keep, so a deleted `.fna` would have passed — the kept checksums are now rewritten to exactly the installed files and a missing one is a hard failure (proved by deleting a file and watching it fail, exit 1). → `results/genome_manifest_notes.md` |
| S5a | **Genomic sweep: bait panel, calibrated pipeline, pilot** — the instrument for S5, measured rather than ported | S4 | completed 2026-09-04 | **38 baits (30 ITPR + 8 RyR control), 121,294 residues, derived from census v3 by seven enforced rules; pipeline piloted on 6 genomes, 24 cells, 43 loci, 0 control failures.** Three calibrations that a port could not supply. **(1) miniprot's `-G`.** A too-small max-intron does not lose a gene, it *splits* one, and a split ITPR reads out of the ledger as `fragment` — so `-G` is a threshold on the call. Measured across an 11-species Ensembl panel: **no ITPR gene has an intron over the 200 kb default** (widest 152,216 bp, human ITPR1) while the RyR control exceeds it twice (RYR2 227,927 bp). The rule scales from the widest measured intron rather than from the intron-per-Gbp ratio, which is largest in the *smallest* genomes and would have asked for a 7.4 Mbp `-G` on the lungfish. **(2) The unlabelled baits were competing as a fourth paralog.** Rescue attribution ranked `vertebrate_basal` (the gar/chimaera/lamprey seeds — evidence that a region is an ITPR, not a hypothesis about *which*) against the real paralogs: in *Todus* an ITPR1 trace at a true 0.235 margin over ITPR2 was reported at 0.027. Three cells moved from `tblastn_trace_ambiguous` to `tblastn_trace` once fixed. **(3) `ATTRIBUTION_REL_MARGIN` does not transfer (D25).** The PIEZO port's 0.333 comes from a family 40–50 % identical where these are 61–68 %; at complete loci whose paralog identity the annotation independently establishes, **7 of 9 sit below it** (range 0.225–0.347, median 0.281). A complete locus is an upper bound on a rescue fragment, so under 0.333 no absence claim could ever be attributed → **0.22**, with the limitation stated. **D14 was never in doubt at genome scale: at all 43 loci only one family's baits aligned at all** — sharper than the protein level, where S1 had all six RyR decoys promoted at 45. **The RyR control paid immediately**: in *Todus mexicanus* (the `zero_hit_proteome` bird) it found one RyR locus where a bird has three, and unannotated — what is missing there is contiguity, not genes. **D4's contiguity bar, set at the median measured ITPR span (142,212 bp) rather than an invented number: 120/309 manifest genomes fail it — 66 % of birds against 11 % of fish, and 68 % of margin species against 12 % of order reps.** The margin set was chosen for what its *proteomes* lack and turns out to be very largely the same genomes whose assemblies cannot hold the gene; the two signals are confounded. **And the paralog recovered from a fragmented assembly is systematically the short one** — below the bar ITPR3 (82 kb median span) is found 2/3 while ITPR1 (186 kb) and ITPR2 (231 kb) are 0/3, so S0's 6.5× span asymmetry returns as a detection bias pointing the same way as the loss signal *(pilot of 6 — direction and mechanism, pending: S5b at 309)*. Screen validated by its own negative controls on every build. → `results/genome_ledger/report.md` |
| S5b | **The full 309-genome sweep** → per-genome ledger (found / lost / assembly-gap) + novel gene models → census v4 | S5a | in_progress | Budget measured, not estimated: **303 genomes, 547 Gbp, 6.7 h compute (miniprot 36 s/Gbp, 8 threads), 3.4 h wall at 2 shards, 164 GB download**; 547 GB of FASTA fits in 1,427 GB free, so keep it for S8/S10. Two giants (*Protopterus* 40 Gbp, *Lissotriton* 23 Gbp) need the chunked path and are outside the driver's 12 GB default cut. Re-calibrate `ATTRIBUTION_REL_MARGIN` on rescue regions of known identity, and test the span-bias result at scale. |
| S20 | **Non-vertebrate sweep** — the family's true range across eukaryotic reference proteomes, and whether the land-plant / dikarya absence is a genome fact or a database fact | S3 | pending | |
| S23 | **Targeted invertebrate + protist + plant/fungal genome sweep** (triggered by S20): the negative claims at genome level, and copy number outside vertebrates | S20 | pending | |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl on census representatives per clade × kingdom, RyR outgroup included | S2, S3, S5b, S20 | pending | |
| S7 | ML phylogeny: IQ-TREE 2 + ModelFinder + 1000 UFBoot, rooted on RyR; resolve which two of ITPR1/2/3 are sisters (AU test on the three rooted topologies) | S6, S23 | pending | |
| S8 | Synteny: shared flanking-gene analysis across the ITPR loci | S2, S5b | pending | |
| S9 | ML selection: PAL2NAL + codeml branch/site models, HyPhy RELAX | S6 | pending | |
| S10 | Annotation-bug molecular validation: assembly-version audit, RNA-seq spanning evidence, miniprot gene models, for the two worst cases S5b/S18 surface | S5b | pending | |
| S11 | Structures: AFDB coverage, TM-align vs the cryo-EM IP3R and RyR references, Foldseek AFDB-wide sweep | S2, S6 | pending | |
| S12 | Expression evidence (SRA junction-spanning reads + atlases) for whichever paralog or lineage the census leaves in doubt | S5b | pending | |
| S13 | Gene-tree/species-tree reconciliation — date the duplications that made ITPR1/2/3 | S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3, S5b, S7, S8, S9, S10, S11, S20, S23, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Loss dynamics** — ancestral-state reconstruction of per-paralog presence/absence, ORF-integrity screen from miniprot frameshift/stop counts, dating any pseudogene fossils | S5b, S13 | high | pending | |
| S16 | **Duplication history** — are ITPR1/2/3 2R ohnologs; are teleost itpr1a/itpr1b from 3R; copy-number landscape and retention asymmetry | S7, S8, S13 | high | pending | |
| S17 | **Constraint & function** — per-site conservation mapped onto the cryo-EM channel; do the SCA15/SCA29/Gillespie, anhidrosis and neuropathy variants sit in the constrained core? Is the IP3-binding core more constrained than the pore? | S6, S9, S11 | high | pending | |
| S18 | **Annotation-quality audit** — how often a real ITPR locus is missing, fragmentary, split, unnamed or filed under the wrong paralog (or as a RyR) across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5b, S15 | high | pending | |
| S19 | **Methods results** — per-method contribution ("what would proteome-only searching have missed?"), assembly contiguity as a confounder of loss claims, bait-panel design sensitivity | S5b, S15, S18 | medium | pending | |
| S21 | **Gene architecture** — exon/intron structure across the genome scope from the miniprot CDS blocks; are database "fragments" real exon boundaries or annotation failures; is the ~58-exon architecture conserved | S5b, S18 | medium | pending | |
| S22 | **Ligand-site evolution** — the IP3-binding core is the one part RyR does not share functionally. Is it under different constraint from the pore, does it differ between paralogs, and does it change in lineages that lost the upstream PLC/IP3 pathway | S9, S17 | medium | pending | |
| S14c | **Manuscript rewrite pass** — one full pass over the draft once every analysis has landed, with the figure audit's lessons applied | S14a, S24 | medium | pending | |

---

## Emergent tasks & new aims

Anything discovered mid-session that deserves its own work goes here rather
than expanding the task in progress.

| Added | From | Task | Status |
|-------|------|------|--------|
| 2026-08-18 | setup | Confirm the external drive is attached and `data_root.txt` points at it before S4/S5 | open |
| 2026-08-18 | setup | ~~40 Viridiplantae + 41 Fungi PF08709 records exist in a family textbooks say plants and fungi lack~~ — **narrowed by S2, not closed.** The records are not scattered: in Viridiplantae *all* 20 ITPR calls are **Chlorophyta** (11 taxa, incl. *Chlamydomonas reinhardtii* with the complete five-signature architecture) and **Streptophyta — the land-plant lineage — has 0 calls from 15 records in 13 taxa**; in Fungi every call is in an early-diverging phylum (Mucoromycota 15, Chytridiomycota 6, Basidiobolomycota 3, Entomophthoromycota 1) and **Dikarya contributes no records to the search space at all**. That is the shape of a loss in the derived lineage of both kingdoms, and it is still a statement about what UniProt holds. → `results/census_v2/lineage_calls.tsv` | open (S20 → S23) |
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
| 2026-09-03 | S2 | **Seven records carry the complete five-signature ITPR architecture in 1,528–1,993 aa** — 700+ residues short of the shortest real family member — and none is flagged `Fragment` by UniProt, because a truncated gene model submitted as a whole protein is not marked as one. All are unnamed locus tags, five from *Hymenochirus boettgeri* (two loci) and two from chironomid midges. The call on them is correct and the records are wrong: a starting list for the correction register, and the reason length stays a recorded column after it stopped being part of the call. → `results/census_v2/short_complete.tsv` | open (S18) |
| 2026-09-03 | S2 | **2,177 records (14.1 %) are `unassigned`** — a partial architecture cannot be called by a rule that reads absence as evidence. They are a fragment population (median 678 aa vs 2,671 for a called ITPR; only 5 % inside the size band), so this is the intended behaviour rather than a shortfall, but it is 14 % of the census and S3's profile sweep is what resolves it. The deepest true members are in here: *Dictyostelium* iplA carries 2 of 5 signatures and is `unassigned` | open (S3) |
| 2026-09-03 | S2 | **The seeded search space holds one bacterial record** (Bacteroidota, no ITPR call). Almost certainly a horizontal-transfer or contamination artefact rather than a real prokaryotic family member, but it is the only prokaryote any of this project's searches has returned and S20's negative claims should name it rather than be surprised by it | open (S20) |
| 2026-09-04 | S3 | **A deep-branch seed is 8.5× more expensive to search with, and the cost is measurable up front.** The *Acanthamoeba* ITPR passes **17.0 % of the 14.4 M-protein database through HMMER's MSV filter** against an expected 2.0 %, so 2.46 M sequences reach the expensive Viterbi/Forward stages; the human and fly seeds behave normally. Its round-1 search took longer than the other two seeds' entire runs. S20 sweeps non-vertebrate proteomes with exactly this kind of seed, so its runtime budget should assume the filter is near-useless for protist queries — and the filter pass rate is printed in every log, so it can be checked after one round rather than discovered after a day | open (S20) |
| 2026-09-03 | S3 | **The one record the two instruments genuinely disagree on is *Tieghemostelium lacteum* A0A152A7I8** — 2,845 aa, ITPR-sized, which S2 called RYR because it carries PF06459 ("Ryanodine Receptor TM 4-6") and which `itpr.hmm` beats `ryr.hmm` on by 303 to 133 bits. PF06459 is the RyR *transmembrane* module — the part the two families share structurally — so this is D14b's own conditional-specificity caveat firing on a real amoebozoan receptor. One record in 11,965, and it is a genuine correction rather than noise. The second conflict, *Symbiodinium* A0A1Q9CN86, is a 9,504-aa protein where only 11 % of the target aligns, and should not be called by either instrument | open (S20) |
| 2026-09-03 | S3 | **1,794 sweep hits carry an ITPR or RYR gene symbol but fall under D22's 200-position gate** (`results/census_v3/subthreshold_family_fragments.tsv`) — pieces of split or truncated gene models, headed by `Ryr3_1` (224), `Itpr2_0` (212), `Itpr1_0` (198). They are real family genes whose annotation has been broken into fragments too short to span a family domain. This is the largest concrete lead the project has for the annotation-quality audit, and it is a *per-proteome* count, so S5's genome sweep can test each one against the DNA | open (S18 → S10) |
| 2026-09-03 | S3 | **604 proteins in vertebrate reference proteomes are absent from the InterPro census entirely** (`results/census_v3/novel_hits.tsv`) — 188 called ITPR, 416 RYR, every one of them a partial-evidence record 209–942 aa long in species like *Eptatretus burgeri* (hagfish), *Pleuronectes platessa* and *Myotis davidii*. A signature census misses fragmentary gene models because a fragment carries too few signatures to be enumerated; a profile finds them because it does not need the annotation. Quantifies what "annotation-derived census" costs | open (S18, S19) |
| 2026-09-03 | S3 | **15 of 758 swept vertebrate taxa have no ITPR record at all**, and their gene sets are small (median 10,042 proteins vs a modal well-annotated proteome). On this evidence they are thin annotations, not losses — but they are the concrete starting list for the genome sweep, and D4 requires both bars before any of them is called an absence. → `results/census_v3/proteomes_without_hits.tsv` | open (S4 → S5) |
| 2026-09-04 | S5a | **The contiguity bar and the margin rules select very largely the same genomes.** 120 of the 309 manifest assemblies have a contig N50 below the median measured ITPR genomic span (142,212 bp) and so cannot carry the gene on one contig — but not at random: **66 % of Aves against 11 % of Actinopteri, and 68 % of margin species against 12 % of order representatives.** S4 chose the margin set for what its *proteomes* lacked; it turns out to be mostly the set whose *assemblies* cannot hold the gene. Every downstream absence claim in birds is confounded until the two are separated, and S19's contiguity floor is now a prerequisite for S15 rather than a methods footnote. → `results/genome_ledger/contiguity_bar.tsv` | open (S19 → S15) |
| 2026-09-04 | S5a | **Recovery from a fragmented assembly tracks gene span, so the detection bias points the same way as the loss signal.** Below the contiguity bar the pilot finds ITPR3 (median span 82 kb) in 2 of 3 genomes and ITPR1 (186 kb) and ITPR2 (231 kb) in 0 of 3; above it, 9 of 9. S0's 6.5× span asymmetry between paralogs of near-identical protein length was filed as a curiosity — it is a systematic bias in which paralog a poor assembly appears to have lost. Six genomes is a direction and a mechanism, not a result; S5b tests it at 309 and S19 must model it. | open (S5b → S19 → S15) |
| 2026-09-04 | S5a | **The paralog labels run out below the well-annotated clades.** Six bait slots could not be filled from the whole census — ITPR1 and ITPR2 in chondrichthyans, ITPR2 in the coelacanth grade, all three in cyclostomes — because no labelled, full-length, complete-architecture record exists there at all. The `vertebrate_basal` baits cover those genomes without a paralog claim, so the sweep finds the loci but cannot name them; S7's phylogeny is what assigns them, and until it does, no per-paralog statement may be made about sharks, lampreys, hagfish or the coelacanth. → `results/s5_baits/bait_build_stats.json` | open (S7) |
| 2026-09-04 | S5a | **Teleost ITPR1 3R duplicates are half-named in RefSeq.** Both *Takifugu rubripes* and *Genypterus blacodes* carry two ITPR1 loci, one annotated `itpr1b` and the other left as an unnamed `LOC` — in *Takifugu*, `LOC101074739`. Two RyR loci in the same genome are unnamed the same way. A name-based census cannot reach these, which is why the sweep carries every locus whose annotation does not already name it correctly into census v4. A concrete, paired worked example for the correction list. | open (S18) |
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


**D25 — A threshold ported from the PIEZO project must be re-measured before
it is used, not merely restated in this project's units.** S5a inherited
`ATTRIBUTION_REL_MARGIN` as that project's 1.5x bit-score ratio expressed as
D7's 0.333. Restating it made it *look* derived while leaving it tuned to a
family whose paralogs are 40-50 % identical, where ITPR1/2/3 are 61-68 % (S1).
Measured here, 7 of 9 complete loci of independently-known paralog identity
fall below it — so the inherited value would have reported every rescue trace
ambiguous and made every absence claim unattributable. D0 says the PIEZO
decisions are not to be re-derived; this is the boundary of that rule. A
*method* ports. A number that encodes how far apart that family's paralogs sit
does not, and neither does one that encodes its gene sizes: the same session
found miniprot's `-G` needed measuring for the same reason. Where a ported
constant cannot be re-measured yet, it is used with its limitation stated in
the code and the per-row value recorded so the cut can be revisited without
re-running.

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

**D14b — At census scale D14 is a positive architecture test, audited
against a label it never sees** (S2, 2026-09-03). A record is called RYR
because it carries a RyR-specific signature (PF02026 / PF06459 / PF21119 /
PF00622), and ITPR because it carries the *complete* IP3-receptor
architecture (PF08709 + PF01365 + PF08454 + PF02815 + PF00520) and none of
them. Absence of RyR evidence is allowed to count **only** when the ITPR
architecture is complete, because a RyR annotated well enough to show all
five shared signatures would also show its own; a partial record stays
`unassigned` and is counted. Length is on every row and enters only as
support for a medium-confidence call. Two things make this auditable rather
than assumed: gene symbols are scored against the architecture call but
never fed into it (`rule_audit.tsv` — 6,191 labelled records, 0
disagreements), and the whole call is re-checked on sequence alone over the
per-phylum core panel with the labelled-bait margin. Note the rule's one
conditional claim: PF00622 (SPRY) sits in ~114,000 UniProt proteins and is
in no sense RyR-specific — it is diagnostic *inside this search space*, and
that is what the audit row tests.

**D20 — Enumeration completeness is cursor exhaustion plus an independent
second count, never the API's advertised total** (S2, 2026-09-03).
InterPro's `count` field for PF08709 reports 12,339 while the same endpoint
serves 12,507 distinct accessions, stably across re-queries; UniProt's
independent count for the same signature is 12,506. A census that stopped
at the advertised number would have dropped 168 proteins silently. Record
all three numbers, treat the cursor chain running out as the completeness
test, and archive every raw page so the parse can be redone offline. The
same decision covers host failover: InterPro answers on `/interpro/api/`
and `/interpro/wwwapi/`, which fail independently — during S2 the
documented host answered 1 request in 12 and the website's host 12 of 12 —
so every request tries both before it sleeps.

**D21 — A domain census is a union of signatures, never a query on the
defining one** (S2, 2026-09-03). Measured: PF08709, the IP3-binding core
that names the family, is absent from 2,911 of the 15,417 enumerated
records, 758 of which this census calls ITPR across 385 taxa — among them
*Dictyostelium* iplA (Q9NA13), a characterised receptor in the project's
own positive panel, which carries PF01365 and PF08454 and none of PF08709 /
PF02815 / PF00520. Only 66.5 % of the search space carries all three seeds.
Conversely a single-signature census over-returns: 49 % of zebrafish
PF08709 records are RyRs (S0). Both directions fail, so the seed set is
declared in `src/utils/family.py:CENSUS_PFAM_IDS` and MIR is deliberately
excluded from it — carried by the O-mannosyltransferases too, it widens the
space without adding evidence, and is kept as an annotation column.

**D24 — Anything a committed artefact is built from is aligned
single-threaded** (S3, 2026-09-03). MAFFT L-INS-i with `--thread -1` is not
reproducible: the same 22 RyR seeds aligned twice on this machine gave 8,510
and 8,468 columns, and the profiles built from them 4,933 and 4,908 match
states, because the iterative refinement stage combines partial results in
whatever order the threads finish. Discovered by rebuilding the profiles
after an unrelated edit and finding the match-state count had moved. A
profile that changes when it is rebuilt cannot be the profile a committed
result was produced with, so `scripts/s3_build_seed.py` pins `--thread 1`
(byte-identical across three runs, 87 s instead of 15 s, once) and records
the SHA-256 of every seed set, alignment and profile in
`seed_build_stats.json` so a future drift is visible in the committed data
rather than only in a count. This applies to S6's alignment upgrade too.

**D10a — For this family, sister-family *level* is not evidence of drift;
only its *rise* is** (S3, 2026-09-04). D10's kill criterion was first coded
as a flat ceiling — kill a jackhmmer run when more than 5 % of the model's
included targets are assigned to the sister family. It fired on every run at
round 1, because a single human ITPR1 sequence searched at E ≤ 1e-5 already
returns **32 % ryanodine receptors before any iteration has happened**. That
is not contamination: ITPR and RyR are genuine homologues sharing the entire
pore, and a search sensitive enough to find *Acanthamoeba* is necessarily
sensitive enough to find RYR1. The level is a fact about the two families'
shared ancestry; only the change in it can be attributed to iterating the
model. K1 now takes round 1 as the baseline and kills on a rise of more than
10 percentage points. Measured: across the accepted runs the sister share is
flat or falling (−2.2 to +0.0 points), which is what a run that has not
drifted looks like — and the threshold was reasoned from what round 1 *is*,
not fitted to the observed numbers.

**D10b — K1 is blind to off-family drift; K3 is what catches it** (S3,
2026-09-04, on the third jackhmmer seed). D10a made K1 watch the *rise* in
sister-family share. The *Acanthamoeba* seed then diverged in the one
direction that rule cannot see: it accreted 18,670 targets neither profile
scores at all, which **diluted** the RyR share rather than raising it — the
K1 statistic fell 27 points while the model grew 5× and its family content
did not change. A run can therefore be badly divergent and read as
maximally clean on K1. Two consequences, both already in code: the ceiling
rule K3 is not a formality but the only rule that fires on this failure
mode, and `jackhmmer_model_composition.tsv` — what each final model is
actually built from, by profile call × evidence class — is reported for
every run rather than folded into a pass/fail. Do not replace K3 with a
"converged or not" test in S20's non-vertebrate sweep.

**D22 — A profile match is family evidence only if it spans a family
domain** (S3, 2026-09-03). Both profiles are full-length channel models, so
a protein sharing one small module with either of them scores against it —
and `ryr.hmm` carries SPRY (PF00622), which D14b already recorded as sitting
in ~114,000 UniProt proteins and being "in no sense RyR-specific". Measured:
without a gate the first run of the vertebrate sweep called **14,981
proteins RYR, 12,874 of them matching under a tenth of the profile** and
named *Tnnc2*, *Rspry1*, *Cabp1*, *Rnf123*, *Ash2l* — EF-hand and SPRY
proteins, not receptors. So the winning profile must span at least **200
match states** before its score counts. The number is measured rather than
tuned, from this project's own S0 domain coordinates: the shortest observed
instance of PF08709 — the IP3-binding core that names the family — is 200
aa, and the longest observed SPRY is 137 aa, so the floor sits in that gap.
It admits any match spanning a family domain, excludes any match spanning
only the shared module, and keeps split gene models (the 200–300 position
band is almost entirely `Itpr1_0` / `Itpr2_1` / `Ryr3_0`, pieces of real
genes). Cost on the calibration set: 86 of 12,020 architecture-called
records lose their profile verdict (0.7 %). Every row carries an `evidence`
class — `architecture` (≥ 50 % of the profile), `partial`, `module` — so the
gate is visible in the table rather than only in the code.

**D23 — Census v3 states two verdicts per record and never silently
reconciles them** (S3, 2026-09-03). The architecture rule (D14b) reads
InterPro's annotation; profile assignment reads the residues. Where they
agree the call is `high`; where only one speaks the row says which; where
they disagree the record is `conflict`, kept and reported. This is what
makes the 3,353 records S2 could not call resolvable — a partial
architecture defeats a rule that reads absence as evidence and does not
defeat a sequence profile — and it is what stops the second instrument from
being used to quietly overwrite the first.

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
