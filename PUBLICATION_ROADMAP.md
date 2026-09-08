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

> **Next session: S7.** S6 is complete: `results/msa_v2/` holds the
> alignment every downstream result stands on — 134 representatives, 1,790
> trimmed columns, 96.6 % of them parsimony-informative, with an audit table
> saying why each tip is in it. Read `results/msa_v2/report.md` first. Four
> things S7 inherits. **The outgroup and the pre-2R grade are in by rule**,
> so the tree can be rooted and the duplications placed. **The sister
> question has a preview and not an answer**: ITPR1 × ITPR2 leads at 0.788
> with a non-overlapping interquartile range — carry it in as the hypothesis
> the AU test is run against, not as a result. **The cyclostome trio is the
> sharpest 2R test available**: lamprey, hagfish and *Myxine* each carry
> three loci, all six fall nearest ITPR1, and the control shows the lean is
> real rather than a property of the metric — but identity cannot say
> whether they are 1:1 orthologs or a lineage-specific expansion, and S7 +
> S8 must. **A tip's group is not evidence**: the R7 novel models carry
> their bait's cell as a working hypothesis, and the tree is what tests it.

| ID | Task (one session each) | Depends | Status | Results (headline) |
|----|-------------------------|---------|--------|--------------------|
| S0 | Literature baseline + scope confirmation: verify every `[lit]` claim in `docs/ip3r_background.md`, build `docs/ip3r_review_2026.md`, smoke-test the app against live APIs | — | completed 2026-08-18 | **19 `[lit]` claims audited vs 51 refs (45 primary): 12 verified, 3 qualified, 2 struck, 1 → `[open]`, 1 → `[db]`.** `docs/ip3r_review_2026.md` written. All `[db]` numbers re-derived exactly. **D14 measured: 49 % (53/109) of zebrafish PF08709 records are RyRs.** Exon/span claim false — ITPR3 spans 76 kb, not "hundreds of kb"; span varies 6.5× across paralogs, length 3 %. Ensembl `/xrefs/symbol/homo_sapiens/` stalls (species-specific, `BRCA2` too); client switched to `/lookup/symbol/` → human ITPR1 **0 → 24 variants**. Smoke 3/3 on human (UniProt 34 / NCBI 44 / Ensembl 41) and on zebrafish; the 8-species panel is 2/3 on Ensembl latency alone. → `results/s0_baseline/report.md` **Extended 2026-08-19 (user request): the review is now illustrated** — 12 generated figures, 24 → 32 pages, rendered by `scripts/s0_review_figures.py` from committed tables only (D13, D19). Measuring the structure for the figures produced three results the text did not have: the pore's two constrictions recovered blind and landing on the GGGVGD filter motif and on Phe2513/Ile2517; the IP₃→gate distance resolved into 103 Å axial vs 120 Å through space; and ***Dictyostelium* iplA carries none of PF08709, PF02815 or PF00520** — a characterised receptor the family's defining signature does not find. |
| S1 | Toolchain install + positive/negative control benchmark (RyR is the sharp decoy) | S0 | completed 2026-08-18 | **Recall 24/25 (96 %); specificity 31/31 (100 %) — but only after S1 had to implement D14.** All 12 binaries resolve with recorded versions (`results/toolchain_manifest.txt`); MAFFT proven invoked (rc=0, 56 seqs → 9,494 cols). **On the first run all six RyR decoys were promoted at 45** — they carry all four diagnostic Pfams, so `pfam+20` also satisfies the D3 gate; specificity was 25/31 and every failure was a RyR. Added the **labelled-bait sister-family test** to `discovery/candidates.py` (D14/D7, margin 0.10, positive test on distances, name never consulted, length band not the call) → all six now capped at 39 with the margin recorded. The margin is 31/31 correct under both identity metrics (true-ITPR +0.065…+0.874 vs RyR −0.868…−0.536). **Caveat: the deepest true member, *Dictyostelium* iplA, has a margin of +0.065 — inside the D7 no-call band**, so profile-based assignment is mandatory for the deep branches. One recall miss: fly `Itpr` at 35, its nearest sibling 0.342 vs the 0.35 breadth threshold (0.420 under `covered_only`). → `results/benchmark_controls/report.md` |
| S2 | Uncapped InterPro enumeration of the family Pfams → census v2, with a positive ITPR/RYR call on every record | S1 | completed 2026-09-03 | **15,421 proteins across 1,488 taxa; ITPR 6,433 / RYR 6,807 / unassigned 2,181 (14.1 %).** All three seeds walked their cursor chains to the end (63 + 67 + 60 pages, 0 restarts). **InterPro's advertised `count` is wrong in both directions** — PF08709 +168, PF01365 +56, PF08454 −176 — while what it *serves* matches UniProt's independent count to within 2 records every time (D20); the first version of this task's own completeness test failed PF08454 on a complete walk because of it. The documented API host answered 1 request in 12 during the run while `/interpro/wwwapi/` answered 12 of 12, so the client now fails over between them. **The call is a positive architecture test audited against a label it never sees (D14b): 6,191 symbol-labelled records, 0 disagreements**, and re-checked on sequence alone over the per-phylum core panel — 27/27 and 26/26 agreement wherever the bait margin decides, with all four sign flips inside D7's no-call band. **A PF08709-only census would have missed 2,914 records, 758 of them called ITPR across 385 taxa** (D21) — including *Dictyostelium* iplA, which the union recovers and correctly leaves `unassigned` at 2/5 signatures. **The plant and fungal records resolve phylogenetically**: every Viridiplantae call is Chlorophyta (Streptophyta 0/15 records, 13 taxa) and every fungal call is in an early-diverging phylum (Dikarya contributes no records at all). Delta vs census v1 ×9.8 with **no enumeration holes** — all 866 absences verdicted. 7 records carry the complete architecture in 1,528–1,993 aa: truncated gene models UniProt does not flag. → `results/census_v2/report.md` |
| S3 | Profile-HMM sweep (itpr.hmm + ryr.hmm) over vertebrate reference proteomes + jackhmmer-to-convergence completeness argument → census v3 | S1, S2 | completed 2026-09-04 | **Census v3: 16,039 records — ITPR 8,000, RYR 7,432, unassigned 605, conflict 2 — with two independent verdicts on every row (D23).** Two profiles built from the S2 census under enforced selection rules (34 ITPR / 22 RyR seeds, 2,684 / 4,930 match states) and swept over **763 vertebrate reference proteomes, 14,414,821 canonical proteins, 6.97 G residues**. **The instrument was calibrated before use**: profile assignment vs S2's architecture call over the whole seeded space agrees **11,875 / 11,876** (seeds excluded), and resolves **2,314 of the 3,361 records the architecture rule could not call**; S2's gene-symbol fallback on 1,219 records is overturned exactly once. **The one genuine disagreement is a correction**: *Tieghemostelium lacteum* A0A152A7I8 (2,845 aa), called RYR by S2 on PF06459 — "Ryanodine Receptor TM 4-6", which is the pore both families share — beaten 303 to 133 bits by `itpr.hmm`. That is D14b's own conditional-specificity caveat firing, once. **The sweep needed a gate the roadmap did not have (D22)**: the first run called **14,981 vertebrate proteins RYR**, 12,874 under a tenth of the profile and named *Tnnc2* / *Rspry1* / *Cabp1* / *Rnf123*, because `ryr.hmm` carries SPRY and `itpr.hmm` has nothing to match it. The floor is **200 match states**, measured from this project's own S0 coordinates (shortest PF08709 200 aa, longest SPRY 137 aa); it costs 86 of 12,020 records their profile verdict and keeps split gene models. **Completeness measured in both directions**: of 2,787 v2 ITPR records from a swept proteome the sweep did not return, **zero were in the database and missed** — all 2,787 are UniProtKB entries the reference proteomes do not contain; and the sweep adds **618 proteins the InterPro census never returned**, all 209–942 aa fragmentary gene models. **jackhmmer, all three seeds** (the *Acanthamoeba* run folded in 2026-09-04 after 7.6 h, against 2.4 h and 2.5 h): `itpr_fly` converged in 10 rounds (D10 clean); `itpr1_human` and `itpr_acanthamoeba` both reached the ceiling (K3). **The family core is seed-independent** — a vertebrate, an insect and an amoebozoan seed recover *identical* family content: ITPR/architecture **1,656 and ITPR/partial 952 in all three models**, RYR/architecture 1,565–1,569. What differs is everything that is not the family: the *Acanthamoeba* model rests on 26,266 targets against human's 7,222, and **18,670 of them are proteins neither profile scores at all** (7 in each of the other two runs). **That drift is invisible to K1** — off-family accretion *dilutes* the sister-family share instead of raising it, and it falls 27 points across this run — so K3 is the rule that caught it (**D10b**). Module-only matches are 10–34 % of these models, which is why they asymptote rather than reach zero. **Three method fixes**: MAFFT is not reproducible at `--thread -1` (**D24**), K1 measured sister-family *level* where only its *rise* is evidence (**D10a**), and K3 punished a run that converged on its last round. → `results/census_v3/report.md` |
| S4 | Genome scope: declared assembly manifest (the denominator) + download tooling | S1 | completed 2026-09-04 | **309 genomes: 161 vertebrate orders ∪ 169 margin species, 552.5 Gbp (≈553 GB FASTA, ≈166 GB as zip, against 1,434 GB free).** Scope confirmed with the user: the *full* margin union, so the bird question is settled per species rather than by sample. **The margin species are derived from the committed census tables, not hand-listed** (unlike the PIEZO port) — four rules with stated thresholds in `s4_manifest_lib.derive_margins`: `zero_hit_proteome` 15, `missing_paralog` (<3 paralogs) 101, `fragment_only` (longest record < `family.MIN_LENGTH_AA` = 2,000 aa) 100, `anchor` 5; 61 genomes carry more than one reason and 11 carry three. All 169 resolved to an assembly — **no margin species is excluded for want of a genome**. **The margin is a bird problem**: 130 of 169 margin species are Aves, against 21 Actinopteri and 10 Mammalia, which is S3's annotation-depth result reappearing as a scope requirement. **D9 is live in the denominator**: 168 RefSeq / 141 GenBank, and **33 genomes carry no gene set at all**, so no annotation claim may be quoted for them. Three assemblies exceed 10 Gbp (*Protopterus annectens* 40.1, *Lissotriton helveticus* 23.2, *Bombina bombina* 10.0) and dominate both download and S5 runtime. 67 reference species have no NCBI order rank (65 percomorphs + 2 sharks) and cannot represent one; the margin rules still reach them. `fetch_genomes.py` tested end to end on 3 small genomes — md5-verified against NCBI's own `md5sum.txt`, `.done` resume, fetch → search → purge, and a failing search correctly keeping its genome. **Two bugs found by testing**: margin species were matched by a last-wins taxid lookup while order reps used the rank function, putting one species in twice under two accessions; and `--verify` skipped any checksummed file the install does not keep, so a deleted `.fna` would have passed — the kept checksums are now rewritten to exactly the installed files and a missing one is a hard failure (proved by deleting a file and watching it fail, exit 1). → `results/genome_manifest_notes.md` |
| S5a | **Genomic sweep: bait panel, calibrated pipeline, pilot** — the instrument for S5, measured rather than ported | S4 | completed 2026-09-04 | **38 baits (30 ITPR + 8 RyR control), 121,294 residues, derived from census v3 by seven enforced rules; pipeline piloted on 6 genomes, 24 cells, 43 loci, 0 control failures.** Three calibrations that a port could not supply. **(1) miniprot's `-G`.** A too-small max-intron does not lose a gene, it *splits* one, and a split ITPR reads out of the ledger as `fragment` — so `-G` is a threshold on the call. Measured across an 11-species Ensembl panel: **no ITPR gene has an intron over the 200 kb default** (widest 152,216 bp, human ITPR1) while the RyR control exceeds it twice (RYR2 227,927 bp). The rule scales from the widest measured intron rather than from the intron-per-Gbp ratio, which is largest in the *smallest* genomes and would have asked for a 7.4 Mbp `-G` on the lungfish. **(2) The unlabelled baits were competing as a fourth paralog.** Rescue attribution ranked `vertebrate_basal` (the gar/chimaera/lamprey seeds — evidence that a region is an ITPR, not a hypothesis about *which*) against the real paralogs: in *Todus* an ITPR1 trace at a true 0.235 margin over ITPR2 was reported at 0.027. Three cells moved from `tblastn_trace_ambiguous` to `tblastn_trace` once fixed. **(3) `ATTRIBUTION_REL_MARGIN` does not transfer (D25).** The PIEZO port's 0.333 comes from a family 40–50 % identical where these are 61–68 %; at complete loci whose paralog identity the annotation independently establishes, **7 of 9 sit below it** (range 0.225–0.347, median 0.281). A complete locus is an upper bound on a rescue fragment, so under 0.333 no absence claim could ever be attributed → **0.22**, with the limitation stated. **D14 was never in doubt at genome scale: at all 43 loci only one family's baits aligned at all** — sharper than the protein level, where S1 had all six RyR decoys promoted at 45. **The RyR control paid immediately**: in *Todus mexicanus* (the `zero_hit_proteome` bird) it found one RyR locus where a bird has three, and unannotated — what is missing there is contiguity, not genes. **D4's contiguity bar, set at the median measured ITPR span (142,212 bp) rather than an invented number: 120/309 manifest genomes fail it — 66 % of birds against 11 % of fish, and 68 % of margin species against 12 % of order reps.** The margin set was chosen for what its *proteomes* lack and turns out to be very largely the same genomes whose assemblies cannot hold the gene; the two signals are confounded. **And the paralog recovered from a fragmented assembly is systematically the short one** — below the bar ITPR3 (82 kb median span) is found 2/3 while ITPR1 (186 kb) and ITPR2 (231 kb) are 0/3, so S0's 6.5× span asymmetry returns as a detection bias pointing the same way as the loss signal *(pilot of 6 — direction and mechanism, pending: S5b at 309)*. Screen validated by its own negative controls on every build. → `results/genome_ledger/report.md` |
| S5b | **The full 309-genome sweep** → per-genome ledger (found / lost / assembly-gap) + novel gene models → census v4 | S5a | completed 2026-09-05 | **309/309 genomes swept, 0 failures, 1,236 cells, 2,144 loci; census v4 = 17,097 records.** **The RyR positive control fired in every one of the 309**, so no genome is excluded on control grounds. **D14 held absolutely: across all 2,144 loci only one family's baits aligned at all** — the ITPR and RyR panels never once contested a locus. That is a sharper separation than the protein level affords (S1 had all six RyR decoys promoted at 45 before the sister test existed) and it is a result in its own right: the family ambiguity is a property of searching protein fragments, not of the two families. **Only 4 cells are `absent`, and all four are cyclostome**: *Petromyzon marinus* and *Myxine glutinosa* each carry ITPR1 and lack ITPR2/ITPR3. Both are above the contiguity bar with a firing control, so the call survives both of D4's gates — but the panel has **no labelled ITPR2 or ITPR3 bait for cyclostomes at all** (4 of its 6 unfilled slots), so this is one findable ITPR in each, consistent with the vertebrate trio arising at 2R after the cyclostome split, and S7 is what turns it into a statement about 2R. **D4's bar is the headline caveat: 120/309 genomes cannot carry the gene on one contig — 66 % of Aves against 11 % of Actinopteri, 68 % of margin species against 12 % of order reps.** Recovery is **98–99 % above the bar and 57–70 % below it**. **S5a's span-bias hypothesis is confirmed at scale** (it wobbled 15 → 8 → 12.5 points during the sweep and settled at **13**): below the bar ITPR3 70 %, ITPR1 61 %, ITPR2 57 %, and ITPR3 is the shortest gene — so the detection bias runs the same direction as the loss signal the margin species were chosen for. **The attribution margin is now validated on the evidence it acts on**: 53 rescue fragments whose paralog an assembly's own annotation names agree with the attribution **53/53**, over margins 0.227–0.445. **0 fall below the 0.22 in force; 31 of 53 would have fallen below the inherited 0.333** and been reported ambiguous. S5a's recalibration was right and the floor sits 0.007 below the observed minimum. **Annotation quality differs by paralog, and the control matters**: uncontrolled the gap looked like ITPR3 54 % vs ITPR2 14 %, which was assembly quality; holding contiguity constant across 487 loci it is **ITPR3 88 % correctly named vs ITPR1 65 %** — a 23-point gap between genes of near-identical protein length in the same genomes. **Census v4 adds 1,058 gene models from 224 genomes**: **318 ITPR models exist only as DNA** and **167 more sit inside an annotated gene carrying no family name**, unreachable by any name-based search. 3 loci are claimed by a cell other than the one their annotation names (1 with a coherent sibling locus) — supplied to S18, not adjudicated. **Two instrument fixes found by the run.** `MIN_LOCUS_IDENTITY` (0.40) — in the two >20 Gbp genomes a 2 Mbp `-G` chained shared-module hits at 23–34 % identity into apparent loci, inflating *Lissotriton* ITPR3 to 10 loci; measured, the 571 annotation-confirmed loci have a **minimum identity of 0.759** and the junk sits below 0.35, so any floor in 0.35–0.50 drops the same 22 and zero real loci. Status was never wrong (the best locus always won) but copy number was. And **miniprot output is now written atomically**: the driver reuses any non-empty GFF, so an interrupted run left a truncated file the next run accepted as complete — undetectable afterwards. Audited: 0 genomes affected. **The budget was wrong in S5a and is corrected**: the sweep is download-bound (~5.5 MB/s aggregate, not scaling with shard count), 164 GB over ~19 h wall against ~2.5 h of compute; S5a projected 3.4 h from compute alone. miniprot is *faster* than budgeted (15 s/Gbp, not 36). Both giants completed on the chunked path (25 and 14 chunks, 1.8 h and 2.5 h). → `results/genome_ledger/report.md` |
| S20a | **Non-vertebrate sweep** — the family's true range across eukaryotic reference proteomes, and whether the land-plant / dikarya absence is a genome fact or a database fact | S3 | completed 2026-09-05 | **6,928 reference proteomes, 63,144,898 proteins, 24.93 G residues swept with S3's own two profiles — 662 carry an ITPR call.** The four eukaryote groups partition Eukaryota with S3's `vertebrata`, so the two sweeps' denominators add; archaea taken whole (634), bacteria genus-stratified (3,537 of 17,981) under a rule chosen to be generous to the hypothesis it tests — one proteome per genus, the *largest*, because a bigger proteome is a more sensitive place to find a homolog. 23 listed proteomes are unpublished in the release FTP tree and are recorded as exclusions rather than aborting the fetch. **One search, two sensitivities**: every `hmmsearch` ran at `-E 10` and the primary E ≤ 1e-5 call is a filter on the same domtblout, so the relaxed set is a superset of the strict one from the same search rather than a second experiment. **Both of S2's headline absences are confirmed with the seeding filter removed: Streptophyta 0/384 (16.3 M proteins) against 15/48 Chlorophyta, and Dikarya 0/1,353 against 28/1,527 fungal proteomes.** Archaea 0/634, bacteria 0/3,537. **But the fungal losses are patchy, not basal** — Glomeromycota 0/27, Mortierellomycota 0/18, Kickxellomycota 0/35 and Microsporidia 0/29 are empty too, so this is repeated independent loss rather than one event. **Every negative carries a positive control inside the same search** (brief step 5): at E ≤ 10 with the four family Pfam models, PF08709 — the IP₃-binding core that names the family — returns **0** substantial matches in land plants against **26** in Chlorophyta, and **0** in Dikarya against **16** in Mucoromycota, while PF02815 (MIR, which every eukaryote carries on other proteins) returns **633** and **4,376** in those same genomes. The instrument demonstrably works there and finds everything except the receptor. **All 99 plant and fungal ITPR records were chased individually** (brief step 3), pooled from both instruments so the UniProtKB entries outside reference proteomes are not quietly dropped: **47 `real_gene`, 52 `fragment`, 0 contaminants, 0 without genome backing**. The contamination test had full power — every record has a cross-kingdom blastp hit, at **19.9–45.8 % identity, median 24.1 %**, nowhere near the 95 % call. **D14 outside the vertebrates: 0 of 704** targets scored by both profiles above the bit floor fall inside the no-call band. **And a trap caught before it became a result**: the sweep reports 142 RYR calls in fungi, algae and protists, all 500–2,000 aa proteins spanning **4–9 %** of a 4,930-state model — read as gene counts they would invent ryanodine receptors across half the eukaryotic tree, so every call is now reported with its model coverage (63 % of ITPR calls are architecture-level against 1 % of RYR ones). The two that *are* architecture-level are genuine and matter: *Salpingoeca rosetta* (choanoflagellate, 5,340 aa, not a seed) and *Capsaspora owczarzaki* (filasterean, 6,625 aa, a `ryr.hmm` seed and flagged circular) carry both families at full length, so the ITPR/RyR duplication predates Metazoa. **jackhmmer (brief step 4) for the two groups the negative claims rest on. viridiplantae: converged in 5 rounds, D10 clean, zero sister content in every round — and all 70 targets in the converged model are Chlorophyta, including all 6 that only iteration found**, so the land-plant absence is not a sensitivity artefact. **Fungi ran the ceiling and K2 fired at round 3**: the included set grew **34.3×** (42 → 1,442, reaching 6,302 by round 10) while the family content went 33 → 35 — the growth is entirely off-family, so **K1 read 0.000 in every round**. That is S3's **D10b** in a new setting, caught this time by K2 rather than by the round ceiling. The accepted model (rounds 1–2) reaches Ascomycota exactly twice and both records are dolichyl-phosphate-mannose mannosyltransferases — the MIR-domain sharer S1's decoy panel was built around, not a receptor, so the Dikarya absence survives its only loophole. Iteration-only targets are **named** in `jackhmmer_iteration_only_s20.tsv` rather than counted, and the report calls them known decoys only when their own protein names say so. **Census v5: 17,882 records** (+785 from this sweep), 8,807 ITPR across 1,401 taxa, lineage on every row including S5b's 488 genomic models joined to S4's manifest. **No v4 call is overturned and S3's 2 conflicts stand** — the first version of the merge resolved them, until it was noticed that S20 scores with the same profile pair that produced one side of each disagreement, so it is the same instrument in a different database rather than a third opinion. **Three silent bugs fixed**: the presence table was overwritten with a fraction of its own denominator by per-group runs; hits were attributed by taxid where 48 taxids carry two reference proteomes each (now measured by a header pass over each proteome file); and **UniProt's `tax_id:` search is hierarchical**, so a 100-term batch matched more than 100 records, the page capped, and requested taxids came back as blank rows — *Arabidopsis* and *Chlamydomonas* among them — now the exact `taxonIds/` endpoint with the batch below its measured 25-row page cap. **The one-search-two-sensitivities design was tested, not assumed** (`s20_test_sensitivity.py`), and the test disproved the assumption it was written to confirm — twice. `--domE` acts on the *conditional* E-value, normalised by how many sequences passed, so a looser `-E` inflates every c-Evalue (1.41× measured) and drops marginal domains; and HMMER prints the sequence E-value to two significant figures, so 42 of 13,770 plant targets whose true E sits just above 1e-5 print as `1e-05` and a `<=` filter admits them. **Neither reaches a result**: over 798 + 13,728 shared targets every assignment and every D22 gate decision agrees, all 42 boundary targets are gate-declined, and the one domain-row difference that moves anything moves an evidence class on a single record (the task uses the relaxed side, which excludes −2.0-bit domains from coverage). The test also refuses to **pass vacuously** — its first run went green on `archaea`, comparing two empty sets. → `results/s20_sweep/report.md` |
| S20b | **The remaining per-group convergence runs** — jackhmmer to convergence over `protista_other` and `metazoa_nonvert` under D10 | S20a | completed 2026-09-05 | Split out because the per-round cost is the **alignment, not the search**, and does not fall with more cores (~45–65 min/round at any thread count here); the two runs took 8.2 h and 8.4 h. **Both hit K3, and the two K3s mean opposite things.** `protista_other` drifted: it grew 721 → 26,148 targets but **never by more than 3.18× in one round so K2 never fired**, and its sister share *fell* (0.28 % → 0.16 %, `sister_rise` negative in 6 of 10 rounds) so K1 moved the wrong way — the sharpest **D10b** case in the project. Its finished model is **772/22,913 (3 %) family**, its largest single contribution **11,344 Apicomplexa proteins** in a clade this sweep called the family in 0/60 proteomes. `metazoa_nonvert` did not drift at all: growth never above **1.21×**, `sister_rise` peaked at **+0.042** against K1's 0.10, and its model is **2,198/2,908 (76 %) records the profiles call ITPR or RYR**. It had simply not finished, at ~35–70 new targets a round. **K3 returns the same verdict for both because it counts rounds rather than content** — right for a ceiling, and the reason the report tabulates the composition beside the verdict. Two instrument corrections fell out: `s3_kill` reports a K3 run's pre-ceiling rounds as `accepted` (correct for K1/K2, where the rounds before the drift are usable), so composition rows now carry the verdict and a `supports_completeness` flag and a disowned run is reported separately rather than mixed into a completeness argument — without which 22,913 accreted targets, and a 3 MB table naming them, would have read as one. **The completeness result itself is clean and does not depend on either verdict: in all four groups the iterated model settles on *exactly* the single pass's ITPR count** — viridiplantae 54 (round 2), protista 729 (round 3), fungi 35 (round 4), metazoa 1,194 (round 5) — so iteration finds no receptor the single pass missed anywhere in the non-vertebrate tree. What the extra rounds bought was 24,153 more non-family targets in the protists and 870 in the fungi. → `results/s20_sweep/report.md` |
| S23a | **The non-vertebrate genomic sweep's instrument** — declared genome manifest (the denominator), gene-span calibration outside the vertebrates, copy-number bait panel, and a pilot | S20a | completed 2026-09-05 | **194 genomes, 100.4 Gbp declared; 80 baits (37 ITPR + 16 RyR + 27 MIR control); 14 anchors piloted — and *four* thresholds plus one selection rule ported from S5 turned out not to transfer.** **The headline is D26: S5's positive control does not exist here.** S5b can write "the RyR control fired in all 309" because every vertebrate has three RyRs; S20 found architecture-level RyR in **2 of 6,928** non-metazoan proteomes, so a land plant's silent RyR is the *correct* answer and witnesses nothing. Carried over unchanged it would have made every negative claim in this task unfalsifiable while looking controlled. The replacement is the **MIR-domain sharer** (PF02815 — the mannosyltransferases S1's decoy panel was built from and S20 used as its in-search control, where PF08709 returns 0 in land plants and PF02815 returns 633), drawn **per control clade** because a chlorophyte enzyme is not a control for *Arabidopsis*. 28 control clades, **22 strong / 5 weak**, graded rather than counted — and the grading paid immediately. **The denominator** is five rules under one principle — *sample most finely where the negative claim is* — all derived from S20's presence table except the anchors: phylum_rep 66, class_rep 87, absence_clade 35, anchor 14, copy_number 36, from 24,596 NCBI eukaryote reference assemblies less 6,216 vertebrate ones. Deriving G3 rather than restating S20a's list found **three absences that summary never named — Bacillariophyta 0/16, Rhodophyta 0/12 and Cestoda 0/11, a metazoan clade**. NCBI files reference assemblies under *strain* taxids, which silently cost four of fourteen anchors — the apicomplexan, microsporidian, chytrid and amoebozoan the claims are named after — until resolution went through the ancestor lineage. **The thresholds.** Gene spans measured from NCBI annotations (not from miniprot, whose `-G` shapes the loci it reports) over 16 genes in 14 bands: **3,739–324,840 bp, median 19,435**, a ~100× range against 6.5× inside the vertebrates, so the contiguity bar is **per group** — metazoa 83 kb, protists and fungi 7–9 kb. S5's single 142,212 bp bar failed 120/309 genomes; here **21 of 194** fall below theirs, and the genomes carrying the negative claims are judged against ~9 kb. The bait length band and the architecture-exception floor both had to be **stratified per band**: measured over the whole population (53 % Arthropoda) they are 2,550–3,100 aa and ~1,095 bits and they discard *Bodo saltans*, a genuine deep euglenozoan receptor; stratified they are **2,450–3,250 aa and 752 bits**. **The pilot found two more, which is what a pilot is for.** (1) `MIN_LOCUS_IDENTITY` = 0.40 does not transfer: S5b measured a wide empty gap in the vertebrates (confirmed loci ≥ 0.759, junk 0.23–0.34), but in *Chlamydomonas reinhardtii* — which carries a documented complete 5/5 receptor — miniprot returned **9 ITPR-family alignments at 24.6–28.4 % identity and the floor discarded all nine**. (2) S5's R3 spread rule was never ported, and porting it as a *filter* made things worse (85 → 78 baits, Chlorophyta down to one) because the shortlist was score-ranked and Chlorophyta's nine records are led by seven *Cymbomonas* ones: **a spread rule has to shape the shortlist, not filter its output.** Fixed, the panel carries both *Cymbomonas* and *Chlamydomonas*. Also fixed: a reused S3 seed was consuming its band's quota instead of adding to it. **The screen earned its place again** — it rejected a control bait UniProt calls a "MIR domain-containing protein" from the tapeworm *Rodentolepis nana* that `ryr.hmm` scores at **1,121 bits against 179**. **Pilot: 14 anchor genomes, 0 failures.** All six positive controls recovered and matched to their own annotation — *Drosophila* Itpr, *C. elegans* itr-1, *Nematostella*, *Strongylocentrotus*, *Dictyostelium* iplA and *Chlamydomonas* — and 8 `no_locus`. **The Dikarya and land-plant absences are now controlled genome facts**: *S. cerevisiae* (4 MIR loci), *Neurospora* (2), *Arabidopsis* (1), *Oryza* (2), *Physcomitrium* (2) and a microsporidian (1) each return the control and no receptor. **One genome is `uncontrolled` and it is the one the weak-control flag predicted**: *Toxoplasma gondii* returns neither a receptor nor a MIR locus, so "Apicomplexa 0/60" cannot yet go to assembly level — without the control this sweep would have reported that absence as a finding. D14 held at every locus: no MIR or RyR locus was ever called ITPR, and *Drosophila*'s ITPR locus has a family margin of 1.0 (the RyR baits put no alignment there at all). → `results/s23_scope/report.md` |
| S23b | **The blocking measurements, and the instrument they change** — the three items S23a's pilot left in front of the full sweep | S23a | completed 2026-09-06 | **All three transferred badly, and one of them was blocking a headline claim.** **(1) The control.** S23a's `uncontrolled` *Toxoplasma* was not a PF02815 problem, it was **fixing one profile in advance for every clade**. Six candidate profiles — all large, deeply conserved, multi-exon eukaryotic families — are now measured over each clade's *own* swept reference proteomes and each clade takes the one that is there. **Apicomplexa takes Myosin_head, present in 36/36 Aconoidasida and 23/23 Conoidasida proteomes at ~1,450 aa**, against PF02815's one protein per class — and *Toxoplasma* goes from `uncontrolled` to **controlled across a kingdom boundary**, which is exactly what "Apicomplexa 0/60" was blocked on. **Rhodophyta vindicates the design independently**: red algae carry myosin in **8 %** of their proteomes, so the fixed-profile approach would have given them a control found in 1 of 12 — measured, they take SMC_N at 12/12. The MIR bait stays in the panel whatever the measurement says, because dropping it would buy a proof-of-search and sell D14's negative control. 27 clades take Myosin_head, 1 takes SMC_N; 56 control baits, 48 strong / 8 weak. A **cross-kingdom tier** came free from the existing panel — every clade's control bait is already searched in every genome, so "a bait from another kingdom aligns across this locus" is a stronger statement than "the assembly is readable" and costs nothing to record. **(2) `MIN_LOCUS_IDENTITY`.** The sweep now **records** every cluster to 0.15 and applies the call floor downstream (D27), because a floor measured from the population it has already filtered is circular — S5b could only measure 0.40 because it had the loci 0.40 excluded. `s23_calibrate_loci.py` measures it against loci whose identity the assembly's **own annotation** establishes, and **refuses to write below 50 genomes / 15 confirmed loci**: a smoke-test over one genome produced 0.25 from two loci, wrote it, and the next sweep read it back and moved three genomes out of `no_locus` on the strength of it. **(3) Span against CDS footprint.** *Drosophila*'s 22 kb *Itpr* sits in a **297,487 bp cluster around an 8,514 bp CDS footprint — 35×**. Harmless for a status call, fatal for a copy count, so copy number is counted on **non-overlapping complete alignments, not on clusters** (D28); the *Drosophila* copy is bounded at 10.8 kb. **Two bugs, both latent.** `s5_classify`'s family-name list had drifted from `family.py` and was missing `itr-1`, so the pilot read *C. elegans*'s correctly-recovered receptor as an annotation naming something else. And the rescue HSP filter read `h["start"]`/`h["end"]`, which `parse_tblastn` does not produce — it never fired because every `no_locus` genome up to *Salpingoeca rosetta* returned zero HSPs, so the generator's predicate was never evaluated; replaced with S5's own `filter_hsps_outside`. Panel rebuilt 80 → **108 baits** (37 ITPR + 16 RyR + 55 control), 277,611 residues. → `results/s23_baits/control_profile_choice.tsv`, `results/s23_scope/report.md` |
| S23c | **The full 194-genome sweep** → per-genome copy-number ledger + the plant/fungal/apicomplexan absence claims at assembly level → census v6 | S23b | completed 2026-09-06 | **194/194 genomes, 100 Gbp, 0 failures. All 35 absence clades hold at assembly level, and all 35 are controlled** — Ascomycota 0/31, Streptophyta 0/25, Basidiomycota 0/17, Magnoliopsida 0/3, **Apicomplexa 0/3**, Microsporidia 0/2, and 29 more. **Zero genomes are `uncontrolled`**: 115 `controlled_cross_kingdom`, 77 `controlled_by_target`, 1 `controlled_partial`, 1 `no_control_bait` — against 1 of 14 uncontrolled in the S23a pilot. **The headline methodological result is a negative one, and it changed the instrument.** The brief asked to re-measure `MIN_LOCUS_IDENTITY`. Two things came out. (1) **The annotation axis barely exists here**: of 917 recorded clusters only **21** sit on a gene whose name says anything (10 name the family, 11 name something else), because outside the vertebrates gene models carry locus tags — *Chlamydomonas* files its receptor as `CHLRE_16g665450v5`. So a second axis was added: every recorded cluster scored against `itpr.hmm`/`ryr.hmm` (D23), re-derived from the archived GFFs. (2) **No threshold on identity separates the two populations** — confirmed loci reach down to 0.193, contradicted ones up to 0.318 — and coverage is no better (Youden J 0.71 vs 0.70, measured, after a draft claim that coverage separates was contradicted by the data and removed). **S5b's inherited 0.40 would have discarded 87 confirmed loci, 56 of them complete gene models — a third of everything the sweep found.** Outside the vertebrates a locus's identity to its nearest bait measures how far away the nearest bait is, not whether it is a gene. Identity is therefore **retired as a call gate** and the profile call carries it (D14/D23), validated against the one axis it does not share: **10/10** on annotation-confirmed loci, declining **10 of 11** contradicted ones. The single exception is *Emiliania huxleyi* `IPR1` at 1,008 bits over 95 % of its bait — almost certainly an IP₃-receptor name the family list does not carry, reported as a conflict rather than reconciled (adding `ipr1` as a substring would collide with every InterPro accession). **Copy number is the deliverable and it is wide**: 0 in 116 genomes, 1 in 43, 2–6 in 32, and then *Dysidea avara* 8, ***Stentor coeruleus* 13** and ***Macrostomum lignano* 18** — which answers S20a's open question directly, since *Macrostomum*'s 62 records resolve to 18 real genes. ***Cymbomonas* 3 copies from 2 clusters**, the other S20a question. D28's rule earned its place: it recovered **7 genes in 7 genomes** (174 copies against 167 clusters) that cluster-counting merges. **D14 held everywhere**: 189 of 195 graded loci have a family margin of 1.0, minimum 0.597, and **0 control loci were ever called ITPR**. Census v6 = 18,065 rows (+183 genomic models; 100 of them `annotated_unnamed`, the dominant class outside the vertebrates). Three bugs fixed, two pre-existing in shared code: `s3_assign` crashed on a hit scoring exactly 0.0 bits with no counterpart (two separate `None` subscripts, never fired until the profiles ran over genomic models); the rescue HSP filter read keys `parse_tblastn` does not produce; and the control now requires a **complete** recovery, symmetric with the family call, admissibility stated as a set rather than a `startswith("controlled")` prefix test that would have admitted `controlled_partial` silently. → `results/s23_scope/report.md` |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl on census representatives per clade × kingdom, RyR outgroup included | S2, S3, S5b, S20a | completed 2026-09-06 | **134 representatives from census v6's 18,065 records by eight coded rules (D8), MAFFT L-INS-i → 11,796 columns (76.3 % gaps, 64.2 min, `--thread 1` per D24), trimAl `-automated1` → 1,790 columns kept (15.2 %, 7.83 % gaps), of which 96.6 % are parsimony-informative.** Median tip coverage 0.96; 1 tip of 134 under half. **Two rules the data forced before the alignment was built.** (1) *A paralog label is vertebrate-only*: ITPR1/2/3 are a 2R product, and four non-vertebrate records carry a type number by annotation transfer (*Tetrabaena*, *Acanthamoeba*, *Branchiostoma*, *Phallusia*) — kept in the audit, grouped by taxonomic grade. (2) *A bait attribution is not a paralog label where it is constant*: the cyclostome band holds 19 ITPR records over 12 loci, 6 carry a label, **every one `ITPR1` and every one from the S5 bait attribution rather than an annotation** — because lamprey, hagfish and *Myxine* each carry **three** full-length loci and the ITPR1 bait wins all of them. R3 now takes the complete copy set of two species per band unlabelled. **Results.** D14's separation **confirmed**: within-paralog identity 0.908, ITPR-to-RyR 0.259, separation **0.650** against S1's 0.579 (compared on the separation, not the absolutes — the two estimators differ). **Sister preview: ITPR1 × ITPR2 leads at 0.788 with an interquartile range (0.773–0.806) that does not overlap either other pair** (ITPR2×ITPR3 0.746, ITPR1×ITPR3 0.736) — a hypothesis for S7's AU test, not an answer. **The cyclostome trio: all 6 loci fall nearest ITPR1, and the control says the signal is real** — the invertebrate, protist and RyR groups lean the same way at a median margin of 0.007, the cyclostome loci at 0.043, 6× that. **One bug D11 caught**: the domain track mapped protein residues by counting ungapped positions in the *trimmed* row, which renumbers every residue after the first discarded column — it put the domains at ~2/3 of their true position and dropped the pore off the end. Fixed by carrying trimAl's own `-colnumbering` output as `column_map.tsv`; the pore now lands on the alignment's highest conservation peak. **One stated deviation**: the brief's "every novel model from S5/S23" is 424 sequences and would fill the tree with bird and teleost gene models; R7 takes one per (clade band, paralog cell) and one per non-vertebrate group instead. → `results/msa_v2/report.md` |
| S7 | ML phylogeny: IQ-TREE 2 + ModelFinder + 1000 UFBoot, rooted on RyR; resolve which two of ITPR1/2/3 are sisters (AU test on the three rooted topologies) | S6, S23c | completed 2026-09-07 | **The sister question has an answer: ITPR2 + ITPR3, with ITPR1 outside.** Unconstrained ML groups them at SH-aLRT 100 / UFBoot 100; the AU test over the three constrained topologies (10,000 RELL replicates) rejects **ITPR1+ITPR2 (p-AU 1.8e-05)** and **ITPR1+ITPR3 (1.65e-05)** and does not reject ITPR2+ITPR3 (0.476) or the ML tree (0.525) — both survivors carry the same pair. **This contradicts S6's identity preview**, which ranked ITPR1x ITPR2 highest, and the review's §7.4 audit found no published support-annotated ML analysis with an RyR outgroup that fixes the pair. Tree: 134 tips x 1,797 trimmed columns, `Q.insect+R7` chosen by a measured two-stage scan (`-m MFP` projected ~20.9 h and was abandoned; both stages committed), logL -215,452.0, threads pinned (D24). **69.5 %** of 131 internal nodes clear SH-aLRT >= 80 *and* UFBoot >= 95. `--bnni` (169 min): **no claim weakened or lost**, 9 of 10 clear both thresholds again; the tenth, the bare ITPR1 core, was already unsupported (47.8/95 -> 47.5/73). RBH on all 5 names the tree declines to place: **5 of 5 upheld**, so those are tree uncertainty, not annotation error. Cyclostomes: all 6 loci sit in **two well-supported cyclostome-only clades**, each holding hagfish *and* lamprey — neither S6's expansion nor 1:1 orthology, but duplications older than the hagfish/lamprey split. 3R pairs are broken **only by other tips of the same paralog**, which is what a duplication older than the species looks like (D33). **Three method defects found and fixed, two of them invalidating a first AU result** (D32 constraint taxon set, D34 best-tree reporting, D35 length ruler + diversity spread) — the last forced a full S6+S7 re-run, which changed one tip (*Volvox carteri* replaces *Tetrabaena socialis*) and **left the sister answer unchanged**. -> `results/phylogeny/report.md` |
| S8 | Synteny: shared flanking-gene analysis across the ITPR loci | S2, S5b | completed 2026-09-07 | **Neighbourhood confirms the paralog cells, and the surviving 2R links run through ITPR1.** 2,144 loci across 309 genomes (from the per-genome `summary.json`, not the ledger, so a teleost cell's *a*/*b* copies and the lamprey's three are separate rows); 274 genomes carry a gene table. Every real pair is scored against **matched random-window control pairs in the same two genomes** (3 per genome, seeded), so annotation depth, naming convention, window and key rule are all held constant. **Within-paralog Jaccard is 216x-413x its own matched null** (ITPR1 0.211, ITPR2 0.219, ITPR3 0.125) with 98-99.8 % of individual pairs beating their own control; **every cross-paralog and cross-family class sits at or below the null** (max mean J 0.0002 over 168,241 ITPR x RyR pairs) - D14 confirmed by an instrument that shares nothing with the sequence evidence. **The 2R paralogon is visible only through gene-family root keys**, and exactly two links survive: **BHLHE40/BHLHE41 (ITPR1-ITPR2, 62 %/85 % of species, 84x background)** and **GRM7/GRM4 (ITPR1-ITPR3, 42 %/53 %, 93x)**; **ITPR2 and ITPR3 share none at any bar from 10 % to 50 %** - so the pair S7's AU test makes sisters is the one pair whose neighbourhoods retain nothing (a deletion record, not a phylogeny; reported as a tension, not an overturn). **ITPR3's neighbourhood is the one that does not travel**: cross-class J 0.062 against 0.160/0.158, a 2.5x gap, while within a class the three span 0.254-0.347 - in human that ground is the MHC region at 6p21. A flank-consensus **paralog caller** (threshold measured, not typed - 0.4 maximises call rate minus random-window false-call rate) is **405/405 correct on 503 annotation-confirmed loci at a 0.008 false-call rate on 726 random windows**, and adds **131 paralog assignments no random neighbourhood could have produced** (87 on loci with no paralog-named annotation, 36 on extra copies in multi-copy cells, 8 in assembly-gap cells). **The question S7 handed to S8 it cannot answer**: all 6 cyclostome loci carry 20 informative flanks but reach an overlap of at most 2 against the gnathostome consensus, and random windows reach 2 - *underpowered*, not negative. 11 constructed negative controls run on every build (one of them caught a hash-order nondeterminism in `flank_consensus.tsv`); every table's SHA-256 is recorded. -> `results/synteny/report.md` |
| S9a | **The codon alignment every selection test stands on** — a validated CDS for every vertebrate tip, PAL2NAL cross-checked against an independent mapping, codon-aware trimming, and the tree-derived selection sets | S6, S7 | completed 2026-09-07 | **57 / 57 vertebrate tips carry a CDS that provably encodes the exact protein S6 aligned** (43 UniProt-route, 14 miniprot locus realignments), 3,253 codons, trimAl keeps **2,459 (75.6 %)**. **Two of the 57 needed a route change to get there, and both are human**: UniProt lists 5 Ensembl transcripts for ITPR1 and 2 for ITPR2 and in each the *first* is not the aligned isoform, so a route returning the first CDS that downloads substitutes a different isoform for the protein the tree was built on — the routes are now candidate generators and the first that **validates** wins (**D36**); human ITPR1 took 9 candidates. Every translation/protein disagreement is masked to `NNN` (24 codons total, 15 of them internal stops), so pal2nal is never handed a pair that disagrees. Selection sets are **S7's extended paralog clades** re-derived from `rooted.nwk` and cross-checked against `paralog_clades.tsv` — ITPR1 19 / ITPR2 13 / ITPR3 19, which nests **7 unlabelled `vertebrate_basal` tips** inside paralog clades (D30 read forwards) and leaves **6 in no clade at all: exactly the six cyclostome loci S7 handed to S8 and S8 reported underpowered**. 14 constructed negative controls on every build; **two found real defects** — the per-paralog subset files were written from a `set`, so their row order was hash-seed dependent and a rebuild changed a phylip file underneath a running codeml job (fixed; re-verified across three `PYTHONHASHSEED` values), and S1's toolchain manifest never probed codeml / yn00 / pal2nal.pl / hyphy, so S9 opened by reporting a tool 'missing' that had been installed all along (now probed; PAML 4.10.10, pal2nal v14, HyPhy 2.5.101). → `results/selection/` |
| S9b | **The model-based selection tests** — codeml one-ratio, pairwise, two-ratio, branch-site model A on each paralog stem, M1a/M2a + M7/M8 site models, and HyPhy RELAX, with every LRT BH-corrected | S9a | completed 2026-09-08 | **40/40 codeml jobs (22.7 h) + 3 RELAX runs; 12 LRTs, BH-corrected. Three independent framings agree that ITPR1 is held roughly twice as tightly as the other two.** One-ratio ω: **ITPR1 0.0238, ITPR3 0.0415, ITPR2 0.0430** (curated-CDS subsets 0.0212 / 0.0371 / 0.0318, so the genome gene models are not driving them). Whole-tree two-ratio, every one significant: ITPR1 foreground **0.0241** vs background 0.0432 (q 2.9e-63), ITPR2 0.0435 vs 0.0317, ITPR3 0.0455 vs 0.0308. RELAX: **ITPR1 k = 9.36 *intensified*** (LRT 11,317), ITPR2 k = 0.908 and ITPR3 k = 0.836 *relaxed* — the same ordering from a statistic that compares whole distributions rather than point estimates. **Three results needed the fitted parameter, not the p-value, before they could be stated.** (1) M8 beats M7 in all three at q ≈ 2e-4, but its extra class sits at **ω = 1.00000**, codeml's boundary, on 0.3–0.7 % of sites with 0/0/1 sites at BEB ≥ 0.95 — a small class of *unconstrained* sites, not positive selection; M2a vs M1a is 2ΔlnL = 0.00 in all three, its positive class carrying a proportion of exactly zero. (2) Branch-site model A is significant on all three stems, but ω₂ is pinned at codeml's **999 bound on ITPR2 and ITPR3** with the likelihood flat above it (ITPR2's restarts give 162 and 999 at the *same* lnL), so only **ITPR1's stem — ω₂ = 5.53 on 11.0 % of sites, stable across restarts, 8 sites at BEB ≥ 0.95 and 5 at ≥ 0.99 — is reported as a result**. (3) **3 of 12 branch-site restarts converged below their own nested null**, one per stem, and at a *different* initial ω each time (ITPR1 ω₀ = 4, ITPR2 0.5, ITPR3 1.5) — no single starting value would have been safe, which is D37's case measured rather than argued. Also **synonymous saturation within a single paralog**: 88.8 % of all within-paralog pairs past dS = 1.5, dN plateauing at 0.1 while pairwise dS runs past 50, so the pairwise matrix is a diagnostic and every ω is tree-based. One silent failure fixed: HyPhy writes non-finite branch estimates as the bare token `inf`, which is not legal JSON — the ITPR1 run succeeded, the parse threw, and the largest effect in the analysis came back as a blank row. → `results/selection/report.md` |
| S10 | Annotation-bug molecular validation: assembly-version audit, RNA-seq spanning evidence, miniprot gene models, for the two worst cases S5b/S18 surface | S5b | pending | |
| S11 | Structures: AFDB coverage, TM-align vs the cryo-EM IP3R and RyR references, Foldseek AFDB-wide sweep | S2, S6 | pending | |
| S12 | Expression evidence (SRA junction-spanning reads + atlases) for whichever paralog or lineage the census leaves in doubt | S5b | pending | |
| S13 | Gene-tree/species-tree reconciliation — date the duplications that made ITPR1/2/3 | S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3, S5b, S7, S8, S9b, S10, S11, S20, S23c, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Loss dynamics** — ancestral-state reconstruction of per-paralog presence/absence, ORF-integrity screen from miniprot frameshift/stop counts, dating any pseudogene fossils | S5b, S13 | high | pending | |
| S16 | **Duplication history** — are ITPR1/2/3 2R ohnologs; are teleost itpr1a/itpr1b from 3R; copy-number landscape and retention asymmetry | S7, S8, S13 | high | pending | |
| S17 | **Constraint & function** — per-site conservation mapped onto the cryo-EM channel; do the SCA15/SCA29/Gillespie, anhidrosis and neuropathy variants sit in the constrained core? Is the IP3-binding core more constrained than the pore? | S6, S9b, S11 | high | pending | |
| S18 | **Annotation-quality audit** — how often a real ITPR locus is missing, fragmentary, split, unnamed or filed under the wrong paralog (or as a RyR) across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5b, S15 | high | pending | |
| S19 | **Methods results** — per-method contribution ("what would proteome-only searching have missed?"), assembly contiguity as a confounder of loss claims, bait-panel design sensitivity | S5b, S15, S18 | medium | pending | |
| S21 | **Gene architecture** — exon/intron structure across the genome scope from the miniprot CDS blocks; are database "fragments" real exon boundaries or annotation failures; is the ~58-exon architecture conserved | S5b, S18 | medium | pending | |
| S22 | **Ligand-site evolution** — the IP3-binding core is the one part RyR does not share functionally. Is it under different constraint from the pore, does it differ between paralogs, and does it change in lineages that lost the upstream PLC/IP3 pathway | S9b, S17 | medium | pending | |
| S14c | **Manuscript rewrite pass** — one full pass over the draft once every analysis has landed, with the figure audit's lessons applied | S14a, S24 | medium | pending | |

---

## Emergent tasks & new aims

Anything discovered mid-session that deserves its own work goes here rather
than expanding the task in progress.

| Added | From | Task | Status |
|-------|------|------|--------|
| 2026-09-08 | S9b | **The ITPR1 stem is the one branch-site result that stands, and it wants following up.** ω₂ = 5.53 on 11.0 % of sites, stable across three restarts, 8 sites at BEB ≥ 0.95 and 5 at ≥ 0.99 — while the rest of the tree sits at ω ≈ 0.03. The sites themselves are in `results/selection/codeml/bs_alt_ITPR1_w0.5/mlc` (BEB block) on the ITPR1 alignment's own coordinates. S17 maps per-site conservation onto the cryo-EM channel; those eight positions should be carried onto the same structure, because "a handful of sites changed fast on the branch that made ITPR1" is a very different claim if they sit in the IP₃-binding core than if they sit in a disordered linker. | open (S17) |
| 2026-09-07 | S9a | **Synonymous saturation inside a single paralog, not only between them.** 94 % of within-ITPR1 and 85 % of within-ITPR2 pairwise comparisons exceed dS = 1.5, and dN plateaus near 0.1 while pairwise dS runs past 50 — a paralog set spanning shark to teleost to mammal has burned its fourfold-degenerate sites. The tree-based models handle this far better than pairwise ML, but no ω in S9 is a point estimate with a small error. Worth re-estimating ω **within a shallow clade** (mammals alone, or per vertebrate class) where dS is still determined, and reporting the two side by side: if the paralog ranking survives at a depth where dS is estimable, it is a much stronger claim than one taken across 450 Myr. | open (S9b or S17) |
| 2026-09-07 | S9a | **S1's toolchain manifest silently omitted a whole task's tools.** codeml, yn00, `pal2nal.pl` and `hyphy` were installed in the `piezo1` env from the start and were never probed, so S9 opened by reporting `pal2nal.pl` missing on a machine that had it. Now probed (PAML 4.10.10 via yn00, which reports a version where codeml does not; pal2nal v14; HyPhy 2.5.101). The general lesson is that the manifest's coverage is itself untested — nothing checks that every external binary the `scripts/` tree shells out to appears in `TOOLS`. A grep-based check over `subprocess.run([...])` call sites would close it. | open (S19) |
| 2026-09-07 | S8 | **S7 asked S8 which side of the vertebrate duplication each cyclostome ITPR lineage attaches to, and flanking-gene synteny cannot answer it.** All 6 loci carry 20 informative flank symbols, so the window is not the limit; the highest overlap any reaches with a gnathostome paralog consensus is 2, and random control windows reach 2. The two species' calls also point at different paralogs. The obstacle is that cross-species flank matching here is by gene *symbol*, and 550 My of rearrangement plus cyclostome annotations that name 27-29 % of coding genes leave no shared vocabulary. Needs a name-independent instrument: Compara/OrthoFinder orthology calls on the flanking proteins, or the 2R paralogon reconstructed from a cyclostome-anchored gene tree. | open |
| 2026-09-07 | S8 | **Two ohnologous flank families survive beside the ITPRs and both connect to ITPR1** — BHLHE40/BHLHE41 (ITPR1–ITPR2, 84x background) and GRM7/GRM4 (ITPR1–ITPR3, 93x) — while ITPR2 and ITPR3 share none. S7's AU test makes ITPR2+ITPR3 the sister pair. Both can be true (flank retention is deletion, not duplication order), but the manuscript's evolution section needs the paralogon stated properly: which 2R block each ITPR sits in, reconstructed from the flanking ohnolog quartets rather than from two surviving pairs. | open |
| 2026-09-07 | S8 | **The ITPR3 neighbourhood decays 2.5x faster across vertebrate classes than ITPR1's or ITPR2's** (cross-class Jaccard 0.062 vs 0.160/0.158) even though ITPR3 is the *best*-recovered paralog in S5 (232 `found_annotated` cells). In human its flanks are MHC-region genes (BAK1, ZBTB9, SYNGAP1, DAXX, TAPBP). Worth testing whether the decay is a property of the region rather than of the gene, by measuring the same cross-class flank retention for non-ITPR genes matched on genomic position. | open |
| 2026-09-07 | S8 | **131 loci now carry a paralog assignment from neighbourhood alone** (`results/synteny/unplaced_loci.tsv`, `null_verdict=supported`), 87 of them at loci whose own annotation names no paralog. These are annotation leads with independent evidence behind them and should feed S18's annotation audit and any census revision. | open |
| 2026-09-06 | S7 (user question) | **The S6 length-fit ruler falls back to a global median for any group with fewer than 5 complete-architecture records, and that inverted a plant tree tip.** `s6_select_reps.length_targets()` uses a group's own median only when the group has >= 5 non-fragment 5/5-signature records; below that it uses the global default, **2,694 aa**, which is set by the 3,749 vertebrate and 930 non-vertebrate metazoan records. Two groups fall through: **Viridiplantae** (2 such records) and **Amoebozoa** (1). `s6_rep_spec.rank_key`'s docstring states the opposite intent — "the target is the group's median so the plant grade is not scored against a vertebrate ruler" — so the fallback silently does the thing the rule says it avoids. Measured effect: every plant candidate below the top two ties on rank components 1-7, so component 8 (`length_fit`) decides the third slot. Against 2,694 aa, *Tetrabaena socialis* (2,680 aa, **3/5** signatures) is 14 aa off and wins; *Volvox carteri* f. *nagariensis* `D8UKJ4` (3,167 aa, **4/5** signatures, S20 verdict `real_gene`, profile margin 0.78) is 473 aa off and loses. Against the plant grade's own ruler — 3,055 aa over its complete-architecture records, 3,080 aa over the R5 pool — Volvox carteri is 87-112 aa off and Tetrabaena 375-400 aa off, and **the pick inverts**. So the tree's third plant tip is an artefact of the fallback, not of the stated rule. Fixing it means re-running S6's selection, the alignment and therefore S7's tree. **RESOLVED 2026-09-06 (D35)**, and the fix needed a second one beside it: with the ruler corrected the third slot went to a *second Chlamydomonas*, because `pick_diverse`'s nested key throttles on a phylum level that has only one value among the candidates. Both fixed; the effect on the whole set is **one tip** — Tetrabaena out, Volvox carteri in — and S6 and S7 were re-run end to end. | closed 2026-09-06 |
| 2026-09-06 | S7 (user question) | **The red-algal absence rests on 12 proteomes but only *one* assembly.** Rhodophyta contributes **no record at all** to census v6, and S20 swept 12 red algal reference proteomes at both E-value sensitivities with zero ITPR in every one — including *Porphyra yezoensis* itself (UP000798662, filed as *Pyropia yezoensis*, 13,100 proteins), *Porphyra umbilicalis*, *Porphyridium purpureum*, *Chondrus crispus*, *Gracilariopsis chorda*, three *Galdieria*, *Cyanidioschyzon merolae* and *Cyanidiococcus yangmingshanensis*. S23 took it to assembly level for **one** genome, *Cyanidioschyzon merolae* GCF_000091205.1, controlled (`controlled_cross_kingdom`, 5 control loci / 2 complete) → "absence holds at assembly level". That is a thinner genome-level test than the phylum deserves, and it is the *most reduced* red alga: Cyanidiophyceae are extremophile specialists with ~5,000 genes, whereas Bangiophyceae like *Porphyra* carry ~13,000, so the one assembly checked is the one where loss is least surprising. A proteome zero is a gene-caller result; if any claim leans on red algae, at least one Bangiophyceae assembly (*Porphyra umbilicalis* or *P. yezoensis*) needs the genomic sweep | open (S23 follow-up) |
| 2026-09-06 | S6 | **The cyclostome grade carries three ITPR loci per genome and the sweep's ITPR1 bait wins every one of them.** Sea lamprey (LOC116947491 / LOC116951169 / LOC116952798), inshore hagfish and *Myxine glutinosa* each carry three full-length loci; across 19 census records in that band, 6 carry a paralog label, all of them `ITPR1`, all from the S5 bait attribution and none from an annotation. A constant attribution is no information about which locus is which, so S6 enters the whole grade unlabelled. **Whether those three are 1:1 orthologs of ITPR1/2/3 (2R, shared with gnathostomes) or a cyclostome-specific expansion is the single sharpest test of the 2R hypothesis this project can run** — and it needs synteny as well as the tree, because three tips at 60-70 % identity will not separate the hypotheses on sequence alone | open (S7 tree + AU test; S8 synteny on the three lamprey loci) |
| 2026-09-06 | S6 | **Four non-vertebrate records carry a vertebrate paralog number by annotation transfer** — *Tetrabaena socialis* and *Acanthamoeba castellanii* as "type 2" in their protein name, *Branchiostoma lanceolatum* and *Phallusia mammillata* as `ITPR1`/`ITPR2` gene symbols. ITPR1/2/3 are a 2R product, so none of those can be a vertebrate paralog by descent. S6 keeps the label in its audit and groups the tips by taxonomic grade instead. Worth a line in the annotation-audit deliverable (S18): a database that ships a 2R paralog number on an amoeba is making a claim it cannot support | open (S18) |
| 2026-08-18 | setup | Confirm the external drive is attached and `data_root.txt` points at it before S4/S5 | open |
| 2026-08-18 | setup | ~~40 Viridiplantae + 41 Fungi PF08709 records exist in a family textbooks say plants and fungi lack~~ — **narrowed by S2, not closed.** The records are not scattered: in Viridiplantae *all* 20 ITPR calls are **Chlorophyta** (11 taxa, incl. *Chlamydomonas reinhardtii* with the complete five-signature architecture) and **Streptophyta — the land-plant lineage — has 0 calls from 15 records in 13 taxa**; in Fungi every call is in an early-diverging phylum (Mucoromycota 15, Chytridiomycota 6, Basidiobolomycota 3, Entomophthoromycota 1) and **Dikarya contributes no records to the search space at all**. That is the shape of a loss in the derived lineage of both kingdoms, and it is still a statement about what UniProt holds. → `results/census_v2/lineage_calls.tsv` | **closed 2026-09-05 by S20**: the answer is loss, and it is a genome fact. Sweeping the proteomes themselves — 384 land-plant and 1,353 Dikarya reference proteomes, no Pfam-annotation filter — returns **zero** in both, while the sister lineages next door return 15/48 chlorophyte and 28/1,527 fungal proteomes. All 99 plant/fungal records chased individually: 47 real genes, 52 fragmentary models, **0 contaminants** (cross-kingdom identity 19.9–45.8 %, median 24.1 %). → `results/s20_sweep/report.md` |
| 2026-08-18 | setup | ~~**Ensembl REST is unreliable right now**~~ — **diagnosed and fixed in S0.** Not general flakiness: `/xrefs/symbol/**homo_sapiens**/{symbol}` stalls indefinitely (no response, no error) for `BRCA2` as well as `ITPR1`, while the same endpoint answers in 0.6 s for `danio_rerio` and `/lookup/symbol/homo_sapiens/` answers normally. `src/databases/ensembl.py:_symbol_to_ids` now uses `lookup/symbol` with `xrefs` as fallback. Evidence: `results/s0_baseline/ensembl_endpoint_probe.tsv` | closed 2026-08-18 |
| 2026-08-18 | S0 | **Ensembl is slow enough to be a scheduling problem, separate from the stall.** Measured: 14 s for a 451-byte `lookup/symbol`, 11 s for a `lookup/id?expand=1`; one gene in one species costs ~95 s end to end. The default 8-species panel × 3 genes therefore needs ~38 min, so `run_headless`'s budget went 300 s → 900 s and a full panel sweep still needs `--species`. The real fix is to parallelise `EnsemblClient.search`'s per-species loop (the other clients already return in seconds) — a client change, not a session's worth of work, but out of S0's scope | open |
| 2026-08-18 | S0 | **Genomic span varies 6.5× across the three human paralogs (ITPR3 76 kb → ITPR2 498 kb) while protein length varies 3 %.** Found while correcting a false `[lit]` claim. Intron-content asymmetry between paralogs of identical architecture is a result, not a footnote — and D16 says the comparison must be paired within genome | open (S21) |
| 2026-08-18 | S0 | **Zebrafish `LOC101884734` (4,900 aa, 7 records) is an unnamed RyR-sized locus** returned by an ITPR-diagnostic Pfam query. First concrete instance of the unnamed-locus problem; keep it as a worked example for the correction list | open (S18) |
| 2026-09-05 | S23a | **The Apicomplexa absence has no working positive control — measured, not predicted.** Both apicomplexan classes carry a *single* PF02815 (MIR) protein across the 60 proteomes S20 swept (3–4 % clade coverage), so the control was flagged `weak` at build time. In the pilot it **failed outright**: *Toxoplasma gondii* ME49 returns no ITPR locus **and no MIR locus**, so the genome is `uncontrolled` and supports no absence claim at all. Without the control this sweep would have reported "no ITPR in *Toxoplasma*" as a finding. Apicomplexa needs a second control of a different kind before "0/60" can be taken to assembly level — a conserved single-copy ortholog panel (BUSCO-style), or the assembly's own annotated gene set used as the proof-of-search  **RESOLVED in S23b:** the fix was not a second profile but choosing the control profile *per clade, by measurement* — Apicomplexa takes Myosin_head at 36/36 and 23/23 swept proteomes, and S23c returned Apicomplexa 0/3 with all three genomes controlled. | closed 2026-09-06 |
| 2026-09-05 | S23a | **`MIN_LOCUS_IDENTITY` = 0.40 is untested outside the vertebrates — a caution, not (as first recorded here) a demonstrated failure.** The *Chlamydomonas* case that prompted this row turned out to be the bait panel: with only a distant *Cymbomonas* bait the genome's nine ITPR-family alignments came back at 24.6–28.4 % identity and the 0.40 floor discarded all nine, which is what that floor is *for*. Once B3's spread rule put *Chlamydomonas reinhardtii*'s own record in the panel the locus was recovered cleanly (`found_annotated`, 1 full locus). So the floor did not misfire — but it also has not been *measured* for this scope, and S5b's justification for it (a wide empty gap: confirmed loci ≥ 0.759, junk 0.23–0.34) was established on vertebrates, where every genome has a same-class bait. Bands here are whole phyla and 13 slots are unbaited, so a genome whose nearest bait is a phylum away could still lose a real locus to it. S23b should re-measure it S5b's way — against loci whose identity an assembly's own annotation confirms — and report how many `no_locus` calls sit just under the floor  **RESOLVED in S23c, and the answer was that the threshold does not work at all**: no value separates the two populations (confirmed down to 0.193, contradicted up to 0.318), and 0.40 would discard 87 confirmed loci, 56 of them complete genes. Retired rather than retuned — D29. | closed 2026-09-06 |
| 2026-09-05 | S23a | **The metazoan `-G` inflates locus spans, which risks *under*-counting copy number — the mirror of S5b's problem.** The metazoan max-intron is 650 kb, derived from the 325 kb *Octopus* ITPR, but a *Drosophila* Itpr spans 22 kb; in the pilot its locus clustered 20 bait alignments across **297 kb, 13× the gene**. The call is unaffected — the winning alignment covers 100 % of the annotated gene's CDS (`frac_cds` 1.0) and the gene covers 3.9 % of the locus — so status and grade are right. The risk is that two real ITPR genes within 650 kb would be absorbed into one locus and counted once, and **copy number is this task's deliverable**. S5b hit the opposite failure (inflated `n_loci`) and fixed it with `MIN_LOCUS_IDENTITY`, which cannot help here because the top alignment is perfect. S23b must measure locus-span / CDS-footprint across the sweep and, where chaining is severe, count distinct non-overlapping full-length alignments rather than loci  **RESOLVED in S23b/S23c:** measured at 35x for *Drosophila* (297,487 bp cluster / 8,514 bp footprint); copy number is now counted on non-overlapping complete alignments (D28), which recovered 7 genes in 7 genomes. | closed 2026-09-06 |
| 2026-09-05 | S23a | **Cestoda 0/11 — a metazoan clade with no ITPR record.** Fell out of deriving S20's absences from the presence table rather than restating the summary's list, along with Bacillariophyta 0/16 and Rhodophyta 0/12. Tapeworms are highly reduced parasites, so a genuine loss inside Bilateria is plausible and would be the first metazoan loss this project has seen; S23b's sweep decides it. Related: the one control bait the chimera screen rejected was a tapeworm protein UniProt calls a "MIR domain-containing protein" and `ryr.hmm` scores at 1,121 bits | open (S23b) |
| 2026-09-05 | S23a | **No Viridiplantae ITPR gene span could be measured**, because the chlorophyte census records carry locus tags NCBI has no gene record for. Every plant genome in S23 is therefore judged against the global median contiguity bar rather than its own group's. `s23_calibration.contiguity_bar()` returns the fallback in its own return value so this is visible, but a measured plant span (from *Chlamydomonas*, once S23b has a gene model for it) would close it | open (S23b) |
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
| 2026-09-05 | S5b | **Cyclostomes carry one findable ITPR, not three.** *Petromyzon marinus* and *Myxine glutinosa* are the only genomes in 309 with an `absent` call, and both are absent for ITPR2 **and** ITPR3 while carrying ITPR1. Both clear D4's contiguity bar and both fired the RyR control, so the call survives the two gates that disqualify every other candidate absence. The caveat is the panel's own: it holds **no labelled ITPR2 or ITPR3 bait for cyclostomes at all**, so the honest statement is one findable ITPR per genome rather than a paralog-resolved absence. Consistent with the ITPR1/2/3 trio arising in 2R after the cyclostome split — which is S7's question, and S13's to date. → `results/genome_ledger/genome_ledger.tsv` | open (S7 → S13 → S16) |
| 2026-09-05 | S5b | **ITPR1 is the worst-named paralog, and only a controlled comparison shows it.** Holding assembly contiguity constant over 487 found loci, an annotated gene names the right paralog for **ITPR3 88 %, ITPR2 87 %, ITPR1 65 %** of the time. Uncontrolled the same data said ITPR3 54 % vs ITPR2 14 %, which was assembly quality wearing annotation quality's clothes — the comparison reversed twice as the sweep grew (n=81 no gap, n=309 a 23-point gap). Any census built on gene symbols inherits this as an apparent difference in copy number between genes of near-identical protein length. | open (S18) |
| 2026-09-05 | S5b | **485 ITPR gene models that no name-based search can reach**: 318 exist only as DNA (no gene model at the locus, or an assembly with no gene set) and 167 sit inside an annotated gene carrying no family name. That is the measured cost of an annotation-derived census, per genome and per paralog, and it is what census v4 adds over v3. → `results/census_v4/genome_models.tsv` | open (S18, S19) |
| 2026-09-05 | S5b | **A large `-G` manufactures phantom loci in giant genomes.** In *Lissotriton* (23 Gbp) and *Protopterus* (40 Gbp), where the `-G` cap binds at 2 Mbp, miniprot chained shared-channel-module hits at 23–34 % identity across megabases into apparent loci — ITPR3 read as 10 loci in the newt. The real gene was always the best-scoring one so no *status* was wrong, but copy number was, and copy number is a result. Fixed by `MIN_LOCUS_IDENTITY` = 0.40, measured from a clean gap (571 annotation-confirmed loci have minimum identity 0.759; the junk is below 0.35). **S20 and S23 sweep genomes of similar size and must carry this floor.** | open (S20, S23) |
| 2026-09-04 | S5a | **The contiguity bar and the margin rules select very largely the same genomes.** 120 of the 309 manifest assemblies have a contig N50 below the median measured ITPR genomic span (142,212 bp) and so cannot carry the gene on one contig — but not at random: **66 % of Aves against 11 % of Actinopteri, and 68 % of margin species against 12 % of order representatives.** S4 chose the margin set for what its *proteomes* lacked; it turns out to be mostly the set whose *assemblies* cannot hold the gene. Every downstream absence claim in birds is confounded until the two are separated, and S19's contiguity floor is now a prerequisite for S15 rather than a methods footnote. → `results/genome_ledger/contiguity_bar.tsv` | open (S19 → S15) |
| 2026-09-04 | S5a | **Recovery from a fragmented assembly tracks gene span, so the detection bias points the same way as the loss signal.** Below the contiguity bar the pilot finds ITPR3 (median span 82 kb) in 2 of 3 genomes and ITPR1 (186 kb) and ITPR2 (231 kb) in 0 of 3; above it, 9 of 9. S0's 6.5× span asymmetry between paralogs of near-identical protein length was filed as a curiosity — it is a systematic bias in which paralog a poor assembly appears to have lost. Six genomes is a direction and a mechanism, not a result; S5b tests it at 309 and S19 must model it. | open (S5b → S19 → S15) |
| 2026-09-04 | S5a | **The paralog labels run out below the well-annotated clades.** Six bait slots could not be filled from the whole census — ITPR1 and ITPR2 in chondrichthyans, ITPR2 in the coelacanth grade, all three in cyclostomes — because no labelled, full-length, complete-architecture record exists there at all. The `vertebrate_basal` baits cover those genomes without a paralog claim, so the sweep finds the loci but cannot name them; S7's phylogeny is what assigns them, and until it does, no per-paralog statement may be made about sharks, lampreys, hagfish or the coelacanth. → `results/s5_baits/bait_build_stats.json` | open (S7) |
| 2026-09-04 | S5a | **Teleost ITPR1 3R duplicates are half-named in RefSeq.** Both *Takifugu rubripes* and *Genypterus blacodes* carry two ITPR1 loci, one annotated `itpr1b` and the other left as an unnamed `LOC` — in *Takifugu*, `LOC101074739`. Two RyR loci in the same genome are unnamed the same way. A name-based census cannot reach these, which is why the sweep carries every locus whose annotation does not already name it correctly into census v4. A concrete, paired worked example for the correction list. | open (S18) |
| 2026-08-18 | setup | AlphaFold DB returned models for 8/8 human ITPR queries in the smoke test — better coverage than the PIEZO family had. Worth checking early whether AFDB covers full-length ITPRs or only fragments, since it changes S11's scope | open (S11) |
| 2026-09-05 | S20 | **The fungal losses are patchy, not basal.** The family is present in Mucoromycota (18/34), Chytridiomycota (6/16), Basidiobolomycota (2/2), Zoopagomycota and Entomophthoromycota — and absent from Glomeromycota (0/27), Mortierellomycota (0/18), Kickxellomycota (0/35) and Microsporidia (0/29) as well as all of Dikarya. That is repeated independent loss across the fungi, not one event at their base, and it is a different claim from the one S2 supported | open (S23) |
| 2026-09-05 | S20 | **Apicomplexa 0/60** — the largest zero-hit protist clade in the sweep, and a group with a well-studied calcium biology. Bacillariophyta 0/16 and Rhodophyta 0/12 sit beside it. Whether these are losses or assembly/annotation gaps needs an assembly search | open (S23) |
| 2026-09-05 | S20 | ***Cymbomonas tetramitiformis* contributes 12 of the 22 complete green-algal genes.** Either a real copy-number expansion in the prasinophytes or a duplicated assembly; a proteome cannot tell those apart | open (S23) |
| 2026-09-05 | S20 | **Both families are already full length in the closest unicellular relatives of animals.** *Salpingoeca rosetta* (choanoflagellate, 5,340 aa, not a profile seed) and *Capsaspora owczarzaki* (filasterean, 6,625 aa, but a `ryr.hmm` seed so circular) each carry an architecture-level ryanodine receptor alongside an ITPR. The ITPR/RyR duplication predates Metazoa, which is a rooting constraint S7 can use | open (S7, S13) |
| 2026-09-05 | S20 | **UniProt's `tax_id:` search field is hierarchical and silently truncating.** A batched `tax_id:A OR tax_id:B …` query matches every descendant of each id, so a 100-term batch can match far more than 100 records, the page caps at `size`, and requested taxids come back missing with a 200 OK. Anything in this project that resolves taxonomy in batches must use `taxonomy/taxonIds/` and verify the returned set. Fixed in `s20_taxa.py`; worth checking if any earlier task did the same | open |

---
| 2026-09-06 | S23c | **`scripts/s20_report_results.py` is 691 lines, over the project's 500-line rule.** A standing violation inherited from S20, found while checking S23's own files. Not fixed in S23c because touching it risks S20's committed report and the task in progress was S23. The natural seam is the same one S23 used: the scope/instrument half against the results half. | open |

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

**D27 — A threshold cannot be calibrated from the population it has already
filtered.** S5b could measure `MIN_LOCUS_IDENTITY` = 0.40 only because its
sweep had *recorded* the loci that 0.40 excludes; a sweep that applies a floor
at search time has thrown away the evidence its own calibration needs, and
whatever it reports afterwards is a description of the filter rather than a
measurement of the threshold. So a threshold that is going to be measured gets
two values: a **recording floor** set far below any plausible answer, and a
**call floor** applied downstream and read back from the committed
calibration (D13 applied to constants). S23b's sweep records at 0.15 and calls
at whatever `locus_calibration.json` says, with
`s23_calibration.call_min_identity()` returning *whether* the number in force
is measured or inherited so no report can present one as the other.

Two corollaries. **(a) The calibration must refuse a sample too thin to be a
measurement**, and this is not hypothetical: a smoke-test run of
`s23_calibrate_loci.py` over one genome derived 0.25 from two confirmed loci,
wrote it, and the next sweep read it back and moved three genomes out of
`no_locus`. The guard is `s20_test_sensitivity.py`'s refusal to pass
vacuously, applied to a calibration. **(b) The evidence a floor is measured
against must not come from the instrument being calibrated.** Here it is the
assembly's own annotation — which is also why the family-name list has to be
right outside the vertebrates, and why `s5_classify`'s drifted local copy,
missing `itr-1`, was putting a correctly-recovered nematode receptor in the
wrong population.

**D29 — A similarity threshold measures distance to the reference panel, not
membership, and outside a well-sampled clade those are different things.**
S5b's `MIN_LOCUS_IDENTITY` = 0.40 works in the vertebrates because every
genome there has a bait from its own class, so a locus's identity to its
nearest bait really does track whether it is the gene. S23c measured the same
statistic over 917 clusters in 194 non-vertebrate genomes and it does not
separate: confirmed loci reach down to **0.193** while contradicted ones reach
up to **0.318**, and the inherited 0.40 would discard **87 confirmed loci, 56
of them complete gene models**. Coverage is no better (Youden J 0.71 against
identity's 0.70) — a draft of this decision claimed coverage separates, and
the measurement contradicted it. The reason is structural rather than
accidental: with bands a whole phylum wide, identity to the nearest bait is a
measure of *how far away the nearest bait is*.

So a family call is not made on a similarity threshold outside the
vertebrates. It is made by the profiles (D14/D23), which is where this project
puts every other family call, and the threshold is retired rather than
retuned. **A retired threshold has to be validated against evidence the
replacement does not share**: here the assembly's own annotation, where the
profile gate agrees with **10 of 10** loci named for the family and declines
**10 of 11** named for something else. That axis is small — 21 of 917 clusters
carry an informative gene name at all, because outside the vertebrates gene
models carry locus tags — and its being small is reported rather than smoothed
over, since it is the reason the second axis was needed.

**D28 — A cluster is not a gene, and where a count is the deliverable the
difference has to be paid for.** Loci are alignment clusters chained at
`LOCUS_GAP` inside a `-G` window, and outside the vertebrates that window is
650 kb: *Drosophila*'s 22 kb *Itpr* sits in a **297,487 bp cluster around an
8,514 bp CDS footprint, 35x**. S5 could tolerate this because its best locus
always won and the cell status was right either way. S23 cannot, because two
real genes inside one such chain would be counted **once** and copy number is
the result. A copy is therefore defined on the evidence rather than on the
clustering: a genomic interval carrying an alignment that covers a *complete*
bait, with two overlapping alignments being one copy (several baits hitting
one gene is the normal case) and two non-overlapping ones being two, however
the clustering happened to group them. The span-to-footprint ratio is reported
per group so a reader can see where it bites.

**D26 — A positive control does not port across the scope it controls; it
has to be re-chosen for the clade the claim is about, and graded.** D25 says
a ported *threshold* must be re-measured. S23a found the same is true of the
*control itself*, and for a sharper reason: a control that is absent for
biological reasons cannot distinguish "the gene is not here" from "the search
did not run". S5 rests every one of its 309 genomes on the RyR positive
control, which is sound there because every vertebrate has three RyRs. S20
measured architecture-level RyR in **2 of 6,928** non-metazoan proteomes, so
carrying that control into the plant and fungal genomes would have made every
absence claim in S23 unfalsifiable while looking controlled. The replacement
is the MIR-domain sharer (PF02815), chosen because it is inside the family's
own signature set and S20 already used it as its in-search control at
proteome level — and it is drawn **per control clade**, because a chlorophyte
mannosyltransferase is not a control for *Arabidopsis*.

Two corollaries, both of which cost this session a rebuild. **(a) A control's
selection rule must not be tuned to what the control usually looks like.** The
first build required 600 aa, which selects the full mannosyltransferases and
silently left five clades — including both apicomplexan classes — with no
control at all, i.e. it deleted the hardest negative claims rather than
reporting them. **(b) A control has to be graded, not merely present.** The
apicomplexan clades carry one PF02815 protein each across 60 swept proteomes;
that control is real but thin, and reporting it identically to a
full-length mannosyltransferase found in 100 % of a class would overstate
every claim resting on it. `control_strength()` records `strong` / `weak`
with the reason.

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

**D30 — A paralog label is only meaningful inside the vertebrates, and a
bait attribution is only a label where it varies** (S6, 2026-09-06). Two
halves of one rule, both measured rather than asserted. *Outside the
vertebrates*: ITPR1/2/3 are a 2R product, so a *Tetrabaena* or
*Acanthamoeba* record whose UniProt protein name reads "receptor type 2"
carries that number by annotation transfer from a vertebrate, not by
descent — four representatives were affected. The label is kept in the
audit with its source; the figure group is the taxonomic grade. *Inside the
vertebrates*: the sweeps record which bait won each locus, and for most
clades that discriminates, because the three baits win different loci. In
the cyclostome band it does not — the sea lamprey, the inshore hagfish and
*Myxine glutinosa* each carry **three** full-length ITPR loci and the ITPR1
bait wins **all** of them. A constant attribution is no information about
which locus is which, and filling an `ITPR1@cyclostomata` grid cell from it
would assert an orthology the data does not contain.
`s6_rep_spec.informative_attribution()` decides this per band by counting
distinct attributed cells, so the exception is found rather than named.
Measured consequence: the whole cyclostome row of the paralog grid is
reported unfilled — which is a result about how the vertebrate base is
annotated, not a gap in the search — and the grade enters the alignment as
complete, unlabelled copy sets instead.

**D35 — A group is never measured against another group's ruler, and a
key level with one value is not a spread** (S6 re-run, 2026-09-06). Two
selection rules were stated correctly and implemented in a way that did
the opposite, and together they cost the tree a real tip. (a)
`length_targets()` used a group's own median length only when the group
held at least five non-fragment *complete-architecture* records, and
otherwise borrowed the **global** median — 2,694 aa, a number set by the
3,749 vertebrate and 930 invertebrate records. Viridiplantae (2 such
records) and Amoebozoa (1) fell through, though `rank_key`'s own
docstring promised "the group's median so the plant grade is not scored
against a vertebrate ruler". The target now comes from the group's own
non-fragment records carrying **at least four** of the five signatures,
which barely moves a well-populated group (Vertebrata 2,671 → 2,671, SAR
2,973 → 2,978) and gives the starved ones a real ruler (Viridiplantae 2
→ 13 records, 3,182 aa; Amoebozoa 1 → 4, 2,886 aa); the tier that
produced each target is committed in `length_targets.tsv`, so a number
computed from four records is visible as such. (b) `pick_diverse()`'s
nested key throttles on **every** prefix level, including one that has a
single value across the candidates — and a level with one value carries
no diversity information, it can only throttle. Every plant candidate is
Chlorophyta, so the phylum prefix admitted one row per wave and by wave 3
the genus quota was loose enough to pass a *second* Chlamydomonas.
Single-valued levels are now skipped. **Measured effect of both, on the
whole 134-tip set: one tip.** *Tetrabaena socialis* (2,680 aa, 3 of 5
signatures) out, *Volvox carteri* f. *nagariensis* (3,167 aa, 4 of 5,
S20 verdict `real_gene`, profile margin 0.78) in. Nothing else moved,
because every other cell's top key level is multi-valued and every other
group already had five clean records of its own. Found by a reader
asking why a well-known genome was missing from the figure — the D11
check ("look at the figure") arriving from outside.

**D36 — A CDS is accepted only if it encodes the aligned protein, and a
route returns the first candidate that *validates*, not the first that
downloads.** (S9a, new.) dN/dS is a statement about codons, so a codon
alignment is only meaningful if every nucleotide sequence in it provably
encodes the exact protein the alignment and the tree were built from. Two
consequences, both of which the first implementation got wrong.

*Which sequence.* A UniProt entry cross-references every Ensembl transcript
of its gene, and the entry's own sequence is one particular isoform. Human
ITPR1 lists **five** and human ITPR2 **two**, and in both the first is not
the isoform S6 aligned. Returning the first CDS that fetches does not fail —
it substitutes a different isoform for the protein the tree was built on,
which is a *wrong* codon alignment rather than a missing one, and it fails
silently on exactly the two records the paper leans on hardest. The routes
are therefore candidate generators, the caller translates each and keeps the
first that matches, and the ones rejected on the way are recorded (human
ITPR1 needed nine).

*What to do with the residue that still disagrees.* Every position where the
translation and the aligned protein differ is masked to `NNN` — not only
internal stops, which is all the ported implementation masked. A codon that
encodes something else is not that protein's codon at that site whatever the
cause, and codeml has to see missing data rather than a residue that is not
there. This is also what lets a locus rerun that reproduces a genome gene
model to within 1 % be used at all: the nine N-terminal residues miniprot
places differently on an 800 kb region than on a 3.2 Gbp genome are masked,
and the other 2,657 are exact.

**D37 — Branch-site model A is restarted from several initial ω by
construction, not repaired afterwards.** (S9a/S9b, new.) A nested
alternative cannot have a lower optimum than its own null, yet codeml
reaches one routinely on alignments of this size — the PIEZO project's first
`bs_alt` converged **below** its own `bs_null` and needed a restart script
written after the fact to recover. Running four initial ω from the start
makes "the best of several optima" the reported number rather than a repair,
and every restart is committed to `bs_restarts.tsv` with a below-null flag,
so a stem the search is *still* stuck on is visible in the data instead of
being reported as a negative result. Two related rules travel with it: the
branch-site LRT is tested against the **50:50 mixture** its boundary null
implies rather than a plain χ²₁, and BH correction runs across the whole LRT
family because S9 asks the same test once per paralog.

**D38 — A significant likelihood-ratio test is not a claim; the fitted
parameter is.** (S9b, new.) Three of S9's model comparisons came back
significant and none of the three meant what its p-value appeared to say,
each for a different reason, and all three were only visible in a number
codeml prints *outside* the test:

* **M8 vs M7** beats its null in all three paralogs at q ≈ 2×10⁻⁴, and the
  class it adds sits at **ω = 1.00000** — codeml's boundary — on 0.3–0.7 %
  of sites. A beta distribution on [0, 1] cannot represent a spike at the
  neutral boundary, so adding one class that lands exactly there improves
  the fit significantly and says nothing about adaptation. What the test
  found is that under 1 % of sites are *unconstrained*.
* **M2a vs M1a** returns 2ΔlnL = 0.00 with the positive class estimated at
  ω = 20–94 and a proportion of **exactly zero**: the alternative has
  collapsed onto its own null, and the large ω is a parameter fitted to no
  sites.
* **Branch-site model A** is significant on all three stems, and on two of
  them ω₂ is pinned at codeml's **999 upper bound** with the likelihood
  flat above it — ITPR2's restarts reach the same lnL at ω₂ = 162 and 999.
  That is an unidentified parameter, not a large effect, and it is exactly
  what §4.2's saturation predicts for a branch of that age.

So every model comparison in this project reports the parameter its claim
is about beside the test statistic, and the report's wording is chosen from
the parameter. The rule generalises the one D22 applies to profile hits: a
score is evidence only when the thing it scores is present.

**D32 — A topology constraint names only the taxa its hypothesis is
about** (S7, 2026-09-06). IQ-TREE's `-g` places freely exactly the taxa a
constraint *omits*; a taxon that is **listed** — even in a top-level
polytomy, which reads as "free" and is not — is thereby fixed outside every
group the constraint declares. S7's three sister constraints listed all 134
tips, and so silently also required the fourteen unlabelled vertebrate tips
that this tree nests *inside* the paralog clades to sit outside them. That
cost is identical across the three hypotheses and has nothing to do with
which two paralogs are sisters, so the AU test rejected all three at
ΔlogL ≈ 1,500 — including the arrangement the ML tree itself holds at
SH-aLRT 100 / UFBoot 100. This is the same failure D14's tree-corrected
membership was written to prevent, arriving through the constraint's
*taxon set* instead of through its labels. The rule is now enforced rather
than remembered: each constraint carries the three paralog cores and the
outgroup and nothing else, and `s7_test_tree.py` T5 fails a constraint that
names a free tip. The outgroup must stay, because on three groups alone
`((A,B),C)` is an unrooted trifurcation that asserts only what every
hypothesis shares.

**D33 — Co-orthologs of a duplication older than the species are expected
*not* to be sisters** (S7, 2026-09-06). S6 put the teleost 3R pairs in the
representative set as a stress test on its own naming, with the stated
expectation that a species' two same-paralog copies come back as sisters
and that anything else means the naming or the alignment is wrong. That
expectation is wrong for a duplication the species inherited: 3R predates
the teleosts, so zebrafish *itpr1a* belongs with other teleosts' *itpr1a*,
not with zebrafish *itpr1b*. Both S7 sets duly failed the test, and both
are broken **only by other tips of the same paralog** — which is the
signature of the shared duplication, not of a naming error. So "not
sisters" is two findings, and the report separates them by what the copies'
smallest containing clade holds: same paralog only → a duplication older
than the species, as expected; anything else → the failure S6 was looking
for. S8's synteny and S9's ohnolog assignments inherit the distinction.

**D34 — A search's reported tree and its best tree are different questions,
and both get stated** (S7, 2026-09-06). A *constrained* search explores a
different path through tree space and on a 134-tip alignment can land on a
better optimum than the free search did. S7's ITPR2+ITPR3 constraint reached
−214,733.9 against the unconstrained search's −214,743.3 — a constraint
costing **negative** likelihood, which is a fact about the search and not
about the hypothesis. The tree reported stays the unconstrained one, because
substituting a tree found under a constraint makes the topology partly an
assumption and every downstream table is built on it; but it is then not the
global optimum, and what has to be said is *which claims differ*. Silently
swapping the better tree in, and silently not mentioning it, are both
wrong. `s7_bnni.py --tree/--label/--out` scores the committed claim set
against any alternative tree so the answer is a table, and the report
renders "yes, it is the best" when no constrained search wins.

**D31 — An alignment has two coordinate systems, and the map between them
is data** (S6, 2026-09-06). Anything that starts from a protein residue —
the domain track here, S11's clinical variants, S22's ligand sites — has to
cross from residue numbering to *trimmed* column numbering, and the obvious
route is wrong in a way that renders: counting ungapped positions in the
trimmed row silently renumbers every residue after the first discarded
column. It put human ITPR1's domains at roughly two thirds of their true
position and dropped the two C-terminal ones off the end of the figure, and
the figure looked plausible. `scripts/s6_msa.py` now runs trimAl a second
time with `-colnumbering` and commits the result as
`results/msa_v2/column_map.tsv`; with the map applied the pore lands on the
alignment's highest conservation peak, which nothing in the mapping knows
about. Caught by D11 — looking at the figure — and by no table check.

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
