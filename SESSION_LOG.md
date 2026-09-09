# SESSION_LOG.md — running notes per session

Chronological record of what was run, what resulted, and what is next.
`FINDINGS.md` holds the biology; this file holds the work.

---

## 2026-08-18 — S-setup: project created

**What ran.** Created this project as a working copy of the PIEZO project's
machinery pointed at the IP3-receptor family.

Ported from `../piezo_genes` and adapted:

- `src/` — the whole app (databases, analysis, discovery, investigation,
  GUI, CLI). The family is now defined in exactly one place,
  `src/utils/family.py`; nothing else in `src/` hard-codes a gene name.
  `interpro.has_piezo_signature` / `batch_fetch_piezo_signatures` became
  `has_family_signature` / `batch_fetch_family_signatures` and read their
  Pfam list from there. Verified: every module imports and the wiring
  resolves to the ITPR constants.
- `src/utils/data_root.py` — added `require_data_root()` and
  `python -m src.utils.data_root --require`, which exits non-zero when the
  external drive is absent. The old behaviour silently fell back to the
  internal disk, which is fine for a cached API page and fatal for a
  400 GB genome download.
- `scripts/dashboard.py` — the PIEZO-specific sweep probes were replaced by
  the generic `results/session_live.json` contract, figures are now
  auto-discovered from `results/*/figures/*.png` as well as listed
  explicitly, and the ledger parser now accepts sub-task ids (S14a), which
  the original silently skipped.
- `scripts/figstyle.py` — palette re-keyed to ITPR1/2/3 with the ryanodine
  receptors as a distinct accent group.
- `scripts/s14_*.py`, `build_findings_page.py` — the manuscript-assembly,
  claim-checking, deposit and findings-page machinery, with the figure map
  and section list set to this project's plan and the claims list emptied
  down to a worked example of each operation.
- `presets/` — six ITPR presets, including a domain scan that excludes the
  ryanodine receptors explicitly.

**What was checked live** (not assumed):

- UniProt: Q14643 / Q14571 / Q14573 are human ITPR1/2/3 at 2,758 / 2,701 /
  2,671 aa; P21817 / Q92736 / Q15413 are RYR1/2/3 at 5,038 / 4,967 /
  4,870 aa.
- InterPro: PF08709 (Ins145_P3_rec), PF01365 (RIH), PF08454 (RIH_assoc),
  PF02815 (MIR) — and human RYR1 carries **all four**, which is why D14
  exists.
- Taxonomic distribution of PF08709 and PF01365 across Metazoa, SAR,
  Discoba, Fungi, Viridiplantae, Amoebozoa, bacteria, archaea and several
  named species. Recorded in `docs/ip3r_background.md` §4 as a planning
  snapshot to be re-derived in S2.
- Zebrafish carries itpr1a, itpr1b, itpr2, itpr3 (2,634–2,819 aa), with
  RyRs returned by the same query at 4,856–5,109 aa.

**What was written.** `PUBLICATION_ROADMAP.md` (27-task ledger, session
protocol, storage rules, and the Decisions log D0–D17 carrying the PIEZO
project's methodological rules forward), `docs/session_briefs.md` (a brief
per ledger row), `docs/ip3r_background.md`, `docs/analysis_catalogue.md`,
`CLAUDE.md`, `INTERFACE.md`, `README.md`, `FINDINGS.md`, `roadmap.md`,
`manuscript/README.md`.

**Smoke tests run** (the point being to prove the port works, not to
produce results — nothing was saved into `results/`):

- `--preset ip3r --species "Homo sapiens"` → **115 variants**: NCBI 43,
  UniProt 31, Ensembl 41.
- `--preset ip3r_discovery --analyze --discover` → 56 variants, 35 with
  sequences aligned, 8 family-signature blocks derived, **0 novel
  candidates** — the correct answer for a query that returns only the three
  known human paralogs, and a check that the evidence gate is not promoting
  known genes.
- Ensembl REST misbehaved: 45.6 s timeout on the first run, then 55–73 s
  successes, then an HTTP 500. Fixed what could be fixed on this side —
  `EnsemblClient._get` now retries 3× with backoff on read timeouts and 5xx,
  and `src/cli.py` raises the per-request timeout from 45 s to 90 s — and
  logged the rest as an Emergent row. NCBI, UniProt and AlphaFold were
  reliable throughout.

**State.** `results/` is empty by design. `data_root.txt` points at
`/Volumes/FANTOM/IP3R_DATA`; the drive was **not attached** at setup, and
`--require` correctly refuses.

**Next session: S0** — verify every `[lit]` claim in the background
document against primary sources, write `docs/ip3r_review_2026.md`,
re-derive the database snapshot, and smoke-test the app against live APIs.
S0 needs no bulk storage, so it can run before the drive is attached.

---

## 2026-08-18 — S0: literature baseline + scope confirmation

**Task.** S0 (topmost `pending`, no dependencies). Completed.

**Protocol.** Dashboard opened + watcher backgrounded. No git remote, so
nothing to pull. `python -m src.utils.data_root --require` → **exit 1**: the
FANTOM drive is not attached. S0 needs no bulk storage, so the session
continued on that basis and the drive remains a hard blocker for S4/S5.
Note: the anaconda `base` python has no `biopython`; the `piezo1` env
(Python 3.11.15, biopython 1.87, requests 2.34.2) does, and everything here
ran under it. **S1 should record that env in `toolchain_manifest.txt`.**

**What ran.**

- `scripts/s0_db_snapshot.py` — re-derives every `[db]` number in
  `docs/ip3r_background.md` (InterPro signature counts, taxonomic
  distribution, UniProt reference + sister-family panel, per-protein Pfam
  architecture, the zebrafish PF08709 query). **All reproduced exactly.**
- `scripts/s0_gene_structure.py` — exon counts, genomic spans and cytobands
  for ITPR1/2/3 from Ensembl `lookup/symbol?expand=1`.
- Literature pass — 19 atomic claims extracted from the 11 `[lit]` blocks,
  checked against 51 references (45 primary) via Europe PMC.
- `scripts/s0_report.py` — renders `results/s0_baseline/report.md` purely
  from the committed tables (D13).
- App smoke tests: `--preset ip3r` and `--preset ip3r_zebrafish`.

**Results.**

- **Literature.** 12 claims verified as written, 3 kept but qualified,
  **2 struck**, 1 retagged `[open]`, 1 upgraded to `[db]`.
  `docs/ip3r_review_2026.md` is the verified baseline;
  `results/s0_baseline/lit_claims.tsv` is the audit trail.
- **The exon/span claim was false.** Measured: ITPR1 62 exons / 354,174 bp,
  ITPR2 57 / 497,888 bp, ITPR3 58 / **76,245 bp**. Range is 57–62 exons, not
  58–60, and ITPR3 does not span "hundreds of kb". Span varies **6.5×** across
  paralogs while protein length varies 3 % → Emergent row for S21.
- **"Both families independently expanded to three paralogs" was retagged
  `[open]`** — the independence of the two triplications is question Q2. It
  was about to be an assumption feeding its own answer.
- **Gillespie syndrome was corrected**: both recessive biallelic *and* de novo
  dominant-negative mechanisms. ITPR3's phenotype was extended: p.Arg2524Cys
  causes a multisystemic immunodeficiency, not only CMT1J. Both → S17.
- **D14 was measured, not asserted.** The zebrafish query
  `taxonomy_id:7955 AND xref:pfam-PF08709`, quoted in the planning document
  as evidence for four IP3 receptors, returns **109 records / 10 gene
  symbols — 53 (49 %) of them ryanodine receptors**, including an unnamed
  4,900 aa `LOC101884734`. In that dataset the length band happened to
  separate the families perfectly; that is annotation luck, not a rule.

**The Ensembl problem — two faults, both diagnosed, both addressed.** The
setup session logged "Ensembl REST is unreliable". That was one word for two
different things.

*Fault 1, a per-species stall.* Ensembl is up (`/info/ping` 0.65 s, release
15.12), but **`/xrefs/symbol/homo_sapiens/{symbol}` stalls indefinitely** —
no response, no error, for `BRCA2` as well as `ITPR1` — while the same
endpoint answers in **0.6 s for `danio_rerio`**. That is exactly why the
human preset lost Ensembl entirely while the zebrafish preset got all three
sources. `_symbol_to_ids` now resolves through `/lookup/symbol/`, and falls
back to `xrefs` **only on a clean 404** — never after a transport error,
which would walk straight back into the stall. Measured effect: human ITPR1
went from **0 variants to 24**.

*Fault 2, latency — only visible once fault 1 was fixed.* Ensembl's speed is
unstable: the same 451-byte `lookup/symbol` call measured **0.61 s, 7.61 s
and 13.91 s** in one session, and the `expand=1` call that carries the
transcript/exon payload costs ~12 s every time. One gene in one species is
~12 s at best, 95 s at worst; the default 8-species panel × 3 genes is 24
sequential pairs, i.e. **5–38 minutes**, straddling `run_headless`'s 300 s
budget. Raised to 900 s. Parallelising the per-species loop is the real fix
→ Emergent.

**Four smoke runs, recorded in `results/s0_baseline/smoke_test.tsv`.** The
decisive one is run 4: `--preset ip3r --species "Homo sapiens"` → UniProt 34,
NCBI 44, **Ensembl 41 in 28.5 s, 3/3 sources**. Run 2 (zebrafish) was also
3/3. The 8-species panel is still a 2/3 run for the latency reason above, and
should be driven with `--species` for now.

A retry budget could never have fixed fault 1, and no amount of retrying
diagnoses fault 2 — worth remembering the next time a source "is flaky".

**Files.** `docs/ip3r_review_2026.md` (new); `docs/ip3r_background.md`
(corrected, struck text left visible); `results/s0_baseline/` (11 tables +
report); `scripts/s0_db_snapshot.py`, `scripts/s0_gene_structure.py`,
`scripts/s0_report.py`; `src/databases/ensembl.py` (fix).

**Next session: S1** — toolchain install + the positive/negative control
benchmark. Read `results/s0_baseline/report.md` first: the 49 % RyR
contamination rate is the reason the decoy panel matters, and the fact that
the length band worked perfectly in that one dataset is the reason S1 must
show the *scorer* separates ITPR from RyR without leaning on it.

---

## 2026-08-18 — S0 addendum: the literature baseline expanded to a full review

**Request.** The S0 review document was judged too thin; it was rebuilt as a
publication-scale review, fully referenced, with a typeset PDF.

**What changed.**

- **Bibliography 51 → 137 references** (117 primary, 20 review). Candidates were
  harvested from Europe PMC across ~60 topic queries covering discovery,
  structure, gating, regulation, cell physiology, paralogues, evolution,
  genetic models, disease and pharmacology — each query run twice, ranked by
  citation count *and* by recency, so the set is neither purely canonical nor
  purely recent. Candidates were curated by hand; **all bibliographic metadata
  was then fetched programmatically by PMID**, so nothing in the reference
  table was transcribed by hand.
- **`docs/ip3r_review_2026.md` is now generated, not written.** Source is
  `docs/review/00_frontmatter.md` … `12_methods_audit.md` (13 files, each under
  the 500-line limit); `scripts/s0_review_build.py` assembles them, renumbers
  the stable `[Rnn]` keys into order of first appearance, renders the
  bibliography from `references.tsv`, and fails the build if a cited key has no
  reference row. The assembled file is 1,022 lines — over the 500-line rule,
  which is why it is a build artefact with modular sources, the same pattern
  `manuscript/` uses.
- **PDF**: `docs/ip3r_review_2026.pdf`, 24 pages A4, pandoc 3.6.4 + xelatex,
  Palatino, running heads, TOC, booktabs tables, hyperlinked PMIDs and DOIs.

**Three build traps, fixed and commented in the script.** (1) Pandoc passes raw
HTML through and the LaTeX writer then *drops* it — every `<sub>`/`<sup>` was
silently flattened, so "IP₃R" typeset as "IP3R". The build now rewrites them
into pandoc's `~x~`/`^x^` syntax while the markdown source keeps the HTML form
so it still renders on a web front end. (2) The YAML title plus the source H1
printed the title twice and gave it its own TOC entry; the build drops the H1.
(3) `linkcolor: [RGB]{...}` cannot survive YAML → LaTeX; the colour is now
defined in `header-includes` and referenced by name.

**Accuracy pass (D11 applied to prose).** Every page of the PDF was rendered and
read. One outright mis-citation was caught and fixed — a paper on end-stage
heart failure had been cited for a claim about antibody disagreement — and four
wordings were pulled back to what their sources actually say (a purification
paper does not establish loss in ataxic mutants; an autoradiography claim
became "highest density in cerebellum"; a smooth-muscle isolation paper is not
a reconstitution replicate; *Drosophila* itpr disruption affects metamorphosis
and ecdysone release, not moulting).

**One inconsistency this surfaced.** `references.tsv` now backs both the claim
audit and the review, so `s0_report.py` was counting 137 references for a
19-claim audit. It now derives the audit bibliography from `lit_claims.tsv`.
That exposed two references (R06, R50) sitting in the audit set but attached to
no claim; both genuinely support claims C01 and C07 and were attached, so the
audit is a consistent 51 and the review a consistent 137.

**Not changed:** no ledger status, no roadmap task. S0 remains `completed`;
this is a deepening of its deliverable, not new scope.

---

## 2026-08-18 — S1: toolchain + control benchmark

**Task.** S1 — install the toolchain, prove MAFFT is really wired in, and
measure the discovery scorer's recall and specificity against ground truth,
with the ryanodine receptors as the sharp decoy.

**Session-protocol notes.** The external drive was **not attached**
(`/Volumes/FANTOM` absent). S1 needs no bulk storage, so it was worked
anyway — as S0 was — and nothing was downloaded to the internal disk. S4/S5
remain blocked. Separately, `python -m src.utils.data_root --require` was
dying with `ModuleNotFoundError: No module named 'Bio'` before it ever
reached the drive check, because `src/utils/__init__` → `results_writer` →
`src.core.__init__` → `search` → `databases` → `Bio.Entrez`. The protocol's
step-3 gate has to work in a bare interpreter, so `SearchOrchestrator` is now
re-exported lazily (PEP 562). The gate now reports the drive.

**Toolchain (step 1).** All 12 binaries resolve, versions recorded in
`results/toolchain_manifest.txt` by `scripts/s1_toolchain.py`: MAFFT v7.526,
HMMER 3.4 (hmmbuild/hmmsearch/jackhmmer), BLAST+ 2.16.0+ (blastp/tblastn/
makeblastdb), miniprot 0.18-r281, trimAl v1.5.rev1, IQ-TREE 2.3.6, NCBI
`datasets` 18.35.0, Foldseek 10.941cd33. Nothing needed installing. BLAST+,
`datasets` and Foldseek are not on the bare PATH — they live in the reused
`piezo1` conda env (python 3.11.15, biopython 1.87), which is now recorded as
Decision **D18**; the manifest marks each tool `PATH` or `env` rather than
letting an env-only tool pass as a PATH tool.

**MAFFT wiring (step 2).** `MafftTracer` wraps `subprocess.run` and records
every mafft call with its return code — necessary because
`alignment.py:_mafft` falls back to the star alignment *silently* on a
non-zero return code. One call, `mafft --auto`, rc=0, 56 sequences → 9,494
columns. Trace committed.

**Panels.** 25 positive controls (ITPR1/2/3 across human, mouse, rat, chicken,
*Xenopus*, zebrafish incl. itpr1a/itpr1b, cow, plus fly `Itpr`, worm `itr-1`,
*Dictyostelium* `iplA`) and 31 decoys (6 RyRs across 3 species, POMT1/2 as
MIR-domain sharers, 10 in-band channels, 11 in-band non-channels, 2 giants),
all fetched live from UniProt and committed as JSON + FASTA. Sea urchin Itpr
did not resolve and is recorded in `panel_positives_missing.txt`.

*One panel bug worth remembering:* the Swiss-Prot preference test was
`"reviewed" in entryType` — and `"UniProtKB unreviewed (TrEMBL)"` contains
`"reviewed"`, so every TrEMBL entry ranked as reviewed and the tie broke on
length. Human RYR2, FLNA and rat Itpr1/2/3 all came back as the wrong (longer,
unreviewed) entry. Fixed and re-fetched; the panel is now canonical Swiss-Prot
wherever one exists.

**The result that changed the code.** On the first run **all six RyR decoys
were promoted at 45** — `pfam+20,cluster+10,breadth+15`. They carry all four
family-diagnostic Pfams, so the domain component fires by construction and
also satisfies the D3 evidence gate. Specificity was 25/31 (81 %) and *every
failure was a ryanodine receptor*. The roadmap anticipated this and required
D14 to be strengthened before S2, so S1 implemented it: a labelled-bait
sister-family test in `discovery/candidates.py` (`sister_paralogs`,
`sister_margin=0.10` per D7, `exclude_sister_family`). A first draft that
skipped sister-family members *by gene symbol* was thrown away — that is
exactly the assumption D14 forbids, and it would do nothing for the unnamed
RyR-sized loci S0 found. The shipped test is positional: identity to the
nearest labelled sister bait vs identity to the nearest known paralog, margin
recorded on every candidate, gene symbol never consulted, length band
contributing only the size component. All six RyRs now cap at 39.

**Numbers.** Recall 24/25 (96 %) — ITPR1 8/8, ITPR2 7/7, ITPR3 7/7,
invertebrate grade 2/3. Specificity 31/31 (100 %). Bait margin 31/31 correct
under both identity metrics; true-ITPR margins +0.065…+0.874, RyR margins
−0.868…−0.536, no overlap.

**Two caveats the benchmark surfaced, both logged as Emergent.**
(1) *Dictyostelium* `iplA` has a margin of +0.065 — inside the D7 10 %
no-call band. The bait margin is a metazoan instrument; the deep branches
need profile-based assignment. (2) The only recall miss, fly `Itpr` at 35,
fails by 0.008: its nearest non-known sibling is worm `itr-1` at 0.342
against a 0.35 breadth threshold. Under `covered_only=True` identity the same
pair scores 0.420 and breadth fires. That metric change would also lift every
RyR from 0.107 to 0.256 identity — into the twilight zone, +20 each — which
is survivable only now the sister test exists. It is a re-weighting, so it is
logged for S6 rather than done silently.

**Files.** `scripts/s1_toolchain.py`, `s1_panels.py`, `s1_lib.py`,
`s1_benchmark.py`, `s1_report.py`; `results/toolchain_manifest.txt`;
`results/benchmark_controls/` (report.md + 9 tables/artefacts).

**Next session: S2** — uncapped InterPro enumeration of the family Pfams into
census v2, with a positive ITPR/RYR call on every record. Carry the RyR bait
panel into every analysis set or the sister test has nothing to measure
against.

---

## 2026-08-19 — Review figures (user request, outside the ledger)

**Task.** Illustrate `docs/ip3r_review_2026.md` with structural diagrams,
sequence alignments and supporting figures. Not a ledger row; recorded here
and against S0, whose output the review is.

**Twelve figures, and the rule they follow (D19).** Nothing is drawn by hand
and nothing is drawn from a live query. Four data scripts fetch once and
commit small tables — `s0_figdata_domains.py` (InterPro domain coordinates),
`s0_figdata_structure.py` (PDB 6DQN), `s0_figdata_align.py` (two MAFFT
alignments of the committed control panels), `s0_figdata_curated.py` (the
three literature tables) — and five figure modules plus a driver render from
those tables and nothing else. Each figure carries a corner tag: *measured*,
*computed*, *schematic* or *curated*.

**The structure was measured rather than described.** PDB 6DQN (human
IP₃R3, IP₃-bound, 3.33 Å — the review's [R24]) was reduced to one
C4-symmetric subunit's Cα trace (0.06 Å RMSD between subunits, so the other
three are drawn by rotation), a pore-radius profile and a measurements file.
The 14 MB mmCIF went to a work directory outside the repo; the script refuses
a `--work-dir` inside it.

Three results came out of that which the text did not have:

1. **The pore's two constrictions were recovered blind.** Nothing in the
   calculation knows the sequence. The narrowest luminal point lands on
   Asn2472/Gly2473 — the residue before the **GGGVGD** selectivity-filter
   motif and its first glycine — and the cytosolic constriction on
   **Phe2513 and Ile2517**, one helical turn apart, the arrangement reported
   for the IP₃R1 gate. The motif check is now an assertion in the script: it
   exits non-zero if the luminal minimum is not on the motif.
2. **"Roughly 100 Å from the gate" is the axial component.** Measured:
   103 Å along the pore axis, 62 Å out from it, **120 Å through space**. §2.4
   now says both.
3. ***Dictyostelium* iplA carries none of PF08709, PF02815 or PF00520.** A
   characterised IP₃ receptor, in this project's own positive panel, that the
   family's *defining* signature does not find. Added to §7.5 and Q4 and
   logged as Emergent for S2/S3.

**Build.** `s0_review_build.py` now numbers figures in order of first
appearance and resolves `{fig:<slug>}` → `Figure N`, so no source file
carries a figure number — the same treatment the citations already get. A
`{fig:}` reference with no figure exits 3; a placed figure whose png/pdf is
missing exits 4. Both guards were tested by breaking them. For LaTeX the png
is swapped for the vector pdf inside a `center` block, so figures stay where
the text puts them. Figures are drawn at `W_REVIEW` (6.38 in, the review's
own text block) rather than `figstyle.W_FULL` (6.7 in, the manuscript's),
because a 5 % placement scale is the exact failure figstyle exists to prevent.

**One correction to the review's content.** SCA29 was going to be drawn as
dominant-negative alongside the Gillespie and *ITPR3* variants. §9.1's
sources do not establish that — they establish dominant and missense, and
that one variant is a *gain* of function. It is drawn as **unresolved**, and
the disease figure distinguishes three resolutions of evidence (`point` /
`domain` / `gene`), which makes visible that only two of nine entries have a
residue behind them.

**Files.** `scripts/s0_figdata_{domains,structure,align,curated}.py`,
`s0_fig_lib.py`, `s0_figs_{structure,sequence,genomics,concepts,clinical}.py`,
`s0_review_figures.py`; `results/s0_baseline/review_figures/` (11 tables);
`docs/figures/` (12 png + 12 pdf); 10 of the 13 review sections edited.
Review: 24 → 32 pages.

**Next session: S2**, unchanged — but S2 now inherits a harder requirement
from the *Dictyostelium* result: the InterPro enumeration cannot be the whole
census, because the defining signature demonstrably misses a real member.


---

## 2026-09-03 — S2: uncapped InterPro enumeration (census v2)

**Task.** Enumerate the family's Pfam signatures to exhaustion, make a
positive ITPR/RYR call on every record, and deliver census v2 with a delta
against what the app's own searches found.

**Result.** 15,421 proteins across 1,488 taxa — ITPR 6,433, RYR 6,807,
unassigned 2,181 (14.1 %). All three seed signatures walked their cursor
chains to the end (63 + 67 + 60 pages, no restarts). Raw pages, both parsed
intermediates, the 47 MB seeded-space FASTA and the full 1,251-sequence
representative set are under the data root; the repo gets the census, the
audit tables and a 48-sequence per-phylum core panel.

**The API is not what its own count says it is.** InterPro advertises a
`count` per signature and it is wrong in both directions — PF08709 +168,
PF01365 +56, PF08454 −176 against what the same endpoint actually serves,
stable across re-queries. In all three cases the *served* set matches
UniProt's independent count to within two records, so the advertised number
is the wrong one. Both directions bite, and the second one bit here: the
first version of this task's completeness test used `served ≥ advertised`
and failed PF08454 on a walk that had run its cursor to the end. The test
is now cursor exhaustion with both counts recorded beside it (**D20**).

**And the documented host was down for the whole run.** `/interpro/api/`
answered 1 request in 12 while `/interpro/wwwapi/` — the host the InterPro
website itself calls, same payloads, same counts — answered 12 of 12.
`list_proteins_with_pfam` now tries both hosts before it sleeps, and got
resume support (`start_url`/`on_page`, `CursorExpired`) so a dead cursor
restarts a signature rather than silently resuming onto a shifted set.

**The call.** Positive architecture test (**D14b**): RYR on any RyR-specific
signature, ITPR on the *complete* five-signature architecture with none of
them. Absence counts as evidence only when the architecture is complete,
because a RyR annotated well enough to show all five shared signatures
would also show its own. Length is on every row and never decides.
Audited two ways: against 6,191 gene symbols the rule never sees (**0
disagreements**), and on sequence alone via D14's bait margin over the core
panel — 27/27 and 26/26 agreement wherever the margin decides, with all
four sign flips inside D7's no-call band, exactly where S1 predicted the
deep branches would land.

**Getting the architecture cheaply was the design move.** Per-protein
InterPro lookups would have been ~15,000 requests against an API answering
one in twelve. UniProt returns the same Pfam matches as a *field*
(`xref_pfam`) — one streamed query gives the whole search space with its
architecture, lineage and length. That also made the two databases
independent views to reconcile rather than pool: they agree on 15,417 of
15,421, and the four InterPro-only rows are *Taenia solium* fragments with
newly-issued accessions where UniProt's cross-reference has not caught up.

**Three results the task was not asked for.**

*A single-signature census fails.* PF08709 is the IP3-binding core, the
signature that names the family, and building the census on it alone would
have missed 2,914 records — 758 of them called ITPR, across 385 taxa. Only
66.5 % of the space carries all three seeds. This closes the question the
review-figure session raised: *Dictyostelium* iplA carries PF01365 and
PF08454 and none of PF08709/PF02815/PF00520, and last session's note said
S2's enumeration would not return it. It does — the union recovers it — and
it is correctly `unassigned` at 2 of 5 signatures, which is the honest
outcome and the argument for S3's profiles (**D21**).

*The plant and fungal records are phylogenetically clean.* Every one of the
20 Viridiplantae ITPR calls is **Chlorophyta** (11 taxa, *Chlamydomonas
reinhardtii* with the complete architecture) and **Streptophyta has zero**
from 15 records in 13 taxa. Every fungal call sits in an early-diverging
phylum — Mucoromycota 15, Chytridiomycota 6, Basidiobolomycota 3,
Entomophthoromycota 1 — and Dikarya contributes no records to the search
space at all. That is the shape of a loss in the derived lineage of both
kingdoms. It is still a statement about what UniProt holds; S20/S23 make it
a statement about genomes.

*The size band caught what the architecture could not.* Seven records carry
the complete five-signature architecture in 1,528–1,993 aa — 700+ residues
short of the shortest real member — and none is flagged `Fragment`, because
a truncated gene model submitted as a whole protein is not marked as one.
The call on them is right and the records are wrong. That is the case for
keeping length as a column after it stopped being part of the call.

**Two figure fixes worth recording.** The first draft of `census_space` and
`census_growth` were stacked bars on a log axis, where segments do not add
up — the ITPR share looked like 80 % of a bar it was 45 % of. Both are now
size-on-log and composition-on-linear as separate panels. And
`census_growth` carries its own caveat in its title, because census v1 is
the app's search bundles rather than a family-wide harvest, so "100 % new"
outside the vertebrates is the searches' scope, not their failure.

**Files.** `scripts/s2_{lib,interpro,uniprot,call,sequences,verify,delta,
figures,report,run}.py`; `results/census_v2/` (14 tables, 5 figures,
`report.md`); `src/databases/interpro.py` (host failover, resume, backoff);
`src/utils/family.py` (`CENSUS_PFAM_IDS` is now the three enumeration seeds,
MIR excluded with the reason).

**Next session: S3** — profile-HMM sweep (`itpr.hmm` + `ryr.hmm`) over
vertebrate reference proteomes, jackhmmer to convergence, census v3. It
inherits two jobs from here: resolve the 2,181 unassigned records, and give
the deep branches the best-profile assignment the bait margin cannot.

---

## S3 — 2026-09-04 — profile-HMM sweep, census v3

**Task.** Build `itpr.hmm` and `ryr.hmm`, sweep the vertebrate reference
proteomes, iterate jackhmmer to convergence, merge into census v3 with an
ITPR-vs-RYR margin on every hit.

**Ran.** `s3_build_seed.py` (34 ITPR + 22 RyR seeds from the S2 census,
MAFFT L-INS-i, hmmbuild → 2,684 and 4,930 match states) · `s3_fetch_
proteomes.py` (763 vertebrate reference proteomes, 14,414,821 canonical
proteins, 6.97 G residues, 8.8 GB) · `s3_calibrate.py` (both profiles over
S2's archived seeded space) · `s3_run_sweep.py` (2 × hmmsearch, 3 ×
jackhmmer) · `s3_census_v3.py` · `s3_figures.py` · `s3_report.py`.

**The instrument was measured before it was used.** Both profiles were run
over S2's whole seeded space and the assignment scored against S2's
architecture call — different evidence entirely, one reading annotation and
one reading residues. **11,875 agree, 1 disagrees** (seeds excluded), and
the profiles resolve **2,314 of the 3,361 records the architecture rule
could not call**. S2's one non-architectural step, the low-confidence
gene-symbol fallback on 1,219 records, is overturned exactly **once**.

**The single disagreement is a real correction.** *Tieghemostelium lacteum*
A0A152A7I8, 2,845 aa, called RYR by S2 on **PF06459** — "Ryanodine Receptor
TM 4-6", which is the *pore*, the part both families share. `itpr.hmm` beats
`ryr.hmm` 303 to 133 bits. D14b wrote down that one of its RyR signatures
was only conditionally specific; this is that caveat firing, once, in
11,876 chances.

**The trap, and D22.** The first sweep called **14,981 vertebrate proteins
RYR**, 12,874 of them matching under a tenth of the profile and named
*Tnnc2*, *Rspry1*, *Cabp1*, *Rnf123*, *Ash2l*. Cause: `ryr.hmm` carries
SPRY, which `itpr.hmm` has nothing to match against, so every EF-hand and
SPRY protein in the proteome won by default. The gate is **200 match
states**, measured from this project's own S0 domain coordinates — the
shortest PF08709 it has ever measured is 200 aa, the longest SPRY 137 aa,
and the floor sits in that gap. It costs 86 of 12,020 architecture-called
records their profile verdict and keeps split gene models, whose 200–300
position band is almost entirely `Itpr1_0` / `Itpr2_1` / `Ryr3_0`.

**Completeness, measured the expensive way.** Of 2,787 census-v2 ITPR
records from a swept proteome that the sweep did not return, every one was
looked up in the 8.8 GB database itself: **zero were in it and missed**. All
2,787 are UniProtKB entries the taxon's reference proteome does not contain,
so they were never searched. And in the other direction the sweep adds
**618 proteins the InterPro census never returned** — all 209–942 aa,
fragmentary gene models a signature query cannot enumerate.

**Census v3: 16,039 records — ITPR 8,000, RYR 7,432, unassigned 605,
conflict 2**, every row carrying both instruments' verdicts and which of
them spoke (**D23**). 1,578 records change call from v2, all of them
`unassigned` → a call.

**Three things the session had to fix in its own method.**

*MAFFT is not reproducible with `--thread -1`* (**D24**). Rebuilding the
profiles after an unrelated edit changed them: the same 22 RyR seeds gave
8,510 and 8,468 columns on consecutive runs, and profiles of 4,933 and
4,908 match states. Found by noticing the match-state count had moved. Now
pinned to `--thread 1` (byte-identical across three runs) with a SHA-256 on
every seed set, alignment and profile. Both hmmsearch stages were re-run
against the rebuilt profiles; the calibration numbers moved by ≤ 3 records.

*D10's kill criterion measured the wrong thing* (**D10a**). K1 was coded as
a flat 5 % ceiling on sister-family content and fired on every run at round
1, because a single ITPR1 sequence at E ≤ 1e-5 already returns **32 % RyRs
before any iteration**. That is shared ancestry, not contamination. K1 now
takes round 1 as the baseline and kills on a **rise** of more than 10
points; the accepted runs are flat to falling (−2.2 to 0.0).

*K3 punished a converged run.* It fired on the ceiling without consulting
jackhmmer's own convergence verdict, so the run that converged **on** round
10 was marked killed. `evaluate()` now takes `converged`.

**jackhmmer.** `itpr_fly` **converged in 10 rounds** (5914 → 695 → 71 → 58
→ 470 → 24 → 1 → 2 → 0 → 0), D10 clean. `itpr1_human` reached the ceiling
at an asymptote of 3–12 new targets a round (K3, 9 of 10 rounds accepted).
The two final models differ by **single targets in every category** despite
seeds ~600 My apart — a completeness statement that does not rest on either
converging. A third of both models is module-only matches, which is why they
asymptote rather than reach zero. The third seed, *Acanthamoeba* L8GF85,
was **still running at session end** and is left to finish: it passes
**17.0 % of the database through HMMER's MSV filter against an expected
2.0 %**, so 2.46 M sequences reach the expensive stages of every round and
one round costs more than either other seed's entire run. Its log is
archived; `s3_run_sweep.py --stage jackhmmer --parse-only --seed-tag
itpr_acanthamoeba` folds it in without searching again. Nothing in census v3
depends on it.

**Files.** `scripts/s3_{hmm_lib,seed_spec,build_seed,assign,kill,fetch_
proteomes,calibrate,run_sweep,census_v3,figures,report}.py`;
`results/hmm_sweep/` (both profiles, seed manifest, calibration, sweep
assignments); `results/census_v3/` (census, novel hits, conflicts, the
two-directional completeness tables, convergence, model composition, 4
figures, `report.md`).

**Next session: S4** — the declared assembly manifest, which is the
denominator every later absence claim is measured against. S3 hands it a
concrete starting list: 15 swept taxa with no ITPR record at all, whose
gene sets are thin enough that annotation is the likelier explanation.

---

## 2026-09-04 (session 5, part 1) — S3 addendum: third jackhmmer seed folded in

**What ran.** The *Acanthamoeba* jackhmmer run left going at the end of the
S3 session finished (10 rounds, 27,242 s = 7.6 h, not converged) and had
been re-parsed from its archived log. Folded in: `s3_census_v3.py` (221 s,
now 3/3 seeds), `s3_figures.py`, `s3_report.py`.

**Result.** D10 verdict `killed (K3)` — the ceiling without convergence.
The three final models were then compared by composition, and the
family-called categories are *identical*: ITPR/architecture 1,656 and
ITPR/partial 952 in all three; RYR/architecture 1,565–1,569. The models
differ only outside the family: `itpr_acanthamoeba` rests on 26,266 targets
against `itpr1_human`'s 7,222, **18,670 of them with no sweep hit from
either profile** (7 in each of the other runs).

**New decision D10b — K1 is blind to off-family drift.** The off-family
accretion *diluted* the sister-family share: K1's statistic **fell 27
points** while the model grew 5×. A maximally divergent run reads as
maximally clean on the rule written to catch divergence. K3 is the only
rule that fires on it. Recorded in the roadmap with the instruction not to
replace K3 with a convergence test in S20.

**Code.** `scripts/s3_report.py` had grown to 530 lines with the rewrite, so
section 4 moved to **`scripts/s3_report_d10.py`** (367 + 207 lines). The
two stale hand-written sentences it contained ("about a third is
module-only", "the two seeds differ by single targets") are now computed
per evidence class from `jackhmmer_model_composition.tsv` — with three
seeds the first was true of two runs and the second of none. Added a
wall-clock column read from the per-seed sweep-stats JSONs, and rescaled
figure panel *b* to the data (the third trace was clipped off-canvas by a
y-limit chosen for the two-seed case).

**Next: S4** — the declared assembly manifest.

---

## 2026-09-04 (session 5, part 2) — S4: the genome manifest

**Task.** S4 — declare the assembly manifest that is the denominator for
every later absence claim, and build the download tooling.

**Scope decision (put to the user).** Order representatives are fixed by the
roadmap; the choice was which margin set joins them. Chosen: the **full
margin union**, so the bird question is answered per species rather than by
sample. Alternatives offered were dropping the 100 fragment-only species, or
capping Aves at a stratified ~25.

**What ran.** `datasets summary genome taxon 7742 --reference` → 6,205
vertebrate reference assemblies; taxonomy for 6,185 taxids in 400-taxid
chunks; both dumps archived under `<data_root>/raw_api/ncbi_datasets/` so
the build reruns offline. Then `s4_build_manifest.py` (32 s cold, 2 s cached).

**Result: 309 genomes**, 161 orders ∪ 169 margin species, 552.5 Gbp —
≈553 GB of FASTA, ≈166 GB as zip, against 1,434 GB free on the drive.

**The margin species are derived, not typed.** This is the one real
departure from the PIEZO port, whose margin list was 28 hand-written species.
Here four rules with stated thresholds run over the committed census tables
(`s4_manifest_lib.derive_margins`): `zero_hit_proteome` 15, `missing_paralog`
(<3 paralogs) 101, `fragment_only` (longest record < `family.MIN_LENGTH_AA`)
100, `anchor` 5. 61 genomes carry two reasons, 11 carry three. The
denominator therefore rebuilds from the census and cannot drift from it.

**130 of the 169 margin species are birds** — S3's annotation-depth result
(Aves median 13,894 proteins vs Mammalia 34,127) reappearing as a scope
requirement. D9 is live in the manifest too: 168 RefSeq / 141 GenBank, and
**33 genomes carry no gene set at all**.

**`fetch_genomes.py` tested end to end** on the three smallest genomes
(*Takifugu rubripes*, *Genypterus blacodes*, *Lepidogalaxias
salamandroides*, ~1.4 GB total, ~30 s each): md5 verification against NCBI's
own `md5sum.txt`, `.done` resume, `--verify`, fetch → search → purge, and a
failing `--search-cmd` correctly keeping its genome.

**Two bugs the testing found.**
1. Margin species were matched to an assembly by a **last-wins taxid
   lookup** while order reps used the rank function; a species with two
   reference assemblies entered the manifest twice under two accessions
   (310 rows / 309 taxids, and zebrafish lost its `anchor` reason). Both
   paths now use `rank_key`.
2. `--verify` **silently skipped** any checksummed file the install does not
   keep (`dataset_catalog.json`), so a deleted or truncated `.fna` would have
   verified clean. The install now rewrites `md5sum.txt` to exactly the
   files it kept, and a missing name is a hard failure — proved by deleting
   `sequence_report.jsonl` and watching `--verify` exit 1.

**Files.** `scripts/s4_manifest_lib.py`, `scripts/s4_build_manifest.py`,
`scripts/s4_notes.py`, `scripts/fetch_genomes.py`;
`results/genome_manifest.tsv` (309 rows), `results/genome_manifest_notes.md`.

**Next session: S5** — tblastn + miniprot over these 309 genomes, producing
the found / lost / assembly-gap ledger. Start with the 130 bird margin
species; the three >10 Gbp assemblies should be scheduled deliberately.

---

## 2026-09-04 — S5a: the genomic sweep's instrument, built and calibrated

**Task.** S5 (genomic tblastn + miniprot sweep → per-genome ledger → census
v4). Scope is S4's 309 genomes / 552.5 Gbp, so it was split per the session
protocol: **S5a** builds and calibrates the instrument and pilots it,
**S5b** runs the full sweep.

**The bait panel** (`s5_bait_spec.py` + `s5_build_baits.py` +
`s5_bait_screen.py`). 38 baits — 30 ITPR + 8 RyR control, 121,294 residues —
**derived from `census_v3.tsv` by seven enforced rules**, not hand-listed, in
the shape S4 used for its margin species. 25 reused S3 seeds + 13 additions
(a Passeriform, a percomorph and a squamate per paralog; the bands that
dominate the scope take two baits each). D5's chimera screen is this
project's own instrument: both S3 profiles under S3's relative margin, plus a
profile-envelope test. Six slots could not be filled from the whole census —
no labelled full-length record exists for ITPR1/ITPR2 in chondrichthyans,
ITPR2 in the coelacanth grade, or any paralog in cyclostomes.

**The screen has its own negative control**, run on every build: three
synthetic failures built from real panel baits, each aimed at the rule
responsible. It earned its place by *failing on first run* — it asked the
envelope to reject a truncation, which it cannot and should not, since a
truncated protein aligns 100 % of *itself* in one clean segment. Truncation
is the length band's job; the envelope's job is fusion.

**Three calibrations a port could not supply.**

1. **miniprot's `-G`** (`s5_intron_calibration.py`). A too-small max-intron
   does not lose a gene, it *splits* one, and a split ITPR reads out of the
   ledger as `fragment`. Measured across 11 Ensembl species: **0/28 ITPR
   genes have an intron over the 200 kb default** (widest 152,216 bp, human
   ITPR1); the RyR control exceeds it twice (RYR2 227,927 bp). The rule
   scales from the widest measured intron, deliberately *not* from the
   intron-per-Gbp ratio — that statistic is largest in the smallest genomes
   and extrapolates to a 7.4 Mbp `-G` on the lungfish.
2. **The unlabelled baits were competing as a fourth paralog.** Rescue
   attribution ranked `vertebrate_basal` (gar/chimaera/lamprey — evidence
   that a region is *an* ITPR, not a hypothesis about which) against the real
   paralogs. In *Todus* an ITPR1 trace at a true 0.235 margin over ITPR2 was
   reported at 0.027. Three cells moved from `tblastn_trace_ambiguous` to
   `tblastn_trace` once the ranking was restricted to canonical paralogs.
3. **`ATTRIBUTION_REL_MARGIN` does not transfer** (new **D25**). Inherited as
   the PIEZO port's 1.5x ratio restated as 0.333, from a family 40–50 %
   identical where these are 61–68 %. At complete loci whose paralog identity
   the annotation independently establishes, **7/9 sit below it** (0.225–
   0.347, median 0.281). A complete locus is an upper bound on a rescue
   fragment, so under 0.333 no absence could ever be attributed → **0.22**,
   limitation stated in the code, per-region margins recorded so S5b can
   re-cut without re-running.

**The pilot.** 6 genomes spanning the classification space (2 teleosts, 3
birds incl. the `zero_hit_proteome` and an unannotated assembly, 1 mammal);
24 cells, 43 loci, **0 RyR control failures**, 0 `absent` calls. The RyR
control paid immediately: in *Todus mexicanus* it found **one** RyR locus
where a bird has three, and unannotated — what is missing there is
contiguity, not genes.

**D14 at genome scale: at all 43 loci only one family's baits aligned at
all.** The ITPR and RyR panels never contested a locus — sharper than the
protein level, where S1 had all six RyR decoys promoted at 45 before the
sister test existed.

**D4's contiguity bar**, set at the median measured ITPR span (142,212 bp)
rather than an invented number: **120/309 manifest genomes fail it — 66 % of
Aves vs 11 % of Actinopteri, 68 % of margin species vs 12 % of order reps.**
And below the bar the pilot recovers ITPR3 (82 kb span) 2/3 while ITPR1
(186 kb) and ITPR2 (231 kb) are 0/3 — the detection bias points the same way
as the loss signal.

**Files.** `scripts/s5_bait_spec.py`, `s5_build_baits.py`,
`s5_bait_screen.py`, `s5_intron_calibration.py`, `s5_calibration.py`,
`s5_sweep_lib.py`, `s5_genome_io.py`, `s5_classify.py`, `s5_rescue.py`,
`s5_run_sweep.py`, `s5_ledger.py`, `s5_calibrate_margin.py`, `s5_budget.py`,
`s5_report.py`; `results/s5_baits/`, `results/genome_ledger/`.

**Next session: S5b** — the full sweep. Budget measured, not estimated: 303
genomes, 547 Gbp, **6.7 h compute** (miniprot 36 s/Gbp at 8 threads), 3.4 h
wall at 2 shards, 164 GB download; 547 GB of FASTA fits in 1,427 GB free, so
keep it for S8/S10. Run the two giants (*Protopterus* 40 Gbp, *Lissotriton*
23 Gbp) separately on the chunked path. Re-calibrate the attribution margin
on rescue regions of known identity, and test the span-bias result at 309.

---

## 2026-09-04 — S5b: the full sweep (running), and hardening it while it runs

**Task.** S5b — the 309-genome sweep → ledger + census v4. User decided to
keep every genome on disk (547 GB against 1,427 GB free) rather than
`--delete-after`, since S8 synteny and S10 validation both want the loci back.

**Running.** Six shards (`s5_run_sweep.py --shard i/6 --threads 3`) plus a
supervisor (`<data_root>/genome_sweep/_overnight.sh`) that logs progress every
10 min and sweeps the two giants on the chunked path once the shards drain.
All detached (PPID 1). Zero failures through 77/307.

**The budget was wrong, and the error was an omission.** S5a projected 3.4 h
wall from compute rates. Measured on the running sweep: **the sweep is
download-bound.** NCBI delivers ~5.5 MB/s in aggregate and that barely moves
with shard count (2 shards ~3.4, 6 shards ~5.5), so it is a bandwidth ceiling
rather than a concurrency problem — 164 GB is ~8 h of transfer against ~2.5 h
of compute. Going 2 → 6 shards did help (42 genomes/hour measured on a settled
window, ~6.4 h to go) but nothing like linearly. miniprot itself is **faster**
than budgeted: 15 s/Gbp at 3 threads, not 36. `s5_budget.py` now carries a
`DOWNLOAD_MB_PER_S` term and reports which side binds.

**Census v4** (`s5_census_v4.py`). Merges the sweep's gene models into v3,
scoring each against `itpr.hmm`/`ryr.hmm` so a v4 row rests on the instrument
that called v3 (D23). It records **how** a database holds a locus rather than
only whether: `genome_only`, `annotated_unnamed`, `annotated_other_paralog`,
`annotated_named`. The PIEZO port filtered on cell status alone and would have
dropped the most useful population here — an annotated gene carrying no family
name, like *Takifugu*'s second ITPR1 3R duplicate beside `itpr1b` as
`LOC101074739`. The disagreement label is named for the observation, not a
verdict, and `sibling_locus_for_annot_paralog` records what makes a case
coherent (the genome also carrying a separate locus for the paralog the
annotation names). First case found: *Nibea albiflora*, a locus annotated
`ITPR2` that the ITPR3 cell claims at 100 % coverage and paralog margin 0.306,
in a genome that carries a separate ITPR2 locus.

**The fragment-level margin calibration S5a deferred is done.** Rescue regions
are now a committed table, and those overlapping a gene the assembly names for
a paralog carry an identity independent of the bait scores. 11 such regions so
far: the attribution agrees **11/11**, over margins 0.227–0.370. **0/11 fall
below the 0.22 threshold S5a set; 7/11 would have fallen below the inherited
0.333** and been reported ambiguous. The recalibration is vindicated on the
fragments it was actually for. (Also fixed a wording bug that would have
printed whatever constant is live as "the inherited threshold".)

**Hardening — the chunked miniprot path, before the giants run on it.**
`s5_test_chunked.py` runs a genome both ways and requires the same loci, baits,
coverage and cell statuses. On *Takifugu*: 10 loci and all 4 cell statuses
identical, worst boundary difference 9 bp (tolerance 50, justified — miniprot
extends a terminal exon from flanking context and there is less at a chunk
edge, while a real offset bug moves loci by megabases). Three findings:
1. **The mRNA-ID collision is real** — 167 rows carry 60 distinct raw IDs over
   6 chunks. The dedup guard in `parse_miniprot_gff` is load-bearing.
2. **The test's own memory measurement was wrong**: `ru_maxrss` is a
   high-water mark, so per-phase `after − before` reported 0.00 GB for the
   chunked run whatever it used.
3. **`CHUNK_BP` was too large.** Measured, miniprot wants **6.3–8.0 GB RSS per
   Gbp**, putting the ported 4 Gbp chunk at ~25 GB on a 34 GB machine → now
   2.5 Gbp (~20 GB). Extra chunks are nearly free: a locus never spans a chunk
   boundary.

**Hardening — resume correctness.** The driver reuses any non-empty
`miniprot.gff`, and miniprot wrote it directly, so an interrupted run left a
truncated file the next run accepted as complete — undetectable afterwards,
because a GFF cut at a line boundary parses perfectly. **Audited the live run:
74 completed GFFs all end on a complete record, and 0 genomes were marked done
after the 6-shard restart while carrying a GFF from before it**, so the earlier
shard kills corrupted nothing (they landed during download, not alignment).
`run_miniprot` and the chunk concatenation now write atomically via a
`.partial` rename. `s5_test_resume.py` proves both properties: a failed run
leaves neither file, and a re-run reproduces the summary **including through
miniprot** — D24's question asked of a second tool, and miniprot is
reproducible at a fixed thread count where MAFFT is not at `--thread -1`.

**Four ledger figures** (`s5_figures.py`), from committed tables only. The
contiguity panel bins rather than smooths: a sliding window straddling D4's bar
reports a contiguity no genome in it has and draws the curve through the very
threshold the panel exists to show. `figstyle.STATUS` gained the `fragment` key
it lacked.

**Next.** Finish the sweep, then `s5_ledger.py` → `s5_calibrate_margin.py` →
`s5_census_v4.py` → `s5_budget.py` → `s5_figures.py` → `s5_report.py`. Decide
the attribution margin last, on the full fragment distribution — changing it
re-derives the rescue attribution only, not the sweep. `s5_report.py` still
carries S5a's pilot framing and needs reworking for 309 genomes.

### S5b addendum — the report generator, and a control that changed the answer

Drafted the full-sweep report as `s5_report.py` (instrument sections) +
`s5_report_results.py` (results), following the `s3_report.py` /
`s3_report_d10.py` split so both stay under 500 lines and share the caller's
table loader.

**The design constraint: the generator must test S5a's pilot claims, not
narrate them.** Six-genome findings are hypotheses at 309. Each results
section computes its statistic, compares it against the pilot's value recorded
in `PILOT`, and renders `confirmed` / `not confirmed` / `contradicted` /
`underpowered` from the comparison, printing the pilot's number beside the new
one. A generator written around the pilot's answers would print them whatever
the sweep found.

Exercised against the partial sweep (81 genomes), it already earned itself
three times:

1. **Caught a latent break**: the S5a section still read
   `n_contested_below_threshold`, a JSON key renamed when the calibration
   started distinguishing inherited from in-force. It would have failed at the
   end of the sweep instead of now.
2. **The span-bias hypothesis renders `Confirmed` but much weaker than the
   pilot** — below the bar ITPR3 67 % vs ITPR1 52 % (pilot: 67 % vs 0 %). The
   verdict is right and the magnitude is honest, which is what printing the
   pilot's number beside it is for.
3. **A control that reversed a finding.** Per-paralog annotation completeness
   looked dramatic uncontrolled — ITPR3 correctly named at 54 % of found loci
   against ITPR2 at 14 %, a 39-point gap. Holding contiguity constant (only
   assemblies whose contigs can carry the gene) it **vanishes**: 38 % / 44 % /
   40 %, within 6 points. The apparent annotation-quality difference was
   assembly quality. The section now reports both blocks, says which one the
   verdict reads, and warns when there are too few above-bar loci to control.

Also of note in the partial ledger: **0 cells are `absent` at 81 genomes**, and
D14 still holds absolutely — at all 566 loci only one family's baits aligned.

### S5b complete — 309/309, 0 failures

**The sweep.** 309 genomes, 1,236 cells, 2,144 loci, zero failures. Six shards
plus a supervisor that picked up the two giants once they drained; both giants
completed on the chunked path (*Protopterus* 40 Gbp / 25 chunks / 1.8 h,
*Lissotriton* 23 Gbp / 14 chunks / 2.5 h). Wall clock ~19 h, download-bound
throughout.

**Two instrument fixes the run itself surfaced.**

1. **`MIN_LOCUS_IDENTITY` = 0.40.** The giants came back with ITPR3 at 10 loci
   (*Lissotriton*) and 8 (*Protopterus*). Inspection: exactly one locus per
   paralog was real (coverage 1.0, identity 0.79–0.91, annotated), and every
   extra sat at 10–27 % coverage and **23–34 % identity**, often overlapping
   unrelated genes (PSME4, TRAPPC9, CDH20, EXOC1) — the shared channel module,
   chained across megabases by a 2 Mbp `-G`. 51 % of loci in assemblies over
   5 Gbp were like this, against 9 % elsewhere. The floor is measured from a
   clean gap: the **571 loci whose paralog an annotation independently
   confirms have a minimum identity of 0.759**, and any floor in 0.35–0.50
   drops the same 22 junk loci and zero real ones. No *status* was ever wrong
   (the real gene always scored best) but `n_loci` was inflated, and copy
   number is a result. All 309 re-derived from cached GFFs after the fix.
2. **Atomic miniprot output.** The driver reuses any non-empty `miniprot.gff`,
   so a run killed mid-write left a truncated file the next run accepted as
   complete — undetectable afterwards, since a GFF cut at a line boundary
   parses perfectly. Audited before fixing: 74 GFFs all ended on a complete
   record and 0 genomes were marked done carrying a pre-restart GFF, so the
   earlier shard kills corrupted nothing.

**Results.**
- **RyR positive control fired in 309/309.** Nothing excluded on control grounds.
- **D14 held absolutely**: across all 2,144 loci only one family's baits
  aligned. The two panels never contested a locus.
- **4 `absent` cells, all cyclostome**: lamprey and hagfish each carry ITPR1
  and lack ITPR2/ITPR3, both above the contiguity bar with a firing control.
- **D4's bar**: 120/309 genomes cannot hold the gene on one contig. Recovery
  98–99 % above it, 57–70 % below.
- **Span bias confirmed at 13 points** (ITPR3 70 %, ITPR1 61 %, ITPR2 57 %
  below the bar). It read 15 → 8 → 12.5 → 13 as the sample composition shifted
  during the sweep, which is exactly why the report generator was built to
  render the verdict from the final table rather than narrate the pilot.
- **Attribution margin validated on fragments**: 53 rescue regions of
  annotation-established identity agree 53/53, margins 0.227–0.445; 0 below
  the 0.22 in force, **31 of 53 below the inherited 0.333**.
- **Annotation quality differs by paralog** — contiguity-controlled over 487
  loci, ITPR3 88 % correctly named vs ITPR1 65 %. Uncontrolled the same data
  said something different and wrong; the comparison reversed twice as n grew.
- **Census v4 = 17,097 records**: +1,058 gene models, of which **318 ITPR
  models exist only as DNA** and 167 more sit in an annotated gene with no
  family name.

**Next: S20**, the non-vertebrate sweep — now the topmost unblocked row.

---

## 2026-09-05 — S20: the non-vertebrate sweep

**Task.** The family's true range across eukaryotic reference proteomes, and
whether the land-plant / Dikarya absence is a genome fact or a database fact.

**Instrument.** S3's two profiles reused unchanged — building a new one would
have made the vertebrate and non-vertebrate numbers incomparable, which is the
whole point of the exercise. Thirteen new modules (`scripts/s20_*.py`), all
under 500 lines.

**Scope.** Six groups: the four non-vertebrate eukaryote groups that partition
Eukaryota with S3's `vertebrata`, plus all 634 archaeal reference proteomes and
a genus-stratified bacterial sample of 3,537 drawn from 17,981. Sampling rule
stated and deliberately generous to the hypothesis it tests — one proteome per
genus, the one with the *most* proteins, because a bigger proteome is a more
sensitive place to find a homolog. **6,928 proteomes, 63,144,898 proteins,
24.93 G residues.** 23 proteomes UniProt lists are not published in the release
FTP tree and 404 permanently; recorded as exclusions rather than aborting the
fetch.

**One search, two sensitivities.** Every `hmmsearch` ran at `-E 10` and the
primary E ≤ 1e-5 call was taken by filtering the same domtblout, since `-E` is
a reporting threshold and does not touch the acceleration filters. The relaxed
claim is therefore a superset of the strict one from the same search, not a
second experiment that might have differed some other way.

**Three silent bugs, all found by checking rather than by failure.**
1. **The presence table was being overwritten with a fraction of its own
   denominator.** The sweeps run one group at a time and `presence()` was
   computing over only the groups of that invocation. Now always over every
   group with an assignment table.
2. **Hits were attributed to proteomes by taxid.** 3 protist, 18 plant and 27
   fungal taxids carry *two* reference proteomes each, so a taxid key credits
   both with either one's hits. Attribution is now measured, by a header-only
   pass over each proteome file, cached under the data root.
3. **UniProt's `tax_id:` search is hierarchical.** `tax_id:3702` matches
   *Arabidopsis* and every strain under it, so a 100-term OR batch matched more
   than 100 records, the page capped at `size`, and requested taxids fell off
   as blank rows — *Arabidopsis* and *Chlamydomonas* among them, which are not
   taxa this task can afford to lose. Switched to the exact `taxonIds/`
   endpoint, whose page cap (25, measured) now sits above the batch size, so a
   short page can only mean genuine absence. The fetch also verifies that
   nothing came back that was not asked for.

**A fourth thing that would have become a result.** The sweep reports 142 RYR
calls in fungi, green algae and protists. They are 500–2,000 aa proteins
matching 4–9 % of a 4,930-state model — the shared module, not the gene. The
report now breaks every call out by model coverage: 63 % of ITPR calls are
architecture-level against 1 % of RYR calls.

**Results.**
- **662/6,928 proteomes carry an ITPR call.** Present across Metazoa, SAR,
  Discoba, Amoebozoa, Haptophyta, Chlorophyta and five early fungal phyla.
- **Streptophyta 0/384** (16.3 M proteins) — S2's seeded-space finding
  CONFIRMED with the seeding filter removed. Chlorophyta 15/48.
- **Dikarya 0/1,353** — CONFIRMED. But the fungal losses are patchy, not
  basal: Glomeromycota 0/27, Mortierellomycota 0/18, Kickxellomycota 0/35 and
  Microsporidia 0/29 are also empty.
- **Archaea 0/634, Bacteria 0/3,537.**
- **The negatives have positive controls inside the same search.** At E ≤ 10
  with the four family Pfam models: PF08709 returns 0 substantial matches in
  land plants against 26 in Chlorophyta, and 0 in Dikarya against 16 in
  Mucoromycota — while PF02815 (MIR, which every eukaryote carries on other
  proteins) returns 633 and 4,376 in those same genomes. The instrument works
  there; it finds everything except the receptor.
- **All 99 plant and fungal records chased individually: 47 `real_gene`, 52
  `fragment`, zero contaminants, zero without genome backing.** Cross-kingdom
  identity to the nearest relative runs 19.9–45.8 %, median 24.1 % — the deep
  homology range, nowhere near the 95 % contamination call.
- **D14 outside the vertebrates: 0 of 704** targets scored by both profiles
  above the bit floor fall inside the no-call band.
- **The two architecture-level RYR records outside Metazoa are genuine**, in
  *Capsaspora owczarzaki* (a `ryr.hmm` seed, so circular, and flagged as such)
  and *Salpingoeca rosetta* (not a seed) — the two closest unicellular
  relatives of animals carry both families at full length.
- **jackhmmer, viridiplantae: converged in 5 rounds, D10 clean**, zero
  sister-family content in every round. All 70 targets in the converged model
  are Chlorophyta, including all 6 that only iteration found — so the
  land-plant absence is not a sensitivity artefact.
- **Census v5: 17,882 records** (+785 from this sweep), 8,807 ITPR across 1,401
  taxa, lineage columns on every row including S5b's 488 genomic gene models,
  joined to S4's manifest. **No v4 call is overturned and the 2 conflicts stand**:
  the first version of the merge let this sweep break them, but S20 scores with
  the same two profiles that produced one side of each disagreement, so it is
  the same instrument re-scoring the same protein in a different database, not
  a third opinion.

- **jackhmmer, fungi: ran the ceiling, D10's K2 fired at round 3.** The
  included set grew **34.3×**, 42 → 1,442, and by round 10 rested on 6,302 —
  while the family content went 33 → 35. The growth is *entirely* off-family,
  so **K1 read 0.000 in every round**: this is S3's D10b in a new setting,
  caught by K2 rather than by the round ceiling. The accepted model (rounds
  1–2) reaches Ascomycota exactly twice, and both records are
  dolichyl-phosphate-mannose mannosyltransferases — the MIR-domain sharer S1's
  decoy panel was built around, not a receptor. Iteration-only targets are now
  named in `jackhmmer_iteration_only_s20.tsv` rather than counted, and the
  report only calls them known decoys when their own protein names say so.

- **jackhmmer, protista_other: 10 rounds, D10's K3.** The third distinct
  outcome from the same coded criterion, and the sharpest D10b case yet. The
  model grew 721 → 26,148 targets, but **never by more than 3.18× in one
  round, so K2 never fired**; and the sister share *fell* from 0.28 % to
  0.16 %, with `sister_rise` **negative in 6 of 10 rounds**, so K1 not only
  missed it but moved the wrong way. Only the ceiling caught it. What it
  accreted is the point: **11,344 Apicomplexa proteins**, a clade in which
  this task's own sweep called the family in 0/60 proteomes. D10's words for
  K3 are that no completeness claim may rest on the run, but `s3_kill` still
  reports its pre-ceiling rounds as `accepted` — correct for K1/K2, where the
  rounds before the drift are usable. The composition rows now carry the
  verdict and a `supports_completeness` flag, and the report presents a
  disowned run separately instead of mixing it into the completeness table.
- **jackhmmer, metazoa_nonvert: 10 rounds, K3 — and the fourth outcome is a
  K3 that means the opposite of the protists'.** Nothing drifted: growth never
  exceeded **1.21×**, `sister_rise` peaked at **+0.042** against K1's 0.10
  limit, and the finished model is **2,198/2,908 (76 %) records the profiles
  call ITPR or RYR**, its largest single contribution 451 arthropod RyRs. It
  simply had not finished — new targets were still trickling in at ~35–70 a
  round. Against the protist K3, whose model is **772/22,913 (3 %) family** and
  whose largest contribution is 11,344 apicomplexan proteins. **The rule
  returns the same verdict for both**, because K3 is about rounds rather than
  content — which is right for a ceiling, and a reason to read the composition
  beside the verdict rather than instead of it. The report now tabulates that
  contrast instead of pooling the two.
- **Family saturation, read off every run, killed or not.** The count of
  family members settles early while the model keeps growing: viridiplantae
  at round 2 (54), protista at round 3 (729 — exactly what the single pass
  found, after which the model grew 1,995 → 26,148 without adding one),
  fungi at round 4 (35), metazoa at round 5 (1,194). **In all four groups the
  saturated count is *exactly* the single pass's ITPR count** — 54, 729, 35,
  1,194 — so iteration finds no receptor the single pass missed anywhere in the
  non-vertebrate tree, and that statement does not depend on the two D10
  verdicts. What the extra rounds bought was 24,153 more non-family targets in
  the protists and 870 in the fungi.

**Ledger split (S20 → S20a / S20b), both now complete.** Everything in S20's completion criteria is delivered, and both convergence runs the negative claims rest on are done. The other two are cost-bound rather than unfinished-in-principle: the per-round cost is the *alignment*, not the search, and it does not fall with more cores — a round over a group with thousands of included targets runs ~45-65 min whatever thread count this machine gives it, so protista needs ~5 h and the invertebrates ~8 h. Both were launched here and are still iterating; S20b collects them.

**Next.** S23,
then S23, which S20 sharpens: the absences to take to genome level are
Streptophyta, Dikarya, Glomeromycota/Mortierellomycota, Apicomplexa (0/60) and
the *Cymbomonas* copy-number question.

---

## 2026-09-05 — S23a: the non-vertebrate sweep's instrument

**Ledger split.** S23 is a genome-scale sweep, so it was split on the S5a/S5b
precedent: **S23a** builds and measures the instrument (scope, calibration,
bait panel, copy-number classifier, pilot), **S23b** runs it. Downstream deps
(S7, S14a) repointed to S23b.

**The headline is a method finding, and it would have voided the task's own
results if it had gone unnoticed: S5's positive control does not transfer,
and neither did three of its thresholds.** New decision **D26**.

### The control (D26)

S5b can write "the RyR positive control fired in every one of the 309, so no
genome is excluded on control grounds" because every vertebrate has three
RyRs. S20 measured architecture-level RyR in **2 of 6,928** non-vertebrate
proteomes. Carrying that control into a land-plant or Dikarya genome would
have made *every negative claim in S23 unfalsifiable while looking
controlled*: a plant with no RyR locus is the correct answer, so its silence
says nothing about whether the search ran.

The replacement is the **MIR-domain sharer** (PF02815 — the protein
O-mannosyltransferases S1's decoy panel was built from and S20 already used
as its in-search control at proteome level, where PF08709 returns 0 matches in
land plants and PF02815 returns 633). It is inside the family's own signature
set, it is drawn **per control clade** (a chlorophyte mannosyltransferase is
not a control for *Arabidopsis*), and it doubles as the sharpest available
decoy for D14.

Two corollaries, each of which cost a rebuild:

- **A control's selection rule must not be tuned to what the control usually
  looks like.** The first build required 600 aa, which selects the full
  mannosyltransferases — and silently left five clades with no control at all,
  **including both apicomplexan classes**, i.e. it deleted the hardest
  negative claim rather than reporting it. The floor is now the PF02815
  model's own length.
- **A control has to be graded.** 22 of 28 control clades get a `strong`
  control (a named, full-length mannosyltransferase across ≥ 25 % of the
  clade's proteomes); 6 are `weak`, with the reason in the row. **The
  apicomplexan clades carry one PF02815 protein each across 60 swept
  proteomes** (3–4 % coverage), so "Apicomplexa 0/60" at genome level rests on
  a 233 aa control found once → emergent task.

**The screen rejected one control on its first run**, which is exactly what a
negative control is for: `A0A0R3TP59` from the tapeworm *Rodentolepis nana*,
which UniProt calls a "MIR domain-containing protein", scores **1,121 bits on
`ryr.hmm` against 179 on `itpr.hmm`**. A MIR protein the profiles call a
family member is a finding, not a control.

### Three thresholds that had to be stratified, all the same error

Each was first measured over the whole in-scope population, and that
population is 53 % Arthropoda. A threshold measured on the best-sampled clade
describes the sampling, not the family — the same correction S20 made when it
stratified its bacterial sample by genus.

| threshold | unstratified | stratified | what the unstratified value deleted |
|---|---|---|---|
| bait length band | 2,550–3,100 aa | **2,450–3,250 aa** | *Bodo saltans*, a 2,356–4,222 aa euglenozoan ITPR |
| architecture-exception floor | ~1,095 bits | **752 bits** | the same records, at 735–782 bits |
| contiguity bar / `-G` | one global bar | **per group** | see below |

**The contiguity bar is per group because the family's span varies ~100×
outside the vertebrates against 6.5× inside.** Measured from NCBI's own gene
annotations (deliberately not from miniprot, whose `-G` shapes the loci it
reports, so a measurement taken that way cannot falsify the setting it
calibrates) over 16 genes in 16 species across 14 bands: spans run **3,739 bp
(*Perkinsus*) to 324,840 bp (*Octopus*)**, and **metazoan genes are ~9× longer
than protist and fungal ones** (83 kb median vs 7–9 kb). S5's single
vertebrate bar was 142,212 bp and 120 of 309 genomes failed it — S5b's
headline caveat. Here only **21 of 194** fall below their group's bar, and the
genomes carrying this task's negative claims are judged against ~9 kb.

`-G` moves the other way and is floored at miniprot's own default: a `-G`
below an unmeasured species' largest intron **splits** its gene, and a split
ITPR reads out of a *copy-number* ledger as `fragment` — worse here than in
S5, because in this task the count is the result.

### The denominator

**194 genomes, 100.4 Gbp** (≈100 GB FASTA, ≈30 GB zip against 805 GB free),
from 24,596 NCBI eukaryote reference assemblies less 6,216 vertebrate ones.
Five rules under one principle — **sample most finely where the negative claim
is** — all derived from S20's committed presence table except the anchors,
which are hand-written deliberately because an anchor is a genome whose answer
is known from the literature.

G1 phylum_rep 66 · G2 class_rep 87 · G3 absence_clade 35 · G4 anchor 14 ·
G5 copy_number 36. Every one of the 35 absence clades has at least one genome.

Two bugs the build found in itself:

- **NCBI files a reference assembly under the *strain* taxid**, not the
  species, so a species-taxid lookup missed four of fourteen anchors —
  *Toxoplasma*, *Encephalitozoon*, *Batrachochytrium* and *Dictyostelium*,
  i.e. the apicomplexan, the microsporidian, the chytrid and the amoebozoan
  the negative claims are named after. Fixed by indexing every ancestor taxid.
- **The unfilled-slot table blamed the census for a threshold's work.** Its
  first version printed "no census record in this band at all" for seven
  bands, four of which have dozens of records that fail the *shape* rules. It
  now names the stage that lost each slot.

### Copy number, not paralog cells

`s23_classify.py` replaces S5's three named cells: loci won by the ITPR baits
are graded `full` / `fragment` / `scrap` and the genome's copy number is its
count of `full` ones, with `merge_split_loci` folding neighbouring same-strand
loci whose bait spans are *complementary* — a split gene inflates the very
number this task reports, which is the error S5 could tolerate (its best locus
always won) and this one cannot. Both counts are kept so the merge is
auditable.

The control verdict gained a case the pilot found immediately:
**`controlled_by_target`**. *Dictyostelium* recovers its iplA gene and matches
it to its own annotation, and the first version still reported it
**UNCONTROLLED**, because no MIR bait exists for Amoebozoa (controls are built
for the clades carrying a negative claim). A control exists to make a
*negative* interpretable; a genome where the family itself was found has
already proved the search reached it. `no_control_bait` is likewise
distinguished from `uncontrolled` — a silent control and an absent one are
different facts.

### Regression check

`s5_bait_screen.py` was parameterised by spec module so one screen serves both
panels (and S5's own negative controls still run on every S23 build).
Re-running `s5_build_baits.py` reproduces the committed S5 panel **byte for
byte apart from its build timestamp**, which was reverted.

### What the pilot found — two more ported thresholds, and a rule I failed to port

The pilot's job is to break the instrument before 100 GB is spent on it, and
it did, twice.

**1. A positive control failed, and the first diagnosis was wrong.**
*Chlamydomonas reinhardtii* carries a documented complete 5/5 receptor
(A0A2K3CTW4, 3,210 aa, 705 bits) and the pilot returned `no_locus`. The
alignments were there — **9 of them at 24.6–28.4 % identity** — and
`MIN_LOCUS_IDENTITY` = 0.40 had discarded all nine, so the obvious reading was
that a fourth S5 threshold does not transfer, and that is what went into the
emergent list.

It was the **bait panel**. The only Chlorophyta bait was a *Cymbomonas*
prasinophyte, and 25–28 % is what a chained module match to a distant bait
looks like — precisely the population that floor exists to remove. Once the
panel carried *Chlamydomonas*' own record the locus came back cleanly
(`found_annotated`, one full locus, annotation-matched). The floor did not
misfire; the row has been corrected to say so.

What survives is a narrower caution. S5b justified 0.40 by a wide empty gap
(confirmed loci ≥ 0.759 against junk at 0.23–0.34) measured where *every*
genome has a same-class bait. Bands here are whole phyla and 13 slots are
unbaited, so a genome whose nearest bait is a phylum away could still lose a
real locus to it. S23b should re-measure it S5b's way and report how many
`no_locus` calls sit just under the floor.

**2. S5's R3 spread rule was not ported, and that is what left *Chlamydomonas*
with no near bait.** S5 required a band's second bait to come from a different
NCBI order than its first. S23's bands are whole **phyla**, far wider than
S5's classes, so the rule matters more — and without it Chlorophyta filled
both its slots with two *Cymbomonas* records of one genus while
*Chlamydomonas*, hundreds of millions of years away, got none.

Porting the rule as a *filter* made things worse, which is worth recording:
applied only at selection time it removed the duplicate *Cymbomonas* without
admitting *Chlamydomonas*, taking the panel from 85 baits to 78 and Chlorophyta
down to one. The reason is that the shortlist handed to the screen was ranked
by profile score alone, and Chlorophyta's nine shape-passing records are led by
**seven *Cymbomonas* ones at 814–1,039 bits against *Chlamydomonas* at 705** —
so a six-deep score-ranked shortlist never reached it. **A spread rule has to
shape the shortlist, not filter its output.** `wanted` now takes the best
candidate per order first, in score order, and fills the remainder from the
rest.

**3. B6 was letting a reused seed consume its band's quota** rather than add
to it, which is the mechanism behind (2): Chlorophyta's single slot went to
S3's *Cymbomonas* seed and no census bait was derived at all. "Reused, not
re-derived" means the seeds are kept; it never meant they replace the
derivation. Fixing it took the panel from 68 baits to 85 (42 ITPR + 16 RyR +
27 MIR control, 221,994 residues).

Two reporting bugs were fixed alongside: the control manifest was writing its
`strength` / `clade_frac` columns nowhere, so the report rendered "0 strong,
0 weak"; and the unfilled-slot test compared against quota+seeds, marking
seed-only bands as gaps.

**Next.** S23b, in this order: (1) re-measure `MIN_LOCUS_IDENTITY` for this
scope against annotation-confirmed loci; (2) find a second, non-MIR control
for Apicomplexa, without which "0/60" cannot go to assembly level; (3) measure
locus-span / CDS-footprint and, where chaining is severe, count non-overlapping
full-length alignments rather than loci; then the full 194-genome sweep, the
absence claims at assembly level, and census v6.

---

## 2026-09-06 — S23b: the blocking measurements, and a control chosen by data

**Task.** S23b as the brief defined it: three blocking items before the full
sweep, then the sweep. The three turned into a session; the sweep is running
and became S23c.

### 1. The control (blocking item 2) — and it was blocking a headline claim

S23a left *Toxoplasma gondii* `uncontrolled`: neither a receptor nor a MIR
locus, so "Apicomplexa 0/60" could not go to assembly level at all. The
diagnosis in the S23a report was "find a second, non-MIR control". That was
half right. **The fault was not PF02815, it was fixing one profile in advance
for every clade.** A control's job is to fire in the clade whose absence is
the claim, so which protein family makes the best control is a property of the
clade — and it is measurable.

`s23_control_profiles.py` declares six candidates, every one a large, deeply
conserved, multi-exon eukaryotic family (MIR, Myosin_head, E1-E2_ATPase,
SMC_N, Kinesin, and AAA as an explicit fallback), so whichever a clade picks
the control exercises spliced alignment over a long gene rather than one small
domain. `s23_control_select.py` runs each over each clade's own swept
reference proteomes — 20 searches over the four archived DBs — and takes the
one that is actually there.

- **Apicomplexa takes Myosin_head: 36/36 Aconoidasida and 23/23 Conoidasida
  swept proteomes, median ~1,450 aa**, against PF02815's one protein per class
  across 60. *Toxoplasma* now returns 2 control loci and is
  `controlled_cross_kingdom`.
- **Rhodophyta vindicates the design independently.** Red algae carry myosin
  in **8 %** of their swept proteomes — the fixed-profile approach would have
  handed them a control found in 1 of 12. Measured, they take SMC_N at 12/12.
- 27 clades take Myosin_head, 1 takes SMC_N. 56 control baits, 48 strong / 8
  weak.
- **The MIR bait stays in the panel regardless.** It is the only control
  inside the family's own signature set, so a MIR locus called ITPR is D14's
  sharpest possible failure. Dropping it for the better proof-of-search would
  buy a control and sell a negative control.
- AAA is marked `fallback` and is *searched only if the ranked candidates
  leave a clade empty*. None did, so it was not run, and
  `control_profile_coverage.tsv` records it as `not_searched` rather than
  leaving a blank a reader would mistake for "searched, found nothing".

**A second tier came free.** The panel already searches every clade's control
bait in every genome, so "a bait from another kingdom also aligns across this
locus" costs nothing to record and says far more than "the assembly is
readable": it shows the search crosses the divergence any receptor here would
have to be found across. 87 of the first 127 genomes are
`controlled_cross_kingdom`.

### 2. `MIN_LOCUS_IDENTITY` (blocking item 1) → D27

The brief asked for a re-measurement. The measurement could not be made at
all as the pipeline stood: the sweep filtered at the floor, so the loci a
calibration needs were never recorded. **S5b could measure 0.40 only because
it had the loci 0.40 excludes.** The sweep now records to
`RECORD_MIN_IDENTITY` = 0.15 and applies the call floor downstream, read back
from the committed calibration; `call_min_identity()` returns *whether* the
number in force is measured or inherited, so no report can present one as the
other.

`s23_calibrate_loci.py` measures it against loci whose identity the
assembly's **own annotation** establishes, classifying each recorded cluster
`confirmed` / `sister` / `contradicted` / `unnamed` / `no_annotation`. It
**refuses to write below 50 genomes / 15 confirmed loci** — and that guard
earned its place inside an hour: a smoke-test over the single re-swept pilot
genome derived 0.25 from two loci, wrote it, and the next sweep read it back
and moved *Neurospora*, *Toxoplasma* and *Oryza* out of `no_locus`.

Making the evidence class right mattered more than expected. `unnamed` had to
be separated from `contradicted` by a positive test on the name: outside the
vertebrates most gene models carry locus tags, and *Chlamydomonas* files its
receptor as `CHLRE_16g665450v5`, *Strongylocentrotus* as `LOC594527`. Treating
an uninformative name as a contradiction would have put two correctly-
recovered genes in the junk population.

### 3. Span against CDS footprint (blocking item 3) → D28

Measured: *Drosophila*'s 22 kb *Itpr* sits in a **297,487 bp cluster around an
8,514 bp CDS footprint — 35×**, because `-G` is 650 kb for the metazoa and 26
alignments chained across it. The status call was right; a copy count would
not be, because two real genes inside one such chain are counted once. Copy
number is now counted on **non-overlapping complete alignments**, not on
clusters. The *Drosophila* copy is bounded at 10.8 kb.

### 4. Two latent bugs, both found by the work above

- `s5_classify.ITPR_NAME_HINTS` had drifted from `family.KNOWN_NAME_SUBSTRINGS`
  and was missing `itr-1`, so the pilot read *C. elegans*'s correctly-recovered
  receptor as an annotation naming something else. Now sourced from
  `family.py` (CLAUDE.md's one-place rule). No vertebrate symbol changes.
- The rescue HSP filter read `h["start"]`/`h["end"]`, which `parse_tblastn`
  does not produce (`sstart`/`send`). It never fired because every `no_locus`
  genome up to *Salpingoeca rosetta* returned zero HSPs, so the generator's
  predicate was never evaluated. Replaced with S5's own `filter_hsps_outside`
  rather than a second reimplementation. *Salpingoeca* now returns
  `tblastn_trace`.
- Also closed: the recording floor had quietly made the control cells easier
  to satisfy than the family call (they were counted with no floor and no
  grade bar), so a genome could have been declared `controlled` on a chained
  fragment. Controls are now held to the same floor and grade bar.

### 5. The sweep

Panel rebuilt 80 → **108 baits** (37 ITPR + 16 RyR + 55 control), 277,611
residues. Launched over all 194 genomes; **127 done, 0 failures** at the time
of writing, ~8 h of miniprot left. All five downstream modules
(`s23_ledger`, `s23_calibrate_loci`, `s23_figures`, `s23_census_v6`,
`s23_report`) run clean against the partial set, and the report renders itself
at three scales so a 127-of-194 run cannot wear the finished sweep's heading.

**Partial results (127 genomes, smallest-first, so biased to small
assemblies):** 34 of 35 absence clades have at least one controlled genome and
**every one holds at assembly level** — Ascomycota 0/31, Basidiomycota 0/17,
Apicomplexa 0/3, Streptophyta 0/6 (one trace-only). Copy number 0/1/2/4 =
104/17/5/1. Two genomes `uncontrolled` (*Intoshia linei*, *Allopauropus
danicus*, both reduced metazoan genomes), correctly excluded from every
absence claim.

**Next.** S23c: resume the sweep, then `s23_calibrate_loci.py` → a `--redo`
reclassify pass to apply the measured floor → ledger → figures → census v6 →
report.

---

## 2026-09-06 (cont.) — S23c: the sweep completes, and the identity floor is retired

The 194-genome sweep finished (**194/194, 0 failures, 100 Gbp**), which
unblocked the whole downstream chain. What was meant to be a mechanical
finish turned into the session's sharpest methodological result.

### The calibration could not be done the way the brief asked

The brief said to measure `MIN_LOCUS_IDENTITY` "against loci whose identity an
assembly's own annotation confirms", as S5b did over 571 loci. Outside the
vertebrates that evidence barely exists: of **917 recorded clusters, 21** sit
on a gene whose name says anything — 10 name the family, 11 name something
else. Most gene models here carry locus tags (*Chlamydomonas*'s receptor is
`CHLRE_16g665450v5`; *Strongylocentrotus*'s is `LOC594527`).

So a **second axis** was added: every recorded cluster's translated model
scored against `itpr.hmm` / `ryr.hmm` (D23), re-derived from the archived
miniprot GFFs so the calibration reruns offline. 226 profile-confirmed.

### And then the measurement said the threshold does not work

- **Identity does not separate.** Confirmed loci reach down to **0.193**;
  contradicted ones reach up to **0.318**.
- **Coverage is no better.** Youden J 0.710 against identity's 0.695. A draft
  of the report claimed coverage separates where identity does not; the
  measurement contradicted it and the claim was removed.
- **S5b's inherited 0.40 would discard 87 confirmed loci, 56 of them complete
  gene models** — a third of everything the sweep found.

The reason is structural: with bands a whole phylum wide, a locus's identity
to its nearest bait measures *how far away the nearest bait is*.

**So identity was retired as a call gate** (set to the recording floor) and
the **profile call** carries it — which is where this project puts every other
family call (D14/D23). A retired threshold has to be validated against
evidence the replacement does not share, so the gate was scored against the
annotation axis: **10/10** on loci the assembly names for the family,
declining **10 of 11** it names for something else. → **D29**.

The one exception is worth naming: *Emiliania huxleyi* gene **`IPR1`**, at
1,008 bits over 95 % of its bait against 329 for `ryr.hmm`. `IPR1` is not in
the family name list, so the annotation axis scored it as naming a different
gene. On the evidence the name list is short, not the gate wrong — but `ipr1`
was **not** added as a name substring, because it collides with every InterPro
accession. It is reported as a conflict.

### Three bugs, two of them pre-existing in shared code

- `s3_assign.assign` crashed twice on a hit whose full-sequence bit score is
  exactly 0.0 with no counterpart from the other profile: `winner` was picked
  by score comparison and then used to subscript a record that is `None`.
  Never fired in S3 or S5; fired immediately when the profiles were first run
  over genomic gene models, where marginal alignments at or below 0 bits are
  normal. The winner is now picked as a *record*. Ties still resolve to
  `unassigned`, so no existing call changes.
- The rescue HSP filter (S23a) read `h["start"]`/`h["end"]`, which
  `parse_tblastn` does not produce. Replaced with S5's own
  `filter_hsps_outside`.
- **The control now requires a complete recovery**, symmetric with the family
  call: a fragmentary control shows the search finds fragments, which is not
  what an absence rests on. Admissibility is stated as an explicit set, not a
  `startswith("controlled")` prefix test — which is exactly how the new
  `controlled_partial` verdict would have been admitted silently.

### Results

**All 35 absence clades hold at assembly level, and all 35 are controlled.**
Ascomycota 0/31, Streptophyta 0/25, Basidiomycota 0/17, Magnoliopsida 0/3,
**Apicomplexa 0/3**, Microsporidia 0/2, and 29 more. **Zero genomes
uncontrolled** — 115 `controlled_cross_kingdom`, 77 `controlled_by_target`, 1
`controlled_partial`, 1 `no_control_bait` — against 1 of 14 uncontrolled in
the pilot.

**Copy number**, the deliverable: 0 in 116 genomes, 1 in 43, 2–6 in 32, then
*Dysidea avara* 8, ***Stentor coeruleus* 13**, ***Macrostomum lignano* 18**.
That answers S20a's open question — *Macrostomum*'s 62 database records
resolve to **18 real genes** — and *Cymbomonas* to **3 copies from 2
clusters**, the other one. D28's copy rule recovered **7 genes in 7 genomes**
(174 against 167) that cluster-counting merges.

**D14 held everywhere**: 189 of 195 graded loci have a family margin of 1.0,
minimum 0.597, and **0 control loci were ever called ITPR**.

Census v6: **18,065 rows** (+183 genomic models). 100 of them are
`annotated_unnamed` — the dominant class outside the vertebrates, and the
population S18 will want.

**Next.** S6 — the alignment upgrade. Two things it inherits: copy number
ranges 0–18, so a representative rule assuming one gene per species is wrong;
and any rule stated in sequence identity has to be re-derived (D29).

---

## 2026-09-06 — S6: the representative alignment

**Task.** S6 — alignment upgrade. MAFFT L-INS-i + trimAl over census
representatives chosen per clade × per kingdom (D8), RyR outgroup included.
→ `results/msa_v2/`.

### What ran

`s6_select_reps.py` (census v6 → 134 representatives), `s6_msa.py`
(MAFFT L-INS-i 64.2 min single-threaded, trimAl, both identity matrices,
conservation, coverage, site classes), `s6_figures.py` (4 figures),
`s6_report.py` + `s6_report_results.py` (`report.md`, `stats.md`),
`s6_test_selection.py` (21 controls, all passing).

### The two rule changes the data forced, before the alignment was built

Both were found by looking at what the selector actually picked, not by
planning:

- **A paralog label is vertebrate-only.** The first dry run grouped
  *Acanthamoeba castellanii* and *Tetrabaena socialis* as `ITPR2`, from a
  UniProt protein name reading "receptor type 2". ITPR1/2/3 are a 2R
  product; those numbers are annotation transfer. `group_of()` now returns
  a paralog only inside the vertebrates, and the raw label survives in the
  audit with `paralog_source`.
- **A bait attribution is not a paralog label where it is constant.**
  Checking why the cyclostome grid cells were thin turned up 19 ITPR
  records over 12 loci with six labels, *all* `ITPR1`, *all* from the S5
  bait attribution. The sea lamprey and both hagfishes carry three
  full-length loci each and the ITPR1 bait wins all six.
  `informative_attribution()` now decides per band by counting distinct
  attributed cells, and R3 takes the complete copy set of two species per
  band rather than one tip per species — three lamprey genes against
  ITPR1/2/3 is a test of 2R, one is a sample. (**D30**)

The second change cost a restart of a 40-minute MAFFT run. It was worth it:
two thirds of the cyclostome evidence would not have been in the alignment.

### Bugs found and fixed

- **The domain track was in the wrong coordinate system** — it counted
  ungapped positions in the *trimmed* human row and called that the residue
  number, which renumbers every residue after the first discarded column.
  Domains landed at ~2/3 of their true position and the two C-terminal ones
  vanished off the end. Fixed by committing trimAl's own `-colnumbering`
  output as `column_map.tsv`. Found by looking at the figure (D11); no
  table check could have caught it. (**D31**)
- Tip labels collided: the three *Myxine* loci are all
  `GCF_964187855.1|ITPR1|…`, so a two-field label made one name for three
  sequences. The duplicate-label guard rejected the set rather than
  silently aligning 132 of 134.
- `pick_diverse` reset its per-key tally each wave, so wave 2 could add two
  more of a key wave 1 had already taken — five SAR slots went to three
  *Triparma* species. Keys are also nested now: `(phylum, genus)` spreads
  on phylum first, which is what stopped three oomycetes taking the slots
  Ciliophora should have had.
- A species-dedup key on the raw census string treated *Prymnesium parvum*
  and "Prymnesium parvum (Toxic golden alga)" as two species.
- `s6_msa.py` now refuses an `aln.fasta` that is not an alignment of the
  current `representatives.fasta`. A long single-threaded MAFFT and an
  edited selector overlap easily and the failure is silent —
  `align_stats.json` would carry the SHA-256 of a file the alignment was
  not built from, which is the drift D24 exists to make visible.

### Results

**The alignment.** 134 representatives, 374,650 residues → **11,796
columns** (76.3 % gaps) → trimAl `-automated1` → **1,790 columns kept
(15.2 %)**, 7.83 % gaps, **96.6 % parsimony-informative** (1,730 columns).
Median tip coverage 0.96; one tip of 134 below half.

**D14 confirmed on the separation.** Within a vertebrate paralog group
0.908, each paralog group to the RyR outgroup 0.259, **separation 0.650**
against S1's 0.579. Compared on the separation rather than the absolute
values, because S1 scored 31 control sequences against their nearest bait
on pairwise alignments and this table averages all pairs of a trimmed MSA —
different estimators of the same quantity.

**The sister question, previewed.** ITPR1 × ITPR2 **0.788**, ITPR2 × ITPR3
0.746, ITPR1 × ITPR3 0.736. The leading pair's interquartile range
(0.773–0.806) does not overlap either other pair's, so the ranking holds
across the middle half of every comparison. It is still not a phylogenetic
estimate; S7's AU test is the answer.

**The cyclostome trio, and its control.** All six cyclostome loci fall
nearest ITPR1 (0.73–0.85). The obvious objection is that ITPR1 may simply
be the slowest-evolving paralog, so everything deep is nearest it — so the
same statistic was run over the groups that are certainly not vertebrate
paralogs. Invertebrates lean ITPR1 31/47 but at a **median margin of
0.007**; protists 0.001; the RyR outgroup 0.002. The cyclostome margin is
**0.043, six times that**. The lean is a fact about those loci. Whether it
means a cyclostome-specific expansion from an ITPR1-like ancestor or 1:1
orthologs with ITPR2/ITPR3 diverging after the cyclostome split, identity
cannot say — S7 and S8 can.

**Next.** S7 — the ML phylogeny. It reads `trimmed.fasta`, roots on the six
RyR tips, and has three things to test that S6 handed it as hypotheses and
not results: ITPR1/2/3 monophyly, the ITPR1 × ITPR2 sister preview, and the
cyclostome trio.

## 2026-09-06/07 — S7, the ML phylogeny (and an S6 re-run it forced)

**Task.** S7, picked as the single `in_progress` ledger row. The tree, the
model, the AU test on the three sister hypotheses, RBH on every naming call
the tree contradicts, figures and report.

### What ran

The session opened onto a partly-built S7: a finished ML tree and analysis
tables, three constrained trees, and an AU test that had never been run.
Finishing it exposed three method defects, two of which had already
invalidated the AU result sitting on disk.

**D32 — the constraint named every tip.** IQ-TREE's `-g` places freely
exactly the taxa a constraint *omits*; a taxon that is listed, even in a
top-level polytomy, is pinned outside every group the constraint declares.
The constraint files named all 134 tips and so forced the unlabelled
vertebrate tips out of the paralog clades the tree nests them in, equally in
all three hypotheses. The AU test rejected every one at ΔlogL ≈ 1,500 —
including the arrangement the ML tree itself holds at 100/100. `T5` now
fails a constraint that names a free tip.

**A prefix collision, found by reading `ps`.** The previous session's `au`
had been interrupted and its three constrained searches orphaned; they were
still running 50 minutes later, writing the same `--prefix` as the three I
had just started. The AU test was reading whichever checkpoint `.treefile`
happened to be on disk. Fixed twice over: every stage now resumes on the
`.iqtree` report rather than the checkpoint treefile, and `claim()` refuses
a prefix a live process owns.

**The RBH step could not fire.** It filtered `naming_conflicts.tsv` to
`reassigned`, of which this tree has none, and reported "the tree
contradicts no census label" while the table it had just read held five.
Widened to the whole table; `compared_against` records whether a row was
tested against the tree's proposal or the census name.

**D35 — and this one came from a reader's question.** Asked why *Volvox
carteri* was missing from the figure, the answer turned out not to be "it
isn't in the census": it is there, called ITPR at high confidence by both
instruments, 4 of 5 signatures, S20 verdict `real_gene`. It lost the S6
representative slot on `length_fit` — against a ruler borrowed from another
kingdom. `length_targets()` used a group's own median only with ≥5
complete-architecture records and otherwise took the *global* median
(2,694 aa, a metazoan number); Viridiplantae had 2 and Amoebozoa 1.
Fixing that alone gave the slot to a *second Chlamydomonas*, because
`pick_diverse()` throttles on a key level that has one value across the
candidates. Both fixed, both now carry negative controls (C7, C8), and C7
was verified to fail against the pre-fix implementation.

**The re-run.** S6 and S7 end to end: selection → MAFFT L-INS-i (61 min,
single-threaded by D24) → trimAl → matrices → model selection (55 min) →
ML search (71 min) → analyse → AU + `--bnni` in parallel (~2.5 h) → RBH →
figures → report. One interruption: `s6_msa.py` was launched with system
`python3` and died at the identity matrices with no biopython — *after*
MAFFT had succeeded. Recovered with `--stats-only`, which reused the
alignment on a SHA-256 match rather than repeating the hour. The stale
`matrices` section it left behind had already been picked up by S6's
figures and report; caught by checking mtimes against the input.

### Results

**The sister question is answered: ITPR2 + ITPR3, ITPR1 outside.** AU over
10,000 RELL replicates rejects ITPR1+ITPR2 (p-AU 1.8e-05) and ITPR1+ITPR3
(1.65e-05); ITPR2+ITPR3 (0.476) and the ML tree (0.525) survive and carry
the same pair. **S6's identity preview picked ITPR1+ITPR2 and is
contradicted** — the pair identity ranked highest is the one likelihood
rejects hardest.

**The tree.** 134 tips × 1,797 columns, `Q.insect+R7`, logL −215,452.0,
69.5 % of 131 internal nodes clearing both thresholds. `--bnni` weakened
and lost nothing; the one claim not clearing both thresholds, the bare
ITPR1 core, was already below them (47.8/95 → 47.5/73), and the tree
prefers a slightly different ITPR1 grouping — the 19-tip extended clade at
100/100.

**All 5 disputed names upheld by RBH**, so those are the tree's uncertainty
about where to hang a tip, not annotation error.

**The cyclostomes are neither of S6's two readings.** All six loci sit in
two well-supported cyclostome-only clades, each holding hagfish *and*
lamprey — duplications older than the hagfish/lamprey split, not a
lineage-specific expansion and not three 1:1 ohnologs. *(pending: S8)*

**S6's 3R stress test failed for the right reason (D33).** Neither
same-species pair is sisters, but both are broken *only by other tips of
the same paralog* — the signature of a duplication older than the species,
which is what teleost 3R is. The naming passes; the expectation was wrong.

**The re-run changed one tip and no conclusion.** *Volvox carteri* replaces
*Tetrabaena socialis*; it pairs with *Chlamydomonas reinhardtii* at 100/100.
The model was re-selected independently and came back `Q.insect+R7` again.
Support improved (64.9 % → 69.5 % of nodes). The sister answer held.

**Next.** S8 — synteny. It inherits two questions this tree could not
settle: which side of the vertebrate duplication each cyclostome lineage
attaches to, and the ITPR1 core's weak support.

## 2026-09-07 — S8: synteny, and a null for every number in it

**Task.** S8 — flanking-gene analysis across the ITPR loci. Deps S2 and S5b
were complete; the row was the topmost unblocked `pending`.

**What was built.** Ten modules, all under the 500-line budget:
`s8_flank_lib.py` (loci, windows, symbol keys, Jaccard), `s8_control.py`
(the matched random-window null), `s8_paralogon.py` (the 2R test and the
consensus caller), `s8_run_synteny.py` (driver), `s8_tables.py`,
`s8_test_flanks.py` (11 negative controls), `s8_priors.py`,
`s8_figures.py`, `s8_report.py` + `s8_report_results.py`. Full run 45–60 s
over 2,144 loci in 309 genomes; 274 of them carry a gene table.

**Three method decisions that changed the answer.**

*Loci come from the per-genome `summary.json`, not the ledger.* The ledger
carries one row per genome × cell and therefore only the best locus. A
teleost ITPR1 cell holds *itpr1a* and *itpr1b*; the sea lamprey's holds
three. Reading the best one would have compared *itpr1a* in one species
against *itpr1b* in the next and reported the mismatch as a synteny result.

*Every Jaccard is scored against a matched control pair in the same two
genomes.* Naming density across this genome set varies ~4× (human 97 % of
coding genes named, sea lamprey 27 %), so a raw Jaccard confounds orthology
with annotation depth. The control holds the two genomes, their naming
conventions, the window and the key rule constant; the paired comparison is
a sign test with both-zero ties dropped and counted.

*The paralogon is invisible to symbol matching by construction.* 2R
ohnologs almost never share a symbol, so cross-paralog Jaccard is zero
whether or not a paralogon exists. The root key (`BHLHE40`, `BHLHE41` →
`BHLHE`) is what makes the pair visible, and the same rule applied to the
random windows is what stops a promiscuous root family manufacturing one.

**Results.** Within-paralog Jaccard 216×–413× its own matched null, with
98–99.8 % of individual pairs beating their own control; every
cross-paralog and cross-family class at or below the null (max mean J
0.0002 over 168,241 ITPR × RyR pairs). Exactly two ohnologous flank
families survive, **both connected to ITPR1** — BHLHE40/41 (ITPR1–ITPR2,
62 %/85 % of species, 84× background) and GRM7/GRM4 (ITPR1–ITPR3,
42 %/53 %, 93×) — and **ITPR2/ITPR3 share none at any bar from 10 % to
50 %**. ITPR3's neighbourhood is the one that does not travel: cross-class
J 0.062 against 0.160/0.158.

**The caller.** A flank consensus paralog caller, leave-one-species-out,
threshold chosen by maximising call rate minus random-window false-call
rate (0.4; accuracy is 1.000 across the whole sweep, so it separates
nothing and is not what is optimised). 405/405 correct on 503
annotation-confirmed loci, 6/726 random windows called. It adds **131
paralog assignments above the null's own maximum**, 87 of them at loci
whose annotation names no paralog.

**What it could not do.** S7 handed S8 the cyclostome question explicitly.
All six loci carry 20 informative flank symbols — the window is not the
limit — but the best overlap any reaches with a gnathostome consensus is 2,
and random windows reach 2. Reported `underpowered`, not negative, and
logged as an emergent task needing a name-independent instrument.

**Two bugs the discipline caught.** `frac_sweep()` originally varied a
module constant that `ConsensusCaller.consensus()` had bound as a default
argument, so the sweep silently reported one threshold six times. And a
two-run diff found `flank_consensus.tsv` differing between processes:
ranked tables built by walking a Python set are hash-seeded unless the sort
carries a final tiebreak. T11 now tests order invariance on every build,
and `synteny_stats.json` records the SHA-256 of every table.

**Next.** S9 — ML selection (dN/dS), deps S6 only, which is complete.

## 2026-09-07 (cont.) — S9a: the codon alignment, and two silent failures it walked into

**Task.** S9 as written is one session's worth of instrument and several
sessions' worth of PAML, so it was split in the ledger the way S5 was:
**S9a** builds the codon alignment every selection test stands on, **S9b**
runs the models. S9a is complete; S9b's suite is running and resumable.

**What the instrument had to guarantee.** dN/dS is a statement about codons,
so every sequence in the alignment has to provably encode the *exact*
protein S6 aligned and S7 built its tree from. 57 vertebrate family tips,
43 through UniProt cross-references (Ensembl → ENA → RefSeq `coded_by`) and
14 through a miniprot locus realignment. **57 / 57 validated.**

**The first silent failure was in the ported route itself, and it hit the
two records the paper leans on hardest.** The PIEZO implementation returns
the first CDS a route successfully downloads. A UniProt entry cross-references
every Ensembl transcript of its gene, and the entry's own sequence is one
particular isoform: human ITPR1 lists **five** transcripts and human ITPR2
**two**, and in each the first is not the isoform S6 aligned. That does not
fail — it substitutes a different isoform for the protein the tree was built
on, which is a *wrong* codon alignment rather than a missing one. First run:
`validation_failed` on human ITPR1 and human ITPR2 and nothing else. The
routes are now candidate generators and the caller keeps the first that
**validates** (D36); human ITPR1 needed nine candidates, ITPR2 three.

**Masking is now total, not partial.** The port masks internal stops. It
leaves any *other* translation/protein disagreement in place, which hands
pal2nal a pair that does not agree — and pal2nal resolves that by dropping
the sequence, quietly. Every disagreement is masked to `NNN` now (24 codons
across the whole set, 15 of them internal stops), and the aligned protein is
written with `X` at those positions so the two files always agree.

That change is also what recovered the 57th tip. *Hymenochirus boettgeri*
ITPR3 spans 795 kb; a locus rerun indexes 800 kb where the sweep indexed
3.2 Gbp, so miniprot placed the **first exon** differently — 9 residues of
2,666, same length, everything else exact. The port refuses anything but an
exact match. Same length and within 1 % is now accepted as the same gene
model with those nine codons masked, and beyond that is still refused
(`s9_test_codon.py` T11 tests all three directions).

**Selection sets come from the tree, not the census.** codeml's branch
models mark a *node*: a foreground that is not a clade does not fail, it
marks a larger one and returns a well-formed ω for a hypothesis nobody
asked. So the three sets are S7's **extended paralog clades**, re-derived
from `rooted.nwk` with S7's own rule and cross-checked against
`paralog_clades.tsv` — a size or membership disagreement is a hard failure.
ITPR1 19, ITPR2 13, ITPR3 19. That nests **7 unlabelled `vertebrate_basal`
tips** inside paralog clades (D30 read forwards) and leaves **6 in no clade
at all — exactly the six cyclostome loci S7 handed to S8 and S8 reported
underpowered.** They stay in the whole-tree analyses as background, because
dropping them would change the branch lengths every other estimate is made
on, and they are in no foreground.

**The second silent failure was mine, and the self-tests did not have it
until it happened.** The per-paralog subset files were written from a
`set`, so their row order is hash-seed dependent: same alignment, same
likelihood, different SHA-256 on every run — and rebuilding to add a column
rewrote `codon_ITPR1.phy` underneath a running codeml job. This is S8's T11
in a new place. Fixed by writing rows in the order they are asked for, and
now tested two ways: T13 that the order follows the input, and **T14 that
three separate interpreters at different `PYTHONHASHSEED` values agree** —
same-process repetition cannot see this class of bug. Confirmed T14 fires on
the set-based version before trusting it.

**S1's toolchain manifest never probed this task's tools.** codeml, yn00,
`pal2nal.pl` and `hyphy` have been in the `piezo1` env all along and were
not in `TOOLS`, so S9 opened by reporting `pal2nal.pl` missing on a machine
that had it. Added, along with `stdin=DEVNULL` and a scratch cwd in
`probe()` — codeml *prompts* for a control file on a terminal, and in a
directory holding a `codeml.ctl` it would start a real analysis instead of
printing its banner. codeml reports no version of its own; yn00 from the
same build does, and the manifest records that the number is borrowed.
PAML 4.10.10, pal2nal v14, HyPhy 2.5.101.

**Numbers.** 3,253 codons; trimAl `-automated1` keeps 2,459 (75.6 %),
chosen on the protein and applied codon-aware. PAL2NAL cross-checked
nucleotide-by-nucleotide against an independent in-house mapping for all 57
sequences. 14 negative controls pass on every build.

**S9b, in flight.** ITPR1 ω = **0.0238**, ITPR2 ω = **0.0430**; the
curated-CDS sensitivity subsets give 0.0212 and 0.0318, so the genome gene
models are not driving the estimates. The result that changes how the rest
must be read: **synonymous saturation is reached inside a single paralog,
not only between them** — 94 % of within-ITPR1 and 85 % of within-ITPR2
pairs exceed dS = 1.5, and dN plateaus near 0.1 while pairwise dS runs past
50. A paralog set spanning shark to teleost to mammal has burned its
fourfold-degenerate sites. The pairwise matrix is therefore a diagnostic and
not an estimate, and every ω reported comes from a tree-based model. Logged
as an emergent task: re-estimate ω inside a shallow clade where dS is still
determined, and report the two side by side.

**Branch-site model A is restarted by construction (D37)** — four initial ω
on every paralog stem, with the spread committed — rather than repaired
after the fact the way the PIEZO project had to.

**Next.** S9b: finish the codeml suite (`python scripts/s9_codeml.py`, fully
resumable), then `s9_relax.py`, `s9_tables.py`, `s9_report.py`,
`s9_figures.py`.

## 2026-09-08 — S9b: the models land, and three significant tests that do not mean what they say

**Task.** The codeml suite finished overnight: **40/40 jobs in 22.7 h** on
7 workers, plus 3 HyPhy RELAX runs. 12 likelihood-ratio tests, BH-corrected
across the family. `scripts/s9_run.py` ran `codeml → relax → tables →
report → figures` unattended and all five stages returned clean.

**The headline is agreed on by three instruments that share no
machinery.** ITPR1 is held roughly twice as tightly as the other two:

| | ITPR1 | ITPR2 | ITPR3 |
|---|---|---|---|
| one-ratio ω | **0.0238** | 0.0430 | 0.0415 |
| two-ratio ω (foreground vs background) | **0.0241** vs 0.0432 | 0.0435 vs 0.0317 | 0.0455 vs 0.0308 |
| RELAX k | **9.36 intensified** | 0.908 relaxed | 0.836 relaxed |

Every two-ratio test is significant (ITPR1 q = 2.9e-63). RELAX compares
whole ω distributions rather than point estimates, so its agreement with
the one-ratio ranking is a check rather than the same number twice.

**Three significant tests needed the fitted parameter before they could be
stated, and none of them meant what its p-value looked like. This became
D38.**

*M8 vs M7* beats its null in all three paralogs at q ≈ 2e-4 — and the class
it adds sits at **ω = 1.00000**, codeml's boundary, on 0.3–0.7 % of sites,
with 0 / 0 / 1 sites reaching BEB ≥ 0.95. A beta distribution on [0, 1]
cannot represent a spike at the neutral boundary, so adding one class that
lands exactly there fits significantly better and says nothing about
adaptation. *M2a vs M1a* is 2ΔlnL = 0.00 in all three: the positive class is
estimated at ω = 20–94 with a proportion of **exactly zero**.

*Branch-site model A* is significant on all three stems, and on two of them
ω₂ is pinned at codeml's **999 upper bound** with the likelihood flat above
it — ITPR2's restarts reach the same lnL at ω₂ = 162 and 999, which is the
definition of an unidentified parameter, and exactly what §4.2's saturation
predicts for a branch that old. **Only the ITPR1 stem is reported as a
result**: ω₂ = 5.53 on 11.0 % of sites, stable across three restarts, 8
sites at BEB ≥ 0.95 and 5 at ≥ 0.99. A pipeline that printed the three
q-values would have reported its strongest signal on the stem whose
parameter is least determined.

**D37 earned its keep, measurably.** 3 of the 12 branch-site restarts
converged *below* their own nested null — one on every stem — and at a
**different initial ω each time** (ITPR1 ω₀ = 4, ITPR2 0.5, ITPR3 1.5). No
single starting value would have been safe, so this is the case for
restarting rather than for choosing better, and it is now a measurement
rather than an argument. The figure shows it at a glance: three points left
of the zero line, one per colour.

**One silent failure, in the direction that hides a result.** HyPhy writes
non-finite per-branch ω estimates as the bare token `inf`, which is not
legal JSON. The ITPR1 RELAX run *succeeded* — its log carries a complete
fit — but `json.loads` threw, `run_one` returned `{}` on the exception, and
ITPR1 came back as `k = None` with 0 branches. Beside two real answers that
reads as a negative result, and it would have buried k = 9.36, the largest
effect in the analysis. Now repaired on read, counted in the table, and a
run whose output cannot be parsed gets a `status` rather than a blank `k`.

**Scheduling notes, for the next long task.** Longest-first ordering put 18
expensive whole-tree jobs ahead of 6 cheap site models, which then waited a
day behind them; a second driver on the cheap ones cleared them in 2.5 h.
That is only safe because `run_job` now takes a **per-job process lock** —
two codeml processes in one directory interleave their fixed output
filenames and still parse into a plausible number (`s7_run.claim`'s incident
in a second tool). RELAX was moved off the critical path entirely once it
was noticed that it depends only on the codon alignment and the tree, both
of which existed twelve hours earlier; it then took 80 minutes.

**Estimating.** My first estimate (6–8 h) was a guess and wrong by ~4×.
Calibrating codeml's own `rub` round counter against finished jobs —
rounds-to-convergence per free parameter, and minutes per round — gave
7.0–8.7 h per whole-tree job against the 10–11 h I had extrapolated from
M2a, and predicted 05:30; it finished 06:37.

**Next.** S10 — annotation-bug molecular validation (deps S5b, complete).
S9's one live lead is logged as emergent: the eight ITPR1-stem BEB sites
should be carried onto the cryo-EM channel in S17, because "a handful of
sites changed fast on the branch that made ITPR1" means something very
different in the IP₃-binding core than in a disordered linker.

---

## 2026-09-08 — S10: annotation-bug molecular validation

**Task.** S10 (deps S5b, complete). Brief: take the two worst annotation
failures the sweep surfaced and prove them at the molecular level.

**Case selection, because "the two worst" has to be a rule.**
`s10_case_spec.py` writes it out and `s10_select.py` applies it to every
locus the sweep recovered (880). The measurement is **annotation loss** —
the share of a gene's coding footprint no single annotated model delivers,
read off S5's own `frac_cds` rather than recomputed. Five eligibility rules;
**E5 is the one that earns its place**: *does this annotation build genes
this long anywhere else in the same genome?* Without it the ranking's top
rows are *Cirrhinus mrigala* (longest annotated gene genome-wide 42 kb) and
*Saguinus oedipus* (138 kb against a 726 kb locus) — genome-wide length
ceilings that would have been written up as bugs at this gene. E5 removes
exactly those 6 loci in 2 genomes and nothing else in the sweep.

**Result: 382 eligible loci, median loss 0, 360 at exactly 0, 7 failures in
3 genomes.** Three modes labelled (`omission` / `truncation` /
`fragmentation`); two are selectable and the rule takes the worst of each,
one case per genome, because taking the top two of the single ranking gives
two omissions and leaves the fragmentation claim unvalidated.

**Case A — *Nibea albiflora* ITPR2, omission.** 56 coding exons, all 56 with
no annotated model; chromosome assembly, 84× headroom; the same annotation
gets ITPR1 and ITPR3 right. Splices into 2,673 codons with **0 internal
stops** against 14.2 expected under neutral drift (computed from the locus's
own codon usage, not quoted). 55/55 introns spliceable; 52/55 exon
boundaries shared by a majority of 35 independently annotated genomes; and
S8's committed flank consensus places it between **SSPN and BHLHE41**,
ITPR2's two most conserved neighbours (197 and 183 of 215 species). The
genome files **8 of 14 family-named models as pseudogenes** — 31.6 % of its
whole gene set is pseudogene — and its "ITPR2"-named model is the **ITPR3**
gene (89.3 % to the ITPR3 locus, its own translated protein). That
adjudicates one of the three naming conflicts S5b deliberately left open.

**Case B — *D. eleginoides* ITPR3, fragmentation.** One 68 kb gene as three
protein-coding models tiling residues 1–67, 52–467, 468–1593, with the 3′
41 % unmodelled. 59/59 canonical introns, 59/59 boundaries shared by a
majority of 40 genomes, 0 internal stops.

**Both species: 0 ITPR protein records in any database, 3 complete genes in
the DNA.**

**The transcript step is an honest negative with its denominator.** Remote
BLAST was abandoned after ~50 min queued and replaced with a local search of
every transcript record NCBI holds for each species — 43 and 10, no TSA.
Faster, archived, offline on re-run. RNA-seq does exist (182 and 19 runs) and
is handed to S12 with the probes committed.

**The control caught a real bug in the counting rule.** Running the probes
against the locus's own genomic DNA — where nothing can span a junction by
construction — returned **6 false spans of 55**: blastn extends a
high-scoring alignment a dozen bases past the junction into the intron, which
clears an 8 nt anchor. Added a probe-coverage requirement derived from what a
probe *is* (contiguous spliced sequence, so a genomic match tops out near
half its length) rather than from the artefact's size. Control now 112 hits,
0 spanning. `n_spanning_anchor_only` is committed so the difference between
the two rules stays visible.

**Two other bugs found by building.** The tiling subject set was census v4's
*non-redundant* models, so a query whose own locus was missing landed on its
nearest paralog 19 Mb away and read as a confident naming disagreement —
fixed with an `at_own_locus` coordinate check plus a complete subject set.
And the ORF step read its expected protein from `novel_models.faa`, which
lacks 3 of the 8 family loci in these two genomes; it now reads the sweep
GFF's `##STA` row, which exists for every model.

**47 negative controls** (`s10_test_evidence.py`) run before any table is
written and the driver refuses to write if they fail.

**Next.** S11 — Structures (deps S5b, complete). Two things it inherits are
in the roadmap's next-session note; the sharpest is that AFDB coverage has to
be measured against census *records*, since a locus with no protein record
cannot have a model, and S10 found two swept species with three complete
genes and zero records between them.

## 2026-09-08 (cont.) — S11: structures, and a family AlphaFold DB does not hold

**Task.** S11 — AFDB coverage, TM-align against the cryo-EM IP3R and RyR
references, per-domain pLDDT, optional Foldseek sweep. Dependencies S2 and
S6 both completed. Completion criteria: `structure_manifest.tsv`, a TM-score
table, per-domain pLDDT, a `report.md` rendered from the tables, and a
calibration figure showing the negative controls. All pass.

**What ran.** Fifteen new modules (`scripts/s11_*.py`, 4,300 lines, every
file inside the 500-line budget), driven by `s11_run.py` through ordered
stages `afdb → panel → tmalign → plddt → foldseek → tables → figures →
report`. The whole pipeline reproduces from cache in 25 s; the cold run was
dominated by TM-align at 5,373 s of wall clock over 435 pairs.

### Results

**AlphaFold DB does not hold this family.** 20.9 % of the census's 8,319
UniProt-shaped ITPR records have a model covering ≥ 95 % of the protein —
but of the **5,861 records at or above the family's own 2,000 aa floor, 13
do (0.2 %)**. Median modelled record 392 aa against 2,674 aa unmodelled.
Coverage is a function of length, not taxonomy, and it is concentrated
exactly on the fragments a structural argument can do least with. Of S6's
134 representatives — the set every alignment, tree and selection result in
this project stands on — **9** have a usable model.

**AFDB answers a canonical accession with an isoform.** Asked for the three
human paralogs it returns `Q14643-4` (2,695 of 2,758 aa), `Q14571-2` —
**181 residues of a 2,701-residue ITPR2** — and only ITPR3 canonical. A
probe taking `payload[0]` reports all three as modelled. `afdb_probe` now
ranks every record in the response, prefers the queried accession, and
records what was served; the 181-residue model enters the panel and is
rejected there by the minimum-chain rule rather than quietly used.

**References resolved by query.** 237 RCSB entries / 414 polymer entities
enumerated on the family's own Pfam signatures, each entity's family decided
by *this project's census* on the UniProt accession RCSB maps it to — no
title is read. Enumerating on the **union** of the four signatures is
load-bearing: RCSB's Pfam annotation of 6DQN carries PF02815/PF08454/
PF01365/PF00520 and **not** PF08709, so a query on the naming signature
alone misses the project's own IP3R reference. Primary references 7LHF
(ITPR1, 2.96 Å), 9YKK (ITPR2, 2.95 Å), 8TKG (ITPR3, 2.50 Å), 9NMO (RYR1,
2.40 Å), 7U9X (RYR2, 2.58 Å). **No full-length RYR3 cryo-EM entry exists.**
Plus a six-state ITPR3 panel and three negative controls derived from S1's
committed decoy panel; the MIR-sharer class is reported unfilled — no
experimental entry.

**D14 confirmed by a sixth instrument that reads only coordinates: 20/20
(100 %) of the callable census-named structures agree, and 0/3 negative
controls receives a family call.** Calibration on the panel itself: same
protein in different conformations 0.78, IP3R × IP3R 0.43, IP3R × RyR 0.39
(0.64 the other way round), controls 0.19 with **0 of 81** control pairs
reaching the 0.50 same-fold bar under either normalisation.

**Per-domain confidence: the IP3-binding core is the best-modelled domain
(median pLDDT 83.9) and the pore the worst (71.0)**, against 69.5 outside
the annotated domains. That is the answer S17 and S22 needed — the part
those tasks are scoped around is the part predicted structures carry best.

**Foldseek over the AFDB Swiss-Prot subset returns the family and nothing
else.** 36 hits confirm census ITPR records, 545 fall below the fold bar,
and all 8 distinct above-bar hits the census has never held are SDF2 and
SDF2L1 across human, mouse, cow, *Arabidopsis* and *Dictyostelium* —
verified from UniProt to carry **PF02815 and nothing else**. Zero novel
structural leads. AFDB holds no RyR model at all, so a sister-family verdict
there is unreachable by construction and the report says so.

### Three things the controls and the self-tests caught

**The negative controls broke the family call, which is what they are for.**
Gated on the relative margin alone — D7's 10 %, inherited from every earlier
stage — the rule calls **all three controls ITPR**: a dynein heavy chain, a
Cav2.1 and a talin each beat their own runner-up by ~30 % of their score
while scoring 0.17–0.27 against everything. A margin between two non-matches
is still a margin. The call now requires the winner to clear TM-align's
published 0.50 same-fold bar before the margin is read (**D39**). With the
gate the controls are declined and 20/20 real structures still agree.

**A killed process was cached as a permanent negative.** An early ad-hoc
TM-align run was interrupted; killing its children left the Python parent
orphaned to init, and it spawned work beside the driver that replaced it for
~20 minutes — S7 recorded the same incident for IQ-TREE. Worse, the
SIGKILLed pairs were written to the cache as zero-score results with empty
error text: **ten permanent false negatives, all against one control**,
caught only by `s11_tmalign_run.self_test`'s requirement that every pair
parse a score. Fixed both ways (**D41**): a non-`ok` result is never
written, and the stage claims a pid-checked lock that a *dead* pid does not
hold. Both now have negative controls in the self-test.

**Two claims I wrote that the data falsified.** I stated that every X-ray
ITPR entry is a ≤ 604 aa binding-core construct — 5GUG and 5X9Z are 2,217 aa
cytosolic-domain crystals at 7.3–7.4 Å, and the report now computes the
ranges instead of asserting them. And I explained the declined
non-vertebrate models as "too short to reach the bar however good they are";
adding the **arithmetic ceiling** (residues over reference length) showed
every one of them had the headroom — ceilings 0.50–0.63 against scores of
0.30–0.45 — so they fall short on similarity, not on length alone. That
column is now in the table (**D40**).

### Instrument notes

- The self-test is 26 constructed checks and was **mutation-tested**: three
  deliberate rule breakages (drop the compound-state pattern, invert the
  TM normalisation, make `largest_chain` return the first chain) and all
  three were caught.
- The pair classifier was duplicated in the figure and the report with
  different label strings, and a whole class silently vanished from one
  panel. Now one function, `s11_tables.pair_class`, with a shared ordering.
- `structures_stats.json` hashes every `.tsv` in the directory rather than
  the ones a given invocation wrote — a `--from tmalign` run had been
  recording 8 of 14.

### Next

S12 — expression evidence (SRA junction-spanning reads) for whichever
paralog or lineage the census leaves in doubt. S10 left it two validated
loci with committed junction probes and measured RNA-seq availability
(182 and 19 runs), which is the natural starting panel.

---

## 2026-09-08 — S12: expression evidence

**Task.** S12 — expression evidence (SRA junction-spanning reads + atlases)
for whichever paralog or lineage the census leaves in doubt. Completed.

### What the census left in doubt, and why it was this

S10 finished by naming the question it could not answer. It had found 7 loci
where the IP3-receptor gene this project recovers from the genome reaches no
annotated gene model, in assemblies that deliver the other 375 of 382 whole,
and it wrote: *"RNA-seq for it does exist and reaching it needs the streaming
aligner S12 builds; these junctions are handed there rather than
half-answered here."* So the scope was inherited, not chosen — and it is
**derived from S10's committed ranking** by four rules in `s12_panel.py`
rather than hand-listed: every eligible locus above 0.5 annotation loss, plus
every other family locus in the same genome as an internal control. That
resolves to exactly S10's 7 failures across 3 species, covering all three
failure modes, with 5 controls beside them.

### What ran

67 public RNA-seq runs (**536,000,000 reads**, 16 studies, 11 tissues),
streamed `fastq-dump -X 4000000 | hisat2` against a per-species reference
holding every recovered family locus, three housekeeping anchors and a
reversed decoy for each. 2 h 30 m wall clock, 0 failures, resumable per run.

### Results

- **7 of 7 loci are transcribed and spliced.** Each is detected in 13–31 of
  its runs and 6–9 tissues; 15,263 reads read through their splice junctions.
- **298 of 314 (94.9 %) of the junctions no annotated model spans are
  crossed by reads.** At full depth that is the same rate as the junctions
  the annotation *does* model (95.2 %).
- **Read coverage falls outside the annotation in the proportion S10
  measured from the coding footprint** — 7/7 within 0.25, four loci at
  1.00 vs 1.00. Same quantity, different evidence (aligned reads vs a GFF).
- ***Nibea albiflora* has no IP3 receptor in any protein database by any
  route** — ITPR2 unannotated, ITPR1 and ITPR3 annotated *as pseudogenes*
  at 0.98 and 0.97 of the coding footprint, all three transcribed in 31–32
  of 32 runs.
- **The deposit cross-check is still `underpowered`, now with a species
  that should have answered it.** *D. mawsoni* has 37,166 mRNA records;
  854 junction probes return 0 spanning hits — and 0 for the annotated RyR
  control. Median deposit length 547 nt against an 8 kb transcript. The
  genomic negative control fires (1,791 hits, 0 spans), so the criterion
  discriminates rather than never firing.

### Three inherited assumptions that did not survive measurement

1. **The validation rule.** Validating a reference by translating it and
   requiring high identity is the wrong instrument: a frameshift costs the
   frame from where it sits, so a *correct* reference scores 1.00 with no
   frameshifts and 0.92 with nine, and any floor across that range rejects
   correct references. Replaced by **colinear block placement**, which has
   no tuned threshold — a block may fail to place only if the model's own
   frameshift count explains it. Reported per locus against its own broken
   versions: reversed order 1.7 %, wrong strand 0 %.
2. **The decoy floor does not transfer.** The PIEZO project measured 0 decoy
   reads across 78 runs. Here: **42 reads in 2 of 469 run × locus
   comparisons**, both on one decoy, confined to ~64 bp of 8,185 — a
   pileup, not porous mapping. The prior renders `contradicted`; no
   detection call moves, because the decoy is a per-run floor and the
   affected locus clears it in the same run.
3. **The closed-set risk is measurable and was measured.** PIEZO argued it
   away in a caveat. Every reference tiled exhaustively with synthetic
   reads and mapped back: **0 of 12,500 cross-mapped**.

### Instrument notes

- The self-test is 11 constructed checks run before anything is written, and
  was **mutation-tested on 5 deliberate rule breakages, all 5 caught**. Two
  slipped through the first attempt: a junction-class lookup by index rather
  than coordinate (a one-intron test cannot distinguish them — T5 now uses a
  two-intron minus-strand gene), and an anchor-rule break my harness had
  applied to the wrong file.
- `miniprot --trans` is load-bearing for the housekeeping anchors: without
  it there is no `##STA` line, every model has an empty protein, and the
  coverage filter rejects them all — which reads as "this genome has no
  housekeeping genes".
- Housekeeping bait accessions are **resolved by query against a declared
  length band**. The first version hard-coded three and got all three wrong,
  including a 427 aa "GAPDH" (the enzyme is 333) and an accession serving
  nothing. A wrong bait still aligns and still produces reads.
- The three anchors are unequal for a measured reason: EEF1A1 has 5–6
  genomic copies, GAPDH 3, RPL13A 1–3, and copy number orders their
  detection exactly. The reference holds one model each, so surplus copies
  multi-map below the MAPQ floor.
- `write_tsv` is now atomic (temp + rename): stages overlap in practice and
  a half-written `reference_table.tsv` read by the quantifier is a short
  table with no error.
- `expression_stats.json` hashes every `.tsv` in the directory, not the ones
  a given invocation wrote — the same fix S11 needed.

### Next

S13 — reconciliation and dating. The loss audit it owes the ledger must
distinguish a real absence from an annotation absence, which S10 and S12
have now separated for these 7 loci; and "present in the census" differs
from "present in the proteome" by more than S3 measured, because a gene can
be annotated as a pseudogene and reach no protein record at all.

---

## 2026-09-08 — S13: gene-tree / species-tree reconciliation

**Task.** S13, the topmost pending row with its dependency (S7) complete.
Data root attached, 681.5 GB free. Nothing bulk was downloaded: every input
was already committed.

### What ran

`scripts/s13_run.py` — ordered stages `species_tree → reconcile → losses →
cyclostome → stats → figures → report`, self-tests first. Whole task runs in
about 5 s; the cost here was in the rules, not the compute.

- **`s13_species_tree.py`** — the accepted species tree as an **input**
  (D15): 31 species, 29 named internal nodes, each with a literature age,
  the spread of published estimates, a **stem age** and its source.
  `--check` validates before use.
- **`s13_reconcile.py`** — 5 topologies (ML, the three AU-scored sister
  constraints, the `--bnni` re-search) × 3 variants (all tips / cyclostome
  loci pruned / unsupported nodes collapsed) = 12 reconciliations, 3
  refused. Plus the minimum-event rooting check over all 112 edges.
- **`s13_losses.py`** — every species × paralog cell asked of the S5 ledger.
- **`s13_cyclostome.py`** — the six cyclostome loci joined to S8's flank
  calls, and the long-branch check.
- **`s13_test_recon.py`** — 14 groups of constructed negative controls.

### What it found

- **The two duplications are not on the same branch.** ITPR1 vs
  ITPR2+ITPR3 sits on the **vertebrate stem** (older than crown Vertebrata,
  published estimates 480–615 Ma, no upper bound this tree can set); ITPR2
  vs ITPR3 sits on the **gnathostome stem**, bracketed **462–563 Ma**.
- **The topology is irrelevant and the taxon sampling is everything.** All
  five gene trees give 14 dup / 53 loss with the cyclostome loci and 12 / 51
  without — including the two AU-*rejected* sister arrangements. Drop the
  six cyclostome tips and the older placement falls back to the gnathostome
  stem.
- **The placement does not rest on the weak node.** The ITPR1 split sits on
  a 17.4/54 gene-tree node; collapsing every node below S7's own bar merges
  two vertebrate-stem duplications into one carried at 100/100 and the
  placement holds. Minimum-event rooting picks a *different* edge (64 events
  against the outgroup rooting's 67) and agrees on the placement.
- **The long-branch objection was measured and does not apply.** The six
  cyclostome tips are 0.96–1.04× the median root-to-tip distance, ranking
  18–55 of 57.
- **0 of 51–53 implied losses survive contact with the genomes.**

### What bit, and what changed because of it

- **The classic LCA duplication rule is only valid on binary trees.** The
  support-collapsed root is a four-way polytomy mapping Cyclostomata /
  Gnathostomata / Cyclostomata / Gnathostomata; no child maps to Vertebrata,
  so the binary rule calls it a *speciation* — the collapse would have been
  reported as the duplications disappearing. Replaced with the non-binary
  rule (Vernot et al. 2008). **D44.**
- **A constrained IQ-TREE search writes no support values**, so the collapse
  variant dissolved all three constrained topologies into one 134-tip
  polytomy and reported it as a collapse of 40 nodes. `collapse_unsupported`
  now refuses a tree with no support labels, and the three cells are
  recorded as refused with their reason.
- **The first loss audit reported four corroborated losses and every one was
  a bait-panel limit.** ITPR2 and ITPR3 read `absent` in both cyclostomes
  while both genomes carry **three ITPR loci apiece**, all filed in the
  ITPR1 cell because S5 has no cyclostome-labelled bait. Added the
  `paralog_unassignable` verdict, a positive test on the genome's spare
  locus count (S5) and the species' tree-unplaced tips (S7). **D45.**
- **T11 failed first because my expected loss counts were wrong, not the
  code.** Two copies in one species duplicate *in that species* (0 losses),
  and a 1:1 speciation across the root implies 2. The test now carries four
  hand-derived cases with the derivation in the docstring.
- Mutation-tested on three deliberate rule breakages, all three caught.
- `s13_report_results.py` hit 514 lines and was split into
  `s13_report_results.py` (§5–§6) + `s13_report_audit.py` (§7–§11), with
  the two matrix helpers moved into `s13_lib.py` so both halves read the
  matrix the same way.

### Next

S15 — loss dynamics. Its character matrix must come from the genome sweep,
not from tips of the tree (S13's own loss count is 0 of 51–53 once audited),
and `paralog_unassignable` is a state it needs: a Dollo count that reads the
cyclostome cells as losses would score two independent losses per cyclostome
that never happened. The 43 `tblastn_trace` cells are the ones the brief
asks S15 to disambiguate by synteny, and S8's consensus caller with its
committed calibration and null already exists for that.

---

## 2026-09-08 — S15a: the loss instrument

Split the S15 ledger row (the protocol's rule for a task too large for one
session). **S15a** is the character matrix and the three things it needs
before a loss can be counted; **S15b** is the counting, and S15a changed
what that task is.

### What ran

`scripts/s15_run.py` — eight stages, `recon → synteny → integrity → tree →
matrix → tables → figures → report`, all green, ~40 s end to end, entirely
offline apart from one archived `datasets` call for 991 ancestor taxon
names. 15 constructed negative controls run before anything is written.
→ `results/loss_dynamics/` (19 tables, 4 figures, `report.md`, 2.3 MB).

### What resulted

**No ITPR paralog is absent from any of the 309 vertebrate genomes.** 0 of
927 genome × paralog cells reaches the loss state; in all 189 assemblies
contiguous enough to carry the gene the minimum is 3.00 gene-equivalents.
783 cells are one placed locus, 90 are truncated by their contig, 7 are
partial, **43 are a gene reassembled across contigs**, and 4 are
`paralog_unassignable` — the cyclostomes, exactly the four S13 flagged,
recovered here by a rule reading the sweep's own locus counts rather than
S13's table.

**The brief's synteny step is answered by a number, not a call.** S8's
consensus caller is 100 % accurate at every key count it acts on (13/13,
14/14, 378/378) — and of 432 rescue regions in the undecided cells, **273
sit on a contig carrying no annotated gene at all**, 151 have too few keys
and **8 reach the four-key floor**. The median region has 0 informative
neighbours and its contig extends 39.9 kb, shorter than the gene. Where the
caller does reach, it agrees with the alignment attribution 6/6. Accurate
and unavailable.

**So the instrument is a reference reassembled across contigs.** Computed
outside every locus the aligner found (the sweep's own cross-paralog
exclusion, inherited not re-invented), one reference at a time, with the
bar measured against a decoy that cost no new search: regions the 38-bait
panel attributes to a paralog the aligner already placed at a locus in the
same genome. Candidates median 0.795 (min 0.174), decoy median 0.030 (max
0.100), Youden J = 1.00, bar 0.137 at the gap's midpoint.

### What went wrong on the way, and what it changed

- **The first calibration did not separate — J = 0.52.** The decoy was every
  region attributed elsewhere, and in a genome where two paralogs are *both*
  shattered a fragment attributed to the other one is a piece of a real
  gene, so the decoy contained genes. Restricting it to accounted-for
  paralogs gives J = 1.00. The middle population is kept and committed:
  `co_trace`, 20 regions, median 0.631 — it *is* the measured size of the
  paralog-attribution problem in a shattered assembly, and it is why S15b's
  primary coding has to be family-level presence per genome. **D46.**
- **The ORF screen's confounder is not contiguity.** Contig N50 barely moves
  lesion density (ρ = −0.077); bait identity moves it a great deal
  (ρ = −0.397, p = 1.5e-67). So the paired within-genome test is run twice,
  once identity-matched — and ITPR2's apparent indel excess **disappears**
  (q = 0.902) while **ITPR3's survives** (39 genomes to 14, q = 0.0032),
  with the RyR control showing no excess. **D47.**
- **`_fmt` wrote floats at four decimal places**, so a p-value of 2.1e-07
  landed in the table as `0.0000` and the report read it back as zero and
  rendered "p < 1e-300". Floats are now `%.6g`, and the report has its own
  `pfmt()` — a p-value at four decimal places hides how strong a claim is,
  not how weak.
- **Two negative controls failed on first run and both were the test's
  fault, not the code's.** T5 used a bait key the fixture did not define;
  T14 asserted that the human+mouse clade must be labelled `Mammalia`, when
  `_suppress_unary` deliberately keeps the *deepest* named node of a unary
  chain and with two mammals sampled that is `Euarchontoglires`. T14 now
  asserts membership and that the label is *some* taxon on the shared path.
- **`compare_with_s13()` first reported 21 of 29 curated clades as
  unrecovered** while every matched node carried the right name: it was
  counting assemblies of species S13 never sampled as intruders. Asked
  correctly — do any assemblies S13 places *outside* the clade fall inside
  it here — the answer is 23 of 23 testable clades recovered, 0
  disagreements.
- **`contig_spans_gene` cannot be keyed by genome.** The ledger's column is
  a property of a cell's best locus, and a last-wins dict over cells put
  most genomes on the wrong side of D4's bar. It now calls
  `s5_calibration.spans_a_gene` on the genome's own contig N50.
- Mutation-tested on three deliberate rule breakages, all three caught, and
  the stdlib t-tail validated against seven published critical values.
- `s15_report_results.py` is at 495 lines — inside the budget but the next
  addition to it needs a split.

### Next

S15b, and S15a has changed it. There are no losses to count, so the
**sensitivity matrix is the deliverable**: which combinations of coding,
evidence bar, branch lengths and contiguity filter *manufacture* a loss.
`loss_candidates.tsv` is built for it, with the rule that stopped each
near-miss in its row. Three constraints come with it — code family-level
presence as primary (D46); do not fit Mk rates to an invariant character;
there are no pseudogene fossils, so the shared-lesion Poisson test must be
reported with its denominator. The lead worth following instead is D47's
ITPR3 indel excess, which nothing in this project explains.

---

## 2026-09-08 — S15b: the loss counts

**Task.** S15b — Dollo parsimony as the primary count, Mk fits over the
coding × evidence × branch-length × contiguity sensitivity matrix, and the
pseudogene-fossil lesion analysis. Status: **completed**.

### What ran

`scripts/s15b_run.py`, ordered stages `dollo → sensitivity → mk → fossils
→ tables → figures → report`, 21 constructed controls before anything is
written, 3 mutations applied and all 3 caught. Everything it reads is
committed (S15a's matrix, S15a's 309-genome tree, S15a's ORF-integrity
tables, S13's node calibrations), so the whole task is offline and a rerun
is deterministic. Full clean run: ~150 s, dominated by the 522 Mk fits.

Twelve new modules, all inside the 500-line budget:
`s15b_lib` (Newick + the three branch-length schemes), `s15b_coding` (the
evidence ladder and the recoder), `s15b_dollo`, `s15b_mk`,
`s15b_sensitivity`, `s15b_fossils`, `s15b_tables`, `s15b_priors`,
`s15b_figures`, `s15b_test_counts`, `s15b_run`, and the three-way report
split `s15b_report` + `s15b_report_results` + `s15b_report_lesions`.

### What resulted

- **Dollo count 0**, family-level and per paralog, across 927 cells in 309
  genomes. The routine that returns it finds a constructed loss on a known
  edge (T1), merges sister losses (T2) and bounds a polytomy's (T4).
- **The sensitivity matrix is the deliverable.** 32 settings × 2 codings ×
  3 branch-length schemes. Family coding manufactures a loss in **2 of 32**
  settings (worst case 1 genome); paralog-resolved in **18 of 32** (up to 45
  loss edges, 37 independent). D46 was a judgement in S15a and is a
  measurement now → **D48**.
- Moving the reconstruction bar across the **whole gap its calibration
  measured** (0.100–0.174) changes no cell in 927. Turning D45 off alone
  manufactures the 4 cyclostome losses at every rung. Tightening the
  evidence alone, both rules on, manufactures **exactly one** cell —
  *Bothrops jararaca* ITPR2, reassembled at 0.796 across 6 contigs.
- **Branch lengths change a Dollo count in 0 of 32 settings**, reported as
  invariant rather than dropped (D48).
- **Mk stated, not fitted**: every likelihood monotone to its boundary on
  every model, axis and scheme — measured on a rate grid, with ARD's *gain*
  axis rising to the edge. All 12 fits at the operating point are refusals.
  The 405 fits that are possible are fits to manufactured losses; their rate
  spans a factor of 495 across the three branch-length schemes.
- **No pseudogene fossils**: 44 of 1,760 scored loci clear the lesion bar
  and **44 of 44 are at full coverage**. The generous three-reading screen
  fires on 7 ITPR loci, all at coverage 1.00 with 1–2 internal stops.
- **D47's lead is a bird result**: ITPR3 25 genomes to 2 in Aves
  (q = 4.5e-5) against 7 to 6 in Actinopteri (p = 1.0) — but 21 of the 27
  bird pairs are below D4's contiguity bar, so the verdict is
  `underpowered` → **D49**.
- **Unasked bonus.** Dollo's unpinned gain nodes land at **Vertebrata** for
  ITPR1 and **Gnathostomata** for ITPR2/ITPR3 — S13's placement, from an
  instrument that reads no gene tree, no alignment and no reconciliation.
  Reported as *orthogonal*, not as corroboration: it rests on the same
  cyclostome cells S13 used.

### What went wrong, and what it changed

- **T3 caught a real design gap on the first run.** Unpinned Dollo scores a
  clade-wide absence as *ancestral*, not as a loss. That is right for ITPR2
  and ITPR3 (S13 places their duplication inside the tree) and wrong for the
  family (S20/S23 found it across the eukaryotes). The gain rule is now
  declared per character in `s15b_sensitivity.PIN_GAIN`, with T3/T3b as its
  two halves.
- **The fossil screen flagged the entire RyR control.** `not_live` asked the
  ITPR character matrix for a state the RyR cell does not have, and got
  "not present" for all 26. The reading is now marked *unavailable* when the
  cell has no matrix row, and T16 is the control.
- **The first ARD likelihood profile was ER in disguise** — profiling along
  the diagonal q01 = q10 is ER by construction and drew the same curve twice
  in two colours. ARD is now profiled on each axis separately, and the gain
  axis rising to the grid's edge is the better statement anyway.
- **A figure used `hash()` for jitter**, which is not stable across runs.
  Replaced with a point's rank in its own row (D24 applied to a jitter).
- Two report bugs, both caught by reading the rendered output against its
  own tables: a filter on `loosest_evidence` printed "D45 off manufactures 0
  losses" beside a table showing 4, and the Mk sentence had its refusal
  count inverted. Both now read the number off `dollo_counts.tsv`.
- `s15b_report_results.py` hit 680 lines and was split into
  `s15b_report_results` (the count) + `s15b_report_lesions` (the half that
  argues against it), 377 and 381 lines.

### Next

S16 — duplication history: are ITPR1/2/3 2R ohnologs, and are teleost
itpr1a/itpr1b from 3R. Three things S15b puts in front of it. The gain-node
result restates the 2R question — if both WGD rounds predate the cyclostome
divergence, ITPR2 and ITPR3 should each have a cyclostome co-ortholog and
neither does. The 309-genome species tree exists and places every swept
assembly, so retention asymmetry no longer has to be stated on S6's 134
representatives. And the per-**locus** copy table must be read out of each
sweep `summary.json`, not out of the ledger, which holds only each cell's
best locus.

---

## 2026-09-08 — S16: duplication history

**Task.** S16, the only unblocked pending row (S7, S8, S13 all complete).
Are ITPR1/2/3 2R ohnologs, are the teleost itpr1a/itpr1b copies 3R, and what
is the copy-number landscape. Data root attached (679 GB free), dashboard
opened and watched, `git pull` clean.

### What ran

Twelve new modules, `scripts/s16_*.py`, driven by `s16_run.py` in eight
stages (`copies → map → paralogon → quartet → teleost → tables → figures →
report`). Full run 42 s warm, ~3 min cold (25 BioMart requests). The method
is ported from the PIEZO project's S16 (`../piezo_genes/scripts/s16_*.py`)
and every result is this family's own.

- **`copies`** — the per-locus table the brief asks for, read out of all 309
  sweep `summary.json` files rather than the ledger: 2,146 loci, 1,841
  copies, every count carrying D4's contiguity flag and S15a's lesion bar
  (read back off S15a's committed verdicts, not recomputed).
- **`map`** — human Compara paralogy from BioMart against a **pinned** dated
  archive (`jun2026.archive.ensembl.org`, Ensembl Genes 116), with the
  archive's own registry committed beside the map so the release is evidence
  and not a sentence. 22,564 undirected pairs touch a neighbourhood gene;
  2,559 of 3,567 S8/window symbols resolve (71.7 %).
- **`paralogon` / `quartet`** — the 2R test at three window sizes × three
  duplication-node vocabularies, against a real-window permutation null
  (D17, 5,000 draws per direction), BH at two scopes, plus a pooled
  one-hypothesis-per-family test, the genome-wide block scan, the quartet
  and the 309-genome replication with its own matched-random-window null.
- **`teleost`** — the five 3R predictions, flank sets joined out of S8's
  committed `flanks.tsv` (which flanked every locus, so this is a join and
  not a re-extraction), DCS against two references, and cross-anchor block
  identity.

### What resulted

- **2R.** ITPR1's neighbourhood is paralogous to ITPR2's in 141/175 genomes
  and to ITPR3's in 89/152, against **24/932 (2.6 %)** matched random
  neighbourhoods in the same genomes. ITPR2 vs ITPR3 is 4/149 — the
  background exactly. **Dating the links split them**: `GRM7 ↔ GRM4` is
  *Vertebrata* (84 genomes), `BHLHE40 ↔ BHLHE41` is *Opisthokonta* and is
  not an ohnolog pair. Both are the two families S8 found by root-key
  symbols; the date is what S16 adds (**D51**).
- **The RyR control is what makes the human test readable.** Pooled at ±10,
  ITPR scores 1 dated link against a null of 0.060 (p 0.039) and RYR 1
  against 0.022 (p 0.021). Neither survives correction across all 135 tests.
  So the single-genome test measures the instrument's ceiling, not a
  difference between the families, and the claim rests on the replication
  (**D52**).
- **3R.** 97.3 % of 73 teleost genomes above D4's bar carry two ITPR1 (mean
  1.97) against 1.04 for ITPR2 and ITPR3 — while the same genomes carry
  **5.82 RyRs**. Pre-3R ray-fins 1/1/1 + 3; extra-WGD lineages 3/2/2 + 8.
  The two copies partition the ancestral block disjointly in 45/49 and 46/49
  genomes against the two references, and **705 of 705** cross-anchor
  assignments agree (p 1.2e-212), corroborated 6/6 by the bait (sequence)
  call.

### What was built, broken and fixed

- The split-model merge **fired on 0 of 2,146 loci** — S5's own clustering
  already chained miniprot's alignments. That zero is only reportable
  because T1 constructs two halves that *each* clear the copy bar and
  requires one copy out; the first version of T1 used halves below the bar
  and would have passed whether the merge worked or not.
- **The quartet test's six significant pairs were circular.** Each was a
  window against the block selected, out of ~23,000, for being the most
  paralogous to that window. Now flagged `selection_circular` in the table
  rather than in the prose, and the headline is the 59 non-circular pairs,
  of which 0 are enriched (**D53**).
- The cross-anchor bait-concordance check was unrunnable because it took
  `anchors[0]`, which is ranked on flank richness and whose two copies won
  the same bait. It now takes the first anchor whose own copies differ; the
  check went from 0/0 to 6/6.
- Two anchor-selection and reporting bugs found by reading the rendered
  output against its own tables: `n_above_bar` was written inside a per-cell
  loop (all four cells share one denominator, so it silently reported the
  last cell's), and the report printed "six anchors from six orders (6)"
  because `headline()` kept only the metric value and not its note.
- `s16_report_results.py` hit 513 lines and was split into
  `s16_report_results` (copies + 2R) and `s16_report_3r` (3R + caveats +
  hand-off), 281 and 259 lines. Every S16 module is now under 500.
- 31 constructed negative controls pass before anything is written.
- **Reproducibility checked rather than claimed** (D24): a second full run
  reproduces all 24 committed tables byte for byte, SHA-256 against
  SHA-256, from the `duplication_stats.json` the first run wrote.

### Emergent

Three rows added or updated at session end. **S8's paralogon headline needs
one qualification**: it reports two retained ohnolog families between the
ITPR neighbourhoods, and dating them leaves **one** — any manuscript
sentence built on the pair must say one (S14a). **Separating 2R round R1
from R2 needs a non-symbol instrument**; S8 and S16 have now hit the same
wall from opposite sides, and the cyclostome loci are the evidence that
would resolve it. And the S5a row on half-named teleost 3R duplicates is
**now quantified**: of the 73 ray-finned genomes carrying two ITPR1 copies,
only 20 name both, 34 name one and leave the other a `LOC`, and 19 name
neither — with `itpr1b` used 45 times against `itpr1a` 15, so the databases
are systematically better at one 3R co-ortholog than the other (S18).

### Next

S17 — constraint & function. S16 hands it an asymmetry that singles out
ITPR1 twice over, from two independent kinds of evidence: its neighbourhood
is the one that retained 2R ohnologs with both other paralogs, and it is the
only paralog whose 3R duplicate was kept. Whether that is one fact or two is
a constraint question, and S9's per-paralog ω estimates are already
committed. `results/duplication/loci.tsv` is the per-locus table S8, S15a
and S16 each had to rebuild from the sweep summaries — anything downstream
needing more than one locus per cell should read it rather than the ledger.

---

## 2026-09-08 — S17: constraint & function

**Task.** The project's mechanistic payoff: which parts of the IP₃ receptor
are evolutionarily intolerant, is the IP₃-binding core under different
constraint from the pore, and is per-site conservation good enough to help
read the family's human variants. Dependencies S6, S9b and S11 all
completed; S17 was the topmost unblocked `pending` row in the analysis
block.

### What ran

`scripts/s17_run.py`, nine ordered stages, 79 s warm on a full rerun
(~15 min cold, of which HyPhy FEL is ~8 min and the two MAFFT passes ~4.5).
Fourteen constructed negative controls run before anything is written.

| stage | what it produced |
|---|---|
| `domains` | `domain_map.tsv`, `functional_sites.tsv` — the architecture in all three human numberings |
| `orthologs` | 264 / 249 / 265 orthologues per paralog from the 309 sweep `summary.json` files, plus the shape screen |
| `conservation` | four layers on every residue of every human paralog, the element tests, the metric controls |
| `variants` | 1,753 ClinVar missense records + 27 UniProt variants, the AUC and the paralog audit |
| `fel` | HyPhy FEL on S9's three codon alignments, 2,459 sites each |
| `paint` | 16 S11 ITPR structures × 2 layers into the B-factor column |
| `tables` / `figures` / `report` | 26 hashed tables, 4 figures, `report.md` |

### What resulted

**Depth first.** msa_v2 carries 19 / 13 / 19 tips per paralog, which cannot
score a column of a 2,700-residue protein. The sweep's `##STA` translations,
joined to `summary.json` by `mp_id`, gave **249–265 full-length orthologues
per paralog** — an instrument that had been on disk since S5b and had not
been used.

**The screen the coverage bar could not do.** Bait coverage is
`aligned_aa / bait_len`, so a model that covers the bait *and* carries two
thousand extra residues passes every quality bar — which is exactly what
S5b recorded a large `-G` manufacturing in the giant genomes. Scoring each
sequence on the fraction of its **own** residues inside reference columns
found one: *Lissotriton helveticus* ITPR2, 4,976 aa, 0.483 against a curated
minimum of 0.972 and a next-lowest of 0.944. The bar went in that gap
(0.713, both edges committed). With it in, MAFFT opened the ITPR2 alignment
to 5,676 columns; without it, 3,380. One chimeric model was inflating the
alignment every ITPR2 per-site score is read off by 68 %.

**The gate is the answer to the first question.** The gate and the
selectivity filter are the most constrained elements on both the JSD and the
composition-free metric, in all three paralogs, and **the gate's five
residues are 100 % identical between all three human copies**. RIH-associated
is the most constrained *domain*.

**The Pfam name for the ligand site does not contain the ligand site.**
Joining S0's per-accession InterPro coordinates to S0's 6DQN measurements
showed that **none of the ten IP₃ contacts lies in PF08709**, the signature
Pfam calls *Inositol 1,4,5-trisphosphate/ryanodine receptor*. They sit in
MIR and RIH — the domains shared with the ryanodine receptor. So the ligand
question was asked of the measured contacts throughout, and they are more
constrained than the rest of the domains carrying them (p = 0.013 / 6e-4 /
4e-3 against their own elements).

**One element was hiding inside another, and finding it changed the
headline.** Unresolved, the Pfam channel domain read as *less* constrained
than the receptor's own linkers — JSD 0.697, p = 0.96 for the one-sided
test. That is a strange thing to report about the pore of an ion channel, so
three explanations were measured rather than argued: a JSD/composition
artefact (ruled out — the composition-free metric agreed), a miniprot
gene-model artefact in the exon-dense TM region (ruled out — the curated
subset agreed), and an unresolved element (confirmed). A **50-residue
luminal loop, located by geometry** on 6DQN's own membrane span rather than
drawn on the conservation profile, is the cause. With it separated the pore
module clears the linker control (p = 0.005), and the loop is the least
conserved element in the receptor on every instrument available.

**The variant test came out against this task's own design.** 1,753 ClinVar
missense records, **0 dropped** by the transcript-numbering check — all
three genes file on a transcript whose translated CDS *is* the UniProt
canonical, which is a result of the check and not a reason it was
unnecessary (the same code on PIEZO2 found 512 of 773 positions
disagreeing). Both of S0's residue-level citations were recovered by the
positive control. **1,546 records (88 %) are VUS.** On one fixed set of 44
pathogenic and 34 benign positions, AUC = **family 0.872 > vert 0.854 >
deep 0.758 > shallow 0.684**: the depth was worth building (the `shallow`
control is the worst layer) but **taxonomic breadth beats within-gene
depth**, and this task's own instrument is the second-best of four. §7.2
says so with the number.

**FEL corroborates §5 with a different instrument.** 5,766 purifying sites
and **1** diversifying across the three paralogs — inside BH's own error
budget of ~288 false rejections, so not evidence of a selected site. Per
element: gate, filter and IP₃ contacts 100 % purifying at median β = 0;
luminal loop median β 0.33–0.40 with 27–33 % purifying against 71–86 %
protein-wide.

### Two things the self-tests and the report caught

**Mutation-tested on four deliberate rule breakages, all four caught by the
test responsible** — a suite that has never been shown to fail is a suite
nobody has checked. Gaps counted as observations in `column_stats` (T2), the
within-protein control selected by a prefix test on the termini (T8), the
containing element placed ahead of the pore elements in `PRIMARY_ORDER` (T6),
and the shape bar put at the lowest curated record rather than the gap
midpoint (T10).

**T8 found a real bug on its first run.** The within-protein control was
selected by a `startswith` test on `nterm`, `cterm` and `linker_` — and
`nterm_trefoil` is a *domain* whose name begins with `nterm`. 225 residues
of the element §2.1 is about were quietly inside the control set, so every
other element was being compared against a set containing one of them. The
membership test is now a prefix for linkers and an exact match for the
termini.

**S6's identity number is a different measurement, and it checks out.** S6
reports between-paralog covered identity at 0.741–0.791 and S17 measures
0.640–0.703. Neither is wrong: **S6 measured on `trimmed.fasta`**. S17
cannot, because trimAl deletes 39 of the luminal loop's 51 residues — the
element the headline rests on. Every row of
`paralog_identity_by_element.tsv` now carries the same pair measured S6's
way, and it reproduces S6's committed matrix exactly (0.753 / 0.773 /
0.812). The verdict is `orthogonal`, not `confirmed`.

### Decisions

**D54** — a coordinate carried from one protein to another must arrive on an
anchor the alignment does not know about (the filter's GGGVGD motif, the
gate's lining residues), and a failed anchor aborts rather than writing a
plausible number; and an element no annotation carries must be located by
measurement rather than drawn, or it is the profile explaining itself.

### Emergent

Five rows. The **luminal loop** is a new object — 50 residues, the least
conserved element in the receptor, paralog identity 0.13–0.31, worst-resolved
in cryo-EM — and two questions follow that are not S17's (is it the luminal
Ca²⁺/ERp44 insert, and is its length variable across the sweep, which S21
would see as a single variable exon). The **family layer beats the deep
layer** as a classifier, so S22 and any published variant resource should
read `family_jsd`. **ITPR2's one pathogenic record** becomes a falsifiable
prediction, since its gate and IP₃ contacts are as constrained as ITPR1's —
the emptiness is ascertainment. **1,546 scored VUS** are a submittable
artefact S14a's deposit rules do not currently mention. And **the gate
cannot distinguish the paralogs**, so a pore-motif shortcut for paralog
assignment is unavailable in principle — one line beside D14.

### Next

S18 — annotation-quality audit (`S5b, S15`, both completed; the topmost
unblocked `pending` row). S17 touches it only indirectly: it read the
annotation solely through S5's gene models, so the correction list is
unaffected. What S17 does hand forward is `constraint_<gene>_<acc>.tsv` for
S22's ligand-site question and S24's supplementary figures, and `painted/`
for any structure figure. The one S17 finding S18 should carry is that the
luminal loop is where a database "fragment" boundary would be least
surprising and least informative — a 50-residue low-complexity insert that
trimAl deletes and cryo-EM cannot resolve is exactly where gene callers
disagree.

---

## 2026-09-08 — S18: annotation-quality audit

Session protocol: dashboard opened and watched, `git pull` clean, data root
`/Volumes/FANTOM/IP3R_DATA` attached with 676.5 GB free. No `in_progress`
row; topmost unblocked `pending` was **S18** (deps S5b and S15, both
completed, priority high).

### What ran

Twelve new modules under `scripts/s18_*.py`, all under the 500-line budget,
plus a small additive refactor of `s10_gff.py`.

- `s10_gff.read_annotation_windows` / `read_miniprot_models` — one pass per
  file for many windows or many models, sharing `_build_genes` with the
  per-window readers so S10 and S18 cannot disagree about what a gene is.
  S10's 47 negative controls re-run and pass unchanged.
- `s18_run.py --only loci` — 2,144 gene-scale loci in 309 assemblies scored
  against 274 archived `genomic.gff.gz` files (4 min per full pass, cached).
- `s18_run.py --only protein` — 11,402 full-length family protein records
  resolved through `s6_lib.SequenceStore` (0 missing) and blastp'd against
  the committed 38-bait panel.
- `zero`, `corrections`, `tables`, `figures`, `report`.
- `s18_test_audit.py` — 45 constructed negative controls, run before
  anything is written; mutation-tested on 9 deliberate rule breakages, all 9
  caught.

### What resulted

**The audit's own premise is contradicted by its control.** ITPR loci are
73.9 % complete and 26.1 % failing; the ryanodine receptors, in the same
assemblies through the same pipelines, are 77.9 % / 22.1 %. No overall
difference survives BH correction (q = 0.13 raw, 0.82 above D4's bar). One
state does separate, and it is the family-specific one: an ITPR locus is
2.7x more likely than a RyR locus to be held *only* by a non-coding feature
(42 vs 16, q = 0.006) — and that does not survive the contiguity control
either (4 vs 5 above D4's bar, q = 1.0), so the excess is confined to
assemblies too broken to carry the gene.

**D9 is the largest effect in the task.** RefSeq gene sets deliver 98.8 % of
these loci complete; submitter-deposited GenBank ones 37.5 %, every state
differing at q < 1e-300. Held above D4's contiguity bar it is 99.4 % against
63.8 % and does not close, so about a third of the archive gap is assembly
quality and the rest is the gene set.

**D4 is the second.** The ITPR failure rate falls from 26.1 % to 6.7 % across
the contiguity bar — two thirds of what looks like an annotation problem is a
contig too short to hold a 2,700-residue gene.

**All 15 of S3's zero-hit reference proteomes are gene-caller failures.**
Every species has a genome in the S4 scope and every genome carries the gene;
0 `genome_also_empty`, 0 `undecidable_no_genome`. Eleven of the fifteen are
birds.

**The protein records are named correctly and cannot be found.** 4 of 11,402
sequence calls disagree with the census; 5 records are named for the sister
family, all non-vertebrate and all under 200 bits; 52 of 8,306 vertebrate
symbols name a paralog the panel assigns elsewhere. But 3,872 records carry a
placeholder gene symbol and 2,395 carry none — **55.0 % of the family's
full-length protein records have no usable gene symbol** — and 66 more are
named for the superfamily, which separates neither family.

297 corrections written, 52 `high` priority, 18 withheld under D6.

### Six things the build caught

1. **A bait lookup keyed on the wrong column.** `bait_labels()` keyed on the
   manifest's `id` (the full FASTA header) while the blast subject ids resolve
   to accessions. Every lookup missed and the sequence call came back
   `no_call` on all 11,402 records — which looks exactly like a family nothing
   can be assigned to, not like a bug.
2. **A name rule that manufactured 66 wrong-family errors.** UniProt's
   commonest name for a non-vertebrate family record is "RyR/IP3R Homology
   associated domain-containing protein", and "Inositol
   1,4,5-trisphosphate/ryanodine receptor" is close behind. Both name *both*
   families; `name_family` resolves them to RYR because it tests the RyR
   hints first. New verdict `family_ambiguous`, and `_names_both()` matches
   the inositol half separately because the shared word "receptor" breaks
   every ITPR hint substring.
3. **A paralog rule that put 6,660 records in the wrong-paralog cell.**
   "Inositol 1,4,5-trisphosphate receptor" with no type number is not a wrong
   paralog. New verdict `paralog_unspecified`.
4. **A cache that survived a rule change.** The measurement cache held the
   finished audit rows, verdicts included, so after (3) the committed table
   still carried the old labels and mouse *Itpr1* read as the annotation
   naming a different paralog, silently. Verdicts moved out of `measure()`
   into `apply_verdicts()`, and the cache now carries a SHA-256 of the reader
   modules. → **D55**.
5. **A model-id lookup that missed 21 loci.** miniprot restarts its
   identifiers at MP000001 *per chunk*, so a chunked genome's concatenated
   GFF repeats every one and the sweep disambiguates the duplicates. Reading
   the raw `ID=` attribute therefore missed every locus in the four giant
   genomes. Replaced with `s5_sweep_lib.parse_miniprot_gff` — the sweep's own
   parser, so there is one implementation rather than two.
6. **A silent fallback behind it.** A locus whose model could not be
   recovered fell back to measuring the annotation against the whole locus
   *span*, which includes every intron: one *Protopterus* ITPR2 got a 2.4 Mb
   denominator under an 8 kb gene, forcing `unannotated` whatever the
   annotation held. There is no fallback now — such a locus is reported
   `cds_unavailable` and left out of the denominator, and after (5) there
   are none. Both are covered by T19, which resolves every locus of a real
   chunked genome.

### Decisions

- **D55** — a cache may hold what a parser found, never what a rule decided.
- **D56** — a threshold this project already has is not re-derived by the
  task that inherits it; it is validated. S18 used
  `s5_classify.ANNOT_CDS_FRAC` (0.50) and spent its calibration measuring
  where 0.50 sits in the distribution it is applied to: over 1,077
  correctly-named, fully-recovered loci a single model covers a median 0.993,
  so the inherited bar is that distribution's 1.3 % point and is conservative.
  Across bars 0.30–0.95 `complete` moves only 76.4 % → 68.1 %.

### Emergent

Four rows added: the non-coding demotion as a validation target for S10's
machinery; the 55 % symbol deficit as a per-method number S19 needs; the
missing RYR3 bait, which 74 records now depend on; and a route for actually
submitting the 297 corrections, which no ledger row covers.

### Next

S19 — methods results (`S5b, S15, S18`, all now completed; the topmost
unblocked `pending` row). S18 hands it two things directly. The per-method
contribution question has a protein-side counterpart it did not have before:
55 % of full-length family records have no usable gene symbol, so a
name-driven search reaches under half of what the databases hold — the
counterpart to S5b's 318 DNA-only models. And S19's contiguity-confounder
step now has a measured version of exactly its question: the ITPR failure
rate across D4's bar (26.1 % → 6.7 %), with the RyR control beside it and the
archive (D9) separated out.

---

## 2026-09-08 — S19: methods results

**Task.** S19 — what each search method was worth, written as results.
Status `completed 2026-09-08`. Whole task rebuilds offline in 43 s
(`python scripts/s19_run.py`).

### What ran

Eight stages over committed artefacts. One derived asset: the swept
accession universe, built by one `grep '^>'` pass over each of the seven
reference-proteome FASTAs (42 GB, ~15 min, cached under the data root, so a
rerun is offline). The bait ablation re-parses the 309 retained
`miniprot.gff` files and caches each genome's alignments in compact form.

- `contribution` + `recovery` — census growth v1→v6, per-channel record
  sets, the head-to-head inside each searched database by length band, the
  per-gene recovery question, and cost against yield.
- `contiguity` + `floor` — the false-negative rate on both control series,
  its bins and strata, the tests, the floor scan and the neighbourhood check.
- `panel` — 19 bait panels × 309 genomes × 4 cells, validated against the
  ledger first.
- `drift` — the seven jackhmmer runs re-derived from the raw logs, each kill
  rule scored as a classifier.
- `inference` — S9/S13/S15a/S17's limits recomputed from their own tables.
- `tables`, `figures`, `report`.

37 tables, 4 figures, `report.md` (564 lines), `methods_stats.json`.

### What resulted

**The control design.** S15b reconstructs no losses anywhere in the scope, so
all 923 assignable cells hold a gene that is there and every ledger cell not
`found` is a false negative of the method: **140/923 (15.2 %)**. The RyR
sister cell is an independent replicate using none of S15a's states:
**42/309 (13.6 %)**, Fisher *p* = 0.58. The four `paralog_unassignable`
cyclostome cells are in neither series.

**Contiguity is the whole of it.** Missed cell median contig N50 23,460 bp
against 3,396,515 bp; odds of finding the gene 8.1× per tenfold (RyR 20.0×);
3/512 ITPR and 0/172 RyR cells missed on chromosome-level assemblies. Below
D4's bar: ITPR3 30.0 %, ITPR1 39.2 %, ITPR2 43.3 % — S5b's span bias as a
false-negative rate. 57 of the 182 misses sit where the paralog's own S8
consensus neighbourhood is also missing.

**D4's a-priori bar survives calibration.** 142,212 bp gives 0.9 % residual
on the ITPR series, 0.0 % on the RyR, retaining 189/309 genomes. Conservative
against 5 % (reached at 100 kb), about right against 1 %. Cost: 38.8 % of the
scope, disproportionately the margin species.

**The panel ablation reproduces the ledger 1,236/1,236.** Four human baits
recover 782 of 783; dropping a clade band costs ≤ 2 cells; dropping the RyR
control changes no ITPR call; dropping one paralog's baits costs 238–241. The
three unlabelled baits change 0 cells and alone recover 0. 68 of 2,179 call
changes are gains, all one mechanism (top-*scoring* vs best-*covering* bait
at the 0.70 bar). One bait alone recovers the gene at any identity above 0.5.

**Profile HMM vs domain annotation, inside one database.** At gene scale the
sweep adds 1 record to 3,135 in the vertebrates and 2 to 1,021 in the
non-vertebrate metazoa; its entire gain is under 1,000 aa — except in the
protists, where it adds 89 gene-scale records. Iteration returned no record
the profile pair calls family that one pass had not (4,960 of 5,130, plus
19,969 non-family). Three maximally-unlike seeds intersect on 4,785 family
records and differ by ≤ 162.

**K1 fires on 0 of 7 jackhmmer runs, including all 3 that drifted.** Scored
against a drift outcome measured on the finished model: K1 sensitivity 0.00,
K2 0.33, K3 1.00 at specificity 0.25. Moving K1's own 0.10 threshold to the
off-family share gives 1.00/1.00. Proposed, not applied.

**940 of 1,232 demonstrated genes (76.3 %) are unreachable from any protein
database**, and it is not a margin-species artefact (74.3 % vs 78.5 %).

### Bugs found and fixed

- **A negative control was overwriting a committed table.** T11 called the
  real `criterion_trace`, which wrote its two constructed rows over
  `kill_criterion_trace.tsv`; the report then read 2 jackhmmer runs where
  there are 7, and nothing failed. Fixed with a `write` flag, and T15 now
  checks the SHA-256 of every committed table before and after the suite.
  Mutation-tested: T15 catches it. → **D60**
- **jackhmmer target names were never intersecting the database universe.**
  The log names targets as `sp|ACC|NAME`, which `acc_key` deliberately leaves
  untouched (a genome model id also carries pipes), so every membership test
  returned empty and iteration silently looked as if it had added nothing.
- **The head-to-head was counting raw profile targets.** 18,501 vertebrate
  targets were scored and 5,130 called; the first version credited the sweep
  with the 13,371 D22's gate declined.
- **A logistic p-value printed as 0.** `1 - Φ(z)` cancels to exactly zero
  above |z| ≈ 6 in double precision; switched to `erfc`.
- **The self-test crashed on a clean tree.** Four checks read tables the
  stages had not written yet; they now report *skipped* and the suite re-runs
  from the `tables` stage.

### Decisions

- **D57** — a search's sensitivity is measurable, not estimable, when the
  family has no losses; and the control must be able to fail.
- **D58** — D4's contiguity bar was chosen a priori and is now calibrated; it
  stands, and its cost is printed beside it.
- **D59** — the kill rule written to catch iterative drift measures the wrong
  axis; the replacement is proposed rather than applied.
- **D60** — a self-test must not be able to damage the artefact it tests.

### Emergent

Five rows: the panel could be four baits (and the same ablation should be run
on S23's 68-bait non-vertebrate panel, where breadth probably is not free);
the one-line D10 fix and where it should be validated next; the three
quarters of genes with no protein record and what a submittable form would
be; the ablation's silence about tblastn rescue; and the 57 misses whose
whole genomic region is absent.

### Next

**S21 — gene architecture.** The instrument is already validated and cached:
S19's panel simulation reproduces the ledger 1,236/1,236 by re-parsing the
309 retained GFFs. Scope any architecture claim above D4's bar, or it will
measure contig lengths rather than exons.

---

## 2026-09-08 — S14a: manuscript assembly

**Task.** The submission package, built by a script from the committed
tables. Topmost pending ledger row with every dependency complete.

### What ran

`python scripts/s14_assemble.py` — ordered stages `figures → claims →
stitch → pdf → deposit` — exits 0:

| stage | result |
|---|---|
| figures | 7 main + 14 Extended Data, 56 panel files, 0 missing |
| claims | 175/175 load-bearing numbers re-verified |
| stitch | 16 sections, 13,822 words, 29 references resolved |
| pdf | 48 pages, A4, 1.2 MB, 0 missing glyphs |
| deposit | 1,767 files, 237 MB, SHA-256 per file |

### What was written

The 16 numbered section files, which are the only hand-written files in the
package. Title: *Retained in every vertebrate, lost repeatedly elsewhere: a
503-genome census of the IP₃ receptor family*. The results are structured on
the three the project actually has — the family is ancestrally eukaryotic
and lost repeatedly outside the animals; no vertebrate paralog is lost
anywhere in 309 genomes; and three quarters of the demonstrated genes are
unreachable from any protein database — with the duplication history and the
constraint map as the two mechanistic sections between them.

Also: `manuscript/reviewer_checklist.md` (the self-audit and 7 open items
for a human), a rewritten `manuscript/README.md`, and `scripts/s14_refs.py`.

### What the build caught

- **The claims ledger failed 5 of 175 rows on its first run**, every one
  because the claim addressed the wrong row rather than because a number was
  wrong: a vertebrate class named `Cyclostomata` where the ledger records
  `Hyperoartia` and `Myxini`; a synteny pair class quoted without its
  `cross_` prefix; two AlphaFold rows that needed the `ALL` group to be a
  single row; and one report phrase quoted loosely.
- **Reading the tables rather than this roadmap caught two stale ledger
  entries.** The alignment is 11,777 columns trimmed to **1,797**, not
  11,796 trimmed to 1,790 as an earlier entry recorded — the numbers moved in
  S7's forced S6+S7 re-run and the ledger text did not follow. The manuscript
  uses the tables.
- **D11 was done rather than asserted: all 56 panels were opened and read
  against their legends before the legends were written**, and two draft
  legends were wrong. *Nibea albiflora* has 55 spliceable introns of which 54
  are GT-AG and one is a minor site, not 55 GT-AG; and the non-vertebrate
  copy-number figure is drawn over the 193 *controlled* genomes, not all 194.
- **The PDF silently dropped the family's own name.** TeX Gyre Termes has no
  subscript glyphs, so `IP₃` typeset as `IP` on a page that otherwise looked
  right, with only a `Missing character` warning in a log nothing read. Fixed
  and generalised (D62); the log is now read and is empty.
- **All four build guards were tested by breaking each on purpose** — a
  missing figure, a missing section, a cited key with no reference row, a
  failing claim — and each sets a non-zero exit.

### Decisions

- **D61** — a manuscript's citations are stable keys resolved at build time
  and its bibliography is rendered from `references.tsv`, never typed.
- **D62** — a typeset build is not finished until its own log has been read
  for missing glyphs.

### Housekeeping

`s14_claims.py` reached 937 lines and was split three ways
(`s14_claims_scope/history/function.py`) to meet the project's 500-line
budget. `s14_lib.py`'s figure maps, deposit directory list and
`BULK_EXCLUSIONS` were rewritten from the plan the port shipped to what this
project actually produced; every regeneration command in `BULK_EXCLUSIONS`
was run to produce the data it regenerates.

### Emergent

Three rows: the discussion cites 29 of 137 curated references and the
remaining claims have not been audited the way S0 audited the baseline; the
15.2 % false-negative rate applies to every per-cell number in the paper and
only some are reported both sides of the contiguity bar; and 28 committed
publication figures — mostly the evidence for method decisions the paper
states in prose — are in the deposit but not the manuscript.

### Next

**S24 — supplementary figures and the figure audit.** `SUPPLEMENTARY_FIGURES`
in `scripts/s14_lib.py` is deliberately empty and the build fails on a
missing figure, so S24 fills that list and re-runs the chain. The figure
audit it also owns should read the panels at *printed* size in the assembled
PDF; S14a read them at source resolution, which is a different check.

---

## 2026-09-08 — S24: supplementary figures, and the figure audit

**Task.** S24 — the six supplementary figures showing the alignments and
structures the main figures rest on, and a figure-by-figure audit of every
main and Extended Data figure against its own legend (D11).

### The supplementary figures

Six, all drawn from committed files only, `results/supplementary/figures/`:
the representative alignment with the columns trimAl kept marked in the
input's own coordinates; the ligand core and the pore module at residue
resolution with every pathogenic position's residue printed for all three
paralogues; the per-paralogue deep alignments the constraint map is computed
on; the trimmed codon alignment; the constraint map painted on all three
cryo-EM channels with the selection layer beside it; and every labelled
variant with its per-element enrichment test. `SUPPLEMENTARY_FIGURES` in
`scripts/s14_lib.py` is filled and the manuscript now builds at 56 pages with
180/180 claims re-verified (five new claim rows, C176–C180).

**The two guards ran first and are hard failures (D63).** trimAl writes no
column map; S6 recovered one with `-colnumbering` and S24 does not trust it,
so every one of the 1,797 trimmed columns is compared against the input
column the map names over all 134 sequences. An off-by-one would still map
every column to a column and would renumber every residue claim downstream
with no other symptom. The residue joins were checked the same way: 1,780
variant residues against their own paralogue's table, 2,699 aligned partners
against the other paralogue's.

**The refusal that shaped Supplementary Fig. 6 (D64).** A structure carries a
human variant position only if every residue it shares with the human
per-residue table carries the same amino acid. Human ITPR2 (9YKK) and ITPR3
(8TKG) pass at 2,168/2,168 and 2,210/2,210. The ITPR1 cryo-EM reference is a
*rat* structure (736/2,300) and AlphaFold DB's human ITPR1 model is the
2,695-residue Q14643-4 isoform (479/2,695) — S11's isoform trap arriving in a
second place — so ITPR1's 55 pathogenic positions are not drawn.

### The audit

Every main and Extended Data figure was opened and read against its legend.
**26 findings: 16 legend corrections, 10 figure fixes.** All are in
`results/supplementary/figure_findings.tsv` with the committed table each
correction was re-derived from. The ones that mattered:

- **Fig. 3** — the legend said four backbone `100/100` labels; the tree draws
  **five**, and they are the two ancestors of each boxed clade, not the
  extended clades.
- **Fig. 2** — the colour key was mis-mapped (four blues, not two; the
  fragmentary class is the palest blue, not the grey), the class count was
  wrong (six of thirteen classes, 302 of 309 genomes), the missing `absent`
  colour was explained by the wrong rule (the ledger *does* hold those four
  cyclostome cells as absent — their classes have one genome each and are not
  drawn), and "Aves carry the most non-blue area" is Lepidosauria on the
  fraction the panel plots.
- **ED Fig. 7c** — "each of the 51–53 implied losses" for a panel that plots
  51 to **102**.
- **ED Fig. 10a** — "the 29-structure panel" for a panel with **30** bars.
- **ED Fig. 1b** — the bar the legend asks the reader to compare against was
  hidden behind the bar it was being compared with; now an open outline.
- **ED Fig. 8d** — two heat maps side by side on independent colour scales,
  so the one manufactured loss in panel a was drawn as dark as the 45 in
  panel b.
- Panel letters were uppercase in Fig. 7 and ED Figs 9, 12 and 13 and
  lowercase everywhere else; Fig. 6 carried a lone panel letter `a`.

The mechanical half of the audit now runs on every build (`s24_audit.py`):
every figure has a legend, every legend a figure, and each Extended Data
figure's legend letters match its panel files. It caught the uppercase Fig. 7
letters on its first run.

### The reproducibility fix (D65)

matplotlib stamps the wall clock into a PDF's `/CreationDate`, so every
figure in this project differed from its own rebuild by two bytes and no
SHA-256 recorded against a figure pdf meant anything. `figstyle.save` now
drops the field; 12 of 12 S24 files rebuild byte-identically. The eleven
figure modules touched this session were re-rendered under the fix; the rest
become reproducible on their next rebuild.

### Testing

17 negative controls run before anything is written, 5 mutation tests, all
caught. **Two initially passed on the broken code** and both are now recorded
in the test: the duplicate-column case has to be built where the two input
columns hold the *same* residues, or the content walk catches it first and
the one-to-one guard is never exercised; and a byte-identity check on a saved
figure passes vacuously whenever both saves land in the same second, so the
property is tested directly instead (the pdf must carry no timestamp).

### Housekeeping

`s24_figs_alignment.py` reached 529 lines and was split
(`s24_figs_inputs.py`), with `binned` moved to `s24_lib` so the two halves
cannot bin a profile differently. The reviewer checklist's open item 2
(supplementary figures not built) is closed; the remaining item is reading
the panels at printed size in the assembled PDF, which this session did for
the six new figures and not for the other 56.

### Next

**S14c — the manuscript rewrite pass.** Its dependencies (S14a, S24) are both
complete. Version the current draft rather than overwriting it; the PIEZO
project froze v1 when its framing changed and that turned out to be worth
doing. S21 and S22 remain pending in the analysis block.

---

## 2026-09-08 — S21: the gene itself, and what a 58-exon architecture is worth

Ledger row S21, gene architecture. Everything it needed was already on the
drive: the sweep's 309 retained `miniprot.gff` files carry one CDS record per
aligned block with its genomic interval, its span in the bait's own residue
numbering and its phase, the genomes are still there so splice dinucleotides
can be read rather than assumed, and the assemblies' gene sets are there too.
No new alignment, no new download, and the whole task rebuilds offline in 52
seconds.

### The result

**The IP₃ receptor is a 58-exon gene and its genomic span is not conserved.**
Over 1,378 genes in 189 genomes above D4's contiguity bar: ITPR1 58 exons in
147 kb, ITPR2 57 in 244 kb, ITPR3 58 in 58 kb, with coding lengths of
8,250 / 8,100 / 8,008 bp. A 4.2-fold span spread at a 3 % coding-length
spread, and the ordering is consistent gene by gene inside genomes rather
than an artefact of averaging (D16 paired sign tests, BH-corrected across the
family; ITPR3 is the shorter gene than ITPR1 in 161 of 182 genomes carrying
both). The RyR control, measured through the identical instrument in the same
assemblies, is 104 exons over 14,910 bp — nearly twice the gene in both, at
the same mean exon length.

**The three paralogues share ~48 of their ~58 intron positions**, against
0.55 expected, in every one of 181–183 genomes tested. **The ryanodine
receptors share one**, in 0 of 188. The sister family that carries every
ITPR-diagnostic Pfam domain — the hazard this project is built around — has
an exon structure with no ancestry in common with this one.

**Database "fragments" are annotation failures, not gene boundaries.** Of the
291 loci S18 called `split` or `fragmentary`, 228 carry at least one annotated
model terminus sitting *inside an exon* of the gene model, where nothing
splices. Three are broken entirely at real junctions.

### The instrument, corroborated twice

99.89 % of 112,254 junctions read as a canonical or minor splice pair off the
genome; 0.78 % carry a frame step. And 94.5 % of 188,146 annotated CDS block
edges from an independent pipeline land *exactly* on a sweep exon boundary,
over 1,171 loci in 164 genomes — S10's two-case check generalised to the
scope, with D9's contrast surviving it (RefSeq 94.9 %, GenBank 91.1 %).

### Three things the brief got wrong, and what replaced them (D66–D68)

**The frameshift-pair merge has nothing to merge.** The brief expected
miniprot to emit two CDS records either side of a frameshift. Measured over
all 149,148 consecutive block pairs in the sweep, no such pair exists: query
spans are contiguous across every one, the smallest genomic gap anywhere is
10 bp, and it reads `GT..AG`. An indel appears *inside* a block. The
calibration therefore refuses to derive a bar — one of its two populations is
empty — and the floor goes where the evidence is (10 bp) rather than at the
declared 30 bp fallback, which would have merged ten junctions the genome
calls splice sites. The merge fires on 0 of 112,254 junctions and T5
constructs a 2 bp pair to prove it can act (**D66**).

**The tandem-duplication test measured paralogy.** "The same bait aligns
twice at disjoint positions" fires in essentially every vertebrate genome,
because the three paralogues are 61–68 % identical and every bait aligns at
all three genes: sensitivity 0.997, **specificity 0.16**. Scoped to the cell's
own loci with `s5_classify.cell_loci` imported unchanged — D14 applied to a
pairwise test — specificity goes to 0.977 against S16's copy call, which the
detector never sees, and the 3R teleost check agrees on 261 of 267 cells
(**D67**). No locus in the sweep encodes the same part of the protein twice.

**The fragment test had the wrong unit.** Scoring every annotated CDS block
edge answers the boundary-concordance question a second time — a model's
internal boundaries are its own splice sites and canonical by construction.
The unit is the annotated *model's terminus*, with the gene's own ends
excluded by rule because a real gene legitimately starts and stops inside an
exon. That changed the answer's shape: 44 % of internal termini land on a
boundary the gene model has, not 81 %.

### The frame, and the anchor that could have failed

Boundaries are in the bait's numbering and the sweep used 38 baits, so
everything travels bait → its cell's human reference (pairwise MAFFT at
`--thread 1`) → a column of S6's committed alignment. The frame is checked
rather than assumed: all 14 residues S0 measured on 6DQN — the ten IP₃
contacts, the two filter and the two gate residues, in each paralogue's own
numbering from S17's `functional_sites.tsv` — land in the **same column** in
all three paralogues. T8 shifts one row by a column and requires the gate to
refuse.

### Testing

21 negative controls run before anything is written, split across
`s21_test_arch.py` and `s21_test_claims.py` to stay inside the file budget.
Six mutation tests, all caught. Three controls failed on their first run and
two of the three were the rule's fault, not the test's: the terminus verdict
had no value for "the annotation ends inside an intron of the model", which
is the two pipelines disagreeing rather than one of them inventing a
boundary, and the Poisson-binomial tail returned 0.9999999999999991 where the
answer is "nothing was measured".

### Housekeeping

Two joins were guarded rather than trusted, and one guard earned its place
immediately: `loci_with_mp()` requires the contig and start of every S16
locus to match the summary it takes `mp_id` from, and found that S16 files
two `vertebrate_basal` loci out of `other_loci` rather than out of a cell.

### Next

**S22 — ligand-site evolution.** Its dependencies (S9b, S17) are complete. The
alignment frame and the anchor test are reusable as they stand, and S17's
per-residue table already carries the measured IP₃ contacts. Expect the
blocking step to be scope: the lineages that lost the upstream PLC/IP₃
pathway are not enumerated anywhere in this project yet.

---

## 2026-09-09 — S22: ligand-site evolution

**Task.** The one module the ryanodine receptors do not share functionally is
the IP₃-binding core. Ask what evolution did to it: against the pore, against
its own contact residues, and in lineages that lost the upstream enzyme.

### What ran

`scripts/s22_*.py`, ten ordered stages behind `s22_run.py`
(`modules → shells → paired → contacts → omega → plc → lineage →
deep_lineage → tables → figures → report`). Everything but `plc` and
`deep_lineage` is offline and rebuilds byte-identically in under a minute;
`plc` is 10 `hmmsearch` runs over 33 GB of reference proteomes (~25 min) and
`deep_lineage` is 668 pairwise MAFFT alignments (~6 min). 32 committed
tables, 4 figures, `results/ligand_site/`.

### The blocking step was scope, and it was not blocking

The previous session expected the lineage list to be the problem. It was
derived rather than read: PI-PLC presence — a protein carrying **both**
halves of the catalytic TIM barrel, PF00387 and PF00388 — swept over all
3,527 eukaryotic reference proteomes with S20's design pointed at a different
profile. 760 of 763 vertebrate proteomes carry one, which is the positive
control for the search. **64 proteomes carry an ITPR and no PI-PLC.**

### Results

1. **The pore is more conserved than the ligand core.** Paired per
   orthologue — one core number and one pore number per sweep orthologue,
   246-262 per paralogue — the pore leads by 0.024 identity in ITPR1
   (223 tips to 32, q = 3.6e-34) and 0.020 in ITPR3 (218 to 42, q = 1.6e-28),
   with ITPR2 flat.
2. **And that reverses on one boundary.** Leave the 50-residue luminal loop
   inside PF00520, which is how InterPro draws it, and all three paralogues
   flip to core > pore at q < 1e-37. Both answers are right about their own
   region; neither is right about "the pore" (D69).
3. **The metrics disagree and the composition-free one wins.** Per-column
   JSD sees no difference between the modules; `frac_modal` sides with the
   paired test. JSD is a divergence from a background amino-acid table, so a
   transmembrane module scores low at equal conservation — S17 measured that
   and this is where it bites.
4. **The constrained unit is the pocket, not the contacts.** The ten measured
   contacts beat the rest of the binding core (q = 0.048 / 0.017 / 0.041) and
   beat the rest of the 15 Å pocket in none of the three. Every shell out to
   15 Å is above the whole-protein mean and there is no step at 4.5 Å.
   FEL agrees from the other side: 100 % of contact sites purifying in all
   three paralogues, the share falling 0.13-0.25 across the shells.
5. **The lineage question has an answer, and it is a bounded null.** Pooled,
   the 64 taxa's ligand core looks relaxed (p = 9.7e-6) — and that is a clade
   artefact: they sit at median pore identity 0.358 against 0.589 for the
   rest, and the paired statistic is itself correlated with divergence.
   Matched to PLC-present records within 0.03 pore identity, all 36 that
   enter the test find controls and the effect is gone: median within-pair
   difference −0.0064 (95 % CI −0.016 to +0.015), 19 to 17, p = 0.87 (D71).

### The instrument

Both modules are defined by measurement rather than taken from Pfam, twice
each, and each definition is checked against something it does not contain —
a failure raises (D69). Every residue within 15 Å of IP₃ is measured
**all-atom** in **six** independent IP₃-bound human ITPR3 depositions, with
S0's ten contacts recovered in the structure S0 used as a hard-failure
positive control (D70). The consensus contact set is **twelve**: Ala276 and
Arg411 are inside 4.5 Å in a majority of the depositions and outside it in
6DQN. The RyR positive control is measured through S22's *own* pairwise
instrument at the *same* divergence as the test group, which is what makes
the null in item 5 bounded rather than empty.

### Testing

44 constructed negative controls run before anything is written, split
across `s22_test_ligand.py` and `s22_test_lineage.py` to stay inside the file
budget, with one entry point. **Ten mutation tests, all ten caught.** Two
mutations were missed on the first pass and both exposed a real weakness
rather than a bad mutation: T1-T3 read the committed `module_map.tsv` instead
of calling the rules, so a builder that widened a module passed every check —
fixed by exercising `M.build()` and adding T3b, which compares the committed
map against what the rules produce now. The transfer refusal fires on nothing
in this data, so T13b makes it fire.

### Next

**S14c — the manuscript rewrite pass.** It is the last unblocked pending row
(S14b is human-gated). S22 hands it three things the draft does not have: the
core-versus-pore boundary problem, which any figure or sentence about "the
pore" now has to declare; the twelve-residue contact set against the ten the
draft quotes; and a lineage result that is a bounded null rather than a
caveat.

### Addendum (same session) — the dashboard was reporting S12 as pending

Asked whether S12 needed more work. It does not: 18 committed tables, its
self-test passes, and its report renders every section. But the **dashboard**
was showing it as `pending`, and the cause is worth recording because it is a
protocol hazard rather than a display bug.

S12's Results cell quotes a shell pipeline — `` `fastq-dump | hisat2` ``.
That pipe is legal inside a markdown table cell but `dashboard.py` split the
row on every `|`, so S12 came out with 6 cells instead of 5, was read as the
*analysis* table's layout (which has an extra Priority column), and its
status was taken from the Results prose. The prose does not contain the word
"completed", so the row fell through to `pending`.

Why that matters: session-protocol step 4 picks the first `in_progress` row,
else the topmost unblocked `pending` one. A session driven from the dashboard
rather than from the ledger text would have started **S12** this morning
instead of S22, and redone finished work. Nothing failed; the number was
simply wrong.

Both halves fixed. The pipe is escaped in the roadmap (`\|`), and
`parse_ledger` now splits on **unescaped** pipes only and locates the status
cell **by pattern rather than by position** — a positional read cannot tell a
missing Priority column from a shifted row. A row with no recognisable status
is now skipped with a message on stderr instead of silently becoming
`pending`. All 33 rows re-parse correctly: only S14b and S14c are pending,
and both legitimately.

---

## 2026-09-09 — S14c: the manuscript rewrite pass

The brief for this row is four questions. Does the paper still lead on its
strongest result; are the figure numbers still right; does every claim row
still pass; and version the previous draft rather than overwriting it. Three
of the four turned up something.

**What the draft did not know.** S14a assembled the package on 2026-09-08 and
S21 and S22 both landed after it. Neither analysis was in the paper at all —
no Results text, no figure, no legend, no claim row, and neither results
directory in `DEPOSIT_DIRS`, so their tables were not even being deposited.
That makes this a rewrite rather than an edit, so the S14a draft is frozen in
`manuscript_v1/` before anything changed, with `FROZEN.md` recording why and
how to recover its figure set from the committed manifest rather than by
storing a second copy of 17 MB.

**The new section.** S21's result belongs next to the origin section, not in
the annotation section where its fragment half sits: the three paralogues
share a median 46–49 of about 58 intron positions in every one of 181–183
genomes at 85–88× chance, and the ryanodine receptors share one, in none of
183–188. Intron position is a character no protein alignment produces, which
is exactly what makes it worth having — the ITPR/RyR separation that every
search in this project is built around is here made by evidence that shares
nothing with the evidence that made it everywhere else. S22 folded into the
machine section and both into the abstract, discussion, limits and methods.

**The figure numbering was wrong, and not only because of the two additions.**
Extended Data Figs 2 and 3 were cited by no sentence in the paper, and the
methods figure was numbered 14 and first cited in the third Results section.
Renumbered to 16 figures in strict order of first mention, and the check that
found it is now code: `s24_audit.citation_order()` requires every Extended
Data figure to have a legend, be cited, and be cited in order. Mutation-tested
four ways — a swapped citation, a removed one, a citation with no legend, and
a forward reference that breaks only the order.

**Three defects in the build itself.** The deposit list omitted two completed
tasks. The reviewer checklist claimed the build reads the LaTeX log for
`Missing character` and it did not — pandoc reports only its own warnings
unless asked — so `s14_pdf.py` now runs `--verbose`, scans the log and exits
non-zero on a dropped glyph; mutation-tested on two glyphs the document font
lacks, and currently zero. And the checklist itself, the last hand-typed
artefact in a package where everything else is generated, was four editions
stale at 175 claims against a ledger of 180, so `s14_claims.py` now fails the
build unless the checklist's headline count is the ledger's.

**One number was wrong in my own new text and the table caught it.** I wrote
that the share of sites under purifying selection falls across the ligand
pocket "out to the furthest shell". It does not: it falls to each paralogue's
*third* shell and the outermost sits slightly above it in all three, so the
pattern is a step onto a floor rather than a gradient. S22's report says so
explicitly; I had read the ledger summary instead of the table. Corrected in
the text and pinned by three claim rows.

**The lead was re-examined and left alone.** The title's two halves — retained
in every vertebrate, lost repeatedly elsewhere — are still the two strongest
results and the two the design was built to be able to make. Both new
analyses are about what the gene and the protein are *like*, not about where
they are, and neither displaces the census. The abstract gained a sentence on
each.

Claims 180 → 276, all passing. 17 sections, 17,807 words, 60 pages, 7 main +
16 Extended Data + 6 Supplementary figures, 1,937 deposited files. The build
is green end to end and every one of its guards has now been broken on
purpose at least once.

Two emergent rows. `figstyle.save()` checks the size a figure *declared*
rather than the size the tight bounding box actually wrote, so three panels
are saved wider than the text block and silently scaled down when placed —
a general fault, not a fault in those three figures. And the new Results
section asserts that the two families' exon structures have no common
ancestry without saying how that came about, which is a parsimony job on the
intron characters of the kind S15b already runs.

**Two tasks added at the user's request (2026-09-09, after S14c closed).**
The observation behind them is that the manuscript reads as a compressed
thesis rather than as a single paper, which is true of the numbers: 17,807
words and 29 references stand on 109,243 words of committed task reports, 68
recorded decisions and a 137-reference review.

**S25 — the thesis.** The long form. Its point is not length but the three
things the paper had to drop: why each instrument is built the way it is (the
Decisions log has never been written as prose), what was measured and
abandoned, and the ~300 constructed negative controls as a body of work. The
constraint that makes it a real task rather than a reformat is the
bibliography: S0's audit rule applies unchanged, so every new reference is
audited on entry and committed as a table. "Far more references" must not
become "far more references nobody checked". The claims ledger extends to it,
because a longer document is a larger surface for drift.

**S26 — the paper series.** The grouping is the deliverable and is derived
from six rules committed before the assignment, not handed down. Five of the
rules are ordinary (one question, its own controls, a declared scope, four to
seven figures, no result primary in two papers); the sixth is the one that
matters — *what does this paper claim if none of the others is ever
published?* — because it is the only honest test of a series against a slice.
A starting proposal of five papers is in the brief, offered to be revised by
the rules rather than instead of them. The S14 machinery generalises from one
package to N, but the claims ledger stays single with a paper column, so a
number quoted in two papers cannot disagree between them.

The two must not contradict each other on how the results group. They are
independent and either may run first; whichever does commits the assignment
table, and the other adopts it or records why a chapter and a paper are not
the same unit.

---

## 2026-09-09 — S25: the thesis, and a bibliography that had to be checked

### What ran

`python scripts/s25_assemble.py` — eight stages, exit 0. The document is
**58,436 words across 15 chapters and 5 appendices, 103 figures, 80
references, 173 typeset pages**, in `thesis/`.

The order of work was: the chapter grouping first, committed before any prose
(the brief's step 1); then the build machinery with every guard; then the
reference audit; then the chapters; then the claims ledger, which was written
last because it can only be written against text that exists.

### The grouping, and the one rule that had to move

Seven rules, in `thesis/chapter_rules.md`, four of them enforced by the build
rather than asserted. T6 is the interesting one and it is enforced in the
strong direction: **every entry under `results/` is either assigned to a
chapter with a rule and a reason, or named in an exclusion list saying why it
is not a result.** An unassigned directory fails the build, which also means a
future task's results cannot be silently left out. 37 assigned across 13
chapters, one excluded (the dashboard's live panel).

T3 was written as "four to eight figures" and widened to four to twelve
**before any prose was written**, with the reason stated: a paper's figure
budget is set by a journal and a chapter's by a reader's attention. Chapter 7
(neighbourhood, reconciliation, duplication) is one argument and needs twelve.
The introduction, methods and discussion are exempt as exposition.

### The reference audit caught 9 of 58

This is the result of the session. The rule is S0's, unchanged: a reference
added and not audited is worse than no reference. What makes it enforceable is
that **nothing bibliographic is typed.** Each new reference is declared by an
identifier alone plus a distinctive phrase its title must carry; the build
resolves the identifier live, admits it only if the phrase is there, and
writes the bibliographic row from the fetched record.

Nine identifiers written from memory resolved to entirely different papers. A
duplication-inference algorithm's returned MrBayes. A reconciliation method's
returned a paper on statistical challenges in real-time PCR. A morphological
likelihood model's returned a paper on species names in phylogenetic
nomenclature. A vertebrate ancestral-genome reconstruction's returned a paper
on structured RNAs. Every one would have entered a bibliography looking
completely normal.

One of the nine was the other failure: the identifier was right and the
*phrase* was wrong, corrected by adjusting the phrase.

The rule then runs in the other direction too. An audited reference that is
never cited is decoration, and that also fails — which is what forced all 58
into the text rather than 26.

### The ledger, closed in both directions (D72, D73)

220 numbers verified. **151 are carried from the manuscript's ledger through
the same engine rather than a fork** — `s14_claims.check` was extracted from
`run()` for the purpose — so a number quoted in both documents is recovered
once and cannot disagree between them. That is S26 step 2's rule applied a
task early.

The new half is D73: **a claim's value must also appear in the chapter that
declares it.** The manuscript's ledger fails when the document states a number
no table produces; it cannot fail when the ledger declares a check the
document never makes. That guard fired on seven rows on its first run — five
numbers the thesis had not actually stated and two claims pointed at the wrong
table.

### Every guard broken, on every build (D74)

The brief asked for each guard to be broken on purpose once.
`scripts/s25_test_guards.py` does it on every build as the first stage: 15
cases, each declaring a fragment the guard's own message has to contain, run
against a sandboxed copy of `thesis/`, with the SHA-256 of every committed
file checked before and after (D60).

Its own first run had four cases returning non-zero for the *wrong* reason —
two path bugs in the harness, one guard reading a module-level constant the
sandbox could not redirect, and one mutation applied to a list the driver had
already copied. That is exactly what the message-fragment requirement is for.

### The one fix outside the thesis (D75)

`s0_report.py` computed "the wider bibliography extends this to N references"
as the size of `references.tsv`. Correct only while the review's bibliography
and the shared table were the same set — and S25 added 58 references to that
table, so the next render of a committed report would have silently restated
137 as 195. The count now comes from the keys `docs/review/*.md` actually
cites. D13 protects a report from a stale number, not from one that was always
measuring something adjacent.

### Testing

- 15/15 build guards fire with their own message; 243 committed files
  unchanged by the suite.
- 220/220 claims verified.
- `s14_assemble.py` still exits 0 (276/276) after the engine extraction, and
  `s0_review_build.py --check` still resolves all 137 review citations.
- Every new module is under the 500-line budget; the longest is 325 lines.

### Next

S14b (human-gated: Zenodo DOI, repo public, preprint) and **S26, the paper
series**, which starts from `thesis/chapter_assignment.tsv` and has one
declared departure to resolve: the thesis puts all of S19 in the chapter that
builds the search, and S26's proposal splits it across two papers. Appendix E
states it in both directions rather than picking.

### Addendum (same session) — the typeset PDF was not reproducible either

Re-running the build after committing showed the thesis PDF changing by two
bytes: an embedded build timestamp, which is the defect D65 fixed for figure
PDFs one level down. Setting `SOURCE_DATE_EPOCH` in the pandoc call removes
it.

It does not make the file byte-identical, and that is stated rather than
claimed away. Two builds of the same document now differ in **exactly 64
bytes of 3.5 million** — two copies of a random 16-byte trailer `/ID` that
this `xdvipdfmx` writes regardless of `FORCE_SOURCE_DATE`. So a checksum
recorded against a typeset PDF still means nothing, and `s14_pdf.py` has the
same defect without even the timestamp fix. Recorded as an emergent item.

### Addendum (same session) — the thesis was rewritten for register

The user read the draft and reported that the prose was confusing: headings
and openings were sentence fragments that gave no context on their own, and
em-dashes were overused. Both are fair, and both are the wrong register for a
scientific document rather than a matter of taste.

**All 27 chapter files were rewritten.** Every heading now names its subject
without needing the one above it, so that it is legible from a table of
contents or a search result. Every section opens by naming what it is about
rather than pointing at it with a pronoun. And there are **zero em-dashes in
the chapter sources**, down from 464, with each one replaced by the
punctuation or the sentence break the sense actually wanted rather than by a
mechanical substitution.

Some examples of what changed. The title went from *a genome-scale census and
its instruments* to *A genome-scale census of the inositol
1,4,5-trisphosphate receptor family, and the methods built to make it*.
*The problem stated precisely* became *Why the separation has to be a positive
test rather than a filter*. *The claim that failed, and what it turned into*
became *The background claim about gene size was false, and it became a
measurement*. *Completeness, asked in the expensive direction* became
*Completeness checked in the expensive direction: which known records the
profiles missed*.

Word count rose from 58,436 to 60,018 and the typeset document from 173 to 177
pages, which is what saying a thing plainly costs.

**The build passed unchanged throughout.** 220 claims, 103 figure placements,
80 citations and all 15 guards survived a full-prose rewrite of the entire
document, checked after every file. That is worth recording for its own sake:
none of the document's checkable content lives in its prose, so the prose
could be replaced wholesale without touching a single number. The one thing
the guards did catch was a citation dropped when a paragraph was reworded, and
the "audited but never cited" rule fired on it immediately.

Recorded as **D76**, which states the three rules and applies them to S26's
papers as well.

### Addendum (same session) — title, authorship, and a second grammar pass

Three corrections after the user read the rewritten draft.

**The title.** It had kept the appended-fragment pattern D76 removed
everywhere else, in the one place a reader sees first. It is now *A
genome-scale census of the inositol 1,4,5-trisphosphate receptor family*, with
the scope moved into a subtitle.

**The authorship.** The thesis is authored by Claude, with the correspondent
named separately as the point of contact. Both now live in `s25_lib` and are
read by the PDF stage rather than typed into it, so the front matter and the
PDF metadata cannot disagree. Recorded as **D77**.

**A second grammar pass**, which is the part worth recording. The first pass
fixed the headings that were obviously confusing and removed every em-dash. A
systematic check found more:

- **30 section headings still had no finite verb**, being participles
  (*A domain-architecture call, audited against...*), gerunds (*Enumerating
  the search space, and...*) or bare noun phrases (*The loss count*). All are
  now complete clauses.
- **Fifteen bold paragraph lead-ins were fragments** (*Its taxonomic range.*,
  *The model-violation guard.*, *Settled.*, *The neighbourhood.*). All now
  carry a predicate.
- **Appendix A's five class labels** were noun phrases with counts attached;
  they are now sentences.
- **One passage listed the figure-audit findings as five verbless fragments.**
  It is now five sentences.

The checks used are worth keeping for S26. Every heading was parsed for a
finite verb, every bold lead-in likewise, and every body sentence scanned the
same way, with the ten hits inspected by hand and nine confirmed as irregular
verbs the check does not know. Dashes were counted by type: zero em-dashes
remain in the chapter sources, and the only en-dashes are a numeric range and
a two-name compound, both of which are correct usage rather than the
parenthetical dash.

60,201 words, 177 pages. The build passed unchanged throughout.

### Addendum (same session) — a chapter on how the project was carried out

The user asked for a description of how Claude Code was used: a short section
in the front matter, and a full chapter at the end covering how the project
organised itself, surveyed the literature, decided what to download and
analyse, ran the analyses unattended, and produced the results, figures and
documents.

**Chapter 16 is that chapter**, in two files, sitting after the general
discussion and before the appendices, which were renumbered to make room. It
covers the one-task session protocol and why it exists, the claim-by-claim
literature audit, the download scope decided by rule and committed as a
manifest, the analysis menu written before the analyses ran, the resumable
stage drivers that made unattended operation possible, where the agent's
judgement actually went (into measuring thresholds), how tables, reports and
figures are generated rather than written, how the three documents are built
and guarded, what the human collaborator contributed, and what the arrangement
is bad at.

**Its numbers are measured rather than recalled**, which is the part worth
recording. `s25_production.py` is a new build stage that derives them from
files under version control: 35 sessions across 9 days, 35 ledger tasks of
which 33 are complete, 82 recorded decisions, 372 analysis scripts, 455
committed tables, 111 committed figures, and 28 rendered reports totalling
109,250 words. Twelve are declared in the claims ledger, so a sentence about
the project's scale fails the build when the scale changes.

Two things were deliberately kept out of the ledger. Total lines of code and
the commit count change on every edit and every commit, so the chapter states
them qualitatively and the committed table holds the exact figures.

The chapter's closing section is about the failure modes of the arrangement
rather than its strengths, because those are the useful part: a fluent agent
produces plausible prose about work it has not checked (the reference audit
caught nine of fifty-eight citations); a rule that fires on nothing looks
exactly like a rule that cannot fire (which is why every zero here has a
constructed control); and a style pass that has not been checked mechanically
only fixes what somebody noticed (which is what the second grammar pass
found). All three are things this project did, caught and recorded.

Recorded as **D79**. 64,377 words, 185 pages, 232 claims, build exit 0.

The ledger then caught its own chapter twice while it was being written: adding
D79 moved the decision count from 82 to 83, and committing the new statistics
table moved the table count. Both sentences were corrected because a check
refused the document, which is the behaviour the chapter describes.

### Addendum (same session) — legends that say why, and a discussion that cites

The user read the assembled thesis again and asked for two changes: every
figure legend to explain the importance of what its figure shows, and §15.1 to
discuss the importance of each new result with references backing the claims.

**All 103 figure legends were expanded.** Each already described the figure
and defended its construction, which is what a legend in this project's
per-task reports does, and none said why the panel was worth drawing. Each now
closes with what it changes for a reader: the coverage figure with the ceiling
it sets on every structural argument, the exon tracks with why the coding
fraction is the point, the ablation panel with what it says about how to spend
a bait-panel budget. The importance sentence goes last, because that is where
it is read last.

**§15.1 went from 409 words to 1,893**, restructured into six subsections, one
per new result, each stating the result and then what it changes about a
standing position. Three of them are arguments rather than observations once
the literature is beside them: the land-plant absence answers a question the
plant literature has left open for two decades, the two duplication splits
fall on the two branches the genome-scale reconstructions assign to 1R and 2R,
and the teleost retention of ITPR1 is the same dosage statement as its human
deletion phenotype made 300 million years apart.

**Eighteen new references, and all eighteen passed the audit first time.**
They were declared by identifier plus a title phrase and resolved live, the
same rule that caught nine wrong identifiers when the bibliography was first
built. The count of added references went 58 to 76, which moved two claims and
four sentences across three chapters and the appendix; each was corrected
because the build refused the document.

Recorded as **D80**. 74,184 words, 194 pages, 232 claims, 15/15 guards, build
exit 0.

## S27 — the editorial pass: every sentence read in context

The user asked for the thesis to be read sentence by sentence, in context, and
held to four standards: factually correct, grammatical, meaningful, and useful
to the reader. They also asked whether it meets the demands of clear
structure, an evidence-based approach and precise communication, checked
against IMRAD.

Every build guard passed before this pass began and passed unchanged after it.
The pass still found five classes of defect, none of them reachable by a
checker, because each is a statement that is well formed, internally
consistent, and wrong or useless in context.

**Repetition, 50 legends of 103.** D80 asked every figure legend to explain
why its figure matters. That was satisfied by appending a sentence to each
legend rather than merging one in, so half of them now stated one point twice,
either inside the legend or against the paragraph beside it. Each half read
well alone, which is why it was invisible while writing. All 50 were rewritten
as single paragraphs.

**One number wrong at source.** `s18_report_results.py` had 764 vertebrate
reference proteomes typed into it, which is the line count of the sweep's
manifest including its header. The sweep's own statistics file says 763, and
so do three other chapters, so the package had disagreed with itself by one
for four editions. The generator now reads the number, and the audit report
was re-rendered.

**Four stale counts.** 75 decisions against 84, 26 sessions against 35, 58
added references against 76, and a chapter-16 anecdote about a correction that
has since happened twice. Two new claims, T83 and T84, now pin Chapter 14's
figures to `production_stats.tsv`.

**Two wrong cross-references** in Chapter 6, both pointing a reader at the
section before the one that answers the question.

**Six results chapters had no citations at all**: 5, 7, 8, 12, 13 and the
search-worth half of 4. Chapter 5 reports the range result, which bears
directly on a plant-literature question two decades old, and cited nothing.
Each now carries the literature its claims bear on, taking cited references
from 111 to 116.

### What was added rather than fixed

The central research questions are now stated explicitly in §1.4 rather than
left implicit in four gap statements. §14.1 carries a **methods-location map**,
because methods distributed by design still have to be findable, and a reader
wanting a procedure rather than an argument had no way in. The front matter
carries a **data and code availability** section and a practical-implications
paragraph closing the abstract. N50 and the Bayesian information criterion are
defined at first use, which were the only two unexplained abbreviations left.

### Three of the five classes are now checks

`scripts/s25_prose.py` is a new build stage: no em-dash in a chapter source,
no legend that says the same thing twice, no legend that restates the
paragraph beside it. Three guard cases were added to break each on purpose, so
the suite is 18 rather than 15. The other two classes, wrong cross-references
and missing citations, remain a reading problem and are recorded as such.

Recorded as **D81**, **D82** and **D83**. 73,083 words, 194 pages, 234 claims,
18/18 guards, 0 prose failures, build exit 0.
