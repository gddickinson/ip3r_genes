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

**Next.** The two remaining jackhmmer runs (protista, metazoa) and
then S23, which S20 sharpens: the absences to take to genome level are
Streptophyta, Dikarya, Glomeromycota/Mortierellomycota, Apicomplexa (0/60) and
the *Cymbomonas* copy-number question.
