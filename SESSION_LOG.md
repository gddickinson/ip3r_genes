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
