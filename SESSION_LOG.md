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
