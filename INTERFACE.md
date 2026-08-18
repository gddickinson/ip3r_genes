# INTERFACE.md — Module Navigation Map

> Read this **before** opening any source files. Each `.py` is < 500 lines
> and has one focused responsibility.

## Top-level

| File | Purpose |
|------|---------|
| `run.py` | CLI/GUI entry — `python run.py [--email you@example.com] [--headless --preset ip3r]`. |
| `PUBLICATION_ROADMAP.md` | **The multi-session publication plan**: session protocol, task ledger (S0–S24), emergent tasks, and the Decisions log — including the methodological rules inherited from the PIEZO project. Read at the start of every session. |
| `docs/session_briefs.md` | Detailed per-task instructions for every ledger row: goal, steps, completion criteria, outputs. |
| `docs/ip3r_background.md` | **The biology baseline.** Every statement tagged `[db]` (verified against a live database, with the query), `[lit]` (literature, verified in S0) or `[open]` (a question this project answers). Includes the ITPR/RYR hazard and the database-scale snapshot. |
| `docs/ip3r_review_2026.md` | **The literature review** — a 24-page, 137-reference review of the family (architecture, gating, regulation, cell physiology, paralogues, evolution, genetic models, disease, pharmacology, open questions), with §12 carrying the S0 claim audit. **Generated — never hand-edit it**; edit `docs/review/*.md` and re-run `scripts/s0_review_build.py`. Cite its sources, never `ip3r_background.md`. |
| `docs/ip3r_review_2026.pdf` | The typeset review (A4, pandoc + xelatex), built by `s0_review_build.py --pdf`. |
| `docs/review/` | The review's **source**: numbered section files `00_frontmatter` … `12_methods_audit`, each under the 500-line limit. Citations are written as stable keys (`[R22]`, `[R22, R25]`) resolved against `results/s0_baseline/references.tsv`; the build renumbers them into order of first appearance. |
| `docs/analysis_catalogue.md` | **What the harvested data can answer** — the analysis menu behind ledger rows S15–S22, each with its inputs, method, deliverable and invalidating caveat, plus the report skeleton mapping analyses to manuscript sections. |
| `CLAUDE.md` | Session-start instructions, conventions, the family definition rule. |
| `INTERFACE.md` | (this file) navigation map. |
| `README.md` | User-facing docs + the project status board. |
| `SESSION_LOG.md` | Running notes per session: what ran, what resulted, what's next. |
| `FINDINGS.md` | **The biological story, task by task** — a dated plain-language entry appended after every completed task; unconfirmed claims marked *(pending: which task confirms it)*. |
| `roadmap.md` | Feature-development history of the app itself (inherited v1.0–v1.7 plus this project's changes). |
| `data_root.txt` | Path to the bulk-data root, read by `src/utils/data_root.py`. Edit one line to move storage. |
| `requirements.txt` | Third-party deps (`biopython`, `requests`). |
| `presets/` | Bundled JSON queries: `ip3r`, `ip3r_zebrafish`, `ip3r_discovery`, `ip3r_paralog_mine` (Compara paralog mine), plus the special-mode presets `ip3r_domain_scan` (`"mode": "domain_scan"`) and `ip3r_all` (`"mode": "exhaustive"` — full harvest + every analysis). Presets may carry `known_paralogs` and a `discovery` block of scorer-threshold overrides, honoured by both the CLI and the GUI Discovery dialog. |
| `cache/` | On-disk JSON cache (auto-created, gitignored). |
| `results/` | Committed analysis output — one directory per task. `s0_baseline/` holds the S0 database snapshot, the literature audit (`lit_claims.tsv`, `references.tsv`), the Ensembl endpoint probe and the rendered `report.md`. |
| `manuscript/` | The submission package, built by `scripts/s14_assemble.py` in S14a. See `manuscript/README.md`; nothing here is hand-edited except the numbered section files. |

## `scripts/` — roadmap-session tooling (not part of the app)

Ported and ready to use from session one:

| File | Purpose |
|------|---------|
| `dashboard.py` | **Project dashboard generator** (stdlib-only) → `dashboard.html` (gitignored): overall + per-task progress parsed from the roadmap ledger, a live panel for the in-progress task driven by `results/session_live.json`, key figures (explicit list first, then auto-discovered from `results/*/figures/*.png`, base64-embedded), and FINDINGS.md rendered. Meta-refreshes every 30 s. Session protocol step 0: `--open` plus a background `--watch`. |
| `figstyle.py` | **The one figure style every publication figure goes through.** Page geometry imported from `s14_lib`, so a figure is drawn at the width it is placed at. The validated palette: `PARALOG` (the fixed ITPR1/2/3 trio), `GROUP` (adds the non-vertebrate grade and the RyR outgroup in accent violet), `STATUS`/`QUALITY` (the ordered blue→grey→red evidence scales), `CLINICAL` (reserved status colours). `panel()`, `despine()`, `hgrid()`, and `save()`, which writes png+pdf at 400 dpi and **raises** rather than writing a figure whose tight bbox runs off-canvas or whose text needs a glyph the font lacks. |
| `s14_lib.py` | Manuscript-assembly shared helpers: page geometry, the section order, the main / Extended Data / Supplementary figure maps (the *plan* — source paths are where the tasks write their figures), the deposit rules, and stdlib TSV/JSON/SHA-256/PNG-size helpers. |
| `s14_figures.py` | Copies every committed figure into `manuscript/figures/` under its publication number, in every format the analysis produced. Nothing is re-plotted, so a manuscript figure cannot differ from the one in its results directory; files not in the manifest are deleted on every build. |
| `s14_claims.py` | **The manuscript's drift guard.** Every load-bearing number is declared with its source table and the op that recovers it (`count` / `cell` / `nunique` / `sum` / `json` / `grep`, with `col__startswith` and `col__in` filters). Ships with an empty `CLAIMS` list and a worked example of each op. |
| `s14_deposit.py` | Walks the deposited tree → `deposit_manifest.tsv` with size + SHA-256 per file, and `deposit_notes.md` listing the bulk classes deliberately excluded **with the command that regenerates each**. |
| `s14_pdf.py` | The typeset PDF: builds a pandoc markdown with each figure placed above its own legend, then xelatex via pandoc. Three traps are handled in code and commented there. |
| `s14_assemble.py` | The driver — ordered stages `figures → claims → stitch → deposit`; `--only` / `--from` / `--list`; non-zero exit on a missing figure, a missing section or a failed claim. |
| `s0_db_snapshot.py` | **S0 step 4** — re-derives every `[db]` number in `docs/ip3r_background.md` (InterPro signature + taxonomy counts, the UniProt reference and sister-family panel, per-protein Pfam architecture, the zebrafish PF08709 query) into `results/s0_baseline/*.tsv`. Re-run it whenever a `[db]` number is about to be quoted. |
| `s0_gene_structure.py` | **S0 step 4** — exon counts, genomic spans and cytobands for ITPR1/2/3 from Ensembl `lookup/symbol?expand=1` (deliberately *not* `xrefs/symbol`, which stalls for `homo_sapiens` — see the S0 report). This is what falsified the "~58–60 exons, hundreds of kb" claim. |
| `s0_report.py` | **S0** — renders `results/s0_baseline/report.md` purely from the committed tables (D13). Nothing in that report is hand-written. The claim-audit reference count comes from `lit_claims.tsv`, not from the size of `references.tsv`, because that table also backs the review. |
| `s0_review_build.py` | **Assembles `docs/ip3r_review_2026.md` from `docs/review/*.md`** and renders the bibliography from `references.tsv`; `--pdf` typesets it via pandoc + xelatex, `--check` validates without writing. A cited key with no reference row is a build error, so the text and the bibliography cannot drift. Also rewrites `<sub>`/`<sup>` into pandoc syntax for LaTeX and drops the duplicate H1. |
| `build_findings_page.py` + `findings_page.css` | Renders `docs/findings_summary.md` to one self-contained HTML page, inlining every linked figure as a downscaled WebP data URI, with a paralog summary card read live from the committed tables (and omitted entirely until they exist, so it cannot show a card of zeroes). |

Each ledger task adds its own `s<n>_*.py` here. The PIEZO project's
equivalents (`../piezo_genes/scripts/`, mapped in its `INTERFACE.md`) are
reference implementations worth reading first.

## `src/utils/family.py` — the family definition

**The single place the project's subject is defined.** Pfam signatures
(`FAMILY_PFAM_IDS`), the census signatures, the size band
(`MIN_LENGTH_AA`/`MAX_LENGTH_AA`, set to exclude the ~5,000 aa ryanodine
receptors), the reference panel (human ITPR1/2/3), the sister-family panel
(RYR1/2/3), the literature keyword, and the known-name substrings. Pointing
this app at another family is an edit of this file plus the presets.

## `src/core/` — data models, search orchestration, caching

| File | Key types / functions |
|------|-----------------------|
| `models.py` | `ProteinVariant`, `GeneRecord`, `SearchQuery`, `SearchResult`, `SearchStatus`. |
| `cache.py` | `DiskCache(root, ttl_s)` — JSON-file cache keyed by `(source, query)`. |
| `search.py` | `SearchOrchestrator(cache, email, …)` — fans queries across enabled DBs on worker threads, marshals results back via `queue.Queue`. Special-cases Foldseek (waits for sequence clients, harvests bait sequences). |

## `src/databases/` — one client per source, common interface

| File | Class | Notes |
|------|-------|-------|
| `base.py` | `DatabaseClient` | Abstract; subclasses implement `search()`. |
| `ncbi.py` | `NCBIClient` | Biopython `Bio.Entrez` against the `protein` index. |
| `ensembl.py` | `EnsemblClient` | `rest.ensembl.org` `/lookup/symbol/{species}/{symbol}` → `/lookup/id?expand=1` → transcripts → translations. **S0 changed the symbol-resolution path**: `/xrefs/symbol/homo_sapiens/…` stalls indefinitely (per-species server fault, `BRCA2` too), while `lookup/symbol` works and `xrefs` still works for other species — so `xrefs` is now only a fallback. Evidence: `results/s0_baseline/ensembl_endpoint_probe.tsv`. |
| `uniprot.py` | `UniProtClient` | `rest.uniprot.org/uniprotkb/search`; surfaces canonical + named isoforms. |
| `alphafold.py` | `AlphaFoldClient` | `alphafold.ebi.ac.uk/api/prediction/{acc}` — pivots off UniProt; reports pLDDT and model URLs. |
| `foldseek.py` | `FoldseekClient` | `search.foldseek.com/api/` — async structure-based remote-homology search; opt-in. |
| `compara.py` | `ComparaClient` | **Phylogenomic paralog mining** via Ensembl Compara gene trees — every within-species paralog of the query genes, regardless of naming. Flags adjacent family loci as possible split gene models. Note: ITPRs and RYRs share a gene tree, so its paralog lists include RyRs by design. |
| `blast.py` | `BlastClient` | **Sequence-bait BLAST** via NCBI's public queue — blastp of the longest fetched sequence per gene against nr; finds homologs however they are named. Async, 2–10 min. |
| `interpro.py` | (functions) | InterPro/Pfam helpers — `list_proteins_with_pfam()` (census mode: `max_results=None` paginates to the API's own count, with retries, `strict` failure, raw-page `dump_dir`, `stats` out-dict), `batch_fetch_family_signatures()`, `has_family_signature()`, `fetch_uniprot_sequence()`. Not in `AVAILABLE_CLIENTS`; used by discovery + domain scan. |

`AVAILABLE_CLIENTS: dict[str, type[DatabaseClient]]` in
`databases/__init__.py` is what the orchestrator looks up by name;
`SEQUENCE_SOURCES` / `STRUCTURE_SOURCES` group them for the UI, and
`BAIT_SOURCES` (`BLAST`, `Foldseek`) marks clients the orchestrator runs
*after* the sequence DBs, feeding them the longest fetched sequence per gene
as bait.

## `src/gui/` — tkinter/ttk widgets

| File | Class | Responsibility |
|------|-------|----------------|
| `app.py` | `MainWindow` | Top-level window, menus, status bar, queue polling, export, bundle save (+ report). Exposes `run(project_root, email)` and `load_variants()`. |
| `search_panel.py` | `SearchPanel` | Gene + species + DB checkboxes + preset loader; routes `"mode": "domain_scan"` presets to `on_special_preset`. |
| `results_view.py` | `ResultsView` | Sortable / filterable `ttk.Treeview` of variants; `add_context_action()` adds right-click items. |
| `details_view.py` | `DetailsView` | Read-only labels + clickable URL + scrollable sequence preview. |
| `tools.py` | `ToolsController` | The **Analysis menu**: sequence analysis, novel-paralog discovery, domain-bait scan, deep-dive investigation, selection test (dN/dS), fold check (ESMFold), presence/absence matrix, and `augment_bundle()`. |
| `dialogs.py` | `AnalysisDialog`, `DiscoveryDialog`, `InvestigateDialog`, `DomainScanDialog`, `SelectionTestDialog`, `FoldCheckDialog` | Modal option dialogs; read `.result` after construction. |
| `task_runner.py` | `TaskRunner` | One-at-a-time worker-thread bridge, results marshalled back via queue + `after()`. |
| `text_window.py` | `TextWindow` | Reusable scrollable monospace output window with Save-as. |
| `phylo_view.py` | `PhyloWindow` | **Alignment & tree viewer**: canvas phylogram (novel tips ★-marked) over a conservation-shaded alignment overview. |

## `src/analysis/` — alignment, tree, clustering, mutations, evolution

| File | Key types / functions |
|------|-----------------------|
| `alignment.py` | `pairwise_align()`, `progressive_msa()` (star / MAFFT-optional). |
| `distance.py` | `identity_matrix(covered_only=False)` → `DistanceTable`; `covered_only=True` scores identity over mutually covered columns only (fragment-aware). |
| `tree.py` | `build_nj_tree()`, `tree_to_newick()`, `tree_to_ascii()`. |
| `clusters.py` | `cluster_variants()` → `ClusterReport` with novel-candidate flags. |
| `mutations.py` | `call_mutations()` → list of `MutationCall` vs a reference. |
| `evolution.py` | `conservation_per_column()`, `variable_regions()`. |
| `selection.py` | **dN/dS** (Nei–Gojobori 1986, pure Python): `ng86()`, `fetch_cds()`, `candidate_partners()`, `selection_test()`. |
| `presence.py` | **Presence/absence retention matrix** → `PresenceReport` (member × species counts + retention calls; auto-flags lineage-specific retention). |
| `motifs.py` | **Family-signature PSSMs**: `derive_family_signatures()` → `FamilySignatureSet.scan()` — the domain-evidence fallback for candidates absent from UniProt/InterPro. |
| `esmfold.py` | **ESMFold segment fold-check** (`api.esmatlas.com`): mean pLDDT for a ≤400 aa segment — for IP3Rs, the C-terminal pore module. |
| `figures.py` | Publication figures (matplotlib): heatmap, NJ tree, conservation, length histogram, sources bar. |
| `pipeline.py` | `analyse(variants, …)` — the full pipeline; `write_analysis()` persists artefacts. |

## `src/discovery/` — novel-paralog hunter

| File | Key types / functions |
|------|-----------------------|
| `candidates.py` | `DiscoveryConfig`, `Candidate`, `DiscoveryReport`, `discover_novel_paralogs(…)`. Composite 8-criterion scorer + **promotion evidence gate** (D3: score ≥ 40 requires ≥ 1 family-specific component — outlier / domain / fold / split — else capped at 39). Defaults read from `utils/family.py`. |
| `domain_scan.py` | `run_domain_scan(…)` — the fourth-paralog hunt: enumerate every protein with a family Pfam ID via InterPro, filter out known paralogs *and the ryanodine receptors*, fetch sequences, run analysis + discovery. |
| `exhaustive.py` | `run_exhaustive_hunt(…)` — the all-in-one `"mode": "exhaustive"` pipeline: Compara mine + ortholog expansion + InterPro enumeration → merged census (`itpr_like_census.csv` + FASTA) → analyse + discovery + presence matrix + highlighted tree + bundle/report. |

## `src/investigation/` — deep-dive on one candidate accession

| File | Key types / functions |
|------|-----------------------|
| `pipeline.py` | `Investigator(accession, options).run(on_progress)` → `InvestigationResult`; `InvestigationOptions` (panel defaults to human ITPR1/2/3 from `utils/family.py`). |
| `synthesis.py` | `assess_signals()`, `verdict_from_signals()`, `render_case_file()`. |
| `domain_arch.py` | Line 1 — full InterPro domain architecture. |
| `uniprot_detail.py` | Line 2 — detailed UniProt cross-refs + TM regions. |
| `structure.py` | Lines 3–4 — AlphaFold pLDDT + optional Foldseek. |
| `phylo_context.py` | Line 5 — panel MSA + NJ placement vs the known paralogs. |
| `synteny_lite.py` | Line 6 — genomic coordinates / neighbourhood. |
| `literature.py` | Line 7 — PubMed search, keyword from `utils/family.py`. |

## `src/utils/`

| File | Functions |
|------|-----------|
| `family.py` | The family definition (see above). |
| `species.py` | `SpeciesInfo`, `SPECIES_LOOKUP`, `DEFAULT_SPECIES_PANEL`, `resolve_species()`, `ensembl_species_slug()`. Covers the model organisms, the teleost 3R lineages and the non-metazoan margins. |
| `data_root.py` | `get_data_root()` (falls back to `<project>/data` with a warning), **`require_data_root()`** (raises — use before anything bulk), `free_bytes()`. `python -m src.utils.data_root --require` is the session-protocol check. |
| `exporters.py` | `write_fasta()`, `write_csv()`, `write_json()`. |
| `results_writer.py` | `make_bundle_dir()`, `write_bundle()`. |
| `report.py` | `write_report()` — publication-style markdown + standalone HTML with embedded figures. |

## `src/cli.py` — headless pipeline

`run_headless(project_root, preset, …, analyze, discover, interpro_lookup,
save_results)`: load preset → build `SearchQuery` (`"mode": "domain_scan"`
presets route to `run_domain_scan`) → `SearchOrchestrator` fans queries →
optional `analyse()` → optional `discover_novel_paralogs()` → optional
`write_bundle()` + `write_report()`.
`run_investigate(project_root, accession, …)` deep-dives one accession.

## Data Flow

```
SearchPanel.on_search ──▶ MainWindow._on_search ──▶ SearchOrchestrator.start
                                                          │
                                            ┌─────────────┼─────────────┐
                                            ▼             ▼             ▼
                                   NCBI / Ensembl / UniProt / Compara / AlphaFold
                                            │             │             │
                                            └────────┬────┴─────────────┘
                                                     ▼
                                          (sequence DBs done)
                                                     ▼
                                gather longest sequences ──▶ BLAST + Foldseek
                                                     ▼
                                                queue.Queue
                                                     ▼
                                  MainWindow._drain_queue (after())
                                                     ▼
                                              ResultsView.extend
                                                     ▼
                                                DetailsView.show
```

Post-search, the **Analysis menu** (`ToolsController` + `TaskRunner`) runs
`analyse()` / `discover_novel_paralogs()` / `run_domain_scan()` /
`Investigator` on a worker thread and shows results in `TextWindow`s.
