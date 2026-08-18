# roadmap.md — app feature history

The publication plan lives in `PUBLICATION_ROADMAP.md`. This file tracks the
**application's** own features.

## Inherited (v1.0–v1.7, built in the PIEZO project)

| Version | What it added |
|---------|---------------|
| v1.0 | Multi-database search (NCBI, Ensembl, UniProt) with a threaded orchestrator, disk cache, sortable results table and details pane. |
| v1.1 | AlphaFold DB client; results bundles (CSV + JSON + FASTA + metadata) and the publication-style markdown/HTML report. |
| v1.2 | Analysis pipeline: MSA (star / MAFFT), identity matrix, NJ tree, clustering, mutation calling, conservation. |
| v1.3 | Novel-paralog discovery scorer with its composite criteria and evidence gate; InterPro/Pfam domain lookups. |
| v1.4 | Domain-bait scan ("is there a fourth paralog?") driven from a preset. |
| v1.5 | Deep-dive investigation pipeline (domain architecture, UniProt detail, structure, phylogenetic placement, synteny, literature) with a synthesised verdict and case file. |
| v1.6 | Compara paralog mining, BLAST and Foldseek bait clients, presence/absence retention matrix, family-signature PSSMs, ESMFold segment fold-check, dN/dS selection test. |
| v1.7 | Exhaustive hunt mode: one preset that harvests, merges, analyses, scores, and writes the bundle, tree and report. |

## This project's changes

| Date | Change |
|------|--------|
| 2026-08-18 | **Family definition centralised** in `src/utils/family.py`: Pfam signatures, size band, reference panel, sister-family panel, literature keyword and known-name substrings. Nothing else in `src/` hard-codes a gene name, so re-pointing the app at another family is an edit of that file plus the presets. |
| 2026-08-18 | `require_data_root()` + `python -m src.utils.data_root --require`: bulk work fails fast when the external drive is absent instead of silently filling the internal disk. |
| 2026-08-18 | Dashboard generalised: no per-task probes, a documented `results/session_live.json` contract, figure auto-discovery, and sub-task ledger ids (S14a) no longer skipped. |

## Known gaps / candidate work

- The discovery scorer has no explicit sister-family penalty. Today an ITPR
  and a RyR are separated by the size band and the Pfam evidence; S1's
  benchmark will say whether that is enough, and if not the fix belongs in
  `src/discovery/candidates.py`, not in a per-task workaround.
- `ComparaClient` returns RyRs among ITPR paralogs by design (they share a
  gene tree). Consumers must filter; the presets document it.
- No client for a dedicated variant database. S17 fetches ClinVar through
  E-utilities in its own script; if a second task needs it, promote it to
  `src/databases/`.
