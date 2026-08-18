# CLAUDE.md — Project-specific instructions

@INTERFACE.md

## ⚠ Publication project — read this first every session

This project is executing a multi-session publication plan. At the start of
EVERY session:

1. Read `PUBLICATION_ROADMAP.md` and follow its **Session protocol**
   exactly: open the dashboard (`python3 scripts/dashboard.py --open` plus a
   background `--watch`), `git pull`, verify the data root
   (`python -m src.utils.data_root --require` — it exits non-zero if the
   external drive is not attached, and **that is a stop, not a warning**),
   pick the single next task from the ledger (first `in_progress`, else
   first unblocked `pending`), announce it, and work only that task to its
   completion criteria.
2. Detailed task instructions live in `docs/session_briefs.md`.
3. At session end, without being asked: update the roadmap ledger +
   Results (+ Emergent tasks for anything new), append to `SESSION_LOG.md`,
   **append a biological-findings entry to `FINDINGS.md`** (plain-language
   summary of what the finished task changed biologically; mark unconfirmed
   claims *(pending: task)*), refresh `README.md`, then `git add -A`, commit
   (`S<n>: <one-line outcome>`), and push if a remote exists.
4. Bulk data (genomes, proteomes, BLAST DBs, structures) go under
   `get_data_root()` (`src/utils/data_root.py`) — never into the repo.

## Project summary

A GUI app for enumerating every protein isoform / ortholog of a gene across
NCBI, Ensembl, UniProt, AlphaFold DB and (opt-in) Foldseek, plus the
publication pipeline built on top of it. **Subject: the IP3 receptor (ITPR)
family** — the ER's ligand-gated calcium-release channel, three vertebrate
paralogs (ITPR1/2/3), ~2,700 residues, ~58 exons, tetrameric.

The biology baseline is `docs/ip3r_background.md`; every statement there is
tagged `[db]` / `[lit]` / `[open]` by how far it can be trusted.

**The family-specific hazard.** Ryanodine receptors (RYR1/2/3, ~5,000 aa)
carry *every* ITPR-diagnostic Pfam domain (PF08709, PF01365, PF08454,
PF02815). They are inside every search this project runs. Separating ITPR
from RYR is a positive test at every stage — best-profile assignment or a
labelled-bait margin, with the length band as support, never as the call
(roadmap Decisions **D14**).

**Provenance.** The app, the tooling and the protocol are ported from the
PIEZO project (`../piezo_genes`), whose `INTERFACE.md` maps a reference
implementation for most tasks here. Port the method; never port a result.

## Conventions

- Every `.py` < 500 lines (split if it would exceed). Read `INTERFACE.md`
  before opening source files.
- The family is defined in **one** place, `src/utils/family.py` — Pfam IDs,
  size band, reference accessions, literature keyword. Nothing else in
  `src/` hard-codes a gene name.
- Each database client subclasses `src/databases/base.py:DatabaseClient` and
  normalises its source's records into `ProteinVariant`. New clients are
  registered in `src/databases/__init__.py:AVAILABLE_CLIENTS`.
- The GUI never makes network calls on the main thread — everything goes
  through `SearchOrchestrator` → worker threads → `queue.Queue` →
  `MainWindow._drain_queue` (polled with `root.after()`).
- Caching is opportunistic: misses are recorded as cache failures but never
  block the request; cache writes are best-effort.
- Analysis scripts live in `scripts/s<n>_*.py`, one task per prefix, and
  render their `report.md` purely from their committed tables so report and
  data cannot drift (Decisions D13).

## Running

```
python run.py
python run.py --email george.dickinson@gmail.com    # NCBI Entrez politeness
python run.py --headless --preset ip3r --save-results
python scripts/dashboard.py --open
```

## Adding a new database

1. Implement `src/databases/<name>.py:NewClient(DatabaseClient)`.
2. Register it in `src/databases/__init__.py:AVAILABLE_CLIENTS`.
3. Add a checkbox in `src/gui/search_panel.py:_build()` (sequence vs
   structure grouping).
4. Update `INTERFACE.md`.
