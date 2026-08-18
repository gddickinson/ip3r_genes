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
