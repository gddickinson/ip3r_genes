# Reviewer checklist — self-audit of the submission package

Built by `python scripts/s14_assemble.py`. Every statement below is either
checked by that build or listed under **Open items** as something it cannot
check.

## What the build guarantees

| Check | How | Status |
|---|---|---|
| Every manuscript figure is the figure its analysis committed | `s14_figures.py` copies, never re-plots; files not in the manifest are deleted on every build | 7 main + 14 Extended Data, 56 panel files, 0 missing |
| Every load-bearing number matches its source table | `s14_claims.py`, 175 declared claims re-read from the committed tables | 175/175 pass → `claims_check.tsv` |
| The text and the bibliography cannot drift | citations are stable keys resolved against `results/s0_baseline/references.tsv` at build time; a cited key with no row is a build error | 29 references, all resolved |
| Every deposited file has a checksum | `s14_deposit.py` | one row per file with size and SHA-256, 237 MB of committed results → `deposit_manifest.tsv` (the file count moves with each build, which also deposits its own outputs) |
| Every excluded bulk class has a regeneration command | `s14_lib.BULK_EXCLUSIONS`; each command was run to produce the data it regenerates | 9 classes → `deposit_notes.md` |
| The PDF contains no dropped glyph | the build converts every super/subscript and symbol the document font lacks and the log is read for `Missing character` | 0 missing characters |
| Every figure is placed at the width it was drawn | `s14_pdf.py` reads each PNG's physical size | 1 figure exceeds the text block by 0.02 in (Extended Data Fig. 14b) |

## What was checked by hand

- **Every one of the 56 figure panels in the package was opened and read
  against its own legend before that legend was written** (roadmap D11). Two
  corrections came out of it: the splice-site count for *Nibea albiflora*
  is 55 spliceable introns of which 54 are GT-AG and one is a minor site, not
  55 GT-AG; and the non-vertebrate copy-number figure is drawn over the 193
  controlled genomes, not all 194, which the legend now states.
- **The claims ledger caught five statements that did not match their
  table** on its first run, all of them in how a row was addressed rather
  than in the number itself (a taxonomic class named `Cyclostomata` where the
  ledger records `Hyperoartia` and `Myxini`; a pair class named without its
  `cross_` prefix; two AlphaFold rows that needed the `ALL` group to be a
  single row; one report phrase quoted loosely).
- **Numbers were taken from the committed tables, not from the roadmap's own
  summaries.** Two disagreed: the alignment is 11,777 columns trimmed to
  1,797, not 11,796 trimmed to 1,790 as an earlier ledger entry recorded, and
  the manuscript uses the tables.

## Known limits, stated in the paper

- The vertebrate scope is one genome per order plus the species the protein
  record left in doubt, not every vertebrate genome.
- 120 of 309 assemblies cannot hold the gene on one contig; every result is
  reported above and below that bar.
- The search's own false-negative rate is 15.2 %, measured against an
  independent replicate rather than estimated.
- Synonymous saturation means no pairwise distance in the paper is
  interpretable; every ω is tree-based.
- AlphaFold covers 0.2 % of the full-length records, so nothing structural
  rests on a predicted model of a whole subunit.

## Open items — these need a human before submission

1. **Affiliation, funding and competing-interests text.** The front matter
   carries a placeholder affiliation and no funding statement.
2. **A figure-by-figure audit against the final page proofs.** Every main and
   Extended Data panel was opened and read against its own legend at source
   resolution (S24; 23 findings, all fixed, recorded in
   `results/supplementary/figure_findings.tsv`), and the mechanical half of
   that audit — every figure has a legend, every legend a figure, panel
   letters match panel files — now runs on every build. Nobody has read them
   at printed size in the assembled PDF.
3. **Reference verification.** The 29 cited references come from a table
   audited in the project's literature baseline, but no one has re-checked
   each DOI and PMID resolves.
4. **Extended Data Figs 10d and 14b are 0.02–0.04 in wider than the text
   block.** Harmless in this build (the PDF scales them by under half a
   percent); worth redrawing before submission.
5. **Journal choice and its formatting.** The package is written to a generic
   Nature-style structure (main figures, Extended Data, Supplementary) and
   has not been fitted to a specific journal's limits.
6. **Deposit and DOI.** Nothing has been uploaded; the manifest is ready.

## How to rebuild

```
python scripts/s14_assemble.py          # figures -> claims -> stitch -> pdf -> deposit
python scripts/s14_assemble.py --list   # the stages
python scripts/s14_claims.py --verbose  # every claim, one line each
```

The build exits non-zero on a missing figure, a missing section, a cited key
with no reference row, or a load-bearing number that no longer matches its
source table.
