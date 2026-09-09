# Reviewer checklist — self-audit of the submission package

Built by `python scripts/s14_assemble.py`. Every statement below is either
checked by that build or listed under **Open items** as something it cannot
check.

## What the build guarantees

| Check | How | Status |
|---|---|---|
| Every manuscript figure is the figure its analysis committed | `s14_figures.py` copies, never re-plots; files not in the manifest are deleted on every build | 7 main + 16 Extended Data + 6 Supplementary, 140 files, 0 missing |
| Every figure has a legend, every legend a figure, and panel letters match panel files | `s24_audit.py`, run on every S24 build | 23 main and Extended Data figures audited, 0 mismatches |
| Every Extended Data figure is cited, and in order of first mention | `s24_audit.citation_order()`, run on every S24 build against the stitched body | all 16 cited; first mentions run 1 → 16 → `figure_citation_order.tsv` |
| Every load-bearing number matches its source table | `s14_claims.py`, 276 declared claims re-read from the committed tables | 276/276 pass → `claims_check.tsv` |
| The text and the bibliography cannot drift | citations are stable keys resolved against `results/s0_baseline/references.tsv` at build time; a cited key with no row is a build error | 29 references, all resolved |
| Every deposited file has a checksum | `s14_deposit.py` | one row per file with size and SHA-256, 1,936 files and 262.5 MB of committed results → `deposit_manifest.tsv` (the file count moves with each build, which also deposits its own outputs) |
| Every excluded bulk class has a regeneration command | `s14_lib.BULK_EXCLUSIONS`; each command was run to produce the data it regenerates | 9 classes → `deposit_notes.md` |
| The PDF contains no dropped glyph | `s14_pdf.py` runs pandoc `--verbose`, scans the LaTeX log for `Missing character` and exits non-zero on any | 0 missing characters |
| Every figure is placed at the width it was drawn | `s14_pdf.py` reads each PNG's physical size and clamps to the text block | 3 panels are drawn wider than the block and are scaled down (see Open items) |

## What was checked by hand

- **Every one of the figure panels in the package was opened and read against
  its own legend before that legend was written** (roadmap D11), first at
  S14a and again at S24, which read all of them at source resolution and
  recorded **26 findings** in `results/supplementary/figure_findings.tsv` —
  16 legends corrected and 10 figures redrawn. The mechanical half of that
  audit now runs on every build.
- **The claims ledger caught five statements that did not match their
  table** on its first run, all of them in how a row was addressed rather
  than in the number itself (a taxonomic class named `Cyclostomata` where the
  ledger records `Hyperoartia` and `Myxini`; a pair class named without its
  `cross_` prefix; two AlphaFold rows that needed the `ALL` group to be a
  single row; one report phrase quoted loosely).
- **Numbers were taken from the committed tables, not from the roadmap's own
  summaries.** Three disagreed: the alignment is 11,777 columns trimmed to
  1,797, not 11,796 trimmed to 1,790 as an earlier ledger entry recorded; and
  the fall in purifying selection across the ligand pocket's shells is
  measured to each paralogue's *lowest* shell, not its furthest — the
  outermost shell sits slightly above the third in all three, so the pattern
  is a step onto a floor and not a gradient. The manuscript uses the tables.

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
- The ligand-core-against-pore comparison runs on the deep orthologue
  alignments, which are vertebrate by construction, and nothing in it
  separates constraint on ligand binding from constraint on domain packing.

## Open items — these need a human before submission

1. **Affiliation, funding and competing-interests text.** The front matter
   carries a placeholder affiliation and no funding statement.
2. **A figure-by-figure audit against the final page proofs.** Every panel
   has been read at source resolution (S24). Nobody has read them at printed
   size in the assembled 60-page PDF.
3. **Reference verification.** The 29 cited references come from a table
   audited in the project's literature baseline, but no one has re-checked
   each DOI and PMID resolves.
4. **Three panels are drawn wider than the 6.7 in text block** and are
   scaled down when placed: Extended Data Fig. 9d by 0.24 in (a 3.5 %
   reduction, so its type sets slightly smaller than its neighbours'),
   Fig. 14d by 0.04 in and Fig. 3b by 0.02 in. The cause is general rather
   than local: `figstyle.save()` checks the size the figure *declared*, not
   the size the tight bounding box actually wrote, so a figure declared at
   the block width can still be saved wider. Worth fixing in `figstyle` and
   redrawing before submission.
5. **Journal choice and its formatting.** The package is written to a generic
   Nature-style structure (main figures, Extended Data, Supplementary) and
   has not been fitted to a specific journal's limits.
6. **Deposit and DOI.** Nothing has been uploaded; the manifest is ready.

## Version history

`manuscript_v1/` holds the S14a draft frozen on 2026-09-08, before S21 and
S22 were written into the paper. `manuscript_v1/FROZEN.md` says what changed
and how to recover that version's figure set.

## How to rebuild

```
python scripts/s14_assemble.py          # figures -> claims -> stitch -> pdf -> deposit
python scripts/s14_assemble.py --list   # the stages
python scripts/s14_claims.py --verbose  # every claim, one line each
```

The build exits non-zero on a missing figure, a missing section, a cited key
with no reference row, a glyph the document font cannot set, or a
load-bearing number that no longer matches its source table.
