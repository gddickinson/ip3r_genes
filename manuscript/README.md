# manuscript/ — the submission package

Built by `python scripts/s14_assemble.py`. **Nothing here is hand-edited
except the numbered section files** and `reviewer_checklist.md`.

## How it builds

```
python scripts/s14_assemble.py            # figures → claims → stitch → pdf → deposit
python scripts/s14_assemble.py --list     # the stages
python scripts/s14_assemble.py --only claims --verbose
```

The build exits non-zero on a missing figure, a missing section, a cited key
with no reference row, a glyph the document font cannot set, or a
load-bearing number that no longer matches its source table. Every guard was
tested by breaking it on purpose.

## What is here

| File | What it is |
|---|---|
| `00_frontmatter.md` … `15_references.md` | the **only** hand-written files; the section order is `scripts/s14_lib.py:SECTION_ORDER` |
| `manuscript.md` | generated: the sections stitched, citations renumbered, bibliography rendered. Never edit it |
| `itpr_family_manuscript.pdf` | generated: the typeset article, 60 pages, every figure placed above its own legend at the width it was drawn |
| `figures/` | generated: every main and Extended Data panel, copied under its publication number in every format the analysis produced |
| `figure_manifest.tsv` | publication number → source path → checksum → drawn size |
| `claims_check.tsv` | every load-bearing number, its source table, the op that recovers it, and whether it still agrees |
| `deposit_manifest.tsv`, `deposit_notes.md` | one row per deposited file with size and SHA-256, plus the bulk classes excluded and the command that regenerates each |
| `reviewer_checklist.md` | the self-audit, and the list of what still needs a human |

## Rules

- **Never edit `manuscript.md`.** It is overwritten on every build.
- **Figures are copied, never re-plotted**, so a manuscript figure cannot
  differ from the one in its own results directory. Files not in the manifest
  are deleted on every build, so a renumbering cannot leave a stale figure.
- **Every load-bearing number needs a claim row** in one of
  `scripts/s14_claims_scope.py`, `s14_claims_history.py` or
  `s14_claims_function.py` (roadmap D12), naming its source table and the
  operation that recovers it.
- **Citations are stable keys** (`[R22]`, `[R22, R25]`) resolved against
  `results/s0_baseline/references.tsv` at build time. A section file never
  carries a citation number, and a cited key with no reference row is a build
  error.
- **Look at the figure** before writing its legend (D11). Every panel in this
  package was opened and read against its legend in S14a and again at source
  resolution in S24, which recorded 26 findings. The mechanical half — every
  figure has a legend, every legend a figure, panel letters match panel files,
  and every Extended Data figure is cited in order of first mention — runs on
  every S24 build.
- If the framing changes substantially, freeze the old draft in
  `manuscript_v1/` rather than overwriting it, and record why there. S14c did
  exactly that; see `manuscript_v1/FROZEN.md`.

## Figure set

7 main figures, 16 Extended Data figures over 57 panel files, and 6
Supplementary figures. **The Extended Data set is numbered in order of first
mention** and `s24_audit.citation_order()` checks that on every build: before
S14c it was not in order, and two of its figures were cited by no sentence at
all.

Supplementary Figures 1–6 are written by **S24**, and `SUPPLEMENTARY_FIGURES`
in `scripts/s14_lib.py` was deliberately empty until then — the build exits
non-zero on a missing figure, and a placeholder for a figure nobody has drawn
would fail for a reason that is not a defect.
