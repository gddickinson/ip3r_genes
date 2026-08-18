# manuscript/ — the submission package

**Nothing here yet.** The package is built in S14a; this file records how,
so the build is not reinvented.

## How it builds

```
python scripts/s14_assemble.py            # figures → claims → stitch → deposit
python scripts/s14_assemble.py --list     # the stages
python scripts/s14_pdf.py                 # typeset PDF (pandoc + xelatex)
```

`s14_assemble.py` exits non-zero if any figure, section or load-bearing
number is missing — which it will, until the analyses have run. That is the
intended behaviour, not a bug to work around.

## Rules

- **Never edit `manuscript.md`.** It is stitched from the numbered section
  files (`00_frontmatter.md` … `15_references.md`, listed in
  `scripts/s14_lib.py:SECTION_ORDER`) and overwritten on every build.
- **Figures are copied, never re-plotted.** `s14_figures.py` copies each
  committed figure from its results directory under its publication number,
  so a manuscript figure cannot differ from the one the analysis produced.
  Files not in the manifest are deleted on every build, so a renumbering
  cannot leave a stale figure behind.
- **Every load-bearing number needs a claim row** in `scripts/s14_claims.py`
  naming its source table and the operation that recovers it (roadmap
  Decisions D12). A re-run that changes a table then fails loudly instead of
  leaving the text stale.
- **Look at the figure** before writing its legend (D11).
- If the framing changes substantially, freeze the old draft in
  `manuscript_v1/` rather than overwriting it, and record why in that
  directory's README.

## Planned figure set

Main figures 1–7, Extended Data 1–11 and Supplementary 1–6 are mapped in
`scripts/s14_lib.py` to the results paths that will produce them. The
working title is *The IP3 receptor family across the eukaryotes*; the real
title is chosen in S14a from what was actually found.
