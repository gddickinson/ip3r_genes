## Data availability

Every table and figure this paper reports is committed in the project
repository under `results/` (the benchmark, the census editions v2, v3, v5 and
v6, the proteome sweeps, the genome sweep and its bait panel, the structural
panel and the iterative-search audit), and `papers/range/deposit_manifest.tsv`
lists each deposited file with its SHA-256. Bulk inputs (reference proteomes,
genome assemblies, structure files and archived API responses) are not
redistributed; each is addressed by a committed accession manifest, and
`papers/range/deposit_notes.md` gives the command that regenerates it from its
public source.

## Code availability

The analysis code is in the repository's `scripts/` directory, one script
prefix per analysis stage, and `python scripts/s26_assemble.py --paper range`
rebuilds this paper from the committed tables, re-verifying every
load-bearing number against the table it comes from.

## Author contributions

George Dickinson conceived and directed the project and reviewed the work.
The analysis code and the text were written by Claude (Anthropic) under his
direction; Chapter 16 of the accompanying thesis describes how the work was
carried out.

## Competing interests

The author declares no competing interests.

## References

Citations are rendered from `results/s0_baseline/references.tsv`, numbered in
order of first appearance.

{references}
