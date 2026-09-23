## Data availability

Every table this paper's numbers and figures are drawn from is committed
under `results/ligand_site/` in the project repository, with the SHA-256 of
each recorded in `ligand_site_stats.json` and in this paper's deposit
manifest. Bulk inputs (reference proteomes, structure files and per-genome
sweep output) are public and are listed with the command that regenerates
each in `deposit_notes.md`.

## Code availability

The analysis is `scripts/s22_run.py` with its modules and negative controls;
this paper is built by `scripts/s26_assemble.py --paper ligand`, which
re-verifies every load-bearing number against its source table.

## Author contributions

George Dickinson directed the project, set its questions and reviewed its
results. The analysis code and the text were written by Claude (Anthropic)
under his direction; Chapter 16 of the project thesis describes how.

## Competing interests

The author declares no competing interests.

## Companion papers

{companions}

## References

Citations are rendered from `results/s0_baseline/references.tsv`, numbered in
order of first appearance.

{references}
