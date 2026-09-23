## Data availability

Every table and figure this paper reports is committed under `results/` in
the project repository, in the directories `msa_v2`, `phylogeny`, `synteny`,
`duplication`, `reconciliation` and `gene_architecture`, with the inference
limits in `methods`. The paper's deposit manifest lists every deposited file
with its SHA-256. Genome assemblies, proteomes and the per-genome sweep output
are not redistributed; they are retrieved from NCBI and UniProt by the
accession lists deposited beside them.

## Code availability

The scripts that produce every table are in `scripts/`, one prefix per
analysis (`s6_`, `s7_`, `s8_`, `s13_`, `s16_`, `s21_`). Each analysis runs its
constructed negative controls before it writes anything.
`python scripts/s26_assemble.py --paper origin` rebuilds this paper from the
committed tables and fails on a number that no longer matches its source.

## Author contributions

George Dickinson directed the project, set its questions and reviewed its
results. The analysis code and the text of this paper were written by Claude
(Anthropic) under his direction; Chapter 16 of the accompanying thesis
describes how the work was carried out.

## Competing interests

The author declares no competing interests.

## Companion papers

{companions}

## References

Citations are rendered from `results/s0_baseline/references.tsv`, numbered in
order of first appearance.

{references}
