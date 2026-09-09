## Data availability

Every committed table, figure, alignment, tree and report in this study is in
the deposited archive, listed file by file with its size and SHA-256 in
`manuscript/deposit_manifest.tsv`, so any table can be verified against its
checksum without re-running an analysis. `manuscript/deposit_notes.md` names
each class of bulk data deliberately excluded and gives the command that
regenerates it from its public source.

All primary data are public and unmodified: NCBI Datasets genome assemblies
and their annotations for the 503 genomes named in the two manifests; UniProt
reference proteomes for the 7,691 proteomes swept; InterPro and Pfam
signatures; Ensembl Compara paralogy from a pinned dated archive, whose
registry is committed beside the map; RCSB PDB entries and AlphaFold DB
models for the structure panel; ClinVar missense records; and the SRA runs
listed in the expression tables. No new sequence data were generated.

## Code availability

The analysis is a single repository of scripted stages: one script prefix per
task, each rendering its own report purely from its committed tables, each
running a suite of constructed negative controls before it writes anything.
`python scripts/s14_assemble.py` rebuilds this manuscript package — figures,
claims, stitched text, typeset PDF and deposit manifest — and exits non-zero
on a missing figure, a missing section or a load-bearing number that no longer
matches its source table.

## Author contributions

G.D. designed the study, wrote the analysis code, performed the analyses and
wrote the manuscript.

## Competing interests

The author declares no competing interests.

