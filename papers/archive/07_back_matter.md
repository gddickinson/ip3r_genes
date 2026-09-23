## Data availability

The 297 corrections are offered as a table, `results/annotation_audit/corrections.tsv`,
one row per proposal with the assembly accession, the annotation source,
the contig coordinates and strand, the current state, name and biotype,
the proposal, the priority with the rule that assigned it, the withheld
flag with its reason, and the path of an archived evidence file a curator
can open. The table is deposited with this paper and is intended for
submission to the RefSeq and GenBank annotation teams and to UniProt.

Every committed table and figure behind this paper is listed file by file,
with its size and SHA-256, in `papers/archive/deposit_manifest.tsv`, and
`papers/archive/deposit_notes.md` names each class of bulk data excluded
with the command that regenerates it. All primary data are public and were
used unmodified: NCBI Datasets genome assemblies and their annotations for
the 309 genomes of the scope, UniProt reference proteomes and protein
records, and the Sequence Read Archive runs listed in
`results/expression/runs_selected.tsv`. No new sequence data were
generated.

## Code availability

The analysis is a single repository of scripted stages, each rendering its
report from committed tables and each running constructed negative controls
before it writes anything. `python scripts/s26_assemble.py --paper archive`
rebuilds this paper from those tables and exits non-zero on a missing
figure, a cited reference with no bibliography row, or a load-bearing
number that no longer matches its source table.

## Author contributions

George Dickinson conceived and directed the project and checked its
results. The analysis code and the text of this paper were written by
Claude (Anthropic) under his direction; Chapter 16 of the accompanying
thesis describes how the work was divided and what was corrected.

## Competing interests

The author declares no competing interests.

## Companion papers

{companions}

## References

{references}
