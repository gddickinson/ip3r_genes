# Genome manifest notes (S4) — the declared denominator

Every absence claim this project makes is a claim about the assemblies listed in `genome_manifest.tsv`, and about nothing else.

## The scope rule

**One best reference assembly per vertebrate order, UNION every species the protein-level census leaves undecided.**

Order representatives are ranked *annotated > RefSeq > assembly level > scaffold N50* (D9: a RefSeq gene set and a submitter GenBank gene set are not comparable evidence, so the annotation source is ranked on and then carried into the manifest rather than discarded).

The margin species are **derived from the committed census tables**, not hand-listed — the rule and its threshold are in `scripts/s4_manifest_lib.py:derive_margins`, so the denominator can be rebuilt from the census and cannot drift from it:

| Rule | Threshold | Species | What the genome sweep decides |
|---|---|---:|---|
| `zero_hit_proteome` | 0 ITPR records in the reference proteome | 15 | whether the gene is absent or merely unannotated |
| `missing_paralog` | < 3 ITPR records | 101 | whether a paralog is lost or the proteome is shallow |
| `fragment_only` | longest record < 2,000 aa (`family.MIN_LENGTH_AA`) | 100 | whether the gene model is broken or the protein is truly short |
| `anchor` | the calibration species | 5 | nothing — it is the control |

A species can fire several rules; the union is **169 species**, of which 169 resolved to an assembly.

`missing_paralog` and `fragment_only` are *questions, not findings*. S3 measured copy number tracking annotation depth — bird reference proteomes hold a median 13,894 proteins against Mammalia's 34,127 — so a low count is exactly what a genome search exists to resolve. That is why these species are in the denominator rather than in a results table.

## What the manifest contains

- **309 genomes**, from 161 vertebrate orders and 169 margin species (a species that is both is one row with both reasons).
- Rows per reason: `order_rep` 161, `missing_paralog` 101, `fragment_only` 100, `zero_hit_proteome` 15, `anchor` 5
- Orders per class: Actinopteri 68, Amphibia 3, Aves 42, Chondrichthyes 12, Cladistia 1, Coelacanthimorpha 1, Crocodylia 1, Dipnoi 1, Hyperoartia 1, Lepidosauria 2, Mammalia 27, Myxini 1, Testudines 1
- Margin species per class: Aves 130, Actinopteri 21, Mammalia 10, Amphibia 5, Lepidosauria 3
- Annotation source: GenBank 141, RefSeq 168; **33 carry no gene set at all** (D9 — an annotation claim may not be quoted for these).
- Species in the reference dump with no NCBI order rank: 67 (not order-eligible; the margin rules can still pull them in).

## Storage

- Total sequence: **552.5 Gbp**.
- Uncompressed FASTA ≈ **553 GB**; download as zip ≈ **166 GB** (at 30% of sequence length).
- Free on the data root at build time: **1,433 GB** — sufficient for the whole manifest resident at once.

### Assemblies over 10 Gbp

These dominate both the download and S5's per-genome runtime.

| Accession | Species | Gbp | Reasons |
|---|---|---:|---|
| GCF_019279795.1 | Protopterus annectens | 40.1 | order_rep:Ceratodontiformes |
| GCF_964261635.1 | Lissotriton helveticus | 23.2 | order_rep:Caudata |
| GCF_027579735.1 | Bombina bombina | 10.0 | order_rep:Anura |

### Reference species without an NCBI order rank

NCBI taxonomy places these in *incertae sedis* buckets with no ranked order, so they cannot represent one. The margin rules still reach them.

- **Actinopteri** (65): Abudefduf saxatilis, Abudefduf troschelii, Acanthochromis polyacanthus, Ambassis agassizii, Ambassis buruensis, Ambassis kopsii, Ambassis urotaenia, Amphiprion clarkii, Amphiprion ocellaris, Amphiprion percula, Amphistichus argenteus, Argyrosomus japonicus …
- **Chondrichthyes** (2): Pristiophorus japonicus, Squatina squatina

## Reproducing this

```
python3 scripts/s4_build_manifest.py            # offline from cache
python3 scripts/s4_build_manifest.py --refresh-assemblies
```

Raw `datasets` dumps are archived under `<data_root>/raw_api/ncbi_datasets/`, so the manifest re-derives offline from the same bytes it was first built from.
