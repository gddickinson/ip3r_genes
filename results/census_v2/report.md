# S2 — Uncapped InterPro enumeration (census v2)

*Rendered from the committed tables by `scripts/s2_report.py` (D13). Do not hand-edit.*

## 1. What was enumerated

Three Pfam signatures define the search space: **PF08709** (Ins145_P3_rec), **PF01365** (RYDR_ITPR), **PF08454** (RIH_assoc). `PF02815` (MIR) is deliberately not among them — it is carried by the O-mannosyltransferases as well as by both receptor families, so it widens the space without adding evidence — and is kept as an annotation column instead.

| signature | InterPro `count` | UniProt count | pages | unique accessions | unique − `count` | cursor restarts | walked to end |
|---|---|---|---|---|---|---|---|
| PF08709 | 12339 | 12506 | 63 | 12507 | 168 | 0 | 1 |
| PF01365 | 13177 | 13232 | 67 | 13233 | 56 | 0 | 1 |
| PF08454 | 12066 | 11888 | 60 | 11890 | -176 | 0 | 1 |

**InterPro's own `count` is not a completeness criterion, and it is wrong in both directions.** PF08709 advertises 12,339 and serves 12,507 (+168); PF01365 advertises 13,177 and serves 13,233 (+56); PF08454 advertises 12,066 and serves 11,890 (-176). The advertised values were stable across re-queries, so this is not a mid-walk edit — and in every case what the API *serves* matches UniProt's independent count for the same signature to within two records (PF08709 12,507 vs 12,506; PF01365 13,233 vs 13,232; PF08454 11,890 vs 11,888), so it is the advertised number that is wrong, not the walk.

Both directions bite. Stopping at an *under*-stated count drops records silently. Trusting an *over*-stated one declares a complete walk incomplete — which is exactly what happened here: the first version of this task's own completeness test failed PF08454 and exited non-zero on a walk that had run its cursor chain to the end. The recorded test is therefore cursor exhaustion (`next: null`), with both counts kept beside it.

### API availability

InterPro answers the same queries on two hosts: the documented `/interpro/api/` and `/interpro/wwwapi/`, which the InterPro website itself calls. During this task the documented host returned HTTP 500 to 11 of 12 probe requests while `wwwapi` served 12 of 12, with identical payloads and the same `count`. `src/databases/interpro.py` now tries both hosts before it sleeps, which is why the walk completed at all.

### What one signature would have missed

The enumeration is a union of three signatures rather than a query on the family's defining one, and this is why:

| signature | records carrying it | records without it | …of those, called ITPR | taxa |
|---|---|---|---|---|
| PF08709 | 12507 | 2914 | 758 | 385 |
| PF01365 | 13233 | 2188 | 262 | 118 |
| PF08454 | 11890 | 3531 | 728 | 366 |

A census built on **PF08709 alone** — the IP3-binding core, the signature that names the family — would have enumerated 12,507 of the 15,421 proteins here and missed 2,914, including 758 this census calls ITPR across 385 taxa. *Dictyostelium* iplA (Q9NA13), a characterised IP3 receptor and a member of this project's own S1 positive panel, is one of the records it would have missed: it carries PF01365 and PF08454 and neither PF08709 nor PF02815 nor PF00520. It **is** in this census, and it is `unassigned` — 2 of 5 architecture signatures is not enough for a positive call, which is the honest outcome and the reason S3 builds profiles.

Only 10,256 records (66.5 %) carry all three seeds; the rest are found by one or two.

## 2. The search space

15,421 distinct proteins carry at least one seed signature, across 1,488 taxa.

The two databases were reconciled rather than pooled, since each is a separate view of the same Pfam matches: **15,417** both, **4** interpro only. The 4 InterPro-only records are all *Taenia solium* fragments of 145–964 aa with newly-issued accessions — InterPro has matched them and UniProt's own Pfam cross-reference has not caught up, which is release skew rather than a disagreement. Nothing was in UniProt's set and missing from InterPro's.

| lineage | records | called ITPR | taxa |
|---|---|---|---|
| Vertebrata | 10936 | 4920 | 764 |
| Metazoa (non-vertebrate) | 3706 | 1131 | 555 |
| SAR | 528 | 268 | 60 |
| Viridiplantae | 75 | 20 | 26 |
| Eukaryota (other) | 57 | 26 | 12 |
| Discoba | 56 | 36 | 25 |
| Fungi | 46 | 25 | 35 |
| Amoebozoa | 12 | 7 | 9 |
| unclassified | 4 | 0 | 1 |
| Bacteria | 1 | 0 | 1 |

### Where the family is, outside the animals

Records are the wrong unit for this question — one well-sequenced alga contributes a dozen — so the table counts **taxa**, and keeps the phyla that put records in the search space without earning a single ITPR call, because those zeros are the interesting rows.

| lineage | phylum | records | taxa | ITPR calls | …high confidence | taxa with a call |
|---|---|---|---|---|---|---|
| Amoebozoa | Evosea | 11 | 8 | 6 | 0 | 6 |
| Amoebozoa | Discosea | 1 | 1 | 1 | 1 | 1 |
| Bacteria | Bacteroidota | 1 | 1 | 0 | 0 | 0 |
| Discoba | Euglenozoa | 36 | 22 | 21 | 0 | 18 |
| Discoba | Heterolobosea | 20 | 3 | 15 | 7 | 3 |
| Eukaryota (other) | (no phylum) | 40 | 7 | 20 | 11 | 6 |
| Eukaryota (other) | Haptophyta | 14 | 3 | 5 | 1 | 3 |
| Eukaryota (other) | Preaxostyla | 3 | 2 | 1 | 0 | 1 |
| Fungi | Mucoromycota | 18 | 16 | 15 | 4 | 15 |
| Fungi | Chytridiomycota | 9 | 6 | 6 | 4 | 5 |
| Fungi | Basidiobolomycota | 5 | 2 | 3 | 3 | 2 |
| Fungi | Entomophthoromycota | 1 | 1 | 1 | 1 | 1 |
| Fungi | Glomeromycota | 10 | 8 | 0 | 0 | 0 |
| Fungi | Kickxellomycota | 2 | 1 | 0 | 0 | 0 |
| Fungi | Zoopagomycota | 1 | 1 | 0 | 0 | 0 |
| SAR | Ciliophora | 343 | 15 | 185 | 37 | 14 |
| SAR | Oomycota | 93 | 26 | 42 | 18 | 20 |
| SAR | (no phylum) | 91 | 18 | 41 | 16 | 11 |
| SAR | Perkinsozoa | 1 | 1 | 0 | 0 | 0 |
| Viridiplantae | Chlorophyta | 60 | 13 | 20 | 2 | 11 |
| Viridiplantae | Streptophyta | 15 | 13 | 0 | 0 | 0 |
| unclassified | (no phylum) | 4 | 1 | 0 | 0 | 0 |

**The plant and fungal records are not scattered — they are phylogenetically clean.** In Viridiplantae every one of the 20 ITPR calls is in **Chlorophyta** (11 taxa, including *Chlamydomonas reinhardtii* with the complete five-signature architecture), and **Streptophyta — the land-plant lineage — has 0 calls** from 15 records in 13 taxa. In Fungi every call sits in an early-diverging phylum (Mucoromycota 15, Chytridiomycota 6, Basidiobolomycota 3, Entomophthoromycota 1), while Dikarya — Ascomycota and Basidiomycota, the yeasts and moulds — contribute **no records to the search space at all**.

That is a much sharper statement than the one this project started from ("plants and fungi lack the family"), and it is the shape a loss looks like: present in the early-diverging lineage of both kingdoms, absent from the derived one. It is also, still, a statement about **what UniProt holds**. Turning it into a statement about genomes is exactly what S20 and S23 are for, and the emergent row that asked this question is narrowed rather than closed.

## 3. The ITPR / RYR call

Every record is called on **domain architecture**: RYR when it carries any of PF02026 (RyR), PF06459 (RR_TM4-6), PF21119 (RYDR_Jsol), PF00622 (SPRY); ITPR when it carries the complete IP3-receptor architecture (PF08709 + PF01365 + PF08454 + PF02815 + PF00520) and none of them. Length is recorded on every row and enters only as support for a medium-confidence call — the size band (2,000–3,600 aa) was chosen to exclude RyRs, so letting it decide would beg the question.

| call | confidence | records |
|---|---|---|
| RYR | high | 5998 |
| ITPR | high | 5187 |
| unassigned | none | 2181 |
| ITPR | low | 770 |
| ITPR | medium | 476 |
| RYR | low | 450 |
| RYR | medium | 359 |

- **ITPR 6,433** (41.7 %), median length 2671 aa, 1,228 taxa
- **RYR 6,807** (44.1 %), median length 4856 aa
- **unassigned 2,181** (14.1 %), median length 678 aa

**1,220 of those calls (7.9 %) rest on a gene symbol alone** — records whose architecture is too partial to decide, carrying a name that is not. D14 forbids the *scorer* from consulting a candidate's own symbol, and that still holds: here the symbol is admitted as explicitly-labelled `low` evidence at census scale, it is recorded in the `reason` column of every row it touched, and the architecture call is kept in its own column so the audit below can score it without symbols anywhere in the input.

### The rule audit

The architecture rule never sees a gene symbol, so symbols are an independent label to score it against. Every symbol-labelled record in the census is used, not a sample.

| test | n | agree | disagree | no call | accuracy where decided |
|---|---|---|---|---|---|
| architecture call vs symbol ITPR | 3249 | 2479 | 0 | 770 | 1.0000 |
| architecture call vs symbol RYR | 2942 | 2492 | 0 | 450 | 1.0000 |
| no symbol-labelled ITPR carries SPRY (PF00622) | 3249 | 3249 | 0 | 0 | 1.0000 |
| RyR (PF02026) → symbol RYR | 2341 | 2341 | 0 | 0 | 1.0000 |
| RR_TM4-6 (PF06459) → symbol RYR | 2274 | 2274 | 0 | 0 | 1.0000 |
| RYDR_Jsol (PF21119) → symbol RYR | 2408 | 2408 | 0 | 0 | 1.0000 |
| SPRY (PF00622) → symbol RYR | 2335 | 2335 | 0 | 0 | 1.0000 |

The RyR rule leans on one conditional claim — that SPRY, which sits in ~114,000 UniProt proteins and is in no sense RyR-specific, *is* diagnostic among proteins that already carry a seed signature. That is the row to check: it is a claim about this search space, not about SPRY.

## 4. The unassigned pile

2,181 records (14.1 %) satisfy neither positive test. They are overwhelmingly fragments: median 678 aa against 2671 aa for a called ITPR, and only 143 of them (6.6 %) fall inside the size band at all. A partial annotation cannot be called either way by a rule that reads absence as evidence, which is the intended behaviour rather than a shortfall — S3's profile sweep is what resolves them.

| source | seed signatures | ITPR signatures (of 5) | length | records |
|---|---|---|---|---|
| both | PF08709 | 1 | out_of_band | 442 |
| both | PF01365;PF08709 | 3 | out_of_band | 402 |
| both | PF01365 | 1 | out_of_band | 323 |
| both | PF08454 | 1 | out_of_band | 242 |
| both | PF08454 | 2 | out_of_band | 209 |
| both | PF08709 | 2 | out_of_band | 131 |
| both | PF01365;PF08454 | 3 | out_of_band | 79 |
| both | PF01365 | 2 | out_of_band | 79 |

### What the size band caught instead

Length never decides a call, which leaves the band free to catch something else. 7 records carry the **complete** five-signature IP3-receptor architecture inside 2,000 residues — 1,528–1,993 aa, against 2,000 for the shortest real family member. An intact architecture in two-thirds of the length is a truncated gene model, and none of them is flagged `Fragment` by UniProt, because a truncated model submitted as a whole protein is not marked as one. All 7 are unnamed locus tags from 3 species.

| accession | gene | species | length (aa) | protein existence |
|---|---|---|---|---|
| A0A8T2J808 | GDO86_004943 | Hymenochirus boettgeri (Congo dwarf clawed frog) | 1528 | Inferred from homology |
| A0A8T2J034 | GDO86_004943 | Hymenochirus boettgeri (Congo dwarf clawed frog) | 1556 | Inferred from homology |
| A0A8T2J539 | GDO86_004943 | Hymenochirus boettgeri (Congo dwarf clawed frog) | 1563 | Inferred from homology |
| A0A9J6C9M0 | PVAND_008274 | Polypedilum vanderplanki (Sleeping chironomid midge) | 1880 | Inferred from homology |
| A0A9N9RQE9 | CHIRRI_LOCUS3757 | Chironomus riparius | 1950 | Inferred from homology |
| A0A8T2IWK9 | GDO86_007249 | Hymenochirus boettgeri (Congo dwarf clawed frog) | 1986 | Inferred from homology |
| A0A8T2IVU7 | GDO86_007249 | Hymenochirus boettgeri (Congo dwarf clawed frog) | 1993 | Inferred from homology |

The call on these is *correct* — they are ITPRs — and the records are still wrong. That is the case for keeping length as a recorded column on every row rather than dropping it once it stopped being part of the call, and it is a starting list for S18's correction register.

## 5. The call checked against sequence

The architecture call is annotation-derived, so it is checked against something annotation-independent: the labelled-bait identity margin from D14, run over the committed per-phylum core panel. Nothing in that test reads a Pfam list or a gene symbol — each panel member is aligned with the six human references and assigned to whichever family it is closer to.

Panel: 48 sequences + 6 labelled baits; MAFFT rc=0, 21,133 columns, 331 s.

Agreement is scored separately inside and outside D7's no-call band (|margin| < 0.10), because a margin of 0.0005 is not a disagreement — it is the test declining to separate two families at that distance.

| identity metric | agree where the margin decides (abs ≥ 0.10) | agree inside the no-call band |
|---|---|---|
| full alignment | **27/27** | 20/21 |
| covered only | **26/26** | 19/22 |

Every row where the two calls differ, under either metric, with its margin:

| accession | species | phylum | architecture | sequence | margin | metric | decisive |
|---|---|---|---|---|---|---|---|
| A0A0M0J4V2 | Chrysochromulina tobinii | Haptophyta | ITPR | RYR | -0.0007 | covered_only | 0 |
| A0A152A7I8 | Tieghemostelium lacteum (Slime mold) (Dictyostelium lacteum) | Evosea | RYR | ITPR | 0.0005 | covered_only | 0 |
| A0A1Q9CN86 | Symbiodinium microadriaticum (Dinoflagellate) (Zooxanthella microadriatica) |  | RYR | ITPR | 0.005 | covered_only | 0 |
| A0A152A7I8 | Tieghemostelium lacteum (Slime mold) (Dictyostelium lacteum) | Evosea | RYR | ITPR | 0.0261 | full_alignment | 0 |

All of them sit inside the no-call band, which is the expected place for the deepest branches: S1 already recorded that *Dictyostelium* iplA — a characterised IP3 receptor — has a bait margin of +0.065, inside the same band. Profile assignment (S3), not a pairwise margin, is what settles these.

## 6. Delta against census v1

Census v1 is what the application returns when it is pointed at the family by name — the committed search bundles plus S1's control panels, 1,571 accessions. It was never meant as a family-wide harvest, so the delta measures what a domain search adds to a name search rather than one census against another.

- v1 1,571 → **v2 15,421** (×9.8), 705 shared
- 866 v1 accessions are absent from v2, each with a verdict:

| verdict | accessions |
|---|---|
| outside UniProt namespace | 806 |
| isoform of a censused canonical | 29 |
| s1_decoy (correctly absent) | 25 |
| no seed signature annotated | 6 |

**No enumeration holes.** Every absence is structural: an accession outside UniProt's namespace was never in the searched space, a `-2` isoform is represented by its canonical parent, S1's decoys are supposed to be absent, and the remainder are short human ITPR fragments (102–408 aa) that carry no seed signature — the real and stated limit of a domain census.

## 7. Representatives

1,251 representatives (ITPR 633, RYR 618): the longest record per species per call, with vertebrates thinned to one per taxonomic order. 48 of them form the per-phylum **core panel** committed as `census_v2_core_panel.faa` — the only FASTA in the repo, since a census-scale one would be re-committed at every later census version.

## 8. What is committed, and where the bulk went

Under the data root: the raw InterPro pages (`raw_api/interpro/<PFAM>/`), both enumerations' parsed intermediates, the UniProt sweep, the seeded-space FASTA (47 MB) and the full representative FASTA. In the repo: the census, the call summary and audit, the seed-contribution and lineage tables, the unassigned profile, the truncated-model list, the sequence check, the delta tables, the representative table, the core panel and 5 figures (`census_growth`, `census_lengths`, `census_lineage`, `census_margin`, `census_space`).

## 9. Caveats

- **This is a UniProt-space census.** RefSeq and Ensembl proteins reachable by name are not in it; that is S3–S5's job, and it is why 806 v1 accessions sit outside it.
- **A domain census cannot see a protein with no domain annotated.** The six unexplained v1 records are exactly that case.
- **A record is not a gene.** One species contributes as many rows as UniProt holds isoforms, redundant TrEMBL entries and alternative models for it, so no copy-number statement can be read off this census — counting paralogs needs the clustering in S6/S7 and the genomic evidence in S5. Where this report counts breadth it counts **taxa**, not records, for that reason.
- The unassigned pile is not noise to be tuned away; it is the fragment population, and it is carried forward with its reasons.
- Every call is architecture-derived and inherits InterPro's annotation. The sequence check in §5 is the independent axis, and it covers the core panel, not the whole census.

