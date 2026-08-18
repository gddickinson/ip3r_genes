# S0 — literature baseline & scope confirmation

*Rendered from the committed tables in `results/s0_baseline/` by `scripts/s0_report.py`. Database snapshot re-derived 2026-08-18.*

## 1. Literature verification

19 atomic claims were extracted from the `[lit]` statements in `docs/ip3r_background.md` and checked against 51 references (45 primary, 6 review).

| verdict | claims | meaning |
|---|---|---|
| `verified` | 12 | primary source(s) support the claim as written |
| `verified_qualified` | 3 | supported, but the wording overstates the evidence; qualified in the review |
| `corrected` | 2 | false or incomplete as written; **struck and replaced** |
| `downgraded_open` | 1 | not established; retagged `[open]` because it is a question this project answers |
| `upgraded_db` | 1 | re-derived from a live database; no longer a literature claim |

### The 4 claims that did not survive as written

**C09 (2. The genes) — `upgraded_db`**  
> Loci: ITPR1 3p26.1, ITPR2 12p11.23, ITPR3 6p21.31

Re-derived from Ensembl (release 15.12) in results/s0_baseline/gene_structure.tsv: all three cytogenetic bands confirmed exactly. No longer a literature claim.

**C11 (2. The genes) — `corrected`**  
> Each is a ~58-60 exon gene spanning hundreds of kb

FALSE as written. Ensembl 15.12 canonical transcripts: ITPR1 62 exons / 354,174 bp; ITPR2 57 exons / 497,888 bp; ITPR3 58 exons / 76,245 bp. The exon range is 57-62, not 58-60, and ITPR3 spans 76 kb - not "hundreds of kb". The 6.5-fold span range across paralogs is itself a finding for S21.

**C14 (3. The sister family) — `downgraded_open`**  
> Both families independently expanded to three vertebrate paralogs

The three-paralog state of each family is a database fact, but the INDEPENDENCE of the two triplications is exactly open question Q2 and has no verified primary source. It is downgraded to [open] and answered by S7/S13/S16. Asserting it as background would make the project's own conclusion an assumption.

**C16 (5. Disease) — `corrected`**  
> ITPR1: specific, often C-terminal/pore-proximal variants cause Gillespie syndrome, which behaves as dominant-negative in the tetramer

Incomplete as written. Gillespie syndrome is caused by BOTH biallelic (recessive) ITPR1 variants and de novo heterozygous variants (Gerber 2016); only the latter act dominant-negatively via the channel domain (McEntagart 2016). Stating only the dominant-negative mechanism omits half the genetics and would misdirect the S17 constraint analysis.

Full per-claim table with references: `lit_claims.tsv`; the bibliography: `references.tsv`; the prose baseline: `docs/ip3r_review_2026.md`.

## 2. Database snapshot re-derived

Every `[db]` number in `docs/ip3r_background.md` was re-queried on 2026-08-18. **All of them reproduced exactly.**

### InterPro protein counts per Pfam signature

| Pfam | name | proteins | note |
|---|---|---|---|
| PF08709 | Ins145_P3_rec | 12,338 | the IP3-binding core (beta-trefoil) |
| PF01365 | RYDR_ITPR (RIH) | 13,177 | RyR and IP3R homology domain |
| PF08454 | RIH_assoc | 12,062 | RIH-associated |
| PF02815 | MIR | 23,453 | also in O-mannosyltransferases - not family-specific |
| PF00520 | Ion_trans | 206,115 | generic voltage-gated-channel pore domain |

### Taxonomic distribution

| Taxon | PF08709 Ins145_P3_rec | PF01365 RIH |
|---|---|---|
| Metazoa | 12,149 | 12,730 |
| SAR (stramenopiles/alveolates/rhizaria) | 309 | 215 |
| Discoba (incl. kinetoplastids) | 48 | 35 |
| Fungi | 41 | 16 |
| Viridiplantae | 40 | 56 |
| Amoebozoa | 7 | 11 |
| Bacteria | 0 | 1 |
| Archaea | 0 | 0 |
| Arabidopsis thaliana | 0 | 0 |
| Saccharomyces cerevisiae | 0 | 0 |
| Paramecium tetraurelia | 21 | 14 |
| Trypanosoma brucei | 1 | 1 |
| Chlamydomonas reinhardtii | 1 | 1 |
| Dictyostelium discoideum | 0 | 1 |

The Q1 anomaly is confirmed and unchanged: **40 Viridiplantae** and **41 Fungi** proteins carry PF08709, the IP3-binding core, while *Arabidopsis thaliana* and *Saccharomyces cerevisiae* have none. S2/S20 identify what they are.

### Reference panel and sister family

| Gene | UniProt | family | length (aa) | status |
|---|---|---|---|---|
| ITPR1 | Q14643 | ITPR | 2758 | reviewed |
| ITPR2 | Q14571 | ITPR | 2701 | reviewed |
| ITPR3 | Q14573 | ITPR | 2671 | reviewed |
| RYR1 | P21817 | RYR | 5038 | reviewed |
| RYR2 | Q92736 | RYR | 4967 | reviewed |
| RYR3 | Q15413 | RYR | 4870 | reviewed |

**The D14 hazard, re-confirmed at the record level.** All 4 of the ITPR-diagnostic Pfam signatures (PF02815, PF08709, PF01365, PF08454) are carried by **all three** human ryanodine receptors as well as all three IP3 receptors. Per-protein detail: `pfam_architecture.tsv`.

## 3. What a single Pfam query actually returns (the D14 hazard, measured)

The background document cites `taxonomy_id:7955 AND xref:pfam-PF08709` as evidence that zebrafish carries four IP3 receptors. Re-running it returns **109 protein records** across **10 gene symbols** — and only 4 of those genes are IP3 receptors:

- IP3R-sized (≤ 3,600 aa): **56 records** across 4 genes — `itpr1a`, `itpr1b`, `itpr2`, `itpr3`
- RyR-sized (> 3,600 aa): **53 records** across 6 genes — `LOC101884734`, `ryr1a`, `ryr1b`, `ryr2a`, `ryr2b`, `ryr3`

**49% of the records returned by the IP3-binding-core Pfam are ryanodine receptors.** This is the single most useful number S0 produced: it is a measured floor on the contamination any Pfam-driven enumeration inherits, in the exact query the planning document quoted as a clean result. It is also why D14 requires a positive ITPR/RYR call rather than a length filter alone — here the length band happens to separate the two cleanly, but that is a fact about zebrafish annotation quality, not a rule.

One of the RyR-sized genes is unnamed (`LOC101884734`, 4,900 aa) — an early example of the unnamed-locus problem S18 audits. Full table: `zebrafish_itpr.tsv`.

## 4. Gene architecture — a `[lit]` claim that failed

| Gene | Ensembl | location | span (bp) | exons (canonical) | band (measured) | band (claimed) |
|---|---|---|---|---|---|---|
| ITPR1 | ENSG00000150995 | chr3:4,493,345-4,847,518 | 354,174 | 62 | 3p26.1 | 3p26.1 |
| ITPR2 | ENSG00000123104 | chr12:26,335,352-26,833,239 | 497,888 | 57 | 12p11.23 | 12p11.23 |
| ITPR3 | ENSG00000096433 | chr6:33,620,331-33,696,575 | 76,245 | 58 | 6p21.31 | 6p21.31 |

All three cytogenetic bands are confirmed, so that claim moves from `[lit]` to `[db]`. The architecture claim does not survive: the canonical exon count is **57–62**, not 58–60, and **ITPR3 spans 76,245 bp** — not "hundreds of kb". Genomic span varies **6.5-fold across the three paralogs** while protein length varies by 3%. That asymmetry is a result in waiting for S21, not a detail.

## 5. App smoke test against live APIs

| Preset | Scope | Client | Source | Status | Hits | s | Note |
|---|---|---|---|---|---|---|---|
| ip3r | 8-species panel | pre-fix | UniProt | ok | 629 | 4.1 | live |
| ip3r | 8-species panel | pre-fix | NCBI | ok | 600 | 8.3 | live |
| ip3r | 8-species panel | pre-fix | Ensembl | failed | 0 | >300 | no response at all; /xrefs/symbol/homo_sapiens stalls |
| ip3r | 8-species panel | pre-fix | TOTAL | partial | 1229 | - | 2/3 sources; results/2026-08-18_091322_itpr1_itpr2_itpr3 |
| ip3r_zebrafish | Danio rerio | pre-fix | UniProt | ok | 87 | 0.0 | cache hit |
| ip3r_zebrafish | Danio rerio | pre-fix | NCBI | ok | 90 | 0.0 | cache hit |
| ip3r_zebrafish | Danio rerio | pre-fix | Ensembl | ok | 33 | 137.5 | danio_rerio is not affected by the stall |
| ip3r_zebrafish | Danio rerio | pre-fix | TOTAL | ok | 210 | - | 3/3 sources; results/2026-08-18_092259_itpr1a_itpr1b_itpr2 |
| ip3r | 8-species panel | post-fix | UniProt | ok | 629 | 4.5 | live, cold cache |
| ip3r | 8-species panel | post-fix | NCBI | ok | 600 | 8.7 | live, cold cache |
| ip3r | 8-species panel | post-fix | Ensembl | failed | 0 | >300 | stall fixed but too slow: 24 sequential gene-species pairs at ~12-95 s each > the then-300 s budget |
| ip3r | 8-species panel | post-fix | TOTAL | partial | 1229 | - | 2/3; results/2026-08-18_123301_itpr1_itpr2_itpr3 - this is what raised run_headless timeout_s 300 -> 900 |
| ip3r | Homo sapiens | post-fix | UniProt | ok | 34 | 2.1 | live, cold cache |
| ip3r | Homo sapiens | post-fix | NCBI | ok | 44 | 4.7 | live, cold cache |
| ip3r | Homo sapiens | post-fix | Ensembl | ok | 41 | 28.5 | human Ensembl data recovered - was 0 before the fix |
| ip3r | Homo sapiens | post-fix | TOTAL | ok | 119 | - | 3/3 sources; results/2026-08-18_123618_s0_smoke_human |

**2 of the 4 runs returned all three sources**, including a human-scoped run after the client fix. Runs 1 and 3 are the instructive failures.

**Ensembl failed the 8-species human preset twice, for two different reasons.** Neither is an outage and neither is a slow network:

| Endpoint | Species | HTTP | s | Verdict |
|---|---|---|---|---|
| /info/ping | - | 200 | 0.65 | up |
| /info/rest | - | 200 | 1.25 | up (release 15.12) |
| /xrefs/symbol/{sp}/ITPR1?object_type=gene | homo_sapiens | timeout | >75 | STALLS - no response, no error |
| /xrefs/symbol/{sp}/ITPR1 | homo_sapiens | timeout | >75 | STALLS - not a query-string artefact |
| /xrefs/symbol/{sp}/BRCA2?object_type=gene | homo_sapiens | timeout | >75 | STALLS - not ITPR-specific |
| /xrefs/symbol/{sp}/itpr1a?object_type=gene | danio_rerio | 200 | 0.61 | OK - not endpoint-wide, so it is per-species |
| /lookup/symbol/{sp}/ITPR1 | homo_sapiens | 200 | 0.61 / 7.61 / 13.91 | OK but HIGHLY VARIABLE - three readings of the same 451-byte call |
| /lookup/symbol/{sp}/ITPR1?expand=1 | homo_sapiens | 200 | 11.70 | OK - full transcript+exon structure; this is the call that costs |
| /lookup/id/ENSG00000150995?expand=1 | homo_sapiens | 200 | 12.94 | OK - same payload by gene id |

**Fault 1 — a per-species stall.** `src/databases/ensembl.py` resolved gene symbols through `/xrefs/symbol/{species}/{symbol}`. That path **stalls indefinitely for `homo_sapiens`** — for `BRCA2` as well as `ITPR1`, so it is neither gene-specific nor a query-string artefact — while the same endpoint answers in 0.6 s for `danio_rerio`, which is why the zebrafish preset got all three sources and the human preset got none. A retry budget could never have fixed this: there is nothing to retry against. `_symbol_to_ids` now resolves through `/lookup/symbol/`, which works, and falls back to `xrefs` only on a clean 404 — never after a transport error, which would walk straight back into the stall. **Human ITPR1 went from 0 variants to 24.**

**Fault 2 — latency, uncovered by fixing fault 1.** Ensembl's speed is unstable: the *same* 451-byte `lookup/symbol` call measured 0.61 s, 7.61 s and 13.91 s within one session, and the `expand=1` call that carries the transcript/exon payload costs ~12 s every time. So one gene in one species costs ~12 s at best and was measured at 95 s at worst. The default panel is 8 species × 3 genes = 24 sequential pairs, i.e. **roughly 5 to 38 minutes** — straddling the old 300 s budget, which is why run 3 still reported 2/3 while run 4 (3 pairs, 28.5 s) sailed through. `run_headless`'s budget is raised 300 s → 900 s; a full panel sweep should use `--species` until `EnsemblClient.search`'s per-species loop is parallelised (Emergent). Do not treat any single latency figure here as a rate — treat the spread as the design constraint.

## 6. Bibliography

51 references, 45 of them primary. Machine-readable: `references.tsv`.

- **R01** Streb H *et al.* (1983) Release of Ca2+ from a nonmitochondrial intracellular store in pancreatic acinar cells by inositol-1,4,5-trisphosphate. *Nature*. PMID [6605482](https://pubmed.ncbi.nlm.nih.gov/6605482/); doi:10.1038/306067a0
- **R02** Furuichi T *et al.* (1989) Primary structure and functional expression of the inositol 1,4,5-trisphosphate-binding protein P400. *Nature*. PMID [2554142](https://pubmed.ncbi.nlm.nih.gov/2554142/); doi:10.1038/342032a0
- **R03** Ferris CD *et al.* (1989) Purified inositol 1,4,5-trisphosphate receptor mediates calcium flux in reconstituted lipid vesicles. *Nature*. PMID [2554143](https://pubmed.ncbi.nlm.nih.gov/2554143/); doi:10.1038/342087a0
- **R04** Mignery GA *et al.* (1989) Putative receptor for inositol 1,4,5-trisphosphate similar to ryanodine receptor. *Nature*. PMID [2554146](https://pubmed.ncbi.nlm.nih.gov/2554146/); doi:10.1038/342192a0
- **R05** Bosanac I *et al.* (2002) Structure of the inositol 1,4,5-trisphosphate receptor binding core in complex with its ligand. *Nature*. PMID [12442173](https://pubmed.ncbi.nlm.nih.gov/12442173/); doi:10.1038/nature01268
- **R06** Bosanac I *et al.* (2005) Crystal structure of the ligand binding suppressor domain of type 1 inositol 1,4,5-trisphosphate receptor. *Mol Cell*. PMID [15664189](https://pubmed.ncbi.nlm.nih.gov/15664189/); doi:10.1016/j.molcel.2004.11.046
- **R07** Berridge MJ *et al.* (1993) Inositol trisphosphate and calcium signalling. *Nature*. PMID [8381210](https://pubmed.ncbi.nlm.nih.gov/8381210/); doi:10.1038/361315a0
- **R08** Yao Y *et al.* (1995) Quantal puffs of intracellular Ca2+ evoked by inositol trisphosphate in Xenopus oocytes. *J Physiol*. PMID [7738847](https://pubmed.ncbi.nlm.nih.gov/7738847/); doi:10.1113/jphysiol.1995.sp020703
- **R09** Rizzuto R *et al.* (1998) Close contacts with the endoplasmic reticulum as determinants of mitochondrial Ca2+ responses. *Science*. PMID [9624056](https://pubmed.ncbi.nlm.nih.gov/9624056/); doi:10.1126/science.280.5370.1763
- **R10** Csordas G *et al.* (2006) Structural and functional features and significance of the physical linkage between ER and mitochondria. *J Cell Biol*. PMID [16982799](https://pubmed.ncbi.nlm.nih.gov/16982799/); doi:10.1083/jcb.200604016
- **R11** Pinton P *et al.* (2008) Calcium and apoptosis: ER-mitochondria Ca2+ transfer in the control of apoptosis. *Oncogene*. PMID [18955969](https://pubmed.ncbi.nlm.nih.gov/18955969/); doi:10.1038/onc.2008.308
- **R12** Bartok A *et al.* (2019) IP3 receptor isoforms differently regulate ER-mitochondrial contacts and local calcium transfer. *Nat Commun*. PMID [31427578](https://pubmed.ncbi.nlm.nih.gov/31427578/); doi:10.1038/s41467-019-11646-3
- **R13** Bezprozvanny I *et al.* (1991) Bell-shaped calcium-response curves of Ins(1,4,5)P3- and calcium-gated channels from endoplasmic reticulum of cerebellum. *Nature*. PMID [1648178](https://pubmed.ncbi.nlm.nih.gov/1648178/); doi:10.1038/351751a0
- **R14** Marchant JS *et al.* (1997) Cooperative activation of IP3 receptors by sequential binding of IP3 and Ca2+ safeguards against spontaneous activity. *Curr Biol*. PMID [9210378](https://pubmed.ncbi.nlm.nih.gov/9210378/); doi:10.1016/s0960-9822(06)00230-5
- **R15** Alzayady KJ *et al.* (2016) Defining the stoichiometry of inositol 1,4,5-trisphosphate binding required to initiate Ca2+ release. *Sci Signal*. PMID [27048566](https://pubmed.ncbi.nlm.nih.gov/27048566/); doi:10.1126/scisignal.aad6281
- **R16** Bezprozvanny I *et al.* (1993) ATP modulates the function of inositol 1,4,5-trisphosphate-gated channels at two sites. *Neuron*. PMID [7686381](https://pubmed.ncbi.nlm.nih.gov/7686381/); doi:10.1016/0896-6273(93)90319-m
- **R17** Supattapone S *et al.* (1988) Cyclic AMP-dependent phosphorylation of a brain inositol trisphosphate receptor decreases its release of calcium. *Proc Natl Acad Sci USA*. PMID [2847175](https://pubmed.ncbi.nlm.nih.gov/2847175/); doi:10.1073/pnas.85.22.8747
- **R18** Chen R *et al.* (2004) Bcl-2 functionally interacts with inositol 1,4,5-trisphosphate receptors to regulate calcium release from the ER in response to inositol 1,4,5-trisphosphate. *J Cell Biol*. PMID [15263017](https://pubmed.ncbi.nlm.nih.gov/15263017/); doi:10.1083/jcb.200402193
- **R19** Ando H *et al.* (2003) IRBIT, a novel inositol 1,4,5-trisphosphate (IP3) receptor-binding protein, is released from the IP3 receptor upon IP3 binding to the receptor. *J Biol Chem*. PMID [12525476](https://pubmed.ncbi.nlm.nih.gov/12525476/); doi:10.1074/jbc.M210119200
- **R20** Higo T *et al.* (2005) Subtype-specific and ER lumenal environment-dependent regulation of inositol 1,4,5-trisphosphate receptor type 1 by ERp44. *Cell*. PMID [15652484](https://pubmed.ncbi.nlm.nih.gov/15652484/); doi:10.1016/j.cell.2004.11.048
- **R21** Maeda N *et al.* (1990) A cerebellar Purkinje cell marker P400 protein is an inositol 1,4,5-trisphosphate (InsP3) receptor protein. Purification and characterization of InsP3 receptor complex. *EMBO J*. PMID [2153079](https://pubmed.ncbi.nlm.nih.gov/2153079/); doi:10.1002/j.1460-2075.1990.tb07386.x
- **R22** Fan G *et al.* (2015) Gating machinery of InsP3R channels revealed by electron cryomicroscopy. *Nature*. PMID [26458101](https://pubmed.ncbi.nlm.nih.gov/26458101/); doi:10.1038/nature15249
- **R23** Fan G *et al.* (2018) Cryo-EM reveals ligand induced allostery underlying InsP3R channel gating. *Cell Res*. PMID [30470765](https://pubmed.ncbi.nlm.nih.gov/30470765/); doi:10.1038/s41422-018-0108-5
- **R24** Paknejad N *et al.* (2018) Structural basis for the regulation of inositol trisphosphate receptors by Ca2+ and IP3. *Nat Struct Mol Biol*. PMID [30013099](https://pubmed.ncbi.nlm.nih.gov/30013099/); doi:10.1038/s41594-018-0089-6
- **R25** Baker MR *et al.* (2021) Cryo-EM structure of type 1 IP3R channel in a lipid bilayer. *Commun Biol*. PMID [34035440](https://pubmed.ncbi.nlm.nih.gov/34035440/); doi:10.1038/s42003-021-02156-4
- **R26** Zalk R *et al.* (2015) Structure of a mammalian ryanodine receptor. *Nature*. PMID [25470061](https://pubmed.ncbi.nlm.nih.gov/25470061/); doi:10.1038/nature13950
- **R27** Yan Z *et al.* (2015) Structure of the rabbit ryanodine receptor RyR1 at near-atomic resolution. *Nature*. PMID [25517095](https://pubmed.ncbi.nlm.nih.gov/25517095/); doi:10.1038/nature14063
- **R28** Efremov RG *et al.* (2015) Architecture and conformational switch mechanism of the ryanodine receptor. *Nature*. PMID [25470059](https://pubmed.ncbi.nlm.nih.gov/25470059/); doi:10.1038/nature13916
- **R29** Nakagawa T *et al.* (1991) The subtypes of the mouse inositol 1,4,5-trisphosphate receptor are expressed in a tissue-specific and developmentally specific manner. *Proc Natl Acad Sci USA*. PMID [1648733](https://pubmed.ncbi.nlm.nih.gov/1648733/); doi:10.1073/pnas.88.14.6244
- **R30** Danoff SK *et al.* (1991) Inositol 1,4,5-trisphosphate receptors: distinct neuronal and nonneuronal forms derived by alternative splicing differ in phosphorylation. *Proc Natl Acad Sci USA*. PMID [1849282](https://pubmed.ncbi.nlm.nih.gov/1849282/); doi:10.1073/pnas.88.7.2951
- **R31** Nucifora FC *et al.* (1995) Molecular cloning of a cDNA for the human inositol 1,4,5-trisphosphate receptor type 1, and the identification of a third alternatively spliced variant. *Brain Res Mol Brain Res*. PMID [7500840](https://pubmed.ncbi.nlm.nih.gov/7500840/); doi:10.1016/0169-328x(95)00194-w
- **R32** Ross CA *et al.* (1992) Three additional inositol 1,4,5-trisphosphate receptors: molecular cloning and differential localization in brain and peripheral tissues. *Proc Natl Acad Sci USA*. PMID [1374893](https://pubmed.ncbi.nlm.nih.gov/1374893/); doi:10.1073/pnas.89.10.4265
- **R33** Taylor CW *et al.* (1999) Expression of inositol trisphosphate receptors. *Cell Calcium*. PMID [10668562](https://pubmed.ncbi.nlm.nih.gov/10668562/); doi:10.1054/ceca.1999.0034
- **R34** Mangla A *et al.* (2020) Type 3 inositol 1,4,5-trisphosphate receptor: A calcium channel for all seasons. *Cell Calcium*. PMID [31790953](https://pubmed.ncbi.nlm.nih.gov/31790953/); doi:10.1016/j.ceca.2019.102132
- **R35** Cai X *et al.* (2012) Ancestral Ca2+ signaling machinery in early animal and fungal evolution. *Mol Biol Evol*. PMID [21680871](https://pubmed.ncbi.nlm.nih.gov/21680871/); doi:10.1093/molbev/msr149
- **R36** Prole DL *et al.* (2011) Identification of intracellular and plasma membrane calcium channel homologues in pathogenic parasites. *PLoS One*. PMID [22022573](https://pubmed.ncbi.nlm.nih.gov/22022573/); doi:10.1371/journal.pone.0026218
- **R37** Mackrill JJ *et al.* (2012) Ryanodine receptor calcium release channels: an evolutionary perspective. *Adv Exp Med Biol*. PMID [22453942](https://pubmed.ncbi.nlm.nih.gov/22453942/); doi:10.1007/978-94-007-2888-2_7
- **R38** van de Leemput J *et al.* (2007) Deletion at ITPR1 underlies ataxia in mice and spinocerebellar ataxia 15 in humans. *PLoS Genet*. PMID [17590087](https://pubmed.ncbi.nlm.nih.gov/17590087/); doi:10.1371/journal.pgen.0030108
- **R39** Hara K *et al.* (2008) Total deletion and a missense mutation of ITPR1 in Japanese SCA15 families. *Neurology*. PMID [18579805](https://pubmed.ncbi.nlm.nih.gov/18579805/); doi:10.1212/01.wnl.0000311912.05593.1e
- **R40** Marelli C *et al.* (2011) SCA15 due to large ITPR1 deletions in a cohort of 333 white families with dominant ataxia. *Arch Neurol*. PMID [21555639](https://pubmed.ncbi.nlm.nih.gov/21555639/); doi:10.1001/archneurol.2011.81
- **R41** Huang L *et al.* (2012) Missense mutations in ITPR1 cause autosomal dominant congenital nonprogressive spinocerebellar ataxia. *Orphanet J Rare Dis*. PMID [22986007](https://pubmed.ncbi.nlm.nih.gov/22986007/); doi:10.1186/1750-1172-7-67
- **R42** Zambonin JL *et al.* (2017) Spinocerebellar ataxia type 29 due to mutations in ITPR1: a case series and review of this emerging congenital ataxia. *Orphanet J Rare Dis*. PMID [28659154](https://pubmed.ncbi.nlm.nih.gov/28659154/); doi:10.1186/s13023-017-0672-7
- **R43** Gerber S *et al.* (2016) Recessive and Dominant De Novo ITPR1 Mutations Cause Gillespie Syndrome. *Am J Hum Genet*. PMID [27108797](https://pubmed.ncbi.nlm.nih.gov/27108797/); doi:10.1016/j.ajhg.2016.03.004
- **R44** McEntagart M *et al.* (2016) A Restricted Repertoire of De Novo Mutations in ITPR1 Cause Gillespie Syndrome with Evidence for Dominant-Negative Effect. *Am J Hum Genet*. PMID [27108798](https://pubmed.ncbi.nlm.nih.gov/27108798/); doi:10.1016/j.ajhg.2016.03.018
- **R45** Klar J *et al.* (2014) Abolished InsP3R2 function inhibits sweat secretion in both humans and mice. *J Clin Invest*. PMID [25329695](https://pubmed.ncbi.nlm.nih.gov/25329695/); doi:10.1172/JCI78173
- **R46** Ronkko J *et al.* (2020) Dominant mutations in ITPR3 cause Charcot-Marie-Tooth disease. *Ann Clin Transl Neurol*. PMID [32949214](https://pubmed.ncbi.nlm.nih.gov/32949214/); doi:10.1002/acn3.51151
- **R47** Beijer D *et al.* (2025) A recurrent missense variant in ITPR3 causes demyelinating Charcot-Marie-Tooth with variable severity. *Brain*. PMID [38938188](https://pubmed.ncbi.nlm.nih.gov/38938188/); doi:10.1093/brain/awae178
- **R48** Molitor A *et al.* (2024) A pleiotropic recurrent dominant ITPR3 variant causes a complex multisystemic disease. *J Exp Med*. PMID [39270020](https://pubmed.ncbi.nlm.nih.gov/39270020/); doi:10.1084/jem.20232178
- **R49** Hytonen MK *et al.* (2025) IP3 receptor depletion in a spontaneous canine model of Charcot-Marie-Tooth disease 1J with amelogenesis imperfecta. *Dis Model Mech*. PMID [39804930](https://pubmed.ncbi.nlm.nih.gov/39804930/); doi:10.1242/dmm.052078
- **R50** Prole DL *et al.* (2019) Structure and Function of IP3 Receptors. *Cold Spring Harb Perspect Biol*. PMID [30745293](https://pubmed.ncbi.nlm.nih.gov/30745293/); doi:10.1101/cshperspect.a035063
- **R51** Cabello-Murgui J *et al.* (2024) ITPR3-associated neuropathy: Report of a further family with adult onset intermediate Charcot-Marie-Tooth disease. *Eur J Neurol*. PMID [39287469](https://pubmed.ncbi.nlm.nih.gov/39287469/); doi:10.1111/ene.16466
