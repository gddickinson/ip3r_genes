# S3 — profile-HMM sweep and the completeness argument

_Rendered from the committed tables on 2026-09-04 07:51 by `scripts/s3_report.py` (D13)._

The census this project can defend rests on two instruments that see different evidence. S2 built the first: a positive **architecture** test over InterPro's annotation of each record (D14b). S3 builds the second: **best-profile assignment** — score every sequence against `itpr.hmm` and `ryr.hmm` and take the winner, but only when it wins by a margin (D7, D14). Where the two agree the call is firm; where only one speaks it says which; where they disagree the record is kept as a conflict rather than resolved by preference.

## 1. The two profiles

| Profile | Seeds | Alignment columns | Match states | Clades |
|---|---:|---:|---:|---|
| `itpr.hmm` | 34 | 6,942 | 2,684 | ITPR1 6, ITPR2 5, ITPR3 5, vertebrate_basal 3, invertebrate 8, non_metazoan 7 |
| `ryr.hmm` | 22 | 8,516 | 4,930 | RYR1 4, RYR2 4, RYR3 4, vertebrate_basal 2, invertebrate 7, non_metazoan 1 |

Aligned with v7.526 (2024/Apr/26) L-INS-i, **single-threaded** (D24). MAFFT's return code is checked and a ragged alignment aborts the build, because S1 established that a silent MAFFT failure degrades to a star alignment with no other symptom — and the thread count is pinned because `--thread -1` is not reproducible: the same RyR seeds aligned twice here gave 8,510 and 8,468 columns and profiles of 4,933 and 4,908 match states. The SHA-256 of every seed set, alignment and profile is recorded in `seed_build_stats.json`, so a rebuild that drifts is visible in the data rather than only in a count.

| Profile | Seed set | Alignment | Profile |
|---|---|---|---|
| `itpr.hmm` | `7c4c1e335e333a5e…` | `246ad551bc789824…` | `fb4b9463dfeda40a…` |
| `ryr.hmm` | `bae02e2d1998b449…` | `84f04c955bb210ee…` | `05386168cd96ebad…` |

Every seed is drawn from the S2 census and carries that census's call, so the profiles are labelled by a rule rather than by a gene name — the same separation S2 kept between its architecture call and the symbols it audited against. The selection rules are stated in `scripts/s3_seed_spec.py` and enforced in code: a seed whose census call, architecture count or length has moved since the manifest was written aborts the build.

5 seeds carry an incomplete architecture and are in the set deliberately, each with its reason:

| Accession | Profile | Species | Arch | Why |
|---|---|---|---|---|
| `A0A1S3H6V5` | itpr | *Lingula anatina* | 4/5 | Lingula — brachiopod; 4/5, lophotroch. breadth |
| `B3RXX7` | itpr | *Trichoplax adhaerens* | 4/5 | Trichoplax — placozoan; 4/5, basal breadth |
| `Q9NA13` | itpr | *Dictyostelium discoideum* | 2/5 | Dictyostelium iplA — characterised receptor, 2/5 signatures; the R3 exception (D21) |
| `A0AAE0LDZ1` | itpr | *Cymbomonas tetramitiformis* | 3/5 | Cymbomonas — Chlorophyta; 3/5, the green algal record the plant question rests on |
| `A0AAV7K8I1` | ryr | *Oopsacas minuta* | 4/5 | Oopsacas — glass sponge; 4/5 shared signatures but 3 RyR-specific ones, so the call is firm |

## 2. Calibrating the instrument before using it

Both profiles were run over S2's archived seeded space — **15,381 records** whose call came from the architecture rule, which reads annotation rather than residues. That makes the agreement between the two a real test and not a restatement. The seed sequences are inside that set by construction, so the number that counts is the one with the 56 seeds removed:

| Scored against | Agree | Disagree | Agreement | Profile declines to call |
|---|---:|---:|---:|---:|
| the architecture call, all records | 11,930 | 1 | 0.9999 | 89 |
| the architecture call, seeds excluded | 11,875 | 1 | 0.9999 | 89 |
| the census call incl. symbol fallback, seeds excluded | 12,612 | 2 | 0.9998 | 570 |

Thresholds: a call needs the winning profile to clear **30 bits**, to span at least **200 match states** (D22), and to beat the loser by more than **10% of its own score**. The margin is relative rather than absolute because bit scores scale with alignable length: a fixed gap would call every full-length protein confidently and no fragment at all. The span floor is measured rather than tuned — in this project's own S0 domain coordinates the shortest observed PF08709, the IP3-binding core that names the family, is 200 aa, and the longest observed SPRY is 137 aa, so the floor sits in the gap between them. What it costs is in the table above: 89 records the architecture rule called lose their profile verdict.

The second instrument earns its place on the records the first could not call: **2,314 of the 3,360 records the architecture rule leaves `unassigned` get a call from the profiles**, because a partial architecture defeats a rule that reads absence as evidence and does not defeat a sequence profile.

S2 filled 1,219 of those from the record's gene symbol at low confidence, which is the one place its call was not purely architectural. The profiles decide those on sequence and overturn **1 of 1,219** — so the symbol fallback was sound, and it is now checked rather than assumed.

**2 records where the two instruments disagree** (`calibration_disagreements.tsv`):

| Accession | Gene | Species | aa | Architecture | Profile | Margin |
|---|---|---|---:|---|---|---:|
| `A0A152A7I8` | DLAC_01018 | *Tieghemostelium lacteum* | 2845 | RYR | ITPR | 55.6% |
| `A0A1Q9CN86` | RYR1 | *Symbiodinium microadriaticum* | 9504 | RYR | ITPR | 55.1% |

## 3. The sweep

* Database: **763 UniProt vertebrate reference proteomes** — 14,414,821 canonical proteins, 6,972,755,129 residues (8.8 GB). One canonical protein per gene: isoform sets are excluded because they add isoforms of genes already counted, and the census counts genes.
* `itpr.hmm`: 6,640 targets at E ≤ 1e-5 (8 min)
* `ryr.hmm`: 17,137 targets at E ≤ 1e-5 (43 min)
* Assignment over the union: **18,501 targets** → ITPR 2,608, RYR 2,522, unassigned 13,371

Before any of that count means anything, **12,609 hits are declined as module-only matches** (D22). Both profiles are full-length channel models, so a protein sharing one small domain with either of them scores against it — and `ryr.hmm` carries SPRY, which D14b had already flagged as sitting in ~114,000 UniProt proteins and being in no sense RyR-specific. Without the gate this sweep called 14,981 vertebrate proteins RYR; the declined pile is headed by *Rspry1*, *Tnnc2*, *Ddx1*, *Cabp1*, *Ccm1*, *Rnf123* — EF-hand and SPRY proteins, not receptors. The declined rows are not duplicated into a table of their own: they are in `hmm_sweep/hmmsearch_assignments.tsv` with `evidence = module`.

Of those, **1,794 carry a gene symbol from this family or its sister** (`subthreshold_family_fragments.tsv`). They are not SPRY proteins that happen to score: they are pieces of split or truncated gene models, sitting under the gate because there is not enough of the protein left to span a family domain. They are kept as an annotation lead rather than discarded.

**618 sweep hits are absent from census v2** (`novel_hits.tsv`) — proteins in a vertebrate reference proteome that the InterPro enumeration never returned. Their calls: ITPR 188, RYR 430.

| Accession | Gene | Species | aa | itpr.hmm | ryr.hmm | Call |
|---|---|---|---:|---:|---:|---|
| `L5LL60` | MDA_GLEAN10006836 | *Myotis davidii* | 942 | 1158.2 | 236.8 | ITPR |
| `A0A4Z2GV61` | ITPR2_2 | *Liparis tanakae* | 917 | 1110.7 | 179.1 | ITPR |
| `A0AAD6AM61` | JOQ06_015010 | *Pogonophryne albipinna* | 747 | 1105.5 | 214.6 | ITPR |
| `A0A9N7YEG8` | PLEPLA_LOCUS10589 | *Pleuronectes platessa* | 669 | 1040.6 | 227.6 | ITPR |
| `A0A8C4RCA7` | — | *Eptatretus burgeri* | 856 | 979.2 | 0.0 | ITPR |
| `A0ABV0UVR8` | ITPR3_2 | *Ilyodon furcidens* | 692 | 975.5 | 214.4 | ITPR |
| `A0A9N7YAW9` | PLEPLA_LOCUS7203 | *Pleuronectes platessa* | 590 | 904.9 | 223.0 | ITPR |
| `A0A6I9PBQ8` | LOC104959425 | *Notothenia coriiceps* | 664 | 900.0 | 199.6 | ITPR |
| `A0A226N2M9` | ASZ78_000997 | *Callipepla squamata* | 567 | 893.1 | 210.8 | ITPR |
| `A0A226NM79` | ASZ78_000425 | *Callipepla squamata* | 545 | 875.6 | 228.1 | ITPR |
| `A0A2P4SWD5` | CIB84_007851 | *Bambusicola thoracicus* | 542 | 874.9 | 227.2 | ITPR |
| `A0A3M0JEH4` | DUI87_24098 | *Hirundo rustica rustica* | 535 | 862.7 | 223.6 | ITPR |

## 4. jackhmmer to convergence, and D10

| Seed | Rounds | New targets per round | Converged | Wall clock | D10 verdict |
|---|---:|---|---|---:|---|
| `itpr1_human` | 10 | 5934 → 681 → 71 → 42 → 452 → 20 → 6 → 12 → 9 → 3 | False | 2.4 h | killed (K3) |
| `itpr_acanthamoeba` | 10 | 5246 → 2562 → 1185 → 8678 → 1075 → 677 → 3300 → 867 → 1312 → 1421 | False | 7.6 h | killed (K3) |
| `itpr_fly` | 10 | 5914 → 695 → 71 → 58 → 470 → 24 → 1 → 2 → 0 → 0 | True | 2.5 h | clean |

D10 requires a coded kill criterion rather than a judgement by eye. For this family the failure mode is specific: an ITPR-seeded run drifts across the shared domain architecture into the ryanodine receptors and reports a converged, confident, wrong answer. Three rules, evaluated per round on that round's own inclusion list (`scripts/s3_kill.py`): **K1** the sister family's share of the model rising more than 10 percentage points above its round-1 value; **K2** the included set growing more than 10× in one round; **K3** hitting the 10-round ceiling without converging. Rounds from the first firing on are reported but excluded from the merge.

**K1 measures the rise, not the level, and that is a result in itself** (D10a). The rule was first coded as a flat 5 % ceiling on sister-family content, and it fired on every run at round 1 — because a single ITPR sequence searched at E ≤ 1e-5 already returns 32%–37% ryanodine receptors before any iteration has happened (itpr1_human, itpr_acanthamoeba, itpr_fly). That is not contamination: the two families are genuine homologues sharing the whole pore, and a search sensitive enough to reach *Acanthamoeba* is necessarily sensitive enough to reach RYR1. The level is a fact about their shared ancestry; only the change in it can be attributed to iterating the model.

* `itpr1_human` (Q14643): killed — reached the 10-round ceiling without converging; no completeness claim may rest on this run; 9 of 10 rounds accepted, 7,226 targets.
* `itpr_acanthamoeba` (L8GF85): killed — reached the 10-round ceiling without converging; no completeness claim may rest on this run; 9 of 10 rounds accepted, 24,902 targets.
* `itpr_fly` (P29993): clean — no kill rule fired; 10 of 10 rounds accepted, 7,235 targets.

### What an iterated model is actually built from

A convergence curve says a search has stopped finding things. It does not say what it found. Every target supporting each run's final model was therefore classified by the sweep's own two-profile verdict and D22 evidence class:

| Composition of the final model | `itpr1_human` | `itpr_acanthamoeba` | `itpr_fly` |
|---|---:|---:|---:|
| unassigned, module | 2,425 (33.6%) | 2,643 (10.1%) | 2,424 (33.5%) |
| ITPR, architecture | 1,656 (22.9%) | 1,656 (6.3%) | 1,656 (22.9%) |
| RYR, architecture | 1,565 (21.7%) | 1,569 (6.0%) | 1,566 (21.7%) |
| ITPR, partial | 952 (13.2%) | 952 (3.6%) | 952 (13.2%) |
| RYR, partial | 617 (8.5%) | 776 (2.9%) | 623 (8.6%) |
| unassigned, no sweep hit | 7 (0.1%) | 18,670 (71.1%) | 7 (0.1%) |

**Module-only matches are 10.1%–33.6% of these models** — the SPRY and EF-hand proteins D22's span gate keeps out of the census. That is what an iterative search at E ≤ 1e-5 accretes if nothing stops it, and it is also why these runs decay to an asymptote of a few new targets a round rather than to zero: the tail being walked is the long tail of proteins sharing one small domain, not the family.

**The family core is the same whichever seed finds it.** The runs started from a human ITPR1, a fly Itpr and an *Acanthamoeba* receptor — a vertebrate, an insect and an amoebozoan — and their final models carry the same family-called content: ITPR/architecture 1,656 in all 3; RYR/architecture 1,565–1,569; ITPR/partial 952 in all 3; RYR/partial 617–776. Starting points that far apart reaching the same core is a completeness statement that does not depend on any one of them converging.

**What differs between them is everything that is not the family.** `itpr_acanthamoeba`'s final model rests on 26,266 targets against `itpr1_human`'s 7,222, and 18,670 of them are proteins neither profile scores at all (7 in `itpr1_human`). **K1 cannot see that drift**: off-family accretion *dilutes* the sister-family share instead of raising it — across this run it falls from 36.6% to 8.9% while the model grows 5.0× — so K1 reads divergence as the opposite. K3 is what caught it, and the composition table above is why the run is reported rather than silently dropped.

## 5. Census v3

**16,039 records** — census v2's 15,421 with a second verdict on every row, plus 618 novel sweep hits.

| Call | Records |
|---|---:|
| ITPR | 8,000 |
| RYR | 7,432 |
| unassigned | 605 |
| conflict | 2 |

| Instruments speaking | Records | Share |
|---|---:|---:|
| both | 12,669 | 79.0 % |
| profile | 2,194 | 13.7 % |
| neither | 605 | 3.8 % |
| architecture | 571 | 3.6 % |

**1,578 records change call from census v2** — every one of them a record the profiles could resolve and the architecture rule could not:

| Transition | Records |
|---|---:|
| unassigned → ITPR | 1,379 |
| unassigned → RYR | 197 |
| RYR → conflict | 2 |

Conflicts (the instruments disagreeing): **2**.

## 6. Completeness, in both directions

A sweep is only a completeness argument if it is checked in the expensive direction as well: how many census-v2 ITPR records from a swept proteome did the profiles *fail* to find? **2,787** — and each one was then looked up in the sweep database itself, one streaming pass, rather than explained away by its gene name:

| Verdict | Records |
|---|---:|
| not in the DB; another entry of the same gene was hit | 1,987 |
| not in the sweep DB — never searched | 800 |

**Not one of them was in the database and missed.** Every apparent miss is a UniProtKB entry that the taxon's *reference proteome* does not contain — a reference proteome is one canonical protein per gene, and the InterPro census enumerated all of UniProtKB — so it was never searched. The profiles' sensitivity failures over 763 vertebrate reference proteomes number **zero** (`missed_by_hmm.tsv`).

**15 of 763 vertebrate reference proteomes carry no ITPR record** (`proteomes_without_hits.tsv`). That is a list of leads for the genome sweep, not a list of losses: a reference proteome is an annotation, and S4/S5 test these against the genomes themselves.

## 7. What this does not settle

* The sweep is **Vertebrata only**. Every count above is a statement about 763 vertebrate reference proteomes, and nothing here bears on the plant, fungal or protist questions S2 opened — those are S20's, using the same code path and manifest format.
* A reference proteome is a *gene set*, not a genome. A gene missing from one may be missing from the annotation, not from the DNA (D9). The zero-hit list is where S5 starts, not what it concludes.
* Profile assignment inherits the seeds' breadth. The deep branches are represented by seven non-metazoan seeds; a lineage more diverged than any of them is outside what this instrument has been shown to call.

