# S1 control benchmark — toolchain, recall, specificity, ITPR/RYR separation

_Generated from the committed tables on 2026-08-18 15:14 · benchmark run 92.8 s · 25 positive controls + 31 decoys · promotion threshold score ≥ 40 · size band 2000–3600 aa._

**Headline.** Recall 24/25 (96%); specificity 31/31 (100%), including 6/6 ryanodine-receptor decoys called correctly. The RyR decoys **failed** this benchmark on first run (all six promoted at 45) and pass only because S1 implemented the labelled-bait sister-family test that roadmap D14 specifies. The margin itself is 31/31 correct.

## 1. Toolchain

| Tool | Version | Resolves via | Needed by |
|---|---|---|---|
| `mafft` | v7.526 | PATH | S1/S6 |
| `hmmbuild` | 3.4 | PATH | S3 |
| `hmmsearch` | 3.4 | PATH | S3 |
| `jackhmmer` | 3.4 | PATH | S3 |
| `blastp` | 2.16.0+ | env | S1/S5 |
| `tblastn` | 2.16.0+ | env | S5 |
| `makeblastdb` | 2.16.0+ | env | S5 |
| `miniprot` | 0.18-r281 | PATH | S5/S10 |
| `trimal` | v1.5.rev1 | PATH | S6 |
| `iqtree2` | 2.3.6 | PATH | S7 |
| `datasets` | 18.35.0 | env | S4 |
| `foldseek` | 10.941cd33 | env | S11 (opt) |

All 12 binaries resolve; full paths and Python package versions in `results/toolchain_manifest.txt`. Tools marked `env` live in the reused `piezo1` conda env rather than on the bare PATH (Decisions D18).

## 2. MAFFT is really invoked (step 2)

`analyse(use_mafft=True)` made 1 call(s) to `mafft --auto …` returning 0; 56 sequences aligned to 9,494 columns. `src/analysis/alignment.py:_mafft` falls back to the star alignment silently on a non-zero return code, so the return code is the evidence, not the fact that a call happened. Trace: `mafft_trace.json`. The MSA-derived family signature set has 8 blocks.

## 3. Positive controls — recall

Each hold-out run removes one paralog name from `known_paralogs`; its orthologs must still surface at ≥ 40. The invertebrate / non-metazoan grade (fly `Itpr`, worm `itr-1`, *Dictyostelium* `iplA`) is scored in the baseline run, because no hold-out can protect a member whose name was never recognised in the first place — that is the real use case.

| Group | Run | Members | Recalled ≥ 40 | Recall |
|---|---|---|---|---|
| ITPR1 | holdout_ITPR1 | 8 | 8 | 100% |
| ITPR2 | holdout_ITPR2 | 7 | 7 | 100% |
| ITPR3 | holdout_ITPR3 | 7 | 7 | 100% |
| invert_grade | baseline | 3 | 2 | 67% |
| **overall** | | **25** | **24** | **96%** |

Every vertebrate ortholog scored 50 on the same three components (`size+15,pfam+20,breadth+15`). Note what does *not* fire: the twilight-zone outlier component never contributes to a vertebrate positive, because ITPR1/2/3 are 61–68 % identical to each other — far above the 40 % novelty ceiling. **Recall for this family rests on the domain component.** Every positive took it by the Pfam route (`pfam+20`); the MSA-signature fallback was never exercised, so a census route that can attach neither leaves a true member at 30 — below the threshold — and the fallback's own reliability is untested here. Per-protein table: `positive_controls.tsv`.

### The one miss

* `P29993` **Itpr** (Drosophila melanogaster) scored 35 — fired `size+15,pfam+20`. The taxonomic-breadth component needs one non-known sibling at ≥ 0.35 identity; its nearest is `Uni_itr-1_Caenorha_Q9Y0A1` at **0.342** under the full-alignment metric the scorer uses — short by 0.008. Under the fragment-aware metric (`covered_only=True`, the one S6 uses) the same pair scores 0.420, and breadth would fire: **True**.

This is a one-species artefact as much as a scorer flaw — in the S2 census the invertebrate grade will be represented by many lineages that corroborate each other. But the honest worst case is now measured: **a lone, 40–60 %-diverged true family member with no sibling in the set scores 35 and is missed.** Table: `recall_failures.tsv`. Metric change → Emergent, not silent.

## 4. Negative controls — specificity

**31/31 decoys stay below 40 in every run — specificity 100%.** `score_max` is the worst score across the baseline and all three hold-outs.

| Decoy | Category | aa | in band | family Pfams | Baseline | Max | Pass | Components (at max) |
|---|---|---|---|---|---|---|---|---|
| `P21817` RYR1 | RyR (sister family) | 5038 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `Q92736` RYR2 | RyR (sister family) | 4967 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `Q15413` RYR3 | RyR (sister family) | 4870 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `E9PZQ0` Ryr1 | RyR (sister family) | 5035 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `E9Q401` Ryr2 | RyR (sister family) | 4966 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `A0A8M3ALN1` ryr3 | RyR (sister family) | 4863 | no | 4 | 39 | 39 | ✓ | pfam+20,cluster+10,breadth+15,SISTER→39 |
| `Q9Y490` TLN1 | in-band non-channel | 2541 | yes | 0 | 39 | 39 | ✓ | size+15,cluster+10,breadth+15,GATED→39 |
| `P26039` Tln1 | in-band non-channel | 2541 | yes | 0 | 39 | 39 | ✓ | size+15,cluster+10,breadth+15,GATED→39 |
| `Q9Y6A1` POMT1 | MIR-domain sharer | 747 | no | 1 | 30 | 30 | ✓ | pfam+20,cluster+10 |
| `Q9UKY4` POMT2 | MIR-domain sharer | 750 | no | 1 | 30 | 30 | ✓ | pfam+20,cluster+10 |
| `O00555` CACNA1A | in-band channel | 2506 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q13936` CACNA1C | in-band channel | 2221 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q15878` CACNA1E | in-band channel | 2313 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `P35498` SCN1A | in-band channel | 2009 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q14524` SCN5A | in-band channel | 2016 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q9BX84` TRPM6 | in-band channel | 2022 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q8TDX9` PKD1L1 | in-band channel | 2849 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q7Z442` PKD1L2 | in-band channel | 2459 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `P21333` FLNA | in-band non-channel | 2647 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q01082` SPTBN1 | in-band non-channel | 2364 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q13402` MYO7A | in-band non-channel | 2215 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q5S007` LRRK2 | in-band non-channel | 2527 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q13315` ATM | in-band non-channel | 3056 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q9BZ29` DOCK9 | in-band non-channel | 2069 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `P46939` UTRN | in-band non-channel | 3433 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q96RL7` VPS13A | in-band non-channel | 3174 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q8NFP9` NBEA | in-band non-channel | 2946 | yes | 0 | 25 | 25 | ✓ | size+15,cluster+10 |
| `Q15858` SCN9A | in-band channel | 1988 | no | 0 | 10 | 10 | ✓ | cluster+10 |
| `Q96QT4` TRPM7 | in-band channel | 1865 | no | 0 | 10 | 10 | ✓ | cluster+10 |
| `P98161` PKD1 | giant (out of band) | 4303 | no | 0 | 10 | 10 | ✓ | cluster+10 |
| `Q14204` DYNC1H1 | giant (out of band) | 4646 | no | 0 | 10 | 10 | ✓ | cluster+10 |

## 5. The ITPR/RYR separation (D14) — what this session had to change

All 6 ryanodine-receptor decoys carry all four family-diagnostic Pfam signatures (PF08709, PF02815, PF01365, PF08454) — the `family Pfams` column above. On the first run of this benchmark that handed each of them the +20 domain component, which also satisfies the D3 evidence gate, and every one scored **45** on `pfam+20,cluster+10,breadth+15` and was promoted. Specificity was 25/31 (81%), and every single failure was a RyR.

The roadmap anticipated exactly this and required the rule to be strengthened *before* S2. `src/discovery/candidates.py` now runs a positive sister-family test on every candidate: identity to the nearest labelled sister bait versus identity to the nearest known paralog, with the D7 margin (10%). A candidate closer to the sister family by more than that margin is assigned to it and capped below the threshold (`SISTER→39`), with the margin recorded on the candidate. Three properties matter:

* it is a **positive test**, not a name filter — the candidate's own gene symbol is never consulted, so the unnamed RyR-sized loci S0 found (e.g. zebrafish `LOC101884734`) are called the same way;
* the **length band is not the call** — it contributes only the size component, exactly as D14 requires;
* the margin is **recorded**, so every exclusion is auditable.

### Margin measurements

| Class | n | identity to ITPR bait | identity to RyR bait | margin | calls correct |
|---|---|---|---|---|---|
| true ITPR (all positives) | 25 | 0.14–0.98 | 0.08–0.11 | +0.07…+0.87 | 25/25 |
| RyR decoys | 6 | 0.10–0.11 | 0.64–0.97 | -0.87…-0.54 | 6/6 |
| other decoys | 25 | 0.03–0.11 | 0.02–0.06 | +0.00…+0.05 | n/a (no truth) |

The two classes do not overlap: the narrowest true-ITPR margin is +0.065 and the narrowest RyR margin is -0.536, a gap of 0.601. The margin is 31/31 correct under the full-alignment metric and 31/31 under the fragment-aware one. Table: `bait_margin.tsv`.

**But the deepest branch sits inside the D7 no-call band.** *Dictyostelium* `iplA`, a true family member, has a margin of +0.065 — smaller than the 10% margin D7 requires for a call. It is not misassigned (the margin points the right way, so the test does not exclude it), but the labelled-bait margin **cannot confidently call the non-metazoan grade**. That is a scope limit on this instrument, and it is why D14's other route — best-profile assignment with `itpr.hmm` vs `ryr.hmm` — is not optional for the deep branches S20 will reach.

### Why the full-alignment identity metric understated the risk

RyR is ~1.8× the length of an ITPR, so full-alignment identity dilutes every RyR-vs-ITPR comparison toward zero: the RyR decoys sit at 0.105–0.110 identity to the nearest ITPR bait, **below** the 0.20 novelty floor, so the twilight-zone component never fired. Under the fragment-aware metric the same comparisons rise to 0.249–0.258 — *inside* the 0.20–0.40 twilight zone, worth a further +20. Any future switch to `covered_only=True` (which the recall failure above argues for) would have promoted every RyR to 65 without the sister test. Both columns are in `bait_margin.tsv`.

## 6. What this benchmark does not establish

* **Decoy provenance.** Decoys enter as named UniProt entries. A decoy arriving via Compara or BLAST would gain the +10 homology-provenance component; read every score in the 30–39 band as a 40+ risk on that route.
* **Unnamed sister-family loci.** The sister test needs at least one *labelled* sister bait in the analysis set to measure against. Every search that expects RyR contamination must therefore carry the RyR reference panel (`src/utils/family.py:SISTER_PANEL`) — S2 onward.
* **Fold, split-annotation and MSA-signature-fallback components** never fired here: no Foldseek run, no split gene models in a UniProt-only panel, and every panel member had a real InterPro record so the fallback was never reached. Their contribution is untested.
* **One species per symbol.** Breadth is measured over a 56-sequence panel, not the census; the recall failure above is partly an artefact of that.

## Files

| File | Contents |
|---|---|
| `toolchain.tsv` / `../toolchain_manifest.txt` | tool versions, paths, env |
| `panel_positives.json` / `.fasta` | the 25 positive controls as fetched |
| `panel_decoys.json` / `.fasta` | the 31 decoys as fetched |
| `panel_positives_missing.txt` | panel entries UniProt could not resolve |
| `domain_hits.tsv` | live InterPro family-Pfam hits per accession |
| `mafft_trace.json` | proof MAFFT ran, with return code |
| `positive_controls.tsv` | per-protein recall table |
| `negative_controls.tsv` | per-protein specificity table |
| `recall_failures.tsv` | why each missed positive was missed |
| `bait_margin.tsv` | the D14 margin, both identity metrics |
| `summary.json` | every headline number this report quotes |
