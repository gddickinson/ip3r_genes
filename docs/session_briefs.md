# Session briefs — detailed instructions per roadmap task

Companion to `PUBLICATION_ROADMAP.md` (protocol + ledger live there). Each
brief: goal, steps, completion criteria, outputs. Commands are starting
points — verify versions and flags at run time and record what was actually
executed in the Results column and `SESSION_LOG.md`.

**Reference implementations.** The PIEZO project (`../piezo_genes`) solved
most of these problems once already; its `INTERFACE.md` maps every script to
what it does, and reading the equivalent script before writing a new one
will save a session. Port the *method*. Do not port constants, thresholds
tuned to that family, or — under any circumstances — results.

**The family-specific hazard, repeated because it will bite:** ryanodine
receptors carry every ITPR-diagnostic Pfam domain. Every brief that touches
a search has an ITPR/RYR separation step. That step is not optional
(Decisions D14).

---

## S0 — Literature baseline + scope confirmation

**Goal.** Turn `docs/ip3r_background.md` from a planning document into a
verified baseline, and prove the app runs end to end against live APIs
before any heavy machinery is built on it.

**Steps.**
1. Confirm storage with the user: `python -m src.utils.data_root --require`.
   If the drive is not attached, say so and stop before anything bulk. The
   rest of S0 needs no bulk storage, so the session can continue.
2. Verify every `[lit]` claim in `docs/ip3r_background.md` against a primary
   source. Use PubMed/Europe PMC; prefer the primary paper over a review for
   any specific claim (a disease association, a structural feature, a
   splice-variant name). Collect 15–25 references.
3. Write `docs/ip3r_review_2026.md`: the verified baseline, one citation per
   claim, sections matching the background document. Anything that fails
   verification is **struck from `ip3r_background.md`**, not softened.
4. Re-derive the `[db]` numbers in the background document (the Pfam counts
   and taxonomic distribution) so the file's own snapshot is current, and
   record the date.
5. Smoke-test the app: `python run.py --headless --preset ip3r
   --save-results` and `--preset ip3r_zebrafish`. Confirm NCBI, Ensembl and
   UniProt all return, and that the bundle writes. Record the counts.
6. Confirm the scope decisions with the user if they are present: vertebrate
   genome scope (default: best reference assembly per vertebrate order, as
   in the PIEZO project — S4 fixes it), and whether the non-vertebrate range
   question (S20) is in scope for the paper. Default: yes, it is the
   headline.

**Done when.** `docs/ip3r_review_2026.md` exists with a citation on every
claim; `ip3r_background.md` holds no unverified `[lit]` tag; the app has run
against live APIs with the result counts recorded.

---

## S1 — Toolchain + control benchmark

**Goal.** Install the local toolchain; prove the pipeline's recall and
specificity on ground truth, so the scorer is a validated instrument rather
than a heuristic.

**Steps.**
1. `brew install mafft hmmer blast miniprot brewsci/bio/trimal
   brewsci/bio/iqtree2` (fall back to `conda install -c bioconda` for
   anything brew lacks). Record exact versions in
   `results/toolchain_manifest.txt`. The PIEZO project used a `piezo1` conda
   env for BLAST+; make or reuse an env and record which.
2. Wire MAFFT: run an analysis with `use_mafft=True` and prove it is
   actually invoked (the PIEZO project used a tracer — see
   `../piezo_genes/scripts/bench_lib.py:MafftTracer`).
3. Positive controls: hold out ITPR1, ITPR2 and ITPR3 in turn from
   `known_paralogs` and confirm each held-out ortholog scores ≥ 40. Report
   recall.
4. Negative controls — this is where the family differs from PIEZO. The
   decoy panel must include **RYR1/RYR2/RYR3**, which share every diagnostic
   domain, plus ~20 other large channels and non-channels in the same size
   range (PKD1, TRPM family, CACNA1 subunits, SCN subunits, LRRC8, and a few
   large non-channel proteins of 2,500–3,500 aa). Decoys must score < 40.
   Report specificity, and report the RyR decoys separately: if RyRs pass
   the scorer, the ITPR/RYR separation rule (D14) needs strengthening before
   S2 rather than after.
5. Write `results/benchmark_controls/report.md` with per-protein tables and
   any scorer re-weighting decisions (re-weighting → Emergent, not silent).

**Done when.** All tools run from PATH with recorded versions; recall and
specificity numbers exist with per-protein tables; the RyR decoys are called
correctly; benchmark report committed.

---

## S2 — Uncapped InterPro enumeration (census v2)

**Goal.** The complete protein lists for the family's Pfam signatures, with
a positive ITPR/RYR call on every record.

**Steps.**
1. Enumerate `PF08709`, `PF01365` and `PF08454` to exhaustion via
   `list_proteins_with_pfam(..., max_results=None)`. Keep the polite sleep,
   handle cursor expiry, dump raw pages to
   `data_root/raw_api/interpro/`. Expect ~12–13 k per signature and
   substantial overlap. `PF02815` (MIR) is *not* family-specific — it is in
   the presets for candidate scoring, not for enumeration.
2. **Separate ITPR from RYR.** For every record, make a positive call:
   presence of a RyR-specific Pfam (`PF02026`, `PF06459`, `PF21119`,
   `PF00622`), length band, and gene-symbol evidence, with a recorded
   confidence and a `reason` column. Records that cannot be called stay
   `unassigned` and are counted — an unassigned pile that grows past a few
   percent is a finding, not a nuisance.
3. Fetch sequences for what the analysis needs (longest per species per
   clade bucket; cap configurable), not for all 15 k.
4. Dedupe against the app-level census; produce
   `results/census_v2/itpr_census.csv` + FASTA + a delta report.

**Done when.** InterPro pagination reaches the API's own count (recorded)
for every signature; every record carries an ITPR/RYR/unassigned call with a
reason; census v2 + delta committed; raw dumps archived under the data root.

---

## S3 — Profile-HMM sweep + convergence argument

**Goal.** Sensitivity beyond BLAST/PSSM, and the "no new hits" completeness
evidence.

**Steps.**
1. Curate a seed alignment of ~30 diverse true ITPRs spanning vertebrates,
   invertebrates and the non-metazoan hits from S2. `hmmbuild itpr.hmm`.
   **Build `ryr.hmm` the same way from a RyR seed set** — assignment by best
   profile is what makes D14 operational.
2. Download UniProt reference proteomes (vertebrates first; tens of GB → data
   root). `hmmsearch --domtblout` at E ≤ 1e-5 with both profiles; assign each
   hit to the better-scoring profile and record the bitscore margin.
3. `jackhmmer` from several seeds (human ITPR1, an invertebrate Itpr, a
   protist ITPR) against the same set; iterate to convergence; record
   per-iteration new-hit counts — the completeness curve. Kill divergent runs
   by the coded criterion (D10).
4. Merge novel hits into census v3 with `source="HMM"`; delta report.

**Done when.** The final jackhmmer iteration yields zero new hits (or a
documented asymptote); domtblout archived; census v3 + the convergence plot
committed; every hit carries its ITPR-vs-RYR margin.

---

## S4 — Genome scope manifest (the denominator)

**Goal.** Declare exactly which genomes the census claims to cover.

**Steps.**
1. Decide scope with the user if present. Default, following the PIEZO
   project: one best assembly per vertebrate order (~160), plus every
   species where the S2/S3 status is uncertain (the margin species). That
   gave 194 assemblies there; expect a similar number.
2. Build `results/genome_manifest.tsv` via the NCBI `datasets` CLI:
   accession, species, taxid, order, assembly level, size, contig/scaffold
   N50, date, annotation source (`GCF_` vs `GCA_`, per D9).
3. Write `scripts/fetch_genomes.py`: md5-verified, resumable, `.done`
   markers, `--delete-after-search` for fetch → search → delete operation.
4. Estimate total GB and check it against free space on the drive. **Do not
   bulk-download before storage is confirmed.**

**Done when.** Manifest committed with per-genome rows and a size estimate;
fetch script tested end to end on two small genomes.

---

## S5 — Genomic sweep + per-genome ledger

**Goal.** Find ITPRs that exist only as DNA, and produce the
found/lost/gap ledger that makes an absence claim scientific.

**Steps.**
1. Build the bait panel (`scripts/s5_build_baits.py`): one correctly
   labelled full-length bait per paralog per major clade, screened for
   chimeras (D5), **plus one RyR bait per genome as an internal positive
   control** — a genome where the RyR bait finds nothing has an assembly or
   pipeline problem, not a biological one.
2. Per manifest genome: `miniprot --gff` spliced alignment of the panel →
   candidate loci; cluster loci; classify each genome × paralog cell as
   `found-annotated / found-unannotated / fragment / absent (assembly OK) /
   assembly-gap`; whole-genome `tblastn` rescue for zero-locus cells, with
   trace attribution against the full panel.
3. Novel unannotated models → protein predictions → census v4.
4. `results/genome_ledger.tsv` + `_wide.tsv` + per-genome matrix and
   by-class figures.
5. Storage discipline: `--delete-after-search` unless the drive has room for
   everything. The ledger is resumable; this task may split (S5a fish, S5b
   tetrapods…) per protocol.

**Done when.** Every manifest genome has a ledger row for each of
ITPR1/ITPR2/ITPR3 with evidence paths; the RyR control cell is `found` in
essentially every genome (exceptions investigated and recorded); novel
models in census v4; figures committed.

---

## S20 — Non-vertebrate sweep (the family's true range)

**Goal.** The headline range result, and the project's sharpest question:
the family is overwhelmingly metazoan in the databases, yet holds a few tens
of plant and fungal records while *Arabidopsis* and *S. cerevisiae* hold
none. Which of those is a fact about genomes and which about databases?

**Steps.**
1. `itpr.hmm` (E ≤ 1e-5, the S3 protocol) over every non-vertebrate
   eukaryotic reference proteome, partitioned by group: metazoa
   (non-vertebrate), fungi, viridiplantae, protists (SAR / Discoba /
   Amoebozoa separately — the S2 snapshot says they differ), plus an
   archaeal and a genus-stratified bacterial sample for the negative claim.
2. Every hit gets the ITPR-vs-RYR call (D14) and a lineage row
   (`s20_taxa`-style taxid → phylum/kingdom table).
3. **Chase the plant and fungal hits individually.** For each, ask: is it a
   real gene in a real assembly, a contaminant, or a mis-annotation? Record
   a verdict per record with its evidence. This is the finding, whichever
   way it goes.
4. jackhmmer from non-vertebrate seeds to convergence, per group.
5. A relaxed multi-profile panel (E ≤ 10) over fungi/archaea/bacteria, so
   the negative claims are "not found at a stated sensitivity" rather than
   "not found".
6. Census v5 with lineage columns on every row (D8), delta report,
   convergence plot, presence-by-phylum figure.

**Done when.** Presence/absence per proteome and per phylum committed; every
plant and fungal record has an individual verdict; the negative claims state
their sensitivity; census v5 + figure committed.

---

## S23a — The non-vertebrate sweep's instrument  *(completed 2026-09-05)*

**Goal.** Build and *measure* the instrument before spending 100 GB on it.
Delivered: the declared 194-genome manifest (`s23_scope.py` G1-G5), the
gene-span calibration (`s23_span_calibration.py`, read back by
`s23_calibration.py`), the 80-bait panel with its MIR-domain positive control
(`s23_bait_spec.py` B1-B7, `s23_controls.py`), the copy-number classifier
(`s23_classify.py`) and a 14-genome anchor pilot. → `results/s23_scope/report.md`

**What it established, and S23b must not undo.** Four S5 thresholds and one
S5 selection rule do not transfer outside the vertebrates: the contiguity bar
and `-G` (now per group), the bait length band and the architecture-exception
floor (now stratified per band), `MIN_LOCUS_IDENTITY` (**still unfixed — see
below**) and R3's taxonomic spread rule. And **D26**: S5's RyR positive
control does not exist in plants or fungi, so the MIR-domain sharer replaces
it, drawn per control clade and graded strong/weak.

---

## S23b — The full non-vertebrate sweep

**Blocking work first**, in this order, all from S23a's pilot:

1. **Re-measure `MIN_LOCUS_IDENTITY` for this scope, and report the margin.**
   S5b's 0.40 came from a wide empty gap in the vertebrates (confirmed loci
   >= 0.759, junk 0.23-0.34) measured where every genome has a same-class
   bait. Bands here are whole phyla and 13 slots are unbaited. (The
   *Chlamydomonas* case that first looked like a failure of this floor was
   the bait panel — see S23a — so this is a caution, not a known defect.)
   Measure it the way S5b did, against loci whose identity an assembly's own
   annotation confirms, and report how many `no_locus` calls sit just under
   the floor.
2. **Find a second, non-MIR control for Apicomplexa.** Both apicomplexan
   classes carry one PF02815 protein each across 60 swept proteomes, and in
   the pilot *Toxoplasma gondii* came back with neither a receptor nor a
   control — `uncontrolled`, supporting no absence claim. Without a working
   control "Apicomplexa 0/60" cannot go to assembly level. A conserved
   single-copy ortholog panel, or the assembly's own annotated gene set as the
   proof-of-search.
3. **Measure locus span against CDS footprint.** The metazoan `-G` is 650 kb
   (from the 325 kb *Octopus* gene) and a *Drosophila* Itpr spans 22 kb; in the
   pilot its locus chained 20 alignments across 297 kb, 13x the gene. The call
   was right (`frac_cds` 1.0) but two real genes within 650 kb would be
   counted once — and **copy number is this task's deliverable**.

**Goal.** Then take the negative claims to genome level. A proteome absence is
an annotation fact; only an assembly search makes it a biological one.

**Steps.**
1. Manifest: best assembly per invertebrate phylum and per class in the big
   phyla, plus the S20 priority species — every lineage whose proteome was
   zero-hit, and specifically *Arabidopsis*, rice, a moss, a
   chlorophyte, *S. cerevisiae*, *Neurospora*, a chytrid and a
   microsporidian. Size guard so a fragmentary assembly is not scored.
2. Bait panel from the S6/msa representatives, one per phylum where
   possible, chimera-screened (D5).
3. **Copy-number** classification (not the 3-paralog cell model — outside
   vertebrates the question is how many ITPRs a genome has): full /
   fragment / scrap loci, conservative `n_full`, statuses found /
   fragment_only / assembly_gap / no_locus.
4. tblastn rescue for zero-locus genomes; per-genome evidence under the data
   root; ledger + copy-number figure + pooled novel models.

**Done when.** Every manifest genome has a copy-number row with evidence
paths; the plant/fungal absence claims are backed by assembly searches at a
stated sensitivity; figure committed.

---

## S6 — Alignment upgrade

**Goal.** A publication-quality MSA underneath every downstream result.

**Steps.** Select representatives per clade × per kingdom by explicit coded
rules (D8) — never longest-per-species — including the RyR outgroup set and
every novel model from S5/S23. MAFFT L-INS-i (or `--auto` if size forces it;
record the choice); trimAl `-automated1`; report length, %gaps, columns
kept; regenerate the identity matrix (fragment-aware,
`identity_matrix(covered_only=True)`) and conservation off the trimmed MSA.

**Done when.** `results/msa_v2/` holds aln, trimmed, stats.md, identity and
conservation tables plus figures; S7 and S9 read it.

---

## S7 — ML phylogeny

**Goal.** The tree that supports or corrects the three-paralog picture — and
settles which two paralogs are sisters, which the literature does not agree
on.

**Steps.** `iqtree2 -s trimmed.fasta -m MFP -B 1000 -alrt 1000 -T AUTO`;
root on the ryanodine receptors; render the publication figure (clade
colours from `figstyle.GROUP`, support at claim nodes only). Explicitly test
ITPR1/2/3 monophyly, and run an **AU test over the three rooted sister
hypotheses** ((1,2),3) / ((1,3),2) / ((2,3),1) with tree-corrected clade
membership. RBH-verify any naming call the tree contradicts.

**Done when.** Treefile, support values, rooted figure and `report.md`
committed; the sister question has an AU-test answer or a documented
failure to resolve.

---

## S8 — Synteny

**Goal.** Orthology backed by genomic neighbourhood, not sequence alone.

**Steps.** Extract ±10 coding genes around every ITPR locus from the S5
`genes_slim` tables; reject uninformative symbols and locus tags; build
per-paralog Jaccard matrices, pair-class statistics against a within-genome
control, and a flank consensus per paralog. Test any pseudogene or
fragment locus against the consensus.

**Done when.** `results/synteny/` holds flanks, locus sets, Jaccard
matrices, pair-class stats, report and figures.

---

## S9 — ML selection

**Goal.** dN/dS across the family, with a translation-validated codon
alignment underneath it.

**Steps.** Acquire CDS for every vertebrate tip in the representative set
(Ensembl REST, UniProt→transcript xref chains, and the S5 genome models —
miniprot models need the `--aln` base-level route to recover exact CDS);
validate every CDS by translating it and requiring an exact match to the
protein in the MSA. PAL2NAL codon alignment cross-checked against an
independent in-house mapping; trimAl columns applied codon-aware. Then:
per-paralog M0, pairwise runmode −2, whole-tree M0 vs two-ratio, branch-site
model A on each paralog stem (restart from several initial ω — model A can
converge below its own null), M1a/M2a and M7/M8 site models, and HyPhy
RELAX.

**Done when.** `cds_status.tsv` shows a validated CDS for every tip or a
recorded reason; ω tables, LRTs with BH correction, report and figure
committed.

---

## S10 — Annotation-bug molecular validation

**Goal.** Take the two worst annotation failures S5/S18 surface and prove
them at the molecular level, so the audit's claims are not just
disagreements between databases.

**Steps.** Assembly-version audit (do the two annotations even sit on the
same build?); count **disjoint coding model blocks**, not genes; miniprot the
full-length protein onto the genomic region; blastp-tile each fragment back
onto it; classify exons falling in unannotated gaps; then classify every
intron by which feature its flanking exons fall in, and search ±90 nt CDS
probes across informative junctions against TSA/EST/nt restricted to the
species, counting only hits contiguous across the junction.

**Done when.** Each case has a committed evidence table, an exon-track
figure, and a `report.md` rendered purely from the tables (D13).

---

## S11 — Structures

**Goal.** Confirm that what the census calls an ITPR folds like one, and put
a structural frame under the constraint work.

**Steps.**
1. **Resolve the reference structures by query, not from memory** — RCSB
   search for the family, then pick the cryo-EM IP3R and RyR entries with a
   recorded resolution and state. Hard-coded PDB IDs from a planning
   document are how the wrong protein ends up in a figure.
2. AFDB coverage of the family: how many census proteins and how many
   representatives have a model at all. Expect it to be poor for
   ~2,700-residue proteins; that absence is itself a result.
3. Assemble the panel: AFDB models + cryo-EM references + **two negative
   controls** (an unrelated channel of similar size, and a large non-channel)
   + any predicted models generated for the project. Record `model_coverage`
   — AFDB entries often model an older, shorter UniProt release.
4. All-vs-all / vs-reference TM-align, resumable; per-domain pLDDT with the
   domain boundaries transferred by pairwise alignment; optional Foldseek
   AFDB-wide sweep with UniProt taxonomy resolution.

**Done when.** `structure_manifest.tsv`, TM-score table, per-domain pLDDT
and a `report.md` rendered from the tables; the calibration figure shows the
negative controls.

---

## S12 — Expression evidence

**Goal.** Show that a gene the census claims is real is transcribed and
spliced — the standard bar for an unannotated or contested model.

**Steps.** Pick a species panel where **all three paralogs** have a
validated CDS, so ITPR1/ITPR2 act as internal positive controls in the same
libraries. Build references: the paralog CDS + housekeeping + a
**reversed-CDS decoy per paralog** (same length and composition, no
homology) to set the spurious-mapping floor, plus exon-junction offsets.
Select SRA runs by **verified BioSample attributes**, not by keyword alone.
Stream `fastq-dump -X | hisat2`, counting reads, junction-spanning reads
(≥ 8 nt each side), covered bases and a coverage vector. Detection needs all
three: junction reads ≥ 2, reads ≥ 5, and reads above that run's own decoy.

**Done when.** Per-run and per-tissue tables, an atlas cross-check, report
and figure committed.

---

## S13 — Reconciliation & dating

**Goal.** Place the duplications that made ITPR1/2/3 on a dated species
tree.

**Steps.** Hand-curate the species tree as an **input** (D15) with a
literature age, spread and source on every calibrated node, and a `--check`
mode that fails on a missing species, an uncalibrated node or a parent
younger than its child. LCA reconciliation (Goodman/Page; losses per Zmasek
& Eddy) across every gene-tree topology from S7 and both cyclostome
variants. Check every implied loss against the S5 ledger: sampling artefact
vs corroborated loss vs loss-with-remnant.

**Done when.** Event tables, duplication placements, a loss audit, the dated
backbone figure and a report are committed, with the answer stated as a
bracket across topologies rather than a single number.

---

## S15 — Loss dynamics

**Goal.** How many independent losses per paralog, where on the tree, and
are the remnants still recognisable?

**Steps.** Disambiguate the sweep's undecided tblastn traces by synteny,
with **measured** accuracy against known loci. Build the genome × paralog
character matrix with bitscore- and synteny-based states side by side. Run
the ORF-integrity screen with the assembly-error confounder controlled: a
paired within-genome test against that genome's own other paralogs, a
disabling-density threshold calibrated on the controls, and a correlation of
density against contig N50. Then Dollo parsimony as the primary count plus
Mk model fits (ER / ARD / and an irreversible model with gain pinned to
zero), over a sensitivity matrix of coding × evidence × branch lengths ×
contiguity filter. For any dead locus, read its lesions in reference-protein
coordinates and test shared disabling lesions against a Poisson null — with
independently decayed pseudogenes from *different* loss events as the
control, because tblastn keeps whatever is most conserved.

**Done when.** Loss events with per-event support, the sensitivity matrix,
the integrity tables and the fossil tables are committed with a report and
figure.

---

## S16 — Duplication history

**Goal.** Are ITPR1/2/3 2R ohnologs, and are teleost itpr1a/itpr1b from 3R?

**Steps.** Build a per-**locus** copy table (read every alignment out of each
sweep summary, including the secondary loci the ledger's one-row-per-cell
view cannot see; merge overlapping same-paralog models so a split model
cannot manufacture a duplication). For 2R: cross-window paralogy from a
pinned Ensembl BioMart archive release against a permutation null of real
genomic windows (D17), with the ITPR genes removed from their own windows; a
genome-wide block scan ranking blocks by shared gene *families*; a quartet
test; and cross-species replication. For 3R: pre-3R outgroup / teleost /
extra-WGD grouping with Acipenseriformes and Salmoniformes held out, flank
sets for **every** copy, double-conserved synteny against both a tetrapod
and a pre-3R ray-finned reference, and cross-anchor block identity.

**Done when.** The 2R and 3R questions each have a stated answer with its
null, its effect size and its caveat; tables, report and figure committed.

---

## S17 — Constraint & function

**Goal.** The project's mechanistic payoff: does the constrained core
explain the clinical variants, and is the IP3-binding core under different
constraint from the pore?

**Steps.**
1. Build deep within-paralog alignments (one locus per genome per paralog —
   3R ohnologs must not enter a within-gene sample) at strict coverage and
   identity floors, plus curated references.
2. Transfer the domain architecture onto each reference's own numbering by
   pairwise alignment, with a primary per-residue assignment (most specific
   element wins) and named linkers.
3. Four conservation layers on the same residues: deep within-paralog,
   vertebrate, family-wide, and the representative-alignment layer kept as
   the **control** for what the deep alignment bought.
4. Variants: ClinVar plus UniProt features. **Map transcripts, do not assume
   them** — fetch each cited transcript's own translated CDS, align it to the
   canonical, and require the reference amino acid to match. Then the ROC-AUC
   test of whether pathogenic positions are more constrained than benign,
   with a protein-wide control.
5. HyPhy FEL per site carried onto the reference protein through a validated
   coordinate transfer (every mapped site must agree in amino acid).
6. Paint constraint into the B-factor column of the structures (unscored
   residues as −1, never 0).

**Done when.** Per-residue constraint tables in canonical numbering, the
per-element comparison against each protein's own distal regions, the
variant test with its AUC and null, the painted structures, report and
figure — all committed.

---

## S18 — Annotation-quality audit

**Goal.** Quantify how badly the family is recorded, and produce a
correction list worth sending to the databases.

**Steps.** Score every gene-scale locus in the genome scope:
complete / split / fragmentary / noncoding / unannotated (same strand only),
a name verdict and a symbol verdict applied identically to genome gene names
and UniProt gene fields, joined to S15's integrity call (the veto, D6).
Break the results down by annotation source (`GCF_` vs `GCA_`, D9). Run the
protein-side audit over the S3 full-length hits: a sequence-based paralog
call against the labelled bait panel with a recorded margin (D7 — and here
the wrong-paralog risk includes calling an ITPR a RyR), and the Pfam-recall
intersection against S2. Resolve every S3 zero-hit proteome against its own
genome. Then write the correction items: assembly, coordinates, current
state, proposal, evidence file, priority.

**Done when.** The locus audit, the by-source breakdown, the protein audit,
the curated case register and the correction list are committed with a
report and figure.

---

## S19 — Methods results

**Goal.** What each method was worth — the paper's honest methods section,
as results.

**Steps.** Census growth and per-method contribution, restricted to the
accessions actually present in each searched database (build the universe by
one `grep '^>'` pass per proteome FASTA — "the enumeration holds a record the
HMM did not return" must be separable from "that record was never in the
searched database"). The per-**gene** recovery question: how many genes are
unreachable from any protein database. Assembly contiguity as a confounder,
using a never-lost paralog as the control whose every non-`found` cell is a
false negative, and a floor scan that calibrates D4's contiguity bar. Exact
bait-panel ablation by re-clustering retained alignments from bait subsets
(validate the reproduction cell-for-cell against the committed ledger before
trusting any ablation). Iterative-search drift from the raw jackhmmer logs,
including the killed runs. Limits of the inference methods, recomputed from
the committed tables.

**Done when.** Every table above is committed, the ablation's validation
passes cell-for-cell, and any coverage cap the workflow imposes is logged
rather than silent.

---

## S21 — Gene architecture

**Goal.** The exon/intron architecture of a ~58-exon gene family, from data
the sweep already produced.

**Steps.** Parse the CDS blocks out of every retained miniprot GFF (each
carries its genomic span, its span in the bait protein's coordinates and its
phase); merge the pair of records miniprot emits either side of a
frameshift, or a broken locus looks exon-rich. Fix one canonical bait per
paralog as a common coordinate frame and map those onto alignment columns so
boundaries are comparable *between* paralogs. Then: exon count / CDS length /
span / intron sizes per paralog and class, all cross-paralog comparisons
paired within genome (D16); conserved intron positions and the
shared-ancestral-intron test (same column *and* same phase, against a
permutation null drawn from columns where both paralogs have residues);
whether each `split`/`fragmentary` annotation's pieces are bounded by real
splice sites; and the tandem-duplication test (one bait aligning twice over
the same part of the protein at disjoint positions), with a detector control
that proves the test fires when two genes really exist.

**Done when.** Gene models, exon blocks, the architecture comparisons, the
boundary tests and the fragment analysis are committed with a report and
figure.

---

## S22 — Ligand-site evolution

**Goal.** The one part of the receptor the ryanodine receptors do not share
functionally is the IP3-binding core. Ask what evolution did to it.

**Steps.** Using S9's codon alignments and S17's per-residue constraint and
domain map: compare constraint and ω in the IP3-binding core against the
pore module and against the protein as a whole, per paralog, paired within
alignment. Ask whether the binding-core residues that contact IP3 (from the
structures, transferred by alignment) are more constrained than the rest of
the core. Then the comparative question: in lineages where the upstream PLC
/ IP3 pathway is reported reduced or absent, is the binding core relaxed
relative to the pore? Report the taxon list and its power honestly — this
test is easy to underpower, and saying so is a result.

**Done when.** The per-element comparison, the contact-residue test and the
lineage test are committed with their nulls, effect sizes and a power
statement.

---

## S14a — Manuscript assembly

**Goal.** The submission package, built by a script, from the committed
tables.

**Steps.**
1. `scripts/s14_figures.py` copies each committed figure into
   `manuscript/figures/` under its publication number, in every format the
   analysis produced. Nothing is re-plotted, so a manuscript figure cannot
   differ from the one in its results directory. Files not in the manifest
   are deleted on every build.
2. Write the numbered section files listed in `s14_lib.SECTION_ORDER`.
   Never edit the stitched `manuscript.md` — it is overwritten.
3. `scripts/s14_claims.py`: a claim row for **every** load-bearing number
   (D12), naming its source table and the operation that recovers it. The
   file ships with an empty list and worked examples of each op.
4. `scripts/s14_deposit.py` walks the deposited tree → manifest with size
   and SHA-256 per file, plus the bulk classes deliberately excluded **with
   the command that regenerates each**.
5. `scripts/s14_pdf.py` typesets it. Traps already handled in that script's
   comments: pandoc parses YAML metadata as markdown, `\(...\)` needs
   `+tex_math_single_backslash`, and Latin Modern silently drops ω/χ/≥/≈.
6. Write `manuscript/reviewer_checklist.md` — the self-audit, with a section
   listing what still needs a human.

**Done when.** `python scripts/s14_assemble.py` runs the whole chain and
exits non-zero on a missing figure, a missing section or a failed claim; the
PDF builds; the checklist is committed.

---

## S24 — Supplementary figures + figure audit

**Goal.** Show the alignments and structures the main figures rest on, and
**look at every figure** (D11).

**Steps.** Render the representative alignment and the columns the tree
actually saw (trimAl writes no column map — recover it by an exact column
walk); the clinically labelled positions at residue level across the
paralogs, with a guard that each stated residue really is the residue at its
column; the within-paralog alignments the constraint map is computed on; and
the trimmed codon alignment. Then the structure figures: the constraint map
painted on the channels, and every labelled variant with its enrichment
test. Finally, open every main and Extended Data figure and read it against
its own legend. Errors of the kind "the legend says four, the figure draws
three" are invisible to every table check in the project.

**Done when.** The supplementary figures and their stats table are
committed, and the audit's findings are fixed in the manuscript text with
each fix recorded in the session log.

---

## S14b — Deposit + release (human-gated)

Not completable autonomously. Zenodo DOI, flipping the repo public (D2),
verifying every reference, affiliation and funding text, and the preprint
upload all need the author. Prepare everything, list the open items in
`manuscript/reviewer_checklist.md`, and hand over.

---

## S14c — Manuscript rewrite pass

One full pass once every analysis has landed: does the paper still lead on
its strongest result? Are the figure numbers still right? Does every claim
row still pass? Version the previous draft rather than overwriting it — the
PIEZO project froze v1 in `manuscript_v1/` when its framing changed, and
that turned out to be worth doing.
