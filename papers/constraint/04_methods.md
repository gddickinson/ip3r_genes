## Methods

### Scope

The rate analyses cover the 57 vertebrate IP₃ receptor tips of the
134-sequence representative alignment, which carries the phylogeny of the
family ({paper:origin}); the six ryanodine receptor outgroup tips are
excluded because synonymous sites do not survive that distance. The
conservation analyses cover the three human paralogues (UniProt Q14643,
Q14571 and Q14573 [R157]) in their own numbering, on four nested layers: a
deep layer of 264, 249 and 265 orthologues of each paralogue, one locus per
genome from a sweep of 309 vertebrate genome assemblies ({paper:retention});
a vertebrate layer of the 57 vertebrate tips; a family-wide layer of all 128
IP₃ receptor tips; and a shallow layer of each paralogue's own
representative tips (19, 13 and 19), kept as the control for what the deep
sets bought. The ryanodine receptors are excluded from every layer, so no
layer scores what the superfamily conserves. The variant analyses cover
every ClinVar missense record on the three genes, 1,753 in all.

### Coding sequences and the codon alignment

Every tip needs a nucleotide sequence that provably encodes its aligned
protein. Candidates came from three routes (Ensembl, UniProt
cross-references and genome gene models from the sweep), and each route was
treated as a generator of candidates rather than an answer, because a UniProt
entry lists every transcript of its gene and for human ITPR1 and ITPR2 the
first listed is not the aligned isoform. Each candidate was translated
against the aligned protein and every disagreement masked to `NNN`, not only
internal stops; 24 codons were masked across the set. Genome gene models were
realigned with miniprot in a mode that keeps frameshift residues [R164] and
accepted only if they reproduced the swept protein. PAL2NAL [R148] built the
codon alignment, which was cross-checked nucleotide by nucleotide against an
independent in-house protein-to-codon mapping, any disagreement aborting the
build. trimAl [R139] columns were chosen on the protein and applied
codon-aware, whole triplets only, keeping 2,459 codons of 3,253. Twelve
constructed negative controls run on every build and check refusal: a
frame-shifted sequence, a sequence encoding a different protein, and a
non-monophyletic foreground must each be refused.

### Selection sets and codon models

The three selection sets are the extended paralogue clades of the maximum
likelihood tree ({paper:origin}), re-derived from the rooted tree and
cross-checked against the committed clade table, because a codeml branch
model given a foreground that is not a clade silently marks a larger one.
Six cyclostome tips that the tree places in no paralogue clade stay in the
whole-tree analyses as background and are in no foreground. With PAML 4
[R145] we fitted one-ratio models per paralogue and on curated sequences
only, pairwise estimates, a whole-tree one-ratio null, a two-ratio model per
paralogue clade, branch-site model A on each paralogue stem [R146], and site
models M1a/M2a and M7/M8 within each paralogue, with Bayes empirical Bayes
posteriors [R147]. Branch-site model A was started from ω = 0.5, 1.5, 2.5
and 4 on every stem, and the best restart is reported; the branch-site test
is read against a 50:50 mixture of χ²₀ and χ²₁ because its null sits on a
parameter boundary. Benjamini–Hochberg correction [R170] runs across all 12
likelihood-ratio tests together. Pairwise synonymous distances above 1.5 are
flagged as saturated. RELAX [R150] in HyPhy [R149] tested each paralogue
clade against the other two, leaving unplaced tips unlabelled rather than in
the reference.

### Per-site selection

FEL [R151] was run on each paralogue's codon alignment and every site
carried onto the human protein by realignment; a site whose amino acid
disagreed after transfer is written with no residue. Because synonymous
sites are saturated, we report the median non-synonymous rate and the share
of sites called purifying at q < 0.05 per element, and take ω only over sites
whose synonymous rate is identifiable. That rate is unidentifiable at 671 of
7,368 scored sites.

### Orthologue sets and conservation

Deep orthologues were taken from the sweep's gene models at bait coverage
≥ 0.8, identity ≥ 0.4 and at least 1,800 aligned residues, one per genome and
paralogue, with lesion-rich loci excluded. Because bait coverage cannot see
a model that is too long, each sequence was also scored on the fraction of
its own residues landing in reference columns, with the bar placed at the
midpoint of the gap in that distribution (0.7133). It removed one chimeric
model, which had been inflating the ITPR2 alignment from 3,380 to 5,676
columns. Alignments were built with MAFFT L-INS-i [R138] at a single thread
for reproducibility. Columns were scored by Jensen–Shannon divergence against
the BLOSUM62 background [R178] after position-based sequence weighting
[R171], with gaps treated as missing data and the occupancy of each score
reported. A composition-free metric, the modal-residue fraction, and a
recomputation on the curated sequences alone are written beside every
element, not only where a result surprised.

### Elements, controls and functional sites

Pfam elements [R155] were taken per accession from InterPro [R156] and need
no transfer. The selectivity filter, the gate and the ten IP₃ contacts
(residues within 4.5 Å of IP₃) were measured on the human ITPR3 structure
6DQN [R24] and transferred to ITPR1 and ITPR2 by pairwise alignment, each
with an anchor test that aborts rather than writing a coordinate: the filter
must arrive on the GGGVGD motif and the gate on the same lining residues. The
luminal loop was defined geometrically on the same structure. Residues
outside every element were assigned to named linkers and termini, which form
the within-protein control; each element is tested against its protein's
linkers by a one-sided Mann–Whitney test. Functional residues were tested
against the whole protein and against the rest of their own elements.
Identity between paralogues was measured over mutually covered columns of
the untrimmed alignment, with the trimmed-alignment value carried beside it.
Fourteen constructed negative controls run before any table is written,
including a mutated filter motif that the anchor must refuse and a domain
that must not enter the linker control; they were mutation-tested on four
deliberate rule breakages.

### Variants

ClinVar missense records were placed only after each cited transcript's own
translated coding sequence was aligned to the canonical sequence and the
reference amino acid required to match; all three genes file on a transcript
whose translation is the canonical, so none was dropped. The two
residue-level variants localised in the literature must be recovered, as a
positive control. Each layer was scored as a classifier by ROC AUC with a
Mann–Whitney test on three contrasts: pathogenic or likely pathogenic
against benign or likely benign, pathogenic against the whole protein (the
control for ascertainment), and the first contrast restricted to the 44
pathogenic and 34 benign positions every layer scores, which is the contrast
the layers are ranked on. Enrichment by element used Fisher's exact test.
Variants of uncertain significance were stratified against the pathogenic and
benign medians of their own gene, thresholds chosen by the labelled
distributions rather than by us.

### Statistics and reproducibility

Element and site comparisons are one-sided Mann–Whitney tests against the
named control; classifier significance is a Mann–Whitney test on the two
labelled sets. Variant classes were collapsed into pathogenic or likely
pathogenic, benign or likely benign, uncertain, conflicting and other, and
only the first two enter a classifier. No test is reported on fewer than
three observations, which is why the two-residue filter-lining and
gate-lining sets carry means only. A p-value is written with its significant
figures and never rounded to zero. MAFFT runs single-threaded because the
L-INS-i result is not reproducible across thread counts, and every alignment,
table and figure is recorded with its SHA-256 so that a rebuild which drifts
is visible in the data.

### Structures

The deep-layer conservation and the inverted FEL rate were written into the
B-factor column of each IP₃ receptor structure in the panel, with unscored
residues set to −1 so that missing data cannot read as low constraint. A
chain mapping under half of itself to its paralogue was refused. A structure
may carry a human variant position only if every residue it shares with the
human table carries the same amino acid.
