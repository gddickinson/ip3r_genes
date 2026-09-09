## Extended Data figure legends

**Extended Data Fig. 1 | Copy number and controlled absence outside the
vertebrates.** (**a**) Complete IP₃ receptor gene models per genome across the
193 controlled non-vertebrate genomes, grouped by kingdom, with the vertebrate
paralogue count marked. Most fungal and plant genomes carry none; most
non-vertebrate metazoan genomes carry one; nine genomes carry more than three
and four carry six or more, the largest being *Macrostomum lignano* at 18. (**b**) Every absence clade taken to assembly level, the 22 with the
largest proteome denominators shown: the open outline is genomes searched, the
grey bar is genomes in which a measured positive control recovered a
comparably long, deeply conserved gene, and the blue bar — invisible because
it is zero everywhere — is genomes with an IP₃ receptor. The grey bar filling
the outline in every row is the result. (**c**) All 99 plant and fungal
records chased individually: verdict per record, and each record's identity to
its nearest protein outside its own kingdom, sorted. Every record sits between
19.9 % and 45.8 %, far below the 95 % that marks an assembly contaminant and
below the 80 % that would flag one as a cross-kingdom outlier; both lines are
drawn.

**Extended Data Fig. 2 | The vertebrate genomic sweep and its contiguity
confounder.** (**a**) What the sweep found for each paralogue and for the
ryanodine receptor positive control across all 309 genomes. (**b**) Recovery
against assembly contiguity, binned on the contiguity bar rather than across
it, and the same comparison as a two-bar contrast either side of it: 98–99 %
of cells are found above the bar against 57–70 % below it, and below the bar
the ordering is *ITPR3* > *ITPR1* > *ITPR2*, which is the order of their gene
spans. (**c**) Loci found per genome for each paralogue, with the ryanodine
receptors beside them: the *ITPR1* second copy in the ray-finned fish is the
only substantial multi-copy signal in the family, against two to six copies
for the control, whose three genes the sweep most often recovers as two or
three loci.

**Extended Data Fig. 3 | The representative alignment every downstream result
stands on.** (**a**) Per-column conservation of the 1,797-column trimmed
alignment, drawn as a rolling mean of 25 columns with the alignment-wide mean
dashed, and human *ITPR1*'s Pfam architecture mapped through the alignment
rather than scaled onto it, the domains coloured by whether they are shared
with the ryanodine receptors or are the generic pore. (**b**) All-against-all
identity over mutually covered columns for the 134 representatives, ordered
and side-barred by group; the three bright blocks at top left are the
vertebrate paralogues, and the whole cross-family field is dark. (**c**)
Per-sequence coverage of the trimmed alignment by group. Median coverage is
0.96 and one tip of 134 falls below half.

**Extended Data Fig. 4 | The sister test, node support and paralogue
placement.** (**a**) Each of the three rooted sister hypotheses as a
constrained maximum-likelihood tree: the log-likelihood cost against the best
tree, and the approximately unbiased p-value over 10,000 RELL replicates.
*ITPR1*+*ITPR2* and *ITPR1*+*ITPR3* each cost 114 log-likelihood units and are
rejected; *ITPR2*+*ITPR3* costs 0.5 and is not. (**b**) SH-aLRT against
ultrafast bootstrap for all 131 internal nodes, with the joint threshold
region shaded and the nodes any claim in the paper rests on ringed; 91 of 131
nodes (69 %) clear both. (**c**) The eleven vertebrate tips the tree declines
to place inside a paralogue clade, the support of the clade that does place
each, and that clade's composition. Reciprocal best hits upheld the database
name for all five names the tree disputes, so these are tree uncertainty
rather than annotation error.

**Extended Data Fig. 5 | Each paralogue has its own genomic neighbourhood.**
(**a**) The human *ITPR1*, *ITPR2* and *ITPR3* neighbourhoods drawn as gene
tracks at ±10 coding genes, with the two gene families shared between two
IP₃ receptor loci linked, and the prevalence of each shared family across
species against its random-window background beneath. (**b**) Mean Jaccard
similarity of flanking-gene symbol sets by pair class, one locus per species,
each against its own matched random-window control (black diamonds).
Within-paralogue classes run 173- to 413-fold above the null; every
cross-paralogue and cross-family class is at the null. (**c**) The consensus
paralogue caller: how many consensus keys a locus shares, for annotation-
confirmed loci against random windows, and the threshold sweep that chose the
operating point by maximising call rate minus false-call rate. (**d**)
Neighbourhood conservation within against across vertebrate classes.
*ITPR3*'s neighbourhood is the one that does not travel.

**Extended Data Fig. 6 | The teleost genome duplication and the ancestral
block.** (**a**) Mean gene copies per genome for the three paralogues and the
ryanodine receptor control in the pre-3R ray-finned outgroups, in the
teleosts, and in the lineages with a further whole-genome duplication; and
the fraction of genomes carrying more than one copy, by vertebrate class.
(**b**) For each two-copy teleost genome, the number of ancestral-block gene
symbols each *ITPR1* copy retains: points off both axes are genomes in which
both copies keep part of one block, and the bar chart counts genomes meeting
that criterion and the stricter one that the two copies' symbol sets are
disjoint, against a tetrapod and a pre-3R ray-finned reference. (**c**) The
cross-anchor assignment margin per reference genome, and the copy counts
recomputed at seven coverage bars.

**Extended Data Fig. 7 | Dating the duplications that made *ITPR1*, *ITPR2*
and *ITPR3*.** (**a**) The hand-curated, literature-calibrated species tree on
a linear time axis, with each reconciliation variant's duplication placement
marked at the node it maps to and the published age spread drawn as a band on
the nodes a placement uses. Every internal node carries a calibration; only
those two are banded, because a band on all 29 would obscure the tree. (**b**) The full matrix — five topologies by three tip
variants — showing the deepest paralogue duplication in each cell, and the
distribution of total events over all 112 rootings of the vertebrate subtree
with the outgroup rooting and the minimum-event rooting marked. (**c**) What
each implied loss turns out to be when asked of the genomes — 51 or 53 in the
ten variants that keep the tree's own resolution, 83 and 102 in the two that
collapse its unsupported nodes: none is corroborated. (**d**) Root-to-tip distances for all 57 vertebrate
tips with the six cyclostome loci marked — two of them 0.0004 substitutions per
site apart and drawn as one line — and each of those loci beside the
independent flanking-gene call for it: four fall inside their own null and two
have too few informative neighbours to call at all.

**Extended Data Fig. 8 | The loss instrument, and what it takes to
manufacture a loss.** (**a**) The calibration: how much of a reference protein
is reassembled outside every placed locus, for candidate genes, for paralogues
that are themselves shattered in the same genome, and for the decoy — a
paralogue whose gene the aligner has already placed elsewhere in that genome.
The two informative populations separate completely; the bar sits at the
midpoint of the gap and both gap edges are committed. Beside it, how far each
undecided cell's reference is spread across contigs, which is what that
coverage is a coverage of. (**b**) Lesion density
against bait identity and against contig N50, showing that the confounder is
the alignment and not the assembly, and the within-genome paired sign tests
that follow from it. (**c**) Why the synteny route could not be used: 273 of
432 rescue regions sit on a contig carrying no annotated gene, while the
caller is accurate at every key count it can act on. (**d**) The sensitivity
matrix: Dollo losses under every combination of coding, evidence bar and the
two decision rules. Zero is drawn as an explicit zero, never as an empty cell.

**Extended Data Fig. 9 | Selection across the three paralogues.** (**a**)
One-ratio ω per paralogue on a logarithmic axis, with the estimate from
curated coding sequences only beside each. (**b**) Whole-tree two-ratio
contrasts, each clade against the rest of the tree, and the RELAX selection
intensity parameter k for the same three contrasts. (**c**) dN against dS for
every within-paralogue pair on log–log axes with the neutral diagonal and the
saturation bar drawn; 85–94 % of pairs are past the bar, which is why every ω
quoted is tree-based. (**d**) Every branch-site restart against its own
nested null: points left of the line reached a lower optimum than the null
they are tested against and are local optima, not results — one per stem, at a
different starting ω each time.

**Extended Data Fig. 10 | Predicted and experimental structures across the
family.** (**a**) The panel: sequence length and resolved or
modelled residues for every reference, state-panel member, negative control and
predicted model. Twenty-nine of the thirty are long enough to score; the
thirtieth is the record AlphaFold DB serves for human *ITPR2*, a 181-residue
isoform, drawn here because that is the coverage result. (**b**) The structural family call: best TM-score against
an IP₃ receptor reference against best against a ryanodine receptor reference,
with TM-align's own random and same-fold bars drawn rather than described;
open symbols are declined for falling below the fold bar, and marker shape is
the role: triangle experimental, circle predicted, square control. Beside it, what a
TM-score means on this panel — conformation is worth 0.22, and no control pair
reaches the fold bar. (**c**) Mean AlphaFold confidence per domain: the
IP₃-binding core is the best-modelled domain and the pore the worst. (**d**)
AlphaFold DB coverage of the census by group, and against record length, where
the modelled mass sits below ~1,300 residues and the peak at a full-length
subunit is almost entirely unmodelled.

**Extended Data Fig. 11 | Constraint by element, at the ligand site, and as a
variant classifier.** (**a**) Per-element constraint in all three paralogues on
a divergence metric and on a composition-free one, with the linker control's
mean over the three proteins drawn as one dashed line per panel; the two
metrics agree, including on the luminal loop.
(**b**) The ten measured IP₃ contacts, the two filter-lining residues and the
two gate-lining residues against two controls — the whole protein and the rest
of their own elements — and, beside it, between-paralogue identity per
element, where the gate is at 1.00 and the luminal loop at 0.13–0.31.
(**c**) Receiver-operating curves for all four constraint layers as
classifiers of pathogenic against benign missense variants on one fixed set of
positions, and where the uncertain variants fall relative to the pathogenic
median on each layer.

**Extended Data Fig. 12 | The annotation audit against its own controls.**
(**a**) The completeness bar validated rather than re-derived: over 1,077 loci
the annotation names correctly, the best single model's share of the gene, with
the inherited 0.50 bar marked at that distribution's 1.3 % point; and every
state recounted across bars from 0.30 to 0.95. (**b**) Annotation state for
each paralogue and for the ryanodine receptor control, and the family against
the control with and without the contiguity restriction: the difference does
not survive correction either way. (**c**) What the protein databases call the
11,402 full-length family records, by gene symbol and by protein name, and the
resolution of the 15 reference proteomes that returned nothing — all 15 are
gene-caller failures.

**Extended Data Fig. 13 | Two annotation failures validated at the exon and by
RNA-seq.** (**a**) Exon tracks at true genomic width for the two validated
cases, with each exon coloured by what the annotation holds at that interval
and the annotated models drawn beneath. *Nibea albiflora ITPR2* has no
annotated gene on any of its 56 exons; *Dissostichus eleginoides ITPR3* has 23
of 60 exons with no annotated gene and is called as three separate models.
(**b**) The denominator: annotation loss across all 382 eligible loci, and the
eligibility rule that does the work — whether the same annotation builds genes
that long anywhere else in the same genome. (**c**) The checks on the
alignment evidence itself: exon-boundary concordance with independently
annotated genomes, splice dinucleotides, and internal stop codons observed
against the number expected under neutral drift. (**d**) Every junction of
every locus in the expression panel, at its position in the spliced coding
sequence, coloured by whether an annotated model spans it and marked where no
read crossed it.

**Extended Data Fig. 14 | What the search methods were worth.** (**a**) The
false-negative rate against contig N50 for the IP₃ receptor cells and for the
independent ryanodine receptor sister series, binned on the contiguity bar;
the residual rate at every possible floor, with the 5 % line drawn; and what
each floor costs in genomes retained. (**b**) The bait-panel ablation as a
change from the full 38-bait panel on a symmetric-logarithmic axis, and the
recall of a single bait against its identity to the target. (**c**) Where a
profile HMM earns its place — almost entirely below 1,000 residues in the two
well-sampled groups, and at every length in the protists; the fungal and plant
panels above 1,000 residues rest on a dozen to two dozen records each — and the
fraction of demonstrated genes that no protein database holds. (**d**) The kill criterion written for iteration drift, measured: the
sister-family share it acts on, the off-family share that actually moves, and
each rule scored as a classifier of the outcome over seven runs.

