## Figure legends

**{fig:census}.** The three IP₃ receptor paralogues are found in nearly every
vertebrate genome searched. Each panel is one paralogue, each row a
vertebrate class, and the bar stacks the evidence for the gene as a fraction
of that class's genomes, from a full-coverage locus through truncated and
fragmentary evidence to the absent call. The red that remains sits in a few
classes, and those are the classes with the poorest assemblies, so this raw
ledger is the problem the later figures resolve rather than a map of loss.

**{fig:false_negatives}.** The search's miss rate falls with assembly
contiguity in two independent control series. (a) The share of
known-present cells the ledger missed is plotted against contig N50, for the
IP₃ receptor series and for the ryanodine receptor series, which uses none of
the state rules. (b) The residual miss rate is scored at every candidate
contiguity floor, with the a-priori bar marked. (c) Each floor removes
genomes, and the panel counts them. Agreement between the two series is what
stops the rate being read as the search agreeing with itself.

**{fig:panel}.** Paralogue coverage decides the panel's recall and
phylogenetic breadth barely matters. (a) Each of 18 ablated panels is drawn
as its change in recovered cells from the full 38-bait panel on a
symmetric-log axis; only the ablations that remove a paralogue's own baits,
and the one that removes every labelled bait, move far. (b) A single bait
recovers its cell at any identity to the gene above 0.5. The simulation
behind both panels reproduces the committed ledger cell for cell.

**{fig:matrix}.** No genome × paralogue cell is absent once every alternative
explanation is tested. (a) All 927 cells are drawn in the state their
evidence supports, with genomes ordered by contig N50 and the contiguity bar
marked. (b) Summing genomic sequence per genome gives the gene-equivalents
each assembly holds; above the bar every genome holds three. The absences of
the raw ledger have become truncations, reassemblies and four unassignable
cyclostome cells.

**{fig:reconstruction}.** A calibrated bar separates a shattered gene from
cross-paralogue similarity. (a) Reassembled coverage is shown for the
candidate regions, for the decoys attributed to a paralogue already placed
elsewhere in the same genome, and for the co-shattered regions between them,
with the gap shaded and the operating point drawn. (b) Every undecided cell
clears the bar however many contigs its gene is spread over. The third
population was not anticipated, and separating it is what made the
calibration separate at all.

**{fig:sensitivity}.** A family-level loss is almost impossible to
manufacture from this data and a paralogue-level one is easy. (a) The Dollo
loss count is given for the family-level coding in every setting of the
evidence ladder crossed with the cyclostome and contiguity rules. (b) The
same grid is shown for the paralogue-resolved coding. The operating point is
ringed in both, and a zero is printed as a zero because an empty cell would
read as a setting nobody measured.
