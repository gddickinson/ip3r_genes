## Extended Data

**{fig:calibration}.** The completeness bar sits in the thin tail of the
population it was validated against. Panel a is the share of the gene that
the best single annotated model delivers, on a logarithmic count axis, over
the 1,077 loci whose own annotation names the correct paralogue on an
assembly able to carry the gene. The inherited bar at 0.5 is drawn where it
falls, below 1.3 % of the population. Panel b re-counts every IP₃ receptor
locus state as the bar moves from 0.3 to 0.95, and the complete band
narrows only slightly across the whole range. A number in a legend cannot
show whether a threshold falls in a gap or in the middle of a population,
and this panel shows that it falls in the gap.

**{fig:cases}.** The validated cases come from a measured tail, and the
annotation's own proteins show what it delivers at each. (a) Annotation
loss over the 382 eligible loci has 360 at zero, on a logarithmic count
axis. Beside it, the fifth eligibility rule plots each locus's span against
the longest gene its own annotation builds and marks in red the six loci it
excludes. (b) The annotation's own proteins, translated from its GFF3 and
genome, are tiled onto each case's recovered protein. No model delivers any
residue of the *Nibea albiflora* *ITPR2* gene, and the model the annotation
names *ITPR2* tiles onto the *ITPR3* locus instead. Three models deliver
1,609 of the 2,687 residues of the *Dissostichus eleginoides* *ITPR3* gene.

**{fig:case_checks}.** Three independent checks agree on both cases. Panel
a ranks the exon boundaries of each case by the share of independently
annotated genomes that draw them, and that share is above 0.9 within the
first ten boundaries of both loci. Panel b counts the introns of each case
by splice dinucleotide, and all but one of them have a canonical pair. In
panel c, no case locus has an internal stop, and neither does any other
family locus of the same genomes, against 4 to 25 stops expected under
neutral drift. No single check is decisive, but they fail in different
circumstances and none fails here.

**{fig:expression}.** Read evidence rises far above its controls, and it
answers the question the transcript deposits cannot. (a) The pooled read
count on every family locus of the three species stands beside the count
on its own reversed decoy on a logarithmic axis, and a cross marks each
locus the annotation loses. Every decoy but one sits at the floor, and the
exception collects a few dozen reads against several thousand on its
locus. (b) Streamed reads cross most unannotated junctions at every locus,
while the species' submitted transcript records cross none. A species with
ten transcript deposits cannot show that a gene is silent, and reporting
that search as underpowered keeps a database's emptiness from reading as
an absence.

**{fig:gap_coverage}.** Two independent instruments agree about the
sequence the annotation loses. Each bar is the fraction of one locus's
read coverage that falls outside any annotated coding block, and the tick
on the bar is the share of the coding footprint that the exon-level
measurement found no model delivering. One instrument aligns a protein to
a genome and the other maps reads to a spliced reference, so their
agreement rules out the loss being an artefact of either.
