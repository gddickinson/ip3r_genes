## Extended Data

**{fig:enumeration}.** The search space was enumerated by three signatures
because no single one reaches the whole family. (**a**) The size of each
signature's protein set and of every overlap between them is drawn, with each
region broken down by what the architecture rule calls it; a census built on
PF08709 alone would miss 2,914 proteins. (**b**) The census is broken down by
lineage with both family calls; the vertebrate bar is the tallest because
sequencing has concentrated there, which is why range is measured in
proteomes rather than records.

**{fig:margin}.** The labelled-bait margin leaves an empty gap between the two
families across the whole census. Each record's margin between its identity
to the nearest labelled IP₃ receptor bait and to the nearest labelled
ryanodine receptor bait is drawn with the no-call band marked. No record
falls between the narrowest correct call and the narrowest correct rejection,
so the 0.10 margin is not a compromise between sensitivity and specificity.

**{fig:census_growth}.** An exhaustive domain enumeration found what a search
by name had not, and length alone could not have made the call. (**a**) Each
taxonomic group's records are split into those a name-based search had
already returned and those the enumeration added, and outside the vertebrates
the name-based share is close to zero. (**b**) The two called families
separate by length, while the records neither rule can call are short
fragments far below both.

**{fig:instruments}.** Two instruments reading different evidence speak for
most of the census, and gene-set size confounds per-proteome copy number.
(**a**) Records are split by whether the architecture rule, the profiles,
both or neither call them, and the profiles' contribution is almost entirely
records the architecture rule could not decide. (**b**) Receptor calls per
vertebrate proteome are plotted against the size of its gene set, and the
fifteen proteomes with no call sit among the smallest gene sets.

**{fig:completeness}.** Iterative searches from distant seeds recover the same
family and differ only in what is not family. (**a**) The three vertebrate-database runs,
seeded from a human, a fly and an amoebozoan receptor, are traced by new
targets per round and by their ryanodine receptor share. (**b**) The same traces are
shown for the runs seeded inside each non-vertebrate eukaryotic group. (**c**) Each
kill rule is scored as a classifier of drift measured on the finished model: the
rule written for this family's hazard fires on none of the drifted runs, and
the same threshold applied to the off-family share separates all seven.

**{fig:separation_outside}.** The two-profile separation holds outside the
vertebrates. Each target scored by either profile in the non-vertebrate sweep
is plotted by its two scores, with the no-call band drawn; of 2,769 targets
scored by both profiles above the floor, one falls in the band.

**{fig:genome_instrument}.** Two genomic thresholds were measured for this
scope because neither transfers from the vertebrates. (**a**) The identity of
every recorded locus cluster, with loci confirmed and contradicted by the
assembly's own annotation shown separately; identity does not separate them,
so a second axis scores every locus against both profiles. (**b**) Locus span
against the coding footprint inside it, by group; the ratio differs between
metazoan and protist genomes, so the intron setting and the contiguity bar
are set per group.

**{fig:structures}.** Structural evidence exists for a small part of the
family. (**a**) Every structure in the panel is
drawn with the length its record claims and the length it resolves. (**b**) Prediction confidence is shown per
domain in subunit order, with the IP₃-binding core the best-modelled domain and the
pore the worst of the named ones.

**{fig:afdb}.** AlphaFold DB holds almost no full-length model of this
family. (**a**) The modelled fraction of the census is drawn per taxonomic
group. (**b**) Modelled and unmodelled records are drawn against record
length; the models sit mostly below about 1,300 residues, while the peak at
about 2,700, a full-length subunit, is almost entirely unmodelled.
