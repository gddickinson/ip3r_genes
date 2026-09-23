## Methods

### Scope

The question is asked of three declared spaces. The proteome space is 763
vertebrate UniProt reference proteomes (14,414,821 canonical proteins) and
6,928 further reference proteomes (63,144,898 proteins) in six groups:
non-vertebrate Metazoa, Fungi, Viridiplantae, other protists, all archaea, and
bacteria sampled as the proteome with the most proteins in each genus
[R157]. The four eukaryotic groups partition Eukaryota with the vertebrate
set, so the denominators add. Canonical proteomes are used because isoform
sets add isoforms of genes already counted. Proteomes listed by UniProt but
absent from the release tree were excluded and recorded rather than retried.
The enumeration space is every UniProtKB protein carrying PF08709, PF01365 or
PF08454, walked through the InterPro API to cursor exhaustion. The genome space
is 194 non-vertebrate assemblies (100.4 Gbp) drawn from 24,596 NCBI eukaryotic
reference assemblies less 6,216 vertebrate ones by the five rules given in
Results, with every rule's selections committed.

### The family call

**The benchmark was run before any search.** The discovery scorer was run
against 25 positive and 31 negative controls fetched from UniProt (Swiss-Prot
preferred), with MAFFT [R138] alignments whose return code was checked, since
a silent MAFFT failure degrades to a star alignment. The labelled-bait margin
compares a candidate's identity to the nearest labelled IP₃ receptor bait with
its identity to the nearest labelled ryanodine receptor bait and assigns it
only when one wins by more than 0.10; identity is also recorded over mutually
covered columns only. The negative controls, the six ryanodine receptors
included, are all rejected.

**The architecture rule reads annotation.** A record is a ryanodine receptor
if it carries PF02026, PF06459, PF21119 or PF00622, and an IP₃ receptor if it
carries all of PF08709, PF01365, PF08454, PF02815 and PF00520 and none of those
four. Length enters only as support. PF00622 (SPRY) occurs in many unrelated
proteins, so its use is conditional on a family seed signature and is audited
under that condition. The rule's accuracy is scored against gene symbols it
never reads.

**The profiles read residues.** Seeds were taken from the census with the
census's call, one per species per family, disjoint between the two sets, and
every selection rule is enforced when the profile is built, so a seed whose
call, architecture or length has drifted aborts the build. Seeds were aligned
with MAFFT L-INS-i single-threaded, because multithreaded runs of the same
seeds gave different alignments, and profiles were built with hmmbuild
[R153]. A target is assigned when the winning profile scores at least 30 bits
over at least 200 match states at E ≤ 1e-5 and beats the other by more than 10
% of its own score. The 200-state span gate is set between the shortest
observed PF08709 (200 residues) and the longest observed SPRY (137); without
it a sensitive channel profile calls troponin and helicases ryanodine
receptors. The profiles were calibrated against the architecture call on the
archived census sequences, with and without the seeds, before any sweep.
Both instruments' verdicts are kept on every row and merged by a stated rule:
agreement gives a high-confidence call, one instrument speaking gives its
call, a disagreement is kept as a conflict, and neither leaves the record
unassigned.

### Completeness

Enumerated IP₃ receptor records that the vertebrate sweep did not report were
each looked up in the 8.8 GB search database by one streaming pass, to
separate records never searched from records searched and missed. Iterative
searches were run with jackhmmer [R154] to convergence or to a ceiling of ten
rounds, from three vertebrate-database seeds and from one seed native to each
non-vertebrate eukaryotic group, and each round's included set was recorded
from the raw log. The kill criterion was evaluated per round on that round's
included set, and rounds from its first firing onward were excluded from
every merge. Drift was then measured on each finished model as the share of
its included targets that neither profile scores, and each rule of the
criterion was scored as a classifier of that outcome.

### Two sensitivities for every negative claim

Every non-vertebrate search was run once at E ≤ 10 and filtered at E ≤ 1e-5,
so that the strict and the relaxed sets come from one search. The relaxed
panel adds the four family Pfam models, downloaded from InterPro, to the two
profiles. That the filtered strict set equals a strict run was tested rather
than assumed: no assignment and no gate decision moves, the few domain rows
that differ all score at or below 0 bits, and the test fails if it compares
fewer than 25 targets. Hits were attributed to proteomes through a measured
accession-to-proteome map, because several taxa carry two reference
proteomes.

### The plant and fungal chase

Every plant and fungal record called IP₃ receptor by either instrument was
assigned one of seven verdicts by ordered rules, each writing the number it
fired on: unresolved, module-only (under the 200-state floor), contamination
suspect (at least 95 % identity over half the query to a protein outside its
kingdom, by BLAST [R158]), cross-kingdom outlier (80–95 %), no genome backing,
fragment (under 2,000 residues), and real gene.

### The genome sweep and its controls

Each assembly was searched with miniprot [R164] using a panel of 108 baits
chosen by seven written rules: every bait's label is its clade band and never a
paralogue, length bands and architecture floors are measured per band, and
every bait passes a chimera screen against both profiles. The maximum intron
length was set per group from sixteen measured gene spans taken from NCBI's
own annotations rather than from miniprot, whose setting would otherwise shape
the measurement that calibrates it. Where the IP₃ receptor baits won nothing,
a family-level tblastn rescue searched the whole assembly. Loci were assigned
to a family by the winning bait, graded full, fragment
or scrap, and counted as genes after merging split models whose bait spans are
complementary and splitting clusters that contain two complete models. The
identity floor for a locus was measured against annotation-confirmed loci;
because identity did not separate confirmed from contradicted loci outside
the vertebrates, every locus was also scored against both profiles.

The positive control for each clade carrying a negative claim was chosen by
running six candidate families (MIR, myosin head, E1–E2 ATPase, SMC
N-terminal, kinesin and AAA) over that clade's swept proteomes and taking the
one present in the most. Each genome receives a verdict: controlled across a
kingdom boundary when a control is recovered and baits from another
kingdom-level group align across it, controlled by target when a receptor is
found, partial when only control fragments are found, or no control bait. An
absence is counted only over controlled genomes, in both the numerator and the
denominator.

### Structures

Candidate structures were enumerated from RCSB on the union of the four family
signatures, because RCSB's annotation of human ITPR3 lacks PF08709.
References were cryo-EM, full length within the family's band, one per
paralogue, with a state panel on the paralogue offering most conformations and
negative controls drawn from the benchmark decoy panel and size-matched.
Every structure was reduced to one chain. AlphaFold DB coverage was measured as
modelled residues over the census length. TM-align [R166] was run on all
pairs, and a structure was called for a family when its best
reference-normalised score cleared 0.5 and beat the other family by a 10 %
relative margin. Model confidence was reported per domain, with boundaries
transferred by pairwise alignment.

### Reproducibility

Every number in this paper is read from a committed table and re-verified
against it on every build of the paper. Raw API responses are archived so that
each parse re-runs offline, and tool versions are recorded in the committed
toolchain manifest.
