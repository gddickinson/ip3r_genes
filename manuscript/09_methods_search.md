## Methods — search and scope

All analyses are scripted, and every table, figure and report in the deposit
is regenerated from committed inputs by a named script. Tool versions are
recorded in `results/toolchain_manifest.txt`: MAFFT v7.526, HMMER 3.4,
BLAST+ 2.16.0+, miniprot 0.18-r281, trimAl v1.5.rev1, IQ-TREE 2.3.6, PAML
4.10.10, pal2nal v14, HyPhy 2.5.101, TM-align 20240303, HISAT2 2.2.3,
sra-tools 3.4.1, NCBI datasets 18.35.0. Where a tool's output depends on
thread count it was run at a pinned thread count and the command line,
version, runtime and SHA-256 of every input and output recorded; MAFFT
L-INS-i and IQ-TREE are not reproducible on this machine at automatic thread
selection, and both were pinned for that reason.

**The family definition** lives in one module: the Pfam signatures PF08709
(IP₃-binding core), PF02815 (MIR), PF01365 (RIH) and PF08454 (RIH-associated),
the 2,000–4,000 aa size band, the human reference and sister-family panels,
and the known-name substrings. Nothing else in the code names a gene.

**Separating the two families is a positive test at every stage.** Because
RYR1 carries all four IP₃ receptor-diagnostic Pfam signatures, no stage
assumes a family: the census call requires the *complete* IP₃ receptor
architecture with none of the ryanodine receptor-specific signatures, with
length as support only; profile assignment requires the winning profile to
beat the loser by more than 10 % of its own bit score; genomic loci are
assigned to a family by alignment score before any paralogue question is
asked; and structures require the winner to clear TM-align's own 0.50
same-fold bar before a family call is made at all. The rule is audited
against evidence it never sees — gene symbols for the census call
(6,191 labelled records, 0 disagreements), the architecture call for profile
assignment (11,875/11,876 agreement outside the seeds).

**Protein-level enumeration.** The four family signatures were walked to
exhaustion through the InterPro API with every raw page archived and the
census re-parsed from the archive, then joined to one streamed UniProt query
carrying full Pfam architecture, length and lineage. A PF08709-only census
would have missed 2,914 records, 758 of them called IP₃ receptor across 385
taxa, so the enumeration uses the union of the four.

**Profile sweeps.** Two profile HMMs were built from census-derived seed
sets under five enforced selection rules — a seed whose census call,
architecture or length has drifted since the manifest was written aborts the
build — and swept over 763 vertebrate reference proteomes (14,414,821
proteins, 6.97 G residues) and, in a second sweep, 6,928 archaeal, bacterial
and non-vertebrate eukaryotic reference proteomes (63,144,898 proteins,
24.93 G residues) whose four eukaryote groups partition Eukaryota against
the vertebrate set. Every search ran at a relaxed reporting threshold with
the primary E ≤ 10⁻⁵ call taken as a filter on the same output, so the two
sensitivities come from one search; the equivalence of that design was
tested rather than assumed, and the test found two real effects (conditional
E-value normalisation, and two-significant-figure printing at the threshold)
neither of which moves a call. A match spanning fewer than 200 profile match
states is counted but not enrolled, a floor measured from this family's own
domain coordinates; without it 12,874 vertebrate proteins carrying only a
SPRY domain were called ryanodine receptor.

**Completeness.** jackhmmer was iterated to convergence from single
sequences in seven runs across vertebrate, invertebrate, protist, plant and
fungal databases, with three coded kill rules evaluated per round on that
round's own inclusion list. In all four non-vertebrate groups the iterated
model settles on exactly the single pass's IP₃ receptor count, so iteration
finds no receptor the single pass missed. Completeness was also measured in
the other direction: each protein-census record the sweep failed to return
was looked up in the 8.8 GB sequence database itself, and none was present
and missed.

**Genome scope.** The vertebrate denominator is 309 assemblies, declared
before the sweep: the best assembly of each of 161 vertebrate orders, ranked
by annotation status, RefSeq category, assembly level and scaffold N50,
united with 169 margin species derived from the committed census tables by
four stated rules (no protein hits in the reference proteome; fewer than
three paralogues; only fragmentary records; literature anchor). The
non-vertebrate denominator is 194 assemblies chosen by five rules from
24,596 NCBI eukaryote reference assemblies less the vertebrates, under one
principle — sample most finely where the negative claim is — with
per-phylum and per-class representatives, one per class in every clade the
proteome sweep found empty, 14 literature anchors and the highest-copy
species per class.

**Genomic sweep.** A 38-bait panel (30 IP₃ receptor, 8 ryanodine receptor as
positive control) derived from the census by seven enforced rules was
aligned to each vertebrate genome with miniprot; empty cells were re-asked
by whole-genome tblastn and attributed family-first. The maximum intron
setting is a threshold on the call, not a performance knob — too small a
value splits a gene and a split gene reads as a fragment — so it was
measured: no IP₃ receptor gene in an 11-species panel has an intron over the
200 kb default (widest 152,216 bp) while the ryanodine receptor control
exceeds it twice. The rescue attribution margin inherited from a related
project (0.333) was re-measured at loci whose paralogue the assembly's own
annotation establishes, where 7 of 9 complete loci fall below it, and reset
to 0.22; the reset was then validated on the evidence it acts on, 53 of 53
rescue fragments agreeing with the annotation over margins 0.227–0.445. The
non-vertebrate sweep uses a 108-bait panel with per-group intron and
contiguity settings measured the same way, and identity to the nearest bait
was retired as a call gate there after it was shown not to separate
confirmed from contradicted loci (Youden J 0.71); the profile call carries
it instead.

**Positive controls for absence.** In the vertebrates the ryanodine
receptors are the control and fired in all 309 genomes. Outside the animals
they are not a control, because architecture-level ryanodine receptor occurs
in 2 of 6,928 non-metazoan proteomes; six candidate control profiles — all
large, deeply conserved, multi-exon eukaryotic families — were therefore
measured over each clade's own swept proteomes and each clade takes the one
that is present in it. A genome that returns neither its own control nor a
cross-kingdom one supports no absence claim.

**Contiguity.** The contiguity bar is the median measured IP₃ receptor gene
span, 142,212 bp: an assembly whose contigs are shorter cannot carry the
gene on one contig. It was chosen from gene geometry with no error rate in
it and calibrated only afterwards, where it gives a 0.9 % residual
false-negative rate on the IP₃ receptor series and 0.0 % on the ryanodine
receptor sister series while retaining 189 of 309 genomes. Every result is
reported above and below it.

