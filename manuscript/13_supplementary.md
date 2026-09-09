## Supplementary information

Supplementary Figures 1–6 show the alignments and structures the main figures
rest on: the representative alignment and the columns the tree actually saw;
the pore module and the IP₃-binding core at residue resolution across the
three paralogues, with each clinically labelled position checked against the
residue at its column; the per-paralogue alignments the constraint map is
computed on; the trimmed codon alignment behind every ω estimate; the
constraint map painted onto the channel; and the structural results panel.
None of them re-aligns or re-renders anything: each is drawn from the same
committed file as the main figure it supports, so a supplementary panel that
disagreed with its main figure would be a bug rather than a difference of
method.

*These figures are produced by the supplementary-figure stage of the analysis
and are not part of this build; the figure manifest lists zero supplementary
files, and the reviewer checklist records this as an open item.*

**A resource, offered with its limit stated.** The per-residue constraint
tables carry all four conservation layers, the site-wise selection rate, the
structural element and the measured-site flags for every residue of the three
human paralogues, and `variants.tsv` joins them onto all 1,753 harvested
ClinVar records. `vus_stratification.tsv` places each of the 1,546 uncertain
variants against its own gene's labelled distributions. That is a
stratification and not a call: it says where a variant sits on an axis that
separates the labelled ones, and nothing about the variant's effect. It is
offered as evidence to combine with others, not as a classifier to act on.

Supplementary Tables are the committed analysis tables themselves, listed with
their sizes and checksums in `manuscript/deposit_manifest.tsv`. The tables a
reader is most likely to want are: the genome manifests (the two declared
denominators), the census (`census_v6.tsv`, one row per record with both
instruments' verdicts), the per-genome ledger, the character matrix, the
per-residue constraint tables for the three human paralogues, the correction
list, and the per-task statistics files, each of which carries the parameters,
the self-test status and the SHA-256 of every table its task wrote.

