## Supplementary information

Supplementary Figures 1–6 show the alignments and structures the main figures
rest on: the representative alignment and the columns the tree actually saw;
the IP₃-binding core and the pore module at residue resolution across the
three paralogues, with each clinically labelled position checked against the
residue at its column; the per-paralogue alignments the constraint map is
computed on; the trimmed codon alignment behind every ω estimate; the
constraint map painted onto the channel; and every labelled variant with its
per-element enrichment test. None of them re-aligns or re-renders anything:
each is drawn from the same committed file as the main figure it supports, so
a supplementary panel that disagreed with its main figure would be a bug
rather than a difference of method.

**Two joins are checked before any of them is drawn, and both are hard
failures.** trimAl writes no map between the untrimmed alignment and the
trimmed one, so the map committed with the alignment is verified by an
exhaustive column walk: all 1,797 trimmed columns against all 134 sequences,
because an off-by-one would renumber every residue claim downstream and would
have no other symptom. And every clinically labelled position is checked
against the residue its own paralogue's per-residue table holds there, and
every aligned partner against the residue the *other* paralogue's table holds
— 1,780 variants and 2,699 aligned partners.

**A structure carries a variant position only if it earns it.** Human *ITPR2*
(9YKK) and *ITPR3* (8TKG) match the human per-residue table at every residue
they share with it; the *ITPR1* cryo-EM reference is a rat structure and the
AlphaFold DB model of human *ITPR1* is the 2,695-residue Q14643-4 isoform, so
neither carries this gene's numbering and *ITPR1*'s pathogenic positions are
not drawn on coordinates that are not theirs. That refusal is the panel, not
a caveat to it.

### Supplementary figure legends

**Supplementary Fig. 1 | The representative alignment, and the columns the
tree actually saw.** (**a**) The 134 representative proteins over the 11,777
columns MAFFT L-INS-i produced, each cell shaded by the fraction of its
column block that is residue rather than gap, rows ordered and side-barred by
group. Beneath it, the 1,797 columns trimAl kept — the alignment the tree,
the selection tests and the constraint map were all computed on. The map
between the two coordinate systems was verified by an exhaustive column walk
before this figure was drawn. (**b**) Column occupancy along the input
alignment. (**c**) The same occupancy as two distributions, kept against cut:
the median kept column is 0.98 occupied and the median cut column 0.03, so
what trimAl removed was the sparse columns and not a region.

**Supplementary Fig. 2 | The ligand core and the pore module at residue
resolution across the three paralogues.** (**a**) Deep-layer constraint along
the β-trefoil and MIR domains of each paralogue in its own numbering, with
every pathogenic position marked. (**b**) The same for the channel region —
selectivity filter, gate, and the geometrically located luminal loop drawn
beside them rather than inside the pore module, since whether it counts as
pore is what reverses the comparison in the main text.
(**c**) Every pathogenic position in those two regions with the residue each
paralogue carries there; a red letter is a paralogue that carries a different
residue. Twenty of the 22 positions are the same residue in all three. Every
letter was checked against that paralogue's own per-residue table.

**Supplementary Fig. 3 | The within-paralogue alignments the constraint map is
computed on.** (**a**) How many of each paralogue's 249–265 sweep orthologues
cover each residue of the human reference. (**b**) The shape screen that bait
coverage cannot do: the fraction of each sequence's own residues that land in
a reference column, with the calibrated bar drawn and the one sequence it
dropped ringed. (**c**) The four conservation layers, with the shallow
msa_v2-only set beside each deep set as the control for what the depth bought.

**Supplementary Fig. 4 | The trimmed codon alignment behind every ω estimate.**
(**a**) Per-codon occupancy of the 3,253-codon PAL2NAL alignment over 57 tips;
trimAl kept 2,459 codons. (**b**) Per-tip coverage of the trimmed alignment,
grouped by the selection set the tree placed the tip in. (**c**) Codons masked
to NNN by validation, by the route the coding sequence came from: the genome
gene models carry the masking, the UniProt records almost none.

**Supplementary Fig. 5 | The constraint map painted onto the channel.**
(**a**–**c**) One subunit of each paralogue's cryo-EM reference, projected on
its two principal axes and coloured by deep-layer constraint from the file
S17 painted. (**d**) The same *ITPR3* structure painted with the selection
layer instead; 59 % of its scored residues sit at the ceiling, having no
non-synonymous substitution anywhere in the tree. Residues with no score are
grey and never the low end of the scale, because the luminal loop is both the
least conserved element and the worst resolved.

**Supplementary Fig. 6 | Every labelled variant, and the per-element
enrichment test.** (**a**, **b**) Pathogenic and benign missense positions on
the two structures whose numbering the check admits. (**c**) The check: the
fraction of each candidate structure's residues whose amino acid matches the
human per-residue table. Human *ITPR2* and *ITPR3* match at every residue;
the rat *ITPR1* reference and the Q14643-4 isoform AlphaFold DB serves for
human *ITPR1* do not, so that gene's pathogenic positions are not placed.
(**d**) Pathogenic positions per element as an odds ratio against the rest of
the protein, Benjamini–Hochberg corrected across all 15 element tests.

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

