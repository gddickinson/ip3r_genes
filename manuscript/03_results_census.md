### Three paralogues, found in the DNA of 309 vertebrate genomes

The vertebrate scope is 309 assemblies — one per vertebrate order plus every
species the protein record left in doubt — 552.5 Gbp, declared before the
search and fixed. A 38-protein bait panel (30 IP₃ receptor, 8 ryanodine
receptor) was aligned to each genome with a spliced aligner, and cells the
alignment left empty were re-asked with a whole-genome translated search.
All 309 completed with no failures; the ryanodine receptor positive control
fired in every one of them, so no genome is excluded on control grounds.
The sweep recorded 1,236 genome × cell results over 2,144 loci (Fig. 2 and
Extended Data Fig. 2).

**The sister family never once contested a locus.** Across all 2,144 loci
only one family's baits aligned at all. That separation is sharper at the
genome than at the protein level, where a benchmark against curated decoys
had promoted all six ryanodine receptors as novel IP₃ receptor paralogues
until an explicit sister test was added; the family ambiguity is a property
of comparing protein fragments, not of the two families.

**Only four cells in 1,236 were called absent, and all four are
cyclostome.** *Petromyzon marinus* and *Myxine glutinosa* each carry
*ITPR1* and no *ITPR2* or *ITPR3*. Both assemblies are contiguous enough to
carry the gene and both have a firing control, so the call survives the
usual gates — but the bait panel has no cyclostome-labelled *ITPR2* or
*ITPR3* bait at all, and each of these genomes in fact carries three
full-length family loci that all fall to the *ITPR1* bait. These are
therefore cells whose paralogue cannot be assigned, not absences, and we
treat them as such throughout.

**The genomes hold genes the databases do not.** Adding the sweep's gene
models to the census contributed 1,058 models from 224 genomes. Of the IP₃
receptor models, **318 exist only as DNA** — no protein record of any kind —
and **167 more sit inside an annotated gene that carries no family name**,
so they are unreachable by any search that starts from a gene symbol.

**How much the search missed, measured rather than estimated.** Because no
loss is reconstructed anywhere in this scope (below), every cell whose gene
is independently known to be present and which the ledger did not grade
*found* is a false negative of the method. That is 140 of 923 cells,
**15.2 %**. The ryanodine receptor sister cell — present in every
vertebrate, swept by the same aligner in the same assemblies, and using none
of the same evidence — gives 42 of 309, **13.6 %**, indistinguishable from
it (Fisher p = 0.58), so the rate is a property of the method and not of the
paralogue calls.

Every miss is an assembly. A missed cell's median contig N50 is 23,460 bp
against 3,396,515 bp for a found one; the odds of finding the gene rise
8.1-fold per tenfold of contig N50 (19.99-fold for the ryanodine receptor
series); and chromosome-level assemblies miss 3 of 512 IP₃ receptor cells
and 0 of 172 ryanodine receptor cells (Extended Data Fig. 3). We therefore
report every downstream result twice, once over all 309 genomes and once
over the 189 whose contigs are longer than the median measured gene span
(142,212 bp). Above that bar the false-negative rate is 0.9 % on the IP₃
receptor series and 0.0 % on the ryanodine receptor series — a bar chosen
from gene geometry with no error rate in it, and calibrated only afterwards.

**Bait-panel design.** Re-running the whole cell-assignment chain over 19
ablated panels reproduces the committed ledger 1,236 of 1,236 cells, so the
ablation is exact rather than modelled. Phylogenetic breadth is nearly free:
four human baits recover 782 of the 783 cells the full 38-bait panel
recovers, and dropping any single clade band costs at most two cells. What
matters is paralogue coverage — dropping one paralogue's baits costs
238–241 cells — and one bait alone recovers the gene at any identity above
0.5.

