## Introduction

Inositol 1,4,5-trisphosphate receptors (IP₃Rs) are the ligand-gated calcium
release channels of the endoplasmic reticulum, the step that turns a
receptor-generated second messenger into a cytosolic calcium signal
[R01, R07]. Vertebrates carry three paralogues, *ITPR1*, *ITPR2* and *ITPR3*,
which arose in the two whole-genome duplications at the base of the
vertebrates ({paper:origin}) and differ in tissue distribution and
regulation [R29, R101]. Their individual loss in mice has distinct
consequences: ataxia and seizures
without type 1 [R114], and failure of exocrine secretion when types 2 and 3
are removed together [R115]. In humans, deletions and point variants of
*ITPR1* cause spinocerebellar ataxia [R38, R40, R42], and loss of type 2
function abolishes sweating [R45].

That each paralogue matters in one mammal does not say that each is
indispensable across the vertebrates. Gene loss is a common outcome of
duplication, and families as old as this one usually show it somewhere
[R193]. A lineage that had dropped one IP₃R paralogue would be a natural
experiment in calcium signalling, and finding one requires only a genome
search. The difficulty is not the search but the negative. A survey that
reports that no lineage has lost a gene is making a claim about everything
it failed to find, and its value depends on how often the same search misses
a gene that is really there. Genome assemblies are uneven across the
vertebrates [R195], annotation is less reliable still [R194], and an IP₃R
gene spans 76 kb to half a megabase of genomic sequence, so a fragmented
assembly can hide one without any biological cause.

Here we ask whether any vertebrate lineage has lost one of its three IP₃
receptor genes, and we built the analysis so that the answer could be a
measurement in either direction. The scope is 309 vertebrate assemblies,
declared before the search. Each is searched by spliced protein-to-genome
alignment [R164] with a bait panel that carries the ryanodine receptors, the
IP₃R sister family, as an internal positive control. Every genome ×
paralogue cell is then restated under a vocabulary in which exactly one
state can be counted as a loss, reachable only after every alternative
explanation has been excluded by a positive test. Cells the aligner could not
decide are resolved by reassembling the reference across contigs, against a
bar calibrated on a gene known to be elsewhere.

Two measurements make the result readable. The first is the search's own
false-negative rate, obtained from the same data: because the count turns
out to be zero, every cell whose gene is independently known to be present
and which the search did not find is a miss with a known right answer. The
ryanodine receptor cell, present in every vertebrate and scored by none of
the same rules, provides a second, independent series. The second is a
sensitivity matrix that names the analytical settings under which a loss
would appear, so that the zero is reported with its escape routes rather than
as a single number.
