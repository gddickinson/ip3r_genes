## 12. How this review was built, and what failed verification

### 12.1 Method

Candidate literature was harvested from Europe PMC across ~60 topic queries
spanning the sections above, each run twice — ranked by citation count and by
recency — to avoid a purely canonical or purely recent bias. Candidates were
curated by hand; bibliographic metadata for every retained reference was then
fetched programmatically from Europe PMC by PMID, so no author list, title,
journal or year in the bibliography was transcribed manually. The final set is
**137 references — 117 primary research articles and 20 reviews** — spanning
1983 to 2025, with publication type taken from the Europe PMC record rather
than assigned by hand. The reference
table is committed as `results/s0_baseline/references.tsv`, and this document
is assembled from its section files by `scripts/s0_review_build.py`, which
renumbers the stable citation keys into order of first appearance and renders
the bibliography from that table. A cited key with no reference row is a build
error, so the text and the bibliography cannot drift apart.

Database figures quoted in §6.1, §7.1, §7.5 and §11 were re-derived for this
review rather than repeated: InterPro protein counts and taxonomic
distributions from the InterPro API, protein lengths and domain architectures
from UniProt, and gene structure from Ensembl release 15.12. The scripts are
`scripts/s0_db_snapshot.py` and `scripts/s0_gene_structure.py`; the outputs are
the committed tables in `results/s0_baseline/`.

**Figures.** Every figure is generated rather than drawn.
`scripts/s0_review_figures.py` renders each one from committed tables in
`results/s0_baseline/review_figures/`, which four data scripts produce:
domain coordinates from InterPro, structural measurements from PDB 6DQN,
alignments with MAFFT over the committed control panels, and the curated
milestone / regulator / variant tables — whose every citation key is
checked against `references.tsv`, so a figure cannot cite something the
bibliography lacks. Figures are numbered by the build in order of first
appearance, like the citations. Each carries a provenance tag —
*measured*, *computed*, *schematic* or *curated* — which is the `[db]` /
`[lit]` discipline of the planning document applied to pictures. Two
structural results are worth flagging as checks rather than claims: the
narrowest luminal point of the pore falls on the GGGVGD selectivity-filter
motif, and the cytosolic constriction on Phe2513/Ile2517, one helical turn
apart — neither the motif nor the published gate residues were given to
the calculation, which sees only coordinates.

**Limitations.** This is a narrative review, not a systematic one: there is no
pre-registered protocol, no PRISMA flow, and no formal inclusion criteria
beyond topical relevance and the preference for primary sources over reviews
for specific claims. Search was restricted to Europe PMC and to
English-language records. Coverage is deliberately deeper on structure,
evolution and genetics than on cell-type-specific physiology, reflecting the
downstream use.

### 12.2 The audit: four claims that did not survive

This review replaced a planning document whose statements were tagged by how
far they could be trusted. Nineteen atomic claims were checked against the
literature. Twelve were verified as written and three were kept with the
wording qualified; the remaining four are recorded here because a review that
silently drops what it could not confirm is less useful than one that says so.

| Claim as originally written | Verdict | Resolution |
|---|---|---|
| Loci *ITPR1* 3p26.1, *ITPR2* 12p11.23, *ITPR3* 6p21.31 | Confirmed by database | Re-derived from Ensembl; all three exact. Promoted from literature claim to database fact (§6.1). |
| "Each is a ~58–60 exon gene spanning hundreds of kb" | **False as written** | Measured: 57–62 canonical exons, and *ITPR3* spans 76 kb, not "hundreds of kb". Replaced by the measured table in §6.1; the 6.5-fold span asymmetry became open question Q7. |
| "Both families independently expanded to three vertebrate paralogues" | **Not established** | The three-paralogue state of each family is a database fact; the *independence* of the two triplications has no primary source. Moved to open question Q2 (§7.4, §11). |
| Gillespie syndrome "behaves as dominant-negative in the tetramer" | **Incomplete** | Both biallelic recessive and de novo dominant-negative mechanisms exist [R43, R44]. Corrected in §9.1. |

Three further claims are retained but qualified: "most agonist-evoked
Ca<sup>2+</sup> signals in non-muscle cells" is a review-level generalisation
with no measured fraction behind it (§1.2, §5.1); the paralogue tissue-emphasis
summary is qualitative and not a quantitative ranking (§6.5); and the *ITPR2*
anhidrosis claim rests on a single family plus a mouse model (§9.2). One
numerical correction was also applied: the tetramer is **1.3 MDa** as reported
in the primary structural work [R22] — four copies of the 2,758-residue human
subunit give ≈1.25 MDa — and should not be quoted as "~1.2 MDa" unattributed.

### 12.3 Reading this review

Cite the sources, not this document. Where a statement here is a database
figure rather than a literature claim it says so and names the release; those
figures should be re-derived before being quoted, because they will move.
