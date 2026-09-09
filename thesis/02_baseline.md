# 2. Auditing what the literature could be trusted to say

## 2.1 Why this project began by auditing its own background document

Every genome-scale project begins with a background document, meaning the set
of things taken as known against which the new measurements will be read. That
document is normally written from reading, and it is normally the least
scrutinised artefact in the project, because it is not a result.

This one was treated as a result. Nineteen atomic claims were extracted from
the `[lit]`-tagged statements of the project's background document and checked
against their primary sources. Every `[db]`-tagged number was re-queried
against the live database it came from. The application that would do the
searching was run against all three public interfaces to establish that it
worked at all. The exercise took a session and changed four claims, one of
which would have made a later chapter's conclusion into an assumption.

The tagging scheme is the point. A statement is `[db]` if it can be re-derived
from a database with a stated query, `[lit]` if it rests on a published
source, and `[open]` if it is a question this project answers. A background
document that does not separate these three cannot be audited, because there
is no way to say what would falsify any given sentence.

## 2.2 Four of the nineteen literature claims did not survive as written

Twelve of the nineteen were verified as written against their primary sources.
Three were verified but qualified, meaning they were supported by wording that
overstated the strength of the evidence and were re-written accordingly. The
remaining four did not survive in the form they were stated.

**One claim was promoted from literature to database.** The cytogenetic
locations of the three human genes, with ITPR1 at 3p26.1, ITPR2 at 12p11.23
and ITPR3 at 6p21.31, were carried as a literature claim. They were re-derived
from Ensembl and confirmed exactly, so they became a `[db]` statement with a
query attached rather than a citation.

**One claim was demoted to an open question, and this is the important one.**
The background document asserted that the IP₃ and ryanodine receptor families
had independently expanded to three vertebrate paralogues each. The
three-paralogue state of each family is a database fact. The independence of
the two triplications is not, because no primary source establishes it. Left
in place, that sentence would have made Chapters 6 and 7 and the
reconciliation that dates the duplications into an elaborate confirmation of
something already assumed. It was retagged `[open]` before any tree was built.

**Two claims were corrected.** The first is discussed in §2.3. The second
concerned Gillespie syndrome, which the background document described as
caused by pore-proximal ITPR1 variants acting dominant-negatively in the
tetramer. That is half the genetics. The syndrome arises from both biallelic
recessive variants and de novo heterozygous ones, and only the latter act
dominant-negatively through the channel domain. Stating only the
dominant-negative mechanism would have misdirected the constraint analysis in
Chapter 11, which asks where pathogenic variants sit relative to the
constrained core, and would have had it looking for one signal where there are
two.

![](figures/review_disease.png)

**{fig:review_disease}.** The curated disease-variant map, with each variant
recorded as `point`, `domain` or `gene` according to how far its cited source
actually localises it. That third column is the honesty of the figure. A
source that reports a deletion of the gene cannot be drawn as a residue, and a
map that drew it as one would create the impression of a resolution the
literature does not have. Chapter 11 returns to these positions with a
harvested variant set of its own and asks whether they sit where the
conservation says they should.

## 2.3 The background claim about gene size was false, and it became a measurement

The background document stated that each ITPR is "a ~58–60 exon gene spanning
hundreds of kb". Ensembl's canonical transcripts say otherwise
({fig:review_gene_architecture}).

The exon counts are 62, 57 and 58, which is a range of 57 to 62 rather than 58
to 60. That is a small correction. The span is not: ITPR1 covers 354,174 bp,
ITPR2 covers 497,888 bp, and ITPR3 covers **76,245 bp**. "Hundreds of kb" is
wrong for one of the three, and the genomic span varies 6.5-fold across the
three paralogues while the protein length varies by about 3 %.

![](figures/review_gene_architecture.png)

**{fig:review_gene_architecture}.** The three human genes drawn to one scale
from Ensembl coordinates. Exon count is nearly identical across the three and
genomic span is not, and the difference between them is almost entirely intron
length. The claim that failed the audit is the one this figure makes visible:
a family whose protein length varies by 3 % and whose gene length varies by a
factor of six and a half.

A failed background claim is normally just a correction. This one became a
measurement. If three paralogues of one gene, encoding proteins of the same
length, occupy genomic spans that differ by a factor of six and a half, then
either the exon structure differs between them or the introns do, and either
answer says something about how the family has evolved. That question is
Chapter 8, and it exists because a sentence taken from reading turned out to
be false.

The correction also has an operational consequence that reaches into Chapter
4. A protein-to-genome aligner has a maximum intron length, and a gene whose
largest intron exceeds it is split into pieces, so a split gene reads out of a
genomic sweep as a fragment. Setting that parameter is not a performance
choice, and it cannot be made from a claim about "hundreds of kb" that is
wrong for a third of the family.

## 2.4 Every database number in the background document reproduced exactly

Every `[db]` number in the background document was re-queried on the day the
audit ran. All of them reproduced exactly, which is the boring outcome and the
one worth recording: the numbers were right, and they are now right with a
query beside them.

Two of those numbers set up the whole project. The IP₃-binding core signature
PF08709 was carried by 12,338 proteins, the RIH domain PF01365 by 13,177, the
RIH-associated domain PF08454 by 12,062, the MIR domain PF02815 by 23,453,
whose excess over the others is the O-mannosyltransferases that share it, and
the generic pore PF00520 by 206,115, which is why the pore is not a family
signature.

The taxonomic distribution carried an anomaly. Forty Viridiplantae proteins
and forty-one fungal proteins carried PF08709, the IP₃-binding core, while
*Arabidopsis thaliana* and *Saccharomyces cerevisiae* had none. Eighty-one
records in kingdoms whose model organisms have no IP₃ receptor is either a
real distribution or an artefact, and no amount of reading resolves which
({fig:review_taxonomic_range}). Chapter 5 chases all of them individually.

![](figures/review_taxonomic_range.png)

**{fig:review_taxonomic_range}.** What the databases said the family's range
was at the start of this project, shown as signature counts by taxon from a
live InterPro query. The plant and fungal columns are the anomaly. They are
also a good illustration of why a count of records is a weak instrument for a
range question: a record enters this figure by carrying a domain annotation,
and nothing in the count says whether the protein is real, whether it is in
that organism's genome, or whether it is an IP₃ receptor at all.

## 2.5 The sister-family hazard was measured in the query that was supposed to be clean

The background document cited a UniProt query on taxonomy 7955 and Pfam
PF08709 as evidence that zebrafish carries four IP₃ receptors. Re-running that
query returns 109 protein records across ten gene symbols, and only four of
those genes are IP₃ receptors. Fifty-three of the 109 records, or **49 %**,
are ryanodine receptors.

This is the single most useful number the audit produced. It is a measured
floor on the contamination that any signature-driven enumeration of this
family inherits, obtained in the exact query that had been quoted as a clean
result. It is why every later chapter separates the two families by positive
evidence rather than by filtering. It also carries a second finding in
passing: one of the six ryanodine-sized genes in that result, a 4,900-residue
locus, carries no gene name at all. That is the first appearance of the
unnamed-locus problem that Chapter 13 eventually audits across 503 genomes.

The length band happens to separate the two families cleanly in this
particular query. That is a fact about zebrafish annotation quality rather
than a rule, and treating it as one would import the quality of an annotation
into the definition of a gene family.

![](figures/review_family_separation.png)

**{fig:review_family_separation}.** Separating the families by evidence rather
than by name, using the control panel Chapter 3 builds. Panel **a** shows
every control scored against a labelled IP₃ receptor bait and a labelled
ryanodine receptor bait, and the two families fall on opposite sides of the
diagonal with nothing between them. Panel **b** shows the same data as a
margin, with the band inside which this project declines to call either way
drawn rather than described. The one positive inside that band is
*Dictyostelium* iplA, a characterised receptor that the family's own defining
signature does not find. Panel **c** shows the zebrafish query and the 53
records in it that are the sister family.

## 2.6 Two sequence measurements become a coordinate system for later chapters

Two further measurements were made for the literature review, and the rest of
this thesis uses them as a coordinate system rather than as a result.

The first is a per-column conservation profile over a family alignment, in
human ITPR1 numbering ({fig:review_conservation}). The second is a set of
alignment windows anchored on sites measured in the structure, meaning the
three stretches that contact the bound IP₃, the selectivity filter and the
gate, rather than on a residue list taken from a paper
({fig:review_alignment_windows}).

![](figures/review_conservation.png)

**{fig:review_conservation}.** Per-column conservation across a 25-sequence
family alignment, mapped onto human ITPR1 numbering, with the domain
architecture beneath it. Chapter 11 rebuilds this at a very different depth,
using 249 to 265 orthologues of each individual paralogue rather than 25
sequences of the whole family, and the two agree about where the peaks are.
That agreement is the useful thing to know about a profile computed from 25
sequences.

![](figures/review_alignment_windows.png)

**{fig:review_alignment_windows}.** Where the two families agree, and where
they stop agreeing. At the pore they are interchangeable, showing GGGVGD
against GGGIGD in the filter and a gate that differs by conservative
substitution. At the ligand site they are not: the arginines and lysines that
grip the trisphosphate are absent from all three ryanodine receptors, and one
window carries a three-residue deletion. The two families share a
superfamily and a binding site that only one of them uses. Chapter 12 turns
that observation into a measurement.

Measuring the structure for these figures produced three things the text did
not have. The channel's two constrictions were recovered blind from the
coordinates and landed on the GGGVGD filter motif and on Phe2513 and Ile2517,
without the geometry being told where the filter was. The distance from the
ligand to the gate resolved into 103 Å along the axis and 120 Å through space,
two numbers that are usually quoted as one. And *Dictyostelium* iplA, a
characterised IP₃-gated channel, carries none of PF08709, PF02815 or PF00520,
making it a real receptor that the family's defining signature does not find.
That last observation is why Chapter 3's profile seeds include it
deliberately.

![](figures/review_regulation.png)

**{fig:review_regulation}.** The curated regulator map, showing what binds the
receptor, where, and with what effect. It is included here because it is the
part of the biology this thesis does not measure. Nothing in a genome-scale
census can say whether a protein interaction is conserved, and a reader should
be able to see the size of what is being left out.

## 2.7 A fault in one public interface, and why it was worth diagnosing

The audit ended by running the search application against all three public
interfaces. Two of four runs returned all three sources, and the two failures
were more useful than the successes.

Ensembl returned nothing at all for the human panel. This was not an outage
and not a slow network. The endpoint the client used to resolve gene symbols,
`/xrefs/symbol/{species}/{symbol}`, **stalls indefinitely for `homo_sapiens`**,
and it does so for *BRCA2* as well as for *ITPR1*, so the fault is not
gene-specific. The same endpoint answers in 0.6 s for `danio_rerio`. That is
why the zebrafish run got three sources and the human run got none, and it is
why a retry budget could never have fixed it, because there is nothing to
retry against. The client now resolves through the database's other symbol
endpoint [R160], which works, and falls back to the old path only on a clean
404, never after a transport error, which would walk straight back into the
stall. Human ITPR1 went from 0 variants to 24.

Fixing that fault exposed a second one underneath it. Ensembl's latency is
unstable: the same 451-byte call measured 0.61 s, 7.61 s and 13.91 s inside
one session, and the call that carries the transcript and exon payload costs
about 12 s every time. A panel of eight species by three genes is 24
sequential pairs, which is somewhere between five and thirty-eight minutes.
The headless run budget was raised accordingly. The instruction that came out
of it is worth restating, because it applies to every latency figure in this
thesis: do not treat a single measurement of a public interface as a rate, and
treat the spread as the design constraint.

## 2.8 What Chapter 3 inherits from this audit

Chapter 3 inherits four things. A background document whose every statement is
tagged by how far it can be trusted, with four claims corrected. A set of
database numbers with queries attached, re-derived and reproducing exactly. A
measured floor on sister-family contamination of 49 % in the query that had
been quoted as clean, which turns the separation of the two receptor families
from a caveat into a design requirement. And a gene-architecture question,
raised by the one background claim that was simply false, which becomes
Chapter 8.
