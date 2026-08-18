## 7. Evolution

### 7.1 One superfamily, not two families

The IP<sub>3</sub> receptors and the ryanodine receptors are a single
structural superfamily. The relationship was visible in the first
IP<sub>3</sub>R sequence, whose reporting title stated it [R04], and it has
been confirmed at every subsequent level of resolution: the receptors share the
MIR, RIH, RIH-associated and Ion_trans domains; their N-terminal regions are
structurally and functionally conserved to the point that domains can be
compared directly [R58]; and the 2015 near-atomic structures of both families
[R22, R26, R27, R28] show the same overall organisation — a vast cytosolic
solenoid cap transducing ligand binding to a C-terminal pore module — at
different scales (RyR subunits are ~4,900–5,000 aa, nearly twice the size).

Three consequences run through the rest of this article, and through any
computational study of the family.

1. **Any similarity search for one family returns the other.** This is not a
   nuisance to be filtered away; it is a statement about descent, and the
   separation must be made by positive evidence — best-profile assignment or a
   labelled-bait margin — rather than by assumption. The scale of the problem is
   easy to underestimate: a single query for the IP<sub>3</sub>-binding-core
   domain in zebrafish returns 109 protein records across 10 gene symbols, of
   which 53 (49%) are ryanodine receptors, including one unnamed
   4,900-residue locus.
2. **RyR is the natural outgroup** for rooting an IP<sub>3</sub>R phylogeny —
   better conditioned than any invertebrate IP<sub>3</sub>R, because it is a
   genuine sister clade rather than a long branch within the ingroup.
3. **Size is a filter, not evidence.** The two families separate cleanly by
   length in well-annotated genomes, but that is a property of annotation
   quality in those genomes, not a phylogenetic argument.

### 7.2 Origins: the split predates animals

Comparative genomics places the machinery early. Analysis of genomes flanking
the animal–fungal divergence — including the apusozoan *Thecamonas trahens*,
from the putative unicellular sister group to Opisthokonta — finds many
components of animal and fungal Ca<sup>2+</sup> signalling already present in
the common ancestor, together with lineage-specific expansions of
Ca<sup>2+</sup> channels in the unicellular ancestors of animals and in basal
fungi [R35]. Homologues of both intracellular release-channel families are
present in parasitic protists [R36]. The broader comparative literature
places Ca<sup>2+</sup> signalling as an ancient and elaborate system rather
than a metazoan invention [R107, R108], with acidic-store channels following
their own evolutionary trajectory [R109]; the RyR side is reviewed separately
[R37].

### 7.3 The invertebrate single-gene state

Most invertebrates carry a single *itpr*. *Drosophila* is the best-developed
model: disruption of the gene affects larval metamorphosis and
ecdysone release [R110], genetic dissection assigns a vital requirement to aminergic neurons
[R111], and hypomorphs lose flight and the associated neuronal rhythmicity
[R112]. The single *Drosophila* receptor has been characterised biophysically
and behaves recognisably like its vertebrate counterparts [R113]. The
*Xenopus* receptor, cloned early, similarly established conservation of
structure and function across vertebrates [R106].

A single-gene invertebrate state carrying out the functions distributed across
three vertebrate paralogues is the strongest available argument that
vertebrate paralogue specialisation is subfunctionalisation of an ancestral
repertoire rather than the acquisition of new capabilities.

### 7.4 The vertebrate expansion — and what is not established

Vertebrates carry three paralogues; so, independently or not, do the ryanodine
receptors. Teleosts add a further layer: zebrafish carries four
IP<sub>3</sub> receptor genes — *itpr1a*, *itpr1b*, *itpr2* and *itpr3* — with
the *itpr1a*/*itpr1b* pair bearing the signature of the teleost-specific
whole-genome duplication (longest isoform per gene 2,635–2,819 aa, re-derived
here).

Three questions about this history are **not settled by the existing
literature**, and they are stated as open rather than glossed:

- Whether *ITPR1/2/3* are ohnologues from the two rounds of vertebrate
  whole-genome duplication has not been demonstrated with synteny-backed,
  phylogeny-tested evidence.
- Whether the IP<sub>3</sub>R and RyR triplications were **independent** events
  is frequently asserted and, as far as this review's search could establish,
  nowhere demonstrated. It is an attractive parallel, not a result.
- Which two of the three vertebrate paralogues are sisters — the rooted
  topology of the family — is not fixed by any published, support-annotated
  maximum-likelihood analysis with an RyR outgroup.

### 7.5 The taxonomic range problem

Textbook accounts hold that land plants and fungi lack IP<sub>3</sub>
receptors, and the model organisms support this: *Arabidopsis thaliana* and
*Saccharomyces cerevisiae* carry no protein annotated with the
IP<sub>3</sub>-binding core. Yet the InterPro protein set for that same
signature contains 40 Viridiplantae and 41 Fungi entries (re-derived
2026-08-18, alongside 12,149 Metazoa of 12,338 total). The lineage-specific
expansion of Ca<sup>2+</sup> channels reported in basal fungi [R35] means a
real-gene explanation cannot be dismissed a priori — but neither can
mis-annotation or contamination. **Whether the famous absences are facts about
genomes or facts about proteome databases is, at present, unresolved**, and it
cannot be settled from database counts alone.
