## Introduction

The IP₃ receptor releases calcium from the endoplasmic reticulum when it binds
inositol 1,4,5-trisphosphate, the second messenger made when phospholipase C
cleaves a membrane phosphoinositide [R01, R07]. It is a tetrameric channel of
about 2,700 residues per subunit, and in vertebrates it exists as three
paralogues, ITPR1, ITPR2 and ITPR3 [R50, R134]. Its closest relatives are the
ryanodine receptors, calcium-release channels of about 5,000 residues that
share its pore and much of its cytoplasmic architecture [R04, R58].

Where the receptor occurs outside the animals has been argued from small
samples. Comparative surveys of the calcium-signalling toolkit found
receptor-like genes in some protists, green algae and early-diverging fungi,
and none in the land plants or the yeasts [R35, R107, R108, R37]. The land
plants are the case that matters most, because plant cells release calcium in
response to IP₃, no plant receptor gene has been identified, and the standing
review of the question asks in its title whether the receptor is real [R196].
If the family is absent there, whatever performs that physiology is built from
other channels [R197, R198].

Two problems make that question harder than it looks. The first is the
ryanodine receptor. Every Pfam signature diagnostic of the IP₃ receptor,
including the IP₃-binding core (PF08709), the RyR and IP₃R homology domain
(PF01365), the RIH-associated domain (PF08454) and the MIR domain (PF02815),
also occurs in ryanodine receptors. A search for one family therefore returns
the other, and a range built from domain annotation would be a range of the
two families mixed. Length separates them in well-annotated proteomes, but
length is exactly what a fragmentary or unannotated record does not carry.

The second problem is that a database returns records and not genes. "The
family occurs in N proteomes" is a statement about what gene callers found,
and "the family is absent from land plants" is a claim about DNA that a
proteome cannot make on its own. An absence becomes a measurement only when
the same search, in the same genome, demonstrably recovers a comparable gene.
A positive control is what licenses a negative claim, and range surveys
rarely carry one.

We therefore asked one question, in which eukaryotic lineages an IP₃ receptor
gene exists, and answered it in four layers, each checked against the one
before. First, every record is assigned to the IP₃ receptors or the ryanodine
receptors by a positive test, meaning evidence that it is one rather than
evidence that it is not the other, and the tests were benchmarked on decoys
before any search ran. Second, the family's profiles were swept over 7,691
reference proteomes spanning Eukaryota, archaea and a genus-stratified
bacterial sample, with each negative claim made at two stated sensitivities.
Third, every plant and fungal record was chased individually for
contamination. Fourth, the absences were taken to 194 genome assemblies, each
carrying a positive control selected by measurement for its own clade. We
report the family's range with that scope attached to every claim, and the
fold of the structures that exist as an independent check on the family call.
