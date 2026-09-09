"""S25 — the references added for the thesis, declared for audit.

The manuscript cites 29 references and the review 137, all of them biology.
The thesis has to cite the *tools*: a chapter that explains why L-INS-i was
run single-threaded, or why a branch-site likelihood ratio is tested against a
50:50 mixture, is citing a method and must say whose.

**Nothing bibliographic is typed here.** Each row declares only an identifier
(a PubMed ID, or a DOI where the work predates PubMed or is not indexed), a
distinctive phrase that must appear in the title the identifier resolves to,
and what the thesis rests on the reference for. `s25_refs.py` fetches the
metadata live, checks the phrase, and writes the reference row from the
fetched record — so a mistyped identifier fails the audit instead of quietly
putting a different paper in the bibliography, and no author, title, year or
journal can be wrong because none of them is entered by hand (D13 applied to
a bibliography).

`kind` is the audit class, following S0's vocabulary:
  tool    software this project actually ran; the version is in
          `results/toolchain_manifest.txt`
  method  an algorithm, model or statistical procedure implemented or applied
  db      a database queried, cited at the release the project used
  lit     a substantive biological or evolutionary claim
"""

from __future__ import annotations

#: ref_id -> dict(pmid | doi, expect, kind, used_for)
#: `expect` is matched case-insensitively against the fetched title after
#: collapsing whitespace and punctuation.
NEW_REFS: dict[str, dict] = {
    # ---------------------------------------------- alignment and trees
    "R138": dict(pmid="23329690", expect="MAFFT multiple sequence alignment",
                 kind="tool",
                 used_for="every alignment in the project; L-INS-i for the "
                          "representative and seed alignments"),
    "R139": dict(pmid="19505945", expect="trimAl", kind="tool",
                 used_for="automated trimming of the representative and "
                          "codon alignments"),
    "R140": dict(pmid="32011700", expect="IQ-TREE 2", kind="tool",
                 used_for="the maximum-likelihood phylogeny and the "
                          "constrained searches behind the AU test"),
    "R141": dict(pmid="28481363", expect="ModelFinder", kind="method",
                 used_for="model selection; the exhaustive scan this project "
                          "measured and abandoned"),
    "R142": dict(pmid="29077904", expect="UFBoot2", kind="method",
                 used_for="ultrafast bootstrap support and the --bnni guard "
                          "against model violation"),
    "R143": dict(pmid="20525638", expect="New algorithms and methods to "
                                         "estimate maximum-likelihood",
                 kind="method",
                 used_for="the SH-aLRT branch test reported beside every "
                          "bootstrap value"),
    "R144": dict(pmid="12079646", expect="approximately unbiased test",
                 kind="method",
                 used_for="the AU test that decides which two paralogues are "
                          "sisters"),
    # ------------------------------------------------------- selection
    "R145": dict(pmid="17483113", expect="PAML 4", kind="tool",
                 used_for="every codon-model fit: one-ratio, two-ratio, "
                          "branch-site and site models"),
    "R146": dict(pmid="16107592", expect="improved branch-site likelihood",
                 kind="method",
                 used_for="branch-site model A on each paralogue stem, and "
                          "its 50:50 mixture null"),
    "R147": dict(pmid="15689528", expect="Bayes empirical Bayes",
                 kind="method",
                 used_for="the per-site posteriors reported beside every "
                          "significant selection test"),
    "R148": dict(pmid="16845082", expect="PAL2NAL", kind="tool",
                 used_for="the protein-guided codon alignment, cross-checked "
                          "against an independent in-house mapping"),
    "R149": dict(pmid="31504749", expect="HyPhy 2.5", kind="tool",
                 used_for="RELAX and FEL, run locally so the version is "
                          "pinned in the toolchain manifest"),
    "R150": dict(pmid="25540451", expect="RELAX", kind="method",
                 used_for="the test for relaxed or intensified selection "
                          "between paralogue clades"),
    "R151": dict(pmid="15703242", expect="Not so different after all",
                 kind="method",
                 used_for="the per-site rates carried onto the receptor's "
                          "own coordinates"),
    "R152": dict(pmid="3444411", expect="Simple methods for estimating the "
                                        "numbers of synonymous",
                 kind="method",
                 used_for="the pure-Python dN/dS in the app's selection "
                          "tool"),
    # -------------------------------------------- profiles and archives
    "R153": dict(pmid="22039361", expect="Accelerated Profile HMM Searches",
                 kind="tool",
                 used_for="hmmbuild, hmmsearch and the two family profiles "
                          "every census edition is scored against"),
    "R154": dict(pmid="20718988", expect="iterative HMM search procedure",
                 kind="method",
                 used_for="the jackhmmer convergence runs and the kill "
                          "criterion evaluated on them"),
    "R155": dict(pmid="33125078", expect="Pfam", kind="db",
                 used_for="the four family signatures that define the family "
                          "in one place"),
    "R156": dict(pmid="33156333", expect="InterPro protein families and "
                                         "domains database",
                 kind="db",
                 used_for="the enumerated search space and every per-protein "
                          "domain architecture"),
    "R157": dict(pmid="36408920", expect="UniProt", kind="db",
                 used_for="the reference proteomes, the records and the "
                          "taxonomy every census row is filed under"),
    "R158": dict(pmid="20003500", expect="BLAST+", kind="tool",
                 used_for="blastp attribution, tblastn rescue and every "
                          "reciprocal best-hit check"),
    "R159": dict(pmid="38969627", expect="NCBI Datasets", kind="db",
                 used_for="the declared assembly manifests and every genome "
                          "and annotation fetched"),
    "R160": dict(pmid="37953337", expect="Ensembl", kind="db",
                 used_for="gene structure, the paralogy map and the "
                          "endpoint probe behind the S0 client change"),
    "R161": dict(pmid="27141089", expect="Ensembl comparative genomics",
                 kind="db",
                 used_for="the gene trees and duplication nodes that date "
                          "each paralogy link"),
    "R162": dict(pmid="25897122", expect="BioMart", kind="db",
                 used_for="the pinned-archive paralogy map behind the 2R "
                          "test"),
    "R163": dict(pmid="21062823", expect="Sequence Read Archive", kind="db",
                 used_for="every RNA-seq run streamed for the expression "
                          "evidence"),
    # ------------------------------------------------ genome alignment
    "R164": dict(pmid="36648328", expect="miniprot", kind="tool",
                 used_for="protein-to-genome alignment for all 503 swept "
                          "genomes"),
    "R165": dict(pmid="31375807", expect="graph-based genome alignment",
                 kind="tool",
                 used_for="read mapping for the junction-level expression "
                          "evidence"),
    # ------------------------------------------------------ structures
    "R166": dict(pmid="15849316", expect="TM-align", kind="tool",
                 used_for="every structural comparison, and the two "
                          "published bars the calibration figure draws"),
    "R167": dict(pmid="34265844", expect="Highly accurate protein structure "
                                         "prediction",
                 kind="method",
                 used_for="the predicted models the structural panel is "
                          "mostly made of"),
    "R168": dict(pmid="34791371", expect="AlphaFold Protein Structure "
                                         "Database",
                 kind="db",
                 used_for="the coverage probe that found the database does "
                          "not hold this family"),
    "R169": dict(pmid="37156916", expect="Fast and accurate protein "
                                         "structure search",
                 kind="tool",
                 used_for="the optional structure-based homology sweep"),
    # -------------------------------------- statistics and reconciliation
    "R170": dict(doi="10.1111/j.2517-6161.1995.tb02031.x",
                 expect="Controlling the false discovery rate",
                 kind="method",
                 used_for="every family of tests in the project is corrected "
                          "together rather than reported one at a time"),
    "R171": dict(pmid="7966282", expect="Position-based sequence weights",
                 kind="method",
                 used_for="sequence weighting before any column statistic, "
                          "without which every teleost-specific residue "
                          "reads as conserved"),
    "R172": dict(doi="10.2307/2412519",
                 expect="Fitting the gene lineage into its species lineage",
                 kind="method",
                 used_for="the LCA reconciliation mapping"),
    "R173": dict(pmid="11590098", expect="infer gene duplication and "
                                         "speciation events",
                 kind="method",
                 used_for="the loss count attached to each reconciled node"),
    "R174": dict(pmid="18808330", expect="non-binary species trees",
                 kind="method",
                 used_for="the duplication test that stays valid on a "
                          "polytomy, which the familiar binary rule does "
                          "not"),
    "R175": dict(doi="10.2307/2412867",
                 expect="Phylogenetic analysis under Dollo's Law",
                 kind="method",
                 used_for="the primary loss count"),
    "R176": dict(pmid="12116640", expect="likelihood approach to estimating "
                                         "phylogeny from discrete "
                                         "morphological character",
                 kind="method",
                 used_for="the Mk fits run beside the parsimony count"),
    "R177": dict(pmid="15405679", expect="Index for rating diagnostic tests",
                 kind="method",
                 used_for="the separation statistic every calibration in the "
                          "project reports beside its threshold"),
    "R178": dict(doi="10.1109/18.61115",
                 expect="Divergence measures based on the Shannon entropy",
                 kind="method",
                 used_for="the per-column conservation metric, reported "
                          "beside a composition-free control"),
    # ------------------------------------- duplication and genome evolution
    "R179": dict(doi="10.1007/978-3-642-86659-3",
                 expect="Evolution by Gene Duplication", kind="lit",
                 used_for="the two-round whole-genome duplication hypothesis "
                          "the paralogue origin is tested against"),
    "R180": dict(pmid="16128622", expect="Two rounds of whole genome "
                                         "duplication",
                 kind="lit",
                 used_for="the 2R account this project's paralogon test is "
                          "read against"),
    "R181": dict(pmid="15496914", expect="Tetraodon nigroviridis",
                 kind="lit",
                 used_for="the teleost-specific genome duplication behind "
                          "the two-copy cells"),
    "R182": dict(pmid="32313176", expect="Deeply conserved synteny",
                 kind="lit",
                 used_for="the current account of what 2R left in the "
                          "vertebrate genome, against which a weak paralogon "
                          "signal has to be read"),
    "R183": dict(pmid="17652425", expect="Reconstruction of the vertebrate "
                                         "ancestral genome",
                 kind="lit",
                 used_for="the ancestral linkage groups the ITPR "
                          "neighbourhoods are placed on"),
    "R184": dict(pmid="9831563", expect="Zebrafish hox clusters",
                 kind="lit",
                 used_for="the first genomic evidence for the teleost "
                          "duplication"),
    "R185": dict(pmid="27762356", expect="Genome evolution in the "
                                         "allotetraploid frog",
                 kind="lit",
                 used_for="the extra-whole-genome-duplication lineages "
                          "carried as copy-number controls"),
    "R186": dict(pmid="29358652", expect="sea lamprey germline genome",
                 kind="lit",
                 used_for="the cyclostome assemblies whose three ITPR loci "
                          "the bait panel cannot label"),
    "R187": dict(pmid="30464347", expect="Amphioxus functional genomics", kind="lit",
                 used_for="the pre-duplication chordate outgroup"),
    # ------------------------------------------------ time calibrations
    "R188": dict(pmid="35932227", expect="TimeTree 5", kind="db",
                 used_for="the median divergence times behind every node age "
                          "in the species tree"),
    "R189": dict(doi="10.26879/424",
                 expect="Constraints on the timescale of animal evolutionary "
                        "history",
                 kind="lit",
                 used_for="the fossil-calibrated bounds quoted beside each "
                          "TimeTree median"),
    "R190": dict(pmid="28890940", expect="jawed vertebrate timetree",
                 kind="lit",
                 used_for="the alternative deep vertebrate ages that make "
                          "the root calibration the widest disagreement in "
                          "the tree"),
    "R191": dict(pmid="28683774", expect="Phylogenetic classification of "
                                         "bony fishes",
                 kind="lit",
                 used_for="the teleost topology and the membership of "
                          "Ovalentaria"),
    "R192": dict(pmid="26444237", expect="comprehensive phylogeny of birds",
                 kind="lit",
                 used_for="the avian backbone, where the sweep's contiguity "
                          "problem is worst"),
    # ------------------------------------------------- loss and annotation
    "R193": dict(pmid="27087500", expect="Evolution by gene loss", kind="lit",
                 used_for="the expectation that a gene family of this age "
                          "should have lost copies somewhere"),
    "R194": dict(pmid="31097009", expect="Next-generation genome annotation",
                 kind="lit",
                 used_for="the state of automated annotation this project's "
                          "audit measures against"),
    "R195": dict(pmid="33911273", expect="Towards complete and error-free "
                                         "genome assemblies",
                 kind="lit",
                 used_for="the assembly-quality standard behind the "
                          "contiguity bar every absence claim passes"),
}
