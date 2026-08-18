"""The gene family this project is built around — one definition, one place.

Everything family-specific in `src/` reads from here: the Pfam signatures a
candidate must carry, the size band a real family member falls in, the
reference proteins used as the phylogenetic comparison panel, and the
literature keyword. Porting the app to another family is an edit of this
file plus the presets, not a search-and-replace across the tree.

**This project: the IP3 receptor (ITPR) family.**

The inositol 1,4,5-trisphosphate receptor is the ER's ligand-gated Ca(2+)
release channel: a ~2,700-residue subunit that assembles as a homotetramer,
binds IP3 at its N-terminal beta-trefoil core and opens a six-TM pore at its
C-terminus. Vertebrates carry three paralogs (ITPR1/2/3); most invertebrates
carry one; the ryanodine receptors (RYR1/2/3, ~5,000 aa) are its sister
family inside the same Ca(2+)-release superfamily and share the pore module
but not the IP3-binding core.

Numbers here that come from a database (accessions, Pfam IDs, lengths) are
verified live by the S1 control benchmark — see `docs/session_briefs.md`.
"""

from __future__ import annotations

#: Display name of the family.
FAMILY_NAME = "IP3R"

#: Long name, used in report prose.
FAMILY_LONG_NAME = "inositol 1,4,5-trisphosphate receptor"

#: The vertebrate paralogs, in the order figures and tables use them.
PARALOGS = ["ITPR1", "ITPR2", "ITPR3"]

#: The sister family. Not part of the census, but every search that reaches
#: far enough will hit it, so it is named rather than left to surprise a
#: session: RyRs share the pore module and are the natural outgroup for
#: rooting the family tree.
SISTER_FAMILY = "RYR"
SISTER_PARALOGS = ["RYR1", "RYR2", "RYR3"]

#: Pfam signatures. A candidate must carry at least one of these; without
#: one, an outlier sequence is just an outlier, not an ITPR paralog.
#:
#:   PF08709  Ins145_P3_rec  IP3-binding core (beta-trefoil) — ITPR + RyR
#:   PF02815  MIR            N-terminal suppressor domain — ITPR + RyR + POMT
#:   PF01365  RYDR_ITPR      RyR and IP3R homology (RIH) domain
#:   PF08454  RIH_assoc      RIH-associated domain
#:   PF00520  Ion_trans      the six-TM pore module (shared with many channels)
#:
#: The first four are family-diagnostic in combination; PF00520 alone is not.
#: `FAMILY_PFAM_IDS` is what the discovery scorer tests, so it holds the
#: diagnostic set; `PORE_PFAM_IDS` is kept separate for the architecture work.
FAMILY_PFAM_IDS = ["PF08709", "PF02815", "PF01365", "PF08454"]
PORE_PFAM_IDS = ["PF00520"]

#: Pfam IDs whose protein lists define the enumerated search space (S2).
CENSUS_PFAM_IDS = ["PF08709", "PF01365", "PF08454", "PF02815"]

#: Size band for a plausible full-length family member, in residues. Human
#: ITPR1/2/3 are 2,758 / 2,701 / 2,671 aa and invertebrate Itpr is ~2,800;
#: the upper bound sits below the ~5,000 aa ryanodine receptors on purpose,
#: so a RyR cannot score as a novel ITPR paralog by size alone.
MIN_LENGTH_AA = 2_000
MAX_LENGTH_AA = 3_600

#: Reference proteins for the investigation panel: the three human paralogs.
#: (accession-label, UniProt accession)
REFERENCE_PANEL = [
    ("Hs_ITPR1_Q14643", "Q14643"),
    ("Hs_ITPR2_Q14571", "Q14571"),
    ("Hs_ITPR3_Q14573", "Q14573"),
]

#: Sister-family references — used to root trees and as negative controls.
SISTER_PANEL = [
    ("Hs_RYR1_P21817", "P21817"),
    ("Hs_RYR2_Q92736", "Q92736"),
    ("Hs_RYR3_Q15413", "Q15413"),
]

#: Literature search keyword (PubMed title/abstract).
LITERATURE_KEYWORD = "IP3 receptor"

#: Substrings that mark a hit as an already-named family member. Lower-cased
#: comparison; covers the vertebrate symbols and the invertebrate ones.
KNOWN_NAME_SUBSTRINGS = [
    "itpr1", "itpr2", "itpr3", "itpr",
    "ip3r", "insp3r", "itr-1", "iplA",
]

#: Gene symbols treated as "known paralogs" by the discovery scorer unless a
#: preset overrides them.
KNOWN_PARALOGS = ["ITPR1", "ITPR2", "ITPR3"]

#: The C-terminal fraction that holds the pore + gating machinery — the
#: segment worth folding when a full-length prediction is out of reach
#: (`src/analysis/esmfold.py`).
POR_SEGMENT_FROM_END_AA = 400
