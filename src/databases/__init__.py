"""Database clients.

Each client subclasses `DatabaseClient` and returns a list of `ProteinVariant`
for a `SearchQuery`. The clients normalize source-specific schemas into the
common `ProteinVariant` shape so the GUI doesn't need to know which DB a row
came from.

Sequence-based clients (BLAST-style discovery):
    ncbi.py       — NCBI Entrez (Biopython Bio.Entrez)
    ensembl.py    — Ensembl REST (rest.ensembl.org)
    uniprot.py    — UniProt REST (rest.uniprot.org)

Structure-based clients (remote-homology / twilight-zone discovery):
    alphafold.py  — AlphaFold Protein Structure DB (predicted 3D models)
    foldseek.py   — Foldseek structural homology search (opt-in, async)
"""

from .base import DatabaseClient
from .ncbi import NCBIClient
from .ensembl import EnsemblClient
from .uniprot import UniProtClient
from .alphafold import AlphaFoldClient
from .foldseek import FoldseekClient
from .compara import ComparaClient
from .blast import BlastClient

# Sources the orchestrator can instantiate by name. Bait-based clients
# (BLAST, Foldseek) are constructed specially — they need bait sequences
# harvested from the sequence clients' results first.
AVAILABLE_CLIENTS: dict[str, type[DatabaseClient]] = {
    "NCBI": NCBIClient,
    "Ensembl": EnsemblClient,
    "UniProt": UniProtClient,
    "Compara": ComparaClient,
    "AlphaFold": AlphaFoldClient,
    "Foldseek": FoldseekClient,
    "BLAST": BlastClient,
}

SEQUENCE_SOURCES = ["NCBI", "Ensembl", "UniProt", "Compara"]
STRUCTURE_SOURCES = ["AlphaFold", "Foldseek"]
BAIT_SOURCES = ["BLAST", "Foldseek"]   # need a bait sequence → run after sequence DBs

__all__ = [
    "DatabaseClient",
    "NCBIClient",
    "EnsemblClient",
    "UniProtClient",
    "ComparaClient",
    "AlphaFoldClient",
    "FoldseekClient",
    "BlastClient",
    "AVAILABLE_CLIENTS",
    "SEQUENCE_SOURCES",
    "STRUCTURE_SOURCES",
    "BAIT_SOURCES",
]
