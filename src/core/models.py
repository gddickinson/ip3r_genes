"""Data models shared across database clients and the GUI.

A `ProteinVariant` is the unit of comparison — one protein-coding record from
one database. Multiple variants for the same biological protein (e.g., the
NCBI RefSeq record and the matching UniProt entry) are kept as separate
variants so the user can see discrepancies in length, naming, etc.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Optional


class SearchStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    OK = "ok"
    EMPTY = "empty"
    ERROR = "error"


@dataclass
class ProteinVariant:
    """One protein record from one source database.

    Fields are deliberately permissive (most Optional) because the three
    databases expose different subsets. The `raw` dict carries the original
    payload for the detail pane.
    """
    source: str               # "NCBI" | "Ensembl" | "UniProt"
    accession: str            # e.g. "XP_068079622.1" / "ENSP000..." / "Q92508"
    gene_symbol: str          # e.g. "ITPR1"
    species: str = ""         # scientific name, e.g. "Homo sapiens"
    taxon_id: Optional[int] = None
    length_aa: Optional[int] = None
    description: str = ""
    transcript_id: str = ""   # mRNA / transcript accession if known
    gene_id: str = ""         # source-specific gene id
    sequence: str = ""        # protein sequence, populated lazily
    url: str = ""             # human-readable detail URL
    raw: dict[str, Any] = field(default_factory=dict)

    def key(self) -> str:
        """Stable id used for dedup within the results table."""
        return f"{self.source}:{self.accession}"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # raw can be large + non-trivially-serializable; let exporter decide
        d.pop("raw", None)
        return d


@dataclass
class GeneRecord:
    """A higher-level grouping of variants by (source, gene)."""
    source: str
    gene_symbol: str
    gene_id: str
    species: str = ""
    taxon_id: Optional[int] = None
    description: str = ""
    variants: list[ProteinVariant] = field(default_factory=list)


@dataclass
class SearchQuery:
    """A single user query."""
    gene_symbols: list[str]        # one or more gene names
    species: str = ""              # blank = any
    taxon_id: Optional[int] = None # if known, more precise than name
    sources: list[str] = field(default_factory=lambda: ["NCBI", "Ensembl", "UniProt"])
    max_results_per_source: int = 50
    include_sequence: bool = False  # fetch sequences (slower)

    def cache_key(self, source: str) -> str:
        genes = "+".join(sorted(s.lower() for s in self.gene_symbols))
        sp = (self.species or "").lower().replace(" ", "_")
        return f"{source.lower()}::{genes}::{sp}::tax{self.taxon_id or 0}::n{self.max_results_per_source}::seq{int(self.include_sequence)}"


@dataclass
class SearchResult:
    """Per-source outcome handed back from worker threads to the GUI."""
    source: str
    query: SearchQuery
    status: SearchStatus
    variants: list[ProteinVariant] = field(default_factory=list)
    message: str = ""              # human-readable error or note
    elapsed_s: float = 0.0
