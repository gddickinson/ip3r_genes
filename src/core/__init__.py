"""Core data models, search orchestration, and caching.

Contents:
    models.py  — ProteinVariant, GeneRecord, SearchQuery, SearchResult dataclasses
    cache.py   — JSON-file disk cache keyed by (database, query)
    search.py  — SearchOrchestrator: fans out queries across enabled DBs in threads

`SearchOrchestrator` is imported lazily (PEP 562). Importing it eagerly
pulls in `src.databases`, hence `Bio.Entrez`, hence biopython — which made
the session-protocol check (`python -m src.utils.data_root --require`)
die with `ModuleNotFoundError: No module named 'Bio'` in any environment
without biopython, reporting a missing package where the real answer was
"the drive is unplugged". The gate has to work in a bare interpreter.
"""

from .models import ProteinVariant, GeneRecord, SearchQuery, SearchResult, SearchStatus
from .cache import DiskCache

__all__ = [
    "ProteinVariant",
    "GeneRecord",
    "SearchQuery",
    "SearchResult",
    "SearchStatus",
    "DiskCache",
    "SearchOrchestrator",
]


def __getattr__(name: str):
    if name == "SearchOrchestrator":
        from .search import SearchOrchestrator
        return SearchOrchestrator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
