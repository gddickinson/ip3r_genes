"""Core data models, search orchestration, and caching.

Contents:
    models.py  — ProteinVariant, GeneRecord, SearchQuery, SearchResult dataclasses
    cache.py   — JSON-file disk cache keyed by (database, query)
    search.py  — SearchOrchestrator: fans out queries across enabled DBs in threads
"""

from .models import ProteinVariant, GeneRecord, SearchQuery, SearchResult, SearchStatus
from .cache import DiskCache
from .search import SearchOrchestrator

__all__ = [
    "ProteinVariant",
    "GeneRecord",
    "SearchQuery",
    "SearchResult",
    "SearchStatus",
    "DiskCache",
    "SearchOrchestrator",
]
