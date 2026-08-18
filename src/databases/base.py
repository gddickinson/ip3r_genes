"""Abstract base class for database clients."""

from __future__ import annotations

from abc import ABC, abstractmethod

from ..core.models import ProteinVariant, SearchQuery


class DatabaseClient(ABC):
    """Common contract for source-specific clients.

    Subclasses must implement `search`. They should:
      * raise on hard failure (orchestrator catches and reports)
      * return [] for "no hits" (not raise)
      * fill in as many ProteinVariant fields as the source supports
    """

    name: str = "Unknown"

    def __init__(self, email: str = "", api_key: str = "", timeout_s: int = 30) -> None:
        self.email = email
        self.api_key = api_key
        self.timeout_s = timeout_s

    @abstractmethod
    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        """Return a flat list of variants matching the query."""
        raise NotImplementedError
