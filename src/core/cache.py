"""Simple JSON-file disk cache for API responses.

One file per `(source, cache_key)`. Files are small JSON blobs containing the
normalized `ProteinVariant` list and a wall-clock timestamp. The cache is
inspectable by hand (`cat cache/ncbi__itpr1__.json`) which is useful when
debugging which database returned what.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict
from pathlib import Path
from typing import Optional

from .models import ProteinVariant


class DiskCache:
    DEFAULT_TTL_S = 60 * 60 * 24 * 7  # 1 week

    def __init__(self, root: str | os.PathLike, ttl_s: int = DEFAULT_TTL_S) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.ttl_s = ttl_s

    def _safe(self, name: str) -> str:
        return "".join(c if c.isalnum() or c in "._-+:" else "_" for c in name)

    def _path(self, key: str) -> Path:
        return self.root / f"{self._safe(key)}.json"

    def get(self, key: str) -> Optional[list[ProteinVariant]]:
        p = self._path(key)
        if not p.exists():
            return None
        try:
            payload = json.loads(p.read_text())
        except (json.JSONDecodeError, OSError):
            return None
        if time.time() - payload.get("timestamp", 0) > self.ttl_s:
            return None
        return [ProteinVariant(**v) for v in payload.get("variants", [])]

    def put(self, key: str, variants: list[ProteinVariant]) -> None:
        p = self._path(key)
        payload = {
            "timestamp": time.time(),
            "variants": [asdict(v) for v in variants],
        }
        try:
            p.write_text(json.dumps(payload, indent=2, default=str))
        except OSError:
            pass  # cache is best-effort

    def clear(self) -> int:
        n = 0
        for f in self.root.glob("*.json"):
            try:
                f.unlink()
                n += 1
            except OSError:
                pass
        return n
