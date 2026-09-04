"""S2 shared helpers — census v2 paths, signatures, TSV I/O, live panel.

The census's search space is defined by three Pfam signatures (the *seeds*).
Everything else in this module exists to turn that raw enumeration into a
table where every row carries a **positive** ITPR-or-RYR call with the
evidence that produced it.

Why these three seeds and not four: `PF02815` (MIR) is carried by the
O-mannosyltransferases as well as by both receptor families, so it widens
the search space without adding family evidence. It is kept as an
*annotation* column — a record's MIR status is recorded, it just never
decides what the record is. `src/utils/family.py:CENSUS_PFAM_IDS` is the
one place that list lives.

The RyR side of the call is the interesting half. Four signatures separate
the ryanodine receptors from the IP3 receptors *inside this search space*:

    PF02026  RyR       the ryanodine-receptor repeat
    PF06459  RR_TM4-6  the RyR transmembrane block
    PF21119  RyR_SPRY  (name resolved live; see `signature_names()`)
    PF00622  SPRY      B30.2/SPRY

`PF00622` is emphatically **not** RyR-specific in general — it sits in
~114,000 UniProt proteins, most of them TRIM ligases and butyrophilins.
It is diagnostic *conditionally*: among proteins that already carry a seed
signature, SPRY means RyR, because no IP3 receptor has one. That
conditional claim is not assumed here — `results/census_v2/rule_audit.tsv`
tests it against every labelled record in the census.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.family import (                    # noqa: E402
    CENSUS_PFAM_IDS, FAMILY_PFAM_IDS, PORE_PFAM_IDS,
)
OUT_DIR = ROOT / "results" / "census_v2"

#: The enumeration seeds — from the family definition, not redeclared here
#: (CLAUDE.md: the family lives in exactly one file).
SEED_PFAMS = list(CENSUS_PFAM_IDS)

#: Recorded on every row, never used to enumerate or to call: MIR (shared
#: with the O-mannosyltransferases) and the generic six-TM pore.
CONTEXT_PFAMS = [p for p in FAMILY_PFAM_IDS if p not in SEED_PFAMS] + \
    list(PORE_PFAM_IDS)

#: Positive RyR evidence *within the seeded search space*.
RYR_PFAMS = ["PF02026", "PF06459", "PF21119", "PF00622"]

#: Human-readable names, filled by `signature_names()` from InterPro and
#: cached, so no report has to hard-code what a Pfam ID means.
NAME_CACHE = OUT_DIR / "signature_names.tsv"

UNIPROT = "https://rest.uniprot.org/uniprotkb"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"

USER_AGENT = "ip3r_genes/S2 (census v2; https://github.com/gddickinson/ip3r_genes)"


# --------------------------------------------------------------- http
def fetch(url: str, timeout_s: int = 90, tries: int = 12,
          accept: str | None = None) -> tuple[int, bytes, dict]:
    """GET with jittered exponential backoff.

    Returns `(status, body, headers)`. Raises only when every attempt
    failed, because the InterPro API answers a steady fraction of requests
    with a 500 that a retry clears — see `results/census_v2/report.md`.
    """
    from src.databases.interpro import host_variants   # host failover

    last = ""
    for attempt in range(tries):
        for candidate in host_variants(url):
            req = urllib.request.Request(
                candidate, headers={"User-Agent": USER_AGENT})
            if accept:
                req.add_header("Accept", accept)
            try:
                with urllib.request.urlopen(req, timeout=timeout_s) as r:
                    return r.status, r.read(), dict(r.headers)
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}"
                if e.code in (400, 404):      # a dead cursor, not an outage
                    return e.code, b"", {}
            except Exception as e:            # noqa: BLE001 - network is broad
                last = f"{type(e).__name__}: {e}"
        # 2, 4, 8 … capped at 60 s, with a deterministic-ish jitter
        wait = min(60, 2 ** (attempt + 1)) + (attempt % 3) * 1.7
        time.sleep(wait)
    raise RuntimeError(f"GET failed after {tries} attempts ({last}): {url}")


def get_json(url: str, **kw) -> dict:
    status, body, _ = fetch(url, **kw)
    if status != 200:
        raise RuntimeError(f"HTTP {status} for {url}")
    return json.loads(body.decode("utf-8"))


# ---------------------------------------------------------- signatures
def signature_names(refresh: bool = False) -> dict[str, str]:
    """`{pfam_id: short name}` for every signature S2 touches.

    Cached to a committed table so the report can name a signature without
    a network call, and so a renamed Pfam entry is a visible diff.
    """
    if NAME_CACHE.exists() and not refresh:
        out = {}
        for line in NAME_CACHE.read_text().splitlines()[1:]:
            if line.strip():
                pid, name, _role = (line.split("\t") + ["", ""])[:3]
                out[pid] = name
        return out
    rows = []
    for role, ids in (("seed", SEED_PFAMS), ("context", CONTEXT_PFAMS),
                      ("ryr", RYR_PFAMS)):
        for pid in ids:
            meta = get_json(f"{INTERPRO}/entry/pfam/{pid}").get("metadata", {})
            name = meta.get("name") or {}
            short = name.get("short") if isinstance(name, dict) else None
            full = name.get("name") if isinstance(name, dict) else str(name)
            rows.append({"pfam_id": pid, "name": short or full or "",
                         "role": role, "long_name": full or ""})
    write_tsv(NAME_CACHE, rows, ["pfam_id", "name", "role", "long_name"])
    return {r["pfam_id"]: r["name"] for r in rows}


# ------------------------------------------------------------- tsv i/o
def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\t".join(cols or []) + "\n")
        return 0
    cols = cols or list(rows[0].keys())
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ")
                              .replace("\n", " ") for c in cols) + "\n")
    return len(rows)


def read_tsv(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    cols = lines[0].split("\t")
    return [dict(zip(cols, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


# --------------------------------------------------------- live panel
def live(task: str, steps: list[tuple[str, bool]]) -> None:
    """Write `results/session_live.json` for the dashboard's live panel."""
    p = ROOT / "results" / "session_live.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        "task": task,
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "steps": [{"label": lbl, "done": bool(done)} for lbl, done in steps],
    }, indent=1))
