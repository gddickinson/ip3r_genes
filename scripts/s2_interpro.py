"""S2 step 1 — enumerate the seed Pfam signatures to exhaustion.

Three signatures define the census's search space (`s2_lib.SEED_PFAMS`).
For each one this walks InterPro's protein list to the API's own `count`,
archiving every raw page under the data root, then *re-parses the census
out of the archive* rather than out of memory. The parse is therefore
reproducible offline (`--parse-only`) and a re-run costs nothing.

Two failure modes are handled rather than hoped away:

* **Runs of HTTP 500.** The InterPro API answers a variable fraction of
  requests with a 500 that a later identical request serves fine; during
  this task's development it produced 13 consecutive failures over 70 s.
  `list_proteins_with_pfam` now backs off exponentially to a 60 s cap, so
  a page survives roughly eight minutes of outage before giving up.
* **Cursor expiry.** The `next` URLs are opaque cursors. A resumed one that
  has expired comes back 400/404, raising `CursorExpired`; the driver then
  restarts that signature from page 1 and records the restart, because a
  census that silently resumes onto a shifted result set is worse than one
  that costs an extra ten minutes.

Completeness is asserted, not assumed — and the API's own `count` turned
out not to be the right assertion. For PF08709 the walk returned 12,507
distinct accessions against an advertised `count` of 12,339 (stable across
re-queries), so a census that stopped when it reached the advertised total
would have dropped 168 proteins without noticing. The recorded completeness
test is therefore **the cursor chain running out** (`next: null`), with the
advertised count and UniProt's independent count for the same signature
both kept beside it as cross-checks.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.databases.interpro import (           # noqa: E402
    CursorExpired, list_proteins_with_pfam, parse_protein_row,
)
from src.utils.data_root import require_data_root   # noqa: E402
from scripts.s2_lib import fetch as http_fetch      # noqa: E402
from scripts.s2_lib import (                   # noqa: E402
    OUT_DIR, ROOT, SEED_PFAMS, live, write_tsv,
)

ARCHIVE = "raw_api/interpro"
COLS = ["accession", "gene", "name", "species", "taxon_id", "length",
        "reviewed", "in_alphafold", "seed_pfams", "n_seed_pfams"]


def sig_dir(pfam_id: str) -> Path:
    return require_data_root(ROOT) / ARCHIVE / pfam_id


def load_state(pfam_id: str) -> dict:
    p = sig_dir(pfam_id) / "state.json"
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            pass
    return {"pages": 0, "next": None, "count": None,
            "complete": False, "restarts": 0}


def save_state(pfam_id: str, state: dict) -> None:
    d = sig_dir(pfam_id)
    d.mkdir(parents=True, exist_ok=True)
    (d / "state.json").write_text(json.dumps(state, indent=1))


def enumerate_signature(pfam_id: str, page_size: int, log) -> dict:
    """Walk one signature to exhaustion, resuming an interrupted walk."""
    d = sig_dir(pfam_id)
    state = load_state(pfam_id)
    if state.get("complete"):
        log(f"  {pfam_id}: already complete "
            f"({state['pages']} pages, count={state['count']})")
        return state

    def on_page(page_no: int, data: dict, next_url: str | None) -> None:
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{pfam_id}_page_{page_no:04d}.json").write_text(json.dumps(data))
        state.update(pages=page_no, next=next_url,
                     count=state.get("count") or data.get("count"))
        save_state(pfam_id, state)
        if page_no % 10 == 0:
            log(f"  {pfam_id}: page {page_no}, {page_no * page_size} rows")

    for attempt in (1, 2):
        stats: dict = {}
        try:
            list_proteins_with_pfam(
                pfam_id, max_results=None, page_size=page_size,
                dump_dir=None, strict=True, max_retries=11,
                stats=stats, on_page=on_page,
                start_url=state.get("next"), start_page=state.get("pages", 0),
            )
        except CursorExpired as e:
            if attempt == 2:
                raise
            log(f"  {pfam_id}: {e} — restarting from page 1")
            for stale in d.glob(f"{pfam_id}_page_*.json"):
                stale.unlink()
            state = {"pages": 0, "next": None, "count": state.get("count"),
                     "complete": False, "restarts": state.get("restarts", 0) + 1}
            save_state(pfam_id, state)
            continue
        state["count"] = state.get("count") or stats.get("count")
        state["complete"] = True
        state["next"] = None          # the cursor chain ran out
        save_state(pfam_id, state)
        return state
    return state


def uniprot_count(pfam_id: str) -> int | None:
    """UniProt's independent count for the same signature.

    A second database's answer to the same membership question — the only
    external check available on whether the walk is complete.
    """
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           f"query=%28xref%3Apfam-{pfam_id}%29&size=1")
    try:
        _status, _body, headers = http_fetch(url, timeout_s=90, tries=4)
    except RuntimeError:
        return None
    for k, v in headers.items():
        if k.lower() == "x-total-results":
            try:
                return int(v)
            except ValueError:
                return None
    return None


def parse_archive(pfam_id: str) -> list[dict]:
    """Rebuild one signature's rows from its archived pages."""
    rows = []
    for page in sorted(sig_dir(pfam_id).glob(f"{pfam_id}_page_*.json")):
        data = json.loads(page.read_text())
        for entry in data.get("results", []):
            r = parse_protein_row(entry.get("metadata", {}), pfam_id)
            if r is not None:
                rows.append(r)
    return rows


def merge(per_sig: dict[str, list]) -> list[dict]:
    """One row per accession, with the seed signatures it carries."""
    merged: dict[str, dict] = {}
    for pfam_id in SEED_PFAMS:
        for r in per_sig.get(pfam_id, []):
            m = merged.setdefault(r.accession, {
                "accession": r.accession, "gene": r.gene, "name": r.name,
                "species": r.species, "taxon_id": r.taxon_id or "",
                "length": r.length, "reviewed": int(r.reviewed),
                "in_alphafold": int(r.in_alphafold), "_seeds": set(),
            })
            m["_seeds"].add(pfam_id)
            # A row can arrive from several signatures; keep the richest.
            if not m["gene"] and r.gene:
                m["gene"] = r.gene
    out = []
    for m in merged.values():
        seeds = sorted(m.pop("_seeds"))
        m["seed_pfams"] = ";".join(seeds)
        m["n_seed_pfams"] = len(seeds)
        out.append(m)
    return sorted(out, key=lambda r: r["accession"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--page-size", type=int, default=200)
    ap.add_argument("--parse-only", action="store_true",
                    help="rebuild the tables from the archive without "
                         "re-walking the API")
    ap.add_argument("--offline", action="store_true",
                    help="also skip the UniProt cross-check (which is one "
                         "request per signature and degrades to a blank "
                         "column if it fails)")
    args = ap.parse_args()

    def log(msg): print(msg, flush=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = [(f"enumerate {p}", False) for p in SEED_PFAMS] + [("merge", False)]

    counts = []
    per_sig = {}
    for i, pfam_id in enumerate(SEED_PFAMS):
        live("S2", steps)
        if not args.parse_only:
            log(f"[s2] enumerating {pfam_id} …")
            enumerate_signature(pfam_id, args.page_size, log)
        rows = parse_archive(pfam_id)
        per_sig[pfam_id] = rows
        st = load_state(pfam_id)
        api_count = st.get("count")
        uni = None if args.offline else uniprot_count(pfam_id)
        unique = len({r.accession for r in rows})
        counts.append({
            "pfam_id": pfam_id, "api_count": api_count if api_count else "",
            "uniprot_count": uni if uni is not None else "",
            "pages": st.get("pages", 0), "parsed": len(rows), "unique": unique,
            "delta_vs_api_count": (unique - int(api_count)) if api_count else "",
            "restarts": st.get("restarts", 0),
            "cursor_exhausted": int(bool(st.get("complete"))
                                    and st.get("next") is None),
        })
        steps[i] = (f"enumerate {pfam_id}", True)
        live("S2", steps)
        log(f"  {pfam_id}: {len(rows)} parsed / api count {api_count}")

    merged = merge(per_sig)
    write_tsv(OUT_DIR / "interpro_counts.tsv", counts,
              ["pfam_id", "api_count", "uniprot_count", "pages", "parsed",
               "unique", "delta_vs_api_count", "restarts",
               "cursor_exhausted"])
    # The merged rows are an intermediate: every column of them is carried
    # into census_v2.tsv, so committing both would put the same 1.7 MB in
    # the repo twice. It sits beside the raw pages it was parsed from.
    write_tsv(require_data_root(ROOT) / ARCHIVE / "interpro_records.tsv",
              merged, COLS)
    steps[-1] = ("merge", True)
    live("S2", steps)
    log(f"[s2] {len(merged)} unique accessions across {len(SEED_PFAMS)} seeds")

    incomplete = [c["pfam_id"] for c in counts if not c["cursor_exhausted"]]
    if incomplete:
        log(f"[s2] INCOMPLETE: {', '.join(incomplete)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
