"""S11 — the RCSB query layer.

Brief step 1 is emphatic: *resolve the reference structures by query, not
from memory*. A PDB id typed out of a planning document is how the wrong
protein ends up in a figure, and this family makes that easy — the IP3R and
RyR entries share every diagnostic Pfam and half their titles.

So the candidate set is **enumerated** from RCSB by the family's own Pfam
signatures (`src/utils/family.py`, the one place the family is defined),
restricted to entries with an experimental method and a recorded
resolution, and every candidate is written to a committed table before
anything is selected. What is asked of RCSB is only *which entries exist*;
which family each one is gets decided by this project's own instruments in
`s11_refs.py`, never by the entry title.

Two endpoints, both cached permanently under `<data_root>/raw_api/s11/`:
the search API for enumeration and the GraphQL data API for entry and
polymer-entity metadata. GraphQL because the alternative is three REST
calls per entity and ~600 of them for this candidate set.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import (PROJECT_ROOT, cache_dir, cached_json, http_get,  # noqa: E402
                     RCSB_SEARCH)

GRAPHQL = "https://data.rcsb.org/graphql"

sys.path.insert(0, str(PROJECT_ROOT))
from src.utils import family  # noqa: E402

#: The signatures the candidate set is enumerated from. The union of all
#: four and not just PF08709: RCSB's own Pfam annotation of human ITPR3
#: (6DQN) carries PF02815/PF08454/PF01365/PF00520 and **not** PF08709, the
#: signature that names the family — the same absence S2 measured across
#: 2,911 of 15,417 enumerated proteins. Enumerating on the naming signature
#: alone would have missed the project's own IP3R reference.
ENUM_PFAMS = list(family.FAMILY_PFAM_IDS)

ENTRY_FIELDS = """
  rcsb_id
  struct { title }
  exptl { method }
  rcsb_entry_info { resolution_combined deposited_polymer_monomer_count
                    polymer_entity_count deposited_model_count }
  rcsb_accession_info { initial_release_date }
  nonpolymer_entities {
    rcsb_nonpolymer_entity_container_identifiers { nonpolymer_comp_id }
    rcsb_nonpolymer_entity { pdbx_description } }
  polymer_entities {
    rcsb_id
    entity_poly { rcsb_sample_sequence_length pdbx_seq_one_letter_code_can }
    rcsb_polymer_entity { pdbx_description }
    rcsb_polymer_entity_container_identifiers {
      reference_sequence_identifiers { database_accession database_name } }
    rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
    rcsb_polymer_entity_annotation { annotation_id type }
    polymer_entity_instances {
      rcsb_polymer_entity_instance_container_identifiers { auth_asym_id } } }
"""


# --------------------------------------------------------------------------
# Search
# --------------------------------------------------------------------------

def _text(attribute: str, value: str, operator: str = "exact_match") -> dict:
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator,
                           "value": value}}


def _exists(attribute: str) -> dict:
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": "exists"}}


def _run_search(query: dict, key: str, rows: int = 1000) -> list[str]:
    payload = {"query": query, "return_type": "entry",
               "request_options": {"paginate": {"start": 0, "rows": rows},
                                   "results_verbosity": "compact"}}
    got = cached_json(RCSB_SEARCH, key, post=payload, timeout=120)
    if not got:
        return []
    return list(got.get("result_set") or [])


def entries_for_pfam(pfam: str) -> list[str]:
    """Every entry whose polymer entity carries this Pfam, with a resolution.

    Deliberately **not** filtered to cryo-EM at this stage. The brief asks
    for cryo-EM references, but a candidate table that only ever contained
    cryo-EM entries could not show that the X-ray structures of this family
    are all fragments of the binding core — which is what makes "cryo-EM"
    the right restriction rather than a preference.
    """
    query = {"type": "group", "logical_operator": "and", "nodes": [
        _text("rcsb_polymer_entity_annotation.annotation_id", pfam),
        _exists("rcsb_entry_info.resolution_combined"),
    ]}
    return _run_search(query, f"rcsb_search_pfam_{pfam}")


def entries_for_uniprot(accession: str) -> list[str]:
    """Every entry with a polymer entity mapped to this UniProt accession."""
    query = {"type": "group", "logical_operator": "and", "nodes": [
        _text("rcsb_polymer_entity_container_identifiers."
              "reference_sequence_identifiers.database_accession", accession),
        _text("rcsb_polymer_entity_container_identifiers."
              "reference_sequence_identifiers.database_name", "UniProt"),
        _exists("rcsb_entry_info.resolution_combined"),
    ]}
    return _run_search(query, f"rcsb_search_up_{accession}")


def family_candidate_entries() -> tuple[list[str], dict[str, list[str]]]:
    """Union of the family Pfams' entries, plus which Pfam found each."""
    found: dict[str, list[str]] = {}
    for pfam in ENUM_PFAMS:
        for pid in entries_for_pfam(pfam):
            found.setdefault(pid, []).append(pfam)
    return sorted(found), found


# --------------------------------------------------------------------------
# GraphQL entry metadata
# --------------------------------------------------------------------------

def _graphql(entry_ids: list[str], key: str) -> dict:
    path = cache_dir() / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    ids = ", ".join(f'"{e}"' for e in entry_ids)
    query = "{ entries(entry_ids: [%s]) {%s} }" % (ids, ENTRY_FIELDS)
    raw = http_get(GRAPHQL, data=json.dumps({"query": query}).encode(),
                   content_type="application/json", timeout=180)
    payload = json.loads(raw.decode())
    if payload.get("errors"):
        raise RuntimeError(f"RCSB GraphQL: {payload['errors'][:1]}")
    path.write_text(json.dumps(payload))
    return payload


def fetch_entries(entry_ids: list[str], batch: int = 40) -> dict[str, dict]:
    """Entry metadata for many ids, batched and cached by batch content."""
    out: dict[str, dict] = {}
    ordered = sorted(set(entry_ids))
    for i in range(0, len(ordered), batch):
        chunk = ordered[i:i + batch]
        key = "rcsb_gql_" + _digest(chunk)
        payload = _graphql(chunk, key)
        for rec in (payload.get("data", {}).get("entries") or []):
            if rec:
                out[rec["rcsb_id"]] = rec
    return out


def _digest(items: list[str]) -> str:
    import hashlib
    return hashlib.sha256("|".join(items).encode()).hexdigest()[:16]


# --------------------------------------------------------------------------
# Flattening: one row per polymer entity
# --------------------------------------------------------------------------

def entity_rows(entry: dict) -> list[dict]:
    """Flatten one entry into one row per polymer entity.

    Chains are sorted so `first_chain` is deterministic (D24) — RCSB
    returns instance order that varies between entries.
    """
    info = entry.get("rcsb_entry_info") or {}
    res = (info.get("resolution_combined") or [None])[0]
    methods = sorted({m.get("method", "") for m in (entry.get("exptl") or [])})
    ligands = sorted({
        (ne.get("rcsb_nonpolymer_entity_container_identifiers") or {})
        .get("nonpolymer_comp_id", "")
        for ne in (entry.get("nonpolymer_entities") or [])} - {""})
    rows = []
    for pe in (entry.get("polymer_entities") or []):
        ids = (pe.get("rcsb_polymer_entity_container_identifiers") or {})
        refs = ids.get("reference_sequence_identifiers") or []
        accs = [r.get("database_accession", "") for r in refs
                if r.get("database_name") == "UniProt"]
        pfams = sorted({a.get("annotation_id", "")
                        for a in (pe.get("rcsb_polymer_entity_annotation") or [])
                        if a.get("type") == "Pfam"})
        chains = sorted(
            (i.get("rcsb_polymer_entity_instance_container_identifiers") or {})
            .get("auth_asym_id", "")
            for i in (pe.get("polymer_entity_instances") or []))
        chains = [c for c in chains if c]
        orgs = sorted({(o.get("ncbi_scientific_name") or "")
                       for o in (pe.get("rcsb_entity_source_organism") or [])})
        poly = pe.get("entity_poly") or {}
        rows.append({
            "entry_id": entry.get("rcsb_id", ""),
            "entity_id": pe.get("rcsb_id", ""),
            "title": (entry.get("struct") or {}).get("title", "") or "",
            "method": ";".join(methods),
            "resolution": res,
            "released": ((entry.get("rcsb_accession_info") or {})
                         .get("initial_release_date") or "")[:10],
            "n_polymer_entities": info.get("polymer_entity_count") or 0,
            "deposited_residues": info.get("deposited_polymer_monomer_count") or 0,
            "description": ((pe.get("rcsb_polymer_entity") or {})
                            .get("pdbx_description") or ""),
            "uniprot": ";".join(accs),
            "organism": ";".join(o for o in orgs if o),
            "sample_length": poly.get("rcsb_sample_sequence_length") or 0,
            "sequence": (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
            "pfams": ";".join(pfams),
            "chains": ";".join(chains),
            "n_chains": len(chains),
            "first_chain": chains[0] if chains else "",
            "ligands": ";".join(ligands),
        })
    return rows
