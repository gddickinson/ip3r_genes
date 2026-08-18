"""InterPro REST client — protein family / domain annotation lookup.

For novel-paralog discovery we use InterPro/Pfam to ask one question per
candidate:  *does this protein carry the family's domain signature?*

The signatures are declared once in `src/utils/family.py`; for the ITPR
family they are PF08709 (Ins145_P3_rec, the IP3-binding core), PF02815
(MIR), PF01365 (RYDR_ITPR / RIH) and PF08454 (RIH_assoc).

A novel-paralog candidate must hit at least one of these. Without the
signature, an outlier sequence is just an outlier — not an ITPR paralog.

API:  https://www.ebi.ac.uk/interpro/api/
Docs: https://github.com/ProteinsWebTeam/interpro7-api
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests

from ..utils.family import FAMILY_PFAM_IDS

INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"


@dataclass
class DomainHit:
    accession: str        # the query UniProt accession
    pfam_id: str          # e.g. PF15917
    pfam_name: str        # human-readable name
    n_hits: int           # how many copies in this protein
    coverage: float = 0.0 # fraction of the protein covered


@dataclass
class ProteinWithDomain:
    """One row from the InterPro 'list proteins carrying Pfam X' endpoint."""
    accession: str
    name: str
    gene: str
    species: str
    taxon_id: Optional[int]
    length: int
    reviewed: bool        # SwissProt (True) vs TrEMBL (False)
    pfam_id: str
    in_alphafold: bool


def has_family_signature(domain_hits: list[DomainHit]) -> bool:
    return any(h.pfam_id in FAMILY_PFAM_IDS for h in domain_hits)


def parse_protein_row(meta: dict, pfam_id: str) -> Optional[ProteinWithDomain]:
    """One `metadata` object from the list endpoint → ProteinWithDomain.

    Shared by live pagination and by scripts that reload archived raw
    pages, so the parse stays identical in both paths."""
    acc = meta.get("accession", "")
    if not acc:
        return None
    org = meta.get("source_organism", {}) or {}
    taxid = org.get("taxId")
    try:
        taxid = int(taxid) if taxid is not None else None
    except (TypeError, ValueError):
        taxid = None
    try:
        length = int(meta.get("length", 0))
    except (TypeError, ValueError):
        length = 0
    return ProteinWithDomain(
        accession=acc,
        name=str(meta.get("name") or ""),
        gene=str(meta.get("gene") or ""),
        species=str(org.get("scientificName") or org.get("fullName") or ""),
        taxon_id=taxid,
        length=length,
        reviewed=(meta.get("source_database") == "reviewed"),
        pfam_id=pfam_id,
        in_alphafold=bool(meta.get("in_alphafold", False)),
    )


def list_proteins_with_pfam(
    pfam_id: str,
    max_results: Optional[int] = 500,
    page_size: int = 200,
    timeout_s: int = 60,
    polite_sleep_s: float = 0.2,
    dump_dir: Optional[Path] = None,
    strict: bool = False,
    max_retries: int = 4,
    stats: Optional[dict] = None,
) -> list[ProteinWithDomain]:
    """Return every UniProt protein InterPro has annotated with `pfam_id`.

    Paginates until `max_results` is reached — `max_results=None` paginates
    to exhaustion (census mode). The endpoint returns:

        GET /protein/UniProt/entry/pfam/{pfam_id}/?page_size=N

    with `next` cursor URLs and a top-level `count`. Each result has
    `metadata` (accession, name, gene, source_organism, length,
    source_database='reviewed'|'unreviewed', in_alphafold).

    Census-mode extras:
      dump_dir  — archive every raw page JSON as `{pfam_id}_page_NNNN.json`.
      strict    — raise RuntimeError on a page that still fails after
                  `max_retries` backed-off attempts (default: stop early
                  and return what was fetched, the legacy behavior).
      stats     — caller-supplied dict, filled with `count` (the API's own
                  total), `pages`, `fetched`, `complete`.
    """
    out: list[ProteinWithDomain] = []
    url: Optional[str] = (
        f"{INTERPRO_BASE}/protein/UniProt/entry/pfam/{pfam_id}/"
        f"?page_size={page_size}"
    )
    page_no = 0
    api_count: Optional[int] = None
    backoffs = [2, 5, 15, 30]
    while url and (max_results is None or len(out) < max_results):
        data = None
        for attempt in range(max_retries + 1):
            try:
                r = requests.get(url, timeout=timeout_s)
                r.raise_for_status()
                data = r.json()
                break
            except (requests.RequestException, ValueError):
                if attempt < max_retries:
                    time.sleep(backoffs[min(attempt, len(backoffs) - 1)])
        if data is None:
            if strict:
                raise RuntimeError(
                    f"InterPro pagination for {pfam_id} failed at page "
                    f"{page_no + 1} after {max_retries + 1} attempts ({url})")
            break
        page_no += 1
        if api_count is None:
            api_count = data.get("count")
        if dump_dir is not None:
            dump_dir.mkdir(parents=True, exist_ok=True)
            (dump_dir / f"{pfam_id}_page_{page_no:04d}.json").write_text(
                json.dumps(data))
        for entry in data.get("results", []):
            row = parse_protein_row(entry.get("metadata", {}), pfam_id)
            if row is None:
                continue
            out.append(row)
            if max_results is not None and len(out) >= max_results:
                break
        url = data.get("next")
        if url:
            time.sleep(polite_sleep_s)
    if stats is not None:
        stats["count"] = api_count
        stats["pages"] = page_no
        stats["fetched"] = len(out)
        stats["complete"] = (api_count is not None and len(out) >= api_count)
    return out


def fetch_uniprot_sequence(accession: str, timeout_s: int = 20) -> str:
    """Pull a single protein sequence by UniProt accession.

    Used to back-fill sequences for domain-scan candidates so the
    downstream MSA / tree have data to work with.
    """
    if not accession:
        return ""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        r = requests.get(url, timeout=timeout_s)
    except requests.RequestException:
        return ""
    if r.status_code != 200:
        return ""
    parts: list[str] = []
    for line in r.text.splitlines():
        if line.startswith(">") or not line.strip():
            continue
        parts.append(line.strip())
    return "".join(parts)


def fetch_domains_for_uniprot(uniprot_acc: str, timeout_s: int = 20) -> list[DomainHit]:
    """Return Pfam domain hits for a UniProt accession via InterPro."""
    if not uniprot_acc:
        return []
    url = f"{INTERPRO_BASE}/entry/pfam/protein/uniprot/{uniprot_acc}"
    try:
        r = requests.get(url, timeout=timeout_s)
    except requests.RequestException:
        return []
    if r.status_code == 204 or r.status_code == 404:
        return []
    if r.status_code != 200:
        return []
    try:
        data = r.json()
    except ValueError:
        return []

    out: list[DomainHit] = []
    for entry in data.get("results", []):
        meta = entry.get("metadata", {})
        pfam_id = meta.get("accession", "")
        pfam_name = meta.get("name", "")
        locations = []
        for proteins in entry.get("proteins", []):
            for loc in proteins.get("entry_protein_locations", []) or []:
                for frag in loc.get("fragments", []) or []:
                    try:
                        locations.append((int(frag["start"]), int(frag["end"])))
                    except (KeyError, TypeError, ValueError):
                        continue
        out.append(DomainHit(
            accession=uniprot_acc,
            pfam_id=pfam_id,
            pfam_name=pfam_name,
            n_hits=len(locations),
            coverage=0.0,  # we don't pull seq length here; left for future use
        ))
    return out


def batch_fetch_family_signatures(
    uniprot_accessions: list[str],
    sleep_s: float = 0.2,
    timeout_s: int = 20,
) -> dict[str, list[DomainHit]]:
    """Look up each UniProt accession and return the dict {acc: domain_hits}.

    Only domain hits in FAMILY_PFAM_IDS are kept — we don't need every Pfam
    domain in the protein, just whether the family signature is present.
    """
    out: dict[str, list[DomainHit]] = {}
    for acc in uniprot_accessions:
        hits = [h for h in fetch_domains_for_uniprot(acc, timeout_s) if h.pfam_id in FAMILY_PFAM_IDS]
        out[acc] = hits
        time.sleep(sleep_s)  # be polite to InterPro
    return out
