"""S11 shared library — structure acquisition, parsing and TM-align.

S11 asks whether what the census *calls* an ITPR **folds** like one, and
puts a structural frame under the constraint work. Four ingredients:

  1. **Predicted models** for the S6 representatives and for the census at
     large. AlphaFold DB is the cheap source; this family is right on its
     length ceiling (AFDB's monomer pipeline stops at 2,700 residues and a
     vertebrate ITPR is ~2,700), so *how much* AFDB holds is itself the
     first result rather than an assumption — see `s11_afdb_probe.py`.
  2. **Experimental references** — cryo-EM IP3R and RyR entries, resolved
     by an RCSB *query* in `s11_refs.py`, never from a planning document.
  3. **Negative controls** — an unrelated channel of similar size and a
     large non-channel, so the TM-score scale the family is read against
     has a measured floor and is not asserted.
  4. **TM-align** to score every model against those references.

Everything cacheable lands under `<data_root>/structures/` and
`<data_root>/raw_api/s11/`; nothing large enters the repo.

Structure parsing lives in `s11_struct_io.py` and TM-align in
`s11_tmalign.py`; this module is acquisition, caching and table access.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils.data_root import get_data_root  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
PDB_CIF = "https://files.rcsb.org/download"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "ip3r_genes/S11 (george.dickinson@gmail.com)"



def structures_dir() -> Path:
    d = get_data_root(PROJECT_ROOT) / "structures" / "s11"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cache_dir() -> Path:
    d = get_data_root(PROJECT_ROOT) / "raw_api" / "s11"
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir() -> Path:
    d = PROJECT_ROOT / "results" / "structures"
    d.mkdir(parents=True, exist_ok=True)
    return d


# --------------------------------------------------------------------------
# HTTP with an on-disk cache (EBI and RCSB both go flaky under load; S9 and
# S10 both needed retries, so S11 assumes the same).
# --------------------------------------------------------------------------

def http_get(url: str, timeout: int = 60, retries: int = 3,
             backoff: float = 2.0, data: bytes | None = None,
             content_type: str = "") -> bytes:
    last: Exception | None = None
    for attempt in range(retries):
        headers = {"User-Agent": USER_AGENT}
        if content_type:
            headers["Content-Type"] = content_type
        req = urllib.request.Request(url, headers=headers, data=data)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 204):
                raise            # a real "not there" — do not burn retries
            last = exc
        except Exception as exc:  # noqa: BLE001 - transport errors vary
            last = exc
        if attempt < retries - 1:
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries}: {url}: {last}")


def cached_json(url: str, key: str, timeout: int = 60,
                post: dict | None = None) -> Optional[dict | list]:
    """GET/POST JSON with a permanent on-disk cache.

    `None` means the server said 404/204 — recorded as the literal `null`
    so a re-run does not re-ask. Everything else raises.
    """
    path = cache_dir() / f"{key}.json"
    if path.exists():
        txt = path.read_text()
        return None if txt == "null" else json.loads(txt)
    body = json.dumps(post).encode() if post is not None else None
    try:
        raw = http_get(url, timeout=timeout, data=body,
                       content_type="application/json" if body else "")
    except urllib.error.HTTPError as exc:
        if exc.code in (404, 204):
            path.write_text("null")
            return None
        raise
    if not raw.strip():
        path.write_text("null")
        return None
    payload = json.loads(raw.decode())
    path.write_text(json.dumps(payload))
    return payload


# --------------------------------------------------------------------------
# AlphaFold DB
# --------------------------------------------------------------------------

@dataclass
class AfdbEntry:
    accession: str                 # what was asked for
    found: bool
    model_accession: str = ""      # what AFDB actually returned
    entry_id: str = ""
    uniprot_start: int = 0
    uniprot_end: int = 0
    model_length: int = 0
    n_records: int = 0
    version: int = 0
    seq_version_date: str = ""
    cif_url: str = ""
    pdb_url: str = ""
    error: str = ""

    @property
    def is_canonical(self) -> bool:
        """Is the served model the accession that was asked for?

        AFDB keys its API on an accession but may answer with an
        **isoform** record (`Q14643-4`). Two of the three human ITPR
        paralogs are served that way, and one of those isoforms is 181
        residues of a 2,701-residue protein — so "AFDB has a model for
        this accession" is not the same statement as "this protein is
        modelled", and S11 never conflates them.
        """
        return bool(self.model_accession) and self.model_accession == self.accession

    def coverage(self, protein_length: int) -> float:
        if not protein_length or not self.model_length:
            return 0.0
        return round(self.model_length / protein_length, 4)


def afdb_probe(accession: str) -> AfdbEntry:
    """Ask AFDB what model it holds for this UniProt accession.

    Every record in the response is considered, not `payload[0]`: the
    canonical record is preferred, and the longest is the fallback, so a
    short isoform can never displace a full-length model that was also on
    offer. What was served is recorded either way.
    """
    try:
        payload = cached_json(f"{AFDB_API}/{accession}", f"afdb_{accession}")
    except Exception as exc:  # noqa: BLE001
        return AfdbEntry(accession, False, error=str(exc)[:120])
    if not payload:
        return AfdbEntry(accession, False, error="404")
    records = list(payload)

    def rank(rec: dict) -> tuple:
        served = rec.get("uniprotAccession") or ""
        span = int(rec.get("uniprotEnd") or 0) - int(rec.get("uniprotStart") or 0) + 1
        return (0 if served == accession else 1, -span, served)

    rec = min(records, key=rank)
    start = int(rec.get("uniprotStart") or 0)
    end = int(rec.get("uniprotEnd") or 0)
    return AfdbEntry(
        accession=accession,
        found=True,
        model_accession=rec.get("uniprotAccession") or "",
        entry_id=rec.get("entryId") or "",
        uniprot_start=start,
        uniprot_end=end,
        model_length=max(0, end - start + 1),
        n_records=len(records),
        version=int(rec.get("latestVersion") or 0),
        seq_version_date=(rec.get("sequenceVersionDate") or "")[:10],
        cif_url=rec.get("cifUrl", ""),
        pdb_url=rec.get("pdbUrl", ""),
    )


def afdb_fetch(entry: AfdbEntry, dest: Path) -> Optional[Path]:
    """Download the AFDB PDB file for a probed entry (pLDDT in B-factors)."""
    if not entry.found or not entry.pdb_url:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.write_bytes(http_get(entry.pdb_url))
    return dest


# --------------------------------------------------------------------------
# PDB / RCSB
# --------------------------------------------------------------------------

def pdb_fetch(pdb_id: str, dest_dir: Path | None = None) -> Path:
    """Download an mmCIF from RCSB.

    mmCIF and not PDB format by necessity: an IP3R or RyR tetramer is far
    past the 99,999-atom PDB ceiling and RCSB does not distribute a legacy
    PDB file for these entries at all.
    """
    dest_dir = dest_dir or (structures_dir() / "reference")
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{pdb_id.lower()}.cif"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.write_bytes(http_get(f"{PDB_CIF}/{pdb_id.upper()}.cif"))
    return dest


def rcsb_search(query: dict, key: str, rows: int = 200) -> list[str]:
    """Run an RCSB search-API query and return the entry ids, cached."""
    payload = {
        "query": query,
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_verbosity": "compact",
        },
    }
    got = cached_json(RCSB_SEARCH, key, post=payload)
    if not got:
        return []
    return list(got.get("result_set") or [])


def rcsb_entry(pdb_id: str) -> Optional[dict]:
    """Entry-level metadata: method, resolution, title, deposit date."""
    return cached_json(f"{RCSB_DATA}/{pdb_id.upper()}",
                       f"rcsb_entry_{pdb_id.lower()}")


def uniprot_record(accession: str) -> Optional[dict]:
    url = (f"{UNIPROT_REST}/{accession}?fields=accession,id,length,"
           "protein_name,gene_names,organism_name,sequence,xref_pdb,"
           "xref_alphafolddb,xref_pfam")
    return cached_json(url, f"uniprot_{accession}")


def protein_facts(accession: str) -> dict:
    """Name, gene and Pfam set for one accession — what turns "these look
    like MIR proteins" into a measurement."""
    rec = uniprot_record(accession)
    if not rec:
        return {"name": "", "gene": "", "pfams": ""}
    desc = rec.get("proteinDescription") or {}
    name = ((desc.get("recommendedName") or {}).get("fullName") or {}
            ).get("value", "")
    if not name:
        subs = desc.get("submissionNames") or []
        name = ((subs[0].get("fullName") or {}).get("value", "")
                if subs else "")
    genes = rec.get("genes") or [{}]
    gene = ((genes[0].get("geneName") or {}).get("value") or "")
    pfams = sorted({x.get("id", "") for x in (rec.get("uniProtKBCrossReferences") or [])
                    if x.get("database") == "Pfam"} - {""})
    return {"name": name, "gene": gene, "pfams": ";".join(pfams)}


# --------------------------------------------------------------------------
# Project-table access
# --------------------------------------------------------------------------

def load_tsv(path: Path) -> list[dict]:
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_representatives(path: Path | None = None) -> list[dict]:
    path = path or (PROJECT_ROOT / "results" / "msa_v2" / "representatives.tsv")
    return load_tsv(path)


def load_census(path: Path | None = None) -> list[dict]:
    path = path or (PROJECT_ROOT / "results" / "census_v6" / "census_v6.tsv")
    return load_tsv(path)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
