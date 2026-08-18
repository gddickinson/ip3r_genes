"""S1 control panels — the ground truth the discovery scorer is tested on.

Two panels, both fetched live from UniProt and cached as JSON + FASTA so the
benchmark reruns offline and the exact panel that produced a number is
committed alongside it.

POSITIVES — real family members that must score >= 40 when their own name is
held out of `known_paralogs`:
  * ITPR1 / ITPR2 / ITPR3 orthologs across a vertebrate species panel, one
    hold-out run each;
  * the invertebrate / non-metazoan single-Itpr grade (fly Itp-r83A, worm
    itr-1, Dictyostelium iplA, sea urchin). These are never name-protected
    by `known_paralogs` at all, so they are scored in the baseline run —
    the real use case is an unnamed true family member surfacing.

DECOYS — must score < 40 in every run. The panel is built to attack specific
scorer components:
  * RYR1/2/3 across species — the sharp decoy. They carry *every*
    ITPR-diagnostic Pfam (roadmap D14), so they take the +20 domain point by
    construction and are only separable on length and on distance.
  * POMT1/POMT2 — share the MIR domain (PF02815) and nothing else.
  * in-band channels and in-band non-channels (2,000-3,600 aa) — take the
    +15 size point by construction.
  * out-of-band giants — length-filter sanity checks.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from src.core.models import ProteinVariant  # noqa: E402

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_POLITE_SLEEP_S = 0.25
_FIELDS = ("accession,id,protein_name,gene_primary,organism_name,length,"
           "sequence,reviewed")

#: (gene symbol, taxon id, group). `group` is the hold-out run a member is
#: scored in: ITPR1/ITPR2/ITPR3 name a hold-out; `invert_grade` is scored in
#: the baseline run because no hold-out can protect an unnamed member.
POSITIVE_SPEC: list[tuple[str, int, str]] = [
    ("ITPR1", 9606, "ITPR1"), ("Itpr1", 10090, "ITPR1"),
    ("Itpr1", 10116, "ITPR1"), ("ITPR1", 9031, "ITPR1"),
    ("itpr1", 8364, "ITPR1"), ("itpr1a", 7955, "ITPR1"),
    ("itpr1b", 7955, "ITPR1"), ("ITPR1", 9913, "ITPR1"),
    ("ITPR2", 9606, "ITPR2"), ("Itpr2", 10090, "ITPR2"),
    ("Itpr2", 10116, "ITPR2"), ("ITPR2", 9031, "ITPR2"),
    ("itpr2", 8364, "ITPR2"), ("itpr2", 7955, "ITPR2"),
    ("ITPR2", 9913, "ITPR2"),
    ("ITPR3", 9606, "ITPR3"), ("Itpr3", 10090, "ITPR3"),
    ("Itpr3", 10116, "ITPR3"), ("ITPR3", 9031, "ITPR3"),
    ("itpr3", 8364, "ITPR3"), ("itpr3", 7955, "ITPR3"),
    ("ITPR3", 9913, "ITPR3"),
    # The single-Itpr grade — true family members wearing a name no
    # vertebrate-symbol filter recognises.
    ("Itp-r83A", 7227, "invert_grade"),   # Drosophila melanogaster
    ("itr-1", 6239, "invert_grade"),      # Caenorhabditis elegans
    ("iplA", 44689, "invert_grade"),      # Dictyostelium discoideum
    ("Itpr", 7668, "invert_grade"),       # Strongylocentrotus purpuratus
]

#: (gene symbol, taxon id, category).
DECOY_SPEC: list[tuple[str, int, str]] = [
    # --- the sharp decoy: shares every diagnostic Pfam (D14) ---
    ("RYR1", 9606, "RyR (sister family)"),
    ("RYR2", 9606, "RyR (sister family)"),
    ("RYR3", 9606, "RyR (sister family)"),
    ("Ryr1", 10090, "RyR (sister family)"),
    ("Ryr2", 10090, "RyR (sister family)"),
    ("ryr3", 7955, "RyR (sister family)"),
    # --- shares the MIR domain only ---
    ("POMT1", 9606, "MIR-domain sharer"),
    ("POMT2", 9606, "MIR-domain sharer"),
    # --- channels inside or near the family size band ---
    ("CACNA1A", 9606, "in-band channel"),
    ("CACNA1C", 9606, "in-band channel"),
    ("CACNA1E", 9606, "in-band channel"),
    ("SCN1A", 9606, "in-band channel"),
    ("SCN5A", 9606, "in-band channel"),
    ("SCN9A", 9606, "in-band channel"),
    ("TRPM6", 9606, "in-band channel"),
    ("TRPM7", 9606, "in-band channel"),
    ("PKD1L1", 9606, "in-band channel"),
    ("PKD1L2", 9606, "in-band channel"),
    # --- non-channels inside the family size band ---
    ("TLN1", 9606, "in-band non-channel"),
    ("FLNA", 9606, "in-band non-channel"),
    ("SPTBN1", 9606, "in-band non-channel"),
    ("MYO7A", 9606, "in-band non-channel"),
    ("LRRK2", 9606, "in-band non-channel"),
    ("ATM", 9606, "in-band non-channel"),
    ("DOCK9", 9606, "in-band non-channel"),
    ("UTRN", 9606, "in-band non-channel"),
    ("VPS13A", 9606, "in-band non-channel"),
    ("NBEA", 9606, "in-band non-channel"),
    ("Tln1", 10090, "in-band non-channel"),
    # --- giants outside the band ---
    ("PKD1", 9606, "giant (out of band)"),
    ("DYNC1H1", 9606, "giant (out of band)"),
]


def _fetch_one(gene: str, taxon_id: int, note: str, kind: str) -> ProteinVariant | None:
    """Best UniProt entry for (gene, taxon): reviewed first, then longest."""
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{taxon_id}",
        "fields": _FIELDS, "size": "25", "format": "json",
    }
    r = requests.get(UNIPROT_SEARCH, params=params, timeout=45)
    r.raise_for_status()
    results = r.json().get("results", [])
    if not results:
        return None

    def _is_reviewed(e: dict) -> bool:
        # "UniProtKB unreviewed (TrEMBL)" CONTAINS "reviewed" — test for the
        # negative form first or every TrEMBL entry ranks as Swiss-Prot.
        t = (e.get("entryType") or "").lower()
        return "unreviewed" not in t and "reviewed" in t

    def rank(e: dict) -> tuple[int, int]:
        return (int(_is_reviewed(e)),
                len((e.get("sequence") or {}).get("value", "")))

    e = max(results, key=rank)
    genes = e.get("genes") or [{}]
    symbol = (genes[0].get("geneName") or {}).get("value") or gene
    seq = (e.get("sequence") or {}).get("value", "")
    acc = e.get("primaryAccession", "")
    reviewed = _is_reviewed(e)
    prot = ((e.get("proteinDescription") or {}).get("recommendedName") or {})
    prot_name = (prot.get("fullName") or {}).get("value", "")
    prefix = "DECOY" if kind == "decoy" else "CONTROL"
    return ProteinVariant(
        source="UniProt",
        accession=acc,
        gene_symbol=symbol,
        species=(e.get("organism") or {}).get("scientificName", ""),
        taxon_id=taxon_id,
        length_aa=len(seq) or None,
        description=f"{prefix} [{note}] {prot_name}",
        sequence=seq,
        url=f"https://www.uniprot.org/uniprotkb/{acc}",
        raw={"panel_note": note, "panel_kind": kind, "queried_symbol": gene,
             "reviewed": reviewed},
    )


def _write_fasta(path: Path, variants: list[ProteinVariant]) -> None:
    with path.open("w") as f:
        for v in variants:
            f.write(f">{v.accession}|{v.gene_symbol}|"
                    f"{v.species.replace(' ', '_')}|{v.raw.get('panel_note','')}\n")
            for i in range(0, len(v.sequence), 60):
                f.write(v.sequence[i:i + 60] + "\n")


def _to_json(variants: list[ProteinVariant]) -> str:
    return json.dumps([{
        "source": v.source, "accession": v.accession,
        "gene_symbol": v.gene_symbol, "species": v.species,
        "taxon_id": v.taxon_id, "length_aa": v.length_aa,
        "description": v.description, "sequence": v.sequence,
        "url": v.url, "raw": v.raw,
    } for v in variants], indent=1)


def fetch_panel(spec: list[tuple[str, int, str]], kind: str, cache_dir: Path,
                log=print) -> list[ProteinVariant]:
    """Fetch (or reload) one panel. `kind` is 'positive' or 'decoy'."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"panel_{kind}s.json"
    if cache.exists():
        records = json.loads(cache.read_text())
        log(f"[panel:{kind}] loaded {len(records)} from cache {cache.name}")
        return [ProteinVariant(**r) for r in records]

    out: list[ProteinVariant] = []
    missing: list[str] = []
    for gene, taxon_id, note in spec:
        try:
            v = _fetch_one(gene, taxon_id, note, kind)
        except requests.RequestException as exc:
            log(f"[panel:{kind}] {gene}/{taxon_id}: fetch failed ({exc})")
            missing.append(f"{gene}/{taxon_id}")
            continue
        if v is None or not v.sequence:
            log(f"[panel:{kind}] {gene}/{taxon_id}: no UniProt entry — skipped")
            missing.append(f"{gene}/{taxon_id}")
            continue
        flag = "SP" if v.raw.get("reviewed") else "TR"
        log(f"[panel:{kind}] {v.gene_symbol:<10} {v.species:<28} "
            f"{v.accession} {v.length_aa:>5} aa [{flag}] {note}")
        out.append(v)
        time.sleep(_POLITE_SLEEP_S)

    cache.write_text(_to_json(out))
    _write_fasta(cache_dir / f"panel_{kind}s.fasta", out)
    if missing:
        (cache_dir / f"panel_{kind}s_missing.txt").write_text(
            "\n".join(missing) + "\n")
    log(f"[panel:{kind}] cached {len(out)} entries "
        f"({len(missing)} unresolved) → {cache.name}")
    return out


if __name__ == "__main__":
    out_dir = PROJECT_ROOT / "results" / "benchmark_controls"
    fetch_panel(POSITIVE_SPEC, "positive", out_dir)
    fetch_panel(DECOY_SPEC, "decoy", out_dir)
