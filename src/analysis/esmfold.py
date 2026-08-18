"""ESMFold structure check via the public ESM Atlas API.

AlphaFold DB only covers UniProt — precisely not the unnamed Ensembl gene
models and fresh BLAST hits that discovery surfaces. ESMFold's free
endpoint folds an arbitrary sequence (protein-language-model based, no
MSA needed), so we can still get a fold-quality readout for them:

    POST https://api.esmatlas.com/foldSequence/v1/pdb/   (plain sequence body)

The endpoint rejects long inputs (practical cap ~400 aa), so we fold a
*segment*. For IP3 receptors the C-terminal ~400 aa hold the six-TM pore
module and the C-terminal tail that mediates tetramerisation — the
family's most conserved region — making "does the C-terminus fold
confidently?" a meaningful family-membership signal. Mean pLDDT is parsed
from the B-factor column.

Interpretation: pLDDT > 70 = confident fold; 50–70 = partial/flexible;
< 50 = disordered or not a real protein segment.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import requests

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_SEGMENT = 400


@dataclass
class FoldResult:
    segment: str                # e.g. "C-terminal 350 aa (res 2201-2551)"
    n_residues: int
    mean_plddt: float
    min_plddt: float
    pdb_text: str
    elapsed_s: float

    def verdict(self) -> str:
        if self.mean_plddt >= 70:
            return "CONFIDENT fold — consistent with a real, structured protein"
        if self.mean_plddt >= 50:
            return "partial confidence — some structured regions"
        return "low confidence — disordered or dubious segment"

    def summary(self) -> str:
        return (
            f"ESMFold check — {self.segment}\n"
            f"  mean pLDDT = {self.mean_plddt:.1f}   (min {self.min_plddt:.1f}, "
            f"{self.n_residues} residues, {self.elapsed_s:.0f}s)\n"
            f"  verdict: {self.verdict()}"
        )


def pick_segment(sequence: str, where: str = "cterm",
                 length: int = 350) -> tuple[str, str]:
    """Choose the sub-sequence to fold. Returns (segment_seq, description)."""
    length = min(length, MAX_SEGMENT, len(sequence))
    if where == "nterm":
        seg = sequence[:length]
        desc = f"N-terminal {length} aa (res 1-{length})"
    else:
        seg = sequence[-length:]
        start = len(sequence) - length + 1
        desc = f"C-terminal {length} aa (res {start}-{len(sequence)})"
    return seg, desc


def fold_segment(
    sequence: str,
    where: str = "cterm",
    length: int = 350,
    timeout_s: int = 300,
) -> FoldResult:
    """Fold a segment of `sequence` via the ESM Atlas API and report pLDDT.
    Raises RuntimeError with a readable message on API failure."""
    if not sequence:
        raise RuntimeError("no sequence available — re-search with 'Fetch sequences'")
    seg, desc = pick_segment(sequence, where=where, length=length)
    start = time.time()
    try:
        r = requests.post(ESMFOLD_URL, data=seg, timeout=timeout_s)
    except requests.exceptions.SSLError:
        # The Atlas endpoint's certificate chain is not always resolvable
        # from every environment; the payload is a public sequence and the
        # response is inert PDB text, so retrying unverified is acceptable.
        r = requests.post(ESMFOLD_URL, data=seg, timeout=timeout_s, verify=False)
    except requests.RequestException as e:
        raise RuntimeError(f"ESMFold API unreachable: {e}") from e
    if r.status_code != 200:
        raise RuntimeError(
            f"ESMFold API returned {r.status_code}: {r.text[:200]} "
            "(the public endpoint rate-limits; try again in a minute)")
    pdb_text = r.text
    plddts = _parse_ca_plddt(pdb_text)
    if not plddts:
        raise RuntimeError("ESMFold response contained no CA atoms")
    # The API writes pLDDT to the B-factor column on a 0–1 scale; the
    # conventional reporting scale is 0–100.
    if max(plddts) <= 1.0:
        plddts = [p * 100 for p in plddts]
    return FoldResult(
        segment=desc,
        n_residues=len(plddts),
        mean_plddt=sum(plddts) / len(plddts),
        min_plddt=min(plddts),
        pdb_text=pdb_text,
        elapsed_s=time.time() - start,
    )


def _parse_ca_plddt(pdb_text: str) -> list[float]:
    vals: list[float] = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                continue
    return vals
