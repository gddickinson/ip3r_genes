"""Selection-pressure analysis — pairwise dN/dS (Nei–Gojobori 1986).

The decisive test for "is this unnamed gene model a real gene?": a
functional protein-coding gene accumulates far fewer non-synonymous than
synonymous substitutions (dN/dS ≪ 1, purifying selection), while a
pseudogene or annotation artifact drifts toward dN/dS ≈ 1 — the standard
evidence that an unnamed gene model is a real, translated gene.

Everything here is API-only and pure Python:

    1. `fetch_cds(protein_or_gene_id)` — Ensembl REST: protein → parent
       transcript → CDS.
    2. `closest_ortholog(gene_id)` — Ensembl Compara orthologues of the
       *candidate itself*; the highest-%id target gives an unsaturated
       comparison partner (distant partners saturate dS).
    3. `ng86(cds_a, cds_b)` — codon-aligns via the existing protein
       `pairwise_align()`, then Nei–Gojobori counting with all-pathway
       averaging and Jukes–Cantor correction.

Interpretation guide (rendered into the report):
    dN/dS < 0.5   strong purifying selection — functional gene
    0.5 – 0.9     weak purifying selection
    ~1            neutral drift — pseudogene-like
    > 1           positive selection (rare; verify alignment)
    dS > 1.5      saturation warning — treat the ratio as qualitative
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from itertools import permutations
from typing import Optional

import requests

from .alignment import pairwise_align

ENSEMBL_REST = "https://rest.ensembl.org"
_PAUSE_S = 0.12

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_BASES = "TCAG"


@dataclass
class DnDsResult:
    dn: float
    ds: float
    ratio: Optional[float]          # None when dS ≈ 0 or saturated-invalid
    n_codons: int                   # aligned, gap-free, stop-free codon pairs
    syn_sites: float
    nonsyn_sites: float
    syn_subs: float
    nonsyn_subs: float
    saturated: bool = False
    notes: list[str] = field(default_factory=list)

    def verdict(self) -> str:
        if self.ratio is None:
            return "indeterminate (dS unresolvable)"
        if self.saturated:
            qual = "qualitative only — dS saturated"
        else:
            qual = ""
        if self.ratio < 0.5:
            v = "STRONG PURIFYING SELECTION — behaves like a functional gene"
        elif self.ratio < 0.9:
            v = "weak purifying selection"
        elif self.ratio <= 1.1:
            v = "neutral drift — pseudogene-like"
        else:
            v = "possible positive selection (verify alignment)"
        return f"{v}{' (' + qual + ')' if qual else ''}"

    def summary(self) -> str:
        r = "n/a" if self.ratio is None else f"{self.ratio:.3f}"
        lines = [
            f"dN/dS = {r}   (dN={self.dn:.4f}, dS={self.ds:.4f}, "
            f"{self.n_codons} aligned codons)",
            f"  sites: {self.nonsyn_sites:.0f} nonsyn / {self.syn_sites:.0f} syn; "
            f"substitutions: {self.nonsyn_subs:.1f} nonsyn / {self.syn_subs:.1f} syn",
            f"  verdict: {self.verdict()}",
        ]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


# ---- Nei–Gojobori machinery ----------------------------------------------

def _codon_sites(codon: str) -> tuple[float, float]:
    """(synonymous, nonsynonymous) site counts for one codon."""
    aa = GENETIC_CODE.get(codon)
    if aa is None or aa == "*":
        return 0.0, 0.0
    syn = 0.0
    for pos in range(3):
        for b in _BASES:
            if b == codon[pos]:
                continue
            mut = codon[:pos] + b + codon[pos + 1:]
            if GENETIC_CODE.get(mut) == aa:
                syn += 1 / 3
    return syn, 3.0 - syn

def _path_subs(c1: str, c2: str) -> Optional[tuple[float, float]]:
    """(syn, nonsyn) substitution counts averaged over all mutational
    pathways between two codons; None if every path crosses a stop."""
    diffs = [i for i in range(3) if c1[i] != c2[i]]
    if not diffs:
        return 0.0, 0.0
    totals: list[tuple[float, float]] = []
    for order in permutations(diffs):
        cur, syn, non, ok = c1, 0.0, 0.0, True
        for pos in order:
            nxt = cur[:pos] + c2[pos] + cur[pos + 1:]
            aa1, aa2 = GENETIC_CODE.get(cur), GENETIC_CODE.get(nxt)
            if aa1 in (None, "*") or aa2 in (None, "*"):
                ok = False
                break
            if aa1 == aa2:
                syn += 1
            else:
                non += 1
            cur = nxt
        if ok:
            totals.append((syn, non))
    if not totals:
        return None
    return (sum(t[0] for t in totals) / len(totals),
            sum(t[1] for t in totals) / len(totals))

def _jukes_cantor(p: float) -> Optional[float]:
    if p < 1e-9:
        return 0.0
    if p >= 0.749:
        return None  # beyond correctable range
    return -0.75 * math.log(1 - (4.0 / 3.0) * p)


def ng86(cds_a: str, cds_b: str) -> DnDsResult:
    """Pairwise NG86 dN/dS. CDS are codon-aligned via their translations
    using the analysis module's existing pairwise protein aligner."""
    prot_a = translate(cds_a)
    prot_b = translate(cds_b)
    aln_a, aln_b, _ = pairwise_align(prot_a, prot_b)

    notes: list[str] = []
    ia = ib = 0
    S = N = Sd = Nd = 0.0
    n_codons = 0
    for col in range(len(aln_a)):
        ca, cb = aln_a[col], aln_b[col]
        if ca == "-" or cb == "-":
            ia += ca != "-"
            ib += cb != "-"
            continue
        codon_a = cds_a[ia * 3: ia * 3 + 3]
        codon_b = cds_b[ib * 3: ib * 3 + 3]
        ia += 1
        ib += 1
        if len(codon_a) < 3 or len(codon_b) < 3:
            continue
        if GENETIC_CODE.get(codon_a, "*") == "*" or GENETIC_CODE.get(codon_b, "*") == "*":
            continue
        subs = _path_subs(codon_a, codon_b)
        if subs is None:
            continue
        sa, na = _codon_sites(codon_a)
        sb, nb = _codon_sites(codon_b)
        S += (sa + sb) / 2
        N += (na + nb) / 2
        Sd += subs[0]
        Nd += subs[1]
        n_codons += 1

    if n_codons < 30:
        notes.append(f"only {n_codons} comparable codons — result unreliable")
    ps = Sd / S if S else 0.0
    pn = Nd / N if N else 0.0
    ds = _jukes_cantor(ps)
    dn = _jukes_cantor(pn)
    saturated = ds is None or (ds is not None and ds > 1.5)
    if ds is None:
        notes.append("dS beyond Jukes–Cantor correctable range (deep saturation)")
        return DnDsResult(dn or 0.0, float("inf"), None, n_codons, S, N, Sd, Nd,
                          saturated=True, notes=notes)
    if dn is None:
        notes.append("dN beyond correctable range — alignment likely wrong")
        return DnDsResult(float("inf"), ds, None, n_codons, S, N, Sd, Nd,
                          saturated=saturated, notes=notes)
    if saturated:
        notes.append("dS > 1.5 — synonymous sites near saturation")
    ratio = None if ds < 1e-6 else dn / ds
    if ratio is None:
        notes.append("dS ≈ 0 — sequences too similar to measure selection")
    return DnDsResult(dn, ds, ratio, n_codons, S, N, Sd, Nd,
                      saturated=saturated, notes=notes)


def translate(cds: str) -> str:
    prot = []
    for i in range(0, len(cds) - 2, 3):
        aa = GENETIC_CODE.get(cds[i:i + 3].upper(), "X")
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)


# ---- Ensembl data access -------------------------------------------------

def _get_json(path: str, **params) -> Optional[dict]:
    params.setdefault("content-type", "application/json")
    r = requests.get(f"{ENSEMBL_REST}{path}", params=params, timeout=30)
    time.sleep(_PAUSE_S)
    if r.status_code in (400, 404):
        return None
    r.raise_for_status()
    return r.json()


def fetch_cds(ensembl_id: str) -> tuple[str, str, str]:
    """CDS for an Ensembl protein/transcript/gene id.
    Returns (cds, species_slug, gene_or_transcript_id); raises ValueError
    when the id can't be resolved."""
    info = _get_json(f"/lookup/id/{ensembl_id}")
    if not info:
        raise ValueError(f"Ensembl id not found: {ensembl_id}")
    kind = info.get("object_type", "")
    species = info.get("species", "")
    if kind == "Translation":
        transcript = info.get("Parent", "")
    elif kind == "Transcript":
        transcript = ensembl_id
    elif kind == "Gene":
        tr = _get_json(f"/lookup/id/{ensembl_id}", expand=1) or {}
        transcripts = tr.get("Transcript", [])
        canonical = [t for t in transcripts if t.get("is_canonical")]
        transcript = (canonical or transcripts or [{}])[0].get("id", "")
    else:
        raise ValueError(f"unsupported Ensembl object type: {kind}")
    if not transcript:
        raise ValueError(f"no transcript for {ensembl_id}")
    seq = _get_json(f"/sequence/id/{transcript}", type="cds")
    if not seq or not seq.get("seq"):
        raise ValueError(f"no CDS for transcript {transcript}")
    return seq["seq"], species, transcript


def candidate_partners(gene_id: str, species: str, k: int = 6) -> list[dict]:
    """Ranked comparison partners for a dN/dS test: the candidate's own
    Compara orthologs, one2one first, then by target %id. Percent-id
    ranking is unreliable for fragment gene models, so callers should try
    several partners and keep the lowest-dS result — dS itself is the
    cleanest closeness measure."""
    data = _get_json(f"/homology/id/{species}/{gene_id}",
                     type="orthologues", sequence="none")
    if not data or not data.get("data"):
        return []
    targets: list[tuple[int, float, dict]] = []
    for entry in data["data"]:
        for h in entry.get("homologies", []):
            t = h.get("target", {})
            if not t.get("protein_id"):
                continue
            one2one = 0 if h.get("type") == "ortholog_one2one" else 1
            targets.append((one2one, -(t.get("perc_id") or 0.0), t))
    targets.sort(key=lambda x: (x[0], x[1]))
    return [t for _, _, t in targets[:k]]


@dataclass
class SelectionReport:
    query_id: str
    partner_id: str
    partner_species: str
    result: DnDsResult
    n_partners_tried: int = 1
    also_tried: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"Selection test: {self.query_id}  vs  {self.partner_id} "
            f"({self.partner_species})",
            self.result.summary(),
        ]
        if self.also_tried:
            lines.append(f"  (kept lowest-dS partner of {self.n_partners_tried} tried: "
                         + ", ".join(self.also_tried) + ")")
        return "\n".join(lines)


def selection_test(ensembl_protein_or_gene_id: str,
                   partner_id: str = "",
                   max_partners: int = 5,
                   log=lambda _m: None) -> SelectionReport:
    """Full automated dN/dS test for one Ensembl id.

    Without `partner_id`, several of the candidate's own Compara orthologs
    are tried and the lowest-dS (i.e. closest, least-saturated) comparison
    is kept.
    """
    log(f"Fetching CDS for {ensembl_protein_or_gene_id}…")
    cds_q, species, _ = fetch_cds(ensembl_protein_or_gene_id)
    info = _get_json(f"/lookup/id/{ensembl_protein_or_gene_id}") or {}
    gene_id = ensembl_protein_or_gene_id
    if info.get("object_type") == "Translation":
        tr = _get_json(f"/lookup/id/{info.get('Parent','')}") or {}
        gene_id = tr.get("Parent", gene_id)

    if partner_id:
        targets = [{"id": partner_id, "protein_id": partner_id, "species": "?"}]
    else:
        log(f"Finding comparison orthologs of {gene_id} ({species})…")
        targets = candidate_partners(gene_id, species, k=max_partners)
        if not targets:
            raise ValueError(f"no Compara orthologs found for {gene_id}")

    best: Optional[SelectionReport] = None
    tried: list[str] = []
    for t in targets:
        pid = t["protein_id"]
        log(f"Trying partner {pid} ({t.get('species','?')})…")
        try:
            cds_p, _, _ = fetch_cds(pid)
        except (ValueError, requests.RequestException):
            continue
        result = ng86(cds_q, cds_p)
        tried.append(pid)
        rep = SelectionReport(
            query_id=ensembl_protein_or_gene_id, partner_id=pid,
            partner_species=str(t.get("species", "?")), result=result,
        )
        def _key(r: SelectionReport) -> tuple:
            # prefer resolvable, unsaturated, then lowest dS
            res = r.result
            return (res.ratio is None, res.saturated,
                    res.ds if math.isfinite(res.ds) else 1e9)
        if best is None or _key(rep) < _key(best):
            best = rep
        if not result.saturated and result.ratio is not None:
            break  # good partner found — no need to keep trying
    if best is None:
        raise ValueError(f"no partner CDS retrievable for {gene_id}")
    best.n_partners_tried = len(tried)
    best.also_tried = [p for p in tried if p != best.partner_id]
    return best
