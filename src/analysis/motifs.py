"""Family-signature motifs — MSA-derived PSSM blocks, dependency-free.

InterPro/Pfam domain evidence only exists for proteins that made it into
UniProt. The most interesting candidates (unnamed Ensembl gene models,
fresh BLAST hits) aren't there — so the discovery scorer's domain
component was structurally blind to them.

This module closes the gap with a data-driven, offline signature scan:

    1. `derive_family_signatures()` — from the *known-paralog* rows of the
       analysis MSA, pick the K most-conserved gapless windows and turn
       each into a log-odds position-specific scoring matrix (PSSM with
       pseudocounts vs. background frequencies). These windows are the
       family's empirical sequence fingerprint (for IP3Rs they land in
       the pore/CTD and other conserved blades).
    2. `FamilySignatureSet.scan(seq)` — slide each PSSM over any raw
       sequence; a signature "hits" when its best window score exceeds
       `hit_fraction` of its own self-score. Coverage = hits / K.

A candidate carrying most of the family's conserved blocks is family-fold
positive even with zero database annotation — a poor man's profile-HMM
built from data already in hand.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Optional

# Rough UniProt/Swiss-Prot background amino-acid frequencies.
BACKGROUND = {
    "A": 0.083, "R": 0.055, "N": 0.041, "D": 0.055, "C": 0.014,
    "Q": 0.039, "E": 0.067, "G": 0.071, "H": 0.023, "I": 0.059,
    "L": 0.097, "K": 0.058, "M": 0.024, "F": 0.039, "P": 0.047,
    "S": 0.066, "T": 0.053, "W": 0.011, "Y": 0.029, "V": 0.069,
}
_PSEUDO = 0.5


@dataclass
class SignaturePSSM:
    start_col: int                  # MSA column where the window starts
    length: int
    matrix: list[dict[str, float]]  # per-position log-odds
    max_score: float                # score of the consensus (self) match
    consensus: str

    def best_window_score(self, seq: str) -> float:
        if len(seq) < self.length:
            return float("-inf")
        best = float("-inf")
        for i in range(len(seq) - self.length + 1):
            s = 0.0
            for j in range(self.length):
                s += self.matrix[j].get(seq[i + j], -2.0)
            if s > best:
                best = s
        return best


@dataclass
class FamilySignatureSet:
    signatures: list[SignaturePSSM]
    n_source_rows: int
    hit_fraction: float = 0.6       # best score ≥ this × max_score = hit

    def scan(self, seq: str) -> tuple[int, int, list[float]]:
        """Returns (n_hits, n_signatures, normalized best score per sig)."""
        norms: list[float] = []
        hits = 0
        for sig in self.signatures:
            best = sig.best_window_score(seq)
            norm = best / sig.max_score if sig.max_score > 0 else 0.0
            norms.append(round(norm, 3))
            if norm >= self.hit_fraction:
                hits += 1
        return hits, len(self.signatures), norms

    def coverage(self, seq: str) -> float:
        if not self.signatures or not seq:
            return 0.0
        hits, k, _ = self.scan(seq)
        return hits / k


def derive_family_signatures(
    aligned_rows: Iterable[tuple[str, str]],
    n_signatures: int = 8,
    window: int = 15,
    max_gap_fraction: float = 0.2,
) -> Optional[FamilySignatureSet]:
    """Build the signature set from (label, aligned_sequence) rows —
    normally the KNOWN-paralog rows of the analysis MSA. Returns None when
    fewer than 3 rows are available (a PSSM from 2 sequences is noise)."""
    rows = [(lbl, aln) for lbl, aln in aligned_rows if aln]
    if len(rows) < 3:
        return None
    length = min(len(a) for _, a in rows)
    n_rows = len(rows)

    # Per-column conservation (frequency of the modal residue, gaps count
    # against) and gap fraction.
    col_scores: list[float] = []
    for c in range(length):
        col = [a[c] for _, a in rows]
        gaps = sum(1 for x in col if x == "-")
        if gaps / n_rows > max_gap_fraction:
            col_scores.append(0.0)
            continue
        residues = [x for x in col if x != "-"]
        modal = max((residues.count(r) for r in set(residues)), default=0)
        col_scores.append(modal / n_rows)

    # Rank candidate windows by mean conservation; greedily keep
    # non-overlapping ones.
    window_scores = []
    for start in range(0, length - window + 1):
        window_scores.append((sum(col_scores[start:start + window]) / window, start))
    window_scores.sort(reverse=True)
    chosen: list[int] = []
    for score, start in window_scores:
        if score <= 0.0 or len(chosen) >= n_signatures:
            break
        if all(abs(start - c) >= window for c in chosen):
            chosen.append(start)
    if not chosen:
        return None

    sigs: list[SignaturePSSM] = []
    for start in sorted(chosen):
        matrix: list[dict[str, float]] = []
        consensus = []
        for c in range(start, start + window):
            col = [a[c] for _, a in rows if a[c] != "-"]
            scores: dict[str, float] = {}
            n = len(col) or 1
            for aa, bg in BACKGROUND.items():
                freq = (col.count(aa) + _PSEUDO * bg * 20) / (n + _PSEUDO * 20)
                scores[aa] = math.log2(freq / bg)
            matrix.append(scores)
            consensus.append(max(scores, key=scores.get))
        cons = "".join(consensus)
        max_score = sum(m[cons[j]] for j, m in enumerate(matrix))
        sigs.append(SignaturePSSM(
            start_col=start, length=window, matrix=matrix,
            max_score=max_score, consensus=cons,
        ))
    return FamilySignatureSet(signatures=sigs, n_source_rows=n_rows)


def signature_report(sset: FamilySignatureSet) -> str:
    lines = [
        f"Family signature set — {len(sset.signatures)} conserved block(s) "
        f"derived from {sset.n_source_rows} known-paralog sequences "
        f"(hit threshold {sset.hit_fraction:.0%} of self-score):",
    ]
    for i, s in enumerate(sset.signatures, 1):
        lines.append(f"  #{i}  MSA cols {s.start_col}-{s.start_col + s.length - 1}"
                     f"  consensus {s.consensus}")
    return "\n".join(lines)
