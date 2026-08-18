"""Pairwise alignment + simple progressive MSA.

Why no MAFFT/MUSCLE shell-out by default: keeping the tool zero-install for
non-bioinformatics users. If MAFFT is on PATH and `use_mafft=True`, we use
it — much higher quality for large alignments. Otherwise we fall back to a
star-style progressive alignment seeded by the longest sequence, which is
adequate for downstream distance / tree work on closely-related proteins.

For Foldseek-style remote homology or thousand-sequence problems, export
FASTA and run MAFFT or MMseqs2 externally — the analysis module isn't
trying to replace a dedicated MSA program.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Iterable

from Bio import Align
from Bio.Align import substitution_matrices


_BLOSUM62 = substitution_matrices.load("BLOSUM62")


def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.substitution_matrix = _BLOSUM62
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.mode = "global"
    return a


@dataclass
class AlignmentRow:
    label: str
    aligned: str   # equal length across all rows in an MSA


def pairwise_align(seq_a: str, seq_b: str) -> tuple[str, str, float]:
    """Return (aligned_a, aligned_b, score)."""
    if not seq_a or not seq_b:
        return seq_a, seq_b, 0.0
    aligner = _aligner()
    alignments = aligner.align(seq_a, seq_b)
    best = next(iter(alignments))
    # Biopython Alignment supports row indexing: best[0] is aligned seq A
    # with gap characters inserted, best[1] is aligned seq B.
    return str(best[0]), str(best[1]), float(best.score)


def progressive_msa(
    labelled_sequences: Iterable[tuple[str, str]],
    use_mafft: bool = False,
) -> list[AlignmentRow]:
    """Build an MSA from (label, sequence) pairs.

    With `use_mafft=True` and `mafft` on PATH, shell out for a real MSA.
    Otherwise: star alignment — pick the longest sequence as the center,
    align every other sequence to it, then pad gaps so every row has
    identical length.
    """
    pairs = [(l, s) for l, s in labelled_sequences if s]
    if not pairs:
        return []
    if len(pairs) == 1:
        return [AlignmentRow(pairs[0][0], pairs[0][1])]
    if use_mafft and shutil.which("mafft"):
        return _mafft(pairs)
    return _star_align(pairs)


def _star_align(pairs: list[tuple[str, str]]) -> list[AlignmentRow]:
    # Seed: longest sequence — most informative as anchor.
    pairs_sorted = sorted(pairs, key=lambda p: len(p[1]), reverse=True)
    center_label, center_seq = pairs_sorted[0]
    aligned_center = center_seq
    aligned_others: list[tuple[str, str]] = []  # (label, aligned-to-current-center)

    for label, seq in pairs_sorted[1:]:
        a_center, a_other, _ = pairwise_align(aligned_center.replace("-", ""), seq)
        # We need to merge the new center-side alignment with the existing
        # aligned_center, inserting gaps everywhere the new pairwise added
        # them on the center side. This is the classic merge step of star
        # alignment.
        aligned_center, aligned_others = _merge_into_msa(
            aligned_center, aligned_others, a_center, label, a_other,
        )

    rows = [AlignmentRow(center_label, aligned_center)]
    rows.extend(AlignmentRow(l, s) for l, s in aligned_others)
    # Pad shorter rows in case of any off-by-one.
    width = max(len(r.aligned) for r in rows)
    for r in rows:
        if len(r.aligned) < width:
            r.aligned = r.aligned + "-" * (width - len(r.aligned))
    return rows


def _merge_into_msa(
    msa_center: str,
    others: list[tuple[str, str]],
    new_center: str,
    new_label: str,
    new_other: str,
) -> tuple[str, list[tuple[str, str]]]:
    """Merge `new_center / new_other` (a pairwise) into the existing MSA.

    Walks both center alignments column-by-column. A gap in one center
    triggers a gap insertion in the other side's columns. The end result is
    an MSA where the center keeps its original residues, gaps are unioned,
    and `new_other` is added as a new row.
    """
    merged_center: list[str] = []
    new_row: list[str] = []
    new_others: list[list[str]] = [[] for _ in others]
    i = j = 0
    while i < len(msa_center) or j < len(new_center):
        # Past-end on either side: insert gaps from that side and advance
        # the other. Without this, both ".elif c_old == '-':" and
        # ".else: c_new == '-':" branches could fail to advance when one
        # index was already past its sequence's end → infinite loop.
        if i >= len(msa_center):
            merged_center.append("-")
            new_row.append(new_other[j] if j < len(new_other) else "-")
            for k in range(len(others)):
                new_others[k].append("-")
            j += 1
            continue
        if j >= len(new_center):
            merged_center.append(msa_center[i])
            new_row.append("-")
            for k, (_lbl, seq) in enumerate(others):
                new_others[k].append(seq[i] if i < len(seq) else "-")
            i += 1
            continue

        c_old = msa_center[i]
        c_new = new_center[j]
        if c_old == c_new:
            merged_center.append(c_old)
            new_row.append(new_other[j] if j < len(new_other) else "-")
            for k, (_lbl, seq) in enumerate(others):
                new_others[k].append(seq[i] if i < len(seq) else "-")
            i += 1; j += 1
        elif c_old == "-":
            merged_center.append("-")
            new_row.append(new_other[j] if j < len(new_other) else "-")
            for k, (_lbl, seq) in enumerate(others):
                new_others[k].append(seq[i] if i < len(seq) else "-")
            i += 1
        else:  # c_new == "-"
            merged_center.append("-")
            new_row.append(new_other[j] if j < len(new_other) else "-")
            for k in range(len(others)):
                new_others[k].append("-")
            j += 1
    out_others = [(label, "".join(seq)) for (label, _), seq in zip(others, new_others)]
    out_others.append((new_label, "".join(new_row)))
    return "".join(merged_center), out_others


def _mafft(pairs: list[tuple[str, str]]) -> list[AlignmentRow]:
    fasta = "".join(f">{label}\n{seq}\n" for label, seq in pairs)
    proc = subprocess.run(
        ["mafft", "--auto", "--quiet", "/dev/stdin"],
        input=fasta, text=True, capture_output=True, check=False,
    )
    if proc.returncode != 0:
        return _star_align(pairs)
    rows: list[AlignmentRow] = []
    cur_label, cur_parts = None, []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if cur_label is not None:
                rows.append(AlignmentRow(cur_label, "".join(cur_parts)))
            cur_label = line[1:].strip()
            cur_parts = []
        else:
            cur_parts.append(line.strip())
    if cur_label is not None:
        rows.append(AlignmentRow(cur_label, "".join(cur_parts)))
    return rows
