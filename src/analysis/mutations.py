"""Mutation calling vs a chosen reference sequence.

Pairwise-aligns each query sequence to a reference, then enumerates
column-by-column differences. The result is a list of `MutationCall` —
substitutions and indels with reference-coordinate positions, suitable for
spotting:

    * Known disease-causing ITPR1 mutations (the SCA15/SCA29 and Gillespie
      syndrome alleles), verified against ClinVar by the constraint task
      in newly fetched orthologs.
    * Lineage-specific changes — e.g. residues conserved in ITPR1/ITPR3
      but altered in ITPR2.
    * Candidate functional sites where all sequences differ from the
      reference.
"""

from __future__ import annotations

from dataclasses import dataclass

from .alignment import pairwise_align


@dataclass
class MutationCall:
    reference_label: str
    query_label: str
    ref_pos: int             # 1-based position in the reference
    ref_aa: str              # '-' for insertion in query relative to ref
    query_aa: str            # '-' for deletion in query relative to ref
    kind: str                # 'substitution' | 'insertion' | 'deletion'

    def short(self) -> str:
        if self.kind == "substitution":
            return f"{self.ref_aa}{self.ref_pos}{self.query_aa}"
        if self.kind == "deletion":
            return f"{self.ref_aa}{self.ref_pos}del"
        return f"{self.ref_pos}ins{self.query_aa}"


def call_mutations(
    reference_label: str, reference_seq: str,
    query_label: str, query_seq: str,
) -> list[MutationCall]:
    if not reference_seq or not query_seq:
        return []
    a_ref, a_qry, _ = pairwise_align(reference_seq, query_seq)
    out: list[MutationCall] = []
    ref_pos = 0
    for ra, qa in zip(a_ref, a_qry):
        if ra != "-":
            ref_pos += 1
        if ra == qa:
            continue
        if ra == "-":
            out.append(MutationCall(reference_label, query_label, ref_pos, ra, qa, "insertion"))
        elif qa == "-":
            out.append(MutationCall(reference_label, query_label, ref_pos, ra, qa, "deletion"))
        else:
            out.append(MutationCall(reference_label, query_label, ref_pos, ra, qa, "substitution"))
    return out
