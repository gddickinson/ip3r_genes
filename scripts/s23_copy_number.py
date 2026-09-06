"""s23_copy_number.py — what counts as one gene, split out of `s23_classify`.

Copy number is S23's deliverable, so the two rules that decide how many genes
a set of alignments represents are a result rather than an implementation
detail, and they pull in opposite directions:

  `merge_split_loci`  folds two neighbouring same-strand clusters whose bait
                      spans are *complementary* — one gene the aligner broke,
                      not two genes that happen to be adjacent.
  `full_copies`       splits one cluster into the distinct complete gene
                      models inside it — because a cluster is not a gene
                      (D28), and with `-G` at 650 kb for the metazoa a cluster
                      is a large object: *Drosophila*'s 22 kb *Itpr* sits in a
                      297,487 bp one.

Both record what they did, so a copy number can be checked rather than taken.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s23_calibration as cal                                   # noqa: E402
from s5_sweep_lib import COV_FOUND, Aln, Locus                  # noqa: E402

#: Locus grades by bait coverage. `COV_FOUND` (0.70) is S5's own "found" bar,
#: reused so a full locus here and a found cell there mean the same thing.
COV_FULL = COV_FOUND
#: Below this a locus is a `scrap` — the shared channel module matching
#: something, not a gene model. Set at the coverage S5b measured its junk
#: population topping out at (30 %), one level of evidence below `fragment`.
COV_SCRAP = 0.30

#: Merge rule. Two same-strand loci within this distance are candidates for
#: being one split gene. Ten times `s5_sweep_lib.LOCUS_GAP`, because the
#: clustering gap is what already failed to join them.
MERGE_GAP_BP = 100_000
#: ...and they are only merged if their bait query spans overlap by no more
#: than this fraction of the shorter one. Complementary halves of a bait are
#: one gene; two loci each covering the same half are two genes.
MERGE_MAX_QUERY_OVERLAP = 0.20


def grade(coverage: float) -> str:
    if coverage >= COV_FULL:
        return "full"
    if coverage >= COV_SCRAP:
        return "fragment"
    return "scrap"


# ------------------------------------------------------------------ merging

def _query_overlap_frac(a, b) -> float:
    lo = max(a.q_start, b.q_start)
    hi = min(a.q_end, b.q_end)
    if hi < lo:
        return 0.0
    shorter = min(a.q_end - a.q_start + 1, b.q_end - b.q_start + 1)
    return (hi - lo + 1) / shorter if shorter > 0 else 1.0


def merge_split_loci(loci: list[Locus], gap: int = MERGE_GAP_BP,
                     max_overlap: float = MERGE_MAX_QUERY_OVERLAP
                     ) -> tuple[list[Locus], list[dict]]:
    """Fold neighbouring loci that look like one split gene.

    Returns (merged loci, merge records). A merge record names both loci and
    the two numbers the decision was made on, so the ledger can show its
    working — the alternative is a copy number nobody can check.
    """
    by_key: dict[tuple, list[Locus]] = {}
    for L in loci:
        by_key.setdefault((L.contig, L.strand), []).append(L)
    out: list[Locus] = []
    merges: list[dict] = []
    for (contig, strand), group in by_key.items():
        group.sort(key=lambda L: L.start)
        cur = None
        for L in group:
            if cur is None:
                cur = Locus(contig, strand, L.start, L.end, list(L.alns))
                continue
            dist = L.start - cur.end
            ov = _query_overlap_frac(cur.best, L.best)
            if 0 <= dist <= gap and ov <= max_overlap:
                merges.append({
                    "contig": contig, "strand": strand,
                    "kept_start": cur.start, "kept_end": cur.end,
                    "merged_start": L.start, "merged_end": L.end,
                    "gap_bp": dist, "query_overlap": round(ov, 3),
                    "kept_bait": cur.best.bait, "merged_bait": L.best.bait,
                    "why": (f"{dist:,} bp apart on one strand and their bait "
                            f"spans overlap {ov:.0%} — complementary parts of "
                            "one gene, counted once")})
                cur.end = max(cur.end, L.end)
                cur.alns.extend(L.alns)
            else:
                out.append(cur)
                cur = Locus(contig, strand, L.start, L.end, list(L.alns))
        if cur is not None:
            out.append(cur)
    return out, merges


# ------------------------------------------------------------ copy counting

def cds_footprint_bp(a: Aln) -> int:
    """The genomic bases the alignment's exons actually occupy."""
    return sum(e - s + 1 for s, e in (a.cds_blocks or [])) or (a.end - a.start + 1)


def _overlap_frac(a: Aln, b: Aln) -> float:
    ov = min(a.end, b.end) - max(a.start, b.start) + 1
    if ov <= 0:
        return 0.0
    shorter = min(a.end - a.start + 1, b.end - b.start + 1)
    return ov / shorter if shorter > 0 else 1.0


def full_copies(loci: list[Locus], family: str,
                cov_full: float = None,
                max_overlap: float = None) -> list[Aln]:
    """The distinct complete gene models the sweep found, genome-wide.

    **Why this is not a count of loci.** A locus is a cluster of alignments
    chained at `LOCUS_GAP`, and outside the vertebrates `-G` is 650 kb for the
    metazoa, so a cluster is a large object: in the S23a pilot *Drosophila*'s
    single 22 kb Itpr gene sat inside a locus spanning 297 kb — 13x the gene —
    because 26 alignments from other baits chained across it. The status call
    was right (the best alignment covered the gene's CDS completely), but two
    real genes inside one such chain would have been counted **once**, and
    copy number is what this task delivers.

    So a copy is defined on the evidence rather than on the cluster: a genomic
    interval carrying an alignment that covers a complete bait. Two alignments
    that overlap are the same copy — several baits hitting one gene is the
    normal case, not two genes — and two that do not are two copies, however
    the clustering happened to group them. Greedy by score, so the copy is
    anchored on the best evidence for it.
    """
    if cov_full is None:
        cov_full = COV_FULL
    if max_overlap is None:
        max_overlap = cal.copy_max_overlap()[0]
    cands = [a for L in loci for a in L.alns
             if a.family == family and a.coverage >= cov_full]
    cands.sort(key=lambda a: (-a.score, a.contig, a.start))
    kept: list[Aln] = []
    for a in cands:
        if any(a.contig == b.contig and _overlap_frac(a, b) > max_overlap
               for b in kept):
            continue
        kept.append(a)
    return sorted(kept, key=lambda a: (a.contig, a.start))
