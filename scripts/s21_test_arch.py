"""Negative controls for the S21 rules, run before anything is written.

The pattern `s5_bait_screen.self_test()` established, applied to a task whose
failure modes are unusually quiet. Every rule here returns a plausible number
when it is wrong, and four of them return a *good-looking* one:

* a reader that keys models the way miniprot names them silently merges the
  blocks of two different genes in the two chunked assemblies, and the merged
  gene has more exons and longer introns;
* computing an intron as `next.start - prev.end` reverses every minus-strand
  gene, which yields negative lengths on some models and plausible ones on
  others;
* inverting the intron-phase convention swaps phases 1 and 2 everywhere,
  changing every shared-intron call without changing a single count; and
* asking whether the same bait aligns twice at disjoint positions measures
  paralogy, not duplication, in a family whose members are 61-68 % identical —
  and returns a specificity of 0.16 while looking like a working detector.

So these are tests on **refusal** and on **reachability**: each rule must fire
on its own violation, and each measurement that came out zero must be shown able
to come out otherwise. `within_locus` duplications and merged frameshift pairs
are both zero over the whole sweep, and a rule that cannot act is not a result.

T1-T9 and T21 are here; T10-T20 — the shared-intron test, the duplication
detector and the annotation half — are in `s21_test_claims.py`, which imports
this module's helpers so the two halves cannot record a failure differently.

Mutation-tested on six deliberate rule breakages — the phase complement, the
intron-phase convention, the strand branch in the gap, the phase requirement in
the shared-intron match, the cell attribution in the duplication detector, and
the reader's collision key — all six caught, by T6, T7, T3, T13, T16 and T1
respectively.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_architecture as A                                      # noqa: E402
import s21_blocks as B                                            # noqa: E402
import s21_frame as FR                                            # noqa: E402
import s21_lib as L                                               # noqa: E402
import s5_sweep_lib as SW                                         # noqa: E402

FAILURES: list[str] = []
SKIPPED: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    print(f"  {'ok  ' if ok else 'FAIL'} {name}"
          + (f" — {detail}" if detail and not ok else ""))


def skip(name: str, why: str) -> None:
    SKIPPED.append(f"{name}: {why}")
    print(f"  skip {name} — {why}")


# --------------------------------------------------------------------------
# constructed inputs
# --------------------------------------------------------------------------
def _blocks(spec: list[tuple[int, int, int, int, int]]) -> list[dict]:
    return [{"start": s, "end": e, "phase": p, "q_start": qs, "q_end": qe,
             "identity": 1.0} for s, e, p, qs, qe in spec]


def _model(strand: str, spec) -> dict:
    bl = _blocks(spec)
    return {"mp_id": "T", "contig": "c", "strand": strand, "bait": "B",
            "identity": 1.0, "frameshifts": 0, "stop_codons": 0,
            "score": 1.0, "blocks": bl,
            "start": min(b["start"] for b in bl),
            "end": max(b["end"] for b in bl)}


class _Region:
    """A LocusRegion stand-in over a literal sequence, 1-based at `start`."""

    def __init__(self, seq: str, start: int = 1):
        self.seq, self.start = seq.upper(), start
        self.end = start + len(seq) - 1
        self.contig = "c"

    def at(self, s: int, e: int) -> str:
        if s < self.start or e > self.end or e < s:
            return ""
        return self.seq[s - self.start:e - self.start + 1]


class _Aln:
    def __init__(self, contig, start, end, strand, bait, q_start, q_end,
                 bait_len=1000, score=1000.0, identity=0.9, clade="ITPR1",
                 family="ITPR"):
        self.contig, self.start, self.end, self.strand = contig, start, end, strand
        self.bait, self.q_start, self.q_end = bait, q_start, q_end
        self.bait_len, self.score, self.identity = bait_len, score, identity
        self.clade, self.family = clade, family
        self.aligned_aa = q_end - q_start + 1
        self.cds_blocks = []
        self.mp_id = f"{contig}:{start}"
        self.frameshifts = 0
        self.stop_codons = 0
        self.translation = ""

    @property
    def coverage(self):
        return min(1.0, self.aligned_aa / self.bait_len)

    @property
    def span_coverage(self):
        return self.coverage


# --------------------------------------------------------------------------
# T1-T2 — the reader
# --------------------------------------------------------------------------
def t1_reader_agrees_with_the_sweeps_own() -> None:
    """The two readers must return the same genomic blocks, chunked genome too.

    `s5_sweep_lib.parse_miniprot_gff` produced the sweep's every number and the
    keys the summaries record; S21 reads the same file for the phases and query
    spans it does not keep. A disagreement anywhere would move exons.
    """
    meta = {r["id"]: {"clade": r["clade"], "family": r["family"],
                      "length": int(float(r["length"]))}
            for r in L.read_tsv(L.BAIT_MANIFEST)}
    for acc, label in (("GCF_000001405.40", "unchunked"),
                       ("GCF_019279795.1", "25-chunk"),
                       ("GCF_964261635.1", "14-chunk")):
        gff = L.sweep_dir(acc) / "miniprot.gff"
        if not gff.exists():
            skip(f"T1 reader agreement ({label})", "sweep GFF not on disk")
            continue
        theirs = {a.mp_id: sorted(a.cds_blocks)
                  for a in SW.parse_miniprot_gff(gff, meta)}
        mine = {k: sorted((b["start"], b["end"]) for b in m["blocks"])
                for k, m in B.read_models(gff).items()}
        check(f"T1 reader agreement ({label})", theirs == mine,
              f"{sum(1 for k in theirs if theirs.get(k) != mine.get(k))} of "
              f"{len(theirs)} models differ; "
              f"{len(set(theirs) ^ set(mine))} keys differ")


def t2_every_summary_mp_id_is_a_key() -> None:
    """The keys the sweep recorded must all be reachable in the GFF."""
    summaries = L.load_summaries()
    missing = 0
    total = 0
    for acc in ("GCF_019279795.1", "GCF_964261635.1", "GCF_027579735.1"):
        s = summaries.get(acc)
        gff = L.sweep_dir(acc) / "miniprot.gff"
        if not s or not gff.exists():
            skip(f"T2 summary keys ({acc})", "summary or GFF not on disk")
            continue
        want = {loc["mp_id"] for c in (s.get("cells") or {}).values()
                for loc in c.get("loci") or []}
        want |= {loc["mp_id"] for loc in s.get("other_loci") or []}
        keys = set(B.read_models(gff))
        total += len(want)
        missing += len(want - keys)
    if total:
        check("T2 every summary mp_id is a reader key", missing == 0,
              f"{missing} of {total} recorded ids are not keys")


# --------------------------------------------------------------------------
# T3-T7 — gene order, splicing, phase
# --------------------------------------------------------------------------
def t3_minus_strand_runs_in_gene_order() -> None:
    m = _model("-", [(2000, 2099, 0, 1, 33), (1000, 1099, 1, 34, 67)])
    arch = B.architecture(m, 10)
    naive = m["blocks"][1]["start"] - m["blocks"][0]["end"] - 1
    check("T3 minus-strand intron runs in gene order",
          arch["n_introns"] == 1 and arch["introns"][0]["length"] == 900
          and naive < 0,
          f"length {arch['introns'][0]['length'] if arch['introns'] else None}, "
          f"naive {naive}")


def t4_minus_strand_splice_pair_reads_gt_ag() -> None:
    """A canonical minus-strand intron must read GT..AG, never CT..AC."""
    # genomic: ... exon | CT ........ AC | exon ...  (revcomp of GT..AG)
    seq = "A" * 10 + "CT" + "N" * 16 + "AC" + "A" * 10
    region = _Region(seq, start=1)
    donor, acceptor = L.splice_pair(region, (11, 30), "-")
    plus = L.splice_pair(region, (11, 30), "+")
    check("T4 minus-strand splice pair is GT..AG",
          (donor, acceptor) == ("GT", "AG") and plus == ("CT", "AC"),
          f"minus {donor}..{acceptor}, plus {plus[0]}..{plus[1]}")


def t5_merge_fires_on_a_frameshift_pair() -> None:
    """The merge must act on a sub-intron gap and leave a real intron alone.

    `n_merges` is 0 over the whole sweep. That is only a result if the rule can
    act, so it is made to act here: two blocks 2 bp apart with contiguous query
    spans are one exon; the same pair 5 kb apart is two.
    """
    close = _model("+", [(1000, 1099, 0, 1, 33), (1102, 1201, 1, 34, 67)])
    far = _model("+", [(1000, 1099, 0, 1, 33), (6000, 6099, 1, 34, 67)])
    a_close = B.architecture(close, 30)
    a_far = B.architecture(far, 30)
    check("T5 merge fires on a 2 bp gap and not on a 4.9 kb one",
          a_close["n_exons"] == 1 and a_close["n_merges"] == 1
          and a_far["n_exons"] == 2 and a_far["n_merges"] == 0,
          f"close {a_close['n_exons']}/{a_close['n_merges']}, "
          f"far {a_far['n_exons']}/{a_far['n_merges']}")


def t6_frame_step_is_detected_and_absent() -> None:
    """A constructed frame step must be found; a clean chain must give zero."""
    clean = _blocks([(1000, 1099, 0, 1, 33), (2000, 2099, 2, 34, 67),
                     (3000, 3100, 1, 68, 101)])
    broken = _blocks([(1000, 1099, 0, 1, 33), (2000, 2099, 0, 34, 67)])
    check("T6 frame steps counted, and zero on a clean chain",
          B.frame_steps(clean) == 0 and B.frame_steps(broken) == 1,
          f"clean {B.frame_steps(clean)}, broken {B.frame_steps(broken)}")


def t7_intron_phase_is_not_inverted() -> None:
    """GFF phase 1 is intron phase 2 and phase 2 is intron phase 1.

    Inverting this swaps phases 1 and 2 in every row and changes every
    shared-intron call while changing no count, which is why it is a test and
    not a comment.
    """
    got = (A.intron_phase(0), A.intron_phase(1), A.intron_phase(2))
    check("T7 intron phase convention", got == (0, 2, 1), f"got {got}")


# --------------------------------------------------------------------------
# T8-T9 — the frame and the scope
# --------------------------------------------------------------------------
def t8_anchor_test_can_fail() -> None:
    """A frame slipped by one residue must be refused, not transferred.

    The real anchors pass, which on its own says nothing: the test has to be
    shown able to reject. One paralogue's aligned row is shifted by inserting a
    gap, which moves every column it contributes, and `require_anchors` must
    then refuse.
    """
    real = FR.anchor_test()
    ok_real = all(r["one_column"] for r in real) and len(real) >= 10
    rows = L.msa_rows()
    saved = L.msa_rows
    shifted = dict(rows)
    acc = L.FRAME_ACC["ITPR2"]
    shifted[acc] = "-" + shifted[acc][:-1] if shifted[acc][-1] == "-" \
        else "-" + shifted[acc]
    L.msa_rows = lambda: shifted
    try:
        bad = FR.anchor_test()
        refused = False
        try:
            FR.require_anchors()
        except SystemExit:
            refused = True
    finally:
        L.msa_rows = saved
    n_bad = sum(1 for r in bad if not r["one_column"])
    check("T8 anchor test passes on the real frame and refuses a shifted one",
          ok_real and refused and n_bad > 0,
          f"real all-one-column {ok_real}, shifted bad {n_bad}, "
          f"refused {refused}")


def t9_each_scope_rule_fires_on_its_own_violation() -> None:
    base = {"cell": "ITPR1", "is_copy": 1, "contig_spans_gene": 1,
            "coverage": 0.99, "bait": "Q14643"}
    frames_ok = {("Q14643", "ITPR1")}
    cases = [
        ({}, ""),
        ({"cell": "vertebrate_basal"}, "A0_cell_vertebrate_basal"),
        ({"is_copy": 0}, "A1_not_a_gene_copy"),
        ({"contig_spans_gene": 0}, "A2_contig_cannot_span_the_gene"),
        ({"coverage": 0.5}, f"A3_coverage_below_{A.COV_ARCH}"),
        ({"bait": "OTHER"}, "A4_no_usable_frame"),
    ]
    bad = []
    for patch, want in cases:
        got = A.scope_exclusion({**base, **patch}, frames_ok)
        if got != want:
            bad.append(f"{patch or 'clean'} -> {got!r} (wanted {want!r})")
    check("T9 every architecture scope rule fires on its own violation",
          not bad, "; ".join(bad))


def _committed_state() -> dict:
    out = {}
    if L.OUT.exists():
        for p in sorted(L.OUT.rglob("*")):
            if p.is_file():
                out[str(p)] = p.stat().st_mtime_ns
    return out


def t21_suite_writes_nothing(before: dict) -> None:
    after = _committed_state()
    changed = [k for k in set(before) | set(after)
               if before.get(k) != after.get(k)]
    check("T21 the suite alters no committed table", not changed,
          f"{len(changed)} files changed: {changed[:3]}")


def main() -> int:
    print("[s21] negative controls")
    before = _committed_state()
    t1_reader_agrees_with_the_sweeps_own()
    t2_every_summary_mp_id_is_a_key()
    t3_minus_strand_runs_in_gene_order()
    t4_minus_strand_splice_pair_reads_gt_ag()
    t5_merge_fires_on_a_frameshift_pair()
    t6_frame_step_is_detected_and_absent()
    t7_intron_phase_is_not_inverted()
    t8_anchor_test_can_fail()
    t9_each_scope_rule_fires_on_its_own_violation()
    import s21_test_claims
    s21_test_claims.run()
    t21_suite_writes_nothing(before)
    if FAILURES:
        print(f"[s21] {len(FAILURES)} FAILURES")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    if SKIPPED:
        print(f"[s21] all negative controls pass; {len(SKIPPED)} skipped")
        for s in SKIPPED:
            print(f"   - {s}")
        return 0
    print("[s21] all negative controls pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
