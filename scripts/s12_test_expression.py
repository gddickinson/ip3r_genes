"""Negative controls for the S12 rules, run before anything is written.

The pattern is `s5_bait_screen.self_test()`'s and `s10_test_evidence.py`'s:
every rule here returns a *plausible* number when it is wrong, so each is
tested against a case constructed to break it.  Three of these were written
after the rule they test had already failed on real data during this
session, and they are marked where that is so.

Run:  python scripts/s12_test_expression.py [-v]
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s12_locus as SL  # noqa: E402
import s12_quantify as Q  # noqa: E402
from s12_lib import normalise_tissue, revcomp, translate  # noqa: E402


class _Region:
    """A synthetic genome window with the same interface as `s10_evidence.Region`."""

    def __init__(self, seq: str, start: int = 1, contig: str = "test"):
        self.seq, self.start, self.contig = seq, start, contig
        self.end = start + len(seq) - 1

    def at(self, start: int, end: int) -> str:
        if start < self.start or end > self.end or end < start:
            return ""
        return self.seq[start - self.start:end - self.start + 1]


def _model(blocks, strand="+", frameshifts=0):
    return {"mp_id": "TEST", "contig": "test", "strand": strand,
            "start": min(b[0] for b in blocks), "end": max(b[1] for b in blocks),
            "frameshifts": frameshifts, "stop_codons": 0, "bait": "b",
            "cds": [{"start": s, "end": e, "phase": p, "identity": 1.0,
                     "q_start": i * 100 + 1, "q_end": (i + 1) * 100}
                    for i, (s, e, p) in enumerate(blocks)]}


def _two_exon_gene():
    """A 2-exon gene on a 600 nt window: exon1 1-90, intron, exon2 201-290.

    The window runs to 600 nt and carries a third exon-sized block at
    401-490 so the minus-strand, two-intron case in T5 has genome to sit on.
    """
    import random
    random.seed(12)
    ex1 = "ATG" + "".join(random.choice("ACGT") for _ in range(87))
    ex2 = "".join(random.choice("ACGT") for _ in range(90))
    intron = "GT" + "".join(random.choice("ACGT") for _ in range(106)) + "AG"
    filler = "".join(random.choice("ACGT") for _ in range(310))
    seq = ex1 + intron + ex2 + filler
    region = _Region(seq)
    model = _model([(1, 90, 0), (201, 290, 0)])
    prot = translate(ex1) + translate(ex2)
    return region, model, prot, ex1 + ex2


CHECKS: list[tuple[str, callable]] = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# --- T1-T4: the reference construction --------------------------------------

@check("T1 splice concatenates blocks in target order, not coordinate order")
def t1():
    region, model, _, spliced = _two_exon_gene()
    seq, rows = SL.splice(model, region)
    assert seq == spliced, "plus-strand splice wrong"
    assert [r["cds_start"] for r in rows] == [1, 91]
    # minus strand: same blocks, reversed target order, revcomp'd
    rev = copy.deepcopy(model)
    rev["strand"] = "-"
    rev["cds"] = list(reversed(rev["cds"]))
    seq2, _ = SL.splice(rev, region)
    assert seq2 == revcomp(region.at(201, 290)) + revcomp(region.at(1, 90)), \
        "minus-strand splice did not reverse-complement in target order"


@check("T2 a block outside the fetched window raises rather than shortening")
def t2():
    region, model, _, _ = _two_exon_gene()
    bad = copy.deepcopy(model)
    bad["cds"].append({"start": 900, "end": 990, "phase": 0, "identity": 1.0,
                       "q_start": 201, "q_end": 300})
    try:
        SL.splice(bad, region)
    except SL.LocusError:
        return
    raise AssertionError("a block past the window silently shortened the "
                         "reference")


@check("T3 colinear placement rejects reversed order and the wrong strand")
def t3():
    region, model, prot, _ = _two_exon_gene()
    placed, tested = SL.colinear_placement(model, region, prot, minlen=10)
    assert (placed, tested) == (2, 2), f"correct model placed {placed}/{tested}"
    rev = copy.deepcopy(model)
    rev["cds"] = list(reversed(rev["cds"]))
    p, t = SL.colinear_placement(rev, region, prot, minlen=10)
    assert p < t, f"reversed block order still placed {p}/{t}"
    wrong = copy.deepcopy(model)
    wrong["strand"] = "-"
    p, t = SL.colinear_placement(wrong, region, prot, minlen=10)
    assert p == 0, f"wrong strand placed {p}/{t}"


@check("T4 the frame-aware translation uses phase, not frame 0")
def t4():
    # exon1 ends mid-codon: 89 nt, so exon2 carries phase 1.
    import random
    random.seed(7)
    cds = "ATG" + "".join(random.choice("ACGT") for _ in range(177))
    region = _Region(cds[:89] + "GT" + "N" * 96 + "AG" + cds[89:] + "A" * 50)
    model = _model([(1, 89, 0), (190, 290, 1)])
    obs, breaks = SL.frame_aware_translation(model, region)
    assert not breaks, f"a continuous frame reported breaks {breaks}"
    assert obs.startswith(translate(cds)[:25]), \
        "phase ignored: the second exon was translated in the wrong frame"


# --- T5-T7: junction tagging ------------------------------------------------

@check("T5 junction offsets and classes are matched by coordinate, not index")
def t5():
    region, model, _, _ = _two_exon_gene()
    _, rows = SL.splice(model, region)
    introns = [{"start": 91, "end": 200, "class": "between_models",
                "donor": "GT", "acceptor": "AG", "splice_class": "canonical",
                "left_models": "gA", "right_models": "gB"}]
    js = SL.junctions(model, rows, introns)
    assert len(js) == 1 and js[0]["cds_offset"] == 90
    assert js[0]["class"] == "between_models" and js[0]["annotated"] == 0

    # Minus strand with **two** introns, which is the smallest case that can
    # tell a coordinate lookup from an index lookup: in target order the
    # first junction is the genomically *last* intron, so an index-to-index
    # mapping tags each junction with the other one's class and a one-intron
    # test cannot see it.  This case was added after a mutation test found
    # exactly that mapping slipping through.
    three = _model([(1, 90, 0), (201, 290, 0), (401, 490, 0)], strand="-")
    three["cds"] = list(reversed(three["cds"]))
    _, rrows = SL.splice(three, region)
    two_introns = [
        {"start": 91, "end": 200, "class": "unannotated", "donor": "GT",
         "acceptor": "AG", "splice_class": "canonical", "left_models": "",
         "right_models": ""},
        {"start": 291, "end": 400, "class": "within_one_model", "donor": "GT",
         "acceptor": "AG", "splice_class": "canonical", "left_models": "gA",
         "right_models": "gA"}]
    js = SL.junctions(three, rrows, two_introns)
    assert len(js) == 2
    assert (js[0]["intron_start"], js[0]["intron_end"]) == (291, 400), \
        "minus-strand junction 1 did not resolve to the genomically last intron"
    assert js[0]["class"] == "within_one_model" and js[0]["annotated"] == 1, \
        "junction 1 took its class from the wrong intron"
    assert (js[1]["intron_start"], js[1]["intron_end"]) == (91, 200)
    assert js[1]["class"] == "unannotated" and js[1]["annotated"] == 0, \
        "junction 2 took its class from the wrong intron"


@check("T6 only within_one_model counts as an annotated junction")
def t6():
    region, model, _, _ = _two_exon_gene()
    _, rows = SL.splice(model, region)
    for klass, want in (("within_one_model", 1), ("between_models", 0),
                        ("model_to_gap", 0), ("unannotated", 0)):
        js = SL.junctions(model, rows, [{
            "start": 91, "end": 200, "class": klass, "donor": "GT",
            "acceptor": "AG", "splice_class": "canonical",
            "left_models": "", "right_models": ""}])
        assert js[0]["annotated"] == want, f"{klass} scored {js[0]['annotated']}"


@check("T7 a junction whose gap matches no intron is kept, not dropped")
def t7():
    region, model, _, _ = _two_exon_gene()
    _, rows = SL.splice(model, region)
    js = SL.junctions(model, rows, [])
    assert len(js) == 1 and js[0]["class"] == "non_colinear", \
        "an unmatched junction disappeared; a dropped junction is " \
        "indistinguishable from one nothing crossed"


# --- T8-T11: the read-counting rules ----------------------------------------

@check("T8 CIGAR blocks: soft clips and insertions consume no reference")
def t8():
    assert Q.aligned_blocks(101, "10S40M") == [(101, 140)]
    assert Q.aligned_blocks(101, "20M5I20M") == [(101, 120), (121, 140)]
    assert Q.aligned_blocks(101, "20M100N20M") == [(101, 120), (221, 240)]


@check("T9 a junction read needs the anchor on BOTH sides of ONE block")
def t9():
    # 8 nt anchor: crossing by 8 counts, by 7 does not.
    assert Q.spanned([(93, 108)], [100], 8) == [100]
    assert Q.spanned([(94, 108)], [100], 8) == []
    assert Q.spanned([(93, 107)], [100], 8) == []
    # two blocks each ending at the junction is not a span: this is the
    # failure mode S10 found, an alignment that stops at the splice point.
    assert Q.spanned([(80, 100), (101, 130)], [100], 8) == []


@check("T10 a read touching the junction from one side only is not a span")
def t10():
    assert Q.spanned([(60, 100)], [100], 8) == []
    assert Q.spanned([(100, 140)], [100], 8) == []


@check("T11 tissue vocabulary prefers the longest match, not dict order")
def t11():
    assert normalise_tissue("head kidney") == "kidney"
    assert normalise_tissue("swim bladder") == "swim_bladder"
    assert normalise_tissue("whole blood") == "blood"
    assert normalise_tissue("not a tissue at all") is None


def main() -> int:
    verbose = "-v" in sys.argv
    failed = []
    for name, fn in CHECKS:
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001
            failed.append((name, exc))
            print(f"FAIL {name}\n     {type(exc).__name__}: {exc}")
        else:
            if verbose:
                print(f"ok   {name}")
    print(f"\ns12 self-test: {len(CHECKS) - len(failed)}/{len(CHECKS)} passed")
    return 1 if failed else 0


def self_test() -> None:
    """Raise unless every check passes (called by the driver before writing)."""
    if main() != 0:
        raise SystemExit("s12 self-test failed — refusing to write")


if __name__ == "__main__":
    sys.exit(main())
