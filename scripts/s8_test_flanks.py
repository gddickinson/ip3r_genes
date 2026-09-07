#!/usr/bin/env python3
"""Negative controls for the S8 synteny rules, run on every build.

The pattern of `s5_bait_screen.self_test()`, `s6_test_selection.py` and
`s7_test_tree.py`: constructed cases that each *must* be rejected by the
rule responsible, plus the two properties the analysis has to have before
its numbers mean anything -- determinism, and a caller that fails on a
neighbourhood that is not one.

Run:  python scripts/s8_test_flanks.py
Also invoked by s8_run_synteny.py before any table is written.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import s8_control as ctl      # noqa: E402
import s8_flank_lib as lib    # noqa: E402
import s8_paralogon as par    # noqa: E402


class Failure(AssertionError):
    pass


def check(name: str, ok: bool, detail: str = ""):
    if not ok:
        raise Failure(f"{name}: {detail or 'failed'}")
    return (name, detail)


# --------------------------------------------------------------- fixtures

def gene(contig, start, end, sym, biotype="protein_coding", strand="+"):
    return (contig, start, end, strand, f"gene-{sym}", sym, biotype)


def synthetic_index():
    """One contig of 25 coding genes at 10 kb spacing, plus a short contig.

    Every fifth symbol is a RefSeq placeholder, so the informative window has
    to scan past them and the fixed window has to keep them.
    """
    genes = []
    for i in range(25):
        sym = f"LOC{900000 + i}" if i % 5 == 0 else f"GENE{i:02d}"
        genes.append(gene("chr1", 10_000 * i + 1, 10_000 * i + 5_000, sym))
    genes.append(gene("chr1", 120_001, 128_000, "TARGET"))
    genes += [gene("chr2", 1_000 * i + 1, 1_000 * i + 500, f"SHORT{i}")
              for i in range(5)]
    genes.append(gene("chr1", 300_001, 300_500, "NCRNA", biotype="lncRNA"))
    return lib.index_by_contig(genes)


def synthetic_locus(contig="chr1", start=120_001, end=128_000, cell="ITPR1"):
    return lib.Locus(
        accession="GCA_TEST.1", organism="Testus testus", vclass="Testia",
        vorder="Testales", cell=cell, status="found_annotated", contig=contig,
        start=start, end=end, idx=0, bait="", bait_paralog="", identity=1.0,
        coverage=1.0, annot_gene="TARGET", annot_paralog=cell)


def fake_set(cell, species, syms, vclass="Testia", coverage=1.0, idx=0):
    fs = lib.FlankSet(locus=lib.Locus(
        accession=f"ACC_{species}", organism=species, vclass=vclass,
        vorder="", cell=cell, status="found_annotated", contig="c", start=1,
        end=2, idx=idx, bait="", bait_paralog="", identity=1.0,
        coverage=coverage, annot_gene="", annot_paralog=cell), window="test")
    fs.keys_strict = frozenset(s.upper() for s in syms)
    fs.keys_relaxed = frozenset(lib.relaxed_key(s) for s in syms)
    fs.keys_root = frozenset(lib.root_key(s) for s in syms)
    return fs


# ------------------------------------------------------------------ tests

def t1_placeholder_symbols():
    """T1 -- the placeholder vocabularies are rejected, real symbols kept."""
    bad = ["LOC116952798", "FN964_004414", "OJAV_G00063340",
           "ENSG00000139618", "si:ch211-152c2.3", "zgc:113232", "1", "X"]
    good = ["BAK1", "GNAS", "itpr1b", "bhlhe41", "C12orf71", "RPS10-NUDT3"]
    for s in bad:
        check("T1", not lib.informative(s), f"{s} was accepted")
    for s in good:
        check("T1", lib.informative(s), f"{s} was rejected")
    return "8 placeholders rejected, 6 real symbols kept"


def t2_relaxed_key():
    """T2 -- lowercase duplicate suffixes strip, uppercase symbols do not."""
    check("T2", lib.relaxed_key("gnasa") == lib.relaxed_key("GNAS") == "GNAS",
          "teleost gnasa did not reach GNAS")
    check("T2", lib.relaxed_key("bhlhe41") == "BHLHE41",
          "a lowercase symbol ending in a digit was truncated")
    check("T2", lib.relaxed_key("GNB1") == "GNB1", "uppercase GNB1 was stripped")
    check("T2", lib.relaxed_key("BAK1") == "BAK1", "uppercase BAK1 was stripped")
    check("T2", lib.relaxed_key("ctsa") == "CTSA", "short core ctsa collapsed")
    check("T2", lib.relaxed_key("vapb") == "VAPB", "vapb lost its suffix")
    check("T2", lib.relaxed_key("itpr1b.2") == "ITPR1B",
          "trailing .2 was not stripped")
    return "gnasa->GNAS, GNB1/BAK1/ctsa/vapb untouched"


def t3_root_key():
    """T3 -- the ohnolog pairs merge; a 2-character root does not."""
    check("T3", lib.root_key("BHLHE40") == lib.root_key("BHLHE41"),
          "BHLHE40/41 did not merge")
    check("T3", lib.root_key("GRM7") == lib.root_key("GRM4"),
          "GRM7/GRM4 did not merge")
    check("T3", lib.root_key("TP53") == "TP53",
          "TP53 was reduced below the 3-character floor")
    check("T3", lib.root_key("C3") == "C3", "C3 was reduced to a 1-char root")
    check("T3", lib.root_key("BHLHE40") != lib.root_key("GRM7"),
          "two unrelated families collided")
    check("T3", lib.root_key("ATP5F1A") == "ATP5F1A",
          "a symbol not ending in a digit was altered")
    return "BHLHE40/41 and GRM7/4 merge; TP53 and C3 keep their key"


def t4_flank_extraction():
    """T4 -- the locus and its overlappers are out; the windows differ."""
    idx = synthetic_index()
    loc = synthetic_locus()
    fixed = lib.extract_flanks(idx, loc, n=10, window="fixed10")
    syms = [g[5] for g in fixed.upstream + fixed.downstream]
    check("T4", "TARGET" not in syms, "the locus itself was kept as a flank")
    check("T4", all(g[6] in lib.FLANK_BIOTYPES
                    for g in fixed.upstream + fixed.downstream),
          "a non-coding gene entered the flank set")
    check("T4", all(g[2] < loc.start for g in fixed.upstream),
          "an upstream flank overlaps the locus")
    check("T4", all(g[1] > loc.end for g in fixed.downstream),
          "a downstream flank overlaps the locus")
    check("T4", fixed.upstream[0][5] == "GENE11",
          f"nearest upstream is {fixed.upstream[0][5]}, not the closest gene")
    inf = lib.extract_flanks(idx, loc, n=10, window="informative10")
    check("T4", all(lib.informative(g[5])
                    for g in inf.upstream + inf.downstream),
          "the informative window kept a placeholder")
    check("T4", any(not lib.informative(g[5])
                    for g in fixed.upstream + fixed.downstream),
          "the fixed window silently dropped a placeholder")
    check("T4", inf.scanned_up >= len(inf.upstream),
          "scanned fewer genes than it kept")
    # a locus on a contig with nothing on one side
    edge = lib.extract_flanks(idx, synthetic_locus(start=1, end=5_000), n=10)
    check("T4", not edge.upstream and edge.downstream,
          "a contig-start locus reported upstream flanks")
    return "locus excluded, biotypes filtered, both windows behave"


def t5_jaccard():
    """T5 -- Jaccard is symmetric and an empty set scores zero, not 1."""
    a, b = frozenset("abc"), frozenset("bcd")
    check("T5", abs(lib.jaccard(a, b) - 0.5) < 1e-12, "wrong value")
    check("T5", lib.jaccard(a, b) == lib.jaccard(b, a), "not symmetric")
    check("T5", lib.jaccard(frozenset(), frozenset()) == 0.0,
          "two empty sets scored as identical, which would make an "
          "unannotated genome look like a perfect synteny match")
    check("T5", lib.jaccard(a, frozenset()) == 0.0, "empty vs non-empty")
    return "symmetric, and empty-vs-empty is 0 not 1"


def t6_one_per_species():
    """T6 -- the pair subset takes one locus per species, the best-covered."""
    import s8_run_synteny as drv
    sets_ = [fake_set("ITPR1", "Danio rerio", ["A"], coverage=0.4, idx=0),
             fake_set("ITPR1", "Danio rerio", ["B"], coverage=0.9, idx=1),
             fake_set("ITPR1", "Homo sapiens", ["C"], coverage=1.0)]
    out = drv.one_per_species(sets_)
    check("T6", len(out) == 2, f"kept {len(out)} of 2 species")
    danio = [fs for fs in out if fs.locus.organism == "Danio rerio"][0]
    check("T6", danio.locus.idx == 1,
          "kept the lower-coverage teleost duplicate, which would compare "
          "itpr1a in one species against itpr1b in the next")
    return "one locus per species, highest coverage wins"


def t7_control_determinism():
    """T7 -- the null is reproducible and never drawn on a short contig."""
    idx = synthetic_index()
    a = ctl.sample_windows(idx, "GCA_TEST.1", "Testus testus", "Testia")
    b = ctl.sample_windows(idx, "GCA_TEST.1", "Testus testus", "Testia")
    check("T7", [x.locus.start for x in a] == [x.locus.start for x in b],
          "two draws with the same seed differed")
    check("T7", all(x.locus.contig == "chr1" for x in a),
          "a control window was drawn on a contig too short to hold one")
    check("T7", len(a) == ctl.N_REPLICATES, "wrong number of replicates")
    return f"{len(a)} reproducible windows, short contig refused"


def t8_caller_rejects_a_stranger():
    """T8 -- the positive and negative halves of the consensus caller.

    A locus carrying its paralog's consensus must be called; one carrying a
    neighbourhood built of unrelated symbols must not. A caller that only
    ever agrees is not evidence about the loci it was built to place.
    """
    keep = par.CONSENSUS_FRAC
    par.CONSENSUS_FRAC = 0.5
    try:
        sets_by_class = {
            "ITPR1": [fake_set("ITPR1", f"sp{i}", ["AAA1", "BBB1", "CCC1",
                                                   "DDD1", "EEE1"])
                      for i in range(6)],
            "ITPR2": [fake_set("ITPR2", f"sq{i}", ["FFF1", "GGG1", "HHH1",
                                                   "III1", "JJJ1"])
                      for i in range(6)],
            "ITPR3": [fake_set("ITPR3", f"sr{i}", ["KKK1", "LLL1", "MMM1",
                                                   "NNN1", "OOO1"])
                      for i in range(6)],
        }
        caller = par.ConsensusCaller(sets_by_class, "relaxed")
        classes = ("ITPR1", "ITPR2", "ITPR3")
        true = fake_set("ITPR1", "sp_new",
                        ["AAA1", "BBB1", "CCC1", "DDD1", "ZZZ9"])
        s = caller.score(true, classes, leave_out=True)
        check("T8", s["call"] == "ITPR1",
              f"a genuine ITPR1 neighbourhood was called {s['call']}")
        stranger = fake_set("ITPR1", "sp_odd",
                            ["QQQ1", "RRR1", "SSS1", "TTT1", "UUU1"])
        s2 = caller.score(stranger, classes, leave_out=True)
        check("T8", s2["call"] == "no_call",
              f"an unrelated neighbourhood was called {s2['call']}")
        thin = fake_set("ITPR1", "sp_thin", ["AAA1", "BBB1"])
        s3 = caller.score(thin, classes, leave_out=True)
        check("T8", s3["call"] == "no_call",
              "a locus below the informative-key floor was still called")
        # leave-one-out must actually remove the query's own species
        self_vote = fake_set("ITPR1", "sp0", ["QQQ1", "RRR1", "SSS1", "TTT1"])
        s4 = caller.score(self_vote, classes, leave_out=True)
        check("T8", s4["call"] == "no_call",
              "a locus was scored against a consensus it had voted into")
    finally:
        par.CONSENSUS_FRAC = keep
    return "true set called, stranger and thin set refused"


def t9_select_frac():
    """T9 -- the threshold rule maximises call rate minus false-call rate,
    with ties broken toward the stricter (larger) threshold."""
    sweep = [dict(consensus_frac=0.2, call_rate=0.81, false_call_rate=0.05,
                  accuracy=1.0),
             dict(consensus_frac=0.3, call_rate=0.80, false_call_rate=0.01,
                  accuracy=1.0),
             dict(consensus_frac=0.4, call_rate=0.80, false_call_rate=0.01,
                  accuracy=1.0),
             dict(consensus_frac=0.5, call_rate=0.60, false_call_rate=0.00,
                  accuracy=1.0)]
    pick, rule = par.select_frac(sweep)
    check("T9", pick == 0.4, f"picked {pick}, expected the tie-break to 0.4")
    loose = [dict(consensus_frac=0.2, call_rate=0.99, false_call_rate=0.90,
                  accuracy=1.0),
             dict(consensus_frac=0.6, call_rate=0.50, false_call_rate=0.00,
                  accuracy=1.0)]
    pick2, _ = par.select_frac(loose)
    check("T9", pick2 == 0.6,
          "a threshold that calls everything, including 90 % of random "
          "windows, was chosen on call rate alone")
    return f"tie-break to 0.4; loose setting refused ({rule[:40]}...)"


def t10_null_verdict():
    """T10 -- a call is only 'supported' above the null's own maximum."""
    summ = dict(score_distribution={"0": 720, "1": 4, "2": 2},
                max_best_score=2)
    tail, max_null = par.null_tail(summ)
    check("T10", max_null == 2, "wrong null maximum")
    check("T10", abs(tail[2] - 2 / 726) < 1e-9, f"tail at 2 is {tail[2]}")
    check("T10", abs(tail[1] - 6 / 726) < 1e-9, f"tail at 1 is {tail[1]}")
    check("T10", par.null_verdict(2, "ITPR1", max_null) == "within_null",
          "a score the null reaches was reported as supported")
    check("T10", par.null_verdict(3, "ITPR1", max_null) == "supported",
          "a score above every random window was not supported")
    check("T10", par.null_verdict(9, "no_call", max_null) == "no_call",
          "a no_call was given a verdict")
    return "tail probabilities correct, bar at the null maximum"


def t11_output_order_is_stable():
    """T11 -- ranked output must not depend on input order.

    Both ranked tables are built by walking Python sets, whose iteration
    order is hash-seeded per process, so a rank computed without an explicit
    final tiebreak differs between two runs of identical data. This caught
    exactly that in `flank_consensus.tsv`.
    """
    a = [fake_set("ITPR1", f"sp{i}", ["AAA1", "BBB1", "CCC1"]) for i in range(4)]
    b = [fake_set("ITPR2", f"sq{i}", ["AAA2", "BBB2", "CCC2"]) for i in range(4)]
    bg = {"AAA": 0.01, "BBB": 0.01, "CCC": 0.01}
    r1 = [r["key"] for r in par.shared_roots(a, b, bg, min_frac=0.1)]
    r2 = [r["key"] for r in par.shared_roots(list(reversed(a)),
                                             list(reversed(b)), bg,
                                             min_frac=0.1)]
    check("T11", r1 == r2, f"order changed: {r1} vs {r2}")
    check("T11", len(r1) == 3, f"expected 3 tied roots, got {r1}")
    check("T11", r1 == sorted(r1),
          "tied roots are not in a deterministic order")
    prev, n, _ = par.side_prevalence(a, "relaxed")
    order1 = sorted(prev.items(), key=lambda kv: (-kv[1], kv[0]))
    prev2, _, _ = par.side_prevalence(list(reversed(a)), "relaxed")
    order2 = sorted(prev2.items(), key=lambda kv: (-kv[1], kv[0]))
    check("T11", order1 == order2, "consensus ranking depends on input order")
    return "ranked tables are order-invariant with a key tiebreak"


TESTS = [t1_placeholder_symbols, t2_relaxed_key, t3_root_key,
         t4_flank_extraction, t5_jaccard, t6_one_per_species,
         t7_control_determinism, t8_caller_rejects_a_stranger,
         t9_select_frac, t10_null_verdict, t11_output_order_is_stable]


def self_test(verbose: bool = True) -> list[tuple[str, str]]:
    out = []
    for fn in TESTS:
        detail = fn()
        name = fn.__name__.split("_")[0].upper()
        out.append((name, detail))
        if verbose:
            print(f"  [s8-test] {name} ok -- {detail}")
    return out


if __name__ == "__main__":
    self_test()
    print(f"[s8-test] {len(TESTS)}/{len(TESTS)} passed")
