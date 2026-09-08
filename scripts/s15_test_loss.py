"""s15_test_loss.py — constructed negative controls for the S15 rules, run
before anything is written (`s5_bait_screen.self_test()`'s pattern).

Every rule in this task returns a plausible number when it is wrong.  A
reconstruction built without the known-locus exclusion reassembles a
"missing" gene out of the genome's other paralogs and reports 0.95.  A
state machine whose order is wrong calls a cyclostome cell `absent`
instead of `paralog_unassignable` and manufactures four losses.  A bar
read off a decoy that contains real genes lands in the middle of the
positives.  So these are checks on *refusal*: each case is constructed to
break one rule, and the rule responsible has to reject it.

Mutation-tested: three deliberate rule breakages, all three caught (see
`SELF_MUTATIONS` for what was broken and which test fired).
"""

from __future__ import annotations

import math

import s15_lib as lib
import s15_calibrate_recon as cal
import s15_matrix as mx
import s15_reconstruct as recon
import s15_species_tree as tree
import s15_states as st

SELF_MUTATIONS = (
    "dropped the known-locus filter from reconstruct()  -> T2 fired",
    "moved R5 (paralog_unassignable) after R7 (absent)   -> T6 fired",
    "used every elsewhere-attributed region as the decoy -> T9 fired",
)


class Fail(Exception):
    pass


def _check(cond, msg):
    if not cond:
        raise Fail(msg)


# ------------------------------------------------------- interval algebra

def t1_union():
    _check(lib.union_length([]) == 0, "empty union is not 0")
    _check(lib.union_length([(1, 10)]) == 10, "closed interval length wrong")
    _check(lib.union_length([(1, 10), (5, 20)]) == 20, "overlap not merged")
    _check(lib.union_length([(1, 10), (11, 20)]) == 20,
           "abutting intervals not merged")
    _check(lib.union_length([(1, 10), (12, 20)]) == 19,
           "a 1-base gap was closed")
    _check(lib.union_length([(20, 5)]) == 16,
           "a reversed interval was not normalised")


# --------------------------------------------------- the known-locus gate

def _hsp(contig, q_lo, q_hi, s_lo, s_hi, bait="B|X|Sp|ITPR1", bits=100.0):
    return dict(bait=bait, contig=contig, q_lo=q_lo, q_hi=q_hi,
                s_lo=s_lo, s_hi=s_hi, bits=bits, evalue=1e-40)


def t2_known_locus_excluded():
    """An HSP inside a locus the aligner already placed must not be used.

    This is the cross-paralog control: without it a shattered ITPR2 cell
    reassembles itself out of the genome's intact ITPR1 gene, at 65 %
    identity, and reports a presence that is another gene.
    """
    bl = {"B|X|Sp|ITPR1": 1000}
    hsps = [_hsp("c1", 1, 1000, 5000, 8000)]
    foot = [("c1", 4000, 9000)]
    got = recon.reconstruct("A", "Sp", "V", "ITPR1", hsps, [], foot, bl, "own")
    _check(got is not None and got.covered_aa == 1000,
           "the unfiltered reconstruction did not run")
    blocked = recon.reconstruct("A", "Sp", "V", "ITPR1", hsps,
                                [("c1", 5500, 6000)], foot, bl, "own")
    _check(blocked is None,
           "an HSP inside a known locus was used in the reconstruction")
    # and the pad has to be applied, not just the interval
    near = recon.reconstruct("A", "Sp", "V", "ITPR1", hsps,
                             [("c1", 8000 + lib.KNOWN_PAD - 10,
                               8000 + lib.KNOWN_PAD)], foot, bl, "own")
    _check(near is None, "the known-locus pad was not applied")


def t3_one_reference_at_a_time():
    """Coverage must never be unioned across baits.

    Three orthologous baits each covering a different third of their own
    reference would otherwise report 100 % coverage of a protein no bait
    achieved.
    """
    bl = {"b1|X|Sp|ITPR1": 300, "b2|X|Sp|ITPR1": 300, "b3|X|Sp|ITPR1": 300}
    hsps = [_hsp("c1", 1, 100, 10, 400, "b1|X|Sp|ITPR1"),
            _hsp("c1", 101, 200, 500, 900, "b2|X|Sp|ITPR1"),
            _hsp("c1", 201, 300, 1000, 1400, "b3|X|Sp|ITPR1")]
    got = recon.reconstruct("A", "Sp", "V", "ITPR1", hsps, [],
                            [("c1", 1, 9999)], bl, "own")
    _check(abs(got.coverage - 100 / 300) < 1e-9,
           f"coverage was unioned across baits: {got.coverage:.3f}")


def t4_exonic_is_subject_side():
    """gene_equiv counts genomic sequence, so two baits hitting the same
    place must not count it twice."""
    bl = {"b1|X|Sp|ITPR1": 1000, "b2|X|Sp|ITPR1": 1000}
    hsps = [_hsp("c1", 1, 500, 1, 1500, "b1|X|Sp|ITPR1"),
            _hsp("c1", 1, 500, 1, 1500, "b2|X|Sp|ITPR1")]
    got = recon.reconstruct("A", "Sp", "V", "ITPR1", hsps, [],
                            [("c1", 1, 9999)], bl, "own")
    _check(got.exonic_nt == 1500,
           f"the same genomic interval was counted twice: {got.exonic_nt}")


def t5_tiling_vs_pileup():
    """unique_frac must tell a tiling from a pile-up."""
    bl = {"B|X|Sp|ITPR1": 300}
    tiled = [_hsp("c1", 1, 100, 1, 300), _hsp("c2", 101, 200, 1, 300),
             _hsp("c3", 201, 300, 1, 300)]
    piled = [_hsp("c1", 1, 300, 1, 900), _hsp("c2", 1, 300, 1, 900),
             _hsp("c3", 1, 300, 1, 900)]
    a = recon.reconstruct("A", "S", "V", "ITPR1", tiled, [],
                          [("c1", 1, 999), ("c2", 1, 999), ("c3", 1, 999)],
                          bl, "own")
    b = recon.reconstruct("A", "S", "V", "ITPR1", piled, [],
                          [("c1", 1, 999), ("c2", 1, 999), ("c3", 1, 999)],
                          bl, "own")
    _check(a.unique_frac == 1.0, f"a clean tiling scored {a.unique_frac}")
    _check(b.unique_frac == 0.0, f"a three-deep pile-up scored {b.unique_frac}")


# ------------------------------------------------------------ state rules

def _cell(status="absent", cov=0.0, n50=50_000_000, cls="ITPR2"):
    return dict(accession="GCA_TEST", organism="Testus testus",
                vclass="Aves", **{"class": cls}, status=status,
                best_coverage=cov, contig_n50=n50)


def _summary(control=True, cells=None, other=None):
    return dict(organism="Testus testus", vclass="Aves",
                control_ok=control, cells=cells or {},
                other_loci=other or [])


def t6_cyclostome_is_not_a_loss():
    """A cell with no locus beside spare family loci is D45, not a loss."""
    cells = {"ITPR1": {"loci": [{}, {}, {}]}, "ITPR2": {"loci": []},
             "ITPR3": {"loci": []}}
    got = st.classify(_cell(), _summary(cells=cells), None, 0.14)
    _check(got["state"] == "paralog_unassignable",
           f"a cell beside 2 spare ITPR loci was called {got['state']}")
    _check(got["n_spare_itpr_loci"] == 2,
           f"spare loci miscounted: {got['n_spare_itpr_loci']}")


def t7_uncontrolled_genome_supports_nothing():
    cells = {"ITPR1": {"loci": []}, "ITPR2": {"loci": []},
             "ITPR3": {"loci": []}}
    got = st.classify(_cell(), _summary(control=False, cells=cells), None, 0.14)
    _check(got["state"] == "no_control",
           f"a genome whose RyR control failed yielded {got['state']}")


def t8_fragmented_assembly_cannot_prove_absence():
    """Below D4's bar the answer is undecidable, never absent."""
    cells = {c: {"loci": []} for c in lib.ITPR_CELLS}
    low = st.classify(_cell(n50=20_000), _summary(cells=cells), None, 0.14)
    _check(low["state"] == "undecidable_contiguity",
           f"a 20 kb-N50 assembly yielded {low['state']}")
    high = st.classify(_cell(n50=50_000_000), _summary(cells=cells), None, 0.14)
    _check(high["state"] == "absent",
           f"a contiguous, controlled, empty, spare-free cell yielded "
           f"{high['state']} — the loss state is unreachable")
    # and a reconstruction over the bar outranks the absence
    frag = st.classify(_cell(n50=50_000_000), _summary(cells=cells),
                       dict(coverage=0.9, gene_equiv=0.9, n_contigs=7), 0.14)
    _check(frag["state"] == "present_fragmented",
           f"a 0.90 reconstruction yielded {frag['state']}")


# --------------------------------------------------------- the decoy rule

def t9_decoy_excludes_co_shattered_paralogs():
    """The negative must not contain genes that are themselves shattered."""
    regions = [
        dict(contig="c1", start=1, end=100, assigned_clade="ITPR2"),
        dict(contig="c2", start=1, end=100, assigned_clade="ITPR1"),
        dict(contig="c3", start=1, end=100, assigned_clade="ITPR3"),
    ]
    status = {"ITPR1": "found_annotated", "ITPR2": "tblastn_trace",
              "ITPR3": "tblastn_trace"}
    got = recon.region_index(regions, "ITPR2", status)
    _check([r[0] for r in got["own_clade"]] == ["c1"], "own clade wrong")
    _check([r[0] for r in got["decoy_accounted"]] == ["c2"],
           f"the decoy is not the accounted-for paralog: {got['decoy_accounted']}")
    _check([r[0] for r in got["co_trace"]] == ["c3"],
           f"a co-shattered paralog leaked out of co_trace: {got['co_trace']}")


def t10_calibration_refuses_to_pass_vacuously():
    rows = [dict(scope="own_clade", bait="b", coverage=0.9, gene_equiv=0.9)
            for _ in range(3)]
    rows += [dict(scope="decoy_accounted", bait="b", coverage=0.01,
                  gene_equiv=0.01)]
    rows += [dict(scope="co_trace", bait="b", coverage=0.5, gene_equiv=0.5)]
    _, stats = cal.calibrate(rows)
    _check(not stats["usable"],
           "a calibration read off 3 positives and 1 negative was accepted")
    try:
        cal.bar(stats)
    except RuntimeError:
        pass
    else:
        raise Fail("bar() handed out a threshold the calibration rejected")


def t11_youden_reports_non_separation():
    y = lib.youden([0.4, 0.5, 0.6], [0.35, 0.55, 0.65])
    _check(y["j"] < 0.7, f"overlapping populations scored J = {y['j']:.2f}")


# ------------------------------------------------------------- statistics

def t12_sign_test_and_bh():
    r = lib.sign_test([1, 1, 1, 1, 1, 1])
    _check(abs(r["p"] - 2 * (1 / 64)) < 1e-12, f"sign-test p wrong: {r['p']}")
    _check(lib.sign_test([0, 0, 0])["n"] == 0, "ties were not dropped")
    _check(lib.sign_test([0, 0, 0])["n_ties"] == 3, "ties were not counted")
    import s15_integrity as integ
    q = integ._bh([0.01, 0.02, 0.03, 0.04])
    _check(all(a <= b + 1e-12 for a, b in zip(q, q[1:])),
           f"BH output is not monotone: {q}")
    _check(abs(q[0] - 0.04) < 1e-9, f"BH first value wrong: {q[0]}")
    _check(all(x >= p - 1e-12 for x, p in zip(q, [.01, .02, .03, .04])),
           "a BH q-value fell below its own p-value")


def t13_spearman_direction():
    rho, p, n = lib.spearman([1, 2, 3, 4, 5, 6, 7, 8],
                             [8, 7, 6, 5, 4, 3, 2, 1])
    _check(abs(rho + 1.0) < 1e-9, f"a perfect inversion scored rho = {rho}")
    _check(p < 0.01, f"a perfect inversion scored p = {p}")
    rho2, _, _ = lib.spearman([1, 2, 3], [1, 1, 1])
    _check(math.isnan(rho2), "a constant series returned a correlation")


# ------------------------------------------------------------- the tree

def t14_tree_is_a_tree():
    """Membership, tip identity, and the naming rule.

    The label of an internal node is the *most exclusive* taxon containing
    exactly its tips, because `_suppress_unary` keeps the deepest named
    node of a unary chain.  With only two mammals sampled, the clade
    separating them from a bird is therefore `Euarchontoglires` and not
    `Mammalia` — both are true of that node and the deeper one carries
    more information.  So the test asserts membership, and that the label
    is *some* taxon on the shared path, never a particular rank.
    """
    man = [dict(accession="A1", organism="Homo sapiens", taxid="9606"),
           dict(accession="A2", organism="Mus musculus", taxid="10090"),
           dict(accession="A3", organism="Gallus gallus", taxid="9031")]
    root, tips, audit = tree.build(man)
    _check(len(tips) == 3, f"{len(tips)} tips from 3 genomes")
    leaves = [n.name for n in root.leaves()]
    _check(sorted(leaves) == ["A1", "A2", "A3"],
           f"tips are not the accessions: {leaves}")
    _check(all(len(n.children) != 1 for n in root.walk()
               if not n.is_tip and not n.children[0].is_tip),
           "an internal unary node survived")
    clades = {frozenset(n.tips): n.label for n in root.walk()
              if not n.is_tip}
    _check(frozenset({"A1", "A2"}) in clades,
           f"the two mammals are not a clade: {list(clades)}")
    _check(clades[frozenset({"A1", "A2"})] in
           ("Mammalia", "Theria", "Eutheria", "Boreoeutheria",
            "Euarchontoglires"),
           f"the mammal clade is labelled "
           f"{clades[frozenset({'A1', 'A2'})]}, which is not on the "
           f"human-mouse path")
    _check(clades[frozenset({"A1", "A2", "A3"})] in ("Amniota", "Tetrapoda"),
           "the root of a human/mouse/chicken sample is not an amniote node")


def t15_implied_copies_counts_spares():
    """A genome whose three loci are all filed in one cell has three."""
    matrix = [dict(accession="A", organism="X", vclass="V", cell=c,
                   state=("present_single_locus" if c == "ITPR1"
                          else "paralog_unassignable"),
                   best_coverage=1.0, recon_gene_equiv=0.0, contig_n50=10 ** 7,
                   contig_spans_gene=1, control_ok=1, n_spare_itpr_loci=2)
              for c in lib.ITPR_CELLS]
    got = mx.implied_copies(matrix)[0]
    _check(abs(got["implied_copies"] - 3.0) < 1e-9,
           f"a 3-locus cyclostome genome implied "
           f"{got['implied_copies']} copies")


TESTS = [t1_union, t2_known_locus_excluded, t3_one_reference_at_a_time,
         t4_exonic_is_subject_side, t5_tiling_vs_pileup,
         t6_cyclostome_is_not_a_loss, t7_uncontrolled_genome_supports_nothing,
         t8_fragmented_assembly_cannot_prove_absence,
         t9_decoy_excludes_co_shattered_paralogs,
         t10_calibration_refuses_to_pass_vacuously,
         t11_youden_reports_non_separation, t12_sign_test_and_bh,
         t13_spearman_direction, t14_tree_is_a_tree,
         t15_implied_copies_counts_spares]


def self_test(verbose: bool = True) -> dict:
    passed, failed = [], []
    for fn in TESTS:
        try:
            fn()
            passed.append(fn.__name__)
            if verbose:
                print(f"  ok    {fn.__name__}")
        except Exception as exc:                      # noqa: BLE001
            failed.append(f"{fn.__name__}: {exc}")
            if verbose:
                print(f"  FAIL  {fn.__name__}: {exc}")
    return dict(n_tests=len(TESTS), n_passed=len(passed),
                n_failed=len(failed), failures=failed,
                mutations_caught=list(SELF_MUTATIONS),
                status="pass" if not failed else "fail")


if __name__ == "__main__":
    import sys
    r = self_test()
    print(f"\n{r['n_passed']}/{r['n_tests']} passed")
    sys.exit(0 if r["status"] == "pass" else 1)
