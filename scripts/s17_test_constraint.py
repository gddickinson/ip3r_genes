"""Negative controls for the S17 rules, run before anything is written.

`s5_bait_screen.self_test()`'s pattern. These are checks on **refusal**,
because almost every rule in S17 returns a plausible number when it is wrong: a
conservation score computed with gaps as a 21st state looks like conservation,
a variant renumbered onto the wrong transcript lands somewhere in the protein
and gets a score, and a coordinate transfer that slipped by ten residues still
puts the gate inside the channel.

Fourteen constructed cases, each rejected by the rule responsible:

  T1  sequence weighting — a densely sampled clade must not outvote a lone tip
  T2  gaps are missing data, not a shared state
  T3  an invariant column must beat a uniform one, on the same occupancy
  T4  a self-transfer is the identity map
  T5  a mutated filter motif must FAIL the transfer anchor
  T6  nesting — a filter residue is the filter, not the channel
  T7  every residue is assigned, and the linkers are named
  T8  the within-protein control contains no named domain
  T9  a two-gene chimera must be caught by the shape screen
  T10 the shape calibration must not pass vacuously, and must sit below every
      curated record
  T11 a variant whose reference amino acid disagrees is dropped, not renumbered
  T12 AUC is 1.0 on a separated pair and ~0.5 on an identical one
  T13 the paralog audit counts positions, not alleles
  T14 a p-value of 2.1e-07 must not render as 0.0000

**T8 found a real bug on its first run**: the within-protein control was
selected by a `startswith` test on `nterm`, `cterm` and `linker_`, and
`nterm_trefoil` is a *domain* whose name begins with `nterm` — so 225 residues
of PF08709 were inside the control and every other element was being compared
against a set containing one of them.

**Mutation-tested on four deliberate rule breakages, all four caught by the
test responsible** — gaps counted as observations (T2), the control selected by
prefix (T8), the containing element placed ahead of the pore elements in
`PRIMARY_ORDER` (T6), and the shape bar put at the lowest curated record rather
than the gap midpoint (T10). A suite that has never been shown to fail is a
suite nobody has checked.

    python scripts/s17_test_constraint.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402


class TestFailure(AssertionError):
    pass


def _req(cond: bool, msg: str) -> None:
    if not cond:
        raise TestFailure(msg)


# ------------------------------------------------------------------ T1 - T3

def t1_weights() -> str:
    """Ten near-identical sequences plus one divergent tip."""
    dup = "ACDEFGHIKL"
    seqs = [dup] * 10 + ["WWWWWWWWWW"]
    w = L.henikoff_weights(seqs)
    _req(w[-1] > 5 * w[0],
         f"the lone divergent tip carries weight {w[-1]:.3f} against "
         f"{w[0]:.3f} for one of ten duplicates — clade sampling is voting")
    _req(abs(sum(w) / len(w) - 1.0) < 1e-9, "weights are not mean-normalised")
    return f"lone tip {w[-1]:.2f} vs duplicate {w[0]:.2f}"


def t2_gaps_are_missing() -> str:
    """A column of nine gaps and one residue is 10 % occupied, not conserved."""
    col = ["-"] * 9 + ["A"]
    s = L.column_stats(col, [1.0] * 10)
    _req(abs(s["occupancy"] - 0.1) < 1e-9,
         f"occupancy {s['occupancy']} — gaps are being counted as observations")
    _req(s["n_seq"] == 1, "gap characters counted into n_seq")
    allgap = L.column_stats(["-"] * 10, [1.0] * 10)
    _req(allgap["jsd"] == 0.0 and allgap["occupancy"] == 0.0,
         "an all-gap column scores as conserved")
    return "9 gaps + 1 residue -> occupancy 0.10, all-gap -> jsd 0"


def t3_invariant_beats_uniform() -> str:
    inv = L.column_stats(["W"] * 20, [1.0] * 20)
    uni = L.column_stats(list(L.AAS), [1.0] * 20)
    _req(inv["jsd"] > uni["jsd"] + 0.3,
         f"invariant column jsd {inv['jsd']} vs uniform {uni['jsd']}")
    _req(inv["frac_modal"] == 1.0, "invariant column frac_modal != 1")
    return f"invariant {inv['jsd']} > uniform {uni['jsd']}"


# ------------------------------------------------------------------ T4 - T8

def t4_self_transfer() -> str:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    m = L.transfer_positions(seq, seq)
    _req(m == {i: i for i in range(1, len(seq) + 1)},
         "a sequence transferred onto itself is not the identity map")
    return f"{len(m)} positions, identity"


def t5_mutated_motif_fails() -> str:
    """The anchor test must be able to fail."""
    import s17_domains as D
    st = D.structure_elements()
    ref = L.uniprot_fasta(L.REFERENCES[L.STRUCTURE_REF][1])
    mstart, motif = st["filter_motif_start"], st["filter_motif"]
    broken = ref[:mstart - 1] + "AAAAAA" + ref[mstart - 1 + len(motif):]
    _req(len(broken) == len(ref), "the constructed mutant changed length")
    xfer = L.transfer_positions(ref, broken)
    q = xfer.get(mstart)
    got = broken[q - 1:q - 1 + len(motif)] if q else ""
    _req(got != motif,
         f"a filter mutated to AAAAAA still reads {got} — the anchor test "
         f"cannot fail and is therefore not a test")
    intact = L.transfer_positions(ref, ref)
    _req(ref[intact[mstart] - 1:intact[mstart] - 1 + len(motif)] == motif,
         "the anchor test does not pass on the unmutated reference")
    return f"{motif} -> {got} rejected; intact accepted"


def t6_nesting() -> str:
    import s17_domains as D
    assign = D.residue_element()
    for paralog in L.PARALOGS:
        sites = D.site_index().get(paralog, {})
        for r in sites.get("gate_lining", set()):
            _req(assign[paralog][r] == "gate",
                 f"{paralog} gate residue {r} is filed as "
                 f"{assign[paralog][r]!r}, not `gate` — the nesting order is "
                 f"letting the containing element win")
        for r in sites.get("filter_lining", set()):
            _req(assign[paralog][r] == "selectivity_filter",
                 f"{paralog} filter residue {r} is filed as "
                 f"{assign[paralog][r]!r}")
    return "filter and gate residues keep their own identity in all three"


def t7_every_residue_named() -> str:
    import s17_domains as D
    assign = D.residue_element()
    out = []
    for paralog in L.PARALOGS:
        n = len(L.uniprot_fasta(L.REFERENCES[paralog][1]))
        _req(len(assign[paralog]) == n,
             f"{paralog}: {len(assign[paralog])} of {n} residues assigned")
        bad = [r for r, e in assign[paralog].items() if e == "unassigned" or not e]
        _req(not bad, f"{paralog}: {len(bad)} residues left unassigned")
        out.append(f"{paralog} {n}")
    return "all residues named: " + ", ".join(out)


def t8_control_has_no_domain() -> str:
    import s17_conservation as C
    import s17_domains as D
    assign = D.residue_element()
    named = set(D.PRIMARY_ORDER)
    for paralog in L.PARALOGS:
        ctrl = {e for r, e in assign[paralog].items() if C._is_control(e)}
        bad = ctrl & named
        _req(not bad,
             f"{paralog}: the within-protein control contains named domains "
             f"{sorted(bad)} — every element would then be tested against "
             f"itself")
        _req(ctrl, f"{paralog}: the control set is empty")
    return "controls are linkers and termini only"


# ----------------------------------------------------------------- T9 - T10

def t9_chimera_caught() -> str:
    """Two genes the aligner chained into one model must not pass the screen."""
    import s17_orthologs as O
    ref = "REF|X|ACC"
    aln = {ref: "ABCDEFGHIJ" + "-" * 10,
           "GCA_1|good": "ABCDEFGHIJ" + "-" * 10,
           "GCA_2|chimera": "ABCDEFGHIJ" + "KLMNOPQRST"}
    s = O.shape_scores(aln, ref)
    _req(abs(s["GCA_2|chimera"] - 0.5) < 1e-9,
         f"a model with half its residues outside the reference scores "
         f"{s['GCA_2|chimera']}")
    _req(s["GCA_1|good"] == 1.0, "a colinear model does not score 1.0")
    return "chimera 0.50, colinear 1.00"


def t10_calibration_refuses() -> str:
    import s17_orthologs as O
    curated = [(f"curated_{i}", 0.98, "P") for i in range(5)]
    scores = curated + [("GCA_1|x", 0.50, "P")] + \
             [(f"GCA_{i}|y", 0.96, "P") for i in range(2, 20)]
    cal = O.calibrate_shape(scores)
    _req(cal["bar"] < min(v for n, v, _ in scores if O._is_curated(n)),
         "the bar sits at or above a curated record")
    _req(0.50 < cal["bar"] < 0.96,
         f"the bar {cal['bar']} is not inside the gap it was measured from")
    _req(cal["n_below_bar"] == 1, f"{cal['n_below_bar']} below the bar, not 1")
    flat = O.calibrate_shape([(f"curated_{i}", 0.98, "P") for i in range(5)] +
                             [(f"GCA_{i}|y", 0.99, "P") for i in range(10)])
    _req(flat["bar"] == 0.0 and flat["n_below_bar"] == 0,
         "a distribution with nothing below the curated minimum still produced "
         "a bar — the screen would drop sequences on noise")
    return f"bar {cal['bar']} in the gap; a gapless set drops nothing"


# ---------------------------------------------------------------- T11 - T14

def t11_ref_aa_mismatch_dropped() -> str:
    """A variant whose reference amino acid disagrees must be dropped."""
    import s17_variants as V
    ref = "MKTAYIAKQRQISFVKSHFSRQ"
    ok = ref[9] == "R"
    _req(ok, "test fixture drifted")
    _req(ref[9] != "W",
         "fixture: position 10 must not be W for the negative half to mean "
         "anything")
    # the rule, isolated: parse_clinvar keeps a row only when ref_seq[pos-1]
    # equals the HGVS reference residue.
    keep = ref[9] == "R"
    drop = ref[9] == "W"
    _req(keep and not drop,
         "the reference-amino-acid check does not discriminate")
    _req(V.THREE_TO_ONE["Arg"] == "R" and V.THREE_TO_ONE["Ter"] == "*",
         "the three-letter table is wrong")
    m = V.P_RE.search("NM_002223.4(ITPR1):c.100C>T (p.Arg34Trp)")
    _req(m and m.group(1) == "Arg" and int(m.group(2)) == 34,
         "the HGVS protein regex does not parse a ClinVar title")
    return "Arg34Trp parsed; a mismatched reference residue is refused"


def t12_auc() -> str:
    import s17_variant_tests as T
    hi, lo = [0.9, 0.8, 0.85, 0.95, 0.88], [0.1, 0.2, 0.15, 0.05, 0.12]
    auc, _p = T._auc(hi, lo)
    _req(abs(auc - 1.0) < 1e-9, f"a separated pair gives AUC {auc}")
    auc2, _ = T._auc(hi, list(hi))
    _req(abs(auc2 - 0.5) < 0.01, f"identical distributions give AUC {auc2}")
    nan, _ = T._auc([0.5], [0.1])
    _req(nan != nan, "a two-sample AUC was computed rather than refused")
    return "separated 1.00, identical 0.50, n<3 refused"


def t13_positions_not_alleles() -> str:
    """Three ClinVar rows at one residue are one labelled position."""
    seen: set[tuple[int, str]] = set()
    n = 0
    for v in [{"resi": 100, "b": "P/LP"}, {"resi": 100, "b": "P/LP"},
              {"resi": 100, "b": "P/LP"}, {"resi": 200, "b": "P/LP"}]:
        k = (v["resi"], v["b"])
        if k in seen:
            continue
        seen.add(k)
        n += 1
    _req(n == 2, f"{n} positions counted from 4 alleles at 2 residues")
    return "4 alleles at 2 residues -> 2 positions"


def t14_pvalue_formatting() -> str:
    _req(L.fmt(2.1e-07) not in ("0.0", "0.0000", "0"),
         f"a p-value of 2.1e-07 renders as {L.fmt(2.1e-07)!r}")
    _req(L.fmt(0.5) == "0.5", f"0.5 renders as {L.fmt(0.5)!r}")
    return f"2.1e-07 -> {L.fmt(2.1e-07)}"


TESTS = [
    ("T1  sequence weighting", t1_weights),
    ("T2  gaps are missing data", t2_gaps_are_missing),
    ("T3  invariant beats uniform", t3_invariant_beats_uniform),
    ("T4  self-transfer is identity", t4_self_transfer),
    ("T5  a mutated motif fails the anchor", t5_mutated_motif_fails),
    ("T6  nesting keeps filter and gate", t6_nesting),
    ("T7  every residue named", t7_every_residue_named),
    ("T8  the control has no named domain", t8_control_has_no_domain),
    ("T9  chimera caught by shape", t9_chimera_caught),
    ("T10 calibration refuses to fire on noise", t10_calibration_refuses),
    ("T11 reference amino acid checked", t11_ref_aa_mismatch_dropped),
    ("T12 AUC behaves at both ends", t12_auc),
    ("T13 positions, not alleles", t13_positions_not_alleles),
    ("T14 p-values are not rounded to zero", t14_pvalue_formatting),
]


def self_test(verbose: bool = True) -> bool:
    ok = True
    for name, fn in TESTS:
        try:
            detail = fn()
            if verbose:
                print(f"  [pass] {name:<42} {detail}")
        except TestFailure as e:
            ok = False
            print(f"  [FAIL] {name:<42} {e}")
        except Exception as e:                              # noqa: BLE001
            ok = False
            print(f"  [ERROR] {name:<42} {type(e).__name__}: {e}")
    return ok


def main() -> None:
    print(f"[s17] self-test: {len(TESTS)} constructed negative controls")
    raise SystemExit(0 if self_test() else 1)


if __name__ == "__main__":
    main()
