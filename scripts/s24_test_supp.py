"""Negative controls for S24's two guards, run before anything is written.

The pattern is `s5_bait_screen.self_test()`'s: every check here is a check on
*refusal*, because both guards return a plausible answer when they are wrong.
A column map that is one column out still maps every column to a column; a
variant residue read off the wrong table is still an amino acid. Neither
failure has any other symptom, and a supplementary figure whose whole purpose
is to let a reader check the joins would be the last place it showed.

Also here: the two properties the *figures* must have and nothing else tests —
that the group vocabulary they colour by is the one S6 committed, and that the
structure-eligibility rule can actually refuse (it does, twice, on this data,
which is why the rule is worth having).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import figstyle as F                                           # noqa: E402
import s24_lib as L                                            # noqa: E402
import s24_figs_structure as S                                 # noqa: E402
import s24_figs_alignment  # noqa: E402,F401 - import must not break
import s24_figs_inputs     # noqa: E402,F401

FAILURES: list[str] = []
CHECKS: list[str] = []

#: Deliberate rule breakages this suite was mutation-tested against, each
#: recorded with the check that caught it. Written down rather than counted
#: in prose so the report cannot quote a number the suite has outgrown.
MUTATIONS = [
    ("the column walk stops comparing columns", "T2"),
    ("the variant residue check is disabled", "T7"),
    ("the numbering rule admits everything", "T11"),
    ("the one-to-one column guard is disabled", "T4"),
    ("the pdf creation timestamp is restored", "T17"),
]


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append(name)
    if not ok:
        FAILURES.append(f"{name}: {detail}" if detail else name)
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}")


def _raises(fn, *a, **k) -> bool:
    try:
        fn(*a, **k)
    except SystemExit:
        return True
    except Exception:                              # noqa: BLE001
        return False
    return False


# ------------------------------------------------- T1-T5: the column map

def t_column_map() -> None:
    aln = {"a": "MKV-QT", "b": "MK--QT", "c": "MKVWQT"}
    trimmed = {"a": "MKQT", "b": "MKQT", "c": "MKQT"}
    good = [{"trimmed_column": "1", "aln_column": "1"},
            {"trimmed_column": "2", "aln_column": "2"},
            {"trimmed_column": "3", "aln_column": "5"},
            {"trimmed_column": "4", "aln_column": "6"}]
    out = L.verify_column_map(aln, trimmed, good)
    check("T1 a correct map is accepted", out["columns_checked"] == 4)

    # T2 the failure with no other symptom: every column still maps to a
    # column, and the alignment still parses. Only the walk catches it.
    off = [{"trimmed_column": "1", "aln_column": "2"},
           {"trimmed_column": "2", "aln_column": "3"},
           {"trimmed_column": "3", "aln_column": "5"},
           {"trimmed_column": "4", "aln_column": "6"}]
    check("T2 an off-by-one map is refused",
          _raises(L.verify_column_map, aln, trimmed, off))

    short = good[:3]
    check("T3 a map of the wrong length is refused",
          _raises(L.verify_column_map, aln, trimmed, short))

    # The duplicate case has to be built where the two input columns hold
    # the *same* residues, or the content walk catches it first and the
    # one-to-one guard is never exercised (a mutation test found exactly
    # that: disabling the guard left T4 green).
    aln2 = {"a": "MKQQT", "b": "MKQQT", "c": "MKQQT"}
    trim2 = {"a": "MKQQ", "b": "MKQQ", "c": "MKQQ"}
    dup = [{"trimmed_column": "1", "aln_column": "1"},
           {"trimmed_column": "2", "aln_column": "2"},
           {"trimmed_column": "3", "aln_column": "3"},
           {"trimmed_column": "4", "aln_column": "3"}]
    check("T4 a map with a repeated input column is refused",
          _raises(L.verify_column_map, aln2, trim2, dup))

    check("T5 mismatched label sets are refused",
          _raises(L.verify_column_map, aln, {"a": "MKQT"}, good))


# --------------------------------------------- T6-T9: residues at columns

def t_residues() -> None:
    index = {"ITPR1": {10: {"aa": "R"}, 20: {"aa": "K"}},
             "ITPR2": {11: {"aa": "R"}}}
    good_v = [{"gene": "ITPR1", "resi": "10", "ref_aa": "R"}]
    good_p = [{"gene": "ITPR1", "resi": "10", "other": "ITPR2",
               "other_resi": "11", "other_aa": "R", "aligned": "True"}]
    out = L.verify_residues(good_v, good_p, index)
    check("T6 a correct residue join is accepted",
          out["variant_residues_checked"] == 1
          and out["aligned_partners_checked"] == 1)

    bad_v = [{"gene": "ITPR1", "resi": "10", "ref_aa": "K"}]
    check("T7 a variant whose reference residue is not there is refused",
          _raises(L.verify_residues, bad_v, [], index))

    missing = [{"gene": "ITPR1", "resi": "99", "ref_aa": "R"}]
    check("T8 a variant at a position the table lacks is refused",
          _raises(L.verify_residues, missing, [], index))

    bad_p = [{"gene": "ITPR1", "resi": "10", "other": "ITPR2",
              "other_resi": "11", "other_aa": "K", "aligned": "True"}]
    check("T9 a partner residue the other table contradicts is refused",
          _raises(L.verify_residues, good_v, bad_p, index))


# ----------------------------------- T10-T12: the structure-numbering rule

def t_numbering() -> None:
    checks = [S.numbering_check(S.PAINTED / p, g)
              for g, p in S.REFERENCE_FILE.items()]
    checks += [S.numbering_check(S.PAINTED / p, g)
               for g, p in S.MODEL_FILE.items()]
    passed = [c for c in checks if c["carries_human_numbering"]]
    failed = [c for c in checks if not c["carries_human_numbering"]]
    # A rule that admitted everything would be no rule, and one that admitted
    # nothing would make the panel unreachable. Both halves are required.
    check("T10 the numbering rule admits at least one structure",
          len(passed) >= 1, f"{len(passed)} passed")
    check("T11 the numbering rule refuses at least one structure",
          len(failed) >= 1, f"{len(failed)} refused")
    check("T12 a refused structure is refused on a measured fraction",
          all(0.0 <= c["frac_agree"] < 1.0 for c in failed))


# ------------------------------------------------- T13-T15: figure inputs

def t_inputs() -> None:
    reps = L.read_tsv(L.MSA_DIR / "representatives.tsv")
    groups = {r["group"] for r in reps}
    check("T13 every representative's group has a committed colour",
          groups <= set(F.GROUP), f"unknown: {sorted(groups - set(F.GROUP))}")

    tab = L.constraint_table("ITPR1")
    elements = {r["element"] for r in tab}
    known = set(L.ELEMENT_ORDER) | {e for e in elements
                                    if e.startswith("linker") or
                                    e in ("nterm", "cterm", "unassigned")}
    check("T14 every element has a place in the reading order or is a linker",
          elements <= known, f"unknown: {sorted(elements - known)}")

    atoms = L.read_ca_pdb(S.PAINTED / S.REFERENCE_FILE["ITPR3"])
    u, v = L.principal_axes(atoms)
    u2, v2 = L.principal_axes(atoms)
    check("T15 the projection is deterministic", u == u2 and v == v2)


def t_reproducible() -> None:
    """T16 — two saves of one figure must be byte-identical, pdf included.

    matplotlib stamps the wall clock into a PDF's /CreationDate, so before
    `figstyle.save` dropped that key every figure in the project differed
    from its own rebuild by two bytes and no SHA-256 recorded against a pdf
    meant anything. The property is cheap to test and silent when it breaks.
    """
    import hashlib
    import tempfile
    import matplotlib.pyplot as plt

    with tempfile.TemporaryDirectory() as tmp:
        digests = []
        for _ in range(2):
            fig, ax = plt.subplots(figsize=(2, 1))
            ax.plot([0, 1], [0, 1])
            paths = F.save(fig, Path(tmp) / "probe")
            plt.close(fig)
            digests.append([hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in paths])
        check("T16 a figure saves byte-identically twice, pdf included",
              digests[0] == digests[1])
        # Byte-identity alone passes vacuously when both saves land in the
        # same second, which is most of the time and is exactly how a test
        # like this goes quiet (a mutation test found it). The absence of
        # the key is the property; check it directly.
        pdf = (Path(tmp) / "probe.pdf").read_bytes()
        check("T17 the saved pdf carries no creation timestamp",
              b"/CreationDate" not in pdf)


def self_test() -> int:
    FAILURES.clear()
    CHECKS.clear()
    print("[s24] negative controls")
    t_column_map()
    t_residues()
    t_numbering()
    t_inputs()
    t_reproducible()
    if FAILURES:
        print(f"  {len(FAILURES)} FAILED")
        for f in FAILURES:
            print(f"    - {f}")
        return 1
    print("  all passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(self_test())
