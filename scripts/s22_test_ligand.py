"""Constructed negative controls for S22, run before anything is written.

Every rule in this task returns a plausible number when it is wrong, and
four of them return a number that looks *better* when it is wrong:

* a ligand core drawn without checking that it holds the contacts still
  produces a clean core-versus-pore p-value;
* a per-tip paired test that admits a tip covering the pore and not the
  core reports the N-terminal truncation of a gene model as an extreme
  ligand-core divergence;
* a permutation that resamples the *values* instead of the label answers a
  question nobody asked and answers it significantly;
* and a lineage panel that files "this species has no reference proteome"
  as "this species has no PLC" manufactures the test set out of missing
  data.

So most of what follows is a test on refusal, and three of them (T5, T14,
T16) are tests on *reachability* — a rule that can only ever return the
answer S22 reports is not a measurement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402
import s22_shells as SH  # noqa: E402

FAILS: list[str] = []
_BUILT: list[dict] | None = None


def built() -> list[dict]:
    """The module map as the rules produce it *now*.

    T1-T3 originally read the committed `module_map.tsv`, which a broken
    builder would leave untouched from an earlier run — so a rule change
    that silently widened a module passed every check.  The rules are
    therefore exercised here and the committed file is compared against
    them (T3b), which is D13's discipline turned into a test.
    """
    global _BUILT
    if _BUILT is None:
        _BUILT = M.build()
    return _BUILT


def check(name: str, ok: bool, detail: str = "") -> None:
    status = "ok  " if ok else "FAIL"
    print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# --- T1-T3  the module definitions refuse a region that is not the module --

def t1_core_must_hold_every_contact() -> None:
    core = M.residues(built(), "ITPR3", "contact_span")
    contacts = {r["resi"] for r in L.functional_sites()
                if r["paralog"] == "ITPR3" and r["site_class"] == "ip3_contact"}
    check("T1 primary ligand core holds all ten measured contacts",
          contacts <= core, f"{len(contacts & core)}/10")
    shrunk = {i for i in core if i > min(contacts)}
    check("T1b a core missing one contact is not the primary core",
          not contacts <= shrunk)


def t1c_the_builder_itself_is_exercised() -> None:
    """The rule is called rather than its output read (see `built`)."""
    try:
        rows = built()
    except M.ModuleRefusal as exc:
        check("T1c the module builder runs on this project's own data",
              False, str(exc))
        return
    prim = [r for r in rows if r["is_primary"]]
    check("T1c the module builder runs and returns both primary modules",
          len(prim) == 6, f"{len(prim)} primary rows")


def t1d_the_definition_check_refuses() -> None:
    for module, args, why in (
            ("ligand_core", (9, 0, 0), "a core missing one contact"),
            ("ligand_core", (10, 1, 0), "a core holding a filter residue"),
            ("pore_module", (0, 2, 1), "a pore missing a gate residue"),
            ("pore_module", (1, 2, 2), "a pore holding an IP3 contact")):
        try:
            M.check_definition(module, *args, "constructed")
            check(f"T1d {why} is refused", False, "it was accepted")
        except M.ModuleRefusal:
            check(f"T1d {why} is refused", True)


def t2_pore_must_hold_filter_and_gate() -> None:
    pore = M.residues(built(), "ITPR3", "channel_minus_luminal")
    sites = {r["resi"] for r in L.functional_sites()
             if r["paralog"] == "ITPR3"
             and r["site_class"] in ("filter_lining", "gate_lining")}
    check("T2 primary pore module holds both filter and both gate residues",
          sites <= pore, f"{len(sites & pore)}/4")
    luminal = next((r for r in L.domain_map()
                    if r["paralog"] == "ITPR3" and r["element"] == "luminal_loop"),
                   None)
    inside = sum(1 for i in range(luminal["start"], luminal["end"] + 1)
                 if i in pore)
    check("T2b the luminal loop is excluded from the primary pore module",
          inside == 0, f"{inside} loop residues inside")


def t3_modules_are_disjoint() -> None:
    rows = built()
    bad = [p for p in L.PARALOGS
           if M.residues(rows, p, M.PRIMARY["ligand_core"])
           & M.residues(rows, p, M.PRIMARY["pore_module"])]
    check("T3 the two primary modules never overlap", not bad, str(bad))
    path = L.OUT_DIR / "module_map.tsv"
    if not path.exists():
        print("  [skip] T3b module_map.tsv not written yet")
        return
    committed = {(r["paralog"], r["definition"]): (r["start"], r["end"],
                                                   r["excluded"])
                 for r in L.read_tsv(path)}
    fresh = {(r["paralog"], r["definition"]): (str(r["start"]), str(r["end"]),
                                               r["excluded"])
             for r in rows}
    check("T3b the committed module map is what the rules produce now",
          committed == fresh,
          str(sorted(k for k in fresh if committed.get(k) != fresh[k])))


# --- T4-T6  the structure reader -------------------------------------------

CIF_HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""


def _cif(lines: list[str], tmp: Path) -> Path:
    tmp.write_text(CIF_HEADER + "\n".join(lines) + "\n#\n")
    return tmp


def _atom(i, atom, alt, comp, resi, chain, x, y, z, model=1):
    return (f"{'HETATM' if comp == 'I3P' else 'ATOM'}   {i} C {atom} {alt} "
            f"{comp} A 1 {resi} ? {x} {y} {z} 1.000 30.0 ? {resi} {comp} "
            f"{chain} {atom} {model}")


def t4_reader_ignores_second_model_and_altloc(tmp: Path) -> None:
    p = _cif([
        _atom(1, "CA", ".", "ARG", 10, "A", 0.0, 0.0, 0.0),
        _atom(2, "CA", "B", "ARG", 11, "A", 0.0, 0.0, 0.0),
        _atom(3, "CA", ".", "ARG", 12, "A", 0.0, 0.0, 0.0, model=2),
        _atom(4, "P", ".", "I3P", 900, "A", 1.0, 0.0, 0.0),
    ], tmp / "t4.cif")
    prot, lig = SH.read_atoms(p)
    resis = {r[1] for r in prot}
    check("T4 second model and altloc B are both skipped",
          resis == {10}, str(sorted(resis)))
    check("T4b the ligand is read as ligand and not as protein", len(lig) == 1)


def t5_all_atom_beats_ca(tmp: Path) -> None:
    """A residue whose CA is far and whose side chain is close is a contact."""
    p = _cif([
        _atom(1, "CA", ".", "ARG", 10, "A", 10.0, 0.0, 0.0),
        _atom(2, "NH1", ".", "ARG", 10, "A", 3.0, 0.0, 0.0),
        _atom(3, "P", ".", "I3P", 900, "A", 0.0, 0.0, 0.0),
    ], tmp / "t5.cif")
    prot, lig = SH.read_atoms(p)
    d = SH.min_distances(prot, lig, cutoff=SH.SEARCH_RADIUS_A)
    got = d[("A", 10)]["d_same"]
    check("T5 an all-atom distance finds a contact a CA trace would miss",
          abs(got - 3.0) < 1e-6 and SH.shell_of(got) == "contact",
          f"d={got:.2f}")


def t6_cross_subunit_is_separated(tmp: Path) -> None:
    p = _cif([
        _atom(1, "CA", ".", "ARG", 10, "A", 0.0, 0.0, 0.0),
        _atom(2, "P", ".", "I3P", 900, "B", 2.0, 0.0, 0.0),
        _atom(3, "P", ".", "I3P", 901, "A", 9.0, 0.0, 0.0),
    ], tmp / "t6.cif")
    prot, lig = SH.read_atoms(p)
    d = SH.min_distances(prot, lig, cutoff=SH.SEARCH_RADIUS_A)["A", 10]
    check("T6 same-subunit and any-subunit distances are kept apart",
          abs(d["d_same"] - 9.0) < 1e-6 and abs(d["d_any"] - 2.0) < 1e-6,
          f"same={d['d_same']:.1f} any={d['d_any']:.1f}")


def t7_search_radius_is_enforced(tmp: Path) -> None:
    p = _cif([
        _atom(1, "CA", ".", "ARG", 10, "A", SH.SEARCH_RADIUS_A + 5.0, 0.0, 0.0),
        _atom(2, "P", ".", "I3P", 900, "A", 0.0, 0.0, 0.0),
    ], tmp / "t7.cif")
    prot, lig = SH.read_atoms(p)
    d = SH.min_distances(prot, lig, cutoff=SH.SEARCH_RADIUS_A)
    check("T7 a residue past the search radius is absent, not binned",
          ("A", 10) not in d or d[("A", 10)]["d_same"] > SH.SEARCH_RADIUS_A)
    check("T7b shell_of names it beyond_radius rather than a shell",
          SH.shell_of(SH.SEARCH_RADIUS_A + 1) == "beyond_radius")


# --- T8-T10  the statistics ------------------------------------------------

def t8_permutation_resamples_the_label() -> None:
    vals = [float(i) for i in range(100)]
    top = [v >= 90 for v in vals]
    hi = L.permutation_label(vals, top, iters=20000)
    check("T8 the label permutation finds an extreme label extreme",
          hi["p"] is not None and hi["p"] < 0.001, f"p={hi['p']}")
    mid = [40 <= v < 50 for v in vals]
    md = L.permutation_label(vals, mid, iters=20000)
    check("T8b a label sitting at the median is not significant",
          md["p"] is not None and md["p"] > 0.2, f"p={md['p']}")
    empty = L.permutation_label(vals, [False] * 100, iters=100)
    check("T8c an empty label is refused rather than scored",
          empty["p"] is None)


def t9_sign_test_drops_and_counts_ties() -> None:
    r = L.sign_test([1.0, 1.0, -1.0, 0.0, 0.0])
    check("T9 ties are dropped from n and counted",
          r["n"] == 3 and r["n_ties"] == 2, str(r))
    r2 = L.sign_test([1.0] * 10)
    check("T9b ten one-way pairs reach p < 0.01", r2["p"] < 0.01,
          f"p={r2['p']}")


def t10_bh_is_monotone_and_passes_none() -> None:
    ps = [0.001, 0.01, None, 0.5, 0.2]
    q = L.benjamini_hochberg(ps)
    check("T10 None passes through and takes no rank", q[2] is None)
    check("T10b no q-value falls below its own p-value",
          all(q[i] >= ps[i] - 1e-12 for i in (0, 1, 3, 4)),
          str([q[i] for i in (0, 1, 3, 4)]))
    order = sorted((i for i in (0, 1, 3, 4)), key=lambda i: ps[i])
    seq = [q[i] for i in order]
    check("T10c q-values are non-decreasing in p order",
          all(b >= a - 1e-12 for a, b in zip(seq, seq[1:])), str(seq))
    check("T10d the correction is applied on the four scored rows, not five",
          abs(q[0] - 0.001 * 4 / 1) < 1e-12, f"q0={q[0]}")


def main() -> int:
    import tempfile
    print("[s22 self-test] constructed negative controls")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        t1_core_must_hold_every_contact()
        t1c_the_builder_itself_is_exercised()
        t1d_the_definition_check_refuses()
        t2_pore_must_hold_filter_and_gate()
        t3_modules_are_disjoint()
        t4_reader_ignores_second_model_and_altloc(tmp)
        t5_all_atom_beats_ca(tmp)
        t6_cross_subunit_is_separated(tmp)
        t7_search_radius_is_enforced(tmp)
    t8_permutation_resamples_the_label()
    t9_sign_test_drops_and_counts_ties()
    t10_bh_is_monotone_and_passes_none()
    import s22_test_lineage as TL
    TL.run(check)
    if FAILS:
        print(f"[s22 self-test] {len(FAILS)} FAILED: {FAILS}")
        return 1
    print("[s22 self-test] all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
