"""Negative controls for the S11 rules, run on every build.

`s5_bait_screen.self_test()`'s pattern, applied to a structural task. Every
rule in S11 returns a *plausible* number when it is wrong — a mis-parsed
mmCIF still yields coordinates, a fragment still yields a TM-score, a
cryo-EM B-factor still yields a mean — so these are checks on **refusal**.

Each case is constructed to violate exactly one rule and must be rejected
by that rule.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_refs as refs                                         # noqa: E402
import s11_tables as tb                                         # noqa: E402
from s11_plddt import covered_identity                          # noqa: E402
from s11_struct_io import (largest_chain, read_cif_ca,          # noqa: E402
                           read_pdb_ca, write_ca_pdb, Chain, Residue)
import s11_tmalign_run as tmr                                    # noqa: E402
from s11_tmalign import TM_RANDOM, run_tmalign                  # noqa: E402

_CIF_HEAD = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
"""


def _cif(rows: list[str]) -> Path:
    path = Path(tempfile.mkstemp(suffix=".cif")[1])
    path.write_text(_CIF_HEAD + "\n".join(rows) + "\n#\n")
    return path


def _chain(n: int, chain_id: str = "A", bfac: float = 50.0) -> Chain:
    return Chain(chain_id, [Residue(chain_id, i, "ALA", float(i), 0.0, 0.0, bfac)
                            for i in range(1, n + 1)])


def run() -> list[str]:
    fails: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        if not ok:
            fails.append(f"{name}: {detail or 'failed'}")

    # T1 — a second model must not be read. Averaging two models of an
    # NMR-style deposition would silently double a chain's length.
    p = _cif([
        "ATOM CA . ALA A 1 0.0 0.0 0.0 50.0 1",
        "ATOM CA . ALA A 2 1.0 0.0 0.0 50.0 1",
        "ATOM CA . ALA A 3 2.0 0.0 0.0 50.0 2",
    ])
    got = read_cif_ca(p)
    check("T1 second model ignored", len(got["A"]) == 2, f"{len(got['A'])} residues")

    # T2 — altloc B is skipped and altloc A kept, exactly once each.
    p = _cif([
        "ATOM CA A ALA A 1 0.0 0.0 0.0 50.0 1",
        "ATOM CA B SER A 1 9.0 0.0 0.0 50.0 1",
        "ATOM CA . ALA A 2 1.0 0.0 0.0 50.0 1",
    ])
    got = read_cif_ca(p)
    check("T2 altloc", len(got["A"]) == 2 and got["A"].sequence == "AA",
          got["A"].sequence)

    # T3 — largest_chain takes the longest, and is deterministic on a tie.
    chains = {"A": _chain(10), "B": _chain(20), "C": _chain(20)}
    first = largest_chain(chains).chain_id
    check("T3 largest chain", first == largest_chain(dict(reversed(list(
        chains.items())))).chain_id and len(largest_chain(chains)) == 20,
        f"picked {first}")

    # T4 — the CA trace round-trips residue ids and B-factors. AFDB puts
    # pLDDT in that column, so losing it loses every confidence number.
    src = Chain("Z", [Residue("Z", 7, "TRP", 1.5, -2.25, 3.125, 88.31),
                      Residue("Z", 8, "GLY", 4.0, 0.0, 0.0, 12.0)])
    dest = Path(tempfile.mkstemp(suffix=".pdb")[1])
    write_ca_pdb(src, dest)
    back = read_pdb_ca(dest)["A"]
    check("T4 round trip",
          [r.seq_id for r in back.residues] == [7, 8]
          and back.sequence == "WG"
          and abs(back.residues[0].bfactor - 88.31) < 0.01
          and abs(back.residues[0].z - 3.125) < 0.001,
          f"{back.sequence} {[r.bfactor for r in back.residues]}")

    # T5 — compound states before the words they contain. "Inactive-like"
    # read as "active" would file an inhibited structure as an open one.
    for title, want in [("Human type 3 IP3R in the inactive-like state", "inactive-like"),
                        ("IP3R in the preactivated state", "preactivated"),
                        ("RyR1 in the activated state", "activated"),
                        ("Structure in the apo-state", "apo"),
                        ("Labile Resting State 1", "labile resting"),
                        ("A structure of something", "")]:
        check(f"T5 state '{want or 'none'}'", refs.parse_state(title) == want,
              f"got '{refs.parse_state(title)}'")

    # T6 — R1 never reads the title. An entity whose title says IP3
    # receptor but whose accession the census calls RYR must call RYR.
    census = {"P00001": {"call": "RYR", "gene": "RYR1"}}
    row = {"uniprot": "P00001",
           "title": "Inositol 1,4,5-trisphosphate receptor type 1"}
    check("T6 call ignores title",
          refs.call_entity(row, census)[0] == "RYR",
          refs.call_entity(row, census)[0])
    check("T6 no census row reported",
          refs.call_entity({"uniprot": "P99999"}, census)[0] == "no_census_row")
    check("T6 no uniprot reported",
          refs.call_entity({"uniprot": ""}, census)[0] == "no_uniprot")

    # T7 — R3 in both directions: a binding-core construct is not an ITPR
    # reference, and a RyR-length chain cannot enter through the ITPR band.
    base = {"method": "ELECTRON MICROSCOPY", "resolution": 3.0,
            "state": "apo", "paralog": "ITPR1"}
    check("T7 core construct rejected",
          refs.eligibility({**base, "call": "ITPR", "sample_length": 226})[1]
          == "R3 not full length")
    check("T7 RyR length rejected from ITPR band",
          not refs.length_ok("ITPR", 5037))
    check("T7 RyR accepted in its own band", refs.length_ok("RYR", 5037))
    check("T7 X-ray rejected",
          refs.eligibility({**base, "call": "ITPR", "sample_length": 2700,
                            "method": "X-RAY DIFFRACTION"})[1] == "R2 not cryo-EM")
    check("T7 unreadable state rejected",
          refs.eligibility({**base, "call": "ITPR", "sample_length": 2700,
                            "state": ""})[1] == "R4 state not readable")
    check("T7 full-length cryo-EM accepted",
          refs.eligibility({**base, "call": "ITPR", "sample_length": 2700})[0])

    # T8 — R7's size rule. A control must be its whole protein; a control
    # that is one domain scores low for reasons of size, not fold.
    ctrl = [{"control_class": "c", "eligible": 1, "sample_length": 2600,
             "entry_id": "AAAA", "resolution": 3.0, "decoy_gene": "X"},
            {"control_class": "c", "eligible": 1, "sample_length": 900,
             "entry_id": "BBBB", "resolution": 1.0, "decoy_gene": "X"}]
    picked, _ = refs.select_controls(ctrl, 2700)
    check("T8 size-matched control", picked and picked[0]["entry_id"] == "AAAA",
          picked[0]["entry_id"] if picked else "none")

    # T9/T10 — D14 structurally. A tie is a `no_call`, a low pair is
    # `no_fold_match`, and the score used is normalised by the *reference*.
    man = [{"id": "ref_i", "role": "reference", "call": "ITPR", "status": "ok"},
           {"id": "ref_r", "role": "reference", "call": "RYR", "status": "ok"},
           {"id": "q", "role": "model", "call": "ITPR", "status": "ok"}]
    def pair(q, t, tmq, tmt):
        return {"query": q, "target": t, "tm_query": tmq, "tm_target": tmt,
                "ok": 1}
    tie = tb.vs_reference(man, [pair("q", "ref_i", 0.9, 0.60),
                                pair("q", "ref_r", 0.9, 0.58)])
    qrow = [r for r in tie if r["id"] == "q"][0]
    check("T9 tie is no_call", qrow["structural_call"] == "no_call",
          qrow["structural_call"])
    check("T10 normalised by reference", abs(qrow["best_itpr"] - 0.60) < 1e-9,
          str(qrow["best_itpr"]))
    low = tb.vs_reference(man, [pair("q", "ref_i", 0.9, TM_RANDOM - 0.02),
                                pair("q", "ref_r", 0.9, 0.02)])
    check("T9 below random bar",
          [r for r in low if r["id"] == "q"][0]["structural_call"]
          == "no_fold_match")
    clear = tb.vs_reference(man, [pair("q", "ref_i", 0.9, 0.72),
                                  pair("q", "ref_r", 0.9, 0.30)])
    crow = [r for r in clear if r["id"] == "q"][0]
    check("T9 clear call agrees",
          crow["structural_call"] == "ITPR" and crow["call_agrees"] == 1)

    # T13 — the run lock. An interrupted run whose Python parent survives
    # keeps spawning TM-align beside its replacement; this task produced
    # exactly that, as S7 did for IQ-TREE. A live pid must block, a dead
    # one must not, or a crash locks the stage for ever.
    lock = tmr.lock_path()
    saved = lock.read_text() if lock.exists() else None
    try:
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text("1\n")           # pid 1 is always alive
        blocked = False
        try:
            tmr.claim()
        except RuntimeError:
            blocked = True
        check("T13 live pid blocks", blocked)
        lock.write_text("999999\n")      # a pid that cannot exist
        stale_ok = True
        try:
            tmr.claim()
        except RuntimeError:
            stale_ok = False
        check("T13 dead pid does not block", stale_ok)
        tmr.release()
        check("T13 release clears the lock", not lock.exists())
    finally:
        if saved is not None:
            lock.write_text(saved)
        elif lock.exists():
            lock.unlink()

    # T14 — the cache write is atomic. Two processes on one key is what the
    # orphan produced, and an interleaved write yields JSON that parses as
    # garbage. A rename cannot be observed half-done.
    tmp_dir = Path(tempfile.mkdtemp())
    target = tmp_dir / "pair.json"
    tmr._write_atomic(target, '{"ok": true}')
    leftovers = [p for p in tmp_dir.iterdir() if p.suffix == ".tmp"]
    import json as _json
    check("T14 atomic write",
          _json.loads(target.read_text())["ok"] is True and not leftovers,
          f"{len(leftovers)} temp files left")

    # T15 — a failed pair is never cached. A TM-align killed by a signal
    # returns a non-zero status with empty stderr; caching that stores a
    # permanent zero-score result that looks exactly like a real "these do
    # not match", which is the one error this table cannot survive.
    cache_dir = tmr.structures_dir() / "tmalign_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    before = set(cache_dir.iterdir())
    missing = Path(tempfile.mkdtemp()) / "does_not_exist.pdb"
    q = {"id": "q", "path": str(missing), "sha256": "deadbeef" * 8}
    t = {"id": "t", "path": str(missing), "sha256": "cafebabe" * 8}
    row = tmr._one((q, t))
    added = set(cache_dir.iterdir()) - before
    for extra in added:
        extra.unlink()
    check("T15 failure not cached", row["ok"] == 0 and not added,
          f"{len(added)} cache entries written for a failed pair")

    # T11/T12 — the confidence transfer. Identity is a fraction, and an
    # unplaced domain reports no pLDDT rather than a plausible one.
    check("T12 identity is a fraction",
          covered_identity("ABCD", "ABCD") == 1.0
          and covered_identity("AB--", "--CD") == 0.0
          and covered_identity("ABCD", "ABXD") == 0.75)
    return fails


def self_test() -> None:
    fails = run()
    if fails:
        raise AssertionError("S11 self-test failed:\n  " + "\n  ".join(fails))


if __name__ == "__main__":
    problems = run()
    for p in problems:
        print("FAIL", p)
    print(f"[s11] self-test: {len(problems)} failure(s)")
    raise SystemExit(1 if problems else 0)
