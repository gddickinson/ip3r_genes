"""Negative controls for the S19 rules, run before anything is written.

The pattern `s5_bait_screen.self_test()` established, applied to a task whose
failure modes are unusually quiet. Every rule in S19 returns a plausible
number when it is wrong: an empty accession universe makes every record look
absent from the database it was found in, a panel simulation that does not
reproduce the ledger still produces a tidy ablation table, and a control
series built from the wrong states measures the search's agreement with
itself.

So these are tests on **refusal** and on **reachability**: each rule must
fire on its own violation, and each measurement must be able to come out the
other way.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402

FAILURES: list[str] = []
SKIPPED: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    FAILURES.append(f"{name}: {detail}") if not ok else None
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f" — {detail}"
                                                    if detail and not ok
                                                    else ""))


def have(name: str, table: str) -> bool:
    """Is the committed table this check needs on disk?

    The suite runs twice: once **before** any stage writes, where only the
    constructed-input checks can run, and once after, where the artefact
    checks can. A missing table is therefore a skip and not a failure — but
    it is printed as a skip, because a check that quietly does not run is
    worse than one that fails.
    """
    if (S.OUT_DIR / table).exists():
        return True
    SKIPPED.append(f"{name} (needs {table})")
    print(f"  skip {name} — {table} not written yet")
    return False


# ------------------------------------------------------------------ T1-T3 lib

def t1_header_accession() -> None:
    """A FASTA header must yield a bare accession, in all three shapes."""
    check("T1a sp| header", S.header_accession(">sp|Q14643|ITPR1_HUMAN x")
          == "Q14643")
    check("T1b tr| header", S.header_accession(">tr|A0A1|A0A1_MOUSE") == "A0A1")
    check("T1c bare header", S.header_accession(">Q14643 something")
          == "Q14643")


def t2_acc_key_leaves_model_ids() -> None:
    """A genome model id must survive `acc_key` intact.

    Stripping at the first `.` would collapse every model of one assembly
    onto its accession prefix, turning 1,058 gene models into a few hundred
    records — a silent 3x error in every contribution count.
    """
    mid = "GCA_000699545.1|ITPR1|JMFS01131396.1:4289-29317+"
    check("T2a model id untouched", S.acc_key(mid) == mid)
    check("T2b isoform stripped", S.acc_key("Q14643-4") == "Q14643")
    check("T2c version stripped", S.acc_key("Q14643.2") == "Q14643")


def t3_universe_refuses_empty(tmp: Path) -> None:
    """A zero-row parse must refuse to cache.

    An empty universe file is indistinguishable from a database that holds
    none of the query, so every head-to-head cell would read "not in the
    database" and the enumeration would appear to have found nothing the
    profile could have.
    """
    fake = tmp / "empty_refprot.fasta"
    fake.write_text("no headers here\njust sequence\n")
    real = S.data_root
    try:
        S.data_root = lambda: tmp                              # type: ignore
        (tmp / "proteomes").mkdir(parents=True, exist_ok=True)
        (tmp / "proteomes" / "empty_refprot.fasta").write_text(
            fake.read_text())
        (tmp / "raw_api" / "s19").mkdir(parents=True, exist_ok=True)
        try:
            S.build_universe("empty", force=True)
            check("T3 empty universe refused", False, "it cached anyway")
        except SystemExit:
            check("T3 empty universe refused", True)
    finally:
        S.data_root = real                                     # type: ignore


# ---------------------------------------------------------------- T4-T6 control

def t4_control_excludes_undecidable() -> None:
    """A `paralog_unassignable` cell must be in neither control series.

    S15a declined to call those four cyclostome cells present or absent.
    Counting them as present would put four guaranteed false negatives into
    the numerator; counting them as absent would assert a loss S15b did not
    reconstruct. The rule must refuse both.
    """
    cells = [
        {"cell": "ITPR2", "s15_state": "paralog_unassignable", "found": 0},
        {"cell": "ITPR1", "s15_state": "present_single_locus", "found": 1},
        {"cell": "RYR", "s15_state": "", "found": 0},
    ]
    out = S.control_cells(cells)
    kinds = {c["cell"]: c["control"] for c in out}
    check("T4a unassignable excluded", "ITPR2" not in kinds)
    check("T4b present cell in the ITPR series",
          kinds.get("ITPR1") == "itpr_present")
    check("T4c RyR cell is its own series",
          kinds.get("RYR") == "ryr_sister")


def t5_false_negative_is_reachable() -> None:
    """A false negative must be constructible, or the rate is not a measure.

    A control definition that can only ever produce zero misses is not a
    measurement of sensitivity. This builds a cell that is demonstrably
    present and demonstrably not found, and requires it to be counted.
    """
    import s19_contiguity as C
    cells = [{"accession": "GCA_TEST", "organism": "Test test",
              "vclass": "Aves", "cell": "ITPR1", "status": "assembly_gap",
              "found": 0, "s15_state": "present_fragmented", "s15_rule": "R3",
              "level": "Scaffold", "annotated_assembly": "N",
              "annotation_source": "", "contig_n50": 12000,
              "scaffold_n50": 20000, "total_length_bp": 1_000_000_000,
              "best_coverage": 0.3, "recon_coverage": 0.8}]
    grade = C.grade(cells[0]["status"])
    ctrl = S.control_cells(cells)
    check("T5a graded as a negative", grade == "assembly_gap", grade)
    check("T5b counted as a false negative",
          bool(ctrl) and not ctrl[0]["found"])


def t6_floor_scan_is_monotone_in_scope() -> None:
    """Raising the floor must never retain more genomes than a lower one."""
    import s19_floor as C
    check("T6b D4's own bar is on the grid",
          S.contiguity_bar() in C.FLOOR_GRID)
    if not have("T6a scope falls as the floor rises", "absence_floor.tsv"):
        return
    rows = S.read_tsv(S.OUT_DIR / "absence_floor.tsv")
    ok = True
    for series in {r["series"] for r in rows}:
        mine = sorted((r for r in rows if r["series"] == series),
                      key=lambda r: float(r["contig_n50_floor"]))
        kept = [int(r["n_genomes_kept"]) for r in mine]
        ok = ok and all(a >= b for a, b in zip(kept, kept[1:]))
    check("T6a scope falls as the floor rises", ok)


# ----------------------------------------------------------------- T7-T9 panel

def t7_panel_reproduces_ledger() -> None:
    """The full panel must reproduce the committed ledger cell for cell.

    This is the test the whole ablation rests on: if the simulation does not
    reproduce the sweep it is a different instrument, and every ablation
    delta is measured against the wrong baseline.
    """
    if not have("T7 full panel == ledger", "panel_validation.tsv"):
        return
    rows = S.read_tsv(S.OUT_DIR / "panel_validation.tsv")
    check("T7 full panel == ledger", not rows,
          f"{len(rows)} cells disagree")


def t8_unlabelled_bait_cannot_fill_a_cell() -> None:
    """A locus with no labelled paralog bait must fill no paralog cell.

    `s5_classify.cell_loci` offers an unresolved-clade locus to the paralog
    whose bait scores highest *there*, so with no labelled bait present there
    is no paralog to offer it to. A simulation that filled the cell anyway
    would credit the unlabelled baits with recall they cannot have.
    """
    import s5_classify as CL
    import s5_sweep_lib as SW
    aln = SW.Aln(contig="c1", start=1, end=1000, strand="+", score=900.0,
                 identity=0.9, bait="basal_bait", clade="vertebrate_basal",
                 family="ITPR", q_start=1, q_end=900, bait_len=900,
                 aligned_aa=900)
    loci = SW.cluster_loci([aln])
    primary, secondary = CL.cell_loci(loci, "ITPR1")
    check("T8a no primary locus", not primary)
    check("T8b no secondary offer without a labelled bait", not secondary)
    labelled = SW.Aln(contig="c1", start=1, end=1000, strand="+", score=800.0,
                      identity=0.9, bait="itpr1_bait", clade="ITPR1",
                      family="ITPR", q_start=1, q_end=800, bait_len=900,
                      aligned_aa=800)
    loci = SW.cluster_loci([aln, labelled])
    _p, secondary = CL.cell_loci(loci, "ITPR1")
    check("T8c offered once a labelled bait aligns there", len(secondary) == 1)


def t9_ryr_locus_never_enters_an_itpr_cell() -> None:
    """D14 must hold in the simulation as it does in the sweep."""
    import s5_classify as CL
    import s5_sweep_lib as SW
    ryr = SW.Aln(contig="c1", start=1, end=5000, strand="+", score=9000.0,
                 identity=0.9, bait="ryr_bait", clade="RYR", family="RYR",
                 q_start=1, q_end=4900, bait_len=4900, aligned_aa=4900)
    itpr = SW.Aln(contig="c1", start=1, end=5000, strand="+", score=100.0,
                  identity=0.4, bait="itpr1_bait", clade="ITPR1",
                  family="ITPR", q_start=1, q_end=700, bait_len=2700,
                  aligned_aa=700)
    loci = SW.cluster_loci([ryr, itpr])
    primary, secondary = CL.cell_loci(loci, "ITPR1")
    check("T9 RyR-won locus offered to no ITPR cell",
          not primary and not secondary)


# ----------------------------------------------------------------- T10-T12 drift

def t10_accepted_rounds_bound_the_target_set() -> None:
    """A killed run must contribute only its accepted rounds.

    `s20_protista_other` ends with 26,148 targets and D10 accepts 22,913 of
    them at round 9. Reading the domtblout instead would hand a completeness
    argument every round of a run the project disowned.
    """
    import s3_kill as K
    rounds = [{"round": 1, "included": ["a", "b"]},
              {"round": 2, "included": ["a", "b", "c"]},
              {"round": 3, "included": ["a", "b", "c", "d", "e"]}]
    check("T10a accepted set stops at the bar",
          K.accepted_targets(rounds, 2) == {"a", "b", "c"})
    check("T10b whole run when nothing was killed",
          K.accepted_targets(rounds, 3) == {"a", "b", "c", "d", "e"})


def t11_criterion_can_fire_and_can_decline() -> None:
    """The proposed off-family rule must do both, on constructed rounds."""
    import s19_drift as D
    drifting = [["r", "g", 1, 0, 100, 90, 5, 5, 0.05, 0.0, "", 0.05, 0.0, 1],
                ["r", "g", 2, 0, 900, 90, 5, 805, 0.006, -0.044, "9.0",
                 0.894, 0.844, 1]]
    stable = [["s", "g", 1, 0, 100, 90, 5, 5, 0.05, 0.0, "", 0.05, 0.0, 1],
              ["s", "g", 2, 0, 110, 98, 5, 7, 0.045, -0.005, "1.1", 0.064,
               0.014, 1]]
    out = {r[0]: r for r in D.criterion_trace(drifting + stable,
                                              log=lambda *_: None,
                                              write=False)}
    check("T11a fires on a drifting run",
          out["r"][6] == 2, str(out["r"][6]))
    check("T11b declines on a stable run",
          out["s"][6] == "never", str(out["s"][6]))
    check("T11c K1 declines on the drifting run too — D10b",
          out["r"][3] == "never", str(out["r"][3]))


def t12_drift_label_is_measured_not_asserted() -> None:
    """The drift outcome must come from the finished model, and separate."""
    import s19_drift as D
    if not have("T12 drift label measured", "drift_outcome.tsv"):
        return
    rows = S.read_tsv(S.OUT_DIR / "drift_outcome.tsv")
    drift = [float(r["offfamily_frac_final"]) for r in rows
             if r["drifted"] == "1"]
    clean = [float(r["offfamily_frac_final"]) for r in rows
             if r["drifted"] == "0"]
    check("T12a both classes reachable", bool(drift) and bool(clean))
    check("T12b the two populations do not touch",
          bool(drift) and bool(clean) and min(drift) > max(clean),
          f"min drifted {min(drift) if drift else '-'} vs max clean "
          f"{max(clean) if clean else '-'}")
    check("T12c the bar sits in the gap",
          bool(drift) and bool(clean)
          and max(clean) < D.DRIFTED_OFFFAMILY_FRAC < min(drift))


# --------------------------------------------------------- T13-T14 contribution

def t13_head_to_head_uses_called_sets() -> None:
    """The comparison must be on family calls, not raw profile targets.

    18,501 vertebrate targets were scored and 5,130 called; scoring the
    enumeration against the raw list would credit the sweep with 13,371
    records D22's gate declined.
    """
    import s19_contribution as C
    scored, called = C._hmm_called("vertebrata")
    check("T13a called is a strict subset of scored",
          called < scored and bool(called))
    if not have("T13b table counts the called set", "head_to_head.tsv"):
        return
    rows = S.read_tsv(S.OUT_DIR / "head_to_head.tsv")
    vert = [r for r in rows if r["database"] == "vertebrata" and r["n_hmm"]]
    total = sum(int(r["n_hmm"]) for r in vert)
    check("T13b table counts the called set", total <= len(called),
          f"{total} vs {len(called)}")


def t14_scope_split_covers_every_genome() -> None:
    """Every manifest genome must land in exactly one scope class."""
    import s19_recovery as C
    accs = [r["accession"] for r in S.read_tsv(S.MANIFEST)]
    classes = {C._scope_of(a) for a in accs}
    check("T14a only the three classes",
          classes <= {"order representative", "margin species", "both"},
          str(classes))
    check("T14b all three are populated", len(classes) == 3, str(classes))


def _committed_state() -> dict:
    """SHA-256 of every committed S19 table, so the suite can prove it wrote
    nothing. The first build's T11 called the real `criterion_trace`, which
    wrote its two constructed rows over `kill_criterion_trace.tsv`, and the
    report then read 2 jackhmmer runs where there are 7. A test that can
    damage the artefact it tests is worse than no test."""
    import hashlib
    out = {}
    for path in sorted(S.OUT_DIR.glob("*.tsv")) + \
            sorted(S.OUT_DIR.glob("*.json")):
        out[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def t15_suite_writes_nothing(before: dict) -> None:
    after = _committed_state()
    changed = sorted(k for k in before if before[k] != after.get(k))
    check("T15 the suite altered no committed table", not changed,
          ", ".join(changed))


def main() -> int:
    import tempfile
    print("[s19] negative controls")
    before = _committed_state()
    t1_header_accession()
    t2_acc_key_leaves_model_ids()
    with tempfile.TemporaryDirectory() as d:
        t3_universe_refuses_empty(Path(d))
    t4_control_excludes_undecidable()
    t5_false_negative_is_reachable()
    t6_floor_scan_is_monotone_in_scope()
    t7_panel_reproduces_ledger()
    t8_unlabelled_bait_cannot_fill_a_cell()
    t9_ryr_locus_never_enters_an_itpr_cell()
    t10_accepted_rounds_bound_the_target_set()
    t11_criterion_can_fire_and_can_decline()
    t12_drift_label_is_measured_not_asserted()
    t13_head_to_head_uses_called_sets()
    t14_scope_split_covers_every_genome()
    t15_suite_writes_nothing(before)
    if FAILURES:
        print(f"[s19] {len(FAILURES)} FAILURES")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    if SKIPPED:
        print(f"[s19] all negative controls pass; {len(SKIPPED)} skipped "
              f"for want of a table not written yet")
        for s in SKIPPED:
            print(f"   - {s}")
        return 0
    print("[s19] all negative controls pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
