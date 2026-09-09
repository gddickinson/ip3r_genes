"""Constructed negative controls for the S18 rules, run before anything writes.

Every rule in this task returns a plausible number when it is wrong, and three
of them did on the way here: a name rule that scored "RyR/IP3R Homology
associated domain-containing protein" as the database calling an IP3 receptor a
ryanodine receptor, a paralog rule that scored 6,660 records as wrong-paralog
because their names decline to name a paralog at all, and a bait lookup keyed
on the wrong column that made every sequence call come back `no_call`. So the
tests here are mostly tests on **refusal** and on **reachability** — a rule
that cannot fire is not a negative result, and a rule that fires on everything
is not a positive one.

`python scripts/s18_test_audit.py` exits non-zero on the first failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_gff                                                    # noqa: E402
import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402
import s18_corrections as C                                       # noqa: E402
import s18_protein_audit as P                                     # noqa: E402
import s18_zero_hits as Z                                         # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f" — {detail}" if detail
                                                     else ""))
    if not ok:
        FAILURES.append(name)


def _gene(gid, start, end, strand, cds, name="", product="", biotype="protein_coding",
          pseudo=False):
    return {"gene_id": gid, "start": start, "end": end, "strand": strand,
            "name": name or gid, "locus_tag": gid, "biotype": biotype,
            "pseudo": pseudo, "description": product, "cds": s10_gff.merge(cds),
            "exons": s10_gff.merge(cds) or [(start, end)], "mrnas": {}}


# --------------------------------------------------------------------------
def t1_same_strand_only():
    print(" t1_same_strand_only")
    locus = [(1000, 2000), (5000, 6000)]
    win = {"genes": [_gene("g-", 900, 6100, "-", [(1000, 2000), (5000, 6000)],
                           name="ITPR1", product="inositol 1,4,5-trisphosphate "
                                                 "receptor type 1")]}
    plus = R.models_over(win, locus, "+")
    minus = R.models_over(win, locus, "-")
    check("T1a a perfectly covering gene on the other strand is not a model",
          plus == [], f"{len(plus)} models")
    check("T1b the same gene on the locus's own strand is", len(minus) == 1)
    st = R.locus_state(plus, locus, "present", R.COMPLETE_FRAC)
    check("T1c and the locus is therefore unannotated",
          st["state"] == "unannotated", st["state"])


def t2_cds_blocks_not_gene_spans():
    print(" t2_cds_blocks_not_gene_spans")
    # A 200 kb ITPR locus with two exons, and a passenger gene whose *span*
    # runs the whole way across it while its coding sequence sits entirely in
    # the intron. Scored on spans it covers the locus completely; scored on
    # CDS blocks it covers none of it. ITPR introns reach 152 kb (S5's intron
    # calibration), so there is room for exactly this.
    locus = [(1000, 2000), (200_000, 201_000)]
    win = {"genes": [_gene("passenger", 1_500, 200_500, "+",
                           [(50_000, 60_000)], name="SOMEGENE")]}
    over = R.models_over(win, locus, "+")
    check("T2a a passenger whose CDS is all intronic overlaps no coding block",
          over == [], f"{len(over)} models")
    st = R.locus_state(over, locus, "present", R.COMPLETE_FRAC)
    check("T2b so the locus reads unannotated, not complete",
          st["state"] == "unannotated", st["state"])


def t3_every_state_is_reachable():
    print(" t3_every_state_is_reachable")
    locus = [(1000, 2000), (5000, 6000)]          # 2002 coding bp
    itpr = "inositol 1,4,5-trisphosphate receptor type 1"
    cases = {
        "complete": [_gene("g1", 900, 6100, "+", [(1000, 2000), (5000, 6000)],
                           name="ITPR1", product=itpr)],
        # Three pieces, not two: at a bar of 0.5 a two-piece split of a
        # two-exon gene gives each piece exactly half, and `complete` fires
        # first by design. The state is reachable, and what it takes to reach
        # it is worth having written down.
        "split": [_gene("g1", 900, 1400, "+", [(1000, 1400)], name="ITPR1"),
                  _gene("g2", 1500, 2100, "+", [(1500, 2000)], name="ITPR1"),
                  _gene("g3", 4900, 5900, "+", [(5000, 5800)], name="ITPR1")],
        "fragmentary": [_gene("g1", 900, 1400, "+", [(1000, 1400)],
                              name="ITPR1")],
        "noncoding": [_gene("g1", 900, 6100, "+", [], name="ITPR1",
                            biotype="other")],
        "unannotated": [],
    }
    for want, genes in cases.items():
        over = R.models_over({"genes": genes}, locus, "+")
        st = R.locus_state(over, locus, "present", R.COMPLETE_FRAC)
        check(f"T3 {want} is reachable", st["state"] == want,
              f"got {st['state']} ({st['state_reason']})")
    over = R.models_over({"genes": cases["complete"]}, locus, "+")
    for gs, want in (("none", "no_gene_set"), ("unavailable", "gff_unavailable")):
        st = R.locus_state(over, locus, gs, R.COMPLETE_FRAC)
        check(f"T3 {want} overrides even a perfect model",
              st["state"] == want, st["state"])


def t4_split_needs_two_real_pieces():
    print(" t4_split_needs_two_real_pieces")
    locus = [(1000, 2000), (5000, 6000)]
    genes = [_gene("big", 900, 5990, "+", [(1000, 2000), (5000, 5980)],
                   name="ITPR1"),
             _gene("sliver", 5995, 6000, "+", [(5995, 6000)], name="ITPR1")]
    over = R.models_over({"genes": genes}, locus, "+")
    st = R.locus_state(over, locus, "present", 0.999)
    check("T4a a 6 bp sliver is not a piece",
          st["n_coding_pieces"] == 1,
          f"{st['state']} pieces={st['n_coding_pieces']}")
    check("T4a2 so the locus is not a split", st["state"] != "split",
          st["state"])
    genes2 = [_gene("a", 900, 2100, "+", [(1000, 2000)], name="ITPR1"),
              _gene("b", 4900, 6100, "+", [(5000, 6000)], name="ITPR1")]
    over2 = R.models_over({"genes": genes2}, locus, "+")
    st2 = R.locus_state(over2, locus, "present", 0.999)
    check("T4b two halves are", st2["state"] == "split"
          and st2["n_coding_pieces"] == 2, st2["state"])


def t5_complete_beats_split():
    print(" t5_complete_beats_split")
    locus = [(1000, 6000)]
    genes = [_gene("full", 900, 6100, "+", [(1000, 6000)], name="ITPR1"),
             _gene("extra", 900, 3000, "+", [(1000, 3000)], name="ITPR1")]
    over = R.models_over({"genes": genes}, locus, "+")
    st = R.locus_state(over, locus, "present", R.COMPLETE_FRAC)
    check("T5 a locus with a full model and a redundant one is complete",
          st["state"] == "complete", st["state"])


def t6_naming_model_prefers_the_named_one():
    print(" t6_naming_model_prefers_the_named_one")
    locus = [(1000, 6000)]
    # A big unnamed *coding* model beside a small correctly-named one: the
    # name has to win on the name, not on size, or a locus whose real gene the
    # annotation calls LOC999 reads as correctly named by accident.
    genes = [_gene("bigcod", 900, 5500, "+", [(1000, 5500)], name="LOC999"),
             _gene("real", 5600, 6100, "+", [(5600, 6000)], name="ITPR1")]
    nm = R.naming_model(R.models_over({"genes": genes}, locus, "+"))
    check("T6a the small correctly-named model wins over the big unnamed one",
          nm.get("gene_id") == "real", nm.get("gene_id"))
    genes2 = [_gene("bigcod", 900, 5500, "+", [(1000, 5500)], name="LOC999"),
              _gene("small", 5600, 6100, "+", [(5600, 6000)], name="LOC998")]
    nm2 = R.naming_model(R.models_over({"genes": genes2}, locus, "+"))
    check("T6b with nothing named, the biggest coding model wins",
          nm2.get("gene_id") == "bigcod", nm2.get("gene_id"))
    genes3 = [_gene("bigcod", 900, 5500, "+", [(1000, 5500)], name="LOC999"),
              _gene("lnc", 900, 6100, "+", [], name="ITPR1",
                    product="inositol 1,4,5-trisphosphate receptor type 1",
                    biotype="other")]
    nm3 = R.naming_model(R.models_over({"genes": genes3}, locus, "+"))
    check("T6c a correctly-named non-coding gene beats an unnamed coding one",
          nm3.get("gene_id") == "lnc" and not nm3.get("coding"),
          str(nm3.get("gene_id")))


def t7_name_verdicts():
    print(" t7_name_verdicts")
    cases = [
        ("inositol 1,4,5-trisphosphate receptor type 1", "ITPR1", "ITPR",
         "correct_paralog"),
        ("inositol 1,4,5-trisphosphate receptor type 2", "ITPR1", "ITPR",
         "correct_family_wrong_paralog"),
        ("inositol 1,4,5-trisphosphate receptor", "ITPR1", "ITPR",
         "paralog_unspecified"),
        ("Ryanodine receptor 2", "ITPR1", "ITPR", "wrong_family"),
        ("RyR/IP3R Homology associated domain-containing protein", "ITPR1",
         "ITPR", "family_ambiguous"),
        ("Inositol 1,4,5-trisphosphate/ryanodine receptor domain-containing "
         "protein", "ITPR1", "ITPR", "family_ambiguous"),
        ("LOC101234567", "ITPR1", "ITPR", "placeholder"),
        ("N323_12252", "ITPR1", "ITPR", "placeholder"),
        ("myosin heavy chain", "ITPR1", "ITPR", "family_unnamed"),
        ("", "ITPR1", "ITPR", "absent"),
        ("Ryanodine receptor 2", "", "RYR", "correct_paralog"),
    ]
    for text, exp, fam, want in cases:
        got, _claim = R.name_verdict(text, exp, fam)
        check(f"T7 '{text[:44]}' -> {want}", got == want, got)


def t8_one_rule_both_sides():
    print(" t8_one_rule_both_sides")
    text = "inositol 1,4,5-trisphosphate receptor type 3"
    a = R.name_verdict(text, "ITPR3", "ITPR")
    b = R.name_verdict(text, "ITPR3", "ITPR")
    check("T8a the verdict is a pure function of its inputs", a == b)
    check("T8b the protein audit calls the same function",
          P.R.name_verdict is R.name_verdict)


def t9_panel_reach_is_not_an_error():
    print(" t9_panel_reach_is_not_an_error")
    reach = {"ITPR": {"ITPR1", "ITPR2", "ITPR3"}, "RYR": {"RYR1", "RYR2"}}
    v = P._panel_guard("correct_family_wrong_paralog", "RYR3", "RYR", reach)
    check("T9a a RYR3 claim the panel cannot reach is not a naming error",
          v == "paralog_not_callable", v)
    v2 = P._panel_guard("correct_family_wrong_paralog", "ITPR2", "ITPR", reach)
    check("T9b an ITPR2 claim the panel can reach still is", 
          v2 == "correct_family_wrong_paralog", v2)
    v3 = P._panel_guard("correct_paralog", "RYR3", "RYR", reach)
    check("T9c the guard touches nothing else", v3 == "correct_paralog", v3)


def t10_family_call_needs_a_margin():
    print(" t10_family_call_needs_a_margin")
    labels = {"a": {"family": "ITPR", "paralog": "ITPR1"},
              "b": {"family": "RYR", "paralog": "RYR1"}}
    close = [("a|X|Y", 1000.0), ("b|X|Y", 950.0)]
    clear = [("a|X|Y", 1000.0), ("b|X|Y", 300.0)]
    row = {"accession": "Q", "call": "ITPR", "gene": "", "protein_name": "",
           "group": "Vertebrata", "pfams": ""}
    r1 = P.audit_record(row, close, labels, callable_by_panel={})
    r2 = P.audit_record(row, clear, labels, callable_by_panel={})
    check("T10a two baits within D7's margin give no family call",
          r1["seq_family"] == "no_call", r1["seq_family"])
    check("T10b a clear winner does", r2["seq_family"] == "ITPR",
          r2["seq_family"])
    check("T10c and the margin is recorded either way",
          r1["family_rel_margin"] == 0.05 and r2["family_rel_margin"] == 0.7,
          f"{r1['family_rel_margin']} / {r2['family_rel_margin']}")


def t11_paralog_is_a_vertebrate_question():
    print(" t11_paralog_is_a_vertebrate_question")
    labels = {"a": {"family": "ITPR", "paralog": "ITPR1"},
              "c": {"family": "ITPR", "paralog": "ITPR2"}}
    hits = [("a|X|Y", 1000.0), ("c|X|Y", 400.0)]
    prot = {"accession": "Q", "call": "ITPR", "group": "SAR", "gene": "",
            "protein_name": "inositol 1,4,5-trisphosphate receptor type 2",
            "pfams": ""}
    r = P.audit_record(prot, hits, labels, callable_by_panel={})
    check("T11a a protist record is not asked its paralog",
          r["seq_paralog"] == "not_asked", r["seq_paralog"])
    check("T11b so its 'type 2' name is not a wrong-paralog error",
          r["name_verdict"] != "correct_family_wrong_paralog",
          r["name_verdict"])
    r2 = P.audit_record({**prot, "group": "Vertebrata"}, hits, labels,
                        callable_by_panel={"ITPR": {"ITPR1", "ITPR2"}})
    check("T11c a vertebrate record is", 
          r2["name_verdict"] == "correct_family_wrong_paralog",
          r2["name_verdict"])


def t12_scope_rules_each_fire():
    print(" t12_scope_rules_each_fire")
    base = {"call": "ITPR", "source": "InterPro", "fragment": "0",
            "length": "2700"}
    check("T12a a full-length InterPro family record is in scope",
          P.in_scope(base)[0])
    for k, v, tag in (("call", "unassigned", "P1"), ("source", "S5_genome", "P2"),
                      ("fragment", "1", "P3"), ("length", "800", "P3")):
        ok, why = P.in_scope({**base, k: v})
        check(f"T12 {tag} excludes {k}={v}", (not ok) and why.startswith(tag),
              why)


def t13_calibration_refuses_to_pass_vacuously():
    print(" t13_calibration_refuses_to_pass_vacuously")
    few = [{"state": "complete", "coverage": 1.0, "contig_spans_gene": "1",
            "namer_coding": "1", "symbol_verdict": "correct_paralog",
            "product_verdict": "correct_paralog", "best_model_frac": 0.99}
           for _ in range(10)]
    cal = R.calibrate_complete(few, L.quantile)
    check("T13a a 10-locus calibration is refused, not reported",
          cal["usable"] == 0, cal.get("reason", ""))
    many = [dict(f, best_model_frac=0.9 + 0.001 * i)
            for i, f in enumerate(few * 8)]
    cal2 = R.calibrate_complete(many, L.quantile)
    check("T13b 80 loci are enough", cal2["usable"] == 1)
    check("T13c and the bar is the inherited one, not a new one",
          cal2["bar"] == R.COMPLETE_FRAC == 0.5, str(cal2["bar"]))
    silent = [dict(f, symbol_verdict="absent", product_verdict="absent")
              for f in many]
    check("T13d a locus the annotation is silent about never calibrates it",
          R.calibration_population(silent) == [])


def t14_d6_veto_is_recorded_not_dropped():
    print(" t14_d6_veto_is_recorded_not_dropped")
    row = {"accession": "GCA_X", "cell": "ITPR1", "locus_idx": 0,
           "gene_set": "present", "coverage": 0.95, "state": "unannotated",
           "best_model_frac": 0.0, "union_coding_frac": 0.0,
           "n_coding_models": 0, "n_noncoding_models": 0,
           "noncoding_biotype": "", "locus_cds_bp": 8000,
           "contig": "c1", "start": 1, "end": 9000, "strand": "+",
           "source": "GenBank", "organism": "Test sp.", "vclass": "Aves",
           "symbol_verdict": "absent", "product_verdict": "absent",
           "symbol_claim": "", "product_claim": "", "namer_symbol": "",
           "namer_product": "", "namer_biotype": "",
           "top_overlapper_biotype": "",
           "contig_spans_gene": "1", "integrity_verdict": "elevated_lesions"}
    out = C.from_loci([row])
    check("T14a the vetoed locus is still written", len(out) == 1)
    check("T14b flagged vetoed with its reason",
          out and out[0]["vetoed"] == 1 and "D6" in out[0]["veto_reason"])
    ok = C.from_loci([{**row, "integrity_verdict": "intact"}])
    check("T14c an intact one is not vetoed", ok[0]["vetoed"] == 0)
    check("T14d and is high priority", ok[0]["priority"] == "high",
          ok[0]["priority"])
    low = C.from_loci([{**row, "integrity_verdict": "intact",
                        "contig_spans_gene": "0"}])
    check("T14e below D4's bar it drops to low", low[0]["priority"] == "low",
          low[0]["priority"])


def t15_a_rename_needs_an_absolute_score():
    print(" t15_a_rename_needs_an_absolute_score")
    base = {"accession": "Q", "species": "Test sp.", "tax_class": "",
            "seq_family": "ITPR", "family_rel_margin": 0.5,
            "name_verdict": "wrong_family", "symbol_verdict": "family_unnamed",
            "protein_name": "Ryanodine receptor 1", "gene": "",
            "itpr_bits": 100.0, "ryr_bits": 40.0}
    weak = C.from_proteins([base])
    strong = C.from_proteins([{**base, "itpr_bits": 900.0}])
    check("T15a a 100-bit family call is a low-priority rename",
          weak[0]["priority"] == "low", weak[0]["priority"])
    check("T15b a 900-bit one is high", strong[0]["priority"] == "high")
    none = C.from_proteins([{**base, "name_verdict": "correct_paralog"}])
    check("T15c a correctly-named record proposes nothing", none == [])


def t16_zero_hit_verdicts_each_fire():
    print(" t16_zero_hit_verdicts_each_fire")
    genomes = {"GCA_1": {"organism": "Testus testus", "accession": "GCA_1"}}
    z = [{"Proteome Id": "UP1", "Organism": "Testus testus", "Organism Id": "1",
          "Protein count": "9000"}]
    found = {"GCA_1": {"vclass": "Aves", "cells": {
        "ITPR1": {"is_control": False, "best_coverage": 0.99,
                  "loci": [{"contig_edge": False}]},
        "RYR": {"is_control": True, "best_coverage": 1.0, "loci": []}}}}
    partial = {"GCA_1": {"vclass": "Aves", "cells": {
        "ITPR1": {"is_control": False, "best_coverage": 0.2,
                  "loci": [{"contig_edge": True}]}}}}
    empty = {"GCA_1": {"vclass": "Aves", "cells": {
        "ITPR1": {"is_control": False, "best_coverage": 0.0, "loci": []}}}}
    for sm, want in ((found, "gene_caller_missed_it"),
                     (partial, "assembly_cannot_carry_it"),
                     (empty, "genome_also_empty")):
        got = Z.resolve(z, genomes, sm)[0]["verdict"]
        check(f"T16 {want} fires", got == want, got)
    got = Z.resolve(z, {}, {})[0]["verdict"]
    check("T16 undecidable_no_genome fires when no assembly is in scope",
          got == "undecidable_no_genome", got)


def t17_multiwindow_readers_agree():
    print(" t17_multiwindow_readers_agree")
    root = L.data_root() / "genomes"
    gff = next((p for p in [root / "GCF_000001405.40" / "genomic.gff.gz"]
                if p.exists()), None)
    if gff is None:
        check("T17 skipped — no annotation on the drive", True)
        return
    wins = [("NC_000003.12", 4_600_000, 4_800_000),
            ("NC_000012.12", 26_000_000, 26_100_000)]
    multi = s10_gff.read_annotation_windows(gff, wins)
    for w in wins:
        one = s10_gff.read_annotation_window(gff, *w)
        check(f"T17 {w[0]} one-pass == per-window",
              sorted(one["by_gene"]) == sorted(multi[w]["by_gene"]),
              f"{len(one['by_gene'])} vs {len(multi[w]['by_gene'])}")


def t19_no_span_fallback():
    print(" t19_no_span_fallback")
    # The bug this test exists for: miniprot restarts its IDs per chunk, so a
    # chunked genome's GFF repeats every one and the sweep disambiguates them.
    # Looking up the raw id missed 21 loci, and substituting the locus *span*
    # for the missing CDS blocks put a 2.4 Mb denominator under an 8 kb gene
    # and forced `unannotated` whatever the annotation held.
    st = R.locus_state([], [], "present", R.COMPLETE_FRAC)
    check("T19a an empty CDS block list is unscorable, not unannotated",
          st["state"] == "cds_unavailable", st["state"])
    st2 = R.locus_state([], [(1, 100)], "cds_unavailable", R.COMPLETE_FRAC)
    check("T19b and the flag alone is enough", st2["state"] == "cds_unavailable",
          st2["state"])
    from s5_sweep_lib import parse_miniprot_gff
    sw = L.sweep_dir()
    chunked = None
    for p in sorted(sw.glob("GC*/summary.json")):
        import json
        d = json.loads(p.read_text())
        if int(d.get("miniprot_chunks") or 1) > 1:
            chunked = (p.parent, d)
            break
    if chunked is None:
        check("T19c skipped — no chunked genome in the archive", True)
        return
    d_dir, summ = chunked
    ids = {loc["mp_id"] for c in (summ.get("cells") or {}).values()
           for loc in (c.get("loci") or [])}
    have = {a.mp_id for a in parse_miniprot_gff(d_dir / "miniprot.gff", {})}
    missing = sorted(ids - have)
    check(f"T19c every locus of a chunked genome resolves ({d_dir.name}, "
          f"{summ.get('miniprot_chunks')} chunks)", not missing,
          f"{len(missing)} unresolved: {missing[:3]}")


def t18_miniprot_multi_reader_agrees():
    print(" t18_miniprot_multi_reader_agrees")
    sw = L.sweep_dir()
    cand = next((p for p in sorted(sw.glob("GC*/miniprot.gff"))
                 if p.stat().st_size > 0), None)
    if cand is None:
        check("T18 skipped — no archived miniprot GFF", True)
        return
    ids = []
    for line in cand.read_text().splitlines():
        f = line.split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            for part in f[8].split(";"):
                if part.startswith("ID="):
                    ids.append(part[3:])
        if len(ids) >= 4:
            break
    multi = s10_gff.read_miniprot_models(cand, ids)
    for i in ids:
        one = s10_gff.read_miniprot_model(cand, i)
        check(f"T18 {i} one-pass == per-model", one == multi[i])


def main() -> int:
    print("[s18] negative controls")
    for fn in (t1_same_strand_only, t2_cds_blocks_not_gene_spans,
               t3_every_state_is_reachable, t4_split_needs_two_real_pieces,
               t5_complete_beats_split, t6_naming_model_prefers_the_named_one,
               t7_name_verdicts, t8_one_rule_both_sides,
               t9_panel_reach_is_not_an_error,
               t10_family_call_needs_a_margin,
               t11_paralog_is_a_vertebrate_question,
               t12_scope_rules_each_fire,
               t13_calibration_refuses_to_pass_vacuously,
               t14_d6_veto_is_recorded_not_dropped,
               t15_a_rename_needs_an_absolute_score,
               t16_zero_hit_verdicts_each_fire,
               t17_multiwindow_readers_agree,
               t18_miniprot_multi_reader_agrees,
               t19_no_span_fallback):
        fn()
    if FAILURES:
        print(f"[s18] {len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        return 1
    print("[s18] all negative controls pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
