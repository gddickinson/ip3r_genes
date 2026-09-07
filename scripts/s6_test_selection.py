"""Negative controls for the S6 selection rules, run on every build.

`s5_bait_screen.self_test()` established the pattern in this project: a
screen that only ever sees real data cannot demonstrate that it rejects
anything. The rules in `s6_rep_spec.py` are all positive tests, and a
positive test that never fires on a case it should reject is indistinguish-
able from no test at all.

Six controls, each aimed at one rule, each constructed from a real census
row so the failure is the *one* thing being tested:

  C1  a fragment                       -> the length gate rejects it
  C2  a RyR-length record as an ITPR   -> the length gate rejects it
  C3  an S20 `contaminant_suspect`     -> R5's verdict gate rejects it
  C4  a protist named "type 2"         -> is not grouped as ITPR2
  C5  the same species twice           -> `new_species` takes it once
  C6  a diversity pool of one genus    -> `pick_diverse` spreads first

Plus two properties of the selection as a whole, which are what D24 needs
from the step *before* the alignment:

  P1  determinism   two runs choose the same set, in the same order, with
                    the same audit — a selector that reshuffles makes the
                    aligned file unreproducible however the aligner is run
  P2  tip labels    unique, and legal as Newick/IQ-TREE tip names

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s6_test_selection.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import s6_rep_spec as spec                                   # noqa: E402
import s6_select_reps as sel                                 # noqa: E402
from s6_lib import binomial, load_census, pick_diverse       # noqa: E402

NEWICK_SAFE = re.compile(r"^[A-Za-z0-9_.]+$")
FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f" — {detail}" if detail
                                                     else ""))
    if not ok:
        FAILURES.append(f"{name}: {detail}")


# ---------------------------------------------------------- the controls

def c1_fragment() -> None:
    ok, why = spec.passes_shape("ITPR", 480, allow_short=False)
    check("C1 fragment rejected by the length gate", not ok, why)
    ok2, why2 = spec.passes_shape("ITPR", 480, allow_short=True)
    check("C1b a 480 aa record is rejected even with the exception open",
          not ok2, why2)


def c2_ryr_length_as_itpr() -> None:
    ok, why = spec.passes_shape("ITPR", 5_100, allow_short=True)
    check("C2 RyR-length record rejected as an ITPR", not ok, why)
    ok2, _ = spec.passes_shape("RYR", 5_100, allow_short=False)
    check("C2b the same length is accepted as the RyR outgroup", ok2)


def c3_contaminant() -> None:
    bad = ["contaminant_suspect", "fragment", "no_genome_backing",
           "module_only", "unresolved"]
    rejected = [v for v in bad if not spec.plant_fungal_ok(v)]
    check("C3 every bad S20 verdict is rejected by R5's gate",
          len(rejected) == len(bad), f"{rejected}")
    check("C3b `real_gene` and a blank verdict stay eligible",
          spec.plant_fungal_ok("real_gene") and spec.plant_fungal_ok(""))


def c4_name_transfer() -> None:
    g = spec.group_of("ITPR", "ITPR2", "Amoebozoa", is_vertebrate=False)
    check("C4 a non-vertebrate 'type 2' is not grouped as ITPR2",
          g == "protist", f"grouped as {g!r}")
    g2 = spec.group_of("ITPR", "ITPR2", "Vertebrata", is_vertebrate=True)
    check("C4b the same label inside the vertebrates is kept",
          g2 == "ITPR2", f"grouped as {g2!r}")
    g3 = spec.group_of("ITPR", "", "Vertebrata", is_vertebrate=True)
    check("C4c an unlabelled vertebrate becomes vertebrate_basal",
          g3 == "vertebrate_basal", f"grouped as {g3!r}")


def c5_same_species_twice() -> None:
    a = "Prymnesium parvum (Toxic golden alga)"
    b = "Prymnesium parvum"
    check("C5 the two spellings of one species collapse",
          binomial(a) == binomial(b) == "Prymnesium parvum")


def c6_diversity_spread() -> None:
    rows = [{"p": "Oomycota", "g": "Achlya", "s": 9},
            {"p": "Oomycota", "g": "Aphanomyces", "s": 8},
            {"p": "Oomycota", "g": "Saprolegnia", "s": 7},
            {"p": "Ciliophora", "g": "Stylonychia", "s": 6},
            {"p": "Perkinsozoa", "g": "Perkinsus", "s": 5}]
    got = pick_diverse(rows, 3, lambda r: (r["p"], r["g"]),
                       lambda r: r["s"])
    phyla = {r["p"] for r in got}
    check("C6 a score-ranked one-clade pool is spread before it is deepened",
          len(phyla) == 3, f"phyla chosen: {sorted(phyla)}")


def c7_single_valued_level_does_not_throttle() -> None:
    """A key level with one value must not eat the deeper level's spread.

    The C6 pool has three phyla, so the phylum level does real work there
    and the bug hid. Here every row is one phylum and the spread has to
    come from the genus level. Before the fix the phylum prefix admitted
    one row per wave, and by wave 3 the genus quota was 3 — loose enough
    to take a second Chlamydomonas while six other genera waited.
    """
    rows = [{"p": "Chlorophyta", "g": "Chlamydomonas", "s": 9},
            {"p": "Chlorophyta", "g": "Cymbomonas", "s": 8},
            {"p": "Chlorophyta", "g": "Chlamydomonas", "s": 7},
            {"p": "Chlorophyta", "g": "Volvox", "s": 6},
            {"p": "Chlorophyta", "g": "Gonium", "s": 5}]
    got = pick_diverse(rows, 3, lambda r: (r["p"], r["g"]),
                       lambda r: r["s"])
    genera = [r["g"] for r in got]
    check("C7 one phylum, many genera: the spread falls to the genus level",
          len(set(genera)) == 3, f"genera chosen: {genera}")


def c8_no_group_borrows_another_groups_ruler() -> None:
    """Every group with a clean record of its own is measured by it.

    The failure was silent and systematic: a group under the record floor
    took the *global* median, which is a metazoan number, and in the plant
    grade that decided a tip. A group may legitimately have nothing of its
    own — that is what the `global` tier is for — but it must then be a
    group with no clean record at all, not merely a sparse one.
    """
    itpr = [r for r in load_census() if r["call"] == "ITPR"]
    _, audit = sel.length_targets(itpr)
    borrowed = [a for a in audit
                if a["tier"] == "global" and int(a["n_arch4plus"]) > 0]
    check("C8 no group with clean records of its own borrows the global "
          "ruler", not borrowed,
          "; ".join(f"{a['group']} ({a['n_arch4plus']} records)"
                    for a in borrowed) or f"{len(audit)} groups, all own")
    thin = [a for a in audit if a["tier"] != "global"
            and int(a["n_records"]) < 3]
    check("C8b a group's own ruler rests on at least 3 of its records",
          not thin,
          "; ".join(f"{a['group']} n={a['n_records']}" for a in thin)
          or "all >= 3")



# --------------------------------------------------------- the properties

def p1_determinism(rows) -> None:
    a = sel.select(rows)
    b = sel.select(rows)
    ka, kb = list(a.chosen), list(b.chosen)
    check("P1 two selections choose the same accessions in the same order",
          ka == kb, f"{len(set(ka) ^ set(kb))} differ")
    fields = ("rule", "cell", "group", "decided_by", "runner_up")
    diff = [k for k in ka
            if any(a.chosen[k][f] != b.chosen[k][f] for f in fields)]
    check("P1b the audit is identical too", not diff, f"{diff[:5]}")
    check("P1c the unfilled-slot list is identical",
          a.unfilled == b.unfilled)
    return a


def p2_labels(picker) -> None:
    labels = [c["label"] for c in picker.chosen.values()]
    bad = [l for l in labels if not NEWICK_SAFE.match(l)]
    check("P2 every tip label is Newick-legal", not bad, f"{bad[:3]}")
    check("P2b every tip label is unique",
          len(set(labels)) == len(labels),
          f"{len(labels) - len(set(labels))} collisions")
    long = [l for l in labels if len(l) > 90]
    check("P2c no tip label is absurdly long", not long, f"{long[:2]}")


def p3_rules_held(picker) -> None:
    """The rules held on the *real* census, not only on constructed cases."""
    chosen = list(picker.chosen.values())
    para_outside = [c for c in chosen
                    if c["census_group"] != "Vertebrata"
                    and c["group"] in spec.PARALOGS]
    check("P3 no non-vertebrate tip is grouped as a vertebrate paralog",
          not para_outside, f"{[c['label'] for c in para_outside][:3]}")
    if "RYR" in {c["group"] for c in chosen}:
        short_ryr = [c for c in chosen if c["group"] == "RYR"
                     and c["length"] < spec.RYR_MIN_AA]
        check("P3b every outgroup tip clears the RyR floor", not short_ryr,
              f"{[c['label'] for c in short_ryr][:3]}")
    grades = {c["group"] for c in chosen}
    check("P3c the tree has an outgroup", "RYR" in grades)
    check("P3d the tree has the pre-2R vertebrate grade",
          "vertebrate_basal" in grades)


def main() -> int:
    print("S6 selection self-test")
    c1_fragment()
    c2_ryr_length_as_itpr()
    c3_contaminant()
    c4_name_transfer()
    c5_same_species_twice()
    c6_diversity_spread()
    c7_single_valued_level_does_not_throttle()
    c8_no_group_borrows_another_groups_ruler()
    rows = load_census()
    picker = p1_determinism(rows)
    p2_labels(picker)
    p3_rules_held(picker)
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILURES")
        return 1
    print("\nall controls passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
