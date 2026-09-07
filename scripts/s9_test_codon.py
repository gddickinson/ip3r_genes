"""Negative controls for the S9 rules, run on every build.

The pattern `s5_bait_screen.self_test()` established: each case is
constructed so that exactly one rule can reject it, and the rule that is
supposed to reject it must. A codon alignment is the one artefact in this
project where a silent error is invisible downstream — codeml will fit a
model to a frame-shifted alignment and return a perfectly plausible ω — so
the checks here are about *refusal*, not about output being produced.

Twelve cases:

  T1  a CDS one base short of its protein is refused (frame)
  T2  a CDS encoding a different protein is refused (identity)
  T3  an internal stop is masked, not refused, and counted
  T4  every translation/protein disagreement is masked, so the accepted
      CDS always translates to the aligned protein modulo `X`
  T5  the in-house protein->codon map is gap-exact
  T6  `mask_protein` writes `X` at exactly the masked codons
  T7  `strip_gap_codons` removes whole triplets and never shifts frame
  T8  `labelled_newick` refuses a non-monophyletic foreground
  T9  `$1` marks the clade and `#1` marks only its stem
  T10 `unroot` keeps every tip and leaves the root with >= 3 children
  T11 a locus rerun that differs by more than the bound is refused, one
      inside it is accepted, and an exact match wins over a near one
  T12 the S7 cross-check fires when the paralog clades disagree
  T13 a subset's rows follow the order they were asked for
  T14 two *processes* build a subset in the same order (the hash-seed
      check — same-process repetition cannot see this class of bug)

Run:  python scripts/s9_test_codon.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_lib import parse_newick  # noqa: E402
from s9_cds_lib import translate, validate_and_mask  # noqa: E402
from s9_codeml_lib import labelled_newick, unroot  # noqa: E402
from s9_codon_aln import (inhouse_codon_map, mask_protein,  # noqa: E402
                          strip_gap_codons, subset_rows)
from s9_miniprot_cds import _best_alignment  # noqa: E402
from s9_sets import check_against_s7  # noqa: E402

CODONS = {"M": "ATG", "K": "AAA", "L": "CTG", "S": "AGC", "T": "ACC",
          "V": "GTG", "R": "CGT", "W": "TGG", "*": "TAA"}


def encode(prot: str) -> str:
    return "".join(CODONS[a] for a in prot)


def _check(name: str, ok: bool, detail: str = "") -> tuple[str, bool, str]:
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + (f"  — {detail}" if detail else ""))
    return name, ok, detail


def self_test() -> list[tuple[str, bool, str]]:
    out: list[tuple[str, bool, str]] = []
    prot = "MKLSTVRWKLSTVRWMKLSTVRWKLSTVRW" * 4      # 120 aa
    cds = encode(prot)

    # T1 — frame. One base missing shifts everything after it.
    masked, st = validate_and_mask(cds[:-1], prot)
    out.append(_check("T1 short CDS refused (frame)",
                      masked is None and st["note"] == "length_mismatch",
                      f"note={st.get('note')}"))

    # T2 — identity. A CDS of the right length for a different protein.
    other = "W" * len(prot)
    masked, st = validate_and_mask(encode(other), prot)
    out.append(_check("T2 wrong protein refused (identity)",
                      masked is None and st["note"] == "identity_below_99pct",
                      f"note={st.get('note')} mism={st.get('n_mismatch')}"))

    # T3 — an internal stop is masked and counted, not refused.
    stopped = list(prot)
    stopped[10] = "X"
    cds3 = list(encode(prot))
    cds3[30:33] = list("TAA")
    masked, st = validate_and_mask("".join(cds3), "".join(stopped))
    out.append(_check("T3 internal stop masked and counted",
                      masked is not None and masked[30:33] == "NNN"
                      and st["n_stop_masked"] == 1,
                      f"stops={st.get('n_stop_masked')}"))

    # T4 — every disagreement is masked, so the kept CDS translates to the
    # protein modulo X. This is the property pal2nal is entitled to assume.
    cds4 = list(encode(prot))
    cds4[3:6] = list(CODONS["W"])          # residue 1: L -> W, not a stop
    masked, st = validate_and_mask("".join(cds4), prot)
    tr = translate(masked or "")
    consistent = masked is not None and all(
        a == b or a == "X" for a, b in zip(tr, prot))
    out.append(_check("T4 all disagreements masked, CDS matches protein",
                      consistent and masked[3:6] == "NNN"
                      and st["n_mismatch"] == 1,
                      f"mism={st.get('n_mismatch')} masked={st.get('n_masked')}"))

    # T5 — the in-house map is gap-exact.
    mapped = inhouse_codon_map("MK-L", encode("MKL"))
    out.append(_check("T5 in-house codon map is gap-exact",
                      mapped == "ATGAAA---CTG", mapped))

    # T6 — mask_protein marks exactly the masked codons.
    mp = mask_protein("MK-L", "ATGNNNCTG")
    out.append(_check("T6 mask_protein writes X at the masked codon",
                      mp == "MX-L", mp))

    # T7 — gap-codon stripping never shifts frame.
    rows = {"a": "ATG---AAA", "b": "ATG---CTG"}
    got = strip_gap_codons(rows)
    out.append(_check("T7 all-gap codon columns dropped as whole triplets",
                      got == {"a": "ATGAAA", "b": "ATGCTG"}, str(got)))

    # T8/T9 — the branch-label guard. codeml marks a node, so a foreground
    # that is not a clade would silently become a bigger one.
    tree = unroot(parse_newick("(((a:1,b:1)80:1,(c:1,d:1)90:1)70:1,e:1);"))
    try:
        labelled_newick(tree, {"a", "c"}, "clade")
        ok8 = False
    except ValueError:
        ok8 = True
    out.append(_check("T8 non-monophyletic foreground refused", ok8))

    clade = labelled_newick(tree, {"a", "b"}, "clade")
    stem = labelled_newick(tree, {"a", "b"}, "stem")
    out.append(_check("T9 $1 marks the clade, #1 marks one stem",
                      clade.count("$1") == 1 and stem.count("#1") == 1
                      and "(a,b)$1" in clade.replace(":1", "")
                      and "(a,b)#1" in stem.replace(":1", ""),
                      f"{clade} | {stem}"))

    # T10 — PAML wants an unrooted tree.
    rooted = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    un = unroot(rooted)
    out.append(_check("T10 unroot keeps every tip, root has >= 3 children",
                      un.leaf_names() == frozenset("abcd")
                      and len(un.children) >= 3,
                      f"{len(un.children)} children"))

    # T11 — the locus-rerun bound.
    exp = "M" * 300
    def aln(sta: str) -> dict:
        return {"sta": sta, "atn": "A" * 900, "ata": sta}
    far = aln("K" * 10 + "M" * 290)      # 10/300 = 3.3 %, over the 1 % bound
    near = aln("K" * 2 + "M" * 298)      # 2/300 = 0.7 %, inside it
    exact = aln(exp)
    ok11 = (_best_alignment([far], exp) is None
            and _best_alignment([near], exp)[1] == 2
            and _best_alignment([near, exact], exp)[1] == 0
            and _best_alignment([aln("M" * 299)], exp) is None)
    out.append(_check("T11 locus rerun: far refused, near masked, exact wins",
                      ok11))

    # T12 — the S7 cross-check has to be able to fail.
    fake = {"ITPR1": (frozenset({"x"}), ""),
            "ITPR2": (frozenset({"y"}), ""),
            "ITPR3": (frozenset({"z"}), "")}
    out.append(_check("T12 S7 paralog-clade cross-check fires on a mismatch",
                      bool(check_against_s7(fake))))

    # T13/T14 — reproducibility of the per-set files. A set of the same
    # strings iterates in the same order all through one process, so
    # repeating the call proves nothing; the order has to be checked
    # against the *input* order, and across two interpreters with
    # different hash seeds.
    labels = [f"lab{i}" for i in range(12)]
    code_of = {l: f"c{i}" for i, l in enumerate(labels)}
    trimmed = {c: "ATG" * 3 for c in code_of.values()}
    got = list(subset_rows(labels, trimmed, code_of))
    rev = list(subset_rows(labels[::-1], trimmed, code_of))
    out.append(_check("T13 subset rows follow the order they were asked for",
                      got == [code_of[l] for l in labels]
                      and rev == got[::-1], str(got[:3])))
    orders = {_subset_order_in_subprocess(seed) for seed in ("0", "1", "12345")}
    out.append(_check("T14 subset order is identical across hash seeds",
                      len(orders) == 1 and "" not in orders,
                      next(iter(orders))[:40] if orders else "no output"))
    return out


def _subset_order_in_subprocess(seed: str) -> str:
    """Build one subset in a fresh interpreter at a given PYTHONHASHSEED."""
    import os
    import subprocess
    code = (
        "import sys;sys.path.insert(0,%r);sys.path.insert(0,%r)\n"
        "from s9_codon_aln import subset_rows\n"
        "labels=[f'lab{i}' for i in range(12)]\n"
        "code_of={l:f'c{i}' for i,l in enumerate(labels)}\n"
        "trimmed={c:'ATG'*3 for c in code_of.values()}\n"
        "print(','.join(subset_rows(labels,trimmed,code_of)))\n"
    ) % (str(ROOT), str(Path(__file__).resolve().parent))
    env = dict(os.environ, PYTHONHASHSEED=seed)
    p = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, env=env, timeout=120)
    return p.stdout.strip()


def main() -> None:
    print("S9 selection-rule self-test")
    results = self_test()
    failed = [n for n, ok, _ in results if not ok]
    print(f"\n{len(results) - len(failed)}/{len(results)} passed")
    if failed:
        raise SystemExit("FAILED: " + ", ".join(failed))


if __name__ == "__main__":
    main()
