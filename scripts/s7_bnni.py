"""S7 — does the `--bnni` re-run change any claim this task makes?

UFBoot is known to be optimistic when the model is violated, and a
134-tip alignment spanning four kingdoms violates any single
substitution model by construction. `--bnni` is UFBoot's own guard
against that: each bootstrap tree gets an extra round of NNI
optimisation, which strips the support that came from the model rather
than from the data.

It is a **check, not a second answer**. The tree the paper reports is
still the main search's; what this module asks is whether the claims
that tree carries survive the guard. So it re-asks the *same* clade
questions — read out of `claim_members.tsv`, never listed here, so
the two runs cannot be scored on different claims — of the `--bnni`
tree, and
writes one row per claim with both runs' verdicts side by side.

Three things can happen to a claim, and the table says which:

  `held`      the clade is present in both trees and clears SH-aLRT ≥ 80
              and UFBoot ≥ 95 in both.
  `weakened`  present in both, but the support the main run reported does
              not survive the guard. This is the outcome the flag exists
              to expose, and a claim that lands here is a claim the paper
              must state with the `--bnni` number, not the main one.
  `lost`      the clade is not in the `--bnni` tree at all — a topology
              change, not a support change.

A claim that is *not* well supported in the main run and becomes so
under `--bnni` is `strengthened`; it is reported for completeness and
leaned on by nothing, because the reported tree is the main one.

It also serves the other question a set of searches raises: when a
*constrained* search finds a tree better than the unconstrained one —
which happened here, the ITPR2+ITPR3 constraint reaching a likelihood
4.8 units above the main search's — the reported tree is not the global
optimum, and what has to be said is *which claims that changes*. Pass
`--tree/--label/--out` to score the claim set against any alternative
tree. A tree with no support labels (a constrained search runs without
`-B`/`-alrt`) is compared on topology alone and its verdicts read
`present` / `absent`, because grading it `weakened` would report a loss
of support that was never measured.

Run:  python3 scripts/s7_bnni.py
      python3 scripts/s7_bnni.py --tree <nwk> --label <name> --out <tsv>
In:   results/phylogeny/itpr_ml_bnni.treefile
      results/phylogeny/claim_nodes.tsv + claim_members.tsv
      results/phylogeny/rooted.nwk
Out:  results/phylogeny/bnni_comparison.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from s6_lib import read_tsv, write_tsv                        # noqa: E402
from s7_lib import (                                          # noqa: E402
    PHYLO_DIR, find_clade, group_labels, load_groups, parse_newick, reroot,
    support_of,
)

BNNI = PHYLO_DIR / "itpr_ml_bnni.treefile"
MAIN = PHYLO_DIR / "rooted.nwk"
CLAIMS = PHYLO_DIR / "claim_nodes.tsv"
MEMBERS = PHYLO_DIR / "claim_members.tsv"
OUT = PHYLO_DIR / "bnni_comparison.tsv"

MIN_ALRT, MIN_UFBOOT = 80.0, 95.0

COLS = ["claim", "n_tips", "main_is_clade", "main_alrt", "main_ufboot",
        "alt_is_clade", "alt_alrt", "alt_ufboot", "verdict"]


def strong(alrt, ufboot) -> bool:
    try:
        return float(alrt) >= MIN_ALRT and float(ufboot) >= MIN_UFBOOT
    except (TypeError, ValueError):
        return False


def claim_members() -> dict[str, set[str]]:
    """Each claim's tip set, as `s7_analyse.py` recorded it.

    Read, never re-derived: recovering a claim's membership here from
    its size would let this module and `s7_analyse.py` disagree about
    what "ITPR2 + ITPR3" means, and the disagreement would surface as a
    `lost` verdict that is really a bookkeeping difference.
    """
    out: dict[str, set[str]] = {}
    for r in read_tsv(MEMBERS):
        out.setdefault(r["claim"], set()).add(r["label"])
    return out


def compare(alt_path: Path, label: str, out: Path) -> int:
    """Score every claim of the main tree against an alternative tree."""
    if not alt_path.exists():
        print(f"no tree at {alt_path}")
        return 1
    if not (MAIN.exists() and CLAIMS.exists() and MEMBERS.exists()):
        raise SystemExit("run `python3 scripts/s7_analyse.py` first")

    claims = read_tsv(CLAIMS)
    main_tree = parse_newick(MAIN.read_text())
    groups = load_groups()
    alt = parse_newick(alt_path.read_text())
    og = group_labels(groups, "RYR") & alt.leaf_names()
    alt = reroot(alt, og)

    if main_tree.leaf_names() != alt.leaf_names():
        raise SystemExit("the two trees do not carry the same tips")

    # A constrained search runs without -B/-alrt, so its tree carries no
    # support labels at all. Grading such a tree `weakened` would report
    # a loss of support that was never measured, so a tree with no
    # support is compared on **topology only** and the verdicts say so.
    scored = any(support_of(n.name)[0] is not None
                 for n in alt.walk() if not n.is_leaf)
    members = claim_members()
    rows: list[dict] = []
    tally: dict[str, int] = {}
    for c in claims:
        mem = members.get(c["claim"])
        # A claim the main tree does not hold as a clade cannot be
        # `lost` or `absent` from the alternative — it was never there.
        # `claim_members.tsv` carries a tip set for every claim, clade
        # or not, so without this the three non-clade rows scored as
        # three disagreements and the summary read "3 absent" when the
        # two trees in fact agree on everything the report claims.
        if mem is not None and c["is_clade"] != "yes":
            rows.append({**{k: "" for k in COLS}, "claim": c["claim"],
                         "n_tips": c["n_tips"], "main_is_clade": "no",
                         "verdict": "not a clade in the main tree"})
            continue
        if mem is None:
            rows.append({**{k: "" for k in COLS}, "claim": c["claim"],
                         "n_tips": c["n_tips"],
                         "main_is_clade": c["is_clade"],
                         "main_alrt": c["alrt"], "main_ufboot": c["ufboot"],
                         "verdict": "not a clade in the main tree"})
            continue
        ok, sup = find_clade(alt, mem)
        ba, bu = support_of(sup) if ok else (None, None)
        if not scored:
            verdict = "present" if ok else "absent"
        else:
            m_strong = strong(c["alrt"], c["ufboot"])
            b_strong = ok and strong(ba, bu)
            if not ok:
                verdict = "lost"
            elif m_strong and b_strong:
                verdict = "held"
            elif m_strong:
                verdict = "weakened"
            elif b_strong:
                verdict = "strengthened"
            else:
                verdict = "unsupported in both"
        tally[verdict] = tally.get(verdict, 0) + 1
        rows.append({
            "claim": c["claim"], "n_tips": c["n_tips"],
            "main_is_clade": c["is_clade"], "main_alrt": c["alrt"],
            "main_ufboot": c["ufboot"],
            "alt_is_clade": "yes" if ok else "no",
            "alt_alrt": "" if ba is None else f"{ba:g}",
            "alt_ufboot": "" if bu is None else f"{bu:g}",
            "verdict": verdict})

    write_tsv(out, COLS, rows)
    print(f"{label}: "
          + ", ".join(f"{n} {v}" for v, n in sorted(tally.items()))
          + f" -> {out}")
    return 0


def main() -> int:
    args = sys.argv[1:]
    tree, label, out = BNNI, "--bnni", OUT
    while args:
        a = args.pop(0)
        if a == "--tree":
            tree = Path(args.pop(0))
        elif a == "--label":
            label = args.pop(0)
        elif a == "--out":
            out = PHYLO_DIR / args.pop(0)
        else:
            print(__doc__)
            return 2
    if tree is BNNI and not BNNI.exists():
        print(f"no --bnni tree yet ({BNNI.name}) — "
              "run `python3 scripts/s7_run.py bnni` first")
        return 1
    return compare(tree, label, out)


if __name__ == "__main__":
    raise SystemExit(main())
