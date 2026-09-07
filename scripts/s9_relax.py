"""S9 step 3b — HyPhy RELAX: is one paralog's selection *relaxed*?

codeml's branch models ask whether a foreground's ω differs from the
background. RELAX asks a different and more interesting question: whether
the whole distribution of ω on the test branches is pulled *towards* ω = 1
(relaxation, k < 1) or *away* from it (intensification, k > 1). For a
family where every ω is far below 1, "ω is higher here" and "selection is
weaker here" are the same statement — and k is the one that says so
directly, with a test attached.

Three runs, one per paralog: that paralog's clade is the test set and the
other two paralogs are the reference. The unlabelled vertebrate tips that
the S7 tree places in no paralog clade are left **unlabelled** rather than
swept into the reference — a branch whose paralog identity is unresolved
is not evidence about either side of the contrast.

RELAX is run locally so the version is pinned in the toolchain manifest,
rather than depending on a Datamonkey job id that cannot be re-run.

    python scripts/s9_relax.py [--workers 4] [--only ITPR3]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_lib import Node, parse_newick  # noqa: E402
from s9_cds_lib import OUT_DIR, PARALOGS, tool_bin  # noqa: E402
from s9_codeml_lib import mrca, unroot  # noqa: E402

HYPHY_DIR = OUT_DIR / "hyphy"
LIVE = ROOT / "results" / "session_live.json"


#: HyPhy writes non-finite branch estimates as the bare tokens `inf`,
#: `-inf` and `nan`, none of which is legal JSON — Python accepts
#: `Infinity`/`NaN` but not those spellings. The partitioned descriptive
#: model estimates a per-branch ω, and a branch with no synonymous change
#: has an infinite one, so this is normal output rather than a broken run.
#: It matters because the failure is silent in the wrong direction: the
#: analysis succeeds, the parse throws, and a caller that returns `{}` on
#: an exception reports "RELAX found nothing" for a run that found
#: something. The repair is counted so it appears in the table.
_BARE_NONFINITE = re.compile(r"(?<=[:\[,])(\s*)(-?)inf\b|(?<=[:\[,])(\s*)nan\b")


def load_hyphy_json(path: Path) -> tuple[dict, int]:
    """(parsed json, number of non-finite literals repaired)."""
    text = path.read_text()

    n = 0

    def sub(m: re.Match) -> str:
        nonlocal n
        n += 1
        if m.group(3) is not None:
            return f"{m.group(3)}NaN"
        return f"{m.group(1)}{m.group(2)}Infinity"

    fixed = _BARE_NONFINITE.sub(sub, text)
    return json.loads(fixed), n


def tip_sets() -> dict[str, str]:
    with open(OUT_DIR / "tip_codes.tsv") as fh:
        return {r["code"]: r["set"]
                for r in csv.DictReader(fh, delimiter="\t")}


def labelled_tree(tree: Node, test: set[str], reference: set[str]) -> str:
    """Newick with HyPhy `{Test}` / `{Reference}` branch sets.

    HyPhy branch sets are per *branch*, so every branch inside a clade
    carries the tag, not only its tips and not only its stem. A tip in
    neither set is left bare, which is how RELAX is told to estimate it
    without counting it as evidence for either side.
    """
    t_node = mrca(tree, test)
    if t_node is None or (t_node.leaf_names() - test):
        raise SystemExit("the test clade is not monophyletic in tree_all.nwk")

    def rec(n: Node, in_test: bool) -> str:
        inside = in_test or n is t_node
        if n.is_leaf:
            tag = ("{Test}" if inside
                   else "{Reference}" if n.name in reference else "")
            return f"{n.name}{tag}"
        inner = ",".join(rec(c, inside) for c in n.children)
        if inside:
            tag = "{Test}"
        else:
            names = n.leaf_names()
            tag = "{Reference}" if names and names <= reference else ""
        return f"({inner}){tag}"

    return rec(tree, False) + ";"


def run_one(para: str, workers: int, timeout_s: int) -> dict:
    HYPHY_DIR.mkdir(parents=True, exist_ok=True)
    out_json = HYPHY_DIR / f"relax_{para}.json"
    if out_json.exists() and out_json.stat().st_size > 0:
        return load_hyphy_json(out_json)[0]
    sets = tip_sets()
    tree = unroot(parse_newick((OUT_DIR / "tree_all.nwk").read_text()))
    test = {c for c, p in sets.items() if p == para}
    reference = {c for c, p in sets.items()
                 if p in PARALOGS and p != para}
    if not test or not reference:
        raise SystemExit(f"no test/reference tips for {para}")
    nwk = HYPHY_DIR / f"tree_relax_{para}.nwk"
    nwk.write_text(labelled_tree(tree, test, reference) + "\n")
    cmd = [tool_bin("hyphy"), f"CPU={workers}", "relax",
           "--alignment", str(OUT_DIR / "codon_trimmed.fasta"),
           "--tree", str(nwk), "--test", "Test", "--reference", "Reference",
           "--code", "Universal", "--output", str(out_json)]
    (HYPHY_DIR / f"relax_{para}.cmd").write_text(" ".join(cmd) + "\n")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=timeout_s)
    (HYPHY_DIR / f"relax_{para}.log").write_text(
        proc.stdout[-60000:] + "\n---stderr---\n" + proc.stderr[-8000:])
    if not out_json.exists():
        raise SystemExit(f"RELAX produced no output for {para}; see "
                         f"{HYPHY_DIR / f'relax_{para}.log'}")
    return load_hyphy_json(out_json)[0]


def summarise(para: str, res: dict, repaired: int = 0) -> dict:
    """One row of `relax_table.tsv`.

    `status` distinguishes a run that finished and reported k from one whose
    output could not be read. Reporting a failed parse as `k = None` beside
    two real answers invites the reader to take it as a negative result.
    """
    tr = res.get("test results", {})
    tested = res.get("tested", {}).get("0", {})
    k = tr.get("relaxation or intensification parameter")
    return {"paralog": para,
            "k": k,
            "p": tr.get("p-value"),
            "LRT": tr.get("LRT"),
            "direction": ("" if k is None else
                          "relaxed" if k < 1 else "intensified"),
            "n_test_branches": sum(1 for v in tested.values() if v == "Test"),
            "n_reference_branches": sum(1 for v in tested.values()
                                        if v == "Reference"),
            "nonfinite_repaired": repaired,
            "status": "ok" if k is not None else "no_test_result"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", default="")
    ap.add_argument("--timeout-h", type=float, default=12.0)
    args = ap.parse_args()

    todo = [p for p in PARALOGS if not args.only or p == args.only]
    rows = []
    for i, para in enumerate(todo, 1):
        LIVE.write_text(json.dumps({
            "task": "S9", "stage": "relax",
            "steps": [{"label": "RELAX", "done": i - 1, "total": len(todo)}],
            "note": para, "ts": time.strftime("%H:%M:%S")}))
        print(f"[{i}/{len(todo)}] RELAX test={para} ...", flush=True)
        t0 = time.time()
        res = run_one(para, args.workers, int(args.timeout_h * 3600))
        _, repaired = load_hyphy_json(HYPHY_DIR / f"relax_{para}.json")
        row = summarise(para, res, repaired)
        row["minutes"] = round((time.time() - t0) / 60, 2)
        rows.append(row)
        print(f"    k={row['k']} p={row['p']} "
              f"({row['n_test_branches']} test / "
              f"{row['n_reference_branches']} reference branches, "
              f"{row['minutes']} min)", flush=True)

    out = OUT_DIR / "relax_table.tsv"
    cols = ["paralog", "k", "p", "LRT", "direction", "n_test_branches",
            "n_reference_branches", "nonfinite_repaired", "status",
            "minutes"]
    existing: list[dict] = []
    if out.exists():
        with open(out) as fh:
            existing = [r for r in csv.DictReader(fh, delimiter="\t")
                        if r["paralog"] not in {r2["paralog"] for r2 in rows}]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in existing + rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
