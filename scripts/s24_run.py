"""The S24 driver — ordered stages `guards -> figures -> audit -> tables ->
report`, with `--only` / `--from` / `--list`.

`guards` runs first and is not optional: `s24_test_supp.self_test()`, then the
exact column walk and the residue join over the real data. Nothing is written
if either fails, because every panel S24 draws is a claim about a join, and a
figure drawn from coordinates that did not check out is worse than no figure.

A failed stage does not stop the ones after it that do not depend on it, and
the report marks an unfinished section *not run yet* — the S9/S10/S12/S15
pattern. Everything read is committed, so the task runs offline and a rerun
is deterministic.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import figstyle as F                                           # noqa: E402
import s24_lib as L                                            # noqa: E402
import s24_audit as A                                          # noqa: E402
import s24_tables as T                                         # noqa: E402
import s24_test_supp as TEST                                   # noqa: E402

STAGES = ["guards", "figures", "audit", "tables", "report"]
STATS = L.OUT_DIR / "supplementary_stats.json"


def load_stats() -> dict:
    import json
    if STATS.exists():
        with open(STATS) as fh:
            return json.load(fh)
    return {}


def stage_guards(stats: dict) -> None:
    if TEST.self_test() != 0:
        raise SystemExit("S24 negative controls failed — nothing written")
    stats["self_test"] = {
        "checks": len(TEST.CHECKS), "passed": len(TEST.CHECKS),
        "mutations_tested": len(TEST.MUTATIONS),
        "mutations_caught": len(TEST.MUTATIONS),
        "mutations": [{"breakage": b, "caught_by": c}
                      for b, c in TEST.MUTATIONS],
    }
    aln = L.read_fasta(L.MSA_DIR / "aln.fasta")
    trimmed = L.read_fasta(L.MSA_DIR / "trimmed.fasta")
    colmap = L.read_tsv(L.MSA_DIR / "column_map.tsv")
    g = L.verify_column_map(aln, trimmed, colmap)
    g.pop("kept_indices", None)
    index = {p: L.residue_index(p) for p in L.PARALOGS}
    g.update(L.verify_residues(
        L.read_tsv(L.CONSTRAINT_DIR / "variants.tsv"),
        L.read_tsv(L.CONSTRAINT_DIR / "paralog_variant_positions.tsv"),
        index))
    stats["guards"] = g
    print(f"[s24] column map verified by an exact walk: "
          f"{g['columns_checked']:,} columns x {g['sequences']} sequences")
    print(f"[s24] residues verified: {g['variant_residues_checked']:,} "
          f"variants, {g['aligned_partners_checked']:,} aligned partners")


def stage_figures(stats: dict) -> None:
    import s24_figs_alignment as FA
    import s24_figs_inputs as FI
    import s24_figs_structure as FS
    F.use()
    L.FIG_DIR.mkdir(parents=True, exist_ok=True)
    for name, fn in (list(FA.FIGURES.items()) + list(FI.FIGURES.items())
                     + list(FS.FIGURES.items())):
        print(f"[s24] figure {name}")
        fn(L.FIG_DIR, stats)


def stage_audit(stats: dict) -> None:
    out = A.write(L.OUT_DIR)
    stats["audit"] = out
    for p in out["problems"]:
        print(f"[s24] AUDIT PROBLEM: {p}")
    print(f"[s24] audited {out['figures_audited']} manuscript figures, "
          f"{out['findings']} findings recorded")


def stage_tables(stats: dict) -> None:
    out = T.write(stats)
    stats["tables"] = out
    stats["headline"] = T.headline(stats)
    if out["missing"]:
        print(f"[s24] MISSING figure files: {out['missing']}")
    print(f"[s24] supp_figure_stats.tsv: {out['rows']} rows")


def stage_report(stats: dict) -> None:
    import s24_report
    s24_report.render(stats)
    print(f"[s24] report -> {L.OUT_DIR / 'report.md'}")


RUNNERS = {"guards": stage_guards, "figures": stage_figures,
           "audit": stage_audit, "tables": stage_tables,
           "report": stage_report}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", choices=STAGES)
    ap.add_argument("--from", dest="start", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for s in STAGES:
            print(s)
        return 0
    todo = (args.only if args.only else
            STAGES[STAGES.index(args.start):] if args.start else STAGES)

    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    stats = load_stats()
    failed = []
    for name in STAGES:
        if name not in todo:
            continue
        try:
            RUNNERS[name](stats)
        except SystemExit:
            raise
        except Exception:                          # noqa: BLE001
            failed.append(name)
            print(f"[s24] stage {name} FAILED")
            traceback.print_exc()
    L.dump_stats(stats, STATS)
    if failed:
        print(f"[s24] failed stages: {', '.join(failed)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
