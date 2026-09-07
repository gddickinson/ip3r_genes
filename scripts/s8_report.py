"""S8 — renders `results/synteny/report.md` purely from the committed
tables (D13).

Scope and instrument live here; the results live in
`s8_report_results.py` (the `s3_report.py` / `s3_report_d10.py` split, so
both halves stay under 500 lines and take this module's loader and
formatter and so cannot read the tables differently).

The report titles itself from the scale it is rendering, for the same
reason `s5_report.py` and `s23_report.py` do.

Run:  python3 scripts/s8_report.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

OUT = Path(__file__).resolve().parents[1] / "results" / "synteny"


def load(name: str) -> list[dict]:
    p = OUT / name
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    p = OUT / name
    return json.loads(p.read_text()) if p.exists() else {}


def num(v, nd=0) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return f"{f:,.{nd}f}"


def table(header: list[str], rows: list[list]) -> list[str]:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out + [""]


def stat_rows(stats: list[dict], window: str, key: str, stratum: str = "all"):
    return {r["pair_class"]: r for r in stats
            if r["window"] == window and r["key"] == key
            and r["stratum"] == stratum}


# ------------------------------------------------------------------ scope

def section_scope(loci: list[dict], meta: dict) -> list[str]:
    n_gen = len({r["accession"] for r in loci})
    with_tab = meta["n_genomes_with_gene_table"]
    by_cell = Counter(r["cell"] for r in loci)
    multi = Counter()
    for r in loci:
        if int(r["locus_idx"]) > 0:
            multi[r["cell"]] += 1
    L = [f"# S8 — Synteny of the ITPR loci across {num(n_gen)} genomes\n",
         f"*Generated {meta['generated']}; "
         f"{meta['flank_n']} flanking genes per side, "
         f"`{meta['primary_window']}` as the primary window and "
         f"`{meta['call_window']}` for the caller, "
         f"{meta['n_control_replicates']} control windows per genome "
         f"(seed {meta['control_seed']}), runtime "
         f"{num(meta['runtime_s'], 1)} s.*\n",
         "**What this task asks.** S5 assigned a paralog to each locus from "
         "sequence: which bait won the alignment. Synteny asks the same "
         "question with evidence the alignment never sees — what the "
         "*neighbouring* genes are called — so that an orthology claim does "
         "not rest twice on the same measurement.\n",
         "## 1. What is in the analysis\n",
         f"Loci come from the per-genome `summary.json` files rather than "
         f"from the ledger, because the ledger carries one row per genome × "
         f"cell and therefore only the best locus. A teleost ITPR1 cell "
         f"holds *itpr1a* and *itpr1b*; the sea lamprey's holds three. "
         f"Reading the best one would compare *itpr1a* in one species "
         f"against *itpr1b* in the next and report the mismatch as a "
         f"synteny result.\n"]
    L += table(["cell", "loci", "of which extra copies in the same cell"],
               [[c, num(by_cell[c]), num(multi[c])]
                for c in ("ITPR1", "ITPR2", "ITPR3", "RYR")])
    L += [f"- Loci with coordinates across {num(n_gen)} swept genomes: "
          f"**{num(len(loci))}**",
          f"- Genomes carrying an annotation gene table: **{num(with_tab)}** "
          f"(the rest are unannotated assemblies; "
          f"**{num(meta['n_loci_without_gene_table'])}** loci have no flanks "
          "and are absent from every statistic below)",
          f"- Genomes contributing control windows: "
          f"**{num(meta['n_control_genomes'])}** — a genome whose longest "
          f"contig carries fewer than {2 * meta['flank_n'] + 1} coding genes "
          "cannot host a fair control window and is excluded from the null "
          "rather than allowed to degrade it",
          f"- Cells excluded by name: **{num(meta['n_excluded_cells'])}** "
          "(the mechanism exists so that an exclusion cannot be made "
          "silently; S7's relabelling rule fired on no tip, so none is due)",
          ""]
    return L


# ------------------------------------------------------------- instrument

def section_instrument(sets_: list[dict], meta: dict) -> list[str]:
    prim, call = meta["primary_window"], meta["call_window"]
    by_w = {}
    for r in sets_:
        w = r["window"]
        by_w.setdefault(w, []).append(int(r["n_informative"]))
    L = ["## 2. The instrument\n",
         "### 2.1 Two windows, because annotation naming density is not "
         "constant\n",
         "Matching flanks across species means matching gene *symbols*, and "
         "a symbol is only usable if the annotation gave the gene one. "
         "Across this genome set that varies about fourfold — 97 % of human "
         "coding genes carry an informative symbol against 27 % of the sea "
         "lamprey's — so a fixed ±10-gene window hands a well-named genome "
         "20 usable keys and a poorly-named one 5. A Jaccard difference "
         "would then be an annotation difference.\n",
         "So both windows are committed and every table carries a `window` "
         "column:\n"]
    L += table(["window", "rule", "mean informative keys per locus"],
               [[f"`{w}`",
                 ("the brief's rule — the 10 nearest coding genes each side"
                  if w == "fixed10" else
                  f"the 10 nearest *informative* symbols each side, scanning "
                  f"at most {__import__('s8_flank_lib').MAX_SCAN} genes"),
                 f"{sum(v) / len(v):.1f}"]
                for w, v in sorted(by_w.items())])
    L += [f"`{prim}` is primary for the pair statistics because it is the "
          f"brief's rule and it makes no assumption about what a symbol is; "
          f"`{call}` is what the consensus caller uses, because a caller "
          "whose evidence depends on the assembly's naming conventions "
          "cannot be calibrated across them.\n",
          "### 2.2 Three key vocabularies\n",
          "A symbol is normalised before it is compared, and the three "
          "levels answer different questions:\n",
          "| key | rule | what it can see |",
          "|---|---|---|",
          "| `strict` | the symbol, uppercased | the same gene, same name |",
          "| `relaxed` | lowercase teleost/amphibian duplicate suffixes "
          "stripped (`gnasa` → `GNAS`, trailing `.2` dropped); uppercase "
          "symbols untouched, so `GNB1` and `BAK1` survive | the same gene "
          "across clades that name it differently |",
          "| `root` | the relaxed key with its trailing digit run removed "
          "(`BHLHE40`, `BHLHE41` → `BHLHE`), floored at 3 characters so "
          "`TP53` and `C3` keep their key | the same gene **family** — which "
          "is the only level at which a 2R ohnolog pair is visible, because "
          "the two copies almost never carry the same symbol |",
          "",
          "The root key over-merges (every zinc finger reaches `ZNF`). That "
          "is tolerable only because it is scored against a null built with "
          "the *same* rule on random windows, so over-merging inflates the "
          "signal and its background together.\n",
          "### 2.3 The control: matched random neighbourhoods\n",
          "A Jaccard of 0.21 means nothing on its own. Two neighbourhoods "
          "in two well-annotated mammals share vocabulary for reasons that "
          "have nothing to do with orthology, and two in a hagfish and a "
          "lamprey share almost none whatever their history.\n",
          f"So every real pair *(locus in genome A, locus in genome B)* is "
          f"scored against matched control pairs *(random coding gene in A, "
          f"random coding gene in B)* — {meta['n_control_replicates']} "
          f"replicates per genome, drawn with the same window rule and the "
          f"same key rule, seeded per accession so the draw is reproducible. "
          "The control therefore holds constant the two genomes, their "
          "annotation depth, their naming conventions, the window and the "
          "normalisation. What is left is the locus.\n",
          "The paired comparison is reported as a **sign test** rather than "
          "a difference of means: Jaccard is bounded, zero-inflated and "
          "nowhere near normal, and pairs where both the locus and its "
          "control score zero are counted as ties and dropped rather than "
          "scored as failures, because counting them against the locus "
          "would make an unannotated genome look like a negative result.\n"]
    return L


def section_selftest(meta: dict) -> list[str]:
    names = meta.get("self_tests", [])
    L = ["### 2.4 The negative controls, run on every build\n",
         "The rules that decide what a flank is, when two symbols are the "
         "same gene, and when a consensus call is evidence are all "
         "constructed-case tested before anything is measured with them "
         "(`scripts/s8_test_flanks.py`, the pattern of "
         "`s5_bait_screen.self_test()` and `s7_test_tree.py`). A failure "
         "aborts the run.\n"]
    detail = {
        "T1": "the placeholder vocabularies — `LOC116952798`, "
              "`FN964_004414`, `ENSG…`, `si:ch211-…` — are rejected and "
              "real symbols kept",
        "T2": "`gnasa` reaches `GNAS`; `GNB1`, `BAK1`, `ctsa`, `vapb` are "
              "left alone",
        "T3": "`BHLHE40`/`BHLHE41` and `GRM7`/`GRM4` merge; `TP53` and `C3` "
              "are not reduced below the character floor",
        "T4": "the locus itself and its overlappers stay out of its own "
              "flank set, non-coding biotypes are filtered, and the two "
              "windows behave differently on the same table",
        "T5": "Jaccard is symmetric and **empty-vs-empty is 0, not 1** — "
              "otherwise two unannotated genomes would score as a perfect "
              "synteny match",
        "T6": "the pair subset takes one locus per species, the "
              "highest-covered, so a teleost *a*/*b* pair cannot be "
              "compared across species",
        "T7": "the control draw is reproducible and refuses a contig too "
              "short to hold a window",
        "T8": "the caller's two halves — a genuine neighbourhood is called, "
              "an unrelated one and a key-poor one are refused, and a locus "
              "cannot be scored against a consensus it voted into",
        "T9": "the threshold rule refuses a setting that calls everything, "
              "including 90 % of random windows",
        "T10": "a call is only `supported` above the highest score any "
               "random window reached",
        "T11": "a ranked table's row order does not depend on which order "
               "the input arrived in — Python's set iteration is hash-seeded "
               "per process, and a rank without a final tiebreak differs "
               "between two runs of identical data",
    }
    L += table(["test", "what must fail"],
               [[f"**{n}**", detail.get(n, "")] for n in names])
    L += ["`synteny_stats.json` also records the SHA-256 of every committed "
          "table, so a rerun that drifts is visible in the data rather than "
          "only in a diff.\n"]
    return L


def main() -> int:
    meta = load_json("synteny_stats.json")
    if not meta:
        raise SystemExit("run scripts/s8_run_synteny.py first")
    loci = load("loci.tsv")
    sets_ = load("locus_sets.tsv")
    L = section_scope(loci, meta)
    L += section_instrument(sets_, meta)
    L += section_selftest(meta)
    import s8_report_results as res
    L += res.render(load=load, load_json=load_json, table=table, num=num,
                    stat_rows=stat_rows, meta=meta, loci=loci)
    out = OUT / "report.md"
    out.write_text("\n".join(L) + "\n")
    print(f"wrote {out} ({len(L)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
