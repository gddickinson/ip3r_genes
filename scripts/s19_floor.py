"""S19 — calibrating D4's contiguity bar, and what a floor costs.

Split out of `s19_contiguity` to keep both inside the 500-line budget. The
first module measures the false-negative rate; this one asks what to *do*
about it — the residual error at every contiguity floor, which clades each
floor removes, and whether the misses that survive the strictest floor are
about the gene at all.

D4's bar was chosen a priori as the median measured ITPR genomic span, from
gene geometry with no error rate in its derivation. `floor_scan` is the first
time it has been scored against one.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402

#: Contig-N50 floors scanned. The grid straddles D4's own bar, which is
#: included exactly so the calibration contains the value it may replace.
FLOOR_GRID = [0, 10_000, 20_000, 50_000, 100_000, 142_212, 200_000,
              500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000,
              20_000_000, 50_000_000]

#: Target residual false-negative rates. 0.05 is conventional; 0.01 is
#: reported beside it because a survey making 1,236 cell calls turns 5 %
#: into 62 wrong ones.
FN_TARGETS = (0.05, 0.01)

#: Below this share of its paralog's consensus neighbourhood present in the
#: assembly, the *region* is what is missing and the cell says nothing about
#: the gene. S8's own flank consensus supplies the symbols.
NEIGHBOURHOOD_PRESENT = 0.50
FLANK_MIN_FRACTION = 0.30
FLANK_MAX_SYMBOLS = 12
#: S8 wrote its consensus under two key vocabularies and two window rules;
#: the pair used here is the one S10's own flank check reads, so a locus
#: cannot be judged against a different neighbourhood in two tasks.
FLANK_KEY_KIND, FLANK_WINDOW = "relaxed", "fixed10"

def floor_scan(cells: list[dict], log=S.log) -> tuple[list[list], dict]:
    """The residual false-negative rate at every contiguity floor.

    This is the calibration D4 never had. Both control series are scanned, and
    the genomes each floor costs are reported beside the rate it buys, because
    a floor that reaches 1 % by keeping 40 assemblies has not made the scope
    more reliable — it has made it smaller.
    """
    series = {n: [c for c in cells if c["control"] == n]
              for n in ("itpr_present", "ryr_sister")}
    rows = []
    for floor in FLOOR_GRID:
        for name, ctrl in series.items():
            kept = [c for c in ctrl if c["contig_n50"] >= floor]
            if not kept:
                continue
            fn = sum(1 for c in kept if c["false_negative"])
            lo, hi = S.wilson(fn, len(kept))
            rows.append([floor, name, len({c["accession"] for c in kept}),
                         round(len({c["accession"] for c in kept}) / 309, 4),
                         len(kept), fn, round(fn / len(kept), 4),
                         round(lo, 4), round(hi, 4),
                         int(floor == S.contiguity_bar())])
    S.write_tsv(S.out_dir() / "absence_floor.tsv",
                ["contig_n50_floor", "series", "n_genomes_kept",
                 "frac_genomes_kept", "n_cells", "n_false_negative",
                 "fn_rate", "wilson_lo", "wilson_hi", "is_d4_bar"], rows)

    rec: dict = {}
    for name in series:
        mine = [r for r in rows if r[1] == name]
        rec[name] = {}
        for target in FN_TARGETS:
            pick = next((r for r in mine if r[8] <= target), None)
            rec[name][str(target)] = {
                "floor": pick[0] if pick else None,
                "genomes_retained": pick[2] if pick else 0,
                "genomes_retained_frac": pick[3] if pick else 0,
                "achieved_fn_rate": pick[6] if pick else None,
                "achieved_upper_ci": pick[8] if pick else None,
            }
        at_bar = next((r for r in mine if r[9] == 1), None)
        rec[name]["at_d4_bar"] = {
            "floor": S.contiguity_bar(),
            "genomes_retained": at_bar[2] if at_bar else 0,
            "fn_rate": at_bar[6] if at_bar else None,
            "upper_ci": at_bar[8] if at_bar else None,
        }
        log(f"floor_scan {name}: at D4's bar FN="
            f"{rec[name]['at_d4_bar']['fn_rate']}, "
            f"for FN<=0.05 floor={rec[name]['0.05']['floor']}")
    return rows, rec


# ------------------------------------------------ 5. what the floor costs

def clade_composition(cells: list[dict], log=S.log) -> list[list]:
    """Which clades a contiguity floor removes.

    PIEZO's S19 found a filtered scope trades one bias for another: the
    best-assembled genomes are disproportionately the clades of interest. The
    same question is asked here rather than assumed to transfer — and here it
    has a sharper form, because the margin species S4 added for their
    *proteome* gaps turn out to be very largely the same genomes whose
    assemblies cannot hold the gene (S5a), so a contiguity floor removes the
    part of the scope the scope was extended for.
    """
    rows = []
    genomes = {}
    for c in cells:
        genomes.setdefault(c["accession"], c)
    for floor in FLOOR_GRID:
        kept = [g for g in genomes.values() if g["contig_n50"] >= floor]
        if not kept:
            continue
        by_class: dict[str, int] = {}
        for g in kept:
            by_class[g["vclass"]] = by_class.get(g["vclass"], 0) + 1
        total = len(kept)
        for vclass in sorted({g["vclass"] for g in genomes.values()}):
            n_all = sum(1 for g in genomes.values() if g["vclass"] == vclass)
            n = by_class.get(vclass, 0)
            rows.append([floor, vclass, n_all, n,
                         round(n / max(1, n_all), 4), round(n / total, 4)])
    S.write_tsv(S.out_dir() / "floor_clade_composition.tsv",
                ["contig_n50_floor", "vclass", "n_genomes_in_scope",
                 "n_retained", "frac_of_class_retained",
                 "share_of_retained_scope"], rows)
    log(f"floor_clade_composition: {len(rows)} rows")
    return rows


# -------------------------------------------------------- 6. residual causes

def _consensus_flanks() -> dict[str, list[str]]:
    path = S.RESULTS / "synteny" / "flank_consensus.tsv"
    if not path.exists():
        return {}
    out: dict[str, list[str]] = {}
    for r in sorted(S.read_tsv(path),
                    key=lambda r: -S.fnum(r.get("fraction"), float, 0.0)):
        if (r.get("key_kind") != FLANK_KEY_KIND
                or r.get("window") != FLANK_WINDOW):
            continue
        if S.fnum(r.get("fraction"), float, 0.0) < FLANK_MIN_FRACTION:
            continue
        mine = out.setdefault(r["cell"], [])
        if len(mine) < FLANK_MAX_SYMBOLS:
            mine.append(r["symbol"].upper())
    return out


def neighbourhood_check(cells: list[dict], log=S.log) -> list[list]:
    """Is the gene missing, or is its whole neighbourhood missing?

    Contig N50 is genome-wide and a regional assembly defect is invisible to
    it, so no contiguity floor can catch one. The consensus flanks come from
    S8's committed table — the symbols that neighbour this paralog across the
    sweep — and a cell whose assembly carries fewer than half of them is one
    where the region is gone and the cell says nothing about the gene.
    """
    import s8_flank_lib as F

    consensus = _consensus_flanks()
    if not consensus:
        log("neighbourhood_check: no S8 flank consensus — skipped")
        return []
    rows = []
    cache: dict[str, object] = {}
    for c in cells:
        if not c["false_negative"] or c["cell"] not in consensus:
            continue
        acc = c["accession"]
        if acc not in cache:
            try:
                cache[acc] = F.load_genes(acc)
            except Exception:                                # pragma: no cover
                cache[acc] = None
        genes = cache[acc]
        symbols = consensus[c["cell"]]
        if not genes:
            rows.append([acc, c["organism"], c["vclass"], c["cell"],
                         c["status"], c["contig_n50"], len(symbols), "", "",
                         "no_gene_table"])
            continue
        present = {(g[5] or "").upper() for g in genes}
        hit = [s for s in symbols if s in present]
        frac = len(hit) / len(symbols)
        rows.append([acc, c["organism"], c["vclass"], c["cell"], c["status"],
                     c["contig_n50"], len(symbols), len(hit), round(frac, 3),
                     "neighbourhood_present" if frac >= NEIGHBOURHOOD_PRESENT
                     else "neighbourhood_missing"])
    S.write_tsv(S.out_dir() / "residual_neighbourhood.tsv",
                ["accession", "organism", "vclass", "cell", "status",
                 "contig_n50", "n_consensus_symbols", "n_present",
                 "frac_present", "verdict"], rows)
    miss = sum(1 for r in rows if r[9] == "neighbourhood_missing")
    log(f"residual_neighbourhood: {miss}/{len(rows)} false negatives sit "
        f"where the paralog's own neighbourhood is missing too")
    return rows


