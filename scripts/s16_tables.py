"""S16's committed stats file and the headline numbers the report reads.

`headline()` exists for the reason D13 exists: a report that computes its own
numbers can drift from the tables it cites. Everything the report states as a
result comes back through here, out of the committed TSVs.
"""

from __future__ import annotations

import collections
import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_paralogon as P                                      # noqa: E402
import s16_teleost as T                                        # noqa: E402

TABLES = [
    "loci.tsv", "merges.tsv", "copy_number.tsv", "copy_number_wide.tsv",
    "copy_number_by_clade.tsv", "copy_sensitivity.tsv", "expansions.tsv",
    "paralogy_map.tsv", "symbol_to_ensembl.tsv", "windows.tsv",
    "paralogy_links.tsv", "paralogon_test.tsv", "paralogon_blocks.tsv",
    "pooled_test.tsv", "quartet_test.tsv", "paralogon_species.tsv",
    "paralogon_species_summary.tsv", "replication_null.tsv",
    "teleost_copies.tsv", "dispersion.tsv", "dcs_flanks.tsv", "dcs_test.tsv",
    "r3_summary.tsv", "block_assignments.tsv",
]

#: The window the 2R answer is read at. +/-10 is S8's own window and the
#: highest-power setting available (the null mean rises with window size, so
#: a wider window buys genes and loses signal); +/-20 and +/-30 are committed
#: beside it as the sensitivity.
PRIMARY_WINDOW = 10
PRIMARY_LEVEL_SET = "2R_window"


def _rows(out: Path, name: str) -> list[dict]:
    p = out / name
    return L.read_tsv(p) if p.exists() else []


def headline(out: Path | None = None) -> dict:
    out = out or L.OUT_DIR
    h: dict = {}

    loci = _rows(out, "loci.tsv")
    h["n_loci"] = len(loci)
    h["n_merges"] = len(_rows(out, "merges.tsv"))
    h["n_copies"] = sum(1 for r in loci if r["is_copy"] == "True")
    h["n_lesion_rich_copies"] = sum(1 for r in loci if r["is_copy"] == "True"
                                    and r["integrity"] == "lesion_rich")

    cn = _rows(out, "copy_number.tsv")
    h["n_genomes"] = len({r["accession"] for r in cn})
    for cell in L.PARALOGS + (L.CONTROL_CELL,):
        sub = [r for r in cn if r["cell"] == cell]
        above = [r for r in sub if r["contig_spans_gene"] == "1"]
        h[f"copies_{cell}"] = sum(int(r["n_copies"]) for r in sub)
        h[f"multi_{cell}_above_bar"] = sum(1 for r in above
                                           if int(r["n_copies"]) > 1)
    # One denominator, computed once: every cell is scored in every genome,
    # so a per-cell `n_above_bar` written inside the loop would silently
    # report the last cell's count for all four.
    h["n_above_bar"] = len({r["accession"] for r in cn
                            if r["contig_spans_gene"] == "1"})

    # -- 2R
    tests = _rows(out, "paralogon_test.tsv")
    prim = [r for r in tests
            if r["level_set"] == PRIMARY_LEVEL_SET
            and int(r["window_n"]) == PRIMARY_WINDOW]
    h["paralogon_primary"] = {
        f"{r['group_a']}_vs_{r['group_b']}": {
            "class": r["pair_class"], "links": int(r["observed_links"]),
            "null_mean": float(r["null_mean"]), "p": float(r["p_permutation"]),
            "q_stratum": float(r["q_stratum"])}
        for r in prim}
    pooled = _rows(out, "pooled_test.tsv")
    h["pooled_primary"] = {
        r["family"]: {"links": int(r["observed_links"]),
                      "null_mean": float(r["null_mean"]),
                      "p": float(r["p_permutation"])}
        for r in pooled
        if r["level_set"] == PRIMARY_LEVEL_SET
        and int(r["window_n"]) == PRIMARY_WINDOW}
    links = _rows(out, "paralogy_links.tsv")
    h["itpr_links"] = sorted({(r["symbol_a"], r["symbol_b"],
                               r["duplication_level"])
                              for r in links
                              if r["pair_class"] == "ITPR_vs_ITPR"})
    h["itpr_links_2R"] = sorted({(r["symbol_a"], r["symbol_b"],
                                  r["duplication_level"])
                                 for r in links
                                 if r["pair_class"] == "ITPR_vs_ITPR"
                                 and r["in_2R_window"] == "1"})
    h["replication"] = {r["pair"]: {
        "n": int(r["n_genomes"]), "hit": int(r["n_with_link"]),
        "frac": float(r["frac_with_link"]),
        "hit_2R": int(r["n_with_2R_link"]),
        "frac_2R": float(r["frac_with_2R_link"]),
        "classes": int(r["n_classes"]),
        "null_frac": float(r["null_frac"]),
        "p": float(r["p_fisher_vs_null"])}
        for r in _rows(out, "paralogon_species_summary.tsv")}
    nul = _rows(out, "replication_null.tsv")
    h["replication_null"] = {"pairs": len(nul),
                             "with_link": sum(1 for r in nul
                                              if int(r["n_links"]) > 0)}

    # -- 3R
    tel = _rows(out, "teleost_copies.tsv")
    h["teleost_groups"] = {}
    for grp in ("pre_3R_outgroup", "teleost", "extra_wgd"):
        sub = [r for r in tel
               if r["group"] == grp and r["contig_spans_gene"] == "1"]
        if not sub:
            continue
        h["teleost_groups"][grp] = {
            "n": len(sub),
            **{p: round(statistics.mean(int(r[f"{p}_copies"]) for r in sub), 3)
               for p in L.PARALOGS},
            "RYR": round(statistics.mean(int(r["RYR_copies"]) for r in sub), 3),
            f"frac_{T.FOCAL}_multi": round(
                sum(1 for r in sub if int(r[f"{T.FOCAL}_copies"]) > 1)
                / len(sub), 4),
        }
    disp = [r for r in _rows(out, "dispersion.tsv")
            if r["cell"] == T.FOCAL and r["group"] == "teleost"]
    h["dispersion"] = disp[0] if disp else {}
    dcs = _rows(out, "dcs_test.tsv")
    h["dcs"] = {
        "n_pairs": len(dcs),
        "both_tetrapod": sum(1 for r in dcs if r["both_share_tetrapod"] == "1"),
        "disjoint_tetrapod": sum(1 for r in dcs
                                 if r["tetrapod_partition_disjoint"] == "1"),
        "both_fish": sum(1 for r in dcs if r["both_share_fish"] == "1"),
        "disjoint_fish": sum(1 for r in dcs
                             if r["fish_partition_disjoint"] == "1"),
    }
    r3rows = _rows(out, "r3_summary.tsv")
    h["r3"] = {r["metric"]: r["value"] for r in r3rows}
    # The note column carries the anchor names; keeping only the value would
    # have the report print "six anchors from six orders (6)".
    h["r3_notes"] = {r["metric"]: r["note"] for r in r3rows}
    return h


def stats(out: Path | None = None, params: dict | None = None,
          self_test: str = "unknown") -> Path:
    out = out or L.out_dir()
    notes = out / "paralogy_map_notes.md"
    release = ""
    if notes.exists():
        for line in notes.read_text().splitlines():
            if line.startswith("- release served by that host"):
                release = line.split("**")[1] if "**" in line else ""
    payload = {
        "task": "S16",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "self_test": self_test,
        "parameters": {
            "copy_coverage_bar": L.COV_FULL,
            "copy_min_aligned_aa": L.MIN_COPY_ALIGNED_AA,
            "copy_min_identity": L.MIN_COPY_IDENTITY,
            "merge_gap_bp": L.MERGE_GAP_BP,
            "merge_max_query_overlap": L.MERGE_MAX_QUERY_OVERLAP,
            "window_sizes": list(P.WINDOW_SIZES),
            "primary_window": PRIMARY_WINDOW,
            "primary_level_set": PRIMARY_LEVEL_SET,
            "levels_2R_core": list(P.LEVELS_2R_CORE),
            "levels_2R_window": list(P.LEVELS_2R_WINDOW),
            "n_permutations": P.N_PERMUTATIONS,
            "seed": P.SEED,
            "biomart_release": release,
            "biomart_host": P.M.ARCHIVE_HOST if hasattr(P, "M") else "",
            "pre_3R_orders": sorted(T.PRE_3R_ORDERS),
            "extra_wgd_orders": sorted(T.EXTRA_WGD_ORDERS),
            "focal_cell": T.FOCAL,
            "control_cell": L.CONTROL_CELL,
            **(params or {}),
        },
        "headline": headline(out),
        "sha256": {name: L.sha256(out / name)
                   for name in TABLES if (out / name).exists()},
    }
    p = out / "duplication_stats.json"
    p.write_text(json.dumps(payload, indent=1, default=str))
    return p


def main() -> None:
    print(stats())


if __name__ == "__main__":
    main()
