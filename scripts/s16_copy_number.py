"""S16 stage `copies` — the per-locus copy table, and what it says.

The brief's first step: a per-**locus** copy table, read out of every sweep
`summary.json` rather than out of the ledger, with split models merged so
that an aligner artefact cannot manufacture a duplication.

Three numbers are kept apart everywhere, because they are not the same
thing and the difference between them is where a duplication claim goes
wrong:

* **loci** — every cluster the sweep placed, secondary ones included;
* **copies** — loci whose own model covers a complete bait, after the merge;
* **intact copies** — those whose ORF is not lesion-rich on S15a's bar.

And every count carries D4's contiguity flag, because a second copy on a
second contig in an assembly whose contigs cannot carry the gene is the one
place a fragment reads as a duplicate.

The **RyR cell is the positive control** throughout. RYR1/2/3 are a
three-member vertebrate family of the same age and the same 2R candidacy as
ITPR1/2/3, and D14's hazard puts them inside every search this project runs
anyway. A rule that finds the teleost duplicates of one family and not the
other is a rule about the instrument.

Outputs (`results/duplication/`): `loci.tsv`, `merges.tsv`,
`copy_number.tsv`, `copy_number_wide.tsv`, `copy_number_by_clade.tsv`,
`copy_sensitivity.tsv`, `expansions.tsv`.
"""

from __future__ import annotations

import collections
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402

LOCI_HEADER = [
    "accession", "organism", "vclass", "vorder", "assembly_level",
    "contig_n50", "contig_spans_gene", "cell", "cell_status", "slot",
    "locus_idx", "locus_clade", "assigned_via", "contig", "start", "end",
    "strand", "bait", "identity", "coverage", "aligned_aa", "q_start",
    "q_end", "bait_len", "frameshifts", "stop_codons", "lesions",
    "lesion_density", "integrity", "family_margin", "paralog_margin",
    "contig_edge", "longest_n_run", "annot_gene", "annot_frac_cds",
    "annot_paralog_matches", "merged_models", "merged_from",
    "is_copy", "excluded_reason",
]

#: The coverage bars the copy count is measured at. The primary is
#: `s16_lib.COV_FULL`; the others exist so the sensitivity of every
#: duplication statement to that bar is a committed table and not a promise.
COV_BARS = (0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)


def build_loci(log=print) -> tuple[list[dict], list[dict]]:
    """Every locus, merged, with an integrity call and a copy verdict."""
    raw = L.iter_loci()
    recs, merges = L.merge_split_models(raw)
    bar, why = L.lesion_bar()
    for r in recs:
        r["integrity"] = L.integrity_call(r, bar)
        r["excluded_reason"] = L.copy_exclusion(r)
        r["is_copy"] = r["excluded_reason"] == ""
    log(f"[copies] {len(raw)} loci -> {len(recs)} after merging "
        f"{len(merges)} split models; lesion bar {bar:.4f} ({why})")
    return recs, merges


def loci_rows(recs: list[dict]) -> list[list]:
    return [[r.get(c, "") for c in LOCI_HEADER] for r in recs]


# ------------------------------------------------------- copy-number table

def copy_number(recs: list[dict], cov: float = L.COV_FULL
                ) -> dict[tuple[str, str], dict]:
    """(accession, cell) -> counts. Every swept genome gets a row per cell,
    including a zero one, because a count is only a count with a
    denominator."""
    gidx = L.genome_index()
    out: dict[tuple[str, str], dict] = {}
    for acc, g in gidx.items():
        for cell in L.PARALOGS + (L.CONTROL_CELL,):
            out[(acc, cell)] = {
                "accession": acc, "organism": g["organism"],
                "vclass": g["vclass"], "vorder": g["vorder"],
                "assembly_level": g["assembly_level"],
                "contig_n50": g["contig_n50"],
                "contig_spans_gene": g["contig_spans_gene"],
                "cell": cell, "n_loci": 0, "n_copies": 0,
                "n_intact_copies": 0, "n_contigs": 0, "copy_contigs": [],
                "max_coverage": 0.0, "annot_named": 0,
            }
    for r in recs:
        key = (r["accession"], r["cell"])
        if key not in out:
            continue
        c = out[key]
        c["n_loci"] += 1
        c["max_coverage"] = max(c["max_coverage"], r["coverage"])
        c["annot_named"] += int(bool(r["annot_gene"]))
        if L.is_copy(r, cov):
            c["n_copies"] += 1
            c["copy_contigs"].append(r["contig"])
            if r["integrity"] == "intact":
                c["n_intact_copies"] += 1
    for c in out.values():
        c["n_contigs"] = len(set(c["copy_contigs"]))
    return out


CN_HEADER = ["accession", "organism", "vclass", "vorder", "assembly_level",
             "contig_n50", "contig_spans_gene", "cell", "n_loci", "n_copies",
             "n_intact_copies", "n_contigs", "same_contig_copies",
             "max_coverage", "annot_named", "copy_contigs"]


def copy_number_rows(cn: dict[tuple[str, str], dict]) -> list[list]:
    rows = []
    for (acc, cell), c in sorted(cn.items()):
        rows.append([acc, c["organism"], c["vclass"], c["vorder"],
                     c["assembly_level"], c["contig_n50"],
                     int(c["contig_spans_gene"]), cell, c["n_loci"],
                     c["n_copies"], c["n_intact_copies"], c["n_contigs"],
                     int(c["n_copies"] > 1 and c["n_contigs"] < c["n_copies"]),
                     round(c["max_coverage"], 4), c["annot_named"],
                     ";".join(sorted(set(c["copy_contigs"])))])
    return rows


def wide_rows(cn: dict[tuple[str, str], dict]) -> tuple[list[str], list[list]]:
    header = (["accession", "organism", "vclass", "vorder",
               "contig_spans_gene"]
              + [f"{p}_copies" for p in L.PARALOGS]
              + ["ITPR_copies", "RYR_copies", "ITPR_intact", "RYR_intact"])
    by_acc: dict[str, dict] = collections.defaultdict(dict)
    for (acc, cell), c in cn.items():
        by_acc[acc][cell] = c
    rows = []
    for acc, cells in sorted(by_acc.items()):
        any_c = next(iter(cells.values()))
        counts = [cells.get(p, {}).get("n_copies", 0) for p in L.PARALOGS]
        ryr = cells.get(L.CONTROL_CELL, {})
        rows.append([acc, any_c["organism"], any_c["vclass"], any_c["vorder"],
                     int(any_c["contig_spans_gene"])] + counts
                    + [sum(counts), ryr.get("n_copies", 0),
                       sum(cells.get(p, {}).get("n_intact_copies", 0)
                           for p in L.PARALOGS),
                       ryr.get("n_intact_copies", 0)])
    return header, rows


# --------------------------------------------------------- retention shape

BY_CLADE_HEADER = ["vclass", "cell", "n_genomes", "n_above_bar",
                   "mean_copies", "mean_copies_above_bar", "median_copies",
                   "max_copies", "genomes_with_2plus", "genomes_with_0",
                   "frac_multi_above_bar"]


def by_clade(cn: dict[tuple[str, str], dict]) -> list[list]:
    """Copies per vertebrate class, with and without D4's contiguity bar.

    Both columns are printed because uncontrolled they say different
    things: a class of scaffold-level assemblies loses copies to the
    assembly, which reads as a lineage that lost the gene (S5's own
    annotation-quality section made the same mistake once and was fixed by
    holding contiguity constant)."""
    by: dict[tuple[str, str], list[dict]] = collections.defaultdict(list)
    for (_, cell), c in cn.items():
        by[(c["vclass"], cell)].append(c)
    rows = []
    for (vclass, cell), cells in sorted(by.items()):
        above = [c for c in cells if c["contig_spans_gene"]]
        counts = [c["n_copies"] for c in cells]
        counts_above = [c["n_copies"] for c in above]
        multi_above = sum(1 for n in counts_above if n > 1)
        rows.append([
            vclass, cell, len(cells), len(above),
            round(statistics.mean(counts), 3) if counts else 0,
            round(statistics.mean(counts_above), 3) if counts_above else 0,
            statistics.median(counts) if counts else 0,
            max(counts) if counts else 0,
            sum(1 for n in counts if n > 1), sum(1 for n in counts if n == 0),
            round(multi_above / len(above), 4) if above else 0.0,
        ])
    return rows


EXP_HEADER = ["accession", "organism", "vclass", "vorder",
              "contig_spans_gene", "cell", "n_copies", "n_contigs",
              "same_contig", "copy_id", "contig", "start", "end", "strand",
              "coverage", "identity", "integrity", "bait", "annot_gene"]


def expansions(recs: list[dict], cn: dict[tuple[str, str], dict],
               cov: float = L.COV_FULL) -> list[list]:
    """Every genome carrying more than one copy of a cell, one row per copy."""
    multi = {k for k, c in cn.items() if c["n_copies"] > 1}
    rows = []
    for r in sorted(recs, key=lambda x: (x["accession"], x["cell"],
                                         x["contig"], x["start"])):
        key = (r["accession"], r["cell"])
        if key not in multi or not L.is_copy(r, cov):
            continue
        c = cn[key]
        idx = sum(1 for x in rows if x[0] == r["accession"] and x[5] == r["cell"])
        rows.append([r["accession"], r["organism"], r["vclass"], r["vorder"],
                     int(r["contig_spans_gene"]), r["cell"], c["n_copies"],
                     c["n_contigs"],
                     int(c["n_contigs"] < c["n_copies"]),
                     f"{r['accession']}#{r['cell']}#copy{idx + 1}",
                     r["contig"], r["start"], r["end"], r["strand"],
                     round(r["coverage"], 4), round(r["identity"], 4),
                     r["integrity"], r["bait"], r["annot_gene"]])
    return rows


SENS_HEADER = ["cov_bar", "cell", "scope", "n_genomes", "total_copies",
               "genomes_with_2plus", "genomes_with_0", "mean_copies"]


def copy_sensitivity(recs: list[dict]) -> list[list]:
    """Every count in this task, recomputed at every coverage bar.

    A duplication claim that only survives one bar is a claim about the
    bar. The grid is committed so a reader can see which statements move
    (S15b's sensitivity matrix, applied to the one threshold S16 owns).
    """
    rows = []
    for cov in COV_BARS:
        cn = copy_number(recs, cov)
        for cell in L.PARALOGS + (L.CONTROL_CELL,):
            for scope, keep in (("all", lambda c: True),
                                ("above_contiguity_bar",
                                 lambda c: c["contig_spans_gene"])):
                cells = [c for (a, ce), c in cn.items()
                         if ce == cell and keep(c)]
                counts = [c["n_copies"] for c in cells]
                rows.append([cov, cell, scope, len(cells), sum(counts),
                             sum(1 for n in counts if n > 1),
                             sum(1 for n in counts if n == 0),
                             round(statistics.mean(counts), 4)
                             if counts else 0.0])
    return rows


MERGE_HEADER = ["accession", "cell", "contig", "strand", "kept_start",
                "kept_end", "merged_start", "merged_end", "gap_bp",
                "query_overlap", "kept_q", "merged_q", "why"]


def run(out: Path | None = None, log=print) -> dict:
    out = out or L.out_dir()
    recs, merges = build_loci(log)
    L.write_tsv(out / "loci.tsv", LOCI_HEADER, loci_rows(recs))
    L.write_tsv(out / "merges.tsv", MERGE_HEADER,
                [[m[c] for c in MERGE_HEADER] for m in merges])
    cn = copy_number(recs)
    L.write_tsv(out / "copy_number.tsv", CN_HEADER, copy_number_rows(cn))
    wh, wr = wide_rows(cn)
    L.write_tsv(out / "copy_number_wide.tsv", wh, wr)
    L.write_tsv(out / "copy_number_by_clade.tsv", BY_CLADE_HEADER, by_clade(cn))
    L.write_tsv(out / "expansions.tsv", EXP_HEADER, expansions(recs, cn))
    L.write_tsv(out / "copy_sensitivity.tsv", SENS_HEADER,
                copy_sensitivity(recs))
    n_copies = sum(1 for r in recs if r["is_copy"])
    for cell in L.PARALOGS + (L.CONTROL_CELL,):
        cells = [c for (a, ce), c in cn.items() if ce == cell]
        multi = sum(1 for c in cells if c["n_copies"] > 1)
        log(f"[copies] {cell:6s} {sum(c['n_copies'] for c in cells):5d} copies "
            f"in {len(cells)} genomes; {multi} genomes carry >1")
    return {"loci": recs, "merges": merges, "copy_number": cn,
            "n_copies": n_copies}


def main() -> None:
    run()


if __name__ == "__main__":
    main()
