"""s23_span_calibration.py — the two genomic thresholds, measured for S23.

S5a measured these for the vertebrates and both turned out to matter more than
a settings file suggests:

  * **miniprot's `-G`.** A max-intron smaller than the gene's largest intron
    does not lose the gene, it *splits* it, and a split ITPR reads out of the
    ledger as `fragment`. So `-G` is a threshold on the ledger's own call.
  * **D4's contiguity bar.** An assembly whose contigs are shorter than the
    gene cannot carry it, so an absence in such a genome is not evidence of
    absence. S5a set the bar at the median measured vertebrate ITPR span
    (142,212 bp) and 120 of 309 genomes failed it — the single largest caveat
    on S5b's results.

Neither number transfers. Vertebrate ITPR genes are 76-498 kb (S0); a
*Drosophila* Itpr is 22 kb. Carrying S5's bar into this sweep would fail
almost every invertebrate, protist and fungal assembly for want of contigs
they do not need, and would turn a working sweep into a page of caveats.

**How it is measured, and why not with miniprot.** The obvious route — run
the panel at a generous `-G` and measure what comes back — is circular: the
aligner's own max-intron shapes the loci it reports, so a measurement taken
that way cannot falsify the setting it is calibrating. This instead reads the
**annotation's** gene span, from NCBI's own gene records, for genes the census
independently calls ITPR. Two consequences:

  * the span is an *upper bound* on the largest intron (a gene spanning 22 kb
    has no 200 kb intron in it), which is exactly the direction `-G` needs —
    it makes the setting decidable in one direction, the useful one;
  * the panel is derived from `census_v5.tsv`, one species per band, by a
    positive test on the record's own gene symbol. The symbol is used only to
    *find* a gene to measure; the family call is the census's, made on
    architecture and profile, and the measurement is the annotation's.

Outputs -> results/s23_baits/span_calibration.{tsv,json}

Usage:
    python3 scripts/s23_span_calibration.py
    python3 scripts/s23_span_calibration.py --refresh
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from src.utils.data_root import get_data_root                 # noqa: E402
import s23_bait_spec as spec                                  # noqa: E402
from s4_manifest_lib import datasets_bin, load_jsonl          # noqa: E402
from s5_build_baits import read_tsv, write_tsv                # noqa: E402
from s5_classify import name_family                           # noqa: E402

CENSUS = PROJECT_ROOT / "results" / "census_v5" / "census_v5.tsv"
OUT_DIR = PROJECT_ROOT / "results" / "s23_baits"

#: Gene symbols that name this family. A positive test on the label, used only
#: to choose which annotated gene to measure. `itr-` is the nematode symbol
#: (*C. elegans* itr-1) and `ipla` the amoebozoan one — S0 found the family's
#: own defining Pfam does not even reach the latter, so a name list built from
#: the vertebrate symbols alone would have measured no protist gene at all.
GENE_PATTERNS = (r"^itpr\d?[a-z]?$", r"^itr-?\d$", r"^ipla$", r"^ip3r\d?$",
                 r"^insp3r\d?$", r"^itprb?$")
_GENE_RE = re.compile("|".join(GENE_PATTERNS), re.IGNORECASE)

#: The second route into the panel: the record's **protein name** says
#: inositol trisphosphate receptor even though its gene field is a locus tag.
#: Outside the model organisms that is the normal case — a symbol list alone
#: reaches 8 bands, and the protein names reach 28, including every protist
#: and algal band where the span is least predictable from the vertebrates.
#: The locus tag is still what NCBI is queried by, because for these
#: assemblies the locus tag *is* the gene symbol.
def names_the_family(row: dict) -> str:
    """Which route admits this record to the panel, or "" if neither."""
    if _GENE_RE.match((row.get("gene") or "").strip()):
        return "gene_symbol"
    if name_family(row.get("protein_name") or "") == "ITPR":
        return "protein_name"
    return ""

#: How many species per band to measure. One is a species; three is a spread.
PER_BAND = 3

#: The quantile of measured spans that sets the contiguity bar. The **median**,
#: matching S5a: half the family's genes are shorter than the bar, so a contig
#: at the bar carries a typical gene. Deliberately not the maximum, which would
#: fail assemblies that can hold every gene but the largest.
BAR_QUANTILE = 0.50

#: Safety factor on the largest measured span when deriving `-G`. A gene's
#: span bounds its largest intron, so the largest span in the panel is already
#: an upper bound for the panel; the factor covers the species not in it.
G_SAFETY = 2.0

#: miniprot's own default, which is what the rule is being checked against.
MINIPROT_DEFAULT_G = 200_000


def cache_dir() -> Path:
    d = get_data_root() / "raw_api" / "ncbi_datasets" / "s23_genes"
    d.mkdir(parents=True, exist_ok=True)
    return d


def query_gene(symbol: str, taxid: str, refresh: bool = False) -> list[dict]:
    """`datasets summary gene symbol` for one species, cached per query."""
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{taxid}_{symbol}")
    path = cache_dir() / f"{safe}.jsonl"
    if not path.exists() or refresh:
        proc = subprocess.run(
            [datasets_bin(), "summary", "gene", "symbol", symbol,
             "--taxon", str(taxid), "--as-json-lines"],
            capture_output=True, text=True)
        path.write_text(proc.stdout if proc.returncode == 0 else "")
    if not path.stat().st_size:
        return []
    return load_jsonl(path)


def spans_from_record(rec: dict) -> list[dict]:
    """Every annotated genomic span of one gene record."""
    gene = rec.get("gene", rec)
    out = []
    for ann in gene.get("annotations") or []:
        for loc in ann.get("genomic_locations") or []:
            rng = loc.get("genomic_range") or {}
            try:
                begin, end = int(rng["begin"]), int(rng["end"])
            except (KeyError, TypeError, ValueError):
                continue
            out.append({
                "symbol": gene.get("symbol", ""),
                "gene_id": gene.get("gene_id", ""),
                "taxid": gene.get("tax_id", ""),
                "organism": gene.get("taxname", ""),
                "assembly": ann.get("assembly_accession", ""),
                "annotation": ann.get("annotation_name", ""),
                "sequence": loc.get("sequence_name", ""),
                "start": min(begin, end), "end": max(begin, end),
                "span_bp": abs(end - begin) + 1,
            })
    return out


def select_panel(census: list[dict]) -> list[dict]:
    """One to `PER_BAND` census records per band whose gene symbol names it.

    Ranked by profile score, so the gene measured in each band is the
    best-supported family member there rather than the first one listed.
    """
    by_band: dict[str, list[dict]] = defaultdict(list)
    for r in census:
        if r.get("call") != spec.FAMILY_ITPR:
            continue
        band = spec.band_of(r)
        route = names_the_family(r) if band else ""
        if not route or not (r.get("gene") or "").strip():
            continue
        if not (r.get("taxon_id") or "").strip():
            continue
        by_band[band].append(dict(r, route=route))
    panel = []
    for band in spec.BAND_ORDER:
        rows = by_band.get(band, [])
        # gene-symbol records first: an `Itpr` is a surer measurement than a
        # locus tag whose protein name happens to say the family.
        rows.sort(key=lambda r: (r["route"] != "gene_symbol",
                                 -float(r.get("itpr_score") or 0)))
        seen = set()
        for r in rows:
            if r["taxon_id"] in seen:
                continue
            seen.add(r["taxon_id"])
            panel.append(dict(r, band=band))
            if len(seen) >= PER_BAND:
                break
    return panel


def derive(spans: list[dict]) -> dict:
    """The two thresholds, globally **and per group**.

    S5a could take one contiguity bar for all 309 vertebrate genomes because
    S0 had measured the family's span variation there at 6.5x. Outside the
    vertebrates it is ~100x — an *Octopus* ITPR spans 325 kb and a *Perkinsus*
    one 3.7 kb — so a single bar would clear a metazoan assembly whose contigs
    are a quarter of the gene it is being asked about. The bar is therefore
    per group, and the global median is only the fallback for a group with no
    measurement of its own.

    `-G` moves the other way: it is taken per group but **never below
    miniprot's own default**, because a `-G` smaller than an unmeasured
    species' largest intron splits its gene, and this panel has 41
    unmeasurable records against 16 measured. Erring high costs some chained
    junk that `s5_sweep_lib.MIN_LOCUS_IDENTITY` already filters; erring low
    costs genes.
    """
    vals = sorted(s["span_bp"] for s in spans)
    if not vals:
        raise SystemExit("no gene span was measured; the calibration cannot "
                         "set a threshold it did not measure")
    med = int(statistics.median(vals))
    largest = vals[-1]

    by_group: dict[str, list[int]] = defaultdict(list)
    for s in spans:
        by_group[s["group"]].append(s["span_bp"])
    groups = {}
    for g, xs in sorted(by_group.items()):
        xs = sorted(xs)
        groups[g] = {
            "n": len(xs), "min_bp": xs[0], "median_bp": int(statistics.median(xs)),
            "max_bp": xs[-1],
            "contiguity_bar_bp": int(statistics.median(xs)),
            "max_intron_bp": max(MINIPROT_DEFAULT_G, int(xs[-1] * G_SAFETY)),
        }
    return {
        "n_spans": len(vals), "n_species": len({s["taxid"] for s in spans}),
        "n_bands": len({s["band"] for s in spans}),
        "min_span_bp": vals[0], "median_span_bp": med,
        "p90_span_bp": vals[int(len(vals) * 0.9)], "max_span_bp": largest,
        "contiguity_bar_bp": med,
        "bar_quantile": BAR_QUANTILE,
        "max_intron_upper_bound_bp": largest,
        "derived_max_intron_bp": max(MINIPROT_DEFAULT_G,
                                     int(largest * G_SAFETY)),
        "miniprot_default_g": MINIPROT_DEFAULT_G,
        "default_g_adequate": bool(largest <= MINIPROT_DEFAULT_G),
        "g_safety_factor": G_SAFETY,
        "by_group": groups,
        "s5_vertebrate_bar_bp": 142_212,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    census = read_tsv(CENSUS)
    panel = select_panel(census)
    print(f"panel: {len(panel)} census records across "
          f"{len({r['band'] for r in panel})} bands")

    rows, missed = [], []
    for r in panel:
        recs = query_gene(r["gene"], r["taxon_id"], args.refresh)
        got = [s for rec in recs for s in spans_from_record(rec)]
        # NCBI returns every gene of that symbol in the taxon; keep the ones
        # whose taxid is the species we asked about.
        got = [s for s in got if str(s["taxid"]) == str(r["taxon_id"])]
        if not got:
            missed.append({"accession": r["accession"], "band": r["band"],
                           "gene": r["gene"], "species": r.get("species", ""),
                           "taxid": r["taxon_id"],
                           "why": "no annotated gene of that symbol at NCBI"})
            continue
        best = max(got, key=lambda s: s["span_bp"])
        rows.append(dict(best, band=r["band"], accession=r["accession"],
                         route=r["route"],
                         group=spec.group_of_band(r["band"]),
                         census_species=r.get("species", ""),
                         itpr_score=r.get("itpr_score", ""),
                         length_aa=r.get("length", ""),
                         n_locations=len(got)))
        print(f"  {r['band']:20s} {best['organism'][:30]:30s} "
              f"{best['symbol']:8s} {best['span_bp']:>10,} bp")

    stats = derive(rows)
    write_tsv(OUT_DIR / "span_calibration.tsv",
              ["band", "group", "organism", "census_species", "symbol",
               "route", "gene_id", "taxid", "accession", "length_aa",
               "itpr_score",
               "assembly", "annotation", "sequence", "start", "end",
               "span_bp", "n_locations"], rows)
    write_tsv(OUT_DIR / "span_calibration_missed.tsv",
              ["accession", "band", "gene", "species", "taxid", "why"], missed)
    (OUT_DIR / "span_calibration.json").write_text(json.dumps(stats, indent=1))

    print(f"\n{stats['n_spans']} genes in {stats['n_species']} species across "
          f"{stats['n_bands']} bands ({len(missed)} unmeasurable)")
    print(f"  spans {stats['min_span_bp']:,} - {stats['max_span_bp']:,} bp, "
          f"median {stats['median_span_bp']:,}")
    print(f"  global contiguity bar -> {stats['contiguity_bar_bp']:,} bp "
          f"(S5's vertebrate bar: {stats['s5_vertebrate_bar_bp']:,} bp)")
    print(f"  largest span {stats['max_span_bp']:,} bp bounds every intron; "
          f"miniprot's default -G {MINIPROT_DEFAULT_G:,} is "
          + ("adequate" if stats["default_g_adequate"] else "NOT adequate")
          + f" -> -G {stats['derived_max_intron_bp']:,}")
    print("  per group (bar / -G):")
    for g, s in stats["by_group"].items():
        print(f"    {g:16s} n={s['n']:2d}  {s['min_bp']:>9,} - "
              f"{s['max_bp']:>9,} bp   bar {s['contiguity_bar_bp']:>9,}   "
              f"-G {s['max_intron_bp']:>9,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
