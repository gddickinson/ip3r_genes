"""s5_rescue.py — whole-genome tblastn rescue and trace attribution.

When miniprot finds no locus for a cell, an absence claim needs a second,
independent search. This runs it, clusters the surviving HSPs into regions,
and re-scores each region against the **full** bait panel so a trace is
credited to whatever actually explains it.

Two attributions matter here and they are not the same question:

  *Family* — a piece of a ryanodine receptor must never be reported as the
  remnant of a missing ITPR. D14, and the reason the RyR baits are in the
  panel at all.

  *Paralog* — a fragment of the genome's own ITPR1 must not be reported as
  the remnant of a missing ITPR3. This is the harder of the two here:
  ITPR1/2/3 are 61-68 % identical (S1), against 40-50 % for the family the
  PIEZO port measured its margin on, so a trace discriminates between ITPR
  paralogs *less* well than that port's threshold assumed. The margin is
  therefore inherited as a floor rather than as a tuned value, it is recorded
  on every region so it can be re-cut, and a region inside the band is
  reported as `tblastn_trace_ambiguous` rather than resolved.

Rescue baits are derived from the committed panel, not listed here, so they
cannot drift from `bait_manifest.tsv`.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s5_sweep_lib import CLASSES, CONTROL_CLASS, RESCUE_E  # noqa: E402

REGION_MERGE = 20_000     # merge rescue HSPs this close into one region
REGION_PAD = 5_000        # flanking bp extracted for re-scoring
MAX_REGIONS = 12          # cap on regions re-scored per cell

#: Relative bit-score margin a region's top paralog must beat the runner-up
#: paralog by before the region counts as evidence for it, in D7's units:
#: `(top - second) / top`.
#:
#: **Measured, not inherited.** It arrived as the PIEZO port's 1.5x bit-score
#: ratio restated (0.333), from a family whose paralogs are 40-50 % identical
#: where these are 61-68 % (S1). `s5_calibrate_margin.py` measured what
#: separation is actually available: across the S5a pilot, at complete
#: full-length loci whose paralog identity the assembly's own annotation
#: establishes and whose call is correct, the paralog margin runs 0.225-0.347
#: with a median of 0.281 — and **7 of 9 sit below 0.333**. A complete locus
#: is an upper bound on a rescue fragment, so a threshold complete evidence
#: fails cannot be met by a fragment: under the inherited value every rescue
#: trace would be reported ambiguous and no absence claim could ever be
#: attributed.
#:
#: 0.22 is the highest threshold that rejects none of those correct calls,
#: rounded down to two places. Its limitation is stated rather than hidden:
#: it is derived from complete loci and applied to fragments, so it is a
#: ceiling on the right answer, not the right answer. S5b re-calibrates it
#: against rescue regions of known identity — a genome carrying two paralogs
#: and genuinely missing the third supplies exactly those. Every region
#: carries its own margin in the ledger, so any re-cut is possible without
#: re-running the sweep.
#: -> results/genome_ledger/margin_calibration.{tsv,json}
ATTRIBUTION_REL_MARGIN = 0.22

#: How many bands a cell's rescue baits should span. A rescue query is a
#: whole-genome tblastn, so it is the expensive search in the pipeline;
#: 3 baits from 3 different clade bands is the spread/cost trade.
RESCUE_BAITS_PER_CLASS = 3


def rescue_baits(bait_meta: dict[str, dict]) -> dict[str, list[str]]:
    """cell class -> bait ids, spread across clade bands (derived, not listed).

    Prefers one bait per band, crown-first, so a rescue query carries a
    mammal, a bird and a fish rather than three mammals.
    """
    from s5_bait_spec import BANDS

    out: dict[str, list[str]] = {}
    for cell in (*CLASSES, CONTROL_CLASS):
        pool = [(bid, m) for bid, m in bait_meta.items()
                if m["clade"] == cell]
        picked: list[str] = []
        for band in BANDS:
            hit = [bid for bid, m in pool if m["band"] == band
                   and bid not in picked]
            if hit:
                picked.append(sorted(hit)[0])
            if len(picked) >= RESCUE_BAITS_PER_CLASS:
                break
        if len(picked) < RESCUE_BAITS_PER_CLASS:      # thin band coverage
            picked += [bid for bid, _ in sorted(pool) if bid not in picked][
                :RESCUE_BAITS_PER_CLASS - len(picked)]
        out[cell] = picked
    return out


def region_margin(region: dict) -> tuple[str | None, float]:
    """(top paralog, its relative bit-score margin over the runner-up paralog).

    Ranked over the **canonical paralog clades only**. The panel's
    `vertebrate_basal` baits carry no paralog assignment — they are the gar,
    chimaera and lamprey seeds, which exist because those bands have no
    labelled record — so they are evidence that a region is an ITPR, not a
    rival hypothesis about which ITPR it is.

    Ranking them alongside ITPR1/2/3 is how the first pilot run got its
    margins: in *Todus mexicanus* an ITPR1 trace scored ITPR1 916.1 against
    `vertebrate_basal` 891.6 and ITPR2 700.6, and reporting the runner-up as
    the deep unlabelled bait collapsed a real 0.235 paralog margin to 0.027.
    Every ITPR1 and ITPR2 trace in that genome was deflated the same way.
    `clade_bits` still records the unlabelled scores — they are useful
    evidence, just not a competing paralog.
    """
    bits = region.get("clade_bits") or {}
    ranked = sorted(((c, b) for c, b in bits.items() if c in CLASSES),
                    key=lambda kv: -kv[1])
    if not ranked:
        return region.get("assigned_clade"), 0.0
    if len(ranked) == 1:
        return ranked[0][0], 1.0
    (top, top_bits), (_, second_bits) = ranked[0], ranked[1]
    return top, ((top_bits - second_bits) / top_bits if top_bits else 0.0)


def region_family(region: dict, bait_meta: dict[str, dict]) -> tuple[str, float]:
    """(winning family, its relative margin) — D14 before the paralog question.

    Collapsing the per-clade bits to per-family first is what stops three
    close ITPR paralog scores from being read as an ambiguous *family* call:
    the family question is usually decided even when the paralog is not.
    """
    bits = region.get("clade_bits") or {}
    per_family: dict[str, float] = {}
    for clade, b in bits.items():
        fam = CONTROL_CLASS if clade == CONTROL_CLASS else "ITPR"
        per_family[fam] = max(per_family.get(fam, 0.0), b)
    if not per_family:
        return "", 0.0
    ranked = sorted(per_family.items(), key=lambda kv: -kv[1])
    if len(ranked) == 1:
        return ranked[0][0], 1.0
    (top, tb), (_, sb) = ranked[0], ranked[1]
    return top, ((tb - sb) / tb if tb else 0.0)


def split_by_margin(regions: list[dict], cell: str, bait_meta: dict,
                    margin: float = ATTRIBUTION_REL_MARGIN
                    ) -> tuple[list, list, list]:
    """Regions for `cell` -> (decisive, ambiguous, attributed elsewhere).

    Family first (D14): a region the RyR baits win is attributed elsewhere
    however its ITPR paralog scores rank, and vice versa for the control cell.
    """
    want_family = CONTROL_CLASS if cell == CONTROL_CLASS else "ITPR"
    decisive, ambiguous, elsewhere = [], [], []
    for r in regions:
        fam, fam_margin = region_family(r, bait_meta)
        top, rel = region_margin(r)
        r["attribution_rel_margin"] = round(rel, 3)
        r["family_call"] = fam
        r["family_rel_margin"] = round(fam_margin, 3)
        if fam != want_family:
            elsewhere.append(r)
        elif top != cell:
            elsewhere.append(r)
        elif rel >= margin:
            decisive.append(r)
        else:
            ambiguous.append(r)
    return decisive, ambiguous, elsewhere


def run_makeblastdb(fna: Path, db_prefix: Path) -> None:
    if Path(str(db_prefix) + ".nin").exists() or Path(str(db_prefix) + ".00.nin").exists():
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        ["makeblastdb", "-in", str(fna), "-dbtype", "nucl",
         "-out", str(db_prefix)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {proc.stderr.strip()[:300]}")


OUTFMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"


def run_tblastn(query_faa: Path, db_prefix: Path, out_tsv: Path,
                threads: int = 8, evalue: float = 1e-3) -> None:
    proc = subprocess.run(
        ["tblastn", "-query", str(query_faa), "-db", str(db_prefix),
         "-evalue", str(evalue), "-num_threads", str(threads),
         "-outfmt", OUTFMT, "-out", str(out_tsv)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"tblastn failed: {proc.stderr.strip()[:300]}")


def parse_tblastn(out_tsv: Path) -> list[dict]:
    hsps = []
    if not out_tsv.exists():
        return hsps
    for line in open(out_tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        hsps.append({"bait": f[0], "contig": f[1], "pident": float(f[2]),
                     "length": int(f[3]), "sstart": int(f[6]), "send": int(f[7]),
                     "evalue": float(f[8]), "bits": float(f[9])})
    return hsps


def cluster_hsps(hsps: list[dict], e_cut: float = RESCUE_E,
                 merge: int = REGION_MERGE) -> list[dict]:
    """Group significant rescue HSPs into candidate regions per contig."""
    sig = [h for h in hsps if h["evalue"] <= e_cut]
    by_contig: dict[str, list[dict]] = {}
    for h in sig:
        by_contig.setdefault(h["contig"], []).append(h)
    regions = []
    for contig, group in by_contig.items():
        for h in group:
            h["_lo"], h["_hi"] = sorted((h["sstart"], h["send"]))
        group.sort(key=lambda h: h["_lo"])
        cur = None
        for h in group:
            if cur and h["_lo"] <= cur["end"] + merge:
                cur["end"] = max(cur["end"], h["_hi"])
                cur["hsps"].append(h)
            else:
                cur = {"contig": contig, "start": h["_lo"], "end": h["_hi"],
                       "hsps": [h]}
                regions.append(cur)
    for r in regions:
        best = min(r["hsps"], key=lambda h: h["evalue"])
        r.update({"n_hsps": len(r["hsps"]), "best_evalue": best["evalue"],
                  "best_pident": best["pident"], "best_bait": best["bait"],
                  "aligned_aa": sum(h["length"] for h in r["hsps"])})
        del r["hsps"]
    regions.sort(key=lambda r: r["best_evalue"])
    return regions


def attribute_region(fetch_fn, region: dict, all_baits: Path, tmp_dir: Path,
                     bait_meta: dict[str, dict]) -> dict:
    """Re-score one rescue region against the full bait panel."""
    seq = fetch_fn(region["contig"], region["start"] - REGION_PAD,
                   region["end"] + REGION_PAD)
    region["assigned_clade"] = None
    region["assigned_bait"] = None
    if not seq:
        return region
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{region['contig']}_{region['start']}".replace("|", "_")
    sub = tmp_dir / f"{tag}.fna"
    sub.write_text(f">{region['contig']}_region\n{seq}\n")
    out = tmp_dir / f"{tag}.tsv"
    proc = subprocess.run(
        ["tblastn", "-query", str(all_baits), "-subject", str(sub),
         "-evalue", "1e-3", "-outfmt", OUTFMT, "-out", str(out)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        region["attribution_error"] = proc.stderr.strip()[:200]
        return region
    per_bait: dict[str, float] = {}
    for h in parse_tblastn(out):
        per_bait[h["bait"]] = per_bait.get(h["bait"], 0.0) + h["bits"]
    if not per_bait:
        return region
    per_clade: dict[str, float] = {}
    for bait, bits in per_bait.items():
        clade = bait_meta.get(bait, {}).get("clade", "unknown")
        per_clade[clade] = max(per_clade.get(clade, 0.0), bits)
    region["assigned_clade"] = max(per_clade, key=per_clade.get)
    region["assigned_bait"] = max(per_bait, key=per_bait.get)
    region["clade_bits"] = {k: round(v, 1) for k, v in
                            sorted(per_clade.items(), key=lambda kv: -kv[1])}
    sub.unlink(missing_ok=True)
    out.unlink(missing_ok=True)
    return region


def annotate_regions(regions: list[dict], gene_index) -> None:
    if gene_index is None:
        return
    for r in regions:
        genes = gene_index.overlapping(r["contig"], r["start"], r["end"])
        r["genes"] = [{"name": g[5] or g[4], "biotype": g[6]} for g in genes[:6]]


def summarize_tblastn(hsps: list[dict], e_cut: float = RESCUE_E) -> dict | None:
    sig = [h for h in hsps if h["evalue"] <= e_cut]
    if not sig:
        return None
    best = min(sig, key=lambda h: h["evalue"])
    return {"n_hsps": len(sig), "best_evalue": best["evalue"],
            "best_contig": best["contig"], "best_pident": best["pident"],
            "aligned_aa": sum(h["length"] for h in sig),
            "contigs": sorted({h["contig"] for h in sig})}
