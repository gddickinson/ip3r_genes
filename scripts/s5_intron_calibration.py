"""s5_intron_calibration.py — measure real ITPR/RYR introns, set miniprot's -G.

miniprot's `-G` (max intron) defaults to 200 kb. A gene whose largest intron
exceeds it is not missed — it is **split**, and a split gene reads out of the
S5 ledger as two `fragment` loci instead of one `found` gene. So `-G` is a
threshold on the ledger's own call, and it has to come from this family's
measured intron sizes rather than from the PIEZO port's.

S0 already measured the *spans* (`results/s0_baseline/gene_structure.tsv`:
ITPR1 354 kb, ITPR2 498 kb, ITPR3 76 kb in human) and found the 6.5x
asymmetry between paralogs of near-identical protein length. A span is not an
intron, though: 498 kb over 57 exons could be evenly spread. This asks the
question the span cannot answer — how big is the *largest single intron* —
for both families across a vertebrate panel, and reports the scaling of that
number with genome size, which is what a formula for unseen genomes needs.

Source: Ensembl REST `lookup/symbol/{species}/{symbol}?expand=1`, whose exon
coordinates are already in the payload S0 used for exon counts. S0 established
that `/xrefs/symbol/homo_sapiens/` stalls while `lookup/symbol` works, so this
uses `lookup/symbol` only.

Outputs -> results/s5_baits/intron_calibration.tsv + intron_calibration.json

Usage:
  python scripts/s5_intron_calibration.py
  python scripts/s5_intron_calibration.py --species homo_sapiens,gallus_gallus
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

OUT_DIR = PROJECT_ROOT / "results" / "s5_baits"
ENSEMBL = "https://rest.ensembl.org"
TIMEOUT = 120

#: The panel. Chosen to span the sweep's clade bands wherever Ensembl carries
#: the species, so the genome-size scaling is measured across a real range
#: (zebrafish 1.4 Gbp -> human 3.1 Gbp -> axolotl-scale genomes are absent
#: from Ensembl main, which is why the formula still has to extrapolate).
SPECIES = [
    ("homo_sapiens", "mammalia"),
    ("mus_musculus", "mammalia"),
    ("gallus_gallus", "aves"),
    ("taeniopygia_guttata", "aves"),
    ("anolis_carolinensis", "reptilia"),
    ("xenopus_tropicalis", "amphibia"),
    ("danio_rerio", "actinopteri"),
    ("oryzias_latipes", "actinopteri"),
    ("latimeria_chalumnae", "sarcopterygian_fish"),
    ("callorhinchus_milii", "chondrichthyes"),
    ("petromyzon_marinus", "cyclostomata"),
]

SYMBOLS = ["ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "RYR3"]


def genome_size(species: str) -> int:
    """Total assembly length, so the intron rule can be expressed as a ratio.

    A max-intron threshold has to extrapolate: the sweep's biggest assemblies
    (Protopterus 40 Gbp, Lissotriton 23 Gbp) are nowhere near Ensembl, and
    intron length scales with genome size. Measuring the ratio here is what
    turns 11 measured genomes into a rule for 309.
    """
    data = get_json(f"{ENSEMBL}/info/assembly/{species}?")
    if not isinstance(data, dict):
        return 0
    total = data.get("base_pairs") or data.get("golden_path")
    return int(total or 0)


def get_json(url: str, tries: int = 3) -> dict | list | None:
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as fh:
                return json.loads(fh.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None                    # species has no such symbol
            time.sleep(2 * (attempt + 1))
        except Exception:                      # noqa: BLE001 — timeouts, resets
            time.sleep(2 * (attempt + 1))
    return None


def introns_of(transcript: dict) -> list[int]:
    """Intron lengths of one transcript, from its exon coordinates."""
    exons = sorted(((int(e["start"]), int(e["end"]))
                    for e in transcript.get("Exon", [])), key=lambda x: x[0])
    return [b1 - a2 - 1 for (_, a2), (b1, _) in zip(exons, exons[1:])
            if b1 - a2 - 1 > 0]


def measure(species: str, symbol: str) -> dict | None:
    url = f"{ENSEMBL}/lookup/symbol/{species}/{symbol}?expand=1"
    data = get_json(url)
    if not isinstance(data, dict) or "Transcript" not in data:
        return None
    coding = [t for t in data["Transcript"]
              if t.get("biotype") == "protein_coding"]
    if not coding:
        return None
    # The canonical transcript where Ensembl marks one, else the transcript
    # with the most exons — the same fallback s0_gene_structure.py uses.
    canonical = next((t for t in coding if t.get("is_canonical")), None)
    if canonical is None:
        canonical = max(coding, key=lambda t: len(t.get("Exon", [])))
    ints = introns_of(canonical)
    # The widest intron over *every* coding transcript, because miniprot has
    # to span whichever exon pair its own alignment picks, not just the
    # canonical one's.
    all_ints = [i for t in coding for i in introns_of(t)]
    return {
        "species": species, "symbol": symbol,
        "gene_id": data.get("id", ""),
        "seq_region": str(data.get("seq_region_name", "")),
        "span_bp": int(data.get("end", 0)) - int(data.get("start", 0)) + 1,
        "n_exons_canonical": len(canonical.get("Exon", [])),
        "n_coding_transcripts": len(coding),
        "max_intron_canonical": max(ints) if ints else 0,
        "median_intron_canonical": (sorted(ints)[len(ints) // 2] if ints else 0),
        "max_intron_any_transcript": max(all_ints) if all_ints else 0,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", help="comma-separated subset of the panel")
    ap.add_argument("--symbols", help="comma-separated subset of the symbols")
    args = ap.parse_args()

    panel = SPECIES
    if args.species:
        keep = {s.strip() for s in args.species.split(",")}
        panel = [(s, b) for s, b in panel if s in keep]
    symbols = ([s.strip() for s in args.symbols.split(",")]
               if args.symbols else SYMBOLS)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows, misses = [], []
    for species, band in panel:
        gbp = genome_size(species)
        print(f"{species}  ({gbp:,} bp)", flush=True)
        for symbol in symbols:
            t0 = time.time()
            rec = measure(species, symbol)
            if rec is None:
                misses.append(f"{species}/{symbol}")
                print(f"  -  {species:<24} {symbol:<6} not found "
                      f"[{time.time() - t0:.0f}s]", flush=True)
                continue
            rec["band"] = band
            rec["family"] = "RYR" if symbol.startswith("RYR") else "ITPR"
            rec["genome_bp"] = gbp
            rec["intron_per_gbp"] = (round(rec["max_intron_any_transcript"]
                                           / (gbp / 1e9), 1) if gbp else 0)
            rows.append(rec)
            print(f"  ok {species:<24} {symbol:<6} "
                  f"span {rec['span_bp']:>9,}  max intron "
                  f"{rec['max_intron_canonical']:>8,} (canonical) / "
                  f"{rec['max_intron_any_transcript']:>8,} (any) "
                  f"[{time.time() - t0:.0f}s]", flush=True)

    cols = ["species", "band", "family", "symbol", "gene_id", "seq_region",
            "genome_bp", "span_bp", "n_exons_canonical", "n_coding_transcripts",
            "median_intron_canonical", "max_intron_canonical",
            "max_intron_any_transcript", "intron_per_gbp"]
    with open(OUT_DIR / "intron_calibration.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["family"], r["symbol"], r["species"])):
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    itpr = [r for r in rows if r["family"] == "ITPR"]
    ryr = [r for r in rows if r["family"] == "RYR"]

    def summarise(group: list[dict], label: str) -> dict:
        if not group:
            return {}
        canon = [r["max_intron_canonical"] for r in group]
        anyt = [r["max_intron_any_transcript"] for r in group]
        worst = max(group, key=lambda r: r["max_intron_any_transcript"])
        over = [r for r in group if r["max_intron_any_transcript"] > 200_000]
        scaled = [r for r in group if r["intron_per_gbp"]]
        steepest = (max(scaled, key=lambda r: r["intron_per_gbp"])
                    if scaled else None)
        print(f"\n{label}: n={len(group)}  "
              f"max intron (canonical) median {sorted(canon)[len(canon) // 2]:,} "
              f"max {max(canon):,}")
        print(f"    over miniprot's 200 kb default (any transcript): "
              f"{len(over)}/{len(group)}")
        print(f"    widest: {worst['species']} {worst['symbol']} "
              f"{worst['max_intron_any_transcript']:,} bp "
              f"({worst['genome_bp'] / 1e9:.1f} Gbp genome)")
        if steepest:
            print(f"    steepest scaling: {steepest['species']} "
                  f"{steepest['symbol']} "
                  f"{steepest['intron_per_gbp']:,.0f} bp of intron per Gbp")
        return {"n": len(group),
                "median_max_intron": sorted(canon)[len(canon) // 2],
                "max_max_intron": max(canon), "max_any_transcript": max(anyt),
                "n_over_200kb": len(over),
                "widest": f"{worst['species']}/{worst['symbol']}",
                "widest_bp": worst["max_intron_any_transcript"],
                "max_intron_per_gbp": (steepest["intron_per_gbp"]
                                       if steepest else 0),
                "steepest": f"{steepest['species']}/{steepest['symbol']}"
                if steepest else ""}

    itpr_s, ryr_s = summarise(itpr, "ITPR"), summarise(ryr, "RYR")

    # ---- the -G rule the sweep uses.
    #
    # Not the per-Gbp ratio: that statistic is largest in the *smallest*
    # genomes measured (Petromyzon, 0.9 Gbp, an 81 kb RYR2 intron -> 92 kb per
    # Gbp) and extrapolating it linearly to a 40 Gbp lungfish asks for a 7.4 Mbp
    # -G, which is an artefact of dividing by a small denominator rather than a
    # measurement of anything.
    #
    # Instead: anchor on the widest intron actually measured, in the genome it
    # was measured in, and scale from there. Both families are in the anchor
    # because the RyR control baits are in the same miniprot run — but note
    # the asymmetry in what a too-small -G costs. A split ITPR reads out of the
    # ledger as `fragment` and corrupts the deliverable; a split RyR still
    # shows RyR present, so the control survives it.
    widest = max(itpr_s.get("widest_bp", 0), ryr_s.get("widest_bp", 0))
    anchor_bp = max((r["genome_bp"] for r in rows
                     if r["max_intron_any_transcript"] == widest), default=0)
    scale = int(round(2 * widest, -5))          # 2x safety at the anchor
    floor = 200_000                             # miniprot's own default
    cap = 2_000_000                             # tractability, not biology
    rule = {"widest_measured_bp": widest, "anchor_genome_bp": anchor_bp,
            "scale_bp": scale, "floor_bp": floor, "cap_bp": cap,
            "safety_factor": 2.0,
            "formula": "min(cap, max(floor, scale_bp * genome_bp / anchor_genome_bp))"}

    def g_for(genome_bp: float) -> int:
        return min(cap, max(floor, int(scale * genome_bp / anchor_bp)))

    # The floor stays at miniprot's own 200 kb default rather than at the
    # anchor's 500 kb: a 1.1 Gbp bird genome's widest measured intron is 93 kb
    # (RYR2) and 39 kb for any ITPR, so 200 kb is already 2x the worst case
    # there. Forcing every small genome up to the human-scale -G would buy no
    # recall and would only make spurious long-range joins easier.
    print(f"\n-G rule: min({cap:,}, max({floor:,}, {scale:,} x genome_bp / "
          f"{anchor_bp:,}))")
    print(f"    scale = 2x the widest intron measured ({widest:,} bp, in a "
          f"{anchor_bp / 1e9:.1f} Gbp genome)")
    for gbp, who in ((1.0, ""), (3.1, "human"), (10.0, "Bombina bombina"),
                     (23.2, "Lissotriton helveticus"),
                     (40.1, "Protopterus annectens")):
        val = g_for(gbp * 1e9)
        print(f"    {gbp:>5.1f} Gbp {who:<24} -> {val:>9,} bp"
              + ("   [capped]" if val == cap else ""))

    summary = {"itpr": itpr_s, "ryr": ryr_s, "max_intron_rule": rule,
               "n_rows": len(rows), "misses": misses}
    (OUT_DIR / "intron_calibration.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(f"\nwrote {OUT_DIR / 'intron_calibration.tsv'} ({len(rows)} rows, "
          f"{len(misses)} not found)")


if __name__ == "__main__":
    main()
