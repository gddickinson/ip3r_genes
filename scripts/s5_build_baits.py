"""s5_build_baits.py — build the S5 genomic-sweep bait panel from the census.

Applies `s5_bait_spec.py`'s seven rules and **enforces** every one of them, in
the shape `s3_build_seed.py` established: a candidate whose census call,
architecture count or length has drifted since the last build does not
silently become a bait, it fails the build. The panel is derived from
`results/census_v3/census_v3.tsv` rather than hand-listed, so it rebuilds when
the census does.

The chimera screen (R4/D5) is this project's own instrument: both S3 profiles
are run over the shortlist, a candidate must be assigned to its own family
under S3's relative margin, and the winning profile's envelope must cover the
bait without a large unaligned interior. That catches the failure a length
band cannot — a fusion, a read-through, a mis-joined gene model.

Sequences come from committed / archived files only: S3's seed FASTAs for the
reused seeds, and S2's archived seeded-space FASTA under the data root for the
additions. No network calls.

Outputs -> results/s5_baits/
  baits.faa            the panel, one record per bait
  bait_manifest.tsv    every bait with its rule provenance
  bait_screen.tsv      every candidate screened, passed or not, with the reason
  bait_build_stats.json  slot fill, SHA-256 of the panel, tool versions

Usage:
  python scripts/s5_build_baits.py
  python scripts/s5_build_baits.py --shortlist 5      # deeper screen per slot
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from src.utils.data_root import require_data_root          # noqa: E402
import s5_bait_spec as spec                                # noqa: E402
import s5_bait_screen as screen_lib                        # noqa: E402
from s3_assign import REL_MARGIN, MIN_SCORE                # noqa: E402
from s3_hmm_lib import acc_key                             # noqa: E402

CENSUS = PROJECT_ROOT / "results" / "census_v3" / "census_v3.tsv"
HMM_DIR = PROJECT_ROOT / "results" / "hmm_sweep"
OUT_DIR = PROJECT_ROOT / "results" / "s5_baits"


# ------------------------------------------------------------------ io

def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as out:
        for k, v in seqs.items():
            out.write(f">{k}\n")
            for i in range(0, len(v), 60):
                out.write(v[i:i + 60] + "\n")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_from_archive(accessions: set[str]) -> dict[str, str]:
    """Pull sequences out of S2's archived seeded-space FASTA (offline)."""
    path = require_data_root() / "raw_api" / "uniprot" / "s2_seed_sweep.fasta"
    if not path.exists():
        raise SystemExit(f"S2 archive missing: {path}\n"
                         "  rerun scripts/s2_sequences.py to rebuild it")
    out: dict[str, str] = {}
    name, buf = None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name in accessions:
                    out[name] = "".join(buf)
                parts = line[1:].split("|")
                name = acc_key(parts[1] if len(parts) > 2 else line[1:].split()[0])
                buf = []
            else:
                buf.append(line.strip())
    if name in accessions:
        out[name] = "".join(buf)
    return out


# ------------------------------------------------------- seeds already ours

def load_s3_seeds() -> tuple[dict[str, str], list[dict]]:
    """R6 — the S3 seeds, with their manifest rows, keyed by bait id."""
    seqs = read_fasta(HMM_DIR / "itpr_seed.faa")
    seqs.update(read_fasta(HMM_DIR / "ryr_seed.faa"))
    rows = read_tsv(HMM_DIR / "seed_manifest.tsv")
    missing = {r["id"] for r in rows} - set(seqs)
    if missing:
        raise SystemExit(f"S3 seed manifest names {len(missing)} sequences the "
                         f"seed FASTAs do not carry: {sorted(missing)[:3]}")
    return seqs, rows


# ---------------------------------------------------------------- selection

def build_candidates(census: list[dict]) -> tuple[dict, dict[str, int]]:
    """(paralog, band) -> ranked candidate rows, plus each paralog's reference
    length (the median of its reviewed records — measured, not typed)."""
    by_paralog_all: dict[str, list[int]] = defaultdict(list)
    for r in census:
        p = spec.paralog_of(r.get("gene", ""), r.get("protein_name", ""))
        if p and str(r.get("reviewed", "")) == "1" \
                and r.get("call") == spec.family_of(p):
            by_paralog_all[p].append(int(r.get("length") or 0))
    reference_len = {p: int(statistics.median(v))
                     for p, v in by_paralog_all.items() if v}

    cand: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in census:
        band = spec.band_of(r)
        if not band:
            continue
        p = spec.paralog_of(r.get("gene", ""), r.get("protein_name", ""))
        if not p:
            continue                                   # R1: unlabelled
        if r.get("call") != spec.family_of(p):
            continue                                   # R1: label vs call
        ok, _why = spec.passes_shape(p, int(r.get("length") or 0),
                                     r.get("n_itpr_arch", ""))
        if not ok:
            continue                                   # R2
        cand[(p, band)].append(r)
    for key, rows in cand.items():
        ref = reference_len.get(key[0], 2700)
        rows.sort(key=lambda r: spec.rank_candidate(r, ref))
    return cand, reference_len


def bait_id(row: dict, paralog: str) -> str:
    species = (row.get("species", "").split("(")[0].strip()
               .replace(" ", "_").replace("/", "_")) or "sp"
    return f"{row['accession']}|{paralog}|{species}"


def fill_slots(cand: dict, seeds_by_slot: dict, shortlist: int) -> list[dict]:
    """One shortlist per unfilled slot, best-ranked first, distinct orders."""
    wanted: list[dict] = []
    for band in spec.BANDS:
        slots = [(p, "ITPR") for p in spec.ITPR_PARALOGS] + [("RYR", "RYR")]
        for paralog, family in slots:
            need = spec.quota(band, family) - len(seeds_by_slot.get((paralog, band), []))
            if need <= 0:
                continue
            pool = (cand.get((paralog, band), []) if family == "ITPR"
                    else [r for p in spec.RYR_PARALOGS
                          for r in cand.get((p, band), [])])
            if family == "RYR":
                pool = sorted(pool, key=lambda r: spec.rank_candidate(r, 4980))
            used_orders = {s.get("order", "") for s in
                           seeds_by_slot.get((paralog, band), [])}
            taken = 0
            for row in pool:
                if taken >= need * shortlist:
                    break
                if row.get("order") and row["order"] in used_orders:
                    continue                     # R3: a distinct order
                wanted.append({"row": row, "slot_paralog": paralog,
                               "band": band, "family": family,
                               "rank": taken})
                used_orders.add(row.get("order", ""))
                taken += 1
    return wanted


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shortlist", type=int, default=3,
                    help="candidates screened per unfilled bait slot")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    census = read_tsv(CENSUS)
    census_by_acc = {acc_key(r["accession"]): r for r in census}
    seed_seqs, seed_rows = load_s3_seeds()

    # --- R6/R7: which S3 seeds enter the panel, and in which slot
    panel: dict[str, str] = {}
    manifest: list[dict] = []
    seeds_by_slot: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in seed_rows:
        clade = row["clade"]
        if clade in spec.OUT_OF_SCOPE_CLADES:
            continue                                       # R7
        c = census_by_acc.get(acc_key(row["accession"]))
        band = spec.band_of(c or {})
        if not band:
            if (c or {}).get("group") == "Vertebrata":
                raise SystemExit(
                    f"R6: S3 seed {row['id']} is a vertebrate the band rules "
                    f"cannot place (class={(c or {}).get('class', '')!r}, "
                    f"order={(c or {}).get('order', '')!r}) — add it to "
                    "s5_bait_spec.BAND_OF_ORDER rather than dropping it.")
            continue                                       # not a vertebrate
        family = "RYR" if row["profile"].upper() == "RYR" else "ITPR"
        if family == "RYR":
            slot_paralog = "RYR"                           # R5: one cell
        elif clade in spec.UNLABELLED_CLADES:
            slot_paralog = clade
        elif clade in spec.ITPR_PARALOGS:
            slot_paralog = clade
        else:
            continue
        rec = {"id": row["id"], "accession": row["accession"],
               "paralog": clade if clade not in spec.UNLABELLED_CLADES else "",
               "clade": ("RYR" if family == "RYR" else slot_paralog),
               "band": band, "family": family, "slot": f"{slot_paralog}@{band}",
               "species": row["species"], "order": (c or {}).get("order", ""),
               "length": row["length"], "arch": row["arch"],
               "reviewed": row["reviewed"], "origin": "s3_seed",
               "note": row.get("note", "")}
        seeds_by_slot[(slot_paralog, band)].append(rec)
        manifest.append(rec)
        panel[row["id"]] = seed_seqs[row["id"]]

    # RyR seeds are per band, one each (R5): keep the best-ranked per band.
    for band in spec.BANDS:
        got = seeds_by_slot.get(("RYR", band), [])
        if len(got) > spec.quota(band, "RYR"):
            got.sort(key=lambda r: (-int(r["reviewed"] or 0),
                                    abs(int(r["length"]) - 4980), r["id"]))
            for drop in got[spec.quota(band, "RYR"):]:
                manifest.remove(drop)
                panel.pop(drop["id"], None)
            seeds_by_slot[("RYR", band)] = got[:spec.quota(band, "RYR")]

    # --- R1-R3: shortlist census candidates for every slot the seeds leave open
    cand, reference_len = build_candidates(census)
    wanted = fill_slots(cand, seeds_by_slot, args.shortlist)
    need_seqs = {acc_key(w["row"]["accession"]) for w in wanted}
    got = extract_from_archive(need_seqs)
    missing = need_seqs - set(got)

    shortlist_seqs = dict(panel)
    for w in wanted:
        acc = acc_key(w["row"]["accession"])
        if acc in got:
            w["id"] = bait_id(w["row"], w["slot_paralog"])
            shortlist_seqs[w["id"]] = got[acc]

    # --- R4: screen the incumbents and the shortlist in one pass, and put the
    # screen itself through its negative control before either is believed.
    meta_by_id = {m["id"]: m for m in manifest}
    with tempfile.TemporaryDirectory(prefix="s5_baits_") as tmp:
        self_test_rows = screen_lib.self_test(panel, meta_by_id, Path(tmp))
        verdicts = screen_lib.screen(shortlist_seqs, Path(tmp))

    screen_rows: list[dict] = []
    for rec in manifest:                       # incumbents must pass too
        ok, why = screen_lib.verdict(
            verdicts.get(rec["id"], dict(screen_lib.NO_HIT)), rec["family"])
        v = verdicts.get(rec["id"], {})
        screen_rows.append({"id": rec["id"], "slot": rec["slot"],
                            "origin": "s3_seed", "rank": 0,
                            "passed": int(ok), "reason": why, **v})
        if not ok:
            raise SystemExit(
                f"R4: S3 seed {rec['id']} fails the S5 chimera screen — {why}\n"
                "  the seed set and the profiles disagree; fix S3 before S5.")

    added_by_slot: dict[tuple[str, str], int] = defaultdict(int)
    for w in sorted(wanted, key=lambda w: (w["band"], w["slot_paralog"], w["rank"])):
        if "id" not in w:
            continue
        row, key = w["row"], (w["slot_paralog"], w["band"])
        have = len(seeds_by_slot.get(key, [])) + added_by_slot[key]
        v = verdicts.get(w["id"], {})
        ok, why = screen_lib.verdict(v, w["family"])
        full = have >= spec.quota(w["band"], w["family"])
        screen_rows.append({"id": w["id"], "slot": f"{w['slot_paralog']}@{w['band']}",
                            "origin": "s5_addition", "rank": w["rank"],
                            "passed": int(ok),
                            "reason": ("slot already filled" if ok and full
                                       else why), **v})
        if not ok or full:
            panel.pop(w["id"], None)
            continue
        added_by_slot[key] += 1
        panel[w["id"]] = shortlist_seqs[w["id"]]
        manifest.append({
            "id": w["id"], "accession": row["accession"],
            "paralog": w["slot_paralog"] if w["family"] == "ITPR"
            else spec.paralog_of(row.get("gene", ""), row.get("protein_name", "")),
            "clade": w["slot_paralog"], "band": w["band"], "family": w["family"],
            "slot": f"{w['slot_paralog']}@{w['band']}",
            "species": row.get("species", ""), "order": row.get("order", ""),
            "length": row.get("length", ""), "arch": row.get("n_itpr_arch", ""),
            "reviewed": row.get("reviewed", ""), "origin": "s5_addition",
            "note": why})

    # --- outputs
    write_fasta(OUT_DIR / "baits.faa", panel)
    write_tsv(OUT_DIR / "screen_self_test.tsv",
              ["case", "rejected_by", "built_from", "description", "length",
               "rejected", "assignment", "rel_margin", "envelope_frac",
               "max_internal_gap", "n_segments", "reason"],
              self_test_rows)
    write_tsv(OUT_DIR / "bait_manifest.tsv",
              ["id", "accession", "family", "clade", "paralog", "band", "slot",
               "species", "order", "length", "arch", "reviewed", "origin",
               "note"],
              sorted(manifest, key=lambda r: (r["family"], r["clade"],
                                              r["band"], r["id"])))
    write_tsv(OUT_DIR / "bait_screen.tsv",
              ["id", "slot", "origin", "rank", "passed", "assignment",
               "rel_margin", "itpr_score", "ryr_score", "envelope_frac",
               "max_internal_gap", "n_segments", "reason"],
              screen_rows)

    filled, unfilled = {}, []
    for band in spec.BANDS:
        for paralog, family in ([(p, "ITPR") for p in spec.ITPR_PARALOGS]
                                + [("RYR", "RYR")]):
            key, want = (paralog, band), spec.quota(band, family)
            have = len(seeds_by_slot.get(key, [])) + added_by_slot[key]
            filled[f"{paralog}@{band}"] = f"{have}/{want}"
            if have < want:
                unfilled.append(f"{paralog}@{band} ({have}/{want})")

    stats = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "census": str(CENSUS.relative_to(PROJECT_ROOT)),
        "n_baits": len(panel), "n_itpr": sum(1 for m in manifest
                                             if m["family"] == "ITPR"),
        "n_ryr": sum(1 for m in manifest if m["family"] == "RYR"),
        "total_residues": sum(len(s) for s in panel.values()),
        "n_s3_seeds": sum(1 for m in manifest if m["origin"] == "s3_seed"),
        "n_s5_additions": sum(1 for m in manifest
                              if m["origin"] == "s5_addition"),
        "reference_lengths": reference_len,
        "slots": filled, "unfilled_slots": unfilled,
        "screened": len(screen_rows),
        "screen_failures": sum(1 for r in screen_rows if not r["passed"]),
        "archive_misses": sorted(missing),
        "thresholds": {"bait_band_aa": list(spec.BAIT_BAND_AA),
                       "ryr_band_aa": list(spec.RYR_BAND_AA),
                       "min_envelope_frac": spec.MIN_ENVELOPE_FRAC,
                       "max_internal_gap_aa": spec.MAX_INTERNAL_GAP_AA,
                       "rel_margin": REL_MARGIN, "min_score": MIN_SCORE},
        "screen_self_test": {r["case"]: bool(r["rejected"])
                             for r in self_test_rows},
        "sha256": {"baits.faa": sha256(OUT_DIR / "baits.faa")},
    }
    (OUT_DIR / "bait_build_stats.json").write_text(json.dumps(stats, indent=1) + "\n")

    print(f"S5 bait panel: {len(panel)} baits "
          f"({stats['n_itpr']} ITPR + {stats['n_ryr']} RyR control), "
          f"{stats['total_residues']:,} residues")
    print(f"  {stats['n_s3_seeds']} reused S3 seeds, "
          f"{stats['n_s5_additions']} S5 additions; "
          f"{stats['screen_failures']}/{len(screen_rows)} screened out")
    if unfilled:
        print(f"  unfilled slots ({len(unfilled)}): {', '.join(unfilled)}")
    if missing:
        print(f"  ! {len(missing)} shortlisted accessions absent from the "
              f"S2 archive: {sorted(missing)[:5]}")
    print(f"wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
