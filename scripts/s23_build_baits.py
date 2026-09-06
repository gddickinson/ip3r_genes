"""s23_build_baits.py — build the S23 bait panel from census v5 by B1-B7.

The rules live in `s23_bait_spec.py` and are *enforced* here, not trusted: a
seed whose census call, architecture count or length has drifted since the
manifest was written aborts the build, exactly as S3's seed builder does.

Sequences come from committed or archived files only — S3's seed FASTAs for
the reused seeds (B6), S2's archived seeded-space FASTA for the census
additions, and the archived reference-proteome DBs for the MIR control baits
(`s23_controls.py`). Nothing here touches the network.

Outputs -> results/s23_baits/
    baits.faa                 the panel
    bait_manifest.tsv         one row per bait with its band, role and origin
    bait_rules.tsv            every rule, its threshold and what it did
    bait_candidates.tsv       every candidate considered, ranked, with verdict
    unfilled_slots.tsv        bands with no bait, and why
    control_manifest.tsv      the MIR controls and the clade each controls for
    screen_self_test.tsv      the screen's own negative controls (every build)
    bait_build_stats.json     counts, the measured band, SHA-256 of the panel

Usage:
    python3 scripts/s23_build_baits.py
    python3 scripts/s23_build_baits.py --shortlist 8
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import s23_bait_spec as spec                                   # noqa: E402
import s23_control_profiles as ctl_profiles                    # noqa: E402
import s23_control_select as controls                          # noqa: E402
import s23_controls as mir_controls                            # noqa: E402
import s5_bait_screen as screen_lib                            # noqa: E402
from s3_hmm_lib import acc_key                                 # noqa: E402
from s5_build_baits import (read_tsv, sha256,                  # noqa: E402
                            write_fasta, write_tsv)
from s23_bait_select import (DEFAULT_SHORTLIST, bait_id,       # noqa: E402
                             build_candidates, fill, load_s3_seeds, measure)

CENSUS = PROJECT_ROOT / "results" / "census_v5" / "census_v5.tsv"
HMM_DIR = PROJECT_ROOT / "results" / "hmm_sweep"
OUT_DIR = PROJECT_ROOT / "results" / "s23_baits"


# ------------------------------------------------------------------- output

RULES = [
    ("B1", "family from the census call, never a gene symbol",
     "call in (ITPR, RYR)"),
    ("B2", "RyR is a positive control only where RyR is expected; every "
           "other genome carries a control drawn from the clade the claim is "
           "about, and **which profile** that control uses is measured per "
           "clade rather than fixed in advance (s23_control_select)",
     "RYR_CONTROL_GROUPS = " + ",".join(sorted(spec.RYR_CONTROL_GROUPS))
     + "; candidates = "
     + ",".join(c["pfam"] for c in ctl_profiles.CANDIDATES)),
    ("B3", "measured length band + one bait per band, two in the deep bands",
     "DEEP_BANDS = " + ",".join(sorted(spec.DEEP_BANDS))),
    ("B4", "chimera screen: S3's margin + envelope coverage, S5's module",
     f"MIN_ENVELOPE_FRAC={spec.MIN_ENVELOPE_FRAC}, "
     f"MAX_INTERNAL_GAP_AA={spec.MAX_INTERNAL_GAP_AA}"),
    ("B5", "architecture exception for a measured incomplete architecture "
           "above the derived score floor",
     f"percentile={spec.ARCH_EXCEPTION_PERCENTILE}"),
    ("B6", "S3 seeds reused where they fall in a band", "seed_manifest.tsv"),
    ("B7", "vertebrate baits left out — S5 owns that panel", "band_of() == ''"),
]

MANIFEST_COLS = ["id", "accession", "band", "group", "role", "paralog",
                 "clade", "family", "length", "species", "taxon_id", "gene",
                 "protein_name", "n_itpr_arch", "arch_exception",
                 "itpr_score", "ryr_score", "rel_margin", "envelope_frac",
                 "origin", "shape_reason", "screen_reason"]


def manifest_row(m: dict) -> dict:
    return {
        "id": m["id"], "accession": m["accession"], "band": m["band"],
        # A control bait's group is the reference-proteome DB it came out of,
        # not the band it controls for: the cross-kingdom control tier
        # (`s23_classify.cross_group_support`) compares a bait's origin
        # against the genome's, and a class-level band name resolves to no
        # group at all.
        "group": m.get("db_group") or spec.group_of_band(m["band"]),
        "role": m["role"],
        # `paralog` and `clade` keep the S5 column names so `s5_sweep_lib`'s
        # bait loader reads this manifest unchanged; outside the vertebrates
        # the paralog is the band (B1).
        "paralog": m["band"], "clade": m["band"], "family": m["role"],
        "length": len(m["sequence"]), "species": m.get("species", ""),
        "taxon_id": m.get("taxon_id", ""), "gene": m.get("gene", ""),
        "protein_name": (m.get("protein_name") or "")[:80],
        "n_itpr_arch": m.get("n_itpr_arch", ""),
        "arch_exception": m.get("arch_exception", 0),
        "itpr_score": m.get("itpr_score", ""), "ryr_score": m.get("ryr_score", ""),
        "rel_margin": m.get("rel_margin", ""),
        "envelope_frac": "", "origin": m.get("origin", ""),
        "shape_reason": m.get("shape_reason", ""),
        "screen_reason": m.get("screen_reason", ""),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shortlist", type=int, default=DEFAULT_SHORTLIST)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--force-controls", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    census = read_tsv(CENSUS)
    in_scope = [r for r in census if spec.band_of(r)]
    itpr_rows = [r for r in in_scope if r.get("call") == spec.FAMILY_ITPR]
    meas = measure(itpr_rows)
    print(f"census v5: {len(census):,} rows, {len(in_scope):,} in scope, "
          f"{len(itpr_rows):,} ITPR")
    print(f"B3 band {meas['band_aa']} aa, B5 floor "
          f"{meas['arch_floor_bits']} bits (from "
          f"{meas['measured_on_records']} complete-architecture records)")

    slots, ref_len, rejected, in_band_counts = build_candidates(census)
    census_by_acc = {acc_key(r["accession"]): r for r in census}
    _, seeds = load_s3_seeds(census_by_acc)
    print(f"B6: {len(seeds)} S3 seeds fall in an S23 band")

    with tempfile.TemporaryDirectory(prefix="s23_baits_") as tmp:
        work = Path(tmp)
        panel, unfilled, screen_rows = fill(slots, seeds, args.shortlist,
                                            work, in_band_counts)

        ctl_rows, ctl_unfilled, ctl_seqs, ctl_stats, ctl_choice = \
            controls.build(args.threads, args.force_controls)
        controls.write_tables(ctl_stats, ctl_choice)
        # **Both controls travel, and they do different jobs.** The measured
        # one above is the proof-of-search: whichever profile a clade actually
        # carries, chosen by coverage. The MIR one below stays because it is
        # the only control that is *in the family's own signature set* — a MIR
        # protein called ITPR is the sharpest possible failure of D14, and the
        # panel is what gives that failure somewhere to show up. Dropping it
        # for the better proof-of-search would buy a control and sell a
        # negative control. Where a clade's measured choice *is* PF02815 the
        # two coincide and the duplicate is dropped.
        mir_rows, mir_unfilled, mir_seqs = mir_controls.build(
            args.threads, args.force_controls)
        have = {r["accession"] for r in ctl_rows}
        for r in mir_rows:
            if r["accession"] in have:
                continue
            ctl_rows.append(dict(r, pfam="PF02815", profile_label="MIR",
                                 reason="D14 decoy control for "
                                        f"{r['control_for']}: " + r["reason"]))
            have.add(r["accession"])
        ctl_seqs = {**mir_seqs, **ctl_seqs}
        ctl_unfilled = [u for u in ctl_unfilled
                        if u["clade"] in {m["clade"] for m in mir_unfilled}]
        ctl_meta: dict[str, dict] = {}
        for r in ctl_rows:
            bid = bait_id(r, r["control_for"], spec.CONTROL_ROLE)
            ctl_meta[bid] = dict(r, id=bid, band=r["control_for"],
                                 role=spec.CONTROL_ROLE,
                                 sequence=ctl_seqs[r["accession"]],
                                 db_group=r["group"],
                                 origin=f"{r['pfam']} over the archived "
                                        f"{r['group']} proteome DB",
                                 shape_reason=r["reason"])
        ctl_verdicts = screen_lib.screen(
            {k: v["sequence"] for k, v in ctl_meta.items()}, work)
        for bid, v in ctl_verdicts.items():
            # B2/D14: a control bait must be assigned to NEITHER family. If a
            # MIR protein is called ITPR or RYR it is not a control, it is a
            # finding — and it must not go in the panel as a control.
            assigned = v["assignment"]
            ctl_meta[bid]["screen_reason"] = (
                f"assigned {assigned} by the screen "
                f"(itpr {v['itpr_score']:.0f} / ryr {v['ryr_score']:.0f} bits)")
            ctl_meta[bid]["screen_ok"] = assigned not in ("ITPR", "RYR")
            screen_rows.append({"id": bid, "accession": ctl_meta[bid]["accession"],
                                "band": ctl_meta[bid]["band"],
                                "role": spec.CONTROL_ROLE,
                                "length": len(ctl_meta[bid]["sequence"]),
                                "passed": int(ctl_meta[bid]["screen_ok"]),
                                "reason": ctl_meta[bid]["screen_reason"], **v})
        controls_kept = [m for m in ctl_meta.values() if m["screen_ok"]]
        controls_rejected = [m for m in ctl_meta.values() if not m["screen_ok"]]

        all_rows = panel + controls_kept
        seqs = {m["id"]: m["sequence"] for m in all_rows}
        meta = {m["id"]: {"family": m["role"], "band": m["band"],
                          "clade": m["band"], "paralog": m["band"],
                          "length": len(m["sequence"])} for m in all_rows}
        self_rows = screen_lib.self_test(
            {k: v for k, v in seqs.items()
             if meta[k]["family"] in ("ITPR", "RYR")},
            meta, work, spec)

    write_fasta(OUT_DIR / "baits.faa", seqs)
    write_tsv(OUT_DIR / "bait_manifest.tsv", MANIFEST_COLS,
              [manifest_row(m) for m in all_rows])
    write_tsv(OUT_DIR / "bait_rules.tsv", ["rule", "statement", "threshold"],
              [{"rule": a, "statement": b, "threshold": c} for a, b, c in RULES])
    write_tsv(OUT_DIR / "bait_candidates.tsv",
              ["accession", "band", "role", "length", "n_itpr_arch",
               "itpr_score", "verdict", "detail"],
              [{"accession": r["accession"], "band": r["band"],
                "role": r["role"], "length": r.get("length", ""),
                "n_itpr_arch": r.get("n_itpr_arch", ""),
                "itpr_score": r.get("itpr_score", ""),
                "verdict": r["verdict"], "detail": r["detail"]}
               for r in rejected])
    write_tsv(OUT_DIR / "unfilled_slots.tsv",
              ["band", "role", "wanted", "filled", "records", "shape_ok",
               "with_sequence", "why"], unfilled)
    write_tsv(OUT_DIR / "control_manifest.tsv",
              ["id", "accession", "pfam", "profile_label", "control_for",
               "control_rank", "swept", "organism", "taxid", "length",
               "score", "model_cov", "candidates", "clade_frac", "strength",
               "strength_note", "group", "protein_name", "kept",
               "screen_reason", "reason"],
              [dict(m, kept=int(m["screen_ok"]),
                    protein_name=m["protein_name"][:80])
               for m in ctl_meta.values()])
    write_tsv(OUT_DIR / "control_unfilled.tsv",
              ["rank", "clade", "swept", "with_itpr", "why"], ctl_unfilled)
    write_tsv(OUT_DIR / "screen.tsv",
              ["id", "accession", "band", "role", "length", "passed",
               "assignment", "rel_margin", "itpr_score", "ryr_score",
               "envelope_frac", "max_internal_gap", "n_segments", "reason"],
              screen_rows)
    write_tsv(OUT_DIR / "screen_self_test.tsv",
              ["case", "rejected_by", "built_from", "description", "rejected",
               "reason", "length", "assignment", "rel_margin", "itpr_score",
               "ryr_score", "envelope_frac", "max_internal_gap"], self_rows)

    n_itpr = sum(1 for m in all_rows if m["role"] == spec.FAMILY_ITPR)
    n_ryr = sum(1 for m in all_rows if m["role"] == spec.FAMILY_RYR)
    stats = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "census": str(CENSUS.relative_to(PROJECT_ROOT)),
        "census_rows": len(census), "in_scope_rows": len(in_scope),
        **meas,
        "baits": len(all_rows), "itpr": n_itpr, "ryr": n_ryr,
        "mir_controls": len(controls_kept),
        "mir_controls_rejected": len(controls_rejected),
        "residues": sum(len(s) for s in seqs.values()),
        "reused_s3_seeds": sum(1 for m in all_rows
                               if str(m.get("origin", "")).startswith("S3")),
        "arch_exceptions": sum(1 for m in all_rows
                               if m.get("arch_exception")),
        "unfilled_slots": len(unfilled),
        "control_clades": len(ctl_meta) + len(ctl_unfilled),
        "control_clades_unfilled": len(ctl_unfilled),
        "self_test_cases": len(self_rows),
        "baits_sha256": sha256(OUT_DIR / "baits.faa"),
    }
    (OUT_DIR / "bait_build_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"\n{len(all_rows)} baits ({n_itpr} ITPR + {n_ryr} RyR + "
          f"{len(controls_kept)} proof-of-search control), "
          f"{stats['residues']:,} residues, "
          f"{stats['reused_s3_seeds']} reused S3 seeds, "
          f"{stats['arch_exceptions']} architecture exceptions")
    if unfilled:
        print(f"{len(unfilled)} unfilled slot(s): "
              + ", ".join(f"{u['band']}/{u['role']}" for u in unfilled))
    if controls_rejected:
        print(f"{len(controls_rejected)} control bait(s) rejected by the "
              "screen — a MIR protein the profiles call a family member is a "
              "finding, not a control")
    return 0


if __name__ == "__main__":
    sys.exit(main())
