"""S2 step 5 — check the architecture call with sequences, not annotation.

`s2_call.py` separates ITPR from RYR on domain architecture. That is one
kind of evidence, and it inherits whatever InterPro's annotation got wrong.
This module tests it against a second, independent kind: the labelled-bait
identity margin S1 built for D14, run over the committed per-phylum core
panel plus the six human baits.

Nothing here consults a gene symbol or a Pfam list. Each panel member is
aligned with the baits and assigned to whichever family it is closer to;
`margin` is (identity to the nearest ITPR bait) − (identity to the nearest
RyR bait), so a positive margin is a sequence-level ITPR call. Agreement
with the architecture call is the result; disagreement names the records
S3's profile sweep has to settle.

Both identity metrics are reported, because S1 found they disagree on
fragments: the full-alignment metric divides by every column, so a genuine
member covering half the alignment scores low, while `covered_only=True`
scores over mutually covered columns. A record whose call flips between
them is doing so for a stated reason.

Agreement is scored separately inside and outside D7's no-call band. A
margin of 0.0005 is not a disagreement — it is the test declining to
separate two families at that distance, and pooling those rows with the
decisive ones would understate the check in one direction and overstate it
in the other.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis.alignment import progressive_msa          # noqa: E402
from src.analysis.distance import identity_matrix           # noqa: E402
from src.discovery.candidates import DiscoveryConfig        # noqa: E402
from src.utils.family import REFERENCE_PANEL, SISTER_PANEL  # noqa: E402
from scripts.s1_lib import MafftTracer                      # noqa: E402
from src.utils.data_root import require_data_root           # noqa: E402
from scripts.s2_lib import OUT_DIR, ROOT, read_tsv, write_tsv  # noqa: E402
from scripts.s2_sequences import load_fasta                 # noqa: E402

#: The six human baits are pulled from the seeded-space sweep rather than
#: from S1's panel files, because the RyR half of the panel lives in
#: `panel_decoys.fasta` and pairing the two would let a file-naming slip
#: decide which family a bait is labelled with.
BAIT_SOURCE = ("raw_api", "uniprot", "s2_seed_sweep.fasta")
COLS = ["accession", "species", "group", "phylum", "arch_call", "confidence",
        "seq_call", "agrees", "decisive", "margin", "best_itpr_bait",
        "id_to_itpr", "best_ryr_bait", "id_to_ryr", "metric"]

#: D7's no-call band, read off the discovery config's own default rather
#: than retyped, so the sequence check and the scorer cannot drift apart.
NO_CALL = DiscoveryConfig.__dataclass_fields__["sister_margin"].default


def bait_pairs(seqs: dict[str, str]) -> list[tuple[str, str]]:
    """The six human references as (label, sequence).

    The `BAIT_ITPR|` / `BAIT_RYR|` prefixes are the *only* family labels in
    this analysis, and they are attached to the human references, never to
    a panel member.
    """
    out = []
    for family, panel in (("ITPR", REFERENCE_PANEL), ("RYR", SISTER_PANEL)):
        for label, acc in panel:
            if acc in seqs:
                out.append((f"BAIT_{family}|{label}", seqs[acc]))
    return out


def margins(dist, covered: bool) -> dict[str, dict]:
    """Labelled-bait margin per panel member (D14, S1's `bait_margins`)."""
    labels = dist.labels
    idx = {l: i for i, l in enumerate(labels)}
    itpr = [l for l in labels if l.startswith("BAIT_ITPR|")]
    ryr = [l for l in labels if l.startswith("BAIT_RYR|")]
    out = {}
    for lbl in labels:
        if lbl.startswith("BAIT_"):
            continue

        def best(baits):
            hits = [(dist.identity[idx[lbl]][idx[b]], b) for b in baits]
            return max(hits) if hits else (float("nan"), "")

        i_id, i_b = best(itpr)
        r_id, r_b = best(ryr)
        out[lbl] = {"best_itpr_bait": i_b, "id_to_itpr": round(i_id, 4),
                    "best_ryr_bait": r_b, "id_to_ryr": round(r_id, 4),
                    "margin": round(i_id - r_id, 4),
                    "seq_call": "ITPR" if i_id > r_id else "RYR",
                    "metric": "covered_only" if covered else "full_alignment"}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", default=str(OUT_DIR / "census_v2_core_panel.faa"))
    args = ap.parse_args()

    panel = list(load_fasta(Path(args.panel)).items())
    bait_fasta = require_data_root(ROOT).joinpath(*BAIT_SOURCE)
    baits = bait_pairs(load_fasta(bait_fasta))
    if len(baits) != 6:
        raise SystemExit(f"expected 6 labelled baits, found {len(baits)} in "
                         f"{bait_fasta}")
    pairs = panel + baits
    print(f"[s2] aligning {len(pairs)} sequences "
          f"({len(panel)} panel + {len(baits)} baits) …", flush=True)

    t0 = time.time()
    with MafftTracer() as tracer:
        msa = progressive_msa(pairs, use_mafft=True)
    width = len(msa[0].aligned) if msa else 0
    if not tracer.calls or tracer.calls[0]["returncode"] != 0:
        raise SystemExit("MAFFT did not run cleanly — a star-alignment "
                         "fallback would make every margin meaningless")
    print(f"[s2] MAFFT rc={tracer.calls[0]['returncode']} "
          f"{len(msa)}×{width} cols in {time.time() - t0:.0f}s")

    census = {r["accession"]: r for r in read_tsv(OUT_DIR / "census_v2.tsv")}
    rows, disagree = [], 0
    for covered in (False, True):
        dist = identity_matrix(msa, covered_only=covered)
        for lbl, m in margins(dist, covered).items():
            c = census.get(lbl, {})
            agrees = int(m["seq_call"] == c.get("call"))
            if not covered and not agrees:
                disagree += 1
            rows.append({
                "accession": lbl, "species": c.get("species", ""),
                "group": c.get("group", ""), "phylum": c.get("phylum", ""),
                "arch_call": c.get("call", ""),
                "confidence": c.get("confidence", ""),
                "agrees": agrees,
                "decisive": int(abs(m["margin"]) >= NO_CALL), **m})
    write_tsv(OUT_DIR / "sequence_check.tsv",
              sorted(rows, key=lambda r: (r["metric"], r["margin"])), COLS)
    (OUT_DIR / "sequence_check_meta.json").write_text(json.dumps({
        "panel": Path(args.panel).name, "n_panel": len(panel),
        "n_baits": len(baits), "mafft_calls": tracer.calls,
        "alignment_columns": width, "seconds": round(time.time() - t0),
        "no_call_band": NO_CALL,
        "agreement": {
            m: {"decisive_n": sum(1 for r in rows
                                  if r["metric"] == m and r["decisive"]),
                "decisive_agree": sum(r["agrees"] for r in rows
                                      if r["metric"] == m and r["decisive"]),
                "band_n": sum(1 for r in rows
                              if r["metric"] == m and not r["decisive"]),
                "band_agree": sum(r["agrees"] for r in rows
                                  if r["metric"] == m and not r["decisive"])}
            for m in ("full_alignment", "covered_only")},
    }, indent=1))

    for metric in ("full_alignment", "covered_only"):
        sub = [r for r in rows if r["metric"] == metric]
        dec = [r for r in sub if r["decisive"]]
        band = [r for r in sub if not r["decisive"]]
        print(f"  {metric:<16s} decisive |margin| >= {NO_CALL:.2f}: "
              f"{sum(r['agrees'] for r in dec)}/{len(dec)} agree   "
              f"no-call band: {sum(r['agrees'] for r in band)}/{len(band)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
