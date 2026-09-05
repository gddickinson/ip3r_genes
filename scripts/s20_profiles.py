"""S20 step 5 — the relaxed panel's profiles.

The negative claims ("no ITPR in Dikarya / land plants / archaea /
bacteria") are the sharpest thing this task says, and a negative made with
one instrument at one threshold is worth very little. Two widenings, both
stated:

**Sensitivity.** Every search runs at E ≤ 10 and the primary call is taken
by filtering the same domtblout at E ≤ 1e-5. `-E` is hmmsearch's *reporting*
threshold — the acceleration filters are unchanged by it — so the strict set
is exactly the set the strict run would have produced, and the relaxed set
is a superset from the same search rather than a second experiment.

**Profiles.** `itpr.hmm` is a 2,684-state full-length channel model. A
bacterial or plant protein carrying only the IP₃-binding core would score
poorly against it at any E-value simply for being 200 aa against 2,684
states. So the relaxed panel adds the family's four **Pfam domain** models,
downloaded from InterPro and committed here:

    PF08709  Ins145_P3_rec   the IP₃-binding core that names the family
    PF02815  MIR             shared with POMT1/2 (S1's decoy panel)
    PF01365  RYDR_ITPR       shared with the ryanodine receptors
    PF08454  RIH_assoc       shared with the ryanodine receptors

Three of the four are shared with something, which is the point: for a
*negative* claim a promiscuous profile is the right instrument. Anything it
finds is then handed to the ITPR-vs-RYR call (D14) rather than counted.

Run:  python3 scripts/s20_profiles.py            # fetch + verify
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from scripts.s20_lib import S20_DIR, HMM_SWEEP_DIR, log, write_json  # noqa: E402
from src.utils.family import FAMILY_PFAM_IDS  # noqa: E402

PFAM_DIR = S20_DIR / "pfam"
HOSTS = ("https://www.ebi.ac.uk/interpro/api",
         "https://www.ebi.ac.uk/interpro/wwwapi")

#: Why each Pfam is in the relaxed panel. Keyed off `family.FAMILY_PFAM_IDS`
#: so the panel cannot drift from the family definition (CLAUDE.md's rule).
PFAM_NOTES = {
    "PF08709": "IP3-binding core — the signature that names the family",
    "PF02815": "MIR — shared with POMT1/2",
    "PF01365": "RYDR_ITPR — shared with the ryanodine receptors",
    "PF08454": "RIH_assoc — shared with the ryanodine receptors",
}


def fetch_pfam_hmm(pfam: str, dest: Path) -> dict:
    """Download one Pfam HMM, gunzipping if InterPro serves it compressed."""
    if dest.exists() and dest.stat().st_size > 0:
        return {"pfam": pfam, "path": str(dest), "status": "cached",
                "bytes": dest.stat().st_size}
    last = ""
    for attempt in range(1, 4):
        for host in HOSTS:
            try:
                r = requests.get(f"{host}/entry/pfam/{pfam}?annotation=hmm",
                                 timeout=180)
                if r.status_code != 200:
                    last = f"HTTP {r.status_code} @{host}"
                    continue
                body = (gzip.decompress(r.content) if r.content[:2] == b"\x1f\x8b"
                        else r.content)
                text = body.decode()
                if not text.startswith("HMMER3"):
                    last = f"not an HMM @{host}"
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_text(text)
                return {"pfam": pfam, "path": str(dest), "status": "fetched",
                        "bytes": dest.stat().st_size}
            except requests.RequestException as exc:
                last = f"{type(exc).__name__} @{host}"
        time.sleep(3 * attempt)
    raise RuntimeError(f"could not fetch {pfam} HMM: {last}")


def hmm_header(path: Path) -> dict:
    """NAME / ACC / LENG out of an HMM file — proof of what was downloaded."""
    out: dict[str, str] = {}
    for line in path.read_text().splitlines()[:20]:
        for key in ("NAME", "ACC", "DESC", "LENG"):
            if line.startswith(key):
                out[key.lower()] = line.split(None, 1)[1].strip()
    return out


def relaxed_panel() -> list[tuple[str, Path]]:
    """(label, hmm path) for every profile in the relaxed panel.

    The two full-length family profiles first — the relaxed claim is made
    with the same instrument as the primary one, widened — then the four
    Pfam domain models.
    """
    panel = [("itpr", HMM_SWEEP_DIR / "itpr.hmm"),
             ("ryr", HMM_SWEEP_DIR / "ryr.hmm")]
    panel += [(pf, PFAM_DIR / f"{pf}.hmm") for pf in FAMILY_PFAM_IDS]
    return panel


def build() -> dict:
    rows = []
    for pf in FAMILY_PFAM_IDS:
        info = fetch_pfam_hmm(pf, PFAM_DIR / f"{pf}.hmm")
        info.update(hmm_header(Path(info["path"])))
        info["why"] = PFAM_NOTES.get(pf, "")
        rows.append(info)
        log("s20_profiles", f"{pf} {info['status']:8s} "
                            f"{info.get('name', '?')} LENG={info.get('leng', '?')}")
    stats = {"panel": [lbl for lbl, _ in relaxed_panel()], "pfam": rows}
    write_json(S20_DIR / "relaxed_panel.json", stats)
    return stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    build()
    for label, path in relaxed_panel():
        if not path.exists():
            log("s20_profiles", f"MISSING {label}: {path}")
            return 1
    log("s20_profiles", f"relaxed panel ready: "
                        f"{', '.join(l for l, _ in relaxed_panel())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
