"""s23_control_profiles.py — choose each clade's positive control by measuring it.

**The problem this exists to fix.** S23a replaced S5's RyR positive control
with the MIR-domain sharer (PF02815), because outside Metazoa a silent RyR is
the correct answer and witnesses nothing (D26). That was the right move for
land plants and Dikarya, where PF02815 is everywhere. It was the wrong move
for Apicomplexa, and the pilot said so: both apicomplexan classes carry **one**
PF02815 protein each across 60 swept proteomes, and *Toxoplasma gondii* came
back `uncontrolled` — no receptor and no control — so "Apicomplexa 0/60" could
not go to assembly level at all.

The fault is not PF02815. It is **fixing one profile in advance for every
clade**. A control's job is to fire in the clade whose absence is the claim, so
which protein family makes the best control is a property of the clade, and it
is measurable: run a panel of candidate profiles over that clade's own swept
reference proteomes and see which one is actually there.

**The candidates.** Every candidate is a large, deeply conserved,
multi-exon eukaryotic protein family — deliberately, because the control has
to exercise the same instrument the receptor search does. A control on a
single-exon 150 aa protein would show that miniprot can align *something* and
say nothing about whether it can recover a 2,700 aa, ~58-exon gene. PF02815
stays in the panel as a candidate and wins wherever it deserves to; it is no
longer assumed.

**What is *not* chosen by measurement**, and must not be: the family call.
A control bait is still screened (`s5_bait_screen`) and must be assigned to
*neither* ITPR nor RyR. A control protein that the profiles call a family
member is a finding, not a control.

Outputs -> results/s23_baits/control_profile_coverage.tsv   every candidate x clade
           results/s23_baits/control_profile_choice.tsv     the winner per clade

Usage:
  python3 scripts/s23_control_profiles.py --fetch          # profiles only
  python3 scripts/s23_control_profiles.py --threads 6      # fetch + search + choose
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import requests                                            # noqa: E402

from src.utils.data_root import require_data_root          # noqa: E402

PFAM_DIR = PROJECT_ROOT / "results" / "s23_scope" / "pfam"
S20_PFAM_DIR = PROJECT_ROOT / "results" / "s20_sweep" / "pfam"
OUT_DIR = PROJECT_ROOT / "results" / "s23_baits"
HOSTS = ("https://www.ebi.ac.uk/interpro/api",
         "https://www.ebi.ac.uk/interpro/wwwapi")

#: The candidate control profiles, each with the reason it is a candidate.
#: The list is short and every entry is a *large multi-domain eukaryotic*
#: family, so whichever one a clade selects, the control exercises spliced
#: alignment over a long multi-exon gene rather than a single small domain.
#: None of them is the family: PF02815 is the only one the receptor shares,
#: and it is exactly the one the screen has to reject as a family member.
CANDIDATES = (
    {"pfam": "PF02815", "label": "MIR",
     "why": "the incumbent — one of the family's own four signatures, carried "
            "by the protein O-mannosyltransferases; the sharpest available "
            "decoy as well as a control"},
    {"pfam": "PF00063", "label": "Myosin_head",
     "why": "the myosin motor domain — a large multi-exon motor present in "
            "every eukaryotic lineage including the apicomplexan glideosome"},
    {"pfam": "PF00122", "label": "E1-E2_ATPase",
     "why": "the P-type ATPase core — ~1,000 aa pumps (SERCA, PMCA, the "
            "aminophospholipid and heavy-metal pumps) universal in eukaryotes"},
    {"pfam": "PF02463", "label": "SMC_N",
     "why": "the SMC condensin/cohesin ATPase — ~1,200 aa, single-copy per "
            "SMC gene and conserved from bacteria to metazoa"},
    {"pfam": "PF00225", "label": "Kinesin",
     "why": "the kinesin motor domain — universal in eukaryotes and, unlike "
            "dynein, retained in the angiosperms"},
    {"pfam": "PF00004", "label": "AAA", "fallback": True,
     "why": "the AAA+ ATPase module — the broadest eukaryotic control "
            "available, present in every proteome, so a clade that fails "
            "even this one has an assembly problem. Marked `fallback`: it is "
            "searched only if the candidates above leave a control clade "
            "with nothing, which is what being a fallback means. A profile "
            "consulted only when the others fail is not the same instrument "
            "as one ranked against them, and searching it regardless would "
            "put it in the ranking table as if it were"},
)

def primary() -> tuple:
    return tuple(c for c in CANDIDATES if not c.get("fallback"))


def fallbacks() -> tuple:
    return tuple(c for c in CANDIDATES if c.get("fallback"))


#: E-value for the candidate searches. Same as `s23_controls.CONTROL_E`: this
#: is bait selection, not the relaxed negative sweep.
CONTROL_E = 1e-5


def profile_path(pfam: str) -> Path:
    """Where a candidate's HMM lives, preferring S20's committed copy."""
    shared = S20_PFAM_DIR / f"{pfam}.hmm"
    if shared.exists() and shared.stat().st_size > 0:
        return shared
    return PFAM_DIR / f"{pfam}.hmm"


def fetch_profile(pfam: str) -> dict:
    """Download one candidate HMM from InterPro (two-host failover, as S20)."""
    dest = profile_path(pfam)
    if dest.exists() and dest.stat().st_size > 0:
        return {"pfam": pfam, "status": "cached", "path": str(dest),
                "bytes": dest.stat().st_size}
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = ""
    for attempt in range(1, 4):
        for host in HOSTS:
            try:
                r = requests.get(f"{host}/entry/pfam/{pfam}?annotation=hmm",
                                 timeout=180)
            except requests.RequestException as exc:       # noqa: PERF203
                last = f"{type(exc).__name__} @{host}"
                continue
            if r.status_code != 200:
                last = f"HTTP {r.status_code} @{host}"
                continue
            body = (gzip.decompress(r.content)
                    if r.content[:2] == b"\x1f\x8b" else r.content)
            text = body.decode()
            if not text.startswith("HMMER3"):
                last = f"not an HMM @{host}"
                continue
            dest.write_text(text)
            name = next((ln.split()[1] for ln in text.splitlines()
                         if ln.startswith("NAME")), "")
            return {"pfam": pfam, "status": "fetched", "path": str(dest),
                    "bytes": dest.stat().st_size, "name": name}
        time.sleep(2 * attempt)
    raise SystemExit(f"could not fetch {pfam}: {last}")


def profile_name(pfam: str) -> str:
    path = profile_path(pfam)
    if not path.exists():
        return ""
    with path.open() as fh:
        for line in fh:
            if line.startswith("NAME"):
                return line.split()[1]
            if line.startswith("HMM "):
                break
    return ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true",
                    help="download the candidate profiles and stop")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    PFAM_DIR.mkdir(parents=True, exist_ok=True)
    for c in CANDIDATES:
        st = fetch_profile(c["pfam"])
        print(f"  {c['pfam']:9s} {profile_name(c['pfam']):16s} "
              f"{st['status']:8s} {st['bytes']:>9,} B")
    if args.fetch:
        return 0

    import s23_control_select as sel
    sel.build_and_write(threads=args.threads, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
