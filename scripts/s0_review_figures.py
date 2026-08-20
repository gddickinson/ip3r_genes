#!/usr/bin/env python3
"""Render every figure in `docs/ip3r_review_2026.md`.

Figures are drawn only from committed tables (Decision D13 applied to
figures): `results/s0_baseline/review_figures/` plus the S0 and S1 tables that
were already committed. Nothing here queries a database or re-reads a
structure — that is the job of the three `s0_figdata_*.py` scripts, which are
run once and whose output is committed:

    python scripts/s0_figdata_domains.py       # InterPro domain coordinates
    python scripts/s0_figdata_structure.py     # cryo-EM measurements
    python scripts/s0_figdata_align.py         # MAFFT alignments

Then:

    python scripts/s0_review_figures.py                # all of them
    python scripts/s0_review_figures.py --only gating  # one, by slug
    python scripts/s0_review_figures.py --list

Every figure writes `docs/figures/<slug>.png` (read by GitHub) and
`docs/figures/<slug>.pdf` (placed by the PDF build, which swaps the extension).
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

import s0_figs_structure as _structure       # noqa: E402
import s0_figs_sequence as _sequence         # noqa: E402
import s0_figs_genomics as _genomics         # noqa: E402
import s0_figs_concepts as _concepts         # noqa: E402
import s0_figs_clinical as _clinical         # noqa: E402

#: slug -> (callable, section it illustrates). The slug is what the section
#: files reference as `figures/<slug>.png` and what `{fig:<slug>}` resolves to,
#: so renaming one here means renaming it in the review text as well.
FIGURES = {
    "discovery_timeline":   (_clinical.discovery_timeline, "§1"),
    "domain_architecture":  (_structure.domain_architecture, "§2.1, §7.1"),
    "channel_structure":    (_structure.channel_structure, "§2.2–2.4, §3"),
    "gating_logic":         (_concepts.gating_logic, "§3.2–3.3"),
    "alignment_windows":    (_sequence.alignment_windows, "§3.1, §7.1"),
    "regulation_map":       (_concepts.regulation_map, "§4"),
    "signal_hierarchy":     (_concepts.signal_hierarchy, "§5.1–5.3"),
    "conservation_profile": (_sequence.conservation_profile, "§2, §6, §9"),
    "gene_architecture":    (_genomics.gene_architecture, "§6.1, Q7"),
    "family_separation":    (_genomics.family_separation, "§7.1, Q4"),
    "taxonomic_range":      (_genomics.taxonomic_range, "§7.5, Q1"),
    "disease_map":          (_clinical.disease_map, "§9"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=None,
                    help="render just this slug (repeatable)")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for slug, (_, sec) in FIGURES.items():
            print(f"  {slug:22s} {sec}")
        return 0

    todo = args.only or list(FIGURES)
    unknown = [s for s in todo if s not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure slug(s): {unknown}")

    failed = []
    for slug in todo:
        fn, sec = FIGURES[slug]
        try:
            paths = fn()
            kb = Path(paths[0]).stat().st_size / 1024
            print(f"  ok  {slug:22s} {sec:14s} {kb:6.0f} KB")
        except Exception:                      # noqa: BLE001 - report them all
            failed.append(slug)
            print(f"  FAIL {slug}")
            traceback.print_exc()
    if failed:
        print(f"\n{len(failed)} figure(s) failed: {', '.join(failed)}",
              file=sys.stderr)
        return 1
    print(f"\n{len(todo)} figure(s) written to "
          f"{(ROOT / 'docs' / 'figures').relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
