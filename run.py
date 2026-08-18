"""Entry point for the Protein Variant Finder.

GUI usage:
    python run.py
    python run.py --email you@example.com
    python run.py --preset ip3r --auto-search

Headless usage (no GUI; useful for CI / batch runs):
    python run.py --headless --preset ip3r --save-results
    python run.py --headless --preset ip3r --species "Homo sapiens" --max 25 --save-results
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description="Protein Variant Finder")
    p.add_argument("--email", default="",
                   help="Email for NCBI Entrez identification (Entrez asks for one).")
    p.add_argument("--preset", default="",
                   help="Preset name (e.g. 'ip3r', 'ip3r_zebrafish'); see presets/.")
    p.add_argument("--auto-search", action="store_true",
                   help="GUI mode: run preset's search on startup.")
    p.add_argument("--headless", action="store_true",
                   help="Run without a GUI — requires --preset.")
    p.add_argument("--save-results", action="store_true",
                   help="Headless mode: write a results bundle to results/<timestamp>_<label>/.")
    p.add_argument("--species", default="",
                   help="Override preset species (e.g. 'Homo sapiens').")
    p.add_argument("--max", type=int, default=0,
                   help="Override preset max_per_source.")
    p.add_argument("--label", default="",
                   help="Label for the results folder name (defaults to gene list).")
    p.add_argument("--analyze", action="store_true",
                   help="Headless mode: run sequence analysis (alignment / tree / clusters / mutations).")
    p.add_argument("--discover", action="store_true",
                   help="Headless mode: rank candidate novel paralogs (e.g. a fourth ITPR).")
    p.add_argument("--interpro", action="store_true",
                   help="Discovery: fetch InterPro Pfam signatures for each UniProt accession (slower).")
    p.add_argument("--known-paralog", action="append", default=[],
                   help="Repeatable. Discovery: gene symbol(s) that count as 'known' (defaults to preset genes).")
    p.add_argument("--investigate", default="",
                   help="Headless mode: deep-dive investigation of one UniProt accession (e.g. A0A0P1B5Q5).")
    p.add_argument("--vs", default="",
                   help="Comma-separated UniProt accessions for the phylogenetic comparison panel (default: the family reference panel in src/utils/family.py — human ITPR1/2/3; add a RYR accession to test which side of the superfamily split a candidate falls on).")
    p.add_argument("--foldseek", action="store_true",
                   help="Investigation: also run Foldseek structural search (slow, async).")
    args = p.parse_args()

    project_root = Path(__file__).resolve().parent
    if args.investigate:
        from src.cli import run_investigate
        run_investigate(
            project_root=project_root,
            accession=args.investigate,
            email=args.email,
            panel_csv=args.vs,
            run_foldseek=args.foldseek,
            save_results=args.save_results,
            label_override=args.label,
        )
        return
    if args.headless:
        if not args.preset:
            p.error("--headless requires --preset")
        from src.cli import run_headless
        run_headless(
            project_root=project_root,
            preset=args.preset,
            email=args.email,
            species_override=args.species,
            max_override=args.max,
            save_results=args.save_results,
            label_override=args.label,
            analyze=args.analyze,
            discover=args.discover,
            known_paralogs=args.known_paralog or None,
            interpro_lookup=args.interpro,
        )
    else:
        from src.gui import run as run_gui
        run_gui(
            project_root=project_root,
            email=args.email,
            preset=args.preset,
            auto_search=args.auto_search,
        )


if __name__ == "__main__":
    main()
