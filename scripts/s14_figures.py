"""S14 stage `figures` — collect the manuscript figure set.

Copies each main and Extended Data figure out of `results/` into
`manuscript/figures/` under its publication number, in every format the
analysis produced (PNG always; PDF where the figure script wrote one), and
writes `manuscript/figure_manifest.tsv` mapping publication number ->
source path -> checksum. Nothing is regenerated: the figures are the ones the
analyses committed, so a manuscript figure cannot silently differ from the
figure in its own results directory.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import s14_lib as lib

FORMATS = (".png", ".pdf")


def _copy_one(src_stem: str, dest_stem: str) -> list[dict]:
    """Copy every available format of one source figure. Returns manifest rows."""
    rows: list[dict] = []
    found_any = False
    for ext in FORMATS:
        src = lib.RESULTS / f"{src_stem}{ext}"
        if not src.exists():
            continue
        found_any = True
        dest = lib.MS_FIGS / f"{dest_stem}{ext}"
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        row = {
            "figure": dest_stem,
            "format": ext.lstrip("."),
            "source": str(src.relative_to(lib.ROOT)),
            "deposited_as": str(dest.relative_to(lib.ROOT)),
            "bytes": src.stat().st_size,
            "sha256": lib.sha256(src),
        }
        if ext == ".png":
            w_in, h_in = lib.png_size_inches(src)
            row["width_in"] = f"{w_in:.3f}"
            row["height_in"] = f"{h_in:.3f}"
            if w_in > lib.W_FULL + 0.02 or h_in > lib.H_MAX + 0.05:
                print(f"  OVERSIZE: {dest_stem} is {w_in:.2f} x {h_in:.2f} in "
                      f"(page allows {lib.W_FULL} x {lib.H_MAX})",
                      file=sys.stderr)
        rows.append(row)
    if not found_any:
        print(f"  MISSING: results/{src_stem}.[png|pdf]", file=sys.stderr)
    return rows


def run() -> int:
    lib.MS_FIGS.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    missing = 0

    for num, slug, src_stem, caption in lib.MAIN_FIGURES:
        got = _copy_one(src_stem, f"Fig{num}_{slug}")
        if not got:
            missing += 1
        for r in got:
            r["kind"] = "main"
            r["number"] = num
            r["caption_stub"] = caption
        rows.extend(got)

    for num, slug, src_stems, caption in lib.EXTENDED_FIGURES:
        for i, src_stem in enumerate(src_stems, start=1):
            suffix = "" if len(src_stems) == 1 else f"{chr(96 + i)}"
            got = _copy_one(src_stem, f"ExtDataFig{num}{suffix}_{slug}")
            if not got:
                missing += 1
            for r in got:
                r["kind"] = "extended_data"
                r["number"] = num
                r["caption_stub"] = caption
            rows.extend(got)

    for num, slug, src_stems, caption in lib.SUPPLEMENTARY_FIGURES:
        for i, src_stem in enumerate(src_stems, start=1):
            suffix = "" if len(src_stems) == 1 else f"{chr(96 + i)}"
            got = _copy_one(src_stem, f"SuppFig{num}{suffix}_{slug}")
            if not got:
                missing += 1
            for r in got:
                r["kind"] = "supplementary"
                r["number"] = num
                r["caption_stub"] = caption
            rows.extend(got)

    # Renumbering a figure leaves the old file behind, and a stale
    # `Fig5_structures.png` beside a live `ExtDataFig5_structures.png` is
    # exactly the kind of thing that reaches a submission. The manifest is the
    # only authority on what belongs here.
    keep = {Path(r["deposited_as"]).name for r in rows}
    for stale in sorted(lib.MS_FIGS.iterdir()):
        if stale.is_file() and stale.name not in keep:
            stale.unlink()
            print(f"  removed stale {stale.name}")

    lib.write_tsv(
        lib.MS / "figure_manifest.tsv", rows,
        ["kind", "number", "figure", "format", "width_in", "height_in",
         "source", "deposited_as", "bytes", "sha256", "caption_stub"],
    )
    n_main = len({r["figure"] for r in rows if r["kind"] == "main"})
    n_ed = len({r["figure"] for r in rows if r["kind"] == "extended_data"})
    n_supp = len({r["figure"] for r in rows if r["kind"] == "supplementary"})
    print(f"[s14 figures] {n_main} main + {n_ed} Extended Data + {n_supp} "
          f"Supplementary figure files -> manuscript/figures/ "
          f"({len(rows)} files, {missing} missing)")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(run())
