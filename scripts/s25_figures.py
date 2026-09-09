"""S25 stage `figures` — collect the thesis's figures, and guard the set.

Nothing is re-plotted. Every figure is copied from the results directory that
committed it, in both formats it was saved in, and the manifest records the
source and the SHA-256 of each copy — so a thesis figure cannot differ from
the one beside the data it was drawn from (D13, D19; rule T7).

Four failures, each a way the set could stop being trustworthy:

* a declared figure whose png or pdf is missing at its source;
* a committed figure in a chapter-assigned results directory that the thesis
  neither places nor lists in `UNPLACED` (rule T6 applied to figures);
* a slug used twice;
* a file in `thesis/figures/` that the manifest does not name — deleted on
  every build, so a renamed figure cannot leave its predecessor behind.
"""

from __future__ import annotations

import shutil
import sys

import s14_lib as lib
import s25_figmap as fm
import s25_lib as L

FIELDS = ["chapter", "number", "slug", "source", "formats", "width_in",
          "height_in", "caption_stub", "sha256_png"]


def numbering() -> dict[str, str]:
    """slug -> 'chapter.n', numbered per chapter in declaration order."""
    seen: dict[int, int] = {}
    out = {}
    for chapter, slug, _src, _cap in fm.FIGURES:
        seen[chapter] = seen.get(chapter, 0) + 1
        out[slug] = f"{chapter}.{seen[chapter]}"
    return out


def _committed_figures() -> set[str]:
    """Every figure stem committed under a chapter-assigned results tree."""
    stems = set()
    for entry in L.ASSIGNMENT:
        d = L.RESULTS / entry / "figures"
        if d.is_dir():
            stems |= {f"results/{entry}/figures/{p.stem}"
                      for p in d.glob("*.png")}
    stems |= {f"docs/figures/{p.stem}" for p in L.DOCS_FIGS.glob("*.png")}
    return stems


def run() -> int:
    L.TH_FIGS.mkdir(parents=True, exist_ok=True)
    nums = numbering()
    rows, problems, keep = [], [], set()

    slugs = [s for _c, s, _p, _cap in fm.FIGURES]
    for slug in sorted({s for s in slugs if slugs.count(s) > 1}):
        problems.append(f"slug {slug!r} is declared more than once")

    placed = set()
    for chapter, slug, src, cap in fm.FIGURES:
        placed.add(src)
        formats = []
        for ext in ("png", "pdf"):
            source = L.ROOT / f"{src}.{ext}"
            if not source.exists():
                problems.append(f"{slug}: missing {src}.{ext}")
                continue
            dest = L.TH_FIGS / f"fig_{nums[slug]}_{slug}.{ext}"
            shutil.copy2(source, dest)
            keep.add(dest.name)
            formats.append(ext)
        png = L.ROOT / f"{src}.png"
        w = h = ""
        if png.exists():
            try:
                wf, hf = lib.png_size_inches(png)
                w, h = f"{wf:.2f}", f"{hf:.2f}"
            except ValueError as exc:                        # noqa: PERF203
                problems.append(f"{slug}: {exc}")
        rows.append({
            "chapter": str(chapter), "number": nums[slug], "slug": slug,
            "source": src, "formats": ",".join(formats),
            "width_in": w, "height_in": h, "caption_stub": cap,
            "sha256_png": lib.sha256(png) if png.exists() else "",
        })

    for stem in sorted(_committed_figures() - placed):
        if stem not in fm.UNPLACED:
            problems.append(
                f"{stem}.png is committed under a chapter's results tree but "
                f"is neither placed in the thesis nor listed in UNPLACED")

    for stale in sorted(p for p in L.TH_FIGS.iterdir()
                        if p.is_file() and p.name not in keep):
        stale.unlink()

    L.write_tsv(L.TH / "figure_manifest.tsv", rows, FIELDS)
    for p in problems:
        print(f"  [FAIL] {p}", file=sys.stderr)
    print(f"[s25 figures] {len(rows)} figures across "
          f"{len({r['chapter'] for r in rows})} chapters, "
          f"{len(keep)} files -> thesis/figures/")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(run())
