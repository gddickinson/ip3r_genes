"""S26 shared declarations: the paper series, its paths, and paper loading.

The single manuscript is one package: one section order, one figure map, one
deposit list. S26 turns it into N packages, and this module holds what those
packages have in common:

* **`SERIES`**: the papers in submission order. A paper may cite only a paper
  that comes before it (`paper_rules.md`, the citation graph), so this order
  is a claim the build checks rather than a list somebody typed once.
* **`paper(pid)`**: each paper's configuration, loaded from
  `s26_paper_<pid>.py`. One module per paper keeps every file inside the
  500-line budget and lets each paper be edited without touching the others.
* **`SECTIONS`**: the section files each paper is stitched from, in order.
  They have the manuscript's shape (front matter, introduction, results,
  discussion, methods, legends, Extended Data, references) because a paper in
  this series is the manuscript's format, not the thesis's.

Page geometry, TSV and SHA-256 helpers come from `s14_lib`, so a paper
measures a figure the same way the manuscript and the thesis do.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import s14_lib as lib

ROOT = lib.ROOT
RESULTS = lib.RESULTS
PAPERS = ROOT / "papers"
W_FULL = lib.W_FULL
H_MAX = lib.H_MAX

read_tsv = lib.read_tsv
write_tsv = lib.write_tsv
sha256 = lib.sha256

#: Submission order. The order is the citation graph read forwards: the range
#: paper defines the family call every other paper relies on, the origin
#: paper places the three paralogues, retention counts them, the archive
#: audit measures how databases record what retention found, and the two
#: channel papers come last because they stand on the orthologue sets.
SERIES = ["range", "origin", "retention", "archive", "constraint", "ligand"]

#: The author line is the manuscript's. The contributions statement is not:
#: the manuscript's says its author wrote the code, and the project's own
#: record (thesis Chapter 16, D79) says Claude did, so every paper says so.
AUTHOR = "George Dickinson"
CORRESPONDENCE = "george.dickinson@gmail.com"

#: Section files, in stitch order. `{references}` goes in the last one.
SECTIONS = [
    "00_frontmatter.md",
    "01_introduction.md",
    "02_results.md",
    "03_discussion.md",
    "04_methods.md",
    "05_figure_legends.md",
    "06_extended_data.md",
    "07_back_matter.md",
]
#: Sections whose `**{fig:...}.**` paragraphs are legends, not body text.
LEGEND_SECTIONS = {"05_figure_legends.md", "06_extended_data.md"}
#: Sections a figure reference counts as a body mention in (R3, ordering).
BODY_SECTIONS = ["01_introduction.md", "02_results.md", "03_discussion.md",
                 "04_methods.md"]

KINDS = {"main": "Fig.", "ed": "Extended Data Fig.",
         "supp": "Supplementary Fig."}
#: P4: the main-figure band.
MAIN_MIN, MAIN_MAX = 4, 7

FIG_REF = re.compile(r"\{fig:([a-z0-9_]+)\}")
PAPER_REF = re.compile(r"\{paper:([a-z]+)\}")
LEGEND_HEAD = re.compile(r"^\*\*\{fig:([a-z0-9_]+)\}\.\*\*")

_CACHE: dict[str, object] = {}


def paper(pid: str):
    """The configuration module of one paper (cached)."""
    if pid not in SERIES:
        raise KeyError(f"{pid!r} is not a paper in the series {SERIES}")
    if pid not in _CACHE:
        _CACHE[pid] = importlib.import_module(f"s26_paper_{pid}")
    return _CACHE[pid]


def number(pid: str) -> int:
    return SERIES.index(pid) + 1


def paper_dir(pid: str) -> Path:
    return PAPERS / pid


def figures(pid: str) -> list[tuple[str, str, list[str], str]]:
    """(slug, kind, source stems, caption stub) for every figure a paper
    places; a single stem may be given as a string."""
    out = []
    for slug, kind, stems, stub in paper(pid).FIGURES:
        out.append((slug, kind, [stems] if isinstance(stems, str)
                    else list(stems), stub))
    return out


def section_text(pid: str, name: str) -> str:
    path = paper_dir(pid) / name
    return path.read_text(encoding="utf-8") if path.exists() else ""


def body_text(pid: str) -> str:
    return "\n\n".join(section_text(pid, n) for n in BODY_SECTIONS)


def results_entries() -> list[str]:
    """Every entry directly under results/, as 'results/<name>'."""
    return sorted(f"results/{p.name}" for p in RESULTS.iterdir()
                  if not p.name.startswith("."))


def stem_entry(stem: str) -> str:
    """The results entry a figure stem (relative to results/) belongs to."""
    return "results/" + stem.split("/", 1)[0]
