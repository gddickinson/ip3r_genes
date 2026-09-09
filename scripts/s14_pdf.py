"""S14 stage `pdf` — typeset the manuscript with its figures embedded.

Builds a pandoc-ready markdown from the numbered section files, inserting each
figure image immediately above its own legend, then renders it with
`pandoc --pdf-engine=xelatex` to `manuscript/itpr_family_manuscript.pdf`.

  python scripts/s14_pdf.py            # build the PDF
  python scripts/s14_pdf.py --keep-tex # also keep the intermediate .md/.tex

The figures used are the ones `s14_figures.py` copied under their publication
numbers, PDF preferred over PNG where the analysis wrote a vector version, so
the typeset figure is the committed figure. Requires pandoc and a LaTeX engine;
if either is missing the stage reports why and exits non-zero without touching
the rest of the package.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import date

import s14_lib as lib
import s14_refs

PDF_OUT = lib.MS / "itpr_family_manuscript.pdf"
MD_BUILD = lib.MS / ".pdf_build.md"
TEX_PREAMBLE = lib.MS / ".pdf_preamble.tex"

#: Unicode superscript characters -> their plain forms, for LaTeX math.
SUPERSCRIPTS = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
    "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+",
}

#: The same for subscripts. TeX Gyre Termes has no subscript glyphs at all, so
#: `IP₃` reached the first build as `IP` with a missing-character warning and
#: nothing on the page — the exact failure mode this module's font choice was
#: made to avoid, one plane over.
SUBSCRIPTS = {
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5",
    "₆": "6", "₇": "7", "₈": "8", "₉": "9", "₋": "-", "₊": "+",
}

#: Symbols Latin Modern has no glyph for, mapped to math mode so the PDF does
#: not silently drop them (xelatex reports these only as warnings).
MATH_MAP = {
    "ω": r"\(\omega\)", "χ": r"\(\chi\)", "α": r"\(\alpha\)",
    "β": r"\(\beta\)", "ρ": r"\(\rho\)", "Δ": r"\(\Delta\)",
    "Σ": r"\(\Sigma\)", "μ": r"\(\mu\)", "≥": r"\(\geq\)",
    "≤": r"\(\leq\)", "≈": r"\(\approx\)", "≠": r"\(\neq\)",
    "×": r"\(\times\)", "±": r"\(\pm\)", "−": r"\(-\)",
    "→": r"\(\rightarrow\)", "·": r"\(\cdot\)", "Å": r"\AA{}",
    "↔": r"\(\leftrightarrow\)", "′": r"\(^{\prime}\)",
}

#: Fonts are set by *file name*, not family name: fontconfig on this machine
#: does not resolve the TeX Gyre families, but kpathsea finds the .otf files in
#: the TeX tree. Termes (a Times clone) carries the Greek and relational glyphs
#: Latin Modern lacks, which is what silently dropped omega and >= before.
HEADER_INCLUDES = r"""
\usepackage{etoolbox}
\usepackage{float}
\usepackage{xcolor}
\definecolor{linkblue}{RGB}{26,71,133}
\setmainfont{texgyretermes}[Extension=.otf, UprightFont=*-regular,
  BoldFont=*-bold, ItalicFont=*-italic, BoldItalicFont=*-bolditalic]
\setsansfont{texgyreheros}[Extension=.otf, UprightFont=*-regular,
  BoldFont=*-bold, ItalicFont=*-italic, BoldItalicFont=*-bolditalic, Scale=0.94]
\setmonofont{texgyrecursor}[Extension=.otf, UprightFont=*-regular,
  BoldFont=*-bold, ItalicFont=*-italic, BoldItalicFont=*-bolditalic, Scale=0.86]
\setmathfont{texgyretermes-math.otf}
\AtBeginEnvironment{longtable}{\footnotesize}
\AtBeginEnvironment{tabular}{\footnotesize}
\setlength{\emergencystretch}{3em}
\renewcommand{\arraystretch}{1.15}
\usepackage{sectsty}
\allsectionsfont{\sffamily}
\setcounter{tocdepth}{3}
"""


def _figure_files() -> dict[tuple[str, int], list[tuple[str, float]]]:
    """publication figure -> [(path, drawn width in inches)], vector preferred.

    The width comes from the PNG the figure was rendered at, so each figure is
    placed at the size it was drawn rather than stretched to the margin — that
    is what keeps 7 pt in one figure the same size as 7 pt in the next.
    """
    manifest = lib.MS / "figure_manifest.tsv"
    if not manifest.exists():
        return {}
    rows = lib.read_tsv(manifest)
    width: dict[tuple[str, int, str], float] = {}
    for row in rows:
        if row["format"] == "png" and row.get("width_in"):
            width[(row["kind"], int(row["number"]), row["figure"])] = \
                float(row["width_in"])
    best: dict[tuple[str, int, str], dict[str, str]] = {}
    for row in rows:
        key = (row["kind"], int(row["number"]), row["figure"])
        seen = best.get(key)
        if seen is None or (seen["format"] == "png" and row["format"] == "pdf"):
            best[key] = row
    grouped: dict[tuple[str, int], list[tuple[str, float]]] = defaultdict(list)
    for key, row in sorted(best.items()):
        drawn = min(width.get(key, lib.W_FULL), lib.W_FULL)
        # A figure that all but fills the text block is set to the block, so
        # the page is used; only a deliberately small figure stays small. The
        # scale factor is therefore never below 1 and never above ~1.2, which
        # is what keeps type sizes comparable from figure to figure.
        placed = lib.W_FULL if drawn >= 0.8 * lib.W_FULL else drawn
        grouped[(key[0], key[1])].append((row["deposited_as"], placed))
    return grouped


def _insert_images(text: str, kind: str, pattern: re.Pattern) -> str:
    """Put each figure's image(s) directly above its legend paragraph."""
    files = _figure_files()
    out: list[str] = []
    for line in text.split("\n"):
        match = pattern.match(line)
        if match:
            number = int(match.group(1))
            paths = files.get((kind, number), [])
            if not paths:
                print(f"  no image for {kind} figure {number}", file=sys.stderr)
            for path, w_in in paths:
                rel = path.split("manuscript/", 1)[-1]
                out.append(f"![]({rel}){{width={w_in:.2f}in}}")
                out.append("")
        out.append(line)
    return "\n".join(out)


def latex_safe(text: str) -> str:
    """Make the section prose survive the LaTeX writer intact."""
    # <sup>18</sup> -> pandoc superscript
    text = re.sub(r"<sup>(.*?)</sup>", r"^\1^", text)
    # runs of unicode superscript digits -> real math superscripts
    charset = "".join(SUPERSCRIPTS)
    text = re.sub(
        f"([{charset}]+)",
        lambda m: (r"\(^{"
                   + "".join(SUPERSCRIPTS[c] for c in m.group(1)) + r"}\)"),
        text,
    )
    # runs of unicode subscript digits -> real math subscripts
    sub_charset = "".join(SUBSCRIPTS)
    text = re.sub(
        f"([{sub_charset}]+)",
        lambda m: (r"\(_{"
                   + "".join(SUBSCRIPTS[c] for c in m.group(1)) + r"}\)"),
        text,
    )
    for char, replacement in MATH_MAP.items():
        text = text.replace(char, replacement)
    return text


def _frontmatter() -> tuple[str, str, str]:
    """Return (title, author_block, remaining_body) from the front matter."""
    raw = (lib.MS / lib.SECTION_FRONTMATTER).read_text(encoding="utf-8")
    title_match = re.search(r"^# (.+)$", raw, flags=re.MULTILINE)
    title = title_match.group(1).strip() if title_match else "Manuscript"
    body = raw[title_match.end():] if title_match else raw
    author_block, _, rest = body.partition("---\n")
    return title, author_block.strip(), rest.strip()


def build_markdown() -> str:
    title, author_block, front_rest = _frontmatter()
    author = "George Dickinson"
    for line in author_block.split("\n"):
        if line.startswith("**Author:**"):
            author = re.sub(r"<sup>.*?</sup>", "",
                            line.split("**Author:**", 1)[1]).strip()

    yaml = [
        "---",
        "title: '" + latex_safe(title).replace("'", "''") + "'",
        f'author: "{author}"',
        f'date: "{date.today().isoformat()}"',
        "documentclass: article",
        "papersize: a4",
        "fontsize: 11pt",
        "geometry: margin=2.2cm",
        "linestretch: 1.12",
        "colorlinks: true",
        # Named colour, defined in the preamble: pandoc parses metadata as
        # markdown, so a "[RGB]{...}" value arrives at xcolor escaped.
        "linkcolor: linkblue",
        "urlcolor: linkblue",
        "citecolor: linkblue",
        "---",
        "",
    ]

    # The title block already prints the author; keep only the affiliation,
    # correspondence and keywords lines from the section file.
    kept = [ln for ln in author_block.split("\n")
            if not ln.startswith("**Author:**")]
    parts: list[str] = ["\n".join(yaml)]
    parts.append(latex_safe("\n".join(kept).strip()))
    parts.append("")
    parts.append(latex_safe(front_rest))
    parts.append("\n\\newpage\n\\tableofcontents\n\\newpage\n")

    main_pat = re.compile(r"^\*\*Fig\. (\d+) \|")
    ed_pat = re.compile(r"^\*\*Extended Data Fig\. (\d+) \|")
    supp_pat = re.compile(r"^\*\*Supplementary Fig\. (\d+) \|")

    for name in lib.SECTION_ORDER:
        if name == lib.SECTION_FRONTMATTER:
            continue
        path = lib.MS / name
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        if name == lib.SECTION_FIGURE_LEGENDS:
            text = _insert_images(text, "main", main_pat)
            parts.append("\n\\newpage\n")
        elif name == lib.SECTION_EXTENDED_DATA:
            text = _insert_images(text, "extended_data", ed_pat)
            parts.append("\n\\newpage\n")
        elif name == lib.SECTION_SUPPLEMENTARY:
            text = _insert_images(text, "supplementary", supp_pat)
            parts.append("\n\\newpage\n")
        elif name.startswith(("09_", "15_")):
            parts.append("\n\\newpage\n")
        parts.append(latex_safe(text))
        parts.append("")

    doc = "\n".join(parts)
    # The same renumbering the stitch stage does, so the typeset document and
    # `manuscript.md` carry identical citation numbers and one bibliography
    # rendered from `references.tsv` (D13). A cited key with no reference row
    # is reported here and leaves the keys in place rather than silently
    # dropping the citation.
    try:
        doc, bib, _ = s14_refs.resolve(doc)
        doc = doc.replace(s14_refs.MARKER, latex_safe(bib))
    except (KeyError, FileNotFoundError) as exc:
        print(f"  CITATIONS: {exc}", file=sys.stderr)
    return doc


def run(keep_intermediate: bool = False) -> int:
    if not shutil.which("pandoc"):
        print("[s14 pdf] pandoc not found — install it to build the PDF",
              file=sys.stderr)
        return 1
    engine = next((e for e in ("xelatex", "lualatex", "pdflatex")
                   if shutil.which(e)), None)
    if engine is None:
        print("[s14 pdf] no LaTeX engine found (xelatex/lualatex/pdflatex)",
              file=sys.stderr)
        return 1

    MD_BUILD.write_text(build_markdown(), encoding="utf-8")
    # The preamble goes in through --include-in-header, never through YAML:
    # pandoc parses metadata values as markdown, which turns `*-regular` into
    # emphasis and `[Extension=...]` into a span before LaTeX ever sees them.
    TEX_PREAMBLE.write_text(HEADER_INCLUDES.strip() + "\n", encoding="utf-8")
    cmd = [
        "pandoc", MD_BUILD.name,
        # tex_math_single_backslash is what makes the \\(...\\) delimiters this
        # script writes parse as math; without it pandoc reads "\\(" as an
        # escaped bracket and \\omega lands in text mode as a missing glyph.
        "--from", ("markdown+pipe_tables+superscript+subscript"
                   "+raw_tex+tex_math_single_backslash"),
        "--pdf-engine", engine,
        "--include-in-header", TEX_PREAMBLE.name,
        "--resource-path", ".",
        "-o", PDF_OUT.name,
    ]
    proc = subprocess.run(cmd, cwd=lib.MS, capture_output=True, text=True,
                          timeout=900)
    if proc.returncode != 0:
        print(proc.stdout[-4000:], file=sys.stderr)
        print(proc.stderr[-4000:], file=sys.stderr)
        print(f"[s14 pdf] pandoc failed (engine {engine})", file=sys.stderr)
        return 1
    if proc.stderr.strip():
        for line in proc.stderr.strip().split("\n")[:10]:
            print(f"  pandoc: {line}")
    if not keep_intermediate:
        MD_BUILD.unlink(missing_ok=True)
        TEX_PREAMBLE.unlink(missing_ok=True)

    size = PDF_OUT.stat().st_size
    print(f"[s14 pdf] {PDF_OUT.relative_to(lib.ROOT)} "
          f"({lib.size_str(size)}, engine {engine})")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--keep-tex", action="store_true",
                    help="keep the intermediate build markdown")
    args = ap.parse_args()
    raise SystemExit(run(keep_intermediate=args.keep_tex))
