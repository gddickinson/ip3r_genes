"""S25 stage `pdf` — typeset the thesis.

The stage reads `thesis/thesis.md`, the file the stitch stage already wrote,
rather than re-concatenating the chapters. That is the whole design: the
manuscript's PDF stage re-resolves citations and figure numbers itself and
takes care to do it the same way twice, and reading the finished document
instead removes the possibility of the two disagreeing — the typeset thesis is
the markdown thesis, converted.

What it adds is the LaTeX half: the page setup, a table of contents, each
figure swapped for its vector PDF at the width it was drawn, a page break
before each chapter, and the five traps `s14_pdf` documents, imported from
that module rather than restated:

1. the preamble goes in through `--include-in-header`, never through YAML,
   because pandoc parses metadata as markdown;
2. fonts are set by file name, since fontconfig on this machine does not
   resolve the TeX Gyre families;
3. `tex_math_single_backslash` is required or the `\\(...\\)` this module
   writes is read as an escaped bracket;
4. sub- and superscript runs are converted to math unconditionally, because
   TeX Gyre Termes has no subscript glyphs and `IP3` would otherwise print as
   `IP`;
5. pandoc runs `--verbose` so the LaTeX log is visible, and **any** dropped
   glyph fails the stage — a missing character is a warning to xelatex and a
   hole in the page to a reader.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import date

import s14_pdf
import s25_lib as L

PDF_OUT = L.TH / "itpr_family_thesis.pdf"
MD_BUILD = L.TH / ".pdf_build.md"
TEX_PREAMBLE = L.TH / ".pdf_preamble.tex"

FIG_IMG = re.compile(r"!\[\]\(figures/([A-Za-z0-9_.]+)\.png\)")


def _widths() -> dict[str, float]:
    """placed filename stem -> width in inches, from the figure manifest."""
    path = L.TH / "figure_manifest.tsv"
    if not path.exists():
        return {}
    out = {}
    for row in L.read_tsv(path):
        stem = f"fig_{row['number']}_{row['slug']}"
        try:
            drawn = min(float(row["width_in"] or L.W_FULL), L.W_FULL)
        except ValueError:
            drawn = L.W_FULL
        out[stem] = L.W_FULL if drawn >= 0.8 * L.W_FULL else drawn
    return out


def build_markdown() -> str:
    body = (L.TH / "thesis.md").read_text(encoding="utf-8")
    body = re.sub(r"<!--.*?-->", "", body, flags=re.DOTALL).lstrip()
    # The YAML block supplies the title; leaving the H1 in prints it twice and
    # puts it in its own table of contents.
    body = re.sub(r"^# .+\n", "", body, count=1)

    widths = _widths()

    def img(m: re.Match) -> str:
        stem = m.group(1)
        w = widths.get(stem, L.W_FULL)
        return ("\\begin{center}\n"
                f"\\includegraphics[width={w:.2f}in,"
                "height=0.78\\textheight,keepaspectratio]"
                f"{{figures/{stem}.pdf}}\n"
                "\\end{center}")

    body = FIG_IMG.sub(img, body)
    # A page break before each chapter heading. The horizontal rules the
    # stitch stage writes between chapters become the break, so a chapter
    # never starts halfway down a page.
    body = re.sub(r"\n---\n\n(# )", r"\n\n\\newpage\n\n\1", body)
    body = s14_pdf.latex_safe(body)

    yaml = [
        "---",
        "title: '" + s14_pdf.latex_safe(L.TITLE).replace("'", "''") + "'",
        "subtitle: '" + s14_pdf.latex_safe(L.SUBTITLE).replace("'", "''") + "'",
        'author: "' + L.AUTHOR + '"',
        f'date: "{date.today().isoformat()}"',
        "documentclass: report",
        "papersize: a4",
        "fontsize: 11pt",
        "geometry: margin=2.2cm",
        "linestretch: 1.15",
        "colorlinks: true",
        "linkcolor: linkblue",
        "urlcolor: linkblue",
        "citecolor: linkblue",
        "toc: true",
        "toc-depth: 2",
        "---",
        "",
    ]
    return "\n".join(yaml) + body


def run(keep_intermediate: bool = False) -> int:
    if not (L.TH / "thesis.md").exists():
        print("[s25 pdf] thesis/thesis.md does not exist — run the stitch "
              "stage first", file=sys.stderr)
        return 1
    if not shutil.which("pandoc"):
        print("[s25 pdf] pandoc not found", file=sys.stderr)
        return 1
    engine = next((e for e in ("xelatex", "lualatex", "pdflatex")
                   if shutil.which(e)), None)
    if engine is None:
        print("[s25 pdf] no LaTeX engine found", file=sys.stderr)
        return 1

    MD_BUILD.write_text(build_markdown(), encoding="utf-8")
    # graphicx is normally pulled in by pandoc when it lowers a markdown
    # image; every figure here is already raw `\includegraphics`, which
    # pandoc passes straight through without noticing it needs the package.
    TEX_PREAMBLE.write_text(
        s14_pdf.HEADER_INCLUDES.strip()
        + "\n\\usepackage{graphicx}\n\\usepackage{longtable,booktabs,array}\n",
        encoding="utf-8")
    cmd = ["pandoc", MD_BUILD.name,
           "--from", ("markdown+pipe_tables+superscript+subscript"
                      "+raw_tex+tex_math_single_backslash"),
           "--pdf-engine", engine,
           "--include-in-header", TEX_PREAMBLE.name,
           "--resource-path", ".", "--verbose",
           "-o", PDF_OUT.name]
    # D65 one level up: `figstyle.save()` drops the creation timestamp from a
    # figure's pdf so a figure rebuilt from the same code on the same data is
    # byte-identical. The *document* pdf had the same defect, which makes any
    # checksum recorded against it meaningless, and SOURCE_DATE_EPOCH fixes
    # it — the build date is now the epoch rather than today.
    #
    # It does not make the file byte-identical, and the residue is stated
    # rather than claimed away: this xdvipdfmx still writes a random 16-byte
    # trailer /ID, twice, so two builds of the same document differ in exactly
    # 64 bytes of 3.5 million and in nothing else. FORCE_SOURCE_DATE=1 is set
    # as well and does not change that on this build.
    env = dict(os.environ, SOURCE_DATE_EPOCH="0",
               FORCE_SOURCE_DATE="1", TZ="UTC")
    proc = subprocess.run(cmd, cwd=L.TH, capture_output=True, text=True,
                          timeout=1800, env=env)
    if proc.returncode != 0:
        print(proc.stdout[-4000:], file=sys.stderr)
        print(proc.stderr[-4000:], file=sys.stderr)
        print(f"[s25 pdf] pandoc failed (engine {engine})", file=sys.stderr)
        return 1

    log = proc.stdout + proc.stderr
    missing = sorted({ln.strip() for ln in log.split("\n")
                      if "Missing character" in ln})
    if missing:
        for ln in missing[:20]:
            print(f"  {ln}", file=sys.stderr)
        print(f"[s25 pdf] {len(missing)} distinct missing characters — the "
              f"document font lacks a glyph the text uses", file=sys.stderr)
        return 1
    if not keep_intermediate:
        MD_BUILD.unlink(missing_ok=True)
        TEX_PREAMBLE.unlink(missing_ok=True)
    import s14_lib as lib
    print(f"[s25 pdf] {PDF_OUT.relative_to(L.ROOT)} "
          f"({lib.size_str(PDF_OUT.stat().st_size)}, engine {engine})")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--keep-tex", action="store_true")
    raise SystemExit(run(keep_intermediate=ap.parse_args().keep_tex))
