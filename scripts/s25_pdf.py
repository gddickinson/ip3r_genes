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
that module rather than restated. S28 added the legend's own typography
(**R4**): the stitch stage already knows which paragraphs are legends, since
each opens `**Figure N.M.**`, so every one is wrapped in `\figlegendbegin`
... `\figlegendend`, a command pair the preamble defines as a small sans
paragraph at unit leading, indented on both sides and set ragged-right (a
justified legend on that narrow a measure stretches its word spacing), and
the figure and its legend are held together inside one float so a page break
cannot come between them. The float is `[!htbp]` rather than `[H]`: pinned in
place, a figure too tall for the rest of its page left half the page empty,
six times in the first S28 read, and `\clearpage` before every chapter keeps a
floated figure inside its own chapter. Commands rather than an environment, because pandoc passes an
unknown command through and parses the markdown after it, whereas everything
inside a raw environment would reach LaTeX unconverted. `check_legends()`
fails the stage if any legend paragraph is left unwrapped.

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
LEGEND_PARA = re.compile(r"^\*\*Figure \d+\.\d+\.\*\*.*", re.S)

#: The legend's typography, and the float that binds it to its figure.
LEGEND_PREAMBLE = r"""
\usepackage{setspace}
\usepackage{xurl}
% A paragraph that cannot be set within the margin (a DOI is one unbreakable
% word) may loosen its spacing rather than run into the margin.
\setlength{\emergencystretch}{3em}
\definecolor{legendink}{RGB}{40,40,36}
\newcommand{\figlegendfont}{\sffamily\small\color{legendink}}
\newcommand{\figblockbegin}{\begin{figure}[!htbp]\centering}
\newcommand{\figblockend}{\end{figure}}
\newcommand{\figlegendbegin}{\par\begingroup\figlegendfont
  \setstretch{1.0}\setlength{\parindent}{0pt}\setlength{\leftskip}{1.4em}
  \setlength{\rightskip}{1.4em plus 2.5em}\vspace{0.3em}}
\newcommand{\figlegendend}{\par\endgroup}
"""


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
        return ("\\figblockbegin\n\n"
                f"\\includegraphics[width={w:.2f}in,"
                "height=0.72\\textheight,keepaspectratio]"
                f"{{figures/{stem}.pdf}}")

    body = FIG_IMG.sub(img, body)
    # A DOI is one word to TeX, so a reference line ending in one cannot
    # break and runs into the margin; set as a URL (xurl) it may break at a
    # slash or a dot. PDF only: thesis.md keeps the plain form.
    body = DOI_RE.sub(lambda m: f"doi:\\url{{{m.group(1)}}}", body)
    body = wrap_legends(body)
    # A page break before each chapter heading. The horizontal rules the
    # stitch stage writes between chapters become the break, so a chapter
    # never starts halfway down a page.
    body = re.sub(r"\n---\n\n(# )", r"\n\n\\clearpage\n\n\1", body)
    # The stitch stage also writes a rule between the two files of a long
    # chapter and after the last one; on the page those are stray lines in
    # the middle of a chapter, so every remaining separator is dropped.
    body = re.sub(r"\n---\n(?=\n|$)", "\n", body)
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
        "geometry: margin=2.0cm",   # s14_lib.W_FULL: figures are drawn for a 17.0 cm block
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


#: An overfull line narrower than this is invisible on the page; wider, it
#: is text in the margin. LaTeX's own report is in points.
OVERFULL_PT = 4.0
DOI_RE = re.compile(r"\bdoi:(10\.[0-9]{4,}/[^\s]*[^\s.,;])")
OVERFULL_RE = re.compile(r"Overfull \\hbox \(([0-9.]+)pt too wide\)")


def _overfull(log: str) -> list[str]:
    """Every overfull horizontal box wider than the bar, with the line."""
    out = []
    lines = log.split("\n")
    for i, ln in enumerate(lines):
        m = OVERFULL_RE.search(ln)
        if m and float(m.group(1)) > OVERFULL_PT:
            context = lines[i + 1].strip() if i + 1 < len(lines) else ""
            out.append(f"{ln.strip()} | {context[:90]}")
    return out


def check_log(log: str) -> list[str]:
    """The typeset-page defects only the LaTeX log records.

    Shared with the paper series (S26), so every document this project
    typesets is held to the same page: no line wider than the margin by more
    than `OVERFULL_PT`, no figure block taller than a page, and no character
    the document font cannot set.
    """
    out = []
    overfull = _overfull(log)
    if overfull:
        out += [f"  {ln}" for ln in overfull[:20]]
        out.append(f"{len(overfull)} overfull line(s) wider than "
                   f"{OVERFULL_PT} pt run into the margin")
    too_tall = [ln for ln in log.split("\n") if "Float too large" in ln]
    if too_tall:
        out.append(f"{len(too_tall)} figure block(s) taller than a page")
    missing = sorted({ln.strip() for ln in log.split("\n")
                      if "Missing character" in ln})
    if missing:
        out += [f"  {ln}" for ln in missing[:20]]
        out.append(f"{len(missing)} distinct missing characters: the "
                   f"document font lacks a glyph the text uses")
    return out


def wrap_legends(body: str, pattern: re.Pattern = LEGEND_PARA) -> str:
    """Wrap every legend paragraph in the command pair, closing the float.

    `pattern` recognises a legend; the paper series (S26) passes its own,
    since a paper's legends open `**Fig. N.**` rather than `**Figure N.M.**`.
    """
    paras = body.split("\n\n")
    # The markers sit on lines of their own: pandoc reads `\cmd **` as the
    # starred form of the command and leaves the bold marker in the output.
    return "\n\n".join(
        f"\\figlegendbegin\n\n{p}\n\n\\figlegendend\n\n\\figblockend"
        if pattern.match(p.strip()) else p for p in paras)


def check_legends(body: str, pattern: re.Pattern = LEGEND_PARA) -> list[str]:
    """R4: every legend paragraph is wrapped, and every float is closed."""
    problems = []
    paras = [p.strip() for p in body.split("\n\n")]
    for i, p in enumerate(paras):
        if pattern.match(p) and (i == 0 or paras[i - 1]
                                     != "\\figlegendbegin"):
            problems.append(f"legend not wrapped in its own typography: "
                            f"'{p[:60]}'")
    n_open = body.count("\\figblockbegin")
    n_close = body.count("\\figblockend")
    if n_open != n_close:
        problems.append(f"{n_open} figure blocks opened, {n_close} closed")
    return problems


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

    md = build_markdown()
    for p in check_legends(md):
        print(f"  [FAIL] {p}", file=sys.stderr)
    if check_legends(md):
        print("[s25 pdf] a legend would be set as body text (R4)",
              file=sys.stderr)
        return 1
    MD_BUILD.write_text(md, encoding="utf-8")
    # graphicx is normally pulled in by pandoc when it lowers a markdown
    # image; every figure here is already raw `\includegraphics`, which
    # pandoc passes straight through without noticing it needs the package.
    TEX_PREAMBLE.write_text(
        s14_pdf.HEADER_INCLUDES.strip()
        + "\n\\usepackage{graphicx}\n\\usepackage{longtable,booktabs,array}\n"
        + LEGEND_PREAMBLE,
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

    problems = check_log(proc.stdout + proc.stderr)
    if problems:
        for ln in problems:
            print(f"[s25 pdf] {ln}", file=sys.stderr)
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
