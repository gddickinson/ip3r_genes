"""S26 stage `pdf`: typeset one paper from its stitched `paper.md`.

Typesets the stitched file as it stands, the thesis's rule, so the markdown
and the PDF cannot disagree about a figure number or a citation. Each
figure's committed vector PDF is placed immediately above its legend, at the
width it was drawn, inside one float with the legend (a page break cannot
separate them). The legend is set in the thesis's legend typography (R4).

The page checks are the thesis's, not the manuscript's, as the 2026-09-23
emergent row asked. `s25_pdf.check_log` is called unchanged: no overfull line
wider than 4 pt, no float taller than a page, no dropped glyph. Text
conversion is `s14_pdf.latex_safe`, so sub- and superscripts and Greek are
set the same way in every document.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from datetime import date

import s14_pdf
import s25_pdf
import s26_lib as L

LEGEND_PARA = re.compile(
    r"^\*\*(?:Fig\.|Extended Data Fig\.|Supplementary Fig\.) \d+\.\*\*.*",
    re.S)
KIND_OF = {"Fig.": "main", "Extended Data Fig.": "ed",
           "Supplementary Fig.": "supp"}
HEAD = re.compile(r"^\*\*(Fig\.|Extended Data Fig\.|Supplementary Fig\.) "
                  r"(\d+)\.\*\*")


def pdf_path(pid: str):
    return L.paper_dir(pid) / f"paper{L.number(pid)}_{pid}.pdf"


def _placements(pid: str) -> dict[tuple[str, int], list[tuple[str, float]]]:
    """(kind, number) -> [(figure pdf path, placed width in inches)]."""
    rows = L.read_tsv(L.paper_dir(pid) / "figure_manifest.tsv")
    widths = {r["figure"]: float(r["width_in"]) for r in rows
              if r["format"] == "png" and r.get("width_in")}
    out: dict[tuple[str, int], list[tuple[str, float]]] = {}
    for r in rows:
        if r["format"] != "pdf":
            continue
        drawn = min(widths.get(r["figure"], L.W_FULL), L.W_FULL)
        placed = L.W_FULL if drawn >= 0.8 * L.W_FULL else drawn
        out.setdefault((r["kind"], int(r["number"])), []).append(
            (r["deposited_as"], placed))
    for v in out.values():
        v.sort()
    return out


def build_markdown(pid: str) -> tuple[str, list[str]]:
    text = (L.paper_dir(pid) / "paper.md").read_text(encoding="utf-8")
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S).lstrip()
    m = re.search(r"^# (.+)$", text, re.M)
    title = m.group(1).strip() if m else L.paper(pid).TITLE
    text = text[m.end():] if m else text
    text = re.sub(r"^\*\*Author:\*\*.*\n", "", text, flags=re.M)
    fails: list[str] = []
    place = _placements(pid)
    paras = text.split("\n\n")
    out = []
    for p in paras:
        h = HEAD.match(p.strip())
        if h:
            key = (KIND_OF[h.group(1)], int(h.group(2)))
            files = place.get(key, [])
            if not files:
                fails.append(f"{pid}: no figure file for {h.group(1)} "
                             f"{h.group(2)}")
            # One height budget per figure, shared by its panels: a
            # three-file figure given 0.66 of the page per file is a float
            # taller than the page (found by the range and archive papers).
            share = 0.66 if len(files) < 2 else round(0.68 / len(files), 3)
            imgs = "\n\n".join(
                f"\\includegraphics[width={w:.2f}in,height={share}"
                f"\\textheight,keepaspectratio]{{{path}}}"
                for path, w in files)
            out.append("\\figblockbegin\n\n" + imgs)
        out.append(p)
    body = "\n\n".join(out)
    body = s25_pdf.DOI_RE.sub(lambda m: f"doi:\\url{{{m.group(1)}}}", body)
    body = s25_pdf.wrap_legends(body, LEGEND_PARA)
    fails += [f"{pid}: {x}" for x in s25_pdf.check_legends(body, LEGEND_PARA)]
    # A page break before the legends, the Extended Data and the back matter,
    # so a figure section never starts halfway down the discussion.
    body = re.sub(r"\n(## (?:Figure legends|Extended Data|Data availability))",
                  r"\n\\clearpage\n\n\1", body)
    body = re.sub(r"\n---\n(?=\n|$)", "\n", body)
    body = s14_pdf.latex_safe(body)
    yaml = [
        "---",
        "title: '" + s14_pdf.latex_safe(title).replace("'", "''") + "'",
        f'author: "{L.AUTHOR}"',
        f'date: "Paper {L.number(pid)} of {len(L.SERIES)}, '
        f'{date.today().isoformat()}"',
        "documentclass: article", "papersize: a4", "fontsize: 11pt",
        "geometry: margin=2.0cm", "linestretch: 1.12", "colorlinks: true",
        "linkcolor: linkblue", "urlcolor: linkblue", "citecolor: linkblue",
        "---", "",
    ]
    return "\n".join(yaml) + body, fails


def page_count(path) -> int | str:
    """Pages, from pdfinfo; '?' when it is not installed (the PDF's own
    object streams are compressed, so counting page objects in the bytes
    returns zero)."""
    if not shutil.which("pdfinfo"):
        return "?"
    out = subprocess.run(["pdfinfo", str(path)], capture_output=True,
                         text=True).stdout
    m = re.search(r"^Pages:\s+(\d+)", out, re.M)
    return int(m.group(1)) if m else "?"


def run(pid: str, keep: bool = False) -> int:
    if not (L.paper_dir(pid) / "paper.md").exists():
        print(f"[s26 pdf {pid}] paper.md missing; run stitch first",
              file=sys.stderr)
        return 1
    if not shutil.which("pandoc"):
        print(f"[s26 pdf {pid}] pandoc not found", file=sys.stderr)
        return 1
    engine = next((e for e in ("xelatex", "lualatex") if shutil.which(e)),
                  None)
    if engine is None:
        print(f"[s26 pdf {pid}] no LaTeX engine found", file=sys.stderr)
        return 1
    md, fails = build_markdown(pid)
    if fails:
        for f in fails:
            print(f"  [FAIL] {f}", file=sys.stderr)
        return 1
    d = L.paper_dir(pid)
    (d / ".pdf_build.md").write_text(md, encoding="utf-8")
    (d / ".pdf_preamble.tex").write_text(
        s14_pdf.HEADER_INCLUDES.strip()
        + "\n\\usepackage{graphicx}\n\\usepackage{longtable,booktabs,array}\n"
        + s25_pdf.LEGEND_PREAMBLE, encoding="utf-8")
    out = pdf_path(pid)
    cmd = ["pandoc", ".pdf_build.md",
           "--from", ("markdown+pipe_tables+superscript+subscript"
                      "+raw_tex+tex_math_single_backslash"),
           "--pdf-engine", engine, "--include-in-header", ".pdf_preamble.tex",
           "--resource-path", ".", "--verbose", "-o", out.name]
    env = dict(os.environ, SOURCE_DATE_EPOCH="0", FORCE_SOURCE_DATE="1",
               TZ="UTC")
    proc = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                          timeout=1800, env=env)
    if proc.returncode != 0:
        print(proc.stderr[-3000:], file=sys.stderr)
        print(f"[s26 pdf {pid}] pandoc failed", file=sys.stderr)
        return 1
    problems = s25_pdf.check_log(proc.stdout + proc.stderr)
    for ln in problems:
        print(f"[s26 pdf {pid}] {ln}", file=sys.stderr)
    if not keep:
        (d / ".pdf_build.md").unlink(missing_ok=True)
        (d / ".pdf_preamble.tex").unlink(missing_ok=True)
    pages = page_count(out)
    print(f"[s26 pdf {pid}] {out.relative_to(L.ROOT)} ({pages} pages, "
          f"engine {engine}): {len(problems)} page problem(s)")
    return 1 if problems else 0
