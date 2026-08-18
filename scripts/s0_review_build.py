#!/usr/bin/env python3
"""Assemble `docs/ip3r_review_2026.md` from its section files, and typeset it.

The review is written as numbered section files under `docs/review/` (each one
kept under the project's 500-line limit). This script is the only thing that
writes `docs/ip3r_review_2026.md`, so the assembled document cannot drift from
its sources, and the bibliography is rendered from the committed reference
table rather than maintained by hand (Decisions D13).

Citations are written in the sections as stable keys — `[R07]`, `[R22, R25]` —
and renumbered here to journal style `[7]`, `[22,25]` **in order of first
appearance**, so the printed bibliography reads 1..N down the page while the
source files stay editable without renumbering anything.

    python scripts/s0_review_build.py            # assemble the markdown
    python scripts/s0_review_build.py --pdf      # ... and typeset the PDF
    python scripts/s0_review_build.py --check    # validate only, write nothing

Exit codes: non-zero if a cited key has no reference row, if a section file is
missing from the manifest, or if the PDF stage was asked for and failed.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SECTIONS_DIR = ROOT / "docs" / "review"
REFS_TSV = ROOT / "results" / "s0_baseline" / "references.tsv"
OUT_MD = ROOT / "docs" / "ip3r_review_2026.md"
PDF_OUT = ROOT / "docs" / "ip3r_review_2026.pdf"
BUILD_MD = ROOT / "docs" / ".review_build.md"

TITLE = "The inositol 1,4,5-trisphosphate receptor family"
SUBTITLE = ("Architecture, gating, evolution and disease — a verified baseline "
            "for a genome-scale census of ITPR1/2/3")

CITE_RE = re.compile(r"\[((?:R\d{2,3})(?:\s*,\s*R\d{2,3})*)\]")


def read_refs() -> dict[str, dict]:
    lines = REFS_TSV.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    refs = {}
    for ln in lines[1:]:
        if not ln.strip():
            continue
        row = dict(zip(header, ln.split("\t")))
        refs[row["ref_id"]] = row
    return refs


def section_files() -> list[Path]:
    files = sorted(SECTIONS_DIR.glob("[0-9][0-9]_*.md"))
    if not files:
        raise SystemExit(f"no section files in {SECTIONS_DIR}")
    return files


def assemble(refs: dict[str, dict]) -> tuple[str, list[str], list[str]]:
    """Concatenate sections, renumber citations, append the bibliography."""
    body_parts = []
    for path in section_files():
        text = path.read_text(encoding="utf-8").rstrip()
        body_parts.append(text)
    body = "\n\n".join(body_parts)

    order: list[str] = []          # ref keys, in order of first appearance
    missing: list[str] = []

    def repl(m: re.Match) -> str:
        keys = [k.strip() for k in m.group(1).split(",")]
        nums = []
        for k in keys:
            if k not in refs:
                if k not in missing:
                    missing.append(k)
                nums.append("?")
                continue
            if k not in order:
                order.append(k)
            nums.append(str(order.index(k) + 1))
        return "[" + ",".join(nums) + "]"

    body = CITE_RE.sub(repl, body)

    biblio = ["", "## References", ""]
    for i, key in enumerate(order, 1):
        r = refs[key]
        authors = r["authors"]
        # journal style: first six authors then et al.
        parts = [a.strip() for a in authors.split(",") if a.strip()]
        if len(parts) > 6:
            authors = ", ".join(parts[:6]) + " *et al.*"
        else:
            authors = ", ".join(parts)
        stop = "" if authors.endswith("*") else "."
        line = (f"{i}. {authors}{stop} {r['title']}. "
                f"*{r['journal']}* **{r['year']}**.")
        if r["pmid"]:
            line += f" PMID [{r['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/)."
        if r["doi"] and r["doi"] != "-":
            line += f" doi:[{r['doi']}](https://doi.org/{r['doi']})."
        biblio.append(line)
    biblio.append("")

    uncited = [k for k in refs if k not in order]
    return body + "\n" + "\n".join(biblio), missing, uncited


def yaml_header() -> str:
    return "\n".join([
        "---",
        f'title: "{TITLE}"',
        f'subtitle: "{SUBTITLE}"',
        f'date: "{date.today().isoformat()}"',
        "documentclass: article",
        "classoption: [11pt]",
        "geometry: [a4paper, margin=2.4cm]",
        "linestretch: 1.08",
        "colorlinks: true",
        "linkcolor: linkblue",
        "urlcolor: linkblue",
        "citecolor: linkblue",
        "numbersections: false",
        "toc: true",
        "toc-depth: 2",
        "mainfont: Palatino",
        "sansfont: Helvetica Neue",
        "monofont: Menlo",
        "header-includes: |",
        "    \\usepackage{xcolor}",
        "    \\definecolor{linkblue}{RGB}{20,70,140}",
        "    \\usepackage{titlesec}",
        "    \\usepackage{fancyhdr}",
        "    \\usepackage{microtype}",
        "    \\pagestyle{fancy}",
        "    \\fancyhf{}",
        "    \\fancyhead[L]{\\footnotesize\\textsc{IP\\textsubscript{3} receptor family}}",
        "    \\fancyhead[R]{\\footnotesize\\thepage}",
        "    \\renewcommand{\\headrulewidth}{0.3pt}",
        "    \\titleformat{\\section}{\\normalfont\\Large\\bfseries\\sffamily}{\\thesection}{0.7em}{}",
        "    \\titleformat{\\subsection}{\\normalfont\\large\\bfseries\\sffamily}{\\thesubsection}{0.7em}{}",
        "    \\titleformat{\\subsubsection}{\\normalfont\\normalsize\\bfseries\\sffamily}{\\thesubsubsection}{0.7em}{}",
        "    \\setlength{\\parskip}{0.35em}",
        "    \\usepackage{longtable,booktabs,array}",
        "    \\setlength{\\emergencystretch}{3em}",
        "---",
        "",
    ])


SUB_RE = re.compile(r"<sub>(.+?)</sub>")
SUP_RE = re.compile(r"<sup>(.+?)</sup>")


def for_latex(markdown: str) -> str:
    """Adapt the GitHub-readable markdown for pandoc/LaTeX.

    The section files use HTML `<sub>`/`<sup>` so that IP3R and Ca2+ render
    correctly when the assembled markdown is read on a web front end. Pandoc
    passes raw HTML straight through and the LaTeX writer then drops it, which
    silently flattens every subscript in the typeset document — so they are
    rewritten here into pandoc's own `~x~` / `^x^` syntax.

    The leading H1 is also removed: the YAML block already supplies the title,
    and leaving both prints it twice and puts it in its own table of contents.
    """
    markdown = SUB_RE.sub(r"~\1~", markdown)
    markdown = SUP_RE.sub(r"^\1^", markdown)
    lines = markdown.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("# "):
            del lines[i]
            break
    return "\n".join(lines)


def build_pdf(markdown: str) -> int:
    if not shutil.which("pandoc"):
        print("pandoc not found — cannot typeset", file=sys.stderr)
        return 1
    engine = next((e for e in ("xelatex", "lualatex", "pdflatex")
                   if shutil.which(e)), None)
    if engine is None:
        print("no LaTeX engine found — cannot typeset", file=sys.stderr)
        return 1
    BUILD_MD.write_text(yaml_header() + for_latex(markdown),
                        encoding="utf-8")
    cmd = ["pandoc", str(BUILD_MD), "-o", str(PDF_OUT),
           f"--pdf-engine={engine}", "--from",
           "markdown+smart+pipe_tables+footnotes+tex_math_dollars"]
    print(f"  {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    if proc.returncode != 0:
        print(proc.stdout[-3000:], file=sys.stderr)
        print(proc.stderr[-4000:], file=sys.stderr)
        return proc.returncode
    size = PDF_OUT.stat().st_size / 1024
    print(f"  wrote {PDF_OUT.relative_to(ROOT)}  ({size:.0f} KB)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdf", action="store_true", help="also typeset the PDF")
    ap.add_argument("--check", action="store_true",
                    help="validate citations only; write nothing")
    ap.add_argument("--keep-build", action="store_true",
                    help="keep the intermediate pandoc markdown")
    args = ap.parse_args()

    refs = read_refs()
    files = section_files()
    print(f"sections: {len(files)}  references available: {len(refs)}")
    markdown, missing, uncited = assemble(refs)

    n_cited = len(set(CITE_RE.findall(
        "\n".join(f.read_text(encoding='utf-8') for f in files))))
    cited_keys = set()
    for f in files:
        for m in CITE_RE.finditer(f.read_text(encoding="utf-8")):
            cited_keys.update(k.strip() for k in m.group(1).split(","))
    print(f"distinct references cited: {len(cited_keys)}")

    if missing:
        print(f"ERROR: cited but not in references.tsv: {missing}",
              file=sys.stderr)
        return 2
    if uncited:
        print(f"NOTE: in references.tsv but never cited ({len(uncited)}): "
              f"{', '.join(sorted(uncited))}")

    if args.check:
        print("check only — nothing written")
        return 0

    OUT_MD.write_text(markdown, encoding="utf-8")
    print(f"  wrote {OUT_MD.relative_to(ROOT)}  "
          f"({len(markdown.splitlines())} lines)")

    rc = 0
    if args.pdf:
        rc = build_pdf(markdown)
        if not args.keep_build and BUILD_MD.exists():
            BUILD_MD.unlink()
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
