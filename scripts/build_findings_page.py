"""Render `docs/findings_summary.md` into one self-contained HTML page.

The markdown is the single source of truth: this script only converts it, so
the page and the repository document cannot disagree. Every figure the markdown
links to is downscaled and inlined as a WebP data URI, so the result is one
file with no external requests — publishable as an Artifact, openable from a
USB stick, and safe to send to a collaborator who has none of the data.

The paralog summary card at the top is *not* in the markdown: it is read live
from the committed tables (`genome_ledger.tsv`, `loss_events.tsv`,
`omega_table.tsv`), so it cannot drift either.

  python scripts/build_findings_page.py                # -> docs/findings_summary.html
  python scripts/build_findings_page.py -o /tmp/x.html
  python scripts/build_findings_page.py --max-width 1600
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs" / "findings_summary.md"
CSS = Path(__file__).resolve().parent / "findings_page.css"
DEFAULT_OUT = ROOT / "docs" / "findings_summary.html"

#: The page's own name and deck. The markdown's H1 is a document heading
#: ("What this project found ..."); a published page needs a name instead, and
#: the scope belongs in the kicker where a reader can size the work at a glance.
PAGE_TITLE = "The ITPR Census"
#: Scope line. Filled in with the real denominators once the census tasks
#: have run — leave a placeholder rather than a wrong number.
KICKER = ("vertebrate genomes &middot; invertebrate genomes "
          "&middot; eukaryotic reference proteomes")
STANDFIRST = ("What a sequence-level count of the IP3-receptor family found: "
              "from the channel's eukaryotic origin to the state of its "
              "database records.")

#: (display name, css swatch variable, one-line blurb). The blurbs are
#: written once the fates are measured; until then they say what the paralog
#: is, not what happened to it.
PARALOGS = [
    ("ITPR1", "p1", "The cerebellar paralog: Purkinje-cell calcium release, "
                    "and the one with a clinical variant set (SCA15/SCA29, "
                    "Gillespie syndrome)."),
    ("ITPR2", "p2", "The broadly expressed paralog; recessive loss-of-function "
                    "causes isolated anhidrosis in humans."),
    ("ITPR3", "p3", "The epithelial/secretory paralog, and the one carrying "
                    "the dominant neuropathy alleles."),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------ figures

def encode_figure(path: Path, max_width: int) -> tuple[str, int, int]:
    """Downscale to `max_width` and return (data URI, width, height)."""
    from PIL import Image

    with Image.open(path) as im:
        im = im.convert("RGB")
        if im.width > max_width:
            height = round(im.height * max_width / im.width)
            im = im.resize((max_width, height), Image.LANCZOS)
        buf = io.BytesIO()
        im.save(buf, format="WEBP", quality=88, method=6)
        payload = base64.b64encode(buf.getvalue()).decode("ascii")
        return f"data:image/webp;base64,{payload}", im.width, im.height


# ----------------------------------------------------------------- markdown

INLINE = [
    (re.compile(r"`([^`]+)`"), lambda m: f"<code>{html.escape(m.group(1))}</code>"),
    (re.compile(r"\*\*(.+?)\*\*"), lambda m: f"<strong>{m.group(1)}</strong>"),
    (re.compile(r"(?<!\*)\*([^*]+)\*(?!\*)"), lambda m: f"<em>{m.group(1)}</em>"),
]


def inline(text: str) -> str:
    """Escape, then apply the inline marks. Code spans escape their own body."""
    out = html.escape(text, quote=False)
    for pattern, repl in INLINE:
        out = pattern.sub(repl, out)
    return out


def _blocks(lines: list[str]):
    """Yield (kind, payload) blocks from the markdown body."""
    buf: list[str] = []
    kind = None

    def flush():
        nonlocal buf, kind
        if buf:
            yield_val = (kind, list(buf))
            buf, kind = [], None
            return yield_val
        buf, kind = [], None
        return None

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            got = flush()
            if got:
                yield got
        elif stripped.startswith("## "):
            got = flush()
            if got:
                yield got
            yield ("h2", [stripped[3:]])
        elif stripped == "---":
            got = flush()
            if got:
                yield got
            yield ("hr", [])
        elif stripped.startswith("!["):
            got = flush()
            if got:
                yield got
            yield ("img", [stripped])
        elif stripped.startswith("|"):
            if kind != "table":
                got = flush()
                if got:
                    yield got
                kind = "table"
            buf.append(stripped)
        elif stripped.startswith("- "):
            if kind != "ul":
                got = flush()
                if got:
                    yield got
                kind = "ul"
            buf.append(stripped[2:])
        elif kind == "ul":
            buf[-1] += " " + stripped          # continuation of a list item
        else:
            kind = kind or "p"
            buf.append(stripped)
        i += 1
    got = flush()
    if got:
        yield got


IMG_RE = re.compile(r"!\[(.*?)\]\((.*?)\)")
NUM_RE = re.compile(r"^(\d+)\.\s+(.*)$")


def render_body(md: str, max_width: int) -> tuple[str, int]:
    lines = md.split("\n")
    parts: list[str] = []
    n_figs = 0
    open_section = False
    pending_eyebrow = None

    for kind, payload in _blocks(lines):
        if kind == "h2":
            if open_section:
                parts.append("</section>")
            open_section = True
            match = NUM_RE.match(payload[0])
            if match:
                num, title = match.groups()
                parts.append(
                    '<section class="finding"><div class="finding-head">'
                    f'<div class="finding-num">{num}</div>'
                    f"<h2>{inline(title)}</h2>"
                    '<div class="col2" id="eyebrow-slot"></div></div>')
            else:
                parts.append('<section class="finding">'
                             '<div class="finding-head">'
                             '<div class="finding-num"></div>'
                             f"<h2>{inline(payload[0])}</h2></div>")
        elif kind == "hr":
            if open_section:
                parts.append("</section>")
                open_section = False
            parts.append("<hr>")
        elif kind == "img":
            match = IMG_RE.match(payload[0])
            if not match:
                continue
            alt, rel = match.groups()
            src = (SRC.parent / rel).resolve()
            if not src.exists():
                print(f"  missing figure: {rel}", file=sys.stderr)
                continue
            uri, w, h = encode_figure(src, max_width)
            n_figs += 1
            try:
                caption = str(src.relative_to(ROOT))
            except ValueError:
                caption = src.name
            parts.append(
                f'<figure><img src="{uri}" width="{w}" height="{h}" '
                f'alt="{html.escape(alt, quote=True)}" loading="lazy">'
                f"<figcaption>{html.escape(caption)}</figcaption></figure>")
        elif kind == "table":
            rows = [r for r in payload if not set(r.replace("|", "").strip())
                    <= set("-: ")]
            cells = [[c.strip() for c in r.strip("|").split("|")] for r in rows]
            head, body = cells[0], cells[1:]
            thead = "".join(f"<th>{inline(c)}</th>" for c in head)
            tbody = "".join(
                "<tr>" + "".join(f"<td>{inline(c)}</td>" for c in row) + "</tr>"
                for row in body)
            parts.append('<div class="table-scroll"><table><thead><tr>'
                         f"{thead}</tr></thead><tbody>{tbody}</tbody>"
                         "</table></div>")
        elif kind == "ul":
            items = "".join(f"<li>{inline(it)}</li>" for it in payload)
            parts.append(f"<ul>{items}</ul>")
        else:
            text = " ".join(payload)
            if text.startswith("*Evidence:") and text.endswith("*"):
                pending_eyebrow = text[1:-1]
                slot = f'<p class="eyebrow">{inline(pending_eyebrow)}</p>'
                for i in range(len(parts) - 1, -1, -1):
                    if '<div class="col2" id="eyebrow-slot"></div>' in parts[i]:
                        parts[i] = parts[i].replace(
                            '<div class="col2" id="eyebrow-slot"></div>',
                            f'<div class="col2">{slot}</div>')
                        break
                continue
            parts.append(f"<p>{inline(text)}</p>")

    if open_section:
        parts.append("</section>")
    # any section that never supplied an eyebrow keeps an empty slot
    body_html = "\n".join(parts).replace(
        '<div class="col2" id="eyebrow-slot"></div>', "")
    return body_html, n_figs


# --------------------------------------------------------------- trio card

def trio_card() -> str:
    """The three-paralog summary card, read live from the committed tables.

    Returns "" until the census exists — a card of zeroes reads as a result
    and is worse than no card. The genome denominator comes from the
    manifest, so it cannot disagree with the ledger it summarises.
    """
    ledger_p = ROOT / "results" / "genome_ledger.tsv"
    if not ledger_p.exists():
        return ""
    ledger = read_tsv(ledger_p)
    manifest_p = ROOT / "results" / "genome_manifest.tsv"
    n_genomes = (len({r["accession"] for r in read_tsv(manifest_p)})
                 if manifest_p.exists()
                 else len({r.get("accession", "") for r in ledger}))

    def maybe(path: Path) -> list[dict[str, str]]:
        return read_tsv(path) if path.exists() else []

    losses = maybe(ROOT / "results" / "loss_dynamics" / "loss_events.tsv")
    omegas = maybe(ROOT / "results" / "selection" / "omega_table.tsv")
    found, loss_n, omega = {}, {}, {}
    for row in ledger:
        if row["status"].startswith("found"):
            found[row["paralog"]] = found.get(row["paralog"], 0) + 1
    for row in losses:
        loss_n[row["paralog"]] = loss_n.get(row["paralog"], 0) + 1
    for row in omegas:
        if row.get("model") == "M0" and row["set"] in {p for p, _, _ in PARALOGS}:
            omega[row["set"]] = f'{float(row["omega"]):.3f}'

    cards = []
    for name, key, blurb in PARALOGS:
        cards.append(
            f'<div><div class="name">'
            f'<span class="swatch" style="background:var(--{key})"></span>'
            f"{name}</div><dl>"
            f"<dt>found in</dt><dd>{found.get(name, 0)} of "
            f"{n_genomes} genomes</dd>"
            + (f"<dt>losses</dt><dd>{loss_n[name]}</dd>"
               if name in loss_n else "")
            + (f'<dt>&omega;</dt><dd>{omega[name]}</dd>'
               if name in omega else "")
            + f'</dl><p class="lede">{html.escape(blurb)}</p></div>')
    return '<div class="trio">' + "".join(cards) + "</div>"


# ------------------------------------------------------------------- driver

def build(out: Path, max_width: int) -> int:
    md = SRC.read_text(encoding="utf-8")
    title_match = re.search(r"^# (.+)$", md, flags=re.M)
    title = title_match.group(1).strip() if title_match else "Findings"
    md = md[title_match.end():] if title_match else md

    # the markdown's opening italic paragraph is provenance, not a deck: it
    # moves to the colophon so the page opens on the finding, not the build date
    stand_match = re.search(r"^\*(.+?)\*$", md, flags=re.M | re.S)
    provenance = ""
    if stand_match and stand_match.start() < md.find("\n## "):
        provenance = inline(" ".join(stand_match.group(1).split()))
        md = md[:stand_match.start()] + md[stand_match.end():]

    body, n_figs = render_body(md, max_width)
    head_html = (
        '<div class="wrap"><header class="masthead">'
        f'<p class="kicker">{KICKER}</p>'
        f"<h1>{html.escape(PAGE_TITLE)}</h1>"
        f'<p class="standfirst">{STANDFIRST}</p></header>'
        + trio_card())
    colophon = (
        '<div class="colophon">Rendered from <code>docs/findings_summary.md</code> '
        "by <code>scripts/build_findings_page.py</code>. Figures are the "
        "committed analysis output, downscaled and inlined; nothing on this "
        "page is re-plotted or re-computed. "
        "Full manuscript: <code>manuscript/itpr_family_manuscript.pdf</code>."
        f"<br>{provenance}</div>")

    page = (f"<title>{html.escape(PAGE_TITLE)}</title>\n<style>\n"
            f"{CSS.read_text(encoding='utf-8')}\n</style>\n"
            f"{head_html}\n{body}\n{colophon}\n</div>\n")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    size = out.stat().st_size
    print(f"[findings] {out} ({size / 1e6:.1f} MB, {n_figs} figures inlined "
          f"at <= {max_width}px)")
    return 0 if size < 16e6 else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--max-width", type=int, default=1500)
    args = ap.parse_args()
    sys.exit(build(args.out, args.max_width))
