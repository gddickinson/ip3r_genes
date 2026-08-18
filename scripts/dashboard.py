"""dashboard.py — self-contained project dashboard (dashboard.html).

Renders, from files already in the repo (no network, stdlib only):
  * an overall project progress bar (roadmap ledger: completed /
    in-progress / pending);
  * one card per ledger task with status and headline results;
  * a live panel for the in-progress task, driven by
    `results/session_live.json`, which the task's own driver script
    writes as it goes:
        {"task": "S5", "workers": 2,
         "steps": [{"label": "GCF_000001405.40", "done": true}, ...]}
    A task with a long-running sweep should write that file after every
    unit of work; anything else gets a card with no bar, which is fine.
  * the key figures (base64-embedded — the page is one portable file):
    the explicit `FIGURES` list first, then anything else found under
    `results/**/figures/*.png`, so a new figure shows up without an edit
    here;
  * FINDINGS.md rendered as the findings summary.

The page meta-refreshes every 30 s; `--watch` keeps regenerating it so
long-running sessions update live. Session protocol: run
`python scripts/dashboard.py --open` and a background
`python scripts/dashboard.py --watch` at session start.

Run:  python scripts/dashboard.py [--open] [--watch [SECONDS]] [--out PATH]
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DEFAULT = ROOT / "dashboard.html"
PROJECT = "IP3 receptor (ITPR) census & discovery"

#: Figures shown at the top of the page, in this order, when they exist.
#: Sessions add their headline figure here; everything else under
#: `results/**/figures/*.png` is picked up automatically below.
FIGURES: list[tuple[str, str]] = [
    # ("results/<dir>/figures/<name>.png", "S<n> — one-line caption"),
]

#: Auto-discovered figures are capped so the page stays a reasonable size.
MAX_AUTO_FIGURES = 12


def discovered_figures() -> list[tuple[str, str]]:
    """Every `results/**/figures/*.png` not already in FIGURES, newest first."""
    listed = {rel for rel, _ in FIGURES}
    found = []
    for path in sorted((ROOT / "results").glob("*/figures/*.png")):
        rel = str(path.relative_to(ROOT))
        if rel in listed:
            continue
        found.append((path.stat().st_mtime, rel,
                      f"{path.parent.parent.name} — {path.stem}"))
    found.sort(reverse=True)
    return [(rel, cap) for _, rel, cap in found[:MAX_AUTO_FIGURES]]


# ------------------------------------------------------------ ledger parse
def parse_ledger() -> list[dict]:
    """Both roadmap ledger tables → ordered task dicts."""
    text = (ROOT / "PUBLICATION_ROADMAP.md").read_text()
    tasks: list[dict] = []
    for line in text.splitlines():
        if not line.startswith("| S"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        # S14a / S14b / S14c are real rows, so the id pattern allows a suffix
        if len(cells) < 4 or not re.fullmatch(r"S\d+[a-z]?", cells[0]):
            continue
        if len(cells) >= 6:          # analysis table has a Priority column
            tid, title, depends, _prio, status, results = cells[:6]
        else:
            tid, title, depends, status, results = cells[:5]
        if status == "—":            # the "(promoted)" placeholder row
            continue
        s = status.lower()
        state = ("completed" if "completed" in s
                 else "in_progress" if "in_progress" in s or "in progress" in s
                 else "pending")
        date = ""
        m = re.search(r"(\d{4}-\d{2}-\d{2})", status)
        if m:
            date = m.group(1)
        clean = re.sub(r"\*\*|`", "", title)
        tasks.append({"id": tid, "title": clean, "depends": depends,
                      "state": state, "date": date,
                      "results": re.sub(r"\*\*|`", "", results)})
    return tasks


# ------------------------------------------------------------- live probes
def live_status(tasks: list[dict]) -> dict | None:
    """Progress of the in-progress task, from `results/session_live.json`.

    The file is written by whichever driver the current task runs, so the
    dashboard needs no per-task knowledge: it only checks that the file is
    talking about the task the ledger says is in progress. A stale file from
    a previous task is ignored rather than shown, because a progress bar for
    the wrong task is worse than no progress bar.
    """
    current = next((t for t in tasks if t["state"] == "in_progress"), None)
    if current is None:
        return None
    p = ROOT / "results" / "session_live.json"
    if p.exists():
        try:
            data = json.loads(p.read_text())
            if data.get("task") == current["id"]:
                data.setdefault("steps", [])
                return data
        except json.JSONDecodeError:
            pass
    return {"task": current["id"], "steps": []}


# ------------------------------------------------------- findings markdown
def md_to_html(text: str) -> str:
    out: list[str] = []
    in_list = in_quote = False

    def close():
        nonlocal in_list, in_quote
        if in_list:
            out.append("</ul>")
            in_list = False
        if in_quote:
            out.append("</blockquote>")
            in_quote = False

    def inline(s: str) -> str:
        s = html.escape(s, quote=False)
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", s)
        s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
        s = re.sub(r"\[([^]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
        return s

    for raw in text.splitlines():
        line = raw.rstrip()
        if re.fullmatch(r"-{3,}", line):
            close()
            out.append("<hr>")
        elif line.startswith("#"):
            close()
            n = len(line) - len(line.lstrip("#"))
            out.append(f"<h{n+1}>{inline(line[n:].strip())}</h{n+1}>")
        elif line.startswith("> "):
            if not in_quote:
                close()
                out.append("<blockquote>")
                in_quote = True
            out.append(inline(line[2:]) + " ")
        elif line.startswith("- "):
            if not in_list:
                close()
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{inline(line[2:])}</li>")
        elif line.startswith("  ") and (in_list or in_quote):
            out.append(inline(line.strip()) + " ")
        elif not line:
            close()
        else:
            close()
            out.append(f"<p>{inline(line)}</p>")
    close()
    return "\n".join(out)


# ------------------------------------------------------------------ pieces
def img_b64(rel: str) -> str | None:
    p = ROOT / rel
    if not p.exists():
        return None
    return base64.b64encode(p.read_bytes()).decode()


def bar(pct_done: float, pct_active: float = 0.0, label: str = "") -> str:
    seg2 = (f'<div class="seg active" style="width:{pct_active:.1f}%">'
            "</div>" if pct_active > 0 else "")
    return (f'<div class="meter" role="img" aria-label="{html.escape(label)}">'
            f'<div class="seg done" style="width:{pct_done:.1f}%"></div>'
            f"{seg2}</div>")


def task_card(t: dict, live: dict | None) -> str:
    chip = {"completed": "chip done", "in_progress": "chip active",
            "pending": "chip"}[t["state"]]
    label = {"completed": f"completed {t['date']}",
             "in_progress": "in progress", "pending": "pending"}[t["state"]]
    extra = ""
    if (t["state"] == "in_progress" and live
            and live.get("task") == t["id"] and live["steps"]):
        done = sum(1 for s in live["steps"] if s["done"])
        total = len(live["steps"])
        extra = (bar(100 * done / total,
                     label=f"{done} of {total} steps")
                 + f'<div class="muted">{done}/{total} steps done</div>')
    results = html.escape(t["results"][:220] + (
        "…" if len(t["results"]) > 220 else ""))
    dep = html.escape(t["depends"] or "—")
    return (f'<div class="card {t["state"]}">'
            f'<div class="cardhead"><span class="tid">{t["id"]}</span>'
            f'<span class="{chip}">{label}</span></div>'
            f'<div class="tasktitle">{html.escape(t["title"])}</div>'
            f'{extra}'
            f'<div class="muted">depends: {dep}</div>'
            + (f'<div class="results">{results}</div>' if results else "")
            + "</div>")


def fmt_eta(seconds: float) -> str:
    seconds = int(seconds)
    if seconds < 90:
        return f"~{seconds}s"
    if seconds < 5400:
        return f"~{round(seconds / 60)}m"
    return f"~{seconds / 3600:.1f}h"


def live_panel(live: dict | None) -> str:
    if not live or not live["steps"]:
        return ""
    done = sum(1 for s in live["steps"] if s["done"])
    total = len(live["steps"])
    units = "".join(
        f'<span class="unit {"ok" if s["done"] else "todo"}">'
        f'{"✓" if s["done"] else "·"} {html.escape(s["label"])}'
        + (f' <b>{fmt_eta(s["eta_s"])}</b>'
           if not s["done"] and s.get("eta_s") else "")
        + "</span>"
        for s in live["steps"])
    pend = [s.get("eta_s", 0) for s in live["steps"] if not s["done"]]
    workers = max(1, live.get("workers", 1))
    eta_html = ""
    if pend and all(pend):
        eta_html = (f' · rough ETA <b>{fmt_eta(sum(pend) / workers)}</b> '
                    f"({workers} parallel workers; from observed "
                    "per-residue rates)")
    return (f'<section><h2>Live — {live["task"]} session progress</h2>'
            + bar(100 * done / total, label=f"{done} of {total}")
            + f'<div class="muted">{done} of {total} pipeline steps '
            f"complete{eta_html}</div>"
            f'<div class="units">{units}</div></section>')


def figures_html() -> str:
    figs = []
    for rel, caption in FIGURES + discovered_figures():
        b64 = img_b64(rel)
        if b64:
            figs.append(
                f'<figure><img src="data:image/png;base64,{b64}" '
                f'alt="{html.escape(caption)}">'
                f"<figcaption>{html.escape(caption)}</figcaption></figure>")
    return "\n".join(figs)


CSS = """
:root { color-scheme: light dark;
  --surface:#fcfcfb; --card:#f4f3f1; --track:#e5e4e0;
  --ink:#0b0b0b; --ink2:#52514e; --line:#dddcd8;
  --blue:#2a78d6; --green:#008300; --green-bg:#e2efe2; --blue-bg:#e3edf9; }
@media (prefers-color-scheme: dark) { :root {
  --surface:#1a1a19; --card:#242423; --track:#3a3a38;
  --ink:#ffffff; --ink2:#c3c2b7; --line:#3a3a38;
  --blue:#3987e5; --green:#35a135; --green-bg:#1f3320; --blue-bg:#1d2c40; } }
* { box-sizing:border-box; }
body { margin:0; padding:24px; background:var(--surface); color:var(--ink);
  font:14px/1.5 -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width:1100px; margin:0 auto; }
h1 { font-size:20px; margin:0 0 4px; }
h2 { font-size:15px; margin:28px 0 8px; border-bottom:1px solid var(--line);
  padding-bottom:4px; }
.muted { color:var(--ink2); font-size:12px; }
.meter { height:10px; background:var(--track); border-radius:4px;
  overflow:hidden; display:flex; gap:2px; margin:6px 0 2px; }
.seg { height:100%; border-radius:4px; }
.seg.done { background:var(--green); }
.seg.active { background:var(--blue); }
.grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(310px,1fr));
  gap:10px; }
.card { background:var(--card); border:1px solid var(--line);
  border-radius:8px; padding:10px 12px; }
.card.in_progress { border-color:var(--blue); }
.cardhead { display:flex; justify-content:space-between; align-items:center; }
.tid { font-weight:700; }
.tasktitle { margin:4px 0 6px; }
.chip { font-size:11px; padding:1px 8px; border-radius:10px;
  background:var(--track); color:var(--ink2); white-space:nowrap; }
.chip.done { background:var(--green-bg); color:var(--green); }
.chip.active { background:var(--blue-bg); color:var(--blue); }
.results { font-size:12px; color:var(--ink2); margin-top:6px;
  border-top:1px dashed var(--line); padding-top:6px; }
.units { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
.unit { font-size:11px; padding:2px 8px; border-radius:10px;
  background:var(--track); color:var(--ink2); }
.unit.ok { background:var(--green-bg); color:var(--green); }
figure { margin:16px 0; }
figure img { max-width:100%; border:1px solid var(--line);
  border-radius:6px; background:#fff; }
figcaption { font-size:12px; color:var(--ink2); margin-top:4px; }
.controls { display:flex; gap:8px; align-items:center; flex-wrap:wrap;
  margin:10px 0 0; }
.controls button, .controls select { font:12px inherit;
  background:var(--card); color:var(--ink); border:1px solid var(--line);
  border-radius:6px; padding:4px 10px; cursor:pointer; }
.controls button:hover { border-color:var(--blue); }
.unit b { font-weight:600; }
#findings { background:var(--card); border:1px solid var(--line);
  border-radius:8px; padding:4px 18px 12px; }
#findings blockquote { border-left:3px solid var(--line); margin:8px 0;
  padding:2px 12px; color:var(--ink2); }
#findings code { background:var(--track); border-radius:3px;
  padding:0 4px; font-size:12px; }
"""


JS = """<script>
(function () {
  const KEY = "ip3rDashCfg";
  let cfg = { rate: 30, paused: false };
  try { Object.assign(cfg, JSON.parse(localStorage.getItem(KEY) || "{}")); }
  catch (e) {}
  let left = cfg.rate;
  const rateSel = document.getElementById("rate");
  const pauseBtn = document.getElementById("pauseBtn");
  const countdown = document.getElementById("countdown");
  rateSel.value = String(cfg.rate);
  if (![...rateSel.options].some(o => o.value === rateSel.value))
    rateSel.value = "30";
  function save() { localStorage.setItem(KEY, JSON.stringify(cfg)); }
  function paint() {
    pauseBtn.textContent = cfg.paused ? "▶ Resume updates"
                                      : "⏸ Pause updates";
    countdown.textContent = cfg.paused ? "auto-refresh paused"
      : "next refresh in " + left + " s";
  }
  document.getElementById("refreshNow").onclick =
    () => location.reload();
  pauseBtn.onclick = () => { cfg.paused = !cfg.paused;
    left = cfg.rate; save(); paint(); };
  rateSel.onchange = () => { cfg.rate = parseInt(rateSel.value, 10);
    left = cfg.rate; save(); paint(); };
  setInterval(() => {
    if (cfg.paused) return;
    left -= 1;
    if (left <= 0) location.reload(); else paint();
  }, 1000);
  paint();
})();
</script>"""


def build(out_path: Path) -> None:
    tasks = parse_ledger()
    live = live_status(tasks)
    n = max(1, len(tasks))          # never divide by zero on an empty ledger
    n_done = sum(1 for t in tasks if t["state"] == "completed")
    n_act = sum(1 for t in tasks if t["state"] == "in_progress")
    findings_md = (ROOT / "FINDINGS.md")
    findings = md_to_html(findings_md.read_text()) if findings_md.exists() \
        else "<p>No FINDINGS.md yet.</p>"
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ITPR census — project dashboard</title>
<style>{CSS}</style></head><body><main>
<h1>IP3 receptor (ITPR) census &amp; discovery — project dashboard</h1>
<div class="muted">generated {stamp} ·
sources: PUBLICATION_ROADMAP.md, FINDINGS.md, results/</div>
<div class="controls">
<button id="refreshNow" type="button">Refresh now</button>
<button id="pauseBtn" type="button"></button>
<label class="muted">every <select id="rate">
<option value="10">10 s</option><option value="30">30 s</option>
<option value="60">1 min</option><option value="300">5 min</option>
</select></label>
<span id="countdown" class="muted"></span>
</div>
<section><h2>Overall progress</h2>
{bar(100 * n_done / n, 100 * n_act * 0.5 / n,
     f"{n_done} of {n} tasks complete")}
<div class="muted">{n_done} of {n} tasks complete
(green) · {n_act} in progress (blue, half-credit) ·
{n - n_done - n_act} pending · tasks are one session each, so
≈ {n - n_done} working sessions remain</div></section>
{live_panel(live)}
<section><h2>Tasks</h2><div class="grid">
{''.join(task_card(t, live) for t in tasks)}
</div></section>
<section><h2>Key figures</h2>
{figures_html()}
</section>
<section><h2>Findings summary</h2><div id="findings">
{findings}
</div></section>
</main>""" + JS + "</body></html>"
    tmp = out_path.with_suffix(".html.tmp")
    tmp.write_text(page)
    tmp.replace(out_path)          # atomic — never a half-written page
    print(f"[dashboard] wrote {out_path} "
          f"({out_path.stat().st_size / 1e6:.1f} MB)", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--open", action="store_true",
                    help="open the page in the default browser")
    ap.add_argument("--watch", nargs="?", const=30, type=int, default=None,
                    metavar="SECONDS", help="rebuild every N seconds")
    args = ap.parse_args()

    build(args.out)
    if args.open:
        subprocess.run(["open", str(args.out)], check=False)
    if args.watch:
        while True:
            time.sleep(args.watch)
            try:
                build(args.out)
            except Exception as exc:      # keep the watcher alive
                print(f"[dashboard] rebuild failed: {exc}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
