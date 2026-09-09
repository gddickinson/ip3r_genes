"""S25 stage `production` — measure how this project was actually made.

Chapter 16 describes how the work was carried out. Like every other chapter,
its load-bearing numbers have to come from somewhere checkable, so they are
measured from the repository itself rather than recalled: the session log, the
task ledger, the decisions log, the scripts, the committed tables and figures,
the rendered reports, and the git history.

Everything here is derived from files under version control, so the table
rebuilds offline and a number in Chapter 16 cannot drift from the repository
it describes. `s25_claims_thesis.py` declares claims against it exactly as the
analysis chapters declare claims against their own tables.

  python scripts/s25_production.py       # write thesis/production_stats.tsv
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import s25_lib as L

FIELDS = ["metric", "value", "unit", "what_it_counts", "source"]


def _lines(paths) -> int:
    return sum(len(p.read_text(encoding="utf-8", errors="replace").splitlines())
               for p in paths)


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=L.ROOT, capture_output=True,
                              text=True, timeout=30).stdout.strip()
    except Exception:                                        # noqa: BLE001
        return ""


def build() -> list[dict]:
    log = (L.ROOT / "SESSION_LOG.md").read_text(encoding="utf-8")
    road = (L.ROOT / "PUBLICATION_ROADMAP.md").read_text(encoding="utf-8")
    scripts = sorted((L.ROOT / "scripts").glob("*.py"))
    src = sorted((L.ROOT / "src").rglob("*.py"))
    reports = sorted(L.RESULTS.glob("*/report.md"))

    sessions = re.findall(r"^## (\d{4}-\d{2}-\d{2})", log, re.M)
    decisions = sorted(set(re.findall(r"^\*\*(D\d+[a-z]?) ", road, re.M)))
    ledger = re.findall(r"^\| (S\d+[a-z]?) \|", road, re.M)
    # The ledger has two shapes: the main block has three columns before the
    # status and the analysis block has four (it carries a priority), so the
    # status is found by name rather than by position.
    completed = [ln for ln in road.splitlines()
                 if re.match(r"^\| S\d+[a-z]? \|", ln) and "| completed " in ln]
    # Results tables only. The thesis's own generated tables live under
    # thesis/ and counting them would make this metric count itself, since
    # this file is one of them.
    tables = _git("ls-files", "results/*.tsv").splitlines()
    figures = _git("ls-files", "results/*/figures/*.png",
                   "docs/figures/*.png").splitlines()
    commits = _git("rev-list", "--count", "HEAD") or "0"
    coauthored = len([b for b in _git("log", "--format=%b").split("\n")
                      if "Co-Authored-By: Claude" in b])

    report_words = sum(len(p.read_text(encoding="utf-8").split())
                       for p in reports)

    rows = [
        ("sessions_logged", len(sessions), "sessions",
         "dated entries in the session log, one per working session",
         "SESSION_LOG.md"),
        ("session_days", len(set(sessions)), "days",
         "calendar days on which at least one session ran",
         "SESSION_LOG.md"),
        ("first_session", sorted(sessions)[0] if sessions else "-", "date",
         "the first logged session", "SESSION_LOG.md"),
        ("last_session", sorted(sessions)[-1] if sessions else "-", "date",
         "the most recent logged session", "SESSION_LOG.md"),
        ("ledger_tasks", len(ledger), "tasks",
         "rows in the task ledger", "PUBLICATION_ROADMAP.md"),
        ("ledger_completed", len(completed), "tasks",
         "ledger rows marked completed", "PUBLICATION_ROADMAP.md"),
        ("decisions_recorded", len(decisions), "decisions",
         "distinct numbered entries in the decisions log",
         "PUBLICATION_ROADMAP.md"),
        ("analysis_scripts", len(scripts), "files",
         "Python files under scripts/, all version-controlled", "scripts/"),
        ("analysis_script_lines", _lines(scripts), "lines",
         "total lines across those files", "scripts/"),
        ("application_lines", _lines(src), "lines",
         "total lines of the search application under src/", "src/"),
        ("committed_result_tables", len(tables), "files",
         "version-controlled TSV tables under results/, which excludes the "
         "thesis's own generated tables", "git ls-files"),
        ("committed_figures", len(figures), "files",
         "version-controlled figure PNGs in results and review directories",
         "git ls-files"),
        ("task_reports", len(reports), "files",
         "rendered per-task reports, none hand-written", "results/*/report.md"),
        ("task_report_words", report_words, "words",
         "total words across those reports", "results/*/report.md"),
        ("commits", int(commits), "commits",
         "commits on the main branch", "git rev-list"),
        ("commits_coauthored", coauthored, "commits",
         "commits carrying the Claude co-author trailer", "git log"),
    ]
    return [dict(zip(FIELDS, (m, str(v), u, w, s))) for m, v, u, w, s in rows]


def run() -> int:
    rows = build()
    L.write_tsv(L.TH / "production_stats.tsv", rows, FIELDS)
    n = len(rows)
    print(f"[s25 production] {n} measurements of how the project was made "
          f"-> thesis/production_stats.tsv")
    if n < 10:
        print("  [FAIL] the measurement set is incomplete", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
