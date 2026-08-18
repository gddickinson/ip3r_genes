"""S1 step 1 — probe the local toolchain and record exact versions.

Every downstream task shells out to these binaries; a manuscript methods
section has to name the version that produced each result. This script is
the only place those versions are collected, and it writes two artefacts:

    results/toolchain_manifest.txt          — human-readable, committed
    results/benchmark_controls/toolchain.tsv — machine-readable, feeds the
                                               S1 report (D13)

Binaries are looked up on PATH first, then in the recorded conda env
(`ENV_BIN`), so a tool that only exists inside the env is still found and
is recorded *as* env-resident rather than silently passing as a PATH tool.

Run:  python3 scripts/s1_toolchain.py
"""

from __future__ import annotations

import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_TXT = PROJECT_ROOT / "results" / "toolchain_manifest.txt"
OUT_TSV = PROJECT_ROOT / "results" / "benchmark_controls" / "toolchain.tsv"

#: The conda env this project runs in. BLAST+, the NCBI `datasets` CLI and
#: Foldseek live here rather than in Homebrew; biopython/matplotlib too.
#: Reused from the PIEZO project rather than cloned — see Decisions D18.
ENV_NAME = "piezo1"
ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")

#: (binary, argv for its version banner, regex capturing the version, task
#:  that first needs it). `None` regex = keep the first non-empty line.
TOOLS: list[tuple[str, list[str], str | None, str]] = [
    ("mafft",       ["--version"],  r"(v[\d.]+)",                   "S1/S6"),
    ("hmmbuild",    ["-h"],         r"HMMER ([\d.]+)",              "S3"),
    ("hmmsearch",   ["-h"],         r"HMMER ([\d.]+)",              "S3"),
    ("jackhmmer",   ["-h"],         r"HMMER ([\d.]+)",              "S3"),
    ("blastp",      ["-version"],   r"blastp: ([\d.+]+)",           "S1/S5"),
    ("tblastn",     ["-version"],   r"tblastn: ([\d.+]+)",          "S5"),
    ("makeblastdb", ["-version"],   r"makeblastdb: ([\d.+]+)",      "S5"),
    ("miniprot",    ["--version"],  None,                           "S5/S10"),
    ("trimal",      ["--version"],  r"trimAl (v[\d.a-z]+)",         "S6"),
    ("iqtree2",     ["--version"],  r"version ([\d.]+)",            "S7"),
    ("datasets",    ["--version"],  r"version: ([\d.]+)",           "S4"),
    ("foldseek",    ["version"],    None,                           "S11 (opt)"),
]

#: Python packages the pipeline imports.
PY_PACKAGES = ["Bio", "requests", "matplotlib", "numpy", "scipy", "pandas"]


def _which(name: str) -> tuple[str, str]:
    """Return (path, location) where location is 'PATH', 'env' or ''."""
    on_path = shutil.which(name)
    if on_path:
        return on_path, "PATH"
    in_env = ENV_BIN / name
    if in_env.exists():
        return str(in_env), "env"
    return "", ""


def probe(name: str, argv: list[str], pattern: str | None) -> dict:
    path, where = _which(name)
    if not path:
        return {"tool": name, "version": "", "status": "MISSING",
                "path": "", "location": ""}
    try:
        p = subprocess.run([path, *argv], capture_output=True, text=True,
                           timeout=60)
    except (OSError, subprocess.SubprocessError) as exc:
        return {"tool": name, "version": "", "status": f"ERROR: {exc}",
                "path": path, "location": where}
    blob = (p.stdout or "") + "\n" + (p.stderr or "")
    version = ""
    if pattern:
        m = re.search(pattern, blob)
        version = m.group(1) if m else ""
    if not version:
        version = next((ln.strip() for ln in blob.splitlines() if ln.strip()), "")
    return {"tool": name, "version": version[:80], "status": "ok",
            "path": path, "location": where}


def py_versions() -> list[tuple[str, str]]:
    out = []
    for mod in PY_PACKAGES:
        try:
            m = __import__(mod)
            out.append((mod, getattr(m, "__version__", "?")))
        except ImportError:
            out.append((mod, "MISSING"))
    return out


def main() -> int:
    rows = [probe(n, argv, pat) for n, argv, pat, _ in TOOLS]
    need = {n: task for n, _, _, task in TOOLS}
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["tool", "version", "status", "location", "needed_by", "path"]
    with OUT_TSV.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            r["needed_by"] = need[r["tool"]]
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    missing = [r["tool"] for r in rows if r["status"] != "ok"]
    lines = [
        "S1 toolchain manifest — IP3R publication project",
        f"generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"host:      {platform.platform()} / {platform.machine()}",
        f"python:    {sys.version.split()[0]}  ({sys.executable})",
        f"conda env: {ENV_NAME}  ({ENV_BIN.parent})",
        "",
        "External binaries",
        "-----------------",
    ]
    for r in rows:
        mark = "ok " if r["status"] == "ok" else "!! "
        lines.append(f"{mark}{r['tool']:<12} {r['version']:<34} "
                     f"[{r['location'] or '-'}] needed by {need[r['tool']]}")
        lines.append(f"   {r['path'] or '(not found)'}")
    lines += ["", "Python packages", "---------------"]
    for mod, ver in py_versions():
        lines.append(f"   {mod:<12} {ver}")
    lines += [
        "",
        "Notes",
        "-----",
        f"* Tools marked [env] are not on the bare PATH; they resolve inside",
        f"  the `{ENV_NAME}` conda env. Run project scripts with",
        f"  {ENV_BIN / 'python'} so both the packages and these binaries are",
        "  available (see Decisions D18).",
        "* Tools marked [PATH] come from Homebrew and work in any shell.",
    ]
    if missing:
        lines.append(f"* MISSING / broken: {', '.join(missing)}")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[toolchain] → {OUT_TXT}\n[toolchain] → {OUT_TSV}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
