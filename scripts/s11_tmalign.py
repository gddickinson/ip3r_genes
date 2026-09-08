"""S11 — TM-align execution and parsing.

TM-align prints **two** TM-scores, normalised by each input's length. For a
full-length predicted model against a *partially resolved* cryo-EM
structure those differ a great deal, and which one is meaningful depends on
the question being asked — so both are kept, and the report says which it
uses for which claim.

The two bars this task reads scores against are TM-align's own published
ones (Zhang & Skolnick 2005; Xu & Zhang 2010), not numbers this project
chose: below `TM_RANDOM` a pair is statistically indistinguishable from two
random structures of the same size; above `TM_FOLD` they share a fold.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

TM_RANDOM = 0.17
TM_FOLD = 0.50


def tmalign_bin() -> str:
    """D18: the structure toolchain is env-resident, not on the bare PATH."""
    env = Path("/opt/anaconda3/envs/piezo1/bin/TMalign")
    return str(env) if env.exists() else "TMalign"


@dataclass
class TmResult:
    query: str
    target: str
    tm_query: float = 0.0     # normalised by query length
    tm_target: float = 0.0    # normalised by target length
    rmsd: float = 0.0
    aligned: int = 0
    seq_id: float = 0.0
    len_query: int = 0
    len_target: int = 0
    ok: bool = False
    error: str = ""

    @property
    def tm_max(self) -> float:
        return max(self.tm_query, self.tm_target)


def run_tmalign(query_pdb: Path, target_pdb: Path, label_q: str = "",
                label_t: str = "", timeout: int = 3600) -> TmResult:
    """Run TM-align on two CA-trace PDBs and parse its report.

    TM-align prints two TM-scores, normalised by each input's length. For a
    full-length model against a *partially resolved* cryo-EM structure
    those differ a lot, and which one is meaningful depends on the
    question — so both are kept and the report says which it uses.
    """
    res = TmResult(query=label_q or query_pdb.stem,
                   target=label_t or target_pdb.stem)
    try:
        proc = subprocess.run(
            [tmalign_bin(), str(query_pdb), str(target_pdb)],
            capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        res.error = "timeout"
        return res
    except FileNotFoundError:
        res.error = "TMalign not found"
        return res
    if proc.returncode != 0:
        res.error = (proc.stderr or proc.stdout)[:200]
        return res
    # TM-align has used both "Chain_1" and "Structure_1" wording across
    # versions (20240303 prints "Chain_1" in the length lines and
    # "Chain_1" in the TM-score lines); match either, or the parse
    # silently yields zeros for every pair.
    len_re = re.compile(r"Length of (?:Chain|Structure)[_ ](\d)\s*:\s*(\d+)")
    tm_re = re.compile(
        r"TM-score\s*=\s*([\d.]+).*normalized by length of "
        r"(?:Chain|Structure)[_ ](\d)")
    for line in proc.stdout.splitlines():
        s = line.strip()
        m = len_re.search(s)
        if m:
            if m.group(1) == "1":
                res.len_query = int(m.group(2))
            else:
                res.len_target = int(m.group(2))
            continue
        m = tm_re.search(s)
        if m:
            if m.group(2) == "1":
                res.tm_query = float(m.group(1))
            else:
                res.tm_target = float(m.group(1))
            continue
        if s.startswith("Aligned length"):
            parts = [p.strip() for p in s.split(",")]
            try:
                res.aligned = int(parts[0].split("=")[1])
                res.rmsd = float(parts[1].split("=")[1])
                res.seq_id = float(parts[2].split("=")[-1])
            except (IndexError, ValueError):
                pass
    res.ok = (res.tm_query > 0 or res.tm_target > 0) and res.len_query > 0
    if not res.ok and not res.error:
        res.error = "no TM-score parsed"
    return res


def tmalign_version() -> str:
    try:
        proc = subprocess.run([tmalign_bin()], capture_output=True, text=True,
                              timeout=60, check=False)
    except Exception:  # noqa: BLE001
        return "unavailable"
    m = re.search(r"TM-align\s*\(Version\s*([\w-]+)\)", proc.stdout)
    return m.group(1) if m else "unknown"


