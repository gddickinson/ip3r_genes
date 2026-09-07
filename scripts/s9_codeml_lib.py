"""S9 — codeml control-file construction, execution and output parsing.

One `CodemlJob` per PAML run. Jobs are resumable (an `mlc` that already
parses is reused), run in their own directory, and are parsed into a
`CodemlResult` carrying lnL, np, kappa and the ω estimate(s).

Tree helpers cover what the branch models need: unrooting (PAML wants an
unrooted tree), and PAML clade/branch labels — `$1` marks every branch
inside a clade, `#1` marks only the stem branch leading to it.

`labelled_newick` **refuses a non-monophyletic foreground**. That is the
guard the whole branch half of S9 rests on: codeml marks a *node*, so a
foreground whose MRCA also contains other tips does not fail — it marks a
larger clade and returns a perfectly well-formed ω for a hypothesis nobody
asked about. S9's sets come from S7's clades (`s9_sets.py`) so this should
never fire; it is here so that if it ever does, it stops the run.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_lib import Node  # noqa: E402
from s9_cds_lib import tool_bin  # noqa: E402

CTL_DEFAULTS = {
    "noisy": "3", "verbose": "0", "runmode": "0", "seqtype": "1",
    "CodonFreq": "2", "clock": "0", "aaDist": "0", "icode": "0",
    "fix_kappa": "0", "kappa": "2", "fix_alpha": "1", "alpha": "0",
    "Malpha": "0", "ncatG": "3", "getSE": "0", "RateAncestor": "0",
    "Small_Diff": ".5e-6", "cleandata": "0", "method": "0",
}


# ---- tree helpers ----------------------------------------------------------

def unroot(tree: Node) -> Node:
    """PAML wants an unrooted tree: give the root three or more children."""
    if len(tree.children) != 2:
        return tree
    a, b = tree.children
    donor = a if not a.is_leaf else (b if not b.is_leaf else None)
    if donor is None:
        return tree
    other = b if donor is a else a
    other.length = (other.length or 0.0) + (donor.length or 0.0)
    tree.children = [other] + list(donor.children)
    return tree


def mrca(tree: Node, labels: set[str]) -> Node | None:
    """Smallest node containing every label in `labels`."""
    best = None
    for n in tree.walk():
        names = n.leaf_names()
        if labels <= names and (best is None or len(names) < len(best.leaf_names())):
            best = n
    return best


def labelled_newick(tree: Node, fg: set[str], mode: str) -> str:
    """Newick with a PAML label on the `fg` clade.

    mode 'clade' -> `$1` (every branch in the clade is foreground);
    mode 'stem'  -> `#1` (only the branch leading to the clade).
    """
    node = mrca(tree, fg)
    if node is None:
        raise ValueError("foreground clade not found in tree")
    # A non-monophyletic foreground would silently walk up to the root and
    # mark the entire tree, which still runs and still produces numbers.
    extra = node.leaf_names() - fg
    if extra:
        raise ValueError(
            f"foreground is not monophyletic: its MRCA also contains "
            f"{len(extra)} other tips ({sorted(extra)[:5]}...)")
    tag = "$1" if mode == "clade" else "#1"

    def rec(n: Node) -> str:
        mark = tag if n is node else ""
        if n.is_leaf:
            return f"{n.name}{mark}"
        inner = ",".join(rec(c) for c in n.children)
        return f"({inner}){mark}"

    return rec(tree) + ";"


def newick_no_lengths(tree: Node) -> str:
    def rec(n: Node) -> str:
        if n.is_leaf:
            return n.name
        return "(" + ",".join(rec(c) for c in n.children) + ")"
    return rec(tree) + ";"


# ---- jobs ------------------------------------------------------------------

@dataclass
class CodemlJob:
    name: str
    seqfile: Path
    treefile: Path
    workdir: Path
    settings: dict = field(default_factory=dict)
    description: str = ""
    pairwise: bool = False      # runmode=-2: no lnL line, matrices instead

    def ctl_text(self) -> str:
        ctl = dict(CTL_DEFAULTS)
        ctl.update({k: str(v) for k, v in self.settings.items()})
        lines = [f"seqfile = {self.seqfile}", f"treefile = {self.treefile}",
                 "outfile = mlc"]
        lines += [f"{k} = {v}" for k, v in ctl.items()]
        return "\n".join(lines) + "\n"

    def mlc(self) -> Path:
        return self.workdir / "mlc"


@dataclass
class CodemlResult:
    name: str
    lnL: float
    np: int
    kappa: float | None = None
    omegas: list[float] = field(default_factory=list)
    tree_length: float | None = None
    site_classes: str = ""
    raw: str = ""

    @property
    def omega(self) -> float | None:
        return self.omegas[0] if self.omegas else None


LNL_RE = re.compile(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)[:\s]*(-?[\d.]+)")


def parse_mlc(path: Path, name: str = "") -> CodemlResult | None:
    if not path.exists():
        return None
    text = path.read_text()
    m = LNL_RE.search(text)
    if not m:
        return None
    res = CodemlResult(name=name or path.parent.name,
                       lnL=float(m.group(3)), np=int(m.group(2)), raw=text)
    k = re.search(r"kappa \(ts/tv\)\s*=\s*([\d.]+)", text)
    if k:
        res.kappa = float(k.group(1))
    tl = re.search(r"tree length\s*=\s*([\d.]+)", text)
    if tl:
        res.tree_length = float(tl.group(1))
    # one-ratio
    w = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", text)
    if w:
        res.omegas = [float(w.group(1))]
    # branch models: "w (dN/dS) for branches:  0.05 0.21"
    wb = re.search(r"w \(dN/dS\) for branches:\s*([\d.\s]+)", text)
    if wb:
        res.omegas = [float(x) for x in wb.group(1).split()]
    # branch-site / site models
    props = re.search(r"proportion\s+([\d.\s]+)\n", text)
    bg = re.search(r"background w\s+([\d.\s]+)\n", text)
    fg = re.search(r"foreground w\s+([\d.\s]+)\n", text)
    if props and fg:
        res.site_classes = (f"p=({props.group(1).split()}) "
                            f"bg=({bg.group(1).split() if bg else []}) "
                            f"fg=({fg.group(1).split()})")
    return res


def pairwise_done(job: CodemlJob) -> bool:
    """runmode=-2 writes distance matrices, not a likelihood line."""
    m = job.workdir / "2ML.dS"
    return m.exists() and m.stat().st_size > 0


def owner_alive(lock: Path) -> int:
    """The live pid holding `lock`, or 0 (missing, malformed or stale)."""
    if not lock.exists():
        return 0
    try:
        pid = int(lock.read_text().split()[0])
    except (ValueError, IndexError, OSError):
        return 0
    if not pid:
        return 0
    try:
        os.kill(pid, 0)          # signal 0 — existence check only
        return pid
    except ProcessLookupError:
        return 0                 # the owner is gone: stale lock
    except PermissionError:
        return pid               # exists, owned by another user


def run_job(job: CodemlJob, timeout_s: int = 36_000) -> CodemlResult | None:
    """Run codeml unless the job's output is already present.

    A job directory is **claimed** for the duration of the run. codeml
    writes fixed filenames (`mlc`, `rst`, `lnf`, `2ML.dS`) into its working
    directory, so two processes on the same job interleave their output in
    silence and whichever finishes last wins some files and loses others —
    and `parse_mlc` will still find an `lnL(ntime` line and report a
    plausible number for a run that never happened that way. This is
    `s7_run.claim()`'s incident, in a second tool: there it was two IQ-TREE
    searches sharing a `--prefix`, here it is two drivers (the site models
    and the whole-tree branch models are hours apart in cost, so running
    them as separate processes is worth doing) whose job lists overlap.

    A stale lock is taken over rather than treated as an error, so an
    interrupted run needs no cleaning by hand.
    """
    job.workdir.mkdir(parents=True, exist_ok=True)
    if job.pairwise:
        if pairwise_done(job):
            return CodemlResult(name=job.name, lnL=0.0, np=0,
                                site_classes="pairwise")
    else:
        existing = parse_mlc(job.mlc(), job.name)
        if existing is not None:
            return existing
    lock = job.workdir / "RUNNING"
    pid = owner_alive(lock)
    if pid:
        print(f"  {job.name}: pid {pid} already owns this job directory — "
              "skipped", flush=True)
        return None
    lock.write_text(f"{os.getpid()}\n")
    (job.workdir / "codeml.ctl").write_text(job.ctl_text())
    try:
        proc = subprocess.run([tool_bin("codeml"), "codeml.ctl"], cwd=job.workdir,
                              capture_output=True, text=True, timeout=timeout_s)
        (job.workdir / "codeml.stdout").write_text(proc.stdout[-20000:])
        (job.workdir / "codeml.stderr").write_text(proc.stderr[-8000:])
    except subprocess.TimeoutExpired:
        (job.workdir / "TIMEOUT").write_text(f"timeout after {timeout_s}s\n")
        return None
    finally:
        lock.unlink(missing_ok=True)
    if job.pairwise:
        return (CodemlResult(name=job.name, lnL=0.0, np=0,
                             site_classes="pairwise")
                if pairwise_done(job) else None)
    return parse_mlc(job.mlc(), job.name)


# ---- likelihood-ratio tests ------------------------------------------------

def lrt(null: CodemlResult, alt: CodemlResult, df: int | None = None) -> dict:
    """2ΔlnL with a chi-square p-value (survival function, pure stdlib)."""
    import math
    d = df if df is not None else max(1, alt.np - null.np)
    stat = max(0.0, 2.0 * (alt.lnL - null.lnL))
    return {"stat": stat, "df": d, "p": _chi2_sf(stat, d)}


def _chi2_sf(x: float, k: int) -> float:
    """Upper tail of chi-square with k df (regularised incomplete gamma Q)."""
    import math
    if x <= 0:
        return 1.0
    a, xx = k / 2.0, x / 2.0
    if xx < a + 1.0:                      # series for P, then Q = 1 - P
        term, total, n = 1.0 / a, 1.0 / a, 0
        while n < 10_000:
            n += 1
            term *= xx / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, min(1.0, 1.0 - total * math.exp(-xx + a * math.log(xx)
                                                        - math.lgamma(a))))
    # continued fraction for Q
    tiny = 1e-300
    b, c, d = xx + 1.0 - a, 1.0 / tiny, 1.0 / (xx + 1.0 - a)
    h = d
    for i in range(1, 10_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return max(0.0, min(1.0, h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))))


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """BH-adjusted q-values, input order preserved."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        val = min(prev, pvals[idx] * n / i)
        q[idx] = val
        prev = val
    return q
