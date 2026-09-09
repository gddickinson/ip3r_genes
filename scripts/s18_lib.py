"""S18's data layer — the loci, the annotations that hold them, and the joins.

S18 audits two records of the same genes against each other:

* what the **assemblies' own annotations** say (`genomic.gff.gz`, still on the
  drive beside every genome the sweep fetched), and
* what the **protein databases** say (census v6's UniProt-sourced rows),

against one common evidence set — the S5 sweep's alignments. Both readings are
loaded here so no later module can invent a second way of reading either.

Three things are deliberately *not* re-derived here:

* the family and paralog a **name** claims comes from `s5_classify.name_family`
  / `name_paralog` unchanged, because the brief requires the same verdict rule
  applied to genome gene names and UniProt gene fields, and two copies of a
  name rule reliably drift;
* the loci come from the archived per-genome `summary.json` and **not** the
  ledger, which holds one row per genome × cell and therefore only the best
  locus (S8, S15a and S16 each hit this); and
* the annotation parse is `s10_gff`'s, through its multi-window reader, so
  what counts as a gene here and in S10 is one implementation.
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_gff                                                    # noqa: E402
import s15_lib                                                    # noqa: E402
from s5_classify import name_family, name_paralog                 # noqa: E402

PROJECT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT / "results"
OUT = RESULTS / "annotation_audit"
FIGS = OUT / "figures"

MANIFEST = RESULTS / "genome_manifest.tsv"
LEDGER = RESULTS / "genome_ledger" / "genome_ledger.tsv"
INTEGRITY = RESULTS / "loss_dynamics" / "integrity_loci.tsv"
CENSUS_V6 = RESULTS / "census_v6" / "census_v6.tsv"
CENSUS_V3 = RESULTS / "census_v3" / "census_v3.tsv"
BAITS = RESULTS / "s5_baits" / "baits.faa"
BAIT_MANIFEST = RESULTS / "s5_baits" / "bait_manifest.tsv"
HMM_ASSIGN = RESULTS / "hmm_sweep" / "hmmsearch_assignments.tsv"
ZERO_HIT = RESULTS / "census_v3" / "proteomes_without_hits.tsv"
PROTEOME_MANIFEST = RESULTS / "hmm_sweep" / "proteome_manifest.tsv"

ITPR_CELLS = ("ITPR1", "ITPR2", "ITPR3")
CONTROL_CELL = "RYR"
ALL_CELLS = ITPR_CELLS + (CONTROL_CELL,)

# The naming signature: the one Pfam that names this family. S2 measured that
# 2,911 of 15,417 seeded-space proteins do not carry it, which is why the
# Pfam-recall intersection is a result and not a formality.
NAMING_PFAM = "PF08709"
FAMILY_PFAMS = ("PF08709", "PF02815", "PF01365", "PF08454")

read_tsv = s15_lib.read_tsv


def write_tsv(path, rows: list[dict], cols: list[str] | None = None):
    """`s15_lib.write_tsv`, with the column list derivable from the rows.

    S18 writes several tables whose columns depend on which verdicts occurred,
    so the header is the union of the rows' keys **in order of first
    appearance** — deterministic, and never a set (a hash-ordered header is
    the nondeterminism `s8_test_flanks` T11 caught in a ranked table).
    """
    if cols is None:
        cols, seen = [], set()
        for r in rows:
            for k in r:
                if k not in seen and not k.startswith("_"):
                    seen.add(k)
                    cols.append(k)
    return s15_lib.write_tsv(path, rows, cols)
sha256 = s15_lib.sha256
data_root = s15_lib.data_root
sweep_dir = s15_lib.sweep_dir
load_summaries = s15_lib.load_summaries


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    """The dashboard live panel, written under this task's own id."""
    path = RESULTS / "session_live.json"
    path.write_text(json.dumps({
        "task": "S18", "stage": stage,
        "steps": [{"label": a, "done": b} for a, b in steps]}, indent=1))


# --------------------------------------------------------------------------
# scope
# --------------------------------------------------------------------------
def genome_rows() -> dict[str, dict]:
    """The declared denominator (S4), keyed by accession."""
    return {r["accession"]: r for r in read_tsv(MANIFEST)}


def annotation_source(acc: str) -> str:
    """D9: RefSeq and submitter-deposited gene sets are not the same evidence.

    Read off the accession prefix rather than the manifest's
    `annotation_source` column, because that column records *who produced* the
    gene set and this one records *which archive serves it* — an assembly can
    carry an NCBI-produced annotation under a `GCA_` accession, and the reader
    a correction is addressed to is the one holding the record.
    """
    return "RefSeq" if acc.startswith("GCF_") else "GenBank"


def gff_path(acc: str) -> Path | None:
    """The assembly's own annotation, as fetched with the genome."""
    p = data_root() / "genomes" / acc / "genomic.gff.gz"
    return p if p.exists() else None


def loci_of(acc: str, summary: dict) -> list[dict]:
    """Every gene-scale locus the sweep placed in this genome, flattened.

    One row per locus and not per cell: the cell's best locus is what the
    ledger keeps, and a second copy the annotation holds differently is
    exactly the thing this audit is looking for.
    """
    out = []
    for cell, c in (summary.get("cells") or {}).items():
        for i, loc in enumerate(c.get("loci") or []):
            d = dict(loc)
            d.update(accession=acc, cell=cell, locus_idx=i,
                     cell_status=c.get("status", ""),
                     is_control=bool(c.get("is_control")),
                     organism=summary.get("organism", ""),
                     vclass=summary.get("vclass", ""),
                     vorder=summary.get("vorder", ""),
                     assembly_level=summary.get("assembly_level", ""),
                     contig_n50=summary.get("contig_n50", 0),
                     annotated=summary.get("annotated", ""))
            out.append(d)
    return out


# --------------------------------------------------------------------------
# the joins
# --------------------------------------------------------------------------
def integrity_index() -> dict[tuple[str, str, int], dict]:
    """S15's ORF-integrity verdict per locus — the D6 veto's evidence."""
    out = {}
    for r in read_tsv(INTEGRITY):
        out[(r["accession"], r["cell"], int(r["locus_idx"]))] = r
    return out


def contiguity_index() -> dict[str, dict]:
    """The ledger row per genome × cell: D4's bar and the sweep's own status."""
    out = {}
    for r in read_tsv(LEDGER):
        out[(r["accession"], r["class"])] = r
    return out


def bait_labels() -> dict[str, dict]:
    """The committed panel's label per bait **accession** — family, paralog.

    Keyed on the accession and not the manifest's `id`, because the panel's
    FASTA headers come in two shapes (`ACC|GENE|Species` and
    `ACC|GENE|Species|PARALOG`) and only the accession is common to both. The
    first version keyed on `id`, every lookup missed, and the sequence call
    came back `no_call` on all 11,402 records — a failure that looks exactly
    like a family nothing can be assigned to.
    """
    return {r["accession"]: r for r in read_tsv(BAIT_MANIFEST)}


def panel_paralogs() -> dict[str, set[str]]:
    """Which paralogs the committed panel can actually call, per family.

    The panel carries no RYR3 bait (S5's slot table records the six unfilled
    slots), so a RYR3 record cannot be called RYR3 by this instrument however
    well it is annotated. Reporting a wrong-paralog verdict there would be a
    fact about the bait panel, so the audit records what the panel can reach
    and scores the paralog question only inside it.
    """
    out: dict[str, set[str]] = {}
    for r in read_tsv(BAIT_MANIFEST):
        if r.get("paralog"):
            out.setdefault(r["family"], set()).add(r["paralog"])
    return out


def bait_id(header: str) -> str:
    """The panel id inside a sweep locus's `bait` field (`ACC|GENE|Species…`)."""
    return (header or "").split("|")[0]


# --------------------------------------------------------------------------
# sequence + FASTA
# --------------------------------------------------------------------------
def iter_fasta(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif name is not None:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path: Path, seqs: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return path


def tool_bin(name: str) -> str:
    """PATH first, then the recorded conda envs (D18)."""
    from shutil import which
    p = which(name)
    if p:
        return p
    for env in ("piezo1", "s12"):
        cand = Path(f"/opt/anaconda3/envs/{env}/bin/{name}")
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(name)


def blastp_vs_panel(query: Path, panel_db: Path, out: Path,
                    threads: int = 4) -> Path:
    """blastp of a query set against the labelled bait panel.

    `-max_target_seqs` is the whole panel, because the paralog call needs the
    best hit in *each* family, not the best hit overall: a record whose top
    hit is an ITPR1 bait still has to be scored against the RyR baits before
    D7's margin means anything.
    """
    if out.exists() and out.stat().st_size:
        return out
    cmd = [tool_bin("blastp"), "-query", str(query), "-db", str(panel_db),
           "-outfmt", "6 qseqid sseqid pident length bitscore evalue qlen "
                      "slen qstart qend sstart send",
           "-max_target_seqs", "500", "-evalue", "1e-5",
           "-num_threads", str(threads), "-out", str(out)]
    subprocess.run(cmd, check=True)
    return out


def make_panel_db(faa: Path, db_dir: Path) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)
    db = db_dir / faa.stem
    if not (db_dir / (faa.stem + ".phr")).exists():
        subprocess.run([tool_bin("makeblastdb"), "-in", str(faa),
                        "-dbtype", "prot", "-out", str(db)],
                       check=True, capture_output=True)
    return db


# --------------------------------------------------------------------------
# small statistics, stdlib
# --------------------------------------------------------------------------
def quantile(xs: list[float], q: float) -> float:
    if not xs:
        return float("nan")
    s = sorted(xs)
    if len(s) == 1:
        return s[0]
    pos = q * (len(s) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def frac(n: int, d: int) -> float:
    return round(n / d, 4) if d else float("nan")


def fisher_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p, stdlib — the by-source comparison's test.

    D9 says a RefSeq gene set and a submitter-deposited one are not comparable
    evidence; showing they differ needs a test, and the counts here are small
    enough in some cells that a chi-square would not do.
    """
    from math import lgamma, exp

    def logc(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

    n = a + b + c + d
    r1, c1 = a + b, a + c
    base = logc(r1, a) + logc(n - r1, c) - logc(n, c1)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    tot = 0.0
    for k in range(lo, hi + 1):
        lp = logc(r1, k) + logc(n - r1, c1 - k) - logc(n, c1)
        if lp <= base + 1e-9:
            tot += exp(lp)
    return min(1.0, tot)


def benjamini_hochberg(ps: list[float]) -> list[float]:
    idx = sorted(range(len(ps)), key=lambda i: ps[i])
    out = [0.0] * len(ps)
    m = len(ps)
    prev = 1.0
    for rank, i in enumerate(reversed(idx), start=1):
        q = min(prev, ps[i] * m / (m - rank + 1))
        out[i] = round(min(1.0, q), 6)
        prev = q
    return out


__all__ = [n for n in dir() if not n.startswith("_")]
