"""s5_genome_io.py — genome-file plumbing for the S5 sweep.

FASTA offset index (faidx-style), region fetch, N-run detection,
NCBI sequence/assembly report parsing, and GFF3 annotation download
(datasets --include gff3) slimmed to a per-gene table.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import zipfile
from bisect import bisect_right
from pathlib import Path

GFF_ZIP_TMP = "_gff_tmp"


# ---------------------------------------------------------------- FASTA index

def build_fai(fna: Path) -> Path:
    """Build (or reuse) a faidx-style index: name len offset linebases linewidth."""
    fai = fna.with_suffix(fna.suffix + ".s5fai")
    if fai.exists() and fai.stat().st_mtime >= fna.stat().st_mtime:
        return fai
    rows = []
    with open(fna, "rb") as fh:
        name, seqlen, offset, linebases, linewidth = None, 0, 0, 0, 0
        pos = 0
        for raw in fh:
            n = len(raw)
            if raw.startswith(b">"):
                if name is not None:
                    rows.append((name, seqlen, offset, linebases, linewidth))
                name = raw[1:].split()[0].decode()
                seqlen, offset, linebases, linewidth = 0, pos + n, 0, 0
            else:
                stripped = len(raw.rstrip(b"\r\n"))
                if linebases == 0:
                    linebases, linewidth = stripped, n
                seqlen += stripped
            pos += n
        if name is not None:
            rows.append((name, seqlen, offset, linebases, linewidth))
    with open(fai, "w") as out:
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    return fai


def read_fai(fai: Path) -> dict[str, tuple[int, int, int, int]]:
    idx = {}
    for line in fai.read_text().splitlines():
        name, ln, off, lb, lw = line.split("\t")
        idx[name] = (int(ln), int(off), int(lb), int(lw))
    return idx


def fetch_region(fna: Path, idx: dict, contig: str, start: int, end: int) -> str:
    """1-based inclusive region fetch."""
    if contig not in idx:
        return ""
    ln, off, lb, lw = idx[contig]
    start, end = max(1, start), min(ln, end)
    if start > end or lb == 0:
        return ""
    first = off + (start - 1) // lb * lw + (start - 1) % lb
    last = off + (end - 1) // lb * lw + (end - 1) % lb + 1
    with open(fna, "rb") as fh:
        fh.seek(first)
        chunk = fh.read(last - first)
    return chunk.replace(b"\n", b"").replace(b"\r", b"").decode()


def longest_n_run(seq: str) -> int:
    best = run = 0
    for ch in seq:
        if ch in "Nn":
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


# ------------------------------------------------------------- NCBI reports

def read_sequence_lengths(genome_dir: Path) -> dict[str, int]:
    """Map every sequence accession (genbank + refseq) to its length."""
    path = genome_dir / "sequence_report.jsonl"
    lengths: dict[str, int] = {}
    if not path.exists():
        return lengths
    with open(path) as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            for key in ("genbankAccession", "refseqAccession"):
                acc = row.get(key)
                if acc and acc != "na":
                    lengths[acc] = int(row.get("length", 0))
    return lengths


def read_assembly_meta(genome_dir: Path) -> dict:
    path = genome_dir / "assembly_data_report.jsonl"
    if not path.exists():
        return {}
    try:
        row = json.loads(path.read_text().splitlines()[0])
    except (json.JSONDecodeError, IndexError):
        return {}
    info = row.get("assemblyInfo", {})
    org = row.get("organism", {})
    return {
        "assembly_name": info.get("assemblyName"),
        "assembly_level": info.get("assemblyLevel"),
        "organism": org.get("organismName"),
        "taxid": org.get("taxId"),
    }


# ------------------------------------------------------------- annotation

def ensure_annotation(accession: str, genome_dir: Path, quiet: bool = False) -> Path | None:
    """Download the assembly's GFF3 (once), store gzipped; return path or None."""
    gz = genome_dir / "genomic.gff.gz"
    marker = genome_dir / ".gff.done"
    if marker.exists():
        return gz if gz.exists() else None
    tmp = genome_dir / GFF_ZIP_TMP
    zip_path = genome_dir / "_gff.zip"
    cmd = ["datasets", "download", "genome", "accession", accession,
           "--include", "gff3", "--filename", str(zip_path), "--no-progressbar"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        if not quiet:
            print(f"  ! gff3 download failed for {accession}: "
                  f"{proc.stderr.strip()[:200]}")
        marker.write_text("failed\n")
        zip_path.unlink(missing_ok=True)
        return None
    if tmp.exists():
        shutil.rmtree(tmp)
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.startswith(("/", "..")) or ".." in Path(name).parts:
                raise RuntimeError(f"unsafe path in zip: {name}")
        zf.extractall(tmp)
    src = tmp / "ncbi_dataset" / "data" / accession / "genomic.gff"
    if not src.exists():
        marker.write_text("no-gff-in-zip\n")
        shutil.rmtree(tmp, ignore_errors=True)
        zip_path.unlink(missing_ok=True)
        return None
    with open(src, "rb") as fin, gzip.open(gz, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    shutil.rmtree(tmp, ignore_errors=True)
    zip_path.unlink(missing_ok=True)
    marker.write_text("ok\n")
    return gz


def _attr(attrs: str, key: str) -> str:
    for part in attrs.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


def slim_genes(gff_gz: Path) -> list[tuple]:
    """Extract gene features → [(contig, start, end, strand, gene_id, name, biotype)]."""
    genes = []
    with gzip.open(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            attrs = f[8]
            genes.append((f[0], int(f[3]), int(f[4]), f[6],
                          _attr(attrs, "ID"),
                          _attr(attrs, "Name") or _attr(attrs, "gene"),
                          _attr(attrs, "gene_biotype") or f[2]))
    return genes


class GeneIndex:
    """Per-contig sorted gene intervals with overlap query."""

    def __init__(self, genes: list[tuple]):
        self.by_contig: dict[str, list[tuple]] = {}
        for g in genes:
            self.by_contig.setdefault(g[0], []).append(g)
        for lst in self.by_contig.values():
            lst.sort(key=lambda g: g[1])
        self.starts = {c: [g[1] for g in lst] for c, lst in self.by_contig.items()}

    def overlapping(self, contig: str, start: int, end: int) -> list[tuple]:
        lst = self.by_contig.get(contig, [])
        if not lst:
            return []
        i = bisect_right(self.starts[contig], end)
        return [g for g in lst[:i] if g[2] >= start]
