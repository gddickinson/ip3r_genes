"""S20 step 3 (evidence) — what can be known about one plant or fungal hit.

The brief's sharpest question is whether the family's few plant and fungal
records are facts about genomes or facts about databases. That is not
answerable from a bit score, so this module gathers the four independent
lines of evidence a verdict can rest on, and gathers them for *every*
candidate whether or not they agree:

1. **Sequence** — pulled from the archived proteome DBs and S2's archived
   seeded-space FASTA, so the whole chase reruns offline.
2. **The record's own paperwork** — UniProt's `fragment`, protein-existence
   level, Pfam architecture, proteome membership and EMBL cross-references,
   in one batched stream query archived under the data root. A record with
   no genome cross-reference is not a gene anyone has seen in an assembly.
3. **The nearest neighbour outside its own kingdom** — blastp against the
   pooled ITPR calls from every other group. This is the contamination
   test, and it is the only one that can distinguish "a green alga has an
   IP₃ receptor" from "a metazoan sequence is sitting in a green alga's
   assembly": a genuine deep homolog is 20–40 % identical to its metazoan
   relatives, while an assembly contaminant is 95–100 % identical to one
   particular animal.
4. **Its placement in its own group** — the best hit *within* its kingdom,
   so a record that is one of many similar plant proteins reads differently
   from one that is unique in the entire green lineage.

Nothing here decides anything; `s20_verdicts.py` applies the rules.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
import urllib.parse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_lib import log, require_data_root  # noqa: E402
from scripts.s3_hmm_lib import iter_fasta, parse_uniprot_header  # noqa: E402

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
DETAIL_FIELDS = ("accession,reviewed,protein_name,gene_names,organism_name,"
                 "organism_id,length,fragment,protein_existence,xref_pfam,"
                 "xref_proteomes,xref_embl,lineage,cc_caution,annotation_score")

BLAST_FIELDS = ("qseqid sseqid pident length qlen slen qstart qend "
                "sstart send evalue bitscore")
BLASTP = "blastp"
MAKEBLASTDB = "makeblastdb"


def _env_bin(name: str) -> str:
    """BLAST+ is env-resident in this project (S1's toolchain manifest)."""
    env = Path("/opt/anaconda3/envs/piezo1/bin") / name
    return str(env) if env.exists() else name


def work_dir() -> Path:
    d = require_data_root() / "hmmer" / "s20" / "verdicts"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------- 1. sequences
def extract(fasta: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    """One streaming pass over a big FASTA → acc → (header, sequence).

    Keyed on the UniProt accession parsed out of the header, not on the
    whole header, because the same protein is written `sp|ACC|NAME` in the
    proteome files and `tr|ACC|NAME` in others.
    """
    found: dict[str, tuple[str, str]] = {}
    if not fasta.exists() or not wanted:
        return found
    for header, seq in iter_fasta(fasta):
        acc = parse_uniprot_header(header)["accession"]
        if acc in wanted and acc not in found:
            found[acc] = (header, seq)
            if len(found) == len(wanted):
                break
    return found


def gather_sequences(wanted: set[str], sources: list[Path]) -> dict[str, tuple[str, str]]:
    """Walk the archived FASTAs in order until every accession is found."""
    out: dict[str, tuple[str, str]] = {}
    for src in sources:
        todo = wanted - set(out)
        if not todo:
            break
        got = extract(src, todo)
        if got:
            log("s20_evidence", f"  {len(got)} sequences from {src.name}")
        out.update(got)
    missing = wanted - set(out)
    if missing:
        log("s20_evidence", f"  {len(missing)} accessions not in any archived "
                            f"FASTA: {sorted(missing)[:5]}")
    return out


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> int:
    with path.open("w") as f:
        for acc, seq in sorted(seqs.items()):
            f.write(f">{acc}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")
    return len(seqs)


# ------------------------------------------------------------- 2. paperwork
def uniprot_detail(accs: list[str], archive: Path,
                   batch: int = 200) -> dict[str, dict]:
    """UniProt's own record for each accession, batched and archived."""
    import requests
    rows: dict[str, dict] = {}
    header: list[str] = []
    parts: list[str] = []
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        # Keyed on the chunk's *contents*, not its position in the list: the
        # candidate set grows as groups finish sweeping, and an index-keyed
        # cache would serve one run's chunk 3 to the next run's chunk 3.
        key = hashlib.sha256("\n".join(chunk).encode()).hexdigest()[:16]
        dest = archive / f"detail_{key}.tsv"
        if dest.exists() and dest.stat().st_size:
            text = dest.read_text()
        else:
            query = " OR ".join(f"accession:{a}" for a in chunk)
            url = (f"{UNIPROT_STREAM}?query={urllib.parse.quote(query)}"
                   f"&format=tsv&fields={DETAIL_FIELDS}")
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            text = r.text
            dest.write_text(text)
        lines = text.splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        parts += lines[1:]
    for line in parts:
        if not line.strip():
            continue
        vals = line.split("\t")
        rec = dict(zip(header, vals + [""] * (len(header) - len(vals))))
        rows[rec.get("Entry", "")] = rec
    log("s20_evidence", f"  UniProt detail for {len(rows)}/{len(accs)} records")
    return rows


# ------------------------------------------------- 3/4. nearest neighbours
def make_db(fasta: Path) -> Path:
    subprocess.run([_env_bin(MAKEBLASTDB), "-in", str(fasta), "-dbtype",
                    "prot", "-out", str(fasta.with_suffix(""))],
                   check=True, capture_output=True, text=True)
    return fasta.with_suffix("")


def blastp(query: Path, db: Path, out: Path, max_target: int = 50,
           threads: int = 6) -> list[dict]:
    """blastp query→db, tabular. Cached: an existing non-empty out is reused."""
    if not (out.exists() and out.stat().st_size):
        tmp = out.with_suffix(".part")
        subprocess.run([_env_bin(BLASTP), "-query", str(query), "-db", str(db),
                        "-outfmt", f"6 {BLAST_FIELDS}", "-evalue", "1e-3",
                        "-max_target_seqs", str(max_target),
                        "-num_threads", str(threads), "-out", str(tmp)],
                       check=True, capture_output=True, text=True)
        tmp.rename(out)
    names = BLAST_FIELDS.split()
    rows = []
    for line in out.read_text().splitlines():
        if not line.strip():
            continue
        vals = line.split("\t")
        row = dict(zip(names, vals))
        for k in ("pident", "evalue", "bitscore"):
            row[k] = float(row[k])
        for k in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send"):
            row[k] = int(row[k])
        row["qcov"] = round((row["qend"] - row["qstart"] + 1) / row["qlen"], 3)
        rows.append(row)
    return rows


def best_by_partition(hits: list[dict], partition: dict[str, str],
                      query_partition: dict[str, str]) -> dict[str, dict]:
    """Per query: best hit inside its own partition and best outside it.

    `partition` maps a subject accession to its coarse group; the same for
    queries. Self-hits are excluded — a record is not its own evidence.
    """
    out: dict[str, dict] = {}
    for h in sorted(hits, key=lambda r: -r["bitscore"]):
        q, s = h["qseqid"], h["sseqid"]
        if q == s:
            continue
        rec = out.setdefault(q, {})
        same = partition.get(s, "?") == query_partition.get(q, "??")
        key = "within" if same else "outside"
        if key not in rec:
            rec[key] = h
    return out
