"""Shared helpers for the S3 profile-HMM sweep.

Used by s3_build_seed.py, s3_fetch_proteomes.py, s3_run_sweep.py,
s3_assign.py, s3_census_v3.py, s3_figures.py and s3_report.py.

Responsibilities: FASTA IO, UniProt FASTA-header parsing, HMMER domtblout
parsing, jackhmmer per-round log parsing (including the per-round included
target lists D10's kill criterion needs), and loading the census v2 tables.
Nothing else — no searching, no calling, no reporting.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

CENSUS_V2_DIR = PROJECT_ROOT / "results" / "census_v2"
HMM_SWEEP_DIR = PROJECT_ROOT / "results" / "hmm_sweep"
CENSUS_V3_DIR = PROJECT_ROOT / "results" / "census_v3"


# ------------------------------------------------------------------ FASTA IO
def read_fasta(path: Path) -> dict[str, str]:
    """Full-header → sequence, order-preserving."""
    out: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    out[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        out[header] = "".join(chunks)
    return out


def write_fasta(path: Path, items: list[tuple[str, str]], width: int = 60) -> None:
    with path.open("w") as f:
        for header, seq in items:
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def iter_fasta(path: Path):
    """Stream (header, sequence) pairs — for the multi-GB sweep DB."""
    header = None
    chunks: list[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


# ------------------------------------------- UniProt proteome FASTA headers
# e.g. "sp|Q14643|ITPR1_HUMAN Inositol 1,4,5-trisphosphate receptor type 1
#       OS=Homo sapiens OX=9606 GN=ITPR1 PE=1 SV=3"
_OS_RE = re.compile(r"\bOS=(.+?)\s+(?:OX|GN|PE|SV)=")
_OX_RE = re.compile(r"\bOX=(\d+)")
_GN_RE = re.compile(r"\bGN=(\S+)")


def parse_uniprot_header(header: str) -> dict:
    """Parse a UniProt FASTA header into accession / gene / species / etc."""
    name_tok = header.split()[0] if header.split() else header
    parts = name_tok.split("|")
    acc = parts[1] if len(parts) >= 3 else name_tok
    desc = header[len(name_tok):].strip()
    os_m = _OS_RE.search(header)
    gn_m = _GN_RE.search(header)
    ox_m = _OX_RE.search(header)
    return {
        "name": name_tok,
        "accession": acc,
        "reviewed": parts[0] == "sp" if len(parts) >= 3 else False,
        "gene": gn_m.group(1) if gn_m else "",
        "species": os_m.group(1) if os_m else "",
        "taxon_id": int(ox_m.group(1)) if ox_m else None,
        "protein_name": desc.split(" OS=")[0] if " OS=" in desc else desc,
    }


# ------------------------------------------------------------- domtblout IO
DOMTBL_FIELDS = [
    "target_name", "target_acc", "tlen", "query_name", "query_acc", "qlen",
    "full_evalue", "full_score", "full_bias", "dom_num", "dom_of",
    "dom_cvalue", "dom_ivalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "acc",
]
_INT_FIELDS = {"tlen", "qlen", "dom_num", "dom_of", "hmm_from", "hmm_to",
               "ali_from", "ali_to", "env_from", "env_to"}
_FLOAT_FIELDS = {"full_evalue", "full_score", "full_bias", "dom_cvalue",
                 "dom_ivalue", "dom_score", "dom_bias", "acc"}


def parse_domtblout(path: Path) -> list[dict]:
    """Parse a --domtblout file → one dict per domain row."""
    rows: list[dict] = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            toks = line.split(None, len(DOMTBL_FIELDS))  # last col = free text
            row = dict(zip(DOMTBL_FIELDS, toks[:len(DOMTBL_FIELDS)]))
            row["description"] = toks[len(DOMTBL_FIELDS)].strip() \
                if len(toks) > len(DOMTBL_FIELDS) else ""
            for k in _INT_FIELDS:
                row[k] = int(row[k])
            for k in _FLOAT_FIELDS:
                row[k] = float(row[k])
            rows.append(row)
    return rows


def best_hits_by_target(rows: list[dict]) -> dict[str, dict]:
    """Collapse domain rows → per-target summary: best (lowest) full-sequence
    E-value, best full-sequence bit score, merged HMM coverage, domain count.

    The bit score is what D14's profile margin is measured on, so it is taken
    from the full-sequence column (identical across a target's domain rows)
    rather than summed over domains.
    """
    out: dict[str, dict] = {}
    for r in rows:
        t = r["target_name"]
        cur = out.get(t)
        if cur is None:
            cur = out[t] = {
                "target_name": t, "tlen": r["tlen"], "qlen": r["qlen"],
                "full_evalue": r["full_evalue"], "full_score": r["full_score"],
                "n_domains": 0, "hmm_spans": [], "ali_spans": [],
                "description": r["description"],
            }
        cur["full_evalue"] = min(cur["full_evalue"], r["full_evalue"])
        cur["full_score"] = max(cur["full_score"], r["full_score"])
        cur["n_domains"] += 1
        cur["hmm_spans"].append((r["hmm_from"], r["hmm_to"]))
        cur["ali_spans"].append((r["ali_from"], r["ali_to"]))
    for cur in out.values():
        cur["hmm_coverage"] = _merged_span(cur.pop("hmm_spans"), cur["qlen"])
        cur["target_coverage"] = _merged_span(cur.pop("ali_spans"), cur["tlen"])
    return out


def _merged_span(spans: list[tuple[int, int]], total: int) -> float:
    covered = 0
    last_end = 0
    for a, b in sorted(spans):
        a = max(a, last_end + 1)
        if b >= a:
            covered += b - a + 1
            last_end = b
    return round(covered / total, 3) if total else 0.0


# ------------------------------------------------------- jackhmmer log parse
# jackhmmer's log has NO "@@ Round: 1" marker: each round ENDS with
# "@@ New targets included: X", and "@@ Round: N" + "@@ Included in MSA: ...
# from M targets" announce the NEXT round (M = cumulative targets after
# round N-1). So the New-targets lines, in order, are the results of
# rounds 1..k.
_NEWT_RE = re.compile(r"^@@ New targets included:\s+(\d+)")
_TOTT_RE = re.compile(r"^@@ Included in MSA:.*from (\d+) targets")
_CONV_RE = re.compile(r"CONVERGED", re.IGNORECASE)
# A per-round score table starts at "Scores for complete sequences", and its
# data rows look like
#     +         0 6365.9  37.5    0 3354.3  19.8   3.0  3  tr|ACC|ID  descr
# where the leading "+" marks a target newly included this round (absent in
# round 1, where everything reported is included, and absent for targets
# carried over). Below a "------ inclusion threshold ------" divider the rows
# are reported but *not* included, so they are not part of what D10 scores.
_INCL_DIVIDER = "inclusion threshold"
_SCORES_HDR = re.compile(r"^Scores for complete sequences")
_TABLE_END = ("Domain annotation", "Internal pipeline", "Query:", "//")
_NAME_COL = 8          # after any leading +/- marker is stripped


def _parse_score_row(stripped: str) -> str | None:
    """Target name from one score-table data row, or None if it is not one."""
    toks = stripped.split()
    if toks and toks[0] in ("+", "-"):
        toks = toks[1:]
    if len(toks) > _NAME_COL and _looks_numeric(toks[0]):
        return toks[_NAME_COL]
    return None


def parse_jackhmmer_log(path: Path) -> dict:
    """Per-round inclusion counts *and* per-round included target names.

    Returns {"rounds": [{"round", "new_targets", "targets_in_msa",
    "included"}], "converged": bool}. `included` is the list of target names
    at or above the inclusion divider in that round's score table — the
    evidence D10's kill criterion is evaluated on. `targets_in_msa` is the
    cumulative count after that round (None for the final round: the log
    only states it at the start of the following round).

    jackhmmer's log has NO "@@ Round: 1" marker: each round ENDS with
    "@@ New targets included: X", and "@@ Round: N" + "@@ Included in MSA:
    ... from M targets" announce the NEXT round. So the New-targets lines,
    in order, are the results of rounds 1..k.
    """
    new_targets: list[int] = []
    cumulative: list[int] = []
    per_round: list[list[str]] = []
    converged = False

    in_table = False
    below_threshold = False
    current: list[str] = []
    with path.open(errors="replace") as f:
        for line in f:
            if _SCORES_HDR.match(line):
                in_table, below_threshold = True, False
                continue
            if in_table:
                stripped = line.strip()
                if not stripped:
                    continue
                if _INCL_DIVIDER in stripped:
                    below_threshold = True
                    continue
                if stripped.startswith(_TABLE_END):
                    in_table = False
                    # fall through: this line may itself be a round marker
                elif not below_threshold:
                    name = _parse_score_row(stripped)
                    if name:
                        current.append(name)
                    continue
                else:
                    continue
            m = _NEWT_RE.match(line)
            if m:
                new_targets.append(int(m.group(1)))
                per_round.append(current)
                current = []
                continue
            m = _TOTT_RE.match(line)
            if m:
                cumulative.append(int(m.group(1)))
                continue
            if _CONV_RE.search(line):
                converged = True

    # A table with no closing "@@ New targets included" line is a round that
    # did not finish — a truncated log, or a run still going. It is NOT
    # appended as a round with 0 new targets: that is exactly what a
    # *converged* final round looks like, and the completeness argument
    # turns on telling those two apart. It is reported separately instead.
    trailing = current if len(per_round) == len(new_targets) else []

    rounds = [
        {"round": i, "new_targets": nt,
         "targets_in_msa": cumulative[i - 1] if i - 1 < len(cumulative) else None,
         "included": per_round[i - 1] if i - 1 < len(per_round) else []}
        for i, nt in enumerate(new_targets, 1)
    ]
    return {"rounds": rounds, "converged": converged,
            "trailing_incomplete_round": len(trailing)}


def _looks_numeric(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------- census IO
def acc_key(acc: str) -> str:
    return acc.split(".")[0].split("-")[0]


def load_census_v2() -> list[dict]:
    """Census v2 rows (results/census_v2/census_v2.tsv) as dicts."""
    with (CENSUS_V2_DIR / "census_v2.tsv").open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def census_index(rows: list[dict]) -> dict[str, dict]:
    return {acc_key(r["accession"]): r for r in rows}


def read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
