"""S11 — the optional Foldseek sweep over AlphaFold DB.

The brief marks this optional, and the scope it is run at is stated rather
than implied. A sweep over the *whole* of AFDB means the UniProt50 index —
about a quarter of a terabyte and many hours — which is not what this task's
completion criteria rest on. What is run here is the **Swiss-Prot subset**
of AFDB: every reviewed protein with a predicted structure, ~550k models.
That is a real structure-first search of the reviewed proteome, and it can
answer the question this task actually has for Foldseek — *is there
anything shaped like an IP3 receptor that the sequence instruments never
returned?* — while a UniProt50 sweep would mostly add unreviewed
near-duplicates of records the census already holds.

Which subset was searched is recorded in the output, so the negative is a
negative about a declared space (the roadmap's scoping rule) and not an
unbounded "nothing exists".

Hits are resolved back to UniProt taxonomy and scored against the census:
a hit the census already calls the family is a confirmation, a hit it calls
RYR is D14 arriving by a third instrument, and a hit it has never seen is
the interesting one — *if* it is a hit to the whole receptor.

**Coverage is computed here, not taken from Foldseek.** Foldseek reports
`qtmscore`/`ttmscore`, but on this data a 191-residue alignment against a
2,300-residue query still comes back at 0.72, so those numbers cannot be
read as whole-chain TM-scores the way TM-align's are. The coverage the
verdict rests on is therefore `alnlen` over the query chain's own length,
taken from the panel manifest. Without it every MIR-domain protein in
Swiss-Prot reads as a novel receptor: the sweep's above-bar non-census hits
are ~190-residue matches to SDF2, SDF2L1 and their plant and amoebozoan
relatives — the MIR domain the family shares, which is the one control
class the PDB could not fill.

**What this sweep structurally cannot find.** AlphaFold DB holds no
ryanodine receptor model at all (they are past its length ceiling), so a
RyR verdict here is unreachable by construction rather than absent by
evidence. The report says so.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import (load_census, protein_facts, results_dir,    # noqa: E402
                     structures_dir, write_tsv)

#: AFDB subsets Foldseek can build. The name is written into the output so
#: a claim can never outrun the space it was measured in.
SUBSETS = {
    "afdb_swissprot": "Alphafold/Swiss-Prot",
    "afdb_proteome": "Alphafold/Proteome",
    "afdb_uniprot50": "Alphafold/UniProt50",
}

FORMAT = ("query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
          "tend,evalue,bits,alntmscore,qtmscore,ttmscore,taxid,taxname")
COLUMNS = FORMAT.split(",")
OUT_COLUMNS = ["subset", "query", "target", "accession", "alntmscore",
               "qtmscore", "ttmscore", "fident", "alnlen", "query_length",
               "query_coverage", "evalue", "bits", "taxid", "taxname",
               "census_call", "census_group", "census_gene", "target_name",
               "target_gene", "target_pfams", "verdict"]

#: TM-align's published same-fold bar, reused so one number means one thing
#: across the task.
TM_FOLD = 0.50

#: An alignment covering less than this much of the query is a *domain*
#: match, however good its TM-score.
MIN_QUERY_COVERAGE = 0.50


def foldseek_bin() -> str:
    env = Path("/opt/anaconda3/envs/piezo1/bin/foldseek")
    return str(env) if env.exists() else "foldseek"


def db_dir() -> Path:
    d = structures_dir() / "foldseek_db"
    d.mkdir(parents=True, exist_ok=True)
    return d


def ensure_db(subset: str) -> Path:
    """Download the subset index if it is not already on the data root."""
    target = db_dir() / subset
    if (target.with_suffix(".dbtype")).exists() or target.exists():
        return target
    subprocess.run([foldseek_bin(), "databases", SUBSETS[subset], str(target),
                    str(db_dir() / "tmp")], check=True)
    return target


def search(query_pdbs: list[Path], subset: str, out_dir: Path,
           evalue: float = 1e-3, max_seqs: int = 2000,
           threads: int = 2) -> Path:
    out = out_dir / f"foldseek_{subset}.tsv"
    if out.exists() and out.stat().st_size > 0:
        return out
    qdir = out_dir / "queries"
    qdir.mkdir(parents=True, exist_ok=True)
    for p in query_pdbs:
        dest = qdir / p.name
        if not dest.exists():
            dest.write_bytes(p.read_bytes())
    db = ensure_db(subset)
    subprocess.run(
        [foldseek_bin(), "easy-search", str(qdir), str(db), str(out),
         str(out_dir / "tmp"), "--format-output", FORMAT,
         "-e", str(evalue), "--max-seqs", str(max_seqs),
         "--threads", str(threads), "--exhaustive-search", "0"],
        check=True)
    return out


def _accession(target: str) -> str:
    """AFDB entries are named `AF-<accession>-F1-model_v4.cif.gz`."""
    parts = target.split("-")
    return parts[1] if len(parts) > 2 and parts[0] == "AF" else target


def annotate(raw: Path, subset: str,
             query_lengths: dict[str, int] | None = None) -> list[dict]:
    census = {r["accession"]: r for r in load_census()}
    query_lengths = query_lengths or {}
    rows = []
    with open(raw) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(COLUMNS):
                continue
            rec = dict(zip(COLUMNS, parts))
            acc = _accession(rec["target"])
            crow = census.get(acc)
            call = crow.get("call", "") if crow else ""
            qlen = query_lengths.get(Path(rec["query"]).stem, 0)
            try:
                alnlen = int(rec["alnlen"])
            except ValueError:
                alnlen = 0
            qcov = round(alnlen / qlen, 4) if qlen else ""
            try:
                tm = float(rec["alntmscore"])
            except ValueError:
                tm = 0.0
            if call == "ITPR":
                verdict = "confirms census ITPR"
            elif call == "RYR":
                verdict = "sister family (D14)"
            elif tm < TM_FOLD:
                verdict = "below the same-fold bar"
            elif qcov != "" and qcov < MIN_QUERY_COVERAGE:
                verdict = "shared domain, not a receptor"
            else:
                verdict = "novel structural lead"
            rows.append({
                "subset": subset, "query": rec["query"], "target": rec["target"],
                "accession": acc, "alntmscore": rec["alntmscore"],
                "qtmscore": rec["qtmscore"], "ttmscore": rec["ttmscore"],
                "fident": rec["fident"], "alnlen": rec["alnlen"],
                "query_length": qlen or "", "query_coverage": qcov,
                "evalue": rec["evalue"], "bits": rec["bits"],
                "taxid": rec["taxid"], "taxname": rec["taxname"],
                "census_call": call or "not_in_census",
                "census_group": crow.get("group", "") if crow else "",
                "census_gene": crow.get("gene", "") if crow else "",
                "verdict": verdict,
            })
    rows.sort(key=lambda r: (r["query"], -float(r["alntmscore"] or 0)))
    # Name only the hits a claim could rest on. Looking up 555 sub-bar hits
    # would be 555 requests to describe noise; the ones above the fold bar
    # that the census has never held are the rows the report names, so
    # those are the ones whose paperwork is fetched.
    named: dict[str, dict] = {}
    for r in rows:
        if r["census_call"] != "not_in_census":
            continue
        if r["verdict"] in ("below the same-fold bar",):
            continue
        acc = r["accession"]
        if acc not in named:
            named[acc] = protein_facts(acc)
        r.update({"target_name": named[acc]["name"],
                  "target_gene": named[acc]["gene"],
                  "target_pfams": named[acc]["pfams"]})
    return rows


def write(rows: list[dict]) -> Path:
    path = results_dir() / "foldseek_hits.tsv"
    write_tsv(path, OUT_COLUMNS, [[r.get(c, "") for c in OUT_COLUMNS]
                                  for r in rows])
    return path


def run(manifest: list[dict], subset: str = "afdb_swissprot",
        threads: int = 2) -> list[dict]:
    queries = [Path(m["path"]) for m in manifest
               if m.get("status") == "ok" and m.get("role") == "reference"
               and m.get("call") == "ITPR"]
    if not queries:
        return []
    work = structures_dir() / "foldseek"
    work.mkdir(parents=True, exist_ok=True)
    raw = search(queries, subset, work, threads=threads)
    lengths = {Path(m["path"]).stem: int(m.get("resolved_residues") or 0)
               for m in manifest if m.get("path")}
    return annotate(raw, subset, lengths)
