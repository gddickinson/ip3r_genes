"""S7 — reciprocal-best-hit verification of the naming calls the tree
contradicts.

The brief: *RBH-verify any naming call the tree contradicts.* A tip that
the tree moves out of the clade its census label claims is either a
mis-annotated record or a real tree error, and neither the alignment nor
the model can tell those apart — both produced the placement in question.
RBH can, because it knows nothing about this alignment, this model or
this tree:

  forward    the candidate protein -> blastp against the **human
             reference proteome** (UP000005640): the best hit's gene name
             says which paralog it is.
  reciprocal that human protein -> blastp against every sequence in this
             project's own vertebrate reference-proteome DB restricted to
             the candidate's species: the candidate (or another isoform
             of the same gene) must come back best.

A call is `confirmed` only when both directions agree. Where the
reciprocal step has no proteome to run against — a tip that is a genomic
gene model rather than a UniProt record — the row says `forward_only`
rather than borrowing the forward hit's confidence.

Targets are read from `naming_conflicts.tsv`, never listed here, so the
set verified is exactly the set the tree disputed — **every** row of it,
`reassigned` and `unconstrained` alike. The first version filtered to
`reassigned`, of which this tree has none, so the step reported "the
tree contradicts no census label" while the table it had just read held
five tips the tree declined to leave in the clade their name claims.
A verification that cannot fire is not a verification.

What `agrees` compares therefore depends on what the tree offered, and
the row records which comparison was made in `compared_against`:

  reassigned    the tree names a different paralog — RBH is asked
                whether that reassignment is right, so the forward hit
                is compared against `tree_says`.
  unconstrained the tree offers no alternative, only a refusal to place
                the tip — so RBH is asked the other question, whether an
                instrument that has never seen this alignment upholds
                the census name. The comparison is against
                `census_group`, and an agreement here means the label is
                sound and the placement is the tree's own uncertainty.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s7_rbh.py
Out:  results/phylogeny/naming_rbh.tsv
Work: <data_root>/blast_db/s7_rbh/
"""

from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.s6_lib import read_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402
from scripts.s7_lib import MSA_DIR, PHYLO_DIR                # noqa: E402
from src.utils.data_root import get_data_root                 # noqa: E402

DATA = get_data_root()
HUMAN_GZ = DATA / "proteomes" / "vertebrata" / "UP000005640_9606.fasta.gz"
VERT_DB = DATA / "proteomes" / "vertebrata_refprot.fasta"
WORK = DATA / "blast_db" / "s7_rbh"
REPS_FASTA = MSA_DIR / "representatives.fasta"
CONFLICTS = PHYLO_DIR / "naming_conflicts.tsv"
OUT = PHYLO_DIR / "naming_rbh.tsv"

BLASTP = shutil.which("blastp") or \
    "/opt/anaconda3/envs/piezo1/bin/blastp"
MAKEDB = shutil.which("makeblastdb") or \
    "/opt/anaconda3/envs/piezo1/bin/makeblastdb"

GENE = re.compile(r"\bGN=(\S+)")
PARALOG = re.compile(r"\b(ITPR[123]|RYR[123])\b", re.I)


def gene_of(header: str) -> str:
    m = GENE.search(header)
    return m.group(1).upper() if m else ""


def paralog_of(header: str) -> str:
    g = gene_of(header)
    m = PARALOG.search(g) or PARALOG.search(header)
    return m.group(1).upper() if m else ""


def ensure_db(fasta: Path, name: str) -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    db = WORK / name
    if (WORK / f"{name}.pin").exists() or (WORK / f"{name}.pdb").exists():
        return db
    subprocess.run([MAKEDB, "-in", str(fasta), "-dbtype", "prot",
                    "-out", str(db)], check=True,
                   stdout=subprocess.DEVNULL)
    return db


def human_fasta() -> Path:
    out = WORK / "human.faa"
    if out.exists():
        return out
    WORK.mkdir(parents=True, exist_ok=True)
    with gzip.open(HUMAN_GZ, "rt") as fh, open(out, "w") as o:
        shutil.copyfileobj(fh, o)
    return out


def species_fasta(species: str) -> Path | None:
    """Every vertebrate-DB sequence whose header names this species.

    One streaming pass over the concatenated reference-proteome DB; the
    `OS=` field is what HMMER and BLAST both preserve, and matching on it
    is exact in this direction (unlike proteome attribution, S20's trap).
    """
    key = species.split(" (")[0].strip()
    if not key or not VERT_DB.exists():
        return None
    out = WORK / ("sp_" + re.sub(r"\W+", "_", key) + ".faa")
    if out.exists():
        return out if out.stat().st_size else None
    WORK.mkdir(parents=True, exist_ok=True)
    keep, n = False, 0
    with open(VERT_DB) as fh, open(out, "w") as o:
        for line in fh:
            if line.startswith(">"):
                keep = f"OS={key}" in line
                n += keep
            if keep:
                o.write(line)
    return out if n else None


def blast_best(query: Path, db: Path) -> tuple[str, float, float]:
    """Best hit: (subject title, bitscore, percent identity)."""
    r = subprocess.run(
        [BLASTP, "-query", str(query), "-db", str(db), "-max_target_seqs",
         "5", "-evalue", "1e-5", "-num_threads", "4",
         "-outfmt", "6 sseqid stitle bitscore pident"],
        capture_output=True, text=True, check=True)
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 4:
            return f[1], float(f[2]), float(f[3])
    return "", 0.0, 0.0


COLS = ["label", "census_group", "tree_rule", "tree_says", "forward_hit",
        "forward_paralog", "forward_pident", "reciprocal_ok", "rbh_call",
        "compared_against", "agrees"]


def _comparison(conflict: dict) -> tuple[str, str]:
    """Which claim RBH is being asked to check, and its expected answer.

    A `reassigned` tip has a tree-proposed paralog to test; an
    `unconstrained` one does not, and testing it against an empty
    `assigned_to` would score every such row `no` for the trivial reason
    that nothing was proposed.
    """
    if conflict["rule"] == "reassigned" and conflict["assigned_to"]:
        return "tree_says", conflict["assigned_to"]
    return "census_group", conflict["census_group"]


def main() -> int:
    conflicts = read_tsv(CONFLICTS) if CONFLICTS.exists() else []
    if not conflicts:
        write_tsv(OUT, COLS, [])
        print("no naming conflicts — the tree contradicts no census "
              f"label; wrote an empty {OUT.name}")
        return 0

    seqs = read_fasta(REPS_FASTA)
    reps = {r["label"]: r for r in read_tsv(MSA_DIR / "representatives.tsv")}
    hdb = ensure_db(human_fasta(), "human")
    rows = []
    for c in conflicts:
        lab = c["label"]
        seq = seqs.get(lab, "").replace("-", "")
        if not seq:
            rows.append({"label": lab, "census_group": c["census_group"],
                         "tree_rule": c["rule"],
                         "tree_says": c["assigned_to"],
                         "compared_against": _comparison(c)[0],
                         "rbh_call": "no_sequence", "agrees": "—"})
            continue
        q = WORK / "q.faa"
        write_fasta({lab: seq}, q)
        title, bits, pid = blast_best(q, hdb)
        fwd = paralog_of(title)

        recip, call = "", "forward_only"
        sp = reps.get(lab, {}).get("species", "")
        spf = species_fasta(sp)
        if spf and title:
            acc = title.split("|")[1] if "|" in title else ""
            hseqs = read_fasta(human_fasta())
            hkey = next((k for k in hseqs if acc and acc in k), "")
            if hkey:
                hq = WORK / "h.faa"
                write_fasta({hkey.split()[0]: hseqs[hkey]}, hq)
                sdb = ensure_db(spf, "sp_" + re.sub(r"\W+", "_", sp)[:40])
                back, _, _ = blast_best(hq, sdb)
                same = (lab.split("_")[-1] in back) or \
                    (paralog_of(back) and paralog_of(back) == fwd)
                recip = "yes" if same else "no"
                call = "confirmed" if same else "one_way"
        against, want = _comparison(c)
        rows.append({
            "label": lab, "census_group": c["census_group"],
            "tree_rule": c["rule"], "tree_says": c["assigned_to"],
            "forward_hit": (title.split(" OS=")[0])[:60],
            "forward_paralog": fwd, "forward_pident": f"{pid:.1f}",
            "reciprocal_ok": recip or "—", "rbh_call": call,
            "compared_against": against,
            "agrees": "yes" if fwd and fwd == want else "no",
        })
        print(f"{lab[:44]:44s} {c['rule']:13s} "
              f"vs {against}={want or '?':6s} rbh={fwd or '?':6s} {call}")

    write_tsv(OUT, COLS, rows)
    agree = sum(1 for r in rows if r.get("agrees") == "yes")
    print(f"\n{agree} of {len(rows)} disputed tips agree with RBH "
          f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
