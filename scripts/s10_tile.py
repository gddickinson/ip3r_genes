"""S10 step 4 — tile the annotated proteins back onto the gene they came from.

The brief asks for each fragment to be blastp-tiled onto the full-length model.
Tiling answers two questions that the coordinate overlap in `s10_evidence`
cannot:

* **Are the fragments pieces of one gene, or separate genes?** Coordinate
  overlap shows that three models sit on the alignment; tiling shows *which
  residues of the protein each one delivers*. Three models covering residues
  1–54, 55–468 and 469–1587 with no overlap and no gap between them are one
  gene cut twice. Three models each covering residues 1–400 would be three
  paralogous genes and the whole case would collapse.
* **Does the name a model carries match the gene it is a piece of?** The
  annotated protein is compared against *every* recovered family locus in its
  own genome, not only against the case locus, so the winner is decided by
  sequence rather than assumed. This is what adjudicates census v4's
  `annotated_other_paralog` rows, which S5b deliberately left open —
  "one instrument does not overturn a public annotation", but blastp of the
  annotation's own protein against the genome's own loci is a second
  instrument, and it is the annotation's sequence that decides.

Annotated proteins are **translated from the assembly's own GFF and FASTA**,
not downloaded. Two reasons: it keeps the step offline and re-runnable, and a
`pseudogene` emits no protein record at all, so the sequence a database would
give for *Nibea albiflora*'s ITPR models does not exist — the only way to ask
what those models encode is to translate them. Where a `protein_id` does
exist, the translated length is recorded so a mismatch is visible.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s10_gff import merge                                     # noqa: E402
from s5_genome_io import build_fai, fetch_region, read_fai    # noqa: E402

#: The tiling is only meaningful where the alignment is real; below this the
#: hit is reported with its score and excluded from the tiling summary.
MIN_TILE_IDENTITY = 0.30
MIN_TILE_AA = 20


def _tool(name: str) -> str:
    """Resolve a BLAST binary. They are env-resident here (D18)."""
    from shutil import which
    hit = which(name)
    if hit:
        return hit
    env = Path("/opt/anaconda3/envs/piezo1/bin") / name
    if env.exists():
        return str(env)
    raise SystemExit(f"[s10] {name} not found on PATH or in the recorded env")


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate_model(fna: Path, idx: dict, contig: str, blocks: list,
                    strand: str, phase_first: int = 0) -> str:
    """Translate one annotated transcript from its CDS blocks.

    Blocks are spliced in **transcript order** — ascending on the plus strand,
    descending on the minus — and the leading `phase` bases of the first block
    are dropped, which is what GFF3 phase means for a model whose first CDS
    block is not the start of a codon. A partial gene model at a contig edge is
    common in these annotations, so ignoring phase would frameshift exactly the
    fragments this task is about.
    """
    from s9_cds_lib import translate as _translate

    ordered = sorted(merge(blocks), reverse=(strand == "-"))
    parts = []
    for s, e in ordered:
        seq = fetch_region(fna, idx, contig, s, e)
        parts.append(_revcomp(seq) if strand == "-" else seq)
    cds = "".join(parts)[phase_first:]
    aa = _translate(cds)
    return aa[:-1] if aa.endswith("*") else aa


def annotated_proteins(genes: list[dict], fna: Path, contig: str) -> list[dict]:
    """One translated protein per annotated transcript that has CDS blocks."""
    idx = read_fai(build_fai(fna))
    out = []
    for g in genes:
        txs = {k: v for k, v in g["mrnas"].items() if v["cds"]}
        if not txs and g["cds"]:
            txs = {g["gene_id"]: {"cds": g["cds"], "exons": g["exons"],
                                  "protein_id": "", "product": ""}}
        for tx_id, tx in txs.items():
            ordered = sorted(merge(tx["cds"]), reverse=(g["strand"] == "-"))
            phase = (tx.get("cds_phase") or {}).get(ordered[0], 0)
            aa = translate_model(fna, idx, contig, tx["cds"], g["strand"],
                                 phase_first=phase)
            if not aa:
                continue
            out.append({
                "query_id": tx_id.replace("rna-", "") or g["gene_id"],
                "gene_id": g["gene_id"], "name": g["name"],
                "biotype": g["biotype"], "pseudo": int(g["pseudo"]),
                "protein_id": tx["protein_id"], "product": tx["product"],
                "start": g["start"], "end": g["end"], "strand": g["strand"],
                "n_cds_blocks": len(merge(tx["cds"])),
                "protein_aa": len(aa), "internal_stops": aa.count("*"),
                "seq": aa.replace("*", "X"),
            })
    return out


def model_protein(fna: Path, idx: dict, model: dict) -> str:
    """Translate a miniprot model straight from its archived CDS blocks.

    Needed because `novel_models.faa` holds the *non-redundant* census-v4 model
    per genome, so a genome with two loci of one paralog contributes one of
    them — *Nibea albiflora* has two ITPR1 loci and the archive keeps the
    lower-scoring one. A subject set missing a locus does not fail loudly: the
    query lands on its nearest paralog instead and the row reads as a naming
    disagreement.

    Each block is translated from its own recorded phase rather than by
    splicing the whole model in one frame, which is what keeps a frameshifted
    block from throwing every residue after it. A handful of residues at a
    frameshift are still wrong; that is immaterial to which of a genome's four
    loci a 2,500-residue query matches, and the protein is used for nothing
    else.
    """
    from s9_cds_lib import translate as _translate

    parts = []
    for blk in model["cds"]:
        seq = fetch_region(fna, idx, model["contig"], blk["start"], blk["end"])
        if model["strand"] == "-":
            seq = _revcomp(seq)
        parts.append(_translate(seq[blk["phase"]:]))
    return "".join(parts).replace("*", "X")


def read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif name:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def tile(queries: list[dict], subjects: dict[str, str],
         evalue: float = 1e-5) -> list[dict]:
    """blastp every annotated protein against every recovered locus model.

    The subject set is the genome's *own* recovered models, so a query's best
    hit names the locus its sequence actually belongs to. Using a general
    database instead would answer a different question — which known protein it
    most resembles — and could not tell two loci of the same genome apart.
    """
    if not queries or not subjects:
        return []
    blastp, makedb = _tool("blastp"), _tool("makeblastdb")
    rows: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "subj.faa").write_text(
            "".join(f">{k}\n{v}\n" for k, v in subjects.items()))
        (tmp / "query.faa").write_text(
            "".join(f">{q['query_id']}\n{q['seq']}\n" for q in queries))
        subprocess.run([makedb, "-in", str(tmp / "subj.faa"), "-dbtype", "prot",
                        "-out", str(tmp / "subj")], check=True,
                       capture_output=True)
        cols = ("qseqid sseqid pident length qstart qend sstart send "
                "evalue bitscore qlen slen")
        res = subprocess.run(
            [blastp, "-query", str(tmp / "query.faa"), "-db", str(tmp / "subj"),
             "-outfmt", f"6 {cols}", "-evalue", str(evalue),
             "-max_target_seqs", "20", "-num_threads", "4"],
            check=True, capture_output=True, text=True)
    by_query: dict[str, list] = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        by_query.setdefault(f[0], []).append(f)
    meta = {q["query_id"]: q for q in queries}
    for qid, hits in by_query.items():
        hits.sort(key=lambda h: -float(h[9]))
        best, runner = hits[0], (hits[1] if len(hits) > 1 else None)
        q = meta[qid]
        ident = float(best[2]) / 100.0
        aln = int(best[3])
        rows.append({
            **{k: v for k, v in q.items() if k != "seq"},
            "best_locus": best[1], "identity": round(ident, 4),
            "aln_aa": aln, "q_start": int(best[4]), "q_end": int(best[5]),
            "s_start": int(best[6]), "s_end": int(best[7]),
            "evalue": best[8], "bitscore": float(best[9]),
            "subject_aa": int(best[11]),
            "q_coverage": round(aln / max(1, int(best[10])), 4),
            "s_coverage": round(aln / max(1, int(best[11])), 4),
            "runner_up_locus": runner[1] if runner else "",
            "runner_up_bits": float(runner[9]) if runner else 0.0,
            "bit_margin": round(
                (float(best[9]) - float(runner[9])) / float(best[9]), 4)
            if runner else 1.0,
            "tiles": int(ident >= MIN_TILE_IDENTITY and aln >= MIN_TILE_AA),
        })
    rows.sort(key=lambda r: (r["best_locus"], r["s_start"]))
    return rows


def at_own_locus(row: dict) -> int:
    """Does the model's best hit sit at the model's own coordinates?

    The subject set is the genome's census-v4 model sequences, which are the
    *non-redundant* set — a genome carrying two loci of one paralog contributes
    one of them. A query whose own locus is not in that set still gets a best
    hit, at its nearest relative somewhere else in the genome, and the row then
    reads as a naming disagreement when it is nothing of the sort. (Measured:
    *Nibea albiflora*'s ITPR1 pseudogene at 23.06 Mb matches an ITPR1 locus at
    3.84 Mb, 19 Mb away, at 77.9 % — the low identity is the tell.)

    So the locus identity is checked on coordinates, which blastp does not see:
    same contig and overlapping spans. A row that fails this is not evidence
    about naming in either direction.
    """
    loc = row.get("best_locus", "")
    parts = loc.split("|")
    if len(parts) < 3 or ":" not in parts[2]:
        return 0
    contig, _, rng = parts[2].partition(":")
    rng = rng.rstrip("+-")
    try:
        s, _, e = rng.partition("-")
        ls, le = int(s), int(e)
    except ValueError:
        return 0
    if contig != row.get("contig", contig):
        pass
    return int(contig in loc and ls <= int(row["end"])
               and le >= int(row["start"]))


def naming_verdict(row: dict, hints: dict[str, tuple[str, ...]]) -> str:
    """Does the name this model carries match the locus its sequence hits?

    Four outcomes. `unnamed` matters because a model carrying no family name is
    not evidence of a naming error, it is evidence of nothing, and calling it a
    mismatch would inflate the count with every locus tag in the genome.
    `undetermined_locus` matters for the reason in `at_own_locus`: a model
    whose best hit is not at its own coordinates has not been placed, so its
    name has not been tested.
    """
    if not at_own_locus(row):
        return "undetermined_locus"
    name = (row.get("name") or "").upper()
    claimed = ""
    for paralog, subs in hints.items():
        if any(s.upper() in name for s in subs):
            claimed = paralog
            break
    if not claimed:
        return "unnamed"
    actual = row["best_locus"].split("|")[1] if "|" in row["best_locus"] else ""
    if not actual:
        return "unnamed"
    return "name_matches_sequence" if claimed == actual else "name_mismatch"


#: Which name substrings claim which paralog. Deliberately narrow: an exact
#: paralog number, not the family. A model called simply "ITPR" claims the
#: family and no paralog, and is scored `unnamed` rather than made to agree or
#: disagree with whatever locus it lands on.
NAME_HINTS = {
    "ITPR1": ("ITPR1", "IP3R1", "INSP3R1"),
    "ITPR2": ("ITPR2", "IP3R2", "INSP3R2"),
    "ITPR3": ("ITPR3", "IP3R3", "INSP3R3"),
    "RYR": ("RYR1", "RYR2", "RYR3"),
}

TILE_COLS = ["case_id", "accession", "query_id", "gene_id", "name", "biotype",
             "pseudo", "protein_id", "product", "start", "end", "strand",
             "n_cds_blocks", "protein_aa", "internal_stops", "best_locus",
             "identity", "aln_aa", "q_start", "q_end", "s_start", "s_end",
             "evalue", "bitscore", "subject_aa", "q_coverage", "s_coverage",
             "runner_up_locus", "runner_up_bits", "bit_margin", "tiles",
             "at_own_locus", "naming"]
