"""S9 — exact CDS reconstruction for the S5 miniprot genome gene models.

The S5 sweep ran `miniprot --gff --trans`, whose `##STA` protein — the
sequence that went into S6's alignment — is *frameshift-skipping*:
residues around a frameshift are dropped. The GFF CDS blocks therefore do
not splice back into that protein, and naive splicing drifts out of frame
partway through. Splicing the sweep GFF is not offered here for that
reason, rather than offered with a warning.

Rerunning miniprot on the locus with `--aln` adds the base-level
alignment rows:

    ##ATN   aligned genomic nucleotides (exons upper, introns lower)
    ##ATA   its translation, one letter per codon start

`##ATA` keeps the frameshift residues that `##STA` drops, so the mapping
is: codon per `##ATA` residue (the next three exonic nucleotides), then
`##STA` -> `##ATA` by difflib, emitting NNN wherever a `##STA` residue has
no clean `##ATA` counterpart. Every reconstruction is accepted only if the
locus rerun reproduces the S5 protein *exactly*; a near-match is a
different gene model and is refused.
"""

from __future__ import annotations

import bisect
import re
import subprocess
import sys
import tempfile
from difflib import SequenceMatcher
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s5_calibration import max_intron_for  # noqa: E402
from s5_genome_io import build_fai, fetch_region, read_fai  # noqa: E402

PAD_TRIES = (20_000, 5_000, 60_000)


# ---- genome-wide GFF: model lookup ----------------------------------------

def index_models(gff: Path) -> list[dict]:
    """Every miniprot model in a sweep GFF: bait, protein, span."""
    models: list[dict] = []
    paf: list[str] = []
    sta = ""
    with open(gff) as fh:
        for line in fh:
            if line.startswith("##PAF"):
                paf = line.rstrip("\n").split("\t")
                sta = ""
            elif line.startswith("##STA"):
                sta = line.rstrip("\n").split("\t")[1]
            elif "\tmRNA\t" in line:
                f = line.rstrip("\n").split("\t")
                mid = re.search(r"ID=([^;]+)", f[8])
                models.append({
                    "id": mid.group(1) if mid else "",
                    "bait": paf[1] if len(paf) > 1 else "",
                    "protein": sta,
                    "contig": f[0], "start": int(f[3]), "end": int(f[4]),
                    "strand": f[6],
                })
    return models


def find_model(models: list[dict], expected_prot: str,
               locus: tuple[str, int, int, str] | None = None) -> dict | None:
    """The model whose protein is `expected_prot` (msa `X` <-> miniprot `*`).

    `locus` narrows the search first. That matters here and did not in the
    PIEZO port: this sweep reports up to four ITPR loci per genome plus the
    RyR control, and two paralogs of the same length in one assembly would
    make a length-only fallback pick the wrong gene.
    """
    pool = models
    if locus is not None:
        contig, start, end, strand = locus
        near = [m for m in models
                if m["contig"] == contig and m["strand"] == strand
                and m["start"] <= end and m["end"] >= start]
        if near:
            pool = near
    for m in pool:
        if m["protein"].replace("*", "X") == expected_prot:
            return m
    same_len = [m for m in pool if len(m["protein"]) == len(expected_prot)]
    return same_len[0] if len(same_len) == 1 else None


# ---- locus rerun with --aln -------------------------------------------------

def _run_miniprot_aln(region_fa: Path, bait_faa: Path, max_intron: int,
                      threads: int = 8) -> list[dict]:
    cmd = ["miniprot", "-t", str(threads), "--gff", "--trans", "--aln",
           "--outs=0.3", "-N", "60", "-G", str(max_intron),
           str(region_fa), str(bait_faa)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"miniprot --aln failed: {proc.stderr[:300]}")
    out: list[dict] = []
    cur: dict = {}
    for line in proc.stdout.split("\n"):
        if line.startswith("##PAF"):
            if cur.get("sta"):
                out.append(cur)
            cur = {}
        elif line.startswith("##STA"):
            cur["sta"] = line.split("\t")[1]
        elif line.startswith("##ATN"):
            cur["atn"] = line.split("\t")[1]
        elif line.startswith("##ATA"):
            cur["ata"] = line.split("\t")[1]
    if cur.get("sta"):
        out.append(cur)
    return out


def _codons_for_ata(atn: str, ata: str) -> tuple[str, list[str]]:
    """(##ATA protein, one codon per residue) — codon = next 3 exonic nt."""
    up = [i for i, c in enumerate(atn) if c.isupper()]
    prot: list[str] = []
    codons: list[str] = []
    for i, c in enumerate(ata):
        if not (c.isalpha() or c == "*"):
            continue
        k = bisect.bisect_left(up, i)
        codon = ("".join(atn[j] for j in up[k:k + 3])
                 if k + 3 <= len(up) else "NNN")
        prot.append(c)
        codons.append(codon if len(codon) == 3 else "NNN")
    return "".join(prot), codons


def _map_sta_to_codons(sta: str, ata_prot: str, ata_codons: list[str]) -> str:
    """CDS for the ##STA protein; NNN where ##STA has no ##ATA counterpart."""
    out = ["NNN"] * len(sta)
    for tag, i1, i2, j1, j2 in SequenceMatcher(
            None, ata_prot, sta, autojunk=False).get_opcodes():
        if tag != "equal":
            continue
        for off in range(i2 - i1):
            out[j1 + off] = ata_codons[i1 + off]
    return "".join(out)


#: A locus rerun may place a *terminal* exon differently from the
#: whole-genome sweep — the k-mer index of an 800 kb region is not the
#: index of a 3.2 Gbp genome, so the seeds available at the ends differ.
#: An alignment of the same length that agrees with the sweep's protein
#: everywhere but a short run is the same gene model; the residues it
#: disagrees on are masked downstream rather than trusted. Anything
#: beyond this is a different gene model and is refused.
MAX_RERUN_DIFF_FRAC = 0.01


def _best_alignment(alns: list[dict], expected_prot: str
                    ) -> tuple[dict, int] | None:
    """The rerun alignment that reproduces the sweep's protein, and by how
    much it misses. Exact first; then same-length within the diff bound."""
    usable = [a for a in alns if "atn" in a and "ata" in a and a.get("sta")]
    for a in usable:
        if a["sta"].replace("*", "X") == expected_prot:
            return a, 0
    scored: list[tuple[int, dict]] = []
    for a in usable:
        sta = a["sta"].replace("*", "X")
        if len(sta) != len(expected_prot):
            continue
        ndiff = sum(1 for x, y in zip(sta, expected_prot) if x != y)
        if ndiff <= MAX_RERUN_DIFF_FRAC * len(expected_prot):
            scored.append((ndiff, a))
    if not scored:
        return None
    ndiff, a = min(scored, key=lambda t: t[0])
    return a, ndiff


def cds_for_model(assembly: str, sweep_root: Path, genome_fna: Path,
                  baits: dict[str, str], expected_prot: str,
                  threads: int = 8,
                  locus: tuple[str, int, int, str] | None = None,
                  ) -> tuple[str | None, str]:
    """Reconstruct the CDS of the S5 model whose protein is `expected_prot`."""
    gff = sweep_root / assembly / "miniprot.gff"
    if not gff.exists():
        return None, f"no_gff:{assembly}"
    model = find_model(index_models(gff), expected_prot, locus)
    if model is None:
        return None, f"model_not_matched:{assembly}"
    bait_seq = baits.get(model["bait"])
    if bait_seq is None:
        return None, f"bait_missing:{model['bait']}"

    fai = read_fai(build_fai(genome_fna))
    genome_bp = sum(v[0] for v in fai.values())
    max_intron = max_intron_for(genome_bp)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "bait.faa").write_text(f">{model['bait']}\n{bait_seq}\n")
        for pad in PAD_TRIES:
            seq = fetch_region(genome_fna, fai, model["contig"],
                               model["start"] - pad, model["end"] + pad)
            if not seq:
                continue
            (tmp / "region.fa").write_text(">region\n" + seq + "\n")
            try:
                alns = _run_miniprot_aln(tmp / "region.fa", tmp / "bait.faa",
                                         max_intron, threads)
            except RuntimeError as exc:
                return None, str(exc)
            best = _best_alignment(alns, expected_prot)
            if best is None:
                continue
            a, ndiff = best
            ata_prot, ata_codons = _codons_for_ata(a["atn"], a["ata"])
            cds = _map_sta_to_codons(a["sta"], ata_prot, ata_codons)
            loc = (f"{model['contig']}:{model['start']}-{model['end']}"
                   f"{model['strand']}")
            tag = "miniprot_aln" if ndiff == 0 else f"miniprot_aln_d{ndiff}"
            return cds, (f"{tag}:{assembly}:{model['id']}:{loc}"
                         f":pad{pad}:bait={model['bait'].split('|')[0]}")
    return None, f"aln_no_matching_protein:{assembly}:{model['id']}"
