"""S21's data layer — the loci, the exon blocks, and the joins.

S21 measures the *gene* where every task before it measured the search, and
almost everything it needs already exists on the drive:

* the sweep's retained `miniprot.gff` per genome carries every alignment's CDS
  blocks with each block's genomic span, its span in the bait's own residue
  coordinates and its phase — the exon/intron structure, already spliced;
* the genome FASTAs the sweep fetched are still there, so every intron's
  splice dinucleotides can be read rather than assumed; and
* S16's `loci.tsv` is the committed per-locus scope, with the coverage and
  identity bars, D4's contiguity flag and S15a's ORF verdict already joined.

Three things are deliberately *not* re-derived here.

* The loci come from the archived per-genome `summary.json` and **not** the
  ledger, which holds one row per genome × cell and therefore only the best
  locus (S8, S15a, S16 and S18 each hit this).
* The locus scope — which loci count as a gene copy — is S16's committed
  `is_copy`, read back rather than recomputed, because two tasks that decide
  independently what a copy is will disagree.
* The alignment frame is S6's committed `aln.fasta`, whose tips include all
  three human paralogues and human RYR2, so a boundary in one paralogue's
  residue numbering reaches another paralogue's through a column and not
  through a second alignment built here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s15_lib                                                    # noqa: E402
import s18_lib                                                    # noqa: E402

PROJECT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT / "results"
OUT = RESULTS / "gene_architecture"
FIGS = OUT / "figures"

S16_LOCI = RESULTS / "duplication" / "loci.tsv"
TELEOST = RESULTS / "duplication" / "teleost_copies.tsv"
CONTIG_CELLS = RESULTS / "methods" / "contiguity_cells.tsv"
LOCUS_AUDIT = RESULTS / "annotation_audit" / "locus_audit.tsv"
MSA = RESULTS / "msa_v2" / "aln.fasta"
BAITS = RESULTS / "s5_baits" / "baits.faa"
BAIT_MANIFEST = RESULTS / "s5_baits" / "bait_manifest.tsv"
MANIFEST = RESULTS / "genome_manifest.tsv"
INTRON_CAL = RESULTS / "s5_baits" / "intron_calibration.tsv"

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
CONTROL_CELL = "RYR"
CELLS = PARALOGS + (CONTROL_CELL,)

#: The residue frame every boundary is expressed in, one per cell. All four are
#: tips of S6's committed alignment, which is what makes a boundary in an
#: ITPR2 locus comparable with one in an ITPR1 locus: the two reach the same
#: column. Human, because that is the accession every other task's reference
#: numbering is already in (S0's structural measurements, S17's variants).
FRAME_ACC = {"ITPR1": "Q14643", "ITPR2": "Q14571", "ITPR3": "Q14573",
             "RYR": "Q92736"}

read_tsv = s15_lib.read_tsv
write_tsv = s18_lib.write_tsv
sha256 = s15_lib.sha256
data_root = s15_lib.data_root


def sweep_dir(acc: str | None = None) -> Path:
    """The sweep's evidence directory, per genome when given an accession."""
    root = s15_lib.sweep_dir()
    return root / acc if acc else root

load_summaries = s15_lib.load_summaries
iter_fasta = s18_lib.iter_fasta
bait_labels = s18_lib.bait_labels
tool_bin = s18_lib.tool_bin


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    (RESULTS / "session_live.json").write_text(json.dumps({
        "task": "S21", "stage": stage,
        "steps": [{"label": a, "done": b} for a, b in steps]}, indent=1))


def log(msg: str) -> None:
    print(f"  {msg}", flush=True)


# --------------------------------------------------------------------------
# scope
# --------------------------------------------------------------------------
INT_COLS = ("locus_idx", "start", "end", "aligned_aa", "q_start", "q_end",
            "bait_len", "frameshifts", "stop_codons", "lesions",
            "merged_models", "contig_n50")
FLOAT_COLS = ("identity", "coverage", "lesion_density", "family_margin",
              "paralog_margin", "annot_frac_cds")
BOOL_COLS = ("contig_spans_gene", "is_copy", "contig_edge",
             "annot_paralog_matches")


def s16_loci() -> list[dict]:
    """S16's committed per-locus table, typed. The scope, not a re-derivation.

    `is_copy` is S16's bar (coverage >= 0.50, identity >= the sweep floor,
    >= 500 aligned residues) and is read back rather than recomputed. S21 adds
    its own, higher coverage bar on top — a model covering half its bait has
    half the exons, so a copy is not automatically an architecture — and that
    bar is measured in `s21_architecture`, never inherited silently.
    """
    rows = read_tsv(S16_LOCI)
    for r in rows:
        for k in INT_COLS:
            r[k] = int(float(r.get(k) or 0))
        for k in FLOAT_COLS:
            r[k] = float(r.get(k) or 0.0)
        for k in BOOL_COLS:
            r[k] = str(r.get(k)) == "True"
    return rows


def contiguity_index() -> dict[tuple[str, str], dict]:
    """S19's per-cell control flag and contiguity — a join, not a rebuild."""
    return {(r["accession"], r["cell"]): r for r in read_tsv(CONTIG_CELLS)}


def audit_index() -> dict[tuple[str, str, int], dict]:
    """S18's annotation state per locus — where the fragment test starts."""
    return {(r["accession"], r["cell"], int(float(r["locus_idx"]))): r
            for r in read_tsv(LOCUS_AUDIT)}


def genome_fna(acc: str) -> Path | None:
    d = data_root() / "genomes" / acc
    if not d.exists():
        return None
    for p in sorted(d.glob("*_genomic.fna")):
        return p
    for p in sorted(d.rglob("*.fna")):
        return p
    return None


def gff_path(acc: str) -> Path | None:
    return s18_lib.gff_path(acc)


# --------------------------------------------------------------------------
# sequences and the alignment frame
# --------------------------------------------------------------------------
def bait_seqs() -> dict[str, str]:
    """The committed bait panel, keyed by accession (the panel's stable id)."""
    out = {}
    for name, seq in iter_fasta(BAITS):
        out[name.split("|")[0]] = seq.replace("-", "").upper()
    return out


def msa_rows() -> dict[str, str]:
    """S6's aligned rows, keyed by the accession at the end of each tip label."""
    out = {}
    for name, seq in iter_fasta(MSA):
        out[name.split("_")[-1]] = seq.upper()
    return out


quantile = s18_lib.quantile


# --------------------------------------------------------------------------
# genome sequence — the independent axis every boundary claim is read against
# --------------------------------------------------------------------------
import s10_evidence as EV                                         # noqa: E402
from s5_genome_io import build_fai, read_fai                       # noqa: E402

#: The splice rule is S10's, imported rather than restated: `splice_pair`
#: builds a minus-strand intron's pair from the reverse complement of both
#: ends and swaps them, and getting that wrong turns every canonical
#: minus-strand intron into `CT..AC` — a systematic non-canonical signal that
#: would discredit exactly the alignments S21 is measuring.
splice_pair = EV.splice_pair
splice_class = EV.splice_class
SPLICE_CANONICAL = EV.SPLICE_CANONICAL

_FAI: dict[str, dict] = {}


def fai_of(acc: str, fna: Path) -> dict:
    """The genome's offset index, built once per process."""
    if acc not in _FAI:
        _FAI[acc] = read_fai(build_fai(fna))
    return _FAI[acc]


class LocusRegion(EV.Region):
    """`s10_evidence.Region` over exactly one locus, with the index cached.

    S10 fetched a locus plus 200 kb of flank because it was asking what else
    the annotation held nearby; S21 asks only about intron edges, which are
    inside the model's own span by construction, so the fetch is the span and
    the index is read once per genome instead of once per locus.
    """

    def __init__(self, acc: str, fna: Path, contig: str, start: int, end: int):
        from s5_genome_io import fetch_region
        idx = fai_of(acc, fna)
        self.contig, self.start = contig, max(1, start)
        self.seq = fetch_region(fna, idx, contig, self.start, end).upper()
        self.end = self.start + len(self.seq) - 1


def loci_with_mp(summaries: dict | None = None) -> list[dict]:
    """S16's committed loci with the sweep's `mp_id` joined onto each.

    S16's table is the scope but carries no miniprot id, and the id is what
    reaches the CDS blocks. The join is on (accession, cell, locus_idx) —
    S16's index is the enumeration of the summary's own locus list, as S18's
    is — and the contig, start and end are **required to match**: an index
    that had drifted would attach one gene's exon structure to another gene's
    row, which is invisible in every downstream number.
    """
    summaries = summaries if summaries is not None else load_summaries()
    rows = s16_loci()
    n_bad = 0
    for r in rows:
        s = summaries.get(r["accession"]) or {}
        cells = s.get("cells") or {}
        # S16 files the two unassignable `vertebrate_basal` loci out of the
        # summary's `other_loci` list rather than a cell, so both places are
        # looked in — the guard below is what found that.
        loci = ((cells.get(r["cell"]) or {}).get("loci")
                if r["cell"] in cells else s.get("other_loci")) or []
        if r["locus_idx"] >= len(loci):
            r["mp_id"] = ""
            n_bad += 1
            continue
        loc = loci[r["locus_idx"]]
        if (loc.get("contig") != r["contig"]
                or int(loc.get("start") or 0) != r["start"]):
            n_bad += 1
            r["mp_id"] = ""
            continue
        r["mp_id"] = loc.get("mp_id", "")
    if n_bad:
        raise SystemExit(f"[s21] {n_bad} S16 loci do not match the summaries — "
                         "the locus index has drifted; refusing to measure")
    return rows
