"""s5_sweep_lib.py — baits, miniprot execution/parsing and locus clustering.

The alignment layer of the S5 genomic sweep. Rescue (tblastn) is in
`s5_rescue.py`; per-cell status classification is in `s5_classify.py`.

Cells are the three vertebrate paralogs plus **RYR as an internal positive
control** (the brief, and D14): a genome where the RyR control finds nothing
has an assembly or pipeline problem rather than a biological result, and the
RyR baits are in the same miniprot run for a second reason — without them a
genome's own RyR locus is claimed by an ITPR bait at partial coverage and
reads out of the ledger as an ITPR `fragment`.

Statuses per genome x class cell (assigned in s5_classify):
  found_annotated      coverage >= COV_FOUND, locus overlaps an annotated gene
  found_unannotated    coverage >= COV_FOUND, annotation exists but no gene there
  found_no_annotation  coverage >= COV_FOUND, assembly carries no annotation
  fragment             best locus coverage < COV_FOUND
  assembly_gap         fragment whose locus abuts a contig edge or an N-gap
  tblastn_trace        no miniprot locus; tblastn rescue HSP at E <= RESCUE_E
  absent               no miniprot locus and no rescue hit
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s5_genome_io import build_fai, read_fai  # noqa: E402
from s5_calibration import (DEFAULT_MAX_INTRON,  # noqa: F401,E402
                            intron_rule, intron_rule_capped,
                            itpr_span_stats, max_intron_for, spans_a_gene)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
S5_BAITS = PROJECT_ROOT / "results" / "s5_baits"

#: The ledger's ITPR cells, in the order figures and tables use them.
CLASSES = ("ITPR1", "ITPR2", "ITPR3")

#: The internal positive-control cell. Classified like a paralog cell but
#: reported as a control: presence, not which RyR.
CONTROL_CLASS = "RYR"

#: Bait clades that carry no paralog assignment. A locus topped by one of
#: these is offered to whichever ITPR paralog scores highest there (see
#: `s5_classify.classify_class`) — but only within its own family.
UNRESOLVED_ITPR_CLADES = ("vertebrate_basal",)

#: Minimum identity for a cluster of alignments to count as a locus at all.
#:
#: **Measured, with a wide empty gap.** Across the 309-genome sweep, the 571
#: loci whose paralog an assembly's own annotation independently confirms have
#: a minimum identity of **0.759** — not one below it. The population it
#: excludes sits at 0.23-0.34: the shared channel module matching unrelated
#: proteins (PSME4, TRAPPC9, CDH20, EXOC1 ...), chained into apparent loci by
#: a large `-G`. Any floor from 0.35 to 0.50 drops exactly those 22 and zero
#: confirmed loci, so this is a gap rather than a tuned cut; 0.40 sits in the
#: middle of it, well above the ~0.25 cross-family noise S1 measured between
#: ITPR and RyR and well below the 0.759 floor of real loci.
#:
#: It matters most in the giant genomes, where `-G` is 2 Mbp: 51 % of loci in
#: assemblies over 5 Gbp were under 30 % coverage, against 9 % elsewhere. The
#: *status* was never wrong — the best locus wins and the real gene always
#: scored best — but `n_loci` was inflated, and copy number is a result.
MIN_LOCUS_IDENTITY = 0.40

COV_FOUND = 0.70          # bait coverage for a "found" call
LOCUS_GAP = 10_000        # merge alignments this close on one strand
EDGE_BP = 10_000          # locus within this of a contig end -> edge
NGAP_RUN = 100            # N-run length that counts as an assembly gap
RESCUE_E = 1e-5


# ---------------------------------------------------------------- baits

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs: dict[str, str], path: Path) -> Path:
    with open(path, "w") as out:
        for k, v in seqs.items():
            out.write(f">{k}\n")
            for i in range(0, len(v), 60):
                out.write(v[i:i + 60] + "\n")
    return path


def load_baits(baits_dir: Path | None = None
               ) -> tuple[dict[str, str], dict[str, dict]]:
    """(id -> sequence, id -> {clade, family, band, length}).

    Only a *committed, screened* panel: the builders enforce the selection
    rules, so falling back to a raw seed set would quietly run the sweep on an
    unscreened panel. `baits_dir` selects which panel — S5's vertebrate one by
    default, S23's non-vertebrate one when that sweep passes its directory.
    """
    root = baits_dir or S5_BAITS
    faa, manifest = root / "baits.faa", root / "bait_manifest.tsv"
    if not faa.exists() or not manifest.exists():
        raise SystemExit(f"bait panel missing under {root}\n"
                         "  run the panel builder for that task first")
    seqs = read_fasta(faa)
    meta: dict[str, dict] = {}
    with open(manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["id"]] = {"clade": row["clade"], "family": row["family"],
                               "band": row.get("band", row["clade"]),
                               "paralog": row["paralog"],
                               "length": int(row["length"] or 0)}
    unknown = set(seqs) - set(meta)
    if unknown:
        raise SystemExit(f"baits.faa carries {len(unknown)} records the "
                         f"manifest does not: {sorted(unknown)[:3]}")
    for k, s in seqs.items():
        meta[k]["length"] = meta[k]["length"] or len(s)
    return seqs, meta


# ---------------------------------------------------------------- miniprot

# Per-chunk reference size for giant genomes. miniprot has no reference-index
# batching, so a 40 Gbp assembly cannot be indexed in one pass.
#
# **Measured, not inherited.** The PIEZO port used 4 Gbp on the strength of a
# 32 GB machine. `s5_test_chunked.py` measured miniprot's actual appetite on
# this project's panel: 11.76 GB peak RSS over a 1.86 Gbp reference at 2
# threads, i.e. **6.32 GB per Gbp**, which projects a 4 Gbp chunk to ~25 GB.
# That fits the 34 GB machine only just, and the giants run unattended. 2.5
# Gbp projects to ~16 GB and leaves real headroom; the cost is more chunks
# (Protopterus 40 Gbp -> ~16 rather than 10), and chunk count is cheap because
# a locus never spans a chunk boundary — chunks hold whole contigs.
#
# Note the cap cannot bind below one contig: a single 2 Gbp lungfish
# chromosome forms a 2 Gbp chunk whatever this says, because splitting a
# contig would require remapping every coordinate.
CHUNK_BP = 2_500_000_000

#: Measured RSS per Gbp of reference, for the headroom check above.
MINIPROT_GB_RSS_PER_GBP = 6.32

def run_miniprot(fna: Path, baits_faa: Path, out_gff: Path, threads: int = 8,
                 max_intron: int = DEFAULT_MAX_INTRON) -> None:
    """Align the panel to one reference, writing `out_gff` **atomically**.

    The output goes to a `.partial` sibling and is renamed only after miniprot
    exits 0. This matters because the driver reuses any existing non-empty
    GFF: a run killed mid-write — an interrupted sweep, an OOM during the
    chunked path, a machine restart — would otherwise leave a truncated file
    that the next run silently accepts as complete, producing a genome whose
    ledger row looks clean and is short of loci. A partial file is invisible
    after the fact, because a GFF cut at a line boundary parses perfectly.
    """
    tmp = out_gff.with_suffix(out_gff.suffix + ".partial")
    cmd = ["miniprot", "-t", str(threads), "--gff", "--trans",
           "--outs=0.3", "-N", "60", "-G", str(max_intron),
           str(fna), str(baits_faa)]
    try:
        with open(tmp, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                  text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"miniprot failed: {proc.stderr.strip()[:300]}")
        tmp.replace(out_gff)
    finally:
        tmp.unlink(missing_ok=True)


def split_fasta_by_contig(fna: Path, out_dir: Path,
                          chunk_bp: int = CHUNK_BP) -> list[Path]:
    """Split a genome into whole-contig chunks (no coordinate remapping).

    Chunks hold complete contigs, so every reported coordinate stays valid
    without translation. The split is decided from the index *first*, so a
    chunk cannot overshoot by a whole contig — lungfish chromosomes run to
    ~2 Gbp each, and closing a chunk only after passing the limit is how the
    PIEZO port produced a 9.26 Gbp chunk that exhausted memory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    plan_file = out_dir / "chunk_plan.json"
    existing = sorted(out_dir.glob("chunk_*.fna"))
    if existing and plan_file.exists():
        try:
            if json.loads(plan_file.read_text()).get("chunk_bp") == chunk_bp:
                return existing
        except json.JSONDecodeError:
            pass
    for stale in list(out_dir.glob("chunk_*.fna")) + list(out_dir.glob("chunk_*.gff")):
        stale.unlink()          # chunk_bp changed — the old plan is invalid

    idx = read_fai(build_fai(fna))
    plan: list[list[str]] = []
    cur: list[str] = []
    cur_bp = 0
    for name, (length, *_rest) in idx.items():
        if cur and cur_bp + length > chunk_bp:
            plan.append(cur)
            cur, cur_bp = [], 0
        cur.append(name)
        cur_bp += length
    if cur:
        plan.append(cur)
    where = {name: i for i, names in enumerate(plan) for name in names}

    handles = [open(out_dir / f"chunk_{i + 1:03d}.fna", "w")
               for i in range(len(plan))]
    try:
        fh = None
        with open(fna) as src:
            for line in src:
                if line.startswith(">"):
                    fh = handles[where[line[1:].split()[0]]]
                if fh:
                    fh.write(line)
    finally:
        for h in handles:
            h.close()
    plan_file.write_text(json.dumps(
        {"chunk_bp": chunk_bp, "n_chunks": len(plan),
         "chunk_contigs": [len(c) for c in plan],
         "chunk_bp_actual": [sum(idx[n][0] for n in c) for c in plan]}))
    return sorted(out_dir.glob("chunk_*.fna"))


def run_miniprot_chunked(fna: Path, baits_faa: Path, out_gff: Path,
                         work_dir: Path, threads: int = 8,
                         max_intron: int = DEFAULT_MAX_INTRON,
                         chunk_bp: int = CHUNK_BP) -> int:
    """miniprot over contig chunks; concatenate the GFFs. Returns chunk count."""
    chunks = split_fasta_by_contig(fna, work_dir, chunk_bp)
    parts = []
    for i, chunk in enumerate(chunks, 1):
        part = work_dir / f"chunk_{i:03d}.gff"
        if not part.exists() or part.stat().st_size == 0:
            run_miniprot(chunk, baits_faa, part, threads, max_intron)
        parts.append(part)
    tmp = out_gff.with_suffix(out_gff.suffix + ".partial")
    try:
        with open(tmp, "w") as out:
            out.write("##gff-version 3\n")
            for part in parts:
                with open(part) as fh:
                    for line in fh:
                        if not line.startswith("##gff-version"):
                            out.write(line)
        tmp.replace(out_gff)          # atomic, as in run_miniprot
    finally:
        tmp.unlink(missing_ok=True)
    return len(chunks)


# ---------------------------------------------------------------- alignments

@dataclass
class Aln:
    contig: str
    start: int
    end: int
    strand: str
    score: float
    identity: float
    bait: str
    clade: str
    family: str
    q_start: int
    q_end: int
    bait_len: int
    frameshifts: int = 0
    stop_codons: int = 0
    mp_id: str = ""
    translation: str = ""
    aligned_aa: int = 0          # union of CDS query spans (true coverage)
    cds_blocks: list = field(default_factory=list)   # genomic (start, end)

    @property
    def coverage(self) -> float:
        """Fraction of the bait actually aligned (CDS-based, not outer span)."""
        if not self.bait_len:
            return 0.0
        # Frameshift bookkeeping can make CDS query spans overlap slightly, so
        # the union can just exceed the bait length — cap at 1.0.
        if self.aligned_aa:
            return min(1.0, self.aligned_aa / self.bait_len)
        return min(1.0, (self.q_end - self.q_start + 1) / self.bait_len)

    @property
    def span_coverage(self) -> float:
        return ((self.q_end - self.q_start + 1) / self.bait_len
                if self.bait_len else 0.0)


def _attr(attrs: str, key: str) -> str:
    for part in attrs.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


def _union_len(spans: list[tuple[int, int]]) -> int:
    if not spans:
        return 0
    spans.sort()
    total, cs, ce = 0, *spans[0]
    for s, e in spans[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    return total + ce - cs + 1


def parse_miniprot_gff(path: Path, bait_meta: dict[str, dict]) -> list[Aln]:
    """Parse mRNA + CDS features and ##STA translations from miniprot.

    Record order per alignment is `##PAF`, `##STA`, `mRNA`, `CDS...`, so a
    ##STA line belongs to the *next* mRNA, not the previous one. miniprot
    numbers alignments from MP000001 per run, so a chunked genome's
    concatenated GFF repeats every ID; CDS lines always follow their own
    mRNA, so they bind to the record just seen and colliding IDs get a
    unique internal key.
    """
    alns: list[Aln] = []
    by_id: dict[str, Aln] = {}
    q_spans: dict[str, list[tuple[int, int]]] = {}
    pending_sta = ""
    cur_raw = cur_key = ""
    for line in open(path):
        if line.startswith("##STA"):
            pending_sta = line.rstrip("\n").split("\t", 1)[1] if "\t" in line else ""
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        attrs = f[8]
        parts = _attr(attrs, "Target").split()          # "bait qstart qend"
        if f[2] == "mRNA":
            bait = parts[0] if parts else ""
            meta = bait_meta.get(bait, {"clade": "unknown", "family": "unknown",
                                        "length": 0})
            aln = Aln(
                contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                score=float(f[5]) if f[5] not in (".", "") else 0.0,
                identity=float(_attr(attrs, "Identity") or 0.0),
                bait=bait, clade=meta["clade"], family=meta["family"],
                q_start=int(parts[1]) if len(parts) > 2 else 0,
                q_end=int(parts[2]) if len(parts) > 2 else 0,
                bait_len=meta["length"] or 0,
                frameshifts=int(_attr(attrs, "Frameshift") or 0),
                stop_codons=int(_attr(attrs, "StopCodon") or 0),
                mp_id=_attr(attrs, "ID"),
                translation=pending_sta,
            )
            pending_sta = ""
            cur_raw = aln.mp_id
            cur_key = cur_raw
            if cur_key in by_id:
                cur_key = f"{cur_raw}#{len(alns)}"
            aln.mp_id = cur_key
            alns.append(aln)
            by_id[cur_key] = aln
        elif f[2] == "CDS":
            parent = _attr(attrs, "Parent")
            key = cur_key if parent == cur_raw else parent
            if len(parts) > 2:
                q_spans.setdefault(key, []).append(
                    (int(parts[1]), int(parts[2])))
            if key in by_id:
                by_id[key].cds_blocks.append((int(f[3]), int(f[4])))
    for mp_id, spans in q_spans.items():
        if mp_id in by_id:
            by_id[mp_id].aligned_aa = _union_len(spans)
    return alns


# ---------------------------------------------------------------- loci

@dataclass
class Locus:
    contig: str
    strand: str
    start: int
    end: int
    alns: list[Aln] = field(default_factory=list)

    @property
    def best(self) -> Aln:
        return max(self.alns, key=lambda a: a.score)

    @property
    def clade(self) -> str:
        return self.best.clade

    @property
    def family(self) -> str:
        """D14 at locus scale: which family's baits win here.

        The whole point of carrying RyR baits in the panel. Taken from the
        top-scoring alignment, so it is a positive test on alignment score
        and never on a name, a length or a bait's label alone.
        """
        return self.best.family

    def best_of(self, clade: str) -> Aln | None:
        hits = [a for a in self.alns if a.clade == clade]
        return max(hits, key=lambda a: a.score) if hits else None

    def best_of_family(self, family: str) -> Aln | None:
        hits = [a for a in self.alns if a.family == family]
        return max(hits, key=lambda a: a.score) if hits else None

    def family_margin(self) -> float:
        """(winning family score - other family score) / winning score.

        D7's relative margin, reported per locus so a cell whose family call
        is close can be read as close rather than as decided. 1.0 means the
        losing family put no alignment here at all, which is the outcome D14
        is hoping for and the one to check rather than assume.
        """
        itpr = self.best_of_family("ITPR")
        ryr = self.best_of_family("RYR")
        if itpr is None or ryr is None:
            return 1.0
        win, lose = max(itpr.score, ryr.score), min(itpr.score, ryr.score)
        return round((win - lose) / win, 4) if win > 0 else 0.0

    def paralog_margin(self, clade: str) -> float:
        """How far `clade`'s best bait beats the next ITPR paralog's, relative.

        The within-family counterpart of `family_margin`, and the harder
        number: ITPR1/2/3 are 61-68 % identical (S1), so a locus can be
        decisively ITPR and still be only marginally ITPR2 rather than ITPR1.
        Recorded at every locus — including loci whose paralog identity the
        annotation already settles — so `s5_rescue.ATTRIBUTION_REL_MARGIN`
        can be calibrated against loci of known identity instead of being
        inherited from a family with better-separated paralogs.
        """
        mine = self.best_of(clade)
        if mine is None:
            return 0.0
        others = [a.score for a in self.alns
                  if a.family == "ITPR" and a.clade in CLASSES
                  and a.clade != clade]
        if not others:
            return 1.0
        return (round((mine.score - max(others)) / mine.score, 4)
                if mine.score > 0 else 0.0)


def filter_loci(loci: list[Locus],
                min_identity: float = MIN_LOCUS_IDENTITY) -> list[Locus]:
    """Drop clusters whose best alignment is below the identity floor.

    Applied to the *locus*, not to each alignment, so a real locus carrying a
    few poorly-aligned exons survives on the strength of its best one.
    """
    return [L for L in loci if L.best.identity >= min_identity]


def cluster_loci(alns: list[Aln], gap: int = LOCUS_GAP) -> list[Locus]:
    by_key: dict[tuple, list[Aln]] = {}
    for a in alns:
        by_key.setdefault((a.contig, a.strand), []).append(a)
    loci: list[Locus] = []
    for (contig, strand), group in by_key.items():
        group.sort(key=lambda a: a.start)
        cur: Locus | None = None
        for a in group:
            if cur and a.start <= cur.end + gap:
                cur.end = max(cur.end, a.end)
                cur.alns.append(a)
            else:
                cur = Locus(contig, strand, a.start, a.end, [a])
                loci.append(cur)
    return loci
