#!/usr/bin/env python3
"""S8 synteny — loci, flank extraction, symbol keys, Jaccard.

Inputs are S5 artefacts only:
  results/genome_ledger/genome_ledger.tsv        genome-level metadata + cell status
  <data_root>/genome_sweep/<acc>/summary.json    **every** locus, with coordinates
  <data_root>/genome_sweep/<acc>/genes_slim.tsv  the assembly's own gene intervals

Loci come from the per-genome summaries rather than the ledger because the
ledger carries one row per genome x cell and therefore only the *best* locus.
A teleost ITPR1 cell holds itpr1a and itpr1b; a lamprey ITPR1 cell holds
three. Reading the best one would compare itpr1a in one species against
itpr1b in the next and call the resulting mismatch a synteny result.

Two window rules, both committed, because annotation naming density varies
~4x across this genome set (human 97 % of coding genes named, sea lamprey
27 %):
  fixed10        the brief's rule -- the 10 nearest coding genes each side.
  informative10  the 10 nearest *informative* symbols each side, scanning at
                 most MAX_SCAN genes. Equalises key-set size across genomes,
                 so a Jaccard difference cannot be an annotation-depth
                 difference.
"""
from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.data_root import get_data_root  # noqa: E402

LEDGER = PROJECT_ROOT / "results" / "genome_ledger" / "genome_ledger.tsv"

# The cells S5 swept. RYR is the sister-family control (D14): it rides the
# same genomes through the same code, so "within-family flanks are shared"
# is demonstrated on a second family and "ITPR shares flanks with RyR" is
# available as a negative.
ITPR_CLASSES = ("ITPR1", "ITPR2", "ITPR3")
ALL_CLASSES = ITPR_CLASSES + ("RYR",)

# Cells excluded from the locus sets, with the reason recorded in the report.
# Empty by construction: S7's relabelling rule fired on no tip
# (results/phylogeny/membership_audit.tsv), so no census cell is known to be
# in the wrong paralog. The mechanism is kept so an exclusion cannot be made
# silently.
EXCLUDED_CELLS: dict[tuple[str, str], str] = {}

# 'gene' covers non-RefSeq GCA annotations whose slim table carries no
# biotype attribute.
FLANK_BIOTYPES = {"protein_coding", "gene"}

MAX_SCAN = 60          # genes scanned per side in informative-window mode
DEFAULT_N = 10


@dataclass
class Locus:
    accession: str
    organism: str
    vclass: str
    vorder: str
    cell: str            # ITPR1 / ITPR2 / ITPR3 / RYR (the bait-class cell)
    status: str          # the cell's S5 status
    contig: str
    start: int
    end: int
    idx: int             # 0-based rank of this locus within its cell
    bait: str
    bait_paralog: str
    identity: float
    coverage: float
    annot_gene: str
    annot_paralog: str

    @property
    def label(self) -> str:
        org = self.organism.replace(" ", "_").replace("|", "_")
        return f"{org}|{self.accession}|{self.cell}.{self.idx}"


@dataclass
class FlankSet:
    locus: Locus
    window: str
    upstream: list = field(default_factory=list)    # nearest first
    downstream: list = field(default_factory=list)  # nearest first
    overlapping: list = field(default_factory=list)
    scanned_up: int = 0
    scanned_down: int = 0
    keys_strict: frozenset = frozenset()
    keys_relaxed: frozenset = frozenset()
    keys_root: frozenset = frozenset()

    @property
    def n_flanks(self) -> int:
        return len(self.upstream) + len(self.downstream)

    @property
    def n_informative(self) -> int:
        return len(self.keys_relaxed)


# ------------------------------------------------------------------ ledger

def load_genome_meta(ledger: Path = LEDGER) -> dict[str, dict]:
    """accession -> the genome-level columns (identical across its 4 rows)."""
    meta: dict[str, dict] = {}
    with open(ledger) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta.setdefault(r["accession"], {
                "organism": r["organism"], "vclass": r["vclass"],
                "vorder": r["vorder"], "annotated": r["annotated_assembly"],
                "assembly_level": r["assembly_level"],
                "contig_n50": r["contig_n50"], "reasons": r["reasons"],
            })
    return meta


def cell_status(ledger: Path = LEDGER) -> dict[tuple[str, str], str]:
    out = {}
    with open(ledger) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["accession"], r["class"])] = r["status"]
    return out


def bait_paralog(bait: str) -> str:
    """The paralog label a bait id carries, or its family when unlabelled."""
    parts = bait.split("|")
    return parts[3] if len(parts) >= 4 else (parts[1] if len(parts) > 1 else "")


def load_loci(ledger: Path = LEDGER) -> tuple[list[Locus], list[dict], list[str]]:
    """Every locus of every swept genome, from the per-genome summaries.

    Returns (loci, excluded rows, accessions with no summary on disk).
    """
    meta = load_genome_meta(ledger)
    status = cell_status(ledger)
    root = get_data_root() / "genome_sweep"
    loci, excluded, missing = [], [], []
    for acc in sorted(meta):
        path = root / acc / "summary.json"
        if not path.exists():
            missing.append(acc)
            continue
        summ = json.loads(path.read_text())
        m = meta[acc]
        for cell, c in sorted(summ.get("cells", {}).items()):
            if cell not in ALL_CLASSES:
                continue
            st = status.get((acc, cell), c.get("status", ""))
            if (acc, cell) in EXCLUDED_CELLS:
                excluded.append({"accession": acc, "organism": m["organism"],
                                 "cell": cell, "status": st,
                                 "exclude_reason": EXCLUDED_CELLS[(acc, cell)]})
                continue
            for i, L in enumerate(c.get("loci", [])):
                if not (L.get("contig") and L.get("start") and L.get("end")):
                    continue
                ann = L.get("annot_gene") or {}
                loci.append(Locus(
                    accession=acc, organism=m["organism"], vclass=m["vclass"],
                    vorder=m["vorder"], cell=cell, status=st,
                    contig=L["contig"], start=int(L["start"]), end=int(L["end"]),
                    idx=i, bait=L.get("bait", ""),
                    bait_paralog=bait_paralog(L.get("bait", "")),
                    identity=float(L.get("identity") or 0.0),
                    coverage=float(L.get("coverage") or 0.0),
                    annot_gene=ann.get("name", "") if isinstance(ann, dict) else "",
                    annot_paralog=L.get("annot_paralog") or "",
                ))
    return loci, excluded, missing


# -------------------------------------------------------------- gene table

def slim_path(accession: str) -> Path:
    return get_data_root() / "genome_sweep" / accession / "genes_slim.tsv"


def load_genes(accession: str) -> list[tuple] | None:
    """genes_slim.tsv -> [(contig, start, end, strand, gene_id, name, biotype)]."""
    path = slim_path(accession)
    if not path.exists():
        return None
    genes = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            try:
                genes.append((f[0], int(f[1]), int(f[2]), f[3], f[4], f[5], f[6]))
            except ValueError:
                continue
    return genes


def index_by_contig(genes: list[tuple],
                    biotypes: set[str] = FLANK_BIOTYPES) -> dict[str, list[tuple]]:
    """Coding genes per contig, sorted by start. Built once per genome."""
    by: dict[str, list[tuple]] = {}
    for g in genes:
        if g[6] in biotypes:
            by.setdefault(g[0], []).append(g)
    for v in by.values():
        v.sort(key=lambda g: (g[1], g[2]))
    return by


def extract_flanks(index: dict[str, list[tuple]], locus: Locus,
                   n: int = DEFAULT_N, window: str = "fixed10",
                   max_scan: int = MAX_SCAN) -> FlankSet:
    """The n nearest flanking genes each side of the locus footprint.

    window='fixed10'        n nearest coding genes (the brief's rule).
    window='informative10'  n nearest genes with an informative symbol,
                            scanning at most max_scan genes per side.
    """
    on_contig = index.get(locus.contig, [])
    up_all = [g for g in on_contig if g[2] < locus.start]
    down_all = [g for g in on_contig if g[1] > locus.end]
    over = [g for g in on_contig if g[2] >= locus.start and g[1] <= locus.end]
    up_all.sort(key=lambda g: g[2], reverse=True)     # nearest first
    down_all.sort(key=lambda g: g[1])                 # nearest first

    def take(pool):
        if window == "fixed10":
            return pool[:n], min(len(pool), n)
        kept, scanned = [], 0
        for g in pool[:max_scan]:
            scanned += 1
            if informative(g[5]):
                kept.append(g)
                if len(kept) >= n:
                    break
        return kept, scanned

    up, s_up = take(up_all)
    down, s_down = take(down_all)
    fs = FlankSet(locus=locus, window=window, upstream=up, downstream=down,
                  overlapping=over, scanned_up=s_up, scanned_down=s_down)
    fs.keys_strict, fs.keys_relaxed, fs.keys_root = flank_keys(up + down)
    return fs


# ----------------------------------------------------------------- symbols

_UNINFORMATIVE_PREFIXES = ("si:", "zgc:", "wu:", "im:", "zmp:", "sb:", "id:",
                           "cabz", "cr3", "bx", "cu4", "fp2", "fq3", "hmgn")


def informative(sym: str) -> bool:
    """Is this symbol usable for cross-species matching?

    Rejects the placeholder vocabularies: RefSeq LOC ids, community locus
    tags (FN964_004414, OJAV_G00063340), Ensembl stable ids, and the
    clone-derived zebrafish prefixes. 27 % of sea-lamprey coding genes
    survive this and 97 % of human ones, which is why the informative
    window exists.
    """
    if not sym or len(sym) < 2:
        return False
    s = sym.lower()
    if s.startswith(_UNINFORMATIVE_PREFIXES):
        return False
    u = sym.upper()
    if u.startswith("LOC") and u[3:].isdigit():
        return False
    if "_" in u:
        tail = u.rsplit("_", 1)[1]
        if tail.isdigit() or (tail[:1] == "G" and tail[1:].isdigit()):
            return False
    if u.startswith("ENS") and any(c.isdigit() for c in u[-6:]):
        return False
    if not any(c.isalpha() for c in sym):
        return False
    return True


def relaxed_key(sym: str) -> str:
    """Canonical cross-clade key.

    Teleost and amphibian symbols are lowercase and may carry a duplicate
    suffix -- gnasa/gnasb, or a trailing '.2' -- that would break matching
    against the tetrapod symbol (GNAS). Only all-lowercase symbols are
    stripped, so uppercase mammal/bird symbols ending in A/B (GNB1, BAK1)
    are never touched.
    """
    s = sym.strip()
    if s != s.lower():
        return s.upper()
    core = s
    while "." in core and core.rsplit(".", 1)[1].isdigit():
        core = core.rsplit(".", 1)[0]
    # >=5 chars before stripping: gnasa->GNAS but vapa/ctsa keep their suffix
    # (short cores collide: CTSA is not CTS, VAPB must still match VAPB)
    if len(core) > 4 and core[-1] in "ab" and not core[-2].isdigit():
        core = core[:-1]
    return core.upper()


def root_key(sym: str) -> str:
    """The named gene *family* a symbol belongs to: its relaxed key with the
    trailing digit run removed.

    BHLHE40 and BHLHE41 -> BHLHE4; SLC25A3 and SLC25A5 -> SLC25A. This is
    what makes an ohnolog pair visible to a set comparison, because the two
    copies of a 2R quartet almost never carry the same symbol.

    It is a heuristic and it over-merges (TP53 -> TP5). That is tolerable
    only because the null distribution it is scored against is built with
    the *same* rule on random windows, so over-merging inflates signal and
    null together. Roots shorter than 3 characters are dropped rather than
    kept as near-universal keys.
    """
    k = relaxed_key(sym)
    core = k.rstrip("0123456789")
    if len(core) < 3 or not any(c.isalpha() for c in core):
        return k
    return core


def flank_keys(flank_genes: list[tuple]) -> tuple[frozenset, frozenset, frozenset]:
    """(strict, relaxed, root) key sets over the informative flank symbols."""
    strict, relaxed, root = set(), set(), set()
    for g in flank_genes:
        sym = g[5]
        if not informative(sym):
            continue
        strict.add(sym.upper())
        relaxed.add(relaxed_key(sym))
        root.add(root_key(sym))
    return frozenset(strict), frozenset(relaxed), frozenset(root)


# ----------------------------------------------------------------- jaccard

def jaccard(a: frozenset, b: frozenset) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def keys_of(fs: FlankSet, kind: str) -> frozenset:
    return {"strict": fs.keys_strict, "relaxed": fs.keys_relaxed,
            "root": fs.keys_root}[kind]


def pairwise_jaccard(sets_a: list[FlankSet], sets_b: list[FlankSet] | None,
                     kind: str = "relaxed",
                     cross_genome_only: bool = True) -> list[tuple]:
    """All pairs, within one list or across two.

    Returns (label_a, label_b, jaccard, n_shared, acc_a, acc_b).
    Same-genome pairs are dropped by default: two loci in one assembly share
    that assembly's naming conventions, which is not the signal being
    measured.
    """
    out = []
    if sets_b is None:
        for i in range(len(sets_a)):
            for j in range(i + 1, len(sets_a)):
                fa, fb = sets_a[i], sets_a[j]
                if cross_genome_only and fa.locus.accession == fb.locus.accession:
                    continue
                ka, kb = keys_of(fa, kind), keys_of(fb, kind)
                out.append((fa.locus.label, fb.locus.label, jaccard(ka, kb),
                            len(ka & kb), fa.locus.accession, fb.locus.accession))
    else:
        for fa in sets_a:
            for fb in sets_b:
                if cross_genome_only and fa.locus.accession == fb.locus.accession:
                    continue
                ka, kb = keys_of(fa, kind), keys_of(fb, kind)
                out.append((fa.locus.label, fb.locus.label, jaccard(ka, kb),
                            len(ka & kb), fa.locus.accession, fb.locus.accession))
    return out


def matrix_from_pairs(labels: list[str], pairs: list[tuple]) -> dict:
    m = {(a, b): j for a, b, j, _, _, _ in pairs}
    m.update({(b, a): j for a, b, j, _, _, _ in pairs})
    for lab in labels:
        m[(lab, lab)] = 1.0
    return m


def consensus(flank_sets: list[FlankSet], kind: str = "relaxed"):
    """Per symbol: how many of these loci carry it, and the fraction."""
    from collections import Counter
    counts: Counter = Counter()
    n = 0
    for fs in flank_sets:
        keys = keys_of(fs, kind)
        if not keys:
            continue
        n += 1
        counts.update(keys)
    return [(sym, c, c / n if n else 0.0) for sym, c in counts.most_common()], n
