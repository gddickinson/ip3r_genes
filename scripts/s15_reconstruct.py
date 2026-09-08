"""s15_reconstruct.py — how much of a reference protein an assembly holds
when it holds none of it as one gene.

The sweep's `tblastn_trace` cells are the ones a loss count cannot ignore
and cannot use: the aligner found no locus, so the ledger records no
coverage, and a character matrix built from the ledger alone would read
them as candidate absences.  What this module measures instead is the
**non-redundant coverage of the reference protein, reassembled across
contigs**: the union of the query intervals of every significant HSP,
divided by the reference's own length.

Three properties make that a positive test rather than a hopeful one.

**It is computed outside every locus the aligner found.**  The HSP set is
filtered through `s15_lib.outside_known`, which reproduces the sweep's own
`filter_hsps_outside` exclusion — every clustered locus in the genome, any
bait, padded by 5 kb.  So a reference cannot be reassembled out of the
genome's other paralogs' genes, which for a family whose paralogs are
61-68 % identical is the failure mode that matters.

**It is split by the full panel's own attribution, into three
populations and not two.**  `attribute_region` re-blasted each rescue
region against all 38 baits and recorded every clade's bits, so each
region carries a family and a paralog call made on evidence the cell's own
bait set did not produce.  `own_clade` is the cell's own call.  The decoy
is *not* simply "attributed elsewhere": in a genome where two paralogs are
both shattered, a fragment attributed to the other one is a piece of a
real gene, and using it as a negative would measure the wrong thing.  The
negative is `decoy_accounted` — regions attributed to a paralog the
aligner **already found at a locus in this genome**, so the fragment
cannot be that gene.  The third population, `co_trace`, is the
mutually-confusable middle, reported rather than folded into either side:
it is the measured size of the paralog-attribution problem in a shattered
assembly.

**It reports the geometry, not only the number.**  A gene reassembled out
of 15 contigs in a 20 kb-N50 assembly is a different object from one
recovered in two pieces, and `n_contigs` / `n_regions` / `max_piece` are
what let the report say which.  `unique_frac` — the share of covered
reference positions reached from exactly one contig — is the one that
tells a tiling from a pile-up: fragments of three mutually 61-68 %
identical paralogs would stack on the conserved blocks, a shattered gene
tiles its reference once.

**And it counts, on the subject side.**  Reference coverage cannot count
copies: a genome holding only ITPR1 recovers most of the ITPR2 reference
too, at 65 % identity.  `exonic_nt` is the union of the *genomic*
intervals the HSPs align to, which is not shared between paralogs because
they are different places in the genome; divided by three times a
reference's length it is `gene_equiv`, gene-equivalents of family exonic
sequence lying outside every locus the aligner found.  That is the
statistic a loss count can be built on, and it is additive across cells.
"""

from __future__ import annotations

import collections
from pathlib import Path

import s15_lib as lib


class Reconstruction:
    """One cell's reassembly of one reference, with its geometry."""

    __slots__ = ("accession", "organism", "vclass", "cell", "bait",
                 "bait_len", "scope", "covered_aa", "n_hsps", "n_regions",
                 "n_contigs", "max_piece_aa", "bits", "q_intervals",
                 "q_depth", "unique_frac", "exonic_nt", "s_intervals")

    def __init__(self, accession, organism, vclass, cell, bait, bait_len,
                 scope):
        self.accession, self.organism, self.vclass = accession, organism, vclass
        self.cell, self.bait, self.bait_len, self.scope = (
            cell, bait, bait_len, scope)
        self.covered_aa = 0
        self.n_hsps = 0
        self.n_regions = 0
        self.n_contigs = 0
        self.max_piece_aa = 0
        self.bits = 0.0
        self.unique_frac = 0.0
        self.exonic_nt = 0
        self.q_intervals: list[tuple[int, int]] = []
        self.q_depth: dict[int, set] = {}
        self.s_intervals: dict[str, list[tuple[int, int]]] = {}

    @property
    def coverage(self) -> float:
        return self.covered_aa / self.bait_len if self.bait_len else 0.0

    @property
    def gene_equiv(self) -> float:
        """Gene-equivalents of family exonic sequence, subject side."""
        return self.exonic_nt / (3.0 * self.bait_len) if self.bait_len else 0.0

    def row(self) -> dict:
        return dict(accession=self.accession, organism=self.organism,
                    vclass=self.vclass, cell=self.cell, scope=self.scope,
                    bait=self.bait, bait_len=self.bait_len,
                    covered_aa=self.covered_aa, coverage=round(self.coverage, 4),
                    n_hsps=self.n_hsps, n_regions=self.n_regions,
                    n_contigs=self.n_contigs, max_piece_aa=self.max_piece_aa,
                    unique_frac=round(self.unique_frac, 4),
                    exonic_nt=self.exonic_nt,
                    gene_equiv=round(self.gene_equiv, 4),
                    bits=round(self.bits, 1))


def parse_hsps(path: Path, e_cut: float = lib.RESCUE_E) -> list[dict]:
    """The sweep's archived tblastn, `qstart qend` kept.

    `s5_rescue.parse_tblastn` drops the query coordinates — it only needed
    the subject side to cluster regions — and the query side is the whole
    measurement here, so the parse is repeated rather than imported.
    """
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        try:
            ev = float(f[8])
        except ValueError:
            continue
        if ev > e_cut:
            continue
        s_lo, s_hi = sorted((int(f[6]), int(f[7])))
        out.append(dict(bait=f[0], contig=f[1], q_lo=int(f[4]), q_hi=int(f[5]),
                        s_lo=s_lo, s_hi=s_hi, bits=float(f[9]), evalue=ev))
    return out


#: the three scopes, in the order the calibration reads them
SCOPES = ("own_clade", "decoy_accounted", "co_trace")


def region_index(regions: list[dict], cell: str,
                 cell_status: dict[str, str]) -> dict[str, list[tuple]]:
    """Region footprints for one cell, split into the three scopes.

    `assigned_clade` is the full-panel call recorded by
    `s5_rescue.attribute_region`; `vertebrate_basal` is deliberately *not*
    folded into own-clade — it is an unlabelled band, and
    `s15_states.py` decides what one means.

    `cell_status` is the genome's own ledger row per cell, and it is what
    makes the negative a negative: a paralog whose gene sits at a locus is
    accounted for, so a fragment attributed to it is not that gene.
    """
    out: dict[str, list[tuple]] = {s: [] for s in SCOPES}
    for r in regions:
        foot = (r["contig"], int(r["start"]), int(r["end"]))
        clade = r.get("assigned_clade") or ""
        if clade == cell:
            out["own_clade"].append(foot)
        elif clade in cell_status:
            st = cell_status[clade]
            if st.startswith("found"):
                out["decoy_accounted"].append(foot)
            elif st.startswith("tblastn_trace"):
                out["co_trace"].append(foot)
    return out


def reconstruct(accession: str, organism: str, vclass: str, cell: str,
                hsps: list[dict], known: list[tuple[str, int, int]],
                footprints: list[tuple], bait_len: dict[str, int],
                scope: str) -> Reconstruction | None:
    """Best single reference reassembled from HSPs inside `footprints`.

    One reference at a time, never a union over baits: a union over three
    orthologous baits would count a residue covered in any of them and
    report a coverage no single protein achieves.  The cell's answer is
    the best-covered reference of its own clade.
    """
    per: dict[str, Reconstruction] = {}
    seen_regions: dict[str, set] = collections.defaultdict(set)
    for h in hsps:
        if not lib.outside_known(h["contig"], h["s_lo"], h["s_hi"], known):
            continue
        hit = None
        for i, (c, s, e) in enumerate(footprints):
            if c == h["contig"] and lib.overlaps(h["s_lo"], h["s_hi"], s, e):
                hit = i
                break
        if hit is None:
            continue
        b = h["bait"]
        blen = bait_len.get(b)
        if not blen:
            continue
        rec = per.get(b)
        if rec is None:
            rec = per[b] = Reconstruction(accession, organism, vclass, cell,
                                          b, blen, scope)
        rec.q_intervals.append((h["q_lo"], h["q_hi"]))
        for q in range(h["q_lo"], h["q_hi"] + 1):
            rec.q_depth.setdefault(q, set()).add(h["contig"])
        rec.s_intervals.setdefault(h["contig"], []).append(
            (h["s_lo"], h["s_hi"]))
        rec.n_hsps += 1
        rec.bits += h["bits"]
        rec.max_piece_aa = max(rec.max_piece_aa, h["q_hi"] - h["q_lo"] + 1)
        seen_regions[b].add(hit)
    if not per:
        return None
    for b, rec in per.items():
        rec.covered_aa = lib.union_length(rec.q_intervals)
        rec.n_regions = len(seen_regions[b])
        rec.n_contigs = len({footprints[i][0] for i in seen_regions[b]})
        if rec.q_depth:
            rec.unique_frac = sum(
                1 for v in rec.q_depth.values() if len(v) == 1
            ) / len(rec.q_depth)
    best = max(per.values(), key=lambda r: (r.covered_aa, r.bits))
    # exonic_nt is a property of the cell's genomic evidence, not of one
    # bait: every bait of the clade aligns to the same places, so the union
    # is taken over all of them and attached to the reported reference.
    merged: dict[str, list[tuple[int, int]]] = {}
    for rec in per.values():
        for c, ivs in rec.s_intervals.items():
            merged.setdefault(c, []).extend(ivs)
    best.exonic_nt = sum(lib.union_length(v) for v in merged.values())
    return best


def cell_reconstructions(acc: str, summary: dict, cell: str,
                         regions: list[dict], bait_len: dict[str, int],
                         cell_status: dict[str, str]
                         ) -> dict[str, Reconstruction | None]:
    """All three scopes for one cell."""
    path = lib.sweep_dir() / acc / f"tblastn_{cell}.tsv"
    hsps = parse_hsps(path)
    if not hsps:
        return {s: None for s in SCOPES}
    known = lib.known_loci(summary)
    foots = region_index(regions, cell, cell_status)
    org, vcl = summary.get("organism", ""), summary.get("vclass", "")
    return {s: reconstruct(acc, org, vcl, cell, hsps, known, foots[s],
                           bait_len, s) for s in SCOPES}


def run(ledger: list[dict], rescue_rows: list[dict],
        summaries: dict[str, dict]) -> list[dict]:
    """One row per (candidate cell x scope) that has archived HSPs."""
    bait_len = lib.bait_lengths()
    by_cell: dict[tuple[str, str], list[dict]] = collections.defaultdict(list)
    for r in rescue_rows:
        by_cell[(r["accession"], r["cell"])].append(r)
    status: dict[str, dict[str, str]] = collections.defaultdict(dict)
    for r in ledger:
        status[r["accession"]][r["class"]] = r["status"]
    rows: list[dict] = []
    todo = [r for r in ledger if r["status"].startswith("tblastn_trace")
            or r["status"] == "absent"]
    for i, r in enumerate(todo):
        acc, cell = r["accession"], r["class"]
        summary = summaries.get(acc)
        if summary is None:
            continue
        recs = cell_reconstructions(acc, summary, cell,
                                    by_cell.get((acc, cell), []), bait_len,
                                    status[acc])
        for scope, rec in recs.items():
            if rec is None:
                rows.append(dict(accession=acc, organism=r["organism"],
                                 vclass=r["vclass"], cell=cell, scope=scope,
                                 bait="", bait_len=0, covered_aa=0,
                                 coverage=0.0, n_hsps=0, n_regions=0,
                                 n_contigs=0, max_piece_aa=0,
                                 unique_frac=0.0, exonic_nt=0,
                                 gene_equiv=0.0, bits=0.0,
                                 cell_status=r["status"]))
            else:
                row = rec.row()
                row["cell_status"] = r["status"]
                rows.append(row)
        if i % 8 == 0:
            lib.live("reconstruct", [(f"cells {i+1}/{len(todo)}", False)])
    return rows


COLS = ["accession", "organism", "vclass", "cell", "cell_status", "scope",
        "bait", "bait_len", "covered_aa", "coverage", "exonic_nt",
        "gene_equiv", "unique_frac", "n_hsps", "n_regions", "n_contigs",
        "max_piece_aa", "bits"]
