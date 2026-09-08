"""S12's scope — which loci get read-level evidence, and why those.

The S12 question is not "where is this gene expressed". It is the one S10
handed forward: **the census recovered a gene from the genome that the
annotation does not deliver — is that gene transcribed?**  S10 could not
answer it, and said so in the report rather than half-answering it: 55 and
25 junction probes against every transcript record NCBI holds for the two
case species returned 0 hits, out of denominators of 43 and 10 records.
That is a statement about the deposit, not about the gene.

So the scope is **derived, not chosen**: every locus in S10's committed
ranking that the annotation demonstrably failed on (loss above
``LOSS_FLOOR`` among the loci S10's five eligibility rules admit), plus,
in each of those genomes, **every other family locus the sweep recovered**.
The second half is what makes the first half readable — the other paralogs
are known-real genes measured in the same libraries, on the same reference,
by the same aligner, so a zero on the failed locus has something to be a
zero against.

Four rules, each a positive test:

  **P1** the locus is eligible under S10's rules (it is in an assembly whose
        annotation could reasonably have delivered the gene) and its
        annotation loss exceeds ``LOSS_FLOOR``.  This is read off
        ``case_ranking.tsv``; nothing here re-derives it.
  **P2** its species has at least ``MIN_RUNS`` public Illumina RNA-seq runs.
        A species with no libraries cannot be asked the question, and
        including it would put an unanswerable row in the results table.
  **P3** every other ITPR locus the sweep recovered in the same genome
        joins the reference as an internal positive control, whatever its
        annotation status.
  **P4** the genome's RyR locus joins it too.  RyR is this project's sister
        family (D14) and its reads must land on RyR: a family-level
        cross-mapping control that costs one sequence.

``MIN_RUNS`` is deliberately low.  It is a floor on *askability*, not on
power — a species with 4 libraries gives a weak answer and the report says
so, but a species with 0 gives none at all and belongs in a different
column.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
RANKING = PROJECT / "results" / "annotation_bugs" / "case_ranking.tsv"
LEDGER = PROJECT / "results" / "genome_ledger" / "genome_ledger.tsv"

LOSS_FLOOR = 0.5   # half the coding footprint reaches no single annotated model
MIN_RUNS = 4       # public Illumina RNA-seq runs needed to ask the question


@dataclass
class Locus:
    """One reference sequence's worth of gene, with its provenance."""
    accession: str
    organism: str
    short: str
    cell: str                # ITPR1 / ITPR2 / ITPR3 / RYR
    contig: str
    start: int
    end: int
    strand: str
    role: str                # "failed" | "control_paralog" | "control_family"
    # `None`, not 0.0: S10 ranked ITPR loci, never the RyR control, so an
    # unmeasured locus and a locus measured at zero loss must not share a
    # cell.  A 0.0 there would put "the annotation delivers this gene
    # completely" against a gene S10 never looked at.
    loss: float | None = None
    mode: str = ""           # omission / truncation / fragmentation
    cell_status: str = ""
    annot_gene: str = ""
    # The biotype matters as much as the name. A GFF3 `pseudogene` may
    # carry a full set of CDS features and still emit no protein, so
    # "the annotation has this gene" and "the databases serve this
    # protein" are different statements (S10's distinction).
    annot_biotype: str = ""
    annot_frac_cds: float | None = None
    coverage: float = 0.0
    identity: float = 0.0

    @property
    def name(self) -> str:
        return f"{self.cell}|{self.short}"

    @property
    def locus_key(self) -> tuple[str, int, int, str]:
        return (self.contig, self.start, self.end, self.strand)


@dataclass
class Species:
    name: str
    short: str
    accession: str
    loci: list[Locus] = field(default_factory=list)
    note: str = ""

    @property
    def failed(self) -> list[Locus]:
        return [l for l in self.loci if l.role == "failed"]


def short_name(organism: str) -> str:
    """Filename-safe key: first letter of the genus + the epithet."""
    parts = organism.replace(".", "").split()
    if len(parts) >= 2:
        return (parts[0][0] + parts[1]).lower()
    return organism.lower().replace(" ", "_")


def _rows(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"missing input: {path}")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def failed_loci(loss_floor: float = LOSS_FLOOR) -> list[dict]:
    """P1 — every eligible locus whose annotation loss clears the floor.

    Read straight off S10's committed ranking.  `eligible` is S10's own
    five-rule verdict; re-deriving it here would let two tasks disagree
    about which loci the annotation could reasonably have delivered.
    """
    out = []
    for r in _rows(RANKING):
        if r.get("eligible", "").strip() != "True":
            continue
        try:
            loss = float(r.get("loss") or 0)
        except ValueError:
            continue
        if loss <= loss_floor:
            continue
        r["_loss"] = loss
        out.append(r)
    out.sort(key=lambda r: (-r["_loss"], r["organism"], r["cell"]))
    return out


def build_panel(run_counts: dict[str, int] | None = None,
                loss_floor: float = LOSS_FLOOR,
                min_runs: int = MIN_RUNS) -> tuple[list[Species], list[dict]]:
    """Apply P1-P4.  Returns (panel, audit rows).

    `run_counts` maps organism -> number of public Illumina RNA-seq runs.
    Passing None skips P2 entirely rather than silently treating an unknown
    count as zero: a species must be *measured* to be excluded for want of
    libraries, and a caller that has not measured yet gets the pre-P2 panel.
    """
    audit: list[dict] = []
    failed = failed_loci(loss_floor)
    by_acc: dict[str, Species] = {}

    for r in failed:
        acc, org = r["accession"], r["organism"]
        n = run_counts.get(org) if run_counts is not None else None
        if run_counts is not None and (n or 0) < min_runs:
            audit.append({"accession": acc, "organism": org, "cell": r["cell"],
                          "rule": "P2", "verdict": "excluded",
                          "reason": f"only {n or 0} RNA-seq runs "
                                    f"(< {min_runs})"})
            continue
        sp = by_acc.setdefault(acc, Species(name=org, short=short_name(org),
                                            accession=acc))
        sp.loci.append(Locus(
            accession=acc, organism=org, short=sp.short, cell=r["cell"],
            contig=r["contig"], start=int(r["start"]), end=int(r["end"]),
            strand=r["strand"], role="failed", loss=r["_loss"],
            mode=r.get("mode", ""), cell_status=r.get("cell_status", ""),
            coverage=float(r.get("coverage") or 0),
            identity=float(r.get("identity") or 0)))
        audit.append({"accession": acc, "organism": org, "cell": r["cell"],
                      "rule": "P1", "verdict": "included",
                      "reason": f"annotation loss {r['_loss']:.3f} "
                                f"({r.get('mode', '')})"})

    # P3/P4 — the rest of the family in the same genome, from the ledger.
    # A control paralog S10 did rank carries its measured loss; the RyR
    # control and any locus S10 declared ineligible carry none.
    ledger = {(r["accession"], r["class"]): r for r in _rows(LEDGER)}
    ranked: dict[tuple[str, str], float] = {}
    for r in _rows(RANKING):
        if r.get("eligible", "").strip() != "True":
            continue
        try:
            ranked[(r["accession"], r["cell"])] = float(r.get("loss") or 0)
        except ValueError:
            continue
    for acc, sp in by_acc.items():
        have = {l.cell for l in sp.loci}
        for (a, cell), r in sorted(ledger.items()):
            if a != acc or cell in have:
                continue
            if not r.get("contig") or not r.get("start"):
                audit.append({"accession": acc, "organism": sp.name,
                              "cell": cell, "rule": "P3/P4",
                              "verdict": "unavailable",
                              "reason": f"no recovered locus "
                                        f"({r.get('status', '?')})"})
                continue
            role = "control_family" if cell == "RYR" else "control_paralog"
            measured = ranked.get((acc, cell))
            sp.loci.append(Locus(
                accession=acc, organism=sp.name, short=sp.short, cell=cell,
                contig=r["contig"], start=int(r["start"]), end=int(r["end"]),
                strand=r["strand"], role=role,
                loss=measured,
                cell_status=r.get("status", ""),
                annot_gene=r.get("annot_gene", ""),
                coverage=float(r.get("best_coverage") or 0),
                identity=float(r.get("best_identity") or 0)))
            audit.append({"accession": acc, "organism": sp.name, "cell": cell,
                          "rule": "P4" if cell == "RYR" else "P3",
                          "verdict": "included",
                          "reason": f"internal control ({r.get('status','')}"
                                    f"{', annotated as ' + r['annot_gene'] if r.get('annot_gene') else ''})"})
        order = {"ITPR1": 0, "ITPR2": 1, "ITPR3": 2, "RYR": 3}
        sp.loci.sort(key=lambda l: order.get(l.cell, 9))
        modes = sorted({l.mode for l in sp.failed if l.mode})
        sp.note = (f"{len(sp.failed)} annotation failure"
                   f"{'s' if len(sp.failed) != 1 else ''}"
                   f"{' (' + ', '.join(modes) + ')' if modes else ''}")

    panel = sorted(by_acc.values(), key=lambda s: (-len(s.failed), s.name))
    return panel, audit


AUDIT_COLS = ["accession", "organism", "cell", "rule", "verdict", "reason"]
