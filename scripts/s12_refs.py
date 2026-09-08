"""S12 step 1 — build the per-species expression reference.

Each species gets one small FASTA holding, as separate sequences:

  * **every family locus the sweep recovered in that genome** — the loci
    S10 showed the annotation loses, and the ones it does not, spliced out
    of the assembly by `s12_locus` and validated against the sweep's own
    protein.  Same construction for the target and for its controls, so a
    difference between them cannot be a difference in how they were built.
  * **housekeeping anchors** — GAPDH / EEF1A1 / RPL13A, aligned to the
    genome with miniprot and spliced by the *same* module.  Two of the
    three assemblies here name no genes at all (*D. mawsoni* files every
    one of its 29,240 genes as "hypothetical protein"), so an anchor pulled
    out of a GFF exists for one species and not the others; aligning a bait
    works everywhere and, more importantly, builds the anchor and the
    target the same way.
  * **a reversed decoy per reference sequence** — the sequence reversed,
    not reverse-complemented.  Reversal preserves length and base
    composition exactly and destroys homology, so whatever a decoy collects
    is this reference's spurious-mapping floor under this aligner.

It also records, per reference sequence, **every exon junction in CDS
coordinates and whether the annotation models it**.  That is the S12
measurement: a read aligned contiguously across a junction the annotation
does not model is transcript evidence for sequence no annotated gene
delivers, which is precisely what S10 could not obtain from the transcript
deposits and handed here.

Outputs (committed):
    results/expression/reference_table.tsv
    results/expression/junctions.tsv
    results/expression/panel_audit.tsv
Bulk (data root):
    <data_root>/expression/refs/<species>.fasta  (+ hisat2 index)

Run:  python scripts/s12_refs.py [--species nalbiflora,...] [--no-housekeeping]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s12_locus  # noqa: E402
import s12_panel  # noqa: E402
from s12_lib import (  # noqa: E402
    CACHE, DATA, GENOMES, OUT_DIR, SWEEP, cached_get, live, read_fasta,
    tool_bin, write_fasta, write_tsv,
)
from s9_miniprot_cds import index_models  # noqa: E402

# Housekeeping baits.  Three, from a teleost, because one anchor that fails
# to align leaves a library with no scale at all; and a spread of expression
# levels (EEF1A1 very high, RPL13A high, GAPDH high but tissue-variable) so
# the anchor is not one gene's quirk.
#
# Each is **resolved by query with a declared expected length**, never by a
# remembered accession (S11's rule for its structural references).  The
# first version of this module hard-coded three accessions and got all
# three wrong in ways nothing downstream would have caught: a 427 aa
# "GAPDH" (the enzyme is ~333), a 1,829 aa "EEF1A1" (~462) and an
# accession that served nothing at all.  A bait of the wrong protein still
# aligns *somewhere* and still produces reads, so the anchor would have
# silently measured a different gene.
HOUSEKEEPING_BAITS = {
    "GAPDH": {"gene": "gapdh", "length": (320, 350)},
    "EEF1A1": {"gene": "eef1a1a", "length": (440, 480)},
    "RPL13A": {"gene": "rpl13a", "length": (180, 220)},
}
HK_ORGANISM = 7955          # Danio rerio
UNIPROT_SEARCH = ("https://rest.uniprot.org/uniprotkb/search?query="
                  "gene_exact:{gene}+AND+organism_id:{taxid}"
                  "&fields=accession,length,protein_name,reviewed"
                  "&format=tsv&size=20")
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

MINIPROT_MAX_INTRON = 200_000   # housekeeping genes are small; the default
MIN_HK_COVERAGE = 0.60          # bait coverage for a usable anchor


class BaitError(RuntimeError):
    """No record matched a housekeeping bait's declared length band."""


def resolve_housekeeping(symbol: str, spec: dict) -> tuple[str, int]:
    """(accession, length) for a housekeeping bait, by query.

    Reviewed entries first, then length band, then shortest — the canonical
    enzyme rather than a long isoform or a fusion.  A symbol whose search
    returns nothing inside the band raises: a housekeeping anchor that is
    quietly the wrong protein is worse than no anchor, because it still
    produces reads.
    """
    lo, hi = spec["length"]
    url = UNIPROT_SEARCH.format(gene=spec["gene"], taxid=HK_ORGANISM)
    txt = cached_get(url, kind="uniprot_search", key=url)
    rows = [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
    cand = []
    for r in rows:
        if len(r) < 4:
            continue
        try:
            length = int(r[1])
        except ValueError:
            continue
        if not lo <= length <= hi:
            continue
        cand.append((0 if r[3].lower().startswith("review") else 1,
                     length, r[0]))
    if not cand:
        raise BaitError(
            f"{symbol}: no {spec['gene']} record in taxid {HK_ORGANISM} "
            f"with length in {lo}-{hi} (searched {len(rows)} records)")
    cand.sort()
    return cand[0][2], cand[0][1]


def uniprot_seq(acc: str) -> str:
    txt = cached_get(UNIPROT_FASTA.format(acc=acc), kind="uniprot", key=acc)
    seq = "".join(l.strip() for l in txt.split("\n")[1:]
                  if not l.startswith(">"))
    if not seq:
        raise BaitError(f"UniProt served no sequence for {acc}")
    return seq


def genome_files(accession: str) -> tuple[Path, Path | None]:
    gdir = GENOMES / accession
    fna = next(iter(sorted(gdir.glob("*_genomic.fna"))), None)
    if fna is None:
        raise SystemExit(f"no genome FASTA under {gdir} — re-fetch it with "
                         f"scripts/fetch_genomes.py --only {accession}")
    gff = gdir / "genomic.gff.gz"
    return fna, (gff if gff.exists() else None)


def locus_models(accession: str) -> dict:
    """`summary.json`'s loci, keyed by coordinates -> the sweep's own record.

    S5's summaries carry the calls (D13); nothing here re-derives which
    miniprot model belongs to which cell.
    """
    path = SWEEP / accession / "summary.json"
    if not path.exists():
        raise SystemExit(f"no sweep summary at {path}")
    summ = json.loads(path.read_text())
    out = {}
    for cell, rec in summ.get("cells", {}).items():
        for loc in rec.get("loci", []):
            key = (cell, loc["contig"], loc["start"], loc["end"], loc["strand"])
            out[key] = loc
    return out


# ---------------------------------------------------------------------------
# housekeeping anchors via miniprot
# ---------------------------------------------------------------------------

def run_housekeeping_miniprot(accession: str, fna: Path, threads: int,
                              log=print) -> Path:
    """Align the housekeeping baits to one genome (cached under the data root)."""
    out = DATA / "housekeeping" / f"{accession}.gff"
    if out.exists() and out.stat().st_size > 0:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    baits = out.parent / "baits.faa"
    if not baits.exists():
        seqs = {}
        for sym, spec in HOUSEKEEPING_BAITS.items():
            acc, length = resolve_housekeeping(sym, spec)
            seq = uniprot_seq(acc)
            if len(seq) != length:
                raise BaitError(f"{sym}: {acc} search said {length} aa, "
                                f"FASTA served {len(seq)}")
            log(f"    bait {sym}: {acc} ({len(seq)} aa)")
            seqs[sym] = seq
        write_fasta(baits, seqs)
    # Same flags as the sweep (`s5_sweep_lib.run_miniprot`), and `--trans`
    # is not optional: without it miniprot writes no `##STA` line, every
    # model comes back with an empty protein, and the coverage filter below
    # silently rejects all of them — which reads as "this genome has no
    # housekeeping genes" rather than "the aligner was asked the wrong way".
    log(f"    miniprot housekeeping -> {accession} (indexes the genome; "
        f"minutes)")
    tmp = out.with_suffix(".partial")
    with open(tmp, "w") as fh:
        p = subprocess.run(
            [tool_bin("miniprot"), "-t", str(threads), "--gff", "--trans",
             "--outs=0.3", "-N", "60", "-G", str(MINIPROT_MAX_INTRON),
             str(fna), str(baits)],
            stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"miniprot failed on {accession}: {p.stderr[:300]}")
    tmp.rename(out)
    return out


def housekeeping_loci(gff: Path, baits: dict[str, str]) -> dict[str, dict]:
    """Best model per housekeeping symbol, with how many copies the genome has.

    The copy count is not bookkeeping: the reference holds **one** model per
    anchor, so reads from a gene's other genomic copies multi-map and are
    cut by the MAPQ floor.  That is why the three anchors are not equally
    informative — measured on these genomes, EEF1A1 has 5-6 copies, GAPDH 3
    and RPL13A 1-3 — and without the number a near-zero count on a
    housekeeping gene reads as a library that failed rather than as an
    anchor whose reads went to its own paralogs.
    """
    best: dict[str, dict] = {}
    copies: dict[str, int] = {}
    for m in index_models(gff):
        sym = m["bait"]
        if sym not in baits:
            continue
        copies[sym] = copies.get(sym, 0) + 1
        cov = len(m["protein"]) / max(1, len(baits[sym]))
        if cov < MIN_HK_COVERAGE:
            continue
        cur = best.get(sym)
        if cur is None or cov > cur["_cov"]:
            m["_cov"] = round(cov, 4)
            best[sym] = m
    for sym, m in best.items():
        m["_copies"] = copies.get(sym, 0)
    return best


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

def build_species(sp, threads: int, with_housekeeping: bool, log=print
                  ) -> tuple[dict[str, str], list[dict], list[dict],
                             list[dict]]:
    fna, gff = genome_files(sp.accession)
    sweep_gff = SWEEP / sp.accession / "miniprot.gff"
    proteins = {m["id"]: m["protein"] for m in index_models(sweep_gff)}
    loci = locus_models(sp.accession)

    refs: dict[str, str] = {}
    rows: list[dict] = []
    junc_rows: list[dict] = []
    valid_rows: list[dict] = []

    def add(name: str, built: dict, role: str, **extra) -> None:
        refs[name] = built["sequence"]
        rows.append({
            "species": sp.short, "organism": sp.name,
            "accession": sp.accession, "seq": name, "role": role,
            "length": built["length"], "n_exons": built["n_exons"],
            "n_junctions": built["n_junctions"],
            "n_junctions_annotated": built["n_junctions_annotated"],
            "n_junctions_unannotated": built["n_junctions_unannotated"],
            "contig": built["contig"], "start": built["start"],
            "end": built["end"], "strand": built["strand"],
            "mp_id": built["mp_id"], "bait": built["bait"],
            "frameshifts": built["frameshifts"],
            "frame_breaks": built["frame_breaks"],
            "blocks_placed": built["blocks_placed"],
            "blocks_tested": built["blocks_tested"],
            "placement": built["placement"],
            "translation_identity": built["translation_identity"],
            "source": f"miniprot:{sp.accession}:{built['mp_id']}",
            **extra})
        c = built.get("controls", {})
        row = {"species": sp.short, "organism": sp.name, "seq": name,
               "role": role, "frameshifts": built["frameshifts"],
               "cds_bp_expected": extra.pop("_cds_bp", ""),
               "cds_bp_built": built["length"]}
        for tag in ("as_built", "reversed_order", "wrong_strand"):
            placed, tested = c.get(tag, (0, 0))
            row[f"{tag}_placed"] = placed
            row[f"{tag}_tested"] = tested
            row[f"{tag}_rate"] = round(placed / tested, 4) if tested else ""
        valid_rows.append(row)
        for j in built["junctions"]:
            junc_rows.append({"species": sp.short, "seq": name,
                              "junction_index": j["junction_index"],
                              "cds_offset": j["cds_offset"],
                              "left_exon": j["left_exon"],
                              "right_exon": j["right_exon"],
                              "intron_start": j["intron_start"],
                              "intron_end": j["intron_end"],
                              "intron_length": j["intron_length"],
                              "class": j["class"], "annotated": j["annotated"],
                              "donor": j["donor"], "acceptor": j["acceptor"],
                              "splice_class": j["splice_class"],
                              "left_models": j["left_models"],
                              "right_models": j["right_models"]})

    for locus in sp.loci:
        key = (locus.cell, locus.contig, locus.start, locus.end, locus.strand)
        rec = loci.get(key)
        if rec is None:
            log(f"    !! {locus.cell}: no sweep locus at {locus.locus_key}")
            continue
        mp_id = rec.get("mp_id", "")
        prot = proteins.get(mp_id, "")
        if not prot:
            log(f"    !! {locus.cell}: no ##STA protein for {mp_id}")
            continue
        try:
            built = s12_locus.build(fna, gff, sweep_gff, mp_id, prot,
                                    expected_cds_bp=int(rec.get("cds_bp") or 0))
        except s12_locus.LocusError as exc:
            log(f"    !! {locus.cell}: {exc}")
            continue
        ag = rec.get("annot_gene") or {}
        add(locus.name, built, locus.role, cell=locus.cell,
            _cds_bp=int(rec.get("cds_bp") or 0),
            annot_biotype=ag.get("biotype", ""),
            annot_frac_cds=ag.get("frac_cds", ""),
            annotation_loss="" if locus.loss is None else locus.loss,
            failure_mode=locus.mode,
            cell_status=locus.cell_status, annot_gene=locus.annot_gene,
            sweep_coverage=locus.coverage, sweep_identity=locus.identity)
        log(f"    {locus.name:22s} {built['length']:>6} nt  "
            f"{built['n_exons']:>3} exons  "
            f"{built['n_junctions_annotated']:>3}/{built['n_junctions']} "
            f"junctions annotated  placed={built['blocks_placed']}/"
            f"{built['blocks_tested']}")

    if with_housekeeping:
        hk_gff = run_housekeeping_miniprot(sp.accession, fna, threads, log)
        baits = read_fasta(DATA / "housekeeping" / "baits.faa")
        for sym, m in sorted(housekeeping_loci(hk_gff, baits).items()):
            try:
                built = s12_locus.build(fna, gff, hk_gff, m["id"],
                                        m["protein"])
            except s12_locus.LocusError as exc:
                log(f"    !! hk {sym}: {exc}")
                continue
            add(f"hk_{sym}|{sp.short}", built, "housekeeping", cell=sym,
                sweep_coverage=m["_cov"], genome_copies=m.get("_copies", ""))
            log(f"    hk_{sym:19s} {built['length']:>6} nt  "
                f"{built['n_exons']:>3} exons  cov={m['_cov']}")

    # composition-matched, homology-free decoys, one per reference sequence
    for r in list(rows):
        name = f"decoy_{r['seq']}"
        refs[name] = refs[r["seq"]][::-1]
        rows.append({"species": sp.short, "organism": sp.name,
                     "accession": sp.accession, "seq": name, "role": "decoy",
                     "length": r["length"], "n_exons": 0, "n_junctions": 0,
                     "n_junctions_annotated": 0, "n_junctions_unannotated": 0,
                     "cell": r.get("cell", ""),
                     "source": f"reversed({r['seq']})",
                     "note": "reversed CDS: identical base composition, "
                             "no homology"})
    return refs, rows, junc_rows, valid_rows


def hisat2_index(short: str, fa: Path, log=print) -> bool:
    idx = DATA / "index" / short
    idx.parent.mkdir(parents=True, exist_ok=True)
    if Path(str(idx) + ".1.ht2").exists():
        return True
    try:
        subprocess.run([tool_bin("hisat2-build"), "-q", str(fa), str(idx)],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        log(f"  !! hisat2-build failed for {short}: {exc}")
        return False
    log(f"  {short}: hisat2 index -> {idx}")
    return True


REF_COLS = ["species", "organism", "accession", "seq", "role", "cell",
            "length", "n_exons", "n_junctions", "n_junctions_annotated",
            "n_junctions_unannotated", "annotation_loss", "failure_mode",
            "cell_status", "annot_gene", "annot_biotype", "annot_frac_cds",
            "genome_copies", "sweep_coverage", "sweep_identity",
            "contig", "start", "end", "strand", "mp_id", "bait",
            "frameshifts", "frame_breaks", "blocks_placed", "blocks_tested",
            "placement", "translation_identity", "source", "note"]
VALID_COLS = ["species", "organism", "seq", "role", "frameshifts",
              "cds_bp_expected", "cds_bp_built",
              "as_built_placed", "as_built_tested", "as_built_rate",
              "reversed_order_placed", "reversed_order_tested",
              "reversed_order_rate", "wrong_strand_placed",
              "wrong_strand_tested", "wrong_strand_rate"]
JUNC_COLS = ["species", "seq", "junction_index", "cds_offset", "left_exon",
             "right_exon", "intron_start", "intron_end", "intron_length",
             "class", "annotated", "donor", "acceptor", "splice_class",
             "left_models", "right_models"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="", help="comma-separated subset")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--no-housekeeping", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    panel, audit = s12_panel.build_panel()
    if args.species:
        want = set(args.species.split(","))
        panel = [s for s in panel if s.short in want]

    all_rows: list[dict] = []
    all_junc: list[dict] = []
    all_valid: list[dict] = []
    for i, sp in enumerate(panel, 1):
        print(f"[refs] {sp.name} ({sp.short}) — {sp.note}", flush=True)
        live("refs", i - 1, len(panel), sp.name)
        refs, rows, junc, valid = build_species(sp, args.threads,
                                                not args.no_housekeeping)
        fa = DATA / "refs" / f"{sp.short}.fasta"
        write_fasta(fa, refs)
        print(f"  {len(refs)} sequences -> {fa}")
        hisat2_index(sp.short, fa)
        all_rows += rows
        all_junc += junc
        all_valid += valid

    write_tsv(OUT_DIR / "reference_table.tsv", REF_COLS, all_rows)
    write_tsv(OUT_DIR / "junctions.tsv", JUNC_COLS, all_junc)
    write_tsv(OUT_DIR / "reference_validation.tsv", VALID_COLS, all_valid)
    write_tsv(OUT_DIR / "panel_audit.tsv", s12_panel.AUDIT_COLS, audit)
    live("refs", len(panel), len(panel), "done")
    n_un = sum(1 for j in all_junc if not j["annotated"])
    print(f"\nwrote reference_table.tsv ({len(all_rows)} rows), "
          f"junctions.tsv ({len(all_junc)} junctions, {n_un} unannotated), "
          f"panel_audit.tsv ({len(audit)} rows)")


if __name__ == "__main__":
    main()
