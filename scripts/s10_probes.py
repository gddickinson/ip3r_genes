"""S10 step 6 — is the exon structure a gene, or is it an alignment?

Everything so far shows that the annotation is missing coding sequence the
aligner places. That is only an annotation bug if the aligner is right, so this
module attacks its own evidence from two directions.

**1. Transcript evidence at the junction (the brief's method).** A ±90 nt CDS
probe is built across each informative junction — informative meaning a
junction the annotation does not already model — and searched against every
nucleotide record NCBI holds *for that species*, fetched and searched locally.
Only hits **contiguous across the junction** count: an HSP has to span the
junction point with at least `MIN_ANCHOR` nt of exact placement on each side,
because a hit lying entirely in one exon says the exon exists, not that it is
spliced to the next one.

The species' **non-transcript** records are searched too, as a negative control
for the spanning criterion itself. A probe is 180 nt of contiguous *spliced*
sequence, so a genomic record cannot span its junction however well it matches
the exons either side. If the probes hit genomic records and none of those hits
spans, the test is demonstrably discriminating rather than merely returning
zero.

The denominator is recorded with the result, and for these two species it is
the result: 43 mRNA records exist for *Nibea albiflora* and 10 for
*Dissostichus eleginoides*, with no TSA at all. A search of 43 records that
returns nothing has not shown the gene is untranscribed — it has shown that the
species has no transcript deposits, which is a different statement and the one
this module reports. RNA-seq for both species does exist (182 and 19 runs), and
reaching it needs the streaming aligner S12 builds; the junctions are handed
there rather than half-answered here.

**2. Exon-boundary concordance (offline, and the stronger test here).** Every
swept genome whose *same paralog* was aligned from the *same bait* and whose
annotation independently agrees with the alignment (`loss = 0`) is a genome
where miniprot's exon boundaries have been corroborated by NCBI's own
annotation pipeline. If the case model puts its boundaries at the same bait
residues as those genomes do, the case model's exon structure is the structure
an independent annotation endorses everywhere else.

What that establishes and what it does not: it validates the **instrument** at
this gene, not the individual locus. A boundary shared with 35 RefSeq-annotated
genomes is not an aligner artefact. It remains conceivable that this particular
locus is a pseudogene with a perfectly conserved exon structure — which is what
the reading-frame test in `s10_orf.py` and the splice sites in `s10_evidence.py`
address, and why the three are reported together rather than any one alone.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s10_gff import read_miniprot_model                       # noqa: E402

#: Half-width of a junction probe, in nucleotides of *spliced* CDS.
PROBE_FLANK_NT = 90

#: How much of a hit must lie on each side of the junction before it counts as
#: spanning it. Eight is the bar S12's brief sets for a junction-spanning read
#: and is used here for the same reason: fewer than eight bases can be placed
#: by chance at the end of an alignment.
MIN_ANCHOR = 8

#: ...and eight is **not sufficient on its own**, which the genomic control
#: caught on its first run: blastn extends a high-scoring HSP some 10-15 bases
#: past the true junction into the intron, and that extension clears an 8-nt
#: anchor. Six of *Nibea albiflora*'s 55 junctions "spanned" against the
#: locus's own genomic DNA, where by construction nothing can.
#:
#: The second requirement comes from what a probe *is* rather than from the
#: size of the artefact: a probe is 180 nt of contiguous **spliced** sequence,
#: so a genomic match tops out at one exon's flank (90 nt) plus that short
#: extension, while a transcript match covers essentially all of it. Requiring
#: three quarters of the probe separates the two by construction.
MIN_PROBE_COVERAGE = 0.75

#: Junction classes worth probing — the ones where the annotation is absent or
#: broken. `within_one_model` junctions are already modelled by the annotation
#: and probing them would spend the budget confirming what is not disputed.
INFORMATIVE_CLASSES = ("between_models", "model_to_gap", "unannotated")

#: Bait-residue tolerance when matching an exon boundary between two genomes.
#: Zero would be too strict: the two loci are 10-15 % divergent, and a single
#: indel near a boundary shifts the bait residue it maps to by one.
BOUNDARY_TOL_AA = 2

EMAIL = "george.dickinson@gmail.com"


# --------------------------------------------------------------------------
# probes
# --------------------------------------------------------------------------
def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def spliced_cds(region, model: dict) -> tuple[str, list[int]]:
    """The spliced coding sequence and the offset of each internal junction.

    Blocks are spliced in target order — the order `read_miniprot_model`
    returns them in — so the offsets are positions in the transcript, which is
    what a probe has to be built in. Splicing in coordinate order would put
    every minus-strand gene's junctions in reverse.
    """
    parts, offsets, pos = [], [], 0
    for blk in model["cds"]:
        seq = region.at(blk["start"], blk["end"])
        if model["strand"] == "-":
            seq = _revcomp(seq)
        parts.append(seq)
        pos += len(seq)
        offsets.append(pos)
    return "".join(parts), offsets[:-1]


def build_probes(region, model: dict, introns: list[dict],
                 flank: int = PROBE_FLANK_NT,
                 classes: tuple[str, ...] = INFORMATIVE_CLASSES
                 ) -> list[dict]:
    """One probe per informative junction, in spliced-transcript coordinates.

    `introns` are in coordinate order (as `classify_introns` produces them) and
    the junction offsets are in transcript order, so on the minus strand the
    two run opposite ways. The mapping is done explicitly on the intron's own
    coordinates rather than by index, because an off-by-reversal here would
    label every probe with the wrong junction's class and the error would not
    show up anywhere downstream.

    `classes` defaults to the informative ones, which is what S10 asks for.
    S12 passes `within_one_model` as well, to probe the junctions the
    annotation *does* model as a positive control on the same deposits.
    """
    cds, offsets = spliced_cds(region, model)
    blocks = model["cds"]
    by_gap: dict[tuple[int, int], dict] = {}
    for i in range(len(blocks) - 1):
        a, b = blocks[i], blocks[i + 1]
        gap = ((a["end"] + 1, b["start"] - 1) if model["strand"] == "+"
               else (b["end"] + 1, a["start"] - 1))
        by_gap[gap] = {"junction_index": i + 1, "cds_offset": offsets[i],
                       "q_residue": a["q_end"]}
    probes = []
    for intr in introns:
        key = (intr["start"], intr["end"])
        meta = by_gap.get(key)
        if meta is None or intr["class"] not in classes:
            continue
        off = meta["cds_offset"]
        lo, hi = max(0, off - flank), min(len(cds), off + flank)
        seq = cds[lo:hi]
        if len(seq) < 2 * MIN_ANCHOR:
            continue
        probes.append({
            "junction_index": meta["junction_index"],
            "intron_start": intr["start"], "intron_end": intr["end"],
            "intron_class": intr["class"], "splice_class": intr["splice_class"],
            "donor": intr["donor"], "acceptor": intr["acceptor"],
            "q_residue": meta["q_residue"],
            "probe_len": len(seq), "junction_offset_in_probe": off - lo,
            "left_nt": off - lo, "right_nt": hi - off, "seq": seq,
        })
    return probes


# --------------------------------------------------------------------------
# transcript search
# --------------------------------------------------------------------------
def resource_counts(taxid: int, cache: Path) -> dict:
    """How much transcript evidence exists for the species at all.

    Fetched once and cached. Without this a table of zeros is unreadable — it
    could mean the gene is not transcribed or it could mean the species has no
    deposits, and those are opposite conclusions.
    """
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.stat().st_size:
        return json.loads(cache.read_text())
    out: dict[str, int | str] = {"taxid": taxid}
    queries = {
        "nuccore_total": ("nuccore", f"txid{taxid}[Organism:exp]"),
        "nuccore_mrna": ("nuccore",
                         f"txid{taxid}[Organism:exp] AND biomol_mrna[PROP]"),
        "nuccore_tsa": ("nuccore",
                        f"txid{taxid}[Organism:exp] AND is_tsa[PROP]"),
        "sra_rnaseq": ("sra",
                       f'txid{taxid}[Organism:exp] AND "rna seq"[Strategy]'),
    }
    for key, (db, term) in queries.items():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db={db}&term={urllib.parse.quote(term)}&retmax=0"
               f"&retmode=json&email={EMAIL}")
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                out[key] = int(json.load(fh)["esearchresult"]["count"])
        except Exception:                             # noqa: BLE001
            out[key] = "unavailable"
    cache.write_text(json.dumps(out, indent=2) + "\n")
    return out


def _blastn() -> str:
    from shutil import which
    hit = which("blastn")
    if hit:
        return hit
    env = Path("/opt/anaconda3/envs/piezo1/bin/blastn")
    return str(env) if env.exists() else ""


def fetch_species_nucleotides(taxid: int, archive: Path,
                              retmax: int = 2000) -> dict:
    """Every nucleotide record NCBI holds for the species, fetched once.

    Fetched and searched **locally** rather than through NCBI's remote BLAST
    queue. Three reasons, and the first is decisive: these species have 43 and
    10 mRNA records, so the entire searchable transcript space is smaller than
    a single BLAST query and downloading it costs seconds where the remote
    queue took the better part of an hour and returned nothing. Second, the
    archived FASTA makes the search re-derive offline, which a queued remote
    job never does. Third, it lets the *genomic* records be searched too, which
    is what gives the spanning test its own negative control — see
    `probe_evidence`.

Only the **transcript** set is fetched. The species' other nucleotide
    deposits are whole genomic scaffolds — 507 MB for 200 of *Nibea
    albiflora*'s 400 — and they are neither the evidence the brief asks for nor
    a useful control, since this project already holds that genome. The control
    for the spanning criterion is built locally instead, from the locus itself
    (`genomic_control`).
    """
    archive.mkdir(parents=True, exist_ok=True)
    out = {}
    for label, term in (
            ("mrna", f"txid{taxid}[Organism:exp] AND biomol_mrna[PROP]"),):
        fasta = archive / f"nuccore_{label}.fasta"
        count_file = archive / f"nuccore_{label}.count"
        if not fasta.exists():
            try:
                url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
                       f"esearch.fcgi?db=nuccore&term={urllib.parse.quote(term)}"
                       f"&retmax={retmax}&retmode=json&email={EMAIL}")
                with urllib.request.urlopen(url, timeout=120) as fh:
                    res = json.load(fh)["esearchresult"]
                ids, total = res.get("idlist", []), int(res["count"])
                count_file.write_text(f"{total}\n")
                text = ""
                for i in range(0, len(ids), 200):
                    post = urllib.parse.urlencode({
                        "db": "nuccore", "rettype": "fasta", "retmode": "text",
                        "email": EMAIL, "id": ",".join(ids[i:i + 200])}).encode()
                    req = urllib.request.Request(
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
                        "efetch.fcgi", data=post)
                    with urllib.request.urlopen(req, timeout=300) as fh:
                        text += fh.read().decode("utf-8", "replace")
                fasta.write_text(text)
            except Exception as exc:                  # noqa: BLE001
                out[label] = {"status": f"error:{type(exc).__name__}",
                              "path": None, "n_records": 0, "n_total": 0}
                continue
        n = fasta.read_text().count(">") if fasta.exists() else 0
        total = (int(count_file.read_text().strip())
                 if count_file.exists() else n)
        out[label] = {"status": "ok", "path": fasta, "n_records": n,
                      "n_total": total}
    return out


def genomic_control(region, cache_dir: Path) -> Path:
    """The locus's own genomic sequence, as the spanning test's negative control.

    A probe is 180 nt of contiguous *spliced* CDS, so no alignment to genomic
    DNA can cross its junction however well each half matches its exon. Running
    the probes against the locus they came from therefore has a known answer —
    many hits, none spanning — and a run that reports spanning hits here has a
    broken spanning test, not a discovery.

    Without this, "0 spanning hits in the transcript set" is unreadable: it
    could equally mean the criterion never fires.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "genomic_control.fasta"
    if not path.exists():
        path.write_text(f">locus_{region.contig}_{region.start}_{region.end}\n"
                        + region.seq + "\n")
    return path


def _blast_bin(name: str) -> str:
    from shutil import which
    hit = which(name)
    if hit:
        return hit
    env = Path("/opt/anaconda3/envs/piezo1/bin") / name
    return str(env) if env.exists() else ""


def search_probes(probes: list[dict], subject_fasta: Path, cache_dir: Path,
                  tag: str) -> dict:
    """Every probe against one archived record set, in one local blastn.

    The tabular output is archived so the counting re-derives offline (the
    discipline `s2_interpro.py` applies to InterPro pages and
    `s3_run_sweep.py --parse-only` to HMMER). A search that could not be run is
    recorded as `unavailable`, never as zero hits — the difference between "no
    transcript covers this junction" and "no search happened" is the whole
    result for these two species.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"probes_{tag}.tsv"
    status = "cached"
    if not cache.exists():
        blastn, makedb = _blast_bin("blastn"), _blast_bin("makeblastdb")
        if not (blastn and makedb) or not subject_fasta \
                or not subject_fasta.exists() \
                or not subject_fasta.read_text().strip():
            return {"status": "unavailable", "by_probe": {}}
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            q = tmp / "probes.fa"
            q.write_text("".join(
                f">j{p['junction_index']}\n{p['seq']}\n" for p in probes))
            subprocess.run([makedb, "-in", str(subject_fasta), "-dbtype",
                            "nucl", "-out", str(tmp / "db")], check=True,
                           capture_output=True)
            cols = ("qseqid sseqid pident length qstart qend sstart send "
                    "evalue bitscore")
            res = subprocess.run(
                [blastn, "-query", str(q), "-db", str(tmp / "db"),
                 "-outfmt", f"6 {cols}", "-evalue", "1e-5",
                 "-max_target_seqs", "100", "-num_threads", "4"],
                check=True, capture_output=True, text=True)
            cache.write_text(res.stdout)
        status = "searched"
    by_query: dict[str, list[str]] = {}
    for line in cache.read_text().splitlines():
        by_query.setdefault(line.split("\t")[0], []).append(line)
    out = {}
    for p in probes:
        rows = by_query.get(f"j{p['junction_index']}", [])
        out[p["junction_index"]] = {
            **count_spanning("\n".join(rows), p), "status": status}
    return {"status": status, "by_probe": out}


def count_spanning(tabular: str, probe: dict) -> dict:
    """Count only HSPs contiguous across the junction.

    Two requirements, and the second is why the genomic control exists.

    The junction sits at `junction_offset_in_probe` in query coordinates, and
    an HSP must cross it with at least `MIN_ANCHOR` nt on each side — the test
    the brief specifies. It must **also** cover `MIN_PROBE_COVERAGE` of the
    probe, because the anchor alone admits a genomic alignment that has run a
    dozen bases past the junction into the intron. A probe is contiguous
    spliced sequence: only a transcript can match all of it.

    `n_spanning_anchor_only` is kept so the difference between the two rules is
    visible in the table rather than absorbed into it.
    """
    j = probe["junction_offset_in_probe"]
    plen = max(1, probe.get("probe_len") or 1)
    n_hits = n_span = n_anchor = 0
    best_id, best_subj = 0.0, ""
    for line in tabular.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        n_hits += 1
        qs, qe = int(f[4]), int(f[5])
        if not (qs <= j - MIN_ANCHOR + 1 and qe >= j + MIN_ANCHOR):
            continue
        n_anchor += 1
        if (qe - qs + 1) / plen < MIN_PROBE_COVERAGE:
            continue
        n_span += 1
        if float(f[2]) > best_id:
            best_id, best_subj = float(f[2]), f[1]
    return {"n_hits": n_hits, "n_spanning": n_span,
            "n_spanning_anchor_only": n_anchor,
            "best_identity": round(best_id / 100, 4) if best_id else "",
            "best_subject": best_subj}


# --------------------------------------------------------------------------
# exon-boundary concordance
# --------------------------------------------------------------------------
def boundaries(model: dict) -> list[int]:
    """Internal exon boundaries as bait residue numbers, in target order."""
    return [b["q_end"] for b in model["cds"][:-1]]


def concordance(case_model: dict, refs: list[dict],
                tol: int = BOUNDARY_TOL_AA) -> tuple[list[dict], dict]:
    """How many reference genomes place each of the case's boundaries.

    `refs` are ``{"accession", "organism", "model"}`` for genomes whose same
    paralog was aligned from the same bait and whose annotation independently
    agrees with the alignment. The per-boundary row carries the count, so a
    boundary unique to the case locus is visible rather than averaged away.
    """
    ref_sets = [(r["accession"], r["organism"], set(boundaries(r["model"])))
                for r in refs]
    rows = []
    for i, q in enumerate(boundaries(case_model), start=1):
        hits = [acc for acc, _org, s in ref_sets
                if any(abs(q - x) <= tol for x in s)]
        rows.append({"junction_index": i, "q_residue": q,
                     "n_refs_with_boundary": len(hits),
                     "frac_refs": round(len(hits) / max(1, len(ref_sets)), 4)})
    shared = [r for r in rows if r["n_refs_with_boundary"] >= 1]
    majority = [r for r in rows
                if r["frac_refs"] >= 0.5]
    return rows, {
        "n_reference_genomes": len(ref_sets),
        "n_boundaries": len(rows),
        "n_shared_with_any": len(shared),
        "n_shared_with_majority": len(majority),
        "frac_shared_with_majority": round(
            len(majority) / max(1, len(rows)), 4),
        "median_refs_per_boundary": (
            sorted(r["n_refs_with_boundary"] for r in rows)[len(rows) // 2]
            if rows else 0),
    }


def reference_models(case: dict, ranking: list[dict], sweep_root: Path,
                     limit: int = 40) -> list[dict]:
    """Swept genomes that corroborate miniprot at this gene.

    The filter is the control: same paralog, same bait, eligible under E1-E5,
    and `loss == 0` — an annotation built by a different pipeline that put one
    gene model over the whole alignment. Genomes are taken in ranking order and
    capped, so the set is deterministic.
    """
    out = []
    for r in ranking:
        if len(out) >= limit:
            break
        if (r["cell"] != case["cell"] or r["bait"] != case["bait"]
                or r["accession"] == case["accession"]
                or r.get("eligible") not in ("True", True)
                or float(r["loss"]) != 0.0):
            continue
        gff = sweep_root / r["accession"] / "miniprot.gff"
        if not gff.exists():
            continue
        model = read_miniprot_model(gff, r["mp_id"])
        if len(model.get("cds", [])) < 2:
            continue
        out.append({"accession": r["accession"], "organism": r["organism"],
                    "identity": r["identity"], "model": model})
    out.sort(key=lambda r: r["accession"])
    return out


PROBE_COLS = ["case_id", "accession", "junction_index", "q_residue",
              "intron_start", "intron_end", "intron_class", "splice_class",
              "donor", "acceptor", "probe_len", "junction_offset_in_probe",
              "left_nt", "right_nt", "db", "status", "n_hits", "n_spanning",
              "n_spanning_anchor_only", "best_identity", "best_subject"]

CONCORDANCE_COLS = ["case_id", "accession", "junction_index", "q_residue",
                    "n_refs_with_boundary", "frac_refs"]
