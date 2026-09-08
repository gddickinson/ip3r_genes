"""S17 stage 2 — deep per-paralog orthologue sets from the S5 sweep models.

Why this stage exists. S6's msa_v2 is taxonomically broad (134 sequences from
protists to vertebrates) but **thin per paralog**: 15 ITPR1, 11 ITPR2, 18
ITPR3. Eleven sequences cannot support a per-site constraint score for a
2,701-residue protein — most columns would be scored from a dozen observations,
several of them mammals, and a single teleost insertion would read as a
conserved feature.

The depth already exists and has never been used. Every one of the 309 swept
genomes retains its `miniprot.gff`, and the `##STA` line before each gene model
is that locus's translation. Joined to `summary.json` by `mp_id`, that is
244-263 paralog-assigned full-length ITPR proteins per paralog — the sweep did
the orthology assignment already, with a bait panel whose family call is a
positive test (D14) and which S7's tree independently validated.

Four rules, each of which exists because of a specific finding upstream:

1. **One locus per genome x paralog.** S16 measured the teleost 3R duplicates:
   the second copy is an ohnolog, not an orthologue. Taking both would put a
   paralog pair inside a within-gene sample and read its divergence as
   tolerated variation. The sweep's cell already holds only that paralog's
   loci, so the rule is "the best locus of the cell", ranked on coverage first
   — a fragment that happens to be 95 % identical over 600 aa is worse evidence
   about a 2,700-residue protein than a full-length model at 60 %.
2. **A quality bar on the model** — bait coverage >= 0.80, identity >= 0.40,
   >= 1,800 aligned residues. The identity floor is the sweep's own
   `MIN_LOCUS_IDENTITY`, read back from `s5_sweep_lib` rather than retyped.
3. **Lesion-rich loci excluded** (S15a). A decaying or mis-assembled ORF's
   residues are drift, not constraint. S15a's screen is one-sided by
   construction, which is exactly the right direction here: it is used only to
   drop sequence from a constraint sample, never to call a gene dead.
4. **A shape screen the coverage bar cannot do.** Bait coverage is
   `aligned_aa / bait_len`, so a model that covers the bait *and* carries two
   thousand extra residues passes it. That is not hypothetical: S5b recorded
   that a large `-G` manufactures phantom loci in giant genomes, chaining
   shared-channel-module hits 23 Mbp apart into one locus. So every sequence is
   scored on the fraction of its **own** residues that land in columns the
   reference occupies, measured on the alignment itself, and the bar is put in
   the gap the distribution shows rather than chosen (`s15_calibrate_recon`'s
   rule, with both edges committed). One sequence fails it — *Lissotriton
   helveticus* ITPR2, 4,976 aa, 0.483 — against a curated minimum of 0.972.
5. **The reference sequence is the curated human protein, not a gene model.**
   All three references enter from UniProt whatever the filters did, so every
   S17 table is reported in the numbering ClinVar and UniProt use.

    python scripts/s17_orthologs.py             # harvest + align (cached)
    python scripts/s17_orthologs.py --force     # realign from scratch
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402

MIN_COVERAGE = 0.80
MIN_ALIGNED_AA = 1800


def min_identity() -> float:
    """The sweep's own locus identity floor, read back rather than retyped.

    Every genome summary records the value the sweep actually ran at, so the
    floor is recovered from the data instead of from a constant that could have
    drifted since (D13 applied to a threshold).
    """
    import s5_sweep_lib
    return float(s5_sweep_lib.MIN_LOCUS_IDENTITY)


def _best_locus(cell: dict, min_id: float) -> dict | None:
    ok = [x for x in cell.get("loci", [])
          if float(x.get("coverage", 0)) >= MIN_COVERAGE
          and float(x.get("identity", 0)) >= min_id
          and int(x.get("aligned_aa", 0)) >= MIN_ALIGNED_AA]
    if not ok:
        return None
    return max(ok, key=lambda x: (float(x["coverage"]), float(x["score"])))


def harvest(out_dir: Path) -> dict[str, dict[str, str]]:
    """Collect one protein per genome x paralog, plus the curated references."""
    min_id = min_identity()
    lesioned = L.lesioned_cells()
    integ = L.integrity_index()
    reps = {r["label"]: r for r in L.read_tsv(L.MSA_DIR / "representatives.tsv")}
    msa = L.read_fasta(L.MSA_DIR / "representatives.fasta")
    groups = L.msa_groups()

    sets: dict[str, dict[str, str]] = {p: {} for p in L.PARALOGS}
    manifest: list[dict] = []
    sta_cache: dict[str, dict[str, str]] = {}

    for summary, paralog, cell in L.iter_sweep_cells():
        acc = summary["accession"]
        locus = _best_locus(cell, min_id)
        if locus is None:
            manifest.append({
                "paralog": paralog, "accession": acc,
                "organism": summary.get("organism", ""),
                "vclass": summary.get("vclass", ""),
                "vorder": summary.get("vorder", ""),
                "mp_id": "", "bait": "", "coverage": cell.get("best_coverage", ""),
                "identity": cell.get("best_identity", ""), "aligned_aa": "",
                "frameshifts": "", "stop_codons": "", "integrity": "",
                "seq_len": 0, "source": "sweep_model", "kept": False,
                "reason": f"below_quality_bar({cell.get('status','')})",
                "name": "",
            })
            continue
        if acc not in sta_cache:
            sta_cache[acc] = L.sta_proteins(acc)
        seq = sta_cache[acc].get(locus.get("mp_id", ""), "")
        row = integ.get((acc, paralog, int(locus.get("locus_idx", 0)) if
                         str(locus.get("locus_idx", "")).isdigit() else 0), {})
        integrity = row.get("verdict", "unscored")
        is_lesioned = (acc, paralog) in lesioned
        keep = bool(seq) and not is_lesioned
        reason = ("ok" if keep else
                  "no_translation" if not seq else "s15a_elevated_lesions")
        org = summary.get("organism", acc).split("(")[0].strip().replace(" ", "_")
        name = f"{paralog}|{org}|{acc}"
        manifest.append({
            "paralog": paralog, "accession": acc,
            "organism": summary.get("organism", ""),
            "vclass": summary.get("vclass", ""),
            "vorder": summary.get("vorder", ""),
            "mp_id": locus.get("mp_id", ""), "bait": locus.get("bait", ""),
            "coverage": locus.get("coverage", ""),
            "identity": locus.get("identity", ""),
            "aligned_aa": locus.get("aligned_aa", ""),
            "frameshifts": locus.get("frameshifts", ""),
            "stop_codons": locus.get("stop_codons", ""),
            "integrity": integrity, "seq_len": len(seq),
            "source": "sweep_model", "kept": keep, "reason": reason,
            "name": name if keep else "",
        })
        if keep:
            sets[paralog][name] = seq.replace("*", "").replace("J", "X")

    # Curated msa_v2 proteins for species the sweep set does not cover. Rule 1
    # applies here too — one per species per paralog, longest wins — because
    # msa_v2 carries both teleost 3R copies for several species and taking both
    # would reintroduce the ohnolog pair rule 1 exists to keep out.
    by_species: dict[tuple[str, str], list[tuple[int, str, str]]] = {}
    for label, seq in msa.items():
        g = groups.get(label, "")
        if g not in L.PARALOGS or "GCF_" in label or "GCA_" in label:
            continue
        species = reps.get(label, {}).get("species") or label
        by_species.setdefault((g, species), []).append((len(seq), label, seq))
    for (g, species), cands in sorted(by_species.items()):
        cands.sort(reverse=True)
        for rank, (n, label, seq) in enumerate(cands):
            keep = rank == 0
            if keep:
                sets[g][label] = seq.replace("*", "").replace("J", "X")
            manifest.append({
                "paralog": g, "accession": "", "organism": species, "vclass": "",
                "vorder": "", "mp_id": "", "bait": "", "coverage": "",
                "identity": "", "aligned_aa": "", "frameshifts": "",
                "stop_codons": "", "integrity": "curated", "seq_len": n,
                "source": "msa_v2", "kept": keep,
                "reason": "ok" if keep else "same_species_shorter_copy",
                "name": label if keep else "",
            })

    for paralog, (_label, acc, desc) in L.REFERENCES.items():
        name = f"REF|{paralog}|{acc}"
        seq = L.uniprot_fasta(acc)
        sets[paralog][name] = seq
        manifest.append({
            "paralog": paralog, "accession": "", "organism": desc, "vclass": "",
            "vorder": "", "mp_id": "", "bait": "", "coverage": "",
            "identity": "", "aligned_aa": "", "frameshifts": "",
            "stop_codons": "", "integrity": "curated", "seq_len": len(seq),
            "source": "uniprot", "kept": True,
            "reason": "coordinate_reference", "name": name,
        })

    L.write_tsv(out_dir / "orthologs_manifest.tsv", manifest)
    for paralog, seqs in sets.items():
        L.write_fasta(out_dir / f"orthologs_{paralog}.fasta", seqs)
    return sets


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def align(out_dir: Path, paralog: str, force: bool = False) -> tuple[Path, dict]:
    """MAFFT the per-paralog set.

    `--thread 1` by decision (D24): S3 measured the same seeds giving profiles
    of 4,933 and 4,908 match states on consecutive `--thread -1` runs, and
    every per-site score S17 publishes is read off this file. `--auto` picks an
    FFT-NS-class strategy at this size; L-INS-i is what S6 used for 134
    sequences and is not affordable for 260, and the trade is deliberate — this
    alignment exists to score columns *within* one paralog, where sequences are
    60-99 % identical, not to resolve the deep relationships S6/S7 needed.

    A ragged output is a hard failure: S1 established that a silent MAFFT
    failure degrades to a star alignment with no other symptom.
    """
    src = out_dir / f"orthologs_{paralog}.fasta"
    dst = out_dir / f"aln_{paralog}.fasta"
    log = out_dir / f"mafft_{paralog}.log"
    if dst.exists() and not force:
        aln = L.read_fasta(dst)
        return dst, {"paralog": paralog, "cached": True, "n_seq": len(aln),
                     "n_col": len(next(iter(aln.values()))) if aln else 0,
                     "sha256_in": _sha(src), "sha256_out": _sha(dst)}
    t0 = time.time()
    with dst.open("w") as fo, log.open("w") as fl:
        r = subprocess.run(["mafft", "--auto", "--anysymbol", "--thread", "1",
                            str(src)], stdout=fo, stderr=fl)
    if r.returncode != 0:
        raise SystemExit(f"mafft failed for {paralog} (rc={r.returncode}); "
                         f"see {log}")
    aln = L.read_fasta(dst)
    widths = {len(s) for s in aln.values()}
    if len(widths) != 1:
        raise SystemExit(f"mafft returned a ragged alignment for {paralog}: "
                         f"{sorted(widths)[:5]}")
    dt = time.time() - t0
    print(f"[s17] aligned {paralog}: {len(aln)} seqs x {widths.pop()} cols "
          f"in {dt:.0f} s")
    return dst, {"paralog": paralog, "cached": False, "n_seq": len(aln),
                 "n_col": len(next(iter(aln.values()))), "runtime_s": round(dt, 1),
                 "sha256_in": _sha(src), "sha256_out": _sha(dst)}


#: A sequence must be curated (a UniProt/Ensembl protein from msa_v2, or a
#: coordinate reference) to count as a positive for the shape calibration: its
#: status as a real protein is external evidence the aligner did not produce.
def _is_curated(name: str) -> bool:
    return name.startswith("REF|") or not ("GCF_" in name or "GCA_" in name)


def shape_scores(aln: dict[str, str], ref_name: str) -> dict[str, float]:
    """Fraction of each sequence's own residues that sit in reference columns.

    A phantom locus — two genes the aligner chained into one model — covers the
    reference once and puts its extra thousands of residues in columns the
    reference does not occupy. Nothing in the harvest can see that, because
    bait coverage is a fraction of the *bait*.
    """
    ref = aln[ref_name]
    in_ref = {i for i, ch in enumerate(ref) if ch not in "-."}
    out = {}
    for name, seq in aln.items():
        own = sum(1 for ch in seq if ch not in "-.")
        hit = sum(1 for i, ch in enumerate(seq) if ch not in "-." and i in in_ref)
        out[name] = round(hit / own, 4) if own else 0.0
    return out


def calibrate_shape(scores: list[tuple[str, float, str]]) -> dict:
    """Put the bar in the distribution's own gap, and commit both edges.

    `s15_calibrate_recon.operating_point()`'s rule: Youden's threshold on a
    separated pair lands on the lowest positive, which is the most permissive
    bar the data allow, so the operating point is the **midpoint of the gap**
    and the two edges go into the table. Pooled across the three paralogs,
    because a per-paralog search would find a "gap" in any set of 250 numbers.
    Curated records are the positives: the bar must sit below every one of them
    or it is not a shape screen, it is a taxon filter.
    """
    vals = sorted(v for _n, v, _p in scores)
    cur_min = min((v for n, v, _p in scores if _is_curated(n)), default=1.0)
    below = [v for v in vals if v < cur_min]
    if not below:
        return {"bar": 0.0, "gap_low": "", "gap_high": "", "n_below_bar": 0,
                "curated_min": cur_min, "n_curated": sum(
                    1 for n, _v, _p in scores if _is_curated(n)),
                "rule": "no sequence scores below the curated minimum; "
                        "nothing to separate, so the screen drops nothing"}
    edges = below + [cur_min]
    gaps = [(b - a, a, b) for a, b in zip(edges, edges[1:])]
    width, lo, hi = max(gaps)
    bar = round((lo + hi) / 2, 4)
    return {"bar": bar, "gap_low": lo, "gap_high": hi, "gap_width": round(width, 4),
            "curated_min": cur_min,
            "n_curated": sum(1 for n, _v, _p in scores if _is_curated(n)),
            "n_below_bar": sum(1 for _n, v, _p in scores if v < bar),
            "rule": "midpoint of the largest gap below the lowest curated record"}


def screen(out_dir: Path, force: bool = False) -> dict:
    """Pass 1 align -> shape screen -> pass 2 align on the survivors."""
    scores: list[tuple[str, float, str]] = []
    for paralog in L.PARALOGS:
        src = out_dir / f"orthologs_{paralog}.fasta"
        dst = out_dir / f"aln_{paralog}_pass1.fasta"
        _run_mafft(src, dst, paralog, force)
        aln = L.read_fasta(dst)
        ref = f"REF|{paralog}|{L.REFERENCES[paralog][1]}"
        for name, v in shape_scores(aln, ref).items():
            scores.append((name, v, paralog))
    cal = calibrate_shape(scores)
    bar = cal["bar"]
    rows, dropped = [], {p: [] for p in L.PARALOGS}
    for name, v, paralog in sorted(scores, key=lambda x: (x[2], x[1])):
        keep = v >= bar
        if not keep:
            dropped[paralog].append(name)
        rows.append({"paralog": paralog, "name": name,
                     "curated": _is_curated(name), "frac_in_ref_columns": v,
                     "bar": bar, "kept": keep})
    L.write_tsv(out_dir / "ortholog_shape.tsv", rows)
    for paralog, names in dropped.items():
        if not names:
            continue
        seqs = L.read_fasta(out_dir / f"orthologs_{paralog}.fasta")
        for n in names:
            print(f"[s17] shape screen drops {paralog} {n} "
                  f"({len(seqs[n])} aa)")
            seqs.pop(n, None)
        L.write_fasta(out_dir / f"orthologs_{paralog}.fasta", seqs)
    cal["dropped"] = {p: v for p, v in dropped.items() if v}
    return cal


def _run_mafft(src: Path, dst: Path, tag: str, force: bool) -> float:
    log = dst.with_suffix(".log")
    if dst.exists() and not force:
        return 0.0
    t0 = time.time()
    with dst.open("w") as fo, log.open("w") as fl:
        r = subprocess.run(["mafft", "--auto", "--anysymbol", "--thread", "1",
                            str(src)], stdout=fo, stderr=fl)
    if r.returncode != 0:
        raise SystemExit(f"mafft failed for {tag} (rc={r.returncode}); see {log}")
    aln = L.read_fasta(dst)
    if len({len(s) for s in aln.values()}) != 1:
        raise SystemExit(f"mafft returned a ragged alignment for {tag}")
    dt = time.time() - t0
    print(f"[s17] aligned {tag}: {len(aln)} seqs x "
          f"{len(next(iter(aln.values())))} cols in {dt:.0f} s")
    return dt


def run(out_dir: Path = L.OUT_DIR, force: bool = False) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    sets = harvest(out_dir)
    for paralog in L.PARALOGS:
        print(f"[s17] {paralog}: {len(sets[paralog])} sequences harvested")
    cal = screen(out_dir, force=force)
    print(f"[s17] shape bar {cal['bar']} "
          f"(gap {cal.get('gap_low')} -> {cal.get('gap_high')}, "
          f"curated minimum {cal['curated_min']}); "
          f"{sum(len(v) for v in cal.get('dropped', {}).values())} dropped")

    mafft_v = subprocess.run(["mafft", "--version"], capture_output=True,
                             text=True).stderr.strip()
    stats = {"mafft_version": mafft_v, "min_coverage": MIN_COVERAGE,
             "min_identity": min_identity(), "min_aligned_aa": MIN_ALIGNED_AA,
             "shape_calibration": cal, "alignments": []}
    for paralog in L.PARALOGS:
        _p, s = align(out_dir, paralog, force=True if cal.get("dropped", {})
                      .get(paralog) else force)
        stats["alignments"].append(s)
    (out_dir / "align_stats.json").write_text(json.dumps(stats, indent=2))
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    run(args.out, force=args.force)


if __name__ == "__main__":
    main()
