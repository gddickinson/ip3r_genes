"""S3 step 1 — the two seed alignments and the two profile HMMs.

Builds `itpr.hmm` and `ryr.hmm` from the committed manifests in
s3_seed_spec.py: pull each seed's sequence out of S2's archived seeded-space
FASTA, check it against the census row the manifest was written against,
align with MAFFT L-INS-i, and hmmbuild.

Every selection rule in s3_seed_spec.py is enforced here rather than
trusted: a seed whose census call, architecture count or length has moved
since the manifest was written aborts the build, and so does a species
contributing more than MAX_SPECIES_PER_PROFILE seeds to one profile.
MAFFT's version and return code are recorded, because S1 established that a
silent MAFFT failure falls back to a star alignment and nothing downstream
would notice.

Outputs (results/hmm_sweep/):
  itpr_seed.faa / ryr_seed.faa      the seed sets as aligned input
  itpr_seed.aln / ryr_seed.aln      the MAFFT L-INS-i alignments
  itpr.hmm / ryr.hmm                the profiles
  seed_manifest.tsv                 every seed, both profiles, with its rule
  seed_build_stats.json             tool versions, sizes, timings

Run:  python3 scripts/s3_build_seed.py
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_hmm_lib import (  # noqa: E402
    HMM_SWEEP_DIR, acc_key, census_index, load_census_v2, parse_uniprot_header,
    read_fasta, write_fasta, write_tsv,
)
from scripts.s3_seed_spec import (  # noqa: E402
    MAX_SPECIES_PER_PROFILE, MIN_SEED_LEN, PROFILES,
)
from src.utils.data_root import require_data_root  # noqa: E402

SEEDED_FASTA = "raw_api/uniprot/s2_seed_sweep.fasta"


def log(msg: str) -> None:
    print(f"[s3_seed] {msg}", flush=True)


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    log("$ " + " ".join(cmd[:6]) + (" …" if len(cmd) > 6 else ""))
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    log(f"  rc={proc.returncode} in {time.time() - t0:.0f}s")
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-3000:])
        raise SystemExit(f"command failed: {cmd[0]}")
    return proc


def load_seeded_sequences() -> dict[str, tuple[str, str]]:
    """accession → (sequence, full UniProt header) from S2's archive."""
    path = require_data_root() / SEEDED_FASTA
    if not path.exists():
        raise SystemExit(f"S2 seeded-space FASTA missing: {path}")
    out: dict[str, tuple[str, str]] = {}
    for header, seq in read_fasta(path).items():
        meta = parse_uniprot_header(header)
        out[acc_key(meta["accession"])] = (seq, header)
    log(f"{len(out)} sequences from S2's seeded-space archive")
    return out


def collect(profile: str, seeds, expect_call: str, seqs, census
            ) -> tuple[list[tuple[str, str]], list[dict]]:
    """Resolve one profile's manifest into (fasta items, manifest rows),
    enforcing every selection rule."""
    items: list[tuple[str, str]] = []
    rows: list[dict] = []
    species_count: Counter[str] = Counter()
    problems: list[str] = []

    for acc, clade, arch_expected, note in seeds:
        row = census.get(acc)
        if row is None:
            problems.append(f"{acc}: not in census v2")
            continue
        entry = seqs.get(acc)
        if entry is None:
            problems.append(f"{acc}: no sequence in the S2 archive")
            continue
        seq, header = entry
        arch = int(row["n_itpr_arch"])
        # R1 — the census call must be the profile's own call, except the
        # documented iplA exception, which is `unassigned` by design.
        if row["call"] != expect_call and not (
                profile == "itpr" and row["call"] == "unassigned"
                and "exception" in note):
            problems.append(f"{acc}: census call {row['call']!r} "
                            f"≠ {expect_call!r} and not the documented exception")
        # R3 — the manifest records what the seed carries; a drift is a stop.
        if arch != arch_expected:
            problems.append(f"{acc}: architecture {arch}/5, manifest says "
                            f"{arch_expected}/5 — census changed under the manifest")
        # R4 — length floor, with the note carrying any exception.
        if len(seq) < MIN_SEED_LEN:
            problems.append(f"{acc}: {len(seq)} aa < {MIN_SEED_LEN}")
        if len(seq) != int(row["length"]):
            problems.append(f"{acc}: sequence {len(seq)} aa ≠ census "
                            f"{row['length']} aa")
        species = row["species"].split(" (")[0]
        species_count[species] += 1
        sid = f"{acc}|{(row['gene'] or 'unnamed')}|{species.replace(' ', '_')}|{clade}"
        items.append((sid, seq))
        rows.append({
            "profile": profile, "id": sid, "accession": acc, "clade": clade,
            "call": row["call"], "gene": row["gene"], "species": species,
            "taxon_id": row["taxon_id"], "phylum": row["phylum"],
            "length": len(seq), "arch": f"{arch}/5",
            "reviewed": row["reviewed"], "note": note,
        })

    # R2 — no genome dominates a profile.
    for species, n in species_count.items():
        if n > MAX_SPECIES_PER_PROFILE:
            problems.append(f"{species}: {n} seeds in {profile} "
                            f"(> {MAX_SPECIES_PER_PROFILE})")
    if problems:
        for p in problems:
            sys.stderr.write(f"  SEED RULE VIOLATION  {p}\n")
        raise SystemExit(f"{len(problems)} seed-rule violations in {profile}")
    return items, rows


# MAFFT L-INS-i is run **single-threaded on purpose**. With `--thread -1`
# it is not reproducible: the same 22 RyR seeds aligned twice on this
# machine gave 8,510 and 8,468 columns, and the profile built from them
# 4,933 and 4,908 match states — because the iterative refinement stage
# combines partial results in whatever order the threads finish. A profile
# that changes when you rebuild it cannot be the profile a committed result
# was produced with, so the 6× slowdown (87 s vs 15 s, once) is the price
# of an artefact that can be regenerated. Two `--thread 1` runs are
# byte-identical.
MAFFT_THREADS = "1"


def mafft_cmd(infile: Path) -> list[str]:
    """L-INS-i, whichever binary name is installed, single-threaded."""
    linsi = shutil.which("mafft-linsi")
    if linsi:
        return [linsi, "--thread", MAFFT_THREADS, str(infile)]
    mafft = shutil.which("mafft")
    if not mafft:
        raise SystemExit("mafft not on PATH")
    return [mafft, "--localpair", "--maxiterate", "1000",
            "--thread", MAFFT_THREADS, str(infile)]


def build_profile(profile: str, items, hmm_name: str) -> dict:
    faa = HMM_SWEEP_DIR / f"{profile}_seed.faa"
    aln = HMM_SWEEP_DIR / f"{profile}_seed.aln"
    hmm = HMM_SWEEP_DIR / f"{profile}.hmm"
    write_fasta(faa, items)

    cmd = mafft_cmd(faa)
    proc = run(cmd)
    aln.write_text(proc.stdout)
    rows = read_fasta(aln)
    if len(rows) != len(items):
        raise SystemExit(f"{profile}: MAFFT returned {len(rows)} of "
                         f"{len(items)} sequences")
    width = len(next(iter(rows.values())))
    if any(len(s) != width for s in rows.values()):
        raise SystemExit(f"{profile}: ragged alignment — MAFFT did not align")
    log(f"{profile}: {len(rows)} seqs × {width} cols")

    hb = run(["hmmbuild", "--amino", "-n", hmm_name, str(hmm), str(aln)])
    digests = {"seed_faa_sha256": _sha256(faa),
               "alignment_sha256": _sha256(aln),
               "hmm_sha256": _sha256(hmm)}
    # hmmbuild's summary line: idx name nseq alen mlen W eff_nseq re/pos
    summary = [ln for ln in hb.stdout.splitlines()
               if ln.strip().startswith("1 ")]
    mlen = int(summary[0].split()[4]) if summary else None
    log(f"{profile}.hmm: {mlen} match states")
    return {"n_seeds": len(items), "alignment_cols": width,
            "match_states": mlen, "mafft_cmd": " ".join(cmd), **digests,
            "hmmbuild_summary": summary[0].split() if summary else [],
            "seed_lengths": sorted(len(s) for _, s in items)}


def _sha256(path: Path) -> str:
    """Digest of a build product, so a rebuild that drifts is visible in the
    committed stats rather than only in a match-state count."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(binary: str, flag: str = "--version") -> str:
    try:
        p = subprocess.run([binary, flag], capture_output=True, text=True)
        return (p.stdout or p.stderr).strip().splitlines()[0]
    except (OSError, IndexError):
        return "unknown"


def main() -> int:
    t0 = time.time()
    HMM_SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    census = census_index(load_census_v2())
    seqs = load_seeded_sequences()

    manifest: list[dict] = []
    stats: dict = {"profiles": {}}
    for profile, spec in PROFILES.items():
        items, rows = collect(profile, spec["seeds"], spec["expect_call"],
                              seqs, census)
        by_clade = Counter(r["clade"] for r in rows)
        log(f"{profile}: {len(items)} seeds {dict(by_clade)}")
        stats["profiles"][profile] = build_profile(
            profile, items, spec["hmm_name"])
        stats["profiles"][profile]["by_clade"] = dict(by_clade)
        manifest.extend(rows)

    write_tsv(HMM_SWEEP_DIR / "seed_manifest.tsv",
              ["profile", "id", "accession", "clade", "call", "gene",
               "species", "taxon_id", "phylum", "length", "arch",
               "reviewed", "note"], manifest)
    stats.update({
        "mafft_version": tool_version("mafft"),
        "hmmbuild_version": tool_version("hmmbuild", "-h").strip(),
        "census_rows": len(census),
        "elapsed_s": round(time.time() - t0, 1),
    })
    (HMM_SWEEP_DIR / "seed_build_stats.json").write_text(
        json.dumps(stats, indent=1))
    log(f"both profiles built in {stats['elapsed_s']}s → {HMM_SWEEP_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
