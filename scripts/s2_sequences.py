"""S2 step 4 — sequences for the census, and the representative set.

Two products, deliberately kept apart:

* **The seeded-space FASTA** — every sequence in the enumeration, ~45 MB.
  Bulk, so it is written under the data root and never into the repo. S3's
  profile sweep and S6's alignment both start from it, and having it on
  disk means neither has to re-walk the API.
* **The representative FASTA** — what an alignment can actually hold. One
  record per species per call (the longest, which for this family means the
  one least likely to be a fragment), with vertebrates thinned to one per
  taxonomic order because 707 vertebrate species of a three-paralog family
  say nothing 100 do not. Every non-vertebrate species is kept, since that
  is where the census's open questions are.

* **The core panel** — the longest ITPR and the longest RYR of every
  phylum in the census, ~66 sequences. This is the only FASTA that goes in
  the repo: a census-scale FASTA would be re-committed at every later
  census version, and the per-phylum panel is what a reader needs to check
  the separation for themselves. It doubles as the panel `s2_verify.py`
  runs the sequence-level D14 test on.

All three selections are recorded in `census_v2_representatives.tsv` with
the rule that picked each row, so the choice is auditable rather than
implicit in a FASTA header.
"""

from __future__ import annotations

import argparse
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.data_root import require_data_root      # noqa: E402
from scripts.s2_lib import (                           # noqa: E402
    OUT_DIR, ROOT, SEED_PFAMS, fetch, read_tsv, write_tsv,
)

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
REP_COLS = ["accession", "call", "confidence", "rule", "core", "gene",
            "species", "taxon_id", "group", "phylum", "class", "order",
            "length"]


def fasta_path() -> Path:
    d = require_data_root(ROOT) / "raw_api" / "uniprot"
    d.mkdir(parents=True, exist_ok=True)
    return d / "s2_seed_sweep.fasta"


def download_fasta(dest: Path) -> int:
    query = " OR ".join(f"(xref:pfam-{p})" for p in SEED_PFAMS)
    url = (f"{STREAM}?query={urllib.parse.quote(query)}&format=fasta")
    status, body, _ = fetch(url, timeout_s=600, tries=6)
    if status != 200:
        raise RuntimeError(f"UniProt FASTA stream returned HTTP {status}")
    dest.write_bytes(body)
    return body.count(b"\n>") + (1 if body.startswith(b">") else 0)


def load_fasta(path: Path) -> dict[str, str]:
    """FASTA → `{accession: sequence}`.

    Handles UniProt's `sp|ACC|NAME` headers and this project's own
    `ACC|GENE|Species|PARALOG` panel headers, which put the accession in
    the *first* field. Reading the wrong field is silent — it yields a
    dict keyed by gene symbols that simply never matches — so the two
    layouts are distinguished explicitly rather than by field count.
    """
    seqs: dict[str, str] = {}
    acc, buf = None, []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                head = line[1:].split()[0]
                parts = head.split("|")
                acc = (parts[1] if len(parts) > 2 and parts[0] in
                       ("sp", "tr") else parts[0])
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def choose(rows: list[dict]) -> list[dict]:
    """Longest per species per call; vertebrates thinned to one per order."""
    best: dict[tuple[str, str], dict] = {}
    for r in rows:
        if r["call"] not in ("ITPR", "RYR"):
            continue
        key = (r["taxon_id"], r["call"])
        cur = best.get(key)
        if cur is None or int(r["length"]) > int(cur["length"]):
            best[key] = r

    keep: list[dict] = []
    vert_best: dict[tuple[str, str], dict] = {}
    for r in best.values():
        if r["group"] == "Vertebrata":
            # One per order per call; an unranked record keeps its own slot
            # rather than being pooled into a single nameless bucket.
            okey = (r["order"] or f"unranked:{r['taxon_id']}", r["call"])
            cur = vert_best.get(okey)
            if cur is None or int(r["length"]) > int(cur["length"]):
                vert_best[okey] = r
        else:
            keep.append(dict(r, rule="longest per species"))
    keep += [dict(r, rule="longest per vertebrate order")
             for r in vert_best.values()]

    # The core panel: the longest of each call in each phylum. Records with
    # no phylum stated (several protist lineages) fall back to their group,
    # so an unranked lineage is still represented rather than dropped.
    core: dict[tuple[str, str], dict] = {}
    for r in keep:
        ckey = (r["phylum"] or f"[{r['group']}]", r["call"])
        cur = core.get(ckey)
        if cur is None or int(r["length"]) > int(cur["length"]):
            core[ckey] = r
    core_ids = {r["accession"] for r in core.values()}
    for r in keep:
        r["core"] = int(r["accession"] in core_ids)
    return sorted(keep, key=lambda r: (r["call"], r["group"], r["species"]))


def write_fasta(path: Path, rows: list[dict], seqs: dict[str, str]) -> int:
    n = 0
    with path.open("w", encoding="utf-8") as fh:
        for r in rows:
            s = seqs.get(r["accession"])
            if not s:
                continue
            label = f"{r['call']}_{r['species'].split('(')[0].strip().replace(' ', '_')}"
            fh.write(f">{r['accession']} {label} len={len(s)} "
                     f"conf={r['confidence']}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()

    dest = fasta_path()
    if not args.skip_download or not dest.exists():
        n = download_fasta(dest)
        print(f"[s2] {n} sequences → {dest} "
              f"({dest.stat().st_size / 1e6:.1f} MB, data root)", flush=True)
    seqs = load_fasta(dest)
    print(f"[s2] loaded {len(seqs)} sequences")

    rows = read_tsv(OUT_DIR / "census_v2.tsv")
    reps = choose(rows)
    missing = [r["accession"] for r in reps if r["accession"] not in seqs]
    write_tsv(OUT_DIR / "census_v2_representatives.tsv",
              [{c: r.get(c, "") for c in REP_COLS} for r in reps], REP_COLS)
    # The full representative set is bulk (~4.6 MB and re-cut at every later
    # census version), so it lives beside the sweep under the data root.
    bulk = dest.parent / "census_v2_representatives.faa"
    written = write_fasta(bulk, reps, seqs)
    core = [r for r in reps if r["core"]]
    n_core = write_fasta(OUT_DIR / "census_v2_core_panel.faa", core, seqs)
    print(f"[s2] core panel: {n_core} sequences "
          f"({len({r['phylum'] or r['group'] for r in core})} phyla) → repo")

    per_call: dict[str, int] = {}
    for r in reps:
        per_call[r["call"]] = per_call.get(r["call"], 0) + 1
    print(f"[s2] representatives: {len(reps)} rows "
          f"({', '.join(f'{k} {v}' for k, v in sorted(per_call.items()))}), "
          f"{written} written, {len(missing)} without a sequence")
    if missing:
        print(f"     missing: {', '.join(missing[:8])}"
              f"{' …' if len(missing) > 8 else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
