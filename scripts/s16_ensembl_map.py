"""S16 stage `map` — the human paralogy map the 2R test runs on.

The 2R question is not "do the three ITPR neighbourhoods share genes". S8
already answered that with a flat **zero**: cross-paralog symbol Jaccard is
0.000 in every direction. After ~500 Myr the neighbours are no longer the
same genes — they are *paralogs* of each other. So the test needs a paralogy
map, and Ensembl Compara is the reference source for one.

The map is pulled from **BioMart, one query per chromosome, against a pinned
dated archive host**. Pinning is the reproducible choice (D24's discipline
applied to a database release): the rolling `www` host follows the release
cycle, so a rerun a month later would silently score the windows against a
different Compara tree. The archive's own registry is fetched and committed
alongside the map, so the release the numbers were made on is evidence in
the results directory rather than a sentence in a report.

BioMart rather than the REST homology endpoint, for the reason the PIEZO
project measured: the REST route needs one xref plus one homology call per
symbol — thousands of calls — while BioMart returns the entire human
paralogy set in 25 requests, and a *complete* map is what makes the
genome-wide block scan possible at all.

Outputs (`results/duplication/`): `paralogy_map.tsv`, `symbol_to_ensembl.tsv`,
`paralogy_map_notes.md`. Raw TSVs cached under
`<data_root>/raw_api/s16/biomart/`; a warm run makes no network call.
"""

from __future__ import annotations

import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402

#: The pinned archive. `www.ensembl.org` is a *fallback only* and the run
#: records which host answered, because a fallback that silently changes the
#: release is the failure this pin exists to prevent.
ARCHIVE_HOST = "https://jun2026.archive.ensembl.org/biomart/martservice"
FALLBACK_HOST = "https://www.ensembl.org/biomart/martservice"
DATASET = "hsapiens_gene_ensembl"
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

QUERY = (
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
    '<Query virtualSchemaName="default" formatter="TSV" header="1" '
    'uniqueRows="1" count="" datasetConfigVersion="0.6">'
    '<Dataset name="{dataset}" interface="default">'
    '<Filter name="chromosome_name" value="{chrom}"/>'
    '<Attribute name="ensembl_gene_id"/>'
    '<Attribute name="external_gene_name"/>'
    '<Attribute name="hsapiens_paralog_ensembl_gene"/>'
    '<Attribute name="hsapiens_paralog_associated_gene_name"/>'
    '<Attribute name="hsapiens_paralog_subtype"/>'
    '</Dataset></Query>'
)

FLANKS = L.PROJECT_ROOT / "results" / "synteny" / "flanks.tsv"


def biomart_dir() -> Path:
    d = L.cache_dir() / "biomart"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _get(url: str, timeout: int = 600) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "ip3r-genes-s16"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read().decode("utf-8", "replace")


def fetch_registry(out: Path, log=print) -> str:
    """The archive's mart registry — which Ensembl release this map is.

    Committed, because "release 116" is a fact about the run and a report
    that merely asserts it cannot be checked.
    """
    path = biomart_dir() / "registry.xml"
    if not path.exists() or path.stat().st_size == 0:
        path.write_text(_get(f"{ARCHIVE_HOST}?type=registry", timeout=120))
    text = path.read_text()
    (out / "biomart_registry.xml").write_text(text)
    release = ""
    for token in text.split():
        if token.startswith('database="ensembl_mart_'):
            release = token.split("_")[-1].strip('"')
            break
    log(f"[map] pinned host {urllib.parse.urlparse(ARCHIVE_HOST).netloc} "
        f"serves Ensembl Genes release {release or '?'}")
    return release


def fetch_chromosome(chrom: str, log=print) -> tuple[Path, str]:
    """One chromosome's paralogy table, cached on disk. (path, host)."""
    path = biomart_dir() / f"paralogy_chr{chrom}.tsv"
    if path.exists() and path.stat().st_size > 0:
        return path, "cache"
    body = urllib.parse.urlencode({"query": QUERY.format(dataset=DATASET,
                                                         chrom=chrom)})
    last: Exception | None = None
    for host in (ARCHIVE_HOST, FALLBACK_HOST):
        for attempt in range(3):
            try:
                text = _get(f"{host}?{body}")
                if text.lstrip().startswith("<"):
                    raise RuntimeError("HTML response (not a BioMart TSV)")
                if "Query ERROR" in text[:2000]:
                    raise RuntimeError(text[:300])
                path.write_text(text)
                netloc = urllib.parse.urlparse(host).netloc
                log(f"[map] chr{chrom}: {text.count(chr(10))} rows from {netloc}")
                return path, netloc
            except Exception as exc:                            # noqa: BLE001
                last = exc
                time.sleep(2 ** attempt)
    raise SystemExit(f"BioMart failed for chr{chrom}: {last}")


def load_biomart(log=print) -> tuple[dict, dict, list, set]:
    """(symbol -> ensg, ensg -> symbol, [(a, b, level), ...], hosts used)."""
    sym2ensg: dict[str, str] = {}
    ensg2sym: dict[str, str] = {}
    pairs: list[tuple[str, str, str]] = []
    hosts: set[str] = set()
    for chrom in CHROMOSOMES:
        path, host = fetch_chromosome(chrom, log=log)
        hosts.add(host)
        with open(path, newline="") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr, None)
            for row in rdr:
                if len(row) < 5:
                    continue
                gid, name, pid, pname, level = row[:5]
                if name:
                    sym2ensg.setdefault(name.upper(), gid)
                    ensg2sym.setdefault(gid, name)
                if pid:
                    if pname:
                        sym2ensg.setdefault(pname.upper(), pid)
                        ensg2sym.setdefault(pid, pname)
                    pairs.append((gid, pid, level))
    log(f"[map] BioMart: {len(ensg2sym)} named genes, {len(pairs)} directed "
        f"paralog relations")
    return sym2ensg, ensg2sym, pairs, hosts


def flank_symbols() -> set[str]:
    """Every informative flanking symbol S8 saw, in relaxed-key form."""
    syms: set[str] = set()
    if not FLANKS.exists():
        return syms
    with open(FLANKS, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("informative") == "1":
                s = (r.get("relaxed_key") or r.get("symbol") or "").strip()
                if s:
                    syms.add(s.upper())
    return syms


def run(out_dir: Path | None = None, extra_symbols: set[str] | None = None,
        log=print) -> dict:
    out = out_dir or L.out_dir()
    release = fetch_registry(out, log=log)
    sym2ensg_all, ensg2sym, pairs, hosts = load_biomart(log=log)

    wanted = set(flank_symbols()) | {s.upper() for s in (extra_symbols or ())}
    wanted = {s for s in wanted if s}
    resolved = {s: sym2ensg_all.get(s, "") for s in wanted}
    mapped = {s: g for s, g in resolved.items() if g}
    log(f"[map] {len(mapped)}/{len(wanted)} S8/window symbols resolved to a "
        f"human Ensembl gene ({len(mapped) / max(1, len(wanted)):.1%})")

    keep_ids = set(mapped.values())
    seen: set[tuple] = set()
    rows = []
    for a, b, level in pairs:
        key = (a, b) if a < b else (b, a)
        if key in seen or (a not in keep_ids and b not in keep_ids):
            continue
        seen.add(key)
        rows.append([key[0], key[1], ensg2sym.get(key[0], key[0]),
                     ensg2sym.get(key[1], key[1]), level])
    L.write_tsv(out / "paralogy_map.tsv",
                ["gene_a", "gene_b", "symbol_a", "symbol_b",
                 "duplication_level"], rows)
    L.write_tsv(out / "symbol_to_ensembl.tsv",
                ["symbol", "ensembl_gene", "status"],
                [[s, g, "wanted" if s in wanted else "background"]
                 for s, g in sorted(sym2ensg_all.items())])
    log(f"[map] {len(rows)} paralog pairs touching {len(keep_ids)} "
        f"neighbourhood genes; {len(sym2ensg_all)} symbols in the table")

    unmapped = sorted(s for s, g in resolved.items() if not g)
    notes = out / "paralogy_map_notes.md"
    with open(notes, "w") as fh:
        fh.write("# S16 paralogy map — provenance and coverage\n\n")
        fh.write(f"- source: Ensembl BioMart, dataset `{DATASET}`, host "
                 f"`{urllib.parse.urlparse(ARCHIVE_HOST).netloc}` — a dated "
                 f"archive, **pinned** so the map cannot change release "
                 f"under a rerun\n")
        fh.write(f"- release served by that host: **Ensembl Genes "
                 f"{release or 'unrecorded'}** (registry committed as "
                 f"`biomart_registry.xml`)\n")
        fh.write(f"- hosts that actually answered: "
                 f"{', '.join(sorted(hosts)) or 'none'}\n")
        fh.write("- attributes: gene id, gene name, paralogue gene id, "
                 "paralogue name, paralogue last common ancestor\n")
        fh.write(f"- genes named genome-wide: **{len(ensg2sym)}**; directed "
                 f"paralog relations: **{len(pairs)}**\n")
        fh.write(f"- symbols queried (S8 flanks + the human ITPR and RyR "
                 f"windows): **{len(wanted)}**, resolved **{len(mapped)}** "
                 f"({len(mapped) / max(1, len(wanted)):.1%})\n")
        fh.write(f"- undirected pairs committed: **{len(rows)}**\n\n")
        fh.write("An unmapped symbol is a lineage-specific name with no human "
                 "one-to-one (`LOC` ids were already dropped by S8's "
                 "`informative` filter). It can only *lower* the measured "
                 "paralogy between two neighbourhoods, so every link count in "
                 "this task is a floor.\n\n")
        fh.write(f"<details><summary>unmapped symbols ({len(unmapped)})"
                 f"</summary>\n\n")
        fh.write(", ".join(f"`{s}`" for s in unmapped[:600]))
        fh.write("\n\n</details>\n")
    log(f"[map] wrote {notes}")
    return {"release": release, "n_pairs": len(rows),
            "n_named": len(ensg2sym), "n_relations": len(pairs),
            "n_wanted": len(wanted), "n_mapped": len(mapped),
            "hosts": sorted(hosts)}


def full_pairs(gene_ids: set[str]) -> dict[str, list[tuple]]:
    """Adjacency for an arbitrary gene set, streamed from the BioMart cache.

    The committed `paralogy_map.tsv` is deliberately restricted to pairs
    touching a real neighbourhood so it stays a reviewable size. The block
    scan and the quartet test need paralogy for blocks found at run time, so
    those read the cache — same data, no extra network, nothing large added
    to the repo.
    """
    adj: dict[str, list[tuple]] = {}
    for chrom in CHROMOSOMES:
        path = biomart_dir() / f"paralogy_chr{chrom}.tsv"
        if not path.exists():
            continue
        with open(path, newline="") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr, None)
            for row in rdr:
                if len(row) < 5 or not row[2]:
                    continue
                a, b, level = row[0], row[2], row[4]
                if a in gene_ids and b in gene_ids:
                    adj.setdefault(a, []).append((b, level))
                    adj.setdefault(b, []).append((a, level))
    return adj


def load(out_dir: Path | None = None) -> dict:
    """Read the committed map back (no network)."""
    out = out_dir or L.OUT_DIR
    sym2ensg, ensg2sym, pairs = {}, {}, {}
    p = out / "symbol_to_ensembl.tsv"
    if p.exists():
        for r in L.read_tsv(p):
            if r["ensembl_gene"]:
                sym2ensg[r["symbol"]] = r["ensembl_gene"]
                ensg2sym.setdefault(r["ensembl_gene"], r["symbol"])
    p = out / "paralogy_map.tsv"
    if p.exists():
        for r in L.read_tsv(p):
            pairs.setdefault(r["gene_a"], []).append(
                (r["gene_b"], r["duplication_level"]))
            pairs.setdefault(r["gene_b"], []).append(
                (r["gene_a"], r["duplication_level"]))
            ensg2sym.setdefault(r["gene_a"], r["symbol_a"])
            ensg2sym.setdefault(r["gene_b"], r["symbol_b"])
    return {"sym2ensg": sym2ensg, "ensg2sym": ensg2sym, "pairs": pairs}


def main() -> None:
    import s16_paralogon
    run(extra_symbols=s16_paralogon.human_window_symbols())


if __name__ == "__main__":
    main()
