"""S22 stage 6 — the upstream pathway, measured rather than assumed.

The brief's third question asks whether the binding core relaxes in
lineages where the upstream PLC/IP3 pathway is reduced or absent.  A
lineage list taken from reading would be an assumption dressed as a scope,
so the list is derived here from the same reference proteomes the family
was swept over.

**What counts as a phosphoinositide-specific phospholipase C.**  A protein
carrying *both* halves of the catalytic TIM barrel — PF00387 (PI-PLC-X) and
PF00388 (PI-PLC-Y) — in one sequence.  Either half alone is not a PI-PLC:
the X box in particular turns up in unrelated proteins, and counting it
would report a pathway where there is none.  Both halves are counted
separately as well, so a proteome scoring one and not the other is visible
as such rather than silently absent.

**The search is not a second experiment.**  It is S20's design applied to a
different profile: one `hmmsearch` per profile per group DB at `-E 10`, the
primary call taken by filtering the same domtblout at E <= 1e-5, and hits
attributed to proteomes by `s20_lib.accession_to_upid` — the measured
attribution, because 48 eukaryotic taxids in this set carry two reference
proteomes each and a taxid key would credit both with either one's hits.

**The vertebrate group is in the sweep even though nobody doubts it has
PLC.**  A co-occurrence statement needs the cell where both are present as
much as the cell where one is missing, and a vertebrate proteome scoring
zero PI-PLC would be a failure of this instrument rather than a biological
result.  It is the positive control for the search.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import s22_lib as L  # noqa: E402
from scripts.s20_lib import (S20_DIR, group_paths, accession_to_upid,  # noqa: E402
                             log as s20log)
from src.utils.data_root import require_data_root  # noqa: E402

#: The two halves of the PI-PLC catalytic barrel, and why each is here.
PLC_PFAM = {
    "PF00387": "PI-PLC-Y — the C-terminal half of the catalytic TIM barrel",
    "PF00388": "PI-PLC-X — the N-terminal half of the catalytic TIM barrel",
}
EVALUE_PRIMARY = 1e-5
EVALUE_SWEEP = 10.0

#: The five eukaryotic reference-proteome sets.  `vertebrata` is S3's sweep
#: DB and the other four are S20's; together they partition Eukaryota, so
#: the denominators of the two sweeps add (S20's rule, reused).
EUK_GROUPS = ("vertebrata", "metazoa_nonvert", "fungi", "viridiplantae",
              "protista_other")
HOSTS = ("https://www.ebi.ac.uk/interpro/api",
         "https://www.ebi.ac.uk/interpro/wwwapi")


def pfam_dir() -> Path:
    d = L.OUT_DIR / "pfam"
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_pfam(pfam: str) -> Path:
    """One Pfam HMM from InterPro, committed beside the results."""
    dest = pfam_dir() / f"{pfam}.hmm"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    import requests
    last = ""
    for _ in range(3):
        for host in HOSTS:
            try:
                r = requests.get(f"{host}/entry/pfam/{pfam}?annotation=hmm",
                                 timeout=180)
                if r.status_code != 200:
                    last = f"HTTP {r.status_code} @{host}"
                    continue
                body = (gzip.decompress(r.content)
                        if r.content[:2] == b"\x1f\x8b" else r.content)
                text = body.decode()
                if not text.startswith("HMMER3"):
                    last = "not a HMMER3 file"
                    continue
                dest.write_text(text)
                return dest
            except Exception as exc:            # noqa: BLE001
                last = str(exc)
        time.sleep(5)
    raise RuntimeError(f"{pfam}: could not fetch HMM — {last}")


def group_db(group: str) -> Path:
    if group == "vertebrata":
        return require_data_root() / "proteomes" / "vertebrata_refprot.fasta"
    return group_paths(group)["db"]


def group_proteome_dir(group: str) -> Path:
    return require_data_root() / "proteomes" / group


def hmm_out_dir() -> Path:
    d = require_data_root() / "hmmer" / "s22"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_search(pfam: str, group: str, threads: int = 8) -> Path:
    """hmmsearch at E <= 10; the strict call filters this same file."""
    out = hmm_out_dir() / f"{pfam}_vs_{group}.domtblout"
    logf = out.with_suffix(".log")
    if out.exists() and out.stat().st_size > 0 and logf.exists() \
            and "[ok]" in logf.read_text()[-200:]:
        L.log(f"{pfam} vs {group}: cached")
        return out
    db = group_db(group)
    if not db.exists():
        raise FileNotFoundError(f"{group}: sweep DB missing at {db}")
    cmd = ["hmmsearch", "--cpu", str(threads), "-E", str(EVALUE_SWEEP),
           "--domtblout", str(out), str(fetch_pfam(pfam)), str(db)]
    L.log("running: " + " ".join(cmd[:6]) + f" ... {group}")
    t0 = time.time()
    with open(logf, "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"hmmsearch failed ({rc}) for {pfam} vs {group}")
    L.log(f"{pfam} vs {group}: {time.time() - t0:.0f}s")
    return out


def targets(domtbl: Path, evalue: float = EVALUE_PRIMARY) -> set[str]:
    """Accessions whose full-sequence E-value clears the primary threshold."""
    keep: set[str] = set()
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 13:
                continue
            try:
                if float(f[6]) <= evalue:
                    keep.add(f[0])
            except ValueError:
                continue
    return keep


def header_accession(name: str) -> str:
    """`sp|ACC|NAME` and `tr|ACC|NAME` both carry the accession in field 2."""
    parts = name.split("|")
    return parts[1] if len(parts) >= 3 else name


def vertebrate_manifest() -> list[dict]:
    rows = L.read_tsv(PROJECT_ROOT / "results" / "hmm_sweep" /
                      "proteome_manifest.tsv")
    return [{"upid": r["Proteome Id"], "organism": r["Organism"],
             "taxid": r["Organism Id"],
             "protein_count": int(r["Protein count"] or 0)} for r in rows]


def group_manifest(group: str) -> list[dict]:
    if group == "vertebrata":
        return vertebrate_manifest()
    rows = L.read_tsv(group_paths(group)["manifest"])
    return [{"upid": r["upid"], "organism": r["organism"],
             "taxid": r["taxid"],
             "protein_count": int(r["protein_count"] or 0)} for r in rows]


def _attribute(group: str, accs: set[str]) -> dict[str, str]:
    """acc -> upid.  S20's measured attribution, or a header pass for S3's DB."""
    if group != "vertebrata":
        return accession_to_upid(group, accs)
    cache = hmm_out_dir() / "acc2upid_vertebrata.tsv"
    if cache.exists():
        cached = {r["accession"]: r["upid"] for r in L.read_tsv(cache)}
        if accs <= set(cached):
            return {a: cached[a] for a in accs}
    out: dict[str, str] = {}
    files = sorted(group_proteome_dir("vertebrata").glob("*.fasta.gz"))
    s20log("s22_plc", f"vertebrata: attributing {len(accs)} hits across "
                      f"{len(files)} proteome files")
    for path in files:
        upid = path.name.split("_")[0]
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    a = header_accession(line[1:].split()[0])
                    if a in accs:
                        out.setdefault(a, upid)
    L.write_tsv(cache, [{"accession": a, "upid": u}
                        for a, u in sorted(out.items())],
                ["accession", "upid"])
    return out


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = [(g, False) for g in EUK_GROUPS]
    L.live("plc", steps)

    prov = [{"pfam": p, "role": note, "path": str(fetch_pfam(p).relative_to(PROJECT_ROOT)),
             "sha256": L.sha256(fetch_pfam(p))} for p, note in PLC_PFAM.items()]
    L.write_tsv(L.OUT_DIR / "plc_profiles.tsv", prov,
                ["pfam", "role", "path", "sha256"])

    rows: list[dict] = []
    for gi, group in enumerate(EUK_GROUPS):
        per_pfam: dict[str, set[str]] = {}
        for pfam in PLC_PFAM:
            per_pfam[pfam] = {header_accession(t)
                              for t in targets(run_search(pfam, group))}
        both = per_pfam["PF00387"] & per_pfam["PF00388"]
        allhits = set().union(*per_pfam.values())
        attr = _attribute(group, allhits)
        counts: dict[str, dict[str, int]] = {}
        for acc in allhits:
            upid = attr.get(acc)
            if upid is None:
                continue
            c = counts.setdefault(upid, {"x": 0, "y": 0, "both": 0})
            if acc in per_pfam["PF00388"]:
                c["x"] += 1
            if acc in per_pfam["PF00387"]:
                c["y"] += 1
            if acc in both:
                c["both"] += 1
        for m in group_manifest(group):
            c = counts.get(m["upid"], {"x": 0, "y": 0, "both": 0})
            rows.append({
                "group": group, "upid": m["upid"], "organism": m["organism"],
                "taxid": m["taxid"], "protein_count": m["protein_count"],
                "n_plc_x_only_domain": c["x"], "n_plc_y_only_domain": c["y"],
                "n_pi_plc": c["both"],
                "plc_status": "present" if c["both"] > 0 else "absent",
            })
        L.log(f"{group}: {sum(1 for r in rows if r['group'] == group and r['n_pi_plc'])}"
              f"/{sum(1 for r in rows if r['group'] == group)} proteomes carry a PI-PLC")
        steps[gi] = (group, True)
        L.live("plc", steps)

    L.write_tsv(L.OUT_DIR / "plc_repertoire.tsv", rows,
                ["group", "upid", "organism", "taxid", "protein_count",
                 "n_plc_x_only_domain", "n_plc_y_only_domain", "n_pi_plc",
                 "plc_status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
