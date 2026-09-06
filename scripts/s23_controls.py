"""s23_controls.py — the positive-control baits for the non-vertebrate sweep.

**The control S5 used does not transfer, and this is the module that says so.**

In S5 every genome carries a RyR, so a genome where the RyR baits find nothing
has an assembly or a pipeline problem rather than a biological result. That
argument is what lets S5b write "the RyR positive control fired in every one
of the 309, so no genome is excluded on control grounds".

Outside the metazoa it collapses. S20 scored 6,928 reference proteomes against
`ryr.hmm` and found architecture-level RyR in exactly two non-metazoan ones
(*Salpingoeca*, *Capsaspora*). A land plant with no RyR locus is the correct
answer, so RyR cannot be that genome's proof-of-search — and without a proof
of search, "no ITPR in *Arabidopsis*" is indistinguishable from "the sweep did
not run properly on *Arabidopsis*". Every negative claim in S23 rests on
closing that gap.

**The control is the MIR-domain sharer.** PF02815 (MIR) is one of the family's
own four diagnostic signatures and the one every eukaryote carries on
something else — the protein O-mannosyltransferase / dolichyl-phosphate-mannose
mannosyltransferase family. S1's decoy panel was built around exactly these
proteins (POMT1/2), and S20 used them as its in-search positive control at
proteome level: at E <= 10 the four family Pfams return **0** PF08709 matches
in land plants and **633** PF02815 ones, **0** in Dikarya and **4,376**. The
instrument demonstrably works there and finds everything except the receptor.
This module takes that same control one level down, to the assembly.

Three properties make it the right control rather than a convenient one:

  * **It is in the family's own signature set.** A control drawn from an
    unrelated gene would show that *some* search works on the assembly; this
    one shows that a search for one of the four domains that define the family
    works on it.
  * **It is drawn from the clade the claim is about.** A chlorophyte
    mannosyltransferase is not a control for *Arabidopsis*, so the bait for
    each control clade comes from that clade — `s23_scope.control_clades()`
    enumerates them from S20's own presence table.
  * **It is simultaneously the sharpest decoy available.** A MIR bait must be
    assigned to *neither* family by the chimera screen, and at genome level a
    MIR locus must never be called ITPR. D14 gets a negative control in every
    genome for free, which S5's panel never had.

Sequences and hits come from committed or archived files only: the PF02815
model is `results/s20_sweep/pfam/PF02815.hmm`, the search runs against the
archived per-group reference-proteome DBs under the data root, and every hit's
organism and taxid are read out of the domtblout's own description field
(UniProt writes `OS=` and `OX=` there), so nothing here needs the network.
"""

from __future__ import annotations

import gzip
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.utils.data_root import require_data_root      # noqa: E402
import s23_bait_spec as spec                           # noqa: E402
import s23_scope as scope                              # noqa: E402

PFAM_MIR = PROJECT_ROOT / "results" / "s20_sweep" / "pfam" / "PF02815.hmm"

#: The S20 group whose archived proteome DB each control clade is searched in.
#: Straight from S20's own partition of Eukaryota — these are the same four
#: DBs the ITPR sweep ran against, so a control and the negative it controls
#: for come out of one database.
GROUP_DBS = ("metazoa_nonvert", "fungi", "viridiplantae", "protista_other")

#: E-value for the control search. Strict, because a control bait has to be a
#: real MIR protein: this is bait selection, not the relaxed sensitivity
#: sweep S20 ran to make its negative claims at.
CONTROL_E = 1e-5

_OS_RE = re.compile(r"\bOS=(.*?)(?:\s+[A-Z][A-Z]=|$)")
_OX_RE = re.compile(r"\bOX=(\d+)")
_GN_RE = re.compile(r"\bGN=(\S+)")


def hmmer_dir() -> Path:
    d = require_data_root() / "hmmer" / "s23_controls"
    d.mkdir(parents=True, exist_ok=True)
    return d


def proteome_db(group: str) -> Path:
    return require_data_root() / "proteomes" / f"{group}_refprot.fasta"


def run_search(group: str, threads: int = 6, force: bool = False) -> Path:
    """PF02815 over one group's archived proteome DB. Cached on the output."""
    out = hmmer_dir() / f"PF02815_{group}.domtblout"
    if out.exists() and out.stat().st_size > 0 and not force:
        return out
    db = proteome_db(group)
    if not db.exists():
        raise SystemExit(
            f"archived proteome DB missing: {db}\n"
            "  rerun scripts/s20_fetch.py --groups " + group)
    log = out.with_suffix(".log")
    tmp = out.with_suffix(".domtblout.partial")
    cmd = ["hmmsearch", "--cpu", str(threads), "-E", str(CONTROL_E),
           "--domtblout", str(tmp), str(PFAM_MIR), str(db)]
    with open(log, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"hmmsearch PF02815 on {group} failed "
                           f"(rc={proc.returncode}); see {log}")
    tmp.replace(out)
    return out


def parse_hits(domtbl: Path, group: str) -> dict[str, dict]:
    """accession -> best row, with the organism and taxid off the description.

    UniProt's FASTA header carries `OS=` and `OX=`, and HMMER preserves
    everything after the first token in the description column — so the hit's
    taxonomy comes out of the search's own output and needs no second lookup.
    That matters here: S20 measured that 48 taxids carry two reference
    proteomes each, so proteome attribution is not a taxid lookup, but the
    reverse direction (this one) is exact.
    """
    best: dict[str, dict] = {}
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split(None, 22)
            if len(f) < 23:
                continue
            target, tlen, score = f[0], int(f[2]), float(f[7])
            desc = f[22]
            parts = target.split("|")
            acc = parts[1] if len(parts) >= 3 else target
            os_m, ox_m, gn_m = (_OS_RE.search(desc), _OX_RE.search(desc),
                                _GN_RE.search(desc))
            row = {"accession": acc, "target": target, "length": tlen,
                   "score": score, "group": group,
                   "organism": (os_m.group(1).strip() if os_m else ""),
                   "taxid": (ox_m.group(1) if ox_m else ""),
                   "gene": (gn_m.group(1) if gn_m else ""),
                   "protein_name": desc.split(" OS=")[0].strip()}
            cur = best.get(acc)
            if cur is None or row["score"] > cur["score"]:
                best[acc] = row
    return best


def clade_coverage(presence: list[dict], taxonomy: dict[str, dict],
                   hits_by_group: dict[str, dict[str, dict]]
                   ) -> dict[tuple, float]:
    """(rank, clade) -> fraction of the clade's swept proteomes with a MIR hit.

    How much of the clade the control actually speaks for. A control drawn
    from the one apicomplexan proteome of 36 that carries a MIR protein is a
    different kind of evidence from one drawn from a class where every
    proteome has one, and `s23_bait_spec.control_strength` uses this to say
    which is which instead of presenting both as "controlled".
    """
    swept: dict[tuple, set] = defaultdict(set)
    for r in presence:
        tax = taxonomy.get(r["taxid"], {})
        for rank in ("class", "phylum"):
            if tax.get(rank):
                swept[(rank, tax[rank])].add(r["taxid"])
    with_hit: dict[tuple, set] = defaultdict(set)
    for hits in hits_by_group.values():
        for row in hits.values():
            tax = taxonomy.get(row["taxid"], {})
            for rank in ("class", "phylum"):
                if tax.get(rank):
                    with_hit[(rank, tax[rank])].add(row["taxid"])
    return {k: (len(with_hit.get(k, ())) / len(v) if v else 0.0)
            for k, v in swept.items()}


def select(clades: list[dict], taxonomy: dict[str, dict],
           hits_by_group: dict[str, dict[str, dict]],
           coverage: dict[tuple, float] | None = None) -> tuple[list, list]:
    """One control bait per control clade. Returns (chosen, unfilled).

    Ranked by profile score among the hits that sit in the MIR bait band and
    whose own protein name says mannosyltransferase. The name filter is a
    selection convenience, not the call: the chimera screen decides family
    membership, and a control bait's job is to be assigned to *neither*
    family. Where no named record qualifies the rule falls back to the best
    unnamed PF02815 protein in the band and records that it did.
    """
    by_clade: dict[tuple, list[dict]] = defaultdict(list)
    for group, hits in hits_by_group.items():
        for row in hits.values():
            tax = taxonomy.get(row["taxid"])
            if not tax:
                continue
            ok, why = spec.passes_shape(spec.CONTROL_MIR, row["length"], "")
            if not ok:
                continue
            for rank in ("class", "phylum"):
                name = (tax.get(rank) or "").strip()
                if name:
                    by_clade[(rank, name)].append(
                        dict(row, rank=rank, clade=name,
                             named=int(spec.is_mir_sharer_name(
                                 row["protein_name"])),
                             lineage_kingdom=tax.get("kingdom", ""),
                             lineage_phylum=tax.get("phylum", "")))
    chosen, unfilled = [], []
    for c in clades:
        pool = by_clade.get((c["rank"], c["clade"]), [])
        if not pool:
            unfilled.append(dict(c, why="no PF02815 protein in the MIR bait "
                                        "band in any swept proteome of this "
                                        "clade"))
            continue
        # Prefer a *strong* control: a named, full-length mannosyltransferase
        # if the clade has one, and only then fall back on length/score. The
        # first build ranked on name-then-score inside a 600 aa floor, which
        # is the same preference expressed as an exclusion — and it excluded
        # both apicomplexan classes rather than reporting a weak control.
        pool.sort(key=lambda r: (
            -(r["named"] and r["length"] >= spec.CONTROL_STRONG_MIN_AA),
            -r["named"], -r["score"], r["accession"]))
        pick = pool[0]
        frac = (coverage or {}).get((c["rank"], c["clade"]), 0.0)
        strength, why = spec.control_strength(pick["named"], pick["length"],
                                              frac)
        chosen.append(dict(
            pick, control_for=c["clade"], control_rank=c["rank"],
            swept=c["swept"], candidates=len(pool),
            clade_frac=round(frac, 4), strength=strength, strength_note=why,
            reason=(f"best PF02815 protein in {c['clade']} "
                    f"({pick['score']:.0f} bits, {pick['length']} aa"
                    + (", named a mannosyltransferase" if pick["named"]
                       else ", unnamed")
                    + f"); the clade S20 swept {c['swept']} proteomes of "
                      f"for 0 ITPR; control is {strength} — {why}")))
    return chosen, unfilled


def _seq_cache() -> Path:
    return hmmer_dir() / "control_sequences.faa"


def _read_cache() -> dict[str, str]:
    path = _seq_cache()
    if not path.exists():
        return {}
    out, name, buf = {}, None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def _write_cache(seqs: dict[str, str]) -> None:
    with _seq_cache().open("w") as fh:
        for k, v in sorted(seqs.items()):
            fh.write(f">{k}\n")
            for i in range(0, len(v), 60):
                fh.write(v[i:i + 60] + "\n")


def fetch_sequences(rows: list[dict]) -> dict[str, str]:
    """Pull each control bait's sequence out of the archived proteome DBs.

    Cached under the data root. Without the cache every rebuild streams the
    four concatenated reference-proteome DBs — 24.5 GB — to recover 28
    sequences, which on an external drive is the single slowest step in the
    panel build and is repeated whenever any unrelated bait rule changes.
    """
    cached = _read_cache()
    want: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        if r["accession"] not in cached:
            want[r["group"]].add(r["accession"])
    out: dict[str, str] = {r["accession"]: cached[r["accession"]]
                           for r in rows if r["accession"] in cached}
    for group, accs in want.items():
        path = proteome_db(group)
        opener = gzip.open if path.suffix == ".gz" else open
        found, name, buf = 0, None, []
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in accs:
                        out[name] = "".join(buf)
                        found += 1
                        if found == len(accs):
                            break      # every bait in this group is in hand
                    parts = line[1:].split("|")
                    name = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if name in accs and name not in out:
                out[name] = "".join(buf)
    if want:
        _write_cache({**cached, **out})
    return out


def build(threads: int = 6, force: bool = False) -> tuple[list, list, dict]:
    """(control rows, unfilled clades, accession -> sequence)."""
    presence = scope.load_presence()
    taxonomy = scope.load_taxonomy()
    clades = scope.control_clades(presence, taxonomy)
    hits = {g: parse_hits(run_search(g, threads, force), g)
            for g in GROUP_DBS}
    coverage = clade_coverage(presence, taxonomy, hits)
    chosen, unfilled = select(clades, taxonomy, hits, coverage)
    seqs = fetch_sequences(chosen)
    missing = [r["accession"] for r in chosen if r["accession"] not in seqs]
    if missing:
        raise SystemExit(
            f"{len(missing)} control baits have no sequence in the archived "
            f"proteome DBs: {missing[:5]}\n"
            "  the domtblout and the FASTA are out of step — rerun "
            "s23_controls with --force")
    return chosen, unfilled, seqs


if __name__ == "__main__":
    rows, unfilled, seqs = build()
    print(f"{len(rows)} control baits, {len(unfilled)} clades unfilled")
    for r in rows:
        print(f"  {r['strength']:6s} {r['control_for']:20s} "
              f"{r['accession']:12s} {r['length']:5d} aa  "
              f"{r['score']:7.1f} bits  clade {r['clade_frac']:.0%}  "
              f"{r['organism'][:28]}")
    for u in unfilled:
        print(f"  ! {u['clade']:22s} {u['why']}")
