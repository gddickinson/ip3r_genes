"""S6 shared plumbing — paths, census loading, sequence resolution.

Three things live here because both halves of S6 (`s6_select_reps.py`
picking the set, `s6_msa.py` aligning it) need them and neither should own
them:

1. **The census loader**, which normalises the one column S23 left in two
   spellings. Census v6 stacks three merges deep, and `tax_group` arrives
   title-case from S20's lineage pass (`Viridiplantae`) on v5 rows and
   lower-case from `s23_census_v6.py` on the 183 genomic models
   (`viridiplantae`). Counted naively that is two groups, and a selection
   rule stated per group would silently give the plant grade two slots.

2. **Sequence resolution across five stores.** A census v6 accession is a
   UniProt accession from the S2 archive, a UniProt accession that exists
   only inside a reference proteome (the S3/S20 sweep rows), or a
   `<assembly>|<cell>|<contig>:<span>` gene-model id from S5 or S23. The
   first four stores are small files that index in a second; the fifth is
   38 GB of proteome FASTA, so it is streamed once for whatever the first
   four could not answer and the answers are cached under the data root.
   `resolve()` therefore reruns offline after the first pass.

3. **The label**, which becomes an MSA row name, an IQ-TREE tip name and a
   figure tick. It has to survive all three, so it is sanitised to
   `[A-Za-z0-9_.]` and built as `<group>_<Species>_<accession-tail>`.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MSA_DIR = ROOT / "results" / "msa_v2"
CENSUS_V6 = ROOT / "results" / "census_v6" / "census_v6.tsv"

#: The stores a census accession can be resolved from, in the order they are
#: tried. Small-and-indexable first; the proteome DBs are the streamed
#: fallback and are handled separately by `stream_proteomes()`.
_LOCAL_FASTAS = [
    ("s2_seeded", "raw_api/uniprot/s2_seed_sweep.fasta"),
    ("census_v2_rep", "raw_api/uniprot/census_v2_representatives.faa"),
    ("s5_models", "census_v4/genome_models.faa"),
    ("s23_models", "census_v6/genome_models_s23.faa"),
]

#: The S3/S20 sweep DBs, and which census `tax_group` sends a lookup to
#: which file. A miss in the guessed one falls through to all of them.
_PROTEOME_DBS = {
    "Vertebrata": "proteomes/vertebrata_refprot.fasta",
    "Metazoa (non-vertebrate)": "proteomes/metazoa_nonvert_refprot.fasta",
    "Viridiplantae": "proteomes/viridiplantae_refprot.fasta",
    "Fungi": "proteomes/fungi_refprot.fasta",
    "SAR": "proteomes/protista_other_refprot.fasta",
    "Discoba": "proteomes/protista_other_refprot.fasta",
    "Amoebozoa": "proteomes/protista_other_refprot.fasta",
    "Eukaryota (other)": "proteomes/protista_other_refprot.fasta",
}

_SAFE = re.compile(r"[^A-Za-z0-9_.]+")


# --------------------------------------------------------------- data root

def data_root() -> Path:
    from src.utils.data_root import require_data_root
    return require_data_root()


# ------------------------------------------------------------------ FASTA

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs: dict[str, str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as out:
        for k, v in seqs.items():
            out.write(f">{k}\n")
            for i in range(0, len(v), 60):
                out.write(v[i:i + 60] + "\n")
    return path


def uniprot_id(header: str) -> str:
    """`sp|Q14643|ITPR1_HUMAN ...` -> `Q14643`; anything else -> first word."""
    first = header.split()[0]
    parts = first.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr"):
        return parts[1]
    return first


# ----------------------------------------------------------------- census

def normalise_group(value: str) -> str:
    """One spelling per group (see module docstring, point 1)."""
    v = (value or "").strip()
    if not v:
        return "Eukaryota (other)"
    canon = {
        "vertebrata": "Vertebrata",
        "metazoa (non-vertebrate)": "Metazoa (non-vertebrate)",
        "metazoa": "Metazoa (non-vertebrate)",
        "viridiplantae": "Viridiplantae",
        "fungi": "Fungi",
        "sar": "SAR",
        "discoba": "Discoba",
        "amoebozoa": "Amoebozoa",
        "eukaryota (other)": "Eukaryota (other)",
        "other": "Eukaryota (other)",
    }
    return canon.get(v.lower(), v)


def load_census(path: Path = CENSUS_V6) -> list[dict]:
    """Census v6 rows with `group` normalised into `grp`."""
    if not path.exists():
        raise SystemExit(f"census not found: {path}\n  run scripts/s23_census_v6.py")
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["grp"] = normalise_group(r.get("tax_group") or r.get("group") or "")
    return rows


def as_int(value, default=0) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def as_float(value, default=0.0) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


# --------------------------------------------------- sequence resolution

class SequenceStore:
    """Resolve census accessions to sequences across the five stores.

    `resolve(accessions, groups)` returns `(seqs, source, missing)`. The
    proteome streaming pass runs only for what the local files could not
    answer, and everything it finds is appended to a cache FASTA under the
    data root, so a second run of the selector never touches the 38 GB.
    """

    def __init__(self, root: Path | None = None):
        self.root = root or data_root()
        self.cache = self.root / "raw_api" / "uniprot" / "s6_resolved.faa"
        self._local: dict[str, tuple[str, str]] = {}   # acc -> (store, seq)
        self._loaded = False

    def _load_local(self) -> None:
        if self._loaded:
            return
        for store, rel in _LOCAL_FASTAS:
            p = self.root / rel
            if not p.exists():
                continue
            for header, seq in _iter_fasta(p):
                acc = uniprot_id(header)
                self._local.setdefault(acc, (store, seq))
        if self.cache.exists():
            for header, seq in _iter_fasta(self.cache):
                self._local.setdefault(uniprot_id(header), ("s6_cache", seq))
        self._loaded = True

    def resolve(self, wanted: dict[str, str], verbose: bool = True
                ) -> tuple[dict[str, str], dict[str, str], list[str]]:
        """`wanted` maps accession -> census group (for DB guessing)."""
        self._load_local()
        seqs: dict[str, str] = {}
        source: dict[str, str] = {}
        missing = []
        for acc in wanted:
            hit = self._local.get(acc)
            if hit:
                source[acc], seqs[acc] = hit
            else:
                missing.append(acc)
        if missing:
            found = self._stream_proteomes(
                {a: wanted[a] for a in missing}, verbose=verbose)
            for acc, seq in found.items():
                seqs[acc] = seq
                source[acc] = "proteome_db"
            self._append_cache(found)
            missing = [a for a in missing if a not in found]
        return seqs, source, missing

    def _stream_proteomes(self, wanted: dict[str, str], verbose: bool
                          ) -> dict[str, str]:
        """One pass per proteome file that any missing accession points at."""
        files: dict[str, set[str]] = {}
        for acc, grp in wanted.items():
            rel = _PROTEOME_DBS.get(grp)
            files.setdefault(rel or "", set()).add(acc)
        # a group with no mapped DB, or a miss in its own DB, is looked for
        # everywhere -- the sweep DBs partition Eukaryota but the census
        # group is a coarser bucket than the sweep group.
        unmapped = files.pop("", set())
        for rel in set(_PROTEOME_DBS.values()):
            files.setdefault(rel, set()).update(unmapped)
        out: dict[str, str] = {}
        for rel, accs in files.items():
            todo = accs - set(out)
            path = self.root / rel
            if not todo or not path.exists():
                continue
            if verbose:
                print(f"  streaming {path.name} for {len(todo)} accessions "
                      f"({path.stat().st_size / 1e9:.1f} GB) ...", flush=True)
            out.update(_grep_fasta(path, todo))
        return out

    def _append_cache(self, found: dict[str, str]) -> None:
        if not found:
            return
        self.cache.parent.mkdir(parents=True, exist_ok=True)
        with open(self.cache, "a") as out:
            for acc, seq in found.items():
                out.write(f">{acc}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")


def _iter_fasta(path: Path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _grep_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    """Streaming extract; stops early once every wanted accession is found."""
    out: dict[str, str] = {}
    keep, acc, chunks = False, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and acc:
                    out[acc] = "".join(chunks)
                    if len(out) == len(wanted):
                        return out
                acc = uniprot_id(line[1:].strip())
                keep, chunks = acc in wanted, []
            elif keep:
                chunks.append(line.strip())
    if keep and acc:
        out[acc] = "".join(chunks)
    return out


# ------------------------------------------------------------------ labels

def sanitize(text: str) -> str:
    return _SAFE.sub("_", (text or "").strip()).strip("_")


def acc_tail(accession: str) -> str:
    """The part of a gene-model id that makes the tip label unique.

    `GCA_001247695.2|copy1|LGRX02003641.1:10549-33606`
        -> `GCA_001247695.2_copy1_10549`

    The assembly and the cell are not enough on their own: the three
    *Myxine glutinosa* loci are all `GCF_964187855.1|ITPR1|…`, because the
    ITPR1 bait won all three, so a two-field label collapsed them into one
    tip name and the duplicate-label guard rejected the whole set. The
    locus start disambiguates without carrying a full coordinate string
    into every figure tick.
    """
    if "|" not in accession:
        return sanitize(accession)
    parts = accession.split("|")
    tail = parts[:2]
    if len(parts) > 2 and ":" in parts[2]:
        span = parts[2].rsplit(":", 1)[1]
        tail.append(span.split("-")[0])
    return sanitize("_".join(tail))


def binomial(species: str) -> str:
    """`Prymnesium parvum (Toxic golden alga)` -> `Prymnesium parvum`.

    The census carries a species string from two sources: UniProt, which
    appends common names, and the genome sweeps, which do not. A dedup key
    on the raw string treats those as two species, which is how *Prymnesium
    parvum* was chosen twice — once as a UniProt record and once as its own
    genome model.
    """
    return " ".join((species or "").split()[:2])


def make_label(group: str, species: str, accession: str) -> str:
    sp = sanitize(species)[:38] or "sp"
    return f"{sanitize(group)}_{sp}_{acc_tail(accession)}"


# -------------------------------------------------------------- diversity

def pick_diverse(rows: list[dict], n: int, key, score) -> list[dict]:
    """Best-scoring rows, spread over `key` in waves.

    `key(row)` may return a scalar or a tuple, and a tuple is treated as a
    **nested** spread: in wave *w* a row is admitted only if every prefix of
    its key is still under quota *w*. With `key = (phylum, genus)` that
    means wave 1 takes one per phylum *and* one per genus, wave 2 allows a
    second of each, and so on.

    The nesting is what a flat key cannot do. Keyed on `(phylum, genus)`
    flatly, the five SAR slots went to three oomycetes of three genera —
    all distinct keys, all one clade — and Ciliophora, the group that
    actually holds the SAR records, got one. Keyed on phylum alone, two of
    three Discoba slots went to two *Naegleria* species.
    """
    ranked = sorted(rows, key=score, reverse=True)
    # A level with only one distinct value among the candidates carries
    # no diversity information — it exists to spread picks across its
    # values, and there is nothing to spread. Left in, it only throttles:
    # every Viridiplantae candidate is Chlorophyta, so the phylum prefix
    # admitted exactly one row per wave and by wave 3 the genus quota was
    # 3, loose enough to let a *second* Chlamydomonas through. The plant
    # grade got two congeners and no Volvox while six other genera waited.
    depth = max((len(key(r)) if isinstance(key(r), tuple) else 1)
                for r in ranked) if ranked else 0
    live = set()
    for i in range(depth):
        vals = {(key(r) if isinstance(key(r), tuple) else (key(r),))[:i + 1]
                for r in ranked}
        if len({v[i] for v in vals}) > 1:
            live.add(i)
    if not live:                 # every level degenerate: rank order only
        live = {depth - 1} if depth else set()
    out: list[dict] = []
    picked: set = set()          # by identity, not by dict equality
    for wave in range(1, len(ranked) + 1):
        # Seeded from what earlier waves already took: a wave that starts
        # from an empty tally lets wave 2 add *two* more of a key that wave
        # 1 already took one of, which is how five SAR slots went to three
        # Triparma species.
        seen: dict = {}
        for done in out:
            kd = key(done)
            kds = kd if isinstance(kd, tuple) else (kd,)
            for i in range(len(kds)):
                if i in live:
                    seen[kds[:i + 1]] = seen.get(kds[:i + 1], 0) + 1
        for r in ranked:
            if len(out) >= n:
                return out
            if id(r) in picked:
                continue
            k = key(r)
            ks = k if isinstance(k, tuple) else (k,)
            prefixes = [ks[:i + 1] for i in range(len(ks)) if i in live]
            if any(seen.get(pre, 0) >= wave for pre in prefixes):
                continue
            for pre in prefixes:
                seen[pre] = seen.get(pre, 0) + 1
            out.append(r)
            picked.add(id(r))
        if len(out) >= n or len(out) == len(ranked):
            break
    return out[:n]


# ------------------------------------------------------------- dashboard

def write_live(task: str, steps: list[tuple[str, bool]]) -> None:
    path = ROOT / "results" / "session_live.json"
    try:
        path.write_text(json.dumps({
            "task": task, "workers": 1,
            "steps": [{"label": l, "done": d} for l, d in steps],
        }, indent=1))
    except OSError:
        pass


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
