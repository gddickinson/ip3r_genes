"""S12 shared library — paths, cached HTTP, SRA metadata, tissue vocabulary.

S12 asks whether the genes S10 proved the annotation loses are *transcribed*.
The evidence is read-level: public RNA-seq runs are streamed from SRA and
aligned against a small per-species reference holding every family locus the
sweep recovered in that genome, a housekeeping anchor and a reversed decoy
per locus.

Everything network-facing is cached under ``<data_root>/raw_api/s12/`` so a
re-run is offline and the run selection is reproducible.

Two environments (D18).  The metadata and reference steps run under the
`piezo1` env; the streaming aligner needs `sra-tools` and `hisat2`, which
live in the `s12` env.  Nothing here may import the project's `src` package:
``src.utils.__init__`` pulls in the GUI/database stack, which the `s12` env
does not have, so the one stdlib-only module needed is loaded from its file.

Modules in the S12 chain:
    s12_panel.py       the scope rules (P1-P4), derived from S10's ranking
    s12_lib.py         (this file) paths, SRA E-utilities, tissue vocabulary
    s12_refs.py        per-species reference FASTA + junction table
    s12_runs.py        run selection -> runs_selected.tsv
    s12_quantify.py    stream + align + count -> per-run, per-junction counts
    s12_tables.py      the committed tables + stats json
    s12_atlas.py       the independent cross-check
    s12_report*.py     report.md, rendered from the committed tables only
    s12_figures.py     four figures
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

_spec = importlib.util.spec_from_file_location(
    "_s12_data_root", PROJECT / "src" / "utils" / "data_root.py")
_dr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_dr)
get_data_root = _dr.get_data_root

OUT_DIR = PROJECT / "results" / "expression"
DATA = get_data_root() / "expression"
CACHE = get_data_root() / "raw_api" / "s12"
GENOMES = get_data_root() / "genomes"
SWEEP = get_data_root() / "genome_sweep"

EMAIL = "george.dickinson@gmail.com"
TOOL = "ip3r_genes_s12"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# The aligner env (D18).  sra-tools and hisat2 are not on the bare PATH and
# not in `piezo1`; resolving them by name alone would make S12 fail in a way
# that reads as "no reads found" rather than "no aligner".
S12_ENV_BIN = Path("/opt/anaconda3/envs/s12/bin")
PIEZO1_ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")
ENV_BINS = (S12_ENV_BIN, PIEZO1_ENV_BIN)


def tool_bin(name: str) -> str:
    """Absolute path to an external tool: PATH first, then the envs (D18).

    Both envs are searched. S12 straddles them — sra-tools and hisat2 live
    in `s12`, while blast and the datasets CLI live in `piezo1` — and a
    lookup that knew about only one would report a present tool as
    unavailable, which in the stats file reads as a stage that could not
    run rather than a path that was not searched.
    """
    found = shutil.which(name)
    if found:
        return found
    for env in ENV_BINS:
        cand = env / name
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(
        f"{name} not found on PATH or in "
        + ", ".join(str(e) for e in ENV_BINS))


def tool_version(name: str) -> str:
    """Version string for the manifest, or 'unavailable'."""
    try:
        exe = tool_bin(name)
    except FileNotFoundError:
        return "unavailable"
    for flag in ("--version", "-version"):
        try:
            p = subprocess.run([exe, flag], capture_output=True, text=True,
                               timeout=60)
        except (OSError, subprocess.SubprocessError):
            continue
        text = (p.stdout or "") + (p.stderr or "")
        for line in text.split("\n"):
            if re.search(r"\d+\.\d+", line):
                return line.strip()[:80]
    return "unknown"


# ---------------------------------------------------------------------------
# tissue vocabulary
#
# Used twice: to *query* SRA and to *verify* that the returned run really
# carries that tissue in its own sample attributes.  A run whose attributes
# do not corroborate the query term is dropped, so a hit on a study title
# alone cannot mislabel a tissue.
# ---------------------------------------------------------------------------

TISSUES: dict[str, tuple[str, ...]] = {
    "brain": ("brain", "telencephalon", "forebrain", "cerebrum", "cortex",
              "hypothalamus", "pituitary"),
    "heart": ("heart", "cardiac", "cardiomyocyte", "ventricle"),
    "muscle": ("muscle", "skeletal muscle", "myotome", "myotomal"),
    "skin": ("skin", "epidermis", "dermis", "scale"),
    "liver": ("liver", "hepatic", "hepatocyte", "hepatopancreas"),
    "kidney": ("kidney", "renal", "pronephros", "mesonephros", "head kidney"),
    "intestine": ("intestine", "gut", "intestinal", "colon", "stomach",
                  "pyloric caeca"),
    "testis": ("testis", "testes", "testicular", "sperm", "milt"),
    "ovary": ("ovary", "ovarian", "oocyte", "egg", "roe"),
    "gill": ("gill", "gills", "branchial"),
    "spleen": ("spleen", "splenic"),
    "eye": ("eye", "retina", "retinal"),
    "blood": ("blood", "erythrocyte", "leukocyte", "peripheral blood"),
    "embryo": ("embryo", "embryonic", "blastula", "gastrula", "larva",
               "larval", "fry", "juvenile"),
    "fin": ("fin", "caudal fin", "pectoral fin"),
    "swim_bladder": ("swim bladder", "swimbladder", "gas bladder"),
}

# attribute keys that may name the biological source
TISSUE_KEYS = ("tissue", "organ", "source_name", "sourcename", "organism_part",
               "body_site", "isolation_source", "cell_type", "sample_type",
               "dev_stage", "developmental_stage", "tissue_type",
               "tissue-type", "sample_name")


def normalise_tissue(text: str) -> str | None:
    """Map free-text sample metadata onto the controlled tissue vocabulary.

    Longest keyword first: 'head kidney' must win over 'kidney' and
    'swim bladder' over 'blood' when both are present in one string, or the
    tissue a run is filed under depends on dict insertion order.
    """
    low = " " + re.sub(r"[^a-z0-9 ]+", " ", text.lower()) + " "
    best: tuple[int, str] | None = None
    for tissue, words in TISSUES.items():
        for w in words:
            if f" {w} " in low and (best is None or len(w) > best[0]):
                best = (len(w), tissue)
    return best[1] if best else None


# ---------------------------------------------------------------------------
# cached HTTP
# ---------------------------------------------------------------------------

def _cache_path(kind: str, key: str) -> Path:
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    d = CACHE / kind
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{h}.gz"


def cached_get(url: str, kind: str, key: str | None = None,
               pause: float = 0.35, retries: int = 4,
               timeout: int = 30) -> str:
    """GET with a permanent on-disk cache. NCBI asks for <= 3 req/s."""
    path = _cache_path(kind, key or url)
    if path.exists():
        with gzip.open(path, "rt") as fh:
            return fh.read()
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": f"{TOOL} ({EMAIL})"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                text = r.read().decode("utf-8", "replace")
            with gzip.open(path, "wt") as fh:
                fh.write(text)
            time.sleep(pause)
            return text
        except urllib.error.HTTPError as exc:
            # 4xx is a permanent answer; retrying it just burns minutes.
            last = exc
            if 400 <= exc.code < 500 and exc.code != 429:
                raise RuntimeError(f"GET {exc.code}: {url}") from exc
            time.sleep(2 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries}: {url} ({last})")


def eutils(endpoint: str, **params) -> str:
    params.setdefault("email", EMAIL)
    params.setdefault("tool", TOOL)
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    return cached_get(url, kind=endpoint.split(".")[0], key=url)


# ---------------------------------------------------------------------------
# SRA search + metadata
# ---------------------------------------------------------------------------

@dataclass
class SraRun:
    run: str
    experiment: str
    study: str
    study_title: str
    sample: str
    organism: str
    strategy: str
    selection: str
    source: str
    layout: str
    platform: str
    model: str
    spots: int
    bases: int
    attributes: dict = field(default_factory=dict)
    tissue: str | None = None
    tissue_source: str = ""

    @property
    def read_len(self) -> int:
        return int(self.bases / self.spots) if self.spots else 0


def esearch_sra(term: str, retmax: int = 200) -> list[str]:
    txt = eutils("esearch.fcgi", db="sra", term=term, retmax=retmax,
                 retmode="json")
    try:
        return json.loads(txt)["esearchresult"].get("idlist", [])
    except (json.JSONDecodeError, KeyError):
        return []


def esearch_count(term: str) -> int:
    txt = eutils("esearch.fcgi", db="sra", term=term, retmax=0, retmode="json")
    try:
        return int(json.loads(txt)["esearchresult"]["count"])
    except (json.JSONDecodeError, KeyError, ValueError):
        return 0


def _text(node, path: str, default: str = "") -> str:
    if node is None:
        return default
    el = node.find(path)
    return (el.text or default) if el is not None and el.text else default


def efetch_sra(uids: list[str], batch: int = 20) -> list[SraRun]:
    """Full SRA XML for a set of UIDs -> one SraRun per RUN element."""
    out: list[SraRun] = []
    for i in range(0, len(uids), batch):
        chunk = uids[i:i + batch]
        xml = eutils("efetch.fcgi", db="sra", id=",".join(chunk),
                     retmode="xml")
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            continue
        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            out.extend(_runs_from_package(pkg))
    return out


def _runs_from_package(pkg) -> list[SraRun]:
    exp = pkg.find("EXPERIMENT")
    if exp is None:
        return []
    lib = exp.find(".//LIBRARY_DESCRIPTOR")
    layout = "PAIRED" if exp.find(".//LIBRARY_LAYOUT/PAIRED") is not None \
        else "SINGLE"
    plat_el = exp.find("PLATFORM")
    platform = model = ""
    if plat_el is not None and len(plat_el):
        platform = plat_el[0].tag
        model = _text(plat_el[0], "INSTRUMENT_MODEL")

    attrs: dict = {}
    for a in pkg.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        tag = _text(a, "TAG").strip().lower().replace(" ", "_")
        val = _text(a, "VALUE").strip()
        if tag and val:
            attrs[tag] = val
    sample = pkg.find(".//SAMPLE")
    organism = _text(sample, ".//SCIENTIFIC_NAME") if sample is not None else ""
    sample_acc = sample.get("accession", "") if sample is not None else ""
    if sample is not None:
        for tag, key in (("TITLE", "_sample_title"),
                         ("DESCRIPTION", "_sample_description")):
            val = _text(sample, tag)
            if val:
                attrs.setdefault(key, val)
    exp_title = _text(exp, "TITLE")
    if exp_title:
        attrs.setdefault("_experiment_title", exp_title)

    study_title = _text(pkg, ".//STUDY/DESCRIPTOR/STUDY_TITLE")
    st = pkg.find(".//STUDY")
    study = st.get("accession", "") if st is not None else ""

    runs = []
    for run in pkg.findall(".//RUN_SET/RUN"):
        acc = run.get("accession", "")
        if not acc:
            continue
        runs.append(SraRun(
            run=acc, experiment=exp.get("accession", ""), study=study,
            study_title=study_title, sample=sample_acc, organism=organism,
            strategy=_text(lib, "LIBRARY_STRATEGY"),
            selection=_text(lib, "LIBRARY_SELECTION"),
            source=_text(lib, "LIBRARY_SOURCE"),
            layout=layout, platform=platform, model=model,
            spots=int(run.get("total_spots") or 0),
            bases=int(run.get("total_bases") or 0),
            attributes=dict(attrs)))
    return runs


# Which attribute the tissue was read from is recorded on every run, because
# "brain" taken from `tissue=brain` and "brain" guessed from a study title
# are not the same evidence and the table has to say which it is.
TITLE_KEYS = ("_sample_title", "_sample_description", "_experiment_title")


def assign_tissue(run: SraRun, want: str | None = None,
                  allow_title: bool = True) -> bool:
    """Corroborate a tissue from the run's own metadata.

    Attribute keys are tried before free-text titles.  If ``want`` is given,
    only that tissue counts: a run returned by a 'brain' query whose
    attributes say 'liver' is a mislabel and is rejected rather than
    silently refiled.
    """
    for key in TISSUE_KEYS:
        val = run.attributes.get(key)
        if not val:
            continue
        t = normalise_tissue(val)
        if t and (want is None or t == want):
            run.tissue, run.tissue_source = t, f"{key}={val[:60]}"
            return True
    if not allow_title:
        return False
    for key in TITLE_KEYS:
        val = run.attributes.get(key, "")
        if not val:
            continue
        t = normalise_tissue(val)
        if t and (want is None or t == want):
            run.tissue, run.tissue_source = t, f"{key.lstrip('_')}={val[:60]}"
            return True
    return False


# ---------------------------------------------------------------------------
# FASTA + small helpers
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line.strip():
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> Path:
    """Write a TSV **atomically** — temp file, then rename.

    S12's stages overlap in practice: the quantifier reads
    `reference_table.tsv` once per run while a later edit may be rewriting
    it, and a reader that catches a half-written file gets a short table
    with no error, which downstream reads as references that do not exist.
    A rename is atomic on the same filesystem, so a reader sees either the
    old file or the new one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    with open(tmp, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    tmp.replace(path)
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def live(stage: str, done: int, total: int, note: str = "") -> None:
    """Feed the dashboard's live panel (best-effort, never fatal)."""
    payload = {"task": "S12", "stage": stage, "done": done, "total": total,
               "note": note, "ts": time.strftime("%H:%M:%S")}
    try:
        (PROJECT / "results" / "session_live.json").write_text(
            json.dumps(payload))
    except OSError:
        pass


CODONS = {}
_B = "TCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for _i, _a in enumerate(_AA):
    CODONS[_B[_i // 16] + _B[(_i // 4) % 4] + _B[_i % 4]] = _a


def translate(nt: str) -> str:
    return "".join(CODONS.get(nt[i:i + 3].upper(), "X")
                   for i in range(0, len(nt) - len(nt) % 3, 3))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
