"""S9 shared helpers — CDS acquisition + validation for the codon alignment.

dN/dS is a statement about codons, so every claim S9 makes rests on a
nucleotide sequence that provably encodes the *exact* protein that went
into S6's alignment. Three routes, one per representative-tip type:

    1. Ensembl protein/transcript ids    -> Ensembl REST /sequence type=cds
       (reuses `src.analysis.selection.fetch_cds`).
    2. UniProt accessions                -> UniProt REST xrefs, tried in
       order: Ensembl translation -> route 1; EMBL protein_id -> ENA
       browser API CDS fasta; RefSeq protein -> Entrez `coded_by`.
    3. S5 miniprot genome models         -> locus realignment, in
       `s9_miniprot_cds.py`. Splicing the sweep GFF's CDS blocks does
       *not* work here and is not offered: the S6 protein comes from
       `miniprot --trans`, whose `##STA` row skips frameshifts, so the
       blocks run out of frame against it partway through.

Every CDS is validated by translating it and requiring the translation to
match the S6 protein (internal stops appear there as `X`; the matching
codons are masked to NNN so pal2nal and codeml see ambiguity, not stops).
A tip whose CDS cannot be validated is *recorded as failed and dropped* —
never repaired by trimming to a length that happens to fit.

Raw fetches are cached under `<data_root>/raw_api/s9_cds/`.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.analysis.selection import GENETIC_CODE, fetch_cds  # noqa: E402
from src.utils.data_root import get_data_root  # noqa: E402

from s3_hmm_lib import read_fasta  # noqa: E402

MSA_DIR = ROOT / "results" / "msa_v2"
REPS_TSV = MSA_DIR / "representatives.tsv"
REPS_FASTA = MSA_DIR / "representatives.fasta"
ALN_FASTA = MSA_DIR / "aln.fasta"
PHYLO_DIR = ROOT / "results" / "phylogeny"
ROOTED_NWK = PHYLO_DIR / "rooted.nwk"
BAITS_FAA = ROOT / "results" / "s5_baits" / "baits.faa"
OUT_DIR = ROOT / "results" / "selection"
CACHE_DIR = get_data_root() / "raw_api" / "s9_cds"

# The codon alignment is a *vertebrate family* alignment. The RyR outgroup
# is deliberately not in it: it roots the S7 protein tree, but at that
# distance dS is saturated and a synonymous-site model estimated across it
# would be measuring alignment error. The paralog stems are unambiguous
# without it — in the unrooted tree of ITPR1/ITPR2/ITPR3 plus the
# unlabelled vertebrate grade, the branch subtending each paralog clade is
# well defined (D30: a paralog label means nothing outside the vertebrates,
# and nothing outside the vertebrates is in this set).
VERTEBRATE_GROUPS = ("ITPR1", "ITPR2", "ITPR3", "vertebrate_basal")
PARALOGS = ("ITPR1", "ITPR2", "ITPR3")

_PAUSE_S = 0.34

#: The selection toolchain — codeml, yn00, pal2nal.pl, hyphy — is
#: env-resident, not on PATH (D18). S1's manifest never probed it, so the
#: first run of this task found `pal2nal.pl` missing on a machine where it
#: has been installed all along.
ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")


def tool_bin(name: str) -> str:
    """Absolute path to an external tool, PATH first then the piezo1 env."""
    found = shutil.which(name)
    if found:
        return found
    env = ENV_BIN / name
    if env.exists():
        return str(env)
    raise SystemExit(f"{name} not found on PATH or in the piezo1 env "
                     "(see results/toolchain_manifest.txt)")


# ---- tiny IO ---------------------------------------------------------------

def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def translate(cds: str) -> str:
    """Full translation — stops emitted as `*`, never truncated at one.

    `src.analysis.selection.translate` stops at the first terminator, which
    is right for a coding sequence and wrong here: the miniprot models
    carry internal stops that have to be *seen* so their codons can be
    masked.
    """
    return "".join(GENETIC_CODE.get(cds[i:i + 3].upper().replace("U", "T"), "X")
                   for i in range(0, len(cds) - len(cds) % 3, 3))


def _http_get(url: str, cache_name: str, retries: int = 3) -> str | None:
    """GET with an on-disk cache; returns text or None."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / cache_name
    if cache.exists() and cache.stat().st_size > 0:
        return cache.read_text()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ip3r-s9"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                text = resp.read().decode()
            cache.write_text(text)
            time.sleep(_PAUSE_S)
            return text
        except Exception:  # noqa: BLE001
            time.sleep(2.0 * (attempt + 1))
    return None


# ---- route 1: Ensembl ------------------------------------------------------

def cds_from_ensembl(ensembl_id: str, retries: int = 4) -> tuple[str | None, str]:
    """Ensembl REST CDS, cached. Retries with backoff — the service returns
    500/503 in bursts and a single failure would silently drop a tip."""
    cache = CACHE_DIR / f"ensembl_{ensembl_id}.json"
    if cache.exists():
        data = json.loads(cache.read_text())
        if data.get("cds"):
            return data["cds"], data.get("note", "cached")
    last = ""
    for attempt in range(retries):
        try:
            cds, species, transcript = fetch_cds(ensembl_id)
            note = f"ensembl:{transcript}:{species}"
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps({"cds": cds, "note": note}))
            return cds, note
        except Exception as exc:  # noqa: BLE001
            last = str(exc)[:120]
            time.sleep(3.0 * (attempt + 1))
    return None, f"ensembl_failed:{last}"


# ---- route 2: UniProt xrefs ------------------------------------------------

def uniprot_acc_for_ensembl(ensembl_id: str) -> str | None:
    """UniProt accession cross-referencing an Ensembl id (reverse lookup)."""
    url = (f"https://rest.uniprot.org/uniprotkb/search?query=xref:{ensembl_id}"
           "&format=json&fields=accession")
    text = _http_get(url, f"uniprot_xref_{ensembl_id}.json")
    if not text:
        return None
    try:
        results = json.loads(text).get("results", [])
    except Exception:  # noqa: BLE001
        return None
    return results[0]["primaryAccession"] if results else None


def ensembl_candidates(ensembl_id: str):
    """Every CDS the Ensembl route can offer for an Ensembl id."""
    cds, note = cds_from_ensembl(ensembl_id)
    if cds:
        yield cds, note
    acc = uniprot_acc_for_ensembl(ensembl_id)
    if acc:
        for cds2, note2 in uniprot_candidates(acc, skip_ensembl=True):
            yield cds2, f"ens_via_uniprot:{acc}->{note2}"


def _uniprot_xrefs(acc: str) -> dict:
    url = (f"https://rest.uniprot.org/uniprotkb/{acc}.json"
           "?fields=xref_ensembl,xref_embl,xref_refseq,sequence")
    text = _http_get(url, f"uniprot_{acc}.json")
    if not text:
        return {}
    try:
        return json.loads(text)
    except Exception:  # noqa: BLE001
        return {}


def uniprot_candidates(acc: str, skip_ensembl: bool = False):
    """Every CDS a UniProt accession cross-references, in preference order.

    **All** of them, not the first that downloads. A UniProt entry lists
    every Ensembl transcript of its gene, and the entry's own sequence is
    one particular isoform: human ITPR1 lists five transcripts and human
    ITPR2 two, and in both the first is not the canonical one S6 aligned.
    Returning the first CDS that *fetches* silently substitutes a different
    isoform for the protein the tree was built on — which is a wrong codon
    alignment, not a missing one. The caller validates by translation and
    keeps the first that matches.
    """
    entry = _uniprot_xrefs(acc)
    xrefs = entry.get("uniProtKBCrossReferences", [])
    if skip_ensembl:
        xrefs = [x for x in xrefs if x.get("database") != "Ensembl"]
    # (a) Ensembl xrefs -> route 1. The translation id lives in a property;
    # the xref's own `id` is the *transcript* — many non-model species have
    # only the latter, so both are tried.
    for x in xrefs:
        if x.get("database") != "Ensembl":
            continue
        ids = [p["value"].split(".")[0] for p in x.get("properties", [])
               if p.get("key") == "ProteinId" and p.get("value", "-") != "-"]
        if x.get("id"):
            ids.append(x["id"].split(".")[0])
        for eid in ids:
            cds, note = cds_from_ensembl(eid)
            if cds:
                yield cds, f"uniprot->{note}"
    # (b) EMBL protein_id -> ENA CDS fasta
    for x in xrefs:
        if x.get("database") != "EMBL":
            continue
        pid = next((p["value"] for p in x.get("properties", [])
                    if p.get("key") == "ProteinId"), None)
        if not pid or pid == "-":
            continue
        text = _http_get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pid}",
                         f"ena_{pid}.fasta")
        if text and text.startswith(">"):
            seq = "".join(text.split("\n")[1:]).strip().upper()
            if seq:
                yield seq, f"ena:{pid}"
    # (c) RefSeq protein -> Entrez coded_by
    for x in xrefs:
        if x.get("database") != "RefSeq":
            continue
        cds, note = cds_from_refseq(x.get("id", ""))
        if cds:
            yield cds, note


def cds_from_refseq(protein_id: str) -> tuple[str | None, str]:
    """GenPept `coded_by` -> nuccore range fetch (Entrez, no API key)."""
    if not protein_id:
        return None, "refseq_no_id"
    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    gp = _http_get(f"{eutils}?db=protein&id={protein_id}&rettype=gp&retmode=text",
                   f"gp_{protein_id}.txt")
    if not gp:
        return None, "refseq_gp_failed"
    m = re.search(r'/coded_by="(complement\()?(?:join\()?([A-Z_0-9.]+):[<>]?(\d+)'
                  r'\.\.[<>]?(\d+)', gp)
    if not m:
        return None, "refseq_no_coded_by"
    comp, nuc, start, end = m.group(1), m.group(2), m.group(3), m.group(4)
    fa = _http_get(f"{eutils}?db=nuccore&id={nuc}&rettype=fasta&retmode=text"
                   f"&seq_start={start}&seq_stop={end}"
                   + ("&strand=2" if comp else ""),
                   f"nuc_{nuc}_{start}_{end}{'_c' if comp else ''}.fasta")
    if not fa or not fa.startswith(">"):
        return None, "refseq_nuc_failed"
    seq = "".join(fa.split("\n")[1:]).strip().upper()
    return (seq or None), f"refseq:{nuc}:{start}-{end}"


# ---- route 3: miniprot genome models ---------------------------------------

def find_genome_fna(acc: str) -> Path | None:
    gdir = get_data_root() / "genomes" / acc
    if not gdir.is_dir():
        return None
    hits = sorted(gdir.rglob("*.fna")) + sorted(gdir.rglob("*.fa"))
    hits = [h for h in hits if "cds_from" not in h.name and "_rna" not in h.name]
    return hits[0] if hits else None


def cds_from_miniprot_model(assembly: str, expected_prot: str,
                            sweep_root: str = "genome_sweep",
                            threads: int = 8,
                            locus: tuple[str, int, int, str] | None = None,
                            mp_id: str | None = None,
                            ) -> tuple[str | None, str]:
    """Exact CDS for an S5 gene model — see `s9_miniprot_cds` for why the
    locus has to be realigned rather than spliced out of the sweep GFF."""
    from s9_miniprot_cds import cds_for_model  # noqa: PLC0415  (heavy import)

    fna = find_genome_fna(assembly)
    if fna is None:
        return None, f"no_genome_fasta:{assembly}"
    baits = read_fasta(BAITS_FAA)
    return cds_for_model(assembly, get_data_root() / sweep_root, fna,
                         baits, expected_prot, threads, locus)


# ---- validation ------------------------------------------------------------

def validate_and_mask(cds: str, expected_prot: str) -> tuple[str | None, dict]:
    """Check `translate(cds)` against the S6 protein; mask internal-stop
    codons to NNN wherever the protein shows X (or the codon is a stop).

    Returns (masked_cds or None, stats). Acceptance: same length after the
    terminator strip and >= 99 % identity at comparable positions.
    """
    cds = cds.upper().replace("U", "T")
    if len(cds) % 3:
        cds = cds[:-(len(cds) % 3)]
    prot = translate(cds)
    # Strip the terminator codon only when that is what reconciles the two
    # lengths. A gene model whose *final* residue is an internal stop (the
    # miniprot models carry several) already matches, and stripping there
    # would fail an otherwise perfect reconstruction on an off-by-one.
    if prot.endswith("*") and len(prot) == len(expected_prot) + 1:
        prot = prot[:-1]
        cds = cds[:-3]
    stats = {"cds_len": len(cds), "trans_len": len(prot),
             "prot_len": len(expected_prot), "n_mismatch": -1,
             "n_stop_masked": 0, "n_masked": 0, "note": ""}
    if len(prot) != len(expected_prot):
        stats["note"] = "length_mismatch"
        return None, stats
    # Every disagreement is masked, not only the expected ones. A CDS that
    # translates to something other than the protein its alignment column
    # holds is not a CDS of that protein at that site, whatever the cause —
    # an internal stop, an `X`, or a terminal exon the locus rerun placed
    # differently from the sweep. Masking makes the pair consistent and
    # makes codeml treat the site as missing; leaving the codon in would
    # hand a synonymous/non-synonymous count a residue that is not there.
    mism = 0
    masked = list(cds)
    n_stop = n_masked = 0
    for i, (a, b) in enumerate(zip(prot, expected_prot)):
        if a == b:
            continue
        masked[i * 3:i * 3 + 3] = "NNN"
        n_masked += 1
        if a == "*":
            n_stop += 1
        elif not (b == "X" or a == "X"):
            mism += 1
    stats["n_mismatch"] = mism
    stats["n_stop_masked"] = n_stop
    stats["n_masked"] = n_masked
    comparable = sum(1 for a, b in zip(prot, expected_prot)
                     if b != "X" and a not in "*X")
    stats["n_comparable"] = comparable
    if comparable and mism / comparable > 0.01:
        stats["note"] = "identity_below_99pct"
        return None, stats
    return "".join(masked), stats


# ---- representative table --------------------------------------------------

def load_reps() -> list[dict]:
    with open(REPS_TSV) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_vertebrate_tips() -> list[dict]:
    """Rows of representatives.tsv restricted to the vertebrate family groups."""
    return [r for r in load_reps() if r["group"] in VERTEBRATE_GROUPS]


# `GCF_052327205.1|ITPR1|NC_140997.1:147905519-148503708+` — the S5/S23
# gene-model accession. It carries no trailing pipe in this project, unlike
# the PIEZO port's, and a regex requiring one silently routes every genome
# model to the UniProt path, where it 400s.
LOCUS_RE = re.compile(r"^(GC[AF]_[0-9.]+)\|([^|]+)\|([^:]+):(\d+)-(\d+)([+-])$")
MPID_RE = re.compile(r"^(GC[AF]_[0-9.]+)\|(MP\d+)$")


def route_for(row: dict) -> tuple[str, dict]:
    """(route_name, kwargs) for one representative row."""
    acc = row["accession"]
    m = MPID_RE.match(acc)
    if m:
        return "miniprot", {"assembly": m.group(1), "mp_id": m.group(2)}
    m = LOCUS_RE.match(acc)
    if m:
        return "miniprot", {"assembly": m.group(1),
                            "locus": (m.group(3), int(m.group(4)),
                                      int(m.group(5)), m.group(6))}
    if acc.startswith("ENS"):
        return "ensembl", {"ensembl_id": acc}
    return "uniprot", {"acc": acc}
