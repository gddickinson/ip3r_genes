"""S17 shared helpers — constraint & function.

S17 asks what the sequence says the protein *does*: which parts of the IP3
receptor are evolutionarily intolerant, whether the ligand-binding core and the
pore are under different constraint, and whether that constraint is good enough
to help interpret the family's human variants.

Everything it needs is already on disk:

* `<data_root>/genome_sweep/<acc>/miniprot.gff` — the `##STA` line before each
  gene model is miniprot's translation of that locus. Joined to `summary.json`
  by `mp_id`, that is 244-263 paralog-assigned full-length ITPR proteins per
  paralog across the 309 swept genomes, which is what makes a per-site resource
  possible. `results/msa_v2/aln.fasta` (S6) is 134 sequences from protists to
  vertebrates but **thin per paralog** — 15/11/18 tips — and twelve sequences
  cannot score a column of a 2,758-residue protein.
* `results/loss_dynamics/integrity_loci.tsv` — S15a's ORF screen. A locus with
  an elevated lesion density is an assembly or a decaying gene, and either way
  its residues are not evidence of constraint.
* `results/s0_baseline/review_figures/domain_coords.tsv` — the Pfam architecture
  **measured by InterPro on each human accession separately**, so the sequence
  elements need no transfer at all (this is the one place S17 is better off
  than the PIEZO port it comes from, where every boundary was a mouse Piezo1
  number carried across).
* `results/s0_baseline/review_figures/structure_meta.json` — the pore, measured
  on PDB 6DQN. Filter, gate, IP3 contacts and the TM span are in **Q14573**
  numbering and *are* transferred, with a positive test on arrival (D54).
* `results/selection/codon_*.fasta` + `tree_*.nwk` — S9's codon alignments, the
  input to the per-site selection tests.
* `results/structures/structure_manifest.tsv` — S11's panel, the structures the
  constraint is painted onto.

Conventions this module fixes for every S17 stage:

* **Paralog membership comes from the S7 tree, never from a census label.**
  `leaf_group()` is re-exported from `s13_lib` so S13/S15/S16 and S17 answer
  that question the same way.
* **Conservation is sequence-weighted.** 258 vertebrate orthologues are not 258
  independent observations — Actinopteri alone are a third of them. Henikoff &
  Henikoff (1994) position-based weights are applied before any column
  statistic, or every teleost-specific residue reads as conserved.
* **Coordinates are reported in a named reference's own numbering**, with the
  alignment column kept alongside, so nothing downstream has to guess.

IO, metrics and coordinate mapping only. Analysis lives in the stage modules.
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.utils.data_root import get_data_root  # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "constraint"
MSA_DIR = PROJECT_ROOT / "results" / "msa_v2"
SEL_DIR = PROJECT_ROOT / "results" / "selection"
STRUCT_DIR = PROJECT_ROOT / "results" / "structures"
S0_FIG_DIR = PROJECT_ROOT / "results" / "s0_baseline" / "review_figures"
INTEGRITY = PROJECT_ROOT / "results" / "loss_dynamics" / "integrity_loci.tsv"

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")

#: The three coordinate systems every S17 table can be read in. Unlike the
#: PIEZO project this method is ported from, all three paralogs have a human
#: gene with curated variants, so all three references are human and every
#: table is directly in the numbering ClinVar and UniProt use.
REFERENCES = {
    "ITPR1": ("ITPR1_Homo_sapiens_Human_Q14643", "Q14643", "human ITPR1"),
    "ITPR2": ("ITPR2_Homo_sapiens_Human_Q14571", "Q14571", "human ITPR2"),
    "ITPR3": ("ITPR3_Homo_sapiens_Human_Q14573", "Q14573", "human ITPR3"),
}

#: The reference the *structure* was measured on. S0 measured 6DQN (human
#: ITPR3, IP3-bound, 3.33 A): the four-fold axis, the pore-radius profile, the
#: filter and gate residues and the IP3 contacts. Every one of those numbers is
#: in Q14573 numbering, so ITPR1 and ITPR2 get them by transfer.
STRUCTURE_REF = "ITPR3"

AAS = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AAS)}

#: BLOSUM62 background amino-acid frequencies — the null distribution the
#: Jensen-Shannon divergence is measured against (Capra & Singh 2007).
BLOSUM62_BG = {
    "A": 0.074, "R": 0.052, "N": 0.045, "D": 0.054, "C": 0.025,
    "Q": 0.034, "E": 0.054, "G": 0.074, "H": 0.026, "I": 0.068,
    "L": 0.099, "K": 0.058, "M": 0.025, "F": 0.047, "P": 0.039,
    "S": 0.057, "T": 0.051, "W": 0.013, "Y": 0.032, "V": 0.073,
}


# --------------------------------------------------------------------- FASTA

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path: Path, seqs: dict[str, str], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\t".join(cols or []) + "\n")
        return
    cols = cols or list(rows[0])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def fmt(x, nd: int = 4) -> str:
    """Numbers at `%.6g`, not a fixed decimal count.

    S15a's rule: a p-value of 2.1e-07 written at four decimal places reads back
    as `0.0000`, which hides how strong a claim is rather than how weak.
    """
    if x is None or x == "":
        return ""
    if isinstance(x, float):
        if x != x:
            return ""
        return f"{x:.6g}" if (abs(x) < 1e-4 and x != 0) else f"{round(x, nd)}"
    return str(x)


# ------------------------------------------------------- msa_v2 + tree groups

PHYLO_DIR = PROJECT_ROOT / "results" / "phylogeny"


def msa_groups() -> dict[str, str]:
    """label -> tree-corrected group, for every msa_v2 representative.

    The census label is the starting point and **S7's extended paralog clades
    overrule it** (D30): the seven unlabelled `vertebrate_basal` tips the tree
    nests inside a paralog clade join that paralog, exactly as `s9_sets.py`
    reads them. A tip S7 recorded as `unconstrained` keeps its census name,
    because there the tree offers no alternative — only a refusal to place it —
    and inventing one would be S17 deciding a question S7 declined.

    Read from the committed tables rather than re-derived from the Newick: S9
    and S13 both read them, and a third derivation that disagreed would be a
    silent fork in what this project means by "ITPR2".
    """
    out = {r["label"]: normalise_group(r["group"])
           for r in read_tsv(MSA_DIR / "representatives.tsv")}
    for r in read_tsv(PHYLO_DIR / "paralog_clades.tsv"):
        for lab in r["added"].split(";"):
            if lab and lab in out:
                out[lab] = r["paralog"]
    return out


def normalise_group(g: str) -> str:
    """S23 left one census column in two spellings; counted naively that is
    two groups (`s6_lib.normalise_group`'s reason)."""
    return g[0].upper() + g[1:] if g and g[0].islower() and g not in (
        "protist", "plant", "fungi", "invert_metazoa", "vertebrate_basal") else g


def leaf_group(label: str) -> str:
    return msa_groups().get(label, "")


def integrity_index() -> dict[tuple[str, str, int], dict]:
    """(accession, cell, locus_idx) -> S15a integrity row."""
    return {(r["accession"], r["cell"], int(r["locus_idx"])): r
            for r in read_tsv(INTEGRITY)}


def label_accession(rep_row: dict) -> str | None:
    """Assembly accession behind an msa_v2 representative, if it is one."""
    src = rep_row.get("genome_accession", "") or rep_row.get("accession", "")
    return src.split("|")[0] if src.startswith(("GCF_", "GCA_")) else None


def lesioned_cells() -> set[tuple[str, str]]:
    """(accession, cell) whose ORF S15a called lesion-rich.

    S15a's verdict is one-sided by construction (`s10_orf.py`'s rule: zero
    stops falsifies a pseudogene call, a handful does not establish one), so
    this set is used only to *exclude* sequence from a constraint sample — it
    is never reported as a loss.
    """
    return {(r["accession"], r["cell"]) for r in read_tsv(INTEGRITY)
            if r["verdict"] == "elevated_lesions"}


# ---------------------------------------------------------- sequence weights

def henikoff_weights(seqs: list[str]) -> list[float]:
    """Henikoff & Henikoff (1994) position-based sequence weights.

    Without this a clade that happens to be densely sampled votes as many times
    as it has members, and every teleost-specific residue in a 258-sequence
    vertebrate alignment reads as conserved. Weights are normalised to mean 1 so
    weighted counts stay on the same scale as raw counts.
    """
    n = len(seqs)
    if n == 0:
        return []
    if n == 1:
        return [1.0]
    ncol = len(seqs[0])
    w = [0.0] * n
    for c in range(ncol):
        counts: dict[str, int] = {}
        for s in seqs:
            a = s[c]
            if a in AA_INDEX:
                counts[a] = counts.get(a, 0) + 1
        r = len(counts)
        if r < 2:
            continue                      # invariant/empty columns carry no info
        for i, s in enumerate(seqs):
            a = s[c]
            if a in counts:
                w[i] += 1.0 / (r * counts[a])
    total = sum(w)
    if total <= 0:
        return [1.0] * n
    return [x * n / total for x in w]


# ------------------------------------------------------- conservation metrics

def _kl(p: list[float], q: list[float]) -> float:
    s = 0.0
    for pi, qi in zip(p, q):
        if pi > 0 and qi > 0:
            s += pi * math.log2(pi / qi)
    return s


def column_stats(column: list[str], weights: list[float]) -> dict:
    """Weighted per-column conservation statistics.

    Returns Jensen-Shannon divergence against the BLOSUM62 background (the
    Capra & Singh 2007 metric, in bits, 0-1), normalised entropy, the modal
    residue and its weighted fraction, and the occupancy all three are
    conditional on. Gaps are excluded from the distribution and reported
    separately rather than treated as a 21st character — a gap is missing data
    here, not a shared state.
    """
    wsum = 0.0
    counts = {a: 0.0 for a in AAS}
    for a, w in zip(column, weights):
        if a in AA_INDEX:
            counts[a] += w
            wsum += w
    total_w = sum(weights) or 1.0
    if wsum <= 0:
        return {"jsd": 0.0, "entropy_norm": 0.0, "modal": "-",
                "frac_modal": 0.0, "occupancy": 0.0, "n_seq": 0}
    p = [counts[a] / wsum for a in AAS]
    q = [BLOSUM62_BG[a] for a in AAS]
    r = [(pi + qi) / 2 for pi, qi in zip(p, q)]
    jsd = 0.5 * _kl(p, r) + 0.5 * _kl(q, r)
    h = -sum(pi * math.log(pi) for pi in p if pi > 0)
    modal = max(AAS, key=lambda a: counts[a])
    return {
        "jsd": round(jsd, 4),
        "entropy_norm": round(1.0 - h / math.log(20), 4),
        "modal": modal,
        "frac_modal": round(counts[modal] / wsum, 4),
        "occupancy": round(wsum / total_w, 4),
        "n_seq": sum(1 for a in column if a in AA_INDEX),
    }


def profile(alignment: dict[str, str]) -> list[dict]:
    """Per-column stats for a whole alignment, sequence-weighted."""
    seqs = list(alignment.values())
    if not seqs:
        return []
    w = henikoff_weights(seqs)
    return [column_stats([s[c] for s in seqs], w) for c in range(len(seqs[0]))]


# --------------------------------------------------- coordinates and transfer

def col_to_residue(aln_row: str) -> dict[int, int]:
    """alignment column (0-based) -> 1-based residue number in that sequence."""
    out, r = {}, 0
    for c, ch in enumerate(aln_row):
        if ch not in "-.":
            r += 1
            out[c] = r
    return out


def residue_to_col(aln_row: str) -> dict[int, int]:
    return {v: k for k, v in col_to_residue(aln_row).items()}


def mafft_pair(seq_a: str, seq_b: str) -> tuple[str, str]:
    """Pairwise MAFFT alignment of two sequences (returns both aligned rows).

    `--thread 1` by decision (D24): S3 measured two MAFFT runs on identical
    input giving profiles of 4,933 and 4,908 match states at `--thread -1`, and
    every coordinate transfer in S17 has to be reproducible.
    """
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fa"
        p.write_text(f">a\n{seq_a}\n>b\n{seq_b}\n")
        r = subprocess.run(["mafft", "--auto", "--anysymbol", "--thread", "1",
                            "--quiet", str(p)], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"mafft failed ({r.returncode}): {r.stderr[:300]}")
        out = r.stdout
    rows, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            rows[name] = []
        elif name:
            rows[name].append(line.strip().upper())
    a, b = "".join(rows["a"]), "".join(rows["b"])
    if len(a) != len(b) or not a:
        raise RuntimeError("mafft returned a ragged pair alignment")
    return a, b


def transfer_positions(seq_from: str, seq_to: str) -> dict[int, int]:
    """Residue map seq_from -> seq_to via pairwise alignment (1-based)."""
    if seq_from == seq_to:
        return {i: i for i in range(1, len(seq_from) + 1)}
    a, b = mafft_pair(seq_from, seq_to)
    out, ra, rb = {}, 0, 0
    for ca, cb in zip(a, b):
        if ca != "-":
            ra += 1
        if cb != "-":
            rb += 1
        if ca != "-" and cb != "-":
            out[ra] = rb
    return out


# ------------------------------------------------------------ cached fetching

def cache_dir() -> Path:
    d = get_data_root() / "raw_api" / "s17"
    d.mkdir(parents=True, exist_ok=True)
    return d


def uniprot_fasta(acc: str) -> str:
    path = cache_dir() / f"{acc}.fasta"
    if not path.exists():
        txt = urllib.request.urlopen(
            f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=120
        ).read().decode()
        path.write_text(txt)
    return "".join(path.read_text().split("\n")[1:]).strip()


def uniprot_json(acc: str, fields: str = "") -> dict:
    key = f"{acc}{'_' + fields.replace(',', '-') if fields else ''}.json"
    path = cache_dir() / key
    if not path.exists():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        if fields:
            url += f"?fields={fields}"
        path.write_text(urllib.request.urlopen(url, timeout=120).read().decode())
    return json.loads(path.read_text())


# --------------------------------------------------------- sweep gene models

def sweep_dir() -> Path:
    return get_data_root() / "genome_sweep"


def sta_proteins(acc: str) -> dict[str, str]:
    """mp_id -> miniprot translation, for every gene model in one genome.

    The GFF emits `##PAF` then `##STA <protein>` then the `mRNA` line carrying
    `ID=MP******`, so the translation is the last `##STA` seen when the mRNA
    line arrives. `--trans` output skips frameshifts and ends at the gene's own
    terminator, so the trailing `*` is stripped.
    """
    gff = sweep_dir() / acc / "miniprot.gff"
    if not gff.exists():
        return {}
    out: dict[str, str] = {}
    pending: str | None = None
    with gff.open() as fh:
        for line in fh:
            if line.startswith("##STA"):
                pending = (line.split("\t", 1)[1] if "\t" in line
                           else line[5:]).strip()
            elif not line.startswith("#") and "\tmRNA\t" in line:
                mp = None
                for kv in line.rstrip("\n").split("\t")[8].split(";"):
                    if kv.startswith("ID="):
                        mp = kv[3:]
                        break
                if mp and pending:
                    out[mp] = pending.rstrip("*")
                pending = None
    return out


def iter_sweep_cells():
    """Yield (summary_dict, paralog, cell_dict) for every swept genome x cell."""
    for d in sorted(sweep_dir().glob("GC*")):
        sj = d / "summary.json"
        if not sj.exists():
            continue
        try:
            s = json.loads(sj.read_text())
        except json.JSONDecodeError:
            continue
        for paralog in PARALOGS:
            cell = s.get("cells", {}).get(paralog)
            if cell:
                yield s, paralog, cell
