"""s19_lib_base.py — S19's shared IO, minus the proteome universe.

Imported through `s19_lib`, which re-exports this module and
`s19_universe` together; nothing outside imports this name directly.

Original docstring follows.


S19 asks what the **search** was worth, not what the genes are. Three things
this project has that make that answerable more sharply than usual:

* **A family with no losses.** S15b's Dollo count is zero: all 923 assignable
  genome x paralog cells hold a present gene, and the four that are not
  assignable are cyclostome cells the panel has no labelled bait for. So
  every ledger cell that is not `found_*` is a **false negative of the
  method**, and the sweep's sensitivity can be read off directly instead of
  being confounded with real absence. `control_cells()` is that definition.
* **A sister family swept in the same pass.** The RyR cell fired in all 309
  genomes, is present in every vertebrate, and was never used to call an
  ITPR. It is a fourth, independent false-negative series measured by the
  same instrument in the same assemblies.
* **The retained alignments.** miniprot aligns each bait independently, so
  every genome's `miniprot.gff` holds a separable record per bait: dropping
  baits and re-clustering reproduces exactly what the sweep would have
  reported had those baits never been in the panel (`s19_panel`).

The one derived asset built here is the **swept accession universe**
(`build_universe`): without it, "the InterPro enumeration holds a record the
profile HMM did not return" cannot be told apart from "that record was never
in the database the profile HMM searched". One `grep '^>'` pass per proteome
FASTA, cached under `<data_root>/raw_api/s19/`.

Table IO, the sweep layout and the genome index are taken from `s16_lib`
rather than re-implemented, so the two tasks cannot drift.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

import s16_lib as L                                             # noqa: E402

read_tsv = L.read_tsv
write_tsv = L.write_tsv
sweep_dir = L.sweep_dir
data_root = L.data_root
fnum = L.fnum
fmt = L.fmt

RESULTS = PROJECT_ROOT / "results"
OUT_DIR = RESULTS / "methods"
HMM_DIR = RESULTS / "hmm_sweep"
S20_DIR = RESULTS / "s20_sweep"
S23_DIR = RESULTS / "s23_scope"
LEDGER = RESULTS / "genome_ledger" / "genome_ledger.tsv"
MANIFEST = RESULTS / "genome_manifest.tsv"
BAIT_MANIFEST = RESULTS / "s5_baits" / "bait_manifest.tsv"
MATRIX = RESULTS / "loss_dynamics" / "character_matrix.tsv"

PARALOGS = L.PARALOGS                     # ("ITPR1", "ITPR2", "ITPR3")
CONTROL_CELL = L.CONTROL_CELL             # "RYR"
CELLS = PARALOGS + (CONTROL_CELL,)

FOUND_STATUSES = ("found_annotated", "found_unannotated", "found_no_annotation")

#: S15a states that mean "the gene is there" — every one of them is reached by
#: a positive test (a placed locus, or a reassembly clearing the calibrated
#: bar). `paralog_unassignable` is deliberately **not** here: it means the
#: genome carries spare ITPR loci the panel cannot label, which is an
#: undecidable cell rather than a demonstrated present gene, and counting it
#: as either a present gene or an absence would be a decision S15a declined.
PRESENT_STATES = ("present_single_locus", "present_truncated",
                  "present_fragmented", "present_partial")
UNDECIDABLE_STATES = ("paralog_unassignable",)

#: D4's contiguity bar, read back from S5's own calibration rather than
#: retyped (the discipline D13 applies to reports, applied to a constant).
def contiguity_bar() -> int:
    import s5_calibration as C
    return int(C.itpr_span_stats()["median"])


def out_dir() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR


def cache_dir() -> Path:
    d = data_root() / "raw_api" / "s19"
    d.mkdir(parents=True, exist_ok=True)
    return d


def log(msg: str) -> None:
    print(f"[s19] {msg}", flush=True)


def write_json(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, default=str))
    return path


def read_json(path: Path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def is_found(status: str) -> bool:
    return status in FOUND_STATUSES


# ------------------------------------------------------------ census / methods

CENSUS = {
    "v2": RESULTS / "census_v2" / "census_v2.tsv",
    "v3": RESULTS / "census_v3" / "census_v3.tsv",
    "v4": RESULTS / "census_v4" / "census_v4.tsv",
    "v5": RESULTS / "census_v5" / "census_v5.tsv",
    "v6": RESULTS / "census_v6" / "census_v6.tsv",
}

#: The four app bundles S2's delta measured census v1 from. Kept as the
#: committed delta summary rather than re-read, because the bundles are the
#: app's own output format and S2 already resolved every accession in them.
V1_DELTA = RESULTS / "census_v2" / "delta_summary.json"


def census_rows(version: str) -> list[dict]:
    return read_tsv(CENSUS[version])


def acc_key(acc: str) -> str:
    """UniProt accession without an isoform or version suffix.

    Genome-derived rows carry a composite model id
    (`GCA_...|ITPR1|contig:start-end`), and stripping those at the first `-`
    or `.` would collapse every model of one assembly onto its accession
    prefix — 1,058 gene models would read as a few hundred "records".
    Composite ids are returned unchanged.
    """
    acc = (acc or "").strip()
    if any(ch in acc for ch in "|:/"):
        return acc
    return acc.split("-")[0].split(".")[0]


#: The search channels S19 scores, in the order they were run.
METHOD_ORDER = ("app_search", "pfam_enumeration", "profile_hmm_vert",
                "jackhmmer_vert", "genome_miniprot_vert",
                "profile_hmm_nonvert", "jackhmmer_nonvert",
                "genome_miniprot_nonvert")

METHOD_LABEL = {
    "app_search": "targeted database search (NCBI/Ensembl/UniProt/Compara)",
    "pfam_enumeration": "InterPro/Pfam exhaustive enumeration",
    "profile_hmm_vert": "profile HMM (hmmsearch, 763 vertebrate proteomes)",
    "jackhmmer_vert": "jackhmmer iterative (vertebrate proteomes)",
    "genome_miniprot_vert": "genome sweep (miniprot, 309 vertebrate assemblies)",
    "profile_hmm_nonvert": "profile HMM (hmmsearch, 6,928 other proteomes)",
    "jackhmmer_nonvert": "jackhmmer iterative (non-vertebrate proteomes)",
    "genome_miniprot_nonvert": "genome sweep (miniprot, 194 other assemblies)",
}

METHOD_UNIVERSE = {
    "app_search": "public_databases",
    "pfam_enumeration": "uniprotkb",
    "profile_hmm_vert": "vertebrata",
    "jackhmmer_vert": "vertebrata",
    "genome_miniprot_vert": "genomes_309",
    "profile_hmm_nonvert": "nonvert_groups",
    "jackhmmer_nonvert": "nonvert_groups",
    "genome_miniprot_nonvert": "genomes_194",
}

NONVERT_GROUPS = ("fungi", "viridiplantae", "metazoa_nonvert",
                  "protista_other", "archaea", "bacteria_genus")
ALL_GROUPS = ("vertebrata",) + NONVERT_GROUPS


def _tsv_accessions(path: Path, column: str = "accession") -> set[str]:
    if not path.exists():
        return set()
    return {a for a in (acc_key(r.get(column, "")) for r in read_tsv(path)) if a}


def method_accessions() -> dict[str, set[str]]:
    """Accession set per search channel.

    Read from each channel's **own** committed output, never from the census's
    `source` column: the census records which channel a record was *first*
    seen by, and the question here is which channels could see it at all.
    """
    sets: dict[str, set[str]] = {}
    sets["pfam_enumeration"] = {acc_key(r["accession"])
                                for r in census_rows("v2")}
    sets["profile_hmm_vert"] = _tsv_accessions(
        HMM_DIR / "hmmsearch_assignments.tsv")

    nv: set[str] = set()
    for g in NONVERT_GROUPS:
        nv |= _tsv_accessions(S20_DIR / f"assignments_{g}.tsv")
    sets["profile_hmm_nonvert"] = nv

    sets["genome_miniprot_vert"] = {r["model_id"] for r in
                                    read_tsv(RESULTS / "census_v4" /
                                             "genome_models.tsv")}
    sets["genome_miniprot_nonvert"] = {r["model_id"] for r in
                                       read_tsv(RESULTS / "census_v6" /
                                                "genome_models_s23.tsv")}
    sets["jackhmmer_vert"] = jackhmmer_targets("vertebrata")
    sets["jackhmmer_nonvert"] = set().union(
        *[jackhmmer_targets(g) for g in NONVERT_GROUPS]) or set()

    d = read_json(V1_DELTA)
    sets["app_search"] = set()          # membership held only as a count
    sets["app_search_n"] = d.get("v1_accessions", 0)      # type: ignore
    return sets


# ------------------------------------------------------------ jackhmmer runs

#: (tag, group, convergence table) for every jackhmmer run the project made.
#: The vertebrate runs are read from **S3's merged `convergence.tsv`**, not
#: from the per-run tables the sweep wrote beside the profiles: only the
#: merged one carries the per-round family composition (`n_own`, `n_sister`,
#: `n_offfamily`) that D10's rules are evaluated on, and the per-run tables
#: hold new-target counts alone.
JACK_RUNS = (
    ("itpr1_human", "vertebrata", RESULTS / "census_v3" / "convergence.tsv"),
    ("itpr_fly", "vertebrata", RESULTS / "census_v3" / "convergence.tsv"),
    ("itpr_acanthamoeba", "vertebrata",
     RESULTS / "census_v3" / "convergence.tsv"),
    ("s20_viridiplantae", "viridiplantae",
     S20_DIR / "jackhmmer_convergence_s20.tsv"),
    ("s20_fungi", "fungi", S20_DIR / "jackhmmer_convergence_s20.tsv"),
    ("s20_protista_other", "protista_other",
     S20_DIR / "jackhmmer_convergence_s20.tsv"),
    ("s20_metazoa_nonvert", "metazoa_nonvert",
     S20_DIR / "jackhmmer_convergence_s20.tsv"),
)


def convergence_rounds(tag: str) -> list[dict]:
    """Per-round rows for one run, from whichever committed table holds it."""
    for run_tag, _group, path in JACK_RUNS:
        if run_tag != tag or not path.exists():
            continue
        rows = [r for r in read_tsv(path)
                if r.get("seed_tag", r.get("tag", "")) == tag]
        return sorted(rows, key=lambda r: int(fnum(r["round"], float, 0)))
    return []


def jackhmmer_log(tag: str) -> Path:
    """The raw log, which is where the per-round included lists live."""
    if tag.startswith("s20_"):
        return data_root() / "hmmer" / "s20" / f"jackhmmer_{tag}.log"
    return data_root() / "hmmer" / f"jackhmmer_{tag}.log"


def jackhmmer_verdicts() -> dict[str, dict]:
    """The committed D10 verdict for every run, from the session that made it.

    Read rather than recomputed: `s3_kill` already evaluated these on the
    per-round included lists, and a second evaluation here that disagreed
    would be a fork in what this project means by "converged".
    """
    out: dict[str, dict] = {}
    for tag, blob in (read_json(RESULTS / "census_v3" /
                                "census_v3_stats.json").get("jackhmmer")
                      or {}).items():
        out[tag] = {**blob, "group": "vertebrata"}
    for group, blob in (read_json(S20_DIR / "jackhmmer_verdicts_s20.json")
                        .get("verdicts") or {}).items():
        out[blob.get("tag", f"s20_{group}")] = {**blob, "group": group}
    return out


def jackhmmer_targets(group: str) -> set[str]:
    """Accessions in the **accepted** rounds of every run over one database.

    A run D10 killed contributes only the rounds before the kill, because a
    diverged round's inclusion list is not a search result — which for
    `s20_protista_other` is the difference between 729 family records and
    22,913 accreted targets. The per-round included lists exist only in the
    raw log, so this parses the log rather than the domtblout, which holds
    every round's hits with no way to tell them apart.
    """
    import s3_hmm_lib as H
    import s3_kill as K
    # Imported here rather than at module scope: `s19_universe` imports this
    # module, so a top-level import back would be circular.
    from s19_universe import header_accession
    verdicts = jackhmmer_verdicts()
    out: set[str] = set()
    for tag, run_group, _path in JACK_RUNS:
        if run_group != group:
            continue
        path = jackhmmer_log(tag)
        if not path.exists():
            continue
        parsed = H.parse_jackhmmer_log(path)
        acc_rounds = int(fnum(verdicts.get(tag, {}).get("accepted_rounds"),
                              float, len(parsed["rounds"])))
        # jackhmmer's log names targets in FASTA-header form
        # (`sp|Q14643|ITPR1_HUMAN`), which `acc_key` deliberately leaves
        # untouched because a genome model id also carries pipes. Parsing
        # them with `header_accession` is what makes these comparable with
        # every other channel's accessions; without it the intersection with
        # any database universe is empty and iteration silently looks like
        # it added nothing.
        out |= {header_accession(">" + a)
                for a in K.accepted_targets(parsed["rounds"], acc_rounds)}
    return out


# ------------------------------------------------------------ genome ledger

def ledger_cells() -> list[dict]:
    """One row per genome x cell, with S15a's state and the manifest joined.

    The join is what makes the false-negative definition possible: the ledger
    says whether the sweep placed a locus, S15a says whether the gene is
    there, and they are different questions answered by different evidence.
    """
    man = {r["accession"]: r for r in read_tsv(MANIFEST)}
    state = {(r["accession"], r["cell"]): r for r in read_tsv(MATRIX)}
    rows = []
    for r in read_tsv(LEDGER):
        m = man.get(r["accession"], {})
        st = state.get((r["accession"], r["class"]), {})
        rows.append({
            **r,
            "cell": r["class"],
            "found": int(is_found(r["status"])),
            "taxid": m.get("taxid", ""),
            "level": m.get("level", r.get("assembly_level", "")),
            "refseq_category": m.get("refseq_category", ""),
            "annotation_source": m.get("annotation_source", ""),
            "release_date": m.get("release_date", ""),
            "total_length_bp": int(fnum(m.get("total_length_bp"), float, 0)),
            "scaffold_n50": int(fnum(m.get("scaffold_n50"), float, 0)),
            "contig_n50": int(fnum(r.get("contig_n50"), float, 0)),
            "s15_state": st.get("state", ""),
            "s15_rule": st.get("rule", ""),
            "recon_coverage": fnum(st.get("recon_coverage"), float, 0.0),
        })
    return rows


def control_cells(cells: list[dict] | None = None) -> list[dict]:
    """Every cell whose gene is independently known present.

    Two sources, kept apart by the `control` column so the report can quote
    either: the three ITPR cells S15a states are present (923 of 927 — S15b
    reconstructs **no** loss anywhere in the scope), and the RyR cell, which
    every vertebrate carries as three genes and which the sweep ran as its
    own positive control. A cell in either set that the ledger did not find is
    a false negative of the search, not a biological absence.
    """
    cells = cells if cells is not None else ledger_cells()
    out = []
    for c in cells:
        if c["cell"] == CONTROL_CELL:
            out.append({**c, "control": "ryr_sister"})
        elif c["s15_state"] in PRESENT_STATES:
            out.append({**c, "control": "itpr_present"})
    return out


# ------------------------------------------------------------ bait panel

def bait_meta() -> dict[str, dict]:
    """bait id -> its manifest row, plus the fields the ablation groups on."""
    out = {}
    for r in read_tsv(BAIT_MANIFEST):
        out[r["id"]] = {
            "accession": r["accession"], "family": r["family"],
            "clade": r["clade"], "paralog": r.get("paralog", ""),
            "band": r["band"], "slot": r.get("slot", ""),
            "species": r["species"].split(" (")[0].replace(" ", "_"),
            "order": r.get("order", ""),
            "length": int(fnum(r["length"], float, 0)),
            "origin": r.get("origin", ""), "note": r.get("note", ""),
        }
    return out


# ------------------------------------------------------------ small stats

def median(xs) -> float:
    s = sorted(xs)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval — these are small counts over a few hundred."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (max(0.0, centre - half), min(1.0, centre + half))
