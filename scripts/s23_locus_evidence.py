"""s23_locus_evidence.py — what independent evidence says about each locus.

Split out of `s23_calibrate_loci` to keep both under 500 lines. This is the
part that decides, for every recorded cluster, what something *other than the
alignment score* says about it — which is the whole basis on which a threshold
on the alignment score can be calibrated at all.

**Two axes, and they are not equally good.**

`_evidence` reads the assembly's own annotation. It is the independent one —
it knows nothing about this project's baits, profiles or family definition —
and outside the vertebrates it is nearly empty: 21 of 917 clusters sit on a
gene whose name says anything, because most gene models here carry locus tags
(*Chlamydomonas* files its receptor as `CHLRE_16g665450v5`). That emptiness is
why the second axis exists, and it is reported rather than smoothed over.

`profile_evidence` scores each cluster's translated model with `itpr.hmm` and
`ryr.hmm` (D23) — a different algorithm on different data, but **not**
independent of the family definition the way an annotation is. Loci are
re-derived from each genome's *archived* miniprot GFF, the discipline S2 and
S3 apply to their raw output: the translations exist only in that file, so
being able to re-parse them offline is what makes the calibration
re-checkable without re-running the sweep.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from src.utils.data_root import require_data_root          # noqa: E402
import s23_bait_spec as spec                               # noqa: E402
import s23_calibration as cal                              # noqa: E402
import s5_sweep_lib as lib                                 # noqa: E402
from s23_run_sweep import ensure_bait_files                # noqa: E402
from s5_census_v4 import score_models                      # noqa: E402
from s5_classify import name_family                        # noqa: E402

ANNOT_CDS_FRAC = 0.50

def load_summaries() -> list[dict]:
    root = require_data_root() / "s23_sweep"
    out = []
    for path in sorted(root.glob("*/summary.json")):
        try:
            out.append(json.loads(path.read_text()))
        except json.JSONDecodeError:
            continue
    return out


#: Gene "names" that are not names — RefSeq `LOC` placeholders, locus tags and
#: the assorted ways an unnamed model is written. They carry **no** evidence
#: about a locus in either direction, and lumping them in with real symbols is
#: what makes an outside-the-vertebrates calibration go wrong: *Chlamydomonas*
#: files its receptor as `CHLRE_16g665450v5` and *Strongylocentrotus* files
#: its as `LOC594527`, so treating an uninformative name as a contradiction
#: would put two correctly-recovered genes in the junk population.
_PLACEHOLDER = re.compile(
    r"""^(?:
          loc\d+                 |   # RefSeq placeholder
          orf\d+                 |
          gene[-_]?\d+           |
          g\d+                   |
          [a-z0-9]+[_.][a-z0-9.]*\d[a-z0-9.]*   # locus tag: CHLRE_16g665450v5
        )$""",
    re.IGNORECASE | re.VERBOSE)

_VAGUE = ("uncharacter", "hypothetical", "unnamed", "predicted protein")


def is_placeholder(name: str) -> bool:
    low = (name or "").strip().lower()
    if not low:
        return True
    return bool(_PLACEHOLDER.match(low)) or any(v in low for v in _VAGUE)


def _evidence(locus: dict) -> tuple[str, str]:
    """(class, the gene name it rests on) for one locus, from the annotation.

    The classes are what the *assembly's own annotation* says about a cluster,
    which is evidence the alignment score did not produce — that is what makes
    it usable to calibrate a threshold on alignment identity.

    `confirmed`     an annotated gene the assembly names for this family covers
                    the alignment's exons.
    `sister`        the covering gene is named for the ryanodine receptors —
                    D14's sharpest possible failure, so it is counted apart
                    from ordinary junk rather than averaged into it.
    `contradicted`  the alignment's exons fall inside an annotated gene with a
                    real symbol for something else. The S23a pilot's junk
                    clusters land on `Myo81F`, `Sdb`, `Prosap`, `mam` — this is
                    the population a floor has to exclude.
    `unnamed`       a gene covers it but its "name" is a placeholder or a locus
                    tag, so the annotation says nothing either way.
    `no_annotation` nothing annotated covers the alignment's exons.
    """
    genes = locus.get("overlapping_genes") or []
    good = [g for g in genes if g.get("frac_cds", 0) >= ANNOT_CDS_FRAC]
    if not good:
        return "no_annotation", ""
    if locus.get("annot_names_family"):
        return "confirmed", (locus.get("annot_gene") or {}).get("name", "")
    for g in good:
        if name_family(g.get("name", "")) == "RYR":
            return "sister", g["name"]
    real = [g for g in good if not is_placeholder(g.get("name", ""))]
    if real:
        return "contradicted", real[0]["name"]
    return "unnamed", (good[0].get("name") or "")


def collect(summaries: list[dict],
            profiles: dict[str, dict] | None = None) -> list[dict]:
    """Every recorded locus, called or not, with its evidence class."""
    profiles = profiles or {}
    rows = []
    for s in summaries:
        for called, pool in ((1, s.get("loci") or []),
                             (0, s.get("below_gate_loci")
                                 or s.get("below_floor_loci") or [])):
            for d in pool:
                klass, gene = _evidence(d)
                key = (f"{s['accession']}|"
                       f"{_addr(d['contig'], d['start'], d['end'])}")
                v = profiles.get(key) or {}
                pcall = v.get("profile_call", "")
                rows.append({
                    "profile_call": pcall,
                    "profile_confidence": v.get("profile_confidence", ""),
                    "itpr_score": v.get("itpr_score", ""),
                    "ryr_score": v.get("ryr_score", ""),
                    "profile_evidence": ("profile_confirmed" if pcall == "ITPR"
                                         else "profile_contradicted"
                                         if pcall == "RYR" else ""),
                    "accession": s["accession"], "organism": s["organism"],
                    "group": s.get("group", ""), "phylum": s.get("phylum", ""),
                    "called": called, "identity": d.get("identity", 0.0),
                    "coverage": d.get("coverage", 0.0),
                    "grade": d.get("grade", ""), "bait": d.get("bait", ""),
                    "band": d.get("band", ""),
                    "span_bp": d.get("span_bp", 0),
                    "cds_footprint_bp": d.get("cds_footprint_bp", 0),
                    "span_inflation": d.get("span_inflation", 0.0),
                    "evidence": klass, "annot_gene": gene,
                    "status": s.get("status", ""),
                    "n_full": s.get("n_full", 0)})
    return rows


# ------------------------------------------------------------ the profiles

def _addr(contig: str, start: int, end: int) -> str:
    return f"{contig}:{start}-{end}"


def profile_evidence(summaries: list[dict]) -> dict[str, dict]:
    """Score every recorded locus's model with itpr.hmm / ryr.hmm.

    The loci are re-derived from each genome's **archived** `miniprot.gff`,
    the same discipline S2 and S3 apply to their raw output: the translations
    exist only in that file, and being able to re-parse them offline is what
    makes this calibration re-checkable without re-running the sweep.
    """
    root = require_data_root() / "s23_sweep"
    _panel, _rescue, _seqs, meta = ensure_bait_files()
    models: dict[str, str] = {}
    for s in summaries:
        gff = root / s["accession"] / "miniprot.gff"
        if not gff.exists() or gff.stat().st_size == 0:
            continue
        alns = lib.parse_miniprot_gff(gff, meta)
        loci = lib.filter_loci(lib.cluster_loci(alns),
                              min_identity=cal.RECORD_MIN_IDENTITY)
        for L in loci:
            if L.family != spec.FAMILY_ITPR:
                continue
            best = L.best_of_family(spec.FAMILY_ITPR) or L.best
            seq = (best.translation or "").replace("*", "")
            if len(seq) < 60:
                continue
            models[f"{s['accession']}|{_addr(L.contig, L.start, L.end)}"] = seq
    if not models:
        return {}
    return score_models(models)
