"""Shared helpers for the S1 control benchmark (see s1_benchmark.py).

Four responsibilities, nothing else:
  * prove that `analyse(use_mafft=True)` really shells out to MAFFT
    (MafftTracer wraps subprocess.run and records every mafft call);
  * look candidates up in a DiscoveryReport by label or dedup key;
  * decode which scorer components fired, from the evidence dict;
  * compute the labelled-bait margin that operationalises D14/D7 —
    identity to the best human ITPR bait minus identity to the best human
    RyR bait, per panel member.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

#: The labelled baits the D14 margin is measured against (UniProt accessions
#: of the human paralogs; the sister panel comes from utils/family.py).
ITPR_BAIT_ACCS = ["Q14643", "Q14571", "Q14573"]
RYR_BAIT_ACCS = ["P21817", "Q92736", "Q15413"]


class MafftTracer:
    """Record every subprocess call to mafft while active — hard evidence
    that the use_mafft=True code path ran the external binary rather than
    falling back to the star alignment (`_mafft` falls back silently on a
    non-zero return code, so the return code is part of the evidence)."""

    def __init__(self) -> None:
        self.calls: list[dict] = []
        self._orig = None

    def __enter__(self) -> "MafftTracer":
        self._orig = subprocess.run

        def traced(cmd, *args, **kwargs):
            result = self._orig(cmd, *args, **kwargs)
            prog = str(cmd[0]) if isinstance(cmd, (list, tuple)) and cmd else str(cmd)
            if "mafft" in prog:
                self.calls.append({
                    "argv": list(cmd) if isinstance(cmd, (list, tuple)) else [prog],
                    "returncode": result.returncode,
                    "stdout_bytes": len(result.stdout or ""),
                })
            return result

        subprocess.run = traced
        return self

    def __exit__(self, *exc) -> None:
        subprocess.run = self._orig


def dedup_key(v) -> tuple[str, str, Optional[int]]:
    """The key `discover_novel_paralogs` collapses candidates on."""
    return (v.gene_symbol.upper(), v.species, v.length_aa)


def score_lookup(report):
    """label → Candidate, with a dedup-key fallback (discovery collapses
    same (symbol, species, length) hits into one representative)."""
    by_label = {c.label: c for c in report.candidates}
    by_key = {dedup_key(c.variant): c for c in report.candidates}

    def lookup(v):
        lbl = f"{v.source}|{v.accession}|{v.gene_symbol}"
        return by_label.get(lbl) or by_key.get(dedup_key(v))

    return lookup


def components(cand) -> str:
    """Compact list of which scorer components fired, read back out of the
    evidence text. Mirrors src/discovery/candidates.py component by
    component; a change there must be reflected here."""
    if cand is None:
        return ""
    ev = cand.evidence
    flags = []
    if "within family band" in ev.get("size", ""):
        flags.append("size+15")
    dom = ev.get("domain", "")
    if dom.startswith("Pfam signature present"):
        flags.append("pfam+20")
    elif dom.startswith("MSA-signature fallback"):
        flags.append("sig+15")
    if "≥" in ev.get("fold", ""):
        flags.append("fold+20")
    if "paralog twilight zone" in ev.get("outlier", ""):
        flags.append("outlier+20")
    if ev.get("cluster", "").startswith("outside known paralog clusters"):
        flags.append("cluster+10")
    if ev.get("breadth", "").startswith("present in"):
        flags.append("breadth+15")
    if "homology" in ev:
        flags.append("homology+10")
    if "annotation" in ev:
        flags.append("split+10")
    if "sister_cap" in ev:
        flags.append("SISTER→39")
    if "gate" in ev:
        flags.append("GATED→39")
    return ",".join(flags)


def sister_margin(cand) -> str:
    """The recorded ITPR-vs-sister margin ('… margin -53%')."""
    if cand is None:
        return ""
    ev = cand.evidence.get("sister_margin", "")
    return ev.split("margin ")[-1].strip() if "margin " in ev else ""



def nearest_known_identity(cand) -> str:
    """The identity-to-nearest-known figure the outlier component used,
    pulled back out of the evidence string ('… = 62% —')."""
    if cand is None:
        return ""
    ev = cand.evidence.get("outlier", "")
    if "= " not in ev:
        return ""
    tail = ev.split("= ")[-1]
    return tail.split("%")[0].strip() if "%" in tail else ""


def bait_margins(distances, label_by_acc: dict[str, str]) -> dict[str, dict]:
    """Per-label ITPR-vs-RyR labelled-bait margin (roadmap D7/D14).

    For every row of the distance matrix: the best identity to any human
    ITPR bait, the best identity to any human RyR bait, the call, and the
    margin between them. A bait's own row excludes itself.
    """
    if distances is None or not distances.labels:
        return {}
    idx = {lbl: i for i, lbl in enumerate(distances.labels)}
    itpr_labels = [label_by_acc[a] for a in ITPR_BAIT_ACCS if a in label_by_acc]
    ryr_labels = [label_by_acc[a] for a in RYR_BAIT_ACCS if a in label_by_acc]

    def best(i: int, bait_labels: list[str]) -> tuple[str, float]:
        hit, score = "", 0.0
        for b in bait_labels:
            j = idx.get(b)
            if j is None or j == i:
                continue
            v = distances.identity[i][j]
            if v > score:
                hit, score = b, v
        return hit, score

    out: dict[str, dict] = {}
    for lbl, i in idx.items():
        i_hit, i_id = best(i, itpr_labels)
        r_hit, r_id = best(i, ryr_labels)
        out[lbl] = {
            "best_itpr_bait": i_hit, "id_to_itpr": i_id,
            "best_ryr_bait": r_hit, "id_to_ryr": r_id,
            "margin": i_id - r_id,
            "call": "ITPR" if i_id > r_id else ("RYR" if r_id > i_id else "tie"),
        }
    return out


def write_tsv(path: Path, rows: list[dict], cols: list[str] | None = None) -> int:
    """Write dict rows as TSV. Column order from `cols`, else first row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\t".join(cols or []) + "\n")
        return 0
    cols = cols or list(rows[0].keys())
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return len(rows)
