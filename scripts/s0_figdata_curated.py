#!/usr/bin/env python3
"""The review figures' hand-curated tables — and the check that keeps them honest.

Three figures show things no database returns: the sequence of discoveries
(§1), the regulatory input map (§4) and where the disease variants sit (§9).
Those are literature content, so they are curated here rather than measured —
but curated with the same guard the rest of the project uses:

* every row carries the **stable citation keys** of the review sections, and
  this script fails if a key has no row in `references.tsv`, so a figure cannot
  cite something the bibliography does not contain;
* **no year is typed in.** Milestone rows carry a key and a label; the year is
  read from the reference table, so the timeline cannot drift from it;
* variant positions are recorded at the resolution the cited source gives.
  `point` means the source names a residue, `domain` means it localises the
  variant to a named region and nothing finer, `gene` means it does not
  localise it at all. Nothing is promoted to a residue number that a source did
  not state — a review whose subject is the difference between what is
  established and what is repeated cannot invent coordinates for a figure.

    python scripts/s0_figdata_curated.py

Outputs (results/s0_baseline/review_figures/): milestones.tsv, regulators.tsv,
disease_sites.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "s0_baseline" / "review_figures"
REFS = ROOT / "results" / "s0_baseline" / "references.tsv"

# --------------------------------------------------------------- milestones
# (ref key, short label, lane). Lane groups the timeline: what the receptor
# *is*, how it *works*, and what it *does wrong*.
MILESTONES = [
    ("R01", "IP$_3$ releases Ca$^{2+}$ from an intracellular store", "found"),
    ("R21", "P400 is the IP$_3$-binding protein", "found"),
    ("R02", "primary structure — it is the receptor", "found"),
    ("R04", "“similar to ryanodine receptor”", "found"),
    ("R03", "purified receptor conducts Ca$^{2+}$", "found"),
    ("R13", "bell-shaped Ca$^{2+}$ dependence", "works"),
    ("R29", "subtypes differ by tissue and stage", "found"),
    ("R32", "three receptors, not one", "found"),
    ("R08", "release is built from quantal puffs", "works"),
    ("R05", "structure of the IP$_3$-binding core", "works"),
    ("R114", "$Itpr1$ null mice are ataxic", "disease"),
    ("R38", "$ITPR1$ deletion causes SCA15", "disease"),
    ("R45", "$ITPR2$ loss abolishes sweating", "disease"),
    ("R22", "cryo-EM: the tetramer at 4.7 Å", "works"),
    ("R15", "all four sites must be occupied", "works"),
    ("R24", "paralogue structures with Ca$^{2+}$ and IP$_3$", "works"),
    ("R46", "$ITPR3$ variants cause CMT1J", "disease"),
    ("R59", "activation and gating across a ligand series", "works"),
    ("R48", "$ITPR3$ variant causes multisystem disease", "disease"),
    ("R47", "recurrent CMT1J variant in nine families", "disease"),
]

# -------------------------------------------------------------- regulators
# (mechanism, name, effect, section, refs)
REGULATORS = [
    ("intrinsic ligands", "IP$_3$", "activates", "3.1", "R05,R57"),
    ("intrinsic ligands", "Ca$^{2+}$ (co-agonist)", "biphasic", "3.2",
     "R13,R14,R61"),
    ("intrinsic ligands", "ATP", "activates", "3.4", "R16,R70"),
    ("phosphorylation", "PKA", "context", "4.1", "R17"),
    ("phosphorylation", "PKC / CaMKII", "context", "4.1", "R71"),
    ("phosphorylation", "PKG", "context", "4.1", "R72"),
    ("protein partners", "Bcl-2 (BH4)", "inhibits", "4.2", "R18,R73"),
    ("protein partners", "Mcl-1", "inhibits", "4.2", "R76"),
    ("protein partners", "Bcl-X$_L$", "context", "4.2", "R75"),
    ("protein partners", "BIRD-2 peptide", "activates", "4.2", "R74,R132"),
    ("protein partners", "IRBIT", "inhibits", "4.4", "R19,R79"),
    ("protein partners", "CaBPs", "activates", "4.5", "R80"),
    ("protein partners", "calmodulin", "context", "4.5", "R81"),
    ("luminal / ER", "ERp44 (redox, pH)", "inhibits", "4.3", "R20"),
    ("turnover", "ubiquitination", "inhibits", "4.6", "R82"),
    ("turnover", "RNF170 (ERAD)", "inhibits", "4.6", "R83"),
    ("disease proteins", "presenilins (FAD)", "activates", "4.7",
     "R84"),
]

# ------------------------------------------------------------ disease sites
# kind: point (a residue the source names) | domain (a region it names) |
# gene (not localised by the source). `where` is a residue for `point` and a
# Pfam accession for `domain`.
DISEASE = [
    ("ITPR1", "SCA15/16 — heterozygous deletion", "gene", "",
     "haploinsufficiency", "R38,R118,R119"),
    # §9.1 establishes SCA29 as dominant and missense, and shows one variant
    # is a *gain* of function — so the direction of effect is explicitly not
    # settled. Recorded as unresolved rather than assimilated to the
    # dominant-negative mechanism of the Gillespie and ITPR3 variants.
    ("ITPR1", "SCA29 — missense", "gene", "", "unresolved", "R41,R42"),
    ("ITPR1", "gain-of-function ataxia", "domain", "PF08709",
     "gain of function", "R120"),
    ("ITPR1", "Gillespie — de novo", "domain", "PF00520",
     "dominant negative", "R44"),
    ("ITPR1", "Gillespie — biallelic", "gene", "", "recessive loss", "R43"),
    ("ITPR2", "isolated anhidrosis — homozygous missense", "domain",
     "PF00520", "recessive loss", "R45"),
    ("ITPR3", "CMT1J — p.Thr1424Met", "point", "1424", "dominant negative",
     "R47"),
    ("ITPR3", "multisystem disorder — p.Arg2524Cys", "point", "2524",
     "dominant negative", "R48"),
    ("ITPR3", "CMT1J — segregating families", "gene", "",
     "dominant negative", "R46,R51"),
]


def read_refs() -> dict[str, dict]:
    lines = REFS.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    return {r["ref_id"]: r for r in
            (dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln.strip())}


def write(path: Path, cols: list[str], rows: list[tuple]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(c) for c in r) + "\n")
    print(f"  wrote {path.name}  ({len(rows)} rows)")


def main() -> int:
    refs = read_refs()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    used = {k for k, _, _ in MILESTONES}
    for table in (REGULATORS, DISEASE):
        for row in table:
            used.update(k.strip() for k in row[-1].split(",") if k.strip())
    missing = sorted(k for k in used if k not in refs)
    if missing:
        print(f"ERROR: figure rows cite keys with no reference: {missing}",
              file=sys.stderr)
        return 2
    print(f"  {len(used)} citation keys, all present in references.tsv")

    write(OUT_DIR / "milestones.tsv",
          ["ref_id", "year", "lane", "label"],
          [(k, refs[k]["year"], lane, lab) for k, lab, lane in MILESTONES])
    write(OUT_DIR / "regulators.tsv",
          ["mechanism", "name", "effect", "section", "refs"], REGULATORS)
    write(OUT_DIR / "disease_sites.tsv",
          ["gene", "label", "kind", "where", "mechanism", "refs"], DISEASE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
