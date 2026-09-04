"""S2 step 3 — a positive ITPR-or-RYR call on every enumerated record.

The census's search space is defined by signatures both families carry, so
roughly half of what the enumeration returns is a ryanodine receptor. This
module makes the separation, and it makes it as a **positive test** on
domain architecture (D14): a record is called RYR because it carries RyR
signatures, and ITPR because it carries the complete IP3-receptor
architecture and none of them. Length is recorded on every row and enters
only as support for a medium-confidence call; it never decides one, because
the band was chosen to exclude RyRs and would therefore beg the question.

**Why absence is allowed to mean something here.** Not carrying a RyR
signature is, on its own, weak evidence — the annotation could simply be
partial. It becomes strong evidence when the record's IP3-receptor
architecture is *complete*: a RyR annotated well enough to show all five
shared signatures would also show its own. The two rules are therefore
paired, and the pairing is audited rather than asserted —
`rule_audit.tsv` scores the architecture call against gene symbols the
architecture rule never sees.

Records that satisfy neither positive test stay `unassigned` and are
profiled in `unassigned_profile.tsv`. Records whose architecture and symbol
disagree are called on the architecture and listed in `conflicts.tsv`;
those are the mis-annotation candidates S18 will want.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.data_root import require_data_root           # noqa: E402
from src.utils.family import MAX_LENGTH_AA, MIN_LENGTH_AA   # noqa: E402
from scripts.s2_lib import (                                # noqa: E402
    CONTEXT_PFAMS, OUT_DIR, RYR_PFAMS, SEED_PFAMS,
    read_tsv, signature_names, write_tsv,
)

#: The architecture a completely annotated IP3 receptor shows: the three
#: seeds plus MIR and the pore. Human ITPR1/2/3 carry exactly this set.
ITPR_ARCH = set(SEED_PFAMS) | set(CONTEXT_PFAMS)
RYR_SET = set(RYR_PFAMS)

CENSUS_COLS = [
    "accession", "call", "confidence", "reason", "arch_call", "symbol_class",
    "conflict", "gene", "protein_name", "species", "taxon_id", "group",
    "kingdom", "phylum", "class", "order", "family", "length", "in_band",
    "fragment",
    "reviewed", "protein_existence", "in_alphafold", "seed_pfams", "pfams",
    "n_itpr_arch", "ryr_pfams", "source",
]


def symbol_class(gene: str) -> str:
    """Gene symbol → 'ITPR' / 'RYR' / '' — evidence, not the call.

    Normalised to lower-case alphanumerics first, so `itpr2.1`, `ITPR2-1`
    and `Itpr2` are one symbol. `itr` is included because the nematode IP3
    receptor is `itr-1`; `ip3r`/`insp3r` because several databases use the
    protein name as the symbol.
    """
    g = "".join(c for c in gene.lower() if c.isalnum())
    if not g:
        return ""
    if g.startswith(("itpr", "ip3r", "insp3r", "itr")):
        return "ITPR"
    if g.startswith("ryr"):
        return "RYR"
    return ""


def architecture_call(pfams: set[str], length: int,
                      arch_known: bool) -> tuple[str, str, list[str]]:
    """(call, confidence, components) from domain architecture + length."""
    comps: list[str] = []
    ryr_hits = sorted(pfams & RYR_SET)
    n_itpr = len(pfams & ITPR_ARCH)
    in_band = MIN_LENGTH_AA <= length <= MAX_LENGTH_AA

    if ryr_hits:
        comps.append("ryr_pfam=" + "+".join(ryr_hits))
        return "RYR", ("high" if len(ryr_hits) >= 2 else "medium"), comps
    if not arch_known:
        # Enumerated by InterPro but absent from the UniProt sweep: we know
        # which seeds it carries, not what else it carries. Absence of a
        # RyR signature cannot be read as evidence here.
        comps.append("arch_unknown")
        return "unassigned", "none", comps
    if ITPR_ARCH <= pfams:
        comps.append("itpr_arch_complete")
        return "ITPR", "high", comps
    if n_itpr >= 3:
        comps.append(f"itpr_arch_partial={n_itpr}/5")
        if in_band:
            comps.append("length_in_band")
            return "ITPR", "medium", comps
        return "unassigned", "none", comps
    comps.append(f"itpr_arch_weak={n_itpr}/5")
    return "unassigned", "none", comps


def final_call(arch: str, arch_conf: str, comps: list[str],
               sym: str) -> tuple[str, str, str, int]:
    """Fold symbol evidence into the architecture call."""
    comps = list(comps)
    conflict = 0
    call, conf = arch, arch_conf
    if sym:
        comps.append(f"symbol={sym}")
    if arch != "unassigned" and sym and sym != arch:
        conflict = 1
        conf = "medium" if conf == "high" else conf
        comps.append("SYMBOL_CONFLICT")
    elif arch != "unassigned" and sym == arch and conf == "medium":
        conf = "high"
    elif arch == "unassigned" and sym:
        call, conf = sym, "low"
    return call, conf, ";".join(comps), conflict


def build(uniprot: list[dict], interpro: list[dict]) -> list[dict]:
    by_acc = {r["accession"]: r for r in uniprot}
    ip_by_acc = {r["accession"]: r for r in interpro}
    rows = []
    for acc in sorted(set(by_acc) | set(ip_by_acc)):
        u = by_acc.get(acc)
        i = ip_by_acc.get(acc)
        source = ("both" if u and i else "uniprot_only" if u else "interpro_only")
        pfams = set((u or {}).get("pfams", "").split(";")) - {""}
        seeds = (i or {}).get("seed_pfams") or ";".join(
            sorted(pfams & set(SEED_PFAMS)))
        try:
            length = int((u or i or {}).get("length") or 0)
        except ValueError:
            length = 0
        gene = (u or {}).get("gene") or (i or {}).get("gene") or ""
        sym = symbol_class(gene)
        arch, arch_conf, comps = architecture_call(pfams, length, u is not None)
        call, conf, reason, conflict = final_call(arch, arch_conf, comps, sym)
        rows.append({
            "accession": acc, "call": call, "confidence": conf,
            "reason": reason, "arch_call": arch, "symbol_class": sym,
            "conflict": conflict, "gene": gene,
            "protein_name": (u or {}).get("protein_name", ""),
            "species": (u or i or {}).get("species", ""),
            "taxon_id": (u or i or {}).get("taxon_id", ""),
            "group": (u or {}).get("group", ""),
            "kingdom": (u or {}).get("kingdom", ""),
            "phylum": (u or {}).get("phylum", ""),
            "class": (u or {}).get("class", ""),
            "order": (u or {}).get("order", ""),
            "family": (u or {}).get("family", ""),
            "length": length,
            "in_band": int(MIN_LENGTH_AA <= length <= MAX_LENGTH_AA),
            "fragment": (u or {}).get("fragment", ""),
            "reviewed": (u or {}).get("reviewed", (i or {}).get("reviewed", "")),
            "protein_existence": (u or {}).get("protein_existence", ""),
            # Only the InterPro rows carry this; S11 needs it and it is
            # free here, so it rides along rather than being re-fetched.
            "in_alphafold": (i or {}).get("in_alphafold", ""),
            "seed_pfams": seeds, "pfams": ";".join(sorted(pfams)),
            "n_itpr_arch": len(pfams & ITPR_ARCH),
            "ryr_pfams": ";".join(sorted(pfams & RYR_SET)),
            "source": source,
        })
    return rows


def audit(rows: list[dict], names: dict[str, str]) -> list[dict]:
    """Score the architecture call against the symbols it never sees."""
    out = []
    labelled = [r for r in rows if r["symbol_class"]]
    for sym in ("ITPR", "RYR"):
        sub = [r for r in labelled if r["symbol_class"] == sym]
        agree = sum(1 for r in sub if r["arch_call"] == sym)
        wrong = sum(1 for r in sub if r["arch_call"] not in (sym, "unassigned"))
        una = sum(1 for r in sub if r["arch_call"] == "unassigned")
        out.append({
            "test": f"architecture call vs symbol {sym}",
            "n": len(sub), "agree": agree, "disagree": wrong,
            "unassigned": una,
            "accuracy_of_decided": (f"{agree / (agree + wrong):.4f}"
                                    if agree + wrong else ""),
        })
    # The conditional claim the RyR rule rests on: SPRY inside this search
    # space means RyR. Tested, not assumed.
    itpr_lab = [r for r in labelled if r["symbol_class"] == "ITPR"]
    spry_itpr = [r for r in itpr_lab if "PF00622" in r["pfams"]]
    out.append({
        # Phrased as the expectation, not as a count, so "agree" reads the
        # same way here as in every other row of the table.
        "test": "no symbol-labelled ITPR carries SPRY (PF00622)",
        "n": len(itpr_lab), "agree": len(itpr_lab) - len(spry_itpr),
        "disagree": len(spry_itpr), "unassigned": 0,
        "accuracy_of_decided": f"{1 - len(spry_itpr) / len(itpr_lab):.4f}"
                               if itpr_lab else "",
    })
    for pid in RYR_PFAMS:
        hit = [r for r in labelled if pid in r["pfams"]]
        out.append({
            "test": f"{names.get(pid, pid)} ({pid}) → symbol RYR",
            "n": len(hit),
            "agree": sum(1 for r in hit if r["symbol_class"] == "RYR"),
            "disagree": sum(1 for r in hit if r["symbol_class"] == "ITPR"),
            "unassigned": 0,
            "accuracy_of_decided": (
                f"{sum(1 for r in hit if r['symbol_class'] == 'RYR') / len(hit):.4f}"
                if hit else ""),
        })
    return out


def seed_contribution(rows: list[dict]) -> list[dict]:
    """What a census built on one signature alone would have missed.

    `PF08709` is the family's *defining* signature — the IP3-binding core —
    and the obvious thing to build a census on. This table is the reason
    the enumeration is a union of three: for each seed, the records that do
    not carry it, and how many of those this census nonetheless calls ITPR.
    """
    out = []
    for seed in SEED_PFAMS:
        absent = [r for r in rows if seed not in r["seed_pfams"]]
        itpr = [r for r in absent if r["call"] == "ITPR"]
        out.append({
            "pfam_id": seed,
            "present": len(rows) - len(absent),
            "absent": len(absent),
            "absent_called_itpr": len(itpr),
            "taxa_lost": len({r["taxon_id"] for r in itpr}),
        })
    combos: dict[tuple, int] = {}
    for r in rows:
        k = tuple(sorted(r["seed_pfams"].split(";")))
        combos[k] = combos.get(k, 0) + 1
    for k, v in sorted(combos.items(), key=lambda kv: -kv[1]):
        out.append({"pfam_id": "+".join(k), "present": v, "absent": "",
                    "absent_called_itpr": "", "taxa_lost": ""})
    return out


def short_complete(rows: list[dict]) -> list[dict]:
    """Records with the *complete* ITPR architecture and an impossible length.

    The call is architecture-based and length never decides it — which is
    the right rule, and it means the band is free to catch something else.
    A protein carrying all five IP3-receptor signatures inside 2,000
    residues has ~700 fewer than the shortest real family member; the
    architecture is intact and the model is not. UniProt's own `Fragment`
    flag does not catch these, because a truncated gene model submitted as
    a whole protein is not marked as a fragment.
    """
    out = [r for r in rows
           if r["call"] == "ITPR" and r["confidence"] == "high"
           and int(r["length"]) < MIN_LENGTH_AA]
    return sorted(out, key=lambda r: int(r["length"]))


def lineage_calls(rows: list[dict]) -> list[dict]:
    """Calls per group × phylum — where the family is and is not.

    The coarse group tally hides the shape of the answer: "46 fungal
    records" is a very different statement from "every fungal call is in an
    early-diverging phylum and Dikarya has none". This is a statement about
    what UniProt holds, not about what genomes hold; S20 and S23 are what
    turn it into the latter.
    """
    agg: dict[tuple, dict] = {}
    for r in rows:
        k = (r["group"] or "unclassified", r["phylum"] or "(no phylum)")
        a = agg.setdefault(k, {"group": k[0], "phylum": k[1], "records": 0,
                               "itpr": 0, "itpr_high": 0, "ryr": 0,
                               "_taxa": set(), "_itaxa": set()})
        a["records"] += 1
        a["_taxa"].add(r["taxon_id"])
        if r["call"] == "ITPR":
            a["itpr"] += 1
            a["_itaxa"].add(r["taxon_id"])
            a["itpr_high"] += int(r["confidence"] == "high")
        elif r["call"] == "RYR":
            a["ryr"] += 1
    out = []
    for a in agg.values():
        a["taxa"] = len(a.pop("_taxa"))
        a["itpr_taxa"] = len(a.pop("_itaxa"))
        out.append(a)
    return sorted(out, key=lambda r: (r["group"], -r["itpr"], -r["records"]))


def summarise(rows: list[dict]) -> list[dict]:
    keys: dict[tuple[str, str], int] = {}
    for r in rows:
        k = (r["call"], r["confidence"])
        keys[k] = keys.get(k, 0) + 1
    return [{"call": c, "confidence": q, "records": n}
            for (c, q), n in sorted(keys.items(), key=lambda kv: -kv[1])]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    names = signature_names()
    root = require_data_root(OUT_DIR.parents[1])
    uniprot = read_tsv(root / "raw_api" / "uniprot" / "uniprot_records.tsv")
    ip_path = root / "raw_api" / "interpro" / "interpro_records.tsv"
    interpro = read_tsv(ip_path) if ip_path.exists() else []
    rows = build(uniprot, interpro)

    write_tsv(OUT_DIR / "census_v2.tsv", rows, CENSUS_COLS)
    write_tsv(OUT_DIR / "call_summary.tsv", summarise(rows),
              ["call", "confidence", "records"])
    write_tsv(OUT_DIR / "short_complete.tsv", short_complete(rows),
              ["accession", "gene", "species", "group", "length", "fragment",
               "protein_existence", "reviewed", "pfams"])
    write_tsv(OUT_DIR / "lineage_calls.tsv", lineage_calls(rows),
              ["group", "phylum", "records", "taxa", "itpr", "itpr_high",
               "itpr_taxa", "ryr"])
    write_tsv(OUT_DIR / "seed_contribution.tsv", seed_contribution(rows),
              ["pfam_id", "present", "absent", "absent_called_itpr",
               "taxa_lost"])
    write_tsv(OUT_DIR / "rule_audit.tsv", audit(rows, names),
              ["test", "n", "agree", "disagree", "unassigned",
               "accuracy_of_decided"])
    write_tsv(OUT_DIR / "conflicts.tsv",
              [r for r in rows if r["conflict"]], CENSUS_COLS)

    una = [r for r in rows if r["call"] == "unassigned"]
    prof: dict[tuple, int] = {}
    for r in una:
        k = (r["source"], r["seed_pfams"], r["n_itpr_arch"],
             "in_band" if r["in_band"] else "out_of_band")
        prof[k] = prof.get(k, 0) + 1
    write_tsv(OUT_DIR / "unassigned_profile.tsv",
              [{"source": s, "seed_pfams": sp, "n_itpr_arch": n,
                "length_band": b, "records": c}
               for (s, sp, n, b), c in sorted(prof.items(), key=lambda kv: -kv[1])],
              ["source", "seed_pfams", "n_itpr_arch", "length_band", "records"])

    n = len(rows)
    for c in ("ITPR", "RYR", "unassigned"):
        k = sum(1 for r in rows if r["call"] == c)
        print(f"  {c:<11s} {k:6d}  ({100 * k / n:.1f} %)")
    print(f"  conflicts   {sum(r['conflict'] for r in rows):6d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
