"""S18's protein half — what the protein databases call a record, tested.

The genome half asks whether an assembly's annotation holds a gene. This half
asks the same two questions of the protein records themselves: does the family
this project calls a record agree with the family its **name** claims, and —
inside the vertebrates, where ITPR1/2/3 exist — does the paralog agree.

Four rules fix the scope, each a positive test:

* **P1** the record is called ITPR or RYR by census v6, this project's own
  merged call (D23), not by its own annotation;
* **P2** it comes from a protein database (`InterPro` or the S3 `HMM` sweep)
  and not from a genome model, because a genome model is the genome half's
  evidence and auditing it here would be this project marking its own work;
* **P3** it is full length — at least `family.MIN_LENGTH_AA` and not flagged a
  fragment — since a 400 aa fragment's paralog is not a naming failure;
* **P4** the **paralog** question is asked only of vertebrates. ITPR1/2/3 are
  a 2R product (S6's rule, D30): a protist record whose UniProt name reads
  "receptor type 2" carries that number by annotation transfer, and scoring it
  as a wrong-paralog error would assert what S7 was run to test. Every record
  is still audited for **family**, where the ITPR/RyR question is real
  everywhere.

The call itself is D7 with S3's own number: the best bit score in each family
against the committed 38-bait panel, assigned to the winner only when it beats
the loser by more than `s3_assign.REL_MARGIN` of its own score. Relative and
not absolute, for `s3_assign`'s reason — bit scores scale with alignable
length, so a fixed gap calls every full-length protein and no short one.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402
from s3_assign import REL_MARGIN                                  # noqa: E402
from src.utils import family as FAM                               # noqa: E402

MIN_FULL_LENGTH = FAM.MIN_LENGTH_AA
PROTEIN_SOURCES = ("InterPro", "HMM")


# --------------------------------------------------------------------------
# scope
# --------------------------------------------------------------------------
def in_scope(row: dict) -> tuple[bool, str]:
    """P1–P3, with the rule that excluded a record named in its own row."""
    if row.get("call") not in ("ITPR", "RYR"):
        return False, "P1 not a family call in census v6"
    if row.get("source") not in PROTEIN_SOURCES:
        return False, f"P2 not a protein-database record ({row.get('source')})"
    if row.get("fragment") == "1":
        return False, "P3 flagged a fragment"
    try:
        n = int(row.get("length") or 0)
    except ValueError:
        return False, "P3 no length"
    if n < MIN_FULL_LENGTH:
        return False, f"P3 {n} aa < {MIN_FULL_LENGTH} (not full length)"
    return True, ""


def paralog_askable(row: dict) -> bool:
    """P4 — the paralog question is a vertebrate question."""
    return (row.get("group") == "Vertebrata"
            or row.get("tax_group") == "Vertebrata")


# --------------------------------------------------------------------------
# the call
# --------------------------------------------------------------------------
def best_by_label(hits: list[tuple[str, float]], labels: dict[str, dict]
                  ) -> tuple[dict[str, float], dict[str, float]]:
    """Best bit score per family and per paralog over one query's hits."""
    fam: dict[str, float] = {}
    par: dict[str, float] = {}
    for sid, bits in hits:
        lab = labels.get(L.bait_id(sid))
        if not lab:
            continue
        f = lab.get("family") or ""
        p = lab.get("paralog") or ""
        if f:
            fam[f] = max(fam.get(f, 0.0), bits)
        if p and p not in ("ITPR", "RYR"):
            par[p] = max(par.get(p, 0.0), bits)
    return fam, par


def _call(scores: dict[str, float], margin: float) -> tuple[str, float, str]:
    """Winner, relative margin, and the runner-up it had to beat."""
    if not scores:
        return "no_call", 0.0, ""
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top, best = ranked[0]
    runner, second = (ranked[1] if len(ranked) > 1 else ("", 0.0))
    rel = round((best - second) / best, 4) if best > 0 else 0.0
    if rel <= margin:
        return "no_call", rel, runner
    return top, rel, runner


def _panel_guard(verdict: str, claim: str, family: str,
                 callable_by_panel: dict[str, set[str]]) -> str:
    """A paralog the panel cannot reach is not a naming error.

    The committed panel carries no RYR3 bait — S5's `unfilled_slots.tsv`
    records the six slots no labelled full-length record could fill — so a
    RYR3 record is called RYR1 or RYR2 by this instrument however well it is
    annotated, and confidently: the 74 RYR3 symbols here are called RYR2 at
    margins up to 0.36. Scoring them as the database naming the wrong paralog
    would report the bait panel's reach as the databases' error rate, which is
    the same mistake `s15_synteny_reach` refused to make when a caller could
    not arrive.
    """
    if verdict != "correct_family_wrong_paralog":
        return verdict
    reach = callable_by_panel.get(family, set())
    return "paralog_not_callable" if claim and claim not in reach else verdict


def audit_record(row: dict, hits: list[tuple[str, float]],
                 labels: dict[str, dict],
                 margin: float = REL_MARGIN,
                 callable_by_panel: dict[str, set[str]] | None = None) -> dict:
    """One record's sequence call, its two name verdicts and their margins."""
    fam_s, par_s = best_by_label(hits, labels)
    fam_call, fam_rel, fam_runner = _call(fam_s, margin)
    askable = paralog_askable(row)
    if fam_call in ("ITPR", "RYR") and askable:
        within = {k: v for k, v in par_s.items() if k.startswith(fam_call)}
        par_call, par_rel, par_runner = _call(within, margin)
    else:
        par_call, par_rel, par_runner = ("not_asked" if not askable
                                         else "no_call"), 0.0, ""

    expected = par_call if par_call.startswith(("ITPR", "RYR")) else ""
    fam_for_verdict = fam_call if fam_call in ("ITPR", "RYR") else \
        (row.get("call") or "")
    sym_v, sym_claim = R.name_verdict(row.get("gene", ""), expected,
                                      fam_for_verdict)
    nm_v, nm_claim = R.name_verdict(row.get("protein_name", ""), expected,
                                    fam_for_verdict)
    reach = callable_by_panel if callable_by_panel is not None \
        else L.panel_paralogs()
    sym_v = _panel_guard(sym_v, sym_claim, fam_for_verdict, reach)
    nm_v = _panel_guard(nm_v, nm_claim, fam_for_verdict, reach)
    pf = set((row.get("pfams") or "").split(";")) - {""}
    return {
        "accession": row["accession"], "species": row.get("species", ""),
        "group": row.get("group", ""), "tax_class": row.get("tax_class", ""),
        "length": row.get("length", ""), "reviewed": row.get("reviewed", ""),
        "census_call": row.get("call", ""),
        "census_confidence": row.get("confidence", ""),
        "gene": row.get("gene", ""), "protein_name": row.get("protein_name", ""),
        "seq_family": fam_call, "family_rel_margin": fam_rel,
        "family_runner_up": fam_runner,
        "itpr_bits": round(fam_s.get("ITPR", 0.0), 1),
        "ryr_bits": round(fam_s.get("RYR", 0.0), 1),
        "paralog_askable": int(askable),
        "seq_paralog": par_call, "paralog_rel_margin": par_rel,
        "paralog_runner_up": par_runner,
        "symbol_verdict": sym_v, "symbol_claim": sym_claim,
        "name_verdict": nm_v, "name_claim": nm_claim,
        "n_bait_hits": len(hits),
        "panel_reaches_claim": int(
            not sym_claim or sym_claim in reach.get(fam_for_verdict, set())
            or not sym_claim.startswith(("ITPR", "RYR"))),
        "has_naming_pfam": int(L.NAMING_PFAM in pf),
        "n_family_pfam": sum(1 for p in L.FAMILY_PFAMS if p in pf),
        "arch_call": row.get("arch_call", ""),
        "profile_call": row.get("profile_call", ""),
        "instruments": row.get("instruments", ""),
    }


# --------------------------------------------------------------------------
# the Pfam-recall intersection
# --------------------------------------------------------------------------
def pfam_recall(rows: list[dict]) -> list[dict]:
    """Would a Pfam query have found the records this project calls family?

    S2 enumerated the seeded space *from* the family signatures, so every row
    it holds carries at least one by construction — the recall question is
    about the **naming** signature, PF08709, the one a reader looking for
    "the IP3 receptor family" would query. S2 measured 2,911 of 15,417 seeded
    proteins without it; this asks the same of the records the project's own
    two instruments call family, and per instrument, because a record found
    only by the profile sweep never had to carry a Pfam at all.
    """
    out = []
    keys = [("call", "census call"), ("group", "taxonomic group"),
            ("instruments", "which instruments called it")]
    for key, label in keys:
        buckets: dict[str, list[dict]] = {}
        for r in rows:
            buckets.setdefault(r.get(key) or "(none)", []).append(r)
        for name, rs in sorted(buckets.items()):
            n = len(rs)
            with_naming = sum(1 for r in rs
                              if L.NAMING_PFAM in (r.get("pfams") or ""))
            complete = sum(1 for r in rs
                           if all(p in (r.get("pfams") or "")
                                  for p in L.FAMILY_PFAMS))
            out.append({
                "axis": label, "bucket": name, "n_records": n,
                "n_with_naming_pfam": with_naming,
                "frac_with_naming_pfam": L.frac(with_naming, n),
                "n_complete_architecture": complete,
                "frac_complete_architecture": L.frac(complete, n),
                "naming_pfam": L.NAMING_PFAM,
            })
    return out


AUDIT_COLS = [
    "accession", "species", "group", "tax_class", "length", "reviewed",
    "census_call", "census_confidence", "gene", "protein_name",
    "seq_family", "family_rel_margin", "family_runner_up", "itpr_bits",
    "ryr_bits", "paralog_askable", "seq_paralog", "paralog_rel_margin",
    "paralog_runner_up", "symbol_verdict", "symbol_claim", "name_verdict",
    "name_claim", "n_bait_hits", "panel_reaches_claim", "has_naming_pfam",
    "n_family_pfam",
    "arch_call", "profile_call", "instruments",
]


# --------------------------------------------------------------------------
# getting the sequences in front of the panel
# --------------------------------------------------------------------------
def resolve_and_blast(rows: list[dict], work: Path, threads: int = 8
                      ) -> tuple[Path, dict[str, str], list[str]]:
    """Sequences for the in-scope records, then blastp against the panel.

    Sequences come through `s6_lib.SequenceStore` unchanged, so S18 resolves an
    accession exactly the way S6 did — the four small archives first, the
    38 GB of proteome DBs streamed only for what they miss, and everything
    found appended to the shared cache. A record whose sequence no store holds
    is returned in `missing` and audited on its name alone rather than dropped,
    because a record the databases serve but this project cannot re-read is
    itself worth counting.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from s6_lib import SequenceStore

    work.mkdir(parents=True, exist_ok=True)
    store = SequenceStore()
    wanted = {r["accession"]: (r.get("tax_group") or r.get("group") or "")
              for r in rows}
    seqs, _source, missing = store.resolve(wanted, verbose=True)

    query = work / "protein_audit_queries.faa"
    if not query.exists() or query.stat().st_size == 0:
        L.write_fasta(query, seqs)
    db = L.make_panel_db(L.BAITS, work / "panel_db")
    out = L.blastp_vs_panel(query, db, work / "protein_audit_hits.tsv",
                            threads=threads)
    return out, seqs, missing


def parse_hits(path: Path) -> dict[str, list[tuple[str, float]]]:
    """qseqid -> [(sseqid, best bitscore)], the best HSP per query/subject."""
    best: dict[str, dict[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            q, s, bits = f[0], f[1], float(f[4])
            d = best.setdefault(q, {})
            d[s] = max(d.get(s, 0.0), bits)
    return {q: sorted(d.items(), key=lambda kv: -kv[1])
            for q, d in best.items()}
