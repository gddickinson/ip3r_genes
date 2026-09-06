"""S6 step 1 — census v6 -> the MSA representative set.

Applies `s6_rep_spec`'s eight rules and writes the audit. Every chosen row
records the rule that chose it, the cell it filled, which component of the
quality key separated it from the runner-up, and the runner-up's accession —
so a reader can ask "why this record and not that one" of any tip in the
tree S7 builds, and get an answer out of a committed table rather than out
of this docstring.

Length targets are **measured, not typed**: the target a record's length-fit
component is scored against is its own group's median over the census's
complete-architecture, non-fragment ITPR records. S23 measured the family's
genomic span varying ~100x outside the vertebrates against 6.5x inside, so a
single target would score the plant and ciliate grades against a vertebrate
ruler and systematically prefer their longest (chimeric) models.

Run:  /opt/anaconda3/envs/piezo1/bin/python scripts/s6_select_reps.py
Out:  results/msa_v2/representatives.{tsv,fasta}
      results/msa_v2/selection_stats.json
      results/msa_v2/unfilled_slots.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import s6_rep_spec as spec                                  # noqa: E402
from s6_lib import (MSA_DIR, SequenceStore, as_float, as_int,  # noqa: E402
                    binomial, load_census, make_label, pick_diverse,
                    read_tsv, write_fasta, write_live, write_tsv)

FIELDS = ["label", "accession", "rule", "cell", "group", "paralog",
          "paralog_raw", "paralog_source", "family_call", "species", "taxon_id", "band",
          "census_group", "kingdom", "phylum", "class", "order", "length",
          "n_itpr_arch", "fragment", "reviewed", "protein_existence",
          "profile_confidence", "itpr_score", "ryr_score", "rel_margin",
          "source", "db_status", "genome_accession", "copy_index",
          "copy_number", "shape_note", "decided_by", "runner_up",
          "n_candidates", "seq_store", "seq_length", "x_count"]


#: R2 — the teleost 3R co-ortholog naming convention. `paralog_of` has
#: already reduced the symbol to `itpr[123][ab]?`, so this only asks
#: whether the suffix is there.
TELEOST_3R_SYMBOL = re.compile(r"^itpr[123][ab]$")


def bait_ids() -> set[str]:
    """R7 — accessions already carrying a role in this project's searches."""
    out: set[str] = set()
    for rel, col in [("results/s5_baits/bait_manifest.tsv", "accession"),
                     ("results/s23_baits/bait_manifest.tsv", "accession"),
                     ("results/hmm_sweep/seed_manifest.tsv", "accession")]:
        p = Path(rel)
        if p.exists():
            out |= {r[col] for r in read_tsv(p) if r.get(col)}
    return out


def length_targets(itpr: list[dict]) -> dict[str, float]:
    """Per-group median length over the census's cleanest ITPR records."""
    by = defaultdict(list)
    for r in itpr:
        if (r.get("fragment") or "").strip():
            continue
        if str(r.get("n_itpr_arch")).strip() != "5":
            continue
        L = as_int(r.get("length"))
        if L:
            by[r["grp"]].append(L)
    allL = [L for v in by.values() for L in v]
    default = statistics.median(allL) if allL else 2700.0
    return {g: statistics.median(v) if len(v) >= 5 else default
            for g, v in by.items()} | {"__default__": default}


class Picker:
    """Applies a rule to a candidate pool and records the audit."""

    def __init__(self, targets: dict[str, float], baits: set[str],
                 informative: set[str] | None = None):
        self.targets, self.baits = targets, baits
        self.informative = informative or set()
        self.chosen: dict[str, dict] = {}      # accession -> audit row
        self.unfilled: list[dict] = []

    def paralog(self, r: dict) -> tuple[str, str, str]:
        """(effective paralog, what the record said, where it said it).

        A bait attribution in a band where every locus wins the same bait is
        not a paralog label (see `s6_rep_spec.informative_attribution`), so
        the effective label is blank there and the tip goes into the tree
        unassigned. What the record said is kept, never deleted.
        """
        raw, src = paralog_of_row(r)
        if src == "s5_cell" and band_of_row(r) not in self.informative:
            return "", raw, src
        return raw, raw, src

    def key(self, row: dict):
        t = self.targets.get(row["grp"], self.targets["__default__"])
        return spec.rank_key(row, t, self.baits)

    def take(self, pool: list[dict], n: int, *, rule: str, cell: str,
             diverse_key=None, allow_short: bool = False,
             new_species: bool = False, unique_key=None) -> list[dict]:
        """Fill up to `n` slots of one cell, recording the audit.

        `new_species` skips records of a species already represented — the
        expansions (R6) are picked first, so without it the phylum rules
        would spend the Porifera slot on a second *Amphimedon* tip.
        """
        taken_species = {binomial(c["species"]) for c in self.chosen.values()}
        cands = []
        for r in pool:
            if r["accession"] in self.chosen:
                continue
            if new_species and binomial(r.get("species", "")) in taken_species:
                continue
            ok, note = spec.passes_shape(
                r["call"], as_int(r.get("length")), allow_short)
            if ok:
                cands.append((r, note))
        if not cands and new_species:            # fall back rather than skip
            return self.take(pool, n, rule=rule, cell=cell,
                             diverse_key=diverse_key, allow_short=allow_short,
                             unique_key=unique_key)
        if not cands:
            self.unfilled.append({"rule": rule, "cell": cell, "wanted": n,
                                  "filled": 0, "pool": len(pool),
                                  "why": "no record passes the length gate"
                                         if pool else "no census record"})
            return []
        rows = [r for r, _ in cands]
        notes = {id(r): note for r, note in cands}
        if new_species:      # and unique within this call, not only against
            best: dict[str, dict] = {}          # what earlier rules took
            for r in sorted(rows, key=self.key, reverse=True):
                best.setdefault(binomial(r.get("species", "")), r)
            rows = list(best.values())
        if unique_key:       # collapse a locus to its best record, so two
            best = {}        # isoforms of one gene cannot become two tips
            for r in sorted(rows, key=self.key, reverse=True):
                best.setdefault(unique_key(r), r)
            rows = list(best.values())
        if diverse_key and n > 1:
            picks = pick_diverse(rows, n, diverse_key, self.key)
        else:
            picks = sorted(rows, key=self.key, reverse=True)[:n]
        ranked = sorted(rows, key=self.key, reverse=True)
        # by accession, not by `list.index`: census rows are dicts and
        # `.index` compares by value, so two rows that happen to agree on
        # every field would report each other's runner-up.
        pos = {r["accession"]: i for i, r in enumerate(ranked)}
        out = []
        for r in picks:
            i = pos[r["accession"]]
            runner = ranked[i + 1] if i + 1 < len(ranked) else None
            self.chosen[r["accession"]] = self._audit(
                r, rule=rule, cell=cell, note=notes[id(r)],
                decided=spec.decisive_component(
                    self.key(r), self.key(runner) if runner else None),
                runner=runner["accession"] if runner else "",
                n_cands=len(rows))
            out.append(r)
        if len(out) < n:
            self.unfilled.append({"rule": rule, "cell": cell, "wanted": n,
                                  "filled": len(out), "pool": len(pool),
                                  "why": "pool exhausted"})
        return out

    def _audit(self, r, *, rule, cell, note, decided, runner, n_cands):
        para, raw, psrc = self.paralog(r)
        vert = r["grp"] == "Vertebrata"
        grp = spec.group_of(r["call"], para, r["grp"], vert)
        return {
            "label": make_label(grp, r.get("species", ""), r["accession"]),
            "accession": r["accession"], "rule": rule, "cell": cell,
            "group": grp, "paralog": para, "paralog_raw": raw,
            "paralog_source": psrc if para else
            (f"{psrc} (not a label in this band)" if raw else "none"),
            "family_call": r["call"], "species": r.get("species", ""),
            "taxon_id": r.get("taxon_id", ""), "band": band_of_row(r),
            "census_group": r["grp"],
            "kingdom": r.get("tax_kingdom") or r.get("kingdom", ""),
            "phylum": r.get("tax_phylum") or r.get("phylum", ""),
            "class": r.get("tax_class") or r.get("class", ""),
            "order": r.get("tax_order") or r.get("order", ""),
            "length": as_int(r.get("length")),
            "n_itpr_arch": r.get("n_itpr_arch", ""),
            "fragment": r.get("fragment", ""), "reviewed": r.get("reviewed", ""),
            "protein_existence": r.get("protein_existence", ""),
            "profile_confidence": r.get("profile_confidence", ""),
            "itpr_score": r.get("itpr_score", ""),
            "ryr_score": r.get("ryr_score", ""),
            "rel_margin": r.get("rel_margin", ""),
            "source": r.get("source", ""), "db_status": r.get("db_status", ""),
            "genome_accession": r.get("genome_accession", ""),
            "copy_index": r.get("cell", ""), "copy_number": "",
            "shape_note": note, "decided_by": decided, "runner_up": runner,
            "n_candidates": n_cands,
        }

    def force(self, row: dict, *, rule: str, cell: str) -> None:
        if row["accession"] in self.chosen:
            self.chosen[row["accession"]]["rule"] += f"+{rule}"
            return
        self.chosen[row["accession"]] = self._audit(
            row, rule=rule, cell=cell, note="", decided="forced",
            runner="", n_cands=1)


def paralog_of_row(r: dict) -> tuple[str, str]:
    """The paralog a record's *label* names, and where the label came from."""
    p = spec.paralog_of(r.get("gene", ""), "")
    if p:
        return p, "gene_symbol"
    if r.get("cell", "") in spec.PARALOGS + ("RYR",):
        return r["cell"], "s5_cell"
    p = spec.paralog_of("", r.get("protein_name", ""))
    if p:
        return p, "protein_name"
    return "", "none"


def band_of_row(r: dict) -> str:
    return spec.band_of({"class": r.get("tax_class") or r.get("class", ""),
                         "order": r.get("tax_order") or r.get("order", "")})


def informative_bands(vert: list[dict]) -> set[str]:
    """Bands where the sweeps' bait attribution discriminates (see spec)."""
    cells = defaultdict(set)
    for r in vert:
        if r.get("cell"):
            cells[band_of_row(r)].add(r["cell"])
    return {b for b, c in cells.items() if spec.informative_attribution(c)}


def locus_key(r: dict) -> str:
    """One gene = one tip. A gene symbol names the locus where there is one;
    a record with none (several hagfish UniProt entries) stands for itself,
    which over-counts isoforms rather than merging two real loci."""
    g = (r.get("gene") or "").strip().lower()
    return g or r["accession"]


def genus(r: dict) -> str:
    return (r.get("species") or "").split()[0] if r.get("species") else "?"


# --------------------------------------------------------------- the rules

def select(rows: list[dict]) -> Picker:
    itpr = [r for r in rows if r["call"] == "ITPR"]
    ryr = [r for r in rows if r["call"] == "RYR"]
    by_acc = {r["accession"]: r for r in rows}
    p = Picker(length_targets(itpr), bait_ids(),
               informative_bands([r for r in itpr
                                  if r["grp"] == "Vertebrata"]))

    # R1 -- the vertebrate paralog grid ------------------------------------
    for acc, label in spec.FORCED_ITPR.items():
        if acc in by_acc:
            p.force(by_acc[acc], rule="R1", cell=f"forced:{label}")
    vert = [r for r in itpr if r["grp"] == "Vertebrata"]
    labelled = defaultdict(list)
    for r in vert:
        para, _, _ = p.paralog(r)
        b = band_of_row(r)
        if b and para in spec.PARALOGS:
            labelled[(b, para)].append(r)
    for b in spec.BANDS:
        for para in spec.PARALOGS:
            want = 2 if b in spec.DEEP_BANDS else 1
            p.take(labelled.get((b, para), []), want, rule="R1",
                   cell=f"{para}@{b}", diverse_key=genus)

    # R2 -- the teleost 3R co-orthologs ------------------------------------
    pairs = defaultdict(lambda: defaultdict(list))
    for r in vert:
        para, src = paralog_of_row(r)
        if src != "gene_symbol" or para not in spec.PARALOGS:
            continue
        sym = (r.get("gene") or "").lower()
        if TELEOST_3R_SYMBOL.match(sym):               # itpr1a / itpr1b
            pairs[(binomial(r["species"]), para)][sym].append(r)
    ready = [(sp_para, d) for sp_para, d in pairs.items() if len(d) == 2]
    ready.sort(key=lambda kv: -max(
        as_float(r.get("itpr_score")) for v in kv[1].values() for r in v))
    for (sp, para), d in ready[:spec.TELEOST_3R_SPECIES]:
        for sym, rs in sorted(d.items()):
            p.take(rs, 1, rule="R2", cell=f"3R:{sym}@{sp.split()[0]}")

    # R3 -- the unlabelled deep vertebrate grade, one tip per locus --------
    # Two stages, and the order is the point: pick the *species* first, then
    # take all of its loci. Spreading over species first would give two
    # loci each of three cyclostomes, and the 2R question needs a complete
    # copy set — three lamprey genes against ITPR1/2/3 is a test, two is a
    # sample.
    labelled_accs = {r["accession"] for v in labelled.values() for r in v}
    for b in spec.UNLABELLED_BANDS:
        pool = [r for r in vert
                if band_of_row(r) == b
                and r["accession"] not in labelled_accs
                and spec.passes_shape(r["call"], as_int(r.get("length")),
                                      True)[0]]
        loci = defaultdict(set)
        for r in pool:
            loci[binomial(r.get("species", ""))].add(locus_key(r))
        species = sorted(loci, key=lambda s: (-len(loci[s]), s)
                         )[:spec.UNLABELLED_SPECIES_PER_BAND]
        for sp in species:
            p.take([r for r in pool if binomial(r.get("species", "")) == sp],
                   spec.UNLABELLED_LOCI_PER_SPECIES, rule="R3",
                   cell=f"basal@{b}:{sp.replace(' ', '_')}",
                   diverse_key=locus_key, allow_short=True,
                   unique_key=locus_key)
        if not species:
            p.unfilled.append({"rule": "R3", "cell": f"basal@{b}",
                               "wanted": spec.UNLABELLED_LOCI_PER_SPECIES,
                               "filled": 0, "pool": len(pool),
                               "why": "no unlabelled record passes the "
                                      "length gate"})

    # R6 -- the copy-number expansions -------------------------------------
    # Runs before the phylum rules: the expansions are the scarcer evidence,
    # and one per phylum so four slots are not spent on two sponges.
    taken, seen_phyla = 0, set()
    for acc, (sp, n, grp, ph) in copy_number_genomes():
        if taken >= spec.COPY_GENOMES:
            break
        if ph in seen_phyla:
            continue
        pool = [r for r in itpr if r.get("genome_accession") == acc]
        if len(pool) < 2:
            continue
        got = p.take(pool, spec.COPY_PER_GENOME, rule="R6",
                     cell=f"copies:{acc}", allow_short=True)
        for r in got:
            p.chosen[r["accession"]]["copy_number"] = str(n)
        if got:
            taken += 1
            seen_phyla.add(ph)

    # R4 -- non-vertebrate metazoa, per phylum -----------------------------
    met = [r for r in itpr if r["grp"] == "Metazoa (non-vertebrate)"]
    by_phylum = defaultdict(list)
    for r in met:
        by_phylum[(r.get("tax_phylum") or r.get("phylum") or "?")].append(r)
    for ph in sorted(by_phylum, key=lambda k: -len(by_phylum[k])):
        want = spec.METAZOAN_DEEP_PHYLA.get(ph, spec.PER_PHYLUM)
        p.take(by_phylum[ph], want, rule="R4", cell=f"phylum:{ph}",
               diverse_key=lambda r: (r.get("tax_class") or "?", genus(r)),
               allow_short=True, new_species=True)

    # R5 -- non-metazoan eukaryotes, S20 verdict-gated ---------------------
    for g in spec.EUKARYOTE_GROUPS:
        pool = [r for r in itpr if r["grp"] == g
                and spec.plant_fungal_ok(r.get("plant_fungal_verdict", ""))]
        want = spec.EUKARYOTE_DEEP_GROUPS.get(g, spec.PER_PHYLUM)
        p.take(pool, want, rule="R5", cell=f"group:{g}",
               diverse_key=lambda r: (r.get("tax_phylum") or "?", genus(r)),
               allow_short=True, new_species=True)

    # R7 -- the novel gene models no database annotates ---------------------
    novel = [r for r in itpr
             if r.get("db_status") == spec.NOVEL_DB_STATUS
             and r["accession"] not in p.chosen]
    p.take([r for r in novel if r["grp"] == "Vertebrata"],
           spec.NOVEL_VERTEBRATE, rule="R7", cell="novel:vertebrate",
           diverse_key=lambda r: (band_of_row(r), r.get("cell") or "?"))
    p.take([r for r in novel if r["grp"] != "Vertebrata"],
           spec.NOVEL_NONVERTEBRATE, rule="R7", cell="novel:nonvertebrate",
           diverse_key=lambda r: (r["grp"], genus(r)), allow_short=True,
           new_species=True)

    # R8 -- the ryanodine-receptor outgroup --------------------------------
    for acc, label in spec.FORCED_RYR.items():
        if acc in by_acc:
            p.force(by_acc[acc], rule="R8", cell=f"outgroup:{label}")
    nonvert_ryr = [r for r in ryr if r["grp"] != "Vertebrata"]
    p.take(nonvert_ryr, spec.RYR_NONVERT, rule="R8", cell="outgroup:nonvert",
           diverse_key=lambda r: (r.get("tax_phylum") or "?", genus(r)))
    return p


def unfilled_context(rows: list[dict], unfilled: list[dict],
                     picker: "Picker") -> dict:
    """For every empty cell, what the census actually holds in that band.

    The report has to be able to say *why* a cell is empty without a human
    typing the reason, and three statements have to be told apart: no
    record at all, no record carrying an annotation label, and no record
    whose label is a paralog assignment rather than a bait attribution. In
    the cyclostome band it is the third — six loci, every one attributed to
    ITPR1 because the ITPR1 bait wins them all.
    """
    out = {}
    itpr = [r for r in rows if r["call"] == "ITPR" and r["grp"] == "Vertebrata"]
    for u in unfilled:
        band = u["cell"].split("@", 1)[-1].split(":")[0]
        pool = [r for r in itpr if band_of_row(r) == band]
        raw = [r for r in pool if paralog_of_row(r)[0]]
        eff = [r for r in pool if picker.paralog(r)[0]]
        out[u["cell"]] = {
            "band": band, "itpr_records_in_band": len(pool),
            "with_any_label": len(raw),
            "with_a_usable_label": len(eff),
            "labels_present": sorted({paralog_of_row(r)[0] for r in raw}),
            "label_sources": sorted({paralog_of_row(r)[1] for r in raw}),
            "loci_in_band": len({locus_key(r) for r in pool}),
        }
    return out


def copy_number_genomes() -> list[tuple[str, tuple[str, int, str, str]]]:
    """R6 — S23's measured copy numbers above the vertebrate paralog count."""
    p = Path("results/s23_scope/copy_number_ledger.tsv")
    if not p.exists():
        return []
    out = []
    for r in read_tsv(p):
        n = as_int(r.get("n_full"))
        if n > spec.COPY_NUMBER_THRESHOLD:
            out.append((r["accession"], (r.get("organism") or "", n,
                                         r.get("group") or "",
                                         r.get("phylum") or r.get("group") or "?")))
    out.sort(key=lambda kv: -kv[1][1])
    return out


# ------------------------------------------------------------------ driver

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="select and report, do not resolve sequences")
    args = ap.parse_args()

    write_live("S6", [("select representatives", False),
                      ("resolve sequences", False), ("MAFFT L-INS-i", False),
                      ("trimAl", False), ("matrices + figures", False)])
    rows = load_census()
    picker = select(rows)
    chosen = list(picker.chosen.values())
    print(f"selected {len(chosen)} representatives from {len(rows)} census rows")
    for rule, n in sorted(Counter(
            c["rule"].split("+")[0] for c in chosen).items()):
        print(f"  {rule}: {n}")
    if picker.unfilled:
        print(f"  {len(picker.unfilled)} unfilled slots "
              f"(see unfilled_slots.tsv)")
    if args.dry_run:
        for c in sorted(chosen, key=lambda c: (c["group"], c["species"])):
            print(f"    {c['group']:16s} {c['species'][:38]:40s} "
                  f"{c['length']:>5} {c['rule']:4s} {c['cell']}")
        return 0

    # ---- sequences -------------------------------------------------------
    write_live("S6", [("select representatives", True),
                      ("resolve sequences", False), ("MAFFT L-INS-i", False),
                      ("trimAl", False), ("matrices + figures", False)])
    store = SequenceStore()
    wanted = {c["accession"]: c["census_group"] for c in chosen}
    seqs, source, missing = store.resolve(wanted)
    if missing:
        print(f"  WARNING: {len(missing)} representatives have no sequence:")
        for a in missing[:20]:
            print(f"    {a}  ({picker.chosen[a]['species']})")
        chosen = [c for c in chosen if c["accession"] not in set(missing)]

    out_seqs = {}
    for c in chosen:
        s = seqs[c["accession"]].upper().replace("*", "X")
        c["seq_store"] = source[c["accession"]]
        c["seq_length"] = len(s)
        c["x_count"] = s.count("X")
        out_seqs[c["label"]] = s

    dup = [l for l, n in Counter(c["label"] for c in chosen).items() if n > 1]
    if dup:
        raise SystemExit(f"duplicate labels, would collide in the MSA: {dup}")

    chosen.sort(key=lambda c: (spec.GROUP_ORDER.index(c["group"]),
                               c["species"]))
    MSA_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(MSA_DIR / "representatives.tsv", FIELDS, chosen)
    write_fasta({c["label"]: out_seqs[c["label"]] for c in chosen},
                MSA_DIR / "representatives.fasta")
    write_tsv(MSA_DIR / "unfilled_slots.tsv",
              ["rule", "cell", "wanted", "filled", "pool", "why"],
              picker.unfilled)
    stats = {
        "census_rows": len(rows), "selected": len(chosen),
        "unfilled_context": unfilled_context(rows, picker.unfilled,
                                             picker),
        "by_rule": dict(Counter(c["rule"].split("+")[0] for c in chosen)),
        "by_group": dict(Counter(c["group"] for c in chosen)),
        "by_seq_store": dict(Counter(c["seq_store"] for c in chosen)),
        "by_decided_by": dict(Counter(c["decided_by"] for c in chosen)),
        "length_targets": {k: round(v, 1)
                           for k, v in picker.targets.items()},
        "short_exceptions": [c["label"] for c in chosen
                             if c["shape_note"] == "short_exception"],
        "unfilled_slots": len(picker.unfilled),
        "missing_sequence": missing,
        "total_residues": sum(c["seq_length"] for c in chosen),
        "longest": max(c["seq_length"] for c in chosen),
    }
    (MSA_DIR / "selection_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"wrote {MSA_DIR/'representatives.fasta'} — {len(chosen)} seqs, "
          f"{stats['total_residues']:,} residues, longest {stats['longest']}")
    for g in spec.GROUP_ORDER:
        n = stats["by_group"].get(g, 0)
        if n:
            print(f"  {g:18s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
