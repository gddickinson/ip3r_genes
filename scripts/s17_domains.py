"""S17 stage 1 — the architecture, in each reference's own numbering.

Every claim S17 makes is of the form "this *part* of the receptor is more
constrained than that part", so the parts have to be located in the coordinate
system each table is reported in. This family gives that in two pieces which
are not equally trustworthy, and the module keeps them apart:

* **The Pfam elements are measured, per accession, and need no transfer.**
  `results/s0_baseline/review_figures/domain_coords.tsv` holds InterPro's
  coordinates for Q14643, Q14571 and Q14573 *separately*. The PIEZO project
  this method is ported from had one published boundary list in a mouse
  protein's numbering and had to carry it onto every reference; here that whole
  class of error does not arise, and the module refuses to invent it.

* **The structural elements are measured on one protein and are transferred.**
  S0 measured the channel on PDB 6DQN — human ITPR3, so Q14573 numbering — and
  what it measured is the selectivity filter, the gate and the ten residues
  that contact IP3. ITPR1 and ITPR2 get those by pairwise alignment, and the
  transfer is a **positive test, not an assumption** (D54): the filter must
  arrive on the same GGGVGD motif, the gate and the contacts on the same amino
  acids. A transferred element whose residues disagree is written with its
  disagreement in the row rather than used.

One thing falls straight out of joining the two, and it changes how the rest of
S17 has to be read: **none of the ten measured IP3 contacts lies in PF08709**,
the signature Pfam calls "Inositol 1,4,5-trisphosphate/ryanodine receptor".
They sit in MIR and RIH — the two domains the family shares with the ryanodine
receptor. So `nterm_trefoil` is not the ligand site, and S17's ligand-site
question is asked of the *measured* contacts, never of the Pfam label.

Two products:

* `domain_map.tsv` — start/end of every element in every reference, with how it
  was obtained and whether its transfer check passed.
* a per-residue **primary** assignment, most specific element wins, because the
  elements nest (the filter and the gate are inside the channel) and a per-site
  table needs each residue to belong to exactly one thing. Sequence between two
  named elements is a **named linker**, never left blank.

    python scripts/s17_domains.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402

#: Pfam signature -> S17's element name. Two RIH copies are told apart by
#: order of appearance, because PF01365 hits the protein twice and calling both
#: "RIH" would pool an N-terminal and a central domain into one test.
PFAM_ELEMENTS = {
    "PF08709": "nterm_trefoil",
    "PF02815": "MIR",
    "PF01365": "RIH",
    "PF08454": "RIH_assoc",
    "PF00520": "channel",
}

#: Most specific first. A residue takes the first element it falls in, so the
#: filter is the filter rather than "channel".
PRIMARY_ORDER = ("selectivity_filter", "gate", "luminal_loop", "nterm_trefoil",
                 "MIR", "RIH_N", "RIH_C", "RIH_assoc", "channel")

#: Elements that carry the mechanism, in the sense S17's question uses: the
#: ligand end and the conducting end. `linker` and the termini are the
#: within-protein controls.
MECHANISM_ELEMENTS = ("nterm_trefoil", "MIR", "RIH_N", "RIH_C", "RIH_assoc",
                      "channel", "selectivity_filter", "gate")


# --------------------------------------------------------- measured elements

def pfam_spans(acc: str) -> list[tuple[str, int, int, str]]:
    """(element, start, end, pfam) for one accession, from S0's committed table."""
    rows = [r for r in L.read_tsv(L.S0_FIG_DIR / "domain_coords.tsv")
            if r["accession"] == acc]
    if not rows:
        raise SystemExit(f"no domain_coords rows for {acc}")
    out, seen_rih = [], 0
    for r in sorted(rows, key=lambda x: int(x["start"])):
        name = PFAM_ELEMENTS.get(r["pfam"])
        if name is None:
            continue
        if name == "RIH":
            seen_rih += 1
            name = "RIH_N" if seen_rih == 1 else "RIH_C"
        out.append((name, int(r["start"]), int(r["end"]), r["pfam"]))
    return out


def luminal_span(m: dict, pfam: list[tuple[str, int, int, str]]
                 ) -> tuple[int, int] | None:
    """The luminal loop, located by geometry rather than drawn by eye.

    The channel element's mean conservation is bimodal, and a boundary chosen
    after seeing that profile would be the profile explaining itself. So the
    loop is defined on the **structure**: S0 committed one subunit's Ca trace
    of 6DQN and the membrane's own axial span, and the loop is the contiguous
    run of channel residues whose Ca sits beyond the membrane on the luminal
    side. Nothing in that derivation knows what any column looks like.
    """
    ca = {int(r["resseq"]): float(r["z"])
          for r in L.read_tsv(L.S0_FIG_DIR / "structure_ca.tsv")}
    lo = m["tm_span_z_A"][0]
    ch = [(s, e) for n, s, e, _p in pfam if n == "channel"]
    if not ch:
        return None
    s0, e0 = ch[0]
    # A residue the map does not resolve is **missing data, not a return to
    # the membrane** — a disordered luminal loop is precisely what fails to
    # resolve, and breaking the run on it would return a three-residue stub.
    # The run therefore closes only on a residue that is resolved and is not
    # luminal, and how many unresolved residues the span contains is reported
    # beside it, because that count is itself evidence about the loop.
    runs, cur = [], None
    for r in range(s0, e0 + 1):
        if r in ca and ca[r] >= lo:
            if cur:
                runs.append(tuple(cur))
                cur = None
            continue
        if r in ca:
            cur = [r, r] if cur is None else [cur[0], r]
        elif cur:
            cur[1] = max(cur[1], r)          # carry the gap, do not close on it
    if cur:
        runs.append(tuple(cur))
    if not runs:
        return None
    best = max(runs, key=lambda x: x[1] - x[0])
    # trim a trailing gap: the span must end on a residue the structure places
    end = max(r for r in range(best[0], best[1] + 1) if r in ca)
    return (best[0], end)


def structure_elements() -> dict:
    """The channel S0 measured on 6DQN, in Q14573 numbering."""
    m = json.loads((L.S0_FIG_DIR / "structure_meta.json").read_text())
    motif = m["filter_motif"]
    fstart = int(m["filter_motif_start"])
    lining = [int(x[3:]) for x in m["filter_lining_residues"]]
    gate = sorted(int(x[3:]) for x in m["gate_lining_residues"])
    return {
        "pdb": m["pdb_id"], "uniprot": m["uniprot"],
        "filter_motif": motif, "filter_motif_start": fstart,
        "filter_span": (min(lining + [fstart]), fstart + len(motif) - 1),
        "gate_span": (gate[0], gate[-1]),
        "filter_lining": lining,
        "gate_lining": gate,
        "ip3_contacts": sorted(int(x[3:])
                               for x in m["ip3_contacts"]["same_subunit"]),
        "resolution_A": m["resolution_A"],
        "membrane_z": tuple(m["tm_span_z_A"]),
        "luminal_span": luminal_span(
            m, pfam_spans(m["uniprot"])),
    }


# ------------------------------------------------------------------- transfer

def _check(seq_ref: str, seq_tgt: str, ref_pos: list[int],
           xfer: dict[int, int]) -> tuple[int, int, str]:
    """How many transferred positions keep their amino acid, and which do not."""
    ok = bad = 0
    notes = []
    for p in ref_pos:
        q = xfer.get(p)
        if q is None:
            bad += 1
            notes.append(f"{seq_ref[p - 1]}{p}->unaligned")
            continue
        if seq_tgt[q - 1] == seq_ref[p - 1]:
            ok += 1
        else:
            bad += 1
            notes.append(f"{seq_ref[p - 1]}{p}->{seq_tgt[q - 1]}{q}")
    return ok, bad, "; ".join(notes)


def build(out_dir: Path = L.OUT_DIR) -> list[dict]:
    st = structure_elements()
    ref_acc = L.REFERENCES[L.STRUCTURE_REF][1]
    if st["uniprot"] != ref_acc:
        raise SystemExit(f"structure_meta is on {st['uniprot']}, not {ref_acc}")
    ref_seq = L.uniprot_fasta(ref_acc)

    rows: list[dict] = []
    for paralog, (_label, acc, desc) in L.REFERENCES.items():
        target = L.uniprot_fasta(acc)

        # --- Pfam: measured on this accession, no transfer, nothing to check --
        for name, s, e, pfam in pfam_spans(acc):
            rows.append({
                "paralog": paralog, "ref_acc": acc, "reference": desc,
                "element": name, "kind": "pfam", "source": f"InterPro {pfam}",
                "start": s, "end": e, "n_residues": e - s + 1,
                "transferred_from": "", "n_check_ok": "", "n_check_fail": "",
                "check": "not_applicable_measured_on_this_accession",
                "anchor": "", "disagreements": "",
            })

        # --- structural: measured on 6DQN/Q14573, transferred and checked -----
        # The check is on the *anchor* — the filter's own GGGVGD motif, and the
        # gate's two lining residues — because that is a property the pairwise
        # alignment knows nothing about, so it can fail. Per-residue identity
        # inside the span is a separate column: a paralog substitution at a
        # correctly-placed site is a result, not a transfer error, and folding
        # the two into one verdict would report the first as the second.
        xfer = ({i: i for i in range(1, len(ref_seq) + 1)} if acc == ref_acc
                else L.transfer_positions(ref_seq, target))
        structural = [
            ("selectivity_filter", st["filter_span"],
             ("motif", st["filter_motif"], st["filter_motif_start"])),
            ("gate", st["gate_span"], ("residues", st["gate_lining"], None)),
        ]
        if st["luminal_span"]:
            # The anchor here is containment, not a motif: a loop has no
            # conserved string to land on, but a transfer that slipped would
            # put it outside the target's own channel domain, which is a
            # property measured independently on that accession.
            structural.append(("luminal_loop", st["luminal_span"],
                               ("inside_channel", None, None)))
        for name, span, anchor in structural:
            lo, hi = xfer.get(span[0]), xfer.get(span[1])
            ok, bad, notes = _check(ref_seq, target,
                                    list(range(span[0], span[1] + 1)), xfer)
            if acc == ref_acc:
                verdict, anchor_note = "measured_here", ""
            elif lo is None or hi is None:
                verdict, anchor_note = "FAILED", "span could not be placed"
            elif anchor[0] == "inside_channel":
                ch = [(s, e) for n, s, e, _p in pfam_spans(acc) if n == "channel"]
                if ch and ch[0][0] <= lo and hi <= ch[0][1]:
                    verdict = "anchor_confirmed"
                    anchor_note = f"inside channel {ch[0][0]}-{ch[0][1]}"
                else:
                    verdict = "FAILED"
                    anchor_note = f"{lo}-{hi} falls outside channel {ch}"
            elif anchor[0] == "motif":
                motif, mstart = anchor[1], anchor[2]
                q = xfer.get(mstart)
                got = target[q - 1:q - 1 + len(motif)] if q else ""
                if got == motif:
                    verdict = "anchor_confirmed"
                    anchor_note = f"{motif} at {q}"
                else:
                    verdict = "FAILED"
                    anchor_note = f"expected {motif}, found {got or 'nothing'}"
            else:
                aok, abad, anotes = _check(ref_seq, target, anchor[1], xfer)
                verdict = "anchor_confirmed" if abad == 0 else "FAILED"
                anchor_note = anotes or "lining residues identical"
            rows.append({
                "paralog": paralog, "ref_acc": acc, "reference": desc,
                "element": name, "kind": "structural",
                "source": f"{st['pdb']} {st['resolution_A']} A",
                "start": lo or "", "end": hi or "",
                "n_residues": (hi - lo + 1) if (lo and hi) else 0,
                "transferred_from": "" if acc == ref_acc else ref_acc,
                "n_check_ok": ok, "n_check_fail": bad,
                "check": verdict, "anchor": anchor_note,
                "disagreements": notes if bad else "",
            })
            if verdict == "FAILED":
                raise SystemExit(
                    f"[s17] {paralog} {name}: transfer anchor failed "
                    f"({anchor_note}) — refusing to write an unvalidated "
                    f"coordinate")

        # the ten measured IP3 contacts, one row, with the same check
        cok, cbad, cnotes = _check(ref_seq, target, st["ip3_contacts"], xfer)
        placed = [xfer[p] for p in st["ip3_contacts"] if p in xfer]
        rows.append({
            "paralog": paralog, "ref_acc": acc, "reference": desc,
            "element": "ip3_contact_set", "kind": "structural",
            "source": f"{st['pdb']} IP3 contacts <= 4.5 A",
            "start": min(placed) if placed else "",
            "end": max(placed) if placed else "",
            "n_residues": len(placed),
            "transferred_from": "" if acc == ref_acc else ref_acc,
            "n_check_ok": cok, "n_check_fail": cbad,
            "check": ("measured_here" if acc == ref_acc else
                      "all_identical" if cbad == 0 else "placed_with_substitutions"),
            "anchor": f"{cok}/{len(st['ip3_contacts'])} contacts identical",
            "disagreements": cnotes,
        })

    L.write_tsv(out_dir / "domain_map.tsv", rows)
    _write_sites(out_dir, st)
    return rows


def _write_sites(out_dir: Path, st: dict) -> None:
    """Every measured functional residue, in all three references' numbering.

    This is the set the ligand-site question is asked of. It is deliberately
    *not* the PF08709 span: none of the ten IP3 contacts falls inside it.
    """
    ref_acc = L.REFERENCES[L.STRUCTURE_REF][1]
    ref_seq = L.uniprot_fasta(ref_acc)
    groups = [("ip3_contact", st["ip3_contacts"]),
              ("filter_lining", st["filter_lining"]),
              ("gate_lining", st["gate_lining"])]
    rows = []
    for paralog, (_l, acc, desc) in L.REFERENCES.items():
        target = L.uniprot_fasta(acc)
        xfer = ({i: i for i in range(1, len(ref_seq) + 1)} if acc == ref_acc
                else L.transfer_positions(ref_seq, target))
        for site_class, positions in groups:
            for p in positions:
                q = xfer.get(p)
                rows.append({
                    "paralog": paralog, "ref_acc": acc, "reference": desc,
                    "site_class": site_class,
                    "source_resi": p, "source_aa": ref_seq[p - 1],
                    "resi": q or "", "aa": target[q - 1] if q else "",
                    "identical": bool(q) and target[q - 1] == ref_seq[p - 1],
                    "transferred": acc != ref_acc,
                })
    L.write_tsv(out_dir / "functional_sites.tsv", rows)


# ---------------------------------------------------------- residue assignment

def residue_element(out_dir: Path = L.OUT_DIR) -> dict[str, dict[int, str]]:
    """paralog -> {residue: primary element}, from the committed map.

    Sequence outside every named element is a **named linker** — `nterm`,
    `cterm`, or `linker_<upstream>_<downstream>` — because "unassigned" would
    pool the N-terminus, five inter-domain stretches and the whole C-terminal
    tail into one bucket and then use it as the within-protein control.
    """
    path = out_dir / "domain_map.tsv"
    if not path.exists():
        build(out_dir)
    lengths = {p: len(L.uniprot_fasta(a)) for p, (_l, a, _d) in L.REFERENCES.items()}

    spans: dict[str, dict[str, tuple[int, int]]] = {}
    for r in L.read_tsv(path):
        if r["start"] == "" or r["element"] == "ip3_contact_set":
            continue
        spans.setdefault(r["paralog"], {})[r["element"]] = (
            int(r["start"]), int(r["end"]))

    out: dict[str, dict[int, str]] = {}
    for paralog, dom in spans.items():
        assign: dict[int, str] = {}
        for name in PRIMARY_ORDER:
            if name not in dom:
                continue
            s, e = dom[name]
            for r in range(s, e + 1):
                assign.setdefault(r, name)
        # name the gaps by their flanking Pfam elements
        ordered = sorted(((s, e, n) for n, (s, e) in dom.items()
                          if n in ("nterm_trefoil", "MIR", "RIH_N", "RIH_C",
                                   "RIH_assoc", "channel")))
        bounds = [(0, 0, "nterm")] + ordered + \
                 [(lengths[paralog] + 1, lengths[paralog] + 1, "cterm")]
        for (s1, e1, n1), (s2, _e2, n2) in zip(bounds, bounds[1:]):
            lo, hi = e1 + 1, s2 - 1
            if hi < lo:
                continue
            if n1 == "nterm":
                name = "nterm"
            elif n2 == "cterm":
                name = "cterm"
            else:
                name = f"linker_{n1}_{n2}"
            for r in range(lo, hi + 1):
                assign.setdefault(r, name)
        out[paralog] = assign
    return out


def site_index(out_dir: Path = L.OUT_DIR) -> dict[str, dict[str, set[int]]]:
    """paralog -> {site_class: residues}, from `functional_sites.tsv`."""
    path = out_dir / "functional_sites.tsv"
    if not path.exists():
        build(out_dir)
    out: dict[str, dict[str, set[int]]] = {}
    for r in L.read_tsv(path):
        if not r["resi"]:
            continue
        out.setdefault(r["paralog"], {}).setdefault(
            r["site_class"], set()).add(int(r["resi"]))
    return out


def run(out_dir: Path = L.OUT_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = build(out_dir)
    for r in rows:
        if r["kind"] == "structural":
            print(f"[s17] {r['paralog']:<6} {r['element']:<18} "
                  f"{r['start']}-{r['end']}  check={r['check']}"
                  + (f"  [{r['disagreements']}]" if r["disagreements"] else ""))
    assign = residue_element(out_dir)
    for p, a in assign.items():
        n = len(L.uniprot_fasta(L.REFERENCES[p][1]))
        print(f"[s17] {p}: {len(a)}/{n} residues assigned to an element")
    return {"n_rows": len(rows)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    run(ap.parse_args().out)


if __name__ == "__main__":
    main()
