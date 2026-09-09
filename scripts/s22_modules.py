"""S22 stage 1 — the two functional modules, defined by measurement.

Neither module this task compares is a Pfam domain, and that is the reason
the stage exists.  The **ligand core** is the part of the receptor that
holds IP3; Pfam calls the same stretch MIR and RIH_N and knows nothing about
the ligand.  The **pore module** is the ion pathway; Pfam calls it PF00520
and includes in it a fifty-residue luminal loop that is the least conserved
element in the whole receptor (S17).  So both are derived here, and both are
derived twice, because a comparison between two hand-drawn regions is worth
what its boundaries are worth.

Primary definitions, each from a measurement this project made:

* `ligand_core` = the minimal contiguous span containing every IP3 contact
  measured on the structure (S0, 6DQN, <= 4.5 A).  It is the smallest region
  that can be called the binding site without importing a boundary from a
  paper.
* `pore_module` = the PF00520 channel span **minus** the luminal loop S17
  located by geometry.  The loop is solvent-facing, unresolved in half the
  depositions, and at JSD 0.49 against the protein's 0.74; leaving it in
  would make the pore look unconstrained for a reason that has nothing to do
  with the pore.

Sensitivity definitions, committed beside them and recomputed everywhere:

* `ligand_core_ibc` = the published IP3-binding core, ITPR1 224-604 [R05],
  transferred to the other two paralogues through S17's alignment.
* `pore_module_all` = PF00520 as InterPro draws it, luminal loop included.

Every definition is checked against something it does not itself contain:
the ligand core must hold all ten contacts and no filter or gate residue,
the pore module must hold both filter residues and both gate residues and no
IP3 contact.  A definition that fails raises.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s17_lib as S17  # noqa: E402

# The published IP3-binding core, in human ITPR1 numbering, from the crystal
# structure of the isolated core [R05].  Declared as an input (D15's habit):
# this project did not measure it and does not re-derive it.
IBC_REFERENCE = ("ITPR1", 224, 604, "R05")

PRIMARY = {"ligand_core": "contact_span", "pore_module": "channel_minus_luminal"}


class ModuleRefusal(RuntimeError):
    pass


def _sites(paralog: str, cls: str) -> list[int]:
    return sorted(r["resi"] for r in L.functional_sites()
                  if r["paralog"] == paralog and r["site_class"] == cls)


def _element(paralog: str, name: str) -> tuple[int, int] | None:
    for r in L.domain_map():
        if r["paralog"] == paralog and r["element"] == name:
            return (r["start"], r["end"])
    return None


def _transfer_ibc() -> dict[str, tuple[int, int]]:
    """ITPR1 224-604 carried to ITPR2 and ITPR3 through a pairwise alignment.

    S17's `transfer_positions` is used unchanged, so the transfer is the one
    every other cross-paralogue coordinate in this project went through.
    """
    src, lo, hi, _ = IBC_REFERENCE
    out = {src: (lo, hi)}
    src_seq = S17.uniprot_fasta(L.REFERENCES[src][1])
    for p in L.PARALOGS:
        if p == src:
            continue
        dst_seq = S17.uniprot_fasta(L.REFERENCES[p][1])
        m = S17.transfer_positions(src_seq, dst_seq)
        lo2 = next((m[i] for i in range(lo, hi + 1) if i in m), None)
        hi2 = next((m[i] for i in range(hi, lo - 1, -1) if i in m), None)
        if lo2 is None or hi2 is None:
            raise ModuleRefusal(f"IBC boundary did not transfer to {p}")
        out[p] = (lo2, hi2)
    return out


def check_definition(module: str, n_contact: int, n_filt: int, n_gate: int,
                     where: str) -> str:
    """A region is the module it claims to be, or it is not written.

    Split out of `build` so the refusal has a name and can be exercised on
    its own: a definition rule that has never been shown to reject anything
    is not a rule (`s5_bait_screen.self_test()`'s point applied to a
    boundary).
    """
    if module == "ligand_core":
        if n_contact == 10 and n_filt == 0 and n_gate == 0:
            return "all_contacts_no_pore"
    elif n_filt == 2 and n_gate == 2 and n_contact == 0:
        return "both_pore_sites_no_contact"
    raise ModuleRefusal(f"{where}: contacts={n_contact} filter={n_filt} "
                        f"gate={n_gate} — definition rejected")


def build() -> list[dict]:
    """One row per paralogue x module definition, with its checks."""
    ibc = _transfer_ibc()
    rows: list[dict] = []
    for p in L.PARALOGS:
        contacts = _sites(p, "ip3_contact")
        filt = _sites(p, "filter_lining")
        gate = _sites(p, "gate_lining")
        if len(contacts) != 10 or len(filt) != 2 or len(gate) != 2:
            raise ModuleRefusal(
                f"{p}: expected 10 contacts / 2 filter / 2 gate, got "
                f"{len(contacts)}/{len(filt)}/{len(gate)}")
        chan = _element(p, "channel")
        lum = _element(p, "luminal_loop")
        if chan is None or lum is None:
            raise ModuleRefusal(f"{p}: channel or luminal_loop missing")

        defs = {
            "contact_span": (min(contacts), max(contacts), "span of the ten measured IP3 contacts"),
            "ibc_literature": (ibc[p][0], ibc[p][1], f"published IP3-binding core [{IBC_REFERENCE[3]}]"),
            "channel_minus_luminal": (chan[0], chan[1], "PF00520 less the luminal loop"),
            "channel_all": (chan[0], chan[1], "PF00520 as InterPro draws it"),
        }
        excl = {"channel_minus_luminal": [lum]}

        for name, (lo, hi, note) in defs.items():
            module = ("ligand_core" if name in ("contact_span", "ibc_literature")
                      else "pore_module")
            holes = excl.get(name, [])
            resis = [i for i in range(lo, hi + 1)
                     if not any(a <= i <= b for a, b in holes)]
            rs = set(resis)
            n_contact = len(rs & set(contacts))
            n_filt = len(rs & set(filt))
            n_gate = len(rs & set(gate))
            check = check_definition(module, n_contact, n_filt, n_gate,
                                     f"{p}/{name}")
            rows.append({
                "paralog": p, "module": module, "definition": name,
                "is_primary": PRIMARY[module] == name,
                "start": lo, "end": hi,
                "excluded": ";".join(f"{a}-{b}" for a, b in holes),
                "n_residues": len(resis),
                "n_ip3_contacts": n_contact,
                "n_filter_lining": n_filt, "n_gate_lining": n_gate,
                "check": check, "note": note,
            })
    return rows


def residues(rows: list[dict], paralog: str, definition: str) -> set[int]:
    for r in rows:
        if r["paralog"] == paralog and r["definition"] == definition:
            holes = []
            if r["excluded"]:
                for chunk in r["excluded"].split(";"):
                    a, b = chunk.split("-")
                    holes.append((int(a), int(b)))
            return {i for i in range(int(r["start"]), int(r["end"]) + 1)
                    if not any(a <= i <= b for a, b in holes)}
    raise KeyError(f"{paralog}/{definition}")


def module_of(rows: list[dict], paralog: str, resi: int,
              primary_only: bool = True) -> str:
    """Which primary module a residue is in — `other` if neither.

    A residue is never in both: the checks above guarantee the two primary
    definitions are disjoint at the sites that matter, and they are drawn
    from opposite ends of a 2,700-residue protein.
    """
    for module, definition in PRIMARY.items():
        if resi in residues(rows, paralog, definition):
            return module
    return "other"


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = build()
    cols = ["paralog", "module", "definition", "is_primary", "start", "end",
            "excluded", "n_residues", "n_ip3_contacts", "n_filter_lining",
            "n_gate_lining", "check", "note"]
    L.write_tsv(L.OUT_DIR / "module_map.tsv", rows, cols)
    for r in rows:
        if r["is_primary"]:
            L.log(f"{r['paralog']:6s} {r['module']:12s} "
                  f"{r['start']}-{r['end']} n={r['n_residues']:4d} "
                  f"{r['check']}")
    # disjointness, stated as a number rather than assumed
    for p in L.PARALOGS:
        a = residues(rows, p, PRIMARY["ligand_core"])
        b = residues(rows, p, PRIMARY["pore_module"])
        if a & b:
            raise ModuleRefusal(f"{p}: primary modules overlap at {len(a & b)}")
    L.log("primary modules disjoint in all three paralogues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
