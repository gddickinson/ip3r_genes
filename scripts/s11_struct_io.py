"""S11 structure IO — a deliberately small mmCIF/PDB reader.

Biopython is not used here on purpose. S11 only ever needs CA coordinates
per chain, the residue numbering and the B-factor column (which is where
AFDB puts pLDDT), and TM-align itself wants plain PDB. What a full parser
would add is failure modes on the very files this family produces.

**The tetramer trap.** Both families are homotetramers, and an IP3R
tetramer is ~11,000 residues. TM-align aligns two chains; handing it a
whole assembly would compare a tetramer against a monomer model and return
a score about quaternary structure. Every reference is reduced to **one
chain** (`largest_chain`) before it is scored, and the chain used is
recorded in the manifest.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


@dataclass
class Residue:
    chain: str
    seq_id: int
    name: str
    x: float
    y: float
    z: float
    bfactor: float = 0.0


@dataclass
class Chain:
    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        return "".join(AA3.get(r.name, "X") for r in self.residues)

    def __len__(self) -> int:
        return len(self.residues)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_pdb_ca(path: Path) -> dict[str, Chain]:
    """CA atoms per chain from a PDB file (first model, first altloc)."""
    chains: dict[str, Chain] = {}
    seen: set[tuple[str, int]] = set()
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            ch = line[21].strip() or "A"
            try:
                seq_id = int(line[22:26])
                x, y, z = (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))
                bf = float(line[60:66] or 0.0)
            except ValueError:
                continue
            if (ch, seq_id) in seen:
                continue
            seen.add((ch, seq_id))
            chains.setdefault(ch, Chain(ch)).residues.append(
                Residue(ch, seq_id, line[17:20].strip(), x, y, z, bf))
    return chains


def read_cif_ca(path: Path) -> dict[str, Chain]:
    """CA atoms per chain from an mmCIF `_atom_site` loop (first model)."""
    chains: dict[str, Chain] = {}
    cols: dict[str, int] = {}
    in_loop = False
    header = False
    first_model: str | None = None
    seen: set[tuple[str, int]] = set()
    with _open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                in_loop, header = True, True
                cols[s.split(".", 1)[1]] = len(cols)
                continue
            if header and not s.startswith("_atom_site."):
                header = False
            if not in_loop:
                continue
            if s.startswith("#") or s.startswith("loop_") or not s:
                if cols:
                    break
                continue
            parts = s.split()
            if len(parts) < len(cols):
                continue

            def get(name: str, default: str = "") -> str:
                idx = cols.get(name)
                return (parts[idx] if idx is not None and idx < len(parts)
                        else default)

            if get("group_PDB") not in ("ATOM", "HETATM"):
                continue
            if get("label_atom_id").strip('"') != "CA":
                continue
            alt = get("label_alt_id", ".")
            if alt not in (".", "?", "A"):
                continue
            model = get("pdbx_PDB_model_num", "1")
            if first_model is None:
                first_model = model
            if model != first_model:
                break
            ch = get("auth_asym_id") or get("label_asym_id") or "A"
            raw_id = get("auth_seq_id") or get("label_seq_id") or "0"
            try:
                seq_id = int(raw_id)
                x, y, z = (float(get("Cartn_x")), float(get("Cartn_y")),
                           float(get("Cartn_z")))
                bf = float(get("B_iso_or_equiv", "0") or 0)
            except ValueError:
                continue
            name = get("label_comp_id")
            if name not in AA3:
                continue
            if (ch, seq_id) in seen:
                continue
            seen.add((ch, seq_id))
            chains.setdefault(ch, Chain(ch)).residues.append(
                Residue(ch, seq_id, name, x, y, z, bf))
    return chains


def read_structure(path: Path) -> dict[str, Chain]:
    suffix = path.name.lower()
    if ".cif" in suffix:
        return read_cif_ca(path)
    return read_pdb_ca(path)


def write_ca_pdb(chain: Chain, dest: Path, chain_id: str = "A") -> Path:
    """Write one chain's CA trace as PDB.

    TM-align only uses CA atoms, and a CA-only file sidesteps the
    99,999-atom PDB limit an IP3R or RyR tetramer blows straight through.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for i, r in enumerate(chain.residues, start=1):
        lines.append(
            f"ATOM  {i % 100000:5d}  CA  {r.name:>3s} {chain_id}"
            f"{r.seq_id % 10000:4d}    "
            f"{r.x:8.3f}{r.y:8.3f}{r.z:8.3f}  1.00{r.bfactor:6.2f}           C"
        )
    lines.append("TER")
    lines.append("END")
    dest.write_text("\n".join(lines) + "\n")
    return dest


def largest_chain(chains: dict[str, Chain]) -> Chain:
    """The chain a homotetramer's fold question is actually about.

    Ties are broken on chain id so the choice is deterministic (D24): four
    identical chains of a C4 assembly differ only by resolved-residue count
    and two can easily be equal.
    """
    return max(chains.values(), key=lambda c: (len(c), [-ord(x) for x in c.chain_id]))


