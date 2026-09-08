"""S13 step 1 — the accepted species tree for the gene-tree taxa.

Reconciliation compares a *gene* tree against a *species* tree, so the
species tree is an input, not a result (**D15**): it must come from outside
this project's data or the comparison is circular. This module holds one,
hand-curated for the 31 vertebrate species S6 sampled, with every internal
node named and every age carrying its source.

**Ages are brackets, not point estimates.** S9 measured pairwise dS between
the paralogs far into saturation, so no molecular clock this project can run
would date these duplications. What reconciliation gives instead is a
*placement*: the duplication sits on a named branch of the species tree, and
that branch's endpoints bracket its age. Every number below is a literature
estimate for a **species-tree node**, never for a duplication.

**Why a `stem_age_ma` column exists.** The sampled tree omits every lineage
this project did not sequence, so a node's parent *in this tree* is usually
older than the split that actually created it: with no marsupial sampled,
Boreoeutheria's parent here is Amniota, and a duplication mapped to
Boreoeutheria would be bracketed [96, 319] Ma when the honest bracket is
[96, 99]. `stem_age_ma` is the crown age of the smallest accepted clade that
properly contains this one — the divergence that gave the clade its stem —
and `--check` requires `age <= stem_age <= parent age`. Where the sampled
parent *is* the true sister split (the whole deep backbone this task's answer
rests on), the two coincide and the source says so.

**Polytomies where the literature does not resolve.** Ovalentaria's internal
arrangement and the percomorph backbone are not fixed; asserting a pectinate
order there would put a shape into the input that the sources do not carry.
Those nodes are left multifurcating, which the reconciliation handles.

**Shared calibrations with the PIEZO project are shared on purpose.** The
deep vertebrate nodes carry the same ages as `../piezo_genes`'s curated table
because both read the same published sources. That is a shared *input*, not a
ported result.

    python scripts/s13_species_tree.py            # build + validate
    python scripts/s13_species_tree.py --check    # validate only
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s6_lib import binomial                                    # noqa: E402

OUT_DIR = ROOT / "results" / "reconciliation"
REPS_TSV = ROOT / "results" / "msa_v2" / "representatives.tsv"

VERT_GROUPS = ("ITPR1", "ITPR2", "ITPR3", "vertebrate_basal")

# --------------------------------------------------------------- calibrations
# node -> (crown age Ma, lo, hi, stem age Ma, source)
#
# Sources, abbreviated in the committed table:
#   TT5     TimeTree 5 (Kumar et al. 2022, MBE 39:msac174) median estimates
#   B15     Benton et al. 2015 (Palaeontologia Electronica 18.1.1FC) fossil
#           calibration min/max bounds
#   IRI17   Irisarri et al. 2017 (Nat Ecol Evol 1:1370) transcriptomic timetree
#   BET17   Betancur-R et al. 2017 (BMC Evol Biol 17:162) teleost classification
#           — topology only, ages from TT5
#   PRUM15  Prum et al. 2015 (Nature 526:569) avian backbone
#
# `stem` is the crown age of the smallest accepted clade properly containing
# this one. "= parent" means the sampled parent is that clade, so the sampling
# costs nothing at this node.
CALIBRATIONS: dict[str, tuple[float, float, float, float, str]] = {
    "Vertebrata":      (563.0, 480.0, 615.0, 615.0, "TT5 median 615; fossil-calibrated estimates 480-550 (B15, IRI17) - the widest disagreement in the tree; stem is the chordate root, taken as the upper bound"),
    "Cyclostomata":    (459.0, 380.0, 500.0, 563.0, "TT5 median hagfish-lamprey; IRI17 ~470; stem = parent (Vertebrata)"),
    "Gnathostomata":   (462.0, 421.0, 468.0, 563.0, "TT5 median 462; B15 crown bounds 420.7-468.4; stem = parent (Vertebrata)"),
    "Chondrichthyes":  (306.0, 190.0, 421.0, 462.0, "TT5 median (Holocephali-Elasmobranchii); palaeontological estimates run to ~420; stem = parent (Gnathostomata)"),
    "Elasmobranchii":  (264.0, 200.0, 310.0, 306.0, "TT5 median (Galeomorphi-Squalomorphi); stem = parent (Chondrichthyes)"),
    "Squalomorphi":    (190.0, 150.0, 235.0, 264.0, "TT5 median (Squaliformes-Echinorhiniformes); Echinorhinidae's position is unstable, hence the spread; stem = parent (Elasmobranchii)"),
    "Osteichthyes":    (431.0, 420.7, 453.0, 462.0, "TT5 median 431; B15 min 420.7; stem = parent (Gnathostomata)"),
    "Clupeocephala":   (240.0, 215.0, 265.0, 258.0, "TT5 median; stem = Osteoglossocephalai crown 258 (no osteoglossomorph sampled)"),
    "Ostariophysi":    (160.0, 140.0, 190.0, 200.0, "TT5 median (Cypriniformes-Gonorynchiformes); stem = Otomorpha crown 200 (no clupeomorph sampled)"),
    "Percomorphaceae": (120.0, 100.0, 140.0, 135.0, "TT5 median; stem = Acanthomorpha crown 135 (no non-percomorph acanthomorph sampled)"),
    "Ovalentaria":     (110.0,  95.0, 125.0, 120.0, "TT5 median; membership per BET17; stem = parent (Percomorphaceae)"),
    "Sarcopterygii":   (413.0, 390.0, 425.0, 431.0, "TT5 median (Actinistia-Tetrapodomorpha); stem = parent (Osteichthyes)"),
    "Tetrapoda":       (352.0, 330.0, 370.0, 413.0, "TT5 median; B15 min 337; stem = parent (Sarcopterygii)"),
    "Pipidae":         (122.0, 100.0, 145.0, 217.0, "TT5 median (Xenopodinae-Pipinae); stem = Batrachia crown 217 (no salamander or non-pipid frog sampled)"),
    "Amniota":         (319.0, 312.0, 330.0, 352.0, "TT5 median; B15 min 318; stem = parent (Tetrapoda)"),
    "Boreoeutheria":   ( 96.0,  88.0, 105.0,  99.0, "TT5 median; stem = Placentalia crown 99 (no atlantogenatan sampled)"),
    "Euarchontoglires": (88.0,  80.0,  96.0,  96.0, "TT5 median; stem = parent (Boreoeutheria)"),
    "Primates":        ( 43.0,  38.0,  50.0,  88.0, "TT5 median (Catarrhini-Platyrrhini); stem = Euarchontoglires crown (no non-primate euarchontan sampled)"),
    "Muridae":         ( 13.0,  10.0,  16.0,  88.0, "TT5 median (Mus-Rattus); stem = Euarchontoglires crown (no non-murid glire sampled)"),
    "Artiodactyla":    ( 62.0,  55.0,  70.0,  96.0, "TT5 median (Ruminantia-Cetacea); stem = Boreoeutheria crown (no non-artiodactyl laurasiathere sampled)"),
    "Sauria":          (279.0, 260.0, 295.0, 319.0, "TT5 median (Lepidosauria-Archelosauria); stem = parent (Amniota)"),
    "Lepidosauria":    (250.0, 238.0, 260.0, 279.0, "TT5 median (Rhynchocephalia-Squamata); stem = parent (Sauria)"),
    "Squamata":        (195.0, 175.0, 215.0, 250.0, "TT5 median (Gekkota-Unidentata); stem = parent (Lepidosauria)"),
    "Unidentata":      (170.0, 150.0, 190.0, 195.0, "TT5 median (Laterata-Toxicofera); stem = parent (Squamata)"),
    "Aves":            (110.0,  95.0, 125.0, 247.0, "TT5 median (Palaeognathae-Neognathae); PRUM15 backbone; stem = Archosauria crown 247 (no crocodilian sampled)"),
    "Neognathae":      (100.0,  88.0, 115.0, 110.0, "TT5 median (Galloanserae-Neoaves); PRUM15; stem = parent (Aves)"),
    "Telluraves":      ( 68.0,  60.0,  80.0, 100.0, "TT5 median (Afroaves-Australaves); PRUM15; stem = Neoaves crown, taken as the parent (no non-telluravian neoavian sampled)"),
    "Australaves":     ( 60.0,  52.0,  70.0,  68.0, "TT5 median (Psittaciformes-Passeriformes); PRUM15; stem = parent (Telluraves)"),
    "Passerida":       ( 25.0,  20.0,  32.0,  60.0, "TT5 median (Muscicapidae-Passerellidae); stem = Australaves crown (no suboscine sampled)"),
}

# ------------------------------------------------------------------ topology
# Nested tuples: (node name, [children]) for internal nodes, "species" for
# tips. Every internal node must have an entry in CALIBRATIONS.
#
# The backbone above Osteichthyes is the part any conclusion of this task
# rests on, and it is uncontroversial. Everything shallower is present so the
# species×paralog grid the loss audit reads has real branch structure; where
# the literature does not resolve an arrangement it is left as a polytomy
# rather than given one.
TOPOLOGY = (
    "Vertebrata", [
        ("Cyclostomata", ["Myxine glutinosa", "Petromyzon marinus"]),
        ("Gnathostomata", [
            ("Chondrichthyes", [
                "Callorhinchus milii",
                ("Elasmobranchii", [
                    "Chiloscyllium punctatum",
                    ("Squalomorphi", ["Somniosus microcephalus",
                                      "Echinorhinus cookei"]),
                ]),
            ]),
            ("Osteichthyes", [
                ("Clupeocephala", [
                    ("Ostariophysi", ["Danio rerio", "Chanos chanos"]),
                    ("Percomorphaceae", [
                        ("Ovalentaria", ["Nothobranchius furzeri",
                                         "Astatotilapia calliptera",
                                         "Salarias fasciatus"]),
                        "Nibea albiflora",
                    ]),
                ]),
                ("Sarcopterygii", [
                    "Latimeria chalumnae",
                    ("Tetrapoda", [
                        ("Pipidae", ["Xenopus tropicalis",
                                     "Hymenochirus boettgeri"]),
                        ("Amniota", [
                            ("Boreoeutheria", [
                                ("Euarchontoglires", [
                                    ("Primates", ["Homo sapiens",
                                                  "Saguinus oedipus"]),
                                    ("Muridae", ["Mus musculus",
                                                 "Rattus norvegicus"]),
                                ]),
                                ("Artiodactyla", ["Bos taurus",
                                                  "Sousa chinensis"]),
                            ]),
                            ("Sauria", [
                                ("Lepidosauria", [
                                    "Sphenodon punctatus",
                                    ("Squamata", [
                                        "Gekko japonicus",
                                        ("Unidentata", ["Podarcis muralis",
                                                        "Pogona vitticeps"]),
                                    ]),
                                ]),
                                ("Aves", [
                                    "Apteryx mantelli",
                                    ("Neognathae", [
                                        "Gallus gallus",
                                        ("Telluraves", [
                                            "Athene cunicularia",
                                            ("Australaves", [
                                                "Strigops habroptila",
                                                ("Passerida", [
                                                    "Copsychus sechellarum",
                                                    "Spizella passerina"]),
                                            ]),
                                        ]),
                                    ]),
                                ]),
                            ]),
                        ]),
                    ]),
                ]),
            ]),
        ]),
    ])


def species_label(name: str) -> str:
    """Newick-safe species label (spaces break Newick parsers)."""
    return name.replace(" ", "_")


def build_newick(node, parent_age: float | None = None) -> str:
    """Render TOPOLOGY as Newick with branch lengths in millions of years."""
    if isinstance(node, str):
        length = parent_age if parent_age is not None else 0.0
        return f"{species_label(node)}:{length:.6g}"
    name, children = node
    age = CALIBRATIONS[name][0]
    inner = ",".join(build_newick(c, age) for c in children)
    length = (parent_age - age) if parent_age is not None else 0.0
    if length < 0:
        raise ValueError(f"negative branch above {name}: parent is younger")
    return f"({inner}){name}:{length:.6g}"


def walk_internal(node, out=None, parent=None, acc=None):
    """[(node name, parent name or None)] in preorder."""
    acc = [] if acc is None else acc
    if isinstance(node, str):
        return acc
    name, children = node
    acc.append((name, parent))
    for c in children:
        walk_internal(c, parent=name, acc=acc)
    return acc


def walk_tips(node, out=None):
    out = [] if out is None else out
    if isinstance(node, str):
        out.append(node)
        return out
    for c in node[1]:
        walk_tips(c, out)
    return out


def n_tips_under(name: str) -> int:
    def rec(node) -> tuple[bool, int]:
        if isinstance(node, str):
            return False, 1
        nm, children = node
        total = 0
        for c in children:
            hit, n = rec(c)
            if hit:
                return True, n
            total += n
        return (nm == name), total
    hit, n = rec(TOPOLOGY)
    return n if hit else 0


def gene_tree_species() -> set[str]:
    """Binomial of every vertebrate tip in the msa_v2 representative set."""
    with REPS_TSV.open() as fh:
        return {binomial(r["species"]) for r in csv.DictReader(fh, delimiter="\t")
                if r["group"] in VERT_GROUPS}


def check() -> list[str]:
    """Every consistency test the tree has to pass, as failure strings."""
    problems: list[str] = []
    tips = walk_tips(TOPOLOGY)
    if len(tips) != len(set(tips)):
        problems.append("duplicate tip in TOPOLOGY")
    needed = gene_tree_species()
    missing = sorted(needed - set(tips))
    extra = sorted(set(tips) - needed)
    if missing:
        problems.append(f"gene-tree species absent from species tree: {missing}")
    if extra:
        problems.append(f"species tree carries taxa with no gene-tree tip: {extra}")
    ages = {}
    for name, parent in walk_internal(TOPOLOGY):
        if name not in CALIBRATIONS:
            problems.append(f"internal node without a calibration: {name}")
            continue
        age, lo, hi, stem, _ = CALIBRATIONS[name]
        ages[name] = (age, stem)
        if not lo <= age <= hi:
            problems.append(f"{name}: point age {age} outside its spread {lo}-{hi}")
        if stem < age:
            problems.append(f"{name}: stem age {stem} younger than crown age {age}")
    for name, parent in walk_internal(TOPOLOGY):
        if parent is None or name not in ages or parent not in ages:
            continue
        if ages[name][1] > ages[parent][0]:
            problems.append(f"{name}: stem age {ages[name][1]} older than its "
                            f"sampled parent {parent} ({ages[parent][0]})")
    for name in CALIBRATIONS:
        if name not in ages:
            problems.append(f"calibration for a node not in TOPOLOGY: {name}")
    try:
        build_newick(TOPOLOGY)
    except ValueError as exc:                      # age inversion
        problems.append(str(exc))
    except KeyError as exc:                        # already reported above
        problems.append(f"missing calibration while rendering: {exc}")
    return problems


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="validate only; write nothing")
    args = ap.parse_args()

    problems = check()
    for p in problems:
        print(f"  ! {p}")
    if problems:
        raise SystemExit("species tree failed validation")
    tips = walk_tips(TOPOLOGY)
    internals = walk_internal(TOPOLOGY)
    print(f"species tree: {len(tips)} tips, {len(internals)} named internal nodes")
    print(f"all {len(gene_tree_species())} gene-tree species present")
    if args.check:
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nwk = build_newick(TOPOLOGY) + ";"
    (OUT_DIR / "species_tree.nwk").write_text(nwk + "\n")
    print(f"wrote {OUT_DIR/'species_tree.nwk'}")

    with (OUT_DIR / "species_tree_calibrations.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["node", "parent", "age_ma", "age_lo", "age_hi",
                    "stem_age_ma", "n_tips", "source"])
        for name, parent in internals:
            age, lo, hi, stem, src = CALIBRATIONS[name]
            w.writerow([name, parent or "", age, lo, hi, stem,
                        n_tips_under(name), src])
    print(f"wrote {OUT_DIR/'species_tree_calibrations.tsv'} ({len(internals)} nodes)")


if __name__ == "__main__":
    main()
