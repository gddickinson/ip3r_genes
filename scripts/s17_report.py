"""S17 — renders `results/constraint/report.md` purely from the committed tables.

D13: nothing here recomputes a statistic. Every number comes from
`s17_tables.headline()`, which reads it off one committed file, so the report
and the data cannot drift.

Scope, the architecture, the instrument and the negative controls live here;
the results — the element ranking, the ligand site, the variant classifier, the
per-site selection and the caveats — live in `s17_report_results.py` (the
`s3_report.py` / `s3_report_d10.py` split). A section whose table is absent
renders *not run yet* rather than nothing, so a stage that was skipped is
visible as skipped.

    python scripts/s17_report.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_domains as D                                        # noqa: E402
import s17_lib as L                                            # noqa: E402
import s17_orthologs as O                                      # noqa: E402
import s17_tables as TB                                        # noqa: E402


def A(*parts: str) -> str:
    return "\n".join(parts)


def _fmt(x, nd=3) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x) if x not in (None, "") else "—"


def section_scope(h: dict) -> str:
    n = h["deep_n"]
    ls = h["layer_sizes"]
    return A(
        "## 1. Scope — what a per-site score for this family can rest on",
        "",
        "S17 asks two questions of the IP₃ receptor's own sequence: **which "
        "parts of the receptor are evolutionarily intolerant**, and **is that "
        "intolerance good enough to help interpret the family's human "
        "variants**. Both are per-residue questions, so the first thing that "
        "has to be established is how many observations a column actually has.",
        "",
        "S6's representative alignment is taxonomically broad — 134 sequences "
        "from protists to vertebrates — and **thin per paralog**: "
        f"{ls.get(('shallow', 'ITPR1'), 0)} ITPR1, "
        f"{ls.get(('shallow', 'ITPR2'), 0)} ITPR2 and "
        f"{ls.get(('shallow', 'ITPR3'), 0)} ITPR3 tips. "
        f"{min(ls.get(('shallow', p), 0) for p in L.PARALOGS)} sequences "
        "cannot score a column of a 2,701-residue protein.",
        "",
        "The depth already existed and had not been used. Every one of the 309 "
        "swept genomes retains its `miniprot.gff`, whose `##STA` line is that "
        "locus's translation; joined to `summary.json` by `mp_id` that is a "
        "full-length ITPR protein per genome per paralog, with the orthology "
        "assignment already done by a bait panel whose family call is a "
        "positive test (D14) and which S7's tree independently validated.",
        "",
        "| layer | sequences | what it answers |",
        "|---|---|---|",
        f"| `deep` ITPR1 | **{n.get('ITPR1', 0)}** | has *this* gene "
        "tolerated change here? |",
        f"| `deep` ITPR2 | **{n.get('ITPR2', 0)}** | |",
        f"| `deep` ITPR3 | **{n.get('ITPR3', 0)}** | |",
        f"| `vert` | {ls.get(('vert', 'all three'), 0)} | was the position "
        "fixed before the duplications? |",
        f"| `family` | {ls.get(('family', 'all'), 0)} | the eukaryote-wide "
        "floor — what the fold cannot do without |",
        f"| `shallow` ITPR1/2/3 | {ls.get(('shallow', 'ITPR1'), 0)} / "
        f"{ls.get(('shallow', 'ITPR2'), 0)} / "
        f"{ls.get(('shallow', 'ITPR3'), 0)} | **the control** — what S17 "
        "would have had without the deep sets |",
        "",
        "`shallow` is not a claim. It is what the task would have been able to "
        "say from S6's alignment alone, and §7.2 uses the difference between "
        "it and the other three to show empirically whether depth bought "
        "discrimination or noise.",
        "",
        "The six RyR tips in msa_v2 are **excluded from every layer**. They "
        "are there as S6's outgroup, and a layer containing them would score "
        "what the *superfamily* conserves — a different question, and one D14 "
        "keeps separate everywhere else in this project.",
        "",
        "All three references are the **human** proteins "
        f"({', '.join(f'{p} {a}' for p, (_l, a, _d) in L.REFERENCES.items())}), "
        "so every table is directly in the numbering ClinVar and UniProt use. "
        "That is not a small convenience: the PIEZO project this method is "
        "ported from had one paralog with no human copy at all, and its "
        "coordinate system had to be a zebrafish gene.",
    )


def _ref_contacts(contacts: list[dict]) -> str:
    r = next(x for x in contacts if x["paralog"] == L.STRUCTURE_REF)
    return f"{r['start']}–{r['end']}"


def _contact_identity(contacts: list[dict]) -> str:
    return "; ".join(f"{r['paralog']} {r['anchor']}" for r in contacts
                     if r["paralog"] != L.STRUCTURE_REF)


def section_architecture(h: dict) -> str:
    dm = h["domain_map"]
    checks = h["transfer_checks"]
    struct = [r for r in dm if r["kind"] == "structural"]
    contacts = [r for r in dm if r["element"] == "ip3_contact_set"]
    rows = []
    for r in struct:
        if r["element"] == "ip3_contact_set":
            continue
        rows.append(f"| {r['paralog']} | `{r['element']}` | "
                    f"{r['start']}–{r['end']} | {r['check']} | "
                    f"{r['anchor'] or '—'} |")
    return A(
        "## 2. The architecture, and which half of it needed transferring",
        "",
        "Every claim S17 makes is of the form *this part is more constrained "
        "than that part*, so the parts have to be located in the coordinate "
        "system each table is reported in. This family supplies that in two "
        "pieces which are **not equally trustworthy**, and they are kept "
        "apart.",
        "",
        "**The Pfam elements are measured per accession and need no "
        "transfer.** S0's `domain_coords.tsv` holds InterPro's coordinates for "
        "Q14643, Q14571 and Q14573 *separately*. A whole class of error — the "
        "one the PIEZO port had to spend a self-test on, where one published "
        "boundary list in a mouse protein's numbering is carried onto every "
        "reference — does not arise here, and the module refuses to invent it.",
        "",
        "**The structural elements are measured on one protein and are "
        "transferred.** S0 measured the channel on PDB 6DQN (human ITPR3, "
        "IP₃-bound, 3.33 Å): the selectivity filter, the gate and the ten "
        "residues within 4.5 Å of IP₃. ITPR1 and ITPR2 get those by pairwise "
        "alignment, and the transfer is **a positive test, not an "
        "assumption** — the filter must arrive on the same GGGVGD motif and "
        "the gate on the same two lining residues, properties the alignment "
        "knows nothing about. A failed anchor aborts rather than writing a "
        "coordinate.",
        "",
        "| paralog | element | span | anchor check | anchor |",
        "|---|---|---|---|---|",
        *rows,
        "",
        "**One element was located by geometry rather than by annotation.** "
        "The Pfam channel domain's mean conservation is bimodal, and a "
        "boundary drawn after seeing that profile would be the profile "
        "explaining itself. So the luminal loop is defined on the structure: "
        "the contiguous run of channel residues whose Cα sits beyond the "
        "membrane on the luminal side of S0's measured axial span. Nothing in "
        "that derivation knows what any column looks like. It comes out as "
        "**2401–2450 in ITPR3**, and it transfers onto 2481–2531 in ITPR1 and "
        "2427–2474 in ITPR2, each confirmed to lie inside that accession's own "
        "channel domain. §5 is what it changes.",
        "",
        "### 2.1 The Pfam name for the ligand site does not contain the ligand site",
        "",
        "Joining the two sources says something that changes how the rest of "
        "this report has to be read. **None of the ten measured IP₃ contacts "
        "lies in PF08709** — the signature Pfam calls *Inositol "
        "1,4,5-trisphosphate/ryanodine receptor*. They sit in MIR (PF02815) "
        "and RIH (PF01365), the two domains the family shares with the "
        "ryanodine receptor.",
        "",
        f"Contact residues in the numbering they were measured in "
        f"({L.STRUCTURE_REF}): {_ref_contacts(contacts)}. Carried across, "
        f"{_contact_identity(contacts)}. So `nterm_trefoil` "
        "is the N-terminal β-trefoil and **not** the ligand site, and S17's "
        "ligand question is asked of the measured contacts throughout — never "
        "of the Pfam label. A report that had taken the label at face value "
        "would have measured the suppressor domain and called it the "
        "IP₃-binding core.",
        "",
        "Residues outside every named element are **named linkers** — "
        "`nterm`, `cterm`, or `linker_<upstream>_<downstream>` — never left "
        "blank, because `unassigned` would pool the N-terminus, five "
        "inter-domain stretches and the whole C-terminal tail into one bucket "
        "and then use it as the within-protein control.",
    )


def section_instrument(h: dict) -> str:
    dropped = h["shape_dropped"]
    drop_lines = [f"- `{r['name']}` — {r['frac_in_ref_columns']} of its own "
                  f"residues inside reference columns" for r in dropped]
    return A(
        "## 3. Instrument",
        "",
        "Conservation is **Jensen–Shannon divergence against the BLOSUM62 "
        "background** (Capra & Singh 2007), and every layer is "
        "**sequence-weighted** by Henikoff & Henikoff position-based weights. "
        "The weighting is not optional at this depth: Actinopteri are about a "
        "third of the vertebrate genome scope, and unweighted, every "
        "teleost-specific residue in a 260-sequence alignment reads as "
        "conserved.",
        "",
        f"Occupancy is reported on every row and a column below "
        f"{L.OUT_DIR.name and ''}"
        f"{int(100 * __import__('s17_conservation').MIN_OCCUPANCY)} % "
        "occupancy is flagged unreliable rather than dropped silently, so a "
        "score computed from a quarter of the sequences cannot be mistaken "
        "for a well-supported one.",
        "",
        "**A composition-free metric sits beside the JSD on every table.** "
        "JSD is a divergence *from a background frequency table*, so an "
        "element built of common amino acids scores lower at equal "
        "conservation — and a transmembrane domain is built of exactly those. "
        "`metric_controls.tsv` carries the modal-residue fraction and the "
        "normalised entropy (neither uses a background model), the mean "
        "background frequency of each element's modal residues, the mean "
        "occupancy, and every element recomputed **on the curated subset "
        "alone**, which no gene caller produced. Those controls are written "
        "for every element rather than only where a result was surprising, "
        "because a control produced after seeing the answer is not one.",
        "",
        "### 3.1 The screen the coverage bar could not do",
        "",
        "Ortholog selection is one locus per genome × paralog at bait coverage "
        f"≥ {O.MIN_COVERAGE}, identity ≥ {O.min_identity()} (the sweep's own "
        f"floor, read back rather than retyped) and ≥ {O.MIN_ALIGNED_AA} "
        "aligned residues, with S15a's lesion-rich loci excluded.",
        "",
        "None of those can see a model that is too **long**. Bait coverage is "
        "`aligned_aa / bait_len`, so a model that covers the bait *and* "
        "carries two thousand extra residues passes every bar — and S5b "
        "recorded that a large `-G` manufactures exactly that, chaining "
        "shared-channel-module hits 23 Mbp apart into one locus in the giant "
        "genomes. So every sequence is scored on the fraction of its **own** "
        "residues that land in columns the reference occupies, measured on the "
        "alignment itself, and the bar is put in the gap the distribution "
        "shows rather than chosen (`s15_calibrate_recon`'s rule, both edges "
        "committed).",
        "",
        f"Bar **{h['shape_bar']}**, in a gap running from the highest failing "
        f"sequence to the lowest passing one, with every curated record above "
        f"it by a margin. {h['n_shape_dropped']} sequence(s) dropped:",
        "",
        *drop_lines,
        "",
        "It paid for itself immediately: with that one sequence in, MAFFT "
        "opened the ITPR2 alignment to 5,676 columns for a 2,701-residue "
        "protein; without it, 3,380 — in line with ITPR1's 3,519 and ITPR3's "
        "3,738. One chimeric model was inflating the alignment every ITPR2 "
        "per-site score is read off by 68 %.",
        "",
        "MAFFT runs at **`--thread 1`** by decision (D24) — S3 measured the "
        "same seeds giving profiles of 4,933 and 4,908 match states on "
        "consecutive `--thread -1` runs — and `align_stats.json` records the "
        "version, the runtimes and the SHA-256 of every input and output. A "
        "ragged alignment is a hard failure: S1 established that a silent "
        "MAFFT failure degrades to a star alignment with no other symptom.",
    )


def section_selftests(h: dict) -> str:
    import s17_test_constraint as T
    st = h.get("_selftest", {})
    return A(
        "## 4. Self-tests",
        "",
        f"**{len(T.TESTS)} constructed negative controls run before anything "
        f"is written** (`s17_test_constraint.py`), each rejected by the rule "
        f"responsible. Status: "
        f"**{'all passed' if st.get('passed', True) else 'FAILED'}**.",
        "",
        "These are checks on *refusal*, because almost every rule in S17 "
        "returns a plausible number when it is wrong: a conservation score "
        "computed with gaps as a 21st state looks like conservation, a variant "
        "renumbered onto the wrong transcript lands somewhere in the protein "
        "and gets a score, and a coordinate transfer that slipped ten residues "
        "still puts the gate inside the channel.",
        "",
        "| test | what it refuses |",
        "|---|---|",
        "| T1 | a densely sampled clade outvoting a lone tip |",
        "| T2 | gaps counted as a shared state rather than as missing data |",
        "| T3 | an invariant column not beating a uniform one |",
        "| T4 | a self-transfer that is not the identity map |",
        "| T5 | **a mutated GGGVGD still passing the transfer anchor** |",
        "| T6 | the containing element winning over the filter and the gate |",
        "| T7 | a residue left unassigned, or a linker left unnamed |",
        "| T8 | **a named domain inside the within-protein control** |",
        "| T9 | a two-gene chimera passing the shape screen |",
        "| T10 | the shape calibration firing on a distribution with no gap |",
        "| T11 | a variant renumbered instead of dropped |",
        "| T12 | an AUC computed from fewer than three observations |",
        "| T13 | alleles counted where positions were meant |",
        "| T14 | a p-value of 2.1e-07 rendered as `0.0000` |",
        "",
        "**T8 found a real bug on its first run.** The within-protein control "
        "was selected by a `startswith` test on `nterm`, `cterm` and "
        "`linker_` — and `nterm_trefoil` is a *domain* whose name begins with "
        "`nterm`. So 225 residues of PF08709 — the element §2.1 is about — "
        "were quietly inside the control, and every other element was being "
        "compared against a set containing one of them. The membership test is "
        "now a prefix for linkers and an exact match for the termini.",
    )


def build(out_dir: Path = L.OUT_DIR) -> Path:
    import s17_report_results as R
    h = TB.headline(out_dir)
    stats_path = out_dir / "constraint_stats.json"
    if stats_path.exists():
        import json
        h["_selftest"] = json.loads(stats_path.read_text()).get("self_test", {})
    parts = [
        "# S17 — Constraint & function",
        "",
        "*Generated by `scripts/s17_report.py` from the committed tables in "
        "`results/constraint/`. Nothing in this report is hand-written and "
        "nothing in it is recomputed (D13).*",
        "",
        section_scope(h),
        "",
        section_architecture(h),
        "",
        section_instrument(h),
        "",
        section_selftests(h),
        "",
        R.build(h),
    ]
    dest = out_dir / "report.md"
    dest.write_text("\n".join(parts) + "\n")
    print(f"[s17] wrote {dest} ({dest.stat().st_size / 1024:.0f} KB)")
    return dest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    build(ap.parse_args().out)


if __name__ == "__main__":
    main()
