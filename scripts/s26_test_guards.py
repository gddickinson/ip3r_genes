"""S26: break every series guard on purpose, and check that it fires (D74).

The thesis's pattern (`s25_test_guards.py`) applied to the paper series.
Each case breaks one rule, runs the stage that enforces it, and requires a
non-zero status **and** a declared fragment of that guard's own message. A
case that fails for some other reason is itself a failure, not a pass.

Every case runs against a copy of `papers/` in a temporary directory, with
`s26_lib.PAPERS` pointed there. Configuration attributes a case changes are
restored afterwards. The suite records the SHA-256 of every committed file
under `papers/` and `scripts/` before and after, and fails if any has moved
(D60): a self-test that damages what it tests is a failure mode this project
has already had.

  python scripts/s26_test_guards.py
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s25_pdf                                                  # noqa: E402
import s26_assign                                               # noqa: E402
import s26_assignment as A                                      # noqa: E402
import s26_claims                                               # noqa: E402
import s26_figures                                              # noqa: E402
import s26_lib as L                                             # noqa: E402
import s26_pdf                                                  # noqa: E402
import s26_prose                                                # noqa: E402
import s26_rules                                                # noqa: E402
import s26_stitch                                               # noqa: E402

OUT_NAME = "guard_check.tsv"
FIELDS = ["guard", "expected_message", "fired", "message_matched", "verdict",
          "said"]
RESULTS: list[dict] = []
#: The paper the section-level cases are run against: the first in the
#: series, so a case never depends on a later paper having been written.
P = L.SERIES[0]


def _digests() -> dict[str, str]:
    out = {}
    for base in (L.PAPERS, L.ROOT / "scripts"):
        for p in sorted(base.rglob("*")):
            if (p.is_file() and p.name != OUT_NAME
                    and "__pycache__" not in p.parts):
                out[str(p)] = L.sha256(p)
    return out


@contextlib.contextmanager
def sandbox():
    tmp = Path(tempfile.mkdtemp(prefix="s26_guards_"))
    shutil.copytree(L.PAPERS, tmp / "papers",
                    ignore=shutil.ignore_patterns("*.pdf"))
    old = L.PAPERS
    L.PAPERS = tmp / "papers"
    try:
        yield L.PAPERS
    finally:
        L.PAPERS = old
        shutil.rmtree(tmp, ignore_errors=True)


@contextlib.contextmanager
def patched(obj, name, value):
    saved = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, saved)


def case(name: str, expect: str, fn) -> None:
    err, out = io.StringIO(), io.StringIO()
    with sandbox():
        try:
            with contextlib.redirect_stderr(err), \
                    contextlib.redirect_stdout(out):
                status = fn()
        except Exception as exc:                             # noqa: BLE001
            status, err = 1, io.StringIO(f"{type(exc).__name__}: {exc}")
    text = err.getvalue() + out.getvalue()
    fired, matched = status != 0, expect.lower() in text.lower()
    said = next((ln.strip() for ln in text.splitlines()
                 if expect.lower() in ln.lower()),
                text.strip().split("\n")[-1] if text.strip() else "(none)")
    RESULTS.append({"guard": name, "expected_message": expect,
                    "fired": int(fired), "message_matched": int(matched),
                    "verdict": "fired" if fired and matched
                    else "DID NOT FIRE", "ok": fired and matched,
                    "said": said[:200]})


def _append(pid: str, section: str, text: str) -> None:
    f = L.paper_dir(pid) / section
    f.write_text(f.read_text(encoding="utf-8") + text, encoding="utf-8")


# ------------------------------------------------------------------ cases

def _assign_unassigned():
    real = L.results_entries
    with patched(L, "results_entries", lambda: real() + ["results/zz_new"]):
        return s26_assign.run()


def _assign_unmatched_file():
    real = s26_assign._files
    with patched(s26_assign, "_files", lambda e: real(e) + (
            ["zz_stray.tsv"] if e == "results/methods" else [])):
        return s26_assign.run()


def _assign_undeclared_departure():
    deps = dict(A.DEPARTURES)
    deps.pop("results/census_v4")
    with patched(A, "DEPARTURES", deps):
        return s26_assign.run()


def _rules_question_and():
    with patched(L.paper(P), "QUESTION", "Where is it, and why?"):
        return s26_rules.run()


def _rules_sibling_control():
    ctl = L.paper(P).CONTROLS + [("a borrowed control",
                                  "results/phylogeny/au_test.tsv")]
    with patched(L.paper(P), "CONTROLS", ctl):
        return s26_rules.run()


def _rules_sibling_standalone():
    """A paper's standalone answer leaning on a number it only cites.

    The paper and the claim are found rather than named, because the first
    paper cites nothing and so carries no sibling's claim to lean on.
    """
    merged, _ = s26_claims.ledger()
    pid, foreign = next((p, cid) for cid, c in merged.items()
                        for p in c["papers"]
                        if c["primary"] not in ("", "-", p))
    cfg = L.paper(pid)
    with patched(cfg, "STANDALONE_CLAIMS",
                 list(cfg.STANDALONE_CLAIMS) + [foreign]):
        return s26_rules.run()


def _rules_cites_later():
    _append(P, "03_discussion.md", f"\n\nSee {{paper:{L.SERIES[-1]}}}.\n")
    return s26_rules.run()


def _rules_cycle():
    rows = [{"citing": "range", "cited": "origin"},
            {"citing": "origin", "cited": "range"}]
    fails = s26_rules._cycles(rows)
    for f in fails:
        print(f, file=sys.stderr)
    return 1 if fails else 0


def _figures_too_few():
    figs = [f for f in L.paper(P).FIGURES if f[1] != "main"] + [
        f for f in L.paper(P).FIGURES if f[1] == "main"][:3]
    with patched(L.paper(P), "FIGURES", figs):
        return s26_figures.run(P)


def _figures_sibling():
    figs = list(L.paper(P).FIGURES) + [
        ("borrowed", "ed", "phylogeny/figures/sister_au", "borrowed")]
    with patched(L.paper(P), "FIGURES", figs):
        return s26_figures.run(P)


def _figures_undescribed():
    figs = list(L.paper(P).FIGURES)
    slug, kind, stems, cap = next(f for f in figs if f[1] == "ed")
    figs.append((slug + "_twin", kind, stems, cap))
    with patched(L.paper(P), "FIGURES", figs):
        return s26_figures.run(P)


def _prose_em_dash():
    _append(P, "02_results.md", "\n\nA sentence — with a dash.\n")
    return s26_prose.run(P)


def _prose_heading_fragment():
    _append(P, "02_results.md", "\n\n### The shape of the census\n\nIt is "
            "measured.\n")
    return s26_prose.run(P)


def _prose_legend_fragment():
    slug = next(f[0] for f in L.paper(P).FIGURES if f[1] == "main")
    f = L.paper_dir(P) / "05_figure_legends.md"
    text = f.read_text(encoding="utf-8")
    head = f"**{{fig:{slug}}}.**"
    i = text.index(head) + len(head)
    f.write_text(text[:i] + " Distribution of lengths by kingdom." + text[i:],
                 encoding="utf-8")
    return s26_prose.run(P)


def _stitch_unknown_key():
    _append(P, "02_results.md", "\n\nThis is cited [R999].\n")
    return s26_stitch.run(P)


def _stitch_unknown_paper():
    _append(P, "02_results.md", "\n\nThis is shown in {paper:nowhere}.\n")
    return s26_stitch.run(P)


def _claims_wrong_value():
    mod = s26_claims._module(P)
    bad = dict(mod.CLAIMS[0], expect="999999") if mod.CLAIMS else None
    if bad is None:
        print("no new claim in the first paper to break", file=sys.stderr)
        return 0
    with patched(mod, "CLAIMS", [bad] + mod.CLAIMS[1:]):
        return s26_claims.run()


def _claims_padded():
    """A true number the paper never states: the value is computed from its
    table first, so the only thing wrong with the claim is its paper."""
    import s14_claims
    import s25_claims
    mod = s26_claims._module(P)
    text = s26_claims.paper_text(P)
    for src in sorted((L.RESULTS / "genome_ledger").glob("*.tsv")):
        pad = dict(id=s26_claims.PREFIX[P] + "ZZ", claim="never stated",
                   source=str(src.relative_to(L.ROOT)), op="count",
                   expect="")
        pad["expect"], _ = s14_claims._evaluate(pad)
        if not any(r in text for r in s25_claims._renderings(pad["expect"])):
            break
    with patched(mod, "CLAIMS", list(mod.CLAIMS) + [pad]):
        return s26_claims.run()


def _claims_twice():
    mod = s26_claims._module(P)
    if not mod.CLAIMS:
        print("no new claim in the first paper to copy", file=sys.stderr)
        return 0
    twin = dict(mod.CLAIMS[0], id=s26_claims.PREFIX[P] + "ZY")
    with patched(mod, "CLAIMS", list(mod.CLAIMS) + [twin]):
        return s26_claims.run()


def _pdf_legend_unwrapped():
    body = "Text.\n\n**Fig. 1.** The figure shows a thing.\n\nMore text."
    problems = s25_pdf.check_legends(body, s26_pdf.LEGEND_PARA)
    for p in problems:
        print(p, file=sys.stderr)
    return 1 if problems else 0


def _pdf_glyph():
    problems = s25_pdf.check_log(
        "Missing character: There is no ∮ in font texgyretermes!\n")
    for p in problems:
        print(p, file=sys.stderr)
    return 1 if problems else 0


CASES = [
    ("assign: a results entry nobody assigned (P5)", "neither assigned",
     _assign_unassigned),
    ("assign: a file in a split directory no group matches (P5s)",
     "matches 0 file groups", _assign_unmatched_file),
    ("assign: a departure from the thesis nobody declared",
     "departure is not declared", _assign_undeclared_departure),
    ("rules: a question with an 'and' in it (P1)", "contains 'and'",
     _rules_question_and),
    ("rules: a control measured in a sibling paper (P2)", "(P2)",
     _rules_sibling_control),
    ("rules: a standalone claim resting on a sibling's table (P6)",
     "does not stand alone", _rules_sibling_standalone),
    ("rules: a paper citing a later paper", "not earlier in the submission",
     _rules_cites_later),
    ("rules: two papers citing each other", "citation cycle", _rules_cycle),
    ("figures: three main figures (P4)", "P4 allows", _figures_too_few),
    ("figures: a figure drawn from a sibling's results (P5)", "(P5)",
     _figures_sibling),
    ("figures: a figure no body sentence mentions (R3)",
     "no body sentence", _figures_undescribed),
    ("prose: an em-dash in a section", "em-dash", _prose_em_dash),
    ("prose: a results heading with no finite verb",
     "results heading has no finite verb", _prose_heading_fragment),
    ("prose: a legend opening with a noun phrase (R1)",
     "opens without a finite verb", _prose_legend_fragment),
    ("stitch: a cited key with no reference row",
     "no row in references.tsv", _stitch_unknown_key),
    ("stitch: a companion citation to no paper", "names no paper",
     _stitch_unknown_paper),
    ("claims: a number the table does not produce", "expected '999999'",
     _claims_wrong_value),
    ("claims: a claim its paper does not state", "does not state",
     _claims_padded),
    ("claims: one number declared under two identifiers", "declared twice",
     _claims_twice),
    ("pdf: a legend set as body text (R4)", "not wrapped",
     _pdf_legend_unwrapped),
    ("pdf: a glyph the document font lacks", "missing characters",
     _pdf_glyph),
]


def run() -> int:
    before = _digests()
    RESULTS.clear()
    for name, expect, fn in CASES:
        case(name, expect, fn)
    after = _digests()
    tampered = sorted(k for k in set(before) | set(after)
                      if before.get(k) != after.get(k))
    L.write_tsv(L.PAPERS / OUT_NAME, RESULTS, FIELDS)
    n_ok = sum(1 for r in RESULTS if r["ok"])
    for r in RESULTS:
        if not r["ok"]:
            print(f"  [FAIL] {r['guard']}: said {r['said']!r}",
                  file=sys.stderr)
    if tampered:
        print(f"  [FAIL] the suite altered {len(tampered)} committed "
              f"file(s): {tampered[:3]}", file=sys.stderr)
    print(f"[s26 guards] {n_ok}/{len(RESULTS)} guards fired with the "
          f"expected message; {len(before)} committed files unchanged"
          if not tampered else f"[s26 guards] {n_ok}/{len(RESULTS)} fired; "
          f"TAMPERED")
    return 0 if n_ok == len(RESULTS) and not tampered else 1


if __name__ == "__main__":
    sys.exit(run())
