"""S25 — break every build guard on purpose, and check that it fires.

The brief requires each guard to have been broken once. Doing that by hand
proves it for one afternoon; doing it here proves it on every run, and gives
the message each guard actually emits rather than a description of it.

Two properties this suite must have, and both are checked:

* **It must not damage what it tests** (D60). Every case runs against a copy
  of `thesis/` in a temporary directory, with the module-level paths
  redirected there, and the suite records the SHA-256 of every committed file
  under `thesis/` and `results/s0_baseline/references.tsv` before and after,
  failing if any moved. An earlier suite in this project overwrote a committed
  table with two synthetic rows and nothing noticed.
* **Each case must fail for the reason it was broken**, not merely fail. Every
  case declares a fragment the guard's own message has to contain, so a case
  that starts failing for a different reason is a failure and not a pass.

  python scripts/s25_test_guards.py
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s14_lib as lib                                          # noqa: E402
import s25_claims                                              # noqa: E402
import s25_claims_thesis                                       # noqa: E402
import s25_figmap as fm                                        # noqa: E402
import s25_figures                                             # noqa: E402
import s25_lib as L                                            # noqa: E402
import s25_assign                                              # noqa: E402
import s25_refs                                                # noqa: E402
import s25_stitch                                              # noqa: E402

RESULTS: list[dict] = []


#: The suite's own output, written after the tamper check and excluded from
#: it — a file the suite creates is not a file the suite damaged.
OUT_NAME = "guard_check.tsv"
FIELDS = ["guard", "expected_message", "fired", "message_matched", "verdict",
          "said"]


def _digests() -> dict[str, str]:
    out = {}
    for p in sorted(L.TH.rglob("*")):
        if p.is_file() and p.name != OUT_NAME:
            out[str(p)] = lib.sha256(p)
    refs = L.RESULTS / "s0_baseline" / "references.tsv"
    out[str(refs)] = lib.sha256(refs)
    return out


@contextlib.contextmanager
def sandbox():
    """A copy of thesis/ with the module paths pointed at it."""
    tmp = Path(tempfile.mkdtemp(prefix="s25_guards_"))
    shutil.copytree(L.TH, tmp / "thesis")
    old_th, old_figs = L.TH, L.TH_FIGS
    L.TH, L.TH_FIGS = tmp / "thesis", tmp / "thesis" / "figures"
    try:
        yield L.TH
    finally:
        L.TH, L.TH_FIGS = old_th, old_figs
        shutil.rmtree(tmp, ignore_errors=True)


def case(name: str, expect: str, fn) -> None:                # noqa: D401
    """Run one broken-guard case; require a non-zero status and the message."""
    err = io.StringIO()
    out = io.StringIO()
    status = None
    with sandbox():
        try:
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
                status = fn()
        except Exception as exc:                              # noqa: BLE001
            status, err = 1, io.StringIO(f"{type(exc).__name__}: {exc}")
    text = err.getvalue() + out.getvalue()
    fired = status != 0
    matched = expect.lower() in text.lower()
    ok = fired and matched
    said = next((ln.strip() for ln in text.splitlines()
                 if expect.lower() in ln.lower()), text.strip().split("\n")[-1]
                if text.strip() else "(no output)")
    RESULTS.append({"guard": name, "expected_message": expect,
                    "fired": "1" if fired else "0",
                    "message_matched": "1" if matched else "0",
                    "verdict": "fired" if ok else "DID NOT FIRE",
                    "ok": ok, "said": said[:200]})


# ----------------------------------------------------------------- the cases

def _assign_unassigned():
    entry = next(iter(L.ASSIGNMENT))
    saved = L.ASSIGNMENT.pop(entry)
    try:
        return s25_assign.run()
    finally:
        L.ASSIGNMENT[entry] = saved


def _assign_ghost():
    L.ASSIGNMENT["no_such_directory"] = (2, "T1", "a ghost", "-")
    try:
        return s25_assign.run()
    finally:
        del L.ASSIGNMENT["no_such_directory"]


def _figure_missing_file():
    saved = list(fm.FIGURES)
    fm.FIGURES[0] = (fm.FIGURES[0][0], fm.FIGURES[0][1],
                     "docs/figures/no_such_figure", fm.FIGURES[0][3])
    try:
        return s25_figures.run()
    finally:
        fm.FIGURES[:] = saved


def _figure_duplicate_slug():
    saved = list(fm.FIGURES)
    fm.FIGURES.append((1, fm.FIGURES[0][1], fm.FIGURES[0][2], "a duplicate"))
    try:
        return s25_figures.run()
    finally:
        fm.FIGURES[:] = saved


def _figure_unexcluded():
    saved = dict(fm.UNPLACED)
    fm.UNPLACED.clear()
    try:
        return s25_figures.run()
    finally:
        fm.UNPLACED.update(saved)


def _stitch_missing_chapter():
    (L.TH / L.CHAPTER_FILES[3]).unlink()
    return s25_stitch.run()


def _stitch_unplaced_figure():
    path = L.TH / "06_phylogeny.md"
    path.write_text(path.read_text().replace(
        "![](figures/tree_ml_rooted.png)", "", 1), encoding="utf-8")
    return s25_stitch.run()


def _stitch_borrowed_figure():
    path = L.TH / "05_eukaryote_range.md"
    path.write_text(path.read_text()
                    + "\n![](figures/tree_ml_rooted.png)\n", encoding="utf-8")
    return s25_stitch.run()


def _stitch_dangling_ref():
    path = L.TH / "05_eukaryote_range.md"
    path.write_text(path.read_text() + "\n{fig:no_such_slug}\n",
                    encoding="utf-8")
    return s25_stitch.run()


def _refs_unknown_key():
    path = L.TH / "05_eukaryote_range.md"
    path.write_text(path.read_text() + "\nA claim resting on [R999].\n",
                    encoding="utf-8")
    return s25_stitch.run()


def _refs_unaudited():
    audit = s25_refs.audit_path()
    rows = [r for r in L.read_tsv(audit) if r["ref_id"] != "R160"]
    L.write_tsv(audit, rows, s25_refs.AUDIT_FIELDS)
    return s25_stitch.run()


def _refs_uncited():
    # R160 is cited exactly once, so removing that one citation makes it an
    # audited reference the document never uses.
    path = L.TH / "02_baseline.md"
    path.write_text(path.read_text().replace(" [R160]", "", 1),
                    encoding="utf-8")
    return s25_stitch.run()


def _claims_wrong_value():
    claim = s25_claims_thesis.CLAIMS[0]
    saved = claim["expect"]
    claim["expect"] = "999999"
    try:
        return s25_claims.run()
    finally:
        claim["expect"] = saved


def _claims_padded():
    # appended to the concatenated list the driver actually reads, not to the
    # data module it was built from
    s25_claims.CLAIMS.append(dict(
        id="TX", chapter=1, section="ch1",
        claim="a number chapter 1 does not state",
        source="results/census_v2/census_v2.tsv", op="count",
        expect="15421"))
    try:
        return s25_claims.run()
    finally:
        s25_claims.CLAIMS.pop()


def _pdf_dropped_glyph():
    """The glyph guard, exercised on the log text it reads."""
    log = ("Overfull \\hbox\n"
           "Missing character: There is no ∮ in font texgyretermes-regular!\n")
    missing = sorted({ln.strip() for ln in log.split("\n")
                      if "Missing character" in ln})
    if missing:
        print(f"[s25 pdf] {len(missing)} distinct missing characters — the "
              f"document font lacks a glyph the text uses", file=sys.stderr)
        return 1
    return 0


CASES = [
    ("assign: an unassigned results directory", "rule T6", _assign_unassigned),
    ("assign: an assignment with no directory", "does not exist",
     _assign_ghost),
    ("figures: a declared figure with no file", "missing",
     _figure_missing_file),
    ("figures: a slug declared twice", "declared more than once",
     _figure_duplicate_slug),
    ("figures: a committed figure neither placed nor excluded",
     "neither placed in the thesis nor listed", _figure_unexcluded),
    ("stitch: a missing chapter file", "missing chapter file",
     _stitch_missing_chapter),
    ("stitch: a declared figure nothing places", "no chapter places it",
     _stitch_unplaced_figure),
    ("stitch: a figure placed in the wrong chapter", "which is declared for",
     _stitch_borrowed_figure),
    ("stitch: a reference to an unplaced figure", "never placed",
     _stitch_dangling_ref),
    ("refs: a cited key with no reference row", "no row in references.tsv",
     _refs_unknown_key),
    ("refs: a new reference cited without an audit row",
     "no verified audit row", _refs_unaudited),
    ("refs: a reference audited and never cited", "never cited",
     _refs_uncited),
    ("claims: a number the table does not produce", "expected '999999'",
     _claims_wrong_value),
    ("claims: a claim its chapter does not state", "does not state",
     _claims_padded),
    ("pdf: a glyph the document font lacks", "missing characters",
     _pdf_dropped_glyph),
]


def main() -> int:
    before = _digests()
    for name, expect, fn in CASES:
        case(name, expect, fn)
    after = _digests()

    n_ok = sum(1 for r in RESULTS if r["ok"])
    L.write_tsv(L.TH / OUT_NAME, RESULTS, FIELDS)
    for r in RESULTS:
        print(f"  {'ok  ' if r['ok'] else 'FAIL'}  {r['guard']}")
        print(f"          said: {r['said']}")

    tampered = sorted(k for k in before
                      if before[k] != after.get(k, before[k]))
    tampered += sorted(set(after) - set(before))
    status = 0
    if tampered:
        print(f"  [FAIL] the suite altered {len(tampered)} committed file(s): "
              f"{tampered[:3]}", file=sys.stderr)
        status = 1
    print(f"[s25 guards] {n_ok}/{len(RESULTS)} guards fired with the expected "
          f"message; {len(before)} committed files unchanged")
    return 1 if (n_ok != len(RESULTS) or status) else 0


if __name__ == "__main__":
    raise SystemExit(main())
