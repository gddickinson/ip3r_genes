"""S25 stage `refs` — audit every reference added for the thesis, then use it.

Two jobs, in one module because they are two halves of one rule.

**The audit.** Every key cited by the thesis that is not in
`thesis/reference_baseline.txt` — the 137 keys S0 assembled and the review
used — is a new reference, and a new reference is worthless unless it has been
checked. So each is declared in `s25_refs_new.py` by identifier only, resolved
live against PubMed (esummary) or Crossref, and admitted **only if the title
the identifier resolves to carries the distinctive phrase the declaration
says it should**. The bibliographic row is then written from the fetched
record, so nothing about it is typed. A failed lookup, or a title that does
not carry its phrase, is a build error: the reference is not added and
`s25_assemble.py` exits non-zero.

That check catches the failure this rule exists to prevent. A transposed
PubMed ID resolves perfectly well — to a different paper — and a bibliography
assembled by hand has no way to notice.

Responses are archived under the data root, so a re-run needs no network.

**The use.** The stitch stage resolves this project's stable keys (`[R22]`,
`[R138, R140]`) into numbers in order of first appearance and renders the
bibliography from the table, exactly as `s14_refs.py` and
`s0_review_build.py` do — a cited key with no reference row is a build error
in all three.

  python scripts/s25_refs.py --audit          # fetch, check, extend the table
  python scripts/s25_refs.py --audit --offline  # re-check from the archive
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import s25_lib as L
import s25_refs_new as NEW

REFS_TSV = L.RESULTS / "s0_baseline" / "references.tsv"


def baseline_path():
    """Resolved at call time: `L.TH` is redirected by the guard suite."""
    return L.TH / "reference_baseline.txt"


def audit_path():
    return L.TH / "reference_audit.tsv"
REF_FIELDS = ["ref_id", "pmid", "doi", "year", "journal", "authors", "title",
              "role"]
AUDIT_FIELDS = ["ref_id", "kind", "source_db", "query_id", "expect_in_title",
                "fetched_title", "fetched_year", "fetched_journal",
                "fetched_first_author", "verdict", "used_for", "checked"]

EMAIL = "george.dickinson@gmail.com"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CROSSREF = "https://api.crossref.org/works/"

CITE = re.compile(r"\[(R\d{2,3}(?:\s*,\s*R\d{2,3})*)\]")
MARKER = "{references}"


def archive_dir() -> Path:
    try:
        sys.path.insert(0, str(L.ROOT))
        from src.utils.data_root import get_data_root
        d = Path(get_data_root()) / "raw_api" / "s25"
    except Exception:                                        # noqa: BLE001
        d = L.ROOT / "cache" / "s25_refs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _norm(s: str) -> str:
    """Lowercase, strip punctuation, collapse whitespace."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", s.lower())).strip()


# --------------------------------------------------------------- fetching

def fetch_pubmed(pmids: list[str], offline: bool) -> dict[str, dict]:
    """One esummary call for every PubMed ID, archived."""
    if not pmids:
        return {}
    cache = archive_dir() / "pubmed_esummary.json"
    if offline or cache.exists():
        raw = json.loads(cache.read_text(encoding="utf-8"))
        if all(p in raw.get("result", {}) for p in pmids):
            return raw["result"]
        if offline:
            missing = [p for p in pmids if p not in raw.get("result", {})]
            raise SystemExit(f"--offline but the archive lacks {missing}")
    url = (ESUMMARY + "?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "retmode": "json",
         "email": EMAIL, "tool": "ip3r_genes_s25"}))
    with urllib.request.urlopen(url, timeout=90) as fh:
        raw = json.loads(fh.read().decode("utf-8"))
    cache.write_text(json.dumps(raw, indent=1), encoding="utf-8")
    return raw["result"]


def fetch_crossref(doi: str, offline: bool) -> dict:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", doi)
    cache = archive_dir() / f"crossref_{slug}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))["message"]
    if offline:
        raise SystemExit(f"--offline but no archived Crossref record for {doi}")
    req = urllib.request.Request(
        CROSSREF + urllib.parse.quote(doi, safe="/"),
        headers={"User-Agent": f"ip3r_genes/s25 (mailto:{EMAIL})"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        raw = json.loads(fh.read().decode("utf-8"))
    cache.write_text(json.dumps(raw, indent=1), encoding="utf-8")
    time.sleep(0.4)
    return raw["message"]


# ---------------------------------------------------------------- parsing

def _from_pubmed(rec: dict) -> dict:
    authors = [a["name"] for a in rec.get("authors", [])
               if a.get("authtype") == "Author"]
    doi = ""
    for aid in rec.get("articleids", []):
        if aid.get("idtype") == "doi":
            doi = aid.get("value", "")
    year = (rec.get("pubdate") or "")[:4]
    return {"title": (rec.get("title") or "").rstrip("."),
            "journal": rec.get("source", ""), "year": year,
            "authors": ", ".join(authors), "doi": doi,
            "pmid": rec.get("uid", "")}


def _from_crossref(msg: dict) -> dict:
    authors = ", ".join(
        f"{a.get('family', '')} {''.join(p[0] for p in a.get('given', '').split() if p)}".strip()
        for a in msg.get("author", []) if a.get("family"))
    title = (msg.get("title") or [""])[0]
    journal = (msg.get("container-title") or [""])[0] or msg.get("publisher", "")
    parts = (msg.get("issued", {}).get("date-parts") or [[""]])[0]
    return {"title": title.rstrip("."), "journal": journal,
            "year": str(parts[0]) if parts else "",
            "authors": authors, "doi": msg.get("DOI", ""), "pmid": ""}


# ------------------------------------------------------------------ audit

def load_baseline() -> set[str]:
    return {ln.strip()
            for ln in baseline_path().read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.startswith("#")}


def audit(offline: bool = False, verbose: bool = False) -> int:
    pmids = [v["pmid"] for v in NEW.NEW_REFS.values() if v.get("pmid")]
    pubmed = fetch_pubmed(pmids, offline)

    existing = {r["ref_id"]: r for r in L.read_tsv(REFS_TSV)}
    baseline = load_baseline()
    rows, n_fail, added = [], 0, 0

    for ref_id, dec in sorted(NEW.NEW_REFS.items()):
        if ref_id in baseline:
            print(f"  [FAIL] {ref_id} is declared as a new reference but is "
                  f"in the frozen baseline", file=sys.stderr)
            n_fail += 1
            continue
        if dec.get("pmid"):
            src, qid = "pubmed", dec["pmid"]
            rec = pubmed.get(qid, {})
            meta = _from_pubmed(rec) if rec and "error" not in rec else None
        else:
            src, qid = "crossref", dec["doi"]
            try:
                meta = _from_crossref(fetch_crossref(qid, offline))
            except Exception as exc:                          # noqa: BLE001
                print(f"  {ref_id}: crossref {qid}: {exc}", file=sys.stderr)
                meta = None

        if meta is None or not meta["title"]:
            verdict = "failed: identifier did not resolve"
            meta = {"title": "", "journal": "", "year": "", "authors": "",
                    "doi": dec.get("doi", ""), "pmid": dec.get("pmid", "")}
        elif _norm(dec["expect"]) not in _norm(meta["title"]):
            verdict = "failed: title does not carry the declared phrase"
        else:
            verdict = "verified"

        if verdict != "verified":
            n_fail += 1
            print(f"  [FAIL] {ref_id} ({src} {qid}): {verdict}\n"
                  f"         expected {dec['expect']!r}\n"
                  f"         got      {meta['title']!r}", file=sys.stderr)
        else:
            role = {"tool": "software", "method": "method", "db": "database",
                    "lit": "primary"}[dec["kind"]]
            existing[ref_id] = {
                "ref_id": ref_id, "pmid": meta["pmid"], "doi": meta["doi"],
                "year": meta["year"], "journal": meta["journal"],
                "authors": meta["authors"], "title": meta["title"],
                "role": role}
            added += 1
            if verbose:
                print(f"  [ok  ] {ref_id} {meta['year']} {meta['journal']}: "
                      f"{meta['title'][:70]}")

        rows.append({
            "ref_id": ref_id, "kind": dec["kind"], "source_db": src,
            "query_id": qid, "expect_in_title": dec["expect"],
            "fetched_title": meta["title"], "fetched_year": meta["year"],
            "fetched_journal": meta["journal"],
            "fetched_first_author": meta["authors"].split(",")[0].strip(),
            "verdict": verdict, "used_for": dec["used_for"],
            "checked": date.today().isoformat()})

    L.write_tsv(audit_path(), rows, AUDIT_FIELDS)
    if not n_fail:
        _write_refs([existing[k]
                     for k in sorted(existing, key=lambda k: int(k[1:]))])
    print(f"[s25 refs] {len(rows) - n_fail}/{len(rows)} new references "
          f"verified against a live record, {added} written into "
          f"references.tsv ({len(existing)} total) "
          f"-> thesis/reference_audit.tsv")
    return 1 if n_fail else 0


def _write_refs(rows: list[dict]) -> None:
    """Write `references.tsv` with the line endings it has always had.

    This table is shared with the literature review and S0's report and has
    been LF-terminated since S0 wrote it. `csv.writer` defaults to CRLF, so
    writing it through the generic helper would rewrite all 137 inherited rows
    as changed while changing nothing in them — and the point of adding 58
    references to a shared table is that the diff shows 58 added rows.
    """
    with REFS_TSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REF_FIELDS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------ citation resolution

def load_refs() -> dict[str, dict[str, str]]:
    return {r["ref_id"]: r for r in L.read_tsv(REFS_TSV)}


def _format(num: int, row: dict[str, str]) -> str:
    parts = [a.strip() for a in row["authors"].split(",") if a.strip()]
    authors = (", ".join(parts[:6]) + " *et al.*" if len(parts) > 6
               else ", ".join(parts))
    stop = "" if authors.endswith("*") else "."
    line = (f"{num}. {authors}{stop} {row['title'].rstrip('.')}. "
            f"*{row['journal']}* **{row['year']}**.")
    if row.get("pmid"):
        line += f" PMID {row['pmid']}."
    if row.get("doi") and row["doi"] != "-":
        line += f" doi:{row['doi']}."
    return line


def resolve(text: str) -> tuple[str, str, dict]:
    """Renumber citations in order of first appearance; render the list.

    Raises KeyError naming every cited key with no reference row, and every
    cited key that is new and carries no verified audit row.
    """
    refs = load_refs()
    baseline = load_baseline()
    audit = audit_path()
    audited = {r["ref_id"] for r in L.read_tsv(audit)
               if r["verdict"] == "verified"} if audit.exists() else set()
    order: list[str] = []
    missing: list[str] = []
    unaudited: list[str] = []

    def repl(m: re.Match) -> str:
        keys = [k.strip() for k in m.group(1).split(",")]
        for k in keys:
            if k not in refs:
                missing.append(k)
            elif k not in baseline and k not in audited:
                unaudited.append(k)
            if k not in order:
                order.append(k)
        return "[" + ",".join(str(order.index(k) + 1) for k in keys) + "]"

    out = CITE.sub(repl, text)
    problems = []
    if missing:
        problems.append("cited keys with no row in references.tsv: "
                        + ", ".join(sorted(set(missing))))
    if unaudited:
        problems.append("cited keys added after the frozen baseline with no "
                        "verified audit row: " + ", ".join(sorted(set(unaudited))))
    # The audit rule runs in both directions. A reference added and not
    # audited launders an assumption into a citation; a reference audited and
    # never cited is decoration, and a bibliography that grows with decoration
    # is the thing this rule exists to prevent, arrived at from the other
    # side. Both are build errors.
    uncited = sorted(audited - set(order), key=lambda k: int(k[1:]))
    if uncited:
        problems.append("references added for the thesis and never cited: "
                        + ", ".join(uncited))
    if problems:
        raise KeyError("; ".join(problems))

    bib = "\n\n".join(_format(i + 1, refs[k]) for i, k in enumerate(order))
    return out, bib, {"n_cited": len(order), "n_available": len(refs),
                      "n_new": len([k for k in order if k not in baseline])}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit", action="store_true")
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    raise SystemExit(audit(offline=a.offline, verbose=a.verbose)
                     if a.audit else 0)
