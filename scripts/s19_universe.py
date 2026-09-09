"""The swept accession universe: which records each database actually held.

Split out of `s19_lib` to keep both inside the 500-line budget, and because
it is the one part of S19 that touches bulk data: a single `grep '^>'` pass
over each of the seven reference-proteome FASTAs, 42 GB in total, cached
under `<data_root>/raw_api/s19/`.

Why it exists at all. Without the universe, "the InterPro enumeration holds a
record the profile HMM did not return" cannot be told apart from "that record
was never in the database the profile HMM searched", and the head-to-head
comparison would be an accounting artefact of two different search spaces
rather than a comparison of two methods.

`s19_lib` re-exports every name here, so nothing else has to know about the
split.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from s19_lib_base import (ALL_GROUPS, acc_key, cache_dir,  # noqa: F401
                          data_root, log)

_OX_RE = re.compile(r"OX=(\d+)")


def header_accession(header: str) -> str:
    """UniProt accession out of a FASTA header line.

    Reference-proteome headers are `>db|ACC|ENTRY_NAME description`; a few
    carry a bare id. Taking the middle pipe field (else the first token)
    covers both without a regex that can silently match nothing.
    """
    tok = header[1:].split(None, 1)[0] if len(header) > 1 else ""
    if "|" in tok:
        parts = tok.split("|")
        tok = parts[1] if len(parts) >= 3 else parts[-1]
    return acc_key(tok)


def _universe_path(group: str) -> Path:
    return cache_dir() / f"universe_{group}.txt"


def build_universe(group: str, force: bool = False) -> Path:
    """Cache the accession list of one concatenated reference-proteome DB.

    A single `grep '^>'` pass; the seven DBs run 0.7-10 GB each, so this is
    the one place S19 touches bulk data. A zero-row parse **refuses to
    cache** rather than writing an empty universe that would make every
    later membership test read "not in the database".
    """
    out = _universe_path(group)
    if out.exists() and out.stat().st_size and not force:
        return out
    fasta = data_root() / "proteomes" / f"{group}_refprot.fasta"
    if not fasta.exists():
        raise SystemExit(f"proteome DB missing: {fasta}")
    tmp = out.with_suffix(".partial")
    log(f"indexing {fasta.name} ({fasta.stat().st_size / 1e9:.1f} GB)")
    n = 0
    with open(tmp, "w") as fh:
        grep = subprocess.Popen(["grep", "^>", str(fasta)], text=True,
                                stdout=subprocess.PIPE)
        for line in grep.stdout:                             # type: ignore
            acc = header_accession(line.rstrip("\n"))
            if not acc:
                continue
            ox = _OX_RE.search(line)
            fh.write(f"{acc}\t{ox.group(1) if ox else ''}\n")
            n += 1
        grep.wait()
    if not n:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"{group}: parsed 0 headers from {fasta} — refusing "
                         "to cache an empty universe")
    tmp.replace(out)
    log(f"{group}: {n:,} accessions -> {out.name}")
    return out


def _universe_ready(group: str) -> Path:
    path = _universe_path(group)
    if not path.exists() or not path.stat().st_size:
        build_universe(group, force=True)
    return path


def universe_intersect(group: str, query: set[str]) -> set[str]:
    """Which of `query` is in that database — streamed, never held in memory."""
    found = set()
    with open(_universe_ready(group)) as fh:
        for line in fh:
            acc = line.split("\t", 1)[0]
            if acc in query:
                found.add(acc)
    return found


def universe_size(group: str) -> int:
    with open(_universe_ready(group)) as fh:
        return sum(1 for _ in fh)


def universe_taxids(group: str) -> set[str]:
    out = set()
    with open(_universe_ready(group)) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 1 and parts[1]:
                out.add(parts[1])
    return out


