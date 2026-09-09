"""S25 stage `prose` — the prose rules that are checked mechanically.

D78 says a style rule is only applied when it is checked mechanically, and
D81 says an appended sentence is not an edit. Both produced defects that no
other guard in this project could see, so the three cheapest checks now run on
every build.

  1. **No em-dash in a chapter source.** D76 removed 464 of them and replaced
     each with the punctuation or the sentence break the sense wanted. An
     en-dash is legitimate in a numeric range and is not counted.
  2. **A figure legend may not say the same thing twice.** Every pair of
     sentences inside a legend is compared on content words and on six-word
     runs. D80 asked every legend to explain why its figure matters, and
     appending that sentence rather than merging it left 50 of 103 legends
     stating one point twice.
  3. **A legend may not restate the paragraph beside it.** Every legend is
     compared against the two paragraphs either side. A legend has to stand
     alone, and paraphrasing its own neighbour is how it stops being worth
     reading.

Thresholds are deliberately loose, because the cost of a false positive is a
sentence somebody rereads and the cost of a false negative is a paragraph that
wastes a reader's time. They were tuned once, against the 103 legends this
document has, to the point where every remaining flag was a real repetition.

  python scripts/s25_prose.py            # check, print a summary
  python scripts/s25_prose.py --verbose  # print every comparison that fired
"""

from __future__ import annotations

import argparse
import re
import sys

import s25_lib as L

#: Words carrying no topic, excluded before two sentences are compared.
STOP = set(
    "the a an and or of to in is are that this it its for as on by with be "
    "been was were which what not no than then so from at each every one two "
    "three but has have had can cannot could would does do their they them "
    "there here how why when where because rather into over under also only "
    "more most less least any all both same other".split())

#: Shared content words that make two sentences in one legend a repetition.
LEGEND_CONTENT = 6
#: Shared six-word runs that make two sentences in one legend a repetition.
LEGEND_RUNS = 2
#: Shared six-word runs that make a legend a restatement of its neighbour.
NEIGHBOUR_RUNS = 4

LEGEND_RE = re.compile(r"\*\*\{fig:([a-z0-9_]+)\}\.\*\*(.*)", re.S)


def _content(s: str) -> set[str]:
    return {w for w in re.findall(r"[a-z]{3,}", s.lower()) if w not in STOP}


def _runs(s: str, n: int = 6) -> set[str]:
    w = re.findall(r"[a-z0-9]+", s.lower())
    return {" ".join(w[i:i + n]) for i in range(len(w) - n + 1)}


def _sentences(text: str) -> list[str]:
    return re.split(r"(?<=[.!?]) +", " ".join(text.split()))


def _is_prose(para: str) -> bool:
    lines = para.split("\n")
    return not any(l.startswith(("#", "|", "```", "![", ">", "    "))
                   for l in lines) and not para.lstrip().startswith("**{fig:")


def _check_file(path, verbose: bool) -> list[str]:
    text = path.read_text(encoding="utf-8")
    fails = []

    n_em = text.count("—")
    if n_em:
        fails.append(f"{path.name}: {n_em} em-dash(es) in a chapter source")

    paras = text.split("\n\n")
    for i, para in enumerate(paras):
        m = LEGEND_RE.match(para.lstrip())
        if not m:
            continue
        slug, body = m.group(1), m.group(2)
        sents = _sentences(body)

        for a in range(len(sents)):
            for b in range(a + 1, len(sents)):
                ca, cb = _content(sents[a]), _content(sents[b])
                shared = ca & cb
                runs = _runs(sents[a]) & _runs(sents[b])
                if (len(ca) >= 4 and len(cb) >= 4
                        and len(shared) >= LEGEND_CONTENT) or \
                        len(runs) >= LEGEND_RUNS:
                    fails.append(
                        f"{path.name}: legend {slug} says the same thing twice "
                        f"({len(shared)} shared words, {len(runs)} shared runs)")
                    if verbose:
                        print(f"    A: {sents[a][:100]}")
                        print(f"    B: {sents[b][:100]}")

        legend_runs = _runs(para)
        for j in (i - 2, i - 1, i + 1, i + 2):
            if not 0 <= j < len(paras) or not _is_prose(paras[j]):
                continue
            shared = legend_runs & _runs(paras[j])
            if len(shared) >= NEIGHBOUR_RUNS:
                fails.append(
                    f"{path.name}: legend {slug} restates the paragraph "
                    f"{'before' if j < i else 'after'} it "
                    f"({len(shared)} shared runs)")
                if verbose:
                    print(f"    shared: {sorted(shared)[0]}")
    return fails


def run(verbose: bool = False) -> int:
    files = sorted(L.TH.glob("[0-9]*.md"))
    fails: list[str] = []
    n_legends = 0
    for path in files:
        n_legends += sum(
            1 for para in path.read_text(encoding="utf-8").split("\n\n")
            if LEGEND_RE.match(para.lstrip()))
        fails += _check_file(path, verbose)
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    print(f"[s25 prose] {len(files)} chapter sources, {n_legends} legends, "
          f"no em-dash and no repetition, {len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    sys.exit(run(verbose=ap.parse_args().verbose))
