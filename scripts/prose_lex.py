"""A finite-verb detector for the prose rules, stdlib only.

D76 rewrote headings that were noun phrases, D78 recorded that a rule not
checked mechanically is a rule that will be broken again, and S28 found the
same fault surviving in the two places neither pass had looked: the first
sentence of a figure legend and the text drawn inside a figure. This module is
the one detector both places now share (`s25_prose.py` for the chapter
sources, `figcheck.py` for the drawn figures), so the two cannot disagree
about what a sentence is.

**What it can see.** A sentence is taken to carry a finite verb when, after
its parenthetical, bracketed and relative-clause material is removed, it
contains one of: an auxiliary or modal (*is, has, can, must*); an irregular
past tense that is not also a participle (*grew, found, drew*, but not
*drawn*, *shown*, *given*); a third-person `-s` form of a verb in the lexicon
that is not sitting in noun position (not after a determiner, number,
preposition or adjective); a regular `-ed` form not introduced by *as*,
*once*, *when* or a comma and not followed by *rather than* or *as*; or a
base form after a plural noun or pronoun.

**What it cannot see, stated so the reading pass knows where to look.** A
reduced relative clause on a regular verb (*The margin measured across the
census*) reads as a sentence, because *measured* is both a past tense and a
participle and only the meaning tells them apart. A `that`-clause is not
stripped, because *that* is a determiner as often as a relativiser. A
compound noun whose head is a lexicon verb (*the record counts*) can pass.
Every one of these is a false negative, never a false positive, which is the
direction to fail in: the check is a floor under the reading pass, not a
replacement for it.

The lexicon is deliberately generous with verbs and deliberately excludes
the words that are nouns more often than verbs in this document (*result,
count, set, read, split*), because a noun in the lexicon lets a fragment
through and a verb missing from it only asks for a synonym.
"""

from __future__ import annotations

import re

# ------------------------------------------------------------------ lexicon

#: Forms that are finite on their own wherever they appear.
AUX = set("""
am is are was were has have had do does did can cannot could will would
shall should may might must
""".split())

#: Irregular past tenses that are not also past participles. `set`, `put`,
#: `cut`, `read`, `split`, `spread`, `hit`, `let`, `shut`, `cost` are left out
#: on purpose: each is also a base form and a common noun here, and the
#: base-form rule below decides them instead.
PAST = set("""
arose became began bent bit bled blew broke brought built bought came caught
chose clung crept dealt drank drew drove ate fell fed felt fought found fled
flew forbade forgot forgave foresaw froze got gave went grew hung heard hid
held kept knew laid led left lent lay lit lost made meant met outgrew overcame
overrode overtook paid quit rang rode rose ran said saw sold sent shook shed
shone shot showed shrank sang sank sat slid slept spoke spent sprang stood
stole struck strove stuck stung swam swept swung took taught tore told thought
threw understood undertook underwent woke wore won withdrew withheld wound
wrote
""".split())

#: Base forms. Inflected as -s / -es / -ies, -ed / -d / -ied, or used bare
#: after a plural subject.
VERBS = set("""
accept account achieve acquire act add address adjust admit adopt advance
affect agree aim align allow alter amount analyse analyze anchor announce
answer appear apply approach argue arise arrive ask assemble assert assess
assign assume attach attempt attribute audit average avoid bear beat become
begin behave belong bend benefit bias bind block blow borrow bound break bring
build bury buy calibrate call cap capture carry catch cause cease change
characterise characterize charge check choose claim classify clear climb close
cluster code collapse collect colour color combine come commit compare
compensate compete compile complete comprise compute conceal concern conclude
confirm conflict confound connect consider consist constitute constrain
construct consume contain continue contradict contrast contribute control
converge convert copy correct correlate corroborate cost count cover create
credit cross cut damage date deal decide decline decrease deduce default
defend define delay delete deliver demand demonstrate deny depend deposit
derive describe deserve design detect determine develop deviate differ dilute
diminish direct disagree disambiguate discard discount discover discriminate
discuss dismiss display dispute dissolve distinguish distort distribute
diverge divide dominate double doubt download draw drift drive drop duplicate
earn echo edit eliminate emerge emit employ enable encode encounter end
enforce enrol enroll ensure enter enumerate equal escape establish estimate
evaluate evolve exceed exclude execute exercise exhaust exist expand expect
explain exploit explore export expose express extend extract fail fall favour
favor feed feel fetch fill filter find finish fire fit fix flag flank flow fold
follow forbid force forget form found frame free freeze gain gate gather
generate get give go govern grade grant group grow guard guess halve hand
handle hang happen harvest head hear help hide highlight hinge hit hold hope
identify ignore illustrate imply import improve include incorporate increase
indicate infer inflate inherit insert inspect install integrate intend
interact interpolate interpret intersect introduce invalidate invent invert
investigate invite involve isolate join judge jump justify keep kill know label
land last launch lay lead learn leave lend let lie lift limit line link list
live load localise localize locate lock look lose maintain make manage
manufacture map mark mask match matter mean measure meet mention merge migrate
mimic mine mirror misread miss mistake model modify motivate move multiply name
need nest normalise normalize note notice number obey observe obtain occupy
occur offer omit open operate order organise organize outnumber outrank
outperform outweigh overlap overrule overstate overturn owe own pack pad paint
pair parse partition pass pay penalise penalize perform permit persist pick
pin place plan plot point pool populate pose position possess precede predict
prefer prepare present preserve press presume pretend prevent print probe
proceed process produce project promote prompt propose protect prove provide
prune publish pull push put qualify quantify question quote raise range rank
rate reach read realign realise realize reassemble reassign rebuild recall
receive recode recognise recognize recommend recompute reconcile reconstruct
record recover recur reduce refer refine reflect refuse regard register reject
relabel relate relax release rely remain remove render renumber repair repeat
replace replicate report represent reproduce request require rerun rescue
resemble reserve reset reside resist resolve respect respond rest restate
restore restrict resume retain retire retrieve return reveal reverse review
revise rewrite rid ride rise risk root rotate round rule run sample satisfy
save say scale scan score screen search see seed seek seem select sell send
separate serve set settle shade shape share shift shorten show shrink shut
sign signal simulate sit skew skip slip slow smooth solve sort span speak
specify spend split spread stack stand start state stay stem step stick stop
store strengthen stress stretch strike strip structure struggle study submit
substitute subtract succeed suffer suffice suggest suit sum summarise
summarize supply support suppose suppress surface surpass survive suspect
sweep swap switch tag take talk target teach tell tend terminate test thank
think threaten throw tie tile tilt tolerate top total touch trace track train
transcribe transfer transform translate transmit trap travel treat trim trip
trust try tune turn type underestimate underlie undermine understand undo
unfold unify unite unpack update use validate value vanish vary verify veto
view visit vote wait walk want warn wash waste watch weaken wear weigh weight
widen win wish withdraw withhold witness work worry wrap write yield
understate overstate contest collide flip splice coincide stratify resample
sustain accumulate port happen encode relax travel declare
reweight refit rescore reindex regroup recount reread rebin underpin outvote
separate saturate corroborate exercise license veto walk
""".split())

#: Determiners, numbers and prepositions: an `-s` word after one of these is
#: a plural noun, and an `-ed` word after one is an adjective.
DETS = set("""
the a an this that these those its their our his her my your each every all
both some any no most more fewer several many few other such same another
one two three four five six seven eight nine ten eleven twelve fifteen twenty
thirty fifty hundred thousand first second third last next own single whole
""".split())
PREPS = set("""
of in on at by for with from to into onto over under across against between
among through during without within beside above below after before than as
per via along inside outside toward towards behind beyond despite like
""".split())
CONJ = set("and or nor but not only".split())
NOUN_POSITION = DETS | PREPS | CONJ | {"whose"}

#: Words that open a clause whose verb is not the main clause's. Relative and
#: interrogative words, and the subordinating conjunctions this document uses.
CLAUSE_OPENERS = set("""
which who whom whose where what how why when because although though
whereas while until unless if once so before after
""".split())
#: The subordinating conjunctions among them: their clause runs to the next
#: comma. A relative or interrogative word's clause ends at its own verb.
SUBORDINATORS = set("because although though whereas while until unless if "
                    "once so before after".split())
#: A word after which a determiner does not open a contact clause.
ADVERBS = set("so then thus hence here there now still also therefore however "
              "yet even only just never always often already not".split())
#: Of those, the ones that can open a sentence as a free relative (*What
#: forces the placement is one pair*); such a sentence needs a second verb.
FREE_RELATIVE = set("which who whom whose where what how why when".split())

#: Determiners that, directly after a noun, open a contact relative clause
#: (*the columns the tree saw*, *the object this thesis counts*).
CONTACT_DETS = set("the a an this these those its their each every no some "
                   "any".split())

#: A token that stands for a plural subject when a bare base form follows.
PLURAL_PRONOUNS = set("they we you i these those both all some none several "
                      "many few most".split())
IRREGULAR_PLURALS = set("loci taxa genera data criteria analyses hypotheses "
                        "indices matrices phenomena bacteria fungi children "
                        "people".split())

QUANTIFIERS = set("two three four five six seven eight nine ten eleven twelve "
                  "fifteen twenty thirty fifty hundred thousand all both some "
                  "many most several few fourteen sixteen eighteen".split())


#: `-ed` after one of these is a participle or an adjective, not a past tense.
PARTICIPLE_CUES = set("as once when if unless than being having by of and or "
                      ",".split()) | (DETS - QUANTIFIERS)
#: `-ed` before one of these is a reduced relative (*the outcome measured on
#: the finished model*) far more often than a past tense in this document.
PARTICIPLE_NEXT = PREPS | {"rather", "twice", "separately", "together"}

_ADJ_SUFFIX = ("ive", "ous", "ful", "less", "wise")


# ---------------------------------------------------------------- cleaning

def clean(text: str) -> str:
    """Strip markup, parentheticals, citations and figure references."""
    t = re.sub(r"`[^`]*`", " ", text)
    t = re.sub(r"\{fig:[a-z0-9_]+\}", "Figure", t)
    t = re.sub(r"\([^()]*\)", " ", t)
    t = re.sub(r"\[[^\[\]]*\]", " ", t)
    t = re.sub(r"\$[^$]*\$", " ", t)
    t = t.replace("**", " ").replace("*", " ")
    t = re.sub(r"[₀-₉⁰-⁹]", "", t)
    return " ".join(t.split())


def first_sentence(text: str) -> str:
    s = re.split(r"(?<=[.!?]) +(?=[A-Z\"'(*])", clean(text).strip(), maxsplit=1)
    return s[0].strip() if s else ""


def _tokens(sentence: str) -> list[str]:
    return re.findall(r"[A-Za-z][A-Za-z'\-]*|\d[\w.,]*\d|\d|[,;:]", sentence)


def _word(tok: str) -> str:
    return re.sub(r"[^a-z']", "", tok.lower())


# --------------------------------------------------------------- detection

def _stems_s(w: str) -> set[str]:
    out = {w[:-1]}
    if w.endswith("es"):
        out.add(w[:-2])
    if w.endswith("ies"):
        out.add(w[:-3] + "y")
    return out


def _stems_ed(w: str) -> set[str]:
    out = {w[:-2], w[:-1]}
    if w.endswith("ied"):
        out.add(w[:-3] + "y")
    if len(w) > 4 and w[-3] == w[-4]:
        out.add(w[:-3])
    return out


def _is_number(tok: str) -> bool:
    return bool(tok) and tok[0].isdigit()


def _plural(tok: str) -> bool:
    w = _word(tok)
    if w in PLURAL_PRONOUNS or w in IRREGULAR_PLURALS or w in QUANTIFIERS:
        return True
    if any(c.isdigit() for c in tok) or (tok.isupper() and len(tok) > 1):
        return True                       # ITPR3, RyR, PF08709, 309: a name
    return (w.endswith("s") and not w.endswith(("ss", "us", "sis", "'s"))
            and len(w) > 3 and w not in NOUN_POSITION and w not in AUX
            and w not in ("thus", "this"))


def _noun_position(toks: list[str], i: int) -> bool:
    """Is the `-s` word at `i` sitting where a noun sits?"""
    if i == 0:
        return True
    prev = toks[i - 1]
    p = _word(prev)
    if p == "that" or (i == 1 and p in ("this", "that", "these", "those")):
        return False                      # the rule that decides; This changes
    if p in NOUN_POSITION or p.endswith(("'s", "s'")) \
            or prev in (",", ";", ":"):
        return True
    if _is_number(prev):
        return not (i >= 2 and toks[i - 2][:1].isupper())   # Chapter 3 builds
    if "-" in prev:
        return True                       # per-site rates
    if i >= 2 and len(p) > 5 and ((p.endswith("ed") and _stems_ed(p) & VERBS)
                                  or (p.endswith("ing") and p[:-3] in VERBS)):
        return True                       # the measured counts
    return i >= 2 and p.endswith(_ADJ_SUFFIX) and p not in AUX


def _plural_subject(toks: list[str], i: int) -> bool:
    """Walk back from a base form to its subject, skipping prepositional
    phrases, and say whether that subject is plural."""
    j, candidate = i - 1, None
    while j >= 0:
        tok = toks[j]
        w = _word(tok)
        if tok in (",", ";", ":"):
            break
        if w in PREPS:
            candidate = None
        elif w in QUANTIFIERS or _is_number(tok):
            if candidate is None:
                candidate = tok
        elif w in DETS or w in CONJ:
            pass
        elif candidate is None:
            if w.endswith("s") and (_stems_s(w) & VERBS) \
                    and not _noun_position(toks, j):
                return False              # `calls` in *which instrument calls*
            candidate = tok
        j -= 1
    return candidate is not None and _plural(candidate)


def _verb_at(toks: list[str], i: int, in_with: bool) -> bool:
    w = _word(toks[i])
    if not w:
        return False
    if w in AUX or w in PAST:
        return True
    if w.endswith("s") and (_stems_s(w) & VERBS) and not _noun_position(toks, i):
        return True
    if w.endswith("ed") and (_stems_ed(w) & VERBS) and i > 0 and not in_with:
        prev = _word(toks[i - 1]) if toks[i - 1] not in (",", ";", ":") else ","
        nxt = _word(toks[i + 1]) if i + 1 < len(toks) else ""
        if prev not in PARTICIPLE_CUES and nxt not in PARTICIPLE_NEXT:
            return True
    return w in VERBS and _plural_subject(toks, i)


PARTICIPLES = set("drawn shown given taken written known seen been done gone "
                  "broken chosen built made held kept read set split spread "
                  "put cut hit let be being".split())


def _participle_like(tok: str) -> bool:
    w = _word(tok)
    return (w.endswith("ed") or w in PARTICIPLES
            or (w in VERBS and w not in AUX))


def _strip_clauses(toks: list[str], verbs: set[int]) -> list[str]:
    """Drop every clause opened by a subordinator, relative word or contact
    determiner, up to the next comma. A sentence-initial free relative is
    kept, and `has_finite_verb` asks it for two verbs instead."""
    out, i = [], 0
    while i < len(toks):
        tok, w = toks[i], _word(toks[i])
        if w in SUBORDINATORS and i > 0:
            while i < len(toks) and toks[i] not in (",", ";", ":"):
                i += 1
            continue
        prev = toks[i - 1] if i else ","
        contact = (w in CONTACT_DETS and tok == tok.lower() and i > 0
                   and prev not in (",", ";", ":")
                   and _word(prev) not in NOUN_POSITION | ADVERBS
                   and not _word(prev).endswith("ing")
                   and (i - 1) not in verbs and not _is_number(prev))
        if contact or (w in CLAUSE_OPENERS and i > 0):
            j = i
            while j < len(toks) and toks[j] not in (",", ";", ":") \
                    and j not in verbs:
                j += 1
            if j < len(toks) and j in verbs:
                j += 1                     # the clause's own verb ...
                while j < len(toks) and _participle_like(toks[j]):
                    j += 1                 # ... and its participle (is read)
            i = j
            continue
        out.append(tok)
        i += 1
    return out


def finite_verbs(sentence: str) -> list[str]:
    """Every token of `sentence` the rules accept as a main-clause finite
    verb, or an empty list."""
    toks = _tokens(clean(sentence))
    first = set(i for i, t in enumerate(toks) if _verb_at(toks, i, False))
    toks = _strip_clauses(toks, first)
    found, in_with = [], False
    for i, tok in enumerate(toks):
        w = _word(tok)
        if tok in (",", ";", ":"):
            in_with = False
        elif w == "with":
            in_with = True
        if _verb_at(toks, i, in_with):
            found.append(tok)
    if toks and _word(toks[0]) in FREE_RELATIVE and len(found) < 2:
        return []
    return found


def has_finite_verb(sentence: str) -> bool:
    return bool(finite_verbs(sentence))


def is_quantity(label: str) -> bool:
    """A label that is a number, a unit, a symbol or a category, not prose."""
    words = _tokens(clean(label))
    return len(words) < 4
