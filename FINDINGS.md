# FINDINGS.md — what this project has discovered

The running biological story, one dated entry per completed roadmap task.
Findings before methods, plain language, no file paths. Claims that are not
yet confirmed are marked *(pending: which task confirms it)*.

`SESSION_LOG.md` is the chronological record of what was *run*; this file is
the record of what is *known*.

---

## 2026-08-18 — Project set up; nothing measured yet

No analysis has run. What exists is the question and the instrument.

**The subject.** The IP3 receptor is the endoplasmic reticulum's
ligand-gated calcium-release channel — the source of most agonist-evoked
calcium signals in non-muscle cells. Vertebrates have three of them
(ITPR1, ITPR2, ITPR3), each about 2,700 residues, assembled as tetramers,
and each associated with a different human disease: cerebellar ataxia and
Gillespie syndrome for ITPR1, an inability to sweat for ITPR2, and a
peripheral neuropathy for ITPR3.

**What is already visible in the databases** (checked 2026-08-18, and to be
re-derived properly in S2 before any of it is quoted):

- The family looks overwhelmingly like an animal family: 12,149 of the
  12,338 proteins carrying its core domain are metazoan.
- Yet there are 40 plant and 41 fungal records in a family that textbooks
  say plants and fungi lack — while *Arabidopsis* and baker's yeast have
  none at all. Something there is either a real and under-appreciated
  branch of the family, or a set of database errors. Both would be worth
  knowing. *(pending: S2, S20, S23)*
- Zebrafish carries four: itpr1a, itpr1b, itpr2 and itpr3. The first two
  look like a teleost duplicate pair. *(pending: S16)*

**The complication that shapes the whole project.** Ryanodine receptors —
the other big calcium-release channel, twice the size — carry every domain
that marks an IP3 receptor. They will appear in every search this project
runs. That is a nuisance and an opportunity: it means the two families can
be counted by the same instrument and compared directly, and RyR gives the
IP3R tree a proper outgroup instead of a guess.

**What the project intends to find out.** Whether the family really is
absent from plants and fungi; where the three vertebrate paralogs came from
and whether they are the same age as the three ryanodine receptors; whether
any of them has ever been lost; whether the conserved core explains the
three diseases; and how much of what the databases say about this family is
wrong.
