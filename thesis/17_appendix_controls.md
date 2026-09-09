# Appendix B. The constructed negative controls

This appendix lists the constructed negative controls of this project. It is
derived from the test modules rather than counted by hand and is committed as
`thesis/control_inventory.tsv`.

There are **437 named constructed checks across 19 suites.**

## B.1 What is counted, and what a unit is

Two conventions are in use across the project and both are counted, because it
grew them at different times and neither is wrong. Some suites group their
checks into named test functions, while others call a helper once per
assertion, naming each check as it runs. **A named check is the finer unit, so
it is preferred where a module has any, and the table records which unit each
row counts**, because the two are not the same size and adding them together
silently would be worse than either.

The count is derived by parsing each module. Any call to one of the project's
check helpers whose first argument is a string literal is one named control,
and any top-level function whose name matches the project's test-function
convention is one grouped control. A control deleted from a module disappears
from this appendix, and one added appears in it.

## B.2 Where the controls are

| suite | unit | controls |
|---|---|---|
| `s6_test_selection.py` | named check | 24 |
| `s7_test_tree.py` | named check | 18 |
| `s8_test_flanks.py` | named check | 49 |
| `s9_test_codon.py` | named check | 14 |
| `s10_test_evidence.py` | named check | 58 |
| `s11_test_structures.py` | named check | 24 |
| `s12_test_expression.py` | named check | 11 |
| `s13_test_recon.py` | test function | 14 |
| `s15_test_loss.py` | test function | 15 |
| `s15b_test_counts.py` | named check | 21 |
| `s16_test_dup.py` | named check | 21 |
| `s17_test_constraint.py` | test function | 14 |
| `s18_test_audit.py` | named check | 42 |
| `s19_test_methods.py` | named check | 33 |
| `s21_test_arch.py` | named check | 9 |
| `s21_test_claims.py` | named check | 12 |
| `s22_test_ligand.py` | named check | 23 |
| `s22_test_lineage.py` | named check | 18 |
| `s24_test_supp.py` | named check | 17 |

Further controls run inside the analysis modules themselves rather than in a
test module, comprising the bait screen's three synthetic failures, the
chunked-genome equivalence test, the resume test and the
one-search-two-sensitivities equivalence test. They are not in the table
because they are not separately addressable checks.

## B.3 What the controls are checks on

Almost none of them checks that a routine returns the right answer. They check
that it **refuses**, and that it **can act**. Chapter 14 gives the argument.
The short form is that in this project almost every rule returns a plausible
number when it is wrong, and several of the thesis's results are zeros that
would be indistinguishable from a rule that cannot fire.

Every suite runs before its task writes anything, and a failure refuses the
build.

## B.4 Six of the controls found something

Six are described in Chapter 14: the chimera screen's control that failed
because the test was wrong rather than the rule; a reassignment rule whose
positive and negative halves both failed and whose rule was wrong; an
order-invariance check that caught real hash-seeded nondeterminism in a
committed table; a within-protein control selected by a name prefix that swept
the element under test into the control set; a self-test that overwrote the
committed table it was testing; and two mutation tests that passed on broken
code until the case was rebuilt to exercise the guard it was aimed at.
