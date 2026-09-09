"""S25 — the negative-control inventory, derived rather than counted by hand.

The manuscript reports the project's constructed negative controls as a
Methods sentence and a number. Appendix B is that number turned into a table,
and it is **derived from the test modules themselves** rather than typed: a
control that is deleted disappears from the appendix, and one that is added
appears in it, without anybody remembering to update a count.

What is counted is a *constructed control* — a top-level function in a task's
test module whose name begins `t<digit>`, which is this project's convention
for one named check on one rule. Helpers, fixtures and the `self_test` /
`check` entry points are excluded by the same rule, so the number is the
number of checks and not the number of functions.

The docstring's first line is carried through as what the check refuses, so
the appendix says what each suite is for rather than only how large it is.

  python scripts/s25_controls.py         # write thesis/control_inventory.tsv
"""

from __future__ import annotations

import ast
import re
import sys

import s25_lib as L

SCRIPTS = L.ROOT / "scripts"

#: `t1_...`, `t14b_...` — one named constructed control.
TEST_FN = re.compile(r"^t\d+[a-z]?_")

FIELDS = ["module", "task", "unit", "n_controls", "controls"]


#: The call that names one constructed control. Every suite in this project
#: reports each check by name as it runs, through one of these helpers.
CHECK_FN = {"check", "_check", "require", "ok", "_t", "_case"}


def module_controls(path) -> tuple[list[str], list[str]]:
    """(named checks, test functions) for one test module.

    Two conventions are in use and both are counted, because the project grew
    them at different times and neither is wrong. Some suites group their
    checks into `t1_...`, `t2_...` functions; others call a `check(name, ok)`
    helper once per assertion. A named check is the finer unit, so it is
    preferred where a module has any, and the `unit` column records which was
    counted so the two are never silently added together.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    fns = [n.name for n in ast.walk(tree)
           if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
           and TEST_FN.match(n.name)]
    checks = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        name = (node.func.id if isinstance(node.func, ast.Name)
                else node.func.attr if isinstance(node.func, ast.Attribute)
                else None)
        if name in CHECK_FN and isinstance(node.args[0], ast.Constant) \
                and isinstance(node.args[0].value, str):
            checks.append(node.args[0].value)
    return checks, fns


def build() -> list[dict]:
    rows = []
    for path in sorted(SCRIPTS.glob("s*test*.py")):
        checks, fns = module_controls(path)
        if checks:
            unit, names = "named check", checks
        elif fns:
            unit, names = "test function", fns
        else:
            continue
        task = path.name.split("_")[0].upper()
        rows.append({
            "module": f"scripts/{path.name}",
            "task": task,
            "unit": unit,
            "n_controls": str(len(names)),
            "controls": "; ".join(names),
        })
    return rows


def run() -> int:
    rows = build()
    L.write_tsv(L.TH / "control_inventory.tsv", rows, FIELDS)
    total = sum(int(r["n_controls"]) for r in rows)
    print(f"[s25 controls] {total} constructed negative controls across "
          f"{len(rows)} modules -> thesis/control_inventory.tsv")
    if not rows:
        print("  [FAIL] no test modules found", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
