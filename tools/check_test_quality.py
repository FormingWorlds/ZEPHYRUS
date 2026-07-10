#!/usr/bin/env python3
"""AST-based test-quality linter for the ZEPHYRUS test suite.

Enforces the ZEPHYRUS test-quality contract:

- Every test file must declare a module-level ``pytestmark`` containing a tier
  marker (``unit`` / ``smoke`` / ``integration`` / ``slow``).
- Test functions must contain at least 2 assertion statements OR a discriminating
  property-based assertion. Single-assert tests are a known weak pattern.
- Forbidden weak assertions when they stand alone as the sole meaningful
  check in the test: ``result is not None``, ``result > 0``,
  ``len(result) > 0``, ``isinstance(result, dict)``, ``result is None``.
  A weak assertion that accompanies a stronger primary assertion (a
  discrimination guard that uses ``val > 0`` as a sign guard alongside a
  ``pytest.approx(...)`` pin, for example) is NOT flagged.
- Every test function must have a docstring.
- ``==`` adjacent to a numeric literal in a test body is a likely float-comparison
  bug (use ``pytest.approx`` instead).
- Optional dependencies (``hypothesis``) imported at module top
  without a preceding ``pytest.importorskip('<name>')`` (the ``pip install
  --no-deps`` CI image otherwise fails collection). ``mors`` is a hard
  dependency of ZEPHYRUS and is intentionally NOT in this set.

Two modes:

* ``--baseline``  Walk the test suite, write per-rule violation counts to
  ``tools/test_quality_baseline.json``. Run this only after a deliberate sweep
  that has reduced violations; commits should not raise the baseline.
* ``--check``     CI mode. Walk the suite, compare current violation counts to
  the baseline. Exit non-zero if any rule's violation count exceeds the
  baseline. Print the offending files + functions.

Optionally:

* ``--reference-pinned-status``  Print the physics source files that lack at
  least one ``@pytest.mark.reference_pinned`` test in the matching
  ``tests/test_<source>.py``. Does not exit non-zero on its own (advisory).
* ``--physics-invariant-status``  Print physics-source tests that assert no
  invariant and are not tagged ``@pytest.mark.physics_invariant``. Advisory.

All exits in ``--check`` mode use exit code 1 on regression; 0 otherwise.

The script reads no configuration outside its own constants and the baseline
file. It is intentionally dependency-free (pure stdlib) so it can run in any
CI environment.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = REPO_ROOT / 'tests'
SRC_DIR = REPO_ROOT / 'src' / 'zephyrus'
BASELINE_PATH = REPO_ROOT / 'tools' / 'test_quality_baseline.json'

TIER_MARKERS = {'unit', 'smoke', 'integration', 'slow'}

# Optional dependencies. Any test module that imports one of these MUST
# precede the import with ``pytest.importorskip('<name>')`` at module
# scope, otherwise CI's ``pip install --no-deps`` build fails collection.
# ``mors`` is a hard runtime dependency (``fwl-mors`` in
# [project.dependencies]) and so is not listed here: importing it
# unconditionally is correct.
OPTIONAL_DEPS = {
    'hypothesis',
}

# Physics source files. Each one must have a 1:1 test file at
# ``tests/test_<source>.py``.
# Each source must have at least one @pytest.mark.physics_invariant test
# and at least one @pytest.mark.reference_pinned test in its companion
# test file.
PHYSICS_SOURCES = {
    'escape.py',
}

# Utility sources are exempt from the physics-invariant / reference-pinned
# requirement but still subject to the anti-happy-path rules.
UTILITY_SOURCES = {
    '__init__.py',
    'constants.py',
    'planets_parameters.py',
}


def _is_weak_assert(node: ast.Assert) -> str | None:
    """Return a label if ``node`` is a forbidden weak standalone assertion."""
    test = node.test
    # `assert x is None` / `assert x is not None`
    if isinstance(test, ast.Compare) and len(test.ops) == 1:
        op = test.ops[0]
        right = test.comparators[0]
        if (
            isinstance(op, (ast.Is, ast.IsNot))
            and isinstance(right, ast.Constant)
            and right.value is None
        ):
            return 'is_none_or_not_none'
        # `assert len(x) > 0` (must come BEFORE the bare `> 0` check,
        # since it is the more specific shape).
        if (
            isinstance(op, ast.Gt)
            and isinstance(test.left, ast.Call)
            and isinstance(test.left.func, ast.Name)
            and test.left.func.id == 'len'
            and isinstance(right, ast.Constant)
            and right.value == 0
        ):
            return 'len_gt_zero'
        # `assert x > 0`
        if isinstance(op, ast.Gt) and isinstance(right, ast.Constant) and right.value == 0:
            return 'gt_zero'
    return None


def _is_isinstance_assert(node: ast.Assert) -> bool:
    test = node.test
    return (
        isinstance(test, ast.Call)
        and isinstance(test.func, ast.Name)
        and test.func.id == 'isinstance'
    )


def _is_numpy_testing(node: ast.AST) -> bool:
    """True if ``node`` represents the ``numpy.testing`` or ``np.testing`` module."""
    if not isinstance(node, ast.Attribute):
        return False
    if node.attr != 'testing':
        return False
    if isinstance(node.value, ast.Name) and node.value.id in ('np', 'numpy'):
        return True
    return False


def _is_exact_zero(value) -> bool:
    """True for the sentinel ``0.0`` / ``-0.0`` float comparand.

    Asserting an exact-zero result is a legitimate physics check and does
    not need ``pytest.approx``: there is no rounding error to absorb.
    Comparing against any other float literal is flagged.
    """
    return isinstance(value, float) and value == 0.0


def _unwrap_unary_minus_float(operand: ast.AST) -> float | None:
    """Return the float value of ``UnaryOp(USub, Constant(float))``, or None.

    The AST parses ``-1.5`` as ``UnaryOp(USub, Constant(1.5))``, NOT as
    ``Constant(-1.5)``. A naive ``isinstance(node, ast.Constant)`` check
    would miss negative float literals; this helper accepts both forms.
    """
    if (
        isinstance(operand, ast.UnaryOp)
        and isinstance(operand.op, ast.USub)
        and isinstance(operand.operand, ast.Constant)
        and isinstance(operand.operand.value, float)
    ):
        return -operand.operand.value
    return None


def _float_literal_value(operand: ast.AST) -> float | None:
    """Return the float value of a positive or negative float-literal node."""
    if isinstance(operand, ast.Constant) and isinstance(operand.value, float):
        return operand.value
    return _unwrap_unary_minus_float(operand)


def _has_float_eq(node: ast.AST) -> bool:
    """Return True if any descendant uses ``==`` against a non-zero float literal.

    Accepts both ``Constant(1.5)`` (positive literal) and
    ``UnaryOp(USub, Constant(1.5))`` (negative literal). The exact-zero
    carve-out applies to both signs.
    """
    for child in ast.walk(node):
        if isinstance(child, ast.Compare):
            for op, right in zip(child.ops, child.comparators):
                if not isinstance(op, ast.Eq):
                    continue
                right_val = _float_literal_value(right)
                if right_val is not None and not _is_exact_zero(right_val):
                    return True
                left_val = _float_literal_value(child.left)
                if left_val is not None and not _is_exact_zero(left_val):
                    return True
    return False


def _module_pytestmark_tier(tree: ast.Module) -> str | None:
    """Return the tier marker declared in a module-level ``pytestmark``, or None."""
    for stmt in tree.body:
        if not isinstance(stmt, ast.Assign):
            continue
        if not (len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Name)):
            continue
        if stmt.targets[0].id != 'pytestmark':
            continue
        marks = stmt.value
        nodes = marks.elts if isinstance(marks, (ast.List, ast.Tuple)) else [marks]
        for n in nodes:
            tier = _tier_of_mark_node(n)
            if tier is not None:
                return tier
    return None


def _tier_of_mark_node(n: ast.AST) -> str | None:
    """Given a ``pytest.mark.<x>`` or ``pytest.mark.<x>(...)`` node, return the tier name."""
    if isinstance(n, ast.Call):
        n = n.func
    if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Attribute):
        if (
            isinstance(n.value.value, ast.Name)
            and n.value.value.id == 'pytest'
            and n.value.attr == 'mark'
            and n.attr in TIER_MARKERS
        ):
            return n.attr
    return None


FuncDef = (ast.FunctionDef, ast.AsyncFunctionDef)


def _iter_test_functions(tree: ast.Module):
    """Yield every ``test_*`` function defined at module scope or as a
    method of a class at module scope.

    Does NOT recurse into function bodies: a `def test_x()` defined
    inside another function body is a local helper, not an independent
    pytest test, and must not be inspected as one. Recursing via
    ``ast.walk`` would treat such helpers as phantom tests.

    Async functions (`async def test_x`) are yielded alongside
    synchronous ones; pytest-asyncio and other plugins run them.
    """
    for stmt in tree.body:
        if isinstance(stmt, FuncDef) and stmt.name.startswith('test_'):
            yield stmt
        elif isinstance(stmt, ast.ClassDef):
            for sub in stmt.body:
                if isinstance(sub, FuncDef) and sub.name.startswith('test_'):
                    yield sub


def _func_markers(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> set[str]:
    """All ``pytest.mark.<x>`` markers on a function definition."""
    out: set[str] = set()
    for dec in fn.decorator_list:
        n = dec
        if isinstance(n, ast.Call):
            n = n.func
        if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Attribute):
            if (
                isinstance(n.value.value, ast.Name)
                and n.value.value.id == 'pytest'
                and n.value.attr == 'mark'
            ):
                out.add(n.attr)
    return out


def _docstring_of(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> str | None:
    if (
        fn.body
        and isinstance(fn.body[0], ast.Expr)
        and isinstance(fn.body[0].value, ast.Constant)
    ):
        v = fn.body[0].value.value
        if isinstance(v, str):
            return v
    return None


class Violations:
    """Aggregate counts and per-violation details for one rule scan."""

    def __init__(self):
        self.counts: dict[str, int] = defaultdict(int)
        self.details: dict[str, list[str]] = defaultdict(list)

    def add(self, rule: str, where: str) -> None:
        self.counts[rule] += 1
        self.details[rule].append(where)

    def to_baseline(self) -> dict[str, int]:
        return dict(self.counts)


def _count_implicit_assertions(node: ast.AST) -> int:
    """Count assertion-equivalents that are not bare ``assert`` statements.

    Recognized patterns:

    - ``with pytest.raises(...)`` blocks. A ``match=`` keyword counts as a
      second implicit assertion: it imposes a separate falsifiable
      constraint on the exception message, distinct from the type check.
    - ``with pytest.warns(...)`` blocks. Same ``match=`` rule.
    - ``with pytest.deprecated_call()`` blocks.
    - ``mock.assert_called_with(...)`` / ``assert_called_once_with(...)`` /
      ``assert_not_called(...)`` etc. method calls on a Mock object.
    - ``pytest.fail(...)`` calls.
    - ``np.testing.assert_*`` family.
    """
    count = 0
    for child in ast.walk(node):
        if isinstance(child, (ast.With, ast.AsyncWith)):
            for item in child.items:
                ctx = item.context_expr
                if isinstance(ctx, ast.Call) and isinstance(ctx.func, ast.Attribute):
                    if (
                        isinstance(ctx.func.value, ast.Name)
                        and ctx.func.value.id == 'pytest'
                        and ctx.func.attr in ('raises', 'warns', 'deprecated_call')
                    ):
                        count += 1
                        # A ``match=`` argument is a separate falsifiable
                        # constraint on the exception or warning message:
                        # ``pytest.raises(ValueError, match='non-negative')``
                        # is strictly stronger than ``pytest.raises(ValueError)``.
                        if any(kw.arg == 'match' for kw in ctx.keywords):
                            count += 1
        if isinstance(child, ast.Call) and isinstance(child.func, ast.Attribute):
            attr = child.func.attr
            if attr.startswith('assert_called') or attr == 'assert_not_called':
                count += 1
            elif attr.startswith('assert_') and _is_numpy_testing(child.func.value):
                count += 1
            elif isinstance(child.func.value, ast.Name) and child.func.value.id == 'pytest':
                if attr == 'fail':
                    count += 1
    return count


def _importorskip_targets(tree: ast.Module) -> set[str]:
    """Set of dependency names protected by a module-scope ``pytest.importorskip``."""
    out: set[str] = set()
    for node in tree.body:
        call = None
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            call = node.value
        elif isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            call = node.value
        if call is None or not isinstance(call.func, ast.Attribute):
            continue
        if (
            isinstance(call.func.value, ast.Name)
            and call.func.value.id == 'pytest'
            and call.func.attr == 'importorskip'
            and call.args
            and isinstance(call.args[0], ast.Constant)
            and isinstance(call.args[0].value, str)
        ):
            out.add(call.args[0].value)
    return out


def _missing_importorskip(tree: ast.Module) -> list[str]:
    """Return optional-dep names imported at module top without a matching importorskip."""
    protected = _importorskip_targets(tree)
    missing: list[str] = []
    seen: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split('.', 1)[0]
                if root in OPTIONAL_DEPS and root not in protected and root not in seen:
                    missing.append(root)
                    seen.add(root)
        elif isinstance(node, ast.ImportFrom):
            if node.module is None:
                continue
            root = node.module.split('.', 1)[0]
            if root in OPTIONAL_DEPS and root not in protected and root not in seen:
                missing.append(root)
                seen.add(root)
    return missing


def check_file(path: Path) -> Violations:
    v = Violations()
    try:
        rel = str(path.relative_to(REPO_ROOT))
    except ValueError:
        rel = str(path)
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError as e:
        v.add('parse_error', f'{rel}: {e}')
        return v

    module_tier = _module_pytestmark_tier(tree)
    if module_tier is None:
        v.add('missing_module_pytestmark', rel)

    for dep in _missing_importorskip(tree):
        v.add('missing_importorskip', f'{rel}: {dep}')

    for node in _iter_test_functions(tree):
        where = f'{rel}::{node.name}'

        if _docstring_of(node) is None:
            v.add('missing_docstring', where)

        if _has_float_eq(node):
            v.add('float_eq_literal', where)

        asserts = [n for n in ast.walk(node) if isinstance(n, ast.Assert)]
        n_assert = len(asserts) + _count_implicit_assertions(node)
        if n_assert == 0:
            v.add('no_assertions', where)
        elif n_assert == 1:
            v.add('single_assert', where)

        # Weak-assertion shapes. Flag only when the weak assertion stands
        # alone as the sole meaningful check in the test.
        if n_assert == 1 and len(asserts) == 1:
            label = _is_weak_assert(asserts[0])
            if label is not None:
                v.add(f'weak_assert_{label}', where)
            if _is_isinstance_assert(asserts[0]):
                v.add('weak_assert_only_isinstance', where)

    return v


def walk_tests() -> Violations:
    total = Violations()
    for p in sorted(TESTS_DIR.rglob('test_*.py')):
        if '__pycache__' in p.parts:
            continue
        v = check_file(p)
        for rule, n in v.counts.items():
            total.counts[rule] += n
        for rule, details in v.details.items():
            total.details[rule].extend(details)
    return total


def _file_has_decorator(path: Path, decorator_name: str) -> bool:
    """True if any function or class in ``path`` carries ``@pytest.mark.<decorator_name>``."""
    try:
        tree = ast.parse(path.read_text())
    except SyntaxError:
        return False
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for dec in node.decorator_list:
            d = dec
            if isinstance(d, ast.Call):
                d = d.func
            if isinstance(d, ast.Attribute) and isinstance(d.value, ast.Attribute):
                if (
                    isinstance(d.value.value, ast.Name)
                    and d.value.value.id == 'pytest'
                    and d.value.attr == 'mark'
                    and d.attr == decorator_name
                ):
                    return True
    return False


def reference_pinned_status() -> list[str]:
    """Return physics source files lacking a reference_pinned test.

    For each source in PHYSICS_SOURCES, check whether the companion test file
    ``tests/test_<source>`` (with .py preserved) exists and contains at least
    one ``@pytest.mark.reference_pinned`` decoration. Sources without such a
    test are reported.
    """
    missing = []
    for source in sorted(PHYSICS_SOURCES):
        # source is e.g. 'escape.py' -> companion is tests/test_escape.py
        stem = source[:-3] if source.endswith('.py') else source
        test_path = TESTS_DIR / f'test_{stem}.py'
        if not test_path.exists():
            missing.append(f'{source} (no test_{stem}.py)')
            continue
        if not _file_has_decorator(test_path, 'reference_pinned'):
            missing.append(f'{source} (test_{stem}.py has no @pytest.mark.reference_pinned)')
    return missing


def physics_invariant_status() -> list[str]:
    """Return physics-source tests without an explicit invariant marker or
    property-based assertion language.

    The heuristic is intentionally simple: a physics-source test must either
    carry the marker OR contain one of the keywords in its body. This is
    advisory, not blocking.
    """
    flagged = []
    keywords = {'approx', 'assert_allclose', 'monoton', 'conserve', 'symmetric', 'positive'}
    for source in sorted(PHYSICS_SOURCES):
        stem = source[:-3] if source.endswith('.py') else source
        test_path = TESTS_DIR / f'test_{stem}.py'
        if not test_path.exists():
            continue
        try:
            tree = ast.parse(test_path.read_text())
        except SyntaxError:
            continue
        rel = str(test_path.relative_to(REPO_ROOT))
        for node in _iter_test_functions(tree):
            markers = _func_markers(node)
            if 'physics_invariant' in markers:
                continue
            body_src = ast.unparse(node)
            if any(kw in body_src for kw in keywords):
                continue
            flagged.append(f'{rel}::{node.name}')
    return flagged


def load_baseline() -> dict[str, int]:
    if not BASELINE_PATH.exists():
        return {}
    return json.loads(BASELINE_PATH.read_text())


def cmd_baseline() -> int:
    v = walk_tests()
    # Guard against accidental regeneration that would raise the baseline.
    old = load_baseline()
    old_total = sum(old.values())
    new_total = sum(v.counts.values())

    allow = os.environ.get('ZEPHYRUS_TEST_QUALITY_ALLOW_REGRESS') == '1'
    if old and new_total > old_total and not allow:
        print(
            f'Refusing to regenerate baseline: new total ({new_total}) exceeds '
            f'old total ({old_total}). The baseline should only ratchet downward.\n'
            f'If this is intentional (e.g. a new rule was added that surfaces '
            f'pre-existing violations), set ZEPHYRUS_TEST_QUALITY_ALLOW_REGRESS=1.',
            file=sys.stderr,
        )
        return 2
    BASELINE_PATH.write_text(json.dumps(v.to_baseline(), indent=2, sort_keys=True) + '\n')
    print(f'Wrote baseline: {BASELINE_PATH.relative_to(REPO_ROOT)}')
    for rule in sorted(v.counts):
        print(f'  {rule}: {v.counts[rule]}')
    if old:
        delta = new_total - old_total
        print(f'  total: {new_total} ({delta:+d} vs previous baseline {old_total})')
    return 0


def cmd_check() -> int:
    baseline = load_baseline()
    v = walk_tests()
    failed = False
    print('Rule                                       Baseline   Current   Status')
    print('-' * 76)
    all_rules = sorted(set(baseline) | set(v.counts))
    for rule in all_rules:
        b = baseline.get(rule, 0)
        c = v.counts.get(rule, 0)
        status = 'OK' if c <= b else 'REGRESSION'
        if c > b:
            failed = True
        print(f'{rule:42} {b:>8} {c:>9}   {status}')
    # Cross-rule total-violation guard.
    total_baseline = sum(baseline.values())
    total_current = sum(v.counts.values())
    print('-' * 76)
    total_status = 'OK' if total_current <= total_baseline else 'REGRESSION'
    print(f'{"TOTAL":42} {total_baseline:>8} {total_current:>9}   {total_status}')
    if total_current > total_baseline:
        failed = True
    if failed:
        print()
        print('New violations vs baseline:')
        for rule in all_rules:
            b = baseline.get(rule, 0)
            c = v.counts.get(rule, 0)
            if c > b:
                offenders = v.details.get(rule, [])
                print(f'\n  {rule} (+{c - b}):')
                for offender in offenders[:5]:
                    print(f'    {offender}')
                if len(offenders) > 5:
                    print(f'    ... and {len(offenders) - 5} more')
        print()
        print('Reduce violations or, after a deliberate sweep, regenerate the baseline:')
        print('  python tools/check_test_quality.py --baseline')
        return 1
    return 0


def cmd_reference_pinned_status() -> int:
    missing = reference_pinned_status()
    if not missing:
        print('All physics sources have a @pytest.mark.reference_pinned test.')
        return 0
    print('Physics sources missing a @pytest.mark.reference_pinned test:')
    for m in missing:
        print(f'  {m}')
    return 0


def cmd_physics_invariant_status() -> int:
    flagged = physics_invariant_status()
    if not flagged:
        print(
            'All physics-source tests either carry @pytest.mark.physics_invariant '
            'or use property-based language.'
        )
        return 0
    print(
        'Physics-source tests without @pytest.mark.physics_invariant and no '
        'property-based assertion language:'
    )
    for f in flagged:
        print(f'  {f}')
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--baseline',
        action='store_true',
        help='Regenerate tools/test_quality_baseline.json from current state.',
    )
    group.add_argument(
        '--check',
        action='store_true',
        help='CI mode: fail if violations exceed baseline.',
    )
    group.add_argument(
        '--reference-pinned-status',
        action='store_true',
        help='Advisory: list physics sources missing a @pytest.mark.reference_pinned test.',
    )
    group.add_argument(
        '--physics-invariant-status',
        action='store_true',
        help='Advisory: list physics-source tests without an invariant marker or '
        'property-based language.',
    )
    args = parser.parse_args()
    if args.baseline:
        return cmd_baseline()
    if args.check:
        return cmd_check()
    if args.reference_pinned_status:
        return cmd_reference_pinned_status()
    if args.physics_invariant_status:
        return cmd_physics_invariant_status()
    return 0


if __name__ == '__main__':
    sys.exit(main())
