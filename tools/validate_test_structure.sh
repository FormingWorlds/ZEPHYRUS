#!/usr/bin/env bash
# Marker-validation gate for the ZEPHYRUS test suite.
#
# Walks every *.py file under tests/, finds every `def test_*` definition,
# and verifies that the function (or its enclosing class, or the
# module-level pytestmark) carries exactly one of:
#
#     @pytest.mark.unit
#     @pytest.mark.smoke
#     @pytest.mark.integration
#     @pytest.mark.slow
#     @pytest.mark.skip
#
# Any unmarked test is printed as <file>:<line> and the script exits 1.
# Run from repository root:
#
#     bash tools/validate_test_structure.sh
#
# CI invokes this before pytest in .github/workflows/tests.yaml.

set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

if [ ! -d tests ]; then
    echo "ERROR: tests/ not found in $REPO_ROOT" >&2
    exit 2
fi

python3 - <<'PY'
import ast
import pathlib
import sys

ROOT = pathlib.Path('tests').resolve()
TIER = {'unit', 'smoke', 'integration', 'slow'}
SKIP = {'skip'}
ALLOWED = TIER | SKIP

failures: list[str] = []
total = 0


def mark_names_from_decorators(decorators):
    out = []
    for d in decorators:
        # @pytest.mark.NAME
        if isinstance(d, ast.Attribute) and isinstance(d.value, ast.Attribute):
            if (
                isinstance(d.value.value, ast.Name)
                and d.value.value.id == 'pytest'
                and d.value.attr == 'mark'
            ):
                out.append(d.attr)
        # @pytest.mark.NAME(...)
        elif isinstance(d, ast.Call) and isinstance(d.func, ast.Attribute):
            f = d.func
            if (
                isinstance(f.value, ast.Attribute)
                and isinstance(f.value.value, ast.Name)
                and f.value.value.id == 'pytest'
                and f.value.attr == 'mark'
            ):
                out.append(f.attr)
    return out


def module_marks(tree):
    marks = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == 'pytestmark':
                    val = node.value
                    items = val.elts if isinstance(val, (ast.List, ast.Tuple)) else [val]
                    for item in items:
                        marks.extend(mark_names_from_decorators([item]))
    return marks


def walk(tree, parent_marks):
    out = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            cls_marks = parent_marks + mark_names_from_decorators(node.decorator_list)
            out.extend(walk(node, cls_marks))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith('test_'):
                fn_marks = mark_names_from_decorators(node.decorator_list)
                out.append((node, parent_marks + fn_marks))
    return out


for path in sorted(ROOT.rglob('*.py')):
    if any(part.startswith('.') for part in path.parts):
        continue
    if path.name == 'conftest.py':
        continue
    try:
        tree = ast.parse(path.read_text(encoding='utf-8'))
    except SyntaxError as exc:
        failures.append(f'{path}:{exc.lineno}: SyntaxError {exc.msg}')
        continue
    mod_marks = module_marks(tree)
    for fn, marks in walk(tree, mod_marks):
        total += 1
        tier_marks = [m for m in marks if m in TIER]
        has_skip = any(m in SKIP for m in marks)
        rel = path.relative_to(ROOT.parent)
        if len(tier_marks) > 1:
            failures.append(
                f'{rel}:{fn.lineno}: {fn.name} has multiple tier markers '
                f'{sorted(set(tier_marks))}; exactly one of '
                'unit / smoke / integration / slow is required.'
            )
        elif has_skip and len(tier_marks) >= 1:
            failures.append(
                f'{rel}:{fn.lineno}: {fn.name} carries both a tier marker '
                f'{sorted(set(tier_marks))} and skip; use exactly one of '
                'unit / smoke / integration / slow / skip.'
            )
        elif len(tier_marks) == 0 and not has_skip:
            failures.append(
                f'{rel}:{fn.lineno}: {fn.name} carries no marker '
                '(need one of unit / smoke / integration / slow / skip)'
            )

if failures:
    print(f'Marker validation FAILED on {len(failures)} of {total} tests.')
    for f in failures:
        print(f'  {f}')
    sys.exit(1)

print(f'Marker validation OK: {total} tests, all carry a marker.')
PY
