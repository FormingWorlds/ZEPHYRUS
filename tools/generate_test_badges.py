"""Generate a shields.io endpoint-badge JSON file for the test count.

The script invokes ``pytest --collect-only -q`` to count tests without
executing them, then writes the count to a JSON file under the ``--out``
directory in the shields.io endpoint-badge schema:

    {"schemaVersion": 1, "label": "<text>", "message": "<count>", "color": "blue"}

Output is a single file, ``tests-total.json`` (label "tests"), holding the
count of ``not skip`` tests. The suite carries no unit/integration marker
split, so the badge reports one total. A marker split can be added later by
extending ``_BADGES`` with the extra expressions.

Usage
-----
    python tools/generate_test_badges.py --out <dir>

The publish workflow points ``--out`` at a scratch directory and copies the
resulting JSON onto the ``badges`` branch, which shields.io reads by raw URL.

Notes
-----
Running the script does not execute the test suite; only collection is
triggered. Pytest exit code 5 ("no tests collected") is treated as a
successful zero count and the badge writes ``"message": "0"``. Any other
non-zero exit is a hard failure, so an incomplete install that breaks
collection fails the run loudly instead of publishing a wrong count.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

_COLLECT_RE = re.compile(r'^(\d+)(?:/\d+)?\s+tests?\s+collected\b', re.MULTILINE)

_BADGES: tuple[tuple[str, str, str], ...] = (('total', 'tests', 'not skip'),)


def count_tests(marker_expr: str) -> int:
    """Run pytest collection and return the number of selected tests.

    Parameters
    ----------
    marker_expr : str
        Pytest marker expression passed via ``-m``.

    Returns
    -------
    int
        Number of tests pytest collected for the given marker. Exit
        code 5 ("no tests collected") is mapped to 0.

    Raises
    ------
    RuntimeError
        If pytest exits with a non-zero code other than 5, or if the
        trailing summary line cannot be parsed from stdout.
    """
    proc = subprocess.run(
        [sys.executable, '-m', 'pytest', '--collect-only', '-q', '-m', marker_expr],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode == 5:
        return 0
    if proc.returncode != 0:
        raise RuntimeError(
            f'pytest --collect-only -m {marker_expr!r} exited with '
            f'code {proc.returncode}\n'
            f'--- stdout ---\n{proc.stdout}\n'
            f'--- stderr ---\n{proc.stderr}'
        )
    match = _COLLECT_RE.search(proc.stdout)
    if match is None:
        raise RuntimeError(
            f'pytest summary line not found for marker {marker_expr!r}\n'
            f'--- stdout ---\n{proc.stdout}'
        )
    return int(match.group(1))


def write_badge(out_dir: Path, name: str, label: str, count: int) -> Path:
    """Write a shields.io endpoint-badge JSON file.

    Parameters
    ----------
    out_dir : Path
        Directory to write the JSON file into. Must already exist.
    name : str
        Suffix used in the filename ``tests-<name>.json``.
    label : str
        Badge label rendered on the left side of the shield.
    count : int
        Badge message count rendered on the right side of the shield.

    Returns
    -------
    Path
        Path of the written JSON file.
    """
    payload = {
        'schemaVersion': 1,
        'label': label,
        'message': str(count),
        'color': 'blue',
    }
    out_path = out_dir / f'tests-{name}.json'
    out_path.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    return out_path


def main() -> int:
    """Entry point.

    Returns
    -------
    int
        Process exit code (always 0 on success; failures raise).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Directory to write tests-<name>.json badge files into.',
    )
    args = parser.parse_args()
    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    for name, label, expr in _BADGES:
        count = count_tests(expr)
        write_badge(out_dir, name, label, count)
        print(f'{label}: {count}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
