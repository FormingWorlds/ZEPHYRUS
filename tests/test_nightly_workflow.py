"""Tests for the data cache steps of ``.github/workflows/nightly.yml``.

The nightly caches the FWL data tree under a key from ``fwl-io cache key``. A
key over other data than the fetch fills, or a fetch into another tree than the
tests read, fails nothing: the nightly stays green and serves a stale tree. These
tests pin the steps to one selection and one data root, and the selection to
the tracks the suite downloads.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import re
from pathlib import Path

import pytest
import yaml

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO = Path(__file__).parents[1]

# MORS track names and the fwl-mors manifest keys they download.
TRACKS = {'Spada': 'star.tracks.spada_2013', 'Baraffe': 'star.tracks.baraffe_2015'}


def _steps():
    """Return the steps of the nightly in order, with the fwl-io cache calls."""
    workflow = yaml.safe_load(
        (REPO / '.github' / 'workflows' / 'nightly.yml').read_text(encoding='utf-8')
    )
    steps = [step for job in workflow['jobs'].values() for step in job['steps']]
    calls = []
    for i, step in enumerate(steps):
        for line in step.get('run', '').splitlines():
            m = re.match(r'\s*fwl-io cache (key|fetch|check)\b(.*)', line)
            if m:
                calls.append((m[1], i, m[2]))
    return steps, calls


def _selection(args):
    """Return the ``--provider`` and ``--key`` options of one call, in order."""
    return re.findall(r'--(provider|key) (\S+)', args)


def test_nightly_keys_fetches_and_checks_one_selection_in_one_tree():
    """Key, restore, fetch and check run in order on one selection and one root."""
    steps, calls = _steps()
    caches = [
        i for i, step in enumerate(steps) if step.get('uses', '').startswith('actions/cache@')
    ]
    tests = next(i for i, step in enumerate(steps) if 'pytest' in step.get('run', ''))
    assert [name for name, _, _ in calls] == ['key', 'fetch', 'check']
    assert len(caches) == 1
    assert calls[0][1] < caches[0] < calls[1][1] <= calls[2][1] < tests

    key_id = steps[calls[0][1]]['id']
    cache = steps[caches[0]]['with']
    assert cache['key'] == f'${{{{ steps.{key_id}.outputs.key }}}}'
    assert cache['restore-keys'] == f'${{{{ steps.{key_id}.outputs.restore-keys }}}}'
    assert {tuple(_selection(args)) for _, _, args in calls} == {
        (('key', 'star.tracks.spada_2013'),)
    }
    root = cache['path']
    assert steps[tests]['env']['FWL_DATA'] == root
    for name, _, args in calls[1:]:
        assert args.split('--data-root ')[1].strip() == root, name


def test_cache_selection_covers_every_track_the_suite_downloads():
    """The tracks the tests download are the datasets the nightly keys and fetches."""
    downloaded = set()
    for path in (REPO / 'tests').rglob('test_*.py'):
        downloaded |= set(
            re.findall(
                r"DownloadEvolutionTracks\(\s*'([^']+)'", path.read_text(encoding='utf-8')
            )
        )
    assert downloaded, 'no DownloadEvolutionTracks call site found; has the fetch moved?'
    assert downloaded <= set(TRACKS), f'unknown track names {sorted(downloaded - set(TRACKS))}'

    _, calls = _steps()
    assert [name for name, _, _ in calls][:1] == ['key'], 'no fwl-io cache key step'
    selected = {value for kind, value in _selection(calls[0][2]) if kind == 'key'}
    assert {TRACKS[name] for name in downloaded} == selected
