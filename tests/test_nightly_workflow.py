"""Tests for the data cache steps of ``.github/workflows/nightly.yml``.

The nightly caches the FWL data tree under a key from ``fwl-io cache key``. A
key over other data than the fetch fills, a fetch into another tree than the
tests read, or a check that is skipped or allowed to fail, fails nothing: the
nightly stays green and serves a stale tree. These tests pin the steps to one
selection and one data root, run unconditionally, and the selection to the
tracks the suite downloads.

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


def _workflow():
    """Return the nightly jobs, their steps in order and the fwl-io cache calls by name."""
    workflow = yaml.safe_load(
        (REPO / '.github' / 'workflows' / 'nightly.yml').read_text(encoding='utf-8')
    )
    jobs = list(workflow['jobs'].values())
    steps = [step for job in jobs for step in job['steps']]
    calls = {}
    for i, step in enumerate(steps):
        for line in step.get('run', '').splitlines():
            m = re.match(r'\s*fwl-io cache (key|fetch|check)\b(.*)', line)
            if m:
                # One call each, on one line, whose exit status the step returns.
                assert m[1] not in calls and not re.search(r'[|;&]|\\$', m[2]), line
                calls[m[1]] = (i, m[2])
    return jobs, steps, calls


def _selection(args):
    """Return the ``--provider`` and ``--key`` options of one call, in order."""
    return re.findall(r'--(provider|key) (\S+)', args)


def test_nightly_keys_fetches_and_checks_one_selection_in_one_tree():
    """Key, restore, fetch and check run every night in order, on one selection and root."""
    jobs, steps, calls = _workflow()
    assert list(calls) == ['key', 'fetch', 'check']
    key_id = steps[calls['key'][0]]['id']
    out = {k: f'${{{{ steps.{key_id}.outputs.{k} }}}}' for k in ('key', 'restore-keys')}
    caches = [
        i
        for i, s in enumerate(steps)
        if s.get('uses', '').startswith('actions/cache@') and s['with'].get('key') == out['key']
    ]
    tests = next(i for i, s in enumerate(steps) if 'pytest' in s.get('run', ''))
    assert len(caches) == 1
    cache = steps[caches[0]]['with']
    assert calls['key'][0] < caches[0] < calls['fetch'][0] <= calls['check'][0] < tests
    # Both keys from the key step and nothing else: lookup-only would skip the save.
    assert cache == {'path': cache['path'], **out}
    # A skipped or failure-tolerant step leaves an unchecked tree saved and the run green.
    for i in (calls['key'][0], caches[0], calls['fetch'][0], calls['check'][0]):
        assert not {'if', 'continue-on-error'} & steps[i].keys(), steps[i].get('name')
    assert 'env' not in steps[calls['fetch'][0]] and 'env' not in steps[calls['check'][0]]
    assert not any('continue-on-error' in job for job in jobs)

    assert re.findall(r'--prefix (\S+)', calls['key'][1]) == ['fwl-data-nightly-']
    assert {tuple(_selection(args)) for _, args in calls.values()} == {
        (('key', 'star.tracks.spada_2013'),)
    }
    root = cache['path']
    assert steps[tests]['env']['FWL_DATA'] == root
    for name in ('fetch', 'check'):
        args = calls[name][1]
        assert args.count('--data-root') == 1, name
        assert re.search(rf'--data-root {re.escape(root)}(\s|$)', args), name


def test_cache_selection_covers_every_track_the_suite_downloads():
    """The tracks the tests download are the datasets the nightly keys and fetches."""
    downloaded = set()
    for path in (REPO / 'tests').rglob('test_*.py'):
        downloaded |= set(
            re.findall(
                r"""DownloadEvolutionTracks\(\s*['"]([^'"]+)['"]""",
                path.read_text(encoding='utf-8'),
            )
        )
    assert downloaded, 'no DownloadEvolutionTracks call site found; has the fetch moved?'
    assert downloaded <= set(TRACKS), f'unknown track names {sorted(downloaded - set(TRACKS))}'

    _, _, calls = _workflow()
    assert 'key' in calls, 'no fwl-io cache key step'
    selected = {value for kind, value in _selection(calls['key'][1]) if kind == 'key'}
    assert {TRACKS[name] for name in downloaded} == selected
