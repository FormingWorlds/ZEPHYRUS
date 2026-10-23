"""Tests for the data cache steps of ``.github/workflows/nightly.yml``.

The nightly caches the FWL data tree under a key from ``fwl-io cache key``. A
key over other data than the fetch fills, a fetch into another tree than the
tests read, or a check that is skipped or allowed to fail, fails nothing: the
nightly stays green and serves a stale tree. These tests pin the cache steps to
their text, up to the actions/cache version and the outer whitespace of a run,
and their selection to the tracks the suite downloads.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip('yaml')

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO = Path(__file__).parents[1]

# MORS track names and the fwl-mors manifest keys they download.
TRACKS = {'Spada': 'star.tracks.spada_2013', 'Baraffe': 'star.tracks.baraffe_2015'}


SEL, ROOT = '--key star.tracks.spada_2013', '${{ github.workspace }}/fwl_data'


def test_nightly_keys_checks_fetches_and_checks_the_spada_grid():
    """The nightly keys and restores the grid, checks an exact hit, then fills and checks it."""
    workflow = yaml.safe_load(
        (REPO / '.github' / 'workflows' / 'nightly.yml').read_text(encoding='utf-8')
    )
    jobs = list(workflow['jobs'].values())
    steps = [step for job in jobs for step in job['steps']]
    # Fixed step bodies: a shell, an extra command or flag, or another condition would let
    # an unchecked tree be saved, or a damaged exact hit pass.
    i = next(i for i, s in enumerate(steps) if s.get('run', '').startswith('fwl-io cache key'))
    key, cache, hit, fill = (
        {k: v.strip() if k == 'run' else v for k, v in s.items() if k != 'name'}
        for s in steps[i : i + 4]
    )
    assert key == {
        'id': 'cachekey',
        'run': f'fwl-io cache key --prefix fwl-data-nightly- {SEL}',
    }
    out = {k: f'${{{{ steps.cachekey.outputs.{k} }}}}' for k in ('key', 'restore-keys')}
    # A major-version tag or a full commit SHA of actions/cache, never a branch such as main.
    assert re.fullmatch(r'actions/cache@(v[0-9]+|[0-9a-f]{40})', cache.pop('uses', ''))
    assert cache == {'id': 'cache-fwl-data', 'with': {'path': ROOT, **out}}
    assert hit == {
        'if': "steps.cache-fwl-data.outputs.cache-hit == 'true'",
        'run': f'fwl-io cache check {SEL} --data-root {ROOT}',
    }
    assert fill == {
        'run': '\n'.join(
            f'fwl-io cache {c} {SEL} --data-root {ROOT}' for c in ('fetch', 'check')
        )
    }
    # The tests read the tree the steps above filled, after them, on every run.
    tests = next(j for j, s in enumerate(steps) if 'FWL_DATA' in s.get('env', {}))
    assert i + 3 < tests and steps[tests]['env']['FWL_DATA'] == ROOT
    assert not {'if', 'continue-on-error'} & steps[tests].keys()
    for scope in (workflow, *jobs):
        assert not {'defaults', 'continue-on-error', 'if'} & scope.keys()
        assert 'FWL_DATA' not in scope.get('env', {})


def test_cache_selection_covers_every_track_the_suite_downloads():
    """The tracks the tests download are the datasets the nightly keys and fetches."""
    downloaded = set()
    for path in (REPO / 'tests').rglob('*.py'):
        downloaded |= set(
            re.findall(
                r"""DownloadEvolutionTracks\(\s*['"]([^'"]+)['"]""",
                path.read_text(encoding='utf-8'),
            )
        )
    assert downloaded, 'no DownloadEvolutionTracks call site found; has the fetch moved?'
    assert downloaded <= set(TRACKS), f'unknown track names {sorted(downloaded - set(TRACKS))}'
    selected = set(re.findall(r'--key (\S+)', SEL))
    assert {TRACKS[name] for name in downloaded} == selected
