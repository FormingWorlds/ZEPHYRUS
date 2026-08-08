"""Tests for ``tools/nightly_data_cache.py``.

The nightly caches the FWL data tree under a key this script resolves. A key
that stops tracking the data does not fail anything: the nightly stays green
and either refetches every run or, worse for ZEPHYRUS, serves a stale grid
forever, because the Spada tracks land in an unversioned directory that
``mors`` skips whenever it is present. These tests pin the key to the
dataset pins in both directions, pin the workflow to the resolved key, and
pin the restore check against a half-unpacked tree.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import importlib.util
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

# Python 3.10 can only resolve an fwl-mors that predates the Spada Zenodo
# record accessor, so the cache key cannot be resolved there at all. On 3.11
# and above a missing accessor is a real failure: it means the dependency
# moved under us, which is exactly what this file exists to notice.
if sys.version_info < (3, 11):
    pytest.skip(
        'this interpreter resolves an fwl-mors that predates the Spada Zenodo '
        'record accessor, so the nightly cache key cannot be resolved here',
        allow_module_level=True,
    )

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO = Path(__file__).parents[1]

# Stand-in the OSF project id is patched to, long enough that it cannot collide
# with a pin value by accident.
MIRROR_SENTINEL = 'osf-mirror-must-not-be-hashed'


def _cache_module():
    """Load the helper, skipping when its one dependency is absent."""
    pytest.importorskip('mors')
    path = REPO / 'tools' / 'nightly_data_cache.py'
    spec = importlib.util.spec_from_file_location('nightly_data_cache', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pin(monkeypatch, *, record):
    """Point the helper at a fabricated Spada Zenodo record."""
    import mors.data

    monkeypatch.setattr(mors.data, 'get_zenodo_record', lambda name: record, raising=True)


def test_cache_key_moves_with_the_spada_pins_and_not_otherwise(monkeypatch):
    """The key tracks the Zenodo record and the unpack directory, and nothing else.

    The grid unpacks to an unversioned directory and mors skips the download
    whenever it exists, so a re-pin is invisible on disk. The key is the only
    thing that can notice, which is why the record has to move it. The OSF
    mirror is deliberately outside the digest, so mirror drift is untracked.
    """
    mod = _cache_module()

    _pin(monkeypatch, record='15729101')
    baseline = mod.resolve_key()

    # Re-pinning the deposit must move the key, or the stale grid is served on.
    _pin(monkeypatch, record='15729102')
    assert mod.resolve_key() != baseline

    # Same pin, same key: a steady-state night has to hit its own entry.
    _pin(monkeypatch, record='15729101')
    assert mod.resolve_key() == baseline

    # The OSF mirror id is deliberately outside the digest: it does not move
    # when the files behind it do, so hashing it would imply coverage.
    import mors.data

    monkeypatch.setattr(mors.data, 'project_id', MIRROR_SENTINEL, raising=False)
    assert mod.resolve_key() == baseline

    # Pin the material itself. The comparison above only notices a digest that
    # reads the patched attribute; one carrying the mirror id as a literal
    # leaves the key at baseline and would pass.
    assert mod._pins() == {'zenodo': '15729101', 'subdir': f'{mod.SUBDIR}/{mod.DATASET}'}

    # An empty digest would leave the key as the constant prefix alone, which
    # exact-hits its own entry on every run and freezes the tree silently.
    assert baseline.startswith(mod.KEY_PREFIX)
    assert re.fullmatch(rf'{re.escape(mod.KEY_PREFIX)}[0-9a-f]{{64}}', baseline)

    # The unpack directory is hashed alongside the record, so a rename of it
    # has to move the key too. Left last: monkeypatch holds the rename in force.
    monkeypatch.setattr(mod, 'SUBDIR', 'stellar_evolution_tracks_v2')
    assert mod.resolve_key() != baseline


def test_cache_key_refuses_to_resolve_without_the_pins(monkeypatch, capsys):
    """A missing record stops the job with a diagnostic rather than a bare key."""
    mod = _cache_module()
    import mors.data

    monkeypatch.setattr(mors.data, 'get_zenodo_record', lambda name: None, raising=True)
    with pytest.raises(mod.ResolutionError, match='no Zenodo record'):
        mod.resolve_key()
    assert mod.main(['key']) == 1
    assert 'no Zenodo record' in capsys.readouterr().err

    monkeypatch.delattr(mors.data, 'get_zenodo_record', raising=True)
    with pytest.raises(mod.ResolutionError, match='get_zenodo_record'):
        mod.resolve_key()


def test_key_command_writes_the_output_line_the_workflow_reads(monkeypatch, tmp_path):
    """The key subcommand emits exactly the GITHUB_OUTPUT line the cache step reads.

    A malformed line leaves the workflow's key expression empty, and the
    cache step rejects an empty key, so the nightly dies at that step
    instead of caching under the resolved key.
    """
    mod = _cache_module()
    _pin(monkeypatch, record='15729101')
    out = tmp_path / 'gh_output'
    monkeypatch.setenv('GITHUB_OUTPUT', str(out))

    assert mod.main(['key']) == 0
    lines = out.read_text(encoding='utf-8').splitlines()
    assert len(lines) == 1
    assert re.fullmatch(rf'key={re.escape(mod.KEY_PREFIX)}[0-9a-f]{{64}}', lines[0])
    assert lines[0].split('=', 1)[1] == mod.resolve_key()


def test_restore_check_requires_an_unpacked_grid(monkeypatch, tmp_path):
    """Presence of the directory is not enough; the grid has to be unpacked.

    A directory that exists but holds nothing, or still holds the archive, is
    exactly what a half-restored cache looks like, so neither may pass.
    """
    mod = _cache_module()
    base = tmp_path / mod.SUBDIR / mod.DATASET

    # A smaller floor keeps the test inside the unit wall-time budget; every
    # assertion below is expressed against the patched constant, so the
    # count-to-floor relationship under test is unchanged.
    monkeypatch.setattr(mod, 'MIN_FILES', 24)

    # Nothing at all.
    count, problems = mod.check_restored(tmp_path)
    assert count == 0
    assert any('does not exist' in p for p in problems)
    assert mod.main(['check', '--data-root', str(tmp_path)]) == 1

    # Present but never unpacked.
    base.mkdir(parents=True)
    count, problems = mod.check_restored(tmp_path)
    assert count == 0
    assert any(mod.GRID_DIR in p for p in problems)

    # Unpacked but short of the floor.
    # Nested, as the real grid is: per-composition subdirectories with no
    # files at the top level, so a walk that does not recurse counts zero.
    grid = base / mod.GRID_DIR
    grid.mkdir()
    for i in range(5):
        comp = grid / f'X0p7{i}_Z0p001_A1p000'
        comp.mkdir()
        (comp / 'track.dat').write_text('x', encoding='utf-8')
    count, problems = mod.check_restored(tmp_path)
    assert count == 5
    assert any('fewer than' in p for p in problems)

    # Populated, but the archive was left behind by an interrupted unpack.
    for i in range(mod.MIN_FILES):
        comp = grid / f'X0p8{i % 20}_Z0p002_A1p875'
        comp.mkdir(exist_ok=True)
        (comp / f'track_{i}.dat').write_text('x', encoding='utf-8')
    (base / 'fs255_grid.tar.gz').write_text('x', encoding='utf-8')
    count, problems = mod.check_restored(tmp_path)
    assert count == mod.MIN_FILES + 5
    assert any('left unextracted' in p for p in problems)
    assert mod.main(['check', '--data-root', str(tmp_path)]) == 1

    # Complete.
    (base / 'fs255_grid.tar.gz').unlink()
    count, problems = mod.check_restored(tmp_path)
    assert problems == []
    assert mod.main(['check', '--data-root', str(tmp_path)]) == 0


def test_check_refuses_to_run_without_a_data_root(monkeypatch):
    """With no root given, check reports it rather than scanning the working directory.

    Path('') is Path('.'), so a guard applied after the Path is built cannot
    fire and would silently count files wherever the job happens to sit.
    """
    mod = _cache_module()
    monkeypatch.delenv('FWL_DATA', raising=False)

    with pytest.raises(mod.ResolutionError, match='no data root'):
        mod._cmd_check(SimpleNamespace(data_root=None))
    assert mod.main(['check']) == 1

    monkeypatch.setenv('FWL_DATA', '')
    with pytest.raises(mod.ResolutionError, match='no data root'):
        mod._cmd_check(SimpleNamespace(data_root=None))


def test_nightly_workflow_derives_its_key_and_declares_no_restore_prefix():
    """The workflow reads the resolved key and offers no prefix fallback."""
    mod = _cache_module()
    workflow = (REPO / '.github' / 'workflows' / 'nightly.yml').read_text(encoding='utf-8')

    assert 'tools/nightly_data_cache.py key' in workflow
    assert 'key: ${{ steps.cachekey.outputs.key }}' in workflow
    # ANY literal, not just the one this replaced: actions/cache never rewrites
    # an entry whose key it hits, so a literal of any value freezes the tree.
    assert not re.search(rf'key:\s*{re.escape(mod.KEY_PREFIX)}\S', workflow)
    # No restore-keys: a prefix fallback would restore the previous grid into
    # the unversioned directory mors checks for, so mors would skip the
    # download and the stale tree would be saved under the new key.
    assert not re.search(r'^\s*restore-keys:', workflow, re.MULTILINE)
    # The check only means something on an exact hit.
    assert "if: steps.cache-fwl-data.outputs.cache-hit == 'true'" in workflow


def test_cache_script_covers_every_dataset_the_suite_downloads():
    """The dataset the key tracks is the dataset the tests actually fetch.

    Nothing but convention links the constant in the helper to the call sites
    here, so a second dataset added later would go untracked by the cache key
    with nothing failing.
    """
    mod = _cache_module()
    downloaded = set()
    for path in (REPO / 'tests').rglob('test_*.py'):
        for name in re.findall(
            r"DownloadEvolutionTracks\(\s*'([^']+)'", path.read_text(encoding='utf-8')
        ):
            downloaded.add(name)

    assert downloaded, 'no DownloadEvolutionTracks call site found; has the fetch moved?'
    assert downloaded == {mod.DATASET}, (
        f'the suite downloads {sorted(downloaded)} but the cache key tracks '
        f'only {mod.DATASET!r}; an untracked dataset freezes in the cache'
    )
