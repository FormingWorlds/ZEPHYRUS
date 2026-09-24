"""Tests for ``tools/nightly_data_cache.py``.

The nightly caches the FWL data tree under a key this script resolves. A key
that stops tracking the data does not fail anything: the nightly stays green
and either refetches every run or serves a stale tree forever, because
``actions/cache`` never rewrites an entry whose key it hits. These tests pin
the key to the Spada entry of the fwl-mors dataset manifest in both
directions, pin the workflow to the resolved key, and pin the restore check
against a half-unpacked tree.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import importlib.util
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

REPO = Path(__file__).parents[1]

SPADA_RECORD = '15729101'
SPADA_MD5 = 'f76987cf3d1da50435547f44a484a97f'
BARAFFE_RECORD = '15729114'


def _cache_module():
    """Load the helper, skipping when its dependencies are absent."""
    pytest.importorskip('mors')
    pytest.importorskip('fwl_io')
    path = REPO / 'tools' / 'nightly_data_cache.py'
    spec = importlib.util.spec_from_file_location('nightly_data_cache', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest(
    monkeypatch,
    tmp_path,
    *,
    spada_record=SPADA_RECORD,
    spada_md5=SPADA_MD5,
    baraffe=BARAFFE_RECORD,
):
    """Write a fwl-mors-style manifest and point ``mors.data.manifest_path`` at it."""
    import mors.data

    root = tmp_path / 'manifest'
    root.mkdir(exist_ok=True)
    manifest = root / 'mors_manifest.toml'
    manifest.write_text(
        '[star.tracks.baraffe_2015]\n'
        'name = "Baraffe"\n'
        f'zenodo = "10.5281/zenodo.{baraffe}"\n'
        'required_by = ["mors"]\n\n'
        '[star.tracks.spada_2013]\n'
        'name = "Spada"\n'
        f'zenodo = "10.5281/zenodo.{spada_record}"\n'
        'extract = "tar"\n'
        'required_by = ["mors"]\n',
        encoding='utf-8',
    )
    (root / 'star.tracks.baraffe_2015.registry.txt').write_text(
        'BHAC15_tracks.dat md5:4b5c14255ca845880a2e940110861df1\n', encoding='utf-8'
    )
    (root / 'star.tracks.spada_2013.registry.txt').write_text(
        f'fs255_grid.tar.gz md5:{spada_md5}\n', encoding='utf-8'
    )
    monkeypatch.setattr(mors.data, 'manifest_path', lambda: manifest, raising=False)
    return manifest


def test_cache_key_moves_with_the_spada_manifest_entry_and_not_otherwise(monkeypatch, tmp_path):
    """The key tracks the Spada record and checksums, and nothing else in the manifest.

    A re-pinned record or a changed archive checksum must move the key, or the
    stale tree keeps being restored. The Baraffe entry sits in the same
    manifest but ZEPHYRUS does not fetch it, so changing it must leave the
    key alone.
    """
    mod = _cache_module()

    _manifest(monkeypatch, tmp_path)
    baseline = mod.resolve_key()

    # A different record moves the key (the version directory changes name).
    _manifest(monkeypatch, tmp_path, spada_record='15729102')
    record_key = mod.resolve_key()
    assert record_key != baseline

    # Same record, different archive checksum: the key still moves.
    _manifest(monkeypatch, tmp_path, spada_md5='0' * 32)
    digest_key = mod.resolve_key()
    assert digest_key not in (baseline, record_key)

    # Same pins, same key: a steady-state night has to hit its own entry.
    _manifest(monkeypatch, tmp_path)
    assert mod.resolve_key() == baseline

    # The Baraffe entry is not consumed here, so it must not move the key.
    _manifest(monkeypatch, tmp_path, baraffe='15729199')
    assert mod.resolve_key() == baseline

    # An empty digest would leave the key as the constant prefix alone, which
    # exact-hits its own entry on every run and freezes the tree silently.
    assert re.fullmatch(rf'{re.escape(mod.KEY_PREFIX)}[0-9a-f]{{64}}', baseline)


def test_cache_key_material_is_the_versioned_directory_and_the_registry(monkeypatch, tmp_path):
    """The digest is built from the directory fwl-io places the grid in and the registry.

    The two-keys comparison above only proves the key moves; pin what it moves
    with, so a digest of a constant or of the wrong dataset cannot pass.
    """
    mod = _cache_module()
    _manifest(monkeypatch, tmp_path)

    fetcher = mod._fetcher(tmp_path)

    assert fetcher.rel_dir == f'star/tracks/spada_2013/r{SPADA_RECORD}'
    assert fetcher.registry == {'fs255_grid.tar.gz': f'md5:{SPADA_MD5}'}
    assert fetcher.target_dir == tmp_path / fetcher.rel_dir


def test_cache_key_refuses_to_resolve_without_the_pins(monkeypatch, tmp_path, capsys):
    """A missing entry, manifest, registry or checksum stops the job with a diagnostic.

    Each is a distinct way the key could quietly track nothing, so each has
    its own message and none may return a bare key.
    """
    mod = _cache_module()
    import mors.data

    # The manifest no longer declares the Spada entry.
    manifest = _manifest(monkeypatch, tmp_path)
    manifest.write_text(
        '[star.tracks.baraffe_2015]\nname = "Baraffe"\nzenodo = "10.5281/zenodo.15729114"\n',
        encoding='utf-8',
    )
    with pytest.raises(
        mod.ResolutionError, match='star.tracks.spada_2013.*star.tracks.baraffe'
    ):
        mod.resolve_key()
    assert mod.main(['key']) == 1
    assert 'declares no' in capsys.readouterr().err

    # The registry file is gone.
    manifest = _manifest(monkeypatch, tmp_path)
    (manifest.parent / 'star.tracks.spada_2013.registry.txt').unlink()
    with pytest.raises(mod.ResolutionError, match='registry'):
        mod.resolve_key()

    # The registry is empty, so there is no checksum to track; fwl-io refuses it.
    manifest = _manifest(monkeypatch, tmp_path)
    (manifest.parent / 'star.tracks.spada_2013.registry.txt').write_text('', encoding='utf-8')
    with pytest.raises(mod.ResolutionError, match='empty registry'):
        mod.resolve_key()

    # The manifest file itself is gone.
    manifest = _manifest(monkeypatch, tmp_path)
    manifest.unlink()
    with pytest.raises(mod.ResolutionError, match='does not exist'):
        mod.resolve_key()

    # fwl-mors no longer exposes a manifest at all.
    monkeypatch.delattr(mors.data, 'manifest_path', raising=True)
    with pytest.raises(mod.ResolutionError, match='manifest_path'):
        mod.resolve_key()


def test_key_command_writes_the_output_line_the_workflow_reads(monkeypatch, tmp_path):
    """The key subcommand emits exactly the GITHUB_OUTPUT line the cache step reads.

    A malformed line leaves the workflow's key expression empty, and the
    cache step rejects an empty key, so the nightly dies at that step
    instead of caching under the resolved key.
    """
    mod = _cache_module()
    _manifest(monkeypatch, tmp_path)
    out = tmp_path / 'gh_output'
    monkeypatch.setenv('GITHUB_OUTPUT', str(out))

    assert mod.main(['key']) == 0
    lines = out.read_text(encoding='utf-8').splitlines()
    assert len(lines) == 1
    assert re.fullmatch(rf'key={re.escape(mod.KEY_PREFIX)}[0-9a-f]{{64}}', lines[0])
    assert lines[0].split('=', 1)[1] == mod.resolve_key()


def test_restore_check_requires_an_unpacked_grid(monkeypatch, tmp_path):
    """Presence of the directory is not enough; the grid has to be unpacked.

    A version directory that exists but holds nothing, or a grid holding a
    handful of files, is what a half-restored cache looks like, so neither may
    pass. The grid is looked for below the versioned directory fwl-io derives
    from the manifest, so a tree at the old unversioned path counts as absent.
    """
    mod = _cache_module()
    _manifest(monkeypatch, tmp_path)
    base = tmp_path / 'star' / 'tracks' / 'spada_2013' / f'r{SPADA_RECORD}'

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

    assert mod.main(['check', '--data-root', str(tmp_path)]) == 1

    # The old unversioned location is not the grid this key describes.
    old = tmp_path / 'stellar_evolution_tracks' / 'Spada' / mod.GRID_DIR
    old.mkdir(parents=True)
    for i in range(mod.MIN_FILES):
        (old / f'track_{i}.dat').write_text('x', encoding='utf-8')
    assert mod.check_restored(tmp_path)[0] == 5

    # Complete.
    for i in range(mod.MIN_FILES):
        comp = grid / f'X0p8{i % 20}_Z0p002_A1p875'
        comp.mkdir(exist_ok=True)
        (comp / f'track_{i}.dat').write_text('x', encoding='utf-8')
    count, problems = mod.check_restored(tmp_path)
    assert count == mod.MIN_FILES + 5
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
    # No restore-keys: a prefix fallback would restore the previous record's
    # tree, report no cache hit and save that stale tree under the new key.
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
