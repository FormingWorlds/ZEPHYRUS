#!/usr/bin/env python3
"""Cache key and restore check for the FWL data tree the nightly caches.

Two subcommands, both used by ``.github/workflows/nightly.yml``::

    python tools/nightly_data_cache.py key
    python tools/nightly_data_cache.py check

ZEPHYRUS reaches one dataset through its ``fwl-mors`` dependency: the Spada
stellar-evolution grid, fetched by ``mors.DownloadEvolutionTracks('Spada')``.
fwl-mors declares that grid in its shipped dataset manifest
(``mors.data.manifest_path()``), which fwl-io reads: a Zenodo version DOI, a
registry of checksums, and the versioned directory
``star/tracks/spada_2013/r<record-id>`` the archive unpacks into.

``key`` prints ``key=<value>`` for ``GITHUB_OUTPUT``, carrying a digest of that
directory and the registry checksums the manifest pins for the dataset, both
resolved through fwl-io. The key therefore moves when the grid is re-pinned and
stays put otherwise. The manifest key of the Spada entry below is a fixed
namespace segment rather than a resolved property: a rename of the entry on the
fwl-mors side has to be mirrored here by hand, and until it is, the script
stops with a diagnostic naming the entries it found.

A moving key is the whole point. ``actions/cache`` writes an entry only on an
exact-key miss and never rewrites one it hits, so a key that never changes is
never rewritten and the tree it holds cannot follow the data.

Only the Spada entry is hashed. The Baraffe entry in the same manifest is not
fetched by ZEPHYRUS, so a change to it must not empty the cache.

``check`` verifies the restored grid is actually unpacked. The registry pins
the archive, not the extracted members, so the check is structural rather than
per-file: the versioned directory and the grid directory exist and the grid
holds a plausible number of files.

Both subcommands fail with a diagnostic rather than degrade: an empty or
partial digest would leave the key as the constant prefix alone, and a
constant key exact-hits its own entry on every run, freezing the cached
tree with nothing reporting it.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import tempfile
from pathlib import Path

# JANUS carries a script of the same name for the same actions/cache defect. It
# hashes every dataset of the manifest where this one hashes only the Spada
# entry. A change here is worth checking against it.
KEY_PREFIX = 'fwl-data-nightly-'
DATASET = 'Spada'
MANIFEST_KEY = 'star.tracks.spada_2013'
GRID_DIR = 'fs255_grid'
# The unpacked grid holds about 1600 files across per-composition
# subdirectories. The floor catches a partial restore; it cannot detect a
# fully-unpacked tree that is simply out of date, which is what the key is for.
MIN_FILES = 1000


class ResolutionError(RuntimeError):
    """The pins that decide the cached layout could not be resolved."""


def _manifest_path() -> Path:
    """Return the dataset manifest shipped inside the installed fwl-mors.

    Returns
    -------
    Path
        Absolute path of ``mors_manifest.toml``.

    Raises
    ------
    ResolutionError
        When fwl-mors or fwl-io is absent, fwl-mors exposes no manifest, or
        the manifest file it names does not exist.
    """
    try:
        import fwl_io  # noqa: F401
        import mors.data
    except ImportError as exc:
        raise ResolutionError(
            f'fwl-mors and fwl-io are not both importable ({exc}), so the {DATASET} '
            'pins this key tracks cannot be resolved. Install ZEPHYRUS with its '
            'dependencies before resolving the cache key.'
        ) from exc

    if not hasattr(mors.data, 'manifest_path'):
        raise ResolutionError(
            'the installed fwl-mors exposes no mors.data.manifest_path(), so the '
            f'{DATASET} pin is no longer where this script looks for it. Point it at '
            'whatever now declares the Zenodo record and checksums.'
        )

    path = Path(mors.data.manifest_path())
    if not path.is_file():
        raise ResolutionError(f'the dataset manifest fwl-mors names does not exist: {path}')
    return path


def _fetcher(data_root: Path):
    """Build the fwl-io fetcher for the Spada entry of the fwl-mors manifest.

    Parameters
    ----------
    data_root : Path
        Root the fetcher resolves its target directory below.

    Returns
    -------
    fwl_io.Fetcher
        Fetcher whose ``rel_dir``, ``target_dir`` and ``registry`` describe
        the dataset. Building it does not touch the network but creates
        ``data_root`` when it is absent.

    Raises
    ------
    ResolutionError
        When the manifest cannot be located or read, declares no Spada entry,
        has no registry file for it, or fwl-io refuses the entry.
    """
    from fwl_io import create_fetcher, load_manifest

    manifest = _manifest_path()
    try:
        datasets = {ds.key: ds for ds in load_manifest(manifest)}
    except ValueError as exc:
        raise ResolutionError(f'fwl-io could not read {manifest}: {exc}') from exc
    ds = datasets.get(MANIFEST_KEY)
    if ds is None:
        raise ResolutionError(
            f'{manifest} declares no {MANIFEST_KEY!r} entry (it declares '
            f'{sorted(datasets) or "no dataset"}), so this key would track nothing. '
            'Check whether the entry was renamed and update MANIFEST_KEY.'
        )
    if not ds.registry_path.is_file():
        raise ResolutionError(
            f'dataset {ds.key!r} declares a registry at {ds.registry_path}, which does '
            'not exist. The checksums are half of what this key tracks, so resolving it '
            'without them would freeze the cache.'
        )
    try:
        return create_fetcher(
            subdir=ds.subdir,
            zenodo=ds.zenodo,
            registry=ds.registry_path,
            data_root=data_root,
            extract=ds.extract,
        )
    except ValueError as exc:
        # fwl-io refuses, among others, an empty registry, which would leave
        # nothing but a directory name for the key to track.
        raise ResolutionError(f'fwl-io rejected the {MANIFEST_KEY!r} entry: {exc}') from exc


def resolve_key() -> str:
    """Return the cache key for the FWL data tree.

    Returns
    -------
    str
        ``fwl-data-nightly-<sha256>``, a digest of the version directory and the
        registry checksums of the Spada manifest entry.

    Raises
    ------
    ResolutionError
        When the pins cannot be resolved. The digest always covers a directory
        and at least one checksum, so it cannot collapse to ``KEY_PREFIX``.
    """
    with tempfile.TemporaryDirectory() as tmp:
        fetcher = _fetcher(Path(tmp))
    material = [f'dir\t{fetcher.rel_dir}']
    material += [f'file\t{name}\t{fetcher.registry[name]}' for name in sorted(fetcher.registry)]
    digest = hashlib.sha256('\n'.join(material).encode('utf-8')).hexdigest()
    return f'{KEY_PREFIX}{digest}'


def check_restored(data_root: Path) -> tuple[int, list[str]]:
    """Report how completely the Spada grid is restored below ``data_root``.

    Parameters
    ----------
    data_root : Path
        Root of the restored FWL data tree.

    Returns
    -------
    tuple
        The file count found under the grid directory, and a list of
        problems; an empty list means the tree looks complete.
    """
    base = _fetcher(data_root).target_dir
    grid = base / GRID_DIR
    problems: list[str] = []

    if not base.is_dir():
        return 0, [f'{base} does not exist']
    if not grid.is_dir():
        problems.append(f'{grid} does not exist, so the grid was never unpacked')

    count = sum(1 for p in grid.rglob('*') if p.is_file()) if grid.is_dir() else 0
    if count < MIN_FILES:
        problems.append(f'{count} files under {grid}, fewer than the {MIN_FILES} expected')

    return count, problems


def _cmd_key(args: argparse.Namespace) -> int:
    key = resolve_key()
    print(f'Cache key: {key}', file=sys.stderr)
    output = os.environ.get('GITHUB_OUTPUT')
    if output:
        with open(output, 'a', encoding='utf-8') as handle:
            handle.write(f'key={key}\n')
    else:
        print(f'key={key}')
    return 0


def _cmd_check(args: argparse.Namespace) -> int:
    # Test the argument before it becomes a Path: Path('') is Path('.'), so a
    # guard on the Path would pass and quietly check the working directory.
    given = args.data_root or os.environ.get('FWL_DATA')
    if not given:
        raise ResolutionError('no data root to check: pass --data-root or set FWL_DATA.')

    if not Path(given).is_dir():
        # Building the fetcher creates a missing root, and a check must not.
        raise ResolutionError(f'the data root {given} does not exist, so nothing was restored.')
    count, problems = check_restored(Path(given))
    print(f'{DATASET} grid: {count} files present')
    if problems:
        for p in problems:
            print(f'  {p}', file=sys.stderr)
        print(
            'The restored tree does not match the key it was stored under. An '
            'exact-key hit is never re-saved, so delete that cache entry to let '
            'the next run store a complete tree.',
            file=sys.stderr,
        )
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run a subcommand and return its exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('key', help='print the cache key for the FWL data tree')
    check = sub.add_parser('check', help='verify the restored grid is unpacked')
    check.add_argument('--data-root', default=None, help='defaults to FWL_DATA')

    args = parser.parse_args(argv)
    handler = {'key': _cmd_key, 'check': _cmd_check}[args.command]
    try:
        return handler(args)
    except ResolutionError as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 -- reported, never swallowed
        # Anything fwl-mors or fwl-io raises reaches here. Name it rather than
        # let a traceback stand in for the diagnostic this script promises.
        print(
            f'error: resolving the {DATASET} pins through fwl-mors and fwl-io failed: '
            f'{exc!r}. If either package changed its API, update this script.',
            file=sys.stderr,
        )
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
