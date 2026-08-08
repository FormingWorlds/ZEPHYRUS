#!/usr/bin/env python3
"""Cache key and restore check for the FWL data tree the nightly caches.

Two subcommands, both used by ``.github/workflows/nightly.yml``::

    python tools/nightly_data_cache.py key
    python tools/nightly_data_cache.py check

ZEPHYRUS reaches one dataset through its ``fwl-mors`` dependency: the Spada
stellar-evolution grid, fetched by ``mors.DownloadEvolutionTracks('Spada')``.
That grid is pinned by a Zenodo record id with an OSF project as its mirror,
both of which live in the installed ``mors.data``, and it unpacks to the
unversioned directory ``stellar_evolution_tracks/Spada``.

``key`` prints ``key=<value>`` for ``GITHUB_OUTPUT``, carrying a digest of the
Zenodo record and the directory the grid unpacks into, so the key moves when
the grid is re-pinned and stays put otherwise. The directory names below are a
fixed namespace segment rather than a resolved property: nothing queries them
from ``mors``, so a mors-side rename of the Spada unpack path has to be
mirrored here by hand.

A moving key is the whole point. ``actions/cache`` writes an entry only on an
exact-key miss and never rewrites one it hits, so a key that never changes is
never rewritten and the tree it holds cannot follow the data.

The OSF project that mirrors the deposit is deliberately not part of the
digest. Its id does not change when the files inside it change, so hashing
it would imply a coverage this key does not have. Mirror drift is untracked,
as it is for the OSF pins in the sibling JANUS module.

The unversioned path is why this matters more here than for a versioned
dataset. ``mors`` skips the download whenever the directory is present, so a
re-pinned grid does not land beside the old one and announce itself: the
stale tree keeps satisfying the check indefinitely. The cache key is the only
thing that can notice.

``check`` verifies the restored grid is actually unpacked. There is no
committed registry for Spada, so the check is structural rather than
per-file: the grid directory exists, holds a plausible number of files, and
no archive is left behind by an interrupted unpack.

Both subcommands fail with a diagnostic rather than degrade: an empty or
partial digest would collide with the workflow's restore-key prefix and
freeze the cached tree with nothing reporting it.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

# JANUS carries a script of the same name for the same actions/cache defect,
# shaped differently because it consumes a manifest-declared dataset with a
# registry and a versioned path. A change here is worth checking against it.
KEY_PREFIX = 'fwl-data-nightly-'
DATASET = 'Spada'
SUBDIR = 'stellar_evolution_tracks'
GRID_DIR = 'fs255_grid'
# The unpacked grid holds about 1600 files across per-composition
# subdirectories. The floor catches a partial restore; it cannot detect a
# fully-unpacked tree that is simply out of date, which is what the key is for.
MIN_FILES = 1000


class ResolutionError(RuntimeError):
    """The pins that decide the cached layout could not be resolved."""


def _pins() -> dict[str, str]:
    """Return the Spada pins read from the installed fwl-mors.

    Returns
    -------
    dict
        Mapping of pin name to value: the Zenodo record and the directory
        the grid unpacks into.

    Raises
    ------
    ResolutionError
        When fwl-mors is absent, no longer exposes these pins, or reports
        no Zenodo record for the dataset.
    """
    try:
        import mors.data
    except ImportError as exc:
        raise ResolutionError(
            'fwl-mors is not importable, so the Spada pins this key tracks cannot '
            'be resolved. Install ZEPHYRUS with its dependencies before resolving '
            'the cache key.'
        ) from exc

    if not hasattr(mors.data, 'get_zenodo_record'):
        raise ResolutionError(
            'the installed fwl-mors exposes no mors.data.get_zenodo_record(), so '
            f'the {DATASET} pin is no longer where this script looks for it. Point '
            'it at whatever now records the Zenodo deposit.'
        )

    record = mors.data.get_zenodo_record(DATASET)
    if not record:
        raise ResolutionError(
            f'the installed fwl-mors reports no Zenodo record for {DATASET!r}, so '
            'this key would track nothing. Check whether the dataset moved into '
            'the fwl-io manifest, which would need a different key entirely.'
        )

    return {'zenodo': str(record), 'subdir': f'{SUBDIR}/{DATASET}'}


def resolve_key() -> str:
    """Return the cache key for the FWL data tree.

    Returns
    -------
    str
        ``fwl-data-nightly-<sha256>``.

    Raises
    ------
    ResolutionError
        When the pins cannot be resolved, or the digest comes out empty and
        would collapse the key onto the restore-key prefix.
    """
    pins = _pins()
    material = [f'{k}\t{pins[k]}' for k in sorted(pins)]
    if not material:
        raise ResolutionError(
            'no pin was resolved, so the key would be the bare restore-key prefix '
            'and the cached tree could never be rewritten.'
        )
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
    base = data_root / SUBDIR / DATASET
    grid = base / GRID_DIR
    problems: list[str] = []

    if not base.is_dir():
        return 0, [f'{base} does not exist']
    if not grid.is_dir():
        problems.append(f'{grid} does not exist, so the grid was never unpacked')

    count = sum(1 for p in grid.rglob('*') if p.is_file()) if grid.is_dir() else 0
    if count < MIN_FILES:
        problems.append(f'{count} files under {grid}, fewer than the {MIN_FILES} expected')

    # mors deletes the tarball after unpacking, so one left behind means the
    # unpack was interrupted and the tree is not what its key describes.
    leftovers = sorted(p.name for p in base.glob('*.tar.gz'))
    if leftovers:
        problems.append(f'archive left unextracted: {", ".join(leftovers)}')

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

    count, problems = check_restored(Path(given))
    print(f'{SUBDIR}/{DATASET}/{GRID_DIR}: {count} files present')
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
        # Anything mors raises reaches here. Name it rather than let a
        # traceback stand in for the diagnostic this script promises.
        print(
            f'error: resolving the {DATASET} pins through fwl-mors failed: {exc!r}. '
            'Check that the installed fwl-mors still exposes the dataset record '
            'this script reads.',
            file=sys.stderr,
        )
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
