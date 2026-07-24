from __future__ import annotations

try:
    from ._version import __version__, __version_tuple__
except ImportError:
    # Fallback for when the package is not installed (e.g., running from
    # source without setuptools-scm having generated _version.py).
    __version__ = '0.0.0.dev0'
    __version_tuple__ = (0, 0, 0, 'dev0')

# Submodules re-exported so `import zephyrus` exposes the package API.
from zephyrus import collision as collision  # noqa: E402
