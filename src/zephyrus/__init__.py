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
from zephyrus import dispatcher as dispatcher  # noqa: E402
from zephyrus import escape as escape  # noqa: E402

# The regime-dispatcher entry points, re-exported at the package top level.
from zephyrus.dispatcher import (  # noqa: E402
    DispatchSettings as DispatchSettings,
)
from zephyrus.dispatcher import (  # noqa: E402
    EscapeInputs as EscapeInputs,
)
from zephyrus.dispatcher import (  # noqa: E402
    EscapeResult as EscapeResult,
)
from zephyrus.dispatcher import (  # noqa: E402
    dispatch as dispatch,
)
