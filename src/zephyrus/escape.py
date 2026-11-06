"""
!!! info "`escape.py`"
    Released entry point for energy-limited atmospheric escape.<br>
    Authors: Emma Postolec, Harrison Nicholls
"""

# The energy-limited physics lives in zephyrus.hydrodynamic beside the other
# hydrodynamic rate prescriptions; this module re-exports the released entry
# point so the import path `from zephyrus.escape import EL_escape` and the
# names historically bound into this namespace keep working for existing
# callers, including the PROTEUS version pin.

from zephyrus.constants import *  # noqa: F403
from zephyrus.hydrodynamic import EL_escape as EL_escape
from zephyrus.planets_parameters import *  # noqa: F403
