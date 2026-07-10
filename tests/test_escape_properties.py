"""Property-based tests for ``src/zephyrus/escape.py``.

Companion to ``tests/test_escape.py`` holding the Hypothesis-driven property
checks. These live in a separate module because Hypothesis is a develop-extra
dependency: confining the ``importorskip`` to this file keeps the closed-form
pins and the error-contract guards in ``tests/test_escape.py`` running even
when Hypothesis is absent, for example under a ``pip install --no-deps``
image. The physical invariants swept here:

- Positivity / boundedness: the EL rate is strictly positive and finite across
  the physically valid efficiency, flux, radius and mass ranges, and vanishes
  exactly at zero XUV flux.
- Monotonicity / linearity: at fixed geometry the rate is exactly linear in
  the XUV flux, so scaling the flux by ``k`` scales the rate by ``k``.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import numpy as np
import pytest

from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Ms, Re

# hypothesis is a develop-extra dependency; skip this whole module if it is
# unavailable rather than failing collection. The closed-form pins and the
# error-contract guards live in tests/test_escape.py and run unconditionally.
hyp = pytest.importorskip('hypothesis')
given = hyp.given
settings = hyp.settings
st = hyp.strategies

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Reference geometry shared with tests/test_escape.py: Rxuv is chosen distinct
# from Rp (Rxuv = 1.2 * Rp) so the two radius-scaling branches differ.
EPSILON = 0.15
RP = Re
RXUV = 1.2 * Re


@pytest.mark.physics_invariant
@given(
    epsilon=st.floats(min_value=0.1, max_value=0.6),
    fxuv=st.floats(min_value=1e-4, max_value=1e4),
    rp=st.floats(min_value=1e6, max_value=1e8),
    rxuv_factor=st.floats(min_value=1.0, max_value=2.0),
    mp=st.floats(min_value=1e23, max_value=1e27),
)
@settings(max_examples=100, deadline=None, derandomize=True)
def test_el_escape_nonnegative_over_valid_inputs(epsilon, fxuv, rp, rxuv_factor, mp):
    """Escape is strictly positive and finite across physically valid inputs.

    Sweeps the efficiency (literature range 0.1-0.6), diluted XUV flux, radii
    and mass across their physical ranges. Every driver is strictly positive,
    so the rate is strictly positive and finite. The zero-flux limit, where
    the rate must vanish exactly, is checked separately with a fixed input.
    """
    rxuv = rxuv_factor * rp
    val = EL_escape(False, 1.0, 0.0, mp, Ms, epsilon, rp, rxuv, fxuv, scaling=2)
    # Physical drivers are all strictly positive, so the rate is too.
    assert val > 0.0
    # Bounded: the rate stays finite across the swept range (no blow-up).
    assert np.isfinite(val)
    # Limit input: zeroing the XUV flux zeroes the escape rate exactly.
    zero_flux = EL_escape(False, 1.0, 0.0, mp, Ms, epsilon, rp, rxuv, 0.0, scaling=2)
    assert zero_flux == pytest.approx(0.0, abs=1e-30)


@pytest.mark.physics_invariant
@given(
    fxuv=st.floats(min_value=1e-3, max_value=1e4),
    k=st.floats(min_value=0.1, max_value=10.0),
)
@settings(max_examples=100, deadline=None, derandomize=True)
def test_el_escape_linear_in_flux_property(fxuv, k):
    """Scaling the XUV flux by ``k`` scales the escape rate by exactly ``k``.

    Property-based check of linearity across the flux range: the ratio pins
    proportionality independent of the absolute scale, so a spurious offset
    term would fail for some sampled ``fxuv``.
    """
    base = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, fxuv, scaling=2)
    scaled = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, k * fxuv, scaling=2)
    assert scaled == pytest.approx(k * base, rel=1e-12)
    assert base > 0.0
