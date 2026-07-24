"""Property-based tests for ``src/zephyrus/collision.py``.

Companion to ``tests/test_collision.py`` holding the Hypothesis-driven
property checks, in a separate module so the ``importorskip`` keeps the
closed-form pins running when Hypothesis is absent (for example under a
``pip install --no-deps`` image). The physical invariants swept here:

- Boundedness: the loss fraction stays in [0, 1] and finite across the
  physically valid mass, radius, density, speed, and angle ranges.
- Monotonicity: at fixed geometry the loss never decreases with contact
  speed and never increases with impact parameter.
- Reduction: at equal bulk densities the density-weighted interacting
  mass of Eqn. B1 equals the interacting volume of Eqn. B2 exactly.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import numpy as np
import pytest

from zephyrus.collision import mass_loss
from zephyrus.constants import G

# hypothesis is a develop-extra dependency; skip this whole module if it is
# unavailable rather than failing collection. The closed-form pins and the
# error-contract guards live in tests/test_collision.py and run
# unconditionally.
hyp = pytest.importorskip('hypothesis')
given = hyp.given
settings = hyp.settings
st = hyp.strategies

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Strategy bounds: masses from large asteroids to super-Earths, radii and
# densities spanning icy to iron bodies, speeds from rest to twice the
# paper's fitted ceiling so the cap branch is exercised, and the full
# head-on to grazing angle range.
_MASS = st.floats(min_value=1e22, max_value=2e25)
_RADIUS = st.floats(min_value=1e6, max_value=2e7)
_RHO = st.floats(min_value=900.0, max_value=13000.0)
_B = st.floats(min_value=0.0, max_value=1.0)
_VFAC = st.floats(min_value=0.0, max_value=6.0)


def _vesc(m_t, m_i, r_t, r_i):
    """Mutual escape speed at contact, matching the paper's definition."""
    return np.sqrt(2.0 * G * (m_t + m_i) / (r_t + r_i))


@pytest.mark.physics_invariant
@given(m_i=_MASS, m_t=_MASS, rho_i=_RHO, rho_t=_RHO, r_i=_RADIUS, r_t=_RADIUS, b=_B, vfac=_VFAC)
@settings(max_examples=200, deadline=None, derandomize=True)
def test_loss_fraction_is_bounded_over_the_physical_domain(
    m_i, m_t, rho_i, rho_t, r_i, r_t, b, vfac
):
    """The loss fraction is finite and in [0, 1] everywhere in the domain.

    Sweeps independent masses, radii, densities, angles, and speeds up to
    twice the fitted ceiling. The cap at 1, the grazing zero, and the
    clamped interacting mass together guarantee boundedness even in the
    corners the fit was never constrained in (a dense small impactor
    head-on), where the raw common-height cap bookkeeping alone would
    misbehave.
    """
    v_c = vfac * _vesc(m_t, m_i, r_t, r_i)
    x = mass_loss(v_c, m_i, m_t, rho_i, rho_t, r_i, r_t, b)
    assert np.isfinite(x)
    assert 0.0 <= x <= 1.0
    # Grazing anchor inside the sweep: at b = 1 the loss is exactly zero.
    assert mass_loss(v_c, m_i, m_t, rho_i, rho_t, r_i, r_t, 1.0) == pytest.approx(
        0.0, abs=1e-15
    )


@pytest.mark.physics_invariant
@given(m_i=_MASS, m_t=_MASS, rho_i=_RHO, rho_t=_RHO, r_i=_RADIUS, r_t=_RADIUS, b=_B)
@settings(max_examples=200, deadline=None, derandomize=True)
def test_loss_never_decreases_with_contact_speed(m_i, m_t, rho_i, rho_t, r_i, r_t, b):
    """Hitting harder never erodes less, up to saturation at the cap.

    At fixed geometry the bracket scales as v_c squared, so the loss is
    non-decreasing in speed; equality is admitted only where both speeds
    saturate the total-erosion cap. The zero-speed anchor is strict: no
    kinetic energy, no erosion.
    """
    v_esc = _vesc(m_t, m_i, r_t, r_i)
    x_slow = mass_loss(1.0 * v_esc, m_i, m_t, rho_i, rho_t, r_i, r_t, b)
    x_fast = mass_loss(2.5 * v_esc, m_i, m_t, rho_i, rho_t, r_i, r_t, b)
    assert x_fast >= x_slow
    assert mass_loss(0.0, m_i, m_t, rho_i, rho_t, r_i, r_t, b) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.physics_invariant
@given(m_i=_MASS, m_t=_MASS, rho=_RHO, r_i=_RADIUS, r_t=_RADIUS)
@settings(max_examples=200, deadline=None, derandomize=True)
def test_loss_never_increases_toward_grazing_at_equal_density(m_i, m_t, rho, r_i, r_t):
    """A more grazing impact never erodes more than a more direct one.

    At equal bulk densities the interacting mass reduces to the
    interacting volume, whose angle factor (1-b)^2 (1+2b) falls
    monotonically over [0, 1], so the loss is non-increasing from head-on
    to grazing; equality is admitted only where the cap saturates both.
    The sweep is restricted to equal densities because the published
    common-height cap bookkeeping is not monotonic in b for density
    contrasts far outside the fitted range (a dense small impactor on a
    low-density target), which the boundedness sweep still covers.
    """
    v_c = 2.0 * _vesc(m_t, m_i, r_t, r_i)
    losses = [mass_loss(v_c, m_i, m_t, rho, rho, r_i, r_t, b) for b in (0.0, 0.4, 0.8, 1.0)]
    assert all(a >= z for a, z in zip(losses, losses[1:]))
    assert losses[-1] == pytest.approx(0.0, abs=1e-15)


@pytest.mark.physics_invariant
@given(
    m_t=_MASS,
    q=st.floats(min_value=1e-4, max_value=1.0),
    rho=_RHO,
    r_i=_RADIUS,
    r_t=_RADIUS,
    b=_B,
    vfac=st.floats(min_value=0.2, max_value=1.5),
)
@settings(max_examples=200, deadline=None, derandomize=True)
def test_equal_densities_match_the_interacting_volume_form(m_t, q, rho, r_i, r_t, b, vfac):
    """Eqn. B1 collapses to Eqn. B2 exactly whenever the densities match.

    The density-weighted interacting mass must reduce to the closed-form
    interacting volume for equal bulk densities at every geometry and
    angle, including head-on with dissimilar radii where the individual
    common-height caps leave [0, V] and only their sum is meaningful.
    The impactor is drawn as a mass ratio of the target (an impactor is
    the lighter body, matching the paper's suite), which bounds
    M_i/M_tot at 1/2 and the raw law at 0.64 * (2.25 / sqrt(2))^0.65 =
    0.867, so the equality check is never a saturated 1 == 1.
    """
    m_i = q * m_t
    v_c = vfac * _vesc(m_t, m_i, r_t, r_i)
    x = mass_loss(v_c, m_i, m_t, rho, rho, r_i, r_t, b)

    f_v = 0.25 * ((r_t + r_i) ** 3 / (r_t**3 + r_i**3)) * (1 - b) ** 2 * (1 + 2 * b)
    bracket = vfac**2 * (m_i / (m_i + m_t)) ** 0.5 * f_v
    x_expected = 0.64 * bracket**0.65
    assert x_expected < 1.0  # the mass-ratio bound keeps the law sub-cap
    assert x == pytest.approx(x_expected, rel=1e-12, abs=1e-15)
