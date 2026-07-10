"""Tests for ``src/zephyrus/escape.py``.

Exercises the energy-limited (EL) atmospheric-escape mass-loss rate and its
tidal correction. The physical invariants under test:

- Conservation / closed form: the EL rate equals
  ``epsilon * pi * R^3 * Fxuv / (G * Mp * K_tide)`` for the selected radius
  scaling, pinned against the Erkaev et al. (2007) formulation
  (``scaling=2``) and the Lehmer & Catling (2017) variant (``scaling=3``).
- Positivity / boundedness: the rate is non-negative for valid inputs and
  the tidal factor ``K_tide`` lies in ``(0, 1)`` when the Hill radius exceeds
  the XUV radius (``ksi = Rhill / Rxuv > 1``), rising toward 1 as the orbit
  widens; the close-in geometries under test keep it well below 1.
- Monotonicity / symmetry: the rate is linear in ``Fxuv``, scales as
  ``1 / Mp``, is larger with the tidal correction than without, and is zero
  when ``Fxuv = 0``.
- Error contract: an unsupported ``scaling`` raises ``ValueError``, and the
  tidal branch raises ``ValueError`` for ``ksi <= 1``, where the atmosphere
  reaches the Roche lobe and the energy-limited approximation no longer holds.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import numpy as np
import pytest

from zephyrus.constants import G, au2m
from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Ms, Re

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Reference geometry shared across the closed-form pins. Rxuv is chosen
# distinct from Rp (Rxuv = 1.2 * Rp) so the scaling=2 (Rp * Rxuv**2) and
# scaling=3 (Rxuv**3) branches give values 20% apart; a regression that
# swaps the default scaling then fails loudly instead of silently.
EPSILON = 0.15
RP = Re
RXUV = 1.2 * Re
FXUV = 14.67  # W m-2, Earth at 10 Myr (Wordsworth et al. 2018, Fig 9)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_el_escape_scaling2_matches_erkaev2007_closed_form():
    """Pin the EL rate for the scaling=2 (Erkaev et al. 2007) branch.

    The default radius scaling uses ``R^3 = Rp * Rxuv**2``. The pinned value
    is hand-evaluated from ``epsilon * pi * Rp * Rxuv**2 * Fxuv / (G * Mp)``
    at the Earth-like reference geometry with no tidal correction.
    """
    val = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    expected = 6479585.361332079  # kg s-1, hand-evaluated at the reference geometry
    # Same closed form and constants as the source, so exact to float rounding.
    assert val == pytest.approx(expected, rel=1e-9)
    # Wrong-scaling discrimination: scaling=3 uses Rxuv**3, not Rp * Rxuv**2.
    # With Rxuv = 1.2 * Rp the two branches are 20% apart, far above tolerance.
    wrong_scaling3 = EPSILON * np.pi * RXUV**3 * FXUV / (G * Me)
    assert abs(val - wrong_scaling3) > 0.01 * expected
    # Sign guard: an energy-limited rate is always non-negative.
    assert val > 0
    # Scale guard: the rate is order 1e6-1e7 kg s-1. A gross unit slip, such
    # as feeding radii in cm instead of m, lands near 1e12 and trips this
    # band; a dropped pi or epsilon factor stays within the band and is caught
    # by the rel=1e-9 pin above, not here.
    assert 1e6 < val < 1e8


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_el_escape_scaling3_matches_lehmer_catling_closed_form():
    """Pin the EL rate for the scaling=3 (Lehmer & Catling 2017) branch.

    The alternative radius scaling uses ``R^3 = Rxuv**3``. The pinned value
    is hand-evaluated from ``epsilon * pi * Rxuv**3 * Fxuv / (G * Mp)`` at the
    same Earth-like reference geometry.
    """
    val = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=3)
    expected = 7775502.433598499  # kg s-1, hand-evaluated at the reference geometry
    assert val == pytest.approx(expected, rel=1e-9)
    # Wrong-scaling discrimination: scaling=2 uses Rp * Rxuv**2 and is 20%
    # smaller here, so a regression that reverts the exponent fails.
    wrong_scaling2 = EPSILON * np.pi * RP * RXUV**2 * FXUV / (G * Me)
    assert abs(val - wrong_scaling2) > 0.01 * expected
    # Sign guard and scale guard as above.
    assert val > 0
    assert 1e6 < val < 1e8


@pytest.mark.parametrize(
    'scaling', [0, 1, 4, -1], ids=['scaling_0', 'scaling_1', 'scaling_4', 'scaling_neg1']
)
def test_el_escape_invalid_scaling_raises(scaling):
    """An unsupported radius-scaling exponent raises ``ValueError``.

    Only ``scaling`` values 2 and 3 are defined; every other integer must
    reach the guard and raise, and no mass-loss rate is returned.
    """
    with pytest.raises(ValueError, match='Invalid radius exponent'):
        EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=scaling)
    # Contrast: a supported scaling returns a finite positive rate at the
    # same geometry, so the raise is specific to the bad exponent.
    ok = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    assert np.isfinite(ok)
    assert ok > 0


def test_el_escape_tidal_raises_below_roche_lobe():
    """The tidal branch rejects sub-Roche-lobe geometries (``ksi <= 1``).

    The Erkaev et al. (2007) factor ``K_tide = (ksi - 1)**2 (2 ksi + 1) /
    (2 ksi**3)`` has a double root at ``ksi = 1``, so the energy-limited rate
    divides by ``K_tide`` and is defined only for ``ksi > 1``. A close-in,
    highly eccentric orbit (``a = 0.02 au``, ``e = 0.90``) shrinks the
    periapsis Hill radius until ``ksi`` is about 0.39, inside the Roche lobe,
    and the call must raise rather than return a suppressed or divergent rate.
    """
    a = 0.02 * au2m
    e = 0.90  # high eccentricity pulls the periapsis Hill radius inside Rxuv
    rhill = a * (1 - e) * (Me / (3 * Ms)) ** (1 / 3)
    ksi = rhill / RXUV
    # Confirm the constructed geometry is genuinely sub-Roche-lobe.
    assert ksi == pytest.approx(0.3910754106364523, rel=1e-9)
    assert ksi < 1.0
    with pytest.raises(ValueError, match='Roche lobe'):
        EL_escape(True, a, e, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    # Boundary case: setting Rxuv equal to the periapsis Hill radius puts
    # ``ksi`` exactly at the singular point 1, which must raise rather than
    # divide by zero.
    rxuv_at_hill = a * (Me / (3 * Ms)) ** (1 / 3)
    with pytest.raises(ValueError, match='Roche lobe'):
        EL_escape(True, a, 0.0, Me, Ms, EPSILON, RP, rxuv_at_hill, FXUV, scaling=2)
    # Contrast: the same close-in orbit with ``ksi > 1`` (circular, Rxuv well
    # inside the Hill radius) returns a finite positive rate, so the raise is
    # specific to ``ksi <= 1``, not to the tidal branch as a whole.
    ok = EL_escape(True, a, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    assert np.isfinite(ok)
    assert ok > 0.0


@pytest.mark.physics_invariant
def test_el_escape_linear_in_xuv_flux():
    """Doubling the XUV flux doubles the escape rate; zero flux gives zero.

    Linearity in ``Fxuv`` is exact for the closed form. The zero-flux limit
    is the edge case: no driver means no escape.
    """
    base = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    doubled = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, 2 * FXUV, scaling=2)
    # Exact proportionality: a spurious additive offset would break this.
    assert doubled == pytest.approx(2 * base, rel=1e-12)
    assert base > 0
    # Zero-flux limit: the rate collapses to exactly zero.
    zero = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, 0.0, scaling=2)
    assert zero == pytest.approx(0.0, abs=1e-30)


@pytest.mark.physics_invariant
def test_el_escape_decreases_with_planet_mass():
    """Escape scales as ``1 / Mp``: a heavier planet loses mass more slowly.

    The gravitational binding term ``G * Mp`` sits in the denominator, so
    doubling the mass halves the rate. The edge case is a very massive
    planet, whose rate stays strictly positive but small.
    """
    light = EL_escape(False, 1.0, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    heavy = EL_escape(False, 1.0, 0.0, 2 * Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    # Exact 1/Mp relation: doubling Mp halves the rate.
    assert heavy == pytest.approx(0.5 * light, rel=1e-12)
    # Monotone decrease, and the heavier planet still escapes at a positive rate.
    assert heavy < light
    massive = EL_escape(False, 1.0, 0.0, 1e3 * Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    assert 0 < massive < heavy


@pytest.mark.physics_invariant
def test_el_escape_tidal_correction_increases_escape():
    """The tidal correction raises the escape rate for a close-in orbit.

    At ``a = 0.02 au`` the Hill radius is only a few XUV radii, so
    ``K_tide`` departs from 1 by tens of percent. Because ``K_tide`` sits in
    the denominator, the tidal rate exceeds the no-tidal rate. The enhancement
    ``1 / K_tide`` is pinned, and the no-tidal value is the discrimination
    guard: a dropped ``K_tide`` would collapse the ratio to 1.
    """
    a = 0.02 * au2m  # close-in so K_tide is well below 1
    no_tidal = EL_escape(False, a, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    tidal = EL_escape(True, a, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    assert tidal > no_tidal
    # Pin the enhancement factor 1 / K_tide against the hand-evaluated K_tide.
    assert tidal / no_tidal == pytest.approx(1.6005072480546971, rel=1e-9)
    # Dropped-K_tide discrimination: the ratio is well above 1, not ~1.
    assert tidal / no_tidal - 1.0 > 0.1
    # Boundedness: for this close-in geometry the Hill radius stays above Rxuv
    # (ksi > 1), so the backed-out K_tide lies strictly in (0, 1); it reaches 1
    # only in the ksi -> infinity limit of an infinitely wide orbit.
    k_tide = no_tidal / tidal
    assert 0.0 < k_tide < 1.0


@pytest.mark.physics_invariant
def test_el_escape_tidal_eccentricity_increases_escape():
    """A more eccentric orbit escapes faster through the tidal correction.

    The tidal factor is evaluated at periapsis, where the Hill radius scales
    as ``a * (1 - e)``. Raising ``e`` shrinks the periapsis Hill radius, which
    lowers ``K_tide``; since ``K_tide`` sits in the denominator, the escape
    rate rises. The circular orbit (``e = 0``) is the boundary case. The wrong
    ``a * (1 + e)`` sign convention would widen the periapsis Hill radius and
    reverse the inequality, so the direction of the change is the guard.
    """
    a = 0.02 * au2m  # close-in so the tidal correction is sizeable
    circular = EL_escape(True, a, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    eccentric = EL_escape(True, a, 0.3, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    # Higher eccentricity gives a smaller periapsis Hill radius, so more escape.
    assert eccentric > circular
    # Pin the enhancement against the hand-evaluated K_tide(e=0) / K_tide(e=0.3).
    assert eccentric / circular == pytest.approx(1.3114173571726595, rel=1e-9)
    # Sign-convention discrimination: the (1 + e) periapsis slip would push the
    # ratio below 1, reversing the inequality asserted above.
    assert eccentric / circular > 1.0


@pytest.mark.physics_invariant
def test_el_escape_tidal_factor_approaches_unity_for_wide_orbit():
    """The tidal correction vanishes as the orbit widens (``K_tide -> 1``).

    This is the analytical limit: at large ``a`` the Hill radius dwarfs the
    XUV radius, so ``K_tide -> 1`` and the tidal rate converges to the
    no-tidal rate. A close-in orbit keeps ``K_tide`` well below 1, so the
    backed-out ``K_tide`` is strictly larger (closer to 1) for the wide orbit.
    """
    a_wide = 1.0 * au2m
    a_close = 0.02 * au2m
    no_tidal = EL_escape(False, a_wide, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    tidal_wide = EL_escape(True, a_wide, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    tidal_close = EL_escape(True, a_close, 0.0, Me, Ms, EPSILON, RP, RXUV, FXUV, scaling=2)
    k_wide = no_tidal / tidal_wide
    k_close = no_tidal / tidal_close
    # Wide orbit: K_tide within 1% of unity, so the tidal branch barely differs.
    assert k_wide == pytest.approx(1.0, abs=1e-2)
    # Close orbit: K_tide well below 1, so the correction is large.
    assert k_close < 0.9
    # Monotone in orbital distance: the wide orbit sits closer to unity.
    assert k_wide > k_close
