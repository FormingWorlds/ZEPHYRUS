"""Tests for ``src/zephyrus/boiloff.py``.

``boiloff.py`` is a physics source: it decides whether an envelope is bound
tightly enough for the energy-limited formula to apply, and how much of it
survives if not. The invariants exercised here are

  - Boundedness: the surviving mass fraction lies in ``(0, 1]`` for every
    admissible geometry, and is clamped rather than allowed to reach zero.
  - Monotonicity: the fraction decreases as the envelope inflates relative to
    its Bondi radius; the Jeans parameter falls as the radius grows.
  - Pinned values with discrimination guards: the Bondi radius and the Jeans
    parameter are pinned against closed forms, with guards against the
    plausible mean-molecular-weight, exponent and unit errors.
  - Error contract: non-positive masses, radii and temperatures raise.

The threshold consistency between the two regime measures is pinned directly:
the mass factor must return unity at exactly the radius where boil-off is
expected to stop.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import pytest

from zephyrus.boiloff import (
    BOILOFF_RADIUS_FRACTION,
    LAMBDA_CRITICAL,
    MU_ENVELOPE,
    boiloff_mass_factor,
    bondi_radius,
    restricted_jeans,
)
from zephyrus.constants import G, kb, m_H

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# A sub-Neptune with an inflated envelope: 3 Earth masses at 640 K. These
# values sit deep in the boil-off regime, where the wrong branch of the
# physics changes the answer by orders of magnitude rather than percent.
MP_SUBNEPTUNE = 1.80e25
TEQ_SUBNEPTUNE = 638.6


@pytest.mark.physics_invariant
def test_bondi_radius_matches_closed_form_and_scales_with_mass():
    """Bondi radius reproduces its closed form and grows linearly with mass."""
    value = bondi_radius(MP_SUBNEPTUNE, TEQ_SUBNEPTUNE)
    expected = G * MP_SUBNEPTUNE * MU_ENVELOPE * m_H / (2.0 * kb * TEQ_SUBNEPTUNE)
    assert value == pytest.approx(expected, rel=1e-12)

    # Mean-molecular-weight guard: computing this for atomic hydrogen instead
    # of a hydrogen-helium mixture shrinks the radius by the factor 2.35, far
    # outside any sensible tolerance, so a silent convention swap fails here.
    wrong_mu = bondi_radius(MP_SUBNEPTUNE, TEQ_SUBNEPTUNE, mu=1.0)
    assert value / wrong_mu == pytest.approx(MU_ENVELOPE, rel=1e-12)

    # Sign and scale guards: the radius is positive and of order 1e8 m for a
    # sub-Neptune, not 1e5 m (a gram-vs-kilogram slip) or 1e11 m.
    assert value > 0.0
    assert 1e7 < value < 1e9

    # Linear in mass, inverse in temperature: a squared or square-rooted mass
    # dependence would break this ratio.
    assert bondi_radius(2.0 * MP_SUBNEPTUNE, TEQ_SUBNEPTUNE) / value == pytest.approx(2.0)
    assert bondi_radius(MP_SUBNEPTUNE, 2.0 * TEQ_SUBNEPTUNE) / value == pytest.approx(0.5)


@pytest.mark.physics_invariant
def test_restricted_jeans_matches_closed_form_and_falls_with_radius():
    """Jeans parameter reproduces its closed form and falls as the planet inflates."""
    radius = 1.0e8
    value = restricted_jeans(MP_SUBNEPTUNE, radius, TEQ_SUBNEPTUNE)
    expected = G * MP_SUBNEPTUNE * m_H / (kb * TEQ_SUBNEPTUNE * radius)
    assert value == pytest.approx(expected, rel=1e-12)

    # Scale guard: an inflated sub-Neptune sits well below the critical value,
    # which is the whole point of the diagnostic.
    assert 0.0 < value < LAMBDA_CRITICAL

    # Inverse in radius: a regression to a squared radius would halve this
    # ratio instead of leaving it at 2.
    doubled = restricted_jeans(MP_SUBNEPTUNE, 2.0 * radius, TEQ_SUBNEPTUNE)
    assert value / doubled == pytest.approx(2.0, rel=1e-12)

    # A compact planet of the same mass is firmly outside the regime, so the
    # diagnostic separates the two cases rather than returning a constant.
    compact = restricted_jeans(MP_SUBNEPTUNE, 6.4e6, TEQ_SUBNEPTUNE)
    assert compact > LAMBDA_CRITICAL


@pytest.mark.physics_invariant
def test_mass_factor_is_unity_exactly_at_the_boil_off_radius():
    """Surviving fraction reaches unity where boil-off is expected to cease."""
    r_bondi = bondi_radius(MP_SUBNEPTUNE, TEQ_SUBNEPTUNE)
    at_threshold = BOILOFF_RADIUS_FRACTION * r_bondi
    assert boiloff_mass_factor(at_threshold, r_bondi) == pytest.approx(1.0, rel=1e-12)

    # Exponent/asymptote guard: the widely quoted small-alpha form 0.2/alpha
    # gives 1.33 at alpha = 0.15 where the exact expression gives 0.8. A
    # regression to that limit would fail by a wide margin.
    alpha = 0.15
    value = boiloff_mass_factor(alpha * r_bondi, r_bondi)
    assert value == pytest.approx(0.8, rel=1e-12)
    assert abs(value - 0.2 / alpha) > 0.5

    # Inside the threshold the envelope is bound and nothing is removed.
    assert boiloff_mass_factor(0.01 * r_bondi, r_bondi) == pytest.approx(1.0)


@pytest.mark.physics_invariant
def test_mass_factor_is_bounded_and_decreases_as_the_envelope_inflates():
    """Surviving fraction stays in (0, 1] and shrinks monotonically with radius."""
    r_bondi = bondi_radius(MP_SUBNEPTUNE, TEQ_SUBNEPTUNE)
    alphas = (0.05, 0.1, 0.3, 1.0, 1.59, 5.0)
    factors = [boiloff_mass_factor(a * r_bondi, r_bondi) for a in alphas]

    for factor in factors:
        assert 0.0 < factor <= 1.0
    for earlier, later in zip(factors, factors[1:]):
        assert later <= earlier

    # An envelope beyond its Bondi radius, as reached by the inflated cases
    # this diagnostic exists for, keeps only a small fraction of its mass:
    # 2 / (10 * 1.59 + 1) = 2 / 16.9.
    assert boiloff_mass_factor(1.59 * r_bondi, r_bondi) == pytest.approx(0.118343, rel=1e-5)

    # The clamp keeps a grossly unbound envelope from being erased entirely,
    # so downstream mass ratios stay finite.
    assert boiloff_mass_factor(1e6 * r_bondi, r_bondi) == pytest.approx(1e-3)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_published_threshold_equivalence_holds_for_atomic_hydrogen():
    """Published equality of the two boil-off thresholds, and its convention.

    Fossati et al. (2017) and Affolter et al. (2023) state that a restricted
    Jeans parameter of 20 is equivalent to the Owen & Wu (2016) condition that
    the planet radius reaches a tenth of the Bondi radius. The identity holds
    when both are evaluated for atomic hydrogen, since
    ``Rp / R_B = 2 / (mu * Lambda)``; it does not survive a swap to the
    hydrogen-helium mean molecular weight, which this test pins explicitly so
    the convention cannot drift unnoticed.
    """
    radius = 1.0e8
    # Choose the mass that places this radius exactly at Lambda = 20.
    mass = LAMBDA_CRITICAL * kb * TEQ_SUBNEPTUNE * radius / (G * m_H)
    assert restricted_jeans(mass, radius, TEQ_SUBNEPTUNE) == pytest.approx(
        LAMBDA_CRITICAL, rel=1e-12
    )

    atomic_bondi = bondi_radius(mass, TEQ_SUBNEPTUNE, mu=1.0)
    assert radius / atomic_bondi == pytest.approx(BOILOFF_RADIUS_FRACTION, rel=1e-12)

    # Convention guard: the same radius sits at 0.0426 of the hydrogen-helium
    # Bondi radius, so the published equality is specific to atomic hydrogen.
    envelope_bondi = bondi_radius(mass, TEQ_SUBNEPTUNE)
    assert radius / envelope_bondi == pytest.approx(
        BOILOFF_RADIUS_FRACTION / MU_ENVELOPE, rel=1e-12
    )
    assert abs(radius / envelope_bondi - BOILOFF_RADIUS_FRACTION) > 0.05


def test_non_positive_inputs_raise_for_every_entry_point():
    """Each function rejects unphysical geometry rather than returning a number."""
    with pytest.raises(ValueError, match='positive mass'):
        bondi_radius(-1.0, TEQ_SUBNEPTUNE)
    with pytest.raises(ValueError, match='positive mass'):
        bondi_radius(MP_SUBNEPTUNE, 0.0)
    with pytest.raises(ValueError, match='positive mass'):
        restricted_jeans(MP_SUBNEPTUNE, 0.0, TEQ_SUBNEPTUNE)
    with pytest.raises(ValueError, match='positive radii'):
        boiloff_mass_factor(1.0e8, -1.0)
