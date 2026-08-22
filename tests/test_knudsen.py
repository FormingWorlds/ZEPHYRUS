"""Tests for ``src/zephyrus/knudsen.py``.

Exercises the neutral collision cross-section ladder and the sonic-point
Knudsen switch. The physical invariants under test:

- Reference pin: the Laricchiuta et al. (2009) collision integrals reproduce
  the measured room-temperature viscosities of N2, O2, CO, and CO2 within
  7 percent through the first Chapman-Enskog approximation.
- Closed form: the sonic-point Knudsen number equals the Maxwell mean free
  path over the analytic sonic-point scale height; the mixture cross section
  is exactly density weighted.
- Monotonicity / boundedness: cross sections shrink with temperature; the
  temperature-independent geometric rung overshoots the collision-integral
  rung at high temperature (its documented bias); the hysteresis window
  widens or tightens the switch threshold on the correct side.
- Error contract: the hydrogen route rejects non-hydrogen species.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.knudsen import (
    KN_BAND,
    effective_threshold,
    kn_sonic,
    lar_sigma_diff,
    mean_free_path,
    sigma_geometric,
    sigma_mixture,
    sigma_species,
    sigma_zk90_hydrogen,
    sonic_scale_height,
    viscosity_pure,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Momentum-transfer cross sections [m^2] evaluated from the transcribed
# Laricchiuta et al. (2009) fit at transcription time and validated then
# against the measured viscosities below. They serve as transcription
# regression pins: a typo in any Table 3/4/5 coefficient moves them.
PINNED_SIGMA = {
    ('N', 'N'): {
        300: 2.56e-19,
        1000: 2.00e-19,
        3000: 1.58e-19,
        10000: 1.16e-19,
        20000: 9.61e-20,
    },
    ('O', 'O'): {300: 2.29e-19, 10000: 1.06e-19},
    ('C', 'C'): {300: 3.11e-19, 10000: 1.41e-19},
    ('C', 'N'): {300: 2.81e-19, 10000: 1.27e-19},
    ('C', 'O'): {300: 2.64e-19, 10000: 1.21e-19},
    ('N', 'O'): {300: 2.41e-19, 10000: 1.10e-19},
    ('N2', 'N2'): {300: 3.66e-19, 10000: 1.73e-19},
    ('O2', 'O2'): {300: 3.63e-19, 10000: 1.72e-19},
    ('CO', 'CO'): {300: 3.87e-19, 10000: 1.81e-19},
    ('CO2', 'CO2'): {300: 5.34e-19, 10000: 2.22e-19},
    ('N', 'CO2'): {300: 3.56e-19, 10000: 1.60e-19},
    ('O', 'CO2'): {300: 3.34e-19, 10000: 1.55e-19},
    ('N2', 'CO2'): {300: 4.38e-19, 10000: 1.97e-19},
}

# Measured dynamic viscosities at 300 K [Pa s] with the molar masses [g/mol]:
# standard handbook values (e.g. the CRC Handbook of Chemistry and Physics).
MEASURED_VISCOSITY = {
    ('N2', 'N2'): (28.014, 17.9e-6),
    ('O2', 'O2'): (31.998, 20.7e-6),
    ('CO', 'CO'): (28.010, 17.8e-6),
    ('CO2', 'CO2'): (44.009, 15.0e-6),
}


def test_laricchiuta_cross_sections_match_transcription_pins():
    """Every tabulated pair reproduces its pinned cross section at every T.

    The pins were evaluated from the published fit when the coefficient
    tables were transcribed, so any later coefficient corruption fails
    here. The temperature trend is the physical guard: collision integrals
    of these attractive-well pairs shrink monotonically with temperature.
    """
    for pair, vals in PINNED_SIGMA.items():
        for T, ref in vals.items():
            assert lar_sigma_diff(pair, T) == pytest.approx(ref, rel=0.01), (pair, T)
    # Monotone decrease with temperature for a representative pair.
    ts = [300.0, 1000.0, 3000.0, 10000.0, 20000.0]
    sigmas = [lar_sigma_diff(('N', 'N'), t) for t in ts]
    assert all(a > b for a, b in zip(sigmas, sigmas[1:]))


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_viscosities_reproduce_measurements_within_7_percent():
    """The collision integrals reproduce four measured viscosities to 7 percent.

    The first Chapman-Enskog approximation converts the Omega^(2,2)*
    integral into the pure-gas dynamic viscosity; the measured 300 K values
    of N2, O2, CO, and CO2 (CRC Handbook) anchor the whole transcription on
    laboratory data. The 7 percent bound is the observed worst case. The
    scale guard: all four are tens of micropascal seconds, so a unit slip
    in the conversion fails by orders of magnitude.
    """
    for pair, (mass, eta_meas) in MEASURED_VISCOSITY.items():
        eta = viscosity_pure(pair, mass, 300.0)
        assert abs(eta / eta_meas - 1.0) < 0.07, pair
        # Scale guard: micropascal-second range.
        assert 1e-5 < eta < 3e-5


def test_hydrogen_route_pins_and_temperature_dependence():
    """The hydrogen cross sections match their published-route values.

    The Zahnle et al. (1990) Eq. (30) inversion of the Zahnle & Kasting
    (1986) Table I diffusion parameter gives sigma(H-H) = 6.4e-20 m^2 at
    1e4 K and sigma(H2-H2) = 2.07e-19 m^2 at 300 K, with a T^-0.25
    dependence by construction. The error contract: any species other than
    H or H2 is rejected.
    """
    assert sigma_zk90_hydrogen('H', 1e4) == pytest.approx(6.4e-20, rel=0.02)
    assert sigma_zk90_hydrogen('H2', 300.0) == pytest.approx(2.07e-19, rel=0.02)
    r = sigma_zk90_hydrogen('H', 1e4) / sigma_zk90_hydrogen('H', 1e2)
    assert r == pytest.approx((1e4 / 1e2) ** -0.25, rel=1e-9)
    with pytest.raises(ValueError, match='hydrogen route'):
        sigma_zk90_hydrogen('He', 1e4)


def test_ladder_provenance_and_geometric_bias():
    """The ladder assigns the right rung and the hard-sphere bias shows.

    N sits on the collision-integral rung, H on the hydrogen route, He (no
    Laricchiuta entry, not hydrogen) on the geometric rung. A hard sphere
    has no temperature dependence, so at 300 K the geometric N cross
    section is comparable to the collision integral, while at 1e4 K it
    overshoots by a factor of a few, the documented bias that pushes the
    Knudsen number low, toward hydrodynamic verdicts.
    """
    _, prov_n = sigma_species('N', 1e4)
    _, prov_h = sigma_species('H', 1e4)
    _, prov_he = sigma_species('He', 1e4)
    assert prov_n == 'laricchiuta'
    assert prov_h == 'zk90-scaled'
    assert prov_he == 'geometric-vdw'
    geo = sigma_geometric('N')
    assert geo / lar_sigma_diff(('N', 'N'), 300.0) == pytest.approx(1.0, abs=0.25)
    assert 2.0 < geo / lar_sigma_diff(('N', 'N'), 1e4) < 4.0  # high-T overshoot
    # A composite molecule with no tabulated radius falls back to its
    # largest constituent element without raising.
    geo_h2o = sigma_geometric('H2O')
    assert geo_h2o == pytest.approx(math.pi * (2.0 * 1.52e-10) ** 2, rel=1e-12)


@pytest.mark.physics_invariant
def test_kn_sonic_equals_mfp_over_scale_height():
    """The switch value is exactly the mean free path over the scale height.

    Closed-form identity of the construction: Kn_sc from the packaged
    formula must equal ``mean_free_path(sigma, n) / sonic_scale_height(r)``
    with the same mixture cross section. The edge case is a rarefied state
    (tiny n) where Kn grows without bound but stays finite and positive.
    """
    vmr = {'N': 1.0}
    n, r, T = 1e14, 1e7, 8000.0
    kn, sigma, prov = kn_sonic(n, r, vmr, T, gamma=1.0)
    assert kn == pytest.approx(mean_free_path(sigma, n) / sonic_scale_height(r, 1.0), rel=1e-12)
    assert prov == {'N': 'laricchiuta'}
    # Rarefied edge: eight decades less dense means eight decades larger Kn.
    kn_thin, _, _ = kn_sonic(n * 1e-8, r, vmr, T)
    assert kn_thin == pytest.approx(kn * 1e8, rel=1e-9)
    assert math.isfinite(kn_thin)


@pytest.mark.physics_invariant
def test_mixture_cross_section_is_density_weighted():
    """The mixture rule interpolates linearly between the pure endpoints.

    Density weighting means a 25/75 N-CO2 atom mix must give exactly
    ``0.25 sigma_N + 0.75 sigma_CO2``, and any mixture must land strictly
    between the pure endpoints. Renormalization: scaling all mole fractions
    by a constant leaves the result unchanged.
    """
    T = 1e4
    s_n, _ = sigma_species('N', T)
    s_co2, _ = sigma_species('CO2', T)
    mix, _ = sigma_mixture({'N': 0.25, 'CO2': 0.75}, T)
    assert mix == pytest.approx(0.25 * s_n + 0.75 * s_co2, rel=1e-12)
    assert min(s_n, s_co2) < mix < max(s_n, s_co2)
    mix_scaled, _ = sigma_mixture({'N': 2.5, 'CO2': 7.5}, T)
    assert mix_scaled == pytest.approx(mix, rel=1e-12)


def test_hysteresis_window_moves_the_threshold_the_right_way():
    """The hysteresis threshold resists leaving the previous regime.

    Sharp threshold with no memory; leaving a hydrodynamic state the
    threshold rises (harder to switch to hydrostatic); leaving a
    hydrostatic state it falls (harder to switch back). An unrecognized
    previous label falls through to the sharp value. The window must
    bracket the sharp threshold from both sides, and the diagnostic band
    edges stay ordered.
    """
    assert effective_threshold(1.0, 1.5, None) == pytest.approx(1.0)
    up = effective_threshold(1.0, 1.5, 'hydrodynamic:EL')
    dn = effective_threshold(1.0, 1.5, 'hydrostatic')
    assert up == pytest.approx(1.5)
    assert dn == pytest.approx(1.0 / 1.5)
    assert dn < 1.0 < up
    assert effective_threshold(1.0, 1.5, 'boiloff') == pytest.approx(1.0)
    # The printed diagnostic band brackets the default threshold of 1.
    assert KN_BAND[0] < 1.0 <= KN_BAND[1]
