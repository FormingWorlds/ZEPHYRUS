"""Tests for ``src/zephyrus/hydrodynamic.py``.

Exercises the energy-limited and radiation-recombination-limited rates and
the machinery that selects between them. The physical invariants under test:

- Reference pins: the Erkaev et al. (2007) Table 1 tidal enhancement
  factors within 1 percent; the Caldiroli et al. (2022) efficiency fit at
  spot points of its validity box; the Lopez (2017) printed wind mean-mass
  pairs; the Murray-Clay et al. (2009) fiducial hot Jupiter (base Jeans
  parameter, sonic radius, and sonic-point Knudsen numbers with their
  Coulomb cross section).
- Monotonicity / symmetry: the RR rate scales as the square root of the
  XUV flux and the EL rate linearly in it.
- Boundedness / error contract: the subcritical floor holds the density at
  the base value, and the efficiency fit rejects its complex-valued region
  by falling back with a flag.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.constants import G, m_p
from zephyrus.hydrodynamic import (
    caldiroli_efficiency,
    el_rate,
    hill_radius_periapsis,
    k_tide,
    rr_chain,
    selection_mechanism,
    wind_mean_masses,
)
from zephyrus.knudsen import mean_free_path, sonic_scale_height
from zephyrus.planets_parameters import Me, Mjup, Ms, Re, Rjup

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Erkaev et al. (2007) Table 1: (xi, printed enhancement factor 1/K).
ERKAEV_TABLE1 = [
    (2.9, 1.97),
    (3.0, 1.92),
    (3.5, 1.70),
    (3.7, 1.65),
    (3.8, 1.63),
    (4.3, 1.53),
    (5.9, 1.34),
]

# Spot evaluations of the Caldiroli et al. (2022) Appendix A.1 fit across
# its validity box at K = 1: (log10 phi [cgs], F_XUV/rho_p [cgs], eta).
# The eta values were evaluated from the published fitting formulas at
# transcription time, independently of this implementation, so a later
# transcription error in either place breaks the agreement.
CALDIROLI_SPOTS = [
    (12.20, 1e3, 8.8e-1),
    (12.20, 1e6, 1.5e-1),
    (12.80, 1e4, 4.5e-1),
    (13.00, 1e4, 1.2e-1),
    (13.10, 1e5, 2.6e-2),
    (13.29, 1e6, 8.7e-4),
]


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_erkaev_table1_enhancement_factors():
    """The tidal factor reproduces the published Table 1 within 1 percent.

    Erkaev et al. (2007) print the enhancement ``1/K`` for seven observed
    hot Jupiters; each must match, and the factor must fall monotonically
    toward 1 as ``xi`` grows (boundedness: ``K`` in (0, 1) throughout).
    """
    for xi, inv_k in ERKAEV_TABLE1:
        assert 1.0 / k_tide(xi) == pytest.approx(inv_k, rel=0.011), xi
        assert 0.0 < k_tide(xi) < 1.0
    factors = [1.0 / k_tide(xi) for xi, _ in ERKAEV_TABLE1]
    assert all(a >= b for a, b in zip(factors, factors[1:]))


def _planet_for(log_phi, f_over_rho_cgs):
    """Invert (log10 phi, F/rho in cgs) into (M_p, R_p, F_xuv) at 2 Earth radii."""
    r_p = 2.0 * Re
    phi_si = 10**log_phi * 1e-4
    m_planet = phi_si * r_p / G
    rho_p = m_planet / (4.0 / 3.0 * math.pi * r_p**3)
    f_cgs = f_over_rho_cgs * (rho_p * 1e-3)
    return m_planet, r_p, f_cgs * 1e-3


@pytest.mark.reference_pinned
def test_caldiroli_fit_spot_values_and_flux_guard():
    """The efficiency fit reproduces spot evaluations and rejects its floor.

    Six points across the fit's validity box reproduce the published
    formulas to 5 percent, spanning three decades of efficiency (the
    collapse toward 1e-3 at high potential is the physical content). Below
    ``F_XUV / rho_p = 1e2`` cgs the formulas turn complex; the fit must
    return None with the flag instead of a complex or extrapolated number.
    """
    for log_phi, f_over_rho, eta_ref in CALDIROLI_SPOTS:
        m_planet, r_p, f_xuv = _planet_for(log_phi, f_over_rho)
        eta, _flags = caldiroli_efficiency(f_xuv, m_planet, r_p, K=1.0)
        assert eta == pytest.approx(eta_ref, rel=0.05), (log_phi, f_over_rho)
    m_planet, r_p, f_xuv = _planet_for(12.5, 0.5)  # F/rho below the 1e2 bound
    eta, flags = caldiroli_efficiency(f_xuv, m_planet, r_p, K=1.0)
    assert eta is None
    assert flags.get('caldiroli_below_flux_bound') is True


@pytest.mark.reference_pinned
def test_wind_mean_masses_reproduce_lopez_pairs():
    """The ionized-wind mean-mass rule reproduces the published pairs.

    Lopez (2017) prints (mu_wind, mu_plus) = (0.62, 1.3) proton masses for
    a 90/10 H/He wind and (3, 6) for steam (fully dissociated 2:1 H:O).
    The generalized rule (electrons counted, heavies singly ionized) must
    recover both, and the per-ion mass must always be twice the per-particle
    mass by construction.
    """
    mu_w, mu_i = wind_mean_masses({'H': 0.9, 'He': 0.1})
    assert mu_w == pytest.approx(0.62, rel=0.06)
    assert mu_i == pytest.approx(1.3, rel=0.02)
    mu_w, mu_i = wind_mean_masses({'H': 2.0 / 3.0, 'O': 1.0 / 3.0})
    assert mu_w == pytest.approx(3.0, rel=0.01)
    assert mu_i == pytest.approx(6.0, rel=0.01)
    assert mu_i == pytest.approx(2.0 * mu_w, rel=1e-12)


@pytest.mark.physics_invariant
def test_rr_subcritical_floor_semantics():
    """A subcritical sonic point floors at the base with the base density.

    A cool 1e4 K hydrogen wind on an Earth-mass planet has its formal sonic
    radius inside the base radius; the chain must flag it, floor the sonic
    radius at the base, hold the density at the base value, and apply no
    barometric suppression (factor exactly 1).
    """
    rr = rr_chain(1.0 * Me, 10.0, 1.0 * Re, 1.0e4, {'H': 1.0})
    assert rr['subcritical'] is True
    assert rr['R_s'] == pytest.approx(1.0 * Re, rel=1e-12)
    assert rr['R_s_calc'] < 1.0 * Re
    assert rr['rho_s'] == pytest.approx(rr['rho_base'], rel=1e-12)
    assert rr['barometric_factor'] == pytest.approx(1.0, rel=1e-12)
    assert selection_mechanism(rr, el_won=False) == 'RR-selected:subcritical-floor'


@pytest.mark.physics_invariant
def test_rr_barometric_factor_and_mechanism_labels():
    """Supercritical winds carry exp(3/2 - lambda_b), labeled by mechanism.

    A heavy C-O wind on a massive planet is strongly bound: the barometric
    factor is exactly ``exp(3/2 - lambda_b)``, and with ``lambda_b`` above
    the reporting split an RR win is labeled barometric suppression, never
    recombination saturation (the category-error guard). The EL-win label
    is independent of the chain state.
    """
    rr = rr_chain(10.0 * Me, 10.0, 2.0 * Re, 1.0e4, {'C': 1.0 / 3.0, 'O': 2.0 / 3.0})
    assert rr['subcritical'] is False
    assert rr['barometric_factor'] == pytest.approx(math.exp(1.5 - rr['lambda_b']), rel=1e-12)
    assert rr['lambda_b'] > 4.0
    assert selection_mechanism(rr, el_won=False) == 'RR-selected:barometric-suppression'
    assert selection_mechanism(rr, el_won=True) == 'EL-selected'
    # A loosely bound hydrogen wind (base Jeans parameter between the
    # supercritical floor at 2 and the reporting split at 4) lands in
    # genuine recombination saturation.
    rr_h = rr_chain(0.7 * Mjup, 5.0, 2.0 * Rjup, 1.0e4, {'H': 1.0})
    assert not rr_h['subcritical']
    assert 2.0 < rr_h['lambda_b'] < 4.0
    assert selection_mechanism(rr_h, el_won=False) == 'RR-selected:recombination-saturation'


@pytest.mark.physics_invariant
def test_flux_scalings_of_both_limits():
    """The RR rate scales as sqrt(F_XUV) and the EL rate linearly in it.

    The square root comes from ionization equilibrium at the base (ion
    density proportional to sqrt of the ionizing flux); linearity is the
    energy-limited budget. A spurious offset or a wrong power fails the
    exact two-point ratios.
    """
    a = rr_chain(5 * Me, 1.0, 1.5 * Re, 1e4, {'H': 1.0})
    b = rr_chain(5 * Me, 100.0, 1.5 * Re, 1e4, {'H': 1.0})
    assert b['mdot_rr'] / a['mdot_rr'] == pytest.approx(10.0, rel=1e-6)
    lo = el_rate(0.1, 1.0, Re, 1.1 * Re, Me, 1.0)
    hi = el_rate(0.1, 100.0, Re, 1.1 * Re, Me, 1.0)
    assert hi / lo == pytest.approx(100.0, rel=1e-12)
    # Zero-flux limits: no driver, no escape, for both chains.
    assert el_rate(0.1, 0.0, Re, 1.1 * Re, Me, 1.0) == pytest.approx(0.0, abs=1e-30)
    assert rr_chain(5 * Me, 0.0, 1.5 * Re, 1e4, {'H': 1.0})['mdot_rr'] == pytest.approx(
        0.0, abs=1e-30
    )


@pytest.mark.reference_pinned
def test_murray_clay_fiducial_hot_jupiter_anchors():
    """The chain reproduces the Murray-Clay et al. (2009) fiducial planet.

    Their planet (0.7 Jupiter masses, 1.4 Jupiter radii, ionized hydrogen
    wind at 1e4 K) has a base Jeans parameter of 5.49 (5 percent tolerance:
    their radius convention rounds R_Jup), a sonic point at
    ``lambda_b / 2`` planetary radii inside their stated 2 to 4, and
    sonic-point Knudsen numbers, evaluated with their proton-proton Coulomb
    cross section ``1e-13 (T / 1e4 K)^-2 cm^2``, near the printed 1e-4 at
    450 erg cm^-2 s^-1 and 1e-5 at 5e5 erg cm^-2 s^-1 (factor-3 tolerance:
    their scale-height convention differs at order unity).
    """
    m_planet, r_p = 0.7 * Mjup, 1.4 * Rjup
    rr = rr_chain(m_planet, 0.45, r_p, 1.0e4, {'H': 1.0})
    assert rr['lambda_b'] == pytest.approx(5.49, rel=0.05)
    assert 2.0 < rr['R_s'] / r_p < 4.0
    assert rr['R_s'] / r_p == pytest.approx(rr['lambda_b'] / 2.0, rel=1e-9)
    sigma_pp = 1e-13 * 1e-4  # cm^2 -> m^2 at 1e4 K
    for f_si, kn_ref in ((0.45, 1e-4), (500.0, 1e-5)):
        rr = rr_chain(m_planet, f_si, r_p, 1.0e4, {'H': 1.0})
        n_sc = rr['rho_s'] / (rr['mu_plus_wind'] * m_p)
        kn = mean_free_path(sigma_pp, n_sc) / sonic_scale_height(rr['R_s'], 1.0)
        assert kn_ref / 3.0 < kn < kn_ref * 3.0, f_si


def test_hill_radius_periapsis_geometry():
    """The periapsis Hill radius scales with a (1 - e) and the mass ratio.

    Eccentricity shrinks the periapsis linearly; the mass dependence is the
    cube root. The circular Earth-Sun value lands near 0.01 au (the
    well-known Hill-sphere scale), the sanity anchor.
    """
    r0 = hill_radius_periapsis(Me, Ms, 1.496e11, 0.0)
    r3 = hill_radius_periapsis(Me, Ms, 1.496e11, 0.3)
    assert r3 == pytest.approx(0.7 * r0, rel=1e-12)
    r8m = hill_radius_periapsis(8 * Me, Ms, 1.496e11, 0.0)
    assert r8m == pytest.approx(2.0 * r0, rel=1e-12)
    # Earth's Hill radius is about 1.5e9 m (0.01 au).
    assert 1.3e9 < r0 < 1.6e9
