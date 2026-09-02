"""Tests for ``src/zephyrus/nozzle.py``.

Exercises the tidally driven L1 nozzle flow of Jackson et al. (2017,
ApJ 835, 145). The physical anchors under test:

- Reference pins: two planets of their Table 2 reproduced with their own
  input prescriptions; two lobe-filling binaries of their Table 1 landing
  on the Figure 5 solid curve; the equal-mass curvature A(1) = 8 their
  Section 2.1 prints.
- Cross-check: the Eq. (14) volume-averaged potential evaluated at the
  Eggleton lobe radius matches the exact corotating Roche potential at a
  numerically solved L1 point.
- Invariants: A(q) symmetric under mass-ratio inversion and bounded in
  (4, 8]; the lobe radius inside the Hill radius at planetary mass ratios;
  the rate invariant under the choice of launch level along an isothermal
  hydrostatic column; the exponent clamped at lobe contact, where the rate
  equals the lobe-filling boundary value.
- Error contract: a non-positive mass ratio raises.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.constants import G, amu, kb
from zephyrus.nozzle import (
    curvature_a,
    eggleton_lobe_radius,
    isothermal_column_density,
    nozzle_candidate,
    volume_averaged_potential,
)
from zephyrus.planets_parameters import Me, Ms, Re
from zephyrus.profiles import interp_at_pressure, isothermal_profile

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

AU = 1.495978707e11  # m
M_JUP = 1.8987e27  # kg, the value Jackson et al. (2017) state
R_JUP = 7.1492e7  # m, likewise
GYR = 3.15576e16  # s


def _jackson_photosphere(M_p, R_p, T_p, mu_kg):
    """Photospheric density by their Section 3 prescription.

    ``rho_ph = tau / (kappa sqrt(2 pi R_p H))`` with the slant optical
    depth tau = 0.56 (Howe & Burrows 2012) and kappa = 1e-2 cm^2/g.
    """
    g = G * M_p / R_p**2
    scale_height = kb * T_p / (mu_kg * g)
    kappa = 1e-3  # m^2/kg
    return 0.56 / (kappa * math.sqrt(2.0 * math.pi * R_p * scale_height))


@pytest.mark.reference_pinned
def test_table2_planets_reproduce_published_rates():
    """Kepler-21 b and CoRoT-24 b of Jackson et al. (2017) Table 2 within 3%.

    Both donors sit well inside their lobes, where the Eq. (14) potential
    approximation is good and the photospheric-radius convention is
    irrelevant, so the printed rates pin the whole formalism (Eqs. 3, 10,
    13, and 14 plus the Eggleton lobe radius). Inputs are their Table 2
    rows evaluated through their Section 3 prescriptions (mu = 1 amu above
    2000 K, 2 amu below). The 3% tolerance covers physical-constant
    conventions, which the CoRoT-24 b exponent of -16 amplifies an order
    of magnitude; the transcription errors this pin exists to catch move
    the rate by factors of several to decades.
    """
    # (Mp [MJup], Rp [RJup], Tp [K], a [au], Ms [Msun], target [kg/s], mu [amu])
    rows = [
        (0.01598, 0.146, 2411.0, 0.04272, 1.41, 2.62e12, 1.0),
        (0.018, 0.33, 1112.0, 0.05600, 0.91, 2.21e8, 2.0),
    ]
    for mp, rp, tp, a_au, ms, target, mu_amu in rows:
        m_p, r_p, m_s, a = mp * M_JUP, rp * R_JUP, ms * 1.989e30, a_au * AU
        mu_kg = mu_amu * amu
        rho = _jackson_photosphere(m_p, r_p, tp, mu_kg)
        rate, detail = nozzle_candidate(m_p, m_s, a, 0.0, rho, r_p, tp, mu_kg)
        assert not detail['saturated']
        assert math.isclose(rate, target, rel_tol=0.03, abs_tol=0.0)


@pytest.mark.reference_pinned
def test_lobe_filling_binaries_land_on_figure5():
    """Two Table 1 donors filling their lobes land on the Figure 5 solid curve.

    Jackson et al. (2017) Figure 5 plots the lobe-filling limit
    (``Phi_ph = Phi_R``, the exponential saturated at 1) for the Ritter
    (1988) Table A1 binaries with an 0.8 Msun accretor; the separation
    follows from inverting the Eggleton fit at the printed photospheric
    radius. The curve reads about 16 Msun/Gyr at M_d = 1.2 Msun and about
    25 at the M_d = 0.25 peak; 30% covers the figure-reading tolerance.
    """
    m_a = 0.8 * 1.989e30
    r_sun = 6.957e8
    # (M_d [Msun], r_ph [Rsun], T_eff [K], mu [amu], rho_ph [kg/m^3], curve [Msun/Gyr])
    rows = [
        (1.2, 1.17, 6480.0, 1.31, 2.5e-4, 16.0),
        (0.25, 0.25, 3410.0, 1.31, 1.6e-2, 25.0),
    ]
    for md, rph, t_eff, mu_amu, rho, curve in rows:
        m_d = md * 1.989e30
        r_ph = rph * r_sun
        q = m_d / m_a
        a = r_ph / (eggleton_lobe_radius(q, 1.0))
        rate, detail = nozzle_candidate(m_d, m_a, a, 0.0, rho, r_ph, t_eff, mu_amu * amu)
        # The lobe radius round-trips through the inverted fit to float
        # precision, so the exponent is zero up to that jitter rather than
        # the flag being exactly raised.
        assert abs(detail['exponent_applied']) < 1e-9
        assert math.isclose(rate * GYR / 1.989e30, curve, rel_tol=0.30, abs_tol=0.0)


@pytest.mark.reference_pinned
def test_curvature_equal_mass_value():
    """A(1) = 8, the equal-mass value Jackson et al. (2017) print.

    At q = 1 the Eq. (10) fit's denominator collapses to b1/4, so the fit
    returns their stated equal-mass curvature without fit error. At small
    q the fit's leading correction is b1 q^(1/3) with b1 = 2 * 3^(2/3),
    the second-order expansion of the L1 position; the in-text asymptotic
    of the paper carries half this coefficient and is not the fit.
    """
    assert math.isclose(curvature_a(1.0), 8.0, rel_tol=1e-12, abs_tol=0.0)
    q = 1e-9
    b1 = 2.0 * 3.0 ** (2.0 / 3.0)
    assert math.isclose(curvature_a(q) - 4.0, b1 * q ** (1.0 / 3.0), rel_tol=5e-3, abs_tol=0.0)


@pytest.mark.physics_invariant
def test_curvature_symmetric_and_bounded():
    """A(q) = A(1/q) and 4 < A <= 8 across twelve decades of mass ratio.

    The symmetry is their stated property of Eq. (10); the bounds follow
    from the denominator's minimum at q = 1.
    """
    for exponent in (-6.0, -3.0, -1.0, -0.3, 0.0, 0.3, 1.0, 3.0, 6.0):
        q = 10.0**exponent
        a_q = curvature_a(q)
        assert math.isclose(a_q, curvature_a(1.0 / q), rel_tol=1e-12, abs_tol=0.0)
        assert 4.0 < a_q <= 8.0


@pytest.mark.physics_invariant
def test_lobe_radius_inside_hill_radius_at_planetary_ratios():
    """The Eggleton lobe radius sits inside the Hill radius for q << 1.

    The volume-equivalent lobe is smaller than the L1 distance it is the
    volume average of, so the geometric Roche screen (which tests the Hill
    radius) fires no earlier than lobe contact does; at small q the ratio
    of the two tends to the constant 0.49 * 3^(1/3), about 0.71, which is
    the limit of the Eggleton fit against the Hill radius.
    """
    a = 0.05 * AU
    for exponent in (-7.0, -5.0, -3.0, -2.0):
        q = 10.0**exponent
        r_hill = a * (q / 3.0) ** (1.0 / 3.0)
        assert 0.0 < eggleton_lobe_radius(q, a) < r_hill
    q = 1e-7
    r_hill = a * (q / 3.0) ** (1.0 / 3.0)
    ratio = eggleton_lobe_radius(q, a) / r_hill
    assert math.isclose(ratio, 0.49 * 3.0 ** (1.0 / 3.0), rel_tol=5e-3, abs_tol=0.0)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_potential_at_lobe_matches_numerical_l1():
    """Eq. (14) at the Eggleton radius matches the exact Roche L1 potential.

    The exact corotating potential along the star-planet axis is maximized
    numerically between the bodies (the L1 saddle) and compared to the
    volume-averaged expansion evaluated at the Eggleton lobe radius, for a
    planetary mass ratio. Jackson et al. (2017) Figure 3 claims agreement
    to better than 2% away from contact; the barrier built from a deep
    photospheric level moves by less than 0.5% between the two.
    """
    m_d, m_a, a = 3.0 * Me, Ms, 0.07 * AU
    m_t = m_d + m_a
    omega2 = G * m_t / a**3
    x_cm = a * m_a / m_t

    def phi_axis(x):
        return -G * m_d / x - G * m_a / (a - x) - 0.5 * omega2 * (x - x_cm) ** 2

    lo = 0.3 * a * (m_d / (3.0 * m_a)) ** (1.0 / 3.0)
    hi = 3.0 * a * (m_d / (3.0 * m_a)) ** (1.0 / 3.0)
    for _ in range(200):
        third = (hi - lo) / 3.0
        if phi_axis(lo + third) < phi_axis(hi - third):
            lo = lo + third
        else:
            hi = hi - third
    phi_l1_exact = phi_axis(0.5 * (lo + hi))

    r_lobe = eggleton_lobe_radius(m_d / m_a, a)
    phi_l1_eq14 = volume_averaged_potential(r_lobe, m_d, m_a, a)
    assert math.isclose(phi_l1_eq14, phi_l1_exact, rel_tol=1e-4, abs_tol=0.0)
    r_ph = 1.6 * 2.2 * Re
    barrier_exact = phi_l1_exact - volume_averaged_potential(r_ph, m_d, m_a, a)
    barrier_eq14 = phi_l1_eq14 - volume_averaged_potential(r_ph, m_d, m_a, a)
    assert math.isclose(barrier_eq14, barrier_exact, rel_tol=5e-3, abs_tol=0.0)


@pytest.mark.physics_invariant
def test_rate_invariant_under_launch_level_choice():
    """The rate does not depend on which isothermal level launches it.

    Along an isothermal hydrostatic column, rho(r) exp(Phi(r)/v_th^2) is
    constant, so the Bernoulli structure of Eq. (3) cancels the level
    choice; the residual is the size of the tidal terms between the two
    levels, far below 1% for a donor deep inside its lobe.
    """
    m_p, r_p, t = Me, Re, 500.0
    prof = isothermal_profile(m_p, r_p, t, {'CO2': 1.0}, 1e7, 1e-3)
    rates = []
    for p_level in (2000.0, 2.0):
        lev = interp_at_pressure(prof, p_level)
        rate, _ = nozzle_candidate(
            m_p, Ms, 0.1 * AU, 0.0, lev['rho'], lev['r'], lev['T'], lev['mmw']
        )
        rates.append(rate)
    assert rates[0] > 0.0
    assert math.isclose(rates[0], rates[1], rel_tol=1e-2, abs_tol=0.0)


@pytest.mark.physics_invariant
def test_exponent_clamps_at_lobe_contact():
    """At and beyond lobe contact the rate is the lobe-filling boundary value.

    A launch level at the lobe radius saturates the exponential at 1
    (their Figure 5 case); a level beyond it, where the Eq. (14) expansion
    has diverged, returns the same clamped rate rather than an unphysical
    amplification.
    """
    m_p, m_s, a = 0.5 * M_JUP, 1.989e30, 0.015 * AU
    rho, t, mu_kg = 1e-5, 1500.0, 2.0 * amu
    r_lobe = eggleton_lobe_radius(m_p / m_s, a)
    at_contact, d_contact = nozzle_candidate(m_p, m_s, a, 0.0, rho, r_lobe, t, mu_kg)
    beyond, d_beyond = nozzle_candidate(m_p, m_s, a, 0.0, rho, 1.2 * r_lobe, t, mu_kg)
    assert d_contact['saturated'] and d_beyond['saturated']
    assert d_contact['exponent_applied'] == 0.0
    v_th = d_contact['v_th']
    boundary_value = rho * math.exp(-0.5) * v_th * d_contact['area_m2']
    assert math.isclose(at_contact, boundary_value, rel_tol=1e-12, abs_tol=0.0)
    assert math.isclose(beyond, at_contact, rel_tol=1e-12, abs_tol=0.0)


def test_curvature_rejects_unphysical_mass_ratio():
    """A non-positive or non-finite mass ratio raises; a valid one returns.

    Zero, negative, and non-finite ratios have no L1 geometry and must
    raise rather than return a complex root or propagate a NaN. A valid
    planetary ratio on the same path stays inside the fit's bounds.
    """
    for bad in (0.0, -1.0, math.inf, math.nan):
        with pytest.raises(ValueError, match='mass ratio'):
            curvature_a(bad)
    ok = curvature_a(1e-5)
    assert math.isfinite(ok)
    assert 4.0 < ok <= 8.0


def test_public_surface_rejects_unphysical_arguments():
    """Every exported function raises on unphysical input, not silently.

    The module is documented public API, so a direct caller must not get a
    negative mass-loss rate from a negative density, a complex cube root
    from a negative mass ratio, or a bare division by zero from a launch
    level at the origin. Each case below returned one of those before the
    guards existed.
    """
    m_p, m_s, a = 3.0 * Me, Ms, 0.05 * AU
    ok = dict(
        M_p=m_p, M_star=m_s, a=a, e=0.0, rho_ph=1e-6, r_ph=2.0 * Re, T=1200.0, mu_kg=2.3 * amu
    )
    for bad, match in (
        (dict(rho_ph=-1.0), 'rho_ph'),
        (dict(e=1.5), 'e must satisfy'),
        (dict(e=-0.1), 'e must satisfy'),
        (dict(T=-100.0), 'T'),
        (dict(r_ph=0.0), 'r_ph'),
        (dict(M_p=math.inf), 'M_p'),
    ):
        with pytest.raises(ValueError, match=match):
            nozzle_candidate(**{**ok, **bad})
    for bad in (0.0, -1.0, math.nan):
        with pytest.raises(ValueError, match='mass ratio'):
            eggleton_lobe_radius(bad, a)
    with pytest.raises(ValueError, match='r_v'):
        volume_averaged_potential(-1.0, m_p, m_s, a)
    with pytest.raises(ValueError, match='r'):
        isothermal_column_density(1e-6, 1e7, 0.0, m_p, 1e4)
    # A valid call on the same path still returns a positive finite rate.
    rate, _ = nozzle_candidate(**ok)
    assert math.isfinite(rate)
    assert rate > 0.0


def test_orbit_average_rejects_an_empty_quadrature():
    """A non-positive phase count raises rather than dividing by zero.

    The orbit average sums Kepler weights over ``n_phase`` midpoint nodes,
    so a count below one has no weight to divide by. A single node is
    legal and, on a circular orbit, exact.
    """
    m_p, m_s, a = 3.0 * Me, Ms, 0.05 * AU
    args = (m_p, m_s, a, 0.0, 1e-6, 2.0 * Re, 1200.0, 2.3 * amu)
    for bad in (0, -1):
        with pytest.raises(ValueError, match='n_phase'):
            nozzle_candidate(*args, n_phase=bad)
    one, _ = nozzle_candidate(*args, n_phase=1)
    many, _ = nozzle_candidate(*args, n_phase=64)
    assert one == pytest.approx(many, rel=1e-12)
