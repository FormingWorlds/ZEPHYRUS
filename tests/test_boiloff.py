"""Tests for ``src/zephyrus/boiloff.py``.

Exercises the bolometrically driven escape branch. The physical invariants
under test:

- Reference pins: the Lambert-W Mach number is exactly 1 when the launch
  level sits at the Bondi radius (the analytical sonic-point limit of Owen
  & Wu 2016), and collapses by more than six decades by ``x = 0.1``, their
  published shutoff; the restricted Jeans parameter obeys the identity
  ``Lambda = 2 R_B / R_p`` for every composition.
- Monotonicity / boundedness: the Mach number falls monotonically as the
  launch level retreats inside the Bondi radius; the candidate rate never
  exceeds any of its caps.
- The luminosity cap divides by the tidally reduced barrier, so it rises by
  the Erkaev et al. (2007) enhancement factor of their printed Table 1 and
  reduces to the untidal form at K = 1.
- Error contract: ``parker_mach`` rejects arguments outside ``(0, 1]``; the
  timescale diagnostic reports "not evaluated" without reservoirs.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.boiloff import (
    LAMBDA_BAND,
    bolometric_candidate,
    lambda_restricted,
    parker_mach,
    tang_timescale_check,
)
from zephyrus.constants import G, amu, kb
from zephyrus.hydrodynamic import k_tide
from zephyrus.planets_parameters import Me, Re
from zephyrus.profiles import interp_at_pressure, isothermal_profile

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _launch(M_p, R_p, T_eq, comp, p=2000.0):
    """Photospheric-type launch level on an isothermal test atmosphere."""
    prof = isothermal_profile(M_p, R_p, T_eq, comp, 1e7, 1e-3)
    return interp_at_pressure(prof, min(p, float(prof.p[0])))


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_parker_mach_sonic_limit_and_shutoff():
    """The Mach number is 1 at the Bondi radius and shuts off by x = 0.1.

    Analytical limit: ``f(1) = e^-1`` and ``W0(-1/e) = -1``, so the launch
    level at the Bondi radius is exactly sonic. Published shutoff: Owen &
    Wu (2016) place the end of boil-off at ``R_p / R_B = 0.1``, where the
    exponential factor has collapsed the rate; the Mach ratio between
    ``x = 1`` and ``x = 0.1`` must exceed six decades. Monotonicity: the
    Mach number falls strictly as ``x`` decreases.
    """
    assert parker_mach(1.0) == pytest.approx(1.0, abs=1e-6)
    xs = [1.0, 0.5, 0.3, 0.2, 0.1]
    ms = [parker_mach(x) for x in xs]
    assert all(a > b for a, b in zip(ms, ms[1:]))
    assert ms[-1] / ms[0] < 1e-6
    # Sign and scale: Mach numbers are subsonic below the Bondi radius.
    assert all(0.0 < m <= 1.0 for m in ms)


def test_parker_mach_rejects_out_of_domain():
    """Arguments outside (0, 1] raise; the caller owns the inflated clamp.

    Zero, negative, and superunity launch ratios are not meaningful inputs
    to the transonic solution and must raise rather than return a complex
    or extrapolated Mach number. A valid argument on the same path returns
    a finite subsonic value.
    """
    for bad in (0.0, -0.5, 1.5):
        with pytest.raises(ValueError, match='parker_mach'):
            parker_mach(bad)
    ok = parker_mach(0.4)
    assert math.isfinite(ok)
    assert 0.0 < ok < 1.0


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_lambda_equals_two_bondi_radii_over_rp():
    """The restricted Jeans parameter obeys Lambda = 2 R_B / R_p exactly.

    For isothermal gas at the same temperature the identity holds for
    every mean molecular mass, which is what makes the Owen & Wu (2016)
    shutoff at ``R_p / R_B = 0.1`` equal to ``Lambda = 20`` for every
    composition. Checked for H2-, steam-, and CO2-like mean masses; the
    literature band around the threshold brackets 20.
    """
    M_p, R_p, T_eq = 5 * Me, 2 * Re, 800.0
    for mu_amu in (2.3, 18.0, 44.0):
        mu = mu_amu * amu
        lam = lambda_restricted(M_p, R_p, T_eq, mu)
        c2 = kb * T_eq / mu
        R_B = G * M_p / (2.0 * c2)
        assert lam == pytest.approx(2.0 * R_B / R_p, rel=1e-12)
    assert LAMBDA_BAND[0] < 20.0 < LAMBDA_BAND[1]
    # Monotone in mu: heavier gas is more tightly bound.
    lam_light = lambda_restricted(M_p, R_p, T_eq, 2.3 * amu)
    lam_heavy = lambda_restricted(M_p, R_p, T_eq, 44.0 * amu)
    assert lam_heavy > lam_light


@pytest.mark.physics_invariant
def test_luminosity_cap_applies_only_past_the_gate():
    """The interior-luminosity cap joins the candidate only after boil-off.

    While the configuration is inflated (``lambda_gate`` below threshold)
    the candidate is min(Parker, Bondi) and no luminosity cap is computed;
    past the gate the cap joins and bounds the returned rate. The rate must
    never exceed any active cap (boundedness).
    """
    M_p, R_p, T_eq = 3 * Me, 3 * Re, 1000.0
    launch = _launch(M_p, R_p, T_eq, {'H2': 1.0})
    rate_a, det_a = bolometric_candidate(
        M_p, R_p, T_eq, 0.01, launch, 1.0, lambda_gate=10.0, lambda_crit=20.0
    )
    assert det_a['active'] is True
    assert det_a['mdot_luminosity'] is None
    assert rate_a <= min(det_a['mdot_parker'], det_a['mdot_bondi']) * (1 + 1e-12)
    rate_b, det_b = bolometric_candidate(
        M_p, R_p, T_eq, 0.01, launch, 1.0, lambda_gate=30.0, lambda_crit=20.0
    )
    assert det_b['active'] is False
    assert det_b['mdot_luminosity'] is not None
    assert rate_b <= det_b['mdot_luminosity'] * (1 + 1e-12)
    # The luminosity cap can only reduce the candidate, never raise it.
    assert rate_b <= rate_a * (1 + 1e-12)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_luminosity_cap_carries_the_tidal_barrier_reduction():
    """The cap divides by the tidally reduced barrier, K from Erkaev Eq. (17).

    The cap is an interior luminosity divided by the work per unit mass
    needed to lift gas out, ``g R_p K``, the same barrier and the same
    reference radius the energy-limited rate divides by. So a tidally
    reduced barrier raises the cap by the enhancement factor ``1 / K``,
    whose value at ``xi = 3`` is the 1.92 printed in Erkaev et al. (2007),
    A&A 472, 329, Table 1. Discrimination: multiplying by ``K`` instead of
    dividing would lower the cap rather than raise it, and ``K = 1`` must
    reproduce the untidal value identically.
    """
    M_p, R_p, T_eq = 3 * Me, 3 * Re, 1000.0
    launch = _launch(M_p, R_p, T_eq, {'H2': 1.0})
    args = (M_p, R_p, T_eq, 0.01, launch, 1.0)
    _rate_flat, det_flat = bolometric_candidate(*args, lambda_gate=30.0, lambda_crit=20.0)
    k = k_tide(3.0)
    _rate_tidal, det_tidal = bolometric_candidate(
        *args, lambda_gate=30.0, lambda_crit=20.0, k_tide=k
    )
    assert det_flat['k_tide'] == pytest.approx(1.0)
    assert det_tidal['mdot_luminosity'] > det_flat['mdot_luminosity']
    assert det_tidal['mdot_luminosity'] / det_flat['mdot_luminosity'] == pytest.approx(
        1.92, rel=0.01
    )
    # K = 1 is the untidal form, and the cap still bounds the rate.
    _rate_one, det_one = bolometric_candidate(
        *args, lambda_gate=30.0, lambda_crit=20.0, k_tide=1.0
    )
    assert det_one['mdot_luminosity'] == pytest.approx(det_flat['mdot_luminosity'], rel=1e-12)
    assert _rate_tidal <= det_tidal['mdot_luminosity'] * (1 + 1e-12)


def test_wind_temperature_and_opacity_scaling():
    """The wind runs at T_eq / 2^(1/4) and the rate scales as 1 / kappa.

    The wind temperature is the Misener et al. (2025) recommendation for
    the isothermal formulas; the Parker rate carries the photospheric
    opacity inversely, so doubling kappa halves the uncapped Parker rate.
    """
    M_p, R_p, T_eq = 4 * Me, 2.5 * Re, 1000.0
    launch = _launch(M_p, R_p, T_eq, {'H2': 1.0})
    _, det = bolometric_candidate(M_p, R_p, T_eq, 0.01, launch, 1.0, 5.0, 20.0)
    assert det['T_wind'] == pytest.approx(1000.0 / 2**0.25, rel=1e-12)
    _, det2 = bolometric_candidate(M_p, R_p, T_eq, 0.02, launch, 1.0, 5.0, 20.0)
    assert det2['mdot_parker'] == pytest.approx(det['mdot_parker'] / 2.0, rel=1e-12)
    # The Bondi cap does not depend on the opacity.
    assert det2['mdot_bondi'] == pytest.approx(det['mdot_bondi'], rel=1e-12)


def test_inflated_launch_level_clamps_with_flag():
    """A launch level beyond the Bondi radius clamps to sonic, flagged.

    A very hot, loosely bound configuration puts the photosphere outside
    the Bondi radius; the branch must clamp ``x`` to 1 (Mach 1) and raise
    ``bondi_inflated`` rather than raise an exception, because such states
    are physically posed inputs on the boil-off side of the gate. The
    launch level is built by hand: a hydrostatic profile bound at every
    level cannot reach past the Bondi radius, which is the very reason the
    clamp exists for externally supplied levels.
    """
    M_p, R_p, T_eq = 0.6 * Me, 1.8 * Re, 1800.0
    mu = 2.3 * amu
    c2 = kb * (T_eq / 2**0.25) / mu
    r_bondi = G * M_p / (2.0 * c2)
    launch = {'r': 1.2 * r_bondi, 'mmw': mu, 'rho': 1e-6}  # plausible photosphere
    assert launch['r'] > r_bondi  # genuinely beyond R_B
    rate, det = bolometric_candidate(M_p, R_p, T_eq, 0.01, launch, 1.0, 3.0, 20.0)
    assert det['flags'].get('bondi_inflated') is True
    assert det['mach'] == pytest.approx(1.0, abs=1e-9)
    assert math.isfinite(rate)
    assert rate > 0.0


def test_tang_timescale_diagnostic_contract():
    """The termination diagnostic evaluates only with reservoirs supplied.

    Without reservoir masses (or with a zero rate) it reports
    ``evaluated: False``; with them, both timescales are pinned against
    hand-evaluated values (``t_mdot = M_env / Mdot = 1e13 s`` and
    ``t_cool = G M_p M_env / (R_p L) = 1.2226e11 s`` for the Earth-like
    inputs), so a swapped pair or a reversed comparison fails: the slow
    rate is terminated (``t_mdot >= t_cool``) and a fast rate on the other
    side of the inequality is not. A heavier envelope at a fixed rate
    cannot flip the verdict (both timescales scale with the envelope mass).
    """
    assert tang_timescale_check(Me, Re, 1.0, 1e5, None) == {'evaluated': False}
    assert tang_timescale_check(Me, Re, 1.0, 0.0, {'H': 1e18}) == {'evaluated': False}
    out = tang_timescale_check(Me, Re, 1.0, 1e5, {'H': 1e18})
    assert out['evaluated'] is True
    assert out['t_mdot_s'] == pytest.approx(1e13, rel=1e-9)
    assert out['t_cool_s'] == pytest.approx(1.2226e11, rel=1e-3)
    # Swap discrimination: the two timescales differ by two decades here.
    assert out['t_mdot_s'] > 50.0 * out['t_cool_s']
    assert out['terminated'] is True
    fast = tang_timescale_check(Me, Re, 1.0, 1e9, {'H': 1e18})
    assert fast['t_mdot_s'] < fast['t_cool_s']
    assert fast['terminated'] is False
    out2 = tang_timescale_check(Me, Re, 1.0, 1e5, {'H': 2e18})
    assert out2['terminated'] == out['terminated']
