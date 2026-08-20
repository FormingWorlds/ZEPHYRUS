"""Tests for ``src/zephyrus/diagnostics.py``.

Exercises the reporting quantities that accompany every dispatch verdict.
The physical invariants under test:

- Closed forms: the Johnson et al. (2013) transonic energy criterion obeys
  its published scalings; the Guo (2024) triple reduces correctly in the
  wide-orbit and Roche limits; the Erkaev et al. (2007) critical exobase
  temperature recovers its Jupiter normalization and vanishes at the Roche
  lobe.
- Monotonicity / boundedness: the along-profile fluid check reports the
  worst local Knudsen number with its truncation declared; the potential
  screens classify the three regimes in the right order.
- Error contract: the self-consistency screen reports "not evaluated"
  without its optional inputs rather than guessing.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import pytest

from zephyrus.constants import G, amu
from zephyrus.diagnostics import (
    CALDIROLI_THRESHOLD_LOG_PHI,
    DAYSIDE_FACTORS,
    MURRAY_CLAY_EXPONENTS,
    SALZ_SCREEN_LOG_PHI,
    along_profile_fluid_check,
    erkaev_tc,
    guo_triple,
    potential_screens,
    q_net_over_qc,
    self_consistency_screen,
)
from zephyrus.planets_parameters import Me, Mjup, Ms, Re, Rjup
from zephyrus.profiles import isothermal_profile

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
def test_johnson_criterion_scalings():
    """The transonic energy criterion carries its published dependences.

    The net power is linear in the XUV flux and quadratic in the XUV
    radius; the critical power is independent of both, so the ratio
    inherits the linearity (a two-point check). A larger cross section
    lowers the critical power (a more collisional gas is easier to drive
    transonic), raising the ratio.
    """
    args = dict(
        eps=0.1,
        R_xuv=1.2 * Re,
        r_sonic=3.0 * Re,
        r_base=1.1 * Re,
        M_p=5 * Me,
        m_mean=1.008 * amu,
        sigma_c=1e-19,
    )
    r1, q_net1, q_c1 = q_net_over_qc(F_xuv=10.0, **args)
    r2, q_net2, q_c2 = q_net_over_qc(F_xuv=20.0, **args)
    assert q_net2 == pytest.approx(2.0 * q_net1, rel=1e-12)
    assert q_c2 == pytest.approx(q_c1, rel=1e-12)
    assert r2 == pytest.approx(2.0 * r1, rel=1e-12)
    args2 = dict(args, sigma_c=2e-19)
    r3, _, q_c3 = q_net_over_qc(F_xuv=10.0, **args2)
    assert q_c3 == pytest.approx(q_c1 / 2.0, rel=1e-12)
    assert r3 > r1


@pytest.mark.physics_invariant
def test_guo_triple_limits():
    """The regime triple reduces correctly in its limits.

    On a wide orbit the Roche correction vanishes, so lambda* approaches
    lambda from below; inside the Roche limit (correction factor at or
    below its root) lambda* reports 0. The exobase value passes through
    untouched.
    """
    mu = 2.3 * amu
    wide = guo_triple(5 * Me, 1.5 * Re, 800.0, mu, Ms, 1.496e11, 0.0, lambda_exo=12.3)
    assert wide['lambda_exo'] == pytest.approx(12.3, rel=1e-12)
    assert wide['lambda_star'] < wide['lambda_rp']
    assert wide['lambda_star'] == pytest.approx(wide['lambda_rp'], rel=1e-2)
    lam_ref = G * (5 * Me) * mu / (1.380649e-23 * 800.0 * 1.5 * Re)
    assert wide['lambda_rp'] == pytest.approx(lam_ref, rel=1e-9)
    # Deep inside the Roche limit the corrected parameter reports zero.
    close = guo_triple(5 * Me, 1.5 * Re, 800.0, mu, Ms, 5e8, 0.0, lambda_exo=12.3)
    assert close['lambda_star'] == pytest.approx(0.0, abs=0.0)


@pytest.mark.reference_pinned
def test_erkaev_critical_temperature_normalization():
    """The critical exobase temperature recovers its Jupiter normalization.

    Erkaev et al. (2007) normalize to 1.45e5 K for Jupiter values; with the
    exobase at the planetary radius and the Hill radius far away the
    correction factor approaches 1, so the critical temperature approaches
    the normalization itself. At the Roche lobe it must vanish (no barrier
    left), and a farther exobase always lowers it.
    """
    t_far = erkaev_tc(Mjup, Rjup, 1.0 * Rjup, 1e3 * Rjup)
    assert t_far == pytest.approx(1.45e5, rel=1e-2)
    assert erkaev_tc(Mjup, Rjup, 2.0 * Rjup, 1e3 * Rjup) < t_far
    assert erkaev_tc(Mjup, Rjup, 3.0 * Rjup, 2.5 * Rjup) == pytest.approx(0.0, abs=0.0)
    # Mass scaling is linear: twice the mass doubles the barrier.
    assert erkaev_tc(2 * Mjup, Rjup, 1.0 * Rjup, 1e3 * Rjup) == pytest.approx(
        2 * t_far, rel=1e-2
    )


def test_along_profile_fluid_check_reports_truncation():
    """The fluid check walks the covered levels and declares its truncation.

    On a bound isothermal profile with the sonic radius above the top, all
    levels are checked, the worst Knudsen number is deep in the fluid
    regime (dense gas), and the truncation flag is set. A sonic radius
    inside the profile checks only the levels below it.
    """
    prof = isothermal_profile(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, 1e7, 1e-4)
    out = along_profile_fluid_check(prof, 5 * Me, r_sonic=10 * Re)
    assert out['levels_checked'] == len(prof.p)
    assert out['truncated_at_profile_top'] is True
    assert out['fluid'] is True
    # The worst level is the rarefied top, still an order below threshold.
    assert out['worst_kn'] < 0.1
    r_mid = float(prof.r[len(prof.p) // 2])
    out2 = along_profile_fluid_check(prof, 5 * Me, r_sonic=r_mid)
    assert out2['levels_checked'] < len(prof.p)
    assert out2['truncated_at_profile_top'] is False


def test_self_consistency_screen_contract():
    """The snapshot screen evaluates only with reservoirs and age supplied.

    Missing inputs report "not evaluated"; with them, a rate that would
    empty the reservoirs faster than the age flags the snapshot, and a
    slow rate does not.
    """
    assert self_consistency_screen(None, 1e5, 1e16) == {'evaluated': False}
    assert self_consistency_screen({'H': 1e18}, 1e5, None) == {'evaluated': False}
    assert self_consistency_screen({'H': 1e18}, 0.0, 1e16) == {'evaluated': False}
    fast = self_consistency_screen({'H': 1e18}, 1e5, 1e16)
    assert fast['evaluated'] is True
    assert fast['inconsistent'] is True  # 1e13 s depletion against 1e16 s age
    slow = self_consistency_screen({'H': 1e18}, 1e-5, 1e16)
    assert slow['inconsistent'] is False


def test_potential_screens_classify_in_order():
    """The threshold-potential screens order the three verdicts correctly.

    An Earth-like potential sits far below both screens (wind side); a
    compact massive planet lands above the upper screen (no-wind side); an
    intermediate case falls between. The two screen bands themselves must
    be ordered and overlapping in the documented way.
    """
    low = potential_screens(Me, Re)
    assert low['salz_verdict'] == 'wind'
    assert low['above_caldiroli'] is False
    high = potential_screens(10 * Mjup, 1.0 * Rjup)
    assert high['salz_verdict'] == 'no-wind'
    assert high['above_caldiroli'] is True
    mid = potential_screens(2.2 * Mjup, 1.05 * Rjup)
    assert mid['salz_verdict'] == 'intermediate'
    # Screen geometry: the efficiency band starts below the wind screen.
    assert CALDIROLI_THRESHOLD_LOG_PHI[0] < SALZ_SCREEN_LOG_PHI[0]
    assert SALZ_SCREEN_LOG_PHI[0] < SALZ_SCREEN_LOG_PHI[1]


def test_documentation_constants_are_complete():
    """The documented comparison constants carry the published values.

    The Murray-Clay et al. (2009) numerical flux exponents (0.6 and 0.9)
    against the analytic ones carried by the implementation (0.5 and 1.0),
    and their dayside reduction factors, are reporting constants consumers
    rely on; pin them so a silent edit fails.
    """
    assert MURRAY_CLAY_EXPONENTS['RR_numerical'] == pytest.approx(0.6, rel=1e-12)
    assert MURRAY_CLAY_EXPONENTS['EL_numerical'] == pytest.approx(0.9, rel=1e-12)
    assert MURRAY_CLAY_EXPONENTS['RR_analytic_inherited'] == pytest.approx(0.5, rel=1e-12)
    assert MURRAY_CLAY_EXPONENTS['EL_analytic_inherited'] == pytest.approx(1.0, rel=1e-12)
    assert DAYSIDE_FACTORS['energy_limited'] == pytest.approx(0.26, rel=1e-12)
    assert DAYSIDE_FACTORS['recombination_limited'] == pytest.approx(0.31, rel=1e-12)
    # The reduction factors are genuine reductions.
    assert all(0.0 < v < 1.0 for v in DAYSIDE_FACTORS.values())
