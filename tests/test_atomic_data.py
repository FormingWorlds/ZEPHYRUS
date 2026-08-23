"""Tests for ``src/zephyrus/atomic_data.py``.

Exercises the transcribed cooling data and the closed-form rate
coefficients. The physical invariants under test:

- Reference pins: spot values of the Nakayama et al. (2022) Appendix C
  transcription, including the corrected O+ level set; the Badnell (2006)
  recombination fit against its printed coefficients, with an explicit
  discrimination against the garbled variant printed by a later source; the
  Murray-Clay et al. (2009) hydrogen case B value.
- Monotonicity / boundedness: recombination coefficients fall with
  temperature; both cooling channels are positive and rise with temperature
  where their level spacings dictate.
- Closed form: the CO2 band reduces to the coronal (collision-limited) form
  far below the critical density, independent of the radiative constants.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.atomic_data import (
    BADNELL_N,
    CO2_KD,
    HNU_15UM,
    THREE_LEVEL,
    alpha_case_b,
    badnell_alpha_rr,
    co2_band_cooling,
    o_finestructure_cooling,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.reference_pinned
def test_three_level_transcription_spot_values():
    """Spot entries reproduce the printed Nakayama et al. (2022) tables.

    The N transauroral and auroral Einstein coefficients, the hydrogen
    Lyman-alpha coefficient, and the O+ level set (the 4S-2D-2P system with
    weights 4, 10, 6, correcting the printed row whose labels carry the
    neutral-O configuration) pin the transcription. The zero-A entries on
    radiatively forbidden couplings are structural, not missing data.
    """
    n = THREE_LEVEL['N']
    assert n['transitions'][(1, 3)][0] == pytest.approx(5.22e-3, rel=1e-12, abs=0.0)
    assert n['transitions'][(2, 3)][0] == pytest.approx(8.47e-2, rel=1e-12, abs=0.0)
    h = THREE_LEVEL['H']
    assert h['transitions'][(1, 3)][0] == pytest.approx(6.26e8, rel=1e-12, abs=0.0)
    op = THREE_LEVEL['O+']
    assert [lv[1] for lv in op['levels']] == [4, 10, 6]
    assert [lv[0] for lv in op['levels']] == ['4S', '2D', '2P']
    # Structural zeros: the H 2s-2p and C+ 4P-2D couplings carry no A value.
    assert THREE_LEVEL['H']['transitions'][(2, 3)][0] == pytest.approx(0.0, abs=0.0)
    assert THREE_LEVEL['C+']['transitions'][(2, 3)][0] == pytest.approx(0.0, abs=0.0)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_badnell_fit_magnitude_slope_and_misprint_guard():
    """The Badnell recombination fit has the right magnitude and slope.

    The nitrogen coefficient at 1e4 K is pinned at 3.761e-13 cm^3 s^-1,
    hand-evaluated from the printed fit form with the shipped coefficients
    (beside the hydrogen case B 2.7e-13, the order anchor), and the
    coefficient falls monotonically with temperature. Two discrimination
    guards: swapping the T0 and T1 coefficients (adjacent tuple entries)
    gives 4.82e-13, outside the pin; and the garbled rendering of the fit
    printed by Chatterjee & Pierrehumbert (2026, their Eq. 35: product
    turned into a sum, one exponent repeated, exponential argument
    inverted) disagrees with the implemented original by more than a
    factor 2 at 1e4 K.
    """
    a4 = badnell_alpha_rr(1.0e4)
    assert a4 == pytest.approx(3.761e-13, rel=0.02, abs=0.0)
    # Swapped-coefficient discrimination: T0 and T1 interchanged gives
    # 4.82e-13, a 28 percent shift, far outside the 2 percent pin.
    assert abs(a4 - 4.82e-13) > 0.2 * a4
    assert 1e-13 < a4 < 1e-12
    assert badnell_alpha_rr(3.0e4) < a4 < badnell_alpha_rr(3.0e3)
    t0, t1, t2, a_fit, b_fit, c_fit = BADNELL_N
    T = 1.0e4
    expo = 1.0 - b_fit - c_fit * math.exp(-T / t2)
    garbled = a_fit / (
        math.sqrt(T / t0) * (1.0 + math.sqrt(T / t0)) ** expo
        + (1.0 + math.sqrt(T / t1)) ** expo
    )
    assert garbled / a4 > 2.0
    # Transcription pin: the six coefficients are the Z = 7, N = 6 row of
    # Badnell's table, not the copy in the secondary that garbles the form.
    # T2 is the digit-transposition trap, printed 6.739e4 and once carried
    # here as 6.379e4, so it is pinned exactly rather than through a rate.
    assert (t0, t1, t2) == (9.467e-2, 2.954e6, 6.739e4)
    assert (a_fit, b_fit, c_fit) == (6.387e-10, 0.7308, 0.2440)
    transposed = (t0, t1, 6.379e4, a_fit, b_fit, c_fit)
    assert badnell_alpha_rr(5.0e4, transposed) / badnell_alpha_rr(5.0e4) > 1.02


def test_case_b_coefficients_and_temperature_scaling():
    """Case B values are pinned and share the hydrogen temperature exponent.

    Hydrogen carries the Murray-Clay et al. (2009) 2.7e-13 cm^3 s^-1 at
    1e4 K; every element scales as (T/1e4 K)^-0.9 exactly (the documented
    extension of hydrogen's exponent to the heavies); an element without a
    tabulated value falls back to the atomic-O coefficient rather than
    raising, because coefficient provenance is flagged upstream.
    """
    assert alpha_case_b('H', 1.0e4) == pytest.approx(2.7e-13, rel=1e-12, abs=0.0)
    for el in ('H', 'He', 'C', 'N', 'O'):
        ratio = alpha_case_b(el, 2.0e4) / alpha_case_b(el, 1.0e4)
        assert ratio == pytest.approx(2.0**-0.9, rel=1e-12, abs=0.0)
    assert alpha_case_b('Xe', 8000.0) == pytest.approx(alpha_case_b('O', 8000.0), rel=1e-12, abs=0.0)
    # Scale guard: all case B values live in the 1e-13 decade at 1e4 K.
    for el in ('He', 'C', 'N', 'O'):
        assert 5e-14 < alpha_case_b(el, 1.0e4) < 5e-13


@pytest.mark.physics_invariant
def test_co2_band_coronal_limit_and_detailed_balance():
    """Far below the critical density the band cooling is collision limited.

    In the coronal limit every collisional excitation radiates, so the full
    expression must reduce to ``h nu k_e n_M n_CO2`` with the excitation
    rate fixed by detailed balance, ``k_e = 2 k_d exp(-667 K / T)``,
    independent of the Einstein coefficient and the escape probability.
    Zero colliders give exactly zero cooling (the error-contract limit).
    """
    T = 300.0
    n_co2, colliders = 1e6, {'O': 1e6}
    q = co2_band_cooling(n_co2, colliders, T)
    a, b = CO2_KD['O']
    kd = a * T**b
    ke = 2.0 * kd * math.exp(-667.0 / T)
    assert q == pytest.approx(HNU_15UM * ke * colliders['O'] * n_co2, rel=1e-3, abs=0.0)
    assert co2_band_cooling(n_co2, {}, T) == pytest.approx(0.0, abs=0.0)
    assert co2_band_cooling(n_co2, {'O': 0.0}, T) == pytest.approx(0.0, abs=0.0)


def test_co2_band_escape_probability_branches():
    """An overlying CO2 column reduces the cooling through photon trapping.

    The zero-column non-LTE ceiling has escape probability 0.5; a small
    column enters the shallow branch and a large column the steep branch,
    each strictly reducing the cooling relative to the ceiling, and the
    large-column case more strongly (monotone in column).
    """
    T, n_co2, colliders = 300.0, 1e10, {'CO2': 1e10, 'O': 1e8}
    q0 = co2_band_cooling(n_co2, colliders, T, col_co2=0.0)
    q_small = co2_band_cooling(n_co2, colliders, T, col_co2=1e14)  # sigma N ~ 0.64
    q_large = co2_band_cooling(n_co2, colliders, T, col_co2=1e16)  # sigma N ~ 64
    assert q0 > 0.0
    assert q_large < q_small
    # Trapping can only reduce the loss relative to the free-escape
    # ceiling, on both branches.
    assert q_small < q0
    assert q_large < q0


@pytest.mark.physics_invariant
def test_o_finestructure_positive_and_activating():
    """The O fine-structure channel is positive and thermally activating.

    The 228 K and 326 K level spacings make the cooling rise steeply from
    150 K to 300 K; it is linear in the atomic-O density (a two-point ratio
    check), and zero density gives zero cooling.
    """
    q150 = o_finestructure_cooling(1e8, 150.0)
    q300 = o_finestructure_cooling(1e8, 300.0)
    assert q150 > 0.0
    assert q300 > q150
    assert o_finestructure_cooling(2e8, 300.0) / q300 == pytest.approx(2.0, rel=1e-12, abs=0.0)
    assert o_finestructure_cooling(0.0, 300.0) == pytest.approx(0.0, abs=0.0)
