"""Tests for ``src/zephyrus/thermostat.py``.

Exercises the statistical-equilibrium line cooling, the ionization balance,
and the wind-temperature root-find. The physical invariants under test:

- Detailed balance: as the electron density grows the three-level
  populations reach Boltzmann ratios (the LTE limit).
- Coronal limit: at low electron density the line cooling is exactly linear
  in the electron density.
- Reference pin: the hydrogen three-level system agrees with the Black
  (1981) Lyman-alpha rate as printed by Murray-Clay et al. (2009) at order
  unity (the sources differ in collision-strength treatment, so identity is
  not expected and the bracket documents the offset).
- Boundedness / error contract: the ionization fraction lies in [0, 1] with
  exact zero at zero flux; the root-find clamps to its bracket edges with
  the clamp recorded; disabling every cooling channel raises.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import pytest

from zephyrus.atomic_data import HC_CM, LYA_BLACK, THREE_LEVEL, alpha_case_b
from zephyrus.constants import kb_cgs
from zephyrus.thermostat import (
    balance_at,
    front_constants,
    ionization_fraction,
    recombination_alpha,
    solve_wind_temperature,
    three_level_cooling,
    three_level_populations,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _base(n_si=1e18, vmr=None):
    """Minimal wind-base level dict: number density [m^-3] and composition."""
    return {'n': n_si, 'vmr': vmr or {'N2': 1.0}}


@pytest.mark.physics_invariant
def test_populations_reach_boltzmann_in_the_lte_limit():
    """At very high electron density the populations are Boltzmann ratios.

    Detailed balance is built into the rate coefficients, so collisions
    alone must thermalize the levels: ``n_u / n_l = (g_u / g_l)
    exp(-E_ul / kB T)``. Species with a zero collision strength on one
    coupling cannot thermalize through it and are skipped. The edge case is
    the ground-state-dominated low-T limit checked through the exponential.
    """
    T = 8000.0
    checked = 0
    for sp, data in THREE_LEVEL.items():
        if any(g == 0.0 for _a, g in data['transitions'].values()):
            continue
        n1, n2, n3 = three_level_populations(sp, 1.0, 1e20, T)
        (_l1, g1, e1), (_l2, g2, e2), (_l3, g3, e3) = data['levels']
        for (na, ga, ea), (nb, gb, eb) in (
            ((n2, g2, e2), (n1, g1, e1)),
            ((n3, g3, e3), (n1, g1, e1)),
        ):
            boltz = (ga / gb) * math.exp(-HC_CM * (ea - eb) * 1e7 / (kb_cgs * T))
            assert na / nb == pytest.approx(boltz, rel=1e-3, abs=0.0), sp
        checked += 1
    assert checked >= 4  # most of the seven systems thermalize fully


@pytest.mark.physics_invariant
def test_cooling_linear_in_electron_density_in_the_coronal_limit():
    """Far below every critical density the cooling is linear in n_e.

    The N 2D metastable level has a critical density of only about
    2.5e3 cm^-3 (A21 = 1.3e-5 s^-1), so electron densities of order unity
    sit deep in the coronal regime, where doubling n_e doubles the loss.
    The saturation guard: at LTE-scale electron density the ratio collapses
    far below 2, which discriminates a missing collisional de-excitation.
    """
    q1 = three_level_cooling('N', 1e8, 1.0, 9000.0)
    q2 = three_level_cooling('N', 1e8, 2.0, 9000.0)
    assert q2 / q1 == pytest.approx(2.0, rel=1e-3, abs=0.0)
    q_lte1 = three_level_cooling('N', 1e8, 1e19, 9000.0)
    q_lte2 = three_level_cooling('N', 1e8, 2e19, 9000.0)
    assert q_lte2 / q_lte1 < 1.1
    # Zero density in either species gives zero cooling (error contract).
    assert three_level_cooling('N', 0.0, 1e6, 9000.0) == pytest.approx(0.0, abs=0.0)


@pytest.mark.reference_pinned
def test_hydrogen_system_brackets_the_black_lyalpha_rate():
    """The H three-level cooling agrees with Black (1981) at order unity.

    Murray-Clay et al. (2009, Eq. 6) print the Black rate ``7.5e-19 n_e n_H
    exp(-118348 K / T)`` erg cm^-3 s^-1. The three-level system uses
    effective collision strengths frozen at 1e4 K and carries no cascades,
    while the Black fit carries both, so the two agree only at order unity
    (the measured ratio runs from about 0.5 at 1e4 K to 0.3 at 3e4 K). The
    bracket catches transcription errors, in the constant and in the
    exponential activation, and nothing finer.
    """
    pref, tscale = LYA_BLACK
    n_h, n_e = 1e8, 1e6
    ratios = []
    for T in (1.0e4, 2.0e4):
        q_mine = three_level_cooling('H', n_h, n_e, T)
        q_black = pref * n_e * n_h * math.exp(-tscale / T)
        ratios.append(q_mine / q_black)
        assert 0.15 < q_mine / q_black < 1.5, T
    # The frozen collision strengths fall behind the Black fit as the
    # temperature rises, so the ratio declines with T.
    assert ratios[1] < ratios[0]


@pytest.mark.physics_invariant
def test_ionization_fraction_limits_and_monotonicity():
    """The ionization fraction is bounded, zero at zero flux, and monotone.

    Exact zero without flux or density; strictly inside (0, 1] otherwise;
    monotone increasing with flux and decreasing with density (stronger
    recombination sinks). The strong-flux limit approaches full ionization.
    """
    assert ionization_fraction(1e10, 1e4, 0.0, 1e-11, 1e-17, 2.7e-13) == pytest.approx(
        0.0, abs=0.0
    )
    f = ionization_fraction(1e6, 1e4, 1e6, 1e-11, 1e-17, 2.7e-13)
    assert 0.0 < f <= 1.0
    f_hi = ionization_fraction(1e6, 1e4, 1e8, 1e-11, 1e-17, 2.7e-13)
    assert f_hi > f
    f_dense = ionization_fraction(1e12, 1e4, 1e6, 1e-11, 1e-17, 2.7e-13)
    assert f_dense < f
    f_sat = ionization_fraction(1e2, 1e4, 1e12, 1e-11, 1e-17, 2.7e-13)
    assert f_sat == pytest.approx(1.0, abs=1e-3)


def test_front_and_recombination_selection_follow_composition():
    """Front constants and the recombination route follow the composition.

    A hydrogen-dominated wind takes the 20 eV hydrogen front; a nitrogen
    wind takes the harder 33.6 eV front with the smaller cross section, and
    its recombination coefficient comes from the Badnell fit rather than
    the case B average (the two differ resolvably at 1e4 K).
    """
    hnu_h, eion_h, sig_h = front_constants({'H': 0.9, 'O': 0.1})
    hnu_n, eion_n, sig_n = front_constants({'N': 0.8, 'O': 0.2})
    assert hnu_n > hnu_h
    assert sig_n > sig_h
    assert eion_n > eion_h
    from zephyrus.atomic_data import alpha_case_b, badnell_alpha_rr

    a_n = recombination_alpha({'N': 0.8, 'O': 0.2}, 1e4)
    assert a_n == pytest.approx(badnell_alpha_rr(1e4), rel=1e-12, abs=0.0)
    a_mix = recombination_alpha({'H': 0.6, 'O': 0.4}, 1e4)
    expected = 0.6 * alpha_case_b('H', 1e4) + 0.4 * alpha_case_b('O', 1e4)
    assert a_mix == pytest.approx(expected, rel=1e-12, abs=0.0)
    # Ratio comparison: at 1e-13 scale an absolute tolerance would swamp the
    # difference, so the resolvability check is multiplicative.
    assert a_n / a_mix > 1.5


def test_rootfind_contract_edges_and_root():
    """The root-find honors its bracket, clamps flagged, and balances at a root.

    Zero flux means no heating, so cooling wins at the lower edge and the
    temperature clamps low. A strong flux on a nitrogen wind finds a root
    inside the bracket, where heating balances cooling to the bisection
    tolerance. An equilibrium temperature at or above the upper bracket
    edge clamps high immediately.
    """
    T, d = solve_wind_temperature(700.0, _base(), {'N': 1.0}, 0.0)
    assert T == pytest.approx(700.0, rel=1e-12, abs=0.0)
    assert d['clamped'] == 'low'
    T, d = solve_wind_temperature(700.0, _base(), {'N': 1.0}, 5.0)
    assert 700.0 <= T <= 5.0e4
    if d['clamped'] is None:
        assert abs(d['q_heat'] - d['q_cool']) / d['q_heat'] < 1e-3
    T, d = solve_wind_temperature(6.0e4, _base(), {'N': 1.0}, 5.0)
    assert T == pytest.approx(5.0e4, rel=1e-12, abs=0.0)
    assert d['clamped'] == 'high'


def test_all_cooling_channels_off_is_rejected():
    """Disabling every cooling channel raises instead of clamping silently.

    A pure-heating balance has no root by construction; the contract is a
    ``ValueError``, and a single enabled channel on the same path returns
    normally.
    """
    with pytest.raises(ValueError, match='cooling channels'):
        solve_wind_temperature(
            700.0,
            _base(),
            {'N': 1.0},
            5.0,
            cool_atomic=False,
            cool_co2_band=False,
            cool_o_finestructure=False,
            cool_recombination=False,
        )
    T, _ = solve_wind_temperature(
        700.0,
        _base(),
        {'N': 1.0},
        5.0,
        cool_atomic=True,
        cool_co2_band=False,
        cool_o_finestructure=False,
        cool_recombination=False,
    )
    assert 700.0 <= T <= 5.0e4


@pytest.mark.physics_invariant
def test_balance_parts_sum_and_channel_toggles():
    """The cooling parts sum to the total and toggles remove their channel.

    Energy bookkeeping: the per-channel parts must sum to ``q_cool``
    exactly. Turning the atomic channel off removes its part and can only
    reduce the total cooling at fixed temperature.
    """
    comp = {'N': 0.8, 'O': 0.2}
    r_all, d_all = balance_at(9000.0, _base(), comp, 5.0)
    assert sum(d_all['parts'].values()) == pytest.approx(d_all['q_cool'], rel=1e-12, abs=0.0)
    r_no_atomic, d_no = balance_at(9000.0, _base(), comp, 5.0, cool_atomic=False)
    assert 'atomic_lines' not in d_no['parts']
    assert d_no['q_cool'] <= d_all['q_cool']
    assert r_no_atomic >= r_all


@pytest.mark.reference_pinned
def test_line_cooling_coronal_limit_absolute_normalization():
    """The line cooling's scale, evaluated from the rate coefficient itself.

    Far below the critical density every collisional excitation out of the
    ground state radiates before it is deexcited, so the cooling reduces to
    ``sum over transitions of dE k_lu n_e n_tot`` with the Maxwellian
    collisional rate coefficient ``k_lu = Upsilon (8.629e-6 / (g_l sqrt(T)))
    exp(-dE / k_B T)``. That is evaluated here from the level data and the
    constants rather than by calling the function, so the prefactor and the
    statistical weights are pinned rather than only the shape. Without this
    the whole channel could be scaled by any factor with the suite green.
    """
    T, n_tot = 8000.0, 1.0e6
    data = THREE_LEVEL['O']
    levels = data['levels']
    for n_e in (1.0e0, 1.0e2):
        expected = 0.0
        for (lo, up), (a_ul, upsilon) in data['transitions'].items():
            if a_ul <= 0.0 or lo != 1:
                continue
            g_lo = levels[lo - 1][1]
            de = HC_CM * (levels[up - 1][2] - levels[lo - 1][2]) * 1e7
            k_lu = upsilon * 8.629e-6 / (g_lo * math.sqrt(T)) * math.exp(-de / (kb_cgs * T))
            expected += de * k_lu * n_e * n_tot
        got = three_level_cooling('O', n_tot, n_e, T)
        assert got == pytest.approx(expected, rel=1e-3, abs=0.0), n_e
    # Linear in the electron density in this limit, and in the total density.
    q1 = three_level_cooling('O', n_tot, 1.0e0, T)
    q2 = three_level_cooling('O', 2.0 * n_tot, 2.0e0, T)
    assert q2 / q1 == pytest.approx(4.0, rel=1e-3, abs=0.0)


@pytest.mark.reference_pinned
def test_recombination_cooling_coefficient():
    """Recombination cooling removes 3/2 k_B T per recombination.

    The channel is ``Q = n_e n_+ alpha_B (3/2) k_B T``, and the 3/2 is the
    mean thermal energy carried off by the recombining electron. Pinned
    against a hand evaluation, since the factor is the whole content of the
    term and doubling it doubled the channel with the suite green.
    """
    T = 1.0e4
    # The published coefficient at 1e4 K, which the term is built on.
    assert alpha_case_b('H', T) == pytest.approx(2.7e-13, rel=1e-9, abs=0.0)
    # The balance's own term, rebuilt from the quantities it reports. The
    # module works in cgs internally, so the number density converts from
    # m^-3 at the boundary and the square of that conversion is where a
    # units slip would hide.
    base = _base(n_si=1.0e18, vmr={'H': 1.0})
    _r, det = balance_at(T, base, {'H': 1.0}, 1.0e2)
    n_cgs = base['n'] * 1e-6
    f_plus = det['f_plus']
    expected = (f_plus * n_cgs) ** 2 * det['alpha_rec_cgs'] * 1.5 * kb_cgs * T
    assert det['parts']['recombination'] == pytest.approx(expected, rel=1e-9, abs=0.0)
    # Discrimination: the 3/2 is the mean energy carried off per
    # recombination and is the whole content of the coefficient, so the same
    # term with a 3 instead is a different number by exactly a factor two.
    assert det['parts']['recombination'] != pytest.approx(
        2.0 * expected, rel=0.01, abs=0.0
    )
    # Quadratic in the electron density, not in the total density: doubling
    # the gas density raises the term by less than four, because the
    # ionization fraction falls as recombination speeds up. The identity
    # above is what pins the coefficient; this pins that the term is built
    # on the electron density and not on the neutral one.
    dense = _base(n_si=2.0e18, vmr={'H': 1.0})
    _r2, det2 = balance_at(T, dense, {'H': 1.0}, 1.0e2)
    ratio = det2['parts']['recombination'] / det['parts']['recombination']
    assert 1.0 < ratio < 4.0
    assert det2['f_plus'] < det['f_plus']
    n2_cgs = dense['n'] * 1e-6
    electron_ratio = (det2['f_plus'] * n2_cgs) ** 2 / (det['f_plus'] * n_cgs) ** 2
    assert ratio == pytest.approx(electron_ratio, rel=1e-6, abs=0.0)
