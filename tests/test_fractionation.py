"""Tests for ``src/zephyrus/fractionation.py``: exact reductions and anchors.

The closure must reproduce, exactly or at the published values, every
special case of the escape-fractionation lineage it generalizes:

- The three-species deuterium system of Gu & Chen (2023): the escaping and
  retained-helium branches, both critical rates, and the helium admixture
  factor on the deuterium threshold (also the ternary limit of Cherubim
  et al. 2024, their Eq. 11).
- The trace-minor relations of Odert et al. (2018, Eq. 5) and Zahnle et
  al. (1990, Eqs. 35, 36, 42), including the adjudication that the earlier
  Zahnle & Kasting (1986) Eq. (36) drag-deficit weighting is NOT
  reproduced.
- The non-trace three-species relations of Zahnle & Kasting (2023,
  Eqs. 19-20).
- The prescribed-flux partition of Chassefiere (1996, Eqs. 1, 6, 7).
- The universal-b closed form and the Hunten et al. (1987) Earth, Mars,
  and Venus numerical anchors.
- Low-flux collapse onto the lightest species and the zero-flux limit.
- The SI shim: per-element rates conserve the bulk rate at machine
  precision, and rock-forming species raise their provenance flag.

The randomized ensemble sweeps (global properties across thresholds,
brute-force active-set uniqueness, fixed-point stability) live in the
smoke-tier companion ``tests/test_fractionation_ensembles.py``.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import numpy as np
import pytest

from zephyrus.constants import kb_cgs
from zephyrus.fractionation import (
    closure_per_species,
    first_threshold,
    solve_closure,
    unfractionated_split,
)
from zephyrus.planets_parameters import Me, Re

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

AMU_G = 1.66053907e-24  # g


def test_input_validation_error_contract():
    """Malformed solver inputs raise; a valid call on the same path returns.

    Negative flux, mole fractions off unit sum, non-positive masses, and an
    asymmetric coefficient matrix are not physically posed inputs and must
    raise ``ValueError`` rather than return a partial solution.
    """
    m = np.array([1.0, 16.0]) * AMU_G
    X = np.array([0.8, 0.2])
    b = np.array([[np.inf, 1e19], [1e19, np.inf]])
    with pytest.raises(ValueError, match='phi'):
        solve_closure(-1.0, X, m, 400.0, 980.0, b)
    with pytest.raises(ValueError, match='sum to 1'):
        solve_closure(1e-10, np.array([0.8, 0.4]), m, 400.0, 980.0, b)
    with pytest.raises(ValueError, match='positive'):
        solve_closure(1e-10, X, -m, 400.0, 980.0, b)
    b_asym = np.array([[np.inf, 1e19], [2e19, np.inf]])
    with pytest.raises(ValueError, match='symmetric'):
        solve_closure(1e-10, X, m, 400.0, 980.0, b_asym)
    flux = solve_closure(1e-10, X, m, 400.0, 980.0, b)
    assert np.all(flux >= 0.0)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_ternary_deuterium_reductions():
    """The H-He-D system reproduces the Gu & Chen (2023) relations exactly.

    With helium escaping, the trace-deuterium flux matches their Eq. (4)
    (equivalently Cherubim et al. 2024, Eq. 11); with helium
    retained, their Eq. (8); the activation thresholds match their two
    critical rates including the (1 + alpha_2 X_He/X_H)^-1 helium factor
    on the deuterium threshold, which must lower it relative to the
    helium-free case (the discrimination guard on the factor's sign).
    """
    m = np.array([1.0, 4.0, 2.0]) * AMU_G  # H, He, D (trace)
    T, g0 = 1000.0, 1000.0
    kT = kb_cgs * T
    b12 = 1.04e18 * T**0.732  # H-He (Zahnle & Kasting 1986 lineage)
    b13 = 7.183e17 * T**0.728  # H-D (Genda & Ikoma 2008)
    b23 = 5.087e17 * T**0.728  # He-D (scaled)
    b = np.array([[np.inf, b12, b13], [b12, np.inf, b23], [b13, b23, np.inf]])
    x3 = 1e-13  # the printed formulas are O(X_trace) truncations
    x2 = 0.15
    x1 = 1 - x2 - x3
    X = np.array([x1, x2, x3])
    a2, a3 = b13 / b12, b13 / b23
    phi_dl_he = b12 * (m[1] - m[0]) * g0 / kT
    phi_dl_d = b13 * (m[2] - m[0]) * g0 / kT
    phi_crit_he = m[0] * x1 * phi_dl_he
    phi_crit_d = m[0] * phi_dl_d / (1 + a2 * x2 / x1)

    # (a) Supercritical: He escaping, trace D follows their Eq. (4).
    for frac in (1.5, 5.0, 50.0):
        flux = solve_closure(frac * phi_crit_he, X, m, T, g0, b)
        f2, f3 = x2 / x1, x3 / x1
        ref = f3 * (flux[0] + a3 * flux[1] + a2 * phi_dl_he * x2 - phi_dl_d) / (1 + a3 * f2)
        assert flux[2] == pytest.approx(ref, rel=1e-11, abs=0.0), frac

    # (b) Subcritical He (retained), D escaping: their Eq. (8).
    for frac in (1.5, 3.0):
        phi = frac * phi_crit_d
        if phi >= phi_crit_he:
            continue
        flux = solve_closure(phi, X, m, T, g0, b)
        assert flux[1] == pytest.approx(0.0, abs=0.0)
        ref = x3 * (flux[0] * (1 + a2 * x2 / x1) - phi_dl_d) / (x1 + a3 * x2)
        assert flux[2] == pytest.approx(ref, rel=1e-11, abs=0.0), frac

    # (c) Both activation thresholds are sharp at the printed critical rates.
    eps = 1e-6
    assert solve_closure(phi_crit_he * (1 - eps), X, m, T, g0, b)[1] == pytest.approx(
        0.0, abs=0.0
    )
    assert solve_closure(phi_crit_he * (1 + eps), X, m, T, g0, b)[1] > 0.0
    assert solve_closure(phi_crit_d * (1 - eps), X, m, T, g0, b)[2] == pytest.approx(
        0.0, abs=0.0
    )
    assert solve_closure(phi_crit_d * (1 + eps), X, m, T, g0, b)[2] > 0.0

    # (d) The He admixture lowers the D threshold by (1 + a2 X_He/X_H)^-1.
    tiny = 1e-14
    x_b = np.array([1 - x3 - tiny, tiny, x3])
    phi_crit_d_no_he = m[0] * phi_dl_d / (1 + a2 * tiny / x_b[0])
    assert solve_closure(phi_crit_d_no_he * (1 + eps), x_b, m, T, g0, b)[2] > 0.0
    assert solve_closure(phi_crit_d_no_he * (1 - eps), x_b, m, T, g0, b)[2] == pytest.approx(
        0.0, abs=0.0
    )
    assert phi_crit_d < phi_crit_d_no_he


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_two_majors_trace_minor_relations():
    """Two escaping majors with trace minors reproduce the printed relations.

    With H and O both escaping, each entrained trace minor obeys Odert et
    al. (2018, Eq. 5), which equals Zahnle et al. (1990, Eq. 35); where the
    printed formula goes negative the closure returns exactly zero (the
    clamp those authors themselves prescribe); at the O limiting flux the
    light minor follows Zahnle et al. (1990, Eq. 36). The adjudication
    guard: the earlier Zahnle & Kasting (1986, Eq. 36) drag-deficit
    weighting must NOT be reproduced (it deviates by more than 5 percent
    somewhere in the sweep), which discriminates the two printed variants.
    """
    m = np.array([1.0, 16.0, 12.0, 36.0]) * AMU_G  # H, O majors; C, Ar traces
    T, g0 = 400.0, 980.0
    kT = kb_cgs * T
    b12 = 4.8e17 * T**0.75
    b_hc, b_har = 8.0e17 * T**0.70, 1.06e18 * T**0.597
    b_oc, b_oar = 9.0e16 * T**0.80, 5.61e16 * T**0.841
    b_car = 8.0e16 * T**0.78
    b = np.array(
        [
            [np.inf, b12, b_hc, b_har],
            [b12, np.inf, b_oc, b_oar],
            [b_hc, b_oc, np.inf, b_car],
            [b_har, b_oar, b_car, np.inf],
        ]
    )
    x_tr = 1e-13
    n_cmp = 0
    zk_dev = 0.0
    for f2 in (0.1, 0.5, 1.0, 2.0):
        x2 = f2 / (1 + f2) * (1 - 2 * x_tr)
        x1 = 1 - x2 - 2 * x_tr
        X = np.array([x1, x2, x_tr, x_tr])
        philim = m[0] * x1 * b12 * (m[1] - m[0]) * g0 / kT / (x1 + x2)
        for frac in (1.3, 3.0, 10.0, 100.0):
            flux = solve_closure(frac * philim, X, m, T, g0, b)
            w1, w2 = flux[0] / x1, flux[1] / x2
            x2_rel = w2 / w1
            f1_flux = flux[0]
            for k, b1k, b2k in ((2, b_hc, b_oc), (3, b_har, b_oar)):
                xk = (flux[k] / x_tr) / w1
                xk_ref = (
                    1
                    - g0 * (m[k] - m[0]) * b1k / (f1_flux * kT)
                    + (b1k / b12) * f2 * (1 - x2_rel)
                    + (b1k / b2k) * f2 * x2_rel
                ) / (1 + (b1k / b2k) * f2)
                if xk_ref > 1e-6:
                    n_cmp += 1
                    assert xk == pytest.approx(xk_ref, rel=1e-11, abs=0.0), (f2, frac, k)
                    mu2 = m[1] / m[0]
                    xk_zk = (
                        1
                        - g0 * (m[k] - m[0]) * b1k / (f1_flux * kT)
                        + (b1k / b12) * f2 * (1 + f2) * (1 - x2_rel) / (mu2 + f2)
                        + (b1k / b2k) * f2 * x2_rel
                    ) / (1 + (b1k / b2k) * f2)
                    if x2_rel < 0.99:
                        zk_dev = max(zk_dev, abs(xk_zk - xk) / xk)
                elif xk_ref < -1e-6:
                    assert flux[k] == pytest.approx(0.0, abs=0.0), (f2, frac, k)
        # At the limiting flux (O marginally retained) the light minor
        # follows Zahnle et al. (1990, Eq. 36).
        flux = solve_closure(philim * (1 - 1e-9), X, m, T, g0, b)
        x3 = (flux[2] / x_tr) / (flux[0] / x1)
        x3_z90 = (
            1
            - (m[2] - m[0]) / (m[1] - m[0]) * (b_hc / b12)
            + (m[1] - m[2]) / (m[1] - m[0]) * (b_hc / b12) * f2
        ) / (1 + (b_hc / b_oc) * f2)
        assert x3 == pytest.approx(x3_z90, rel=2e-8, abs=0.0), f2
    assert n_cmp >= 8  # enough entrained comparisons to be meaningful
    assert zk_dev > 0.05  # the 1986 variant is resolvably not reproduced


@pytest.mark.reference_pinned
def test_first_entrainment_with_two_retained_heavies():
    """The first-entrainment threshold matches Zahnle et al. (1990, Eq. 42).

    H2 escaping through retained CO2 and N2: their Eq. (42) gives the
    per-background thresholds; whichever is smaller names the species that
    entrains first, sharply, while the other stays retained on both sides.
    """
    m = np.array([2.0, 44.0, 28.0]) * AMU_G
    T, g0 = 400.0, 373.0
    kT = kb_cgs * T
    b12 = 2.3e17 * T**0.75
    b13 = 2.65e17 * T**0.75
    b23 = 1e17 * T**0.75
    b = np.array([[np.inf, b12, b13], [b12, np.inf, b23], [b13, b23, np.inf]])
    f2, f3 = 1.0, 0.5
    x1 = 1 / (1 + f2 + f3)
    X = np.array([x1, f2 * x1, f3 * x1])
    phi12 = b12 * (m[1] - m[0]) * g0 / kT / (1 + f2 + (b12 / b13) * f3)
    phi13 = b13 * (m[2] - m[0]) * g0 / kT / (1 + f3 + (b13 / b12) * f2)
    phi_thresh = m[0] * min(phi12, phi13)
    eps = 1e-6
    flux_lo = solve_closure(phi_thresh * (1 - eps), X, m, T, g0, b)
    flux_hi = solve_closure(phi_thresh * (1 + eps), X, m, T, g0, b)
    idx = 1 if phi12 < phi13 else 2
    assert flux_lo[idx] == pytest.approx(0.0, abs=0.0)
    assert flux_hi[idx] > 0.0
    assert flux_lo[3 - idx] == pytest.approx(0.0, abs=0.0)
    # The library's own first_threshold agrees with the printed expression.
    assert first_threshold(X, m, T, g0, b) == pytest.approx(phi_thresh, rel=1e-9, abs=0.0)


@pytest.mark.reference_pinned
def test_zk23_nontrace_ternary_relations():
    """Non-trace H, O, CO2 reproduce Zahnle & Kasting (2023, Eqs. 19-20).

    Their Eq. (19): the oxygen-crossover flux of hydrogen escaping alone
    over a retained CO2 background, checked in closed form to machine
    precision and through the solver's own bisected activation threshold
    (whose residual floor is the solver's active-set tolerance, bounded at
    1e-9). Their Eq. (20): with H and O both escaping and CO2 static, the
    printed difference relation holds at every flux inside the two-species
    band, not only at an endpoint.
    """
    m = np.array([1.0, 16.0, 44.0]) * AMU_G
    T, g0 = 1000.0, 870.0
    kT = kb_cgs * T
    b12 = 4.8e17 * T**0.75
    b14 = 6.0e19 * (T / 1000.0) ** 0.75  # their Table 2 H-CO2 row
    b24 = 5.0e16 * T**0.75  # plausible heavy-heavy value; Eqs. 19-20 are identities in it
    b = np.array([[np.inf, b12, b14], [b12, np.inf, b24], [b14, b24, np.inf]])

    for f4 in (0.1, 0.5):
        for r21 in (0.5, 2.0):
            f2 = (1 - f4) * r21 / (1 + r21)
            f1 = 1 - f4 - f2
            X = np.array([f1, f2, f4])
            zk19 = (m[1] - m[0]) * g0 * b12 / kT / (1 + f4 * (b12 / b14 - 1))
            phi_star = f1 * (m[1] - m[0]) * g0 / kT / ((f1 + f2) / b12 + f4 / b14)
            assert phi_star / f1 == pytest.approx(zk19, rel=1e-13, abs=0.0), (f4, r21)
            ph1 = first_threshold(X, m, T, g0, b)
            th = _activation_threshold(1, X, m, T, g0, b, ph1 * 1e-4, ph1 * 1e3)
            assert (th / m[0]) / f1 == pytest.approx(zk19, rel=1e-9, abs=0.0), (f4, r21)
            eps = 1e-6
            _, _, lo_act = solve_closure(th * (1 - eps), X, m, T, g0, b, return_diag=True)
            _, _, hi_act = solve_closure(th * (1 + eps), X, m, T, g0, b, return_diag=True)
            assert lo_act == frozenset({0})
            assert hi_act == frozenset({0, 1})  # O entrains here, not CO2

    for f4 in (0.2, 0.6):
        f2 = (1 - f4) / 2
        f1 = 1 - f4 - f2
        X = np.array([f1, f2, f4])
        ph1 = first_threshold(X, m, T, g0, b)
        th_o = _activation_threshold(1, X, m, T, g0, b, ph1 * 1e-4, ph1 * 1e3)
        th_c = _activation_threshold(2, X, m, T, g0, b, ph1 * 1e-4, ph1 * 1e6)
        for frac in (0.2, 0.5, 0.9):
            phi = th_o * (th_c / th_o) ** frac
            flux, _c, act = solve_closure(phi, X, m, T, g0, b, return_diag=True)
            assert act == frozenset({0, 1}), (f4, frac)
            lhs = flux[0] * (1 + f2 / f1 + (f4 / f1) * (b12 / b14)) - flux[1] * (
                1 + f1 / f2 + (f4 / f2) * (b12 / b24)
            )
            rhs = g0 * (m[1] - m[0]) * b12 / kT
            assert lhs == pytest.approx(rhs, rel=1e-12, abs=0.0), (f4, frac)


def _activation_threshold(k, X, m, T, g0, b, lo, hi, niter=100):
    """Mass flux at which the solver first admits species k, by bisection."""
    for _ in range(niter):
        mid = np.sqrt(lo * hi)
        if k in solve_closure(mid, X, m, T, g0, b, return_diag=True)[2]:
            hi = mid
        else:
            lo = mid
    return np.sqrt(lo * hi)


@pytest.mark.reference_pinned
def test_chassefiere_prescribed_flux_partition():
    """The binary partition reproduces Chassefiere (1996, Eqs. 1, 6, 7).

    His crossover mass ``m_c = m_1 + kT F_1 / (b g X_1)`` evaluated at the
    threshold flux equals the heavy mass exactly (his dropout test and the
    closure's threshold are the same statement), and above threshold the
    solved split satisfies his Eqs. (1) and (6) to machine precision at
    every draw. His small-mass approximation Eq. (12) is deliberately not
    asserted: it deviates by up to a factor of a few, as its own stated
    approximation requires.
    """
    rng = np.random.default_rng(112)
    for _ in range(60):
        m1 = AMU_G * 1.0
        m2 = m1 * rng.uniform(2, 40)
        r21 = rng.choice([0.2, 0.5, 1.0, 2.0])
        x2 = r21 / (1 + r21)
        x1 = 1 - x2
        T = rng.uniform(300, 2000)
        g0 = rng.uniform(300, 1500)
        b12 = 4.8e17 * T**0.75
        b = np.array([[np.inf, b12], [b12, np.inf]])
        kT = kb_cgs * T
        X, mm = np.array([x1, x2]), np.array([m1, m2])
        phi_c = b12 * x1 * (m2 - m1) * m1 * g0 / kT
        f1_c = solve_closure(phi_c, X, mm, T, g0, b)[0]
        assert m1 + kT * f1_c / (b12 * g0 * x1) == pytest.approx(m2, rel=1e-12, abs=0.0)
        for frac in (1.2, 2.0, 10.0, 100.0):
            f1, f2 = solve_closure(frac * phi_c, X, mm, T, g0, b)
            f1_ref = frac * phi_c / m1  # his Eq. (5)
            mc = m1 + kT * f1 / (b12 * g0 * x1)  # his Eq. (7)
            pred6 = 1.0 / (1 + (x2 / x1) * (m2 / m1) * (mc - m2) / (mc - m1))
            assert f1 / f1_ref == pytest.approx(pred6, rel=1e-12, abs=0.0)
            pred1 = (x2 / x1) * f1 * (mc - m2) / (mc - m1)  # his Eq. (1)
            assert f2 == pytest.approx(pred1, rel=1e-12, abs=0.0)


@pytest.mark.physics_invariant
def test_universal_b_closed_form():
    """With one common coefficient the fluxes take the analytic closed form.

    When every pair shares one b the fully entrained solution is
    ``Phi_j = X_j [Phi_tot - (b g0 / kT)(m_j - m_bar)]``: the flux excess
    over the mean follows the mass deviation linearly. Exact to 1e-10,
    which any indexing or drag-bookkeeping error breaks.
    """
    rng = np.random.default_rng(105)
    n = 6
    m = np.sort(rng.uniform(1, 50, n)) * AMU_G
    X = rng.dirichlet(np.ones(n))
    T, g0 = 800.0, 1500.0
    kT = kb_cgs * T
    bval = 2e17 * T**0.75
    b = np.full((n, n), bval)
    np.fill_diagonal(b, np.inf)
    mbar = np.sum(m * X)
    var = np.sum(X * (m - mbar) ** 2)
    beta = bval * g0 / kT
    w_needed = beta * (np.max(m) - mbar) * 2
    phi = mbar * w_needed * 4
    flux = solve_closure(phi, X, m, T, g0, b)
    phi_tot = (phi + beta * var) / mbar
    ref = X * (phi_tot - beta * (m - mbar))
    np.testing.assert_allclose(flux, ref, rtol=0, atol=1e-10 * np.max(ref))
    # Conservation on top of the closed form.
    assert np.sum(m * flux) == pytest.approx(phi, rel=1e-12, abs=0.0)


@pytest.mark.reference_pinned
def test_hunten_anchors_earth_mars_venus():
    """The printed Hunten et al. (1987) worked anchors are reproduced.

    Their Earth (crossover mass 140 amu at reference flux 8.1e13
    cm^-2 s^-1), Mars (130 amu, 2.9e13), and Venus (2e11 giving crossover
    mass 1.35 amu) numbers follow from the threshold relation with their
    stated b, g, and T = 400 K; the solver's own activation threshold on a
    trace heavy over pure light gas reproduces each within 1 percent.
    """
    kt400 = kb_cgs * 400.0
    f1_earth = (140 - 1) * AMU_G * 2e19 * 980 / kt400
    assert f1_earth == pytest.approx(8.1e13, rel=0.02, abs=0.0)
    f1_mars = (130 - 1) * AMU_G * 2e19 * 373 / kt400
    assert f1_mars == pytest.approx(2.9e13, rel=0.02, abs=0.0)
    mc_venus = 1 + kt400 * 2e11 / (2.2e19 * 850 * AMU_G)
    assert mc_venus == pytest.approx(1.35, abs=0.02)
    for m2_amu, g0, b12, f1_ref in (
        (140.0, 980.0, 2e19, 8.19e13),
        (130.0, 373.0, 2e19, 2.894e13),
        (mc_venus, 850.0, 2.2e19, 2e11),
    ):
        m = np.array([1.0, m2_amu]) * AMU_G
        X = np.array([1 - 1e-10, 1e-10])
        b = np.array([[np.inf, b12], [b12, np.inf]])
        phi_star = m[0] * X[0] * b12 * (m[1] - m[0]) * g0 / kt400 * 400.0 / 400.0
        eps = 1e-6
        assert solve_closure(phi_star * (1 - eps), X, m, 400.0, g0, b)[1] == pytest.approx(
            0.0, abs=0.0
        )
        assert solve_closure(phi_star * (1 + eps), X, m, 400.0, g0, b)[1] > 0.0
        assert phi_star / m[0] == pytest.approx(f1_ref, rel=0.01, abs=0.0)


@pytest.mark.physics_invariant
def test_low_flux_collapse_and_zero_limit():
    """Low flux collapses onto the lightest species; zero flux returns zero.

    Down to 1e-12 of the first threshold only the lightest species is
    active and carries the whole mass flux exactly; at ``phi = 0`` the
    fluxes are zero with an empty active set and no exception (the
    continuity limit of the multiplier is defined there).
    """
    rng = np.random.default_rng(109)
    for _ in range(4):
        n = int(rng.integers(2, 8))
        m = np.sort(rng.uniform(1, 60, n)) * AMU_G
        X = rng.dirichlet(np.ones(n))
        T = rng.uniform(300, 2000)
        g0 = rng.uniform(200, 3000)
        logb = rng.uniform(17.0, 20.0, (n, n))
        b = 10.0 ** (0.5 * (logb + logb.T)) * T**0.75
        np.fill_diagonal(b, np.inf)
        phi1 = first_threshold(X, m, T, g0, b)
        light = int(np.argmin(m))
        for frac in (1e-3, 1e-6, 1e-12):
            phi = frac * phi1
            flux, _c, act = solve_closure(phi, X, m, T, g0, b, return_diag=True)
            assert act == frozenset({light})
            assert flux[light] == pytest.approx(phi / m[light], rel=1e-12, abs=0.0)
            assert np.sum(m * flux) == pytest.approx(phi, rel=1e-12, abs=0.0)
        flux0, _c0, act0 = solve_closure(0.0, X, m, T, g0, b, return_diag=True)
        np.testing.assert_array_equal(flux0, 0.0)
        assert act0 == frozenset()


@pytest.mark.physics_invariant
def test_per_species_shim_conserves_mass():
    """The SI shim returns non-negative per-element rates summing to the bulk.

    A mixed H-He-O wind at a super-Earth base: the per-element rates must
    sum to the bulk rate at machine precision, all be non-negative, and
    keep hydrogen in the active set (the lightest species always escapes
    when anything does).
    """
    per, diag, _flags = closure_per_species(
        1.0e6, {'H': 0.85, 'He': 0.10, 'O': 0.05}, 8000.0, 5 * Me, 2 * Re
    )
    assert sum(per.values()) == pytest.approx(1.0e6, rel=1e-9, abs=0.0)
    assert diag['mass_conservation_rel'] < 1e-9
    assert all(v >= 0.0 for v in per.values())
    assert 'H' in diag['active_set']
    # Fractionation direction: hydrogen escapes preferentially relative to
    # its base mass fraction when heavies are near their thresholds.
    m_frac_h = 0.85 * 1.008 / (0.85 * 1.008 + 0.10 * 4.0026 + 0.05 * 15.999)
    assert per['H'] / 1.0e6 >= m_frac_h - 1e-9


def test_shim_low_flux_and_rock_former_flag():
    """The shim collapses onto hydrogen at low flux and flags rock formers.

    A tiny bulk rate leaves every heavy retained: oxygen carries exactly
    zero and hydrogen exactly the bulk rate. A silicon-bearing wind raises
    the rock-former provenance flag naming the species.
    """
    per, diag, _ = closure_per_species(1.0e-3, {'H': 0.5, 'O': 0.5}, 8000.0, 5 * Me, 2 * Re)
    assert per['O'] == pytest.approx(0.0, abs=0.0)
    assert per['H'] == pytest.approx(1.0e-3, rel=1e-9, abs=0.0)
    assert diag['retained'] == ['O']
    _per, _diag, flags = closure_per_species(
        1.0e6, {'H': 0.9, 'Si': 0.1}, 8000.0, 5 * Me, 2 * Re
    )
    assert flags.get('rock_former_bij') == ['Si']


def test_unfractionated_split_protocol():
    """The unfractionated split follows reservoirs, else base composition.

    With reservoir masses the split is by reservoir mass fraction exactly;
    without them it falls back to the atomized base composition's mass
    fractions with the substitution flagged, and either way the split sums
    to the bulk rate.
    """
    per, flags = unfractionated_split(4.0, {'H': 3.0e18, 'O': 1.0e18}, {'H': 1.0})
    assert per['H'] == pytest.approx(3.0, rel=1e-12, abs=0.0)
    assert per['O'] == pytest.approx(1.0, rel=1e-12, abs=0.0)
    assert flags == {}
    per2, flags2 = unfractionated_split(4.0, None, {'H': 0.5, 'O': 0.5})
    assert flags2.get('split_from_base_composition') is True
    assert sum(per2.values()) == pytest.approx(4.0, rel=1e-12, abs=0.0)
    # Mass weighting: oxygen outweighs hydrogen at equal mole fractions.
    assert per2['O'] > per2['H']
    # A reservoir that has run dry has no proportions to split by, which is
    # the state the end of an evolutionary track reaches. It takes the same
    # fallback as no reservoir at all, flagged, rather than dividing by zero.
    per3, flags3 = unfractionated_split(4.0, {'H': 0.0, 'O': 0.0}, {'H': 0.5, 'O': 0.5})
    assert flags3.get('split_from_base_composition') is True
    assert per3 == pytest.approx(per2, rel=1e-12, abs=0.0)
    assert sum(per3.values()) == pytest.approx(4.0, rel=1e-12, abs=0.0)
    # Negative reservoir masses are malformed input, not a dry reservoir.
    with pytest.raises(ValueError, match='non-negative'):
        unfractionated_split(4.0, {'H': -1.0, 'O': 2.0}, {'H': 0.5, 'O': 0.5})
