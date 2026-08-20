"""Ensemble sweeps for ``src/zephyrus/fractionation.py``.

Companion to ``tests/test_fractionation.py`` holding the randomized
ensemble verifications, which run for seconds rather than milliseconds and
so carry the smoke tier: they exercise the real solver across hundreds of
configurations per test. The properties under test:

- The two-species limit against the closed-form binary partition of
  Cherubim & Wordsworth (2024, Eqs. 7-9), including flux continuity at the
  crossover and exact mass conservation, over 200 random draws.
- Global properties across activation thresholds on random systems: mass
  conservation, non-negativity (including exactly at bisected thresholds),
  componentwise monotonicity, monotone active-set growth, two-sided
  continuity at every detected threshold, and piecewise linearity of the
  velocity scales within active-set segments.
- Active-set correctness and uniqueness by brute force: every candidate
  active set is enumerated, exactly one satisfies both Karush-Kuhn-Tucker
  conditions, and it matches the solver.
- Stability of the composition fixed point under well-mixed-layer
  relaxation, with the analytic universal-b Jacobian cross-checked against
  finite differences and general-b fixed points linearized numerically.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

from itertools import combinations

import numpy as np
import pytest

from zephyrus.constants import kb_cgs
from zephyrus.fractionation import first_threshold, solve_closure, solve_fixed_active

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

AMU_G = 1.66053907e-24  # g


def isofate_binary(phi, x1, x2, m1, m2, T, g0, b12):
    """Closed-form two-species partition (Cherubim & Wordsworth 2024, Eqs. 7-9)."""
    kT = kb_cgs * T
    h1 = kT / (m1 * g0)
    h2 = kT / (m2 * g0)
    phi_c = b12 * x1 * (m2 - m1) / h1
    if phi < phi_c:
        return phi / m1, 0.0, phi_c
    mbar = m1 * x1 + m2 * x2
    f1 = (x1 * phi + x1 * x2 * (m2 - m1) * b12 / h2) / mbar
    f2 = (x2 * phi + x1 * x2 * (m1 - m2) * b12 / h1) / mbar
    return f1, f2, phi_c


def solve_bruteforce(phi, X, m, T, g0, b, tol=1e-9):
    """Independent active-set selection by exhaustive enumeration (test oracle).

    Tries every nonempty candidate active set, keeps those satisfying both
    Karush-Kuhn-Tucker conditions (strictly positive drifts on the set, the
    retention condition non-positive off it), and returns the accepted
    ``(active_set, w, C)`` triples.
    """
    kT = kb_cgs * T
    n = len(X)
    accepted = []
    for r in range(1, n + 1):
        for cand in combinations(range(n), r):
            try:
                w, c_mult = solve_fixed_active(phi, X, m, T, g0, b, set(cand))
            except np.linalg.LinAlgError:
                continue
            wref = max(abs(w[j]) for j in cand)
            if any(w[j] <= tol * wref for j in cand):
                continue
            ok = True
            for k in range(n):
                if k in cand:
                    continue
                r_k = sum(X[i] * w[i] / b[i, k] for i in cand) - (m[k] * g0 / kT - c_mult)
                if r_k > tol * abs(m[k] * g0 / kT):
                    ok = False
                    break
            if ok:
                accepted.append((frozenset(cand), w, c_mult))
    return accepted


def _random_system(rng, nmin=3, nmax=12, iso_pair=False):
    """Random physical system: masses, composition, T, g, coefficient matrix."""
    n = int(rng.integers(nmin, nmax + 1))
    m = np.sort(rng.uniform(1, 60, n)) * AMU_G
    if iso_pair and n >= 3:
        m[1] = m[0] * (1 + rng.uniform(0.001, 0.05))  # isotope-close pair
        m = np.sort(m)
    X = rng.dirichlet(np.ones(n))
    T = rng.uniform(300, 2000)
    g0 = rng.uniform(200, 3000)
    logb = rng.uniform(17.0, 20.0, (n, n))  # three-decade coefficient spread
    b = 10.0 ** (0.5 * (logb + logb.T)) * T**0.75
    np.fill_diagonal(b, np.inf)
    return n, m, X, T, g0, b


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_binary_limit_matches_closed_form_over_random_draws():
    """200 random binaries match the closed-form partition on both branches.

    Below the crossover the light species carries everything; above it both
    escape with the printed split; the two branches join continuously at
    the crossover; fluxes are non-negative exactly at and beside the
    threshold; and mass is conserved on both branches. Any drag-term or
    threshold error breaks one of the five properties somewhere in the
    draw.
    """
    rng = np.random.default_rng(101)
    worst = 0.0
    for _ in range(200):
        m1 = AMU_G * rng.uniform(1, 4)
        m2 = m1 * rng.uniform(1.5, 40)
        x2 = rng.uniform(0.01, 0.9)
        x1 = 1 - x2
        T = rng.uniform(200, 2000)
        g0 = rng.uniform(100, 3000)
        b12 = rng.uniform(0.5, 5) * 1e17 * T**0.75
        b = np.array([[np.inf, b12], [b12, np.inf]])
        kT = kb_cgs * T
        phi_c = b12 * x1 * (m2 - m1) * m1 * g0 / kT
        X, mm = np.array([x1, x2]), np.array([m1, m2])
        for frac in (0.2, 0.999999, 1.000001, 1.7, 30.0):
            phi = frac * phi_c
            flux = solve_closure(phi, X, mm, T, g0, b)
            f1, f2, _ = isofate_binary(phi, x1, x2, m1, m2, T, g0, b12)
            err = max(abs(flux[0] - f1), abs(flux[1] - f2)) / max(f1, 1e-300)
            worst = max(worst, err)
        assert worst <= 1e-12
        flux_lo = solve_closure(phi_c * (1 - 1e-9), X, mm, T, g0, b)
        flux_hi = solve_closure(phi_c * (1 + 1e-9), X, mm, T, g0, b)
        jump = np.max(np.abs(flux_hi - flux_lo)) / (phi_c / m1)
        assert jump <= 1e-6
        for frac in (1.0, 1 - 1e-13, 1 + 1e-13):
            assert np.all(solve_closure(phi_c * frac, X, mm, T, g0, b) >= 0.0)
        for frac in (0.3, 3.0):
            phi = frac * phi_c
            flux = solve_closure(phi, X, mm, T, g0, b)
            assert m1 * flux[0] + m2 * flux[1] == pytest.approx(phi, rel=1e-12)


@pytest.mark.physics_invariant
def test_global_properties_across_thresholds():
    """Conservation, positivity, monotonicity, and continuity over random scans.

    Twenty random systems (up to twelve species, isotope-close pairs every
    third draw), each scanned over 400 fluxes spanning the full activation
    ladder: mass conservation and non-negativity at every point;
    componentwise monotonicity of the velocity scales in the flux; monotone
    growth of the active set; two-sided continuity at every bisected
    threshold (including non-negativity and conservation exactly there);
    and piecewise linearity of the velocity scales within each active-set
    segment (the solution is piecewise linear in the flux by construction).
    """
    rng = np.random.default_rng(107)
    n_trans = 0
    for trial in range(20):
        _n, m, X, T, g0, b = _random_system(rng, iso_pair=(trial % 3 == 0))
        kT = kb_cgs * T
        bmax = np.max(b[np.isfinite(b)])
        phimax = 3 * np.sum(m * X) * (bmax * g0 / kT) * (m[-1] - m[0])
        phis = np.geomspace(phimax * 1e-5, phimax, 400)
        prev = None
        seg = []
        for phi in phis:
            flux, _c, act = solve_closure(phi, X, m, T, g0, b, return_diag=True)
            w = np.where(X > 0, flux / np.where(X > 0, X, 1.0), 0.0)
            assert abs(np.sum(m * flux) - phi) / phi <= 1e-10
            assert np.all(flux >= 0.0)
            if prev is not None:
                pphi, pw, pact, pflux = prev
                assert not np.any(w - pw < -1e-9 * max(np.max(w), 1e-300))
                assert pact <= act  # the active set only grows with phi
                dref = max(np.max(np.abs(flux)), 1e-300)
                assert np.max(np.abs(flux - pflux)) / dref <= 0.15
                if act != pact:
                    n_trans += 1
                    lo, hi = pphi, phi
                    for _ in range(60):
                        mid = np.sqrt(lo * hi)
                        _, _, amid = solve_closure(mid, X, m, T, g0, b, return_diag=True)
                        if amid == pact:
                            lo = mid
                        else:
                            hi = mid
                    th = np.sqrt(lo * hi)
                    p_lo = solve_closure(th * (1 - 1e-9), X, m, T, g0, b)
                    p_hi = solve_closure(th * (1 + 1e-9), X, m, T, g0, b)
                    sref = max(np.max(np.abs(p_hi)), 1e-300)
                    assert np.max(np.abs(p_hi - p_lo)) / sref <= 1e-5
                    p_at = solve_closure(th, X, m, T, g0, b)
                    assert np.all(p_at >= 0.0)
                    assert abs(np.sum(m * p_at) - th) / th <= 1e-10
                    seg = []
                seg.append((phi, w.copy()))
                if len(seg) >= 3:
                    (p1, w1), (p2, w2), (p3, w3) = seg[-3], seg[-2], seg[-1]
                    w2_lin = w1 + (w3 - w1) * (p2 - p1) / (p3 - p1)
                    assert np.max(np.abs(w2 - w2_lin)) <= 1e-8 * max(np.max(np.abs(w2)), 1e-300)
            else:
                seg = [(phi, w.copy())]
            prev = (phi, w, act, flux)
    assert n_trans >= 20  # the scans genuinely cross activation thresholds


@pytest.mark.physics_invariant
def test_active_set_unique_and_matches_bruteforce():
    """Exactly one candidate active set is KKT-consistent, and it is the solver's.

    For random near-threshold systems every nonempty candidate active set
    is enumerated; exactly one must satisfy both Karush-Kuhn-Tucker
    conditions, and the solver must return that set with the same fluxes
    and multiplier. Partial active sets (the retention branch) must occur
    in the ensemble, or the check would not exercise dropout at all.
    """
    rng = np.random.default_rng(108)
    n_partial, n_total = 0, 0
    for _trial in range(25):
        n, m, X, T, g0, b = _random_system(rng, nmin=3, nmax=5)
        phi1 = first_threshold(X, m, T, g0, b)
        for frac in (0.3, 1.5, 4.0, 12.0, 40.0):
            phi = frac * phi1
            flux, c_mult, act = solve_closure(phi, X, m, T, g0, b, return_diag=True)
            accepted = solve_bruteforce(phi, X, m, T, g0, b)
            n_total += 1
            assert len(accepted) == 1, (n_total, act)
            a_set, w, c_bf = accepted[0]
            assert a_set == act
            wref = max(np.max(np.abs(w)), 1e-300)
            assert np.max(np.abs(X * w - flux)) <= 1e-8 * wref * np.max(X)
            assert c_bf == pytest.approx(c_mult, rel=1e-9)
            if len(a_set) < n:
                n_partial += 1
    assert n_partial > 0  # the retention branch is genuinely exercised


@pytest.mark.physics_invariant
def test_composition_fixed_point_stability_universal_b():
    """Well-mixed-layer fixed points are stable for the universal-b closure.

    Relaxing a layer's composition toward the supply composition at fixed
    total flux, the analytic tangent-space Jacobian at a fully entrained
    universal-b fixed point must have strictly negative real eigenvalues in
    every draw (a runaway fractionation instability would show as a
    positive one), and the analytic Jacobian must match finite differences
    on the first draws (the transcription guard).
    """
    rng = np.random.default_rng(110)
    n_ok = n_run = 0
    for trial in range(600):
        n = int(rng.integers(2, 9))
        m = np.sort(rng.uniform(1, 50, n)) * AMU_G
        X = rng.dirichlet(np.ones(n) * rng.uniform(0.5, 3))
        beta = 10.0 ** rng.uniform(-9, -6)
        mbar = np.sum(X * m)
        margin = rng.uniform(1.02, 20.0)
        a_scale = beta * max(np.max(m) - mbar, 1e-30) * margin + beta * mbar
        w = a_scale - beta * (m - mbar)
        if np.any(w <= 0):
            continue
        a_vec = (beta * (m - mbar) ** 2 - a_scale * m) / mbar
        jac = (
            -np.diag(w)
            - (beta / a_scale) * np.outer(X * (m - mbar), a_vec)
            - beta * np.outer(X, m)
        )
        basis = np.zeros((n, n - 1))
        for j in range(n - 1):
            basis[j, j] = 1.0
            basis[n - 1, j] = -1.0
        jr = np.linalg.lstsq(basis, jac @ basis, rcond=None)[0]
        lam = np.max(np.linalg.eigvals(jr).real)
        n_run += 1
        if lam < 0:
            n_ok += 1
        if trial < 5:
            phi_here = mbar * a_scale - beta * np.sum(X * (m - mbar) ** 2)
            xi = X * w / a_scale
            xi = xi / np.sum(xi)

            def relax(q):
                xq = q / np.sum(q)
                mb = np.sum(xq * m)
                s2 = np.sum(xq * (m - mb) ** 2)
                aq = (phi_here + beta * s2) / mb
                flux_q = xq * (aq - beta * (m - mb))
                return xi * np.sum(flux_q) - flux_q

            jfd = np.zeros((n, n))
            h = 1e-8
            f0 = relax(X)
            for j in range(n):
                qp = X.copy()
                qp[j] += h
                jfd[:, j] = (relax(qp) - f0) / h
            assert np.max(np.abs((jfd - jac) @ basis)) <= 1e-4 * np.max(np.abs(jac @ basis))
    assert n_run > 500
    assert n_ok == n_run  # every admissible fixed point is stable


@pytest.mark.physics_invariant
def test_composition_fixed_point_stability_general_b():
    """General-b composition fixed points are locally stable when entrained.

    Layer relaxation at fixed total mass flux and fixed supply composition,
    general coefficient matrices: numerically constructed fixed points must
    have negative-real-part linearizations in nearly every converged case
    (a small number of draws may fail to converge and are skipped, but
    enough must converge for the check to mean something).
    """
    rng = np.random.default_rng(111)
    n_ok = n_conv = 0
    for _trial in range(12):
        n = int(rng.integers(2, 6))
        m = np.sort(rng.uniform(1, 40, n)) * AMU_G
        T = rng.uniform(300, 2000)
        g0 = rng.uniform(200, 3000)
        b = rng.uniform(0.3, 5, (n, n)) * 1e17 * T**0.75
        b = 0.5 * (b + b.T)
        np.fill_diagonal(b, np.inf)
        kT = kb_cgs * T
        xi = rng.dirichlet(np.ones(n) * 2)
        X = xi.copy()
        phi = None
        converged = False
        for _it in range(600):
            mbar = np.sum(m * X)
            beta_max = np.max(b[np.isfinite(b)]) * g0 / kT
            phi = mbar * beta_max * (np.max(m) - np.min(m)) * 8
            flux = solve_closure(phi, X, m, T, g0, b)
            if np.any(flux <= 0):
                break
            frac = flux / np.sum(flux)
            xn = X * (xi / frac) ** 0.5
            xn /= np.sum(xn)
            if np.max(np.abs(xn - X)) < 1e-13:
                X = xn
                converged = True
                break
            X = xn
        if not converged:
            continue
        n_conv += 1

        def relax(q):
            xq = q / np.sum(q)
            flux_q = solve_closure(phi, xq, m, T, g0, b)
            return xi * np.sum(flux_q) - flux_q

        jac = np.zeros((n, n))
        h = 1e-7
        f0 = relax(X)
        for j in range(n):
            qp = X.copy()
            qp[j] += h
            jac[:, j] = (relax(qp) - f0) / h
        basis = np.zeros((n, n - 1))
        for j in range(n - 1):
            basis[j, j] = 1.0
            basis[n - 1, j] = -1.0
        jr = np.linalg.lstsq(basis, jac @ basis, rcond=None)[0]
        if np.max(np.linalg.eigvals(jr).real) < 0:
            n_ok += 1
    assert n_conv >= 10
    assert n_ok >= n_conv - 1  # at most one marginal linearization tolerated
