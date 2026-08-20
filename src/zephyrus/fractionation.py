"""
!!! info "`fractionation.py`"
    Simultaneous N-species fractionation closure for a hydrodynamic wind.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

import numpy as np

from zephyrus.composition import ELEMENT_AMU
from zephyrus.constants import G, kb_cgs
from zephyrus.diffusion import ROCK_FORMERS, bmatrix, build_rows, masses_g

# The closure generalizes the two-species fractionation of Hunten, Pepin &
# Wallace (1987, Icarus 69, 532) to N species escaping simultaneously
# through mutual binary diffusion, the constant-composition closure of the
# subsonic multispecies wind system of Zahnle et al. (1990, Icarus 84,
# 502), with active-set dropout: at a given total mass flux the heavy
# species partition into an escaping (active) set and a retained set, and
# the retained species exert drag without escaping. The solved system, per
# active species j (w_j = Phi_j / X_j the species velocity scale):
#
#   sum_{i active} X_i (w_i - w_j) / b_ij
#       - w_j sum_{k retained} X_k / b_jk = m_j g0 / kT - C,
#
# plus the mass constraint sum_j m_j X_j Phi_j... i.e.
# sum_{j active} m_j X_j w_j = phi, with C the common inverse scale height
# of the escaping gas (the Lagrange multiplier of the constraint). The
# Karush-Kuhn-Tucker conditions select the active set: active species have
# strictly positive w, and retained species satisfy the retention
# inequality (the drift the escaping gas would impose on them does not
# exceed their gravitational settling).
#
# Exact reductions verified in the companion test suite: the two-species
# limit of Hunten et al. (1987) in the form of Cherubim & Wordsworth
# (2024, ApJ 967, 139, Eqs. 7-9); the three-species deuterium system of Gu
# & Chen (2023, Eqs. 4, 8, 9, 12); the trace-minor relations of Odert et
# al. (2018, Icarus 307, 327, Eq. 5) and Zahnle et al. (1990, Eqs. 35, 36,
# 42); the non-trace three-species relations of Zahnle & Kasting (2023,
# GeCoA 361, 228, Eqs. 19-20); the prescribed-flux partition of
# Chassefiere (1996, Icarus 124, 537, Eqs. 1, 6, 7); the universal-b
# closed form; and the Hunten et al. (1987) Earth, Mars, and Venus
# numerical anchors.
#
# Solver units are cgs (the convention of the source literature and the
# diffusion library): phi in g cm^-2 s^-1, m in g, g0 in cm s^-2, b in
# cm^-1 s^-1, fluxes in cm^-2 s^-1. The public per-species interface
# converts from and to SI at the boundary.


def _validate_inputs(phi, X, m, T, g0, b):
    X, m, b = np.asarray(X, float), np.asarray(m, float), np.asarray(b, float)
    n = len(X)
    if phi < 0:
        raise ValueError('phi must be >= 0')
    if np.any(X < 0) or abs(X.sum() - 1.0) > 1e-6:
        raise ValueError('X must be non-negative and sum to 1')
    if np.any(m <= 0) or T <= 0 or g0 <= 0:
        raise ValueError('m, T, g0 must be positive')
    off = ~np.eye(n, dtype=bool)
    if not np.all(np.isfinite(b[off])) or np.any(b[off] <= 0):
        raise ValueError('off-diagonal b entries must be finite and positive')
    if not np.array_equal(b[off], b.T[off]):
        raise ValueError('b must be symmetric')
    return X, m, b


def solve_fixed_active(phi, X, m, T, g0, b, active):
    """Solve the linear closure on a fixed active set.

    Unknowns: ``w_j`` for j in ``active`` and the common inverse scale
    height ``C``. Returns ``(w_full, C)`` with ``w = 0`` for retained
    species. The single-active case is solved analytically; the general
    case with two-sided diagonal equilibration, which controls the spread
    of roughly 25 decades the matrix entries can span between light-species
    drag terms and heavy-species mass terms.
    """
    kT = kb_cgs * T
    act = sorted(active)
    ret = [k for k in range(len(X)) if k not in active]
    na = len(act)
    if na == 1:
        j = act[0]
        w = np.zeros(len(X))
        w[j] = phi / (m[j] * X[j])
        C = m[j] * g0 / kT + w[j] * sum(X[k] / b[j, k] for k in ret)
        return w, C
    mat = np.zeros((na + 1, na + 1))
    rhs = np.zeros(na + 1)
    for row, j in enumerate(act):
        diag = 0.0
        for col, i in enumerate(act):
            if i == j:
                continue
            mat[row, col] += X[i] / b[i, j]
            diag += X[i] / b[i, j]
        for k in ret:
            diag += X[k] / b[j, k]
        mat[row, row] -= diag
        mat[row, na] = 1.0  # the +C column
        rhs[row] = m[j] * g0 / kT
    for col, j in enumerate(act):
        mat[na, col] = m[j] * X[j]
    rhs[na] = phi
    # Two-sided equilibration: D_r M D_c y = D_r rhs, solution = D_c y.
    dr = 1.0 / np.max(np.abs(mat), axis=1)
    ms = mat * dr[:, None]
    dc = 1.0 / np.max(np.abs(ms), axis=0)
    ms = ms * dc[None, :]
    sol = dc * np.linalg.solve(ms, rhs * dr)
    w = np.zeros(len(X))
    for col, j in enumerate(act):
        w[j] = sol[col]
    return w, sol[na]


def solve_closure(phi, X, m, T, g0, b, return_diag=False):
    """Active-set solution of the N-species fractionation closure.

    Parameters
    ----------
    phi : float
        Total mass flux at the base [g cm^-2 s^-1], non-negative.
    X : array
        Mole fractions, non-negative, summing to 1.
    m : array
        Particle masses [g].
    T : float
        Temperature [K].
    g0 : float
        Gravity at the base [cm s^-2].
    b : array
        Symmetric matrix of binary diffusion parameters [cm^-1 s^-1],
        diagonal ignored (conventionally np.inf).
    return_diag : bool
        When True, also return the multiplier ``C`` and the active set.

    Returns
    -------
    Phi : array
        Number fluxes [cm^-2 s^-1], guaranteed non-negative and satisfying
        ``sum m_i Phi_i = phi``. For ``phi = 0`` the fluxes are zero, the
        active set empty, and ``C`` the continuous limit
        ``min_j m_j g0 / kT``.
    """
    X, m, b = _validate_inputs(phi, X, m, T, g0, b)
    kT = kb_cgs * T
    n = len(X)
    if phi == 0.0:
        flux = np.zeros(n)
        if return_diag:
            return flux, np.min(m) * g0 / kT, frozenset()
        return flux
    active = set(range(n))
    wscale = phi / np.min(m)  # crude magnitude scale for the tolerances
    for _ in range(4 * n + 8):
        w, C = solve_fixed_active(phi, X, m, T, g0, b, active)
        neg = {j for j in active if w[j] < -1e-12 * wscale}
        if neg:
            active -= neg
            if not active:
                raise RuntimeError('empty active set')
            continue
        # Retention check for the inactive species.
        viol, worst = None, 0.0
        for k in range(n):
            if k in active:
                continue
            r_k = sum(X[i] * w[i] / b[i, k] for i in active) - (m[k] * g0 / kT - C)
            if r_k > 1e-12 * abs(m[k] * g0 / kT) and r_k > worst:
                viol, worst = k, r_k
        if viol is None:
            # Clamp within-tolerance negative drifts (magnitudes below the
            # solver's own tolerance) so the returned fluxes honor w >= 0.
            np.maximum(w, 0.0, out=w)
            flux = X * w
            flux[list(set(range(n)) - active)] = 0.0
            if return_diag:
                return flux, C, frozenset(active)
            return flux
        active.add(viol)
    raise RuntimeError('active-set iteration did not converge')


def first_threshold(X, m, T, g0, b):
    """Mass flux at which the first heavy species entrains, in g cm^-2 s^-1.

    Below this flux only the lightest species escapes; the value is the
    smallest crossover over the heavier species of the multi-background
    drag balance.
    """
    kT = kb_cgs * T
    light = int(np.argmin(m))
    best = np.inf
    for k in range(len(X)):
        if k == light:
            continue
        denom = X[light] / b[light, k] + sum(
            X[kk] / b[light, kk] for kk in range(len(X)) if kk != light
        )
        w_star = (m[k] - m[light]) * g0 / kT / denom
        best = min(best, m[light] * X[light] * w_star)
    return best


def closure_per_species(
    mdot: float, element_fractions: dict, T_wind: float, M_p: float, r_base: float
) -> tuple[dict, dict, dict]:
    """Per-element escape rates [kg/s] from the closure at the wind base.

    SI boundary around the cgs solver: the bulk rate ``mdot`` [kg/s]
    converts to the base mass flux, the closure partitions it over the
    atomized composition with the binary-diffusion library coefficients,
    and the number fluxes convert back to per-element mass rates that sum
    to the bulk rate. Rock-forming species (Na, Mg, Si, Fe) raise the
    ``rock_former_bij`` flag: their coefficients sit in the widest
    provenance class, and Na and Mg ionize where rock vapor exists while
    these are neutral-gas coefficients.

    Returns ``(per_element, diagnostics, flags)`` with the diagnostics
    carrying the active set, the retained species, the multiplier, the
    relative mass-conservation residual, and the per-pair coefficient
    provenance.
    """
    flags: dict = {}
    species = sorted(element_fractions, key=lambda el: ELEMENT_AMU[el])
    X = np.array([element_fractions[el] for el in species], float)
    X = X / X.sum()
    m_g = masses_g(species)
    rows = build_rows(species)
    b = bmatrix(species, T_wind, rows=rows)
    prov = {f'{r.i}-{r.j}': (r.provenance, r.uncertainty) for r in rows}
    if any(el in ROCK_FORMERS for el in species):
        flags['rock_former_bij'] = [el for el in species if el in ROCK_FORMERS]

    g0_cgs = (G * M_p / r_base**2) * 1e2  # m/s^2 -> cm/s^2
    phi_cgs = (mdot / (4.0 * math.pi * r_base**2)) * 0.1  # kg/m^2/s -> g/cm^2/s

    flux, C, active = solve_closure(phi_cgs, X, m_g, T_wind, g0_cgs, b, return_diag=True)
    area_cm2 = 4.0 * math.pi * (r_base * 1e2) ** 2
    per_element = {
        el: float(flux[k]) * area_cm2 * float(m_g[k]) * 1e-3 for k, el in enumerate(species)
    }

    total = sum(per_element.values())
    conservation = abs(total - mdot) / mdot if mdot > 0 else 0.0
    diag = dict(
        active_set=sorted(species[k] for k in active),
        retained=[el for el in species if el not in {species[k] for k in active}],
        C_inv_scale_height_cgs=float(C),
        mass_conservation_rel=float(conservation),
        b_provenance=prov,
    )
    return per_element, diag, flags


def unfractionated_split(
    mdot: float, reservoirs: dict | None, element_fractions: dict
) -> tuple[dict, dict]:
    """Split a bulk rate over elements without fractionation.

    The protocol of the energy-limited path: reservoir mass fractions when
    reservoir masses are supplied; otherwise the mass fractions of the
    atomized wind-base composition, with the ``split_from_base_composition``
    flag recording the substitution. The split conserves the bulk rate
    exactly.
    """
    flags: dict = {}
    if reservoirs:
        tot = sum(reservoirs.values())
        fracs = {el: mass / tot for el, mass in reservoirs.items()}
    else:
        flags['split_from_base_composition'] = True
        mass = {el: x * ELEMENT_AMU[el] for el, x in element_fractions.items()}
        tot = sum(mass.values())
        fracs = {el: mm / tot for el, mm in mass.items()}
    return {el: mdot * f for el, f in fracs.items()}, flags
