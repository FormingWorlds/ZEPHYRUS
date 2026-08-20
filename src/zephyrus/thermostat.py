"""
!!! info "`thermostat.py`"
    Wind-temperature thermostat: local heating against radiative cooling.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

import numpy as np

from zephyrus.atomic_data import (
    E_ION_H_EV,
    E_ION_N_EV,
    EV_ERG,
    HC_CM,
    HNU0_H_EV,
    HNU_I_N_EV,
    SIGMA_NU0_H,
    SIGMA_NU_N,
    THREE_LEVEL,
    alpha_case_b,
    badnell_alpha_rr,
    co2_band_cooling,
    o_finestructure_cooling,
)
from zephyrus.constants import kb_cgs

# The thermostat sets the hydrodynamic wind temperature by a local balance
# of photoionization heating against radiative cooling, evaluated at the
# wind-base level (all rates cgs internally):
#
# - Atomic lines: three-level statistical equilibrium for H, C, C+, N, N+,
#   O, and O+ under electron impact with cool-to-space losses, the
#   machinery of Chatterjee & Pierrehumbert (2026, arXiv:2412.05188,
#   their Eqs. 26-30) on the Nakayama et al. (2022) level data; the H
#   system carries Lyman-alpha.
# - The CO2 15 micron band and the atomic O fine structure (Johnstone et
#   al. 2018; see atomic_data for their stated validity limits).
# - Recombination cooling: the continuum part only, Q_RR = n_e n_+
#   alpha_RR(T) (3/2) kB T (Chatterjee & Pierrehumbert 2026, Section 5.2);
#   the exact free-bound emission integral of Tucker & Gould (1966) is not
#   implemented and this form is the stated stand-in.
# - Heating: monochromatic-front photoionization, Q = n_0 sigma
#   (F_XUV / h nu) (h nu - E_ion), with the small thermal correction of
#   order kB T_e against the excess energy dropped. The local ionization
#   fraction comes from photoionization-recombination balance.
#
# The blind spot of a local balance is the temperature structure through
# the sonic region, which no single-level evaluation captures; rates
# computed with this wind temperature inherit that limitation.

T_BRACKET_HIGH = 5.0e4  # K; the lower bracket edge is the equilibrium temperature

_ION_OF = {'C': 'C+', 'N': 'N+', 'O': 'O+'}


def three_level_populations(species: str, n_tot: float, n_e: float, T: float):
    """Steady-state level populations (n1, n2, n3) of a three-level system.

    Electron-impact excitation and de-excitation plus radiative decay, with
    no radiative excitation (every emitted photon escapes: cool-to-space).
    Rate coefficients follow the effective-collision-strength form
    ``k_lu = gamma (8.629e-6 / (g_l sqrt(T))) exp(-E_lu / kB T)`` with
    de-excitation by detailed balance (Nakayama et al. 2022, Eqs. 13-14).
    Densities in cm^-3.
    """
    data = THREE_LEVEL[species]
    (_t1, g1, e1), (_t2, g2, e2), (_t3, g3, e3) = data['levels']
    energy = {
        (1, 2): HC_CM * (e2 - e1) * 1e7,  # J -> erg
        (1, 3): HC_CM * (e3 - e1) * 1e7,
        (2, 3): HC_CM * (e3 - e2) * 1e7,
    }
    g = {1: g1, 2: g2, 3: g3}
    c_up, c_dn, a_rad = {}, {}, {}
    for (lo, up), (a_ul, gam) in data['transitions'].items():
        e = energy[(lo, up)]
        k_lu = gam * 8.629e-6 / (g[lo] * math.sqrt(T)) * math.exp(-e / (kb_cgs * T))
        k_ul = gam * 8.629e-6 / (g[up] * math.sqrt(T))
        c_up[(lo, up)] = n_e * k_lu
        c_dn[(up, lo)] = n_e * k_ul
        a_rad[(up, lo)] = a_ul
    # Balance equations for levels 2 and 3 with the closure n1+n2+n3 = n_tot.
    a22 = -(c_dn[(2, 1)] + a_rad[(2, 1)] + c_up[(2, 3)])
    a23 = c_dn[(3, 2)] + a_rad[(3, 2)]
    a32 = c_up[(2, 3)]
    a33 = -(c_dn[(3, 1)] + c_dn[(3, 2)] + a_rad[(3, 1)] + a_rad[(3, 2)])
    b2, b3 = c_up[(1, 2)], c_up[(1, 3)]
    det = a22 * a33 - a23 * a32
    if det == 0.0:
        return n_tot, 0.0, 0.0
    r2 = (-b2 * a33 + b3 * a23) / det
    r3 = (-a22 * b3 + a32 * b2) / det
    n1 = n_tot / (1.0 + r2 + r3)
    return n1, r2 * n1, r3 * n1


def three_level_cooling(species: str, n_tot: float, n_e: float, T: float) -> float:
    """Line cooling of one species, erg s^-1 cm^-3, cool-to-space."""
    if n_tot <= 0.0 or n_e <= 0.0:
        return 0.0
    data = THREE_LEVEL[species]
    levels = data['levels']
    n = dict(zip((1, 2, 3), three_level_populations(species, n_tot, n_e, T)))
    q = 0.0
    for (lo, up), (a_ul, _gam) in data['transitions'].items():
        if a_ul <= 0.0:
            continue
        e = HC_CM * (levels[up - 1][2] - levels[lo - 1][2]) * 1e7
        q += n[up] * a_ul * e
    return q


def ionization_fraction(
    n_cgs: float, T: float, F_cgs: float, hnu_erg: float, sigma_cm2: float, alpha_cm3s: float
) -> float:
    """Local photoionization-recombination balance ionization fraction.

    Solves ``f^2 / (1 - f) = R`` with ``R = sigma F / (h nu alpha n)``, the
    quadratic root in [0, 1], written as ``f = 2 / (1 + sqrt(1 + 4/R))``:
    the textbook form ``(-R + sqrt(R^2 + 4R)) / 2`` cancels catastrophically
    at large R and can return fractions above 1. Zero flux or zero density
    gives zero.
    """
    if F_cgs <= 0.0 or n_cgs <= 0.0:
        return 0.0
    R = sigma_cm2 * F_cgs / (hnu_erg * alpha_cm3s * n_cgs)
    return 2.0 / (1.0 + math.sqrt(1.0 + 4.0 / R))


def front_constants(element_fractions: dict) -> tuple[float, float, float]:
    """Monochromatic-front constants (h nu [erg], E_ion [erg], sigma [cm^2]).

    Hydrogen front for winds with an atomized hydrogen fraction of one half
    or more, the nitrogen-like front otherwise. Treating the base gas as
    atomized in front of the ionizing continuum is an approximation the
    front constants inherit.
    """
    if element_fractions.get('H', 0.0) >= 0.5:
        return HNU0_H_EV * EV_ERG, E_ION_H_EV * EV_ERG, SIGMA_NU0_H
    return HNU_I_N_EV * EV_ERG, E_ION_N_EV * EV_ERG, SIGMA_NU_N


def recombination_alpha(element_fractions: dict, T: float) -> float:
    """Composition recombination coefficient, cm^3 s^-1.

    Nitrogen-dominated winds use the Badnell fit for N; other compositions
    take the mole-fraction-weighted case B coefficients.
    """
    if element_fractions.get('N', 0.0) >= 0.5:
        return badnell_alpha_rr(T)
    return sum(x * alpha_case_b(el, T) for el, x in element_fractions.items())


def balance_at(
    T: float,
    base: dict,
    element_fractions: dict,
    F_xuv: float,
    cool_atomic: bool = True,
    cool_co2_band: bool = True,
    cool_o_finestructure: bool = True,
    cool_recombination: bool = True,
) -> tuple[float, dict]:
    """Heating minus cooling [erg s^-1 cm^-3] at temperature T.

    Evaluated on the wind-base level ``base`` (a level dict with ``n``
    [m^-3] and ``vmr``; densities converted to cgs internally), for the
    atomized composition ``element_fractions`` under the XUV flux ``F_xuv``
    [W m^-2]. Returns ``(residual, detail)`` with the detail carrying the
    ionization fraction, both totals, and the per-channel parts.
    """
    n_cgs = base['n'] * 1e-6
    F_cgs = F_xuv * 1e3
    hnu, e_ion, sigma = front_constants(element_fractions)
    alpha = recombination_alpha(element_fractions, T)
    f_plus = ionization_fraction(n_cgs, T, F_cgs, hnu, sigma, alpha)
    n_e = f_plus * n_cgs

    q_heat = (1.0 - f_plus) * n_cgs * sigma * (F_cgs / hnu) * max(hnu - e_ion, 0.0)

    q_cool = 0.0
    parts = {}
    if cool_atomic:
        q_lines = 0.0
        for el, x in element_fractions.items():
            if el not in THREE_LEVEL:
                continue
            n_el = x * n_cgs
            q_lines += three_level_cooling(el, (1.0 - f_plus) * n_el, n_e, T)
            ion = _ION_OF.get(el)
            if ion:
                q_lines += three_level_cooling(ion, f_plus * n_el, n_e, T)
        parts['atomic_lines'] = q_lines
        q_cool += q_lines
    if cool_co2_band:
        vmr = base['vmr']
        n_co2 = vmr.get('CO2', 0.0) * n_cgs
        colliders = {
            sp: vmr.get(sp, 0.0) * n_cgs for sp in ('O', 'O2', 'N2', 'CO2', 'He', 'Ar')
        }
        q_co2 = co2_band_cooling(n_co2, colliders, T) if n_co2 > 0.0 else 0.0
        parts['co2_band'] = q_co2
        q_cool += q_co2
    if cool_o_finestructure:
        n_o = element_fractions.get('O', 0.0) * (1.0 - f_plus) * n_cgs
        q_o = o_finestructure_cooling(n_o, T) if n_o > 0.0 else 0.0
        parts['o_finestructure'] = q_o
        q_cool += q_o
    if cool_recombination:
        q_rr = n_e * (f_plus * n_cgs) * alpha * 1.5 * kb_cgs * T
        parts['recombination'] = q_rr
        q_cool += q_rr

    return q_heat - q_cool, dict(
        f_plus=f_plus, q_heat=q_heat, q_cool=q_cool, parts=parts, alpha_rec_cgs=alpha
    )


def solve_wind_temperature(
    T_eq: float,
    base: dict,
    element_fractions: dict,
    F_xuv: float,
    cool_atomic: bool = True,
    cool_co2_band: bool = True,
    cool_o_finestructure: bool = True,
    cool_recombination: bool = True,
    n_scan: int = 120,
) -> tuple[float, dict]:
    """Thermostat root-find on the bracket [T_eq, 5e4 K].

    Scans the balance upward in temperature and bisects the first downward
    crossing of heating through cooling, which selects the lowest stable
    root when several exist. When no root lies inside the bracket the
    temperature clamps to the nearer edge with the ``clamped`` field set
    ('low' when cooling already wins at T_eq, 'high' when heating still
    wins at 5e4 K, where the missing physics is the ionization and line
    inventory beyond the modeled channels). Returns ``(T_wind, detail)``.

    Raises
    ------
    ValueError
        If every cooling channel is disabled: a pure-heating balance has no
        root by construction and would silently clamp high.
    """
    if not (cool_atomic or cool_co2_band or cool_o_finestructure or cool_recombination):
        raise ValueError('all cooling channels disabled; at least one must stay on')
    channels = dict(
        cool_atomic=cool_atomic,
        cool_co2_band=cool_co2_band,
        cool_o_finestructure=cool_o_finestructure,
        cool_recombination=cool_recombination,
    )
    lo, hi = T_eq, T_BRACKET_HIGH
    if lo >= hi:
        _, d = balance_at(hi, base, element_fractions, F_xuv, **channels)
        return hi, dict(clamped='high', **d)
    ts = np.geomspace(lo, hi, n_scan)
    res = [balance_at(float(t), base, element_fractions, F_xuv, **channels)[0] for t in ts]
    idx = None
    for i in range(len(ts) - 1):
        if res[i] > 0.0 >= res[i + 1]:
            idx = i
            break
    if idx is None:
        if res[0] <= 0.0:
            t_w, clamp = lo, 'low'  # cooling already wins at the equilibrium temperature
        else:
            t_w, clamp = hi, 'high'  # heating still wins at the top of the bracket
        _, d = balance_at(t_w, base, element_fractions, F_xuv, **channels)
        return t_w, dict(clamped=clamp, **d)
    a, b = float(ts[idx]), float(ts[idx + 1])
    for _ in range(60):
        m = 0.5 * (a + b)
        r, _ = balance_at(m, base, element_fractions, F_xuv, **channels)
        if r > 0.0:
            a = m
        else:
            b = m
    t_w = 0.5 * (a + b)
    _, d = balance_at(t_w, base, element_fractions, F_xuv, **channels)
    return t_w, dict(clamped=None, **d)
