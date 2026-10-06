"""
!!! info "`atomic_data.py`"
    Atomic and molecular cooling data and closed-form rate coefficients.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

from zephyrus.constants import c as C_LIGHT
from zephyrus.constants import h_planck

# hc in J cm: converts a wavenumber in cm^-1 to an energy in J.
HC_CM = h_planck * C_LIGHT * 100.0

# Three-level systems for atomic line cooling: levels 1 to 3 of each species
# from Nakayama, Ikoma & Terada (2022, ApJ 937, 72), Appendix C Tables 2 to
# 5. levels: (term, statistical weight g, excitation energy [cm^-1]);
# transitions keyed (lower, upper) with 1-based indices:
# (Einstein A [s^-1], effective collision strength at 1e4 K).
# The printed O+ row carries neutral-O term labels beside the O+ weights 4,
# 10, and 6, so its level set is entered as 4S-2D-2P with those weights.
THREE_LEVEL = {
    'H': {
        'levels': [('1s2S', 2, 0.0), ('2s2S', 2, 82258.96), ('2p2P', 6, 82259.17)],
        'transitions': {
            (1, 2): (2.50e-6, 2.42e-1),
            (1, 3): (6.26e8, 5.00e-1),
            (2, 3): (0.0, 0.0),
        },
    },
    'C': {
        'levels': [('3P', 9, 29.59122), ('1D', 5, 10192.67), ('1S', 1, 21648.04)],
        'transitions': {
            (1, 2): (2.43e-4, 1.21),
            (1, 3): (2.13e-3, 7.39e-2),
            (2, 3): (6.38e-1, 3.90e-1),
        },
    },
    'C+': {
        'levels': [('2P', 6, 42.26666), ('4P', 12, 43035.75), ('2D', 10, 74931.60)],
        'transitions': {(1, 2): (4.57e1, 6.57), (1, 3): (2.90e7, 2.92), (2, 3): (0.0, 1.94)},
    },
    'N': {
        'levels': [('4S', 4, 0.0), ('2D', 10, 19227.95), ('2P', 6, 28838.51)],
        'transitions': {
            (1, 2): (1.30e-5, 5.61e-1),
            (1, 3): (5.22e-3, 1.64e-1),
            (2, 3): (8.47e-2, 4.37e-1),
        },
    },
    'N+': {
        'levels': [('3P', 9, 85.22956), ('1D', 5, 16455.11), ('1S', 1, 33218.58)],
        'transitions': {
            (1, 2): (3.90e-3, 1.38),
            (1, 3): (3.20e-2, 8.00e-1),
            (2, 3): (1.14e0, 5.12),
        },
    },
    'O': {
        'levels': [('3P', 9, 76.83111), ('1D', 5, 15868.34), ('1S', 1, 33792.22)],
        'transitions': {
            (1, 2): (8.57e-3, 2.93e-1),
            (1, 3): (7.87e-2, 3.23e-2),
            (2, 3): (1.26e0, 8.83e-3),
        },
    },
    'O+': {
        'levels': [('4S', 4, 0.0), ('2D', 10, 27826.09), ('2P', 6, 42125.60)],
        'transitions': {
            (1, 2): (7.68e-5, 1.33),
            (1, 3): (4.51e-2, 4.06e-1),
            (2, 3): (9.68e-2, 1.70),
        },
    },
}

# Radiative recombination: the Badnell (2006, ApJS 167, 334) fit,
#   alpha_RR = A / [ sqrt(T/T0) (1 + sqrt(T/T0))^(1-B') (1 + sqrt(T/T1))^(1+B') ],
#   B' = B + C exp(-T2/T),
# from the original Eqs. (1)-(2), with the nitrogen row Z = 7, N = 6 of
# Badnell's table. Chatterjee & Pierrehumbert (2026, ApJ 998, 236) Eq. 35
# misprints both the form and T2 (6.379e4 for 6.739e4); neither is used.
# (T0, T1, T2 [K], A [cm^3/s], B, C) for nitrogen.
BADNELL_N = (9.467e-2, 2.954e6, 6.739e4, 6.387e-10, 0.7308, 0.2440)


def badnell_alpha_rr(T: float, coeffs: tuple = BADNELL_N) -> float:
    """Radiative recombination coefficient, cm^3 s^-1 (Badnell 2006, Eqs. 1-2)."""
    T0, T1, T2, A, B, C = coeffs
    bp = B + C * math.exp(-T2 / T)
    s0 = math.sqrt(T / T0)
    s1 = math.sqrt(T / T1)
    return A / (s0 * (1.0 + s0) ** (1.0 - bp) * (1.0 + s1) ** (1.0 + bp))


# Case B at 1e4 K, cm^3 s^-1: H from Murray-Clay et al. (2009, ApJ 693, 23)
# Eq. (7); heavies as AMDPP case A totals minus the ground-state partial,
# as compiled for Malina Ovesen's radiation-recombination module.
CASE_B_1E4K = {
    'H': 2.7e-13,
    'He': 4.37e-13 - 1.56e-13,
    'C': 4.72e-13 - 2.32e-13,
    'N': 3.76e-13 - 1.15e-13,
    'O': 2.72e-13 - 1.31e-13,
}


def alpha_case_b(element: str, T: float) -> float:
    """Case B recombination coefficient, cm^3 s^-1, with temperature scaling.

    Applies the Murray-Clay et al. (2009) hydrogen temperature dependence
    ``(T / 1e4 K)^-0.9`` to every element: the heavies' archival values
    exist only at 1e4 K, so extending hydrogen's exponent to them is a
    documented approximation that matters once the thermostat moves the
    wind temperature away from 1e4 K. Elements without a tabulated value
    take the atomic-O coefficient (the nearest heavy; callers record the
    coefficient provenance upstream).
    """
    a0 = CASE_B_1E4K.get(element)
    if a0 is None:
        a0 = CASE_B_1E4K['O']
    return a0 * (T / 1.0e4) ** -0.9


# CO2 15 micron band: Johnstone et al. (2018, A&A 617, A107) Eqs. (34)-(38),
# Table 1 colliders, cgs. Their exp(-667/T_n) takes the 667 cm^-1 wavenumber
# for a temperature; detailed balance here uses h nu / k_B = 959.7 K.
HNU_15UM = 1.325e-13  # erg, the 15 micron quantum
T_15UM = 959.7  # K, that quantum over the Boltzmann constant
A10_CO2 = 0.46  # s^-1, Einstein coefficient of the bending mode
SIGMA_CO2_15UM = 6.43e-15  # cm^2, band column parameter for the escape probability
CO2_KD = {
    'O': (5.10e-11, -0.59),
    'O2': (4.97e-22, 2.83),
    'N2': (6.43e-21, 2.30),
    'CO2': (4.21e-17, 0.85),
    'He': (4.73e-19, 2.19),
    'Ar': (8.13e-24, 3.19),
}


def co2_band_cooling(n_co2: float, colliders: dict, T: float, col_co2: float = 0.0) -> float:
    """CO2 15 micron band volumetric cooling, erg s^-1 cm^-3.

    Johnstone et al. (2018) Eqs. (34)-(38) with no stellar infrared pumping
    (their S_IR term is zero in the escaping-region application).
    ``n_co2`` [cm^-3]; ``colliders`` maps species to number densities
    [cm^-3]; ``col_co2`` is the overlying CO2 column [cm^-2] for the escape
    probability, whose zero-column limit is the non-LTE ceiling 0.5.
    Excitation rates follow from detailed balance,
    ``k_e = 2 k_d exp(-h nu / k_B T)``.

    The abundances are molecular because this is a base-region coolant, on
    deexcitation fits measured over roughly 150 to 500 K. Above CO2
    dissociation the band should be absent rather than small, which this
    function does not enforce.
    """
    sn = SIGMA_CO2_15UM * col_co2
    if sn > 2.0:
        eps = 0.7202 * sn**-0.613
    elif sn > 0.0:
        # The fit rises through 0.5 below sn = 3.4e-4 and diverges; half the
        # photons escaping is the non-LTE ceiling of this two-level band,
        # and capping there keeps the zero-column case continuous.
        eps = min(0.4732 * sn**-0.0069, 0.5)
    else:
        eps = 0.5
    ke_sum = kd_sum = 0.0
    for sp, n in colliders.items():
        if sp not in CO2_KD or n <= 0.0:
            continue
        a, b = CO2_KD[sp]
        kd = a * T**b
        ke = 2.0 * kd * math.exp(-T_15UM / T)
        kd_sum += kd * n
        ke_sum += ke * n
    if ke_sum == 0.0:
        return 0.0
    n_star = ke_sum * n_co2 / (ke_sum + kd_sum + A10_CO2 * eps)
    return HNU_15UM * A10_CO2 * eps * n_star


def o_finestructure_cooling(n_o: float, T: float) -> float:
    """Atomic O fine-structure cooling at 63 and 147 micron, erg s^-1 cm^-3.

    The closed LTE form of Johnstone et al. (2018) Eqs. (41)-(43), tracing
    to Bates (1951); ``n_o`` is the atomic oxygen density [cm^-3].
    """
    den = 1.0 + 0.6 * math.exp(-228.0 / T) + 0.2 * math.exp(-326.0 / T)
    q63 = 1.67e-18 * math.exp(-228.0 / T) * n_o / den
    q147 = 4.59e-20 * math.exp(-326.0 / T) * n_o / den
    return q63 + q147


# Monochromatic front constants: hydrogen after Murray-Clay et al. (2009),
# sigma_nu0 = 6e-18 (h nu0 / 13.6 eV)^-3 cm^2 at 20 eV; the nitrogen-like
# front of Chatterjee & Pierrehumbert (2026), N I potential from NIST.
EV_ERG = 1.602176634e-12
HNU0_H_EV = 20.0
E_ION_H_EV = 13.6
SIGMA_NU0_H = 6.0e-18 * (HNU0_H_EV / 13.6) ** -3  # cm^2
HNU_I_N_EV = 33.6
E_ION_N_EV = 14.53
SIGMA_NU_N = 1.0e-17  # cm^2

# Black (1981) Lyman-alpha cooling as printed by Murray-Clay et al. (2009)
# Eq. (6), 7.5e-19 n_e n_H exp(-118348 K / T) in cgs: a cross-check
# constant only, as the H three-level system above carries Lyman-alpha.
LYA_BLACK = (7.5e-19, 118348.0)
