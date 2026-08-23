"""
!!! info "`diagnostics.py`"
    Regime diagnostics reported beside every dispatch verdict.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

import numpy as np

from zephyrus.constants import G, kb, rate_floor
from zephyrus.hydrodynamic import k_tide
from zephyrus.knudsen import mean_free_path, sigma_mixture
from zephyrus.planets_parameters import Mjup, Rjup

# Everything in this module is reporting: the quantities let a reader
# translate a regime verdict into the criteria other escape taxonomies use,
# and quantify how close each call sat to its boundaries. Nothing here is
# read back by the dispatch control flow, and the container has no off
# switch: the regime boundaries carry genuine physical uncertainty, and
# printing the translation quantities beside every verdict is the
# mitigation.

# Murray-Clay et al. (2009) fit their numerical models with flux exponents
# 0.6 (radiation-recombination limited) and 0.9 (energy limited); the
# analytic chains implemented here carry the idealized 0.5 and 1.0.
# Rate-level comparisons against their models must budget the difference.
MURRAY_CLAY_EXPONENTS = {
    'RR_numerical': 0.6,
    'RR_analytic_inherited': 0.5,
    'EL_numerical': 0.9,
    'EL_analytic_inherited': 1.0,
}
# Their regime-dependent dayside-heating reduction factors relative to
# full-surface redistribution, for sensitivity analyses.
DAYSIDE_FACTORS = {'energy_limited': 0.26, 'recombination_limited': 0.31}

# Threshold gravitational potentials, log10(-phi) in cgs (erg/g), with
# phi = -G M_p / R_p, the convention both sources use. The Caldiroli et al.
# (2022) band marks where the evaporation efficiency collapses. The second
# screen separates wind-forming from hydrostatic thermospheres and is
# Salz et al. (2016, A&A 585, L2), whose photoionization hydrodynamics
# simulations find the energy-limited concept valid below 13.11, because
# the radiative input is efficiently spent driving the wind, and stable
# thermospheres above about 13.6, because the whole input is re-emitted in
# hydrogen Lyman alpha (above roughly 1.1 R_p) and free-free emission
# (below it). Between the two the wind weakens as the heating efficiency
# falls. Their grid is hydrogen-dominated thermospheres of hot gas planets,
# from super-Earth-sized to massive hot Jupiters, so the screen is out of
# its own scope on a heavy secondary atmosphere and is reported, never
# applied.
CALDIROLI_THRESHOLD_LOG_PHI = (12.9, 13.2)
SALZ_SCREEN_LOG_PHI = (13.11, 13.6)

# The published name of the rate floor, defined once in constants.
RATE_FLOOR_KG_S = rate_floor


def q_net_over_qc(
    eps: float,
    F_xuv: float,
    R_xuv: float,
    r_sonic: float,
    r_base: float,
    M_p: float,
    m_mean: float,
    sigma_c: float,
    gamma: float = 1.0,
    kn_m: float = 1.0,
) -> tuple[float, float, float]:
    """Transonic energy criterion of Johnson et al. (2013, ApJL 768, L4).

    Their Eq. (10) critical power for sustaining a transonic outflow,
    ``Q_c = 4 pi r_* (gamma / (c_c sigma_c Kn_m)) sqrt(2 U(r_*) / m)
    U(r_0)`` with ``c_c = sqrt(2)`` and ``U(r) = G M m / r``, against the
    intercepted, efficiency-degraded XUV power
    ``Q_net = eps pi R_XUV^2 F_XUV``. A ratio well below 1 says the heating
    cannot drive the flow transonic regardless of what a rate formula
    returns. All SI.
    """
    q_net = eps * math.pi * R_xuv**2 * F_xuv
    u_star = G * M_p * m_mean / r_sonic
    u_0 = G * M_p * m_mean / r_base
    q_c = (
        4.0
        * math.pi
        * r_sonic
        * gamma
        / (math.sqrt(2.0) * sigma_c * kn_m)
        * math.sqrt(2.0 * u_star / m_mean)
        * u_0
    )
    return q_net / q_c, q_net, q_c


def guo_triple(
    M_p: float,
    R_p: float,
    T_eq: float,
    mu_kg: float,
    M_star: float,
    a: float,
    e: float,
    lambda_exo: float,
) -> dict:
    """The (lambda_exo, lambda, lambda*) regime triple of Guo (2024).

    Guo (2024, Nat. Astron. 8, 920) classifies escape regimes with the Jeans
    parameter at the planetary radius and its Roche-corrected companion
    ``lambda* = lambda K(xi)``; reporting the triple beside the exobase
    value lets a reader translate the dispatch verdict into that taxonomy.
    The orbital distance is taken at periapsis for consistency with the
    tidal factor elsewhere in the package (a distinction that vanishes on
    circular orbits).
    """
    lam = G * M_p * mu_kg / (kb * T_eq * R_p)
    d = a * (1.0 - e)
    f = (M_p / (3.0 * M_star)) ** (1.0 / 3.0) * d / R_p
    lam_star = lam * k_tide(f) if f > 1.0 else 0.0
    return dict(
        lambda_exo=lambda_exo,
        lambda_rp=lam,
        lambda_star=lam_star,
        thresholds='thermally driven lambda < ~3; tidal lambda* < 3; XUV lambda* > 6',
    )


def erkaev_tc(M_p: float, R_p: float, r_exo: float, R_hill: float) -> float:
    """Tidally corrected critical exobase temperature, in K.

    Erkaev et al. (2007, Eq. 23): the exobase temperature above which the
    thermosphere blows off, ``T_c = T_Jup (M_p / M_Jup)(R_Jup / R_p)
    K(x_Rl / x) / x`` with ``x = r_exo / R_p``, ``x_Rl = R_Hill / R_p``,
    and their Jupiter normalization 1.45e5 K. Returns 0 when the exobase
    reaches the Roche lobe (the barrier is gone).
    """
    x = r_exo / R_p
    x_rl = R_hill / R_p
    if x_rl / x <= 1.0:
        return 0.0
    return 1.45e5 * (M_p / Mjup) * (Rjup / R_p) * k_tide(x_rl / x) / x


def along_profile_fluid_check(
    profile, M_p: float, r_sonic: float, kn_threshold: float = 1.0
) -> dict:
    """Fluid condition along the profile up to the sonic surface.

    Owen & Jackson (2012, MNRAS 425, 2931) require the fluid condition to
    hold everywhere below the sonic surface, not only at it. The check
    walks the supplied profile levels below the sonic radius and records
    the worst local Knudsen number (Maxwell mean free path over the local
    scale height); the truncation at the profile top is declared, because
    the sonic point normally lies above the modeled atmosphere.
    """
    worst = 0.0
    checked = 0
    for i in range(len(profile.p)):
        r = float(profile.r[i])
        if r > r_sonic:
            break
        T = float(profile.T[i])
        n = float(profile.p[i]) / (kb * T)
        vmr = {sp: float(np.asarray(v)[i]) for sp, v in profile.vmr.items()}
        sigma, _ = sigma_mixture(vmr, T)
        H = kb * T * r**2 / (G * M_p * float(profile.mmw[i]))
        kn = mean_free_path(sigma, n) / H
        worst = max(worst, kn)
        checked += 1
    truncated = bool(len(profile.p) and float(profile.r[-1]) < r_sonic)
    return dict(
        levels_checked=checked,
        worst_kn=worst,
        fluid=(worst < kn_threshold),
        truncated_at_profile_top=truncated,
    )


def self_consistency_screen(reservoirs: dict | None, mdot: float, age: float | None) -> dict:
    """Static-snapshot screen: depletion timescale against the snapshot age.

    A dispatched rate that would empty the supplied reservoirs in less than
    the system age flags the snapshot as inconsistent with its own history
    (the atmosphere could not have survived to be observed in this state).
    Reports ``evaluated: False`` when age or reservoirs are absent.
    """
    if not reservoirs or age is None or mdot <= 0.0:
        return {'evaluated': False}
    t_dep = sum(reservoirs.values()) / mdot
    return {'evaluated': True, 't_deplete_s': t_dep, 'age_s': age, 'inconsistent': t_dep < age}


def rate_floor_screen(mdot: float) -> dict:
    """Numerical-content screen: the rate against one proton per year.

    A strongly bound heavy atmosphere returns rates many decades below
    anything with physical content, and a regime label attached to such a
    rate is decided by the ordering of two meaningless numbers. One proton
    crossing the surface per Julian year is the smallest rate worth
    reading. Reporting only: the module never applies the floor, because
    what counts as negligible belongs to the caller, and clearing the floor
    does not make a rate matter (for that, use
    :func:`self_consistency_screen`).
    """
    return {'floor_kg_s': RATE_FLOOR_KG_S, 'above_floor': mdot > RATE_FLOOR_KG_S}


def potential_screens(M_p: float, R_p: float) -> dict:
    """Threshold-potential screens, log10(-phi) in cgs (erg/g).

    Reports where the configuration sits against the Caldiroli et al.
    (2022) efficiency-collapse band and the Salz et al. (2016) screen
    separating wind-forming from hydrostatically stable thermospheres (see
    the module constants for what each threshold means and for the
    hydrogen-dominated scope of the second).
    """
    log_phi = math.log10(G * M_p / R_p * 1e4)
    return dict(
        log_minus_phi_cgs=log_phi,
        caldiroli_threshold=CALDIROLI_THRESHOLD_LOG_PHI,
        above_caldiroli=log_phi > CALDIROLI_THRESHOLD_LOG_PHI[0],
        salz_screen=SALZ_SCREEN_LOG_PHI,
        salz_verdict=(
            'wind'
            if log_phi < SALZ_SCREEN_LOG_PHI[0]
            else 'no-wind'
            if log_phi > SALZ_SCREEN_LOG_PHI[1]
            else 'intermediate'
        ),
        salz_attribution=(
            'Salz et al. (2016, A&A 585, L2): energy-limited escape valid below '
            '13.11, hydrodynamically stable thermospheres above about 13.6, the '
            'wind weakening in between; derived for hydrogen-dominated '
            'thermospheres, so out of scope on a heavy secondary atmosphere'
        ),
    )
