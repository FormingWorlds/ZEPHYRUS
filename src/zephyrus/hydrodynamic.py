"""
!!! info "`hydrodynamic.py`"
    Hydrodynamic escape: energy-limited and radiation-recombination-limited
    rates, and the selection between them.<br>
    Authors: Malina Ovesen, Mara Attia
"""

from __future__ import annotations

import math

from zephyrus.atomic_data import HNU0_H_EV, HNU_I_N_EV, alpha_case_b
from zephyrus.composition import ELEMENT_AMU
from zephyrus.constants import G, ev2joule, kb, m_p
from zephyrus.profiles import SIGMA_NU0

# The branch computes both hydrodynamic limits and takes their minimum:
#
# - Energy limited (EL): Erkaev et al. (2007, A&A 472, 329, their Eq. 21),
#   Mdot = eps pi F_XUV R_p R_XUV^2 / (G M_p K(xi)), with the tidal factor
#   K(xi) of their Eq. (17) at xi = R_Hill/R_p and the periapsis Hill
#   radius. The factor pi encodes full-surface redistribution of the
#   intercepted power.
# - Radiation-recombination limited (RR): the analytic chain of Murray-Clay
#   et al. (2009, ApJ 693, 23, Section 3.2) in the form derived by Malina
#   Ovesen from Lopez (2017, MNRAS 472, 245, Eqs. 4-6): ionization
#   equilibrium at the wind base sets the base ion density proportional to
#   sqrt(F_XUV), and an isothermal Parker wind carries it to the sonic
#   point with the barometric factor exp(3/2 - lambda_b), the exact
#   isothermal value. min(EL, RR) selects RR two physically distinct ways:
#   genuine recombination saturation (the sqrt(F) regime at modest
#   lambda_b) and barometric suppression at large lambda_b, where the label
#   "recombination limited" would be a category error; the selection
#   mechanism is reported so the two are never conflated. Where the two
#   candidates cross in flux is sensitive to the wind temperature: the RR
#   chain carries it through the sound speed, the barometric exponent, and
#   the recombination coefficient, so a thermostat-driven wind temperature
#   can move the EL/RR crossover by an order of magnitude against the
#   canonical fixed 1e4 K evaluation.
# - Efficiency: fixed, or the Caldiroli et al. (2022, A&A 663, A122,
#   Appendix A.1) fit, defined against their R_p^3 geometry and therefore
#   converted by (R_p/R_XUV)^2 before use in the Erkaev form.

RHO_UNIT_CGS = 1e-3  # kg m^-3 -> g cm^-3
FLUX_UNIT_CGS = 1e3  # W m^-2 -> erg s^-1 cm^-2


def hill_radius_periapsis(M_p: float, M_star: float, a: float, e: float) -> float:
    """Periapsis Hill radius a (1 - e) (M_p / 3 M_star)^(1/3), in m."""
    return a * (1.0 - e) * (M_p / (3.0 * M_star)) ** (1.0 / 3.0)


def k_tide(xi: float) -> float:
    """Erkaev et al. (2007) Eq. (17) tidal factor; valid for xi > 1 (K(1) = 0)."""
    return 1.0 - 3.0 / (2.0 * xi) + 1.0 / (2.0 * xi**3)


def el_rate(eps: float, F_xuv: float, R_p: float, R_xuv: float, M_p: float, K: float) -> float:
    """Energy-limited rate eps pi F R_p R_xuv^2 / (G M_p K), in kg/s."""
    return eps * math.pi * F_xuv * R_p * R_xuv**2 / (G * M_p * K)


def caldiroli_efficiency(F_xuv: float, M_p: float, R_p: float, K: float) -> tuple:
    """Evaporation-efficiency fit of Caldiroli et al. (2022, Appendix A.1).

    Their fit is a function of the tidally corrected gravitational
    potential ``phi = K G M_p / R_p`` and the flux-to-density ratio
    ``F_XUV / rho_p``, both in cgs internally. It returns the efficiency
    defined against their ``R_p^3`` rate geometry; the caller converts by
    ``(R_p / R_XUV)^2`` before using it in the Erkaev form. Below their
    validity bound ``F_XUV / rho_p = 1e2`` (cgs) the fitting formulas turn
    complex, so that region is rejected here: the return is ``(None,
    flags)`` with ``caldiroli_below_flux_bound`` set and the caller falls
    back to the fixed efficiency. Outside their fitted box the value is
    still returned, flagged ``caldiroli_out_of_box``.
    """
    flags = {}
    rho_p = M_p / (4.0 / 3.0 * math.pi * R_p**3)
    f_cgs = F_xuv * FLUX_UNIT_CGS
    rho_cgs = rho_p * RHO_UNIT_CGS
    f2 = (f_cgs / rho_cgs) / 1e2
    if f2 < 1.0:
        flags['caldiroli_below_flux_bound'] = True
        return None, flags
    phi_red = K * (G * M_p / R_p) * 1e4  # erg/g
    if not (10**12.17 <= phi_red <= 10**13.29) or f2 > 1e4:
        flags['caldiroli_out_of_box'] = True
    lf2 = math.log10(f2)
    a_coef = 1.682 * lf2**0.2802 - 5.488 if lf2 > 0 else -5.488
    alpha = 0.02489 * f2**-0.0860 - 0.01007 * f2**-0.9543
    eta0 = -0.03973 * lf2**2.173 - 0.01359 if lf2 > 0 else -0.01359
    beta = -0.01799 * f2**0.1723 - 3.3875 * f2**0.0140
    sigma = 1.0 / (1.0 + (phi_red / 10**13.22) ** beta)
    log_eta = a_coef * phi_red**alpha * sigma + eta0 * (1.0 - sigma)
    return 10**log_eta, flags


def wind_mean_masses(element_fractions: dict) -> tuple[float, float]:
    """(mu_wind, mu_plus) of an ionized wind, in proton masses.

    Generalizes the printed mean-mass pairs of Lopez (2017): with hydrogen
    fully ionized, heavier atoms singly ionized, and the electrons counted
    among the particles, the mean mass per particle is half the mean atomic
    mass and the mean mass per ion is the mean atomic mass itself. The rule
    reproduces Lopez's printed H/He pair (0.62, 1.3) and steam pair (3, 6).
    """
    mbar = sum(x * ELEMENT_AMU[el] for el, x in element_fractions.items())
    return mbar / 2.0, mbar


def rr_chain(
    M_p: float, F_xuv: float, R_base: float, T_wind: float, element_fractions: dict
) -> dict:
    """The radiation-recombination-limited chain at the wind base.

    Evaluates the Murray-Clay et al. (2009) analytic chain at wind
    temperature ``T_wind`` for the atomized base composition: sound speed
    and sonic radius ``R_s = G M_p / (2 c_s^2)``, the base Jeans parameter
    ``lambda_b``, the base ion density from photoionization-recombination
    balance (proportional to ``sqrt(F_xuv)``), the barometric factor
    ``exp(3/2 - lambda_b)`` to the sonic point, and the rate
    ``4 pi rho_s c_s R_s^2``.

    When the computed sonic radius falls below the base (a subcritical
    configuration), the sonic radius is floored at the base and the density
    there is the base density (the barometric factor is not applied below
    the base); the ``subcritical`` flag reports it and the caller carries
    it on the result.

    The ionizing front follows the composition: the 20 eV hydrogen front
    for winds with an atomized hydrogen fraction of one half or more, the
    33.6 eV nitrogen-like front otherwise. The composition recombination
    coefficient is the mole-fraction-weighted case B set with its
    documented temperature scaling.

    Returns a dict with ``c_s``, ``R_s``, ``R_s_calc``, ``lambda_b``,
    ``rho_base``, ``rho_s``, ``n_plus_base``, ``n_0_base``,
    ``f_plus_base``, ``mdot_rr`` [kg/s], ``subcritical``,
    ``barometric_factor``, ``mu_wind``, ``mu_plus_wind``, ``hnu0_eV``.
    """
    mu_wind, mu_plus = wind_mean_masses(element_fractions)
    c_s = math.sqrt(kb * T_wind / (mu_wind * m_p))
    r_s_calc = G * M_p / (2.0 * c_s**2)
    subcritical = r_s_calc < R_base
    r_s = max(r_s_calc, R_base)
    lambda_b = G * M_p / (R_base * c_s**2)

    x_h = element_fractions.get('H', 0.0)
    hnu0 = (HNU0_H_EV if x_h >= 0.5 else HNU_I_N_EV) * ev2joule

    # Composition-weighted case B coefficient, cm^3/s -> m^3/s.
    alpha_b = sum(x * alpha_case_b(el, T_wind) for el, x in element_fractions.items()) * 1e-6

    # Base ion density from photoionization-recombination balance with the
    # neutral density at unit optical depth over a scale height substituted,
    # so the photoionization cross section cancels:
    # n_+^2 = F G M / (h nu0 alpha_B c_s^2 R_base^2).
    n_plus_base = (
        math.sqrt(F_xuv * G * M_p / (hnu0 * alpha_b * c_s**2 * R_base**2)) if F_xuv > 0 else 0.0
    )
    rho_base = n_plus_base * mu_plus * m_p
    # Neutral base density from unit optical depth over a scale height:
    # n_0 = G M / (sigma_nu0 c_s^2 R_base^2).
    n_0_base = G * M_p / (SIGMA_NU0 * c_s**2 * R_base**2)
    f_plus = n_plus_base / (n_plus_base + n_0_base) if (n_plus_base + n_0_base) > 0 else 0.0

    if subcritical:
        baro = 1.0
        rho_s = rho_base
    else:
        baro = math.exp(1.5 - lambda_b)
        rho_s = rho_base * baro
    mdot_rr = 4.0 * math.pi * rho_s * c_s * r_s**2
    return dict(
        c_s=c_s,
        R_s=r_s,
        R_s_calc=r_s_calc,
        lambda_b=lambda_b,
        rho_base=rho_base,
        rho_s=rho_s,
        n_plus_base=n_plus_base,
        n_0_base=n_0_base,
        f_plus_base=f_plus,
        mdot_rr=mdot_rr,
        subcritical=subcritical,
        barometric_factor=baro,
        mu_wind=mu_wind,
        mu_plus_wind=mu_plus,
        hnu0_eV=hnu0 / ev2joule,
    )


def selection_mechanism(rr: dict, el_won: bool) -> str:
    """Which mechanism min(EL, RR) actually selected; diagnostic only.

    An RR win means one of two physically different things: genuine
    recombination saturation (the sqrt(F) limitation at modest base Jeans
    parameter) or barometric suppression (large ``lambda_b``: the
    isothermal wind exponentially throttled between base and sonic point,
    which has nothing to do with recombination). The split at
    ``lambda_b = 4`` is a reporting convention, stated as such. This string
    never gates anything.
    """
    if el_won:
        return 'EL-selected'
    if rr['subcritical']:
        return 'RR-selected:subcritical-floor'
    if rr['lambda_b'] >= 4.0:
        return 'RR-selected:barometric-suppression'
    return 'RR-selected:recombination-saturation'
