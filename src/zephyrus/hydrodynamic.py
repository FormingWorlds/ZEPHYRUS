"""
!!! info "`hydrodynamic.py`"
    Hydrodynamic escape: energy-limited and radiation-recombination-limited
    rates, and the selection between them.<br>
    Authors: Emma Postolec, Harrison Nicholls, Malina Ovesen, Mara Attia
"""

from __future__ import annotations

import math

from zephyrus.atomic_data import (
    HNU0_H_EV,
    HNU_I_N_EV,
    SIGMA_NU0_H,
    SIGMA_NU_N,
    alpha_case_b,
)
from zephyrus.composition import ELEMENT_AMU
from zephyrus.constants import G, ev2joule, kb, m_p

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
#   genuine recombination saturation, and barometric suppression at large
#   lambda_b, where the label "recombination limited" would be a category
#   error. The flux scaling does not separate them, since the base ion
#   density follows sqrt(F_XUV) at every lambda_b in this chain; the
#   barometric factor does, and it is reported beside the rate. The
#   selection diagnostic names which candidate won, not why. Where the two
#   candidates cross in flux is sensitive to the wind temperature: the RR
#   chain carries it through the sound speed, the barometric exponent, and
#   the recombination coefficient, so a thermostat-driven wind temperature
#   can move the EL/RR crossover by an order of magnitude against the
#   canonical fixed 1e4 K evaluation.
# - Efficiency: fixed, or the Caldiroli et al. (2022, A&A 663, A122,
#   Appendix A.1) fit, defined against their R_p^3 geometry and therefore
#   converted by (R_p/R_XUV)^2 before use in the Erkaev form.
#
# EL_escape is the released standalone entry point for the energy-limited
# rate (scaling selection, tidal branch, and input validation in one
# self-contained function); el_rate is the bare kernel the regime dispatch
# assembles with its own tidal factor. The two are kept as separate code
# paths on purpose, so the cross-implementation test between them guards
# the scaling and tidal plumbing. zephyrus.escape re-exports EL_escape for
# compatibility with the released import path.

RHO_UNIT_CGS = 1e-3  # kg m^-3 -> g cm^-3
FLUX_UNIT_CGS = 1e3  # W m^-2 -> erg s^-1 cm^-2


def hill_radius_periapsis(M_p: float, M_star: float, a: float, e: float) -> float:
    """Periapsis Hill radius a (1 - e) (M_p / 3 M_star)^(1/3), in m."""
    return a * (1.0 - e) * (M_p / (3.0 * M_star)) ** (1.0 / 3.0)


def k_tide(xi: float) -> float:
    """Erkaev et al. (2007) Eq. (17) tidal factor, for xi > 1.

    The factor is ``(xi - 1)^2 (2 xi + 1) / (2 xi^3)``, which has a double
    root at ``xi = 1`` and rises toward 1 as ``xi`` grows. The energy-limited
    rate divides by it, so the rate diverges as the atmosphere approaches its
    Roche lobe: the factor is 1.5e-6 at xi = 1.001 and 1.2e-2 at xi = 1.1,
    inflating the rate 6.7e5-fold and 83-fold. At and below the root the
    polynomial turns back upward and returns values above 1, which would
    reduce the rate rather than raise it, so the domain is enforced rather
    than extrapolated: a caller at xi <= 1 has a planet filling its lobe and
    needs the overflow machinery, not this factor.
    """
    if not xi > 1.0:
        raise ValueError(
            f'k_tide is defined for xi > 1 and has a double root at 1, got {xi!r}'
        )
    return 1.0 - 3.0 / (2.0 * xi) + 1.0 / (2.0 * xi**3)


def el_rate(eps: float, F_xuv: float, R_p: float, R_xuv: float, M_p: float, K: float) -> float:
    """Energy-limited rate eps pi F R_p R_xuv^2 / (G M_p K), in kg/s."""
    return eps * math.pi * F_xuv * R_p * R_xuv**2 / (G * M_p * K)


def EL_escape(
    tidal_contribution: bool,
    a: float,
    e: float,
    Mp: float,
    Ms: float,
    epsilon: float,
    Rp: float,
    Rxuv: float,
    Fxuv: float,
    scaling: int = 2,
):
    r"""
    Compute the mass-loss rate for Energy-Limited (EL) atmospheric escape.

    The mass-loss rate is given by

    $$
    \dot{M}_\mathrm{EL} = \frac{\epsilon\,\pi\,R^3\,F_\mathrm{XUV}}
                               {G\,M_p\,K_\mathrm{tide}}
    $$

    where $R^3$ is either $R_p R_\mathrm{XUV}^2$ or $R_\mathrm{XUV}^3$
    depending on ``scaling``, and $K_\mathrm{tide}$ is the tidal
    correction factor of Erkaev et al. (2007) when ``tidal_contribution``
    is True, else 1.

    Parameters
    ----------
    tidal_contribution : bool
        If True, include the tidal correction factor $K_\mathrm{tide}$
        (Erkaev et al. 2007). Its argument is
        $\xi \equiv R_\mathrm{Hill}/R$, where $R$ is the radius that
        appears linearly in the $R^3$ term selected by ``scaling``: $R_p$
        for ``scaling=2`` (the convention of Erkaev et al. 2007, whose
        own $\xi$ is the Roche-lobe distance over the planetary radius)
        and $R_\mathrm{XUV}$ for ``scaling=3`` (the single-radius form,
        where $R_\mathrm{XUV}$ is the only radius in the problem). The
        factor is valid for $\xi > 1$, where $0 < K_\mathrm{tide} < 1$
        and the correction enhances escape; it rises monotonically
        toward 1 as $\xi \to \infty$. A ``ValueError`` is raised for
        $\xi \le 1$, where the atmosphere reaches the Roche lobe and the
        energy-limited approximation no longer applies. If False,
        $K_\mathrm{tide} = 1$ (no tidal effects).
    a : float
        Planetary semi-major axis [m]. Only used when
        ``tidal_contribution`` is True.
    e : float
        Orbital eccentricity (dimensionless). Only used when
        ``tidal_contribution`` is True.
    Mp : float
        Planetary mass [kg].
    Ms : float
        Stellar mass [kg]. Only used when
        ``tidal_contribution`` is True.
    epsilon : float
        Escape efficiency factor (dimensionless). Typical literature
        range is $0.1 < \epsilon < 0.6$, but hydrodynamic simulations
        find the effective efficiency falls far below that band for
        strongly bound planets: above a threshold gravitational
        potential, $\log_{10}(G M_p K_\mathrm{tide}/R_p) \approx 12.9$
        to $13.2$ in cgs units (erg g$^{-1}$), it drops to of order
        $10^{-2}$ for compact hot Jupiters (Caldiroli et al. 2022).
    Rp : float
        Planetary radius [m]. Used as a linear factor when
        ``scaling=2``.
    Rxuv : float
        Planetary radius at which the atmosphere becomes optically
        thick to XUV radiation [m]. In PROTEUS this level is placed at
        a fixed pressure, by default 20 mbar following Baumeister et
        al. (2023); that is an optical-photosphere-type level, distinct
        from the roughly nanobar level where the XUV heating is
        actually deposited and the wind is launched (Lopez 2017,
        $P_\mathrm{base} = \mu m_\mathrm{H} g / \sigma_{\nu_0}$).
    Fxuv : float
        XUV flux received by the planet from the host star, in
        W m$^{-2}$.
    scaling : int, optional
        Planet radius scaling exponent. ``2`` (default) uses
        $R_p R_\mathrm{XUV}^2$; ``3`` uses $R_\mathrm{XUV}^3$. Any other
        value raises ``ValueError``.

    Returns
    -------
    escape_EL : float
        Mass-loss rate for energy-limited escape, in kg s$^{-1}$.

    Raises
    ------
    ValueError
        If ``scaling`` is not ``2`` or ``3``, or if
        ``tidal_contribution`` is True and $\xi \le 1$ (the atmosphere
        reaches the Roche lobe, outside the energy-limited regime),
        with $\xi$ built on the radius selected by ``scaling``.

    References
    ----------
    The default radius scaling (``scaling=2``, ``Rp * Rxuv**2``) is the
    energy-limited XUV cross-section form of Watson et al. (1981) and
    Lammer et al. (2003), Equation 6, written as a mass-loss rate by
    Erkaev et al. (2007), Equation 21. The alternative radius scaling
    (``scaling=3``, ``Rxuv**3``) is the single-radius simplification of
    Lopez, Fortney & Miller (2012), Equation 2, Lopez & Fortney (2013),
    Equation 1, and Lehmer & Catling (2017), Equation 1. The tidal
    reduction factor ``K_tide`` is Erkaev et al. (2007), Equation 17.

    1. Watson, A. J., Donahue, T. M., & Walker, J. C. G. (1981).
       The dynamics of a rapidly escaping atmosphere: applications to
       the evolution of Earth and Venus. *Icarus*, 48(2), 150-166.
    2. Lammer, H., Selsis, F., Ribas, I., et al. (2003). Atmospheric
       loss of exoplanets resulting from stellar X-ray and
       extreme-ultraviolet heating. *ApJ*, 598(2), L121-L124.
    3. Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007).
       Roche lobe effects on the atmospheric loss from "Hot Jupiters".
       *A&A*, 472(1), 329-334.
    4. Lopez, E. D., Fortney, J. J., & Miller, N. (2012).
       How thermal evolution and mass-loss sculpt populations of
       super-Earths and sub-Neptunes. *ApJ*, 761(1), 59.
    5. Lopez, E. D., & Fortney, J. J. (2013). The role of core mass
       in controlling evaporation: the Kepler radius distribution and
       the Kepler-36 density dichotomy. *ApJ*, 776(1), 2.
    6. Lehmer, O. R., & Catling, D. C. (2017). Rocky worlds
       limited to ~1.8 Earth radii by atmospheric escape during a
       star's extreme UV saturation. *ApJ*, 845(2), 130.
    7. Lopez, E. D. (2017). Born dry in the photoevaporation desert:
       Kepler's ultra-short-period planets formed water-poor.
       *MNRAS*, 472(1), 245-253.
    8. Baumeister, P., Tosi, N., Brachmann, C., Grenfell, J. L., &
       Noack, L. (2023). Redox state and interior structure control on
       the long-term habitability of stagnant-lid planets.
       *A&A*, 675, A122.
    9. Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I.,
       & Rauscher, E. (2022). Irradiation-driven escape of primordial
       planetary atmospheres II. Evaporation efficiency of sub-Neptunes
       through hot Jupiters. *A&A*, 663, A122.
    """
    # Radius term, and the radius the tidal factor is measured from: the
    # one that appears linearly in R^3, since that is the radius the
    # potential barrier in the denominator refers to.
    match scaling:
        case 2:
            R_cubed = Rp * Rxuv**2
            R_tide = Rp
        case 3:
            R_cubed = Rxuv**3
            R_tide = Rxuv
        case _:
            raise ValueError(f'Invalid radius exponent: {scaling}')

    # Tidal contribution
    if tidal_contribution:
        # ksi is the ratio of the periapsis Hill radius to the radius the
        # scaling selects. K_tide = (ksi-1)^2 (2 ksi + 1) / (2 ksi^3) is
        # non-negative for all ksi > 0 with a double root at ksi = 1, so the
        # energy-limited rate (which divides by K_tide) diverges as ksi -> 1
        # and is only valid for ksi > 1, where the atmosphere sits inside the
        # Roche lobe.
        Rhill = a * (1 - e) * (Mp / (3 * Ms)) ** (1 / 3)
        ksi = Rhill / R_tide
        if ksi <= 1:
            raise ValueError(
                'Tidal energy-limited escape requires the periapsis Hill '
                'radius to exceed the escape-level radius '
                f'(ksi = Rhill/R > 1); got ksi = {ksi:.4g}. At ksi <= 1 the '
                'atmosphere reaches the Roche lobe and the energy-limited '
                'approximation no longer applies.'
            )
        K_tide = 1 - (3 / (2 * ksi)) + (1 / (2 * (ksi**3)))
    else:
        K_tide = 1

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * math.pi * R_cubed * Fxuv) / (G * Mp * K_tide)

    return escape_EL


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
    """(mu_wind, mu_plus) of an ionized wind, in atomic mass units.

    Generalizes the mean-mass pairs of Lopez (2017): with hydrogen fully
    ionized, heavier atoms singly ionized, and the electrons counted among
    the particles, the mean mass per particle is half the mean atomic mass
    and the mean mass per ion is the mean atomic mass itself.

    Two conventions to keep straight. The returned values are in atomic mass
    units, and the call sites multiply by the proton mass where Lopez writes
    the hydrogen atom mass; the three candidate units span 0.36 percent on
    the sound speed, and the proton mass sits 0.04 percent from Lopez's own.
    And Lopez's printed pairs are not both reachable: the rule makes the
    per-ion mass exactly twice the per-particle mass, so the printed steam
    pair (3, 6) is recovered while the printed H/He pair (0.62, 1.3) is
    internally inconsistent by 4.6 percent, 1.3 halving to 0.65. The rule
    follows the per-ion value and the tests pin that reading.
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
    33.6 eV nitrogen-like front otherwise, and the photoionization cross
    section follows the same front rather than staying at hydrogen's. The
    composition recombination coefficient is the mole-fraction-weighted case
    B set with its documented temperature scaling.

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

    # The photon energy and the cross section belong to one front and must
    # be taken from the same one. Taking the energy from the composition and
    # the cross section from hydrogen put a nitrogen-like wind on a section
    # 5.3 times too small, which raised its neutral base density by that
    # factor and understated the reported base ionization fraction.
    x_h = element_fractions.get('H', 0.0)
    hydrogen_front = x_h >= 0.5
    hnu0 = (HNU0_H_EV if hydrogen_front else HNU_I_N_EV) * ev2joule
    sigma_nu0 = (SIGMA_NU0_H if hydrogen_front else SIGMA_NU_N) * 1e-4  # cm^2 -> m^2

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
    n_0_base = G * M_p / (sigma_nu0 * c_s**2 * R_base**2)
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
    """Which candidate min(EL, RR) selected; diagnostic only.

    Three outcomes: the energy-limited rate won, the
    recombination-limited rate won, or it won with the sonic radius
    floored at the wind base (the subcritical configuration of
    :func:`rr_chain`, where the returned value is a floored one rather
    than a transonic wind).

    Why an RR win came out small is a separate question, and this string
    does not answer it. The quantity that does is the barometric factor
    ``exp(3/2 - lambda_b)`` returned beside the rate: near 1 the
    sonic-point density is the base density and the recombination-limited
    base ionization sets the rate, while several decades below 1 the rate
    is small mostly because the isothermal wind cannot carry material
    from the base to the sonic point, which has nothing to do with
    recombination. The flux scaling cannot separate the two, because the
    base ion density follows sqrt(F_XUV) at every ``lambda_b`` here.
    """
    if el_won:
        return 'EL-selected'
    if rr['subcritical']:
        return 'RR-selected:subcritical-floor'
    return 'RR-selected'
