"""
!!! info "`escape.py`"
    Main functions to compute atmospheric escape.<br>
    Authors: Emma Postolec, Harrison Nicholls
"""

import numpy as np

from zephyrus.constants import *
from zephyrus.planets_parameters import *

########################################################### Energy-Limited escape (EL) ###########################################################


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
    escape_EL = (epsilon * np.pi * R_cubed * Fxuv) / (G * Mp * K_tide)

    return escape_EL
