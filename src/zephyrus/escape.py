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
        (Erkaev et al. 2007). It is valid for
        $\xi \equiv R_\mathrm{Hill}/R_\mathrm{XUV} > 1$, where
        $0 < K_\mathrm{tide} < 1$ and the correction enhances escape; the
        factor rises monotonically toward 1 as $\xi \to \infty$. A
        ``ValueError`` is raised for $\xi \le 1$, where the atmosphere
        reaches the Roche lobe and the energy-limited approximation no
        longer applies. If False, $K_\mathrm{tide} = 1$ (no tidal
        effects).
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
        range is $0.1 < \epsilon < 0.6$.
    Rp : float
        Planetary radius [m]. Used as a linear factor when
        ``scaling=2``.
    Rxuv : float
        Planetary radius at which the atmosphere becomes optically
        thick to XUV radiation [m]. Defined at 20 mbar in
        Baumeister et al. (2023).
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
        ``tidal_contribution`` is True and
        $\xi \equiv R_\mathrm{Hill}/R_\mathrm{XUV} \le 1$ (the atmosphere
        reaches the Roche lobe, outside the energy-limited regime).

    References
    ----------
    Based on the formulation of Lopez, Fortney & Miller (2012),
    Equations 2-4. The alternative radius scaling (``scaling=3``)
    follows Lehmer & Catling (2017), Equation 1. The tidal correction
    factor is from Erkaev et al. (2007), Equation 21.

    1. Lopez, E. D., Fortney, J. J., & Miller, N. (2012).
       How thermal evolution and mass-loss sculpt populations of
       super-Earths and sub-Neptunes. *ApJ*, 761(1), 59.
    2. Lehmer, O. R., & Catling, D. C. (2017). Rocky worlds
       limited to ~1.8 Earth radii by atmospheric escape during a
       star's extreme UV saturation. *ApJ*, 845(2), 130.
    3. Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007).
       Roche lobe effects on the atmospheric loss from "Hot Jupiters".
       *A&A*, 472(1), 329-334.
    """
    # Tidal contribution
    if tidal_contribution:
        # ksi = Rhill/Rxuv is the ratio of the periapsis Hill radius to the
        # XUV radius. K_tide = (ksi-1)^2 (2 ksi + 1) / (2 ksi^3) is non-negative
        # for all ksi > 0 with a double root at ksi = 1, so the energy-limited
        # rate (which divides by K_tide) diverges as ksi -> 1 and is only valid
        # for ksi > 1, where the atmosphere sits inside the Roche lobe.
        Rhill = a * (1 - e) * (Mp / (3 * Ms)) ** (1 / 3)
        ksi = Rhill / Rxuv
        if ksi <= 1:
            raise ValueError(
                'Tidal energy-limited escape requires the periapsis Hill '
                'radius to exceed the XUV radius '
                f'(ksi = Rhill/Rxuv > 1); got ksi = {ksi:.4g}. At ksi <= 1 the '
                'atmosphere reaches the Roche lobe and the energy-limited '
                'approximation no longer applies.'
            )
        K_tide = 1 - (3 / (2 * ksi)) + (1 / (2 * (ksi**3)))
    else:
        K_tide = 1

    # Radius
    match scaling:
        case 2:
            R_cubed = Rp * Rxuv**2
        case 3:
            R_cubed = Rxuv**3
        case _:
            raise ValueError(f'Invalid radius exponent: {scaling}')

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * np.pi * R_cubed * Fxuv) / (G * Mp * K_tide)

    return escape_EL
