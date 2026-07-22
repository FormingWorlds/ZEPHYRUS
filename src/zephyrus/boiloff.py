# Boil-off regime diagnostics for weakly bound hydrogen-helium envelopes
# Authors: Tim Lichtenberg
from __future__ import annotations

from zephyrus.constants import G, kb, m_H

# Planet radius, in units of the Bondi radius, at which boil-off ceases
# (Owen & Wu 2016). `boiloff_mass_factor` is normalised to return unity here.
BOILOFF_RADIUS_FRACTION = 0.1

# Mean molecular weight of a hydrogen-helium envelope [dimensionless]
MU_ENVELOPE = 2.35

# Restricted Jeans parameter below which escape is driven by the envelope's
# own thermal energy rather than by stellar XUV heating. Published values span
# 15 to 35 and depend on stellar type and orbital separation.
LAMBDA_CRITICAL = 20.0


def bondi_radius(Mp: float, Teq: float, mu: float = MU_ENVELOPE) -> float:
    """Radius at which the isothermal sound speed equals the escape speed.

    Parameters
    ----------
    Mp : float
        Planetary mass [kg].
    Teq : float
        Planetary equilibrium temperature [K].
    mu : float
        Mean molecular weight of the envelope [dimensionless]. Defaults to a
        hydrogen-helium mixture.

    Returns
    -------
    float
        Bondi radius [m].

    Raises
    ------
    ValueError
        If any argument is non-positive.

    Notes
    -----
    ``R_B = G Mp mu m_H / (2 kb Teq)``, the isothermal sonic radius of
    Owen & Wu (2016), Equation 1. An envelope extending beyond a tenth of this
    radius is losing mass through a thermally driven wind.
    """
    if Mp <= 0.0 or Teq <= 0.0 or mu <= 0.0:
        raise ValueError(
            f'Bondi radius needs positive mass, temperature and molecular weight; '
            f'got Mp={Mp:.4g} kg, Teq={Teq:.4g} K, mu={mu:.4g}'
        )
    return G * Mp * mu * m_H / (2.0 * kb * Teq)


def restricted_jeans(Mp: float, Rp: float, Teq: float) -> float:
    """Ratio of gravitational to thermal energy at the planetary radius.

    Parameters
    ----------
    Mp : float
        Planetary mass [kg].
    Rp : float
        Planetary radius [m].
    Teq : float
        Planetary equilibrium temperature [K].

    Returns
    -------
    float
        Restricted Jeans escape parameter [dimensionless].

    Raises
    ------
    ValueError
        If any argument is non-positive.

    Notes
    -----
    ``Lambda = G Mp m_H / (kb Teq Rp)`` of Fossati et al. (2017), Equation 10,
    evaluated for atomic hydrogen so it depends only on observable bulk
    properties. Values below :data:`LAMBDA_CRITICAL` mark the regime where the
    outflow is powered by the envelope's thermal energy and the energy-limited
    formula is no longer an upper bound on the mass-loss rate.
    """
    if Mp <= 0.0 or Rp <= 0.0 or Teq <= 0.0:
        raise ValueError(
            f'Restricted Jeans parameter needs positive mass, radius and '
            f'temperature; got Mp={Mp:.4g} kg, Rp={Rp:.4g} m, Teq={Teq:.4g} K'
        )
    return G * Mp * m_H / (kb * Teq * Rp)


def boiloff_mass_factor(
    Rp: float,
    R_bondi: float,
    kappa: float = 1.0 / BOILOFF_RADIUS_FRACTION,
    f_min: float = 1e-3,
) -> float:
    """Fraction of the envelope that survives boil-off.

    Parameters
    ----------
    Rp : float
        Planetary radius at the onset of boil-off [m].
    R_bondi : float
        Bondi radius [m], from :func:`bondi_radius`.
    kappa : float
        Reciprocal of the radius, in Bondi units, at which boil-off ceases.
        Defaults to the reciprocal of :data:`BOILOFF_RADIUS_FRACTION`.
    f_min : float
        Lower clamp on the returned fraction [dimensionless], so a grossly
        unbound envelope leaves a residue rather than exactly nothing.

    Returns
    -------
    float
        Surviving envelope mass fraction, in ``(0, 1]``.

    Raises
    ------
    ValueError
        If any argument is non-positive.

    Notes
    -----
    Owen & Wu (2016), Equation 16, with their order-unity central-concentration
    constants set to one:

    ``f = (A_i + 1) / (kappa alpha A_e + 1) -> 2 / (kappa alpha + 1)``

    for ``alpha = Rp / R_B``. The expression returns exactly unity at
    ``alpha = 1 / kappa``, so the correction switches off precisely at the
    radius where boil-off is expected to stop, and falls towards zero for
    envelopes far outside their Bondi radius. Use this exact form rather than
    the asymptotic ``0.2 / alpha`` quoted for small ``alpha``, which exceeds
    unity before that radius is reached.
    """
    if Rp <= 0.0 or R_bondi <= 0.0 or kappa <= 0.0 or f_min <= 0.0:
        raise ValueError(
            f'Boil-off mass factor needs positive radii and parameters; got '
            f'Rp={Rp:.4g} m, R_bondi={R_bondi:.4g} m, kappa={kappa:.4g}, '
            f'f_min={f_min:.4g}'
        )
    alpha = Rp / R_bondi
    return min(1.0, max(f_min, 2.0 / (kappa * alpha + 1.0)))
