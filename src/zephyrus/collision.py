"""
!!! info "`collision.py`"
    Fractional atmospheric mass loss of the target planet in a giant impact.<br>
    Author(s): Anna Grace Ulses
"""

from __future__ import annotations

import numpy as np

from zephyrus.constants import G


def mass_loss(
    v_c: float,
    M_i: float,
    M_t: float,
    rho_i: float,
    rho_t: float,
    R_i: float,
    R_t: float,
    b: float,
) -> float:
    r"""Fractional atmospheric mass loss of the target in a giant impact.

    Implements the scaling law of Kegerreis et al. (2020), their Eqn. 1:

    $X \approx 0.64 \left[ \left(\frac{v_c}{v_{esc}}\right)^2
    \left(\frac{M_i}{M_{tot}}\right)^{1/2}
    \left(\frac{\rho_i}{\rho_t}\right)^{1/2} f_M(b) \right]^{0.65}$

    capped at 1 for total erosion, where subscript (i) is the impactor,
    (t) the target, and $M_{tot} = M_i + M_t$. The mutual escape speed is
    $v_{esc} = \sqrt{2 G (M_t + M_i) / (R_t + R_i)}$, and $f_M(b)$ is the
    fractional interacting mass of their Eqn. B1, built from the
    density-weighted spherical caps of common height
    $d = (R_t + R_i)(1 - b)$. The common-height caps are a linearised
    bookkeeping, so outside the fitted geometry (a much denser, much
    smaller impactor near head-on) the raw $f_M$ can leave $[0, 1]$, and
    is clamped to it here, and can vary non-monotonically with $b$; at
    equal bulk densities $f_M$ reduces exactly to the interacting volume
    $f_V$ of their Eqn. B2.

    Conventions the caller must honour (Kegerreis et al. 2020, Sect. 2):
    $v_c$ is the speed at first contact, not at infinity; the masses and
    radii exclude any atmosphere, with the radii taken at its base; and
    $b \equiv \sin\beta$ for impact angle $\beta$ (0 head-on, 1 grazing).

    The fit is constrained for target masses of roughly 0.3 to 3 Earth
    masses, impactors down to about 0.05 Earth masses, bulk densities of
    about half to double Earth's, speeds of 1 to 3 $v_{esc}$, any angle,
    and thin atmospheres of order 1 percent of the planet mass. The
    median deviation of the simulations from the law is 9 percent,
    rising to about 20 percent for slow, head-on impacts.

    Parameters
    ----------
    v_c : float
        Collision speed at first contact between impactor and target [m/s].
    M_i : float
        Mass of the impactor, excluding any atmosphere [kg].
    M_t : float
        Mass of the target, excluding any atmosphere [kg].
    rho_i : float
        Bulk density of the impactor, excluding any atmosphere [kg/m^3].
    rho_t : float
        Bulk density of the target, excluding any atmosphere [kg/m^3].
    R_i : float
        Radius of the impactor, at the base of any atmosphere [m].
    R_t : float
        Radius of the target, at the base of any atmosphere [m].
    b : float
        Dimensionless impact parameter, the sine of the impact angle,
        in [0, 1]: 0 is head-on, 1 is fully grazing.

    Returns
    -------
    float
        Fractional mass loss of the target body's atmosphere, in [0, 1].

    Raises
    ------
    ValueError
        If ``b`` lies outside [0, 1], if any mass, radius, or density is
        not strictly positive and finite, or if ``v_c`` is negative or
        not finite. Inputs are scalar; arrays are not supported.

    References
    ----------
    1. Kegerreis J.A., Eke V.R., Catling D.C., Massey R.J., Teodoro
       L.F.A., Zahnle K.J. (2020). Atmospheric Erosion by Giant Impacts
       onto Terrestrial Planets: A Scaling Law for any Speed, Angle,
       Mass, and Density. ApJL 901, L31. doi:10.3847/2041-8213/abb5fb
    """
    if not 0.0 <= b <= 1.0:
        raise ValueError(f'Impact parameter b must be in [0, 1], got {b!r}')
    for name, value in (
        ('M_i', M_i),
        ('M_t', M_t),
        ('rho_i', rho_i),
        ('rho_t', rho_t),
        ('R_i', R_i),
        ('R_t', R_t),
    ):
        if not (value > 0.0 and np.isfinite(value)):
            raise ValueError(f'{name} must be strictly positive and finite, got {value!r}')
    if not (v_c >= 0.0 and np.isfinite(v_c)):
        raise ValueError(f'Collision speed v_c must be non-negative and finite, got {v_c!r}')

    # Mutual escape speed of the pair at contact
    v_esc = np.sqrt((2.0 * G * (M_t + M_i)) / (R_t + R_i))

    # Fractional interacting mass f_M (Kegerreis et al. 2020, Eqn. B1):
    # density-weighted spherical caps of common height d, clamped to [0, 1]
    # because the linearised caps can leave the interval outside the fitted
    # geometry. At equal bulk densities this reduces exactly to the
    # interacting volume f_V of their Eqn. B2.
    d = (R_t + R_i) * (1.0 - b)
    v_t_cap = np.pi / 3.0 * d**2 * (3.0 * R_t - d)
    v_i_cap = np.pi / 3.0 * d**2 * (3.0 * R_i - d)
    v_t_full = 4.0 / 3.0 * np.pi * R_t**3
    v_i_full = 4.0 / 3.0 * np.pi * R_i**3
    f_m = (rho_t * v_t_cap + rho_i * v_i_cap) / (rho_t * v_t_full + rho_i * v_i_full)
    f_m = min(max(f_m, 0.0), 1.0)

    m_tot = M_i + M_t
    bracket = (v_c / v_esc) ** 2 * (M_i / m_tot) ** 0.5 * (rho_i / rho_t) ** 0.5 * f_m

    # Fractional atmospheric mass loss of the target, capped at 1 for
    # total erosion (Kegerreis et al. 2020, Eqn. 1)
    x_loss = 0.64 * bracket**0.65
    return min(max(x_loss, 0.0), 1.0)
