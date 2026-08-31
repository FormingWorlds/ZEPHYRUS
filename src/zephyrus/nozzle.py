"""
!!! info "`nozzle.py`"
    Roche-lobe overflow: the tidally driven nozzle flow through L1.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

from zephyrus.constants import G, kb

# Provenance of the branch, all closed form, from one primary:
#
# - Rate: Jackson et al. (2017, ApJ 835, 145) Eq. (3), isothermal mass
#   transfer through the inner Lagrange point from a donor with an
#   extended atmosphere, in the lineage of Ritter (1988, A&A 202, 93)
#   rebuilt to hold at arbitrary mass ratio,
#   Mdot = e^(-1/2) rho_ph exp(-(Phi_L1 - Phi_ph)/v_th^2) v_th
#          * 2 pi v_th^2 / (Omega^2 sqrt(A (A - 1))),
#   with v_th = sqrt(kB T / mu) the isothermal sound speed and
#   Omega^2 = G (M_d + M_a) / a^3 the orbital frequency. The three factors
#   are the density at L1 (a Bernoulli integral from the photosphere, their
#   Eqs. 11 to 13, source of the e^(-1/2)), the transonic speed there, and
#   the elliptical nozzle area around L1 (their Eqs. 8 and 9).
# - Nozzle curvature: A(q) from their Eq. (10) fit,
#   A = 4 + b1 / (b2 + q^(1/3) + q^(-1/3)), b1 = 2 * 3^(2/3),
#   b2 = b1/4 - 2, symmetric under q to 1/q, accurate to 0.3% for all
#   mass ratios (their Figure 2), which is the specific advance over the
#   Ritter (1988) fits that hold only for donor-accretor ratios of roughly
#   0.05 to 25 and fail at planetary values.
# - Lobe radius: the Eggleton (1983, ApJ 268, 368) fit as printed in their
#   Section 2.1, r_R = a 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3))),
#   accurate to 1% for all q.
# - Potentials: their Eq. (14) volume-averaged Roche potential, evaluated
#   at the lobe radius for Phi_L1 and at the photospheric radius for
#   Phi_ph. The expansion converges inside the lobe and not outside; at
#   and beyond lobe contact the exponent is clamped at zero, which is the
#   paper's own lobe-filling case (their Figure 5 solid curves), and the
#   clamped value is a boundary value rather than a trusted rate.
# - Applicability: the overflow description holds where the isothermal
#   sonic radius R_sonic = G M_d / (2 v_th^2) lies at or beyond the L1
#   distance, so that no spherical transonic wind fits inside the lobe and
#   the L1 nozzle is the flow's constriction (their Section 4 and
#   Figure 9). Inward of that the gas chokes at its own sonic surface
#   first and the wind branches of this package are the right
#   description; the candidate rate is still reported there, but the
#   dispatcher does not let it compete. Without this criterion the nozzle
#   area, which grows as the cube of the separation, hands a loosely
#   bound envelope an unbounded rate at separations where the planet is
#   nowhere near its lobe.
# - Stated limitations carried from the primary: the flow is isothermal
#   (their Section 3.1 names the neglected thermal structure), the orbit
#   circular and the rotation synchronous (an eccentric caller is
#   evaluated at periapsis, our convention rather than theirs), and the
#   rate can overestimate the transfer where the escaping gas keeps its
#   orbital angular momentum and disk-stellar torque balance regulates the
#   flow instead (their Eq. 24 and Figure 6); the torque-balance rate
#   needs the stellar tidal dissipation and is not computed here.
#
# The temperature is the model's dominant uncertainty by the authors' own
# statement; which temperature enters v_th is the caller's
# `nozzle_temperature` setting, resolved in the dispatcher.

_B1 = 2.0 * 3.0 ** (2.0 / 3.0)
_B2 = _B1 / 4.0 - 2.0


def curvature_a(q: float) -> float:
    """Dimensionless L1 curvature A(q) of Jackson et al. (2017) Eq. (10).

    ``q`` is the donor-over-accretor mass ratio. The fit is symmetric under
    ``q`` to ``1/q``, equals 8 at equal masses, and tends to 4 at extreme
    ratios; it reproduces the numerical root of their Eq. (7) to 0.3% for
    all ``q`` (their Figure 2).
    """
    if q <= 0.0 or not math.isfinite(q):
        raise ValueError(f'q must be a positive finite mass ratio, got {q!r}')
    return 4.0 + _B1 / (_B2 + q ** (1.0 / 3.0) + q ** (-1.0 / 3.0))


def eggleton_lobe_radius(q: float, separation: float) -> float:
    """Volume-equivalent Roche lobe radius of the donor, in m.

    The Eggleton (1983, ApJ 268, 368) fit as printed in Jackson et al.
    (2017) Section 2.1, accurate to 1% for all mass ratios ``q`` (donor
    over accretor); ``separation`` is the orbital separation [m].
    """
    q13 = q ** (1.0 / 3.0)
    return separation * 0.49 * q13**2 / (0.6 * q13**2 + math.log(1.0 + q13))


def volume_averaged_potential(r_v: float, M_d: float, M_a: float, separation: float) -> float:
    """Volume-averaged Roche potential at volume-equivalent radius ``r_v``.

    Jackson et al. (2017) Eq. (14): the potential of the equipotential
    surface enclosing the same volume as a sphere of radius ``r_v`` around
    the donor ``M_d``, with accretor ``M_a`` at ``separation`` [m], in
    J/kg. The expansion converges inside the donor's Roche lobe and not
    outside it, and disagrees with a direct numerical evaluation by a few
    percent as ``r_v`` approaches the lobe radius, which the primary
    quantifies as about a factor of two in the final rate.
    """
    m_t = M_d + M_a
    x = r_v / separation
    bracket = (
        1.0
        + (m_t / M_d) * x**3 / 3.0
        + (4.0 / 45.0) * ((m_t**2 + 9.0 * M_a**2 + 3.0 * M_a * m_t) / M_d**2) * x**6
    )
    return (
        -(G * M_a / separation + G * M_a**2 / (2.0 * separation * m_t))
        - (G * M_d / r_v) * bracket
    )


def nozzle_candidate(
    M_p: float,
    M_star: float,
    a: float,
    e: float,
    rho_ph: float,
    r_ph: float,
    T: float,
    mu_kg: float,
) -> tuple[float, dict]:
    """Roche-lobe overflow rate through the L1 nozzle, Jackson et al. (2017) Eq. (3).

    Parameters
    ----------
    M_p, M_star : float
        Donor (planet) and accretor (star) masses [kg].
    a, e : float
        Semi-major axis [m] and eccentricity. The geometry is evaluated at
        the periapsis separation ``a (1 - e)``, matching the Roche screen's
        periapsis Hill radius; the primary treats a circular, synchronously
        rotating donor, so the periapsis evaluation is this module's
        convention and an upper bound on the instantaneous rate elsewhere
        on the orbit.
    rho_ph, r_ph : float
        Density [kg m^-3] and radius [m] of the launch level. The profile
        radius stands in for the volume-equivalent photospheric radius
        without the primary's Appendix distortion conversion, a
        few-percent radius convention worth about 1.6x in rate per percent
        near lobe contact and nothing for a donor well inside its lobe.
        The Bernoulli structure makes rho_ph exp(Phi_ph / v_th^2)
        level-invariant along an isothermal column, so the level choice
        itself largely cancels; the temperature is the sensitivity.
    T, mu_kg : float
        Temperature [K] and mean particle mass [kg] evaluating the
        isothermal sound speed and the exponential barrier.

    Returns
    -------
    (rate, detail)
        The candidate rate [kg/s] and a detail dict: the sound speed, the
        mass ratio and curvature, the lobe radius, both potentials and
        their difference, the applied exponent, the nozzle area, and
        ``saturated`` (the photospheric potential reached the L1 value, so
        the exponent was clamped at zero and the rate is the lobe-filling
        boundary value).
    """
    a_peri = a * (1.0 - e)
    v_th = math.sqrt(kb * T / mu_kg)
    omega2 = G * (M_p + M_star) / a_peri**3
    q = M_p / M_star
    a_curv = curvature_a(q)
    r_lobe = eggleton_lobe_radius(q, a_peri)
    phi_l1 = volume_averaged_potential(r_lobe, M_p, M_star, a_peri)
    phi_ph = volume_averaged_potential(r_ph, M_p, M_star, a_peri)
    delta_phi = phi_l1 - phi_ph
    exponent = -delta_phi / v_th**2
    # The saturation test is geometric, not potential-ordered: outside the
    # lobe the Eq. (14) expansion diverges downward, so a level beyond
    # r_lobe reports a spuriously deep Phi_ph and a large positive barrier
    # where the physical barrier is gone. At or beyond contact the
    # exponential is clamped at 1 and the rate is the lobe-filling
    # boundary value, a lower bound on the transfer (the density at the
    # lobe itself exceeds the launch level's).
    saturated = r_ph >= r_lobe or exponent >= 0.0
    if saturated:
        exponent = 0.0
    area = 2.0 * math.pi * v_th**2 / (omega2 * math.sqrt(a_curv * (a_curv - 1.0)))
    rate = rho_ph * math.exp(-0.5 + exponent) * v_th * area
    # The Figure 9 crossover quantity: the overflow description applies
    # where this sonic radius reaches the L1 distance (see module notes).
    r_sonic = G * M_p / (2.0 * v_th**2)
    return rate, dict(
        T_K=T,
        mu_kg=mu_kg,
        v_th=v_th,
        R_sonic=r_sonic,
        q=q,
        A=a_curv,
        a_periapsis=a_peri,
        r_lobe=r_lobe,
        phi_L1=phi_l1,
        phi_ph=phi_ph,
        delta_phi=delta_phi,
        exponent_applied=exponent,
        area_m2=area,
        saturated=saturated,
    )
