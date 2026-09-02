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
# - L1 distance: the small-mass-ratio expansion of the L1 root, checked
#   against the exact stationary point of the corotating axial potential
#   (see ``l1_distance``). The Hill radius is its leading order and sits
#   0.5% outside it at planetary mass ratios, 5% at q = 1e-2.
# - Applicability: the overflow description holds where the isothermal
#   sonic radius R_sonic = G M_d / (2 v_th^2) lies at or beyond the L1
#   distance, so that no spherical transonic wind fits inside the lobe and
#   the L1 nozzle is the flow's constriction. Where the sonic radius lies
#   inside the L1 distance the gas chokes at its own sonic surface first
#   and the wind branches of this package are the right description.
#   Without this criterion the nozzle area, which grows as the cube of the
#   separation, hands a loosely bound envelope an unbounded rate at
#   separations where the planet is nowhere near its lobe. The primary
#   draws the comparison qualitatively, in their Section 4 and Figure 9,
#   to ask which of the two pictures a planet belongs in; making it a gate
#   is this module's sharpening of it and not a rule they state.
# - Orbit average: the returned rate is the time average over the orbit,
#   Kepler-weighted through the eccentric anomaly, and the detail dict
#   also carries it duty-cycled over the applicable arc, which is what a
#   dispatcher competes. A secular caller integrates over many orbital
#   periods and needs the mass carried per unit time. Each phase is
#   evaluated with the circular formula at its own separation; the primary
#   has no eccentric treatment, so the quasi-static evaluation and the
#   duty cycle are both ours. At e = 0 the average is the instantaneous
#   rate exactly.
# - Stated limitations carried from the primary: the flow is isothermal
#   (their Section 3.1 names the neglected thermal structure), the orbit
#   circular and the rotation synchronous (an eccentric caller is averaged
#   over its orbit as above, our convention rather than theirs, and the
#   rotation is synchronous at no single phase of such an orbit), and the
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


def _positive(name: str, value: float) -> float:
    """Return ``value``, or raise if it is not a positive finite number.

    The module is a public surface with its own reference page, so a caller
    reaching it directly gets a message rather than a negative mass-loss
    rate, a complex cube root, or a bare division by zero.
    """
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{name} must be positive and finite, got {value!r}')
    return value


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
    if q <= 0.0 or not math.isfinite(q):
        raise ValueError(f'q must be a positive finite mass ratio, got {q!r}')
    _positive('separation', separation)
    q13 = q ** (1.0 / 3.0)
    return separation * 0.49 * q13**2 / (0.6 * q13**2 + math.log(1.0 + q13))


def l1_distance(q: float, separation: float) -> float:
    """Distance from the donor's center to the inner Lagrange point, in m.

    The small-mass-ratio expansion of the L1 root, ``x_L1/a = eps -
    eps^2/3 - eps^3/9`` with ``eps = (q/3)^(1/3)``, which is the Hill
    radius at leading order and falls inside it beyond that. Checked
    against the exact stationary point of the corotating axial potential:
    the relative error is 8e-7 at ``q = 1e-5`` and 7e-4 at ``q = 1e-2``,
    against 0.5% and 5.3% for the Hill radius itself. ``separation`` is
    the orbital separation [m].
    """
    if q <= 0.0 or not math.isfinite(q):
        raise ValueError(f'q must be a positive finite mass ratio, got {q!r}')
    _positive('separation', separation)
    eps = (q / 3.0) ** (1.0 / 3.0)
    return separation * (eps - eps**2 / 3.0 - eps**3 / 9.0)


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
    _positive('r_v', r_v)
    _positive('M_d', M_d)
    _positive('M_a', M_a)
    _positive('separation', separation)
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


def isothermal_column_density(
    rho_ref: float, r_ref: float, r: float, M_p: float, v_th: float
) -> float:
    """Density at ``r`` on an isothermal hydrostatic column through ``r_ref``.

    ``rho(r) = rho_ref exp[(G M_p / v_th^2)(1/r - 1/r_ref)]``, the closed
    solution of hydrostatic balance in a point-mass potential at constant
    sound speed. This is the column the Bernoulli argument behind the
    launch-level convention assumes: along it the product
    ``rho exp(Phi / v_th^2)`` is constant, so the nozzle rate does not
    depend on which level is called the launch level. It is a device for
    placing that level consistently with the sound speed evaluating the
    barrier, not a claim about the structure below the anchor, which for a
    wind anchor is far hotter than the atmosphere really is there.
    """
    _positive('r', r)
    _positive('r_ref', r_ref)
    _positive('M_p', M_p)
    _positive('v_th', v_th)
    if not math.isfinite(rho_ref) or rho_ref < 0.0:
        raise ValueError(f'rho_ref must be non-negative and finite, got {rho_ref!r}')
    return rho_ref * math.exp((G * M_p / v_th**2) * (1.0 / r - 1.0 / r_ref))


def _phase_state(
    M_p: float,
    M_star: float,
    sep: float,
    rho_ph: float,
    r_ph: float,
    v_th: float,
    q: float,
    a_curv: float,
) -> dict:
    """Rate, geometry, and heat demand at one orbital separation."""
    omega2 = G * (M_p + M_star) / sep**3
    r_lobe = eggleton_lobe_radius(q, sep)
    phi_l1 = volume_averaged_potential(r_lobe, M_p, M_star, sep)
    phi_ph = volume_averaged_potential(r_ph, M_p, M_star, sep)
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
    # Heat the isothermal flow demands per unit mass, which is what the
    # radiation field has to supply for the uncapped model to hold. For a
    # steady flow dh + d(v^2/2) + dPhi = dq, and an isothermal ideal gas
    # has dh = 0, so integrating from a launch level at rest to the
    # transonic point at L1 gives the barrier the rate actually applied
    # plus v_th^2/2. Built from the applied exponent, so it is the
    # clamped barrier at saturation and never the divergent one; the
    # acceleration term survives there and the barrier does not.
    heat = max(-exponent, 0.0) * v_th**2 + 0.5 * v_th**2
    return dict(
        separation=sep,
        r_lobe=r_lobe,
        phi_L1=phi_l1,
        phi_ph=phi_ph,
        delta_phi=delta_phi,
        exponent_applied=exponent,
        area_m2=area,
        saturated=saturated,
        rate=rate,
        power=rate * heat,
        R_L1=l1_distance(q, sep),
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
    n_phase: int = 64,
) -> tuple[float, dict]:
    """Orbit-averaged Roche-lobe overflow rate through the L1 nozzle.

    The rate is Jackson et al. (2017) Eq. (3) at each orbital separation,
    averaged in time over the orbit and duty-cycled over the arc where the
    overflow description applies.

    Parameters
    ----------
    M_p, M_star : float
        Donor (planet) and accretor (star) masses [kg].
    a, e : float
        Semi-major axis [m] and eccentricity.
    rho_ph, r_ph : float
        Density [kg m^-3] and radius [m] of the launch level. The profile
        radius stands in for the volume-equivalent photospheric radius
        without the primary's Appendix distortion conversion, a
        few-percent radius convention worth about 1.6x in rate per percent
        near lobe contact and nothing for a donor well inside its lobe.
        The Bernoulli structure makes rho_ph exp(Phi_ph / v_th^2)
        level-invariant along an isothermal column, so the level choice
        largely cancels there; on a non-isothermal column it does not, and
        the temperature is the leading sensitivity either way.
    T, mu_kg : float
        Temperature [K] and mean particle mass [kg] evaluating the
        isothermal sound speed and the exponential barrier.
    n_phase : int
        Midpoint nodes in eccentric anomaly for the orbit average. The
        quadrature converges to machine precision well below the default:
        measured relative change 4.3e-7 from 16 to 32 nodes and 3e-14 from
        32 to 64 at e = 0.5 on two states. At ``e = 0`` every node holds
        the same value and the average is the instantaneous rate exactly.

    Returns
    -------
    (rate, detail)
        The orbit-averaged Eq. (3) rate [kg/s], unguarded, so that the
        closed form stays directly comparable with the primary's own
        published rates, and a detail dict. The rate a caller should
        compete is ``detail['rate_applicable_kg_s']``, the same average
        duty-cycled over the arc where the overflow description applies.
        In the detail dict, Phase-independent
        entries are the sound speed, the sonic radius, the mass ratio, and
        the curvature; the geometry entries (lobe radius, both potentials,
        the applied exponent, the nozzle area, ``saturated``) are reported
        at periapsis, which is the tightest geometry of the orbit; and the
        orbit entries are the applicable and saturated orbit fractions,
        the periapsis and apoapsis rates, and the averaged lift power,
        which is reported both duty-cycled (pairing with the rate a caller
        competes) and over the full orbit (pairing with the returned
        unguarded rate).
    """
    if n_phase < 1:
        raise ValueError(f'n_phase must be at least 1, got {n_phase!r}')
    for name, value in (
        ('M_p', M_p),
        ('M_star', M_star),
        ('a', a),
        ('r_ph', r_ph),
        ('T', T),
        ('mu_kg', mu_kg),
    ):
        _positive(name, value)
    if not math.isfinite(rho_ph) or rho_ph < 0.0:
        raise ValueError(f'rho_ph must be non-negative and finite, got {rho_ph!r}')
    if not math.isfinite(e) or not 0.0 <= e < 1.0:
        raise ValueError(f'e must satisfy 0 <= e < 1, got {e!r}')
    v_th = math.sqrt(kb * T / mu_kg)
    q = M_p / M_star
    a_curv = curvature_a(q)
    # The Figure 9 crossover quantity: the overflow description applies
    # where this sonic radius reaches the L1 distance (see module notes).
    r_sonic = G * M_p / (2.0 * v_th**2)

    # Orbit average. A secular caller integrates over many orbital periods,
    # so what it needs is the mass carried per unit time rather than the
    # instantaneous rate at one phase. The corotating Roche geometry is
    # defined for a circular synchronous donor, so each phase is evaluated
    # with the circular formula at that separation and the result averaged
    # in time. That quasi-static reading is this module's construction and
    # not the primary's, which has no eccentric treatment. Time weighting
    # is Kepler's, dt proportional to (1 - e cos E) dE, and the separation
    # at that anomaly carries the same factor.
    w_sum = rate_sum = power_sum = duty_rate_sum = duty_power_sum = 0.0
    applicable_sum = saturated_sum = 0.0
    for i in range(n_phase):
        ecc_anomaly = (i + 0.5) * 2.0 * math.pi / n_phase
        w = 1.0 - e * math.cos(ecc_anomaly)
        st = _phase_state(M_p, M_star, a * w, rho_ph, r_ph, v_th, q, a_curv)
        w_sum += w
        rate_sum += w * st['rate']
        power_sum += w * st['power']
        # The applicable arc surrounds periapsis, because the L1 distance
        # grows with separation while the sonic radius does not. Off that
        # arc the gas chokes at its own sonic surface first and the nozzle
        # carries nothing, so the duty-cycled average is what a dispatcher
        # should compete. What the duty cycle leaves out is the wind the
        # planet drives on the rest of the orbit, which one dispatched
        # rate cannot also carry. The returned rate is the unguarded
        # Eq. (3) average, so the closed form stays comparable with the
        # primary's own published rates; the gate is the caller's.
        if r_sonic >= st['R_L1']:
            duty_rate_sum += w * st['rate']
            duty_power_sum += w * st['power']
            applicable_sum += w
        if st['saturated']:
            saturated_sum += w

    peri = _phase_state(M_p, M_star, a * (1.0 - e), rho_ph, r_ph, v_th, q, a_curv)
    apo = _phase_state(M_p, M_star, a * (1.0 + e), rho_ph, r_ph, v_th, q, a_curv)
    return rate_sum / w_sum, dict(
        T_K=T,
        mu_kg=mu_kg,
        v_th=v_th,
        R_sonic=r_sonic,
        q=q,
        A=a_curv,
        n_phase=n_phase,
        a_periapsis=a * (1.0 - e),
        r_lobe=peri['r_lobe'],
        phi_L1=peri['phi_L1'],
        phi_ph=peri['phi_ph'],
        delta_phi=peri['delta_phi'],
        exponent_applied=peri['exponent_applied'],
        area_m2=peri['area_m2'],
        saturated=peri['saturated'],
        R_L1=peri['R_L1'],
        R_sonic_over_R_L1=r_sonic / peri['R_L1'],
        R_sonic_over_R_L1_apoapsis=r_sonic / apo['R_L1'],
        rate_periapsis_kg_s=peri['rate'],
        rate_apoapsis_kg_s=apo['rate'],
        rate_applicable_kg_s=duty_rate_sum / w_sum,
        applicable=applicable_sum > 0.0,
        applicable_orbit_fraction=applicable_sum / w_sum,
        saturated_orbit_fraction=saturated_sum / w_sum,
        power_lift_W=duty_power_sum / w_sum,
        power_lift_full_orbit_W=power_sum / w_sum,
    )
