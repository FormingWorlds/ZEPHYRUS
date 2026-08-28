"""
!!! info "`boiloff.py`"
    Bolometrically driven escape: boil-off and its capped residual.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

from scipy.special import lambertw

from zephyrus.constants import G, kb

# Provenance of the branch, all closed form:
#
# - Rate: Owen & Wu (2016, ApJ 817, 107) isothermal transonic Parker wind
#   with the photospheric Mach number in exact Lambert-W form,
#   Mdot = 4 pi G M_p Mach / (kappa c_s), Mach = sqrt(-W0(-f(x))),
#   f(x) = x^-4 exp(3 - 4/x), x = R_launch/R_B, R_B = G M_p / (2 c_s^2).
#   The exact form is used rather than their large-1/x asymptote, whose
#   order-unity prefactor is absorbed. The 1/kappa dependence enters
#   through the photospheric opacity input.
# - Wind temperature: T_eq / 2^(1/4), the explicit recommendation of
#   Misener et al. (2025, ApJ 980, 152) for the isothermal formulas.
# - Activation: the restricted Jeans parameter
#   Lambda = G M_p mu / (kB T_eq R_p) (Fossati et al. 2017, A&A 598, A90),
#   built with the composition mean molecular mass. For isothermal gas
#   Lambda = 2 R_B / R_p identically, so the Owen & Wu shutoff at
#   R_p/R_B = 0.1 is Lambda = 20 for every composition; the transfer of
#   that hydrogen-calibrated shutoff to other envelopes is an assumption
#   the activation band (15 to 35 across the literature) makes visible.
# - Bondi cap: Gupta & Schlichting (2020, MNRAS 493, 792, their Eq. 10),
#   Mdot_B = 4 pi R_B^2 c_s rho_launch exp(-G M_p / (c_s^2 R_launch)),
#   sonic-point area with the launch-level density; the launch level is
#   identified with the radiative-convective boundary, a documented
#   approximation on a static profile.
# - Luminosity cap, applied only past the activation gate:
#   Mdot_E = L / (g R_p K) with L = 4 pi R_p^2 F_int (Gupta & Schlichting
#   2019, MNRAS 487, 24, their Eq. 9). Capping the residual bolometric
#   channel by the interior luminosity sidesteps the open dispute over how
#   long core-powered mass loss survives after boil-off (Tang et al. 2024,
#   ApJ 976, 221, argue it is brief; Gupta & Schlichting argue it lasts).
#   The barrier the luminosity has to lift the gas over carries the tidal
#   reduction K(xi) of Erkaev et al. (2007, A&A 472, 329, their Eq. 17),
#   xi = R_Hill/R_p, so the cap and the energy-limited rate it competes
#   against measure the same barrier from the same reference radius; K = 1
#   recovers the untidal form and is what a caller with tides off gets.
# - Termination diagnostic: the Tang et al. (2024) Eq. (8) timescale
#   comparison, run as a diagnostic beside the rate's own exponential
#   shutoff, never as a gate.

LAMBDA_BAND = (15.0, 35.0)  # literature spread of the activation threshold


def lambda_restricted(M_p: float, R_p: float, T_eq: float, mu_kg: float) -> float:
    """Restricted Jeans parameter Lambda = G M_p mu / (kB T_eq R_p).

    Dimensionless; built with the composition mean molecular mass ``mu_kg``
    [kg] at the launch level (Fossati et al. 2017). For isothermal gas this
    equals ``2 R_B / R_p`` with the Bondi radius at ``T_eq``.
    """
    return G * M_p * mu_kg / (kb * T_eq * R_p)


def parker_mach(x: float) -> float:
    """Photospheric Mach number of the transonic isothermal Parker wind.

    Exact Lambert-W form of Owen & Wu (2016): ``Mach = sqrt(-W0(-f))`` with
    ``f = x^-4 exp(3 - 4/x)`` and ``x = R_launch / R_B <= 1``. At ``x = 1``
    the launch level sits at the sonic (Bondi) radius and the Mach number
    is 1; for small ``x`` it shuts off exponentially.

    Raises
    ------
    ValueError
        If ``x`` is outside ``(0, 1]``; the caller clamps inflated
        configurations to 1 with a flag before calling.
    """
    if not 0.0 < x <= 1.0:
        raise ValueError('parker_mach needs 0 < x <= 1')
    f = x**-4 * math.exp(3.0 - 4.0 / x)
    # f -> 1/e as x -> 1 (sonic point at the Bondi radius, Mach 1); clamp at
    # the W0 branch point against floating-point overshoot.
    if f >= math.exp(-1.0) * (1.0 - 1e-12):
        return 1.0
    w = lambertw(-f, 0)
    return math.sqrt(max(-w.real, 0.0))


def bolometric_candidate(
    M_p: float,
    R_p: float,
    T_eq: float,
    kappa_photo: float,
    launch: dict,
    F_int: float,
    lambda_gate: float,
    lambda_crit: float,
    k_tide: float = 1.0,
) -> tuple[float, dict]:
    """The bolometrically driven candidate mass-loss rate, in kg/s.

    Computed at every call: while ``lambda_gate < lambda_crit`` the
    atmosphere is inflated enough to boil off and the candidate is
    ``min(Parker rate, Bondi cap)``; past the gate the same machinery stays
    alive as a residual, additionally capped by the interior luminosity
    (see the module provenance notes).

    Parameters
    ----------
    M_p, R_p : float
        Planet mass [kg] and radius [m].
    T_eq : float
        Equilibrium temperature [K]; the wind runs at ``T_eq / 2^(1/4)``.
    kappa_photo : float
        Photospheric opacity [m^2 kg^-1]; the Parker rate scales as its
        inverse.
    launch : dict
        The photospheric working level (from
        :func:`zephyrus.profiles.photospheric_level`): uses ``p``, ``r``,
        ``mmw`` (molecular, not atomized), and ``rho``.
    F_int : float
        Interior heat flux [W m^-2], for the luminosity cap.
    lambda_gate : float
        The restricted Jeans parameter of the configuration.
    lambda_crit : float
        The activation threshold (20 by default upstream; band 15 to 35).
    k_tide : float
        Erkaev tidal reduction factor of the escape barrier, evaluated at
        ``xi = R_Hill / R_p``. Divides the luminosity cap, which is the
        only term that measures a barrier. The default of 1 is the untidal
        form; a non-positive value means the barrier has vanished and the
        cap is dropped.

    Returns
    -------
    (rate, detail)
        The candidate rate [kg/s] and a detail dict carrying the wind
        temperature, sound speed, Bondi radius, Mach number, each cap, the
        tidal factor the cap used, the activation state, and flags
        (``bondi_inflated`` when the launch level sits above the Bondi
        radius).
    """
    T_w = T_eq / 2.0**0.25
    mu = launch['mmw']
    c_s = math.sqrt(kb * T_w / mu)
    R_B = G * M_p / (2.0 * c_s**2)
    R_launch = launch['r']
    x = R_launch / R_B
    flags = {}
    if x > 1.0:
        flags['bondi_inflated'] = True  # photosphere above the sonic radius
        x = 1.0
    mach = parker_mach(x)
    mdot_parker = 4.0 * math.pi * G * M_p * mach / (kappa_photo * c_s)

    rho_launch = launch['rho']
    mdot_bondi = (
        4.0 * math.pi * R_B**2 * c_s * rho_launch * math.exp(-G * M_p / (c_s**2 * R_launch))
    )

    # Optical depth of the launch level to its own opacity, in the
    # plane-parallel form tau = kappa P / g. The Parker rate is derived from
    # a photosphere, so this reports whether the prescribed level and the
    # supplied opacity describe the same surface: tau far from 1 means they
    # do not, and the rate is being evaluated off the definition it came
    # from. Reporting only. The level is prescribed rather than solved for
    # because the activation threshold above is calibrated at a level of its
    # own, so solving here would put the gate and the rate on two surfaces.
    g_launch = G * M_p / R_launch**2
    tau_launch = kappa_photo * launch['p'] / g_launch

    active = lambda_gate < lambda_crit
    caps = [mdot_parker, mdot_bondi]
    mdot_lum = None
    if not active:
        L = 4.0 * math.pi * R_p**2 * F_int
        g = G * M_p / R_p**2
        barrier = g * R_p * k_tide  # J/kg to lift gas out, tides included
        mdot_lum = L / barrier if barrier > 0.0 else math.inf
        caps.append(mdot_lum)
    rate = min(caps)
    if mdot_lum is not None and rate == mdot_lum:
        # The interior luminosity is the binding term. Worth a flag rather
        # than an inference from the branch being past its gate: the cap
        # switches on at the gate, so a state that crosses the activation
        # threshold drops discontinuously (a factor 6.7e3 on a two Earth-mass
        # hydrogen envelope) while keeping the same label.
        flags['luminosity_capped'] = True
    return rate, dict(
        T_wind=T_w,
        c_s=c_s,
        R_B=R_B,
        x=min(R_launch / R_B, 1.0),
        mach=mach,
        mdot_parker=mdot_parker,
        mdot_bondi=mdot_bondi,
        mdot_luminosity=mdot_lum,
        k_tide=k_tide,
        active=active,
        R_sonic=R_B,
        tau_launch=tau_launch,
        p_launch=float(launch['p']),
        flags=flags,
    )


def tang_timescale_check(
    M_p: float, R_p: float, F_int: float, mdot: float, reservoirs: dict | None
) -> dict:
    """Boil-off termination diagnostic after Tang et al. (2024, Eq. 8).

    Boil-off has ended once the mass-loss timescale ``t_Mdot = M_env/Mdot``
    reaches the cooling timescale ``t_cool = G M_p M_env / (R_p L)`` with
    ``L = 4 pi R_p^2 F_int`` (their Eq. 9 with the envelope concentration
    factor at 1 and the radiative-convective boundary at ``R_p``, both
    stated approximations). Returns ``{'evaluated': False}`` when no
    reservoir masses are supplied; otherwise the two timescales and the
    verdict. Diagnostic only: it never gates the rate.
    """
    if not reservoirs:
        return {'evaluated': False}
    M_env = sum(reservoirs.values())
    if M_env <= 0.0 or mdot <= 0.0:
        return {'evaluated': False}
    L = 4.0 * math.pi * R_p**2 * F_int
    t_mdot = M_env / mdot
    t_cool = G * M_p * M_env / (R_p * L)
    return {
        'evaluated': True,
        't_mdot_s': t_mdot,
        't_cool_s': t_cool,
        'terminated': t_mdot >= t_cool,
    }
