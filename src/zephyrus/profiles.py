"""
!!! info "`profiles.py`"
    Atmosphere-profile container, interpolation, and escape working levels.<br>
    Author: Mara Attia
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from zephyrus.composition import atomize, species_mass_amu
from zephyrus.constants import G, amu, kb

# Mixing ratios below zero by less than this are treated as solver noise; a
# chemistry or transport solver can return a small negative mole fraction
# where a species is absent, and a coupled run must not die on it.
VMR_NOISE_FLOOR = 1.0e-12

# Murray-Clay et al. (2009, ApJ 693, 23) photoionization cross section at
# their representative 20 eV photon energy: sigma_nu0 = 6e-18 (h nu / 13.6
# eV)^-3 cm^2. Converted to m^2. Used by the Lopez (2017) wind-base pressure.
SIGMA_NU0 = 6.0e-18 * (20.0 / 13.6) ** -3 * 1e-4  # [m^2]


@dataclass
class Profile:
    """One-dimensional atmosphere profile, ordered base to top.

    Attributes
    ----------
    p : ndarray
        Pressure per level [Pa], strictly decreasing with index.
    r : ndarray
        Radius per level [m], strictly increasing with index.
    T : ndarray
        Temperature per level [K]; realistic profiles are non-monotone.
    vmr : dict
        Species name to an array of volume mixing ratios per level.
    mmw : ndarray
        Mean molecular mass per level [kg per particle].
    kzz : ndarray or None
        Eddy diffusion coefficient per level [m^2 s^-1], optional.
    """

    p: np.ndarray
    r: np.ndarray
    T: np.ndarray
    vmr: dict
    mmw: np.ndarray
    kzz: np.ndarray | None = None

    def validate(self) -> None:
        """Raise ``ValueError`` on a malformed profile.

        Pressure must decrease and radius increase strictly with index;
        temperature is unconstrained beyond positivity. Every mixing-ratio
        array must share the level count, be finite, and be non-negative.

        Finiteness is checked before the sign comparisons, because a
        comparison against NaN is false and a NaN would otherwise pass every
        positivity test and surface far downstream as an error naming some
        unrelated quantity. Mixing ratios below zero by less than
        ``VMR_NOISE_FLOOR`` are solver noise and pass; the consumers ignore
        non-positive weights and renormalize over the rest.
        """
        p, r, T, mmw = map(np.asarray, (self.p, self.r, self.T, self.mmw))
        if not (len(p) == len(r) == len(T) == len(mmw)):
            raise ValueError('profile arrays must share one length')
        if len(p) < 3:
            raise ValueError('profile needs at least 3 levels')
        for name, arr in (('p', p), ('r', r), ('T', T), ('mmw', mmw)):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f'{name} carries a non-finite value')
        if not (np.all(np.diff(p) < 0) and np.all(np.diff(r) > 0)):
            raise ValueError('p must decrease and r increase strictly with index')
        if np.any(p <= 0) or np.any(T <= 0) or np.any(mmw <= 0):
            raise ValueError('p, T, mmw must be positive')
        for sp, x in self.vmr.items():
            arr = np.asarray(x, dtype=float)
            if len(arr) != len(p):
                raise ValueError(f'vmr[{sp}] length mismatch')
            if not np.all(np.isfinite(arr)):
                raise ValueError(f'vmr[{sp}] carries a non-finite value')
            if np.any(arr < -VMR_NOISE_FLOOR):
                raise ValueError(f'vmr[{sp}] is negative beyond solver noise')
        if self.vmr:
            total = sum(np.clip(np.asarray(x, dtype=float), 0.0, None) for x in self.vmr.values())
            if np.any(np.asarray(total) <= 0.0):
                raise ValueError('every level needs at least one species present')


def isothermal_profile(
    M_p: float,
    R_p: float,
    T: float,
    composition: dict[str, float],
    p_surf: float,
    p_top: float,
    n_levels: int = 120,
) -> Profile:
    """Isothermal hydrostatic profile for a fixed molecular composition.

    Integrates ``dr = -(k T r^2 / (G M mu)) d ln p`` outward from ``R_p`` at
    constant temperature and composition. The integration truncates, before
    reaching ``p_top``, at the level where the local Jeans parameter
    ``G M mu / (k T r)`` drops below 2.2, because an isothermal structure is
    unbound beyond that point and a hydrostatic profile there would be
    meaningless.

    Parameters
    ----------
    M_p, R_p : float
        Planet mass [kg] and base radius [m].
    T : float
        Temperature [K].
    composition : dict
        Species name to mole fraction (renormalized internally).
    p_surf, p_top : float
        Base and requested top pressure [Pa].
    n_levels : int
        Number of levels in log pressure.

    Returns
    -------
    Profile

    Raises
    ------
    ValueError
        If the structure is unbound at the surface (fewer than three bound
        levels), which is not a physically posed hydrostatic input.
    """
    tot = sum(composition.values())
    mu = sum(x * species_mass_amu(sp) for sp, x in composition.items()) / tot * amu
    lnp = np.linspace(np.log(p_surf), np.log(p_top), n_levels)
    r = np.empty(n_levels)
    r[0] = R_p
    last = n_levels - 1
    for i in range(n_levels - 1):
        H = kb * T * r[i] ** 2 / (G * M_p * mu)
        r[i + 1] = r[i] - H * (lnp[i + 1] - lnp[i])
        if G * M_p * mu / (kb * T * r[i + 1]) < 2.2:
            last = i + 1
            break
    if last < 2:
        raise ValueError('isothermal profile unbound at the surface; check inputs')
    lnp, r = lnp[: last + 1], r[: last + 1]
    n = last + 1
    vmr = {sp: np.full(n, x / tot) for sp, x in composition.items()}
    return Profile(
        p=np.exp(lnp), r=r, T=np.full(n, float(T)), vmr=vmr, mmw=np.full(n, mu), kzz=None
    )


def interp_at_pressure(profile: Profile, p_target: float) -> dict:
    """Level state interpolated at a target pressure, linear in log pressure.

    Returns a dict with the interpolated radius ``r`` [m], temperature ``T``
    [K], mean molecular mass ``mmw`` [kg], number density ``n`` [m^-3], mass
    density ``rho`` [kg m^-3], per-species ``vmr``, ``kzz`` (or None), and
    the pressure ``p`` [Pa] itself.
    """
    lp = np.log(profile.p)
    lt = np.log(p_target)
    # p decreases with index; np.interp needs increasing x.
    x = lp[::-1]
    r = np.interp(lt, x, profile.r[::-1])
    T = np.interp(lt, x, profile.T[::-1])
    mmw = np.interp(lt, x, profile.mmw[::-1])
    vmr = {sp: float(np.interp(lt, x, np.asarray(v)[::-1])) for sp, v in profile.vmr.items()}
    kzz = float(np.interp(lt, x, profile.kzz[::-1])) if profile.kzz is not None else None
    n = p_target / (kb * T)
    return dict(
        p=float(p_target),
        r=float(r),
        T=float(T),
        mmw=float(mmw),
        n=float(n),
        rho=float(n * mmw),
        vmr=vmr,
        kzz=kzz,
    )


def pressure_at_radius(profile: Profile, r_target: float) -> float:
    """Pressure interpolated at a target radius, log-linear in pressure.

    Clamps to the endpoint pressures outside the covered radius range.
    """
    lp = np.log(profile.p)
    return float(np.exp(np.interp(r_target, profile.r, lp)))


def photospheric_level(profile: Profile, p_photo: float = 2000.0) -> tuple[dict, dict]:
    """The photospheric working level for the energy-limited geometric factor.

    The level is placed at ``p_photo`` (default 20 mbar, the
    optical-photosphere-type level of Baumeister et al. 2023, A&A 675,
    A122). When the profile does not span that pressure, the nearest end
    level is used and the ``photo_clamped`` flag raised.

    Returns
    -------
    (level, flags)
        The interpolated level dict and a flags dict.
    """
    if profile.p[0] < p_photo:
        return interp_at_pressure(profile, float(profile.p[0])), {'photo_clamped': True}
    if profile.p[-1] > p_photo:
        return interp_at_pressure(profile, float(profile.p[-1])), {'photo_clamped': True}
    return interp_at_pressure(profile, p_photo), {}


def lopez_base_pressure(mu_kg: float, g: float) -> float:
    """Lopez (2017) XUV wind-base pressure, in Pa.

    ``P_base = mu g / sigma_nu0`` with the Murray-Clay et al. (2009)
    photoionization cross section at 20 eV, about a nanobar for a hot
    Jupiter: the pressure of the tau = 1 level for XUV photons, where the
    heating is deposited and the wind is launched (Lopez 2017, MNRAS 472,
    245, their Section 2). ``mu_kg`` is the local mean particle mass [kg]
    and ``g`` the local gravity [m s^-2].
    """
    return mu_kg * g / SIGMA_NU0


def wind_base_level(
    profile: Profile,
    M_p: float,
    method: str = 'lopez',
    fixed_pressure: float = 5.0,
    boreas_scalars: dict | None = None,
) -> tuple[dict, dict]:
    """Locate the XUV wind base on the profile.

    Three methods:

    - ``'lopez'`` (default): fixed-point iteration of the Lopez (2017)
      base pressure ``P_base = mu(P) g(r(P)) / sigma_nu0`` on the profile
      (two to four passes converge). When the physical base pressure lies
      above the profile top (``P_base < p_top``), the level clamps to the
      topmost level, the ``base_clamped`` flag is raised, and the clamp
      distance in pressure decades is recorded; callers wanting the base on
      an extended upper structure evaluate it there instead.
    - ``'fixed_pressure'``: the level at ``fixed_pressure`` [Pa].
    - ``'boreas'``: the base radius from the BOREAS mass-loss solver
      (optional dependency), translated to a profile pressure; falls back
      to ``'lopez'`` with the ``base_method_fallback`` flag when BOREAS is
      absent or does not converge. Requires ``boreas_scalars`` with keys
      ``R_p`` [m], ``T_eq`` [K], and ``F_xuv`` [W m^-2].

    Returns
    -------
    (level, flags)
        The interpolated level dict and a flags dict. The level carries
        ``p_physical``, the base pressure the method asked for before any
        clamp, which equals ``p`` whenever the profile covers it. The flags
        report clamps and fallbacks only.
    """
    flags: dict = {}
    if method == 'fixed_pressure':
        p_target = min(max(fixed_pressure, float(profile.p[-1])), float(profile.p[0]))
        if p_target != fixed_pressure:
            flags['base_clamped'] = True
            flags['base_clamp_decades'] = abs(float(np.log10(fixed_pressure / p_target)))
        level = interp_at_pressure(profile, p_target)
        level['p_physical'] = float(fixed_pressure)
        return level, flags

    if method == 'boreas':
        p_boreas = _boreas_base_pressure(profile, M_p, boreas_scalars)
        if p_boreas is not None:
            p_target = min(max(p_boreas, float(profile.p[-1])), float(profile.p[0]))
            if p_target != p_boreas:
                flags['base_clamped'] = True
                flags['base_clamp_decades'] = abs(float(np.log10(p_boreas / p_target)))
            level = interp_at_pressure(profile, p_target)
            level['p_physical'] = float(p_boreas)
            return level, flags
        flags['base_method_fallback'] = 'lopez'

    # Lopez (2017) fixed point: P depends on mu and g, which depend on the
    # level P selects; iterate from the profile top downward.
    lev = interp_at_pressure(profile, float(profile.p[-1]))
    p_target = None
    for _ in range(6):
        g = G * M_p / lev['r'] ** 2
        p_new = lopez_base_pressure(lev['mmw'], g)
        p_new = min(max(p_new, float(profile.p[-1])), float(profile.p[0]))
        if p_target is not None and abs(np.log(p_new / p_target)) < 1e-6:
            p_target = p_new
            break
        p_target = p_new
        lev = interp_at_pressure(profile, p_target)
    g = G * M_p / lev['r'] ** 2
    p_phys = lopez_base_pressure(lev['mmw'], g)
    if p_phys < profile.p[-1] * (1.0 - 1e-9):
        # The profile top is deeper than the physical base level.
        flags['base_clamped'] = True
        flags['base_clamp_decades'] = float(np.log10(profile.p[-1] / p_phys))
        p_target = float(profile.p[-1])
    level = interp_at_pressure(profile, p_target)
    level['p_physical'] = float(p_phys)
    return level, flags


def _boreas_base_pressure(profile: Profile, M_p: float, scalars: dict | None) -> float | None:
    """XUV base pressure from the BOREAS solver, or None on any failure.

    Feeds BOREAS the same scalars its PROTEUS wrapper uses (equilibrium
    temperature, XUV flux, planet mass and radius, and the mass mixing
    ratios of the gases it models, taken here at the profile top), runs its
    bulk mass-loss solve, and converts the returned XUV radius to a profile
    pressure. Every failure mode (missing dependency, missing scalars,
    unsupported composition, non-convergence) returns None so the caller
    can fall back.
    """
    if not scalars:
        return None
    try:
        import boreas  # type: ignore

        params = boreas.ModelParams()
        params.albedo = 0.0
        params.Teq = float(scalars['T_eq'])
        params.FXUV = float(scalars['F_xuv']) * 1e3  # W m^-2 -> erg cm^-2 s^-1
        params.rplanet = float(scalars['R_p']) * 1e2  # m -> cm
        params.mplanet = float(M_p) * 1e3  # kg -> g
        # Mass mixing ratios of the BOREAS-supported gases at the profile top.
        supported = set(params.kappa.keys())
        top = {sp: float(np.asarray(v)[-1]) for sp, v in profile.vmr.items()}
        masses = {sp: species_mass_amu(sp) for sp in top}
        norm = sum(top[sp] * masses[sp] for sp in top if sp in supported)
        if norm <= 0.0:
            return None
        for sp in top:
            if sp in supported:
                setattr(params, f'X_{sp}', top[sp] * masses[sp] / norm)
        params._recompute_composites()
        params._init_opacities()
        result = boreas.MassLoss(params).compute_mass_loss_parameters(
            [params.mplanet], [params.rplanet], [params.Teq]
        )[0]
        if result.get('regime') == 'SKIPPED' or 'RXUV' not in result:
            return None
        r_xuv = float(result['RXUV']) * 1e-2  # cm -> m
        return pressure_at_radius(profile, r_xuv)
    except Exception:
        return None


def atomized_element_fractions(level: dict) -> dict[str, float]:
    """Element mole fractions of the atomized composition at a level dict."""
    return atomize(level['vmr'])
