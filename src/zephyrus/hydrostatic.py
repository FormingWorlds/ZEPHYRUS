"""
!!! info "`hydrostatic.py`"
    Hydrostatic escape: extended upper structure, per-species Jeans escape,
    and the diffusion-limited supply cap.<br>
    Authors: Ioana Balint, Viesturs Strelcs, Mara Attia
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from zephyrus.composition import ELEMENT_AMU, parse_formula, species_mass_amu
from zephyrus.constants import G, amu, kb
from zephyrus.constants import rate_floor as RATE_FLOOR_KG_S
from zephyrus.diffusion import b_mixture
from zephyrus.knudsen import sigma_mixture

# The branch evaluates escape where the gas is too rarefied to sustain a
# hydrodynamic wind, per species:
#
# - Upper structure: the Bates temperature profile
#   T(zeta) = T_exo - (T_exo - T_top) exp(-gamma zeta), in the form Yelle
#   (2024, Icarus 416, 116099, their Eq. 19) uses, anchored at the topmost
#   supplied profile level and integrated hydrostatically, in closed form,
#   in zeta = ln(p_top/p), with composition and mean mass frozen at the
#   anchor. The exobase quantities are read on this inflated structure,
#   never on photospheric values (Johnson et al. 2013, ApJL 768, L4).
# - Exobase: the lowest point where the Maxwell mean free path
#   1/(sqrt(2) sigma n) reaches the local scale height (the convention of
#   Volkov et al. 2011), with the mixture cross section of the Knudsen
#   switch.
# - Jeans escape per species: the effusion flux
#   w_J = sqrt(kT / (2 pi m)) (1 + lambda) exp(-lambda) (Yelle 2024,
#   Eq. 20), multiplied by the flat kinetic enhancement C(lambda) measured
#   in direct simulation Monte Carlo runs: about 1.7 at lambda = 6 falling
#   to about 1.4 at lambda = 15 (Volkov et al. 2011, ApJL 729, L24). Their
#   bulk-velocity correction is not applied on top, since it expresses the
#   same departure from equilibrium.
# - Diffusion-limited supply: Yelle (2024) Eqs. (9)-(11) discretized on the
#   extension: the modified mixing ratio X-tilde grows by the exponential
#   of the integrated (1 - m-tilde/m_bar) D/(D + K) factor, with the
#   thermal diffusion factor alpha = -0.25 for light species (Yelle 2024,
#   after Banks & Kockarts), and the limiting flux is the inverse of the
#   resistance integral g. Binary coefficients come from the diffusion
#   library fallback order; Blanc's law combines pairs into the mixture value.
# - Combination: the harmonic mean of the Jeans and diffusion-limited
#   fluxes, Phi = Phi_J Phi_l / (Phi_J + Phi_l) (Yelle 2024, Eq. 14), both
#   referred to the anchor area (their Eq. 15). The dominant species has no
#   supply limit (it supplies itself) and takes the Jeans flux alone.
# - Escape temperatures: T_esc,neutral = G M m / (2 kB r), the
#   lambda = 2 criterion, and the plasma escape temperature at half that
#   value because the ambipolar field shares the ion's binding with the
#   electron (Chatterjee & Pierrehumbert 2026, ApJ 998, 236, their
#   Eq. 34); a hydrostatic exobase hotter than half the gating escape
#   temperature is unstable (their Figure 10 criterion).
#
# Hydrostatic heavy-element rates are lower limits, since nonthermal
# channels are not modeled; ``hydrostatic_lower_limit`` travels with every
# result.

ALPHA_THERMAL = -0.25  # thermal diffusion factor (Yelle 2024, after Banks & Kockarts)

# Species whose Jeans rate at their profile-top fraction is below the rate
# floor skip the diffusion integrals and keep that rate, an estimate.


def volkov_flat_factor(lam: float) -> float:
    """Kinetic enhancement C(lambda) on the Jeans flux, dimensionless.

    Direct simulation Monte Carlo runs exceed the Jeans flux by a factor
    1.7 near lambda = 6, falling to 1.4 by lambda = 15 (Volkov et al.
    2011); linear between, held at the endpoint values outside, where the
    caller flags the extrapolation on either side. Below lambda = 6 the
    held value likely understates the flux, and trace light species on a
    heavy background reach that side.
    """
    if lam <= 6.0:
        return 1.7
    if lam >= 15.0:
        return 1.4
    return 1.7 + (1.4 - 1.7) * (lam - 6.0) / (15.0 - 6.0)


def jeans_effusion_velocity(T: float, m: float, lam: float) -> float:
    """Jeans effusion velocity sqrt(kT/(2 pi m)) (1 + lambda) exp(-lambda), m/s."""
    return math.sqrt(kb * T / (2.0 * math.pi * m)) * (1.0 + lam) * math.exp(-lam)


def bates_extension(
    profile,
    M_p: float,
    T_exo: float,
    gamma: float = 0.75,
    zeta_max: float = 40.0,
    n_levels: int = 400,
) -> dict:
    """Bates upper structure above the topmost profile level.

    The hydrostatic relation on the Bates temperature profile in
    ``zeta = ln(p_top / p)``, with composition and mean mass frozen at the
    anchor, integrates in closed form: ``d(1/r)/dzeta = -kB T / (G M mu)``
    and the Bates profile integrates analytically, so the structure carries
    no step error at any resolution. It ends, flagged ``unbound``, where the
    local Jeans parameter falls to 2: an isothermal-tail structure is
    unbound beyond that point and the geometry belongs to the overflow and
    boil-off machinery, not to a hydrostatic exosphere.

    Returns a dict of arrays on ``n_levels`` points spanning the extension
    (``zeta``, ``p`` [Pa], ``r`` [m], ``T`` [K], ``n`` [m^-3]), the end
    point ``zeta_end``, the frozen ``mu`` [kg], the normalized species
    ``vmr``, the anchor values ``p0`` and ``r0``, the Bates parameters
    (``T_exo``, ``T_top``, ``gamma``, ``M_p``), and ``unbound``.
    ``bates_state`` evaluates the same structure at any ``zeta``.
    """
    p0 = float(profile.p[-1])
    r0 = float(profile.r[-1])
    t_top = float(profile.T[-1])
    mu = float(profile.mmw[-1])
    # An exobase colder than the anchor would bind the exobase more strongly
    # than its anchor, inverting the construction. Floored and flagged, not
    # raised, since a coupled run's profile top can warm past the prescription.
    floored = T_exo < t_top
    if floored:
        T_exo = t_top
    # Every species at the anchor is carried, however thin: a trace light
    # species can dominate the exospheric loss.
    vmr = {
        sp: float(np.asarray(v)[-1])
        for sp, v in profile.vmr.items()
        if float(np.asarray(v)[-1]) > 0.0
    }
    tot = sum(vmr.values())
    if tot <= 0.0:
        raise ValueError('no species present at the profile top')
    vmr = {sp: x / tot for sp, x in vmr.items()}
    ext = dict(mu=mu, vmr=vmr, p0=p0, r0=r0, T_exo=T_exo, T_top=t_top, gamma=gamma, M_p=M_p)

    def excess_lambda(z: float) -> float:
        _r, t, _n, inv_r = _bates_point(ext, z)
        return G * M_p * mu * inv_r / (kb * t) - 2.0

    if excess_lambda(zeta_max) >= 0.0:
        zeta_end, unbound = zeta_max, False
    elif excess_lambda(0.0) <= 0.0:
        zeta_end, unbound = 0.0, True
    else:
        # The Jeans parameter falls monotonically with zeta (1/r falls and
        # T does not), so the crossing is unique.
        zeta_end = brentq(excess_lambda, 0.0, zeta_max, xtol=1e-12, rtol=1e-12)
        unbound = True
    zeta = np.linspace(0.0, zeta_end, n_levels) if zeta_end > 0.0 else np.zeros(1)
    r, T, n, _inv = bates_state(ext, zeta)
    ext.update(
        zeta=zeta,
        p=p0 * np.exp(-zeta),
        r=r,
        T=T,
        n=n,
        zeta_end=float(zeta_end),
        unbound=unbound,
        t_exo_floored=floored,
    )
    return ext


def _bates_point(ext: dict, zeta: float) -> tuple[float, float, float, float]:
    """Scalar form of ``bates_state``: ``(r, T, n, 1/r)`` at one ``zeta``."""
    r, T, n, inv_r = bates_state(ext, np.asarray(zeta, dtype=float))
    return float(r), float(T), float(n), float(inv_r)


def bates_state(ext: dict, zeta):
    """Radius [m], temperature [K], number density [m^-3], and 1/r at ``zeta``.

    ``T = T_exo - (T_exo - T_top) exp(-gamma zeta)`` integrates to
    ``T_exo zeta - (T_exo - T_top)(1 - exp(-gamma zeta)) / gamma``, which
    gives ``1/r`` exactly. Beyond the point where ``1/r`` reaches zero the
    radius is infinite; ``bates_extension`` stops well before it.
    """
    zeta = np.asarray(zeta, dtype=float)
    t_exo, t_top, gamma = ext['T_exo'], ext['T_top'], ext['gamma']
    decay = np.exp(-gamma * zeta)
    T = t_exo - (t_exo - t_top) * decay
    t_int = t_exo * zeta - (t_exo - t_top) * (1.0 - decay) / gamma
    inv_r = 1.0 / ext['r0'] - kb * t_int / (G * ext['M_p'] * ext['mu'])
    with np.errstate(divide='ignore'):
        r = np.where(inv_r > 0.0, 1.0 / np.where(inv_r > 0.0, inv_r, 1.0), np.inf)
    n = ext['p0'] * np.exp(-zeta) / (kb * T)
    return r, T, n, inv_r


def find_exobase(ext: dict, M_p: float, n_scan: int = 201) -> tuple[float, dict]:
    """Exobase on the extension, as ``(zeta_exo, flags)``.

    The level where the Maxwell mean free path ``1/(sqrt(2) sigma n)`` with
    the mixture cross section equals the local scale height
    ``kB T r^2 / (G M mu)``, located by a bracketing root find on the
    closed-form structure rather than snapped to a grid level, so the
    exobase state does not depend on any quadrature resolution. ``n_scan``
    points bracket the first crossing. When the extension never reaches
    that point its end is used, flagged ``exobase_not_reached``; an exobase
    at the anchor itself returns ``zeta_exo = 0``, flagged
    ``exobase_at_anchor``.
    """
    flags: dict = {}
    zeta_end = float(ext['zeta_end'])

    def excess(z: float) -> float:
        r, t, n, _inv = _bates_point(ext, z)
        sigma, _prov = sigma_mixture(ext['vmr'], t)
        mfp = 1.0 / (math.sqrt(2.0) * sigma * n)
        return mfp / (kb * t * r**2 / (G * M_p * ext['mu'])) - 1.0

    if excess(0.0) >= 0.0:
        flags['exobase_at_anchor'] = True
        return 0.0, flags
    if zeta_end <= 0.0:
        flags['exobase_not_reached'] = True
        return 0.0, flags
    grid = np.linspace(0.0, zeta_end, max(2, int(n_scan)))
    lo = 0.0
    for z in grid[1:]:
        if excess(float(z)) >= 0.0:
            return float(brentq(excess, lo, float(z), xtol=1e-12, rtol=1e-12)), flags
        lo = float(z)
    flags['exobase_not_reached'] = True
    return zeta_end, flags


def hydrostatic_rates(
    profile,
    M_p: float,
    T_exo: float,
    gamma_bates: float = 0.75,
    kzz_default: float = 3.0e2,
    n_levels: int = 400,
) -> tuple[dict, dict]:
    """Per-species hydrostatic escape mapped onto per-element rates [kg/s].

    Builds the Bates extension at the prescribed exobase temperature,
    locates the exobase, and combines the kinetic-corrected Jeans flux with
    the diffusion-limited supply by the harmonic mean, species by species
    (the module notes give the construction and its provenance). The
    species are the ones the supplied profile carries at its top level (the
    profile chemistry decides how atomized the exobase gas is); their rates
    map onto elements stoichiometrically at output.

    The eddy diffusion coefficient comes from the profile's top level when
    a ``kzz`` column is present, else ``kzz_default`` [m^2/s]. Species
    whose supply-free Jeans rate already falls below one proton mass per
    year skip the supply integrals (``pruned`` in the per-species detail)
    and return that rate at their profile-top mole fraction. It is an
    estimate, not a bound, since diffusive enrichment can raise a light
    species' full rate above it (2.0 to 2.6 times on trace hydrogen and
    helium profiles). The shortcut exists because the integrals dominate
    the cost of the branch on many-species profiles.

    ``n_levels`` sets the quadrature resolution of the supply integrals,
    which are first-order accurate in the log-pressure step;
    ``hydrostatic_rates_refined`` refines it to a target. The exobase does
    not depend on it.

    Returns ``(per_element, detail)``; the detail dict carries the exobase
    state, both escape temperatures, the ``dominant`` species that supplies
    itself without a diffusion cap, the per-species terms, coefficient
    provenance, and flags (including ``hydrostatic_lower_limit``, which is
    always on: nonthermal loss channels are absent).
    """
    structure = _exobase_structure(profile, M_p, T_exo, gamma_bates)
    return _rates_on_structure(profile, M_p, structure, kzz_default, n_levels)


def _exobase_structure(profile, M_p: float, T_exo: float, gamma_bates: float) -> tuple:
    """The extension, the exobase on it, and the structure flags."""
    ext = bates_extension(profile, M_p, T_exo, gamma=gamma_bates, n_levels=2)
    zeta_x, flags = find_exobase(ext, M_p)
    if ext['t_exo_floored']:
        flags['t_exo_floored_to_profile_top'] = True
    if ext['unbound']:
        flags['extension_unbound'] = True
    return ext, zeta_x, flags


def _rates_on_structure(
    profile, M_p: float, structure: tuple, kzz_default: float, n_levels: int
) -> tuple[dict, dict]:
    """Per-species rates on a fixed exobase structure, at one quadrature grid."""
    ext, zeta_x, structure_flags = structure
    flags = dict(structure_flags)
    r_x, t_x, n_x, _inv = _bates_point(ext, zeta_x)
    r_0 = float(ext['r0'])
    comp = ext['vmr']
    m_bar = ext['mu']
    dominant = max(comp, key=comp.get)

    if profile.kzz is not None:
        k_eddy = float(profile.kzz[-1])
    else:
        k_eddy = kzz_default

    # Quadrature nodes from the anchor to the exobase itself, so the grid
    # sets only the integration error and never where the exobase sits.
    if zeta_x > 0.0:
        zeta = np.linspace(0.0, zeta_x, max(2, int(n_levels)))
    else:
        zeta = np.zeros(1)
    _r, T, n, _ = bates_state(ext, zeta)
    dtdz_all = ext['gamma'] * (ext['T_exo'] - ext['T_top']) * np.exp(-ext['gamma'] * zeta)

    per_species_rate: dict = {}
    detail_species: dict = {}
    b_prov: dict = {}
    for sp, x0 in comp.items():
        m_i = species_mass_amu(sp) * amu
        lam_x = G * M_p * m_i / (kb * t_x * r_x)
        w_j = jeans_effusion_velocity(t_x, m_i, lam_x)
        c_enh = volkov_flat_factor(lam_x)
        if lam_x > 15.0 or lam_x < 6.0:
            flags['volkov_extrapolated'] = True

        area = 4.0 * math.pi * r_0**2
        phi_jeans_unlimited = (r_x / r_0) ** 2 * c_enh * w_j * x0 * n_x
        pruned = False
        if sp == dominant:
            x_tilde_x = x0
            phi_l = math.inf
        elif area * m_i * phi_jeans_unlimited < RATE_FLOOR_KG_S:
            # Not a bound: diffusive enrichment can raise a light species'
            # full rate above this one by a few times.
            x_tilde_x = x0
            phi_l = math.inf
            pruned = True
        else:
            int1 = 0.0
            g_int = 0.0
            x_tilde = x0
            for k in range(len(zeta) - 1):
                dz = float(zeta[k + 1] - zeta[k])
                t_k = float(T[k])
                n_k = float(n[k])
                bmk, prov = b_mixture(sp, comp, t_k)
                b_prov.update(prov)
                d_mol = bmk / n_k
                dtdz = float(dtdz_all[k])
                m_tilde = m_i + ALPHA_THERMAL * dtdz * (m_bar / t_k)
                int1 += (1.0 - m_tilde / m_bar) * (d_mol / (d_mol + k_eddy)) * dz
                x_tilde = x0 * math.exp(min(int1, 700.0))
                g_int += (
                    kb * t_k * r_0**2 / (x_tilde * n_k * G * M_p * m_bar * (d_mol + k_eddy))
                ) * dz
            x_tilde_x = x_tilde
            phi_l = 1.0 / g_int if g_int > 0.0 else math.inf

        phi_jeans = (r_x / r_0) ** 2 * c_enh * w_j * x_tilde_x * n_x
        if math.isinf(phi_l):
            phi = phi_jeans
        else:
            phi = phi_jeans * phi_l / (phi_jeans + phi_l)
        per_species_rate[sp] = area * m_i * phi
        detail_species[sp] = dict(
            lambda_exo=lam_x,
            w_jeans=w_j,
            volkov_C=c_enh,
            X_tilde_exo=x_tilde_x,
            phi_per_area_r0=phi,
            phi_jeans=phi_jeans,
            phi_diffusion=phi_l,
            dl_bypass=(sp == dominant),
            pruned=pruned,
        )

    per_element: dict = {}
    for sp, rate in per_species_rate.items():
        counts = parse_formula(sp)
        m_sp = species_mass_amu(sp)
        for el, cnt in counts.items():
            per_element[el] = per_element.get(el, 0.0) + rate * (cnt * ELEMENT_AMU[el] / m_sp)

    flags['hydrostatic_lower_limit'] = True

    t_esc_neutral = G * M_p * m_bar / (2.0 * kb * r_x)
    t_esc_plasma = 0.5 * t_esc_neutral
    lam_exo_bulk = G * M_p * m_bar / (kb * t_x * r_x)

    detail = dict(
        r_exo=r_x,
        zeta_exo=zeta_x,
        T_exo=t_x,
        n_exo=n_x,
        r_anchor=r_0,
        composition=comp,
        dominant=dominant,
        m_bar=m_bar,
        species=detail_species,
        per_species_rate=per_species_rate,
        b_provenance=b_prov,
        T_esc_neutral=t_esc_neutral,
        T_esc_plasma=t_esc_plasma,
        lambda_exo_bulk=lam_exo_bulk,
        K_eddy=k_eddy,
        n_levels=n_levels,
        flags=flags,
    )
    return per_element, detail


def hydrostatic_rates_refined(
    profile,
    M_p: float,
    T_exo: float,
    gamma_bates: float = 0.75,
    kzz_default: float = 3.0e2,
    n_levels_min: int = 200,
    n_levels_max: int = 3200,
    rtol: float = 1.0e-2,
) -> tuple[dict, dict]:
    """Hydrostatic rates refined until the quadrature stops moving them.

    The exobase is located once, on the closed-form structure, so only the
    supply integrals depend on the grid. They are first-order accurate in
    the log-pressure step, so the change between a grid and its refinement
    estimates the remaining error. The resolution doubles from
    ``n_levels_min`` until two successive relative changes in the bulk rate
    fall below ``rtol``, so a single chance agreement between two grids does
    not end the refinement, or until ``n_levels_max`` is reached.
    The finest grid's rates are returned, never an extrapolation.

    The detail dict gains a ``convergence`` entry recording the levels used,
    the last relative change in the bulk rate and the worst one across the
    species, the grid history, and whether the target was met. A call that
    exhausts ``n_levels_max`` without meeting it reports ``converged`` False
    rather than raising: the rate is still the best available and the
    caller can see how far from the target it is.
    """
    structure = _exobase_structure(profile, M_p, T_exo, gamma_bates)
    n = max(2, int(n_levels_min))
    n_max = max(n, int(n_levels_max))
    per_element, detail = _rates_on_structure(profile, M_p, structure, kzz_default, n)
    history = [(n, sum(per_element.values()))]
    rel_bulk = math.inf
    rel_species = math.inf
    below = 0
    converged = False
    while n < n_max:
        n_fine = min(2 * n, n_max)
        per_fine, detail_fine = _rates_on_structure(
            profile, M_p, structure, kzz_default, n_fine
        )
        rel_bulk = _relative_change(sum(per_element.values()), sum(per_fine.values()))
        rel_species = max(
            (
                _relative_change(per_element.get(el, 0.0), per_fine.get(el, 0.0))
                for el in set(per_element) | set(per_fine)
            ),
            default=0.0,
        )
        per_element, detail, n = per_fine, detail_fine, n_fine
        history.append((n, sum(per_element.values())))
        below = below + 1 if rel_bulk <= rtol else 0
        if below >= 2:
            converged = True
            break
    detail['convergence'] = dict(
        n_levels=n,
        n_levels_min=int(n_levels_min),
        n_levels_max=n_max,
        rtol=float(rtol),
        rel_change_bulk=float(rel_bulk),
        rel_change_worst_species=float(rel_species),
        converged=bool(converged),
        history=history,
    )
    return per_element, detail


def _relative_change(coarse: float, fine: float) -> float:
    """Relative difference between two grids, zero when both vanish."""
    scale = max(abs(coarse), abs(fine))
    if scale == 0.0:
        return 0.0
    return abs(fine - coarse) / scale


def gate_unstable(
    T_exo: float, detail: dict, gate: str, f_plus_exo: float
) -> tuple[bool, bool]:
    """Escape-temperature gate on the hydrostatic branch.

    The hydrostatic equilibrium is unstable when the exobase temperature
    exceeds half the gating escape temperature (the stability criterion of
    Chatterjee & Pierrehumbert 2026, their Figure 10). ``gate`` selects the
    neutral or the plasma escape temperature; both are always computed, and
    the returned ``contested`` marks points where the two conventions
    disagree, which is where the (unmodeled) ion physics decides the branch
    assignment and both branch rates belong in the diagnostics.
    ``f_plus_exo`` is carried for the caller's record and does not gate.
    """
    un_neutral = T_exo > 0.5 * detail['T_esc_neutral']
    un_plasma = T_exo > 0.5 * detail['T_esc_plasma']
    unstable = un_plasma if gate == 'plasma' else un_neutral
    return unstable, (un_neutral != un_plasma)
