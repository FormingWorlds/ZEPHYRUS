"""
!!! info "`hydrostatic.py`"
    Hydrostatic escape: extended upper structure, per-species Jeans escape,
    and the diffusion-limited supply cap.<br>
    Authors: Ioana Balint, Viesturs Strelcs, Mara Attia
"""

from __future__ import annotations

import math

import numpy as np

from zephyrus.composition import ELEMENT_AMU, parse_formula, species_mass_amu
from zephyrus.constants import G, amu, kb
from zephyrus.diffusion import b_mixture
from zephyrus.knudsen import sigma_mixture

# The branch evaluates escape where the gas is too rarefied to sustain a
# hydrodynamic wind, per species:
#
# - Upper structure: the Bates temperature profile
#   T(zeta) = T_exo - (T_exo - T_top) exp(-gamma zeta), in the form Yelle
#   (2024, Icarus 416, 116099, their Eq. 19) uses, anchored at the topmost
#   supplied profile level and integrated hydrostatically in
#   zeta = ln(p_top/p). Composition and mean mass are frozen at the anchor
#   on the extension. Evaluating the exobase quantities on this extended,
#   inflated structure rather than on photospheric values is essential:
#   the exobase Jeans parameter can differ from the photospheric one by an
#   order of magnitude, and using the latter biases rates toward false
#   retention by up to three decades (Johnson et al. 2013, ApJL 768, L4).
# - Exobase: the first level where the Maxwell mean free path
#   1/(sqrt(2) sigma n) reaches the local scale height (the convention of
#   Volkov et al. 2011), with the mixture cross section of the Knudsen
#   switch.
# - Jeans escape per species: the effusion flux
#   w_J = sqrt(kT / (2 pi m)) (1 + lambda) exp(-lambda) (Yelle 2024,
#   Eq. 20), multiplied by the flat kinetic enhancement C(lambda) measured
#   in direct simulation Monte Carlo runs: about 1.7 at lambda = 6 falling
#   to about 1.4 at lambda = 15 (Volkov et al. 2011, ApJL 729, L24). Their
#   companion bulk-velocity correction is deliberately not applied on top:
#   the two express the same departure from equilibrium and applying both
#   double-counts. Beyond lambda = 15 the factor is held at 1.4, a flagged
#   extrapolation.
# - Diffusion-limited supply: Yelle (2024) Eqs. (9)-(11) discretized on the
#   extension: the modified mixing ratio X-tilde grows by the exponential
#   of the integrated (1 - m-tilde/m_bar) D/(D + K) factor, with the
#   thermal diffusion factor alpha = -0.25 for light species (Yelle 2024,
#   after Banks & Kockarts), and the limiting flux is the inverse of the
#   resistance integral g. Binary coefficients come from the diffusion
#   library ladder; Blanc's law combines pairs into the mixture value.
# - Combination: the harmonic mean of the Jeans and diffusion-limited
#   fluxes, Phi = Phi_J Phi_l / (Phi_J + Phi_l) (Yelle 2024, Eq. 14), both
#   referred to the anchor area (their Eq. 15). The dominant species has no
#   supply limit (it supplies itself) and takes the Jeans flux alone.
# - Escape temperatures: T_esc,neutral = G M m / (2 kB r), the
#   lambda = 2 criterion, and the plasma escape temperature at half that
#   value because the ambipolar field shares the ion's binding with the
#   electron (Chatterjee & Pierrehumbert 2026, arXiv:2412.05188, their
#   Eq. 34); a hydrostatic exobase hotter than half the gating escape
#   temperature is unstable (their Figure 10 criterion) and callers
#   re-route such points to the hydrodynamic branch.
#
# Hydrostatic heavy-element rates are lower limits: the non-thermal
# channels (ion outflow, photochemical ejection, sputtering) that dominate
# heavy-species loss in this regime are not modeled; the
# ``hydrostatic_lower_limit`` flag travels with every result.

ALPHA_THERMAL = -0.25  # thermal diffusion factor (Yelle 2024, after Banks & Kockarts)

# Rates below one proton mass per Julian year are numerical artifacts on
# any planetary reservoir; species whose supply-free Jeans rate already
# sits below this floor skip the diffusion integrals (their harmonic-mean
# rate could only be smaller).
RATE_FLOOR_KG_S = 1.67262192369e-27 / 3.15576e7


def volkov_flat_factor(lam: float) -> float:
    """Kinetic enhancement C(lambda) on the Jeans flux, dimensionless.

    Direct simulation Monte Carlo runs exceed the Jeans flux by a factor
    1.7 near lambda = 6, falling to 1.4 by lambda = 15 (Volkov et al.
    2011); linear between, held at the endpoint values outside, where the
    caller flags the extrapolation.
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

    Integrates the hydrostatic relation on the Bates temperature profile in
    ``zeta = ln(p_top / p)`` with composition and mean mass frozen at the
    anchor. The integration stops, flagged ``unbound``, where the local
    Jeans parameter drops below 2: an isothermal-tail structure is unbound
    beyond that point and the geometry belongs to the overflow and boil-off
    machinery, not to a hydrostatic exosphere.

    Returns a dict of arrays over the extension (``zeta``, ``p``, ``r``,
    ``T``, ``n``) plus the frozen ``mu`` [kg], the normalized species
    ``vmr``, the anchor values ``p0`` and ``r0``, and ``unbound``.
    """
    p0 = float(profile.p[-1])
    r0 = float(profile.r[-1])
    t_top = float(profile.T[-1])
    mu = float(profile.mmw[-1])
    vmr = {
        sp: float(np.asarray(v)[-1])
        for sp, v in profile.vmr.items()
        if float(np.asarray(v)[-1]) > 1e-8
    }
    tot = sum(vmr.values())
    vmr = {sp: x / tot for sp, x in vmr.items()}
    zeta = np.linspace(0.0, zeta_max, n_levels)
    T = T_exo - (T_exo - t_top) * np.exp(-gamma * zeta)
    r = np.empty(n_levels)
    r[0] = r0
    unbound = False
    last = n_levels - 1
    for i in range(n_levels - 1):
        H = kb * T[i] * r[i] ** 2 / (G * M_p * mu)
        r[i + 1] = r[i] + H * (zeta[i + 1] - zeta[i])
        lam_next = G * M_p * mu / (kb * T[i + 1] * r[i + 1])
        if lam_next < 2.0:
            unbound = True
            last = i + 1
            break
    zeta, T, r = zeta[: last + 1], T[: last + 1], r[: last + 1]
    p = p0 * np.exp(-zeta)
    n = p / (kb * T)
    return dict(zeta=zeta, p=p, r=r, T=T, n=n, mu=mu, vmr=vmr, p0=p0, r0=r0, unbound=unbound)


def find_exobase(ext: dict, M_p: float) -> tuple[int, dict]:
    """Exobase index on the extension: mean free path equals scale height.

    Maxwell mean free path with the mixture cross section, against the
    local scale height ``kB T r^2 / (G M mu)``. When the extension never
    reaches that point the top level is used, flagged
    ``exobase_not_reached``; an exobase at the anchor itself keeps one
    integration interval so the supply integrals exist, flagged
    ``exobase_at_anchor``.
    """
    flags: dict = {}
    idx = None
    for i in range(len(ext['zeta'])):
        sigma, _prov = sigma_mixture(ext['vmr'], float(ext['T'][i]))
        mfp = 1.0 / (math.sqrt(2.0) * sigma * ext['n'][i])
        H = kb * ext['T'][i] * ext['r'][i] ** 2 / (G * M_p * ext['mu'])
        if mfp >= H:
            idx = i
            break
    if idx is None:
        idx = len(ext['zeta']) - 1
        flags['exobase_not_reached'] = True
    if idx == 0:
        flags['exobase_at_anchor'] = True
        idx = 1
    return idx, flags


def hydrostatic_rates(
    profile,
    M_p: float,
    T_exo: float,
    gamma_bates: float = 0.75,
    kzz_default: float = 3.0e2,
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
    year skip the supply integrals (``pruned`` in the per-species detail):
    the harmonic mean could only be smaller, and the cost of the integrals
    dominates the branch on many-species profiles.

    Returns ``(per_element, detail)``; the detail dict carries the exobase
    state, both escape temperatures, the ``dominant`` species that supplies
    itself without a diffusion cap, the per-species terms, coefficient
    provenance, and flags (including ``hydrostatic_lower_limit``, which is
    always on: non-thermal loss channels are absent).
    """
    ext = bates_extension(profile, M_p, T_exo, gamma=gamma_bates)
    i_x, flags = find_exobase(ext, M_p)
    if ext.get('unbound'):
        flags['extension_unbound'] = True
    r_x, t_x, n_x = float(ext['r'][i_x]), float(ext['T'][i_x]), float(ext['n'][i_x])
    r_0 = float(ext['r0'])
    comp = ext['vmr']
    m_bar = ext['mu']
    dominant = max(comp, key=comp.get)

    if profile.kzz is not None:
        k_eddy = float(profile.kzz[-1])
    else:
        k_eddy = kzz_default

    zeta = ext['zeta'][: i_x + 1]
    T = ext['T'][: i_x + 1]
    n = ext['n'][: i_x + 1]

    per_species_rate: dict = {}
    detail_species: dict = {}
    b_prov: dict = {}
    for sp, x0 in comp.items():
        m_i = species_mass_amu(sp) * amu
        lam_x = G * M_p * m_i / (kb * t_x * r_x)
        w_j = jeans_effusion_velocity(t_x, m_i, lam_x)
        c_enh = volkov_flat_factor(lam_x)
        if lam_x > 15.0:
            flags['volkov_extrapolated'] = True

        area = 4.0 * math.pi * r_0**2
        phi_jeans_unlimited = (r_x / r_0) ** 2 * c_enh * w_j * x0 * n_x
        pruned = False
        if sp == dominant:
            x_tilde_x = x0
            phi_l = math.inf
        elif area * m_i * phi_jeans_unlimited < RATE_FLOOR_KG_S:
            # The supply-free rate is already numerically negligible; the
            # harmonic mean with any supply limit is smaller still.
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
                dtdz = (float(T[k + 1]) - t_k) / dz
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
        flags=flags,
    )
    return per_element, detail


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
