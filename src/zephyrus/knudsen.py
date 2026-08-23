"""
!!! info "`knudsen.py`"
    Neutral collision cross sections and the sonic-point Knudsen switch.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math

from zephyrus.composition import BONDI_VDW_RADIUS_A, parse_formula
from zephyrus.constants import kb, kb_cgs

SQRT2 = math.sqrt(2.0)

# Diagnostic band on the switch threshold: the transition Knudsen number is
# heating-geometry physics, not a free parameter. Direct simulation Monte
# Carlo runs place it near 0.1 for a sharp heating layer and near 1 for
# distributed heating (Johnson et al. 2013, ApJL 768, L4); the upper edge
# extends the band to 3. Printed alongside every switch verdict, never
# configurable.
KN_BAND = (0.1, 3.0)

# ---------------------------------------------------------------------------
# Rung 1: Laricchiuta et al. (2009, Eur. Phys. J. D 54, 607) phenomenological
# collision integrals. Their Eq. (2) is a double-sigmoid fit to the reduced
# collision integral in x = ln(kT/eps0), Eq. (3) gives the fit coefficients
# a_i as polynomials in the pair parameter beta (coefficients below from
# their electronic-appendix Table 3, neutral-neutral case, m = 6), and
# Eq. (4) sets the dimensional scale sigma^2 = (x0 r_e)^2 with
# x0 = xi1 beta^xi2 (their Table 4). The pair parameters (beta, eps0 in meV,
# r_e in Angstrom) are their Table 5. The momentum-transfer cross section is
# sigma_diff = pi sigma^2 Omega^(1,1)*, with the factor pi converting the
# reduced integral to a cross section. The implementation reproduces the
# measured room-temperature viscosities of N2, O2, CO, and CO2 to within
# 7 percent (see the companion tests).
# ---------------------------------------------------------------------------

# Table 3 (neutral-neutral, m = 6): rows are (c0, c1, c2) of
# a_i(beta) = c0 + c1 beta + c2 beta^2, for i = 1..7.
_TABLE3_M6 = {
    'omega11': [
        (7.884756e-1, -2.438494e-2, 0.0),
        (-2.952759e-1, -1.744149e-3, 0.0),
        (5.020892e-1, 4.316985e-2, 0.0),
        (-9.042460e-1, -4.017103e-2, 0.0),
        (-3.373058e0, 2.458538e-1, -4.850047e-3),
        (4.161981e0, 2.202737e-1, -1.718010e-2),
        (2.462523e0, 3.231308e-1, -2.281072e-2),
    ],
    'omega22': [
        (7.898524e-1, -2.114115e-2, 0.0),
        (-2.998325e-1, -1.243977e-3, 0.0),
        (7.077103e-1, 3.583907e-2, 0.0),
        (-8.946857e-1, -2.473947e-2, 0.0),
        (-2.958969e0, 2.303358e-1, -5.226562e-3),
        (4.348412e0, 1.920321e-1, -1.496557e-2),
        (2.205440e0, 2.567027e-1, -1.861359e-2),
    ],
}
_XI1_M6, _XI2_M6 = 0.8002, 0.049256  # Table 4, m = 6

# Table 5 pair parameters (beta, eps0 [meV], r_e [Angstrom]) for the pairs
# relevant to secondary and CO2/N2/O2-bearing atmospheres.
LARICCHIUTA_PAIRS = {
    ('N', 'N'): (6.61, 6.432, 3.583),
    ('O', 'O'): (6.90, 5.763, 3.423),
    ('C', 'C'): (6.69, 7.861, 3.832),
    ('C', 'N'): (6.65, 6.884, 3.717),
    ('C', 'O'): (6.78, 6.125, 3.653),
    ('N', 'O'): (6.72, 5.989, 3.507),
    ('N2', 'N2'): (8.07, 11.443, 3.829),
    ('O2', 'O2'): (8.14, 11.972, 3.780),
    ('CO', 'CO'): (8.00, 12.264, 3.889),
    ('CO2', 'CO2'): (7.75, 19.911, 4.119),
    ('N', 'CO2'): (6.90, 10.461, 3.892),
    ('O', 'CO2'): (7.19, 9.270, 3.842),
    ('N2', 'CO2'): (7.90, 14.772, 3.986),
}

# ---------------------------------------------------------------------------
# Rung 2: hydrogen, which Laricchiuta et al. do not tabulate. Zahnle et al.
# (1990, Icarus 84, 502) Eq. (30) inverts the binary diffusion parameter
# into a collision cross section, sigma_c = (3 sqrt(pi) / (16 b11))
# sqrt(2 k T / mu11) in cgs, with mu11 = m/2 the like-pair reduced mass.
# The H2-H2 self-diffusion parameter b11 = 4.96e17 T^0.75 cm^-1 s^-1 is
# recovered from the in-H2 column of Zahnle & Kasting (1986, Icarus 68, 462)
# Table I, and the H-H value scales it by 1.91, the mean of that table's
# printed in-H over in-H2 column ratios. The route gives sigma(H-H) =
# 6.4e-20 m^2 at 1e4 K with a T^-0.25 dependence by construction.
# ---------------------------------------------------------------------------
_B11_H2H2 = 4.96e17  # cm^-1 s^-1 prefactor of b11 = 4.96e17 T^0.75
_H_COLUMN_SCALE = 1.91  # ZK86 Table I in-H2 -> in-H column scaling
_M_H2_G = 2.016 * 1.66053907e-24  # g
_M_H_G = 1.008 * 1.66053907e-24  # g


def _lar_ai(beta: float, which: str) -> list[float]:
    return [c0 + c1 * beta + c2 * beta * beta for (c0, c1, c2) in _TABLE3_M6[which]]


def lar_omega_star(pair: tuple, T: float, which: str = 'omega11') -> float:
    """Reduced collision integral Omega^(l,s)* (Laricchiuta et al. 2009, Eqs. 2-3)."""
    beta, eps0_mev, _re = LARICCHIUTA_PAIRS[pair]
    a1, a2, a3, a4, a5, a6, a7 = _lar_ai(beta, which)
    x = math.log(kb * T / (eps0_mev * 1e-3 * 1.602176634e-19))
    s1 = 1.0 / (1.0 + math.exp(-2.0 * (x - a3) / a4))
    s2 = 1.0 / (1.0 + math.exp(-2.0 * (x - a6) / a7))
    return math.exp((a1 + a2 * x) * s1 + a5 * s2)


def lar_sigma2_omega(pair: tuple, T: float, which: str = 'omega11') -> float:
    """Dimensional sigma^2 Omega^(l,s)* in Angstrom^2 (Laricchiuta et al. 2009, Eq. 4)."""
    beta, _eps0, re_a = LARICCHIUTA_PAIRS[pair]
    x0 = _XI1_M6 * beta**_XI2_M6
    return (x0 * re_a) ** 2 * lar_omega_star(pair, T, which)


def lar_sigma_diff(pair: tuple, T: float) -> float:
    """Momentum-transfer cross section pi sigma^2 Omega^(1,1)*, in m^2."""
    return math.pi * lar_sigma2_omega(pair, T, 'omega11') * 1e-20


def viscosity_pure(pair: tuple, mass_gmol: float, T: float) -> float:
    """First-approximation Chapman-Enskog viscosity of a pure gas, in Pa s.

    ``eta = 2.6693e-5 sqrt(M T) / (sigma^2 Omega^(2,2)*)`` poise, with the
    molar mass in g/mol and the collision integral in Angstrom^2 (the
    standard first Chapman-Enskog approximation; see e.g. Hirschfelder,
    Curtiss & Bird 1954). Exposed as the validation route: it anchors the
    Laricchiuta transcription on measured viscosities.
    """
    s2o22 = lar_sigma2_omega(pair, T, 'omega22')
    eta_poise = 2.6693e-5 * math.sqrt(mass_gmol * T) / s2o22
    return eta_poise * 0.1


def sigma_zk90_hydrogen(species: str, T: float) -> float:
    """H-H or H2-H2 momentum-transfer cross section in m^2 (Zahnle et al. 1990, Eq. 30).

    ``sigma_c = (3 sqrt(pi) / (16 b11)) sqrt(2 k T / mu11)`` in cgs, with
    ``mu11 = m / 2`` and ``b11 = 4.96e17 T^0.75 cm^-1 s^-1`` for H2-H2,
    scaled by 1.91 for H-H. Temperature dependence T^-0.25 by construction.
    """
    if species == 'H2':
        b11 = _B11_H2H2 * T**0.75
        mu11 = _M_H2_G / 2.0
    elif species == 'H':
        b11 = _H_COLUMN_SCALE * _B11_H2H2 * T**0.75
        mu11 = _M_H_G / 2.0
    else:
        raise ValueError(f'hydrogen route only covers H and H2, got {species!r}')
    sigma_cm2 = (3.0 * math.sqrt(math.pi) / (16.0 * b11)) * math.sqrt(2.0 * kb_cgs * T / mu11)
    return sigma_cm2 * 1e-4


def sigma_geometric(species: str) -> float:
    """Last-resort geometric hard-sphere cross section pi (2 r_vdW)^2, in m^2.

    Uses the Bondi (1964) van der Waals radius; for a composite molecule
    without a tabulated radius, the largest constituent-element radius sets
    the scale. A hard sphere has no temperature dependence, so against the
    shrinking collision integrals this rung is roughly right at room
    temperature but overshoots by a factor of a few at 1e4 K (2.6 for
    atomic N), which biases the Knudsen number low and the switch toward
    hydrodynamic verdicts. The provenance class records which species sit
    on this rung so the bias stays visible.
    """
    base = species.split('_')[0]
    if base in BONDI_VDW_RADIUS_A:
        r = BONDI_VDW_RADIUS_A[base]
    else:
        r = max(BONDI_VDW_RADIUS_A.get(el, 1.5) for el in parse_formula(base))
    return math.pi * (2.0 * r * 1e-10) ** 2


def sigma_species(species: str, T: float) -> tuple[float, str]:
    """Cross section for one species with its provenance class.

    The ladder, most trusted rung first: the Laricchiuta et al. (2009)
    like-pair collision integral where tabulated; the Zahnle et al. (1990)
    diffusion-inversion route for H and H2; the geometric Bondi-radius
    hard sphere as last resort. Returns ``(sigma [m^2], provenance)`` with
    provenance one of ``'laricchiuta'``, ``'zk90-scaled'``,
    ``'geometric-vdw'``.
    """
    base = species.split('_')[0]
    if (base, base) in LARICCHIUTA_PAIRS:
        return lar_sigma_diff((base, base), T), 'laricchiuta'
    if base in ('H', 'H2'):
        return sigma_zk90_hydrogen(base, T), 'zk90-scaled'
    return sigma_geometric(base), 'geometric-vdw'


def sigma_mixture(vmr: dict[str, float], T: float) -> tuple[float, dict]:
    """Density-weighted effective cross section of a mixture.

    ``sigma_C = sum_k n_k sigma_k / n``, the mixture weighting of
    Chatterjee & Pierrehumbert (2026, ApJ 998, 236, their Eq. 25).
    ``vmr`` maps species to mole fractions (renormalized internally).
    Returns ``(sigma_C [m^2], provenance dict per species)``.
    """
    tot = sum(vmr.values())
    sig = 0.0
    prov: dict[str, str] = {}
    for sp, x in vmr.items():
        if x <= 0.0:
            continue
        s, p = sigma_species(sp, T)
        sig += (x / tot) * s
        prov[sp] = p
    return sig, prov


# ---------------------------------------------------------------------------
# The sonic-point Knudsen switch. Chatterjee & Pierrehumbert (2026,
# ApJ 998, 236) build the sonic-point Knudsen number from the Maxwell
# mean free path 1/(sqrt(2) sigma n) against the analytic sonic-point
# density scale height of their Eq. (17),
#     H_sc = (1 + gamma) r_sc / (4 + sqrt(2) sqrt(5 - 3 gamma)),
# giving their Eq. (18)
#     Kn_sc = (4 + sqrt(2) sqrt(5 - 3 gamma))
#             / (sqrt(2) (1 + gamma) sigma_C n_sc r_sc).
# A flow with Kn_sc at or below the threshold is collisional at its sonic
# point and sustains a hydrodynamic wind; above it, the gas decouples before
# reaching sonic conditions and escape is hydrostatic (Jeans-like). The
# threshold's physical band is KN_BAND above.
# ---------------------------------------------------------------------------


def mean_free_path(sigma: float, n: float) -> float:
    """Maxwell mean free path 1/(sqrt(2) sigma n), in m."""
    return 1.0 / (SQRT2 * sigma * n)


def sonic_scale_height(r_sc: float, gamma: float = 1.0) -> float:
    """Analytic sonic-point density scale height, in m (CP26 Eq. 17)."""
    return (1.0 + gamma) * r_sc / (4.0 + SQRT2 * math.sqrt(5.0 - 3.0 * gamma))


def kn_sonic(
    n_sc: float, r_sc: float, vmr: dict[str, float], T_sc: float, gamma: float = 1.0
) -> tuple[float, float, dict]:
    """Sonic-point Knudsen number (CP26 Eq. 18) for a mixture.

    ``n_sc`` is the heavy-particle number density at the sonic point
    [m^-3], ``r_sc`` the sonic radius [m], ``vmr`` the (atomized) mixture
    composition, ``T_sc`` the temperature the cross sections are evaluated
    at [K], and ``gamma`` the polytropic index (1 for an isothermal wind).
    Returns ``(Kn_sc, sigma_C [m^2], provenance dict)``.
    """
    sigma, prov = sigma_mixture(vmr, T_sc)
    kn = (4.0 + SQRT2 * math.sqrt(5.0 - 3.0 * gamma)) / (
        SQRT2 * (1.0 + gamma) * sigma * n_sc * r_sc
    )
    return kn, sigma, prov


def effective_threshold(kn_crit: float, kn_hysteresis: float, prev_regime: str | None) -> float:
    """The switch threshold, with a hysteresis window for evolutionary use.

    With no previous regime label the threshold is ``kn_crit`` (sharp
    switch). When a previous label is supplied, the threshold moves to
    ``kn_crit * kn_hysteresis`` while leaving a hydrodynamic state and to
    ``kn_crit / kn_hysteresis`` while leaving a hydrostatic state, so a
    time-stepping track cannot chatter between branches on numerical noise.
    """
    if prev_regime is None:
        return kn_crit
    if prev_regime.startswith('hydrodynamic'):
        return kn_crit * kn_hysteresis
    if prev_regime == 'hydrostatic':
        return kn_crit / kn_hysteresis
    return kn_crit
