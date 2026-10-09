"""
!!! info "`collision.py`"
    Fractional atmospheric mass loss of the target planet in a giant impact.<br>
    Author(s): Anna Grace Ulses
"""

from __future__ import annotations

from typing import Any

import numpy as np

from zephyrus.constants import G

# Fitted coefficients for the Roche et al. (2026) giant impact model.
# Source: Roche et al. (2026), arXiv:2610.06077; Zenodo doi:10.5281/zenodo.23192423.
_ROCHE2026_COEFFICIENTS: dict[str, float] = {
    # Near-field mass parameters (fit_params_NF_mass.txt, 23 values):
    'q11': -0.0011404070287601,
    'q12': 0.0005853554706874,
    'q13': 0.0004229310650559,
    'q14': 2.0,
    'q15': -0.0001681893595867,
    'q21': -2.990394733907284,
    'q22': -4.340509285222407,
    'q31': 3.680308721661617,
    'q32': -3.9046955628287696,
    'q33': 1.881415631996061,
    'q34': 2.0,
    'q35': 2.328882600831519,
    'q36': -6.797028501938182,
    'q38': -0.0411818029322558,
    'q41': 161.05720316915694,
    'q42': -7.690384546730376,
    'q43': 7.577016664697488,
    'q44': 1.0308432839505288,
    'q45': 1.033504078552506,
    'q46': -2.3556779061398494,
    'q47': -0.0012879570817664,
    'zeta5': -161.02157740225127,
    'zeta6': 8.28566684431704e-08,
    # Near-field loss parameters (fit_params_NF_loss.txt, 17 values):
    'k11': -535.6899898220761,
    'k12': 843.4184481359747,
    'k13': -847.0795889608318,
    'k14': 494.21158521587927,
    'k15': 0.0002703723957616,
    'k16': -0.0187233083484831,
    'k17': -0.3845105615920718,
    'k21': 0.0029091121113151,
    'k22': 0.0013483201004243,
    'k23': -0.0013547525798418,
    'k24': -0.0029748236386274,
    'k25': 0.0002703723957616,
    'k31': 4811.153182686202,
    'k32': -0.0364090718700348,
    'k33': -0.2382117497656167,
    'k34': -4020.819856042687,
    'k35': 0.0002703723957616,
    # Far-field loss parameters (fit_params_FF_loss.txt, 21 values):
    's11': 1020.3918108627572,
    's12': 9.294519886102522,
    's13': -0.11166497483473946,
    's14': -1030.1005018606777,
    's15': 3.262206184294514e-05,
    's16': -0.0501934230566666,
    's21': -3.8940780439779576,
    's22': -1.1017541187112085,
    's23': -0.0004728717375858641,
    's24': 4.993174134062099,
    's25': 0.000581419775592155,
    's26': 0.003429221121849461,
    's31': -3208.3991776498683,
    's32': 3208.848548705081,
    's33': 0.00038822588810121883,
    's41': -322.486418835046,
    's42': -1.9230068604833475,
    's43': 1.415457700191505,
    's44': 322.9055347383739,
    's45': 0.0001799965131025202,
    's46': -0.1320268872697112,
}

# Table C unlisted constant values (Roche et al. 2026, Table C1-C3).
_ROCHE2026_TABLE_C_CONSTANTS: dict[str, float] = {
    'q37': 1.0,
    'q48': 1.0,
    'q16': 0.0,
    'q17': 0.0,
    'q18': 0.0,
    'q23': 0.0,
    'q24': 0.0,
    'q25': 0.0,
    'q26': 0.0,
    'q27': 0.0,
    'q28': 0.0,
    'k26': 0.0,
    'k27': 0.0,
    'k36': 0.0,
    'k37': 0.0,
    's34': 0.0,
    's35': 0.0,
    's36': 0.0,
}

_ROCHE2026_PARAMS: dict[str, float] = {
    **_ROCHE2026_COEFFICIENTS,
    **_ROCHE2026_TABLE_C_CONSTANTS,
}


def _roche2026_fit(
    b: float | np.ndarray,
    gamma: float | np.ndarray,
    v_c_v_esc: float | np.ndarray,
    M_t_earth: float | np.ndarray,
    M_i_earth: float | np.ndarray,
    Q_R_prime_MJ: float | np.ndarray,
    f_atm: float | np.ndarray,
    R_ratio: float | np.ndarray,
) -> tuple[Any, Any, Any, Any]:
    r"""Evaluate the Roche et al. (2026) giant impact atmospheric mass loss model.

    Implements the scaling law of Roche et al. (2026) for fractional atmospheric
    mass loss ($X_{\text{atm}} = X_{\text{NF}} + X_{\text{FF}}$) using the
    reduced/dimensionless variables of the paper.

    Parameters
    ----------
    b : float or numpy.ndarray
        Impact parameter $\sin\beta$ in [0, 1].
    gamma : float or numpy.ndarray
        Refractory mass ratio $M_i^r / (M_i^r + M_t^r)$ in (0, 1).
    v_c_v_esc : float or numpy.ndarray
        Velocity ratio $v_c / v_{\text{esc}}$ at surface contact.
    M_t_earth : float or numpy.ndarray
        Target refractory mass in Earth masses ($M_t^r / M_{\oplus}$).
    M_i_earth : float or numpy.ndarray
        Impactor refractory mass in Earth masses ($M_i^r / M_{\oplus}$).
    Q_R_prime_MJ : float or numpy.ndarray
        Modified specific impact energy $Q'_R$ in MJ/kg.
    f_atm : float or numpy.ndarray
        Target atmospheric mass fraction $(M_t^{\text{tot}} - M_t^r) / M_t^{\text{tot}}$.
    R_ratio : float or numpy.ndarray
        Refractory radius ratio $R_i^r / R_t^r$.

    Returns
    -------
    f_NF : float or numpy.ndarray
        Near-field atmospheric mass fraction in [0, 1].
    X_NF : float or numpy.ndarray
        Near-field atmospheric mass loss fraction in [0, f_NF].
    X_FF : float or numpy.ndarray
        Far-field atmospheric mass loss fraction in [0, 1 - f_NF].
    X_atm : float or numpy.ndarray
        Total atmospheric mass loss fraction in [0, 1].

    References
    ----------
    1. Roche M.J., Lock S.J., Carter P.J., Leinhardt Z.M. (2026).
       Giant impacts preferentially remove low-mass atmospheres:
       a generalised scaling law for impact-driven atmospheric loss.
       arXiv:2610.06077 (accepted, ApJL); Zenodo doi:10.5281/zenodo.23192423.
    """
    p = _ROCHE2026_PARAMS

    # Near-field envelope fraction zeta_i terms (Eq. 7, general form)
    def _calc_zeta(i: int) -> Any:
        return (
            p[f'q{i}1']
            + p[f'q{i}2'] * b
            + p[f'q{i}3'] * b ** p[f'q{i}4']
            + p[f'q{i}5'] * f_atm
            + p[f'q{i}6'] * f_atm**2
            + p[f'q{i}7'] * M_t_earth ** p[f'q{i}8']
        )

    z1 = _calc_zeta(1)
    z2 = _calc_zeta(2)
    z3 = _calc_zeta(3)
    z4 = _calc_zeta(4)

    # Near-field atmospheric mass fraction f_NF (Eq. 6)
    f_nf = z4 / (1.0 + p['zeta6'] * np.exp(z3 * (R_ratio - z2))) ** z1 + p['zeta5']
    f_nf = np.clip(f_nf, 0.0, 1.0)
    f_ff = 1.0 - f_nf

    # Near-field loss function xi_i (Eq. 8, Eq. 9, general form)
    def _calc_xi_i(i: int, v_rel: Any) -> Any:
        return (
            p[f'k{i}1']
            + p[f'k{i}2'] * gamma
            + p[f'k{i}3'] * (gamma + 0.05) ** 2
            + p[f'k{i}4'] * v_rel ** p[f'k{i}5']
            + p[f'k{i}6'] * (M_t_earth + 1.0)
            + p[f'k{i}7'] * np.log10(f_atm)
        )

    def _calc_xi(v_rel: Any) -> Any:
        x1 = _calc_xi_i(1, v_rel)
        x2 = _calc_xi_i(2, v_rel)
        x3 = _calc_xi_i(3, v_rel)
        return f_nf * (x1 - x2 * (b + x3) ** 2)

    x_nf = _calc_xi(v_c_v_esc)
    lo = np.maximum(0.0, _calc_xi(1.0))
    x_nf = np.clip(x_nf, lo, f_nf)

    # Far-field loss function psi_i (Eq. 10, Eq. 11, general form)
    def _calc_psi_i(i: int) -> Any:
        return (
            p[f's{i}1']
            + p[f's{i}2'] * (gamma + 0.05) ** p[f's{i}3']
            + p[f's{i}4'] * (M_t_earth + 0.05) ** p[f's{i}5']
            + p[f's{i}6'] * np.log10(f_atm)
        )

    p1 = _calc_psi_i(1)
    p2 = _calc_psi_i(2)
    p3 = _calc_psi_i(3)
    p4 = _calc_psi_i(4)

    m_ratio = M_i_earth / M_t_earth
    x_ff = np.clip(
        f_ff * (p1 * np.exp(-p2 * (Q_R_prime_MJ * (1.0 + m_ratio) * (1.0 - b) ** p4)) + p3),
        0.0,
        f_ff,
    )
    x_atm = np.clip(x_nf + x_ff, 0.0, 1.0)
    return f_nf, x_nf, x_ff, x_atm


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
