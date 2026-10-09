"""
!!! info "`collision.py`"
    Fractional atmospheric mass loss of the target planet in a giant impact.<br>
    Author(s): Anna Grace Ulses
"""

from __future__ import annotations

import types
from dataclasses import dataclass
from typing import Any

import numpy as np

from zephyrus.constants import G, c
from zephyrus.planets_parameters import Me


@dataclass
class ImpactLossResult:
    r"""Atmospheric mass loss outcome from a giant impact.

    Parameters
    ----------
    law : str
        Name of the scaling law used ('kegerreis2020' or 'roche2026').
    fraction : float
        Fraction of the target atmosphere lost, in [0, 1].
    flags : tuple of str
        Names of parameters that fall outside the calibrated range.
        Empty when all parameters lie within range. Active stability
        clamps are recorded in diagnostics['clamped'].
    diagnostics : dict of str to Any
        Diagnostic quantities computed during evaluation.
        For 'kegerreis2020': 'v_esc', 'v_ratio', 'gamma', 'mass_ratio', 'f_M'.
        For 'roche2026': 'f_atm', 'M_t_earth', 'gamma', 'b', 'R_ratio',
        'v_ratio', 'v_esc', 'Q_R_prime' (in MJ/kg), 'clamped', 'f_NF',
        'X_NF', 'X_FF'.
    """

    law: str
    fraction: float
    flags: tuple[str, ...]
    diagnostics: dict[str, Any]


# Calibrated parameter ranges for the Roche et al. (2026) scaling law:
# Section 4.1 (p. 7) for f_atm, M_t_earth, gamma, b, and v_ratio;
# Section 3.2 (p. 5) and the authors' fitting data (scaling_law.csv) for R_ratio
# (union of initialized SPH planet suite [0.001, 1.0] and impact suite [0.498, 1.015]).
ROCHE2026_FITTED_RANGE: types.MappingProxyType[str, tuple[float, float]] = (
    types.MappingProxyType(
        {
            'f_atm': (0.01, 0.20),
            'M_t_earth': (0.35, 5.0),
            'gamma': (0.1, 0.5),
            'b': (0.0, 0.9),
            'R_ratio': (0.001, 1.015),
            'v_ratio': (1.0, 3.0),
        }
    )
)
"""Empirical parameter domain for the Roche et al. (2026) scaling law.

Keys:
- ``f_atm``: initial envelope mass fraction [dimensionless]
- ``M_t_earth``: refractory target mass in Earth masses [M_E]
- ``gamma``: impactor mass fraction M_i / (M_t^r + M_i) [dimensionless]
- ``b``: dimensionless impact parameter at mantle contact [dimensionless]
- ``R_ratio``: impactor-to-target refractory radius ratio R_i / R_t [dimensionless]
- ``v_ratio``: impact speed ratio v_c / v_esc at mantle contact [dimensionless]

Flags are recorded when an evaluation parameter falls outside its bounds by
more than the 1% relative tolerance (_ROCHE2026_RANGE_RTOL = 0.01). The lower
v_ratio bound (v_ratio < 1.0) is not flagged because the near-field velocity
floor models sub-escape collisions.
"""

# Relative tolerance on empirical boundary flags to cover grid-point rounding.
_ROCHE2026_RANGE_RTOL: float = 0.01

# Numerical stability ranges used for clamping inputs during fit evaluation.
_ROCHE2026_STABLE_RANGE: dict[str, tuple[float, float]] = {
    'f_atm': (1.0e-6, 0.4),
    'M_t_earth': (1.0e-3, 10.0),
    'gamma': (1.0e-3, 1.0),
}


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

# Terms fixed in Table C1-C3 and absent from the data files.
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
    *,
    b: float | np.ndarray,
    gamma: float | np.ndarray,
    v_c_v_esc: float | np.ndarray,
    M_t_earth: float | np.ndarray,
    mass_ratio: float | np.ndarray,
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
    mass_ratio : float or numpy.ndarray
        Refractory mass ratio $M_i^r / M_t^r$ (dimensionless, unclamped by design).
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
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
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

    one_minus_b = np.maximum(0.0, 1.0 - b)
    with np.errstate(invalid='ignore', divide='ignore'):
        safe_base = np.where(one_minus_b > 0.0, one_minus_b, 1.0)
        geom_factor = np.where(one_minus_b > 0.0, safe_base**p4, 0.0)

    arg = Q_R_prime_MJ * (1.0 + mass_ratio) * geom_factor
    exp_arg = np.minimum(-p2 * arg, 700.0)
    with np.errstate(over='ignore', invalid='ignore'):
        raw_ff = f_ff * (p1 * np.exp(exp_arg) + p3)
        x_ff = np.where(f_ff == 0.0, 0.0, np.clip(raw_ff, 0.0, f_ff))
    if np.ndim(x_ff) == 0:
        x_ff = float(x_ff)
    x_atm = np.clip(x_nf + x_ff, 0.0, 1.0)
    if np.ndim(x_atm) == 0:
        x_atm = float(x_atm)
    return f_nf, x_nf, x_ff, x_atm


def _interacting_mass_fraction_kegerreis(
    R_t: float,
    R_i: float,
    rho_t: float,
    rho_i: float,
    b: float,
) -> float:
    """Fractional interacting mass f_M (Kegerreis et al. 2020, Eqn. B1).

    Density-weighted spherical caps of common height d, clamped to [0, 1]
    because the linearised caps can leave the interval outside the fitted
    geometry. At equal bulk densities this reduces exactly to the
    interacting volume f_V of their Eqn. B2.
    """
    d = (R_t + R_i) * (1.0 - b)
    v_t_cap = np.pi / 3.0 * d**2 * (3.0 * R_t - d)
    v_i_cap = np.pi / 3.0 * d**2 * (3.0 * R_i - d)
    v_t_full = 4.0 / 3.0 * np.pi * R_t**3
    v_i_full = 4.0 / 3.0 * np.pi * R_i**3
    f_m = (rho_t * v_t_cap + rho_i * v_i_cap) / (rho_t * v_t_full + rho_i * v_i_full)
    return float(np.clip(f_m, 0.0, 1.0))


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
    if not (0.0 <= v_c < c and np.isfinite(v_c)):
        raise ValueError(
            f'Collision speed v_c must be non-negative, sub-luminal (< c), and finite, got {v_c!r}'
        )

    # Mutual escape speed of the pair at contact
    v_esc = np.sqrt((2.0 * G * (M_t + M_i)) / (R_t + R_i))

    # Fractional interacting mass f_M (Kegerreis et al. 2020, Eqn. B1).
    f_m = _interacting_mass_fraction_kegerreis(R_t, R_i, rho_t, rho_i, b)

    m_tot = M_i + M_t
    bracket = (v_c / v_esc) ** 2 * (M_i / m_tot) ** 0.5 * (rho_i / rho_t) ** 0.5 * f_m

    # Fractional atmospheric mass loss of the target, capped at 1 for
    # total erosion (Kegerreis et al. 2020, Eqn. 1)
    x_loss = 0.64 * bracket**0.65
    return min(max(x_loss, 0.0), 1.0)


def _as_floats(**kwargs: Any) -> dict[str, float]:
    """Convert scalar arguments to float, rejecting sequences, strings, and overflow."""
    floats: dict[str, float] = {}
    for name, val in kwargs.items():
        if (
            isinstance(val, (str, bytes))
            or getattr(val, 'ndim', 0) > 0
            or isinstance(val, (list, tuple))
        ):
            raise TypeError(f'{name} must be a scalar numeric value, got {type(val).__name__}')
        try:
            floats[name] = float(val)
        except OverflowError as exc:
            raise ValueError(f'{name} exceeds floating-point range') from exc
    return floats


def _check_strictly_positive(**kwargs: float) -> None:
    """Ensure values are strictly positive and finite."""
    for name, val in kwargs.items():
        if not (val > 0.0 and np.isfinite(val)):
            raise ValueError(f'{name} must be strictly positive and finite, got {val!r}')


def mutual_escape_speed(
    M_1: float,
    M_2: float,
    R_1: float,
    R_2: float,
) -> float:
    r"""Mutual escape speed of two spherical bodies at contact [m/s].

    Parameters
    ----------
    M_1 : float
        Mass of the first body [kg].
    M_2 : float
        Mass of the second body [kg].
    R_1 : float
        Radius of the first body [m].
    R_2 : float
        Radius of the second body [m].

    Returns
    -------
    float
        Mutual escape speed $\sqrt{2 G (M_1 + M_2) / (R_1 + R_2)}$ [m/s].

    Raises
    ------
    TypeError
        If any input is a string, bytes, or non-scalar sequence/array.
    ValueError
        If any mass or radius is not strictly positive and finite, or if
        floating-point conversion overflows.
    """
    vals = _as_floats(M_1=M_1, M_2=M_2, R_1=R_1, R_2=R_2)
    _check_strictly_positive(**vals)
    return float(np.sqrt(2.0 * G * (vals['M_1'] + vals['M_2']) / (vals['R_1'] + vals['R_2'])))


def specific_impact_energy(
    v_c: float,
    M_i: float,
    M_t_tot: float,
    R_i: float,
    R_t: float,
    b: float,
) -> float:
    r"""Modified specific impact energy Q'_R in MJ/kg.

    Parameters
    ----------
    v_c : float
        Impact velocity at surface contact [m/s].
    M_i : float
        Impactor mass [kg].
    M_t_tot : float
        Target mass including its atmosphere ($M_\mathrm{t}^\mathrm{tot}$ in
        Roche et al. 2026, Eq. 1 and $Q'_\mathrm{R}$) [kg].
    R_i : float
        Impactor refractory radius [m].
    R_t : float
        Target refractory radius [m].
    b : float
        Impact parameter $\sin\beta$ in [0, 1].

    Returns
    -------
    float
        Modified specific impact energy $Q'_R$ in MJ/kg.

    Raises
    ------
    TypeError
        If any input is a string, bytes, or non-scalar sequence/array.
    ValueError
        If b is outside [0, 1], if masses or radii are not strictly positive
        and finite, or if v_c is negative, non-finite, or >= c.

    References
    ----------
    1. Roche, M. J., Lock, S. J., Dou, J., Carter, P. J., Kegerreis, J. A.,
       & Leinhardt, Z. M. (2025), "Atmospheric Loss during Giant Impacts:
       Mechanisms and Scaling of Near- and Far-field Loss", The Planetary
       Science Journal, 6, 149, doi:10.3847/PSJ/add929, Eqns. 13 to 18.
    2. Leinhardt, Z. M., & Stewart, S. T. (2012), "Collisions between
       gravity-dominated bodies. I. Outcome regimes and scaling laws",
       ApJ 745, 79, doi:10.1088/0004-637X/745/1/79.
    """
    vals = _as_floats(v_c=v_c, M_i=M_i, M_t_tot=M_t_tot, R_i=R_i, R_t=R_t, b=b)
    _check_strictly_positive(
        M_i=vals['M_i'],
        M_t_tot=vals['M_t_tot'],
        R_i=vals['R_i'],
        R_t=vals['R_t'],
    )
    v_c = vals['v_c']
    M_i = vals['M_i']
    M_t_tot = vals['M_t_tot']
    R_i = vals['R_i']
    R_t = vals['R_t']
    b = vals['b']

    if not (0.0 <= b <= 1.0 and np.isfinite(b)):
        raise ValueError(f'Impact parameter b must be in [0, 1], got {b!r}')
    if not (0.0 <= v_c < c and np.isfinite(v_c)):
        raise ValueError(
            f'Collision speed v_c must be non-negative, sub-luminal (< c), and finite, got {v_c!r}'
        )

    impact_param = (R_t + R_i) * b
    if impact_param + R_i <= R_t:
        alpha = 1.0
    else:
        interact_len = R_t + R_i - impact_param
        alpha = (3.0 * R_i * interact_len**2 - interact_len**3) / (4.0 * R_i**3)

    m_tot = M_t_tot + M_i
    mu_alpha = (alpha * M_i) / (1.0 + (alpha * M_i) / M_t_tot)
    q_r_prime = mu_alpha * v_c**2 / (2.0 * m_tot) / 1.0e6
    return float(q_r_prime)


def _eval_roche2026(
    v_c: float,
    M_i: float,
    M_t: float,
    R_i: float,
    R_t: float,
    b: float,
    f_atm: float,
) -> tuple[float, tuple[str, ...], dict[str, Any]]:
    """Validate, clamp, and evaluate the Roche et al. (2026) scaling law."""
    vals = _as_floats(v_c=v_c, M_i=M_i, M_t=M_t, R_i=R_i, R_t=R_t, b=b, f_atm=f_atm)
    _check_strictly_positive(
        M_i=vals['M_i'],
        M_t=vals['M_t'],
        R_i=vals['R_i'],
        R_t=vals['R_t'],
    )
    v_c = vals['v_c']
    M_i = vals['M_i']
    M_t = vals['M_t']
    R_i = vals['R_i']
    R_t = vals['R_t']
    b = vals['b']
    f_atm = vals['f_atm']

    if not (0.0 <= b <= 1.0 and np.isfinite(b)):
        raise ValueError(f'Impact parameter b must be in [0, 1], got {b!r}')
    if not (0.0 <= v_c < c and np.isfinite(v_c)):
        raise ValueError(
            f'Collision speed v_c must be non-negative, sub-luminal (< c), and finite, got {v_c!r}'
        )
    if not (0.0 <= f_atm < 1.0 and np.isfinite(f_atm)):
        raise ValueError(f'f_atm must be in [0, 1) and finite, got {f_atm!r}')

    m_t_tot = M_t / (1.0 - f_atm)
    if not np.isfinite(m_t_tot):
        raise ValueError(f'Total target mass M_t / (1 - f_atm) must be finite, got {m_t_tot!r}')

    gamma = M_i / (M_i + M_t)
    r_ratio = R_i / R_t
    m_t_earth = M_t / Me
    v_esc = mutual_escape_speed(m_t_tot, M_i, R_t, R_i)
    v_ratio = v_c / v_esc
    q_r_prime = specific_impact_energy(v_c, M_i, m_t_tot, R_i, R_t, b)

    clamped: dict[str, float] = {}
    diag: dict[str, Any] = {
        'f_atm': f_atm,
        'M_t_earth': m_t_earth,
        'gamma': gamma,
        'b': b,
        'R_ratio': r_ratio,
        'v_ratio': v_ratio,
        'v_esc': v_esc,
        'Q_R_prime': q_r_prime,
        'clamped': clamped,
        'f_NF': 0.0,
        'X_NF': 0.0,
        'X_FF': 0.0,
    }

    if f_atm == 0.0:
        return 0.0, (), diag

    eval_params = {'f_atm': f_atm, 'M_t_earth': m_t_earth, 'gamma': gamma}
    for name, val in eval_params.items():
        c_lo, c_hi = _ROCHE2026_STABLE_RANGE[name]
        if val < c_lo:
            clamped[name] = c_lo
            eval_params[name] = c_lo
        elif val > c_hi:
            clamped[name] = c_hi
            eval_params[name] = c_hi

    flags: list[str] = []
    param_vals = (
        ('f_atm', f_atm),
        ('M_t_earth', m_t_earth),
        ('gamma', gamma),
        ('b', b),
        ('R_ratio', r_ratio),
        ('v_ratio', v_ratio),
    )
    for name, val in param_vals:
        f_lo, f_hi = ROCHE2026_FITTED_RANGE[name]
        if name == 'v_ratio':
            if val > f_hi * (1.0 + _ROCHE2026_RANGE_RTOL):
                flags.append(name)
        elif val < f_lo * (1.0 - _ROCHE2026_RANGE_RTOL) or val > f_hi * (
            1.0 + _ROCHE2026_RANGE_RTOL
        ):
            flags.append(name)

    f_nf, x_nf, x_ff, x_atm = _roche2026_fit(
        b=b,
        gamma=eval_params['gamma'],
        v_c_v_esc=v_ratio,
        M_t_earth=eval_params['M_t_earth'],
        mass_ratio=M_i / M_t,
        Q_R_prime_MJ=q_r_prime,
        f_atm=eval_params['f_atm'],
        R_ratio=r_ratio,
    )

    if not (
        np.isfinite(x_atm) and np.isfinite(x_nf) and np.isfinite(x_ff) and np.isfinite(f_nf)
    ):
        raise ValueError(f'Roche scaling law produced non-finite result: x_atm={x_atm!r}')

    diag['f_NF'] = float(f_nf)
    diag['X_NF'] = float(x_nf)
    diag['X_FF'] = float(x_ff)
    return float(x_atm), tuple(flags), diag


def mass_loss_roche2026(
    v_c: float,
    M_i: float,
    M_t: float,
    R_i: float,
    R_t: float,
    b: float,
    f_atm: float,
) -> float:
    r"""Fractional atmospheric mass loss of the target in a giant impact.

    Implements the scaling law of Roche et al. (2026).

    Parameters
    ----------
    v_c : float
        Collision speed at first contact [m/s].
    M_i : float
        Impactor refractory mass, excluding atmosphere [kg].
    M_t : float
        Target refractory mass, excluding atmosphere [kg].
    R_i : float
        Impactor refractory radius, at base of atmosphere [m].
    R_t : float
        Target refractory radius, at base of atmosphere [m].
    b : float
        Dimensionless impact parameter $\sin\beta$ in [0, 1].
    f_atm : float
        Target atmospheric mass fraction $M_{\text{atm}} / (M_t + M_{\text{atm}})$ in [0, 1).

    Returns
    -------
    float
        Fractional atmospheric mass loss of the target in [0, 1].

    Raises
    ------
    TypeError
        If any input is a string, bytes, or non-scalar sequence/array.
    ValueError
        If inputs violate physical domain constraints.

    Notes
    -----
    Stability clamps apply to inputs outside the numerical stability bounds.
    This function returns only the loss fraction and does not report
    diagnostic flags; use ``impact_loss`` if validity flags or diagnostics
    are needed. Inputs are scalar; arrays are not supported.

    References
    ----------
    1. Roche M.J., Lock S.J., Carter P.J., Leinhardt Z.M. (2026).
       Giant impacts preferentially remove low-mass atmospheres:
       a generalised scaling law for impact-driven atmospheric loss.
       arXiv:2610.06077 (accepted, ApJL); Zenodo doi:10.5281/zenodo.23192423.
    """
    frac, _, _ = _eval_roche2026(v_c, M_i, M_t, R_i, R_t, b, f_atm)
    return frac


def impact_loss(
    law: str,
    *,
    v_c: float,
    M_i: float,
    M_t: float,
    R_i: float,
    R_t: float,
    b: float,
    rho_i: float | None = None,
    rho_t: float | None = None,
    f_atm: float | None = None,
) -> ImpactLossResult:
    r"""Evaluate atmospheric erosion scaling law for a giant impact.

    Parameters
    ----------
    law : str
        Scaling law name ('kegerreis2020' or 'roche2026').
    v_c : float
        Collision speed at first contact [m/s].
    M_i : float
        Impactor mass, excluding atmosphere [kg].
    M_t : float
        Target mass, excluding atmosphere [kg].
    R_i : float
        Impactor radius, at base of atmosphere [m].
    R_t : float
        Target radius, at base of atmosphere [m].
    b : float
        Dimensionless impact parameter $\sin\beta$ in [0, 1].
    rho_i : float, optional
        Impactor bulk density [kg/m^3] (required for 'kegerreis2020').
    rho_t : float, optional
        Target bulk density [kg/m^3] (required for 'kegerreis2020').
    f_atm : float, optional
        Target atmospheric mass fraction in [0, 1) (required for 'roche2026').

    Returns
    -------
    ImpactLossResult
        Container holding law name, loss fraction, flags, and diagnostics.

    Raises
    ------
    TypeError
        If any input is a string, bytes, or non-scalar sequence/array.
    ValueError
        If law is unsupported or required law-specific inputs are missing.

    Notes
    -----
    Arguments that the selected scaling law does not use are ignored (e.g.
    bulk densities for 'roche2026', or f_atm for 'kegerreis2020'). An airless
    target (f_atm = 0) returns no diagnostic flags.
    """
    if law == 'kegerreis2020':
        if rho_i is None or rho_t is None:
            missing = [k for k, v in (('rho_i', rho_i), ('rho_t', rho_t)) if v is None]
            raise ValueError(f'kegerreis2020 requires {", ".join(missing)}')
        vals = _as_floats(
            v_c=v_c,
            M_i=M_i,
            M_t=M_t,
            rho_i=rho_i,
            rho_t=rho_t,
            R_i=R_i,
            R_t=R_t,
            b=b,
        )
        frac = mass_loss(
            vals['v_c'],
            vals['M_i'],
            vals['M_t'],
            vals['rho_i'],
            vals['rho_t'],
            vals['R_i'],
            vals['R_t'],
            vals['b'],
        )
        v_esc = mutual_escape_speed(vals['M_t'], vals['M_i'], vals['R_t'], vals['R_i'])
        v_ratio = vals['v_c'] / v_esc
        gamma = vals['M_i'] / (vals['M_t'] + vals['M_i'])
        f_m = _interacting_mass_fraction_kegerreis(
            vals['R_t'], vals['R_i'], vals['rho_t'], vals['rho_i'], vals['b']
        )
        diag = {
            'v_esc': float(v_esc),
            'v_ratio': float(v_ratio),
            'gamma': float(gamma),
            'mass_ratio': float(vals['M_i'] / vals['M_t']),
            'f_M': float(f_m),
        }
        return ImpactLossResult(law=law, fraction=frac, flags=(), diagnostics=diag)

    if law == 'roche2026':
        if f_atm is None:
            raise ValueError('roche2026 requires f_atm')
        frac, flags, diag = _eval_roche2026(v_c, M_i, M_t, R_i, R_t, b, f_atm)
        return ImpactLossResult(law=law, fraction=frac, flags=flags, diagnostics=diag)

    raise ValueError(f"Unknown scaling law {law!r}, expected 'kegerreis2020' or 'roche2026'")
