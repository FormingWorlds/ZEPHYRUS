"""
!!! info "`diffusion.py`"
    Binary diffusion parameters b = n D with provenance, and mixture rules.<br>
    Authors: Mara Attia, Viesturs Strelcs
"""

from __future__ import annotations

import math

import numpy as np

from zephyrus.composition import BONDI_VDW_RADIUS_A, ELEMENT_AMU, species_mass_amu
from zephyrus.constants import amu, kb

# The library convention is b = n D in cm^-1 s^-1 (the form of Zahnle &
# Kasting 1986, Icarus 68, 462), fitted as b = b1000 (T/1000 K)^s. Public
# helpers that feed SI code convert at the boundary (1 cm^-1 s^-1 =
# 100 m^-1 s^-1).

KB_SI = kb  # J/K
AMU_G = amu * 1e3  # g

# Molecular masses in amu, needed as scaling anchors, plus atomic deuterium.
MOLECULE_AMU = {
    'H2': 2.01588,
    'D': 2.0141,
    'O2': 31.998,
    'N2': 28.014,
    'CO2': 44.0095,
    'H2O': 18.0153,
    'CH4': 16.043,
}
ALL_MASS = dict(ELEMENT_AMU, **MOLECULE_AMU)

# ---------------------------------------------------------------------------
# Kinetic diameters, pm.
# ---------------------------------------------------------------------------
# Printed by Zahnle & Kasting (2023, GeCoA 361, 228) in the last row of their
# Table 2; the standard published kinetic diameters for the molecules, with
# the atoms H and O assigned the value of the neighboring noble gas.
D_ZK23 = {
    'H2': 289.0,
    'He': 260.0,
    'Ne': 275.0,
    'Ar': 340.0,
    'CO2': 330.0,
    'O2': 346.0,
    'H2O': 265.0,
    'CH4': 380.0,
    'N2': 364.0,
    'H': 260.0,
    'D': 265.0,
    'O': 275.0,
}
# Kr and Xe from the same standard compilation. Consistency check: the van
# der Waals scaling rule below independently gives 365 and 391 pm against
# these 360 and 396, agreeing to 1.4 percent.
D_STANDARD_EXTRA = {'Kr': 360.0, 'Xe': 396.0}

# Elements with no printed kinetic diameter anywhere (C, N, S, and the rock
# formers) get one by scaling the printed atomic-O entry with the ratio of
# Bondi (1964) van der Waals radii. Taking the atomic diameter to be the van
# der Waals size has published precedent in exactly this application (Ito &
# Ikoma 2021, MNRAS 502, 750, their Eq. 36). The rule is checkable: applied
# to He, Ne, and Ar it reproduces the printed diameters to 3 percent.
_VDW_SCALED = ('C', 'N', 'S', 'Na', 'Mg', 'Si', 'Fe')


def diameters() -> dict[str, float]:
    """Kinetic diameters in pm: printed values plus the van der Waals rule."""
    d = dict(D_ZK23)
    d.update(D_STANDARD_EXTRA)
    for s in _VDW_SCALED:
        d[s] = D_ZK23['O'] * BONDI_VDW_RADIUS_A[s] / BONDI_VDW_RADIUS_A['O']
    return d


# Rock-forming species carry two standing warnings that no coefficient
# improves away: every pair involving Na, Mg, Si, or Fe is a scaling on an
# estimated diameter (no measured coefficient exists in any compilation for
# these pairs), and at the temperatures where rock vapor exists Na and Mg
# ionize readily while these are neutral-gas coefficients. The Fe radius
# spans about 20 percent across published compilations and the Mg radius is
# a flagged outlier of Bondi's own table, so those two rungs are soft beyond
# their provenance class (see the notes in composition.BONDI_VDW_RADIUS_A).
ROCK_FORMERS = ('Na', 'Mg', 'Si', 'Fe')

# ---------------------------------------------------------------------------
# Source 1: Zahnle & Kasting (2023, GeCoA 361, 228) Table 2.
# (i, j): (b at 1000 K in cm^-1 s^-1, class, source column as printed).
# Classes: 'M' traces to Marrero & Mason (1972) measurements; 'E' is the
# authors' own scaling estimate from the named analog pairs. All rows carry
# the fitted exponent 0.75.
# ---------------------------------------------------------------------------
ZK23_EXPONENT = 0.75
ZK23_TABLE2 = {
    ('H', 'H'): (1.3e20, 'E', 'H-H2, H2-D2, H2-Ne'),
    ('H', 'O'): (9.0e19, 'E', 'H2-O2, H2-Ne, O-He, H-Ar'),
    ('H', 'O2'): (6.5e19, 'E', 'H2-O2, H2-Ar, H-Ar'),
    ('H', 'CO2'): (6.0e19, 'E', 'H2-CO2, He-CO2, H-Ar'),
    ('O', 'O'): (3.0e19, 'E', 'O-Ar, O-O2, Ne-N2'),
    ('CO2', 'O'): (1.6e19, 'E', 'H2O-CO2, Ne-CO2, O-Ar'),
    ('H2O', 'CO2'): (1.56e19, 'M', 'Marrero and Mason (1972)'),
    ('H2O', 'O2'): (1.59e19, 'M', 'Marrero and Mason (1972)'),
    ('CO2', 'O2'): (1.0e19, 'M', 'Marrero and Mason (1972)'),
    ('H', 'D'): (1.1e20, 'E', 'b11, H2-D2'),
    ('H', 'He'): (1.6e20, 'M', 'Marrero and Mason (1972)'),
    ('H', 'Ne'): (9.3e19, 'E', 'Ne-H2, H-Ar, He-Ne'),
    ('H', 'N2'): (6.5e19, 'E', 'H2-N2'),
    ('H', 'Ar'): (6.5e19, 'M', 'Marrero and Mason (1972)'),
    ('O', 'D'): (6.5e19, 'E', 'b11, H2-D2'),
    ('O', 'He'): (6.0e19, 'M', 'Marrero and Mason (1972)'),
    ('O', 'Ne'): (3.0e19, 'E', 'O-He, O-Ar, H2-Ne'),
    ('O', 'Ar'): (1.8e19, 'M', 'Marrero and Mason (1972)'),
    ('O', 'N2'): (2.0e19, 'E', 'CH4-N2, air-H2O, O-O2'),
    # Printed as 4.3e20; entered as 4.3e19. The printed value sits a full
    # decade above the H2-D2 anchor its own source column names, while
    # 4.3e19 is consistent with it and with the neighboring D rows, so the
    # printed exponent is treated as a typographical error.
    ('CO2', 'D'): (4.3e19, 'E', 'b11, H2-D2 [printed 4.3e20, misprint]'),
    ('CO2', 'He'): (3.56e19, 'M', 'Marrero and Mason (1972)'),
    ('CO2', 'Ne'): (1.62e19, 'M', 'Marrero and Mason (1972)'),
    ('CO2', 'Ar'): (1.0e19, 'M', 'Marrero and Mason (1972)'),
    ('CO2', 'N2'): (1.04e19, 'M', 'Marrero and Mason (1972)'),
}

# ---------------------------------------------------------------------------
# Source 2: Sasaki & Nakazawa (1988, EPSL 89, 323) Table 1, which tabulates
# f_ij = P D_ij for the noble gases against H2 and He at 100, 1000, and
# 10000 K, based on Marrero & Mason (1972). Their Eq. (6) is D_ij = f_ij/P =
# f_ij/(n kB T), so b = n D = f/(kB T), and reading f in SI units
# (Pa m^2 s^-1) gives b [cm^-1 s^-1] = f/(kB T)/100. That unit reading is
# verified, not assumed: it reproduces five in-H2 entries of Zahnle & Kasting
# (1986) Table I (He, Ne, Ar, Kr, Xe) to between 0.02 and 3 percent, two
# independent compilations both citing Marrero & Mason; no other unit choice
# comes within two orders of magnitude (see the companion tests).
# ---------------------------------------------------------------------------
SN88_TABLE1 = {
    ('He', 'H2'): {100: 2.4, 1000: 1.3e2, 10000: 7.6e3},
    ('Ne', 'H2'): {100: 1.7, 1000: 9.4e1, 10000: 5.1e3},
    ('Ne', 'He'): {100: 1.8, 1000: 8.7e1, 10000: 4.8e3},
    ('Ar', 'H2'): {100: 1.0, 1000: 7.1e1, 10000: 3.9e3},
    ('Ar', 'He'): {100: 1.2, 1000: 6.2e1, 10000: 3.6e3},
    ('Kr', 'H2'): {100: 9.2e-1, 1000: 6.1e1, 10000: 3.5e3},
    ('Kr', 'He'): {100: 1.0, 1000: 5.3e1, 10000: 3.2e3},
    ('Xe', 'H2'): {100: 8.3e-1, 1000: 5.0e1, 10000: 2.6e3},
    ('Xe', 'He'): {100: 7.9e-1, 1000: 4.4e1, 10000: 2.7e3},
}


def b_from_sn88(f_value: float, T_of_table: float) -> float:
    """Convert a Sasaki & Nakazawa (1988) Table 1 entry to b = n D in cm^-1 s^-1."""
    return f_value / (KB_SI * T_of_table) / 100.0


def sn88_fit(species: str, partner: str, t_lo: float = 1000, t_hi: float = 10000):
    """Fit b = A T^s to two Sasaki & Nakazawa temperature points.

    Defaults to the 1000 and 10000 K points, both above the 242 K validity
    floor their footnote attaches to the coldest column. Returns ``(A, s)``.
    """
    row = SN88_TABLE1[(species, partner)]
    b_lo = b_from_sn88(row[t_lo], t_lo)
    b_hi = b_from_sn88(row[t_hi], t_hi)
    s = math.log(b_hi / b_lo) / math.log(t_hi / t_lo)
    return b_lo / t_lo**s, s


# ---------------------------------------------------------------------------
# Cross-check set: Zahnle & Kasting (1986, Icarus 68, 462) Table I, b = A T^s
# in cm^-1 s^-1. Not a source of library rows; six pairs appear both here and
# in the 2023 compilation, and their agreement bounds the transcriptions.
# ---------------------------------------------------------------------------
ZK86_TABLE1 = {
    ('He', 'H2'): (5.23e17, 0.75),
    ('He', 'H'): (1.04e18, 0.732),
    ('He', 'O'): (3.44e17, 0.75),
    ('O', 'H2'): (3.0e17, 0.75),
    ('O', 'H'): (4.8e17, 0.75),
    ('Ne', 'H2'): (4.37e17, 0.731),
    ('Ne', 'H'): (7.9e17, 0.731),
    ('Ne', 'O'): (1.5e17, 0.75),
    ('Ar', 'H2'): (2.81e17, 0.75),
    ('Ar', 'H'): (1.06e18, 0.597),
    ('Ar', 'O'): (5.61e16, 0.841),
    ('Kr', 'H2'): (2.3e17, 0.76),
    ('Kr', 'H'): (4.1e17, 0.76),
    ('Kr', 'O'): (4.3e16, 0.841),
    ('Xe', 'H2'): (2.7e17, 0.712),
    ('Xe', 'H'): (4.9e17, 0.712),
    ('Xe', 'O'): (3.5e16, 0.841),
    ('N2', 'H2'): (2.65e17, 0.75),
    ('N2', 'H'): (6.5e17, 0.70),
    ('N2', 'O'): (9.7e16, 0.774),
    ('CO2', 'H2'): (2.3e17, 0.75),
    ('CO2', 'H'): (8.4e17, 0.60),
    ('CO2', 'O'): (7.86e16, 0.776),
    ('H2O', 'H2'): (2.7e17, 0.75),
    ('H2O', 'H'): (6.6e17, 0.70),
    ('H2O', 'O'): (1.06e17, 0.774),
}


def b_zk86(species: str, partner: str, T: float) -> float:
    """Zahnle & Kasting (1986) Table I fit b = A T^s, in cm^-1 s^-1."""
    a, s = ZK86_TABLE1[(species, partner)]
    return a * T**s


# ---------------------------------------------------------------------------
# Minor species in molecular backgrounds: fits b = A T^s in cm^-1 s^-1
# carried as tabulated by the diffusion-limited escape branch (compiled by
# Viesturs Strelcs); individual rows trace to the standard compilations of
# measured binary diffusion coefficients (Marrero & Mason 1972; Zahnle &
# Kasting 1986). Keys are unordered pairs.
# ---------------------------------------------------------------------------
MOLECULAR_BACKGROUND = {
    ('H', 'O2'): (4.75e17, 0.711),
    ('H', 'O'): (5.7e17, 0.708),
    ('H', 'He'): (8.84e17, 0.706),
    ('CO2', 'H'): (8.4e17, 0.6),
    ('H', 'N2'): (4.87e17, 0.698),
    ('H2', 'N2'): (2.8e17, 0.740),
    ('CH4', 'N2'): (7.34e16, 0.75),
    ('H2', 'O2'): (3.06e17, 0.732),
    ('CH4', 'O2'): (7.51e16, 0.759),
    ('CO2', 'H2'): (2.23e17, 0.75),
    ('CH4', 'H2'): (2.3e17, 0.765),
    ('He', 'O'): (3.44e17, 0.749),
    ('Ar', 'O'): (5.51e16, 0.841),
    ('CO', 'O2'): (8.3e16, 0.724),
    ('Ar', 'O2'): (7.17e16, 0.736),
    ('CO2', 'O2'): (5.77e16, 0.749),
    ('O', 'O2'): (9.69e16, 0.774),
    ('He', 'O2'): (3.21e17, 0.71),
    ('CO', 'N2'): (9.28e16, 0.71),
    ('Ar', 'N2'): (6.64e16, 0.752),
    ('CO2', 'N2'): (6.58e16, 0.752),
    ('N2', 'O'): (9.69e16, 0.774),
    ('He', 'N2'): (2.94e17, 0.718),
}

# ---------------------------------------------------------------------------
# The scaling rule and the library assembly.
# ---------------------------------------------------------------------------

# One fractional 1-sigma width per provenance class. 'measured': the 2023 and
# 1986 compilations agree to 0.2 to 4 percent on their shared measured rows;
# 10 percent is deliberately wider and is not a published figure.
# 'estimated': Zahnle & Kasting (1986) demonstrate a 30 percent perturbation
# for their estimated class and state most of their coefficients are
# estimates. 'scaled': the same operation the 2023 authors perform, whose
# in-sample error is of that size. 'scaled*': as 'scaled' plus an estimated
# kinetic diameter on at least one species.
SIGMA_CLASS = {'measured': 0.10, 'estimated': 0.30, 'scaled': 0.30, 'scaled*': 0.30}
CLASS_OF_ZK23 = {'M': 'measured', 'E': 'estimated'}

N_CLOSEST = 3


def eq10(b_anchor: float, target: tuple, anchor: tuple, mass: dict, diam: dict) -> float:
    """Zahnle & Kasting (2023) Eq. (10): scale b between pairs.

    ``b_ij = b_kl sqrt((1/m_i + 1/m_j) / (1/m_k + 1/m_l))
    ((d_k + d_l) / (d_i + d_j))^2``, the reduced-mass and hard-sphere
    diameter scaling their own estimated rows are built with. Masses in amu
    and diameters in pm (the units cancel).
    """
    i, j = target
    k, l = anchor  # noqa: E741
    rm_t = 1.0 / mass[i] + 1.0 / mass[j]
    rm_a = 1.0 / mass[k] + 1.0 / mass[l]
    return b_anchor * math.sqrt(rm_t / rm_a) * ((diam[k] + diam[l]) / (diam[i] + diam[j])) ** 2


def _zk23_lookup(pair: tuple):
    if pair in ZK23_TABLE2:
        return ZK23_TABLE2[pair]
    return ZK23_TABLE2.get(pair[::-1])


def _anchors_for(target: tuple, mass: dict, diam: dict):
    """Measured Table 2 rows to scale a missing pair from.

    Prefer rows sharing exactly one species with the target (the 2023
    authors' own practice); if none exists, fall back to the three rows
    closest in summed kinetic diameter, which minimizes the stretch of the
    hard-sphere factor. The central value downstream is the geometric mean
    over the retained anchors.
    """
    i, j = target
    eligible = [
        (pair, b1000)
        for pair, (b1000, cls, _src) in ZK23_TABLE2.items()
        if cls == 'M' and pair[0] != pair[1]
    ]
    out = [
        (pair, eq10(b1000, target, pair, mass, diam))
        for pair, b1000 in eligible
        if len({i, j} & set(pair)) == 1 and all(s in mass and s in diam for s in pair)
    ]
    if out:
        return out, 'shared'
    dsum_t = diam[i] + diam[j]
    ranked = sorted(
        (p for p, _b in eligible if all(s in mass and s in diam for s in p)),
        key=lambda p: abs((diam[p[0]] + diam[p[1]]) - dsum_t),
    )
    out = [(p, eq10(_zk23_lookup(p)[0], target, p, mass, diam)) for p in ranked[:N_CLOSEST]]
    return out, 'closest'


class Row:
    """One unordered species pair of the library, b(T) = b1000 (T/1000 K)^s."""

    __slots__ = ('i', 'j', 'b1000', 'exponent', 'provenance', 'uncertainty', 'source')

    def __init__(self, i, j, b1000, exponent, provenance, uncertainty, source):
        self.i, self.j = i, j
        self.b1000 = b1000
        self.exponent = exponent
        self.provenance = provenance
        self.uncertainty = uncertainty
        self.source = source

    @property
    def key(self):
        return tuple(sorted((self.i, self.j)))

    def b(self, T: float) -> float:
        """The pair's b = n D at temperature T, in cm^-1 s^-1."""
        return self.b1000 * (T / 1000.0) ** self.exponent


def build_rows(species: list) -> list[Row]:
    """Assemble one Row per unordered pair of ``species``.

    Source order: the Zahnle & Kasting (2023) Table 2 row where one is
    printed (measured or estimated, at its printed value); the Sasaki &
    Nakazawa (1988) noble-gas rows; otherwise the Eq. (10) scaling from the
    measured Table 2 anchors, with the provenance class recording that the
    row is scaled and whether it rests on an estimated diameter.
    """
    mass, diam = ALL_MASS, diameters()
    rows = []
    for a in range(len(species)):
        for c in range(a + 1, len(species)):
            i, j = species[a], species[c]
            hit = _zk23_lookup((i, j))
            if hit is not None:
                b1000, cls, src = hit
                rows.append(
                    Row(i, j, b1000, ZK23_EXPONENT, f'ZK23 T2 [{cls}]', CLASS_OF_ZK23[cls], src)
                )
                continue
            sn = None
            for a_, b_ in ((i, j), (j, i)):
                if (a_, b_) in SN88_TABLE1 and b_ in ('H2', 'He'):
                    sn = (a_, b_)
            if sn is not None:
                a_fit, s_fit = sn88_fit(*sn)
                rows.append(
                    Row(
                        i,
                        j,
                        a_fit * 1000.0**s_fit,
                        ZK23_EXPONENT,
                        'SN88 T1 [M]',
                        'measured',
                        f'Sasaki & Nakazawa (1988) Table 1; fitted exponent '
                        f'{s_fit:.3f}, entered at 0.75',
                    )
                )
                continue
            scaled, submode = _anchors_for((i, j), mass, diam)
            vals = np.array([v for _p, v in scaled])
            b1000 = float(np.exp(np.mean(np.log(vals))))
            est = tuple(s for s in (i, j) if s in _VDW_SCALED)
            rows.append(
                Row(
                    i,
                    j,
                    b1000,
                    ZK23_EXPONENT,
                    f'ZK23 Eq.(10) [{submode}]',
                    'scaled*' if est else 'scaled',
                    'anchors ' + ', '.join(f'{p[0]}-{p[1]}' for p, _ in scaled),
                )
            )
    return rows


def bmatrix(species: list, T: float, rows: list[Row] | None = None) -> np.ndarray:
    """Symmetric b matrix in cm^-1 s^-1, with np.inf on the diagonal.

    The infinite diagonal is the convention the fractionation closure
    expects (a species does not diffuse against itself).
    """
    if rows is None:
        rows = build_rows(species)
    idx = {s: k for k, s in enumerate(species)}
    n = len(species)
    b = np.full((n, n), np.inf)
    for r in rows:
        if r.i in idx and r.j in idx:
            b[idx[r.i], idx[r.j]] = b[idx[r.j], idx[r.i]] = r.b(T)
    off = ~np.eye(n, dtype=bool)
    if not np.all(np.isfinite(b[off])):
        missing = [
            (species[a], species[c])
            for a in range(n)
            for c in range(n)
            if a != c and not np.isfinite(b[a, c])
        ]
        raise ValueError(f'diffusion library incomplete for pairs {missing}')
    return b


def masses_g(species: list) -> np.ndarray:
    """Particle masses in g for the closure, from the shared mass table."""
    missing = [s for s in species if s not in ALL_MASS]
    if missing:
        raise KeyError(f'no mass tabulated for {missing}')
    return np.array([ALL_MASS[s] for s in species]) * AMU_G


# ---------------------------------------------------------------------------
# The pair ladder for arbitrary (possibly molecular) species, used by the
# hydrostatic branch, and Blanc's law for mixtures.
# ---------------------------------------------------------------------------

_PAIR_CACHE: dict = {}


def b_pair(sp_i: str, sp_j: str, T: float) -> tuple[float, str]:
    """Binary diffusion parameter b = n D for one pair, in SI [m^-1 s^-1].

    The ladder, most trusted rung first: the atomic library above (printed
    rows, then the Eq. 10 scaling); the molecular-background table; as a
    last resort, the library value of the nearest-mass covered species,
    with the substitution recorded in the provenance string. Trailing
    state annotations in species names are stripped.

    Returns ``(b [m^-1 s^-1], provenance string)``.
    """
    key = tuple(sorted((sp_i.split('_')[0], sp_j.split('_')[0])))
    if key in _PAIR_CACHE:
        b1000, s, prov = _PAIR_CACHE[key]
        return b1000 * (T / 1000.0) ** s * 100.0, prov
    a, b = key
    if a == b:
        hit = _zk23_lookup(key)
        if hit is not None:
            b1000, cls, _src = hit
            prov = f'ZK23 T2 [{cls}]'
            _PAIR_CACHE[key] = (b1000, ZK23_EXPONENT, prov)
            return b1000 * (T / 1000.0) ** ZK23_EXPONENT * 100.0, prov
    elif a in ALL_MASS and b in ALL_MASS:
        try:
            r = build_rows([a, b])[0]
            _PAIR_CACHE[key] = (r.b1000, r.exponent, r.provenance)
            return r.b(T) * 100.0, r.provenance
        except (ValueError, KeyError):
            pass
    fit = MOLECULAR_BACKGROUND.get(key) or MOLECULAR_BACKGROUND.get(key[::-1])
    if fit is not None:
        a_fit, s_fit = fit
        b1000 = a_fit * 1000.0**s_fit
        prov = 'molecular-background table'
        _PAIR_CACHE[key] = (b1000, s_fit, prov)
        return b1000 * (T / 1000.0) ** s_fit * 100.0, prov

    def _proxy(sp):
        if sp in ALL_MASS:
            return sp
        m = species_mass_amu(sp)
        return min(ALL_MASS, key=lambda s2: abs(ALL_MASS[s2] - m))

    pa, pb = _proxy(a), _proxy(b)
    if pa == pb:
        pb = 'O' if pa != 'O' else 'N'
    r = build_rows([pa, pb])[0]
    prov = f'proxy {a}->{pa}, {b}->{pb} [{r.provenance}]'
    _PAIR_CACHE[key] = (r.b1000, r.exponent, prov)
    return r.b(T) * 100.0, prov


def b_mixture(sp: str, comp: dict[str, float], T: float) -> tuple[float, dict]:
    """Blanc's-law mixture diffusion parameter for one species, SI [m^-1 s^-1].

    ``b_mix = (1 - X_sp) / sum_j (X_j / b_sp,j)`` over the other species of
    the mixture ``comp`` (mole fractions). Returns ``(b_mix, provenance
    dict per pair)``; a species alone in its mixture has nothing to diffuse
    against and returns infinity.
    """
    x_sp = comp.get(sp, 0.0)
    s = 0.0
    prov: dict[str, str] = {}
    for oj, xj in comp.items():
        if oj == sp or xj <= 0.0:
            continue
        b, p = b_pair(sp, oj, T)
        s += xj / b
        prov[f'{sp}-{oj}'] = p
    if s == 0.0:
        return math.inf, prov
    return (1.0 - x_sp) / s, prov
