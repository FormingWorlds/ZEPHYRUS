"""
!!! info "`composition.py`"
    Element masses, chemical-formula parsing, and composition handling.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import re

from zephyrus.constants import amu

# Element masses in atomic mass units. Reactive elements carry the standard
# atomic weights (IUPAC/CIAAW). Neon and argon carry the escape-relevant
# isotope masses (20Ne and 36Ar) rather than the terrestrial elemental
# averages: a primordial or solar-composition inventory is dominated by those
# isotopes, and the binary diffusion coefficients printed by Zahnle & Kasting
# (1986, Icarus 68, 462; 2023, GeCoA 361, 228) are for them. Terrestrial Ar is
# 40Ar-dominated (radiogenic), an 11 percent mass difference. Kr and Xe carry
# the standard atomic weights, which sit within 0.1 and 1.1 percent of the
# 84Kr and 130Xe isotopes those compilations tabulate.
ELEMENT_AMU = {
    'H': 1.008,
    'He': 4.0026,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'Ne': 19.992,
    'Na': 22.990,
    'Mg': 24.305,
    'Si': 28.085,
    'P': 30.974,
    'S': 32.06,
    'Cl': 35.45,
    'Ar': 35.968,
    'K': 39.098,
    'Ca': 40.078,
    'Ti': 47.867,
    'Fe': 55.845,
    'Kr': 83.798,
    'Xe': 131.293,
}

# Van der Waals radii in Angstrom, from Bondi (1964, J. Phys. Chem. 68, 441),
# Tables I and XIV as printed, except Fe: Bondi prints no transition metals,
# so the iron radius is Alvarez (2013, Dalton Trans. 42, 8617), with published
# values spanning roughly 2.0 to 2.44 Angstrom across compilations. The Mg
# value is a flagged outlier of Bondi's own table (derived from the critical
# volume and marked tentative there); Batsanov (2001) gives 2.10 to 2.27 and
# Alvarez (2013) 2.51, so quantities scaled from the Mg radius carry a 36 to
# 55 percent softness beyond their provenance class. Used as the last-resort
# geometric rung of the collision cross-section ladder and by the
# kinetic-diameter scaling rule of the binary-diffusion library.
BONDI_VDW_RADIUS_A = {
    'H': 1.20,
    'He': 1.40,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'Ne': 1.54,
    'Na': 2.27,
    'Mg': 1.73,
    'Si': 2.10,
    'S': 1.80,
    'Ar': 1.88,
    'Fe': 2.44,
    'Kr': 2.02,
    'Xe': 2.16,
}

_FORMULA_TOKEN = re.compile(r'([A-Z][a-z]?)(\d*)')


def parse_formula(name: str) -> dict[str, int]:
    """Element counts of a simple molecular formula string.

    Parses formulas of the kind atmospheric chemistry codes emit, for
    example ``'CO2' -> {'C': 1, 'O': 2}``. Trailing annotations after an
    underscore are stripped, because species names such as ``'O_1'``,
    ``'N_2D'``, or ``'S8_l_s'`` denote excited states or condensates of the
    same stoichiometry.

    Parameters
    ----------
    name : str
        Molecular formula, e.g. ``'H2O'``, ``'CO2'``, ``'He'``.

    Returns
    -------
    dict
        Element symbol to integer count.

    Raises
    ------
    ValueError
        If the string cannot be parsed as a formula or contains an element
        outside ``ELEMENT_AMU``.
    """
    base = name.split('_')[0]
    out: dict[str, int] = {}
    pos = 0
    for m in _FORMULA_TOKEN.finditer(base):
        if m.start() != pos:
            raise ValueError(f'cannot parse formula {name!r}')
        pos = m.end()
        el, cnt = m.group(1), int(m.group(2) or 1)
        if el not in ELEMENT_AMU:
            raise ValueError(f'unknown element {el!r} in formula {name!r}')
        out[el] = out.get(el, 0) + cnt
    if pos != len(base) or not out:
        raise ValueError(f'cannot parse formula {name!r}')
    return out


def species_mass_amu(name: str) -> float:
    """Molecular mass of a formula string, in atomic mass units."""
    return sum(ELEMENT_AMU[el] * n for el, n in parse_formula(name).items())


def atomize(vmr: dict[str, float]) -> dict[str, float]:
    """Element mole fractions of an atomized composition.

    Reduces molecular volume mixing ratios to the mole fractions of their
    constituent atoms, which is the composition a fully dissociated gas
    would have. Used where escape physics operates on atoms: at the XUV
    wind base, where molecules are photodissociated well below the level
    where the wind is launched (Murray-Clay et al. 2009, ApJ 693, 23).

    Parameters
    ----------
    vmr : dict
        Species name to mole fraction. Need not sum to 1; the result is
        renormalized.

    Returns
    -------
    dict
        Element symbol to mole fraction, summing to 1.

    Raises
    ------
    ValueError
        If no species carries a positive mole fraction.
    """
    counts: dict[str, float] = {}
    for sp, x in vmr.items():
        if x <= 0.0:
            continue
        for el, n in parse_formula(sp).items():
            counts[el] = counts.get(el, 0.0) + x * n
    tot = sum(counts.values())
    if tot <= 0.0:
        raise ValueError('empty composition: no species has a positive mole fraction')
    return {el: c / tot for el, c in counts.items()}


def mean_particle_mass(element_fractions: dict[str, float]) -> float:
    """Mean particle mass of an atomized composition, in kg.

    ``element_fractions`` maps element symbols to mole fractions summing
    to 1, as returned by :func:`atomize`.
    """
    return sum(x * ELEMENT_AMU[el] for el, x in element_fractions.items()) * amu
