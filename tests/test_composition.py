"""Tests for ``src/zephyrus/composition.py``.

``composition.py`` is a utility source (element masses, formula parsing, and
composition handling), so it is exempt from the physics-invariant
requirement, but a wrong mass or a mis-parsed formula silently corrupts every
downstream escape rate. These tests pin the shipped element masses, exercise
the formula parser on well-formed and malformed input, and assert the
conservation properties of the atomization step (mole fractions renormalize
to 1, element counts follow stoichiometry) with discrimination guards.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import pytest

from zephyrus.composition import (
    ELEMENT_AMU,
    atomize,
    mean_particle_mass,
    parse_formula,
    species_mass_amu,
)
from zephyrus.constants import amu

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_element_masses_pin_isotope_conventions():
    """Ne and Ar carry the escape-relevant isotope masses, not the averages.

    The table stores 20Ne (19.992 amu) and 36Ar (35.968 amu) because the
    binary-diffusion compilations of Zahnle & Kasting (1986, 2023) tabulate
    those isotopes. The discrimination guards are the terrestrial elemental
    averages (Ne 20.180, Ar 39.948 amu), which must resolvably differ; the
    reactive elements carry the standard atomic weights.
    """
    assert ELEMENT_AMU['Ne'] == pytest.approx(19.992, rel=1e-9)
    assert ELEMENT_AMU['Ar'] == pytest.approx(35.968, rel=1e-9)
    # Discrimination: the elemental-average masses differ by 1 to 11 percent.
    assert ELEMENT_AMU['Ne'] != pytest.approx(20.180, rel=1e-3)
    assert ELEMENT_AMU['Ar'] != pytest.approx(39.948, rel=1e-2)
    # Standard atomic weights for the reactive elements.
    assert ELEMENT_AMU['H'] == pytest.approx(1.008, rel=1e-9)
    assert ELEMENT_AMU['O'] == pytest.approx(15.999, rel=1e-9)


def test_parse_formula_counts_and_strips_annotations():
    """The parser returns stoichiometric counts and strips state suffixes.

    ``CO2`` splits into one carbon and two oxygens; the excited-state name
    ``N_2D`` reduces to plain atomic nitrogen; the condensate ``S8_l_s``
    keeps its eight sulfur atoms. A two-letter element must not be read as
    two one-letter elements (He is helium, not H plus e).
    """
    assert parse_formula('CO2') == {'C': 1, 'O': 2}
    assert parse_formula('N_2D') == {'N': 1}
    assert parse_formula('S8_l_s') == {'S': 8}
    assert parse_formula('He') == {'He': 1}
    # Case boundary: 'CO' is carbon monoxide, not cobalt.
    assert parse_formula('CO') == {'C': 1, 'O': 1}


def test_parse_formula_rejects_malformed_input():
    """Unparseable strings and unknown elements raise ``ValueError``.

    The error contract: lowercase-first tokens, empty strings, and symbols
    outside the mass table must raise rather than silently return a partial
    parse. A valid formula on the same call path returns normally, so the
    raise is specific to the malformed input.
    """
    with pytest.raises(ValueError, match='cannot parse'):
        parse_formula('h2o')
    with pytest.raises(ValueError, match='cannot parse'):
        parse_formula('')
    with pytest.raises(ValueError, match='unknown element'):
        parse_formula('U')
    # Contrast: the same path parses a valid formula.
    assert parse_formula('H2O') == {'H': 2, 'O': 1}


def test_species_mass_amu_matches_stoichiometric_sum():
    """Molecular masses equal the stoichiometric sums of the element masses.

    Water must weigh two hydrogens plus one oxygen (18.015 amu); the guard
    is the reversed formula weight of OH (17.007 amu), which differs by one
    hydrogen and would expose a dropped count.
    """
    m_h2o = species_mass_amu('H2O')
    assert m_h2o == pytest.approx(2 * 1.008 + 15.999, rel=1e-12)
    # Scale guard: water is 18 amu, not 17 (OH) or 19.
    assert 17.5 < m_h2o < 18.5
    assert species_mass_amu('CO2') == pytest.approx(12.011 + 2 * 15.999, rel=1e-12)


def test_atomize_conserves_stoichiometry_and_normalizes():
    """Atomization follows stoichiometry and returns fractions summing to 1.

    Pure CO2 atomizes to 1/3 carbon and 2/3 oxygen. A 50/50 H2-H2O mix has
    hydrogen and oxygen in a 4:1 mole ratio (two hydrogens from each
    molecule against one oxygen from every second molecule). Zero and
    negative entries are ignored, and the result always renormalizes to 1.
    """
    a = atomize({'CO2': 1.0})
    assert a['C'] == pytest.approx(1.0 / 3.0, rel=1e-12)
    assert a['O'] == pytest.approx(2.0 / 3.0, rel=1e-12)
    b = atomize({'H2': 0.5, 'H2O': 0.5, 'CO2': 0.0})
    assert sum(b.values()) == pytest.approx(1.0, rel=1e-12)
    # 2*0.5 + 2*0.5 = 2 hydrogens against 0.5 oxygens: a 4:1 ratio.
    assert b['H'] / b['O'] == pytest.approx(4.0, rel=1e-12)
    assert 'C' not in b  # zero-fraction species contribute nothing


def test_atomize_rejects_empty_composition():
    """An all-zero composition raises rather than dividing by zero.

    The error contract: atomize must not return NaN fractions. A one-species
    composition on the same path returns the trivial answer.
    """
    with pytest.raises(ValueError, match='empty composition'):
        atomize({'H2': 0.0, 'He': -1.0})
    assert atomize({'He': 0.2}) == {'He': 1.0}


def test_mean_particle_mass_reproduces_known_mixtures():
    """The mean particle mass interpolates linearly between the constituents.

    A 90/10 H/He atomic mix has ``0.9 * 1.008 + 0.1 * 4.0026 = 1.307 amu``,
    the textbook mean particle mass of a fully dissociated (but not ionized)
    solar-like gas. The guards bracket against pure hydrogen (1.008) and the
    molecular value (2.3 amu for H2/He), either of which would signal a
    convention slip.
    """
    m = mean_particle_mass({'H': 0.9, 'He': 0.1})
    assert m / amu == pytest.approx(0.9 * 1.008 + 0.1 * 4.0026, rel=1e-12)
    # Convention guards: atomic mean, not pure H and not the molecular mean.
    assert m / amu > 1.2
    assert m / amu < 2.0
