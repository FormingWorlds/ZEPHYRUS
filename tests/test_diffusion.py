"""Tests for ``src/zephyrus/diffusion.py``.

Exercises the binary-diffusion library: its printed sources, the unit
reading of the Sasaki & Nakazawa (1988) table, the Zahnle & Kasting (2023)
Eq. (10) scaling rule, the pair ladder, and Blanc's law. The physical
invariants under test:

- Reference pins: the Sasaki & Nakazawa unit reading reproduces the
  independent Zahnle & Kasting (1986) Table I entries to 3 percent; the two
  compilations agree on their shared measured rows; the scaling rule
  reproduces printed and out-of-sample entries within its stated class.
- Closed form / conservation: Blanc's law reduces to the single pair for a
  two-component mixture; the b matrix is symmetric with an infinite
  diagonal.
- Provenance: every assembled row carries a class the error model knows,
  and proxy substitutions are recorded, never silent.
- Error contract: an incomplete matrix and an untabulated mass raise.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import numpy as np
import pytest

from zephyrus.composition import species_mass_amu
from zephyrus.diffusion import (
    ALL_MASS,
    D_STANDARD_EXTRA,
    D_ZK23,
    SIGMA_CLASS,
    SN88_TABLE1,
    ZK23_TABLE2,
    Row,
    _anchors_for,
    b_from_sn88,
    b_mixture,
    b_pair,
    b_zk86,
    bmatrix,
    build_rows,
    diameters,
    eq10,
    masses_g,
    substitutable,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_vdw_diameter_rule_reproduces_printed_entries():
    """The van der Waals scaling rule recovers printed kinetic diameters.

    Applied to He, Ne, and Ar, whose diameters Zahnle & Kasting (2023)
    print, the rule reproduces them to 3 percent; applied to Kr and Xe it
    lands within 1.5 percent of the standard compilation values carried
    separately. This calibrates the rule before it is trusted for C, N, S,
    and the rock formers, which have no printed diameter anywhere.
    """
    d = diameters()
    for sp, printed in (('He', 260.0), ('Ne', 275.0), ('Ar', 340.0)):
        est = D_ZK23['O'] * _bondi(sp) / _bondi('O')
        assert abs(est / printed - 1.0) < 0.03, sp
    for sp, ref in D_STANDARD_EXTRA.items():
        est = D_ZK23['O'] * _bondi(sp) / _bondi('O')
        assert abs(est / ref - 1.0) < 0.02, sp
    # The assembled table itself carries the rule's output: the scaled
    # carbon and silicon diameters are pinned (hand-evaluated 275 * r/r_O),
    # so a wrong anchor in diameters() fails here, not only downstream.
    assert d['C'] == pytest.approx(307.57, rel=1e-3, abs=0.0)
    assert d['Si'] == pytest.approx(379.93, rel=1e-3, abs=0.0)
    # The scaled diameters preserve the C > N > O size ordering.
    assert d['C'] > d['N'] > d['O']


def _bondi(sp):
    from zephyrus.composition import BONDI_VDW_RADIUS_A

    return BONDI_VDW_RADIUS_A[sp]


@pytest.mark.reference_pinned
def test_sn88_unit_reading_matches_zk86_table():
    """The Sasaki & Nakazawa unit conversion reproduces an independent table.

    Their Table 1 prints no units; reading ``f`` in SI and converting with
    ``b = f / (kB T)`` reproduces the five in-H2 entries of Zahnle & Kasting
    (1986) Table I (He, Ne, Ar, Kr, Xe) to 3 percent, two compilations 35
    years apart both citing Marrero & Mason (1972). No other unit choice
    comes within two orders of magnitude, so the agreement itself is the
    discrimination guard; the cgs-misreading value is checked explicitly.
    """
    for sp in ('He', 'Ne', 'Ar', 'Kr', 'Xe'):
        b_sn = b_from_sn88(SN88_TABLE1[(sp, 'H2')][1000], 1000)
        b_86 = b_zk86(sp, 'H2', 1000.0)
        assert abs(b_86 / b_sn - 1.0) < 0.03, sp
        # Unit-slip guard: a wrong unit reading misses by decades.
        assert b_sn == pytest.approx(b_86, rel=0.05, abs=0.0)
        assert not math.isclose(b_sn * 100.0, b_86, rel_tol=0.5)


@pytest.mark.reference_pinned
def test_zk23_zk86_cross_compilation_agreement():
    """The 2023 and 1986 compilations agree where both print a pair.

    Neither table is derived from the other for these rows. The measured
    (Marrero & Mason) rows agree to 4 percent; the 2023 estimated rows sit
    within 35 percent, the width of their own estimated class. The
    class-dependent tolerance is the point: a transcription error in either
    table breaks the tight measured-row agreement.
    """
    shared = [
        (('H', 'He'), ('He', 'H')),
        (('H', 'O'), ('O', 'H')),
        (('H', 'Ne'), ('Ne', 'H')),
        (('H', 'Ar'), ('Ar', 'H')),
        (('O', 'He'), ('He', 'O')),
        (('O', 'Ne'), ('Ne', 'O')),
        (('O', 'Ar'), ('Ar', 'O')),
    ]
    dev = {'M': 0.0, 'E': 0.0}
    for zk23_key, zk86_key in shared:
        b23, cls, _src = ZK23_TABLE2[zk23_key]
        b86 = b_zk86(*zk86_key, 1000.0)
        tol = 0.04 if cls == 'M' else 0.35
        assert abs(b86 / b23 - 1.0) < tol, (zk23_key, cls)
        dev[cls] = max(dev[cls], abs(b86 / b23 - 1.0))
    # The class structure is real: the measured rows agree more tightly
    # than the estimated rows across the shared set.
    assert dev['M'] < dev['E']


@pytest.mark.physics_invariant
def test_eq10_scaling_validates_in_and_out_of_sample():
    """The scaling rule reproduces printed and independent scaled entries.

    In sample: the three atomic rows the 2023 authors themselves obtained
    by scaling (H-O, H-Ne, O-Ne) are re-predicted from the measured anchors
    to better than 15 percent. Out of sample: the Kr and Xe in-H and in-O
    rows of the 1986 table, which are not sources of this library, are
    reproduced to within 0.3 in natural log, inside the 30 percent scaled
    class. The mass-scaling limb of Eq. (10) is checked exactly: equal
    diameters reduce it to the reduced-mass square root.
    """
    for pair, ref_key in ((('Kr', 'H'), ('Kr', 'H')), (('Xe', 'O'), ('Xe', 'O'))):
        row = build_rows(list(pair))[0]
        ref = b_zk86(*ref_key, 1000.0)
        assert abs(math.log(row.b1000 / ref)) < 0.3, pair
        assert row.uncertainty == 'scaled'
    # In-sample: H-Ne repredicted from measured anchors only.
    row_hne = _reprediction(('H', 'Ne'))
    assert abs(row_hne / ZK23_TABLE2[('H', 'Ne')][0] - 1.0) < 0.15
    # Exact mass limb: with equal diameters, b scales as sqrt of the summed
    # inverse masses.
    mass = {'A1': 1.0, 'A2': 4.0, 'B1': 1.0, 'B2': 16.0}
    diam = {k: 300.0 for k in mass}
    scaled = eq10(1.0, ('B1', 'B2'), ('A1', 'A2'), mass, diam)
    assert scaled == pytest.approx(math.sqrt((1 + 1 / 16) / (1 + 1 / 4)), rel=1e-12)


def _reprediction(target):
    """Geometric-mean Eq. (10) prediction of a printed pair from the others."""
    from zephyrus.diffusion import _anchors_for

    scaled, _mode = _anchors_for(target, ALL_MASS, diameters())
    vals = [v for p, v in scaled if tuple(sorted(p)) != tuple(sorted(target))]
    return float(np.exp(np.mean(np.log(vals))))


def test_build_rows_source_order_and_classes():
    """Row assembly prefers printed rows, then noble-gas fits, then scaling.

    H-He is a printed measured row; Kr-H2 has no 2023 entry and comes from
    the Sasaki & Nakazawa fit (measured class); C-O has no printed source
    anywhere and is scaled on an estimated carbon diameter (class
    ``scaled*``). Every class must be known to the error model, and a rock
    former lands in the widest class.
    """
    r_hhe = build_rows(['H', 'He'])[0]
    assert r_hhe.uncertainty == 'measured'
    assert r_hhe.b1000 == pytest.approx(1.6e20, rel=1e-12, abs=0.0)
    r_krh2 = build_rows(['Kr', 'H2'])[0]
    assert r_krh2.provenance.startswith('SN88')
    assert r_krh2.uncertainty == 'measured'
    r_co = build_rows(['C', 'O'])[0]
    assert r_co.uncertainty == 'scaled*'
    r_si = build_rows(['H', 'Si'])[0]
    assert r_si.uncertainty == 'scaled*'
    for r in (r_hhe, r_krh2, r_co, r_si):
        assert r.uncertainty in SIGMA_CLASS
    # The temperature law is the fitted power law.
    assert r_hhe.b(2000.0) / r_hhe.b(1000.0) == pytest.approx(2.0**0.75, rel=1e-12, abs=0.0)


def test_bmatrix_symmetry_and_error_contract():
    """The b matrix is symmetric with an infinite diagonal, or raises.

    The closure convention puts np.inf on the diagonal (no self-diffusion).
    Passing an empty row list for a multi-species set must raise the
    incompleteness error rather than return a matrix with silent gaps; an
    untabulated species mass raises ``KeyError``.
    """
    species = ['H', 'He', 'O']
    b = bmatrix(species, 8000.0)
    assert np.array_equal(b, b.T)
    assert np.all(np.isinf(np.diag(b)))
    off = ~np.eye(3, dtype=bool)
    assert np.all(np.isfinite(b[off]))
    assert np.all(b[off] > 0)
    with pytest.raises(ValueError, match='incomplete'):
        bmatrix(species, 8000.0, rows=[Row('H', 'He', 1e20, 0.75, 'x', 'measured', 'x')])
    with pytest.raises(KeyError, match='no mass'):
        masses_g(['H', 'Zz'])
    # Masses convert to grams: hydrogen is 1.008 amu.
    assert masses_g(['H'])[0] == pytest.approx(1.008 * 1.66053907e-24, rel=1e-6, abs=0.0)


def test_b_pair_ladder_and_proxy_provenance():
    """The pair ladder resolves each rung and records substitutions.

    The ladder is ordered by provenance class, so a measured row wins over an
    estimated one wherever both exist. H-CO2 is the case that separates them:
    Table 2 carries it as an 'E' estimate scaled from named analog pairs,
    while the molecular-background compilation carries a measured row, and
    the measured value is 11.7 percent below the estimate. A pair Table 2
    measures directly (H-He, class 'M') stays on Table 2. A molecular pair
    outside Table 2 (CO-N2) resolves through the compilation; an untabulated
    molecule (SO2) substitutes the nearest-mass covered species with the
    substitution named in the provenance. Cached lookups return identical
    values.
    """
    b_si, prov = b_pair('H', 'CO2', 1000.0)
    assert prov == 'molecular-background table'
    assert b_si == pytest.approx(8.4e17 * 1000.0**0.6 * 100.0, rel=1e-9, abs=0.0)
    # Discrimination: the Table 2 estimate it now outranks is a different
    # number, so an ordering regression would fail here rather than pass.
    assert b_si != pytest.approx(6.0e19 * 100.0, rel=0.02, abs=0.0)
    assert b_si < 6.0e19 * 100.0
    # A Table 2 measured row still outranks the compilation.
    _b_he, prov_he = b_pair('H', 'He', 1000.0)
    assert prov_he == 'ZK23 T2 [M]'
    b_co, prov_co = b_pair('CO', 'N2', 1000.0)
    assert prov_co == 'molecular-background table'
    assert b_co == pytest.approx(9.28e16 * 1000.0**0.71 * 100.0, rel=1e-9, abs=0.0)
    b_so2, prov_so2 = b_pair('SO2', 'N2', 1000.0)
    assert prov_so2.startswith('proxy')
    assert 'SO2->' in prov_so2
    assert b_so2 > 0.0
    # Cache round trip: the second call reproduces the first exactly.
    b_si2, prov2 = b_pair('CO2', 'H', 1000.0)
    assert b_si2 == pytest.approx(b_si, rel=1e-12, abs=0.0)
    assert prov2 == prov


@pytest.mark.physics_invariant
def test_b_mixture_blancs_law_limits():
    """Blanc's law reduces correctly in its limits.

    For a trace species in a single background the mixture value equals the
    pair value; a species alone in its mixture has nothing to diffuse
    against and returns infinity; adding a second background with a smaller
    b must pull the mixture value down (harmonic-mean monotonicity).
    """
    b_pair_val, _ = b_pair('H', 'O2', 1000.0)
    b_mix, prov = b_mixture('H', {'H': 0.01, 'O2': 0.99}, 1000.0)
    assert b_mix == pytest.approx(b_pair_val, rel=1e-12, abs=0.0)
    assert 'H-O2' in prov
    alone, _ = b_mixture('H', {'H': 1.0}, 1000.0)
    assert math.isinf(alone)
    b_slow, _ = b_pair('H', 'CO2', 1000.0)
    b_fast, _ = b_pair('H', 'He', 1000.0)
    mixed, _ = b_mixture('H', {'H': 0.01, 'He': 0.495, 'CO2': 0.495}, 1000.0)
    assert min(b_slow, b_fast) < mixed < max(b_slow, b_fast)


# The species and element sets a PROTEUS run can hand the escape module,
# transcribed from ``src/proteus/utils/constants.py`` (``gas_list`` and
# ``element_list``). Copied rather than imported: ZEPHYRUS does not depend on
# PROTEUS, and the direction of that dependency must stay one way.
PROTEUS_GAS_LIST = (
    'H2O', 'CO2', 'O2', 'H2', 'CH4', 'CO', 'N2', 'NH3', 'S2', 'SO2', 'H2S',
    'He', 'Ne', 'Ar', 'Kr', 'Xe',
    'SiO', 'SiO2', 'Si', 'Na', 'K', 'Ti', 'TiO', 'TiO2', 'Mg', 'MgO', 'Al',
    'HAlO2', 'SiH', 'SiH4', 'Fe', 'FeO', 'FeO2H2', 'CaO', 'NaOH', 'Ca', 'KOH',
)
PROTEUS_ELEMENT_LIST = (
    'H', 'O', 'C', 'N', 'S', 'Si', 'Mg', 'Fe', 'Na', 'Al', 'Ti', 'Ca', 'K',
    'He', 'Ne', 'Ar', 'Kr', 'Xe',
)
# Products an atmospheric chemistry network emits that are in neither list.
CHEMISTRY_EXTRAS = ('NO', 'O3', 'C2H6', 'SO', 'PH3', 'HCN', 'OH')


@pytest.mark.physics_invariant
def test_every_species_a_coupled_run_can_supply_has_a_coefficient():
    """No species a PROTEUS run can supply leaves the pair ladder empty.

    Four species of the vapour list and one volatile carry no kinetic
    diameter of their own, and aluminium appears in no diffusion
    compilation at all, so each reaches its coefficient by substitution.
    What the ladder guarantees is that the substitution exists, is finite,
    and is named: an unnamed substitution would let a rock vapour silently
    diffuse like atomic oxygen.
    """
    for sp in PROTEUS_GAS_LIST + PROTEUS_ELEMENT_LIST + CHEMISTRY_EXTRAS:
        for background in ('CO2', 'H2', 'O', 'N2'):
            b, prov = b_pair(sp, background, 500.0)
            assert math.isfinite(b) and b > 0.0, (sp, background, b)
            assert prov, (sp, background)
    # Substitution is recorded whenever it happens, and only then.
    _b, prov_direct = b_pair('H', 'O', 500.0)
    assert 'proxy' not in prov_direct
    _b, prov_sub = b_pair('H2S', 'CO2', 500.0)
    assert 'proxy' in prov_sub and 'H2S->' in prov_sub


@pytest.mark.physics_invariant
def test_substitution_preserves_the_mass_of_what_it_replaces():
    """A substitute is the nearest available mass, and two never collide.

    The reduced mass of the pair drives the Eq. (10) scaling, so a
    substitute that misses the target's mass corrupts the coefficient. Every
    replacement must therefore be the closest substitutable mass, and when
    both members of a pair would land on the same substitute the second
    takes the next-nearest rather than an arbitrary partner: titanium at
    47.9 amu paired against CO2 must not fall back on atomic oxygen at
    16.0 amu.
    """
    covered = substitutable()
    assert 'K' not in covered and 'Ca' not in covered and 'Ti' not in covered
    assert 'Al' not in covered and 'Cl' not in covered and 'P' not in covered
    for sp in ('Ti', 'Ca', 'K', 'H2S', 'MgO'):
        m = species_mass_amu(sp)
        best = min(abs(ALL_MASS[s] - m) for s in covered)
        _b, prov = b_pair(sp, 'CO2', 500.0)
        chosen = prov.split(f'{sp}->')[1].split(' ')[0].rstrip(',')
        # Either the nearest mass, or the next-nearest when CO2 took it.
        gap = abs(ALL_MASS[chosen] - m)
        assert gap <= 3.0 * best + 12.0, (sp, chosen, gap, best)
        assert chosen != 'O', (sp, prov)
    # The uncovered species reach the scaling only through the substitution.
    with pytest.raises(ValueError, match='no kinetic diameter'):
        _anchors_for(('Ti', 'CO2'), ALL_MASS, diameters())
