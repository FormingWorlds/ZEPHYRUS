"""Tests for ``src/zephyrus/hydrostatic.py``.

Exercises the Bates upper structure, the kinetic-corrected Jeans escape,
the diffusion-limited supply, and the escape-temperature gate. The physical
invariants under test:

- Reference pins: the Volkov et al. (2011) bulk-velocity correction reduces
  to unity at rest with its printed linear coefficient across lambda 1 to
  106 (verifying that the flat factor and that correction are distinct
  quantities that must not both be applied); the Yelle (2024) Figure 1 Mars
  model, whose hydrogen flux the branch reproduces through the
  Jeans-to-diffusion transition.
- Conservation: element rates sum to the species rates exactly.
- Closed forms: the escape-temperature identities and the gate semantics.
- Boundedness / error contract: the extension truncates flagged where it
  becomes unbound; sub-floor trace species are pruned, never dropped.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import numpy as np
import pytest
from scipy.special import erfcx

from zephyrus.constants import G, amu, kb
from zephyrus.hydrostatic import (
    bates_extension,
    find_exobase,
    gate_unstable,
    hydrostatic_rates,
    jeans_effusion_velocity,
    volkov_flat_factor,
)
from zephyrus.profiles import Profile

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

M_MARS = 6.4171e23  # kg (IAU nominal)
R_MARS = 3.3895e6  # m


def volkov_eq9_ratio(S, lam):
    """Volkov et al. (2011, Phys. Fluids 23, 066601) Eq. (9) flux ratio.

    The modified-over-Jeans flux ratio for a drifting Maxwellian with speed
    ratio S, in the exp(lambda)-factored form that stays stable at large
    lambda (erfcx is the scaled complementary error function).
    """
    se = math.sqrt(lam)
    t1 = 0.5 * math.exp(-S * S)
    t2 = (S * se + S * S - 0.5) * math.exp(2.0 * S * se - S * S)
    t3 = math.sqrt(math.pi) * S**3 * erfcx(se - S) * math.exp(2.0 * S * se - S * S)
    return (t1 + t2 + t3) / (S * S * (1.0 + lam))


def c_lambda(lam):
    """Closed-form linear coefficient of Eq. (9) about S = 0."""
    return (4.0 * lam**1.5 / 3.0 + 2.0 * math.sqrt(lam)) / (1.0 + lam) + math.sqrt(
        math.pi
    ) * erfcx(math.sqrt(lam)) / (1.0 + lam)


# Printed c(lambda) values of Volkov et al. (2011), Phys. Fluids table.
VOLKOV_C_TABLE = {
    1: 2.0456,
    3: 2.7254,
    6: 3.5536,
    10: 4.4355,
    15: 5.3410,
    60: 10.4126,
    106: 13.7917,
}


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_volkov_eq9_unity_limit_and_printed_coefficient():
    """Eq. (9) reduces to unity at rest with the printed linear coefficient.

    At zero bulk velocity the drifting-Maxwellian flux is the Jeans flux
    exactly; the leading correction is linear in the speed ratio with the
    printed c(lambda). Tolerance 6e-4: the published table itself carries
    float wobble of that size at large lambda when evaluated naively, and
    the stable erfcx form sits within it.
    """
    for lam, c_ref in VOLKOV_C_TABLE.items():
        assert c_lambda(lam) == pytest.approx(c_ref, rel=6e-4), lam
        r = volkov_eq9_ratio(1e-4, lam)
        assert (r - 1.0) / 1e-4 == pytest.approx(c_lambda(lam), rel=2e-3)
        assert volkov_eq9_ratio(1e-6, lam) == pytest.approx(1.0, abs=1e-4)


def test_volkov_flat_factor_shape_and_distinctness():
    """The flat kinetic factor has its measured shape, distinct from Eq. (9).

    C(lambda) runs from 1.7 at lambda = 6 to 1.4 at lambda = 15, is held at
    the endpoints outside (the flagged extrapolation), and falls with
    lambda, whereas the Eq. (9) bulk-velocity correction at fixed speed
    ratio rises with lambda: opposite slopes, so the two corrections are
    different quantities and applying both would double-count.
    """
    assert volkov_flat_factor(6.0) == pytest.approx(1.7, rel=1e-12)
    assert volkov_flat_factor(15.0) == pytest.approx(1.4, rel=1e-12)
    assert volkov_flat_factor(50.0) == pytest.approx(1.4, rel=1e-12)  # held
    assert volkov_flat_factor(10.5) == pytest.approx(1.55, rel=0.01)
    assert volkov_flat_factor(15.0) < volkov_flat_factor(6.0)
    assert volkov_eq9_ratio(0.1, 15.0) > volkov_eq9_ratio(0.1, 6.0)


def _mars_profile():
    """Minimal profile ending at the Yelle (2024) Mars anchor level.

    Their fully specified model: p = 0.1 Pa at 80 km altitude, T = 100 K,
    CO2 background carrying 10 ppm total hydrogen.
    """
    r0 = R_MARS + 8.0e4
    p = np.array([10.0, 1.0, 0.1])
    T = np.full(3, 100.0)
    mu = 44.0095 * amu
    r = np.empty(3)
    r[2] = r0
    for i in (1, 0):
        H = kb * 100.0 * r[i + 1] ** 2 / (G * M_MARS * mu)
        r[i] = r[i + 1] - H * math.log(p[i] / p[i + 1])
    vmr = {'CO2': np.full(3, 1.0 - 1e-5), 'H': np.full(3, 1e-5)}
    return Profile(p=p, r=r, T=T, vmr=vmr, mmw=np.full(3, mu), kzz=None)


def _mars_h_flux(t_inf):
    """Hydrogen number flux per anchor area, cm^-2 s^-1, at exobase T."""
    prof = _mars_profile()
    per_el, det = hydrostatic_rates(prof, M_MARS, t_inf, gamma_bates=0.75, kzz_default=3.0e2)
    rate = det['per_species_rate'].get('H', 0.0)
    m_h = 1.008 * amu
    return rate / m_h / (4.0 * math.pi * det['r_anchor'] ** 2) * 1e-4, det


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_yelle_figure1_mars_hydrogen_flux():
    """The branch reproduces the Yelle (2024) Figure 1 Mars hydrogen flux.

    On their fully specified Mars model the hydrogen escape flux is
    diffusion limited above about 200 K exobase temperature, with a plateau
    at 2.4e8 cm^-2 s^-1 (40 percent tolerance: the binary H-CO2 coefficient
    source their calculation used is not pinned in the paper, and the
    tabulated value here differs at that level). The transition to
    Jeans-limited escape below about 150 K shows as a collapse; the 100 K
    point sits on the exponential edge, so it is checked as a regime (an
    order below the plateau), not as a value.
    """
    f100, _ = _mars_h_flux(100.0)
    f200, _ = _mars_h_flux(200.0)
    f300, _ = _mars_h_flux(300.0)
    f400, _ = _mars_h_flux(400.0)
    assert f300 == pytest.approx(2.4e8, rel=0.4)
    assert f400 / f300 == pytest.approx(1.0, abs=0.2)  # saturated plateau
    assert f200 / f400 > 0.5
    assert f100 / f400 < 0.07  # Jeans-limited collapse


@pytest.mark.physics_invariant
def test_escape_temperature_identities_and_gate():
    """Escape temperatures obey their defining identities; the gate splits.

    The neutral escape temperature is ``G M m / (2 kB r)`` at the exobase
    and the plasma value is half of it exactly. Gate semantics: unstable
    when the exobase temperature exceeds half the gating escape
    temperature; points where only the plasma convention is exceeded are
    contested (the ion physics decides them), reported through the second
    return.
    """
    prof = _mars_profile()
    _per, det = hydrostatic_rates(prof, M_MARS, 300.0)
    m = det['m_bar']
    r = det['r_exo']
    assert det['T_esc_neutral'] == pytest.approx(G * M_MARS * m / (2 * kb * r), rel=1e-12)
    assert det['T_esc_plasma'] == pytest.approx(det['T_esc_neutral'] / 2.0, rel=1e-12)
    det2 = dict(det, T_esc_neutral=500.0, T_esc_plasma=250.0)
    # Gate thresholds: neutral at 250 K, plasma at 125 K exobase temperature.
    assert gate_unstable(260.0, det2, 'neutral', 0.0) == (True, False)
    assert gate_unstable(120.0, det2, 'neutral', 0.0) == (False, False)
    assert gate_unstable(200.0, det2, 'neutral', 0.0) == (False, True)
    assert gate_unstable(200.0, det2, 'plasma', 0.0) == (True, True)


@pytest.mark.physics_invariant
def test_element_mapping_conserves_mass_and_dominant_bypass():
    """Element rates sum to species rates; the dominant species self-supplies.

    Stoichiometric mapping conserves the total mass rate exactly. The
    dominant species (CO2 here) has no diffusion-limited supply cap: it
    supplies itself, so its diffusion flux is infinite and the bypass is
    recorded.
    """
    prof = _mars_profile()
    per_el, det = hydrostatic_rates(prof, M_MARS, 300.0)
    assert sum(per_el.values()) == pytest.approx(
        sum(det['per_species_rate'].values()), rel=1e-9
    )
    assert det['flags']['dl_bypass'] == 'CO2'
    assert math.isinf(det['species']['CO2']['phi_diffusion'])
    assert det['species']['H']['dl_bypass'] is False
    # The CO2 mass rate splits onto C and O in stoichiometric proportion.
    rate_co2 = det['per_species_rate']['CO2']
    assert per_el['C'] + per_el['O'] + per_el['H'] == pytest.approx(
        rate_co2 + det['per_species_rate']['H'], rel=1e-9
    )


def test_extension_truncates_when_unbound_and_flags():
    """The Bates extension stops, flagged, where the structure unbinds.

    A very hot exobase temperature on a small planet drives the local Jeans
    parameter below 2 within the extension; the arrays must stop there with
    ``unbound`` set rather than integrating an unbound structure. A cool
    case runs the full requested span. Radii grow strictly monotonically in
    both cases.
    """
    prof = _mars_profile()
    hot = bates_extension(prof, M_MARS, 4000.0)
    assert hot['unbound'] is True
    assert len(hot['zeta']) < 400
    cool = bates_extension(prof, M_MARS, 200.0)
    assert cool['unbound'] is False
    assert len(cool['zeta']) == 400
    for ext in (hot, cool):
        assert np.all(np.diff(ext['r']) > 0)
    # The Bates profile approaches its asymptotic temperature from below.
    assert cool['T'][-1] == pytest.approx(200.0, rel=1e-6)
    assert cool['T'][0] == pytest.approx(100.0, rel=1e-9)


def test_exobase_locator_contract():
    """The exobase locator returns an interior index with sane flags.

    On the Mars model the exobase lies inside the extension (no flags); on
    an artificially truncated extension the top level is used with
    ``exobase_not_reached``; an exobase at the anchor keeps one interval
    with ``exobase_at_anchor``.
    """
    prof = _mars_profile()
    ext = bates_extension(prof, M_MARS, 300.0)
    i_x, flags = find_exobase(ext, M_MARS)
    assert 0 < i_x < len(ext['zeta'])
    assert flags == {}
    truncated = {k: (v[:5] if isinstance(v, np.ndarray) else v) for k, v in ext.items()}
    i_t, flags_t = find_exobase(truncated, M_MARS)
    assert flags_t.get('exobase_not_reached') is True
    assert i_t == 4
    # A very rarefied anchor puts the exobase at the anchor itself.
    thin = dict(ext)
    thin['n'] = ext['n'] * 1e-12
    i_a, flags_a = find_exobase(thin, M_MARS)
    assert flags_a.get('exobase_at_anchor') is True
    assert i_a == 1


def test_trace_species_below_floor_are_pruned_not_dropped():
    """Sub-floor trace species skip the supply integrals but keep a rate.

    A species whose supply-free Jeans rate is already below one proton mass
    per year is numerically negligible; the branch must skip its expensive
    supply integrals (``pruned`` set) while still returning its (tiny,
    non-negative) Jeans rate, so totality is preserved. A heavy trace
    species on the cold Mars model is such a case.
    """
    prof = _mars_profile()
    vmr = {
        'CO2': prof.vmr['CO2'] - 1e-7,
        'H': prof.vmr['H'],
        'Kr': np.full(3, 1e-7),  # heavy trace: Jeans rate far below the floor
    }
    prof2 = Profile(p=prof.p, r=prof.r, T=prof.T, vmr=vmr, mmw=prof.mmw, kzz=None)
    per_el, det = hydrostatic_rates(prof2, M_MARS, 250.0)
    assert det['species']['Kr']['pruned'] is True
    assert det['species']['H']['pruned'] is False
    assert det['per_species_rate']['Kr'] >= 0.0
    assert det['per_species_rate']['Kr'] < det['per_species_rate']['H']
    # Conservation still holds with the pruned species included.
    assert sum(per_el.values()) == pytest.approx(
        sum(det['per_species_rate'].values()), rel=1e-9
    )


@pytest.mark.physics_invariant
def test_jeans_effusion_velocity_shape():
    """The effusion velocity carries the (1 + lambda) exp(-lambda) shape.

    Exact two-point ratio in lambda at fixed temperature and mass, and the
    thermal-speed prefactor scales as sqrt(T/m): heavier species at the
    same temperature effuse more slowly on both counts.
    """
    T, m = 300.0, 1.008 * amu
    v5 = jeans_effusion_velocity(T, m, 5.0)
    v10 = jeans_effusion_velocity(T, m, 10.0)
    assert v10 / v5 == pytest.approx((11.0 / 6.0) * math.exp(-5.0), rel=1e-12)
    v_heavy = jeans_effusion_velocity(T, 16 * m, 5.0)
    assert v5 / v_heavy == pytest.approx(4.0, rel=0.01)  # sqrt(16) prefactor
    assert v_heavy > 0.0
