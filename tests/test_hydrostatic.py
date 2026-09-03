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

from zephyrus.composition import species_mass_amu
from zephyrus.constants import G, amu, kb
from zephyrus.hydrostatic import (
    bates_extension,
    find_exobase,
    gate_unstable,
    hydrostatic_rates,
    hydrostatic_rates_refined,
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
def test_volkov_flat_factor_against_published_correction():
    """The source's flat factor is distinct from the published Eq. (9).

    Two anchors in one test, both against Volkov et al. (2011). First, the
    test-local Eq. (9) oracle is certified against the paper: it reduces to
    unity at zero bulk velocity, and its leading correction is linear in
    the speed ratio with the printed c(lambda) table (tolerance 6e-4, the
    float wobble the published table itself carries at large lambda).
    Second, the source's ``volkov_flat_factor`` carries the measured shape
    (1.7 at lambda = 6 to 1.4 at lambda = 15, held at the endpoints as the
    flagged extrapolation) and falls with lambda while the certified
    Eq. (9) correction rises with it: opposite slopes, so the two
    corrections are different quantities and applying both to the branch
    would double-count.
    """
    for lam, c_ref in VOLKOV_C_TABLE.items():
        assert c_lambda(lam) == pytest.approx(c_ref, rel=6e-4, abs=0.0), lam
        r = volkov_eq9_ratio(1e-4, lam)
        assert (r - 1.0) / 1e-4 == pytest.approx(c_lambda(lam), rel=2e-3, abs=0.0)
        assert volkov_eq9_ratio(1e-6, lam) == pytest.approx(1.0, abs=1e-4)
    assert volkov_flat_factor(6.0) == pytest.approx(1.7, rel=1e-12, abs=0.0)
    assert volkov_flat_factor(15.0) == pytest.approx(1.4, rel=1e-12, abs=0.0)
    assert volkov_flat_factor(50.0) == pytest.approx(1.4, rel=1e-12, abs=0.0)  # held
    assert volkov_flat_factor(10.5) == pytest.approx(1.55, rel=0.01, abs=0.0)
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


def _co2_hydrogen_profile(x_h2):
    """The Mars anchor profile with the hydrogen abundance set explicitly."""
    r0 = R_MARS + 8.0e4
    p = np.array([10.0, 1.0, 0.1])
    T = np.full(3, 100.0)
    mu = 44.0095 * amu
    r = np.empty(3)
    r[2] = r0
    for i in (1, 0):
        H = kb * 100.0 * r[i + 1] ** 2 / (G * M_MARS * mu)
        r[i] = r[i + 1] - H * math.log(p[i] / p[i + 1])
    vmr = {'CO2': np.full(3, 1.0 - x_h2)}
    if x_h2 > 0.0:
        vmr['H'] = np.full(3, x_h2)
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
    at 2.4e8 cm^-2 s^-1. The 40 percent tolerance is deliberate and its
    direction understood: the binary H-CO2 coefficient source their
    calculation used is not pinned in the paper, the tabulated coefficient
    here is smaller, and this implementation sits systematically below the
    anchor (measured plateau 1.70e8, a ratio of 0.71 against the 0.60
    floor), so a future coefficient revision in either direction moves
    this pin and should be re-tuned consciously rather than by widening
    the tolerance. The transition to Jeans-limited escape below about
    150 K shows as a collapse; the 100 K point sits on the exponential
    edge, so it is checked as a regime (an order below the plateau), not
    as a value.
    """
    f100, _ = _mars_h_flux(100.0)
    f200, _ = _mars_h_flux(200.0)
    f300, _ = _mars_h_flux(300.0)
    f400, _ = _mars_h_flux(400.0)
    assert f300 == pytest.approx(2.4e8, rel=0.4, abs=0.0)
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
    assert det['T_esc_neutral'] == pytest.approx(G * M_MARS * m / (2 * kb * r), rel=1e-12, abs=0.0)
    assert det['T_esc_plasma'] == pytest.approx(det['T_esc_neutral'] / 2.0, rel=1e-12, abs=0.0)
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
    recorded against that species and no other.
    """
    prof = _mars_profile()
    per_el, det = hydrostatic_rates(prof, M_MARS, 300.0)
    assert sum(per_el.values()) == pytest.approx(
        sum(det['per_species_rate'].values()), rel=1e-9
    )
    # The dominant species is the most abundant one at the exobase, and it
    # is the only one whose supply cap is bypassed.
    assert det['dominant'] == 'CO2'
    assert det['species']['CO2']['dl_bypass'] is True
    assert math.isinf(det['species']['CO2']['phi_diffusion'])
    assert det['species']['H']['dl_bypass'] is False
    assert math.isfinite(det['species']['H']['phi_diffusion'])
    # The CO2 mass rate splits onto C and O in stoichiometric proportion.
    rate_co2 = det['per_species_rate']['CO2']
    assert per_el['C'] + per_el['O'] + per_el['H'] == pytest.approx(
        rate_co2 + det['per_species_rate']['H'], rel=1e-9
    , abs=0.0)


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
    assert cool['T'][-1] == pytest.approx(200.0, rel=1e-6, abs=0.0)
    assert cool['T'][0] == pytest.approx(100.0, rel=1e-9, abs=0.0)


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
    assert v10 / v5 == pytest.approx((11.0 / 6.0) * math.exp(-5.0), rel=1e-12, abs=0.0)
    v_heavy = jeans_effusion_velocity(T, 16 * m, 5.0)
    assert v5 / v_heavy == pytest.approx(4.0, rel=0.01, abs=0.0)  # sqrt(16) prefactor
    assert v_heavy > 0.0


@pytest.mark.physics_invariant
def test_trace_species_survive_into_the_exobase_anchor():
    """A trace light species keeps its rate however thin it is.

    The escaping flux of a minor species is linear in its mixing ratio, so
    on a heavy background whose own rate is twenty decades lower the bulk
    rate must follow the trace hydrogen down without a step. A lower cut on
    the anchor mixing ratio would instead delete the species that carries
    the whole rate, and the rate would fall to the background's the moment
    the abundance crossed it.
    """
    rates, keys = [], []
    fractions = (1e-6, 1e-8, 1e-9, 1e-12, 1e-15)
    for x_h2 in fractions:
        prof = _co2_hydrogen_profile(x_h2)
        per_el, _det = hydrostatic_rates(prof, M_MARS, 1000.0)
        rates.append(sum(per_el.values()))
        keys.append(set(per_el))
    # Linear in abundance across nine decades, so no threshold sits inside.
    for x, rate in zip(fractions, rates):
        assert rate / x == pytest.approx(rates[0] / fractions[0], rel=1e-3, abs=0.0), (x, rate)
    # Hydrogen is reported at every abundance, never dropped from the split.
    assert all('H' in k for k in keys)
    # Removing it entirely is the only way to lose it, and then the rate
    # collapses to the heavy background, which is what the step looked like.
    bare = _co2_hydrogen_profile(0.0)
    per_el_bare, _ = hydrostatic_rates(bare, M_MARS, 1000.0)
    assert 'H' not in per_el_bare
    assert sum(per_el_bare.values()) < 1e-4 * rates[-1]


@pytest.mark.physics_invariant
def test_supply_quadrature_refines_to_its_target():
    """The supply integrals are refined until they stop moving the rate.

    The integrals are first-order accurate in the log-pressure step, so the
    change between a grid and its refinement estimates what is left to
    converge, and one fixed count cannot report its own error. Doubling from
    a coarse start must therefore drive the change below the target and say
    so, the refined answer must be the finest grid's rather than an
    extrapolation, and a ceiling reached without meeting the target must be
    reported rather than passed off as converged.
    """
    prof = _co2_hydrogen_profile(0.01)
    per, det = hydrostatic_rates_refined(prof, M_MARS, 1000.0, rtol=1e-2)
    conv = det['convergence']
    assert conv['converged'] is True
    assert conv['rel_change_bulk'] <= 1e-2
    assert conv['n_levels'] >= conv['n_levels_min']
    # The returned rates are the ones the finest grid produced.
    per_at_n, _ = hydrostatic_rates(prof, M_MARS, 1000.0, n_levels=conv['n_levels'])
    assert sum(per.values()) == pytest.approx(sum(per_at_n.values()), rel=1e-12, abs=0.0)
    # First order in the step: the coarse grid is measurably off, and the
    # refinement moves the answer toward the fine one monotonically.
    coarse, _ = hydrostatic_rates(prof, M_MARS, 1000.0, n_levels=100)
    fine, _ = hydrostatic_rates(prof, M_MARS, 1000.0, n_levels=3200)
    assert abs(sum(coarse.values()) - sum(fine.values())) / sum(fine.values()) > 1e-2
    assert sum(coarse.values()) > sum(per.values()) >= sum(fine.values())
    # A ceiling below the target is reported, not hidden.
    _p, det_capped = hydrostatic_rates_refined(
        prof, M_MARS, 1000.0, n_levels_min=50, n_levels_max=100, rtol=1e-12
    )
    assert det_capped['convergence']['converged'] is False
    assert det_capped['convergence']['n_levels'] == 100
    assert det_capped['convergence']['rel_change_bulk'] > 1e-12


@pytest.mark.physics_invariant
def test_exobase_temperature_floors_at_the_profile_top():
    """A prescribed exobase temperature never builds a falling thermosphere.

    The extension is the inflated structure the exobase quantities must be
    read from, so it cannot end colder than the level it extends from: that
    puts the exobase more strongly bound than its own anchor, inverting the
    construction and biasing the branch toward retention. The prescribed
    temperature is a stand-in for physics the branch does not solve, and in a
    coupled run the profile top warms over secular time and can pass it, so
    the value floors at the anchor and the call is flagged rather than
    raising and stopping the run.
    """
    prof = _co2_hydrogen_profile(0.01)
    t_top = float(prof.T[-1])
    above = bates_extension(prof, M_MARS, 4.0 * t_top)
    assert above['t_exo_floored'] is False
    assert float(above['T'][-1]) > t_top
    below = bates_extension(prof, M_MARS, 0.25 * t_top)
    assert below['t_exo_floored'] is True
    # Floored, not falling: the extension is isothermal at the anchor value.
    assert float(below['T'][-1]) == pytest.approx(t_top, rel=1e-12, abs=0.0)
    assert float(below['T'][0]) == pytest.approx(t_top, rel=1e-12, abs=0.0)
    # The exobase is no more bound than the level it extends from.
    per, det = hydrostatic_rates(prof, M_MARS, 0.25 * t_top)
    assert det['flags'].get('t_exo_floored_to_profile_top') is True
    lam_anchor = G * M_MARS * det['m_bar'] / (kb * t_top * float(prof.r[-1]))
    lam_exo = G * M_MARS * det['m_bar'] / (kb * det['T_exo'] * det['r_exo'])
    assert lam_exo <= lam_anchor
    # At the anchor temperature exactly, nothing is flagged.
    _p, det_eq = hydrostatic_rates(prof, M_MARS, t_top)
    assert 't_exo_floored_to_profile_top' not in det_eq['flags']


@pytest.mark.physics_invariant
def test_per_element_shares_follow_the_species_that_escape():
    """The split names which element leaves, not just how much in total.

    Element rates summing to the bulk rate is a weak claim: the dispatcher
    renormalizes the split onto the bulk rate, so the sum matches by
    construction and a permuted mapping would conserve mass while moving the
    wrong elements out of the planet. What has to hold is the identity of the
    shares. On a carbon dioxide host carrying one percent hydrogen, hydrogen
    is the only species light enough to escape and carries the whole rate,
    twenty decades above the carbon and oxygen the heavy background supplies,
    and the CO2 that does leave splits onto carbon and oxygen in
    stoichiometric mass proportion.
    """
    prof = _co2_hydrogen_profile(0.01)
    per_el, det = hydrostatic_rates(prof, M_MARS, 1000.0)
    assert set(per_el) == {'H', 'C', 'O'}
    # Hydrogen carries the rate, and by twenty decades, so a permutation onto
    # carbon or oxygen cannot pass as a rounding difference.
    assert per_el['H'] > 1.0e19 * per_el['O']
    assert per_el['H'] > 1.0e19 * per_el['C']
    assert per_el['H'] == pytest.approx(det['per_species_rate']['H'], rel=1e-12, abs=0.0)
    # Oxygen above carbon, in the ratio the CO2 formula fixes: two oxygens of
    # 15.999 against one carbon of 12.011.
    assert per_el['O'] / per_el['C'] == pytest.approx(
        2.0 * species_mass_amu('O') / species_mass_amu('C'), rel=1e-9, abs=0.0
    )
    rate_co2 = det['per_species_rate']['CO2']
    m_co2 = species_mass_amu('CO2')
    assert per_el['C'] == pytest.approx(
        rate_co2 * species_mass_amu('C') / m_co2, rel=1e-9, abs=0.0
    )
    assert per_el['O'] == pytest.approx(
        rate_co2 * 2.0 * species_mass_amu('O') / m_co2, rel=1e-9, abs=0.0
    )


@pytest.mark.reference_pinned
def test_yelle_harmonic_mean_and_area_referral():
    """The Eq. 14 combination and the Eq. 15 area referral, at their scale.

    Yelle (2024, Icarus 416, 116099) combines the effusion flux and the
    diffusion-limited supply by the harmonic mean of Eq. (14),
    ``phi = phi_J phi_l / (phi_J + phi_l)``, and refers the flux from the
    exobase back to the anchor radius through the ``(r_x / r_0)^2`` factor of
    Eq. (15). Both are invisible on a state where the two fluxes differ by
    decades, since there the harmonic mean equals the smaller one and the
    referral is a fixed rescaling: a plain minimum reproduces the combination
    to one part in 1e5 on the default exobase temperature. This test therefore
    runs at the exobase temperature where the two fluxes cross, where the
    harmonic mean is 0.61 of the minimum and a substitution cannot hide.
    """
    prof = _co2_hydrogen_profile(0.01)
    # 130 K puts the hydrogen effusion flux within a factor 1.6 of its supply.
    _per, det = hydrostatic_rates(prof, M_MARS, 130.0)
    d = det['species']['H']
    phi_j, phi_l, phi = d['phi_jeans'], d['phi_diffusion'], d['phi_per_area_r0']
    assert 0.5 < phi_j / phi_l < 3.0, 'the crossover state has drifted'
    assert phi == pytest.approx(phi_j * phi_l / (phi_j + phi_l), rel=1e-12, abs=0.0)
    # Discrimination: a plain minimum is 1.6 times the harmonic mean here.
    assert phi < 0.75 * min(phi_j, phi_l)
    # The harmonic mean is below both arguments, always, which is the property
    # that makes it a supply cap rather than a blend.
    assert phi < phi_j and phi < phi_l
    # Area referral: the flux is quoted at the anchor, so it carries the ratio
    # of the exobase and anchor areas explicitly.
    r_x, r_0 = det['r_exo'], det['r_anchor']
    assert r_x > r_0
    referral = (r_x / r_0) ** 2
    assert referral > 1.05, 'the referral factor is too close to 1 to discriminate'
    # Rebuild the effusion flux from the reported exobase quantities: the
    # referral is the only factor between the local flux and the quoted one.
    # The mixing ratio is the diffusively enriched one at the exobase, not the
    # anchor value, which for a light trace species differs by decades.
    local = d['volkov_C'] * d['w_jeans'] * d['X_tilde_exo'] * det['n_exo']
    assert phi_j == pytest.approx(referral * local, rel=1e-9, abs=0.0)
