"""Tests for ``src/zephyrus/dispatcher.py``.

Exercises the assembled regime dispatcher end to end on synthetic
atmospheres, so the whole file carries the smoke tier (real code path, no
mocks). The properties under test:

- Totality: 200 random physically posed inputs return exactly one label, a
  finite non-negative bulk rate, per-species rates summing to it, and a
  populated diagnostics container, with no exception.
- Routing: an inflated light envelope dispatches to boil-off; a bound heavy
  atmosphere under weak XUV to hydrostatic; a light envelope under strong
  XUV to a hydrodynamic sub-label; a Roche-filling geometry to overflow.
- The Roche screen renames a state without changing its rate, so the
  dispatched rate is continuous across the overflow boundary, and its
  subflag separates the two published overflow geometries by the extent of
  the atmosphere rather than of the photosphere.
- Cross-implementation pin: the dispatcher's energy-limited candidate
  equals the released ``EL_escape`` at the same inputs.
- Boxedness: sabotaging every diagnostics producer with garbage stubs
  leaves the verdict, the rate, and the split unchanged, so no diagnostic
  feeds control flow.
- Error contract: malformed settings and physical states raise.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import math

import numpy as np
import pytest

from zephyrus.composition import species_mass_amu
from zephyrus.constants import G, amu, kb
from zephyrus.dispatcher import (
    REGIME_LABELS,
    DispatchSettings,
    EscapeInputs,
    dispatch,
)
from zephyrus.escape import EL_escape
from zephyrus.hydrodynamic import caldiroli_efficiency
from zephyrus.nozzle import isothermal_column_density, nozzle_candidate
from zephyrus.planets_parameters import Me, Ms, Re
from zephyrus.profiles import Profile, isothermal_profile, photospheric_level

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

AU = 1.496e11  # m

COMPOSITIONS = [
    {'H2': 1.0},
    {'H2': 0.9, 'He': 0.1},
    {'H2O': 1.0},
    {'N2': 1.0},
    {'CO2': 1.0},
    {'CO2': 0.7, 'N2': 0.3},
    {'H2': 0.5, 'H2O': 0.5},
    {'N2': 0.8, 'O2': 0.2},
]


def _inputs(M_p, R_p, T_eq, comp, F_xuv, a=0.1 * AU, e=0.0, p_surf=1e7, p_top=1e-5, **kw):
    """Physically posed dispatch inputs on an isothermal test atmosphere."""
    prof = isothermal_profile(M_p, R_p, T_eq, comp, p_surf, p_top)
    defaults = dict(
        M_p=M_p,
        R_p=R_p,
        M_star=Ms,
        a=a,
        e=e,
        T_eq=T_eq,
        F_xuv=F_xuv,
        F_bol=1e4,
        F_int=0.5,
        kappa_photo=0.01,
        profile=prof,
    )
    defaults.update(kw)
    return EscapeInputs(**defaults)


def _inverted_profile(M_p, R_p, T_base, T_top, comp, p_surf=1e7, p_top=1e-5, n=160):
    """Hydrostatic profile whose temperature rises with altitude.

    The isothermal builder makes the profile temperature identical to
    ``T_eq`` at every level, so no test on it can tell a profile
    temperature from an equilibrium one, and no test on it can see the
    launch-level cancellation fail. This integrates the same hydrostatic
    relation on a temperature that ramps linearly in log pressure, which is
    the shape a real upper atmosphere has and the shape that breaks the
    isothermal Bernoulli argument.
    """
    tot = sum(comp.values())
    mu = sum(x * species_mass_amu(sp) for sp, x in comp.items()) / tot * amu
    lnp = np.linspace(math.log(p_surf), math.log(p_top), n)
    frac = (lnp - lnp[0]) / (lnp[-1] - lnp[0])
    T = T_base + (T_top - T_base) * frac
    r = np.empty(n)
    r[0] = R_p
    last = n - 1
    for i in range(n - 1):
        H = kb * T[i] * r[i] ** 2 / (G * M_p * mu)
        r[i + 1] = r[i] - H * (lnp[i + 1] - lnp[i])
        if G * M_p * mu / (kb * T[i + 1] * r[i + 1]) < 2.2:
            # The same boundedness stop the isothermal builder uses: past
            # it the scale height grows faster than the radius and the
            # integration runs away instead of describing an atmosphere.
            last = i
            break
    lnp, r, T = lnp[: last + 1], r[: last + 1], T[: last + 1]
    n = last + 1
    vmr = {sp: np.full(n, x / tot) for sp, x in comp.items()}
    return Profile(p=np.exp(lnp), r=r, T=T, vmr=vmr, mmw=np.full(n, mu), kzz=None)


def _random_inputs(rng):
    """One random draw over masses, radii, temperatures, and compositions."""
    M_p = rng.uniform(0.5, 20.0) * Me
    R_p = rng.uniform(0.7, 4.0) * Re
    T_eq = rng.uniform(200.0, 2500.0)
    comp = COMPOSITIONS[rng.integers(len(COMPOSITIONS))]
    p_surf = 10 ** rng.uniform(4.0, 9.0)
    p_top = 10 ** rng.uniform(-6.0, -2.0)
    try:
        prof = isothermal_profile(M_p, R_p, T_eq, comp, p_surf, p_top)
    except ValueError:
        return None  # unbound at the surface: not a physically posed input
    return EscapeInputs(
        M_p=M_p,
        R_p=R_p,
        M_star=Ms * rng.uniform(0.3, 1.5),
        a=10 ** rng.uniform(math.log10(0.01), 0.0) * AU,
        e=rng.uniform(0.0, 0.5),
        T_eq=T_eq,
        F_xuv=0.0 if rng.random() < 0.05 else 10 ** rng.uniform(-3.0, 4.0),
        F_bol=10 ** rng.uniform(2.0, 6.0),
        F_int=10 ** rng.uniform(-2.0, 2.0),
        kappa_photo=10 ** rng.uniform(-3.0, -1.0),
        profile=prof,
        prev_regime=None
        if rng.random() < 0.7
        else ['boiloff', 'hydrodynamic:EL', 'hydrostatic'][rng.integers(3)],
        reservoirs=None
        if rng.random() < 0.5
        else {'H': 10 ** rng.uniform(15.0, 20.0), 'O': 10 ** rng.uniform(15.0, 20.0)},
        age=None if rng.random() < 0.5 else 10 ** rng.uniform(14.0, 17.5),
    )


@pytest.mark.physics_invariant
def test_totality_over_random_physical_inputs():
    """Every physically posed input returns one consistent, finite result.

    200 random draws across compositions, masses, radii, fluxes (including
    exactly zero), eccentricities, and profile depths: exactly one known
    label, a finite non-negative bulk rate, finite non-negative per-species
    rates summing to the bulk rate, and a populated diagnostics container,
    with no exception anywhere. This is the conservation and boundedness
    contract of the whole dispatcher.
    """
    rng = np.random.default_rng(42)
    n_ok = 0
    seen = set()
    while n_ok < 200:
        inp = _random_inputs(rng)
        if inp is None:
            continue
        res = dispatch(inp)
        assert res.regime in REGIME_LABELS, res.regime
        assert math.isfinite(res.mdot)
        assert res.mdot >= 0.0
        assert all(math.isfinite(v) and v >= 0.0 for v in res.per_species.values())
        tot = sum(res.per_species.values())
        if res.mdot > 0.0:
            assert tot == pytest.approx(res.mdot, rel=1e-6, abs=0.0)
        else:
            # A zero bulk rate is still a conservation statement, and it is
            # the one a relative tolerance cannot make: nothing may leave.
            # Roughly a quarter of these draws land here, on states so
            # strongly bound that every branch underflows, and skipping them
            # left the split unchecked exactly where it is cheapest to break.
            assert res.mdot == 0.0
            assert tot == 0.0
            assert all(v == 0.0 for v in res.per_species.values())
        assert isinstance(res.diagnostics, dict)
        assert 'knudsen' in res.diagnostics
        seen.add(res.regime)
        n_ok += 1
    # The sweep must genuinely exercise more than one branch.
    assert len(seen) >= 3, seen


def test_routing_boiloff_for_inflated_light_envelope():
    """A hot, loosely bound H2 envelope dispatches to boil-off.

    The restricted Jeans parameter sits below the activation threshold, so
    the bolometric branch owns the point: label ``boiloff``, a positive
    rate even at zero XUV flux (the driver is bolometric, not XUV), and an
    unfractionated split (fractionation stays off on this branch).
    """
    inp = _inputs(2 * Me, 2.2 * Re, 1800.0, {'H2': 1.0}, F_xuv=0.0, a=0.05 * AU)
    res = dispatch(inp)
    assert res.diagnostics['lambda_gate'] < 20.0
    assert res.regime == 'boiloff'
    assert res.mdot > 0.0
    assert 'closure' not in res.diagnostics  # no fractionation on this branch


def test_routing_hydrostatic_for_bound_heavy_atmosphere():
    """A massive CO2 atmosphere under weak XUV dispatches to hydrostatic.

    The wind's sonic point is rarefied (Knudsen number far above the
    switch), the exobase is cool against both escape temperatures, and the
    rates are per-species with the lower-limit flag: non-thermal channels
    are absent from the hydrostatic branch.
    """
    inp = _inputs(10 * Me, 1.8 * Re, 800.0, {'CO2': 1.0}, F_xuv=1.0, a=0.5 * AU)
    res = dispatch(inp)
    assert res.regime == 'hydrostatic'
    assert res.diagnostics['knudsen']['kn_sc'] > res.diagnostics['knudsen']['threshold_applied']
    assert res.flags.get('hydrostatic_lower_limit') is True
    assert res.diagnostics['hydrostatic']['gate_unstable'] is False
    # Cool against both conventions: the point is not contested either.
    assert 'contested_ion' not in res.flags


@pytest.mark.reference_pinned
def test_routing_hydrodynamic_and_el_candidate_matches_el_escape():
    """A strongly irradiated light envelope goes hydrodynamic, EL-consistent.

    The dispatcher's energy-limited candidate must equal the released
    ``EL_escape`` evaluated with the same efficiency, radii, flux, and
    tidal factor (the cross-implementation pin tying the dispatcher to the
    package's public energy-limited contract). The label carries the
    min(EL, RR) winner as its sub-label.
    """
    inp = _inputs(5 * Me, 1.8 * Re, 1100.0, {'H2': 0.9, 'He': 0.1}, F_xuv=200.0, a=0.05 * AU)
    res = dispatch(inp)
    assert res.regime.startswith('hydrodynamic')
    hydro = res.diagnostics['hydrodynamic']
    # Reconstruct the EL candidate through the public entry point.
    eps = hydro['efficiency']
    r_xuv = _photo_radius(inp)
    ref = EL_escape(True, inp.a, inp.e, inp.M_p, inp.M_star, eps, inp.R_p, r_xuv, inp.F_xuv, 2)
    assert hydro['mdot_el'] == pytest.approx(ref, rel=1e-9, abs=0.0)
    assert res.mdot == pytest.approx(min(hydro['mdot_el'], hydro['mdot_rr']), rel=1e-9, abs=0.0)
    winner = 'EL' if hydro['mdot_el'] <= hydro['mdot_rr'] else 'RR'
    assert res.regime == f'hydrodynamic:{winner}'


def _photo_radius(inp):
    """The 20 mbar photospheric radius of an input's profile."""
    from zephyrus.profiles import photospheric_level

    lev, _ = photospheric_level(inp.profile, inp.settings.P_photo)
    return lev['r']


def test_routing_roche_overflow_inside_the_hill_sphere():
    """A planet whose Hill sphere sits inside its radius overflows, flagged.

    At 0.003 au the periapsis Hill radius of a 5 Earth-mass planet drops
    below its own radius: the label is ``roche_overflow`` with the
    dynamical subflag, and the rate comes from the Bondi-capped bolometric
    machinery at the overflow geometry (finite and non-negative).
    """
    inp = _inputs(5 * Me, 1.5 * Re, 1500.0, {'H2': 1.0}, F_xuv=100.0, a=0.003 * AU)
    res = dispatch(inp)
    assert res.diagnostics['roche']['xi_ktide'] < 1.0
    assert res.regime == 'roche_overflow'
    assert res.flags.get('roche_subflag') == 'dynamical'
    assert math.isfinite(res.mdot)
    assert res.mdot >= 0.0


def test_roche_screen_renames_without_changing_the_rate():
    """Crossing the overflow boundary changes the label and not the rate.

    The screen's boundary is a geometric criterion on the winning branch's
    flow radius, so the two sides of the boundary hold the same branch and
    the dispatched rate must be continuous across it. The family here is a
    luminosity-capped bolometric residual whose sonic radius crosses the
    Hill radius as the orbit widens, chosen because every other candidate
    stays subdominant across the bracket: the XUV flux is negligible and
    the nozzle candidate sits outside its applicability criterion, so the
    rename is the only thing that changes. Bisecting in orbital distance
    brackets the label change; the rates on either side agree to machine
    precision and both equal the capped residual. Discrimination:
    substituting the Bondi-capped bolometric rate at the overflow geometry,
    which is what a rate-changing screen returns, is more than a decade
    larger, because that form bypasses the luminosity cap.
    """
    comp = {'H2': 0.9, 'He': 0.1}

    def at(a):
        return dispatch(_inputs(3 * Me, 2.0 * Re, 1000.0, comp, F_xuv=0.1, a=a))

    lo, hi = 0.078 * AU, 0.3 * AU
    inner = at(lo)
    assert inner.regime == 'roche_overflow'
    assert inner.diagnostics['nozzle']['applicable'] is False
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if at(mid).regime == inner.regime:
            lo = mid
        else:
            hi = mid
    below, above = at(lo), at(hi)
    assert below.regime == 'roche_overflow'
    assert above.regime == 'boiloff'
    assert below.mdot == pytest.approx(above.mdot, rel=1e-9, abs=0.0)
    bolo = below.diagnostics['bolometric']
    assert below.mdot == pytest.approx(bolo['mdot_luminosity'], rel=1e-9, abs=0.0)
    assert below.diagnostics['roche']['rate_branch'] == above.regime
    overflow_geometry_rate = min(bolo['mdot_parker'], bolo['mdot_bondi'])
    assert overflow_geometry_rate > 10.0 * below.mdot


def test_roche_subflag_separates_the_two_geometries():
    """The subflag reads the atmosphere's extent, not the photosphere's.

    Owen & Jackson (2012) separate dynamical overflow, where the atmosphere
    itself reaches the lobe, from the narrow band where only the would-be
    sonic surface does. A Mars-mass CO2 planet at 0.028 au with a 2000 K
    exobase has its extended structure outside the Hill radius while its
    photosphere sits far inside, so the subflag is dynamical even though the
    Hill sphere still encloses the planet many times over; an inflated
    hydrogen envelope whose Parker sonic radius alone passes the Hill radius
    gets the other subflag, with its atmosphere well inside.
    """
    dyn = dispatch(
        _inputs(
            0.107 * Me,
            0.53 * Re,
            1600.0,
            {'CO2': 1.0},
            F_xuv=1e-4,
            a=0.028 * AU,
            settings=DispatchSettings(T_exo_value=2000.0),
        )
    )
    roche = dyn.diagnostics['roche']
    assert dyn.regime == 'roche_overflow'
    assert dyn.flags.get('roche_subflag') == 'dynamical'
    assert roche['xi_ktide'] > 1.0  # not the trivial planet-inside-its-lobe case
    assert roche['r_atmosphere'] > roche['R_hill_periapsis']

    # The second subflag: an atmosphere inside its own Roche lobe whose
    # would-be sonic surface sits outside the Hill radius, which is the
    # narrow band Owen & Jackson (2012) describe. The comparator is the
    # lobe rather than the Hill radius, because the lobe is the critical
    # surface and sits about 0.70 of the way out to the Hill radius.
    bound = dispatch(
        _inputs(
            0.5 * Me,
            2.0 * Re,
            1000.0,
            {'CO2': 1.0},
            F_xuv=1e-2,
            a=0.03 * AU,
            settings=DispatchSettings(T_exo_value=2000.0),
        )
    )
    roche_b = bound.diagnostics['roche']
    assert bound.regime == 'roche_overflow'
    assert bound.flags.get('roche_subflag') == 'no_transonic'
    assert roche_b['rate_branch'].startswith('hydrodynamic')
    assert roche_b['xi_flow'] <= 1.0  # the sonic surface is outside the Hill radius
    assert roche_b['r_atmosphere'] < bound.diagnostics['nozzle']['r_lobe']  # the gas is not
    assert roche_b['xi_ktide'] > 1.0

    # The third: the label won by the rate crossing alone, with neither
    # the atmosphere nor the flow radius reaching out.
    neither = dispatch(
        _inputs(
            0.107 * Me,
            0.53 * Re,
            1000.0,
            {'CO2': 1.0},
            F_xuv=1e-4,
            a=0.02 * AU,
            settings=DispatchSettings(T_exo_value=2000.0),
        )
    )
    roche_n = neither.diagnostics['roche']
    assert neither.regime == 'roche_overflow'
    assert roche_n['rate_branch'] == 'roche_overflow'
    assert neither.flags.get('roche_subflag') == 'neither'
    assert roche_n['r_atmosphere'] < neither.diagnostics['nozzle']['r_lobe']
    assert roche_n['xi_flow'] > 1.0


def test_nozzle_win_relabels_with_the_transfer_rate():
    """A nozzle win labels ``roche_overflow`` and carries a real transfer rate.

    A puffy sub-Neptune whose photosphere sits within a few thermal units
    of its lobe dispatches the Jackson et al. (2017) L1 rate: the label is
    ``roche_overflow``, the branch is ``roche_overflow`` too, the subflag reads
    ``dynamical`` because this envelope's extended structure itself reaches
    past the Hill sphere and the geometric reading takes precedence, the
    split is unfractionated (the nozzle is a bulk photospheric flow, not
    the closure's wind base), and ``near_roche`` stays down because it
    warns about the tidal inflation of a bound rate, which this is not.
    """
    res = dispatch(
        _inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU)
    )
    noz = res.diagnostics['nozzle']
    assert res.regime == 'roche_overflow'
    assert res.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    assert res.flags.get('roche_subflag') == 'dynamical'
    assert (
        res.diagnostics['roche']['r_atmosphere'] > res.diagnostics['roche']['R_hill_periapsis']
    )
    assert res.mdot == pytest.approx(noz['rate_kg_s'], rel=1e-12, abs=0.0)
    assert noz['applicable'] is True
    assert 'near_roche' not in res.flags
    assert 'closure' not in res.diagnostics
    assert sum(res.per_species.values()) == pytest.approx(res.mdot, rel=1e-9, abs=0.0)


def test_nozzle_competes_only_inside_its_domain():
    """The nozzle candidate needs applicability and a non-empty rate to win.

    Two refusals, one per condition. A bound CO2 planet has a nozzle rate
    below the one-proton-per-Julian-year floor, so a crossing against the
    similarly empty hydrodynamic rate decides nothing and the verdict
    stands. A residual-driven sub-Neptune at a wide orbit has a nozzle
    candidate above its own dispatched rate, but the isothermal sonic
    radius sits inside the L1 distance (Jackson et al. 2017, their
    Figure 9), so a spherical wind chokes first, the candidate reports
    without competing, and the bolometric verdict stands.
    """
    empty = dispatch(_inputs(Me, Re, 700.0, {'CO2': 1.0}, F_xuv=10.0, a=0.1 * AU))
    noz = empty.diagnostics['nozzle']
    assert noz['applicable'] is True
    assert 0.0 <= noz['rate_kg_s'] < empty.diagnostics['rate_floor']['floor_kg_s']
    assert empty.regime == 'hydrodynamic:EL'

    outside = dispatch(
        _inputs(3 * Me, 2.0 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=0.1, a=0.12 * AU)
    )
    noz_o = outside.diagnostics['nozzle']
    assert noz_o['applicable'] is False
    assert noz_o['R_sonic_over_R_L1'] < 1.0
    # The unguarded Eq. (3) average would have won; the duty-cycled rate
    # the dispatcher competes is zero, because no arc of the orbit is
    # inside the overflow description.
    assert noz_o['rate_full_orbit_kg_s'] > outside.mdot
    assert noz_o['rate_kg_s'] == 0.0
    assert noz_o['applicable_orbit_fraction'] == 0.0
    assert outside.regime == 'boiloff'


def test_nozzle_crossing_is_continuous():
    """The dispatched rate is continuous across the nozzle rate crossing.

    Where the nozzle candidate overtakes the standing branch rate the label
    changes because two rates cross, so the dispatched rate on either side
    of the bisected boundary is the same number: the boundary is a rate
    crossing, unlike the criterion boundaries (the activation gate, the
    Knudsen switch, the applicability edge), whose jumps are results.
    """
    prof_settings = DispatchSettings(T_exo_value=3000.0)

    def at(a):
        return dispatch(
            _inputs(
                0.107 * Me,
                0.53 * Re,
                600.0,
                {'CO2': 1.0},
                F_xuv=1e-4,
                a=a,
                settings=prof_settings,
            )
        )

    lo, hi = 0.010 * AU, 0.014 * AU
    inner = at(lo)
    assert inner.regime == 'roche_overflow'
    assert inner.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    # This donor's extended structure passes its own Roche lobe while
    # staying inside the Hill sphere, so the geometric subflag reads
    # dynamical whichever candidate carries the rate.
    assert inner.flags.get('roche_subflag') == 'dynamical'
    assert inner.diagnostics['roche']['r_atmosphere'] >= inner.diagnostics['nozzle']['r_lobe']
    assert at(hi).regime == 'hydrostatic'
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if at(mid).regime == inner.regime:
            lo = mid
        else:
            hi = mid
    below, above = at(lo), at(hi)
    assert below.regime == 'roche_overflow'
    assert above.regime == 'hydrostatic'
    assert below.mdot == pytest.approx(above.mdot, rel=1e-6, abs=0.0)


def test_nozzle_saturation_flag_marks_lobe_contact():
    """A launch level at or beyond the lobe wins at the boundary value, flagged.

    A loose hydrogen envelope whose photospheric level sits far outside its
    shrunken lobe dispatches the saturated nozzle rate with
    ``nozzle_saturated`` raised; the unsaturated nozzle win of the puffy
    sub-Neptune family must not raise it.
    """
    saturated = dispatch(
        _inputs(
            2 * Me,
            2.5 * Re,
            1500.0,
            {'H2': 0.9, 'He': 0.1},
            F_xuv=13.4,
            a=0.015 * AU,
            p_top=1e3,
        )
    )
    assert saturated.regime == 'roche_overflow'
    assert saturated.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    assert saturated.flags.get('nozzle_saturated') is True
    assert saturated.diagnostics['nozzle']['saturated'] is True
    assert math.isfinite(saturated.mdot)

    unsaturated = dispatch(
        _inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU)
    )
    assert unsaturated.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    assert 'nozzle_saturated' not in unsaturated.flags


def test_nozzle_orbit_average_on_eccentric_wins_only():
    """An eccentric nozzle win dispatches a duty-cycled orbit average.

    The primary treats a circular, synchronously rotating donor, so both
    the quasi-static evaluation at each separation and the duty cycle over
    the arc where the overflow description applies are this module's
    conventions, and the flags make them visible on the result they
    produced. A circular win must raise neither, and its average must be
    the instantaneous rate exactly, so the convention costs nothing where
    it does not apply.
    """
    ecc = dispatch(
        _inputs(
            3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU, e=0.3
        )
    )
    noz = ecc.diagnostics['nozzle']
    assert ecc.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    assert ecc.flags.get('nozzle_orbit_averaged') is True
    # The periapsis separation is the one the geometry is built on, and it
    # is not the semi-major axis: a rate taken at ``a`` instead would sit
    # 31% below the periapsis value on this state.
    assert noz['a_periapsis'] == pytest.approx(0.7 * 0.07 * AU, rel=1e-12)
    # The applicable arc surrounds periapsis, because the L1 distance grows
    # with separation while the sonic radius does not.
    assert noz['R_sonic_over_R_L1'] > 1.0 > noz['R_sonic_over_R_L1_apoapsis']
    assert 0.0 < noz['applicable_orbit_fraction'] < 1.0
    assert ecc.flags.get('nozzle_partial_orbit') is True
    # Duty-cycled below the full-orbit average, which is itself below the
    # periapsis rate: this donor is unsaturated, so the rate falls with
    # separation and periapsis is the richest phase.
    assert noz['rate_kg_s'] < noz['rate_full_orbit_kg_s'] < noz['rate_periapsis_kg_s']
    assert ecc.mdot == pytest.approx(noz['rate_kg_s'], rel=1e-12)

    circ = dispatch(
        _inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU)
    )
    noz_c = circ.diagnostics['nozzle']
    assert circ.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    assert 'nozzle_orbit_averaged' not in circ.flags
    assert 'nozzle_partial_orbit' not in circ.flags
    # At most one by construction, so this is the whole-orbit claim without
    # an equality against a float literal.
    assert noz_c['applicable_orbit_fraction'] >= 1.0
    # Every phase of a circular orbit is the same phase, so the quadrature
    # returns the instantaneous rate and the convention is free.
    assert noz_c['rate_kg_s'] == pytest.approx(noz_c['rate_periapsis_kg_s'], rel=1e-12)
    assert noz_c['rate_apoapsis_kg_s'] == pytest.approx(noz_c['rate_periapsis_kg_s'], rel=1e-12)


@pytest.mark.physics_invariant
def test_wind_launch_level_cancels_along_the_wind_column():
    """In wind mode the launch level cancels along the wind's own column.

    The launch-level convention rests on the Bernoulli invariance of
    ``rho exp(Phi / v_th^2)``, which holds only when the density and the
    sound speed belong to one column. The wind setting launches from the
    wind base at the wind's temperature, so moving the level along the
    isothermal column through that anchor must leave the rate alone. The
    guard matters because taking the density from the profile's own
    (far colder) structure instead moves the rate by more than two decades
    over the same range of levels, which is what the invariance claim
    would otherwise be hiding.

    The base is placed by hand rather than by the Lopez default, which on
    this planet puts the wind base at 1.29 lobe radii: outside the lobe
    the exponent is clamped and the rate is the lobe-filling boundary
    value, which is linear in the launch density and invariant along no
    column at all.
    """
    state = _inputs(
        3 * Me,
        2.2 * Re,
        1000.0,
        {'H2': 0.9, 'He': 0.1},
        F_xuv=13.4,
        a=0.07 * AU,
        settings=DispatchSettings(
            nozzle_temperature='wind', base_method='fixed_pressure', P_base_fixed=5.0
        ),
    )
    noz = dispatch(state).diagnostics['nozzle']
    assert not noz['saturated']
    r_ref, rho_ref, v_th = noz['r_launch'], noz['rho_launch'], noz['v_th']
    reference = noz['rate_full_orbit_kg_s']

    # The discriminator is the same closed form at the wrong sound speed:
    # the profile's own, which is what the launch state carried before the
    # column was made consistent with the barrier.
    v_cold = math.sqrt(kb * state.T_eq / noz['mu_kg'])
    on_column, cold_column = [], []
    for factor in (0.6, 0.8, 1.5, 2.5):
        r = r_ref * factor
        for bucket, speed in ((on_column, v_th), (cold_column, v_cold)):
            rho = isothermal_column_density(rho_ref, r_ref, r, state.M_p, speed)
            bucket.append(
                nozzle_candidate(
                    state.M_p,
                    state.M_star,
                    state.a,
                    state.e,
                    rho,
                    r,
                    noz['T_K'],
                    noz['mu_kg'],
                )[0]
            )
    for value in on_column:
        assert value == pytest.approx(reference, rel=0.01)
    assert max(cold_column) / min(cold_column) > 5.0


def test_nozzle_temperature_setting_selects_and_validates():
    """The nozzle temperature setting moves the diagnostic and validates.

    The default evaluates the sound speed at the photospheric level (the
    primary's construction); ``wind`` evaluates it at the thermostat's wind
    state; anything else is rejected by the settings validator.
    """
    state = dict(F_xuv=13.4, a=0.07 * AU)
    photo = dispatch(_inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, **state))
    assert photo.diagnostics['nozzle']['temperature_mode'] == 'photospheric'
    assert photo.diagnostics['nozzle']['T_K'] == pytest.approx(1000.0, rel=1e-9, abs=0.0)
    wind = dispatch(
        _inputs(
            3 * Me,
            2.2 * Re,
            1000.0,
            {'H2': 0.9, 'He': 0.1},
            settings=DispatchSettings(nozzle_temperature='wind'),
            **state,
        )
    )
    assert wind.diagnostics['nozzle']['temperature_mode'] == 'wind'
    assert wind.diagnostics['nozzle']['T_K'] == pytest.approx(
        wind.diagnostics['hydrodynamic']['T_wind'], rel=1e-9, abs=0.0
    )
    with pytest.raises(ValueError):
        DispatchSettings(nozzle_temperature='photosphere').validate()


def test_nozzle_power_diagnostic_reports_the_lift_cost():
    """The lift power travels beside the two luminosities on every call.

    The nozzle carries no energy cap, so the diagnostic is what shows where
    the isothermal assumption is strained: the power to lift the dispatched
    flow to L1, against the interior luminosity and the intercepted
    instellation, present and finite whether the candidate won or lost.
    """
    for res in (
        dispatch(
            _inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU)
        ),
        dispatch(_inputs(Me, Re, 700.0, {'CO2': 1.0}, F_xuv=10.0, a=0.1 * AU)),
    ):
        noz = res.diagnostics['nozzle']
        for key in ('power_lift_W', 'L_int_W', 'L_bol_intercepted_W'):
            assert math.isfinite(noz[key])
            assert noz[key] >= 0.0
        assert noz['L_int_W'] > 0.0
        assert noz['L_bol_intercepted_W'] > 0.0


def test_diagnostics_are_boxed(monkeypatch):
    """Sabotaging every diagnostics producer changes no dispatch outcome.

    The container is reporting only: no control flow reads it back. The test
    proves it by replacing every diagnostics-side producer with a stub
    returning garbage of the right shape and asserting the regime, the bulk
    rate, and every per-species rate are unchanged. It runs on one state per
    branch, each dispatching a rate of order unity or above, because a state
    whose rate sits near the denormal floor would compare equal to anything
    under any absolute tolerance and the assertion would not discriminate.
    """
    states = {
        'hydrostatic': _inputs(
            0.107 * Me, 0.53 * Re, 440.0, {'CO2': 0.99, 'H2': 0.01}, F_xuv=0.01, a=0.2 * AU
        ),
        'hydrodynamic:EL': _inputs(Me, Re, 1000.0, {'N2': 1.0}, F_xuv=100.0),
        'hydrodynamic:RR': _inputs(Me, Re, 1000.0, {'CO2': 1.0}, F_xuv=5000.0),
        'boiloff': _inputs(
            Me, 1.5 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=10.0, a=0.0775 * AU
        ),
        'roche_overflow': _inputs(
            3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU
        ),
    }
    reference = {k: dispatch(v) for k, v in states.items()}
    for expected, res in zip(states, reference.values()):
        assert res.regime == expected, (expected, res.regime)
        assert res.mdot > 1.0, (expected, res.mdot)
    nan = float('nan')
    monkeypatch.setattr('zephyrus.diagnostics.q_net_over_qc', lambda *a, **k: (nan, nan, nan))
    monkeypatch.setattr('zephyrus.diagnostics.guo_triple', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.erkaev_tc', lambda *a, **k: nan)
    monkeypatch.setattr('zephyrus.diagnostics.potential_screens', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.along_profile_fluid_check', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.self_consistency_screen', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.rate_floor_screen', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.boiloff.tang_timescale_check', lambda *a, **k: {})
    for name, inp in states.items():
        sabotaged = dispatch(inp)
        ref = reference[name]
        assert sabotaged.regime == ref.regime, name
        assert sabotaged.mdot == pytest.approx(ref.mdot, rel=1e-12, abs=0.0), name
        assert set(sabotaged.per_species) == set(ref.per_species), name
        for el, v in ref.per_species.items():
            assert sabotaged.per_species[el] == pytest.approx(v, rel=1e-12, abs=0.0), (name, el)
        assert sabotaged.flags == ref.flags, name
        # The sabotage genuinely reached the container.
        assert sabotaged.diagnostics['guo_triple'] == {}


def test_hysteresis_memory_moves_the_threshold():
    """A previous regime label activates the hysteresis window, flagged.

    With no memory the sharp threshold applies; with a hydrodynamic
    previous label the applied threshold rises by the window factor, and
    the flag records that the memory was consumed.
    """
    inp = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0)
    sharp = dispatch(inp)
    assert 'hysteresis_active' not in sharp.flags
    assert sharp.diagnostics['knudsen']['threshold_applied'] == pytest.approx(1.0)
    inp2 = _inputs(
        5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, prev_regime='hydrodynamic:EL'
    )
    remembered = dispatch(inp2)
    assert remembered.flags.get('hysteresis_active') is True
    assert remembered.diagnostics['knudsen']['threshold_applied'] == pytest.approx(1.5)


def test_base_out_of_range_extend_mode():
    """A profile too shallow for the wind base extends instead of clamping.

    On a profile truncated at 0.01 Pa the physical Lopez base lies above
    the top: the default policy clamps (flagged, with the distance in
    decades); the extend policy evaluates the base on the extended upper
    structure instead, replacing the clamp flag with ``base_extended`` and
    placing the base at a lower pressure than the profile top. Either way
    the diagnostics report the pressure the base method asked for before
    any clamp, which is the only place that quantity is available.
    """
    inp_clamp = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, p_top=1e-2)
    res_clamp = dispatch(inp_clamp)
    assert res_clamp.flags.get('base_clamped') is True
    assert res_clamp.flags['base_clamp_decades'] > 0.0
    # The clamped base sits at the profile top, above its own target, and
    # the recorded distance is the gap between the two. A diagnostic that
    # echoed the clamped level instead would give a zero distance here.
    base_clamp = res_clamp.diagnostics['base_level']
    assert base_clamp['p_Pa'] == pytest.approx(1e-2, rel=1e-9, abs=0.0)
    assert base_clamp['p_physical_Pa'] < base_clamp['p_Pa']
    decades = np.log10(base_clamp['p_Pa'] / base_clamp['p_physical_Pa'])
    assert decades == pytest.approx(res_clamp.flags['base_clamp_decades'], rel=1e-9, abs=0.0)
    settings = DispatchSettings(base_out_of_range='extend')
    inp_ext = _inputs(
        5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, p_top=1e-2, settings=settings
    )
    res_ext = dispatch(inp_ext)
    assert res_ext.flags.get('base_extended') is True
    assert 'base_clamped' not in res_ext.flags
    assert res_ext.diagnostics['base_level']['p_Pa'] < 1e-2
    # On the extension the base reaches its target, so the two agree.
    base_ext = res_ext.diagnostics['base_level']
    assert base_ext['p_physical_Pa'] == pytest.approx(base_ext['p_Pa'], rel=1e-9, abs=0.0)


def test_fractionation_toggle_and_split_protocol():
    """The hydrodynamic split follows the fractionation toggle.

    With fractionation on, a hydrodynamic verdict carries the closure
    diagnostics and per-element rates that need not follow the base mass
    fractions; with it off, the split is unfractionated (flagged when it
    falls back to the base composition), and both sum to the bulk rate.
    """
    on = dispatch(
        _inputs(5 * Me, 1.8 * Re, 1100.0, {'H2': 0.9, 'He': 0.1}, F_xuv=200.0, a=0.05 * AU)
    )
    assert on.regime.startswith('hydrodynamic')
    assert 'closure' in on.diagnostics
    assert sum(on.per_species.values()) == pytest.approx(on.mdot, rel=1e-6, abs=0.0)
    settings = DispatchSettings(fractionate=False)
    off = dispatch(
        _inputs(
            5 * Me,
            1.8 * Re,
            1100.0,
            {'H2': 0.9, 'He': 0.1},
            F_xuv=200.0,
            a=0.05 * AU,
            settings=settings,
        )
    )
    assert off.regime.startswith('hydrodynamic')
    assert 'closure' not in off.diagnostics
    assert off.flags.get('split_from_base_composition') is True
    assert sum(off.per_species.values()) == pytest.approx(off.mdot, rel=1e-6, abs=0.0)


def test_settings_and_inputs_error_contract():
    """Malformed settings and physical states raise before any physics runs.

    Unknown option strings, all cooling channels off, non-positive masses,
    and an eccentricity outside [0, 1) all raise ``ValueError``; the same
    state with the defect repaired dispatches normally.
    """
    with pytest.raises(ValueError, match='base_method'):
        DispatchSettings(base_method='nonsense').validate()
    with pytest.raises(ValueError, match='cooling'):
        DispatchSettings(
            cool_atomic=False,
            cool_co2_band=False,
            cool_o_finestructure=False,
            cool_recombination=False,
        ).validate()
    # Numeric knobs are bounded too. Unbounded, an out-of-domain value either
    # surfaced as a bare math domain error from inside a branch or, for a
    # negative efficiency, changed the regime label with nothing said.
    for field, value, expected in (
        ('gamma_wind', 2.0, 'gamma_wind'),
        ('T_exo_value', -100.0, 'T_exo_value'),
        ('efficiency', -0.5, 'efficiency'),
        ('efficiency', 5.0, 'efficiency'),
        ('kn_crit', -1.0, 'kn_crit'),
        ('kzz', -300.0, 'kzz'),
        ('lambda_crit', -5.0, 'lambda_crit'),
        ('gamma_bates', 0.0, 'gamma_bates'),
        ('kn_hysteresis', 0.5, 'kn_hysteresis'),
        ('P_photo', float('nan'), 'P_photo'),
    ):
        with pytest.raises(ValueError, match=expected):
            DispatchSettings(**{field: value}).validate()
    # The physical edges of each range stay legal: an isothermal wind, the
    # monatomic index that is the domain limit of the sonic scale height, a
    # disabled hysteresis window, and full efficiency.
    for legal in (
        {'gamma_wind': 1.0},
        {'gamma_wind': 5.0 / 3.0},
        {'kn_hysteresis': 1.0},
        {'efficiency': 1.0},
        {'efficiency': 0.6},
    ):
        DispatchSettings(**legal).validate()
    good = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0)
    bad_e = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, e=1.0)
    with pytest.raises(ValueError, match='e must be'):
        dispatch(bad_e)
    bad_m = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0)
    bad_m.M_p = -1.0
    with pytest.raises(ValueError, match='M_p'):
        dispatch(bad_m)
    res = dispatch(good)
    assert res.regime in REGIME_LABELS


def test_caldiroli_efficiency_mode_applies_and_falls_back():
    """The fitted-efficiency mode applies inside its box and falls back below.

    Above the fit's flux-to-density validity bound the dispatched
    energy-limited candidate carries the converted fitted efficiency, which
    must differ resolvably from the fixed default; below the bound the fit
    is rejected and the dispatcher falls back to the fixed efficiency with
    both flags recorded.
    """
    settings = DispatchSettings(efficiency_mode='caldiroli')
    strong = dispatch(
        _inputs(
            5 * Me,
            1.8 * Re,
            1100.0,
            {'H2': 0.9, 'He': 0.1},
            F_xuv=200.0,
            a=0.05 * AU,
            settings=settings,
        )
    )
    eff = strong.diagnostics['hydrodynamic']['efficiency']
    assert 'efficiency_fallback_fixed' not in strong.flags
    assert eff != pytest.approx(0.1, rel=0.05, abs=0.0)
    assert 0.0 < eff < 1.0
    weak = dispatch(
        _inputs(
            5 * Me,
            1.8 * Re,
            1100.0,
            {'H2': 0.9, 'He': 0.1},
            F_xuv=0.1,
            a=0.05 * AU,
            settings=DispatchSettings(efficiency_mode='caldiroli'),
        )
    )
    assert weak.flags.get('caldiroli_below_flux_bound') is True
    assert weak.flags.get('efficiency_fallback_fixed') is True
    assert weak.diagnostics['hydrodynamic']['efficiency'] == pytest.approx(
        0.1, rel=1e-12, abs=0.0
    )


def test_t_exo_thermostat_mode_estimates_and_reports_itself():
    """The thermostat exobase mode estimates a temperature and records itself.

    On the bound CO2 case the estimator returns a temperature inside the
    thermostat bracket (above the equilibrium temperature, below the upper
    bracket edge), the diagnostics record which mode produced it, and the
    dispatch completes with a consistent per-species sum. The mode is a
    property of the call and not a warning about it, so it belongs in the
    diagnostics and not in the flags dictionary.
    """
    settings = DispatchSettings(T_exo_mode='thermostat')
    res = dispatch(
        _inputs(
            10 * Me, 1.8 * Re, 800.0, {'CO2': 1.0}, F_xuv=1.0, a=0.5 * AU, settings=settings
        )
    )
    assert res.diagnostics['hydrostatic']['T_exo_mode'] == 'thermostat'
    assert 'T_exo_thermostat' not in res.flags
    prescribed = dispatch(
        _inputs(10 * Me, 1.8 * Re, 800.0, {'CO2': 1.0}, F_xuv=1.0, a=0.5 * AU)
    )
    assert prescribed.diagnostics['hydrostatic']['T_exo_mode'] == 'prescribed'
    assert (
        prescribed.diagnostics['hydrostatic']['T_exo']
        != res.diagnostics['hydrostatic']['T_exo']
    )
    assert 800.0 <= res.diagnostics['hydrostatic']['T_exo'] <= 5.0e4
    if res.mdot > 0.0:
        assert sum(res.per_species.values()) == pytest.approx(res.mdot, rel=1e-6, abs=0.0)


def test_extend_mode_truncated_extension_keeps_the_clamp():
    """When the extension itself unbinds, the clamp stands, flagged.

    A loosely bound hot H2 envelope truncates its own upper structure
    before reaching the physical base pressure: the extend policy cannot
    place the base there, so the clamped level stands and the truncation
    is recorded alongside the clamp flags. The dispatch still returns a
    consistent result (totality).
    """
    settings = DispatchSettings(base_out_of_range='extend')
    inp = _inputs(
        1 * Me,
        1.0 * Re,
        2000.0,
        {'H2': 1.0},
        F_xuv=10.0,
        p_top=1e-6,  # requested; the bound structure truncates far deeper
        settings=settings,
    )
    res = dispatch(inp)
    assert res.flags.get('base_extension_truncated') is True
    assert res.flags.get('base_clamped') is True
    assert res.regime in REGIME_LABELS
    assert math.isfinite(res.mdot)


def test_stale_input_and_boreas_fallback_flags(monkeypatch):
    """Data-quality and base-method fallbacks surface as flags, not errors.

    An unconverged upstream atmosphere marks the result ``stale_input``
    without changing the contract; requesting the BOREAS base method with
    the dependency absent falls back to the Lopez base, flagged, and the
    negative-flux error contract still raises.
    """
    import sys

    monkeypatch.setitem(sys.modules, 'boreas', None)  # forces ImportError
    settings = DispatchSettings(base_method='boreas')
    inp = _inputs(
        5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, settings=settings, atm_converged=False
    )
    res = dispatch(inp)
    assert res.flags.get('stale_input') is True
    assert res.flags.get('base_method_fallback') == 'lopez'
    assert res.regime in REGIME_LABELS
    bad = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0)
    bad.F_xuv = -1.0
    with pytest.raises(ValueError, match='F_xuv'):
        dispatch(bad)


def test_settings_option_raises_cover_every_knob():
    """Every enumerated settings knob rejects an unknown value.

    The four option strings each raise with a message naming the knob, so
    a typo in a configuration surfaces at validation rather than as a
    silent default. The default settings validate silently.
    """
    with pytest.raises(ValueError, match='base_out_of_range'):
        DispatchSettings(base_out_of_range='nonsense').validate()
    with pytest.raises(ValueError, match='gate'):
        DispatchSettings(gate='nonsense').validate()
    with pytest.raises(ValueError, match='efficiency_mode'):
        DispatchSettings(efficiency_mode='nonsense').validate()
    with pytest.raises(ValueError, match='T_exo_mode'):
        DispatchSettings(T_exo_mode='nonsense').validate()
    DispatchSettings().validate()  # the defaults are a valid configuration


def test_flags_describe_the_branch_that_produced_the_rate():
    """A warning never survives onto a verdict its rate did not come from.

    The flags dictionary is read as a warning set about the returned result,
    so a caution about the wind temperature or the sonic radius must not ride
    along on a bolometric or hydrostatic verdict, whose rate those quantities
    did not set. The hydrodynamic candidates are always computed, because the
    diagnostics report them at every dispatch, which is what makes the
    scoping necessary rather than automatic.
    """
    hydro = dispatch(
        _inputs(Me, Re, 1000.0, {'N2': 0.8, 'O2': 0.2}, F_xuv=100.0, a=0.0775 * AU)
    )
    assert hydro.regime.startswith('hydrodynamic')
    assert hydro.flags.get('subcritical_sonic') is True
    boiloff = dispatch(
        _inputs(Me, 1.5 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=10.0, a=0.0775 * AU)
    )
    assert boiloff.regime == 'boiloff'
    # The same chain runs on this state and still reports subcritical in the
    # diagnostics, but the dispatched rate is the bolometric one.
    assert boiloff.diagnostics['hydrodynamic']['rr_chain']['R_s'] > 0.0
    for leaked in (
        'subcritical_sonic',
        'thermostat_clamped',
        'efficiency_fallback_fixed',
        'hydrostatic_lower_limit',
    ):
        assert leaked not in boiloff.flags, leaked
    # The hydrostatic lower-limit caution appears on, and only on, a
    # hydrostatic verdict.
    static = dispatch(_inputs(Me, Re, 1000.0, {'CO2': 1.0}, F_xuv=0.1, a=0.0775 * AU))
    assert static.regime == 'hydrostatic'
    assert static.flags.get('hydrostatic_lower_limit') is True
    assert 'subcritical_sonic' not in static.flags


def test_one_exobase_temperature_per_call():
    """The base extension and the branch stand on one upper structure.

    Under the extend policy the wind base is re-evaluated on a Bates
    extension, and the hydrostatic branch stands on one too. Both must be
    built at the same exobase temperature: resolving it twice, once from the
    settings and once from the equilibrium temperature, gave one call two
    thermospheres, and under the thermostat they differed by an order of
    magnitude in temperature, which sets the base density the wind rate is
    built from.
    """
    m_p, r_p, t_eq = 10 * Me, 1.8 * Re, 800.0
    # A profile whose top lies above the Lopez base, so the policy engages.
    settings = DispatchSettings(T_exo_mode='thermostat', base_out_of_range='extend')
    res = dispatch(
        _inputs(
            m_p,
            r_p,
            t_eq,
            {'CO2': 1.0},
            F_xuv=10.0,
            a=0.5 * AU,
            p_top=1.0,
            settings=settings,
        )
    )
    prof = isothermal_profile(m_p, r_p, t_eq, {'CO2': 1.0}, 1e7, 1.0)
    assert res.flags.get('base_extended') is True
    t_branch = res.diagnostics['hydrostatic']['T_exo']
    t_base = res.diagnostics['base_level']['T_K']
    t_top = float(prof.T[-1])
    # The thermostat must have moved off the equilibrium temperature, or the
    # test cannot tell the two resolutions apart.
    assert t_branch > 5.0 * t_eq
    # The base sits on the same extension: between the anchor and the exobase,
    # and nowhere near the equilibrium temperature the old path used.
    assert t_top < t_base <= t_branch
    assert t_base > 0.5 * t_branch
    assert t_base != pytest.approx(t_eq, rel=0.1, abs=0.0)


# Each row is (flag, a state that must raise it, a state that must not). Every
# warning the result can carry belongs here: a flag nothing asserts can be
# deleted without the suite noticing, which makes it a comment rather than part
# of the contract. The negative state is what stops a flag that fires on
# everything from passing as discrimination.
FLAG_CASES = (
    (
        'near_roche',
        dict(
            M_p=3 * Me, R_p=2 * Re, T_eq=1000.0, comp={'H2': 0.9, 'He': 0.1}, F_xuv=0.1, a=0.10
        ),
        dict(
            M_p=3 * Me, R_p=2 * Re, T_eq=1000.0, comp={'H2': 0.9, 'He': 0.1}, F_xuv=0.1, a=0.30
        ),
    ),
    (
        'k_tide_undefined',
        dict(M_p=Me, R_p=2 * Re, T_eq=2000.0, comp={'H2': 0.9, 'He': 0.1}, F_xuv=1.0, a=0.004),
        dict(M_p=Me, R_p=2 * Re, T_eq=2000.0, comp={'H2': 0.9, 'He': 0.1}, F_xuv=1.0, a=0.10),
    ),
    (
        'roche_overflow',
        dict(
            M_p=3 * Me,
            R_p=2 * Re,
            T_eq=1000.0,
            comp={'H2': 0.9, 'He': 0.1},
            F_xuv=0.1,
            a=0.0775,
        ),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
    (
        'bolometric_residual',
        dict(
            M_p=3 * Me,
            R_p=2 * Re,
            T_eq=1000.0,
            comp={'H2': 0.9, 'He': 0.1},
            F_xuv=0.1,
            a=0.0775,
        ),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
    (
        'thermostat_clamped',
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'N2': 0.8, 'O2': 0.2}, F_xuv=700.0, a=0.0775),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
    (
        'caldiroli_out_of_box',
        dict(
            M_p=Me,
            R_p=Re,
            T_eq=1000.0,
            comp={'CO2': 1.0},
            F_xuv=10.0,
            a=0.0775,
            settings=DispatchSettings(efficiency_mode='caldiroli'),
        ),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
    (
        'volkov_extrapolated',
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=0.1, a=0.0775),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
    (
        'contested_ion',
        dict(
            M_p=0.107 * Me,
            R_p=0.53 * Re,
            T_eq=440.0,
            comp={'CO2': 0.99, 'H2': 0.01},
            F_xuv=0.01,
            a=0.2,
            settings=DispatchSettings(T_exo_value=6000.0),
        ),
        dict(
            M_p=0.107 * Me,
            R_p=0.53 * Re,
            T_eq=440.0,
            comp={'CO2': 0.99, 'H2': 0.01},
            F_xuv=0.01,
            a=0.2,
        ),
    ),
    (
        'extension_unbound',
        dict(
            M_p=0.107 * Me,
            R_p=0.53 * Re,
            T_eq=440.0,
            comp={'CO2': 0.99, 'H2': 0.01},
            F_xuv=0.01,
            a=0.2,
            settings=DispatchSettings(T_exo_value=6000.0),
        ),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=0.1, a=0.0775),
    ),
    (
        'gate_rerouted',
        dict(
            M_p=0.107 * Me,
            R_p=0.53 * Re,
            T_eq=440.0,
            comp={'CO2': 0.99, 'H2': 0.01},
            F_xuv=0.01,
            a=0.2,
            settings=DispatchSettings(T_exo_value=15000.0),
        ),
        dict(
            M_p=0.107 * Me,
            R_p=0.53 * Re,
            T_eq=440.0,
            comp={'CO2': 0.99, 'H2': 0.01},
            F_xuv=0.01,
            a=0.2,
            settings=DispatchSettings(T_exo_value=3000.0),
        ),
    ),
    (
        'base_clamp_decades',
        dict(
            M_p=10 * Me,
            R_p=1.8 * Re,
            T_eq=800.0,
            comp={'CO2': 1.0},
            F_xuv=1.0,
            a=0.5,
            p_top=1.0,
        ),
        dict(M_p=Me, R_p=Re, T_eq=1000.0, comp={'CO2': 1.0}, F_xuv=10.0, a=0.0775),
    ),
)


@pytest.mark.physics_invariant
def test_every_warning_flag_has_a_state_that_raises_it_and_one_that_does_not():
    """Each flag fires on a state that warrants it and stays off otherwise.

    A flag no test asserts is a comment: it can be deleted and the suite stays
    green, so nothing holds the module to raising it. Each row pins one flag
    against a state that must raise it and a state that must not, and the
    second half is what keeps a flag that fires on everything from passing as
    a working warning.
    """
    for flag, on, off in FLAG_CASES:
        res_on = dispatch(_inputs(**_flag_state(on)))
        assert flag in res_on.flags, f'{flag} did not fire on its positive state'
        res_off = dispatch(_inputs(**_flag_state(off)))
        assert flag not in res_off.flags, f'{flag} fired on its negative state'


def _flag_state(spec):
    """Expand a FLAG_CASES row into _inputs keyword arguments."""
    kw = dict(spec)
    comp = kw.pop('comp')
    m_p, r_p, t_eq = kw.pop('M_p'), kw.pop('R_p'), kw.pop('T_eq')
    kw['a'] = kw['a'] * AU
    return dict(M_p=m_p, R_p=r_p, T_eq=t_eq, comp=comp, **kw)


@pytest.mark.physics_invariant
@pytest.mark.physics_invariant
def test_launch_level_is_a_width_on_a_realistic_column():
    """Off an isothermal column the launch level is a width, not a cancellation.

    The transfer rate is invariant to the launch level only along a column
    whose sound speed matches the density that column carries. Real
    profiles are not isothermal, so the cancellation leaves a residual, and
    the residual is reported rather than assumed small. The isothermal
    control is the contrast: the same sweep on the builder every other test
    uses moves the rate by a few percent, while a modest inversion moves it
    by a factor.
    """
    comp = {'H2': 0.9, 'He': 0.1}
    inverted = _inverted_profile(3 * Me, 2.2 * Re, 1000.0, 1600.0, comp)
    # The profile is warmer than the equilibrium temperature everywhere
    # above its base, so a rate reading T_eq instead of the profile is
    # detectable here and is not on an isothermal fixture.
    assert inverted.T[-1] > inverted.T[0]
    spreads = {}
    for label, profile in (('inverted', inverted), ('isothermal', None)):
        rates = []
        for p_photo in (2.0e4, 2.0e3, 2.0e2, 2.0e1):
            extra = {'profile': profile} if profile is not None else {}
            state = _inputs(
                3 * Me,
                2.2 * Re,
                1000.0,
                comp,
                F_xuv=13.4,
                a=0.07 * AU,
                settings=DispatchSettings(P_photo=p_photo),
                **extra,
            )
            rates.append(dispatch(state).diagnostics['nozzle']['rate_full_orbit_kg_s'])
        spreads[label] = max(rates) / min(rates)
    assert spreads['isothermal'] < 1.2
    assert spreads['inverted'] > 2.0
    noz = dispatch(
        _inputs(3 * Me, 2.2 * Re, 1000.0, comp, F_xuv=13.4, a=0.07 * AU, profile=inverted)
    ).diagnostics['nozzle']
    assert noz['T_K'] > 1000.0  # the profile's photospheric level, not T_eq


def test_nozzle_win_keeps_the_losing_branch_flow_radius():
    """A nozzle win reports the wind's flow radius, not the lobe radius.

    The lobe radius is a fixed fraction of the Hill radius at every
    planetary mass ratio, so substituting it would pin the reported ratio
    at a constant and discard the one geometric fact the screen still
    carries on this branch. The lobe radius has its own key.
    """
    res = dispatch(
        _inputs(3 * Me, 2.2 * Re, 1000.0, {'H2': 0.9, 'He': 0.1}, F_xuv=13.4, a=0.07 * AU)
    )
    roche, noz = res.diagnostics['roche'], res.diagnostics['nozzle']
    assert roche['rate_branch'] == 'roche_overflow'
    assert roche['flow_radius'] != pytest.approx(noz['r_lobe'], rel=1e-6)
    assert roche['xi_flow'] == pytest.approx(
        roche['R_hill_periapsis'] / roche['flow_radius'], rel=1e-12
    )
    # The constant the substitution would have produced, for contrast.
    assert roche['R_hill_periapsis'] / noz['r_lobe'] == pytest.approx(1.418, rel=0.01)


def test_nozzle_win_splits_by_reservoir_not_by_the_losing_branch():
    """A nozzle win carries a bulk split, not the split its rival computed.

    The transfer is a bulk flow through L1 with no per-species physics, so
    the elements leave in their reservoir proportions. The state below is
    the discriminator: its hydrostatic rival gives hydrogen the entire
    flux by Jeans selection, so inheriting that split instead of replacing
    it would be visible as hydrogen dominating a carbon dioxide planet.
    """
    settings = DispatchSettings(T_exo_value=2000.0)
    won = dispatch(
        _inputs(
            0.107 * Me,
            0.53 * Re,
            600.0,
            {'CO2': 0.99, 'H2': 0.01},
            F_xuv=1e-4,
            a=0.008 * AU,
            settings=settings,
        )
    )
    assert won.diagnostics['roche']['rate_branch'] == 'roche_overflow'
    total = sum(won.per_species.values())
    assert total == pytest.approx(won.mdot, rel=1e-12)
    shares = {k: v / total for k, v in won.per_species.items()}
    assert shares['O'] > shares['C'] > shares['H']
    assert shares['H'] < 1e-2
    # The same composition on the branch the nozzle displaced, for contrast.
    rival = dispatch(
        _inputs(
            0.107 * Me,
            0.53 * Re,
            600.0,
            {'CO2': 0.99, 'H2': 0.01},
            F_xuv=1e-4,
            a=0.2 * AU,
            settings=settings,
        )
    )
    assert rival.regime == 'hydrostatic'
    rival_total = sum(rival.per_species.values())
    assert rival.per_species['H'] / rival_total > 0.99


def test_dispatched_split_names_the_element_that_leaves():
    """The dispatched split is checked by identity, not only by its sum.

    The dispatcher renormalizes the per-species rates onto the bulk rate, so
    a sums-to-mdot assertion cannot fail however the shares are assigned: a
    permuted mapping conserves total mass while moving the wrong elements out
    of the planet, which is what the PROTEUS side debits reservoirs by. Both
    branches that produce a split are pinned by which element dominates.
    """
    # Hydrostatic: only hydrogen is light enough to leave the Mars-mass host,
    # and it carries the rate by nineteen decades over the heavy background.
    static = dispatch(
        _inputs(0.107 * Me, 0.53 * Re, 440.0, {'CO2': 0.99, 'H2': 0.01}, F_xuv=0.01, a=0.2 * AU)
    )
    assert static.regime == 'hydrostatic'
    assert set(static.per_species) == {'H', 'C', 'O'}
    assert max(static.per_species, key=static.per_species.get) == 'H'
    assert static.per_species['H'] > 1.0e19 * static.per_species['O']
    assert static.per_species['O'] > static.per_species['C']  # two O per C, heavier
    # Hydrodynamic: a carbon dioxide wind carries oxygen over carbon, in the
    # ratio the closure returns rather than one this test recomputes, but the
    # ordering and the absence of hydrogen are its own statement.
    wind = dispatch(_inputs(Me, Re, 1000.0, {'CO2': 1.0}, F_xuv=5.0e3, a=0.0775 * AU))
    assert wind.regime.startswith('hydrodynamic')
    assert set(wind.per_species) == {'C', 'O'}
    assert wind.per_species['O'] > wind.per_species['C']
    assert 1.5 < wind.per_species['O'] / wind.per_species['C'] < 4.0
    # And the sum still holds, which is the weaker claim of the two.
    for res in (static, wind):
        assert sum(res.per_species.values()) == pytest.approx(res.mdot, rel=1e-9, abs=0.0)


@pytest.mark.reference_pinned
def test_caldiroli_efficiency_geometry_conversion():
    """The fitted efficiency is converted to the geometry it is used in.

    Caldiroli et al. (2022) fit their efficiency against a rate written on an
    ``R_p^3`` geometry, while the dispatcher's energy-limited rate is the
    Erkaev form on ``R_p R_XUV^2`` (``scaling=2``). Decision 12 therefore
    converts the fitted value by ``(R_p / R_XUV)^2`` before using it. The
    conversion is a pure geometric factor, so nothing about the rate's shape
    reveals whether it was applied: dropping it entirely left the suite green.
    This pins it against the two radii the same call reports.
    """
    settings = DispatchSettings(efficiency_mode='caldiroli')
    inp = _inputs(Me, Re, 1000.0, {'CO2': 1.0}, F_xuv=10.0, a=0.0775 * AU, settings=settings)
    res = dispatch(inp)
    hy = res.diagnostics['hydrodynamic']
    raw, _flags = caldiroli_efficiency(10.0, Me, Re, hy['K_tide'])
    assert raw is not None
    # The XUV radius is the photospheric level the settings select, which is
    # the radius the Erkaev form cubes; it is recomputed here from the same
    # profile rather than read back, so the test does not depend on the
    # module reporting it.
    photo, _pf = photospheric_level(inp.profile, settings.P_photo)
    r_xuv = photo['r']
    factor = (Re / r_xuv) ** 2
    assert hy['efficiency'] == pytest.approx(raw * factor, rel=1e-9, abs=0.0)
    # Discrimination: the factor is not 1 on this state, so an unconverted
    # efficiency is a different number.
    assert factor != pytest.approx(1.0, rel=1e-3, abs=0.0)
    assert hy['efficiency'] != pytest.approx(raw, rel=1e-3, abs=0.0)
    # The XUV radius is above the planetary radius, so the conversion always
    # reduces the efficiency: the fitted value refers to a larger absorbing
    # area than the Erkaev geometry charges for.
    assert r_xuv > Re
    assert hy['efficiency'] < raw
