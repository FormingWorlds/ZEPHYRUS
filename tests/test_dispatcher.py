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

from zephyrus.dispatcher import (
    REGIME_LABELS,
    DispatchSettings,
    EscapeInputs,
    dispatch,
)
from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Ms, Re
from zephyrus.profiles import isothermal_profile

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
            assert tot == pytest.approx(res.mdot, rel=1e-6)
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
    assert hydro['mdot_el'] == pytest.approx(ref, rel=1e-9)
    assert res.mdot == pytest.approx(min(hydro['mdot_el'], hydro['mdot_rr']), rel=1e-9)
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


def test_diagnostics_are_boxed(monkeypatch):
    """Sabotaging every diagnostics producer changes no dispatch outcome.

    The container is reporting only: no control flow reads it back. The
    test proves it by replacing every diagnostics-side producer with a
    stub returning garbage of the right shape and asserting the regime,
    the bulk rate, and the per-species split are unchanged; a regression
    in which any diagnostic feeds a routing decision fails loudly here.
    """
    inp = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0)
    reference = dispatch(inp)
    nan = float('nan')
    monkeypatch.setattr('zephyrus.diagnostics.q_net_over_qc', lambda *a, **k: (nan, nan, nan))
    monkeypatch.setattr('zephyrus.diagnostics.guo_triple', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.erkaev_tc', lambda *a, **k: nan)
    monkeypatch.setattr('zephyrus.diagnostics.potential_screens', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.along_profile_fluid_check', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.diagnostics.self_consistency_screen', lambda *a, **k: {})
    monkeypatch.setattr('zephyrus.boiloff.tang_timescale_check', lambda *a, **k: {})
    sabotaged = dispatch(inp)
    assert sabotaged.regime == reference.regime
    assert sabotaged.mdot == pytest.approx(reference.mdot, rel=1e-12)
    for el, v in reference.per_species.items():
        assert sabotaged.per_species[el] == pytest.approx(v, rel=1e-12)
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
    placing the base at a lower pressure than the profile top.
    """
    inp_clamp = _inputs(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, p_top=1e-2)
    res_clamp = dispatch(inp_clamp)
    assert res_clamp.flags.get('base_clamped') is True
    assert res_clamp.flags['base_clamp_decades'] > 0.0
    settings = DispatchSettings(base_out_of_range='extend')
    inp_ext = _inputs(
        5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, F_xuv=1.0, p_top=1e-2, settings=settings
    )
    res_ext = dispatch(inp_ext)
    assert res_ext.flags.get('base_extended') is True
    assert 'base_clamped' not in res_ext.flags
    assert res_ext.diagnostics['base_level']['p_Pa'] < 1e-2


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
    assert sum(on.per_species.values()) == pytest.approx(on.mdot, rel=1e-6)
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
    assert sum(off.per_species.values()) == pytest.approx(off.mdot, rel=1e-6)


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
    assert eff != pytest.approx(0.1, rel=0.05)
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
    assert weak.diagnostics['hydrodynamic']['efficiency'] == pytest.approx(0.1, rel=1e-12)


def test_t_exo_thermostat_mode_estimates_and_flags():
    """The thermostat exobase mode estimates a temperature and flags itself.

    On the bound CO2 case the estimator returns a temperature inside the
    thermostat bracket (above the equilibrium temperature, below the upper
    bracket edge), the flag records the mode, and the dispatch completes
    with a consistent per-species sum.
    """
    settings = DispatchSettings(T_exo_mode='thermostat')
    res = dispatch(
        _inputs(
            10 * Me, 1.8 * Re, 800.0, {'CO2': 1.0}, F_xuv=1.0, a=0.5 * AU, settings=settings
        )
    )
    assert res.flags.get('T_exo_thermostat') is True
    assert 800.0 <= res.diagnostics['hydrostatic']['T_exo'] <= 5.0e4
    if res.mdot > 0.0:
        assert sum(res.per_species.values()) == pytest.approx(res.mdot, rel=1e-6)


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
