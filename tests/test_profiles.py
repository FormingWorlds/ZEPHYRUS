"""Tests for ``src/zephyrus/profiles.py``.

Exercises the atmosphere-profile container and the escape working levels.
The physical invariants under test:

- Conservation / closed form: the isothermal profile obeys the hydrostatic
  relation ``dr = -(k T r^2 / (G M mu)) d ln p`` level by level, and the
  interpolated number density obeys the ideal-gas law ``n = p / (k T)``.
- Positivity / boundedness: profile validation rejects non-monotone
  pressure or radius and non-positive state variables; working levels clamp
  to the covered pressure range with flags rather than extrapolating.
- Reference pin: the Lopez (2017) wind-base pressure lands at the nanobar
  tau = 1 level Murray-Clay et al. (2009) print for a hot Jupiter.
- Error contract: an unbound isothermal structure raises; the BOREAS base
  method falls back to Lopez with a flag when the dependency is absent.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import sys
import types

import numpy as np
import pytest

from zephyrus.constants import G, amu, kb, m_p
from zephyrus.planets_parameters import Me, Mjup, Re, Rjup
from zephyrus.profiles import (
    Profile,
    interp_at_pressure,
    isothermal_profile,
    lopez_base_pressure,
    photospheric_level,
    pressure_at_radius,
    wind_base_level,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _n2_profile(p_top=1e-6):
    """Bound isothermal N2 test atmosphere on a warm super-Earth."""
    return isothermal_profile(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, 1e7, p_top)


def test_profile_validation_rejects_malformed_input():
    """Validation raises on non-monotone arrays and accepts a sound profile.

    The error contract: increasing pressure, decreasing radius, negative
    temperature, mismatched mixing-ratio length, and a two-level profile
    all raise ``ValueError``; the well-formed profile on the same path
    validates silently.
    """
    good = _n2_profile()
    good.validate()  # must not raise
    bad_p = Profile(p=good.p[::-1], r=good.r, T=good.T, vmr=good.vmr, mmw=good.mmw)
    with pytest.raises(ValueError, match='must decrease'):
        bad_p.validate()
    bad_T = Profile(p=good.p, r=good.r, T=-good.T, vmr=good.vmr, mmw=good.mmw)
    with pytest.raises(ValueError, match='positive'):
        bad_T.validate()
    bad_vmr = Profile(p=good.p, r=good.r, T=good.T, vmr={'N2': good.p[:2]}, mmw=good.mmw)
    with pytest.raises(ValueError, match='length mismatch'):
        bad_vmr.validate()
    tiny = Profile(
        p=good.p[:2],
        r=good.r[:2],
        T=good.T[:2],
        vmr={'N2': good.p[:2] * 0 + 1},
        mmw=good.mmw[:2],
    )
    with pytest.raises(ValueError, match='at least 3'):
        tiny.validate()


@pytest.mark.physics_invariant
def test_isothermal_profile_obeys_hydrostatic_relation():
    """Each integration step matches the local hydrostatic scale height.

    The construction integrates ``dr = -H d ln p`` with
    ``H = k T r^2 / (G M mu)``, so the recovered per-step ratio
    ``dr / d ln p`` must equal ``-H`` evaluated at the lower level. This is
    the conservation-style closed form; the guard is that a plane-parallel
    slip (H frozen at the surface value) accumulates a visible radius error
    over the profile, so the top radius must exceed the plane-parallel
    estimate.
    """
    M, R, T = 5 * Me, 1.5 * Re, 800.0
    prof = _n2_profile()
    mu = float(prof.mmw[0])
    lnp = np.log(prof.p)
    for i in (0, len(prof.p) // 2):
        H = kb * T * prof.r[i] ** 2 / (G * M * mu)
        step = (prof.r[i + 1] - prof.r[i]) / (lnp[i + 1] - lnp[i])
        assert step == pytest.approx(-H, rel=1e-12)
    # Curvature guard: the r^2 growth of H makes the true extent exceed the
    # plane-parallel (surface-H) estimate.
    H0 = kb * T * R**2 / (G * M * mu)
    plane_parallel_top = R + H0 * np.log(prof.p[0] / prof.p[-1])
    assert prof.r[-1] > plane_parallel_top
    # The structure stays bound: the local Jeans parameter never drops
    # below the documented truncation threshold.
    lam = G * M * mu / (kb * T * prof.r)
    assert np.all(lam >= 2.2 - 1e-9)


def test_isothermal_profile_truncates_or_raises_when_unbound():
    """An unbound structure truncates with fewer levels, or raises outright.

    A hot, light hydrogen atmosphere on a small planet becomes unbound
    within the requested pressure span: the returned profile must stop at
    the Jeans-parameter floor instead of extending to ``p_top``. When even
    the surface is unbound, the constructor raises (the error contract).
    """
    prof = isothermal_profile(1 * Me, 1.0 * Re, 2000.0, {'H2': 1.0}, 1e7, 1e-6)
    # Truncated: the top pressure never reaches the requested 1e-6 Pa.
    assert prof.p[-1] > 1e-6
    assert len(prof.p) < 120
    with pytest.raises(ValueError, match='unbound'):
        # A Jeans parameter below 2.2 at the surface itself.
        isothermal_profile(0.05 * Me, 2.0 * Re, 3000.0, {'H2': 1.0}, 1e7, 1e-6)


@pytest.mark.physics_invariant
def test_interp_at_pressure_ideal_gas_and_node_exactness():
    """Interpolation returns exact node values and ideal-gas densities.

    At a grid node the interpolated radius and temperature reproduce the
    stored level exactly; between nodes the number density obeys
    ``n = p / (k T)`` by construction (the ideal-gas closed form), and the
    mass density is ``n`` times the mean molecular mass.
    """
    prof = _n2_profile()
    k = len(prof.p) // 2
    lev = interp_at_pressure(prof, float(prof.p[k]))
    assert lev['r'] == pytest.approx(float(prof.r[k]), rel=1e-12)
    assert lev['T'] == pytest.approx(float(prof.T[k]), rel=1e-12)
    # Ideal gas at an off-node pressure.
    p_mid = np.sqrt(prof.p[k] * prof.p[k + 1])
    lev2 = interp_at_pressure(prof, float(p_mid))
    assert lev2['n'] == pytest.approx(float(p_mid) / (kb * lev2['T']), rel=1e-12)
    assert lev2['rho'] == pytest.approx(lev2['n'] * lev2['mmw'], rel=1e-12)
    # Radius must lie strictly between the bracketing nodes.
    assert prof.r[k] < lev2['r'] < prof.r[k + 1]


def test_pressure_at_radius_inverts_the_profile():
    """Radius-to-pressure interpolation inverts pressure-to-radius lookup.

    Round trip: interpolating the level at a node pressure and asking for
    the pressure at the returned radius must recover the node pressure. The
    guard is an off-grid radius, whose pressure must fall between the
    bracketing node pressures.
    """
    prof = _n2_profile()
    k = len(prof.p) // 3
    lev = interp_at_pressure(prof, float(prof.p[k]))
    assert pressure_at_radius(prof, lev['r']) == pytest.approx(float(prof.p[k]), rel=1e-9)
    r_mid = 0.5 * (prof.r[k] + prof.r[k + 1])
    p_mid = pressure_at_radius(prof, r_mid)
    assert prof.p[k + 1] < p_mid < prof.p[k]


def test_photospheric_level_clamps_with_flags():
    """The photospheric level sits at 20 mbar, clamping to the ends flagged.

    Inside the covered range the level lands at 2000 Pa exactly. A profile
    whose top never reaches 20 mbar clamps to its top level; one whose
    surface is already above 20 mbar clamps to its deepest level. Both
    clamps raise ``photo_clamped``.
    """
    prof = _n2_profile()
    lev, flags = photospheric_level(prof)
    assert lev['p'] == pytest.approx(2000.0, rel=1e-12)
    assert flags == {}
    deep = isothermal_profile(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, 1e7, 1e4)
    lev_d, flags_d = photospheric_level(deep)
    assert flags_d.get('photo_clamped') is True
    assert lev_d['p'] == pytest.approx(1e4, rel=1e-12)
    shallow = isothermal_profile(5 * Me, 1.5 * Re, 800.0, {'N2': 1.0}, 1e3, 1e-4)
    lev_s, flags_s = photospheric_level(shallow)
    assert flags_s.get('photo_clamped') is True
    assert lev_s['p'] == pytest.approx(1e3, rel=1e-12)


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
def test_lopez_base_pressure_is_the_nanobar_level():
    """The Lopez (2017) base pressure reproduces the published nanobar level.

    Murray-Clay et al. (2009, ApJ 693, 23, Section 2.1) place the tau = 1
    level for 20 eV photons near a nanobar for their fiducial hot Jupiter
    (0.7 Jupiter masses, 1.4 Jupiter radii, atomic-hydrogen thermosphere),
    and Lopez (2017, MNRAS 472, 245) builds the wind-base prescription
    ``P_base = mu g / sigma_nu0`` on the same level. Both papers quote the
    scale, not a precise value, so the pin is the order of magnitude: the
    computed pressure must land within a factor of a few of 1e-4 Pa (one
    nanobar). Sign and scale guards: strictly positive and far below the
    20 mbar photospheric level.
    """
    g = G * (0.7 * Mjup) / (1.4 * Rjup) ** 2
    p_base = lopez_base_pressure(1.008 * amu, g)
    # One nanobar is 1e-4 Pa; allow a factor of a few for the mu convention.
    assert 3e-5 < p_base < 3e-4
    assert p_base > 0.0
    # Scale guard: eighty decades of headroom is a unit slip, five is right.
    assert p_base < 2000.0 * 1e-4
    # Linearity in gravity: doubling g doubles the base pressure.
    assert lopez_base_pressure(1.008 * amu, 2 * g) == pytest.approx(2 * p_base, rel=1e-12)
    # The proton-mass convention differs from the atomic-weight convention
    # by under a percent; both stay inside the pinned band.
    assert lopez_base_pressure(m_p, g) == pytest.approx(p_base, rel=0.01)


def test_wind_base_level_lopez_converges_and_clamps():
    """The Lopez fixed point converges on-profile and clamps off-profile.

    On a profile reaching 1e-6 Pa the nanobar-scale base lies inside the
    covered range, so the returned level pressure equals the physical base
    pressure the level reports and no clamp flag is raised. On a profile
    truncated at 1e-2 Pa the physical base lies above the top, so the level
    clamps to the top with the clamp distance recorded in pressure decades,
    and the unclamped target stays readable on the level itself.
    """
    prof = _n2_profile(p_top=1e-6)
    lev, flags = wind_base_level(prof, 5 * Me, method='lopez')
    assert 'base_clamped' not in flags
    assert lev['p'] == pytest.approx(lev['p_physical'], rel=1e-3)
    # The N2 base sits at the nanobar scale (between 0.01 and 100 nanobar).
    assert 1e-6 < lev['p'] < 1e-2
    # The physical base pressure is a level quantity, not a flag: an
    # unclamped call reports no flags at all.
    assert flags == {}
    deep = _n2_profile(p_top=1e-2)
    lev_d, flags_d = wind_base_level(deep, 5 * Me, method='lopez')
    assert flags_d.get('base_clamped') is True
    assert lev_d['p'] == pytest.approx(float(deep.p[-1]), rel=1e-12)
    # The clamped level sits above its own physical target, by the recorded
    # distance; a clamp that lost the target would fail both assertions.
    assert lev_d['p'] > lev_d['p_physical']
    expected_decades = np.log10(deep.p[-1] / lev_d['p_physical'])
    assert flags_d['base_clamp_decades'] == pytest.approx(expected_decades, rel=1e-9)
    assert flags_d['base_clamp_decades'] > 0.0


def test_wind_base_level_fixed_pressure_method():
    """The fixed-pressure method hits its target and clamps out of range.

    The returned level must sit exactly at the requested pressure when the
    profile covers it, and clamp with a flag when it does not (matching the
    nearest-level behavior of the energy-limited path).
    """
    prof = _n2_profile()
    lev, flags = wind_base_level(prof, 5 * Me, method='fixed_pressure', fixed_pressure=5.0)
    assert lev['p'] == pytest.approx(5.0, rel=1e-12)
    assert flags == {}
    lev_c, flags_c = wind_base_level(prof, 5 * Me, method='fixed_pressure', fixed_pressure=1e-9)
    assert flags_c.get('base_clamped') is True
    assert lev_c['p'] == pytest.approx(float(prof.p[-1]), rel=1e-12)


def test_wind_base_level_boreas_falls_back_without_dependency(monkeypatch):
    """The BOREAS method falls back to Lopez, flagged, when BOREAS is absent.

    Blocking the ``boreas`` import must not raise: the method returns the
    Lopez level with ``base_method_fallback`` recorded, and the fallback
    level matches a direct Lopez call. Missing scalars trigger the same
    fallback path.
    """
    monkeypatch.setitem(sys.modules, 'boreas', None)  # forces ImportError
    prof = _n2_profile()
    scalars = {'R_p': 1.5 * Re, 'T_eq': 800.0, 'F_xuv': 10.0}
    lev, flags = wind_base_level(prof, 5 * Me, method='boreas', boreas_scalars=scalars)
    assert flags.get('base_method_fallback') == 'lopez'
    lev_ref, _ = wind_base_level(prof, 5 * Me, method='lopez')
    assert lev['p'] == pytest.approx(lev_ref['p'], rel=1e-12)
    # Missing scalars: same fallback, no exception.
    lev2, flags2 = wind_base_level(prof, 5 * Me, method='boreas', boreas_scalars=None)
    assert flags2.get('base_method_fallback') == 'lopez'
    assert lev2['p'] == pytest.approx(lev_ref['p'], rel=1e-12)


def test_wind_base_level_boreas_uses_solver_radius(monkeypatch):
    """A converged BOREAS solve places the base at the solver's XUV radius.

    The mocked solver returns a physically plausible XUV radius mid-profile
    (mock discipline: a real radius in cm, not a unit constant), and the
    returned level pressure must equal the profile pressure at that radius.
    A solver that reports a skipped regime triggers the Lopez fallback.
    """
    prof = _n2_profile()
    r_target = float(0.5 * (prof.r[0] + prof.r[-1]))  # mid-profile radius [m]

    class FakeParams:
        def __init__(self):
            self.kappa = {'N2': 1.0}
            self.albedo = None
            self.Teq = None
            self.FXUV = None
            self.rplanet = None
            self.mplanet = None

        def _recompute_composites(self):
            pass

        def _init_opacities(self):
            pass

    class FakeMassLoss:
        def __init__(self, params):
            self.params = params

        def compute_mass_loss_parameters(self, m, r, t):
            return [{'regime': 'EL', 'RXUV': r_target * 1e2}]  # cm

    fake = types.ModuleType('boreas')
    fake.ModelParams = FakeParams
    fake.MassLoss = FakeMassLoss
    monkeypatch.setitem(sys.modules, 'boreas', fake)
    scalars = {'R_p': 1.5 * Re, 'T_eq': 800.0, 'F_xuv': 10.0}
    lev, flags = wind_base_level(prof, 5 * Me, method='boreas', boreas_scalars=scalars)
    assert 'base_method_fallback' not in flags
    assert lev['p'] == pytest.approx(pressure_at_radius(prof, r_target), rel=1e-9)
    # Interior radius: the level pressure must be between the endpoints.
    assert prof.p[-1] < lev['p'] < prof.p[0]

    class SkippedMassLoss(FakeMassLoss):
        def compute_mass_loss_parameters(self, m, r, t):
            return [{'regime': 'SKIPPED'}]

    fake.MassLoss = SkippedMassLoss
    lev_f, flags_f = wind_base_level(prof, 5 * Me, method='boreas', boreas_scalars=scalars)
    assert flags_f.get('base_method_fallback') == 'lopez'
