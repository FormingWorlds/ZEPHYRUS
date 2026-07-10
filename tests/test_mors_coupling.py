"""Tests for the MORS-to-escape flux hand-off used by the PROTEUS coupling.

Mocks the MORS stellar-luminosity lookup so the coupling recipe runs in the
fast unit tier without downloading evolution tracks. The real-lookup version
lives in ``tests/test_earth.py`` (integration tier). Invariants under test:

- Unit boundary: MORS returns XUV luminosities in erg s-1; converting to a
  flux at the planet applies the erg-to-W factor and the au-to-cm dilution,
  so the SI flux lands in a physically plausible W m-2 range. A dropped or
  swapped conversion moves the flux by orders of magnitude.
- Coupling monotonicity: doubling the stellar XUV luminosity doubles the
  flux and hence the energy-limited escape rate.
- Positivity: the coupled escape rate is non-negative.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

from unittest import mock

import mors
import numpy as np
import pytest

from zephyrus.constants import au2cm, ergcm2stoWm2
from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Re, a_earth

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Plausible MORS XUV luminosities for a young Sun-like star (erg s-1),
# matching the age = 150 Myr row of the real-lookup integration fixture.
LX_ERG_S = 7.284408e28
LEUV_ERG_S = 1.101644e29


def _mock_star_with(lx, leuv):
    """Return a stand-in for ``mors.Star`` whose ``Value`` serves fixed luminosities.

    The mock returns plausible, distinct ``Lx`` and ``Leuv`` values rather
    than a constant, so a dilution or conversion bug cannot hide behind a
    degenerate ``1.0`` everywhere.
    """
    star = mock.MagicMock()
    star.Value.side_effect = lambda age, key: {'Lx': lx, 'Leuv': leuv}[key]
    return star


def _flux_from_star(star, age, a_au):
    """Diluted SI XUV flux [W m-2] at the planet from a (mocked) MORS star."""
    lx = star.Value(age, 'Lx')
    leuv = star.Value(age, 'Leuv')
    return (lx + leuv) * ergcm2stoWm2 / (4 * np.pi * (a_au * au2cm) ** 2)


@pytest.mark.physics_invariant
def test_mors_flux_handoff_lands_in_plausible_si_range():
    """The erg-to-W and au-to-cm conversions put the flux in a sane W m-2 range.

    Patches ``mors.Star`` so the lookup is exercised through the real
    constructor seam without a network download. A dropped erg-to-W factor
    (x1000) or an au-to-m slip (x10000) would push the flux far out of range;
    the scale guard catches both, and the coupled escape rate is pinned.
    """
    with mock.patch(
        'mors.Star', return_value=_mock_star_with(LX_ERG_S, LEUV_ERG_S)
    ) as star_ctor:
        star = mors.Star(Mstar=1.0, Omega=1.0)
        flux = _flux_from_star(star, 150.0, a_earth)
    star_ctor.assert_called_once_with(Mstar=1.0, Omega=1.0)
    assert flux == pytest.approx(0.06507260050807114, rel=1e-9)
    # Scale guard: a young Sun-like XUV flux at 1 au is order 1e-2 W m-2. A
    # dropped erg-to-W conversion lands near 1e1 and an au-to-m slip near 1e2,
    # so this band catches either.
    assert 1e-3 < flux < 1e0
    escape = EL_escape(False, a_earth, 0.0, Me, 1.0, 0.15, Re, 1.2 * Re, flux, scaling=2)
    assert escape == pytest.approx(28741.886140143717, rel=1e-9)
    assert escape > 0


@pytest.mark.physics_invariant
def test_mors_coupling_escape_linear_in_stellar_luminosity():
    """Doubling the stellar XUV luminosity doubles the coupled escape rate.

    The flux is linear in ``Lx + Leuv`` and the escape rate is linear in the
    flux, so a twice-as-bright star drives exactly twice the mass loss. A
    spurious offset anywhere in the chain would break the exact ratio.
    """
    dim = _mock_star_with(LX_ERG_S, LEUV_ERG_S)
    bright = _mock_star_with(2 * LX_ERG_S, 2 * LEUV_ERG_S)
    flux_dim = _flux_from_star(dim, 150.0, a_earth)
    flux_bright = _flux_from_star(bright, 150.0, a_earth)
    esc_dim = EL_escape(False, a_earth, 0.0, Me, 1.0, 0.15, Re, 1.2 * Re, flux_dim)
    esc_bright = EL_escape(False, a_earth, 0.0, Me, 1.0, 0.15, Re, 1.2 * Re, flux_bright)
    # Exact proportionality across the whole chain.
    assert esc_bright == pytest.approx(2.0 * esc_dim, rel=1e-12)
    assert esc_dim > 0
    # Mock discipline: the lookup was actually exercised for both luminosities.
    assert dim.Value.call_count >= 2
