"""MORS-coupled escape regression for an Earth-analogue around a Sun-like star.

Cross-cutting integration test spanning the two seams that ``EL_escape``
depends on in a real PROTEUS run: the MORS stellar-XUV lookup and the
erg-to-W / au-to-cm flux conversion that feeds the energy-limited mass-loss
formula. Unlike the mocked unit coverage in ``tests/test_mors_coupling.py``,
this exercises the real MORS ``Spada`` tracks, so it carries the integration
tier and runs nightly.

Invariants exercised:

- Reference regression: at three stellar ages spanning the young, middle-aged
  and old Sun, the returned ``(Lx, Leuv, Fxuv, escape)`` tuple is pinned to
  hand-verified values, so a drift in the MORS coupling or the conversion
  chain fails loudly.
- Positivity: the diluted XUV flux and the escape rate are strictly positive
  for a physically valid Earth-analogue.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import mors
import numpy as np
import pytest

from zephyrus.constants import au2cm, au2m, ergcm2stoWm2
from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Re, a_earth, e_earth

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]

# Reference tuples: (age [Myr], (Lx [erg s-1], Leuv [erg s-1],
# Fxuv [W m-2], escape [kg s-1])) for an Earth-analogue on a 1 au orbit
# around a 1 Msun, Omega = 1 Omega_sun star on the MORS Spada tracks.
# The ages span the young, middle-aged and old Sun.
TEST_DATA = (
    (150.0, (7.284408e28, 1.101644e29, 0.06507259, 19959.640)),
    (5813.0, (8.851172e27, 3.046169e28, 0.01397853, 4287.620)),
    (10020.0, (7.109151e27, 3.345416e28, 0.01442316, 4423.998)),
)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'age_star,expected',
    TEST_DATA,
    ids=['young_150Myr', 'middle_5.8Gyr', 'old_10Gyr'],
)
def test_earth_analogue_mors_coupled_escape(age_star, expected):
    """Pin the MORS-coupled escape rate for an Earth-analogue at three ages.

    Loads the real MORS Spada tracks, reads the stellar XUV luminosities at
    the given age, dilutes them to a flux at 1 au, and drives the
    energy-limited escape formula. The full ``(Lx, Leuv, Fxuv, escape)``
    tuple is pinned to hand-verified reference values; the flux and escape
    rate are additionally checked to be strictly positive so a sign or
    conversion regression cannot hide behind the tuple comparison.
    """
    # Load the stellar-evolution tracks from MORS (real lookup, nightly tier).
    mors.DownloadEvolutionTracks('Spada')
    star = mors.Star(Mstar=1.0, Omega=1.0)

    lx = star.Value(age_star, 'Lx')
    leuv = star.Value(age_star, 'Leuv')
    # erg s-1 -> W m-2 at the planet: erg-to-W factor with the orbital
    # distance in centimetres so the 1/(4 pi a^2) dilution is consistent.
    fxuv_star_si = (lx + leuv) * ergcm2stoWm2 / (4 * np.pi * (a_earth * au2cm) ** 2)
    # Earth geometry: Rp = Rxuv = Re, so the scaling=2 default reduces to Re**3.
    escape = EL_escape(False, a_earth * au2m, e_earth, Me, 1.0, 0.15, Re, Re, fxuv_star_si)

    ret = (lx, leuv, fxuv_star_si, escape)

    # rtol=1e-5 because the reference tuple is pinned to six/seven sig figs.
    np.testing.assert_allclose(ret, expected, rtol=1e-5, atol=0)
    # Positivity guards: a young Sun-like star drives a strictly positive
    # flux and mass-loss rate; a sign flip in either would pass the tuple
    # comparison only if the reference were also wrong.
    assert fxuv_star_si > 0
    assert escape > 0
