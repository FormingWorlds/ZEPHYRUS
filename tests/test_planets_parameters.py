"""Tests for ``src/zephyrus/planets_parameters.py``.

``planets_parameters.py`` is a utility source (star-planet system constants),
so it is exempt from the physics-invariant requirement. A wrong parameter is a
silent physics error, so these tests pin the shipped values, bracket their
physical magnitude, and assert the derived-parameter relations: the TOI-561
star and planet quantities are defined as fixed multiples of the solar and
Earth values, and those multipliers are the discrimination guards. The
Earth mean density recovered from the pinned mass and radius is checked
against the known rocky-planet value.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import numpy as np
import pytest

from zephyrus.constants import au2m
from zephyrus.planets_parameters import (
    L_TOI561,
    M_TOI561,
    R_TOI561,
    Fxuv_earth_10Myr,
    Fxuv_earth_today,
    Ls,
    M_TOI561b,
    Me,
    Mjup,
    Ms,
    R_TOI561b,
    Re,
    Rjup,
    Rs,
    a_earth,
    e_earth,
    e_TOI561b,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_earth_mass_radius_recover_rocky_density():
    """Earth mass and radius reproduce a rocky-planet mean density.

    The pinned ``Me`` and ``Re`` are the discrimination guard for each other:
    the mean density ``Me / (4/3 pi Re^3)`` must land near the terrestrial
    ``5.5e3 kg m-3``. A decimal slip in either constant moves the density by
    an order of magnitude and fails the band, which a bare value pin would
    miss.
    """
    assert Me == pytest.approx(5.9722e24, rel=1e-9, abs=0.0)
    assert Re == pytest.approx(6.378e6, rel=1e-9, abs=0.0)
    rho = Me / (4 / 3 * np.pi * Re**3)
    # Terrestrial mean density, ~5.5e3 kg m-3; a decimal slip lands outside.
    assert 5.0e3 < rho < 6.0e3
    # Re is the equatorial radius, distinct from the ~6.357e6 polar radius.
    assert Re > 6.36e6


def test_sun_to_earth_mass_ratio():
    """The solar-to-Earth mass ratio is the known ``~3.33e5``.

    ``Ms`` and ``Me`` cross-check each other through the ratio. The guard
    brackets it tightly around ``332946`` so a mantissa or exponent slip in
    either mass fails.
    """
    assert Ms == pytest.approx(1.98847e30, rel=1e-9, abs=0.0)
    assert Rs == pytest.approx(6.957e8, rel=1e-9, abs=0.0)
    ratio = Ms / Me
    # Known solar/Earth mass ratio, order 3.33e5.
    assert 3.32e5 < ratio < 3.34e5


def test_earth_orbit_parameters_are_physical():
    """Earth's semi-major axis and eccentricity sit in their physical ranges.

    ``a_earth`` is the astronomical unit expressed in au (so exactly ``1``),
    and ``e_earth`` is a near-circular eccentricity that must lie in
    ``[0, 1)``. The eccentricity upper-bound check is the edge case: a bound
    orbit cannot have ``e >= 1``.
    """
    # a_earth is 1 au by construction; compare as a number without float ==.
    assert a_earth == pytest.approx(1.0, rel=1e-12, abs=0.0)
    # Bound-orbit eccentricity: 0 <= e < 1, and Earth's is small.
    assert 0.0 <= e_earth < 1.0
    assert e_earth == pytest.approx(0.017, rel=1e-6, abs=0.0)


def test_toi561_star_scaled_from_solar_values():
    """TOI-561 stellar radius, mass and luminosity are fixed solar multiples.

    The source defines ``R_TOI561 = 0.832 Rs``, ``M_TOI561 = 0.805 Ms`` and
    ``L_TOI561 = 0.522 Ls`` (Weiss et al. 2021). Dividing out the solar value
    must recover those multipliers exactly, so a change to either the solar
    anchor or the multiplier is caught.
    """
    assert R_TOI561 / Rs == pytest.approx(0.832, rel=1e-9, abs=0.0)
    assert M_TOI561 / Ms == pytest.approx(0.805, rel=1e-9, abs=0.0)
    assert L_TOI561 / Ls == pytest.approx(0.522, rel=1e-9, abs=0.0)
    # A sub-solar star: each quantity is strictly below its solar anchor.
    assert R_TOI561 < Rs
    assert M_TOI561 < Ms


def test_toi561b_planet_scaled_from_earth_values():
    """TOI-561b radius and mass are fixed Earth multiples of a rocky super-Earth.

    The source defines ``R_TOI561b = 1.37 Re`` and ``M_TOI561b = 2.24 Me``
    (Brinkman et al. 2023). The multipliers are recovered by dividing out the
    Earth anchor, and the pinned mass and radius reproduce a rocky mean
    density. The circular-orbit assumption ``e_TOI561b = 0`` is the boundary
    eccentricity.
    """
    assert R_TOI561b / Re == pytest.approx(1.37, rel=1e-9, abs=0.0)
    assert M_TOI561b / Me == pytest.approx(2.24, rel=1e-9, abs=0.0)
    rho = M_TOI561b / (4 / 3 * np.pi * R_TOI561b**3)
    # Rocky super-Earth mean density, a few 1e3 kg m-3.
    assert 3.0e3 < rho < 6.0e3
    # Circular orbit boundary case.
    assert e_TOI561b == pytest.approx(0.0, abs=1e-12)


def test_jupiter_nominal_values_recover_gas_giant_density():
    """Jupiter's mass and radius give a gas-giant density.

    ``Rjup`` is the IAU 2015 Resolution B3 nominal equatorial radius and
    ``Mjup`` an adopted mass. The
    mean density they imply, about ``1.24 g cm-3``, is the discrimination
    guard: a cm-vs-m slip in the radius moves it by six decades, and swapping
    in the Earth values moves it above ``5 g cm-3``. The mass ratio to Earth
    (about 318) brackets both pins at once.
    """
    assert Mjup == pytest.approx(1.8982e27, rel=1e-9, abs=0.0)
    assert Rjup == pytest.approx(7.1492e7, rel=1e-9, abs=0.0)
    rho = Mjup / (4.0 / 3.0 * np.pi * Rjup**3)
    # Gas-giant mean density in SI: about 1240 kg m-3.
    assert rho == pytest.approx(1240.0, rel=0.02, abs=0.0)
    assert 300.0 < Mjup / Me < 330.0


def test_fxuv_earth_10myr_is_the_sun_track_of_xuv_model_a():
    """``Fxuv_earth_10Myr`` is the Sun track of XUV model A at 1 au.

    Wordsworth et al. (2018, Sect. 3.3, Fig. 9) set ``L_XUV = 1e-3 L`` until
    ``t_sat = 50 Myr`` for the Sun, with ``L`` from Baraffe et al. (2015).
    The 1 M_sun BHAC15 track gives ``log L/L_sun = -0.337`` at
    ``log t = 7.000``; the tolerance covers the 0.001 dex rounding of that
    table. The value 14.67 W m-2 of the Proxima track is more than 20 times
    higher and fails the pin. The flux must also exceed today's value.
    """
    expected = 1e-3 * Ls * 10**-0.337 / (4 * np.pi * au2m**2)
    assert Fxuv_earth_10Myr == pytest.approx(expected, rel=2e-3, abs=0.0)
    assert 1e2 < Fxuv_earth_10Myr / Fxuv_earth_today < 1e3
