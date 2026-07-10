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

from zephyrus.planets_parameters import (
    L_TOI561,
    M_TOI561,
    R_TOI561,
    Ls,
    M_TOI561b,
    Me,
    Ms,
    R_TOI561b,
    Re,
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
    assert Me == pytest.approx(5.9722e24, rel=1e-9)
    assert Re == pytest.approx(6.378e6, rel=1e-9)
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
    assert Ms == pytest.approx(1.98847e30, rel=1e-9)
    assert Rs == pytest.approx(6.957e8, rel=1e-9)
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
    assert a_earth == pytest.approx(1.0, rel=1e-12)
    # Bound-orbit eccentricity: 0 <= e < 1, and Earth's is small.
    assert 0.0 <= e_earth < 1.0
    assert e_earth == pytest.approx(0.017, rel=1e-6)


def test_toi561_star_scaled_from_solar_values():
    """TOI-561 stellar radius, mass and luminosity are fixed solar multiples.

    The source defines ``R_TOI561 = 0.832 Rs``, ``M_TOI561 = 0.805 Ms`` and
    ``L_TOI561 = 0.522 Ls`` (Weiss et al. 2021). Dividing out the solar value
    must recover those multipliers exactly, so a change to either the solar
    anchor or the multiplier is caught.
    """
    assert R_TOI561 / Rs == pytest.approx(0.832, rel=1e-9)
    assert M_TOI561 / Ms == pytest.approx(0.805, rel=1e-9)
    assert L_TOI561 / Ls == pytest.approx(0.522, rel=1e-9)
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
    assert R_TOI561b / Re == pytest.approx(1.37, rel=1e-9)
    assert M_TOI561b / Me == pytest.approx(2.24, rel=1e-9)
    rho = M_TOI561b / (4 / 3 * np.pi * R_TOI561b**3)
    # Rocky super-Earth mean density, a few 1e3 kg m-3.
    assert 3.0e3 < rho < 6.0e3
    # Circular orbit boundary case.
    assert e_TOI561b == pytest.approx(0.0, abs=1e-12)
