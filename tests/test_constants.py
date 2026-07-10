"""Tests for ``src/zephyrus/constants.py``.

``constants.py`` is a utility source (physical constants and unit
conversions), so it is exempt from the physics-invariant requirement, but a
wrong constant is a silent physics error: every downstream rate inherits it.
These tests pin the shipped values and, more importantly, assert the
cross-consistency relations that a single-line typo would break (the SI and
cgs gravitational constants must differ by exactly ``1e3``; the au-to-cm and
au-to-m factors by exactly ``1e2``; the composite ``erg cm-2 s-1`` to
``W m-2`` factor must equal the erg-to-joule factor times ``1e4``). Each pin
carries a discrimination guard so a plausible wrong value fails loudly.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions.
"""

import pytest

from zephyrus.constants import (
    G,
    G_cgs,
    au2cm,
    au2m,
    c,
    erg2joule,
    ergcm2stoWm2,
    ergpersecondtowatt,
    kb,
    s2yr,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def test_gravitational_constant_si_cgs_consistency():
    """The SI and cgs gravitational constants differ by exactly ``1e3``.

    ``G`` in ``m^3 kg^-1 s^-2`` and ``G_cgs`` in ``cm^3 g^-1 s^-2`` must be
    related by the ``m^3/cm^3 = 1e6`` volume factor over the ``kg/g = 1e3``
    mass factor, i.e. ``G_cgs = 1e3 * G``. The escape formula divides by the
    SI ``G``; a cgs value leaking into an SI expression would inflate every
    rate by ``1e3``, so the ratio is the discrimination guard.
    """
    assert G == pytest.approx(6.6743e-11, rel=1e-12)
    # cgs/SI consistency: the ratio must be exactly 1e3.
    assert G_cgs / G == pytest.approx(1e3, rel=1e-9)
    # Discrimination: the SI value is order 1e-11, the cgs value order 1e-8;
    # they must not be interchangeable, so a swap moves the rate by 1e3.
    assert G < 1e-9 < G_cgs


def test_speed_of_light_is_the_exact_si_value():
    """The speed of light is the exact defined SI value, not a rounded ``3e8``.

    Since 1983 the metre is defined so that ``c = 299792458 m s-1`` exactly.
    The boundary case a careless value would take is ``3e8``; the guard
    asserts the shipped constant sits more than ``2e5`` away from it.
    """
    assert c == pytest.approx(299792458.0, rel=1e-12)
    # Discrimination against the common 3e8 rounding.
    assert abs(c - 3e8) > 2e5


def test_unit_conversions_are_mutually_consistent():
    """The length and energy-flux conversions obey their defining ratios.

    ``au2cm`` must be ``1e2`` times ``au2m`` (metre-to-centimetre), and the
    composite ``erg cm-2 s-1`` to ``W m-2`` factor must equal the erg-to-joule
    factor scaled by the ``cm^-2`` to ``m^-2`` area factor of ``1e4``. The
    per-second erg-to-watt factor must equal the plain erg-to-joule factor.
    These relations are what the MORS-to-escape flux hand-off relies on.
    """
    assert erg2joule == pytest.approx(1e-7, rel=1e-12)
    # Length: centimetres are exactly 1e2 metres.
    assert au2cm / au2m == pytest.approx(1e2, rel=1e-9)
    # Energy flux: erg s-1 cm-2 -> W m-2 folds erg->J (1e-7) with cm^-2->m^-2 (1e4).
    assert ergcm2stoWm2 == pytest.approx(erg2joule * 1e4, rel=1e-9)
    # Power: erg s-1 -> W is the plain erg->J factor.
    assert ergpersecondtowatt == pytest.approx(erg2joule, rel=1e-12)


def test_seconds_to_year_uses_the_365_day_convention():
    """``s2yr`` inverts a 365-day year, distinct from the 365.25-day Julian year.

    The shipped factor is ``1 / (3600 * 24 * 365)``. The discrimination guard
    is the Julian-year value ``1 / (3600 * 24 * 365.25)``: the two differ by
    about ``0.07%``, well above the pin tolerance, so the test documents which
    convention is in force and would fail if it were switched.
    """
    assert s2yr == pytest.approx(1.0 / (3600 * 24 * 365), rel=1e-12)
    julian = 1.0 / (3600 * 24 * 365.25)
    # The 365-day and Julian conventions are resolvably different.
    assert s2yr != pytest.approx(julian, rel=1e-4)
    assert s2yr > julian  # a shorter year makes each second a larger fraction


def test_boltzmann_constant_magnitude_and_value():
    """The Boltzmann constant is the shipped three-significant-figure value.

    ``kb`` is stored as ``1.38e-23 J K-1``. The guard brackets its order of
    magnitude so a decimal-place slip (``1.38e-22`` or ``1.38e-24``) fails,
    and pins the value so a mantissa typo is caught.
    """
    assert kb == pytest.approx(1.38e-23, rel=1e-12)
    # Order-of-magnitude bracket: a single decimal-place slip lands outside.
    assert 1e-23 < kb < 2e-23
