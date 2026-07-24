"""Tests for ``src/zephyrus/collision.py``.

Exercises the giant-impact atmospheric mass-loss scaling law of Kegerreis
et al. (2020), ApJL 901, L31: closed-form pins of their Eqn. 1 with
wrong-formula discrimination guards, the density-weighted interacting mass
of Eqn. B1 against the interacting-volume simplification of Eqn. B2, the
total-erosion cap and the grazing limit, the input-validation error
contract, and reference pins against the paper's published simulation
results (their Table 2). See ``docs/How-to/run_tests.md`` for the tier and
marker conventions and ``docs/Validation/collision.md`` for the anchors.
"""

import numpy as np
import pytest

from zephyrus.collision import mass_loss
from zephyrus.constants import G

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Earth-like reference bodies. The bulk density follows from the mass and
# radius, so the mass-ratio and density-ratio terms stay self-consistent.
M_E = 5.972e24  # [kg]
R_E = 6.371e6  # [m]
RHO_E = M_E / (4.0 / 3.0 * np.pi * R_E**3)  # 5513.3 kg/m^3


def _mutual_vesc(M_t, M_i, R_t, R_i):
    """Mutual escape speed of the pair at contact, the paper's Sect. 2."""
    return np.sqrt(2.0 * G * (M_t + M_i) / (R_t + R_i))


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_scaling_law_pins_the_kegerreis_closed_form():
    """Equal twins at the escape speed reproduce Eqn. 1 exactly.

    For two identical Earths head-on at v_c = v_esc every ratio in the
    bracket is unity except the mass ratio M_i/M_tot = 1/2, so Eqn. 1 of
    Kegerreis et al. (2020) collapses to X = 0.64 * 0.5**0.325 = 0.510911,
    a closed form with no geometry left in it. The second pin at
    v_c = 1.5 v_esc keeps the velocity exponent alive. Guards: the
    pre-correction M_i/M_t denominator gives 0.640000 at v_esc and
    1.084 (clamped 1.0) at 1.5 v_esc; a wrong outer exponent of 0.5
    gives 0.807261 at 1.5 v_esc; a wrong velocity exponent of 1 gives
    0.664974. All sit far outside the pin tolerance.
    """
    v_esc = _mutual_vesc(M_E, M_E, R_E, R_E)

    x_1 = mass_loss(v_esc, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 0.0)
    assert x_1 == pytest.approx(0.510911, rel=1e-4)
    # Wrong-denominator guard: M_i/M_t = 1 removes the only non-unit ratio.
    assert abs(x_1 - 0.64) > 0.12

    x_15 = mass_loss(1.5 * v_esc, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 0.0)
    assert x_15 == pytest.approx(0.865494, rel=1e-4)
    # Wrong-exponent guards resolved well above the tolerance.
    assert abs(x_15 - 0.807261) > 0.05  # outer exponent 0.5 instead of 0.65
    assert abs(x_15 - 0.664974) > 0.19  # velocity exponent 1 instead of 2

    # Faster impacts erode more: the two pins are ordered.
    assert x_15 > x_1

    # Absolute anchor: for Earth twins the mutual escape speed equals
    # Earth's own escape speed, 11185.7 m/s. Feeding that as a literal
    # pins the velocity ratio through the code's G and v_esc formula
    # rather than through the test helper, so a wrong gravitational
    # constant or a dropped factor in v_esc shifts this value even
    # though every relative pin above would still pass (a 1 percent G
    # error moves X by 0.65 percent, resolved by the tolerance).
    x_abs = mass_loss(1.11857e4, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 0.0)
    assert x_abs == pytest.approx(0.510911, rel=1e-3)


@pytest.mark.physics_invariant
def test_total_erosion_is_capped_and_grazing_removes_nothing():
    """The loss fraction is capped at 1 and vanishes in the grazing limit.

    Eqn. 1 is "capped at 1 for total erosion": at v_c = 3 v_esc the raw
    power law evaluates to 2.131 for equal twins head-on, so the cap is
    the difference between a usable fraction and unphysical over-removal.
    At b = 1 the interacting mass is exactly zero, so nothing is lost
    however fast the impact; at v_c = 0 nothing is lost however heavy
    the impactor.
    """
    v_esc = _mutual_vesc(M_E, M_E, R_E, R_E)

    x_fast = mass_loss(3.0 * v_esc, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 0.0)
    assert x_fast == pytest.approx(1.0, abs=1e-12)
    # The cap is doing real work here: the raw fit value is 2.131.
    raw = 0.64 * (9.0 * 0.5**0.5) ** 0.65
    assert raw > 2.0

    # Grazing limit: b = 1 gives zero interacting mass, hence zero loss.
    x_graze = mass_loss(3.0 * v_esc, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 1.0)
    assert x_graze == pytest.approx(0.0, abs=1e-12)

    # Zero contact speed: no kinetic energy, no erosion.
    x_still = mass_loss(0.0, M_E, M_E, RHO_E, RHO_E, R_E, R_E, 0.0)
    assert x_still == pytest.approx(0.0, abs=1e-12)


@pytest.mark.physics_invariant
def test_interacting_mass_is_density_weighted_per_eqn_b1():
    """An iron-rich impactor weights the interacting mass, not just volume.

    A half-radius iron impactor (7900 kg/m^3, mass following from its
    density) on a rocky Earth at b = 0.3 and v_c = v_esc pins Eqn. B1 at
    X = 0.281711. The interacting-volume simplification of Eqn. B2, which
    drops the density weighting, gives 0.276082 for the same bodies, a
    0.0056 separation that is two orders of magnitude above the pin
    tolerance, so a regression to the unweighted form cannot pass.
    """
    r_i = 0.5 * R_E
    rho_i = 7900.0  # iron-rich impactor
    m_i = rho_i * 4.0 / 3.0 * np.pi * r_i**3  # 1.0697e24 kg, self-consistent
    v_esc = _mutual_vesc(M_E, m_i, R_E, r_i)

    x_b1 = mass_loss(v_esc, m_i, M_E, rho_i, 5514.0, r_i, R_E, 0.3)
    assert x_b1 == pytest.approx(0.281711, rel=1e-4)

    # Eqn. B2 counterpart (density weighting dropped) for the same bodies.
    f_v = 0.25 * ((R_E + r_i) ** 3 / (R_E**3 + r_i**3)) * (1 - 0.3) ** 2 * (1 + 2 * 0.3)
    x_b2 = 0.64 * ((m_i / (m_i + M_E)) ** 0.5 * (rho_i / 5514.0) ** 0.5 * f_v) ** 0.65
    assert abs(x_b1 - x_b2) > 5.0e-3  # B1 and B2 are resolved apart here


@pytest.mark.physics_invariant
def test_equal_densities_reduce_eqn_b1_to_the_interacting_volume():
    """At equal bulk densities the interacting mass equals the volume form.

    Eqn. B2 is the equal-density limit of Eqn. B1, so for two bodies of
    the same density but different radii the implementation must agree
    with the closed-form interacting-volume expression to float
    precision. This holds even at b = 0, where the common-height cap
    bookkeeping assigns the smaller body a negative cap volume that
    cancels exactly in the sum.
    """
    r_i = 0.4 * R_E
    m_t = RHO_E * 4.0 / 3.0 * np.pi * R_E**3
    m_i = RHO_E * 4.0 / 3.0 * np.pi * r_i**3
    v_esc = _mutual_vesc(m_t, m_i, R_E, r_i)

    for b in (0.0, 0.3, 0.7):
        x = mass_loss(v_esc, m_i, m_t, RHO_E, RHO_E, r_i, R_E, b)
        f_v = 0.25 * ((R_E + r_i) ** 3 / (R_E**3 + r_i**3)) * (1 - b) ** 2 * (1 + 2 * b)
        x_expected = 0.64 * ((m_i / (m_i + m_t)) ** 0.5 * f_v) ** 0.65
        assert x == pytest.approx(x_expected, rel=1e-12)

    # The b sweep is ordered: more grazing, less interacting mass.
    x_head = mass_loss(v_esc, m_i, m_t, RHO_E, RHO_E, r_i, R_E, 0.0)
    x_graze = mass_loss(v_esc, m_i, m_t, RHO_E, RHO_E, r_i, R_E, 0.7)
    assert x_head > x_graze


def test_unphysical_inputs_are_rejected():
    """Out-of-range inputs raise instead of returning plausible nonsense.

    The (1-b)^2 symmetry of the cap height means b = 1.2 would silently
    return a positive loss fraction, so the domain must be enforced.
    Non-positive masses, radii, or densities feed square roots and
    ratios that produce NaN or complex intermediates, and a negative
    contact speed has no physical meaning.
    """
    v = 1.0e4
    good = dict(M_i=M_E, M_t=M_E, rho_i=RHO_E, rho_t=RHO_E, R_i=R_E, R_t=R_E)

    for bad_b in (-0.1, 1.2):
        with pytest.raises(ValueError, match=r'\[0, 1\]'):
            mass_loss(v, b=bad_b, **good)

    for field in ('M_i', 'M_t', 'rho_i', 'rho_t', 'R_i', 'R_t'):
        for bad in (0.0, -1.0, float('nan'), float('inf')):
            broken = dict(good, **{field: bad})
            with pytest.raises(ValueError, match='strictly positive'):
                mass_loss(v, b=0.5, **broken)

    # A NaN or infinite contact speed must raise, not propagate a NaN loss
    # fraction into the caller's mass accounting.
    for bad_v in (-1.0, float('nan'), float('inf')):
        with pytest.raises(ValueError, match='non-negative'):
            mass_loss(bad_v, b=0.5, **good)
    with pytest.raises(ValueError, match=r'\[0, 1\]'):
        mass_loss(v, b=float('nan'), **good)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_scaling_law_reproduces_kegerreis_table2_simulations():
    """The law lands on the paper's own simulation results.

    Three fast, grazing scenarios from Table 2 of Kegerreis et al. (2020)
    (b = 0.7, v_c = 3 v_esc, the regime they report fits tightest), with
    the body radii from their Table 1 and bulk densities following from
    mass and radius. The law reproduces the simulated loss fractions to
    within 4 percent here; the 20 percent tolerance covers the paper's
    stated scatter (9 percent median, about 20 percent worst). The three
    scenarios share the impactor:target mass ratio 10^-0.5, so their
    near-equal simulated losses also exercise the paper's finding that
    the loss is independent of the total mass at fixed mass ratio.
    """
    # (M_t, M_i) in Earth masses; radii in Earth radii (their Table 1);
    # X_sim from their Table 2 (first suite, b = 0.7, v_c = 3 v_esc).
    rows = [
        (10**0.0, 10**-0.5, 0.992, 0.733, 0.520),
        (10**0.25, 10**-0.25, 1.153, 0.856, 0.528),
        (10**-0.5, 10**-1.0, 0.733, 0.538, 0.527),
    ]
    fits = []
    for m_t_e, m_i_e, r_t_e, r_i_e, x_sim in rows:
        m_t, m_i = m_t_e * M_E, m_i_e * M_E
        r_t, r_i = r_t_e * R_E, r_i_e * R_E
        rho_t = m_t / (4.0 / 3.0 * np.pi * r_t**3)
        rho_i = m_i / (4.0 / 3.0 * np.pi * r_i**3)
        v_esc = _mutual_vesc(m_t, m_i, r_t, r_i)
        x_fit = mass_loss(3.0 * v_esc, m_i, m_t, rho_i, rho_t, r_i, r_t, 0.7)
        fits.append(x_fit)
        # 20% tolerance: the paper's stated simulation-to-law scatter.
        assert x_fit == pytest.approx(x_sim, rel=0.20)
        assert 0.0 < x_fit < 1.0

    # Mass-ratio universality: same ratio, near-same loss across a factor
    # of ~5.6 in total mass (the spread here is under 2 percent).
    assert max(fits) == pytest.approx(min(fits), rel=0.05)
