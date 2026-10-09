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

import csv
import io

import numpy as np
import pytest

from zephyrus.collision import (
    _ROCHE2026_COEFFICIENTS,
    _ROCHE2026_TABLE_C_CONSTANTS,
    _roche2026_fit,
    impact_loss,
    mass_loss,
    mass_loss_roche2026,
    mutual_escape_speed,
    specific_impact_energy,
)
from zephyrus.constants import G, M_earth

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


# 12 oracle rows from Roche et al. (2026), arXiv:2610.06077,
# Zenodo doi:10.5281/zenodo.23192423 (authors' scaling_law.csv and impacts file).
_ROCHE2026_ORACLE_CSV = """set,f_atm,M_t_r_earth,M_i_r_earth,M_t_tot_earth,R_t_r_earth,R_i_r_earth,R_ratio,b,gamma,v_c_kms,v_c_v_esc,Q_R_prime_MJkg,f_NF_calc,X_NF_calc,X_FF_calc,X_atm_calc,X_atm_data
A,0.0100205171612713,0.9970246031043678,0.24927721643800169,1.007116430580397,1.0173459965204401,0.67502812051147343,0.66351872698199699,0.29999999999999999,0.20000000000000001,19.27,2,28.541887209782693,0.250519976769084,0.250519976769084,0.3116934236236179,0.56221340039270196,0.61917971501489433
A,0.0100205171612713,0.9970246031043678,0.99691090161648799,1.007116430580397,1.0173459965204401,1.0188975090364123,1.0015250588504592,0.69999999999999996,0.5,22.190000000000001,2,21.77022600953692,0.24891987161603879,0.24891987161603879,0.18184828108652931,0.43076815270256807,0.39370144669541107
A,0.050015747296285197,0.99773538218479818,0.99690269734953563,1.0502651800229117,1.0075866967156619,1.0129610601709269,1.0053338967979459,0,0.5,11.26,1,15.837510342372353,0.6089156069657804,0.24789325220066549,0.027994212660734701,0.2758874648614002,0.20827731001722169
A,0.050015747296285197,0.99773538218479818,0.66531267383424186,1.0502651800229117,1.0075866967156619,0.90174017377539961,0.89495045608950519,0.69999999999999996,0.40000000000000002,15.9,1.5,10.158606138408899,0.30793262783942049,0.16604357767121219,0.029451161572784601,0.19549473924399691,0.20386220995872259
A,0.099920378560586498,0.99771839452616762,0.66531267383424186,1.1084779287976874,1.0037951643625012,0.90174017377539961,0.89833086050786526,0,0.40000000000000002,16.190000000000001,1.5,30.71933776950743,0.60186215686999844,0.38744829857182511,0.059468433528020899,0.44691673209984611,0.4355300137034947
A,0.099920378560586498,0.99771839452616762,0.24927468276732531,1.1084779287976874,1.0037951643625012,0.66951756861323253,0.66698624618144631,0.69999999999999996,0.20000000000000001,15.109999999999999,1.5,6.1880740691416349,0.3272244993677304,0.10392914305129421,0,0.10392914305129421,0.1373774800820031
A,0.2001693807888473,0.99751888605804562,0.99690269734953563,1.2471626643174358,0.99775132424063095,1.0129610601709269,1.0152440147767998,0,0.5,35.450000000000003,3,155.1196221397058,0.76441261046764453,0.76441261046764453,0.23558738953235539,1,0.99999011592834042
A,0.099920378560586498,0.99771839452616762,0.42785599630970522,1.1084779287976874,1.0037951643625012,0.79021890219892132,0.7872312302886818,0.90000000000000002,0.29999999999999999,10.35,1,0.52572239112921526,0.25129319677481021,0,0,0,0.0134884653175223
B,0.050025663901190003,0.34910830517613212,0.23273887011742145,,,,0.8857295894764512,0,0.40000000000000002,,2,25.362807690472678,0.54064633762052949,0.54064633762052949,0.13745119987884191,0.67809753749937141,0.74070392055042422
B,0.049968172845281197,1.9931581262969305,0.49828953157423261,,,,0.67715332501481318,0.69999999999999996,0.20000000000000001,,1.5,9.9733501988963873,0.2365032673150437,0.098076173519816595,0,0.098076173519816595,0.13405858691159919
B,0.050028123964488198,4.9830138432695898,4.9830138432695898,,,,1.0120861984241254,0.29999999999999999,0.5,,1.5,101.63755183729236,0.47484704644023168,0.30125581504975818,0.069771494411668994,0.3710273094614272,0.43950189716330318
B,0.050028123964488198,4.9830138432695898,2.135577361401253,,,,0.80908953470156231,0.90000000000000002,0.29999999999999999,,1.5,3.6938912385513007,0.16266222125253191,0.058404634425177702,0,0.058404634425177702,0.072038965041878802"""


def _get_roche2026_oracle_rows() -> list[dict[str, str]]:
    """Parse the embedded 12 oracle rows from Roche et al. (2026)."""
    return list(csv.DictReader(io.StringIO(_ROCHE2026_ORACLE_CSV.strip())))


def test_roche2026_coefficients_count_and_constants():
    """Verify count of 61 fitted coefficients and exact Table C constants."""
    assert len(_ROCHE2026_COEFFICIENTS) == 61

    # Verify all 61 fitted coefficients are present, float, and non-zero.
    for _name, val in _ROCHE2026_COEFFICIENTS.items():
        assert isinstance(val, float)
        assert val != 0.0
        assert np.isfinite(val)

    # Verify Table C unlisted constant values.
    assert _ROCHE2026_TABLE_C_CONSTANTS['q37'] == pytest.approx(1.0, abs=1e-15)
    assert _ROCHE2026_TABLE_C_CONSTANTS['q48'] == pytest.approx(1.0, abs=1e-15)
    zero_keys = [
        'q16',
        'q17',
        'q18',
        'q23',
        'q24',
        'q25',
        'q26',
        'q27',
        'q28',
        'k26',
        'k27',
        'k36',
        'k37',
        's34',
        's35',
        's36',
    ]
    for key in zero_keys:
        assert _ROCHE2026_TABLE_C_CONSTANTS[key] == 0.0


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_oracle_reproduction():
    """Verify _roche2026_fit matches the 12 oracle rows within tolerance.

    Tests against 12 reference rows from Roche et al. (2026),
    arXiv:2610.06077, Zenodo doi:10.5281/zenodo.23192423.
    """
    rows = _get_roche2026_oracle_rows()
    assert len(rows) == 12

    # Set A: target mass 1 M_E, tight tolerances.
    # Set B: targets 0.35 to 4.98 M_E, looser far-field tolerance.
    for r in rows:
        s = r['set']
        tol_fnf = 1e-12
        tol_xnf = 1e-12
        tol_xff = 1e-10 if s == 'A' else 1e-3
        tol_xat = 1e-10 if s == 'A' else 1e-3

        fnf, xnf, xff, xat = _roche2026_fit(
            float(r['b']),
            float(r['gamma']),
            float(r['v_c_v_esc']),
            float(r['M_t_r_earth']),
            float(r['M_i_r_earth']),
            float(r['Q_R_prime_MJkg']),
            float(r['f_atm']),
            float(r['R_ratio']),
        )
        assert abs(fnf - float(r['f_NF_calc'])) <= tol_fnf
        assert abs(xnf - float(r['X_NF_calc'])) <= tol_xnf
        assert abs(xff - float(r['X_FF_calc'])) <= tol_xff
        assert abs(xat - float(r['X_atm_calc'])) <= tol_xat

    # Vectorized array input check across all 12 rows.
    b_arr = np.array([float(r['b']) for r in rows])
    g_arr = np.array([float(r['gamma']) for r in rows])
    vc_arr = np.array([float(r['v_c_v_esc']) for r in rows])
    mt_arr = np.array([float(r['M_t_r_earth']) for r in rows])
    mi_arr = np.array([float(r['M_i_r_earth']) for r in rows])
    qr_arr = np.array([float(r['Q_R_prime_MJkg']) for r in rows])
    fa_arr = np.array([float(r['f_atm']) for r in rows])
    rr_arr = np.array([float(r['R_ratio']) for r in rows])

    fnf_v, xnf_v, xff_v, xat_v = _roche2026_fit(
        b_arr, g_arr, vc_arr, mt_arr, mi_arr, qr_arr, fa_arr, rr_arr
    )
    assert np.all(abs(fnf_v - np.array([float(r['f_NF_calc']) for r in rows])) <= 1e-12)
    assert np.all(abs(xnf_v - np.array([float(r['X_NF_calc']) for r in rows])) <= 1e-12)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_specific_impact_energy_calculation():
    """Verify Q'_R calculated from SI inputs against oracle column."""
    rows = [r for r in _get_roche2026_oracle_rows() if r['set'] == 'A']
    assert len(rows) == 8

    for r in rows:
        rt = float(r['R_t_r_earth']) * R_E
        ri = float(r['R_i_r_earth']) * R_E
        mt = float(r['M_t_tot_earth']) * M_E
        mi = float(r['M_i_r_earth']) * M_E
        b = float(r['b'])
        vc = float(r['v_c_kms']) * 1e3

        q_r_prime = specific_impact_energy(vc, mi, mt, ri, rt, b)
        expected_qr = float(r['Q_R_prime_MJkg'])
        assert abs(q_r_prime / expected_qr - 1.0) <= 1e-12


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_mutual_escape_speed_calculation():
    """Verify mutual escape speed with M_t^tot vs impact contact velocity ratio."""
    rows = [r for r in _get_roche2026_oracle_rows() if r['set'] == 'A']
    assert len(rows) == 8

    for r in rows:
        rt = float(r['R_t_r_earth']) * R_E
        ri = float(r['R_i_r_earth']) * R_E
        mt_tot = float(r['M_t_tot_earth']) * M_E
        mi = float(r['M_i_r_earth']) * M_E
        vc = float(r['v_c_kms']) * 1e3
        vc_ratio = float(r['v_c_v_esc'])

        v_esc = mutual_escape_speed(mt_tot, mi, rt, ri)
        v_expected = vc / vc_ratio
        assert abs(v_esc / v_expected - 1.0) <= 1e-3

        # Discrimination guard: M_t^r gives 2 to 10 percent error (> 1e-3).
        mt_r = float(r['M_t_r_earth']) * M_E
        v_esc_wrong = mutual_escape_speed(mt_r, mi, rt, ri)
        assert abs(v_esc_wrong / v_expected - 1.0) > 1e-3


@pytest.mark.physics_invariant
def test_roche2026_velocity_floor_property():
    """Verify near-field loss floor holds at impact velocity below escape velocity."""
    rows = _get_roche2026_oracle_rows()[:2]
    assert len(rows) == 2
    for r in rows:
        _, xnf_floor, _, _ = _roche2026_fit(
            float(r['b']),
            float(r['gamma']),
            0.5,
            float(r['M_t_r_earth']),
            float(r['M_i_r_earth']),
            float(r['Q_R_prime_MJkg']),
            float(r['f_atm']),
            float(r['R_ratio']),
        )
        _, xnf_unity, _, _ = _roche2026_fit(
            float(r['b']),
            float(r['gamma']),
            1.0,
            float(r['M_t_r_earth']),
            float(r['M_i_r_earth']),
            float(r['Q_R_prime_MJkg']),
            float(r['f_atm']),
            float(r['R_ratio']),
        )
        assert xnf_floor == pytest.approx(xnf_unity, rel=1e-12, abs=1e-12)


@pytest.mark.physics_invariant
def test_roche2026_zero_atmosphere_fraction():
    """Verify f_atm = 0 returns exactly zero loss fraction without evaluation."""
    x = mass_loss_roche2026(
        v_c=2.0e4,
        M_i=1.0e24,
        M_t=5.972e24,
        R_i=3.0e6,
        R_t=6.371e6,
        b=0.5,
        f_atm=0.0,
    )
    assert x == 0.0

    res = impact_loss(
        'roche2026',
        v_c=2.0e4,
        M_i=1.0e24,
        M_t=5.972e24,
        R_i=3.0e6,
        R_t=6.371e6,
        b=0.5,
        f_atm=0.0,
    )
    assert res.fraction == 0.0
    assert len(res.flags) == 0
    assert res.diagnostics['f_NF'] == 0.0
    assert res.diagnostics['X_NF'] == 0.0
    assert res.diagnostics['X_FF'] == 0.0


@pytest.mark.physics_invariant
def test_roche2026_monotonic_trends():
    """Verify loss fraction decreases with f_atm and b, and increases with v_c."""
    m_t = 5.972e24
    m_i = (0.3 / 0.7) * m_t
    r_t = 6.371e6
    r_i = r_t * (0.3 / 0.7) ** (1.0 / 3.0)
    v_esc = mutual_escape_speed(m_t / (1.0 - 0.01), m_i, r_t, r_i)

    # 1. Monotonic decrease with f_atm
    fa_sweep = [
        mass_loss_roche2026(1.5 * v_esc, m_i, m_t, r_i, r_t, 0.5, fa)
        for fa in (1e-4, 1e-3, 0.01, 0.05, 0.2)
    ]
    assert len(fa_sweep) == 5
    for x1, x2 in zip(fa_sweep[:-1], fa_sweep[1:], strict=True):
        assert x1 > x2

    # 2. Monotonic decrease with b
    b_sweep = [
        mass_loss_roche2026(1.5 * v_esc, m_i, m_t, r_i, r_t, b_val, 0.01)
        for b_val in (0.0, 0.3, 0.5, 0.7, 0.9)
    ]
    assert len(b_sweep) == 5
    for x1, x2 in zip(b_sweep[:-1], b_sweep[1:], strict=True):
        assert x1 > x2

    # 3. Monotonic increase with v_c
    vc_sweep = [
        mass_loss_roche2026(v_fac * v_esc, m_i, m_t, r_i, r_t, 0.5, 0.01)
        for v_fac in (1.0, 1.5, 2.0, 2.5)
    ]
    assert len(vc_sweep) == 4
    for x1, x2 in zip(vc_sweep[:-1], vc_sweep[1:], strict=True):
        assert x2 > x1


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_grazing_continuity_and_value():
    """Verify finite non-zero loss and continuity at grazing impact b -> 1."""
    m_t = 5.972e24
    m_i = (0.3 / 0.7) * m_t
    r_t = 6.371e6
    r_i = r_t * (0.3 / 0.7) ** (1.0 / 3.0)
    v_esc = mutual_escape_speed(m_t / (1.0 - 0.01), m_i, r_t, r_i)

    x_near = mass_loss_roche2026(1.5 * v_esc, m_i, m_t, r_i, r_t, 0.999999, 0.01)
    x_one = mass_loss_roche2026(1.5 * v_esc, m_i, m_t, r_i, r_t, 1.0, 0.01)

    assert np.isfinite(x_near)
    assert np.isfinite(x_one)
    assert x_near == pytest.approx(0.03147, rel=1e-3)
    assert x_one == pytest.approx(0.03147, rel=1e-3)
    assert x_one == pytest.approx(x_near, rel=1e-4)


@pytest.mark.physics_invariant
def test_roche2026_input_contract_and_clamps():
    """Verify input domain validation, stability clamps, and diagnostic flags."""
    m_t = 5.972e24
    m_i = 1.0e24
    r_t = 6.371e6
    r_i = 3.0e6
    v_c = 2.0e4

    # Domain errors
    for b_bad in (-0.1, 1.1, np.nan):
        with pytest.raises(ValueError, match='Impact parameter b'):
            mass_loss_roche2026(v_c, m_i, m_t, r_i, r_t, b_bad, 0.01)

    for m_bad in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match='M_t must be strictly positive'):
            mass_loss_roche2026(v_c, m_i, m_bad, r_i, r_t, 0.5, 0.01)

    for fa_bad in (-0.01, 1.0, 1.5, np.nan):
        with pytest.raises(ValueError, match='f_atm must be in'):
            mass_loss_roche2026(v_c, m_i, m_t, r_i, r_t, 0.5, fa_bad)

    for vc_bad in (-10.0, np.nan, np.inf):
        with pytest.raises(ValueError, match='Collision speed v_c'):
            mass_loss_roche2026(vc_bad, m_i, m_t, r_i, r_t, 0.5, 0.01)

    # Stability clamps and flags via impact_loss
    # f_atm < 1e-6 clamped to 1e-6
    res_fa_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=1e-7
    )
    assert 'f_atm' in res_fa_lo.flags
    res_fa_1e6 = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=1e-6
    )
    assert res_fa_lo.fraction == pytest.approx(res_fa_1e6.fraction, rel=1e-12)

    # f_atm > 0.4 clamped to 0.4
    res_fa_hi = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.5
    )
    assert 'f_atm' in res_fa_hi.flags
    res_fa_04 = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.4
    )
    assert res_fa_hi.fraction == pytest.approx(res_fa_04.fraction, rel=1e-12)

    # Target mass clamped outside [1e-3, 10] M_E
    res_mt_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=1e-4 * M_earth, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'M_t' in res_mt_lo.flags
    res_mt_hi = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=15.0 * M_earth, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'M_t' in res_mt_hi.flags

    # Impactor mass ratio clamped below 1e-3
    res_gamma_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=1e-4 * m_t, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'gamma' in res_gamma_lo.flags

    # Flags for out-of-fitted-range conditions (b > 0.9, gamma > 0.5, v_ratio > 3)
    res_flags = impact_loss(
        'roche2026',
        v_c=4.0 * mutual_escape_speed(m_t / 0.99, 2.0 * m_t, r_t, r_i),
        M_i=2.0 * m_t,
        M_t=m_t,
        R_i=r_i,
        R_t=r_t,
        b=0.95,
        f_atm=0.01,
    )
    assert 'b' in res_flags.flags
    assert 'gamma' in res_flags.flags
    assert 'v_ratio' in res_flags.flags


@pytest.mark.physics_invariant
def test_impact_loss_kegerreis_dispatcher_equivalence():
    """Verify impact_loss with kegerreis2020 matches mass_loss bit-identically."""
    vc = 2.5e4
    mi = 1.0e24
    mt = 6.0e24
    rho_i = 3000.0
    rho_t = 5500.0
    ri = 4.0e6
    rt = 6.4e6
    b = 0.6

    res = impact_loss(
        'kegerreis2020',
        v_c=vc,
        M_i=mi,
        M_t=mt,
        R_i=ri,
        R_t=rt,
        b=b,
        rho_i=rho_i,
        rho_t=rho_t,
    )
    expected = mass_loss(vc, mi, mt, rho_i, rho_t, ri, rt, b)
    assert res.fraction == pytest.approx(expected, rel=1e-15, abs=1e-15)
    assert res.law == 'kegerreis2020'
    assert res.flags == ()
    assert 'v_esc' in res.diagnostics
    assert 'f_M' in res.diagnostics

    # Missing arguments and unknown law
    with pytest.raises(ValueError, match='kegerreis2020 requires'):
        impact_loss('kegerreis2020', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)

    with pytest.raises(ValueError, match='roche2026 requires f_atm'):
        impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)

    with pytest.raises(ValueError, match='Unknown scaling law'):
        impact_loss('unknown_law', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)


@pytest.mark.physics_invariant
def test_roche2026_m_earth_sensitivity():
    """Verify mass scaling sensitivity between 5.972e24 and 5.9724e24 is below 1e-5."""
    rows = _get_roche2026_oracle_rows()
    r = rows[1]
    b = float(r['b'])
    g = float(r['gamma'])
    vc = float(r['v_c_v_esc'])
    mt = float(r['M_t_r_earth'])
    mi = float(r['M_i_r_earth'])
    qr = float(r['Q_R_prime_MJkg'])
    fa = float(r['f_atm'])
    rr = float(r['R_ratio'])

    _, _, _, x_nom = _roche2026_fit(b, g, vc, mt, mi, qr, fa, rr)
    factor = 5.972e24 / 5.9724e24
    _, _, _, x_alt = _roche2026_fit(b, g, vc, mt * factor, mi * factor, qr, fa, rr)
    diff = abs(x_alt - x_nom)
    assert diff < 1.0e-5
    assert abs(x_nom - float(r['X_atm_calc'])) <= 1e-10
