"""Tests for ``src/zephyrus/collision.py``.

Exercises the giant-impact atmospheric mass-loss scaling law of Kegerreis
et al. (2020), ApJL 901, L31: closed-form pins of their Eqn. 1 with
wrong-formula discrimination guards, the density-weighted interacting mass
of Eqn. B1 against the interacting-volume simplification of Eqn. B2, the
total-erosion cap and the grazing limit, the input-validation error
contract, and reference pins against the paper's published simulation
results (their Table 2).

Also exercises the Roche et al. (2026) scaling law and its physical invariants:
- Boundedness: loss fraction bounded in [0, 1] across parameter domain.
- Airless target: zero loss fraction for f_atm = 0.
- Monotonicity: loss fraction increases monotonically with collision speed v_c.
- Continuity: continuous at grazing impact parameter b -> 1.
- Velocity floor: near-field loss floor holds at v_c <= v_esc.
- Validity range: out-of-fitted-range diagnostics flags with 1% relative tolerance.
- Numerical stability: input clamping to prevent overflow outside stable bounds.

See ``docs/How-to/run_tests.md`` for the tier and marker conventions and
``docs/Validation/collision.md`` for the anchors.
"""

import csv
import io

import numpy as np
import pytest

from zephyrus.collision import (
    _ROCHE2026_COEFFICIENTS,
    _ROCHE2026_RANGE_RTOL,
    _ROCHE2026_STABLE_RANGE,
    _ROCHE2026_TABLE_C_CONSTANTS,
    ROCHE2026_FITTED_RANGE,
    _roche2026_fit,
    impact_loss,
    mass_loss,
    mass_loss_roche2026,
    mutual_escape_speed,
    specific_impact_energy,
)
from zephyrus.constants import G
from zephyrus.planets_parameters import Me

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

# Earth-like reference bodies for Kegerreis scaling law pins.
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

    # Absolute anchor: for Earth twins mutual escape speed equals 11185.7 m/s.
    # Feeding this literal directly pins G and v_esc without helper cancellation.
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
# Note for Set B: the impactor masses come from the nominal gamma, which
# explains the 1.02e-4 residual between published columns and fit output.
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
    # Set B: targets 0.35 to 4.98 M_E, far-field tolerance 2e-4.
    for r in rows:
        s = r['set']
        tol_fnf = 1e-12
        tol_xnf = 1e-12
        tol_xff = 1e-10 if s == 'A' else 2e-4
        tol_xat = 1e-10 if s == 'A' else 2e-4

        fnf, xnf, xff, xat = _roche2026_fit(
            b=float(r['b']),
            gamma=float(r['gamma']),
            v_c_v_esc=float(r['v_c_v_esc']),
            M_t_earth=float(r['M_t_r_earth']),
            M_i_earth=float(r['M_i_r_earth']),
            Q_R_prime_MJ=float(r['Q_R_prime_MJkg']),
            f_atm=float(r['f_atm']),
            R_ratio=float(r['R_ratio']),
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
        b=b_arr,
        gamma=g_arr,
        v_c_v_esc=vc_arr,
        M_t_earth=mt_arr,
        M_i_earth=mi_arr,
        Q_R_prime_MJ=qr_arr,
        f_atm=fa_arr,
        R_ratio=rr_arr,
    )
    assert np.all(abs(fnf_v - np.array([float(r['f_NF_calc']) for r in rows])) <= 1e-12)
    assert np.all(abs(xnf_v - np.array([float(r['X_NF_calc']) for r in rows])) <= 1e-12)
    tol_xff_arr = np.array([1e-10 if r['set'] == 'A' else 2e-4 for r in rows])
    tol_xat_arr = np.array([1e-10 if r['set'] == 'A' else 2e-4 for r in rows])
    assert np.all(abs(xff_v - np.array([float(r['X_FF_calc']) for r in rows])) <= tol_xff_arr)
    assert np.all(abs(xat_v - np.array([float(r['X_atm_calc']) for r in rows])) <= tol_xat_arr)

    # Public API entry point verification on reference oracle row 1.
    r0 = rows[0]
    res0 = impact_loss(
        'roche2026',
        v_c=float(r0['v_c_kms']) * 1e3,
        M_i=float(r0['M_i_r_earth']) * Me,
        M_t=float(r0['M_t_r_earth']) * Me,
        R_i=float(r0['R_i_r_earth']) * R_E,
        R_t=float(r0['R_t_r_earth']) * R_E,
        b=float(r0['b']),
        f_atm=float(r0['f_atm']),
    )
    assert res0.fraction == pytest.approx(float(r0['X_atm_calc']), abs=2e-5)
    assert res0.flags == ()
    assert res0.diagnostics['f_NF'] == pytest.approx(float(r0['f_NF_calc']), abs=1e-10)
    assert res0.diagnostics['X_NF'] == pytest.approx(float(r0['X_NF_calc']), abs=1e-4)
    assert res0.diagnostics['X_FF'] == pytest.approx(float(r0['X_FF_calc']), abs=1e-4)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_specific_impact_energy_calculation():
    """Verify Q'_R calculated from SI inputs against oracle column."""
    rows = [r for r in _get_roche2026_oracle_rows() if r['set'] == 'A']
    assert len(rows) == 8

    for r in rows:
        rt = float(r['R_t_r_earth']) * R_E
        ri = float(r['R_i_r_earth']) * R_E
        mt = float(r['M_t_tot_earth']) * Me
        mi = float(r['M_i_r_earth']) * Me
        b = float(r['b'])
        vc = float(r['v_c_kms']) * 1e3

        q_r_prime = specific_impact_energy(vc, mi, mt, ri, rt, b)
        expected_qr = float(r['Q_R_prime_MJkg'])
        assert abs(q_r_prime / expected_qr - 1.0) <= 1e-12

    # Verify M_t_tot keyword parameter name
    q_kw = specific_impact_energy(vc, mi, M_t_tot=mt, R_i=ri, R_t=rt, b=b)
    assert q_kw == q_r_prime


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_mutual_escape_speed_calculation():
    """Verify mutual escape speed with M_t^tot vs impact contact velocity ratio."""
    rows = [r for r in _get_roche2026_oracle_rows() if r['set'] == 'A']
    assert len(rows) == 8

    for r in rows:
        rt = float(r['R_t_r_earth']) * R_E
        ri = float(r['R_i_r_earth']) * R_E
        mt_tot = float(r['M_t_tot_earth']) * Me
        mi = float(r['M_i_r_earth']) * Me
        vc = float(r['v_c_kms']) * 1e3
        vc_ratio = float(r['v_c_v_esc'])

        v_esc = mutual_escape_speed(mt_tot, mi, rt, ri)
        v_expected = vc / vc_ratio
        assert abs(v_esc / v_expected - 1.0) <= 1e-3

        # Discrimination guard: M_t^r gives 0.2 to 6 percent error (> 1e-3).
        mt_r = float(r['M_t_r_earth']) * Me
        v_esc_wrong = mutual_escape_speed(mt_r, mi, rt, ri)
        assert abs(v_esc_wrong / v_expected - 1.0) > 1e-3


def test_mutual_escape_speed_validation():
    """Verify error raises on invalid radius or mass inputs."""
    with pytest.raises(ValueError, match='M_1 must be strictly positive'):
        mutual_escape_speed(0.0, 1e24, 1e6, 1e6)
    with pytest.raises(ValueError, match='M_2 must be strictly positive'):
        mutual_escape_speed(1e24, -1e24, 1e6, 1e6)
    with pytest.raises(ValueError, match='R_1 must be strictly positive'):
        mutual_escape_speed(1e24, 1e24, 0.0, 1e6)
    with pytest.raises(ValueError, match='R_2 must be strictly positive'):
        mutual_escape_speed(1e24, 1e24, 1e6, -1e6)
    with pytest.raises(ValueError, match='M_1 must be strictly positive and finite'):
        mutual_escape_speed(np.inf, 1e24, 1e6, 1e6)
    with pytest.raises(ValueError, match='R_1 must be strictly positive and finite'):
        mutual_escape_speed(1e24, 1e24, np.nan, 1e6)
    with pytest.raises(TypeError, match='must be a scalar'):
        mutual_escape_speed(np.array([1e24, 2e24]), 1e24, 1e6, 1e6)


def test_specific_impact_energy_validation_and_limits():
    """Verify error raises on invalid inputs and geometric limits."""
    with pytest.raises(ValueError, match='Impact parameter b'):
        specific_impact_energy(1e4, 1e24, 1e24, 1e6, 1e6, -0.1)
    with pytest.raises(ValueError, match='Impact parameter b'):
        specific_impact_energy(1e4, 1e24, 1e24, 1e6, 1e6, 1.1)
    with pytest.raises(ValueError, match='M_i must be strictly positive and finite'):
        specific_impact_energy(1e4, 0.0, 1e24, 1e6, 1e6, 0.5)
    with pytest.raises(ValueError, match='M_t must be strictly positive and finite'):
        specific_impact_energy(1e4, 1e24, -1e24, 1e6, 1e6, 0.5)
    with pytest.raises(ValueError, match='R_i must be strictly positive and finite'):
        specific_impact_energy(1e4, 1e24, 1e24, 0.0, 1e6, 0.5)
    with pytest.raises(ValueError, match='R_t must be strictly positive and finite'):
        specific_impact_energy(1e4, 1e24, 1e24, 1e6, -1e6, 0.5)
    with pytest.raises(ValueError, match='Collision speed v_c'):
        specific_impact_energy(-1.0, 1e24, 1e24, 1e6, 1e6, 0.5)
    with pytest.raises(TypeError, match='must be a scalar'):
        specific_impact_energy(np.array([1e4, 2e4]), 1e24, 1e24, 1e6, 1e6, 0.5)

    # Complete capture branch (impact_param + R_i <= R_t)
    qr_headon = specific_impact_energy(1e4, 1e23, 1e25, 1e6, 1e7, 0.0)
    assert qr_headon > 0.0

    # Grazing miss branch (impact_param >= R_t + R_i)
    qr_miss = specific_impact_energy(1e4, 1e24, 1e24, 1e6, 1e6, 1.0)
    assert qr_miss == pytest.approx(0.0, abs=1e-15)


@pytest.mark.physics_invariant
def test_roche2026_velocity_floor_property():
    """Verify near-field loss floor holds at impact velocity below escape velocity."""
    rows = _get_roche2026_oracle_rows()[:2]
    assert len(rows) == 2
    for r in rows:
        _, xnf_floor, _, _ = _roche2026_fit(
            b=float(r['b']),
            gamma=float(r['gamma']),
            v_c_v_esc=0.5,
            M_t_earth=float(r['M_t_r_earth']),
            M_i_earth=float(r['M_i_r_earth']),
            Q_R_prime_MJ=float(r['Q_R_prime_MJkg']),
            f_atm=float(r['f_atm']),
            R_ratio=float(r['R_ratio']),
        )
        _, xnf_unity, _, _ = _roche2026_fit(
            b=float(r['b']),
            gamma=float(r['gamma']),
            v_c_v_esc=1.0,
            M_t_earth=float(r['M_t_r_earth']),
            M_i_earth=float(r['M_i_r_earth']),
            Q_R_prime_MJ=float(r['Q_R_prime_MJkg']),
            f_atm=float(r['f_atm']),
            R_ratio=float(r['R_ratio']),
        )
        assert xnf_floor == pytest.approx(xnf_unity, rel=1e-12, abs=1e-12)


@pytest.mark.physics_invariant
def test_roche2026_zero_atmosphere_fraction():
    """Verify f_atm = 0 returns exactly zero loss fraction without evaluation."""
    x = mass_loss_roche2026(
        v_c=2.0e4,
        M_i=1.0e24,
        M_t=Me,
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
        M_t=Me,
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

    # Diagnostics computed from physical inputs (Item 15)
    v_esc_expected = mutual_escape_speed(Me, 1.0e24, 6.371e6, 3.0e6)
    assert res.diagnostics['v_esc'] == v_esc_expected
    assert res.diagnostics['v_ratio'] == 2.0e4 / v_esc_expected
    assert res.diagnostics['Q_R_prime'] == specific_impact_energy(
        2.0e4, 1.0e24, Me, 3.0e6, 6.371e6, 0.5
    )


@pytest.mark.physics_invariant
def test_roche2026_collision_speed_monotonic_grid():
    """Verify loss fraction increases monotonically with v_c over parameter grid."""
    masses = [0.5 * Me, 1.0 * Me, 3.0 * Me]
    gammas = [0.1, 0.25, 0.4]
    bs = [0.1, 0.5, 0.8]
    fatms = [0.02, 0.05, 0.15]
    v_ratios = [1.0, 1.5, 2.0, 2.5, 3.0]

    for mt in masses:
        rt = 6.371e6 * (mt / Me) ** (1.0 / 3.0)
        for g in gammas:
            mi = mt * (g / (1.0 - g))
            ri = rt * (mi / mt) ** (1.0 / 3.0)
            for b in bs:
                for fa in fatms:
                    vesc = mutual_escape_speed(mt / (1.0 - fa), mi, rt, ri)
                    xs = [
                        mass_loss_roche2026(vr * vesc, mi, mt, ri, rt, b, fa) for vr in v_ratios
                    ]
                    for x1, x2 in zip(xs[:-1], xs[1:], strict=True):
                        assert x2 >= x1
                        if x1 < 1.0:
                            assert x2 > x1


@pytest.mark.physics_invariant
def test_roche2026_grazing_continuity_and_value():
    """Verify finite non-zero loss and continuity at grazing impact b -> 1 for Earth-mass target with gamma = 0.3, f_atm = 0.01, v/v_esc = 1.5."""
    m_t = Me
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
    m_t = Me
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
    assert set(ROCHE2026_FITTED_RANGE.keys()) == {
        'f_atm',
        'M_t_earth',
        'gamma',
        'b',
        'R_ratio',
        'v_ratio',
    }
    assert _ROCHE2026_RANGE_RTOL == pytest.approx(0.01)
    assert set(_ROCHE2026_STABLE_RANGE.keys()) == {'f_atm', 'M_t_earth', 'gamma'}

    # f_atm < 1e-6 clamped to 1e-6 in fit
    res_fa_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=1e-7
    )
    assert 'f_atm' in res_fa_lo.flags
    assert res_fa_lo.diagnostics['clamped']['f_atm'] == pytest.approx(1.0e-6)
    assert 0.0 <= res_fa_lo.fraction <= 1.0
    res_fa_1e6 = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=1e-6
    )
    assert res_fa_lo.fraction == pytest.approx(res_fa_1e6.fraction, rel=1e-5)

    # f_atm > 0.4 clamped to 0.4 in fit
    res_fa_hi = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.5
    )
    assert 'f_atm' in res_fa_hi.flags
    assert res_fa_hi.diagnostics['clamped']['f_atm'] == pytest.approx(0.4)
    assert 0.0 <= res_fa_hi.fraction <= 1.0

    # Target mass clamped outside [1e-3, 10] M_E (Item 13 value assertions)
    res_mt_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=1e-4 * Me, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'M_t_earth' in res_mt_lo.flags
    assert res_mt_lo.diagnostics['clamped']['M_t_earth'] == pytest.approx(1.0e-3)
    assert res_mt_lo.diagnostics['M_t_earth'] == pytest.approx(1e-4)
    _, _, _, x_expected_lo = _roche2026_fit(
        b=0.5,
        gamma=res_mt_lo.diagnostics['gamma'],
        v_c_v_esc=res_mt_lo.diagnostics['v_ratio'],
        M_t_earth=1.0e-3,
        M_i_earth=1.0e-3 * (m_i / (1e-4 * Me)),
        Q_R_prime_MJ=res_mt_lo.diagnostics['Q_R_prime'],
        f_atm=0.01,
        R_ratio=r_i / r_t,
    )
    assert res_mt_lo.fraction == pytest.approx(x_expected_lo, abs=1e-12)

    res_mt_hi = impact_loss(
        'roche2026', v_c=v_c, M_i=m_i, M_t=15.0 * Me, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'M_t_earth' in res_mt_hi.flags
    assert res_mt_hi.diagnostics['clamped']['M_t_earth'] == pytest.approx(10.0)
    assert res_mt_hi.diagnostics['M_t_earth'] == pytest.approx(15.0)
    _, _, _, x_expected_hi = _roche2026_fit(
        b=0.5,
        gamma=res_mt_hi.diagnostics['gamma'],
        v_c_v_esc=res_mt_hi.diagnostics['v_ratio'],
        M_t_earth=10.0,
        M_i_earth=10.0 * (m_i / (15.0 * Me)),
        Q_R_prime_MJ=res_mt_hi.diagnostics['Q_R_prime'],
        f_atm=0.01,
        R_ratio=r_i / r_t,
    )
    assert res_mt_hi.fraction == pytest.approx(x_expected_hi, abs=1e-12)

    # Impactor mass ratio clamped below 1e-3 (Item 13 value assertion)
    res_gamma_lo = impact_loss(
        'roche2026', v_c=v_c, M_i=1e-4 * m_t, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'gamma' in res_gamma_lo.flags
    assert res_gamma_lo.diagnostics['clamped']['gamma'] == pytest.approx(1.0e-3)
    _, _, _, x_expected_gamma = _roche2026_fit(
        b=0.5,
        gamma=1.0e-3,
        v_c_v_esc=res_gamma_lo.diagnostics['v_ratio'],
        M_t_earth=1.0,
        M_i_earth=1.0 * (1e-4 * m_t / m_t),
        Q_R_prime_MJ=res_gamma_lo.diagnostics['Q_R_prime'],
        f_atm=0.01,
        R_ratio=r_i / r_t,
    )
    assert res_gamma_lo.fraction == pytest.approx(x_expected_gamma, abs=1e-12)

    # gamma = 0.9 is not clamped (Item 13)
    res_gamma_09 = impact_loss(
        'roche2026', v_c=v_c, M_i=9.0 * m_t, M_t=m_t, R_i=r_i, R_t=r_t, b=0.5, f_atm=0.01
    )
    assert 'gamma' in res_gamma_09.flags
    assert 'gamma' not in res_gamma_09.diagnostics['clamped']
    _, _, _, x_expected_g09 = _roche2026_fit(
        b=0.5,
        gamma=0.9,
        v_c_v_esc=res_gamma_09.diagnostics['v_ratio'],
        M_t_earth=1.0,
        M_i_earth=9.0,
        Q_R_prime_MJ=res_gamma_09.diagnostics['Q_R_prime'],
        f_atm=0.01,
        R_ratio=r_i / r_t,
    )
    assert res_gamma_09.fraction == pytest.approx(x_expected_g09, abs=1e-12)

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
    assert res_flags.diagnostics['clamped'] == {}
    assert res_flags.diagnostics['b'] == pytest.approx(0.95)
    assert res_flags.diagnostics['gamma'] == pytest.approx(2.0 / 3.0)
    assert res_flags.diagnostics['f_atm'] == pytest.approx(0.01)
    assert res_flags.diagnostics['M_t_earth'] == pytest.approx(1.0)
    assert res_flags.diagnostics['v_ratio'] > 3.0

    # Speed below escape speed is not flagged (held by near-field velocity floor)
    res_sub_vesc = impact_loss(
        'roche2026',
        v_c=0.5 * mutual_escape_speed(m_t / 0.99, m_i, r_t, r_i),
        M_i=m_i,
        M_t=m_t,
        R_i=r_i,
        R_t=r_t,
        b=0.5,
        f_atm=0.01,
    )
    assert 'v_ratio' not in res_sub_vesc.flags
    assert res_sub_vesc.diagnostics['v_ratio'] == pytest.approx(0.5)

    # v/v_esc = 5 is not clamped and differs from value at 3 (Item 13)
    vesc_eval = mutual_escape_speed(
        m_t / (1.0 - 0.01), 0.1 * m_t, r_t, r_t * 0.1 ** (1.0 / 3.0)
    )
    res_v3 = impact_loss(
        'roche2026',
        v_c=3.0 * vesc_eval,
        M_i=0.1 * m_t,
        M_t=m_t,
        R_i=r_t * 0.1 ** (1.0 / 3.0),
        R_t=r_t,
        b=0.7,
        f_atm=0.01,
    )
    res_v5 = impact_loss(
        'roche2026',
        v_c=5.0 * vesc_eval,
        M_i=0.1 * m_t,
        M_t=m_t,
        R_i=r_t * 0.1 ** (1.0 / 3.0),
        R_t=r_t,
        b=0.7,
        f_atm=0.01,
    )
    assert 'v_ratio' in res_v5.flags
    assert 'v_ratio' not in res_v5.diagnostics['clamped']
    assert res_v5.fraction > res_v3.fraction

    # Physical quantities derived from physical inputs
    res_phys = impact_loss(
        'roche2026',
        v_c=v_c,
        M_i=m_i,
        M_t=15.0 * Me,
        R_i=r_i,
        R_t=r_t,
        b=0.5,
        f_atm=0.5,
    )
    assert res_phys.diagnostics['v_esc'] == mutual_escape_speed(15.0 * Me / 0.5, m_i, r_t, r_i)
    assert 'M_t_earth' in res_phys.flags
    assert 'f_atm' in res_phys.flags
    assert res_phys.diagnostics['clamped']['M_t_earth'] == pytest.approx(10.0)
    assert res_phys.diagnostics['clamped']['f_atm'] == pytest.approx(0.4)


@pytest.mark.physics_invariant
def test_impact_loss_kegerreis_dispatcher_equivalence():
    """Verify impact_loss with kegerreis2020 matches mass_loss bit-identically."""
    cases = [
        (2.5e4, 1.0e24, 6.0e24, 3000.0, 5500.0, 4.0e6, 6.4e6, 0.6),
        (2.0e4, 1.0e24, 6.0e24, 3000.0, 5500.0, 4.0e6, 6.4e6, 1.0),
        (2.0e4, 1.0e24, 6.0e24, 3000.0, 5500.0, 4.0e6, 6.4e6, 0.0),
        (1.0e6, 2.0e24, 5.0e24, 3500.0, 5000.0, 5.0e6, 6.0e6, 0.2),
        (3.0e4, 5.0e23, 6.0e24, 7000.0, 3000.0, 2.5e6, 6.4e6, 0.4),
        (2.2e4, 1.5e24, 6.0e24, 4500.0, 4500.0, 4.2e6, 6.4e6, 0.5),
    ]

    for vc, mi, mt, rho_i, rho_t, ri, rt, b in cases:
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
        assert res.fraction == expected
        assert res.law == 'kegerreis2020'
        assert res.flags == ()
        assert 'v_esc' in res.diagnostics
        assert 'f_M' in res.diagnostics

        # Hand calculations check for Kegerreis diagnostics (Item 15)
        v_esc_hand = np.sqrt(2.0 * G * (mt + mi) / (rt + ri))
        assert res.diagnostics['v_esc'] == pytest.approx(v_esc_hand, rel=1e-12)
        assert res.diagnostics['v_ratio'] == pytest.approx(vc / v_esc_hand, rel=1e-12)
        assert res.diagnostics['gamma'] == pytest.approx(mi / (mi + mt), rel=1e-12)
        assert res.diagnostics['mass_ratio'] == pytest.approx(mi / mt, rel=1e-12)
        assert 0.0 <= res.diagnostics['f_M'] <= 1.0

    # Missing arguments and unknown law
    vc, mi, mt, ri, rt, b = 2.5e4, 1.0e24, 6.0e24, 4.0e6, 6.4e6, 0.6
    with pytest.raises(ValueError, match='kegerreis2020 requires'):
        impact_loss('kegerreis2020', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)

    with pytest.raises(ValueError, match='roche2026 requires f_atm'):
        impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)

    with pytest.raises(ValueError, match='Unknown scaling law'):
        impact_loss('unknown_law', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b)


@pytest.mark.physics_invariant
def test_roche2026_m_earth_sensitivity(monkeypatch):
    """Verify Earth-mass definition sensitivity between 5.972e24 and 5.9724e24 is below 1e-5."""
    import zephyrus.collision as zc

    vc = 2.0e4
    mi = 1.0e24
    mt = 5.972e24
    ri = 3.0e6
    rt = 6.371e6
    b = 0.5
    fa = 0.05

    monkeypatch.setattr(zc, 'Me', 5.972e24)
    res_nom = impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=fa)

    monkeypatch.setattr(zc, 'Me', 5.9724e24)
    res_alt = impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=fa)

    diff = abs(res_alt.fraction - res_nom.fraction)
    # Measured difference is 8.22e-6, safely bounded by 1e-5.
    assert diff == pytest.approx(8.2225e-6, rel=1e-2)
    assert diff < 1.0e-5


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_roche2026_oracle_set_a_impact_loss():
    """Verify all 8 Set A rows and all 12 oracle rows through public impact_loss."""
    rows = _get_roche2026_oracle_rows()
    set_a = [r for r in rows if r['set'] == 'A']
    assert len(set_a) == 8

    for i, r in enumerate(set_a):
        vc = float(r['v_c_kms']) * 1e3
        res = impact_loss(
            'roche2026',
            v_c=vc,
            M_i=float(r['M_i_r_earth']) * Me,
            M_t=float(r['M_t_r_earth']) * Me,
            R_i=float(r['R_i_r_earth']) * R_E,
            R_t=float(r['R_t_r_earth']) * R_E,
            b=float(r['b']),
            f_atm=float(r['f_atm']),
        )
        assert res.fraction == pytest.approx(float(r['X_atm_calc']), abs=2e-4)
        v_expected = vc / float(r['v_c_v_esc'])
        assert abs(res.diagnostics['v_esc'] / v_expected - 1.0) <= 1e-3

        # Row 6 has R_ratio = 1.015244 > 1.01, flagged with ('R_ratio',) per Ruling 13.
        # All other Set A rows have no flags.
        if i == 6:
            assert res.flags == ('R_ratio',)
        else:
            assert res.flags == ()


@pytest.mark.physics_invariant
def test_roche2026_fitted_range_boundaries_and_r_ratio():
    """Verify boundary flags with 1% relative tolerance and R_ratio flag per Ruling 13."""
    rt = R_E
    ri = 0.5 * R_E
    mt = Me
    mi = (0.25 / 0.75) * mt
    fa = 0.05
    b = 0.5
    vesc = mutual_escape_speed(mt / (1.0 - fa), mi, rt, ri)
    vc = 2.0 * vesc

    # f_atm: (0.01, 0.2)
    assert (
        'f_atm'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=0.01 * 0.995
        ).flags
    )
    assert (
        'f_atm'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=0.01 * 0.985
        ).flags
    )
    assert (
        'f_atm'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=0.2 * 1.005
        ).flags
    )
    assert (
        'f_atm'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=0.2 * 1.015
        ).flags
    )

    # M_t_earth: (0.35, 5.0)
    assert (
        'M_t_earth'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=0.35 * 0.995 * Me, R_i=ri, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'M_t_earth'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=0.35 * 0.985 * Me, R_i=ri, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'M_t_earth'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=5.0 * 1.005 * Me, R_i=ri, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'M_t_earth'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=5.0 * 1.015 * Me, R_i=ri, R_t=rt, b=b, f_atm=fa
        ).flags
    )

    # gamma: (0.1, 0.5)
    g_lo_in = 0.1 * 0.995
    assert (
        'gamma'
        not in impact_loss(
            'roche2026',
            v_c=vc,
            M_i=(g_lo_in / (1.0 - g_lo_in)) * mt,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )
    g_lo_out = 0.1 * 0.985
    assert (
        'gamma'
        in impact_loss(
            'roche2026',
            v_c=vc,
            M_i=(g_lo_out / (1.0 - g_lo_out)) * mt,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )
    g_hi_in = 0.5 * 1.005
    assert (
        'gamma'
        not in impact_loss(
            'roche2026',
            v_c=vc,
            M_i=(g_hi_in / (1.0 - g_hi_in)) * mt,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )
    g_hi_out = 0.5 * 1.015
    assert (
        'gamma'
        in impact_loss(
            'roche2026',
            v_c=vc,
            M_i=(g_hi_out / (1.0 - g_hi_out)) * mt,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )

    # b: (0.0, 0.9)
    assert (
        'b'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=0.9 * 1.005, f_atm=fa
        ).flags
    )
    assert (
        'b'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=0.9 * 1.015, f_atm=fa
        ).flags
    )

    # R_ratio: (0.001, 1.0) per Ruling 13
    assert (
        'R_ratio'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=rt * 0.001 * 0.995, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'R_ratio'
        in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=rt * 0.001 * 0.985, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'R_ratio'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=rt * 1.0, R_t=rt, b=b, f_atm=fa
        ).flags
    )
    assert (
        'R_ratio'
        not in impact_loss(
            'roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=rt * 1.005, R_t=rt, b=b, f_atm=fa
        ).flags
    )

    # Ruling 13 test: R_i/R_t = 1.2 sets flag and X equals value without flag logic
    r_12 = impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=rt * 1.2, R_t=rt, b=b, f_atm=fa)
    assert 'R_ratio' in r_12.flags
    _, _, _, x_direct = _roche2026_fit(
        b=b,
        gamma=0.25,
        v_c_v_esc=r_12.diagnostics['v_ratio'],
        M_t_earth=1.0,
        M_i_earth=mi / Me,
        Q_R_prime_MJ=r_12.diagnostics['Q_R_prime'],
        f_atm=fa,
        R_ratio=1.2,
    )
    assert r_12.fraction == pytest.approx(x_direct, abs=1e-12)

    # v_ratio: (1.0, 3.0), sub-escape speeds unflagged
    vesc_curr = mutual_escape_speed(mt / (1.0 - fa), mi, rt, ri)
    assert (
        'v_ratio'
        not in impact_loss(
            'roche2026',
            v_c=3.0 * 1.005 * vesc_curr,
            M_i=mi,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )
    assert (
        'v_ratio'
        in impact_loss(
            'roche2026',
            v_c=3.0 * 1.015 * vesc_curr,
            M_i=mi,
            M_t=mt,
            R_i=ri,
            R_t=rt,
            b=b,
            f_atm=fa,
        ).flags
    )


def test_roche2026_validation_contracts(monkeypatch):
    """Verify input validation on mass_loss_roche2026 and impact_loss."""
    import zephyrus.collision as zc

    vc = 2.0e4
    mi = 1.0e24
    mt = Me
    ri = 3.0e6
    rt = 6.371e6
    b = 0.5
    fa = 0.05

    for entry_point in (
        lambda **kw: mass_loss_roche2026(
            v_c=kw.get('v_c', vc),
            M_i=kw.get('M_i', mi),
            M_t=kw.get('M_t', mt),
            R_i=kw.get('R_i', ri),
            R_t=kw.get('R_t', rt),
            b=kw.get('b', b),
            f_atm=kw.get('f_atm', fa),
        ),
        lambda **kw: impact_loss(
            'roche2026',
            v_c=kw.get('v_c', vc),
            M_i=kw.get('M_i', mi),
            M_t=kw.get('M_t', mt),
            R_i=kw.get('R_i', ri),
            R_t=kw.get('R_t', rt),
            b=kw.get('b', b),
            f_atm=kw.get('f_atm', fa),
        ),
    ):
        for bad_m in (0.0, -1.0, np.nan, np.inf):
            with pytest.raises(ValueError, match='M_i must be strictly positive'):
                entry_point(M_i=bad_m)
            with pytest.raises(ValueError, match='M_t must be strictly positive'):
                entry_point(M_t=bad_m)

        for bad_r in (0.0, -1.0, np.nan, np.inf):
            with pytest.raises(ValueError, match='R_i must be strictly positive'):
                entry_point(R_i=bad_r)
            with pytest.raises(ValueError, match='R_t must be strictly positive'):
                entry_point(R_t=bad_r)

        for bad_vc in (-1.0, np.nan, np.inf):
            with pytest.raises(ValueError, match='Collision speed v_c'):
                entry_point(v_c=bad_vc)

        with pytest.raises(TypeError, match='must be a scalar'):
            entry_point(v_c=np.array([1.0, 2.0]))

    # Non-finite scaling law output raises ValueError
    monkeypatch.setattr(zc, '_roche2026_fit', lambda **kw: (np.nan, np.nan, np.nan, np.nan))
    with pytest.raises(ValueError, match='produced non-finite result'):
        impact_loss('roche2026', v_c=vc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=fa)


def test_roche2026_float32_conversion():
    """Verify float32 scalar inputs match float64 result exactly."""
    v_c_f32 = np.float32(2.0e4)
    m_i_f32 = np.float32(1.0e24)
    m_t_f32 = np.float32(5.9722e24)
    r_i_f32 = np.float32(3.0e6)
    r_t_f32 = np.float32(6.371e6)
    b_f32 = np.float32(0.5)
    fa_f32 = np.float32(0.05)

    x32 = mass_loss_roche2026(v_c_f32, m_i_f32, m_t_f32, r_i_f32, r_t_f32, b_f32, fa_f32)
    x64 = mass_loss_roche2026(
        float(v_c_f32),
        float(m_i_f32),
        float(m_t_f32),
        float(r_i_f32),
        float(r_t_f32),
        float(b_f32),
        float(fa_f32),
    )
    assert np.isfinite(x32)
    assert 0.0 <= x32 <= 1.0
    assert x32 == x64


def test_roche2026_broadcast_shape():
    """Verify _roche2026_fit supports non-trivial broadcasting shapes."""
    b_arr = np.linspace(0.1, 0.8, 3)[:, None]
    fa_arr = np.linspace(0.02, 0.15, 4)[None, :]
    fnf, xnf, xff, xatm = _roche2026_fit(
        b=b_arr,
        gamma=0.3,
        v_c_v_esc=1.5,
        M_t_earth=1.0,
        M_i_earth=0.3 / 0.7,
        Q_R_prime_MJ=10.0,
        f_atm=fa_arr,
        R_ratio=0.8,
    )
    assert fnf.shape == (3, 4)
    assert xnf.shape == (3, 4)
    assert xff.shape == (3, 4)
    assert xatm.shape == (3, 4)


def test_roche2026_raw_fnf_greater_than_one():
    """Verify raw f_NF > 1 inside fitted box clamps to 1.0."""
    b = 0.1625
    mt = 1.22 * Me
    fa = 0.1086
    rr = 1.474
    rt = 6.371e6 * (1.22) ** (1.0 / 3.0)
    ri = rt * rr
    mi = 0.3 * mt
    vesc = mutual_escape_speed(mt / (1.0 - fa), mi, rt, ri)
    res = impact_loss(
        'roche2026', v_c=1.5 * vesc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=b, f_atm=fa
    )
    assert res.diagnostics['f_NF'] == pytest.approx(1.0, abs=1e-15)
    assert np.isfinite(res.fraction)
    assert 0.0 <= res.fraction <= 1.0


def test_roche2026_grazing_psi4_negative():
    """Verify grazing b = 1 with psi_4 < 0 returns finite X_NF and zero X_FF."""
    mt = 0.35 * Me
    mi = mt
    fa = 0.2
    rt = 6.371e6 * (0.35) ** (1.0 / 3.0)
    ri = rt
    vesc = mutual_escape_speed(mt / (1.0 - fa), mi, rt, ri)
    res = impact_loss(
        'roche2026', v_c=1.5 * vesc, M_i=mi, M_t=mt, R_i=ri, R_t=rt, b=1.0, f_atm=fa
    )
    assert np.isfinite(res.fraction)
    assert res.fraction == res.diagnostics['X_NF']
    assert res.diagnostics['X_FF'] == 0.0
