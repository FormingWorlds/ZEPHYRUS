"""
!!! info "`dispatcher.py`"
    The escape-regime dispatcher: one call, one regime, one rate.<br>
    Author: Mara Attia
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from zephyrus import boiloff as bl
from zephyrus import diagnostics as dg
from zephyrus import hydrodynamic as hy
from zephyrus import hydrostatic as hs
from zephyrus import knudsen as kn
from zephyrus import nozzle as nz
from zephyrus import thermostat as th
from zephyrus.composition import atomize, mean_particle_mass
from zephyrus.constants import kb, m_p
from zephyrus.fractionation import closure_per_species, unfractionated_split
from zephyrus.profiles import Profile, photospheric_level, wind_base_level

# The dispatcher assembles the escape branches of this package into one
# total prescription: every physically posed input state returns exactly
# one regime label, one bulk mass-loss rate, per-species rates summing to
# it, flags, and a diagnostics container. Exceptions are reserved for
# malformed input. The fixed evaluation order:
#
# 1. The bolometrically driven candidate is computed at every call. When
#    the restricted Jeans parameter sits below its threshold the atmosphere
#    is inflated enough to boil off and that candidate is the rate
#    (Owen & Wu 2016); past the threshold the same machinery survives as a
#    luminosity-capped residual (Gupta & Schlichting 2019) that can still
#    win the final comparison. XUV-driven escape needs a base to launch
#    from, and a bolometrically boiling atmosphere has not built one yet,
#    which is why this test precedes everything (Owen & Schlichting 2024).
# 2. The hydrodynamic candidate: the wind base is located by the
#    configured method, the thermostat sets the wind temperature by local
#    heating-cooling balance, and the candidate is min(EL, RR) with the
#    winner naming the sub-label.
# 3. The sonic-point Knudsen switch decides whether that wind is
#    collisional enough to exist. It lives on the hydrodynamic branch
#    only, never above the boil-off test. A confirmed hydrodynamic label
#    applies the fractionation closure; otherwise the point re-routes to
#    the hydrostatic branch.
# 4. The hydrostatic branch evaluates per-species Jeans escape with the
#    diffusion-limited supply on the extended upper structure. Its
#    escape-temperature gate re-routes thermally unstable exospheres back
#    to the hydrodynamic rate; points where the neutral and plasma gate
#    conventions disagree are flagged contested with both rates recorded.
# 5. The Roche screen, per branch, tests the active flow radius (sonic
#    radius, max(R_XUV, R_s), or exobase radius) against the periapsis
#    Hill radius before the label is finalized. An overflowing point is
#    renamed ``roche_overflow`` and keeps the rate its own branch
#    computed: the screen renames a state and never changes its rate.
#    Near misses raise ``near_roche``.
# 6. The final rate is the largest of the surviving branch rate, the
#    bolometric residual, and the tidally driven L1 nozzle rate
#    (Jackson et al. 2017), labeled by the winner. A nozzle win labels
#    ``roche_overflow`` with a real transfer rate, so that boundary is a
#    rate crossing and the dispatched rate is continuous across it; the
#    step 5 rename keeps its bound-flow lower-limit meaning, and
#    ``diagnostics['roche']['rate_branch']`` says which reading applies.
#
# Diagnostics are boxed: nothing in this module branches on anything the
# diagnostics container carries, and the container has no off switch.

REGIME_LABELS = (
    'boiloff',
    'hydrodynamic:EL',
    'hydrodynamic:RR',
    'hydrostatic',
    'roche_overflow',
)
# The label 'impact' is reserved for impact-driven escape (see collision.py
# for the erosion scaling law; the dispatcher does not yet route to it).

_TINY = 1e-300


@dataclass
class DispatchSettings:
    """Dispatch options; every default is the documented reference choice."""

    base_method: str = 'lopez'  # 'lopez' | 'fixed_pressure' | 'boreas'
    base_out_of_range: str = 'clamp'  # 'clamp' | 'extend'
    P_photo: float = 2000.0  # Pa; photospheric level for the EL geometry
    P_base_fixed: float = 5.0  # Pa; only for base_method = 'fixed_pressure'
    kn_crit: float = 1.0  # sonic-point Knudsen threshold
    kn_hysteresis: float = 1.5  # window factor, consumed only with prev_regime
    gate: str = 'neutral'  # 'neutral' | 'plasma' hydrostatic gate convention
    efficiency: float = 0.1
    efficiency_mode: str = 'fixed'  # 'fixed' | 'caldiroli'
    T_exo_mode: str = 'prescribed'  # 'prescribed' | 'thermostat'
    T_exo_value: float = 1000.0  # K; the prescribed exobase temperature
    cool_atomic: bool = True
    cool_co2_band: bool = True
    cool_o_finestructure: bool = True
    cool_recombination: bool = True
    fractionate: bool = True
    tidal: bool = True
    nozzle_temperature: str = 'photospheric'  # 'photospheric' | 'wind'
    lambda_crit: float = 20.0  # boil-off activation threshold (band 15 to 35)
    gamma_bates: float = 0.75  # Bates profile shape parameter
    kzz: float = 3.0e2  # m^2/s eddy diffusion when the profile carries none
    gamma_wind: float = 1.0  # polytropic index at the sonic point (isothermal)
    hydrostatic_levels_min: int = 200  # first quadrature grid of the supply integrals
    hydrostatic_levels_max: int = 3200  # refinement ceiling
    hydrostatic_rtol: float = 1.0e-2  # target relative change in the bulk rate

    def validate(self) -> None:
        """Raise ``ValueError`` on an unsupported option combination."""
        if self.base_method not in ('lopez', 'fixed_pressure', 'boreas'):
            raise ValueError("base_method must be 'lopez', 'fixed_pressure', or 'boreas'")
        if self.base_out_of_range not in ('clamp', 'extend'):
            raise ValueError("base_out_of_range must be 'clamp' or 'extend'")
        if self.gate not in ('neutral', 'plasma'):
            raise ValueError("gate must be 'neutral' or 'plasma'")
        if self.efficiency_mode not in ('fixed', 'caldiroli'):
            raise ValueError("efficiency_mode must be 'fixed' or 'caldiroli'")
        if self.T_exo_mode not in ('prescribed', 'thermostat'):
            raise ValueError("T_exo_mode must be 'prescribed' or 'thermostat'")
        if self.nozzle_temperature not in ('photospheric', 'wind'):
            raise ValueError("nozzle_temperature must be 'photospheric' or 'wind'")
        if not (
            self.cool_atomic
            or self.cool_co2_band
            or self.cool_o_finestructure
            or self.cool_recombination
        ):
            raise ValueError('all cooling channels disabled; at least one must stay on')
        # Numeric bounds. Outside them the closed forms leave their domains,
        # and what a caller saw was a bare math domain error from inside the
        # branch or, worse, a silently different regime label.
        for name, value in (
            ('P_photo', self.P_photo),
            ('P_base_fixed', self.P_base_fixed),
            ('kn_crit', self.kn_crit),
            ('T_exo_value', self.T_exo_value),
            ('lambda_crit', self.lambda_crit),
            ('kzz', self.kzz),
            ('gamma_bates', self.gamma_bates),
            ('hydrostatic_rtol', self.hydrostatic_rtol),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f'{name} must be a positive finite value, got {value!r}')
        if self.kn_hysteresis < 1.0 or not math.isfinite(self.kn_hysteresis):
            raise ValueError(
                'kn_hysteresis is a window factor at or above 1 (1 disables the '
                f'window), got {self.kn_hysteresis!r}'
            )
        if not 0.0 < self.efficiency <= 1.0:
            raise ValueError(
                f'efficiency is a fraction of the deposited power, got {self.efficiency!r}'
            )
        # The sonic-point scale height of Chatterjee & Pierrehumbert Eq. (17)
        # carries sqrt(5 - 3 gamma), which leaves the reals above the monatomic
        # 5/3. Below 1 the polytrope is no longer a wind solution.
        if (
            self.hydrostatic_levels_min < 2
            or self.hydrostatic_levels_max < self.hydrostatic_levels_min
        ):
            raise ValueError(
                'hydrostatic_levels_min must be at least 2 and no greater than '
                f'hydrostatic_levels_max, got {self.hydrostatic_levels_min!r} and '
                f'{self.hydrostatic_levels_max!r}'
            )
        if not 1.0 <= self.gamma_wind <= 5.0 / 3.0:
            raise ValueError(
                'gamma_wind must lie in [1, 5/3], the domain of the sonic-point '
                f'scale height, got {self.gamma_wind!r}'
            )


@dataclass
class EscapeInputs:
    """One dispatch call's physical state. SI at every boundary."""

    M_p: float  # kg, planet (interior) mass
    R_p: float  # m, planet (interior) radius
    M_star: float  # kg
    a: float  # m, semi-major axis
    e: float  # eccentricity
    T_eq: float  # K, equilibrium temperature
    F_xuv: float  # W m^-2
    F_bol: float  # W m^-2, bolometric instellation (consumed only by the
    #               nozzle power diagnostic; no branch rate reads it)
    F_int: float  # W m^-2, interior heat flux (luminosity cap)
    kappa_photo: float  # m^2 kg^-1, photospheric opacity
    profile: Profile
    settings: DispatchSettings = field(default_factory=DispatchSettings)
    prev_regime: str | None = None  # hysteresis memory for evolutionary use
    atm_converged: bool | None = None  # data-quality passthrough
    age: float | None = None  # s; consumed only by the snapshot screen
    reservoirs: dict | None = None  # element -> kg; screens and the split
    dt: float | None = None  # s; carried for the caller's supply cap only

    def validate(self) -> None:
        """Raise ``ValueError`` on a malformed physical state."""
        for name in ('M_p', 'R_p', 'M_star', 'a', 'T_eq', 'F_bol', 'F_int', 'kappa_photo'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive')
        if self.F_xuv < 0:
            raise ValueError('F_xuv must be >= 0')
        if not (0.0 <= self.e < 1.0):
            raise ValueError('e must be in [0, 1)')
        self.settings.validate()
        self.profile.validate()


@dataclass
class EscapeResult:
    """One dispatch call's outcome."""

    regime: str  # one of REGIME_LABELS
    mdot: float  # kg/s bulk rate, >= 0
    per_species: dict  # element -> kg/s, non-negative, summing to mdot
    flags: dict  # dispatch flags (clamps, screens, fallbacks)
    diagnostics: dict  # boxed reporting container; never gates anything


def dispatch(inputs: EscapeInputs) -> EscapeResult:
    """Dispatch one atmospheric state to its escape regime and rate.

    Runs the fixed evaluation order documented in the module notes and
    returns an :class:`EscapeResult`. Raises only on malformed input;
    every physically posed state returns a labeled, finite, non-negative
    result whose per-species rates sum to the bulk rate.
    """
    inputs.validate()
    st = inputs.settings
    flags: dict = {}
    diag: dict = {
        'documentation': {
            'murray_clay_exponents': dg.MURRAY_CLAY_EXPONENTS,
            'dayside_factors': dg.DAYSIDE_FACTORS,
            'kn_band': kn.KN_BAND,
            'lambda_crit_band': bl.LAMBDA_BAND,
        }
    }
    if inputs.atm_converged is False:
        flags['stale_input'] = True

    r_hill = hy.hill_radius_periapsis(inputs.M_p, inputs.M_star, inputs.a, inputs.e)
    photo, f = photospheric_level(inputs.profile, st.P_photo)
    flags.update(f)
    r_xuv = photo['r']

    # The tidal reduction of the escape barrier, needed by both rate
    # branches that measure one: the energy-limited rate and the
    # luminosity cap on the bolometric residual.
    xi_ktide = r_hill / inputs.R_p
    # The tidal factor has a double root at xi = 1 and the rates divide by it,
    # so it inflates them steeply as the lobe closes: 83-fold at xi = 1.1 and
    # 6.7e5-fold at xi = 1.001. At and below the root the barrier is gone and
    # the factor is undefined, so the rates are computed without it, which is
    # the smaller of the two readings. Such a state is already relabeled by the
    # Roche screen below; the flag says the reduction was dropped rather than
    # applied, and the inflation the factor is contributing is reported beside
    # the rate at every geometry so that a rate set by the divergence rather
    # than by the physics is visible as such.
    if st.tidal and xi_ktide <= 1.0:
        k_factor = 1.0
        flags['k_tide_undefined'] = True
    elif st.tidal:
        k_factor = hy.k_tide(xi_ktide)
    else:
        k_factor = 1.0

    # Step 1: bolometric candidate, computed at every point.
    lam_gate = bl.lambda_restricted(inputs.M_p, inputs.R_p, inputs.T_eq, photo['mmw'])
    bolo_rate, bolo = bl.bolometric_candidate(
        inputs.M_p,
        inputs.R_p,
        inputs.T_eq,
        inputs.kappa_photo,
        photo,
        inputs.F_int,
        lam_gate,
        st.lambda_crit,
        k_tide=k_factor,
    )
    flags.update(bolo['flags'])
    diag['lambda_gate'] = lam_gate
    diag['bolometric'] = {k: v for k, v in bolo.items() if k != 'flags'}
    diag['bolometric']['rate_kg_s'] = bolo_rate

    # Step 2: hydrodynamic candidate (always computed; it is cheap).
    channels = dict(
        cool_atomic=st.cool_atomic,
        cool_co2_band=st.cool_co2_band,
        cool_o_finestructure=st.cool_o_finestructure,
        cool_recombination=st.cool_recombination,
    )
    # The exobase temperature is resolved once and used by both the upper
    # structure the hydrostatic branch stands on and, under the extend
    # policy, the one the wind base is re-evaluated on. Resolving it twice
    # built those two structures at two different temperatures.
    t_exo = _resolve_t_exo(inputs, channels)
    base, f = _resolve_wind_base(inputs, t_exo)
    flags.update(f)
    elements = atomize(base['vmr'])
    t_wind, thermo = th.solve_wind_temperature(
        inputs.T_eq, base, elements, inputs.F_xuv, **channels
    )
    # Warnings about the hydrodynamic candidates are held aside and merged
    # only if one of them wins the route. A warning about the wind
    # temperature or the sonic radius describes a rate that a bolometric or
    # hydrostatic verdict did not dispatch, and the flags dictionary is read
    # as a warning set about the result. What the losing candidate did is
    # still in diag['hydrodynamic'].
    hydro_flags: dict = {}
    if thermo.get('clamped'):
        hydro_flags['thermostat_clamped'] = thermo['clamped']
    rr = hy.rr_chain(inputs.M_p, inputs.F_xuv, base['r'], t_wind, elements)
    if rr['subcritical']:
        hydro_flags['subcritical_sonic'] = True

    eps = st.efficiency
    if st.efficiency_mode == 'caldiroli':
        eta_eff, cf = hy.caldiroli_efficiency(inputs.F_xuv, inputs.M_p, inputs.R_p, k_factor)
        hydro_flags.update(cf)
        if eta_eff is not None:
            # Their efficiency is defined against an R_p^3 rate geometry.
            eps = eta_eff * (inputs.R_p / r_xuv) ** 2
        else:
            hydro_flags['efficiency_fallback_fixed'] = True
    mdot_el = hy.el_rate(eps, inputs.F_xuv, inputs.R_p, r_xuv, inputs.M_p, k_factor)
    mdot_rr = rr['mdot_rr']
    el_won = mdot_el <= mdot_rr
    mdot_hydro = min(mdot_el, mdot_rr)
    hydro_label = 'hydrodynamic:EL' if el_won else 'hydrodynamic:RR'
    diag['hydrodynamic'] = dict(
        mdot_el=mdot_el,
        mdot_rr=mdot_rr,
        efficiency=eps,
        K_tide=k_factor,
        T_wind=t_wind,
        selection_mechanism=hy.selection_mechanism(rr, el_won),
        rr_chain={
            k: rr[k]
            for k in (
                'c_s',
                'R_s',
                'R_s_calc',
                'lambda_b',
                'rho_s',
                'f_plus_base',
                'barometric_factor',
                'mu_wind',
                'mu_plus_wind',
                'hnu0_eV',
            )
        },
    )
    diag['thermostat'] = thermo

    # The tidally driven L1 nozzle candidate (Jackson et al. 2017 Eq. 3),
    # computed at every point: it joins the final comparison on both sides
    # of the activation gate, and its power comparison is an always-on
    # diagnostic. The temperature setting decides which state the flow is
    # launched from: the photospheric level (the primary's own
    # construction, a bolometrically maintained flow) or the wind base at
    # the thermostat's wind state (the upper envelope their Figure 9
    # explores). Both settings launch from one level with one temperature
    # and one mean mass, which is what the Bernoulli cancellation behind
    # the launch-level convention requires; the wind setting rebuilds the
    # launch density from the ideal gas law at the base pressure rather
    # than carrying the photosphere's cold density into a hot sound speed.
    # The radius is still the profile's, so the hot structure is not
    # solved, only its thermodynamic state, and that is a stated limit.
    # The flow is uncapped, faithful to the primary; the lift power
    # reported beside the interior and intercepted stellar luminosities
    # shows where that assumption is strained.
    if st.nozzle_temperature == 'wind':
        t_nozzle, mu_nozzle = t_wind, rr['mu_wind'] * m_p
    else:
        t_nozzle, mu_nozzle = photo['T'], photo['mmw']
    nozzle_full, noz = nz.nozzle_candidate(
        inputs.M_p,
        inputs.M_star,
        inputs.a,
        inputs.e,
        rho_ph=photo['rho'],
        r_ph=photo['r'],
        T=t_nozzle,
        mu_kg=mu_nozzle,
    )
    # The dispatched candidate is the average duty-cycled over the arc
    # where the overflow description applies; the unguarded average is
    # kept beside it so the closed form stays comparable with the
    # primary's published rates. The applicability edge is a criterion
    # boundary like the activation gate, not a rate crossing, and the jump
    # across it is a result to measure rather than hide.
    nozzle_rate = noz['rate_applicable_kg_s']
    nozzle_applicable = noz['applicable']
    noz['rate_kg_s'] = nozzle_rate
    noz['rate_full_orbit_kg_s'] = nozzle_full
    noz['temperature_mode'] = st.nozzle_temperature
    # The power comparison: what the isothermal flow demands against what
    # the planet has. Built from the barrier the rate applied plus the
    # acceleration to the sonic speed, so it stays finite and meaningful
    # at saturation, where the barrier is gone and the acceleration is not.
    noz['L_int_W'] = 4.0 * math.pi * inputs.R_p**2 * inputs.F_int
    noz['L_bol_intercepted_W'] = math.pi * inputs.R_p**2 * inputs.F_bol
    diag['nozzle'] = dict(noz)

    # Step 3: the sonic-point Knudsen switch.
    n_sc = rr['rho_s'] / (rr['mu_plus_wind'] * m_p)  # heavy-particle density
    if n_sc > _TINY:
        kn_sc, sigma_c, sigma_prov = kn.kn_sonic(
            n_sc, rr['R_s'], elements, t_wind, gamma=st.gamma_wind
        )
    else:
        kn_sc, sigma_c, sigma_prov = math.inf, math.nan, {}
    threshold = kn.effective_threshold(st.kn_crit, st.kn_hysteresis, inputs.prev_regime)
    diag['knudsen'] = dict(
        kn_sc=kn_sc,
        threshold_applied=threshold,
        sigma_c=sigma_c,
        provenance=sigma_prov,
        counterfactual_labels={
            edge: ('hydrodynamic' if kn_sc <= edge else 'hydrostatic') for edge in kn.KN_BAND
        },
    )
    if inputs.prev_regime is not None:
        flags['hysteresis_active'] = True

    # Step 4: hydrostatic branch (always evaluated: its exobase quantities
    # feed the diagnostics at every dispatch).
    hs_per_element, hsd = hs.hydrostatic_rates_refined(
        inputs.profile,
        inputs.M_p,
        t_exo,
        gamma_bates=st.gamma_bates,
        kzz_default=st.kzz,
        n_levels_min=st.hydrostatic_levels_min,
        n_levels_max=st.hydrostatic_levels_max,
        rtol=st.hydrostatic_rtol,
    )
    hs_flags = hsd.pop('flags')
    mdot_hs = sum(hs_per_element.values())
    hnu, _e_ion, sigma_front = th.front_constants(elements)
    f_plus_exo = th.ionization_fraction(
        hsd['n_exo'] * 1e-6,
        hsd['T_exo'],
        inputs.F_xuv * 1e3,
        hnu,
        sigma_front,
        th.recombination_alpha(elements, hsd['T_exo']),
    )
    unstable, contested = hs.gate_unstable(hsd['T_exo'], hsd, st.gate, f_plus_exo)
    diag['hydrostatic'] = dict(
        rate_kg_s=mdot_hs,
        convergence=hsd['convergence'],
        T_exo=hsd['T_exo'],
        T_exo_mode=st.T_exo_mode,
        r_exo=hsd['r_exo'],
        f_plus_exo=f_plus_exo,
        T_esc_neutral=hsd['T_esc_neutral'],
        T_esc_plasma=hsd['T_esc_plasma'],
        gate=st.gate,
        gate_unstable=unstable,
        detail=hsd,
    )
    if contested:
        flags['contested_ion'] = True
        diag['contested_ion'] = dict(
            hydrostatic_rate=mdot_hs,
            hydrodynamic_rate=mdot_hydro,
            note=(
                'the neutral and plasma escape-temperature conventions disagree '
                'here, so the branch assignment depends on ion physics this '
                'version does not model; both rates are recorded'
            ),
        )

    # Route. ``branch`` names the physics that produced the rate and decides
    # the split; ``label`` is what the caller reads back and the Roche screen
    # can overwrite it. The two are the same on every state whose flow stays
    # inside the Hill sphere.
    per_species = None
    if lam_gate < st.lambda_crit:
        branch = 'boiloff'
        rate = bolo_rate
        flow_radius = bolo['R_sonic']
    else:
        if kn_sc <= threshold:
            branch = hydro_label
            rate = mdot_hydro
            flow_radius = max(r_xuv, rr['R_s'])
            flags.update(hydro_flags)
        elif unstable:
            branch = hydro_label
            rate = mdot_hydro
            flow_radius = max(r_xuv, rr['R_s'])
            flags['gate_rerouted'] = True
            flags.update(hydro_flags)
        else:
            branch = 'hydrostatic'
            rate = mdot_hs
            per_species = dict(hs_per_element)
            flags.update(hs_flags)
            flow_radius = hsd['r_exo']
        # Step 5: the bolometric residual stays a candidate past the gate.
        if bolo_rate > rate:
            branch = 'boiloff'
            rate = bolo_rate
            per_species = None
            flags['bolometric_residual'] = True
            flow_radius = bolo['R_sonic']
            # The residual displaces whichever candidate had won, so the
            # warnings about that candidate stop describing the result.
            for key in tuple(hydro_flags) + ('hydrostatic_lower_limit',):
                flags.pop(key, None)
    # The nozzle candidate competes last, on both sides of the activation
    # gate, wherever the overflow description applies: where the
    # photosphere approaches the lobe, the tidally driven transfer through
    # L1 outruns every bound-flow estimate, and the label boundary it
    # creates is a rate crossing, continuous by construction. A candidate
    # below the one-proton-per-Julian-year floor does not compete: this
    # label is a rate crossing, not a geometric verdict, and a crossing
    # between two numerically empty numbers would rename the deeply bound
    # corner on no physical content (the floor otherwise stays reported
    # and never applied, and a geometric verdict still ignores it).
    if nozzle_applicable and nozzle_rate > rate and nozzle_rate > dg.RATE_FLOOR_KG_S:
        branch = 'roche_nozzle'
        rate = nozzle_rate
        per_species = None
        flow_radius = noz['r_lobe']
        for key in tuple(hydro_flags) + ('hydrostatic_lower_limit', 'bolometric_residual'):
            flags.pop(key, None)
        if noz['saturated']:
            # The photospheric potential reached the L1 value, so the rate
            # is the lobe-filling boundary value of the model rather than
            # an interior point of it.
            flags['nozzle_saturated'] = True
        if inputs.e > 0.0:
            # The rate is a time average over the orbit, evaluated with
            # the circular formula at each separation. Under saturation it
            # scales as the cube of the separation, so periapsis is a
            # lower bound there and an upper bound while the barrier is
            # unclamped; the average is what a secular caller needs either
            # way.
            flags['nozzle_orbit_averaged'] = True
        if noz['applicable_orbit_fraction'] < 1.0:
            # The overflow description holds only on an arc around
            # periapsis. The rate is duty-cycled over that arc, so it
            # omits the wind the planet drives on the rest of the orbit.
            flags['nozzle_partial_orbit'] = True
    label = 'roche_overflow' if branch == 'roche_nozzle' else branch

    # Step 6: the Roche screen on the active flow radius. The screen renames
    # the state and never touches the rate. Its boundary is a rate
    # comparison, since the branch whose flow radius gets tested is the one
    # that won the final comparison, so reporting the winning branch's own
    # rate keeps the dispatched rate continuous across the boundary;
    # substituting another branch's formula would not. When the rename fires
    # on a bound branch, the rate beside the label is the bound-flow
    # estimate, a lower limit on what tides would do; when the nozzle
    # candidate won above, the rate is the tidally driven transfer itself
    # and the subflag reads ``nozzle``.
    xi_flow = r_hill / flow_radius if flow_radius > 0 else math.inf
    # The outer extent of the atmosphere itself, modeled plus extended,
    # which is what separates the two overflow geometries. It is reported
    # and used for that separation, and deliberately not used to trigger
    # the screen: what the screen asks is whether the escaping flow stays
    # bound, and widening its trigger would move the label boundary itself.
    r_atm = max(float(inputs.profile.r[-1]), hsd['r_exo'])
    diag['roche'] = dict(
        R_hill_periapsis=r_hill,
        flow_radius=flow_radius,
        xi_flow=xi_flow,
        xi_ktide=xi_ktide,
        k_tide=k_factor,
        tidal_inflation=1.0 / k_factor,
        r_atmosphere=r_atm,
        rate_branch=branch,
    )
    if xi_flow <= 1.0 or xi_ktide <= 1.0:
        label = 'roche_overflow'
        flags['roche_overflow'] = True
        # Dynamical overflow when the atmosphere itself reaches the lobe;
        # no transonic solution when only the flow radius does, which is the
        # narrow band Owen & Jackson (2012) describe.
        flags['roche_subflag'] = (
            'dynamical' if (xi_ktide <= 1.0 or r_hill <= r_atm) else 'no_transonic'
        )
    elif branch == 'roche_nozzle':
        # The label arrived through the rate crossing rather than the
        # geometric trigger. The subflag still reads the geometry: an
        # atmosphere that itself reaches the lobe is dynamical overflow
        # whichever candidate carries the rate; ``nozzle`` marks the
        # remaining case, a photosphere close enough to the lobe for the
        # L1 transfer to outrun the bound branches while the structure
        # sits inside the Hill sphere. ``near_roche`` is a warning about
        # the tidal inflation of a bound rate, which this rate is not.
        flags['roche_overflow'] = True
        flags['roche_subflag'] = (
            'dynamical' if (xi_ktide <= 1.0 or r_hill <= r_atm) else 'nozzle'
        )
    elif xi_flow < 1.5:
        flags['near_roche'] = True

    # Per-species split, by the branch that produced the rate.
    if per_species is None:
        if branch.startswith('hydrodynamic') and st.fractionate and rate > 0.0:
            per_species, cdiag, cflags = closure_per_species(
                rate, elements, t_wind, inputs.M_p, base['r']
            )
            flags.update(cflags)
            diag['closure'] = cdiag
        else:
            per_species, sflags = unfractionated_split(rate, inputs.reservoirs, elements)
            flags.update(sflags)

    # Boxed diagnostics (nothing above this line reads them back).
    m_bar = mean_particle_mass(elements)
    ratio, q_net, q_c = dg.q_net_over_qc(
        eps,
        inputs.F_xuv,
        r_xuv,
        rr['R_s'],
        base['r'],
        inputs.M_p,
        m_bar,
        sigma_c if sigma_c == sigma_c else 1e-19,
    )
    diag['johnson_q'] = dict(q_net_over_qc=ratio, q_net_W=q_net, q_c_W=q_c)
    diag['guo_triple'] = dg.guo_triple(
        inputs.M_p,
        inputs.R_p,
        inputs.T_eq,
        photo['mmw'],
        inputs.M_star,
        inputs.a,
        inputs.e,
        hsd['lambda_exo_bulk'],
    )
    diag['potential_screens'] = dg.potential_screens(inputs.M_p, inputs.R_p)
    diag['erkaev_tc_K'] = dg.erkaev_tc(inputs.M_p, inputs.R_p, hsd['r_exo'], r_hill)
    diag['fluid_check'] = dg.along_profile_fluid_check(
        inputs.profile, inputs.M_p, rr['R_s'], st.kn_crit
    )
    diag['tang_timescale'] = bl.tang_timescale_check(
        inputs.M_p, inputs.R_p, inputs.F_int, rate, inputs.reservoirs
    )
    diag['self_consistency'] = dg.self_consistency_screen(inputs.reservoirs, rate, inputs.age)
    diag['rate_floor'] = dg.rate_floor_screen(rate)
    diag['base_level'] = dict(
        p_Pa=base['p'],
        p_physical_Pa=base.get('p_physical'),
        r_m=base['r'],
        T_K=base['T'],
        clamp_decades=flags.get('base_clamp_decades'),
    )

    # Output contract: non-negative bulk rate, per-species rates summing to
    # it exactly.
    rate = max(rate, 0.0)
    per_species = {el: max(v, 0.0) for el, v in per_species.items()}
    tot = sum(per_species.values())
    if tot > 0.0 and rate > 0.0 and abs(tot - rate) / rate > 1e-9:
        per_species = {el: v * rate / tot for el, v in per_species.items()}
    return EscapeResult(
        regime=label, mdot=rate, per_species=per_species, flags=flags, diagnostics=diag
    )


def _resolve_wind_base(inputs: EscapeInputs, t_exo: float) -> tuple[dict, dict]:
    """The wind-base level, and what happens when the profile cannot reach it.

    Locates the base by the configured method; when the physical base
    pressure lies above the profile top and the policy is ``'extend'``,
    the level is re-evaluated on the Bates extension at ``t_exo``, which is
    the exobase temperature the hydrostatic branch is given, so the two
    structures are one structure. Flagged ``base_extended``; under
    ``'clamp'`` (default) the clamped top level and its recorded clamp
    distance stand.
    """
    st = inputs.settings
    boreas_scalars = None
    if st.base_method == 'boreas':
        boreas_scalars = {'R_p': inputs.R_p, 'T_eq': inputs.T_eq, 'F_xuv': inputs.F_xuv}
    base, flags = wind_base_level(
        inputs.profile,
        inputs.M_p,
        method=st.base_method,
        fixed_pressure=st.P_base_fixed,
        boreas_scalars=boreas_scalars,
    )
    if flags.get('base_clamped') and st.base_out_of_range == 'extend':
        p_target = base.get('p_physical')
        if p_target is not None:
            ext = hs.bates_extension(inputs.profile, inputs.M_p, t_exo, gamma=st.gamma_bates)
            p_ext = np.asarray(ext['p'])
            if p_target >= p_ext[-1]:
                lz = np.log(p_ext[::-1])
                lt = math.log(p_target)
                r = float(np.interp(lt, lz, ext['r'][::-1]))
                T = float(np.interp(lt, lz, ext['T'][::-1]))
                n = p_target / (kb * T)
                base = dict(
                    p=float(p_target),
                    p_physical=float(p_target),
                    r=r,
                    T=T,
                    mmw=float(ext['mu']),
                    n=float(n),
                    rho=float(n * ext['mu']),
                    vmr=dict(ext['vmr']),
                    kzz=None,
                )
                flags.pop('base_clamped', None)
                flags.pop('base_clamp_decades', None)
                flags['base_extended'] = True
            else:
                # The extension itself truncates (unbound) above the target
                # level; the clamp to the profile top stands, flagged.
                flags['base_extension_truncated'] = True
    return base, flags


def _resolve_t_exo(inputs: EscapeInputs, channels: dict) -> float:
    """The hydrostatic exobase temperature per the configured mode.

    The default is the prescribed value: the exobase temperature is the
    branch's dominant sensitivity and is owned by the caller. The optional
    thermostat mode evaluates the local heating-cooling balance at the
    profile top; a conduction-free local balance biases the estimate high
    (heating scales with density, the cooling channels with its square),
    so the mode is an estimator, not a default.
    """
    st = inputs.settings
    if st.T_exo_mode != 'thermostat':
        return st.T_exo_value
    top = {
        'n': float(inputs.profile.p[-1]) / (kb * float(inputs.profile.T[-1])),
        'vmr': {sp: float(np.asarray(v)[-1]) for sp, v in inputs.profile.vmr.items()},
    }
    t_exo, _detail = th.solve_wind_temperature(
        inputs.T_eq, top, atomize(top['vmr']), inputs.F_xuv, **channels
    )
    return t_exo
