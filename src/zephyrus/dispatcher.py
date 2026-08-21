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
#    Hill radius before the label is finalized; an overflowing point is
#    labeled ``roche_overflow`` with the Bondi-capped bolometric rate at
#    the overflow geometry, and near misses raise ``near_roche``.
# 6. The final rate is the larger of the surviving branch rate and the
#    bolometric residual, labeled by the winner.
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
    lambda_crit: float = 20.0  # boil-off activation threshold (band 15 to 35)
    gamma_bates: float = 0.75  # Bates profile shape parameter
    kzz: float = 3.0e2  # m^2/s eddy diffusion when the profile carries none
    gamma_wind: float = 1.0  # polytropic index at the sonic point (isothermal)

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
        if not (
            self.cool_atomic
            or self.cool_co2_band
            or self.cool_o_finestructure
            or self.cool_recombination
        ):
            raise ValueError('all cooling channels disabled; at least one must stay on')


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
    F_bol: float  # W m^-2, bolometric instellation (carried with the state;
    #               not consumed by any branch in this version)
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
    )
    flags.update(bolo['flags'])
    diag['lambda_gate'] = lam_gate
    diag['bolometric'] = {k: v for k, v in bolo.items() if k != 'flags'}
    diag['bolometric']['rate_kg_s'] = bolo_rate

    # Step 2: hydrodynamic candidate (always computed; it is cheap).
    base, f = _resolve_wind_base(inputs)
    flags.update(f)
    elements = atomize(base['vmr'])
    channels = dict(
        cool_atomic=st.cool_atomic,
        cool_co2_band=st.cool_co2_band,
        cool_o_finestructure=st.cool_o_finestructure,
        cool_recombination=st.cool_recombination,
    )
    t_wind, thermo = th.solve_wind_temperature(
        inputs.T_eq, base, elements, inputs.F_xuv, **channels
    )
    if thermo.get('clamped'):
        flags['thermostat_clamped'] = thermo['clamped']
    rr = hy.rr_chain(inputs.M_p, inputs.F_xuv, base['r'], t_wind, elements)
    if rr['subcritical']:
        flags['subcritical_sonic'] = True

    xi_ktide = r_hill / inputs.R_p
    k_factor = hy.k_tide(xi_ktide) if (st.tidal and xi_ktide > 1.0) else 1.0
    eps = st.efficiency
    if st.efficiency_mode == 'caldiroli':
        eta_eff, cf = hy.caldiroli_efficiency(inputs.F_xuv, inputs.M_p, inputs.R_p, k_factor)
        flags.update(cf)
        if eta_eff is not None:
            # Their efficiency is defined against an R_p^3 rate geometry.
            eps = eta_eff * (inputs.R_p / r_xuv) ** 2
        else:
            flags['efficiency_fallback_fixed'] = True
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
    t_exo = _resolve_t_exo(inputs, channels)
    if st.T_exo_mode == 'thermostat':
        flags['T_exo_thermostat'] = True
    hs_per_element, hsd = hs.hydrostatic_rates(
        inputs.profile, inputs.M_p, t_exo, gamma_bates=st.gamma_bates, kzz_default=st.kzz
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
        T_exo=hsd['T_exo'],
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

    # Route.
    per_species = None
    if lam_gate < st.lambda_crit:
        label = 'boiloff'
        rate = bolo_rate
        flow_radius = bolo['R_sonic']
    else:
        if kn_sc <= threshold:
            label = hydro_label
            rate = mdot_hydro
            flow_radius = max(r_xuv, rr['R_s'])
        elif unstable:
            label = hydro_label
            rate = mdot_hydro
            flow_radius = max(r_xuv, rr['R_s'])
            flags['gate_rerouted'] = True
        else:
            label = 'hydrostatic'
            rate = mdot_hs
            per_species = dict(hs_per_element)
            flags.update(hs_flags)
            flow_radius = hsd['r_exo']
        # Step 6: the bolometric residual stays a candidate past the gate.
        if bolo_rate > rate:
            label = 'boiloff'
            rate = bolo_rate
            per_species = None
            flags['bolometric_residual'] = True
            flow_radius = bolo['R_sonic']

    # Step 5: the Roche screen on the active flow radius.
    xi_flow = r_hill / flow_radius if flow_radius > 0 else math.inf
    diag['roche'] = dict(
        R_hill_periapsis=r_hill, flow_radius=flow_radius, xi_flow=xi_flow, xi_ktide=xi_ktide
    )
    if xi_flow <= 1.0 or xi_ktide <= 1.0:
        label = 'roche_overflow'
        flags['roche_overflow'] = True
        flags['roche_subflag'] = (
            'dynamical' if (xi_ktide <= 1.0 or r_hill <= photo['r']) else 'no_transonic'
        )
        # The overflow rate is the Bondi-capped bolometric machinery at the
        # overflow geometry.
        rate = min(bolo['mdot_parker'], bolo['mdot_bondi'])
        per_species = None
    elif xi_flow < 1.5:
        flags['near_roche'] = True

    # Per-species split, by label.
    if per_species is None:
        if label.startswith('hydrodynamic') and st.fractionate and rate > 0.0:
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


def _resolve_wind_base(inputs: EscapeInputs) -> tuple[dict, dict]:
    """The wind-base level with the out-of-range policy applied.

    Locates the base by the configured method; when the physical base
    pressure lies above the profile top and the policy is ``'extend'``,
    the level is re-evaluated on the Bates extension (the same upper
    structure the hydrostatic branch uses), flagged ``base_extended``;
    under ``'clamp'`` (default) the clamped top level and its recorded
    clamp distance stand.
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
            t_exo = st.T_exo_value if st.T_exo_mode == 'prescribed' else inputs.T_eq
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
