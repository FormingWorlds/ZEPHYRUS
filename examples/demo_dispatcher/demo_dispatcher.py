"""Escape regimes of small planets, dispatched with `zephyrus.dispatch`.

Runs the regime framework over a set of synthetic atmospheres: one verdict
read field by field, a flux sweep that crosses two regime boundaries, the
boil-off and Roche-overflow labels, the diagnostics container, the four
knobs that move a boundary, and one planet dispatched along a stellar XUV
history. Prints a table per step and writes one figure.

Every function returns its results, so the steps can be imported and reused
one at a time; nothing runs on import. From the repository root:

    mkdir -p output && python examples/demo_dispatcher/demo_dispatcher.py

Companion tutorial: docs/Tutorials/dispatch.md.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import mors
import numpy as np

from zephyrus.composition import atomize, species_mass_amu
from zephyrus.constants import au2m
from zephyrus.diagnostics import RATE_FLOOR_KG_S
from zephyrus.dispatcher import DispatchSettings, EscapeInputs, dispatch
from zephyrus.planets_parameters import Ls, Me, Me_atm, Ms, Re
from zephyrus.profiles import isothermal_profile

# ----------------------------------------------------------------- setup

SIGMA_SB = 5.670374419e-8  # Stefan-Boltzmann constant             [W m-2 K-4]
SEC_PER_YR = 3.15576e7  # Julian year                              [s]

# The reference orbit places a solar-luminosity star's zero-albedo,
# full-redistribution equilibrium temperature at 1000 K, so the orbit, the
# temperature, and the bolometric flux below are mutually consistent.
A_REF = 0.0775 * au2m  # reference semi-major axis                 [m]
A_TRACK = 0.2 * au2m  # semi-major axis for the stellar history    [m]

KAPPA_PHOTO = 0.01  # photospheric opacity, about 0.1 cm2 g-1      [m2 kg-1]
F_INT = 1.0  # interior heat flux                                  [W m-2]
P_SURF = 1.0e7  # profile base pressure, 100 bar                   [Pa]
P_TOP = 1.0e-5  # profile top pressure, 0.1 nanobar                [Pa]
AGE_REF = 1.0e8 * SEC_PER_YR  # snapshot age, 100 Myr              [s]

# One proton crossing the surface per year: the smallest rate that can mean
# anything physically. A returned rate below this is numerical noise, not
# escape. It is a reporting convention for the caller, so the module computes
# whatever the physics gives it and reports the comparison beside every
# verdict, in diagnostics['rate_floor'], rather than applying it. The
# constant is imported from there so this file keeps no second copy.
RATE_FLOOR = RATE_FLOOR_KG_S  # [kg s-1]

COMPOSITIONS = {
    'CO2': {'CO2': 1.0},
    'N2-O2': {'N2': 0.8, 'O2': 0.2},
    'H/He': {'H2': 0.9, 'He': 0.1},
    'CO2 + 1% H2': {'CO2': 0.99, 'H2': 0.01},
}

SWEEP_FLUXES = np.logspace(-2, math.log10(5.0e3), 30)  # [W m-2]


def brand_colors() -> dict:
    """Apply the PROTEUS figure theme when it is installed, and return the
    colors this example draws with.

    The theme is an optional convenience: the hex values are the same either
    way, so the figure carries the same identity without it. Escape is the
    chemistry and escape domain, and the regime colors come from the shared
    categorical cycle, with the hottest regime taking the one red mark.
    """
    palette = {
        'escape': '#1B6FA8',  # domain color of chemistry and escape
        'boiloff': '#E23D28',  # the one hot mark on the figure
        'hydrodynamic:EL': '#1B6FA8',
        'hydrodynamic:RR': '#4FA3D9',
        'hydrostatic': '#7A8894',
        'roche_overflow': '#593E74',
        'solar': '#C8860F',
        'rule': '#3E4A55',
    }
    try:
        import proteus_mpl

        proteus_mpl.use()
        palette['escape'] = proteus_mpl.DOMAINS['chemistry']
        palette['hydrodynamic:EL'] = proteus_mpl.DOMAINS['chemistry']
        palette['hydrodynamic:RR'] = proteus_mpl.COLORS['azure']
        palette['hydrostatic'] = proteus_mpl.COLORS['fog']
        palette['boiloff'] = proteus_mpl.COLORS['magma']
        # The violet slot of the shared cycle, not the tidal module color:
        # these are regime categories, not any module's output.
        palette['roche_overflow'] = proteus_mpl.CYCLE[5]
        palette['solar'] = proteus_mpl.COLORS['solar_deep']
        palette['rule'] = proteus_mpl.COLORS['ink']
    except ImportError:
        pass
    return palette


# Marker per label, so the figure never encodes a regime in color alone.
REGIME_MARKERS = {
    'boiloff': 'D',
    'hydrodynamic:EL': 'o',
    'hydrodynamic:RR': 's',
    'hydrostatic': '^',
    'roche_overflow': 'X',
}


def t_eq_at(a_m: float, luminosity: float = Ls) -> float:
    """Equilibrium temperature [K] for zero albedo and full redistribution."""
    return (luminosity / (16.0 * math.pi * SIGMA_SB * a_m**2)) ** 0.25


def f_bol_at(a_m: float, luminosity: float = Ls) -> float:
    """Bolometric instellation [W m-2] at orbital distance `a_m` [m]."""
    return luminosity / (4.0 * math.pi * a_m**2)


def element_reservoirs(composition: dict, mass_total: float = Me_atm) -> dict:
    """Element inventories [kg] for an atmosphere of `mass_total` [kg].

    The dispatcher consumes reservoirs two ways: the snapshot consistency
    screen compares their total against the dispatched rate, and a
    non-fractionating branch splits its rate over them by mass.
    """
    elements = atomize(composition)
    norm = sum(x * species_mass_amu(el) for el, x in elements.items())
    return {el: mass_total * x * species_mass_amu(el) / norm for el, x in elements.items()}


def build_state(
    composition: str,
    m_earth: float,
    r_earth: float,
    f_xuv: float,
    a_m: float = A_REF,
    settings: DispatchSettings | None = None,
    reservoir_mass: float = Me_atm,
) -> EscapeInputs:
    """One dispatch input state: scalars plus an atmosphere profile.

    Unlike the energy-limited entry point, which takes scalars alone, the
    framework needs the atmospheric structure: the wind base, the exobase,
    and the collisionality of the flow are all properties of the profile,
    not of the surface. The profile here is isothermal at the equilibrium
    temperature, which is enough to exercise every branch; a coupled run
    supplies the atmosphere model's own profile instead.
    """
    comp = COMPOSITIONS[composition]
    m_p, r_p = m_earth * Me, r_earth * Re
    t_eq = t_eq_at(a_m)
    return EscapeInputs(
        M_p=m_p,
        R_p=r_p,
        M_star=Ms,
        a=a_m,
        e=0.0,
        T_eq=t_eq,
        F_xuv=f_xuv,
        F_bol=f_bol_at(a_m),
        F_int=F_INT,
        kappa_photo=KAPPA_PHOTO,
        # The top pressure sits below the nanobar level the XUV wind
        # launches from, so the base never clamps to the profile top.
        profile=isothermal_profile(m_p, r_p, t_eq, comp, P_SURF, P_TOP),
        settings=settings or DispatchSettings(),
        age=AGE_REF,
        reservoirs=element_reservoirs(comp, reservoir_mass),
    )


# ------------------------------------------------------- step 1: one call


def one_verdict(composition: str = 'CO2', f_xuv: float = 10.0) -> dict:
    """Dispatch one state and report every field of the result."""
    result = dispatch(build_state(composition, 1.0, 1.0, f_xuv))
    total = sum(result.per_species.values())
    print(f'\n=== One verdict: {composition}, 1 Me, 1 Re, F_xuv = {f_xuv:g} W m-2 ===')
    print(f'  regime      {result.regime}')
    print(f'  rate        {result.mdot:.4e} kg/s   ({result.mdot * SEC_PER_YR:.3e} kg/yr)')
    for element, rate in sorted(result.per_species.items()):
        print(f'  {element:<11} {rate:.4e} kg/s')
    print(f'  sum         {total:.4e} kg/s, closure error {abs(total - result.mdot):.2e}')
    print(f'  flags       {result.flags or "none"}')
    groups = ', '.join(sorted(result.diagnostics))
    print(f'  diagnostics {len(result.diagnostics)} groups: {groups}')
    return dict(result=result, closure_error=abs(total - result.mdot))


# ---------------------------------------------------- step 2: flux sweeps


def flux_sweep(composition: str, m_earth: float = 1.0, r_earth: float = 1.0) -> list[dict]:
    """Dispatch one planet across the XUV flux range and record each verdict.

    Sweeping the flux at fixed orbit stands in for stellar age: a young star
    delivers orders of magnitude more XUV than the same star does later.
    """
    rows = []
    for f_xuv in SWEEP_FLUXES:
        result = dispatch(build_state(composition, m_earth, r_earth, float(f_xuv)))
        knudsen = result.diagnostics['knudsen']
        rows.append(
            dict(
                F_xuv=float(f_xuv),
                regime=result.regime,
                mdot=result.mdot,
                above_floor=result.diagnostics['rate_floor']['above_floor'],
                kn_sc=knudsen['kn_sc'],
                counterfactual=knudsen['counterfactual_labels'],
                T_wind=result.diagnostics['hydrodynamic']['T_wind'],
                mechanism=result.diagnostics['hydrodynamic']['selection_mechanism'],
                flags=dict(result.flags),
            )
        )
    return rows


def boundary_flux(composition: str, lo: float, hi: float, tol: float = 1e-3) -> float:
    """Locate a label change in flux by bisection, to a tolerance in decades."""
    label_lo = dispatch(build_state(composition, 1.0, 1.0, lo)).regime
    while math.log10(hi / lo) > tol:
        mid = math.sqrt(lo * hi)
        if dispatch(build_state(composition, 1.0, 1.0, mid)).regime == label_lo:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def report_sweep(composition: str, rows: list[dict]) -> list[tuple]:
    """Print a sweep and the fluxes where its label changes."""
    print(f'\n=== Flux sweep: {composition}, 1 Me, 1 Re ===')
    print(
        f'  {"F_xuv":>9}  {"regime":<17} {"rate [kg/s]":>12}  {"Kn_sc":>9}  '
        f'{"label at Kn 0.1":<14} {"label at Kn 3":<14}'
    )
    for row in rows[::3]:
        counter = row['counterfactual']
        print(
            f'  {row["F_xuv"]:9.3g}  {row["regime"]:<17} {row["mdot"]:12.3e}  '
            f'{row["kn_sc"]:9.3g}  {counter[0.1]:<14} {counter[3.0]:<14}'
        )
    changes = []
    for before, after in zip(rows[:-1], rows[1:]):
        if before['regime'] != after['regime']:
            flux = boundary_flux(composition, before['F_xuv'], after['F_xuv'])
            changes.append((before['regime'], after['regime'], flux))
            print(f'  boundary: {before["regime"]} to {after["regime"]} at {flux:.4g} W m-2')
    raised = sorted({flag for row in rows for flag in row['flags']})
    print(f'  flags raised anywhere in the sweep: {raised or ["none"]}')
    below = [row['F_xuv'] for row in rows if not row['above_floor']]
    if below:
        print(
            f'  rates below the one proton per year floor ({RATE_FLOOR:.2g} kg/s) up to '
            f'F_xuv = {max(below):.3g} W m-2: no escape worth reporting there'
        )
    return changes


# ------------------------------------------- step 3: the other two labels


def extreme_labels() -> list[dict]:
    """The boil-off and Roche-overflow labels, one call each."""
    print('\n=== The other two labels ===')
    out = []
    cases = [
        ('boil-off: an inflated hydrogen envelope', 'H/He', 1.0, 1.5, 10.0),
        ('overflow: a puffy envelope filling its Hill sphere', 'H/He', 3.0, 2.0, 0.1),
    ]
    for note, composition, m_earth, r_earth, f_xuv in cases:
        result = dispatch(build_state(composition, m_earth, r_earth, f_xuv))
        diagnostics = result.diagnostics
        print(f'  {note}')
        print(
            f'    regime {result.regime}, rate {result.mdot:.3e} kg/s, '
            f'Lambda {diagnostics["lambda_gate"]:.3g} '
            f'(activation band {diagnostics["documentation"]["lambda_crit_band"]})'
        )
        roche = diagnostics['roche']
        print(
            f'    flow radius {roche["flow_radius"]:.3e} m against Hill radius '
            f'{roche["R_hill_periapsis"]:.3e} m, ratio {roche["xi_flow"]:.3f}'
        )
        print(
            f'    atmosphere reaches {roche["r_atmosphere"] / roche["R_hill_periapsis"]:.3f} '
            f"Hill radii; the rate is the {roche['rate_branch']} branch's"
        )
        print(f'    flags {sorted(result.flags)}')
        out.append(dict(case=note, regime=result.regime, mdot=result.mdot))
    return out


# ------------------------------------------------- step 4: the diagnostics


def read_diagnostics(composition: str = 'CO2', f_xuv: float = 10.0) -> dict:
    """Walk the diagnostics container of one verdict, group by group."""
    result = dispatch(build_state(composition, 1.0, 1.0, f_xuv))
    d = result.diagnostics
    print(f'\n=== Diagnostics: {composition}, F_xuv = {f_xuv:g} W m-2, {result.regime} ===')

    knudsen = d['knudsen']
    print('  Which branch, and how close was the switch')
    print(
        f'    Kn_sc {knudsen["kn_sc"]:.4g} against threshold '
        f'{knudsen["threshold_applied"]:.3g}; at the band edges the label would be '
        f'{knudsen["counterfactual_labels"]}'
    )

    hydro = d['hydrodynamic']
    print('  What set the rate')
    print(
        f'    energy limited {hydro["mdot_el"]:.3e}, recombination limited '
        f'{hydro["mdot_rr"]:.3e} kg/s, selected by {hydro["selection_mechanism"]}'
    )
    print(
        f'    wind temperature {hydro["T_wind"]:.0f} K, efficiency '
        f'{hydro["efficiency"]:.3g}, tidal factor {hydro["K_tide"]:.4f}'
    )

    johnson = d['johnson_q']
    print('  Could the heating drive the flow at all')
    print(
        f'    absorbed power over the critical power {johnson["q_net_over_qc"]:.3g} '
        f'(below 1 says no transonic flow, whatever a rate formula returns)'
    )

    guo, screens = d['guo_triple'], d['potential_screens']
    print('  Translations into other taxonomies')
    print(
        f'    Jeans parameters: exobase {guo["lambda_exo"]:.3g}, at the radius '
        f'{guo["lambda_rp"]:.3g}, tidally corrected {guo["lambda_star"]:.3g}'
    )
    print(
        f'    log potential {screens["log_minus_phi_cgs"]:.3f}; efficiency collapse band '
        f'{screens["caldiroli_threshold"]}, wind screen says {screens["salz_verdict"]}'
    )

    fluid, consistency = d['fluid_check'], d['self_consistency']
    print('  Is the snapshot self-consistent')
    print(
        f'    worst Knudsen below the sonic surface {fluid["worst_kn"]:.3g} over '
        f'{fluid["levels_checked"]} levels, fluid {fluid["fluid"]}, '
        f'truncated at the profile top {fluid["truncated_at_profile_top"]}'
    )
    if consistency['evaluated']:
        print(
            f'    reservoirs empty in {consistency["t_deplete_s"] / SEC_PER_YR:.3e} yr '
            f'against an age of {consistency["age_s"] / SEC_PER_YR:.3e} yr, '
            f'inconsistent {consistency["inconsistent"]}'
        )

    base = d['base_level']
    print('  Where the wind was launched')
    print(
        f'    base at {base["p_Pa"]:.3e} Pa and {base["r_m"]:.4e} m, physical target '
        f'{base["p_physical_Pa"]:.3e} Pa, clamp {base["clamp_decades"]}'
    )
    print(f'  Coefficient provenance: {knudsen["provenance"]}')
    return d


# ----------------------------------------------------------- step 5: knobs


def knob_collisionality(composition: str = 'CO2') -> dict:
    """Locate the wind boundary at each edge of the collisionality band.

    The threshold's physical band spans a factor of 30, because kinetic
    simulations place the fluid-to-kinetic transition near 0.1 for heating
    deposited in a sharp layer and near 1 for distributed heating. Bisecting
    the boundary at each edge turns that band into the quantity a reader
    needs: the range of fluxes over which the label is not decided.
    """
    print('\n=== Knob: the collisionality threshold across its band ===')
    out = {}
    for kn_crit in (0.1, 1.0, 3.0):
        settings = DispatchSettings(kn_crit=kn_crit)

        def label(f_xuv: float, settings: DispatchSettings = settings) -> str:
            return dispatch(build_state(composition, 1.0, 1.0, f_xuv, settings=settings)).regime

        lo, hi = 1.0e-2, 1.0e2
        label_lo = label(lo)
        while math.log10(hi / lo) > 1e-3:
            mid = math.sqrt(lo * hi)
            if label(mid) == label_lo:
                lo = mid
            else:
                hi = mid
        out[kn_crit] = math.sqrt(lo * hi)
        print(f'  kn_crit {kn_crit:>4}: wind sets in at F_xuv = {out[kn_crit]:.4g} W m-2')
    spread = max(out.values()) / min(out.values())
    print(f'  The boundary spans a factor {spread:.1f} across the band.')
    print('  That spread is the width of the boundary, not a parameter to tune.')
    return out


def knob_exobase_temperature() -> list[dict]:
    """Exobase temperature on a hydrostatic verdict with a light species.

    A pure heavy atmosphere returns hydrostatic rates far below the floor,
    so the sensitivity is shown where the branch does physical work: a
    small planet whose carbon dioxide carries one percent hydrogen. The
    bulk rate barely moves, because hydrogen is limited by how fast
    diffusion resupplies it, while the heavy species carry the exponential
    dependence of the Jeans flux.
    """
    print('\n=== Knob: the prescribed exobase temperature ===')
    print('  Mars-mass planet, CO2 with 1% H2, F_xuv = 0.01 W m-2')
    out, hydrogen = [], None
    for t_exo in (1000.0, 1500.0, 2000.0, 3000.0, 4000.0):
        settings = DispatchSettings(T_exo_value=t_exo)
        result = dispatch(build_state('CO2 + 1% H2', 0.107, 0.53, 0.01, settings=settings))
        if t_exo == 1000.0:
            hydrogen = result.diagnostics['hydrostatic']['detail']['species']['H2']
        out.append(
            dict(
                T_exo=t_exo,
                regime=result.regime,
                mdot=result.mdot,
                per_species=dict(result.per_species),
            )
        )
        print(
            f'  T_exo {t_exo:6.0f} K: {result.regime:<12} rate {result.mdot:.4e} kg/s, '
            f'H {result.per_species.get("H", 0.0):.3e}, '
            f'C {result.per_species.get("C", 0.0):.3e}'
        )
    heavy_span = out[-1]['per_species']['C'] / out[0]['per_species']['C']
    bulk_span = out[-1]['mdot'] / out[0]['mdot']
    print(
        f'  Over that range the bulk rate moves by a factor {bulk_span:.2f} while the '
        f'carbon rate moves by a factor {heavy_span:.1e}.'
    )
    print(
        f'  At 1000 K hydrogen sits at Jeans parameter {hydrogen["lambda_exo"]:.3g} with a '
        f'Jeans flux of {hydrogen["phi_jeans"]:.3e} against a diffusion-limited supply of '
        f'{hydrogen["phi_diffusion"]:.3e} (per unit area), so the supply is what binds and '
        f'the exobase temperature barely matters. Carbon and oxygen are Jeans limited, '
        f'which is where the exponential sensitivity went.'
    )
    return out


def knob_fractionation(composition: str = 'CO2', f_xuv: float = 10.0) -> dict:
    """The closure split against the reservoir mass-fraction split."""
    print('\n=== Knob: fractionation on and off ===')
    out = {}
    for fractionate in (True, False):
        settings = DispatchSettings(fractionate=fractionate)
        result = dispatch(build_state(composition, 1.0, 1.0, f_xuv, settings=settings))
        total = sum(result.per_species.values())
        shares = {el: rate / total for el, rate in sorted(result.per_species.items())}
        out[fractionate] = shares
        label = 'closure' if fractionate else 'reservoir mass fractions'
        print(f'  {label:<24} ' + ', '.join(f'{el} {x:.4f}' for el, x in shares.items()))
    print('  A small shift here; the light species is the one enriched.')
    return out


def knob_efficiency(composition: str = 'CO2', f_xuv: float = 10.0) -> list[dict]:
    """The energy-limited efficiency across its literature range, and the fit."""
    print('\n=== Knob: the energy-limited efficiency ===')
    out = []
    for efficiency in (0.1, 0.15, 0.3, 0.6):
        settings = DispatchSettings(efficiency=efficiency)
        result = dispatch(build_state(composition, 1.0, 1.0, f_xuv, settings=settings))
        hydro = result.diagnostics['hydrodynamic']
        out.append(dict(efficiency=efficiency, regime=result.regime, mdot=result.mdot))
        print(
            f'  epsilon {efficiency:4.2f}: {result.regime:<17} rate {result.mdot:.4e} kg/s '
            f'(energy limited {hydro["mdot_el"]:.3e}, recombination limited '
            f'{hydro["mdot_rr"]:.3e})'
        )
    fitted = dispatch(
        build_state(
            composition, 1.0, 1.0, f_xuv, settings=DispatchSettings(efficiency_mode='caldiroli')
        )
    )
    print(
        f'  fitted mode: efficiency {fitted.diagnostics["hydrodynamic"]["efficiency"]:.3g}, '
        f'flags {sorted(fitted.flags)}'
    )
    print('  The guard fires because a 1 Earth-mass planet sits below the fitted potential')
    print('  range, so the fitted value there is an extrapolation.')
    return out


def hysteresis_window(composition: str = 'CO2') -> list[dict]:
    """Show the previous label deciding the verdict inside the window."""
    print('\n=== Evolutionary use: the hysteresis window ===')
    out = []
    for f_xuv in (0.72, 0.78):
        for previous in (None, 'hydrostatic', 'hydrodynamic:EL'):
            state = build_state(composition, 1.0, 1.0, f_xuv)
            state.prev_regime = previous
            result = dispatch(state)
            knudsen = result.diagnostics['knudsen']
            out.append(dict(F_xuv=f_xuv, prev=previous, regime=result.regime))
            print(
                f'  F_xuv {f_xuv:.3f}, previously {str(previous):<17} -> '
                f'{result.regime:<17} (Kn_sc {knudsen["kn_sc"]:.4g}, threshold '
                f'{knudsen["threshold_applied"]:.4g})'
            )
    print('  A time-stepping track cannot chatter across the threshold on numerical noise.')
    return out


# -------------------------------------------------- step 6: a stellar track


def stellar_track(star=None, n_samples: int = 40) -> list[dict]:
    """Dispatch one frozen atmosphere along a stellar XUV history.

    The profile does not evolve, so this is a sequence of static snapshots
    rather than an evolutionary calculation: a real planet's structure
    responds to the loss and to the star. What the sequence does show is
    that the regime label belongs to the state and not to the planet, and
    the snapshot consistency screen reports where a frozen state stops
    being compatible with its own age.
    """
    star = star or mors.Star(Mstar=1.0, Omega=1.0)
    age_myr = np.asarray(star.Tracks['Age'])
    l_xuv = np.asarray(star.Tracks['Lx']) + np.asarray(star.Tracks['Leuv'])
    f_xuv = l_xuv * 1e-7 / (4.0 * math.pi * A_TRACK**2)  # erg/s to W, then flux
    l_bol = np.asarray(star.Tracks['Lbol']) * 1e-7  # erg/s to W

    composition = 'CO2 + 1% H2'
    comp = COMPOSITIONS[composition]
    m_p, r_p = Me, Re
    t_eq_track = np.array([t_eq_at(A_TRACK, float(L)) for L in l_bol])
    # One profile, built at the median equilibrium temperature of the track.
    profile = isothermal_profile(m_p, r_p, float(np.median(t_eq_track)), comp, P_SURF, P_TOP)
    reservoirs = element_reservoirs(comp, 100.0 * Me_atm)

    print(f'\n=== A stellar history: {composition}, 1 Me, 1 Re at 0.2 au ===')
    print(
        f'  {"age [Myr]":>10}  {"F_xuv":>8}  {"regime":<17} {"rate [kg/s]":>12}  '
        f'{"Kn_sc":>9}  snapshot'
    )
    indices = np.unique(np.linspace(0, age_myr.size - 1, n_samples).astype(int))
    rows, previous = [], None
    for i in indices:
        state = EscapeInputs(
            M_p=m_p,
            R_p=r_p,
            M_star=Ms,
            a=A_TRACK,
            e=0.0,
            T_eq=float(t_eq_track[i]),
            F_xuv=float(f_xuv[i]),
            F_bol=f_bol_at(A_TRACK, float(l_bol[i])),
            F_int=F_INT,
            kappa_photo=KAPPA_PHOTO,
            profile=profile,
            settings=DispatchSettings(),
            prev_regime=previous,
            age=float(age_myr[i]) * 1e6 * SEC_PER_YR,
            reservoirs=dict(reservoirs),
        )
        result = dispatch(state)
        previous = result.regime
        consistency = result.diagnostics['self_consistency']
        rows.append(
            dict(
                age_Myr=float(age_myr[i]),
                F_xuv=float(f_xuv[i]),
                regime=result.regime,
                mdot=result.mdot,
                kn_sc=result.diagnostics['knudsen']['kn_sc'],
                inconsistent=consistency.get('inconsistent'),
                per_species=dict(result.per_species),
            )
        )
    for row in rows[::6]:
        verdict = 'inconsistent' if row['inconsistent'] else 'consistent'
        print(
            f'  {row["age_Myr"]:10.1f}  {row["F_xuv"]:8.3g}  {row["regime"]:<17} '
            f'{row["mdot"]:12.3e}  {row["kn_sc"]:9.3g}  {verdict}'
        )
    for before, after in zip(rows[:-1], rows[1:]):
        if before['regime'] != after['regime']:
            print(
                f'  label changes from {before["regime"]} to {after["regime"]} between '
                f'{before["age_Myr"]:.0f} and {after["age_Myr"]:.0f} Myr, and the rate '
                f'drops from {before["mdot"]:.3e} to {after["mdot"]:.3e} kg/s'
            )
    return rows


# --------------------------------------------------------------- the figure


def boundary_band(composition: str, m_earth: float = 1.0, r_earth: float = 1.0) -> tuple:
    """The wind boundary at both edges of the collisionality band.

    Returns the flux where the wind sets in for a threshold of 3 and for a
    threshold of 0.1, which bracket the default of 1. The pair is the width
    the criterion implies, and it is what the figure shades.
    """
    edges = []
    for kn_crit in (3.0, 0.1):
        settings = DispatchSettings(kn_crit=kn_crit)
        lo, hi = 1.0e-2, 1.0e2
        label_lo = dispatch(
            build_state(composition, m_earth, r_earth, lo, settings=settings)
        ).regime
        while math.log10(hi / lo) > 1e-3:
            mid = math.sqrt(lo * hi)
            state = build_state(composition, m_earth, r_earth, mid, settings=settings)
            if dispatch(state).regime == label_lo:
                lo = mid
            else:
                hi = mid
        edges.append(math.sqrt(lo * hi))
    return tuple(edges)


def make_figure(sweeps: dict, boundaries: dict, bands: dict, outpath: str) -> None:
    """Two panels: the same sweep for two compositions, colored by regime."""
    palette = brand_colors()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    y_floor = 1.0e3  # rates below this round to zero on these planets

    for ax, (composition, rows) in zip(axes, sweeps.items()):
        # The dispatched rate, as a line in the escape domain color.
        ax.plot(
            [r['F_xuv'] for r in rows],
            [max(r['mdot'], y_floor) for r in rows],
            color=palette['escape'],
            linewidth=1.2,
            alpha=0.5,
            zorder=0,
        )
        for label in ('hydrostatic', 'hydrodynamic:EL', 'hydrodynamic:RR'):
            group = [r for r in rows if r['regime'] == label]
            if not group:
                continue
            # Filled markers where the rate is on scale, open markers where
            # it has been raised to the axis floor to stay visible.
            for on_scale in (True, False):
                subset = [r for r in group if (r['mdot'] >= y_floor) is on_scale]
                if not subset:
                    continue
                ax.plot(
                    [r['F_xuv'] for r in subset],
                    [max(r['mdot'], y_floor) for r in subset],
                    linestyle='none',
                    marker=REGIME_MARKERS[label],
                    markersize=8,
                    markerfacecolor=palette[label] if on_scale else 'none',
                    markeredgecolor=palette[label],
                    label=label
                    if on_scale or not any(r['mdot'] >= y_floor for r in group)
                    else None,
                )
        # The collisionality criterion spans a factor of 30 in threshold, so
        # the boundary it sets is a band. Shade the band before the line.
        band = bands.get(composition)
        if band:
            ax.axvspan(min(band), max(band), color=palette['hydrostatic'], alpha=0.18, lw=0)
            ax.annotate(
                'switch band',
                xy=(math.sqrt(band[0] * band[1]), 1.0e9),
                ha='center',
                fontsize=11,
                color=palette['rule'],
                fontfamily='Spline Sans Mono',
            )
        for _before, _after, flux in boundaries.get(composition, []):
            ax.axvline(flux, color=palette['rule'], linestyle='--', linewidth=1.0)
            ax.annotate(
                f'{flux:.3g} W m$^{{-2}}$',
                xy=(flux, 1.2e10),
                xytext=(3, 0),
                textcoords='offset points',
                fontsize=11,
                color=palette['rule'],
                fontfamily='Spline Sans Mono',
            )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(3.0e2, 3.0e10)
        ax.set_xlabel('XUV flux [W m$^{-2}$]')
        ax.set_title(_panel_title(composition))

    axes[0].set_ylabel('Mass loss rate [kg s$^{-1}$]')
    axes[0].annotate(
        'open markers: rate rounds to zero',
        xy=(0.03, 0.14),
        xycoords='axes fraction',
        fontsize=11,
        fontfamily='Spline Sans Mono',
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.11))
    fig.subplots_adjust(bottom=0.25, top=0.9, left=0.09, right=0.98, wspace=0.08)
    fig.savefig(outpath, bbox_inches=None)
    fig.savefig(outpath.replace('.pdf', '.png'), bbox_inches=None)
    plt.close(fig)
    print(f'\nwrote {outpath} and its PNG companion')


def make_track_figure(rows: list[dict], outpath: str) -> None:
    """Two stacked panels summarizing one atmosphere along a stellar history.

    The upper panel carries the bulk rate against age, marked by regime,
    with the ages the consistency screen rejects shaded. The lower panel
    carries the per-element rates, which is where the change of branch
    shows its consequence: the wind takes the heavy elements with the
    hydrogen, and the exosphere takes hydrogen alone.
    """
    palette = brand_colors()
    fig, (ax_rate, ax_species) = plt.subplots(
        2, 1, figsize=(8.5, 7.0), sharex=True, height_ratios=[1.0, 1.0]
    )
    ages = [row['age_Myr'] for row in rows]

    # The span the snapshot screen rejects, drawn behind everything.
    bad = [row['age_Myr'] for row in rows if row['inconsistent']]
    if bad:
        for ax in (ax_rate, ax_species):
            ax.axvspan(min(bad), max(bad), color=palette['hydrostatic'], alpha=0.15, lw=0)
    if bad:
        ax_rate.annotate(
            'the snapshot screen rejects these ages',
            xy=(math.sqrt(min(bad) * max(bad)), 3.0e7),
            ha='center',
            fontsize=11,
            color=palette['rule'],
            fontfamily='Spline Sans Mono',
        )

    ax_rate.plot(
        ages, [row['mdot'] for row in rows], color=palette['escape'], lw=1.2, alpha=0.5
    )
    for label in ('hydrodynamic:EL', 'hydrostatic'):
        group = [row for row in rows if row['regime'] == label]
        if not group:
            continue
        ax_rate.plot(
            [row['age_Myr'] for row in group],
            [row['mdot'] for row in group],
            linestyle='none',
            marker=REGIME_MARKERS[label],
            markersize=7,
            markerfacecolor=palette[label],
            markeredgecolor=palette[label],
            label=label,
        )
    # The age at which the label changes, from the sequence itself.
    for before, after in zip(rows[:-1], rows[1:]):
        if before['regime'] != after['regime']:
            crossing = math.sqrt(before['age_Myr'] * after['age_Myr'])
            for ax in (ax_rate, ax_species):
                ax.axvline(crossing, color=palette['rule'], linestyle='--', lw=1.0)
            ax_rate.annotate(
                f'wind ends near {crossing:.0f} Myr',
                xy=(crossing, 1.0e2),
                xytext=(-6, 0),
                textcoords='offset points',
                ha='right',
                fontsize=11,
                color=palette['rule'],
                fontfamily='Spline Sans Mono',
            )

    elements = sorted({el for row in rows for el in row['per_species']})
    species_colors = {'H': palette['escape'], 'C': palette['boiloff'], 'O': palette['solar']}
    for element in elements:
        ax_species.plot(
            ages,
            [max(row['per_species'].get(element, 0.0), 1e-30) for row in rows],
            color=species_colors.get(element, palette['rule']),
            label=element,
        )

    ax_rate.set_yscale('log')
    ax_rate.set_ylabel('Bulk rate [kg s$^{-1}$]')
    ax_rate.set_ylim(1.0e-2, 3.0e8)
    ax_rate.legend(loc='lower left', title='regime')
    ax_species.set_xscale('log')
    ax_species.set_yscale('log')
    ax_species.set_ylim(1.0e-25, 3.0e7)
    ax_species.set_xlabel('Stellar age [Myr]')
    ax_species.set_ylabel('Element rate [kg s$^{-1}$]')
    ax_species.legend(loc='lower left', title='element', ncol=3)
    ax_species.annotate(
        'heavy elements leave with the wind, and stop with it',
        xy=(0.03, 0.62),
        xycoords='axes fraction',
        fontsize=11,
        color=palette['rule'],
        fontfamily='Spline Sans Mono',
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.09, hspace=0.08)
    fig.savefig(outpath, bbox_inches=None)
    fig.savefig(outpath.replace('.pdf', '.png'), bbox_inches=None)
    plt.close(fig)
    print(f'wrote {outpath} and its PNG companion')


def _panel_title(composition: str) -> str:
    """Composition name with typeset subscripts for a panel title."""
    return {
        'CO2': 'CO$_2$',
        'N2-O2': 'N$_2$ and O$_2$',
        'CO2 + 1% H2': 'CO$_2$ with 1% H$_2$',
        'H/He': 'H$_2$ and He',
    }.get(composition, composition)


# ----------------------------------------------------------------- driver


def main(outdir: str = 'output') -> dict:
    """Run every step, print the tables, and write the figure."""
    results = {}
    results['verdict'] = one_verdict()
    sweeps, boundaries, bands = {}, {}, {}
    for composition in ('CO2', 'N2-O2'):
        rows = flux_sweep(composition)
        sweeps[composition] = rows
        boundaries[composition] = report_sweep(composition, rows)
        bands[composition] = boundary_band(composition)
        print(
            f'  the wind boundary spans {min(bands[composition]):.3g} to '
            f'{max(bands[composition]):.3g} W m-2 across the collisionality band'
        )
    results['sweeps'] = sweeps
    results['boundaries'] = boundaries
    results['bands'] = bands
    results['extremes'] = extreme_labels()
    results['diagnostics'] = read_diagnostics()
    results['knob_kn'] = knob_collisionality()
    results['knob_t_exo'] = knob_exobase_temperature()
    results['knob_fractionation'] = knob_fractionation()
    results['knob_efficiency'] = knob_efficiency()
    results['hysteresis'] = hysteresis_window()
    results['track'] = stellar_track()
    make_figure(sweeps, boundaries, bands, f'{outdir}/demo_dispatcher_regimes.pdf')
    make_track_figure(results['track'], f'{outdir}/demo_dispatcher_track.pdf')
    return results


if __name__ == '__main__':
    main()
