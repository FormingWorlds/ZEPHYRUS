"""Dispatching an escape regime: the script behind docs/Tutorials/dispatch.md.

Runs top to bottom, one section per tutorial step, and prints each step's
output under a header. The tutorial page includes each section verbatim,
between the ``# --8<--`` marker lines, so the page and this file cannot
drift apart. Writes two figures to ``output/``. From the repository root:

    python examples/demo_dispatcher/demo_dispatcher.py

The last step reads the MORS stellar tracks (``mors download all``).
"""

# --8<-- [start:setup]
import dataclasses
import math
import os

import mors
import numpy as np
from scipy.constants import Stefan_Boltzmann as SIGMA_SB

from zephyrus import plot
from zephyrus.composition import atomize, species_mass_amu
from zephyrus.constants import au2m, rate_floor, yr_julian
from zephyrus.dispatcher import DispatchSettings, EscapeInputs, dispatch
from zephyrus.planets_parameters import Ls, Me, Me_atm, Ms, Re
from zephyrus.profiles import isothermal_profile

os.makedirs('output', exist_ok=True)
plot.apply_theme()  # PROTEUS figure style, when proteus-mpl is installed
# --8<-- [end:setup]

# --8<-- [start:step1]
print('=== Step 1: build a planetary state ===')

A_REF = 0.0775 * au2m  # orbit where a Sun-like star gives T_eq = 1000 K  [m]
P_SURF, P_TOP = 1.0e7, 1.0e-5  # profile from 100 bar up to 0.1 nbar      [Pa]

COMPOSITIONS = {  # mole fractions, the same at every level
    'CO2': {'CO2': 1.0},
    'N2-O2': {'N2': 0.8, 'O2': 0.2},
    'H/He': {'H2': 0.9, 'He': 0.1},
    'CO2 + 1% H2': {'CO2': 0.99, 'H2': 0.01},
}


def build_state(composition, m_earth, r_earth, f_xuv, settings=None, a=A_REF, inventory=Me_atm):
    """Planet of m_earth Earth masses and r_earth Earth radii, isothermal atmosphere."""
    gas = COMPOSITIONS[composition]
    M_p, R_p = m_earth * Me, r_earth * Re
    T_eq = (Ls / (16.0 * math.pi * SIGMA_SB * a**2)) ** 0.25  # zero albedo  [K]
    elements = atomize(gas)  # mole fraction of each element
    norm = sum(x * species_mass_amu(el) for el, x in elements.items())
    return EscapeInputs(
        M_p=M_p,  # planet mass                       [kg]
        R_p=R_p,  # planet radius                     [m]
        M_star=Ms,  # stellar mass                    [kg]
        a=a,  # semi-major axis                       [m]
        e=0.0,  # eccentricity
        T_eq=T_eq,  # equilibrium temperature         [K]
        F_xuv=f_xuv,  # XUV flux at the planet        [W m-2]
        F_bol=Ls / (4.0 * math.pi * a**2),  # bolometric flux  [W m-2]
        F_int=1.0,  # interior heat flux              [W m-2]
        kappa_photo=0.01,  # photospheric opacity     [m2 kg-1]
        profile=isothermal_profile(M_p, R_p, T_eq, gas, P_SURF, P_TOP),
        settings=settings or DispatchSettings(),
        age=1.0e8 * yr_julian,  # age of the snapshot, 100 Myr  [s]
        reservoirs={  # element inventories             [kg]
            el: inventory * x * species_mass_amu(el) / norm for el, x in elements.items()
        },
    )


state = build_state('CO2', 1.0, 1.0, 10.0)
profile = state.profile
print(f'T_eq        {state.T_eq:.1f} K')
print(f'F_bol       {state.F_bol:.4g} W m-2')
print(f'levels      {profile.p.size}, from {profile.p[0]:.3g} to {profile.p[-1]:.3g} Pa')
print(f'radius      {profile.r[0] / Re:.3f} to {profile.r[-1] / Re:.3f} Earth radii')
print(f'reservoirs  {", ".join(f"{el} {m:.4g} kg" for el, m in state.reservoirs.items())}')
# --8<-- [end:step1]

# --8<-- [start:step2]
print('=== Step 2: dispatch one state and read the result ===')

result = dispatch(state)
total = sum(result.per_species.values())

print(f'regime       {result.regime}')
print(f'mdot         {result.mdot:.4e} kg/s')
for element, rate in result.per_species.items():
    print(f'  {element:<10} {rate:.4e} kg/s')
print(f'closure      {abs(total / result.mdot - 1.0) < 1e-12}')
print(f'flags        {result.flags}')
print(f'diagnostics  {len(result.diagnostics)} groups')
# --8<-- [end:step2]

# --8<-- [start:step3_sweep]
print('=== Step 3: sweep the XUV flux across two regime boundaries ===')

fluxes = np.logspace(-2, math.log10(5.0e3), 30)  # [W m-2]
sweeps = {
    composition: [dispatch(build_state(composition, 1.0, 1.0, f)) for f in fluxes]
    for composition in ('CO2', 'N2-O2')
}

print(
    f'{"F_xuv [W m-2]":>13}  {"regime":<16} {"mdot [kg/s]":>11}  {"Kn_sc":>8}  '
    'label at Kn_crit 0.1 | 3'
)
for f_xuv, out in list(zip(fluxes, sweeps['CO2']))[::3]:
    kn = out.diagnostics['knudsen']
    edges = kn['counterfactual_labels']
    print(
        f'{f_xuv:13.3g}  {out.regime:<16} {out.mdot:11.3e}  {kn["kn_sc"]:8.3g}  '
        f'{edges[0.1]} | {edges[3.0]}'
    )
# --8<-- [end:step3_sweep]


# --8<-- [start:step3_boundaries]
def boundary_flux(composition, lo, hi, settings=None):
    """Flux [W m-2] between lo and hi where the label changes, to 0.001 dex."""

    def label(f_xuv):
        return dispatch(build_state(composition, 1.0, 1.0, f_xuv, settings)).regime

    first = label(lo)
    while math.log10(hi / lo) > 1.0e-3:  # halve the bracket in log flux
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if label(mid) == first else (lo, mid)
    return math.sqrt(lo * hi)


boundaries = {
    'CO2': [boundary_flux('CO2', 0.1, 10.0), boundary_flux('CO2', 100.0, 2000.0)],
    'N2-O2': [boundary_flux('N2-O2', 0.01, 10.0)],
}
bands = {  # the wind boundary at the two edges of the Knudsen criterion
    composition: tuple(
        boundary_flux(composition, 0.01, 100.0, DispatchSettings(kn_crit=kn_crit))
        for kn_crit in (3.0, 0.1)
    )
    for composition in ('CO2', 'N2-O2')
}

for composition in ('CO2', 'N2-O2'):
    lo, hi = bands[composition]
    lowest = sweeps[composition][0]
    print(composition)
    print(f'  boundaries        {", ".join(f"{b:.3g}" for b in boundaries[composition])} W m-2')
    print(f'  wind onset band   {lo:.3g} to {hi:.3g} W m-2, a factor {hi / lo:.2g}')
    print(
        f'  rate at 0.01      {lowest.mdot:.2e} kg/s, above the floor: '
        f'{lowest.diagnostics["rate_floor"]["above_floor"]}'
    )
print(f'rate floor          {rate_floor:.2e} kg/s')
# --8<-- [end:step3_boundaries]

# --8<-- [start:step3_figure]
panels = [
    dict(
        title=title,
        F_xuv=fluxes,
        mdot=[out.mdot for out in sweeps[composition]],
        regime=[out.regime for out in sweeps[composition]],
        boundaries=boundaries[composition],
        band=bands[composition],
    )
    for composition, title in (('CO2', 'CO$_2$'), ('N2-O2', 'N$_2$ and O$_2$'))
]
fig = plot.plot_regime_sweep(panels)
fig.savefig('output/dispatcher_regime_sweep.png')
fig.savefig('output/dispatcher_regime_sweep.pdf')
# --8<-- [end:step3_figure]

# --8<-- [start:step4_boiloff]
print('=== Step 4: boil-off and Roche overflow ===')

inflated = dispatch(build_state('H/He', 1.0, 1.5, 10.0))

print(f'regime     {inflated.regime}')
print(f'mdot       {inflated.mdot:.3e} kg/s')
print(f'Lambda     {inflated.diagnostics["lambda_gate"]:.3g}')
print(f'threshold  20, band {inflated.diagnostics["documentation"]["lambda_crit_band"]}')
print(f'flags      {sorted(inflated.flags)}')
# --8<-- [end:step4_boiloff]

# --8<-- [start:step4_overflow]
admitted = DispatchSettings(residual_mode='luminosity_capped')
puffy = dispatch(build_state('H/He', 3.0, 1.7, 0.1, settings=admitted))
roche = puffy.diagnostics['roche']

print(f'regime            {puffy.regime}')
print(f'mdot              {puffy.mdot:.3e} kg/s, from the {roche["rate_branch"]} rate')
print(f'flow radius       {roche["flow_radius"] / roche["R_hill_periapsis"]:.3f} Hill radii')
print(
    f'atmosphere        {roche["r_atmosphere"] / puffy.diagnostics["nozzle"]["r_lobe"]:.3f} '
    f'Roche-lobe radii'
)
print(f'roche_subflag     {puffy.flags["roche_subflag"]}')
# --8<-- [end:step4_overflow]

# --8<-- [start:step4_default]
default = dispatch(build_state('H/He', 3.0, 1.7, 0.1))
bolometric = default.diagnostics['bolometric']
nozzle = default.diagnostics['nozzle']

print(f'regime            {default.regime}')
print(f'mdot              {default.mdot:.3e} kg/s')
print(f'flow radius       {1.0 / default.diagnostics["roche"]["xi_flow"]:.3f} Hill radii')
print(f'flags             {sorted(default.flags)}')
print(
    f'bolometric        {bolometric["rate_kg_s"]:.3e} kg/s, competes: {bolometric["competes"]}'
)
print(
    f'L1 transfer       {nozzle["rate_full_orbit_kg_s"]:.3e} kg/s, competes: '
    f'{nozzle["rate_kg_s"] > 0.0} (R_sonic / R_L1 = {nozzle["R_sonic_over_R_L1"]:.3f})'
)
# --8<-- [end:step4_default]

# --8<-- [start:step5]
print('=== Step 5: read the diagnostics ===')

d = result.diagnostics  # the carbon dioxide planet of step 2
kn, hydro = d['knudsen'], d['hydrodynamic']

print('Which branch, and how close was the switch?')
print(f'  Kn_sc {kn["kn_sc"]:.3g} against threshold {kn["threshold_applied"]:.3g}')
print(f'  label at the band edges: {kn["counterfactual_labels"]}')
print('What set the rate?')
print(f'  EL {hydro["mdot_el"]:.3e}, RR {hydro["mdot_rr"]:.3e}, PL {hydro["mdot_pl"]:.3e} kg/s')
print(
    f'  {hydro["selection_mechanism"]}, T_wind {hydro["T_wind"]:.0f} K, '
    f'barometric factor {hydro["rr_chain"]["barometric_factor"]:.2g}'
)
print('Could the heating drive a transonic flow?')
print(f'  Q_net / Q_c {d["johnson_q"]["q_net_over_qc"]:.3g}')
print('How does the verdict read in other classifications?')
guo = d['guo_triple']
print(
    f'  lambda_exo {guo["lambda_exo"]:.4g}, lambda_rp {guo["lambda_rp"]:.4g}, '
    f'lambda_star {guo["lambda_star"]:.4g}'
)
print(
    f'  log10(-phi) {d["potential_screens"]["log_minus_phi_cgs"]:.4g} (cgs), '
    f'Salz screen: {d["potential_screens"]["salz_verdict"]}'
)
print('Is the snapshot self-consistent?')
consistency = d['self_consistency']
print(f'  worst Kn below the sonic point {d["fluid_check"]["worst_kn"]:.3g}')
print(
    f'  inventory lasts {consistency["t_deplete_s"] / yr_julian:.3g} yr, '
    f'age {consistency["age_s"] / yr_julian:.3g} yr, '
    f'inconsistent: {consistency["inconsistent"]}'
)
print('Where was the wind launched, and where did the cross sections come from?')
print(
    f'  wind base {d["base_level"]["p_Pa"]:.3g} Pa, '
    f'clamped: {d["base_level"]["clamp_decades"] is not None}'
)
print(f'  cross sections: {kn["provenance"]}')
# --8<-- [end:step5]

# --8<-- [start:step6_kn]
print('=== Step 6: settings that move a boundary ===')

for kn_crit in (0.1, 1.0, 3.0):
    out = dispatch(build_state('CO2', 1.0, 1.0, 1.0, DispatchSettings(kn_crit=kn_crit)))
    print(f'kn_crit {kn_crit:<4} at 1 W m-2: {out.regime}')
# --8<-- [end:step6_kn]

# --8<-- [start:step6_texo]
print(f'{"T_exo [K]":>9}  {"regime":<12} {"mdot [kg/s]":>11}  {"H [kg/s]":>9}  {"C [kg/s]":>9}')
for t_exo in (1000.0, 2000.0, 4000.0):
    settings = DispatchSettings(T_exo_value=t_exo)
    out = dispatch(build_state('CO2 + 1% H2', 0.107, 0.53, 0.01, settings))
    rates = out.per_species
    print(
        f'{t_exo:9.0f}  {out.regime:<12} {out.mdot:11.4g}  {rates["H"]:9.4g}  {rates["C"]:9.3g}'
    )

h2 = out.diagnostics['hydrostatic']['detail']['species']['H2']  # at 4000 K
print(
    f'H2 at 4000 K: lambda_exo {h2["lambda_exo"]:.3g}, '
    f'Jeans flux / diffusion supply {h2["phi_jeans"] / h2["phi_diffusion"]:.3g}'
)
# --8<-- [end:step6_texo]

# --8<-- [start:step6_efficiency]
print(f'{"epsilon":>7}  {"F_xuv":>6}  {"regime":<16} {"EL":>9}  {"RR":>9}  {"PL":>9}  [kg/s]')
for eps in (0.1, 0.3, 0.6):
    for f_xuv in (10.0, 1.0e3):
        out = dispatch(build_state('CO2', 1.0, 1.0, f_xuv, DispatchSettings(efficiency=eps)))
        hy = out.diagnostics['hydrodynamic']
        print(
            f'{eps:7.1f}  {f_xuv:6.0f}  {out.regime:<16} {hy["mdot_el"]:9.3g}  '
            f'{hy["mdot_rr"]:9.3g}  {hy["mdot_pl"]:9.3g}'
        )
print(f'photon-limit efficiency {hy["efficiency_photon_limit"]:.3g}')
# --8<-- [end:step6_efficiency]

# --8<-- [start:step6_fractionation]
for fractionate in (True, False):
    out = dispatch(
        build_state('CO2', 1.0, 1.0, 10.0, DispatchSettings(fractionate=fractionate))
    )
    shares = {el: rate / out.mdot for el, rate in out.per_species.items()}
    print(
        f'fractionate {fractionate!s:<5}  '
        + ', '.join(f'{el} {s:.4f}' for el, s in shares.items())
    )
# --8<-- [end:step6_fractionation]

# --8<-- [start:step6_hysteresis]
print(f'{"F_xuv":>5}  {"previous label":<16} {"regime":<16} {"Kn_sc":>6}  threshold')
for f_xuv in (0.72, 0.78):
    for previous in (None, 'hydrostatic', 'hydrodynamic:EL'):
        hysteresis_state = build_state('CO2', 1.0, 1.0, f_xuv)
        hysteresis_state.prev_regime = previous
        out = dispatch(hysteresis_state)
        kn = out.diagnostics['knudsen']
        print(
            f'{f_xuv:5.2f}  {previous!s:<16} {out.regime:<16} {kn["kn_sc"]:6.3f}  '
            f'{kn["threshold_applied"]:.3g}'
        )
# --8<-- [end:step6_hysteresis]

# --8<-- [start:step7_track]
print('=== Step 7: follow a stellar XUV history ===')

A_TRACK = 0.2 * au2m  # [m]
star = mors.Star(Mstar=1.0, Omega=1.0)  # 1 solar mass, spinning at today's solar rate at 1 Myr
age_myr = np.asarray(star.Tracks['Age'])
l_xuv = (np.asarray(star.Tracks['Lx']) + np.asarray(star.Tracks['Leuv'])) * 1.0e-7  # erg/s to W
l_bol = np.asarray(star.Tracks['Lbol']) * 1.0e-7  # erg/s to W

# One frozen atmosphere with 100 Earth atmospheres of gas, at 40 ages along the track.
planet = build_state('CO2 + 1% H2', 1.0, 1.0, 0.0, a=A_TRACK, inventory=100.0 * Me_atm)
track, previous = [], None
for i in np.unique(np.linspace(0, age_myr.size - 1, 40).astype(int)):
    snapshot = dataclasses.replace(
        planet,
        T_eq=(l_bol[i] / (16.0 * math.pi * SIGMA_SB * A_TRACK**2)) ** 0.25,
        F_xuv=l_xuv[i] / (4.0 * math.pi * A_TRACK**2),
        F_bol=l_bol[i] / (4.0 * math.pi * A_TRACK**2),
        age=age_myr[i] * 1.0e6 * yr_julian,
        prev_regime=previous,
    )
    out = dispatch(snapshot)
    previous = out.regime
    track.append(dict(age_myr=age_myr[i], F_xuv=snapshot.F_xuv, result=out))

print(
    f'{"age [Myr]":>9}  {"F_xuv":>6}  {"regime":<16} {"mdot [kg/s]":>11}  {"Kn_sc":>8}  '
    'consistent'
)
for row in track[::6]:
    out = row['result']
    print(
        f'{row["age_myr"]:9.1f}  {row["F_xuv"]:6.3g}  {out.regime:<16} {out.mdot:11.3e}  '
        f'{out.diagnostics["knudsen"]["kn_sc"]:8.3g}  '
        f'{not out.diagnostics["self_consistency"]["inconsistent"]}'
    )
for before, after in zip(track[:-1], track[1:]):
    if before['result'].regime != after['result'].regime:
        print(
            f'label changes between {before["age_myr"]:.0f} and {after["age_myr"]:.0f} Myr: '
            f'{before["result"].mdot:.3e} to {after["result"].mdot:.3e} kg/s'
        )
# --8<-- [end:step7_track]

# --8<-- [start:step7_figure]
fig = plot.plot_track(
    age_myr=[row['age_myr'] for row in track],
    mdot=[row['result'].mdot for row in track],
    regimes=[row['result'].regime for row in track],
    per_species=[row['result'].per_species for row in track],
    inconsistent=[
        row['result'].diagnostics['self_consistency']['inconsistent'] for row in track
    ],
)
fig.savefig('output/dispatcher_track.png')
fig.savefig('output/dispatcher_track.pdf')
# --8<-- [end:step7_figure]
