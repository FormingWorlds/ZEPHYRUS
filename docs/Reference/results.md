# Dispatch results

The [parameter reference](parameters.md) documents what goes into `zephyrus.dispatch`. This page documents what comes out: the result fields, every flag the framework can raise, and every group in the diagnostics container. The physics behind each quantity is on the [escape regimes](../Explanations/regimes.md) page, and the [tutorial](../Tutorials/dispatch.md) shows the reading order in practice. For a flag that fired on a result you did not expect, the [troubleshooting guide](../How-to/troubleshooting.md) is faster than this page.

---

## The result

`dispatch` returns an `EscapeResult` with five fields.

| Field | Type | What it guarantees |
|---|---|---|
| `regime` | str | One of the five labels below. Always set. |
| `mdot` | float | Bulk mass-loss rate in kg s⁻¹, finite and non-negative. |
| `per_species` | dict | Element symbol to rate in kg s⁻¹, non-negative, summing to `mdot` at machine precision. |
| `flags` | dict | Everything that clamped, fell back, or was screened. Empty when there is nothing to report. |
| `diagnostics` | dict | Reporting container. Nothing in the dispatch logic reads it, and it has no off switch. |

Every physically posed state returns a result. A `ValueError` means the state or the settings are malformed, not that the physics failed: a non-positive mass, radius, stellar mass, orbital distance, equilibrium temperature, bolometric flux, interior flux, or opacity, a negative XUV flux, an eccentricity outside $[0, 1)$, a profile whose pressure does not decrease or whose radius does not increase with index, fewer than three profile levels, mixing-ratio arrays of the wrong length, an unsupported option string, or all four cooling channels disabled at once.

## The five labels

| Label | Physics | Rate |
|---|---|---|
| `boiloff` | Bolometrically driven outflow from an atmosphere too weakly bound to hold itself, which the restricted Jeans parameter tests. The launch level can sit inside or outside the sonic radius; when it is outside, the Mach number is clamped to 1 and `bondi_inflated` says so. | Closed-form transonic Parker wind, Bondi-capped; past the activation gate the same machinery is luminosity-capped and dispatched only when the `residual` setting admits it. |
| `hydrodynamic:EL` | A collisional XUV-driven wind whose rate is set by the energy budget. | The smaller of the two hydrodynamic limits, the energy-limited one winning. |
| `hydrodynamic:RR` | The same wind, with the recombination-limited rate winning the minimum. | Read `rr_chain['barometric_factor']` to see whether the rate is set by the recombination-limited base ionization (factor near 1) or by the wind failing to reach the sonic point (factor decades below 1). |
| `hydrostatic` | Too rarefied for a wind; per-species Jeans escape from the exobase, capped by diffusive resupply. | Natively per-species; heavy-element rates are lower limits. |
| `roche_overflow` | Either the active flow radius reaches the periapsis Hill radius, so the flow is not bound to the planet, or the tidally driven transfer through L1 outruns every bound candidate. | Two readings, separated by `diagnostics['roche']['rate_branch']`. When the geometric screen renamed a bound branch, the rate is that branch's own, a bound-flow estimate, and a lower limit on the tidal transfer only where the nozzle candidate sat outside its applicability criterion; where it was applicable and lost, `diagnostics['nozzle']['rate_kg_s']` is below the dispatched rate. When `rate_branch` reads `roche_overflow`, the rate is the L1 transfer rate of Jackson et al. (2017) itself, dispatched because it won the final comparison; that boundary is a rate crossing, so the dispatched rate is continuous across it. |

A sixth label, `impact`, is reserved for the [giant-impact channel](../Explanations/impacts.md), which the caller invokes directly rather than through the dispatcher.

---

## Flags

The flags dictionary is a warning set: every entry means something happened that you should know about, and a result with nothing to report carries an empty dictionary. Most entries are `True`; five carry a value instead, because the warning is more useful with its magnitude attached (`base_clamp_decades`, `base_method_fallback`, `roche_subflag`, `rock_former_bij`, `thermostat_clamped`). Quantities that merely describe a call, rather than warning about it, live in the diagnostics instead.

The `effect` column says whether the returned rate already reflects the flag or whether the flag reports only.

### The working levels

| Flag | Value | Meaning | Effect |
|---|---|---|---|
| `base_clamped` | `True` | The physical wind-base pressure lies above the profile top, so the base was clamped to the topmost level. | The rate reflects it |
| `base_clamp_decades` | float | How far that clamp moved the base, in pressure decades. | Reporting only |
| `base_extended` | `True` | With `base_out_of_range = 'extend'`, the base was evaluated on the extended upper structure instead of clamping. | The rate reflects it |
| `base_extension_truncated` | `True` | The extension became unbound before reaching the base pressure, so the clamp stands after all. | The rate reflects it |
| `base_method_fallback` | `'lopez'` | The requested base method was unavailable or did not converge, so the default was used. | The rate reflects it |
| `photo_clamped` | `True` | The profile does not span the photospheric level, so the nearest end level was used for the energy-limited geometry. | The rate reflects it |

### The branches

| Flag | Value | Meaning | Effect |
|---|---|---|---|
| `bondi_inflated` | `True` | The launch level sits above the Bondi radius, so the Mach number was capped at one. | The rate reflects it |
| `thermostat_clamped` | `'high'` or `'low'` | The heating against cooling balance had no root inside the bracket, so the wind temperature clamped to the nearer edge. A high clamp at a dense base is collisional quenching of the line coolants, not a failure. | The rate reflects it A clamped value is the edge of the bracket, not a solution: at the high edge the balance usually does have a root near 1e5 K, outside the model rather than inside it, because the coolants are neutral three-level systems and the gas there is fully ionized. Every quantity built on the wind temperature inherits that, including the sonic radius, the recombination-limited rate, and the sonic-point Knudsen number. |
| `subcritical_sonic` | `True` | The isothermal sonic radius came out below the wind base (base Jeans parameter under 2), so it was floored at the base and the barometric factor dropped. A recombination-limited win under this flag is a floored placeholder rather than a rate. | The rate reflects it |
| `caldiroli_out_of_box` | `True` | The fitted efficiency was evaluated outside the range of gravitational potential and flux it was fitted on. The value is returned as an extrapolation. | The rate reflects it |
| `caldiroli_below_flux_bound` | `True` | Below the validity bound of the efficiency fit, where its formulas turn complex. | Rejected |
| `efficiency_fallback_fixed` | `True` | The fitted efficiency was unavailable, so the fixed setting was used. | The rate reflects it |
| `bolometric_residual` | `True` | The luminosity-capped bolometric residual, admitted by `residual = 'luminosity_capped'`, beat the branch that won above and took the label. Never raised under the default setting, which reports the candidate without dispatching it. | The rate reflects it |
| `gate_rerouted` | `True` | The exobase was too hot to stay hydrostatic, so the state was routed back to the hydrodynamic rate. | The rate reflects it |
| `contested_ion` | `True` | The neutral and plasma escape-temperature conventions disagree about the branch. Both rates are recorded in `diagnostics['contested_ion']`. | Reporting only |
| `hysteresis_active` | `True` | A previous regime label was supplied, so the hysteresis window was available. It says the memory was offered, not that it changed anything: read `diagnostics['knudsen']['threshold_applied']` against `kn_crit` to see whether the window actually moved the threshold on this call. | The rate reflects it |

### The hydrostatic branch

| Flag | Value | Meaning | Effect |
|---|---|---|---|
| `hydrostatic_lower_limit` | `True` | Always set on this branch: the nonthermal channels that dominate heavy-species loss from real exospheres are not modeled, so heavy-element rates are lower limits. | Reporting only |
| `volkov_extrapolated` | `True` | A species sits outside the Jeans-parameter range 6 to 15 over which the kinetic enhancement factor was measured, so the factor was held at the nearer endpoint. The two sides differ in cost: above 15 the factor is falling toward 1 and holding 1.4 overstates the flux slightly, while below 6 it is rising and holding 1.7 understates it, and the low side is the one a trace light species on a heavy background actually reaches. | The rate reflects it |
| `exobase_not_reached` | `True` | The extended structure never reaches the level where the mean free path equals the scale height, so its top level was used as the exobase. | The rate reflects it |
| `exobase_at_anchor` | `True` | The exobase landed on the profile top itself; one integration interval was kept so the supply integrals exist. | The rate reflects it |
| `extension_unbound` | `True` | The extended upper structure became unbound before the integration finished. | The rate reflects it |

### The screens and the split

| Flag | Value | Meaning | Effect |
|---|---|---|---|
| `roche_overflow` | `True` | The active flow radius reaches the periapsis Hill radius, or the nozzle candidate won the final comparison. The label changes; the rate is the geometric case's branch rate or the nozzle transfer rate, per `rate_branch`. | Reporting only |
| `roche_subflag` | `'dynamical'`, `'no_transonic'`, or `'neither'` | The geometry under either route into the label, read against the Roche lobe rather than the Hill radius because the lobe is the critical surface and sits about 0.70 of the way out to it. `dynamical`: the atmosphere itself reaches the lobe, its outer extent being `r_atmosphere` against `r_lobe` in `diagnostics['nozzle']`. `no_transonic`: it does not, but the flow radius passes the Hill radius. `neither`: the label came from the rate crossing alone. An atmosphere that spills reads `dynamical` whichever candidate carries the rate; which candidate that was is `rate_branch`. | Reporting only |
| `nozzle_saturated` | `True` | The nozzle won with its launch level at or beyond the lobe radius, where the potential expansion has diverged and the exponential is clamped at 1. The rate is the lobe-filling boundary value of the model, a lower bound on the transfer, rather than an interior point of it. | The rate reflects it |
| `nozzle_orbit_averaged` | `True` | The nozzle won on an eccentric orbit, so the dispatched rate is a time average over that orbit, each phase evaluated with the circular formula at its own separation. The primary has no eccentric treatment, so the quasi-static evaluation is this module's convention. Periapsis is the richest phase while the barrier is unclamped and the poorest once `nozzle_saturated` is raised, where the rate goes as the cube of the separation. | The rate reflects it |
| `nozzle_partial_orbit` | `True` | The overflow description holds only on an arc around periapsis, so the average is duty-cycled over that arc and `applicable_orbit_fraction` is below one. The rate therefore omits the wind the planet drives on the rest of the orbit, which one dispatched rate cannot also carry. | The rate reflects it |
| `near_roche` | `True` | The Hill radius is less than 1.5 flow radii, that is, the flow reaches beyond two thirds of the way to the lobe. The tidal factor is steep there; read `diagnostics['roche']['tidal_inflation']` for how much of the rate it is setting. Not raised on a state already labeled `roche_overflow`, since the warning is about the tidal inflation of a bound rate. | Reporting only |
| `t_exo_floored_to_profile_top` | `True` | The prescribed exobase temperature was below the profile's top level, which would build a thermosphere that cools with height and an exobase more strongly bound than its own anchor. The temperature floors at the anchor. In a coupled run the profile top can warm past a fixed prescription over secular time, so this is reported rather than raised. | Reporting only |
| `luminosity_capped` | `True` | The interior luminosity is the term setting the bolometric rate. The cap applies only once the activation gate has closed, so a state crossing that threshold drops discontinuously, by a factor of thousands, while keeping the same label; this flag is how you see which term won. Reachable only when the residual is admitted, since the cap applies past the gate and only an admitted residual is dispatched there. | Reporting only |
| `k_tide_undefined` | `True` | The Hill radius is inside the planetary radius, so the tidal barrier is gone and the factor has no value. The rates are computed without the tidal reduction, which is the smaller reading, and the Roche screen relabels the state. | Reporting only |

Every flag describes the branch whose rate was dispatched. The hydrodynamic candidates are computed on every call, because the diagnostics report them at every dispatch, so their cautions (`subcritical_sonic`, `thermostat_clamped`, `efficiency_fallback_fixed`, and the fitted-efficiency guards) appear only when a hydrodynamic branch won, and the bolometric residual clears them again when it displaces that winner. What the losing candidate did is in `diagnostics['hydrodynamic']` either way.
| `split_from_base_composition` | `True` | The per-species split used the atomized wind-base composition because no reservoirs were supplied. | The split reflects it |
| `rock_former_bij` | list | Rock-forming species (Na, Mg, Si, Fe) are present in the closure, whose binary-diffusion coefficients for them sit in the widest provenance class. | Reporting only |
| `stale_input` | `True` | The caller passed `atm_converged=False`, so the profile is from a non-converged atmosphere solve. | Reporting only |

---

## Diagnostics

Nineteen groups on a typical call. Nothing in the dispatch control flow reads any of them, and there is no option to switch them off: the regime boundaries carry genuine physical uncertainty, and reporting the translation quantities beside every verdict is how the framework handles that.

| Group | Key contents | What it answers |
|---|---|---|
| `knudsen` | `kn_sc`, `threshold_applied`, `sigma_c`, `provenance`, `counterfactual_labels` | Which side of the collisionality switch the state fell on, how close it sat, what the label would have been at both edges of the criterion band, and where the cross sections came from. |
| `hydrodynamic` | `mdot_el`, `mdot_rr`, `efficiency`, `K_tide`, `T_wind`, `selection_mechanism`, `rr_chain` | Both wind candidates, which one won, the temperature the thermostat returned, and the full recombination-limited chain (sound speed, sonic radius, base Jeans parameter, densities, barometric factor). `selection_mechanism` is one of `EL-selected`, `RR-selected`, or `RR-selected:subcritical-floor`, and says which candidate won, not why it was small; the barometric factor is what answers that. |
| `hydrostatic` | `rate_kg_s`, `convergence`, `T_exo`, `T_exo_mode`, `r_exo`, `f_plus_exo`, `T_esc_neutral`, `T_esc_plasma`, `gate`, `gate_unstable`, `detail` | The exobase state, both escape temperatures with the local ionization fraction, which convention gated the branch, the `detail['dominant']` species that supplies itself without a diffusion cap, and per-species Jeans and diffusion fluxes in `detail['species']`. The `convergence` entry records the quadrature levels the supply integrals ended on, the last relative change in the bulk rate and the worst one across the species, and whether the target was met; a call that hits the ceiling first reports `converged` false and the rate it reached. |
| `bolometric` | `T_wind`, `c_s`, `R_sonic`, `x`, `mach`, `mdot_parker`, `mdot_bondi`, `mdot_luminosity`, `k_tide`, `active`, `residual_mode`, `competes`, `p_launch`, `tau_launch`, `rate_kg_s` | The bolometric candidate in full: each cap separately, the tidal factor the luminosity cap divided by, whether the branch was active, and whether the candidate took part in the final comparison. `residual_mode` echoes the setting, and `competes` is true while the activation gate is open and, past it, only when the residual is admitted, so a candidate rate above the dispatched one beside `competes` false was never a contender rather than a loser. `tau_launch` is the plane-parallel optical depth $\kappa P / g$ of the launch level to the opacity you supplied. The Parker rate is derived from a photosphere, so a `tau_launch` far from 1 says the prescribed level and the opacity do not describe the same surface, and the rate is being read off its own definition at the wrong place. It is reporting only, and the level stays prescribed rather than solved for because the activation threshold that gates the branch is calibrated at a level of its own; solving here would put the gate and the rate on two different surfaces. |
| `lambda_gate` | float | The restricted Jeans parameter that decides boil-off activation. |
| `thermostat` | heating and cooling terms, `clamped` | How the wind temperature was reached, channel by channel. |
| `closure` | `active_set`, `retained`, `inv_H_bar_cgs`, `mass_conservation_rel`, `b_provenance` | How the fractionation closure partitioned a wind, present only where a hydrodynamic branch produced the rate and fractionation is on, including under the `roche_overflow` label when that branch won. `inv_H_bar_cgs` is the inverse of the density scale height every escaping gas shares, in cm⁻¹, which the closure solves for alongside the drifts. |
| `roche` | `R_hill_periapsis`, `flow_radius`, `xi_flow`, `xi_ktide`, `k_tide`, `tidal_inflation`, `r_atmosphere`, `rate_branch` | The overflow screen in full: which radius was tested and against what, how far the atmosphere itself reaches, which branch the rate came from, always one of the regime labels (`roche_overflow` itself when the L1 transfer won, and the producing branch when the screen renamed a bound state), and the tidal factor with the factor by which it is raising the rate. A large `tidal_inflation` means the rate is set by the divergence of the factor at the lobe rather than by the heating. |
| `nozzle` | `rate_kg_s`, `rate_full_orbit_kg_s`, `rate_periapsis_kg_s`, `rate_apoapsis_kg_s`, `applicable`, `applicable_orbit_fraction`, `saturated_orbit_fraction`, `R_sonic`, `R_L1`, `R_sonic_over_R_L1`, `R_sonic_over_R_L1_apoapsis`, `n_phase`, `temperature_mode`, `r_launch`, `rho_launch`, `T_K`, `mu_kg`, `v_th`, `q`, `A`, `a_periapsis`, `r_lobe`, `phi_L1`, `phi_ph`, `delta_phi`, `exponent_applied`, `area_m2`, `saturated`, `power_lift_W`, `power_lift_full_orbit_W`, `L_int_W`, `L_bol_intercepted_W` | The Jackson et al. (2017) L1 transfer candidate at every call, whether it won or not. `rate_kg_s` is what the dispatcher competes: the orbit average duty-cycled over the arc where the overflow description applies, which is where the isothermal sonic radius reaches the L1 distance. `rate_full_orbit_kg_s` is the same average without that gate, so it stays comparable with the primary's own published rates, and the periapsis and apoapsis rates bracket the orbit. The geometry entries (the lobe radius, both potentials, the applied exponent, the nozzle area, `saturated`) are reported at periapsis, the tightest geometry of the orbit. The flow carries no energy cap, so the lift power against `L_int_W` and `L_bol_intercepted_W` is what shows where the isothermal assumption is strained: `power_lift_W` pairs with the competed rate and `power_lift_full_orbit_W` with the unguarded one. Both are built from the barrier the rate applied plus the acceleration to the sonic speed, which is the heat an isothermal flow demands, so they stay finite and non-trivial at saturation where the barrier is gone and the acceleration is not. |
| `johnson_q` | `q_net_over_qc`, `q_net_W`, `q_c_W` | Whether the absorbed power can drive a transonic outflow at all, independently of any rate formula. |
| `guo_triple` | `lambda_exo`, `lambda_rp`, `lambda_star`, `thresholds` | The verdict translated into the Jeans-parameter taxonomy. |
| `potential_screens` | `log_minus_phi_cgs`, `caldiroli_threshold`, `above_caldiroli`, `salz_screen`, `salz_verdict`, `salz_attribution` | The verdict translated into threshold-potential taxonomies: the specific binding energy in erg g⁻¹ against the efficiency-collapse band of Caldiroli et al. (2022) and against the wind-versus-thermosphere thresholds of Salz et al. (2016), 13.11 and 13.6. The Salz simulations are hydrogen-dominated, which the reported attribution states, so treat the verdict as out of scope on a heavy secondary atmosphere. |
| `erkaev_tc_K` | float | The tidally corrected critical exobase temperature above which the thermosphere blows off. |
| `fluid_check` | `levels_checked`, `worst_kn`, `fluid`, `truncated_at_profile_top` | Whether the fluid condition holds everywhere below the sonic surface, not only at it, with the truncation declared. |
| `tang_timescale` | boil-off termination timescales | A consistency check on the bolometric rate's own exponential shutoff. |
| `self_consistency` | `evaluated`, `t_deplete_s`, `age_s`, `inconsistent` | Whether the dispatched rate would have destroyed the supplied inventory within the supplied age. Reports `evaluated: False` without an age or reservoirs. |
| `rate_floor` | `floor_kg_s`, `above_floor` | Whether the dispatched rate has any numerical content, against one proton per Julian year. Reporting only: the module never applies the floor. |
| `base_level` | `p_Pa`, `p_physical_Pa`, `r_m`, `T_K`, `clamp_decades` | Where the wind was launched, and the pressure the base method asked for before any clamp. |
| `documentation` | criterion bands and published exponents | The bands themselves, so a stored result is self-describing: the collisionality band, the boil-off activation band, the numerical against analytic flux exponents of the wind limits, and the dayside-heating reduction factors. |
| `contested_ion` | both branch rates and a note | Present only when the two escape-temperature conventions disagree, which turns on ion physics this version does not model. |

---

## Two conventions worth adopting

**A rate floor.** The framework computes what the physics gives it, including rates like $10^{-123}$ kg s⁻¹ from a strongly bound heavy atmosphere. One proton crossing the planet's surface per year, about $5.3 \times 10^{-35}$ kg s⁻¹, is the smallest rate with physical content; below that, report no escape. The convention belongs to the caller, and the module does not apply it, but it does report it: `diagnostics['rate_floor']` carries the number and whether this call cleared it, so a caller need not keep its own copy. Read it before trusting a regime label on a slow state, because a label decided by the ordering of two rates far below the floor is decided by nothing.

**A relevance test, separately.** Clearing the floor does not make a rate matter: a hundred decades above it can still be grams per year. Use `diagnostics['self_consistency']`, which divides the supplied inventory by the dispatched rate and compares against the supplied age, as the yardstick for whether a rate is worth carrying.
