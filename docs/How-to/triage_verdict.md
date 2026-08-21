# Triage a dispatch verdict

A regime verdict looks wrong, or a flag fired and you want to know whether it matters. This page is the shortest path from a symptom to the thing to look at. Every flag and diagnostics key is defined in the [dispatch results reference](../Reference/results.md), and the physics is on the [escape regimes](../Explanations/regimes.md) page; this is the triage layer between them.

Start with the two questions that dispose of most surprises:

1. **Is the rate physically meaningful at all?** One proton per year through the surface, about 5.3e-35 kg s⁻¹, is the smallest rate with content. Below it, the label is describing an outflow that does not exist.
2. **Was the state near a boundary?** `diagnostics['knudsen']['counterfactual_labels']` gives the label at both edges of the collisionality criterion. When those disagree with each other, the verdict is a choice the criterion made, not a measurement, and the rest of the triage is about which choice.

---

## The rate is zero, or absurdly small

**Symptom.** A `hydrostatic` verdict returns 1e-70 kg s⁻¹ or smaller.

**Cause.** Nothing is wrong. Jeans escape depends exponentially on the Jeans parameter at the exobase, and a heavy species on a strongly bound planet sits at a Jeans parameter of hundreds. Read `diagnostics['hydrostatic']['detail']['species'][name]['lambda_exo']`: above roughly 30, the exponential has already taken the rate out of physical relevance.

**What to do.** Report no escape rather than the number. If you expected loss here, the missing physics is probably not thermal: the nonthermal channels that actually remove heavy species from a real exosphere are absent everywhere in ZEPHYRUS, which is what `hydrostatic_lower_limit` says on every hydrostatic result. A trace of a light species changes the picture completely, because it escapes at the diffusion-limited supply rate rather than its own Jeans rate.

## The bulk rate barely responds to the exobase temperature

**Symptom.** A hydrostatic verdict whose rate is insensitive to `T_exo_value`, against the documented exponential sensitivity.

**Cause.** The species carrying the rate is supply limited, not Jeans limited. The species that carries no cap at all is `diagnostics['hydrostatic']['detail']['dominant']`, which supplies itself. Compare `phi_jeans` against `phi_diffusion` in `diagnostics['hydrostatic']['detail']['species'][name]`: when the diffusion-limited supply is the smaller of the two, the harmonic mean of the two sits near the supply and the exobase temperature has little left to do.

**What to do.** Nothing, but do not report the insensitivity as a general property of the branch: the heavy species in the same result will be moving by orders of magnitude over the same temperature range.

## `thermostat_clamped` fired

**Symptom.** The wind temperature sits at a bracket edge, usually 5e4 K.

**Cause.** The local heating against cooling balance had no root in the bracket. At a dense wind base this is physical: electron densities far above the critical densities of the forbidden lines quench the line coolants collisionally, so nothing balances the heating and the wind runs hot.

**What to do.** Read `diagnostics['hydrodynamic']['T_wind']` and decide whether that temperature is one you are willing to carry. It propagates: the sound speed, the barometric exponent, and the recombination coefficient all depend on it, so it moves the boundary between the two hydrodynamic sub-labels. The thermostat evaluates one level and does not model the temperature structure through the sonic region, which is the limitation behind the clamp.

## `subcritical_sonic` fired, and the label is `hydrodynamic:RR`

**Symptom.** A recombination-limited verdict carrying this flag.

**Cause.** The computed sonic radius fell below the wind base, so it was floored at the base. The recombination-limited rate then wins the minimum through barometric suppression between base and sonic point, not through recombination saturation.

**What to do.** Check `diagnostics['hydrodynamic']['selection_mechanism']`, which names the mechanism. Do not describe such a point as recombination limited in text or in a figure legend: it is the same label reached by different physics, and the distinction is the reason the mechanism is reported.

## `base_clamped` fired

**Symptom.** The wind base sits at the profile top, with a clamp distance of several decades.

**Cause.** The profile does not extend to the pressure where XUV photons are absorbed, near a nanobar. An isothermal hydrogen envelope, in particular, becomes unbound and truncates well below that.

**What to do.** For a production profile, set the top pressure below 1 nanobar and the clamp never engages. Otherwise decide by branch: on `boiloff` the clamp is harmless, because that branch launches from the photospheric level. On a hydrodynamic verdict it is not, because the base density sets the sonic-point density and therefore the collisionality switch itself, so a clamped base moves the boundary. The `base_out_of_range = 'extend'` setting evaluates the base on the extended upper structure instead; `base_extension_truncated` means even that did not reach it.

## `roche_overflow` on a state that should not be overflowing

**Symptom.** The overflow label on a small or quiescent planet, sometimes at a rate that is negligible in absolute terms, and sometimes flipping to another label under a small change in an input or a setting.

**Cause.** The label names the branch that won the rate comparison, and the screen tests that branch's own flow radius. When the bolometric residual wins, the radius tested is its sonic radius, which for a warm low-gravity atmosphere can be larger than the Hill radius even when the residual rate is tiny. So the label reports a geometry, and the rate it carries can be small.

**What to do.** Read three things: `flags['bolometric_residual']` (did the residual take the label), `diagnostics['roche']` (`flow_radius` against `R_hill_periapsis`), and `diagnostics['bolometric']` (which cap set the rate). If the residual won at a rate below the floor, treat the point as no escape and the geometry as a note rather than a result. Points inside 1.5 Hill radii raise `near_roche` instead of the label, and the tidal factor is steep there, so a rate from that region carries the tidal correction's sensitivity with it.

## `contested_ion` fired

**Symptom.** The flag, plus a `diagnostics['contested_ion']` group holding two rates.

**Cause.** The neutral escape temperature and the plasma escape temperature, which is half of it because an ambipolar field shares each ion's binding energy with its electron, disagree about whether the exosphere can stay hydrostatic. The physics that would decide it is ion outflow, which this version does not model.

**What to do.** Report both rates, or report the point as contested. Do not resolve it with the `gate` setting: that setting chooses which convention gates the branch, and choosing one does not make the other wrong. In the nitrogen and oxygen corner the two conventions can differ by four orders of magnitude in rate.

## The label flips between neighboring points of a smooth sweep

**Symptom.** A parameter grid or a time series whose label alternates.

**Cause.** Either the state is sitting inside the collisionality band, where the criterion does not decide the label, or the two candidate rates are within a hair of each other and the minimum keeps changing hands.

**What to do.** For an evolutionary track, supply `prev_regime` from the previous step. That opens a hysteresis window around the threshold, so the previous label wins inside it and a track cannot chatter on numerical noise; `hysteresis_active` confirms it was in use and `diagnostics['knudsen']['threshold_applied']` shows the shifted threshold. For a static grid, quote the boundary with the width its criterion implies rather than as a line: re-dispatch at `kn_crit` of 0.1 and 3 and report the range. The [tutorial](../Tutorials/dispatch.md) measures such a band as a factor of 4.3 in boundary flux for one planet.

## The per-species rates look like a bulk split

**Symptom.** Element shares that match their reservoir mass fractions.

**Cause.** Either fractionation is off, or the verdict is not a confirmed wind. The closure runs only on a hydrodynamic verdict with `fractionate` on; every other branch splits by reservoir mass fractions, except the hydrostatic branch, which is natively per-species.

**What to do.** Check the label and the `fractionate` setting, and look for `split_from_base_composition`, which says the split fell back to the wind-base composition because no reservoirs were supplied. A three percent departure from the mass fractions is a real result for a well-coupled heavy wind, not a sign the closure failed to run; `diagnostics['closure']` carries the active set and the mass residual if you need to confirm it did.

## The rate is large and the planet should not have survived

**Symptom.** A plausible-looking rate on a grid point whose atmosphere would be long gone.

**Cause.** Nothing in the module knows the state's history. `diagnostics['self_consistency']` is the check: it divides the supplied inventory by the dispatched rate and compares against the supplied age.

**What to do.** When `inconsistent` is `True`, the grid point is describing a state that cannot have persisted to the age it claims. That is a statement about the grid, not the rate. Report the flagged fraction of a static grid rather than dropping the points silently.

## The fitted efficiency returns something implausible

**Symptom.** `efficiency_mode = 'caldiroli'` returning a value near one, with `caldiroli_out_of_box`.

**Cause.** The fit was made on sub-Neptunes through hot Jupiters. A one Earth-mass planet sits below the gravitational potential range it covers, so the value is an extrapolation, and the flag says so rather than silently refusing.

**What to do.** Below the fit's flux bound the value is rejected outright and `caldiroli_below_flux_bound` plus `efficiency_fallback_fixed` tell you the fixed setting was used instead. Outside the potential box, prefer a swept fixed efficiency over the extrapolated fit, and report the sweep range rather than one value.

---

## When none of the above applies

Print the whole container for the offending call and read it in the order the [tutorial](../Tutorials/dispatch.md) uses: which branch and how close, what set the rate, whether the heating could drive a flow at all, how the verdict translates into other taxonomies, and whether the snapshot is self-consistent. If the verdict still looks wrong after that, the input state is the next suspect: check the profile spans the pressures the branches need, that the composition is what you meant, and that every scalar is in SI.
