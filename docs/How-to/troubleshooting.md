# Troubleshooting the dispatcher

A regime verdict looks wrong, or a flag fired and you want to know whether it matters. It is the shortest path from a symptom to the thing to look at. Every flag and diagnostics key is defined in the [dispatch results reference](../Reference/results.md), and the physics is on the [escape regimes](../Explanations/regimes.md) page; this page sits between them and starts from what you observed.

Start with the two questions that dispose of most surprises:

1. **Is the rate physically meaningful at all?** One proton per year through the surface, about $5.3 \times 10^{-35}$ kg s⁻¹, is the smallest rate with content. Below it, the label is describing an outflow that does not exist.
2. **Was the state near a boundary?** `diagnostics['knudsen']['counterfactual_labels']` gives the label at both edges of the collisionality criterion. When those disagree with each other, the verdict is a choice the criterion made, not a measurement, and the rest of this page is about which choice.

---

## The rate is zero, or absurdly small

**Symptom.** A `hydrostatic` verdict returns $10^{-70}$ kg s⁻¹ or smaller.

**Cause.** Nothing is wrong. Jeans escape depends exponentially on the Jeans parameter at the exobase, and a heavy species on a strongly bound planet sits at a Jeans parameter of hundreds. Read `diagnostics['hydrostatic']['detail']['species'][name]['lambda_exo']`: above roughly 30, the exponential has already taken the rate out of physical relevance.

**What to do.** Report no escape rather than the number. If you expected loss here, the missing physics is probably not thermal: the nonthermal channels that actually remove heavy species from a real exosphere are absent everywhere in ZEPHYRUS, which is what `hydrostatic_lower_limit` says on every hydrostatic result. A trace of a light species changes the picture completely, because it escapes at the diffusion-limited supply rate rather than its own Jeans rate.

## The bulk rate barely responds to the exobase temperature

**Symptom.** A hydrostatic verdict whose rate is insensitive to `T_exo_value`, against the documented exponential sensitivity.

**Cause.** The species carrying the rate is supply limited, not Jeans limited. The species that carries no cap at all is `diagnostics['hydrostatic']['detail']['dominant']`, which supplies itself. Compare `phi_jeans` against `phi_diffusion` in `diagnostics['hydrostatic']['detail']['species'][name]`: when the diffusion-limited supply is the smaller of the two, the harmonic mean of the two sits near the supply and the exobase temperature has little left to do.

**What to do.** Nothing, but do not report the insensitivity as a general property of the branch: the heavy species in the same result will be moving by orders of magnitude over the same temperature range.

## `thermostat_clamped` fired

**Symptom.** The wind temperature sits at a bracket edge, usually $5 \times 10^{4}$ K.

**Cause.** The local heating against cooling balance had no root in the bracket. At a dense wind base this is physical: electron densities far above the critical densities of the forbidden lines quench the line coolants collisionally, so nothing balances the heating and the wind runs hot.

**What to do.** Read `diagnostics['hydrodynamic']['T_wind']` and decide whether that temperature is one you are willing to carry. It propagates: the sound speed, the barometric exponent, and the recombination coefficient all depend on it, so it moves the boundary between the two hydrodynamic sub-labels. The thermostat evaluates one level and does not model the temperature structure through the sonic region, which is the limitation behind the clamp.

## `subcritical_sonic` fired

**Symptom.** The flag on a hydrodynamic verdict, sometimes with the label `hydrodynamic:RR`.

**Cause.** The isothermal sonic radius $R_\mathrm{s} = G M_\mathrm{p} / (2 c_\mathrm{s}^2)$ came out below the wind base, which happens when the base Jeans parameter $\lambda_\mathrm{b} = G M_\mathrm{p} / (R_\mathrm{base} c_\mathrm{s}^2)$ drops below 2. The gas at the base is then already unbound in the isothermal sense: there is no subsonic region above the base for a wind to accelerate through, so there is no transonic point to anchor a rate on. Below $\lambda_\mathrm{b} = 3/2$ it is worse, because the barometric factor $e^{3/2 - \lambda_\mathrm{b}}$ exceeds 1 and the formula would carry density outward rather than thinning it. The branch responds by flooring the sonic radius at the base and dropping the barometric factor, which keeps the rate finite and monotone, and flags that it did.

Physically, a state in that corner is not an XUV wind launched at an ionization front. It is the blow-off corner, where the whole upper atmosphere streams away, the transonic point sits below the base, and the rate is set by the energy supplied rather than by radiation. The framework has a branch for that, tested before this one on the restricted Jeans parameter, and most states that raise this flag are labeled `boiloff` or `roche_overflow` for that reason.

**What to do.** Read `diagnostics['hydrodynamic']['rr_chain']['lambda_b']` and `diagnostics['hydrodynamic']['selection_mechanism']` together, then split by which candidate won:

- `EL-selected`: the flag is informational. The energy-limited rate is an energy budget and never used the sonic point, so it survives the condition intact.
- `RR-selected:subcritical-floor`: the number is a floored placeholder rather than a rate, because the minimum compared a real energy-limited rate against a formula evaluated outside its own domain. Compare against the bolometric candidate in `diagnostics['bolometric']` before using it, and expect the label to be one of the extreme ones.

## The label is `hydrodynamic:RR` and you want to know why

**Symptom.** A recombination-limited verdict, and the question of whether recombination is what limited it.

**Cause.** The minimum over the two hydrodynamic candidates can select the recombination-limited rate for two physically different reasons. One is genuine saturation: the base ion density follows photoionization against recombination balance, so the rate grows as $\sqrt{F_\mathrm{XUV}}$ while the energy-limited rate grows linearly, and the former wins at high flux. The other is barometric suppression: at large $\lambda_\mathrm{b}$ the isothermal wind is exponentially throttled between the base and the sonic point, and the smallness has nothing to do with radiation.

**What to do.** Read `diagnostics['hydrodynamic']['rr_chain']['barometric_factor']`, which is $e^{3/2 - \lambda_\mathrm{b}}$, the fraction of the base density that survives to the sonic point. Near 1 the sonic-point density is the base density and the rate is set by the recombination-limited base ionization, which is the mechanism the label names. Several decades below 1 the rate is small mostly because the isothermal wind cannot carry material that far, and attributing the smallness to recombination misplaces it.

Do not expect `selection_mechanism` to answer this. It reports which candidate won, `EL-selected`, `RR-selected`, or `RR-selected:subcritical-floor`, and nothing about why. An earlier version split RR wins on the base Jeans parameter; that split was retired because the number behind it had no source and it placed the canonical recombination-limited case of the literature, near $\lambda_\mathrm{b} = 5.5$, on the suppression side.

The flux scaling will not separate them either: the base ion density follows $\sqrt{F_\mathrm{XUV}}$ at every $\lambda_\mathrm{b}$ in this chain, so a fitted exponent sits near 1/2 whichever case holds, departing from it only where the thermostat is still moving the wind temperature quickly with flux. Use the factor, and do not describe a strongly suppressed point as recombination limited in text or in a figure legend.

## `base_clamped` fired

**Symptom.** The wind base sits at the profile top, with a clamp distance of several decades.

**Cause.** The profile does not extend to the pressure where XUV photons are absorbed, near a nanobar. An isothermal hydrogen envelope, in particular, becomes unbound and truncates well below that.

**What to do.** The base level is $\mu g / \sigma_{\nu_0}$, so where it sits depends on the state: tens of nanobars on an Earth-mass carbon dioxide planet, and below a nanobar on a low-gravity hydrogen envelope, where a profile stopping at a nanobar still misses it. Extend the profile past the `p_physical` the diagnostics report for the state at hand rather than trusting one number for all of them. Otherwise decide by branch: on `boiloff` the clamp is harmless, because that branch launches from the photospheric level. On a hydrodynamic verdict it is not, because the base density sets the sonic-point density and therefore the collisionality switch itself, so a clamped base moves the boundary. The `base_out_of_range = 'extend'` setting evaluates the base on the extended upper structure instead; `base_extension_truncated` means even that did not reach it.

## `roche_overflow` on a state that should not be overflowing

**Symptom.** The overflow label on a small or quiescent planet, sometimes at a rate that is negligible in absolute terms, and sometimes flipping to another label under a small change in an input or a setting.

**Cause.** The label reports a geometry, not a rate. The screen tests the flow radius of the branch that won the rate comparison, and when the bolometric residual wins, which it can do only where `residual = 'luminosity_capped'` admits it, that radius is its sonic radius $R_\mathrm{B} = \Lambda R_\mathrm{p} 2^{1/4} / 2$, which grows with the Jeans parameter. A tightly bound heavy atmosphere therefore puts it several Hill radii out while the atmosphere itself sits deep inside, which is the opposite of overflowing. The flipping under a small change is the other half: the label boundary is a comparison between two rate candidates, so wherever the two are within a few percent of each other, a small input change moves the label. The rate does not move with it, because the screen never changes the rate.

**What to do.** Read four things: `diagnostics['roche']['rate_branch']` (which branch the rate came from; `roche_overflow` means the L1 transfer itself was dispatched, and any other value under that label means the geometric screen renamed a state whose rate that branch produced), `r_atmosphere` against `R_hill_periapsis` in the same group (does the atmosphere itself reach the lobe, which is what `roche_subflag` reports as `dynamical`), `diagnostics['rate_floor']['above_floor']` (is the rate a number at all), and `flags['bolometric_residual']`.

The four cases that produces:

- `rate_branch` reading `roche_overflow`: the rate is the tidally driven transfer through the inner Lagrange point (Jackson et al. 2017), dispatched because it beat every bound candidate inside its applicability criterion. It is a real transfer rate rather than a rename; `diagnostics['nozzle']` carries the barrier, the saturation state, and the lift power against the available luminosities. The subflag beside it reads `dynamical` when the atmosphere's own extent passes the lobe and `neither` when it does not, so it describes the geometry and not the mechanism. On an eccentric orbit the rate is an orbit average duty-cycled over the arc where the description applies; read `applicable_orbit_fraction`, and treat a value below one as a rate that omits the wind driven on the rest of the orbit.
- Subflag `dynamical` on a bound branch and a rate above the floor: the atmosphere reaches its Roche lobe while the nozzle candidate lost or sat outside its criterion, the rate is the bound-flow estimate, and the real rate is higher by whatever the tidal flow would carry. Treat it as a lower limit.
- Subflag `no_transonic` with `r_atmosphere` well inside the Hill radius: only the would-be sonic surface passes the lobe. On a heavy bound atmosphere this is the tightly-bound case above, not an overflow; the [regimes page](../Explanations/regimes.md) explains why the two look alike to the screen.
- `above_floor` false: the label is decided by the ordering of two rates with no numerical content. Report no escape and treat the geometry as a note.
- A rate above the floor but far below anything that matters. The floor marks what is distinguishable from zero in floating point, one proton per Julian year, which is many decades below a rate that could change an atmosphere. Around one nozzle win in twenty is dispatched below 1e-20 kg s⁻¹, which clears the floor and means nothing physically. Read `diagnostics['self_consistency']` beside it: what settles whether a rate matters is the depletion timescale against the age you passed in, not the floor.

Points whose flow radius exceeds two thirds of the Hill radius raise `near_roche` instead of the label, and the tidal factor is steep there, so a rate from that region carries the tidal correction's sensitivity with it.

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

**What to do.** Check the label and the `fractionate` setting, and look for `split_from_base_composition`, which says the split fell back to the wind-base composition because no reservoirs were supplied. A 3% departure from the mass fractions is a real result for a well-coupled heavy wind, not a sign the closure failed to run; `diagnostics['closure']` carries the active set and the mass residual if you need to confirm it did.

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
