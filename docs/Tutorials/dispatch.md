# Dispatching an escape regime

The [first run](first_run.md) applied one prescription unconditionally: scalars in, an energy-limited rate out. This tutorial uses the full framework instead. You hand `zephyrus.dispatch` one planetary state, it decides which escape physics that state is in, and it returns the rate that physics gives, the per-species split, flags, and the diagnostics that say how close the state sat to each regime boundary.

By the end you will have crossed two regime boundaries on the same planet, produced the figure below, read a diagnostics container field by field, measured how far a boundary moves across the width of its own criterion, and dispatched one atmosphere along a stellar XUV history.

If you have not installed ZEPHYRUS yet, follow the [installation guide](../How-to/installation.md) and run `mors download all` for the stellar tracks. The [escape regimes](../Explanations/regimes.md) page defines the physics behind every step here; this page is about driving it.

!!! info "What you'll do"
    - Build a planetary state: scalars plus an atmosphere profile
    - Read one verdict field by field, and check the per-species sum
    - Sweep the XUV flux until the regime label changes, twice
    - Meet the boil-off and overflow labels
    - Read the diagnostics container and learn which entries to look at first
    - Move a boundary with each of the four knobs that move one
    - Dispatch a frozen atmosphere along a stellar history

---

## The full script

The complete example is `examples/demo_dispatcher/demo_dispatcher.py`. It writes one figure, so create the output directory first:

```sh
mkdir -p output
python examples/demo_dispatcher/demo_dispatcher.py
```

It prints a table per step and runs in a few seconds, most of which is the boundary bisections. Every step is a function that returns its results, so you can import one at a time:

```python
from examples.demo_dispatcher.demo_dispatcher import flux_sweep, read_diagnostics
```

The rest of this page walks the script step by step.

---

## Step 1: build a state

The energy-limited entry point takes scalars. The framework needs the atmospheric structure as well, because the quantities that decide the regime are properties of the structure and not of the surface: the pressure level where XUV photons are absorbed and the wind is launched, the level where the gas stops colliding often enough to behave as a fluid, and the exobase where individual particles start escaping ballistically.

So a state is scalars plus a `profiles.Profile`:

```python
from zephyrus.dispatcher import DispatchSettings, EscapeInputs, dispatch
from zephyrus.profiles import isothermal_profile

profile = isothermal_profile(M_p, R_p, T_eq, {'CO2': 1.0}, 1.0e7, 1.0e-5)

state = EscapeInputs(
    M_p=M_p,             # planet mass                              [kg]
    R_p=R_p,             # planet radius                            [m]
    M_star=Ms,           # stellar mass                             [kg]
    a=a,                 # semi-major axis                          [m]
    e=0.0,               # eccentricity
    T_eq=T_eq,           # equilibrium temperature                  [K]
    F_xuv=10.0,          # XUV flux at the planet                   [W m-2]
    F_bol=f_bol,         # bolometric instellation                  [W m-2]
    F_int=1.0,           # interior heat flux                       [W m-2]
    kappa_photo=0.01,    # photospheric opacity                     [m2 kg-1]
    profile=profile,
    settings=DispatchSettings(),
    age=age,             # optional, for the consistency screen     [s]
    reservoirs=reservoirs,  # optional, element inventories          [kg]
)
```

In a coupled run the atmosphere module supplies the profile. Standalone, `isothermal_profile` integrates a hydrostatic isothermal structure of fixed composition, which is enough to exercise every branch. Three choices in the script are worth stating, because they are the ones that change results:

- `p_top = 1e-5` Pa, which is 0.1 nanobar. The XUV wind launches near a nanobar, so a profile that stops deeper than that cannot reach its own wind base and the base clamps to the profile top instead, flagged. Setting the top below a nanobar keeps the clamp out of the way.
- `kappa_photo = 0.01` m² kg⁻¹, about 0.1 cm² g⁻¹. The boil-off rate scales as its inverse, so it matters whenever the bolometric branch is in play.
- The orbit is 0.0775 au around a solar-luminosity star, which puts the zero-albedo, full-redistribution equilibrium temperature at 1000 K. Deriving $T_\mathrm{eq}$ from the orbit rather than setting both by hand keeps the state self-consistent.

!!! warning "SI at every boundary"
    Every input is SI: kilograms, meters, seconds, kelvin, W m⁻², m² kg⁻¹. MORS returns cgs luminosities, so the flux conversions in step 7 are explicit.

---

## Step 2: one call, one verdict

```python
result = dispatch(state)
```

For the CO₂ planet at 10 W m⁻², the script prints:

```text
regime      hydrodynamic:EL
rate        2.3477e+06 kg/s   (7.409e+13 kg/yr)
C           6.6097e+05 kg/s
O           1.6867e+06 kg/s
sum         2.3477e+06 kg/s, closure error 0.00e+00
flags       none
diagnostics 17 groups
```

Five fields, and each has a contract. `regime` is one of five labels. `mdot` is a non-negative bulk rate in kg s⁻¹. `per_species` gives element rates that sum to `mdot` at machine precision, which is the property a coupled run relies on when it debits element inventories, and which is worth asserting in your own code. `flags` records every clamp, fallback, and screen that fired; an empty dictionary means nothing needed reporting. `diagnostics` is the container of step 5.

Every physically posed state returns a result. Exceptions are reserved for malformed input, so a `ValueError` from `dispatch` means the state itself is wrong (a negative mass, an eccentricity of 1, a profile whose pressure does not decrease outward), not that the physics failed.

---

## Step 3: cross a boundary

The interesting thing about a dispatcher is where it changes its mind. Sweep the XUV flux on a fixed planet and the label moves, because the flux at a fixed orbit stands in for stellar age: a young star delivers orders of magnitude more XUV than the same star does later.

```python
for f_xuv in np.logspace(-2, np.log10(5e3), 30):
    result = dispatch(build_state('CO2', 1.0, 1.0, float(f_xuv)))
    print(f_xuv, result.regime, result.mdot,
          result.diagnostics['knudsen']['kn_sc'])
```

On the CO₂ planet the sweep crosses two boundaries:

```text
    F_xuv  regime             rate [kg/s]      Kn_sc  at Kn 0.1      at Kn 3
     0.01  hydrostatic         4.989e-123    4.5e+18  hydrostatic    hydrostatic
    0.587  hydrostatic         4.989e-123       3.92  hydrostatic    hydrostatic
     2.28  hydrodynamic:EL      5.356e+05      0.118  hydrostatic    hydrodynamic
     34.5  hydrodynamic:EL      8.090e+06     0.0127  hydrodynamic   hydrodynamic
      520  hydrodynamic:EL      1.222e+08    0.00265  hydrodynamic   hydrodynamic
 2.02e+03  hydrodynamic:RR      2.572e+08    0.00131  hydrodynamic   hydrodynamic
boundary: hydrostatic to hydrodynamic:EL at 0.77 W m-2
boundary: hydrodynamic:EL to hydrodynamic:RR at 567 W m-2
```

The first crossing is the collisionality switch. Below 0.77 W m⁻² the heating is too weak to keep the gas collisional where a wind would go sonic, so the sonic-point Knudsen number `kn_sc` exceeds 1 and escape is per-particle from the exobase. Above it a fluid wind exists and the rate is the smaller of the two hydrodynamic limits. The second crossing at 567 W m⁻² is inside the wind: the recombination-limited rate drops below the energy-limited one and names the label.

Run the same sweep on a nitrogen and oxygen atmosphere of the same mass and radius and both facts change: the wind sets in at 0.103 W m⁻², a factor 7.5 lower, and the recombination-limited rate never wins anywhere in the swept range. Composition moves the boundaries, not just the rates.

![Regime labels across an XUV flux sweep for two compositions](../assets/dispatcher_regime_sweep.png)

Two things in that figure deserve attention.

**The shaded band.** The collisionality threshold has a default of 1 and a physical range of 0.1 to 3, because kinetic simulations place the fluid-to-kinetic transition near 0.1 when the heating is deposited in a sharp layer and near 1 when it is distributed[^johnson]. That is heating-geometry physics, not a tuning parameter, and it makes the boundary a band: on the CO₂ planet the wind sets in anywhere between 0.61 and 2.6 W m⁻², a factor of 4.3, depending on which edge of the criterion you adopt. The dashed line is the default; the band is the honest width. Every verdict reports the label it would have carried at both edges, in `diagnostics['knudsen']['counterfactual_labels']`, so you never have to rerun a sweep to find this out.

**The open markers.** The hydrostatic rates on these two planets are 1e-123 and 1e-71 kg s⁻¹. Those are not small rates, they are zero with numerical noise attached: carbon dioxide at one Earth mass sits at an exobase Jeans parameter of 301, and $e^{-301}$ is not a number with physical content. Two yardsticks keep this straight:

- One proton crossing the surface per year, 5.3e-35 kg s⁻¹, is the smallest rate that can mean anything. Below it, report no escape. The module computes what the physics gives and leaves this convention to you.
- Above that floor, ask whether the rate matters, by comparing it against the inventory and the age. The diagnostics already do this: `diagnostics['self_consistency']` divides the reservoirs by the rate and compares against the age you supplied.

A rate can clear the floor by a hundred decades and still be irrelevant. A Mars-mass planet losing 2e-9 kg s⁻¹ loses 66 grams a year.

---

## Step 4: boil-off and overflow

Two labels do not appear in that sweep, and neither depends on the XUV flux.

An inflated hydrogen envelope flows out on the planet's own thermal energy before stellar XUV matters. The test is the restricted Jeans parameter $\Lambda$, the ratio of gravitational binding to thermal energy at the photospheric level, and the state is labeled `boiloff` while it sits below the threshold:

```text
regime boiloff, rate 1.348e+14 kg/s, Lambda 11.1 (activation band 15 to 35)
flow radius 6.314e+07 m against Hill radius 1.160e+08 m, ratio 1.837
flags ['base_clamp_decades', 'base_clamped', 'subcritical_sonic']
```

At 1 Earth mass and 1.5 Earth radii with a hydrogen and helium envelope, $\Lambda = 11.1$, well below the threshold of 20, and the rate is 1e14 kg s⁻¹ at every flux in the sweep. That planet is not long-lived, which is the point: boil-off is the regime of the first few million years.

The clamp flags are expected here, and harmless: an isothermal hydrogen envelope becomes unbound before it reaches a nanobar, so the profile stops early and the wind base clamps to its top. The boil-off branch launches from the photospheric level, not the wind base, so the clamp does not touch the rate. Flags tell you what happened; deciding whether it matters is your job, and the [triage guide](../How-to/triage_verdict.md) is a shortcut for the common cases.

Push the same envelope to 3 Earth masses and 2 Earth radii and the flow stops being bound to the planet at all:

```text
regime roche_overflow, rate 1.038e+09 kg/s, Lambda 25
flow radius 1.894e+08 m against Hill radius 1.673e+08 m, ratio 0.883
```

The flow radius of the winning branch now exceeds the Hill radius, so the atmosphere spills over the gravitational boundary rather than escaping through any of the other regimes, and the label says so. When you see this label, read `diagnostics['roche']` for the two radii and `diagnostics['bolometric']` for the rate that won: the label names the branch that won the rate comparison, and on a marginal case that branch can be the bolometric residual carrying a rate that is negligible in absolute terms.

---

## Step 5: read the diagnostics

A call returns sixteen or seventeen diagnostic groups, depending on which branch ran, and none of them gates anything: the dispatch control flow never reads them back, and there is no switch to turn them off. The regime boundaries carry real physical uncertainty, and reporting the translation quantities beside every verdict is how the framework handles that instead of hiding it.

The [dispatch results reference](../Reference/results.md) documents every key. What follows is the order to read them in, which is the order the questions occur to you.

**Which branch, and how close was the switch?**

```python
kn = result.diagnostics['knudsen']
kn['kn_sc'], kn['threshold_applied'], kn['counterfactual_labels']
# 0.03003, 1, {0.1: 'hydrodynamic', 3.0: 'hydrodynamic'}
```

A Knudsen number of 0.03 against a threshold of 1 is a wind by a factor of 30, and the counterfactuals confirm the verdict survives both edges of the band. A value within a factor of a few of the threshold, or counterfactuals that disagree, means the label is a choice and not a measurement.

**What set the rate?**

```python
hy = result.diagnostics['hydrodynamic']
hy['mdot_el'], hy['mdot_rr'], hy['selection_mechanism'], hy['T_wind']
# 2.348e+06, 1.148e+07, 'EL-selected', 8431.0
```

Both candidates are always computed, so you can see the margin. `selection_mechanism` matters more than it looks: the minimum can select the recombination-limited rate two physically different ways, through genuine recombination saturation or through plain barometric suppression between the wind base and the sonic point, and calling the second one recombination-limited would be a category error. The mechanism string distinguishes them. `T_wind` is the temperature a local heating against cooling balance returned, 8431 K here rather than the canonical 10 000 K, and it feeds the sound speed, the barometric exponent, and the recombination coefficient, so it moves the crossover flux of step 3.

**Could the heating drive the flow at all?**

```python
result.diagnostics['johnson_q']['q_net_over_qc']   # 10.1
```

The transonic energy criterion[^johnson] compares the absorbed, efficiency-degraded power against the power needed to sustain a transonic outflow. A ratio below 1 says the heating cannot drive the flow sonic no matter what a rate formula returns. At 10.1 there is an order of magnitude in hand.

**How does this verdict translate into other taxonomies?**

```python
result.diagnostics['guo_triple']        # lambda_exo 301, lambda_rp 331, lambda_star 304
result.diagnostics['potential_screens'] # log potential 11.796, verdict 'wind'
result.diagnostics['erkaev_tc_K']       # tidally corrected critical exobase temperature
```

Different papers classify escape with different quantities. Reporting the Jeans-parameter triple and the threshold-potential screens beside the label lets a reader who thinks in one taxonomy check the verdict in theirs.

**Is the snapshot self-consistent?**

```python
result.diagnostics['fluid_check']       # worst_kn 0.124 over 120 levels, fluid True
result.diagnostics['self_consistency']  # empties in 6.95e+04 yr against 1.00e+08 yr
result.diagnostics['tang_timescale']    # boil-off termination timescales
```

The fluid condition has to hold everywhere below the sonic surface, not only at it, so the check walks the profile levels and reports the worst local Knudsen number, declaring the truncation at the profile top. The consistency screen is the sharp one here: at 2.3e6 kg s⁻¹ this planet empties one Earth atmosphere in 70 000 years, against the 100 Myr age supplied with the state. The state is not wrong, but it cannot have persisted, and a static grid point that fails this screen is telling you the grid, not the code, needs a second look.

**Where was the wind launched, and where did the coefficients come from?**

```python
result.diagnostics['base_level']        # p 3.312e-03 Pa, physical target the same, no clamp
result.diagnostics['knudsen']['provenance']  # {'C': 'laricchiuta', 'O': 'laricchiuta'}
```

Collision cross sections come from a provenance ladder, from tabulated collision integrals down to a geometric hard sphere whose bias is documented. The provenance travels with the result, so a rate that rests on the last rung says so.

---

## Step 6: turn the knobs

Four settings move a boundary rather than just a rate. All of them are in the [parameter reference](../Reference/parameters.md); what follows is what each one does to the answer.

**The collisionality threshold, across its band.** Bisect the wind boundary at each edge instead of guessing:

```text
kn_crit  0.1: wind sets in at F_xuv = 2.639 W m-2
kn_crit  1.0: wind sets in at F_xuv = 0.7695 W m-2
kn_crit  3.0: wind sets in at F_xuv = 0.6118 W m-2
```

A factor of 4.3 in boundary position, from a criterion whose own range is a factor of 30. That number is a result, not an error bar to hide: quote a regime boundary with it.

**The exobase temperature.** The hydrostatic branch prescribes it, default 1000 K, and the Jeans flux depends on it exponentially, so it is the branch's dominant sensitivity. The lesson is sharper on a planet where hydrostatic escape does physical work: a Mars-mass planet whose carbon dioxide carries one percent hydrogen.

```text
T_exo    500 K: hydrostatic  rate 3.2861e+02 kg/s, H 3.286e+02, C 1.091e-24
T_exo   1000 K: hydrostatic  rate 3.5941e+02 kg/s, H 3.594e+02, C 7.672e-10
T_exo   2000 K: hydrostatic  rate 3.6912e+02 kg/s, H 3.690e+02, C 3.475e-02
```

Over a factor of four in temperature the bulk rate moves by 12 percent while the carbon rate moves by 22 orders of magnitude. The reason is in the per-species detail: hydrogen sits at an exobase Jeans parameter of 1.59 with a Jeans flux of 5.4e15 against a diffusion-limited supply of 2.7e14, so its escape is set by how fast diffusion resupplies it through the heavy background and the exobase temperature barely enters. Carbon and oxygen are Jeans limited and carry the whole exponential. One case, both halves of the harmonic mean that combines them[^yelle], and a warning against reading a bulk rate as though one mechanism produced it.

**Fractionation.** A confirmed wind partitions over species through the closure; anything else splits by reservoir mass fractions.

```text
closure                  C 0.2815, O 0.7185
reservoir mass fractions C 0.2729, O 0.7271
```

A three percent shift, with the lighter element enriched, which is the right size for a well-coupled heavy wind and a reminder that a per-species output is sometimes a bulk split wearing a per-species shape. Where fractionation is strong, this comparison is the whole story of how the atmospheric mean molecular weight evolves.

**The efficiency.** Sweeping the energy-limited efficiency across its literature range moves the rate linearly and can move the sub-label:

```text
epsilon 0.10: hydrodynamic:EL   rate 2.3477e+06 (EL 2.348e+06, RR 1.148e+07)
epsilon 0.60: hydrodynamic:RR   rate 1.1476e+07 (EL 1.409e+07, RR 1.148e+07)
```

At 0.6 the energy-limited candidate overtakes the recombination-limited one and the label changes without the physics of the wind changing at all: the minimum switched hands, nothing else. The fitted-efficiency option returns 0.791 for this planet with `caldiroli_out_of_box` raised, because a one Earth-mass planet sits below the gravitational potential range the fit was made on[^caldiroli]. That is the guard working. Take the flag seriously rather than the number.

**One more, for evolutionary use.** Supply the previous label and a hysteresis window opens around the threshold, so a time-stepping track cannot chatter between branches on numerical noise:

```text
F_xuv 0.747, previously None              -> hydrostatic       (Kn_sc 1.124, threshold 1)
F_xuv 0.747, previously hydrodynamic:EL   -> hydrodynamic:EL   (Kn_sc 1.124, threshold 1.5)
F_xuv 0.793, previously None              -> hydrodynamic:EL   (Kn_sc 0.898, threshold 1)
F_xuv 0.793, previously hydrostatic       -> hydrostatic       (Kn_sc 0.898, threshold 0.667)
```

Inside the window the previous label wins, in both directions, and the applied threshold in the diagnostics tells you when the memory was in use.

---

## Step 7: dispatch along a stellar history

Nothing so far needed a star with a history. Load one, derive the XUV and bolometric fluxes from its track, and dispatch the same frozen atmosphere at each age:

```python
star = mors.Star(Mstar=1.0, Omega=1.0)
f_xuv = (star.Tracks['Lx'] + star.Tracks['Leuv']) * 1e-7 / (4 * np.pi * a**2)
l_bol = star.Tracks['Lbol'] * 1e-7
```

!!! warning "This is a sequence of snapshots, not an evolution"
    The profile does not change along the track. A real planet's structure responds to the loss, to its own cooling, and to the star, so treat what follows as the same atmosphere asked the same question at many stellar ages. The framework is not wired into the coupled loop yet; see [coupling to PROTEUS](../Explanations/proteus.md) for the status.

For a one Earth-mass planet at 0.2 au, carbon dioxide with one percent hydrogen, and an inventory of 100 Earth atmospheres:

```text
 age [Myr]     F_xuv  regime             rate [kg/s]      Kn_sc  snapshot
       1.0      58.4  hydrodynamic:EL      1.273e+07    0.00941  consistent
      45.0      1.62  hydrodynamic:EL      3.542e+05      0.265  consistent
     310.1       1.4  hydrodynamic:EL      3.045e+05      0.388  inconsistent
    1138.9     0.852  hydrostatic          4.695e-02       4.81  consistent
    3039.5     0.495  hydrostatic          4.695e-02   6.73e+03  consistent
    9439.6     0.324  hydrostatic          4.695e-02   4.83e+07  consistent
label changes from hydrodynamic:EL to hydrostatic between 767 and 939 Myr,
and the rate drops from 2.273e+05 to 4.695e-02 kg/s
```

Four things in that table are worth more than the rest of this page.

The label belongs to the state, not to the planet. Nothing about the planet changed; the star quieted down, the Knudsen number climbed smoothly through its threshold, and the physics of the loss changed character.

The rate drops by almost seven orders of magnitude at the crossing. The switch is deliberately sharp, with no blend function between branches, so the size of that jump is a measurement you can quote rather than an artifact a smoothing function hides. It is also the honest scale of the disagreement between the two prescriptions at the same physical state.

What escapes changes with the branch. In the wind phase carbon and oxygen leave with the hydrogen; in the exosphere phase hydrogen leaves alone, and the heavy elements are pinned at 1e-140 kg s⁻¹. A planet crossing this boundary stops losing its atmosphere and starts losing only its hydrogen.

The hydrostatic rate is flat, at 4.695e-02 kg s⁻¹ for ten billion years, because the branch has no XUV physics in it: prescribed exobase temperature, frozen profile, and diffusion-limited supply that does not know about the star. That flat line is a visible reminder of what the branch does not model. It is also why hydrostatic heavy-element rates carry a lower-limit flag: the nonthermal channels that actually remove heavy species from a real exosphere are absent.

And the consistency screen fires in the middle of the track, not at the ends. Between about 40 and 800 Myr the dispatched rate would have emptied the supplied inventory faster than the star aged, so those snapshots are not compatible with their own ages. After the crossing the screen goes quiet again, but read that carefully: the state became consistent because escape effectively stopped, not because the inventory survived.

---

## Things to try

- **Move the planet.** Run the track at 0.1 and at 0.4 au. The crossing age moves; find where the planet never leaves the wind regime within the age of the star.
- **Change the composition at fixed mass and radius.** Add water instead of hydrogen, or drop the hydrogen entirely, and watch both the boundary position and the post-crossing rate respond.
- **Sweep the exobase temperature on a heavy atmosphere.** Confirm for yourself that a pure carbon dioxide planet returns rates far below the floor at every temperature in the range, and that the boundary between "computed" and "meaningful" is where the floor sits.
- **Print the whole container.** `pprint(result.diagnostics)` on one call, and read the groups this page skipped: `bolometric`, `thermostat`, `closure`, and `documentation`, which carries the criterion bands themselves so a stored result is self-describing.
- **Break a state on purpose.** Give the profile an increasing pressure or a negative mass and confirm the `ValueError`, so you know what a malformed state looks like as opposed to an extreme one.

## Where to go next

- [Escape regimes](../Explanations/regimes.md): every branch, threshold, and equation behind this page.
- [Dispatch results](../Reference/results.md): every flag and every diagnostics key.
- [Triage a verdict](../How-to/triage_verdict.md): a flag fired, or the answer looks wrong.
- [Model parameters](../Reference/parameters.md): every setting and input, with defaults.
- [Limitations](../Explanations/limitations.md): what none of this models.

---

[^johnson]: Johnson, R. E., Volkov, A. N., & Erwin, J. T. (2013). Molecular-Kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow. *The Astrophysical Journal Letters, 768*(1), L4. https://doi.org/10.1088/2041-8205/768/1/L4

[^yelle]: Yelle, R. V. (2024). Diffusion limited escape of hydrogen from Mars. *Icarus, 416*, 116099.

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763
