# Limitations

Every prescription in ZEPHYRUS is a deliberate simplification of a richer physical problem. This page collects what each entry point does not model and what that implies for results; the physics each one does model is defined on the [energy-limited escape](energy_limited.md), [escape regimes](regimes.md), [fractionation](fractionation.md), and [giant impacts](impacts.md) pages.

---

## The energy-limited default (`EL_escape`)

The released entry point applies one prescription unconditionally, so every limitation of that prescription passes through to coupled PROTEUS runs until the regime framework is wired in:

- **No regime awareness.** `EL_escape` returns an energy-limited rate whether or not the state sustains a collisional XUV wind. Outside that regime (weak XUV flux, compact atmosphere, boil-off conditions, Roche-filling geometries) the returned rate can be wrong by orders of magnitude in either direction. The [regime framework](regimes.md) classifies the state first; `EL_escape` does not.
- **Radiative cooling is not accounted for.** Line emission, molecular bands, and recombination divert absorbed XUV power away from driving the outflow, which is why the effective efficiency collapses for strongly bound planets (of order $10^{-2}$ above the threshold potential of Caldiroli et al. 2022; see the [energy-limited page](energy_limited.md)). `EL_escape` treats the efficiency $\epsilon$ as a constant input, so $\epsilon = 1$ is a nonphysical upper limit and even canonical values overestimate the loss for compact planets. The regime framework improves on this two ways (a radiatively cooled wind temperature, and the radiation-recombination cap), but its energy-limited efficiency remains an input as well.
- **Bulk removal.** The rate is split over species by reservoir mass fractions, with no preferential loss of light species. The [fractionation closure](fractionation.md) resolves the partition when the regime framework confirms a wind; `EL_escape` alone cannot. For close-in planets where fractionation matters, bulk-removal rates are a lower bound on how fast the atmospheric mean molecular weight grows.
- **$\epsilon$ is constant in time.** The efficiency in reality evolves with mass, radius, and flux; fixed-$\epsilon$ histories can overestimate late-time loss. The fitted-efficiency option of the regime framework captures the potential dependence, not a separate time dependence.

## The regime framework (`dispatch`)

The framework removes the regime-awareness limitation and carries its own, each flagged on the results it affects:

- **The exobase temperature is prescribed.** The hydrostatic branch's rate depends exponentially on $T_\mathrm{exo}$, which the caller supplies (default 1000 K). This is the branch's dominant sensitivity. The optional local-balance estimator is biased high by construction (heating scales with density, the cooling channels with its square, and conduction is absent) and is deliberately not the default.
- **The thermostat evaluates one level.** The wind temperature comes from a local balance at the wind base; the temperature structure through the sonic region is not modeled. The sensitivity propagates: the wind temperature sets the sonic-point density and therefore feeds the collisionality switch itself.
- **Hydrostatic heavy-element rates are lower limits.** The nonthermal channels that dominate heavy-species loss from real exospheres (ion outflow, photochemical ejection, sputtering, charge exchange, ion pickup) are absent. Every hydrostatic result carries a flag saying so, and states where the neutral and plasma escape-temperature conventions disagree are flagged as contested with both branch rates recorded, because the unmodeled ion physics decides them.
- **Thresholds are calibrated elsewhere.** The boil-off activation threshold is calibrated on hydrogen-rich envelopes and transferred to other compositions through the mean molecular mass; the collisionality threshold carries a factor-30 physical band from the heating geometry. Both bands are reported beside every verdict rather than hidden, but a band is not a resolution.
- **Known extrapolations are held, flagged.** The kinetic enhancement on the Jeans flux is measured up to a Jeans parameter of 15 and held constant beyond; the CO$_2$ band's deexcitation rates are measured over roughly 150 to 500 K; the geometric cross-section rung has a documented high-temperature bias. Each engagement is flagged or provenance-classed.
- **The switch inherits the base-pressure choice.** The sonic-point density scales with the wind-base density, so the base-method setting moves where the switch fires; the setting exposes that dependence rather than resolving it.
- **Overflowing states get a bound-flow rate.** The Roche screen names a state whose flow reaches the Hill sphere and reports the rate its branch computed, which is a lower limit: the tidally driven flow through the inner Lagrange point that such a planet actually drives is not modeled, and neither is the accompanying orbital evolution. The subflag and the reported extent of the atmosphere separate the geometry that genuinely overflows from an atmosphere sitting deep inside its lobe with only its sonic surface outside.
- **Impacts are not dispatched.** The giant-impact channel has a reserved label but is invoked directly by the caller, outside the continuous classification.

## The impact channel (`collision.mass_loss`)

The channel is a single fitted power law, not an impact simulation, and inherits the scope of the simulation suite behind it (see the [giant impacts page](impacts.md) for the fitted domain):

- **Thin atmospheres only.** The fit covers atmospheres of order 1% of the planet mass; a substantially thicker envelope cushions the impactor and the eroded fraction is no longer described by the law.
- **Target-side loss only.** Any atmosphere the impactor carries, and any volatile delivery into the merged body, is outside the function. The underlying simulations show a slow grazing collision with an atmosphere-hosting impactor can leave the target with about 85% of the two bodies' combined atmospheres, so treating the impactor as bare is a caller-side assumption.
- **No mantle or core erosion.** Violent impacts also strip silicate and metal mass; the law tracks the atmospheric fraction only.
- **Chaotic-regime scatter.** Slow, head-on collisions produce chaotic fall-back; the fit carries about 20% scatter there against 9% overall.
- **Fit-domain extrapolation is unflagged.** The function evaluates the law for any physically valid inputs and does not warn when they leave the fitted ranges; staying inside them is the caller's responsibility.

## Not modeled by any entry point

- **Nonthermal escape.** Ion pickup, charge exchange, photochemical escape, sputtering, and unmagnetized ion outflow are absent everywhere. For present-day Earth and Venus these dominate the total loss, at rates of order $10^3$ g s$^{-1}$, many orders of magnitude below the rates ZEPHYRUS produces during the early high-XUV phase it targets; they matter for evolved, weakly irradiated states.
- **Magnetic fields.** No entry point knows about planetary or stellar magnetic fields, which can channel, throttle, or enhance the loss.
- **Photochemistry.** Hazes, aerosols, and photochemically produced species are not tracked; the composition the escape sees is the one the atmosphere model supplies.

## Upstream uncertainties

The XUV flux $F_\mathrm{XUV}$ entering every XUV-driven rate carries large intrinsic uncertainty from the stellar model: saturation timescales vary from about 10 to 300 Myr for Sun-like stars and up to a Gyr for fully convective M dwarfs depending on initial rotation, integrated XUV histories differ by factors of a few to ten between standard prescriptions, and the observational anchors are themselves absorbed by the interstellar medium. Integrated mass-loss histories inherit these factors on top of everything above.

## Practical implications

- Avoid $\epsilon > 0.3$ for rocky planets without a specific reason; $\epsilon \approx 0.15$ is the conservative baseline, and for strongly bound planets consider the fitted-efficiency option of the regime framework.
- For close-in planets around M dwarfs, where fractionation is expected, bulk-removal rates bound the growth of the atmospheric mean molecular weight from below; use the regime framework with fractionation on when the composition history matters.
- Treat single-prescription mass-loss histories as scenario calculations rather than predictions: the regime diagnostics reported beside every framework verdict are the tool for judging how sensitive a given history is to the boundary placements.
