# Limitations

Every prescription in ZEPHYRUS is a deliberate simplification of a richer physical problem. This page collects what each entry point does not model and what that implies for results. The physics each one does model is defined on the [escape regimes](regimes.md), [fractionation](fractionation.md), [energy-limited escape](energy_limited.md), and [giant impacts](impacts.md) pages.

## The regime dispatcher (`dispatch`)

The dispatcher decides which escape physics applies before computing a rate. It carries its own limitations, each flagged on the results it affects.

### Temperatures

- **The exobase temperature is prescribed.** The hydrostatic rate depends exponentially on $T_\mathrm{exo}$, which the caller supplies (default 1000 K). This is the branch's dominant sensitivity. The optional local-balance estimator is biased high by construction (heating scales with density, the cooling terms with its square, and conduction is absent), and it is not the default.
- **The thermostat evaluates one level.** The wind temperature comes from a local balance at the wind base, so the temperature structure through the sonic region is not modelled. The wind temperature sets the sonic-point density and therefore feeds the collisionality switch itself.

### Criteria and their thresholds

- **Thresholds are calibrated elsewhere.** The boil-off threshold is calibrated on hydrogen-rich envelopes and transferred to other compositions through the mean molecular mass. The collisionality threshold carries a factor-30 physical band from the heating geometry. Both bands are reported beside every verdict, but a band is not a resolution.
- **The collisionality switch uses neutral cross sections.** The wind can be substantially ionised, and ion-atom and electron-atom collisions have larger cross sections. Chatterjee & Pierrehumbert (2026) [^cp26] choose the neutral onset deliberately and call it reasonable for a weakly ionised, advection-dominated flow and conservative for rapid mass loss. Including the ions would lower the Knudsen number, so every hydrostatic verdict on an ionised wind is one the fuller physics could overturn, and no hydrodynamic verdict is.
- **The switch inherits the wind-base choice.** The sonic-point density scales with the wind-base density, so the `base_method` setting moves where the switch fires. The setting exposes that dependence rather than resolving it.
- **The handover from boil-off to the XUV wind is decided on rate alone.** Owen & Schlichting (2024) [^owensch] give the mass loss to the XUV wind when its photons are absorbed inside the Bondi radius and it drives the stronger outflow. The dispatcher compares the two rates and does not test where the photons are absorbed. On an envelope inflated enough that the XUV is absorbed beyond its Bondi radius, it could dispatch an XUV rate their map assigns to the bolometric wind.
- **The front switches at half atomised hydrogen.** The thermostat and the recombination coefficient take the hydrogen front at an atomised hydrogen fraction of one half or more and the nitrogen-like front below it, so a composition crossing that value steps the wind temperature and the Knudsen number without a physical change.

### Rates

- **Hydrostatic heavy-element rates are lower limits.** The nonthermal channels that dominate heavy-species loss from real exospheres (ion outflow, photochemical ejection, sputtering, charge exchange, ion pickup) are absent. Every hydrostatic result carries a flag saying so. States where the neutral and plasma escape temperatures disagree are flagged as contested, with both branch rates recorded, because the ion physics that would decide them is not modelled.
- **An unstable exosphere has no rate of its own.** Past the exosphere stability test the outflow is driven by the exosphere's own heat, which no branch models. The dispatcher takes the larger of the XUV wind rate and the Jeans rate, and neither describes that outflow. Volkov et al. (2011) [^volkova] find the kinetic rate near the free-evaporation flux of the heated base below a Jeans parameter of about 2, not a Jeans flux, and the wind rate is energy-limited, which is zero without XUV. The rate there is a placeholder that rises with the exobase temperature, not a prediction.
- **Known extrapolations are held and flagged.** The kinetic enhancement of the Jeans flux is simulated for Jeans parameters of 6 to 15 and held constant outside that range. Below 6 the held value likely understates the flux, and trace light species on a heavy background reach that side. The CO$_2$ band's deexcitation rates are measured over roughly 150 to 500 K, and the band's real limit is CO$_2$ dissociation, which makes it a coolant of the region near the base. The geometric cross-section fallback has a documented high-temperature bias. Each of these is flagged or carries a provenance class when it is used.
- **Some cross sections rest on assumed radii.** Six elements (Al, P, Cl, K, Ca, Ti) have no tabulated radius yet and fall back on an assumed one, with a provenance class of its own ([knudsen validation page](../Validation/knudsen.md)).

### The photon-limited rate

The photon limit is an estimate, not a bound.

- It counts one particle per intercepted ionising photon (Owen & Alvarez 2016, Eq. (10) [^owenalvarez]). That limits primary ionisations rather than mass loss: neutrals heated by collisions escape too, and X-ray photoelectrons ionise more than once, so a real flow can exceed it.
- The source never simulates the regime. The one three-dimensional plane-parallel calculation with numbers to compare sits 1.3 to 2.4 times above the count taken on the wind-base disk (Debrecht et al. 2019, MNRAS 483, 1481 [^debrecht]).
- Here the count uses the energy-limited rate's 20 mbar disk, so that the two share one photon budget. The wind base can sit several times further out (1.04 to 9 times in radius on the states tested), so on this disk the limit sits below the source's own Eq. (10) by the square of that ratio, a factor the energy-limited rate shares. Where the photon limit binds, the rate may therefore be low.
- The mass per ionisation is the wind's mean atomic mass per ion, which on a hydrogen front counts helium as ionised although 20 eV photons do not ionise it. The count per photon is then 5% low at a helium mixing ratio of 0.1 and about 9% at solar abundance, a convention the recombination-limited base density shares.
- The rate scales inversely with the photon energy per ionisation, taken as the front energy (20 eV for a hydrogen front). The plausible range of 13.6 to 40 eV spans a factor of 2.9.

`photon_limit = False` recovers the rate without the photon limit for comparison.

### The fractionation closure

The closure solves the system of Zahnle et al. (1990) [^z90] in the subsonic, constant-composition limit and shares their isothermal and neutral-gas assumptions. Three of its simplifications affect results:

- **A hard zero below threshold.** A retained species carries zero flux. The continuous forms of Zahnle & Kasting (1986, their Eqs. (14) and (16)) [^zk86] and of Zahnle et al. (1990, their Eq. (39)) keep a small nonzero escape rate at and just below a species' threshold flux, so the code underestimates the fluxes of species close to their thresholds. Hunten et al. (1987) [^hunten] have the same sharp threshold.
- **One level, no radial integration.** The closure is algebra at the wind base. Zahnle et al. (1990) integrate the minor species' composition with height (their Eqs. (37) to (39)) and solve the limiting-flux problem as a transonic wind. Under constant composition every term of the drift balance, Eq. (1) of the [fractionation page](fractionation.md), scales with the same power of radius for co-escaping species, so one level suffices for them, while retained species enter only through their mole fractions at the base.
- **Atoms, not molecules.** Zahnle et al. (1990) work with molecular H$_2$, CO$_2$, and N$_2$, while the closure runs on the atomised composition at the wind base (see [where it applies](fractionation.md#where-it-applies)).

### Roche overflow

- **The transfer rate is isothermal and has no energy limit.** The rate through L1 (Jackson et al. 2017, their Eq. (3) [^jackson17]) assumes an isothermal flow from the launch level to L1 on a circular, synchronously rotating orbit, and its authors name the temperature as the dominant uncertainty, which the `nozzle_temperature` setting chooses. The flow carries no energy limit, as in the source, and the diagnostics report the heat it demands beside the interior and intercepted stellar luminosities, so a consumer can see where the isothermal assumption is strained.
- **The launch level cancels only along an isothermal column.** $\rho\, e^{\Psi / v_\mathrm{th}^2}$ is constant along an isothermal hydrostatic column, so there the choice of launch level does not change the rate. A coupled run supplies a profile that is not isothermal, and on it the launch level is a width rather than a cancellation.
- **Two radius conventions are omitted.** The launch radius is the profile radius, standing in for the volume-equivalent photospheric radius without the distortion conversion of the source's Appendix, which their own input chain needs because it starts from a measured transit radius. It is worth about 1.6 times in rate per 1% of radius near contact and nothing well inside the lobe. The volume-averaged potentials carry an error of a few per cent near the lobe, which the source quantifies as about a factor of two in the rate.
- **Two things the source computes are omitted.** The torque-balance transfer rate of their Eq. (24) needs the stellar tidal dissipation and can sit orders of magnitude below the nozzle rate where the torques balance. Under that reading the dispatched transfer rate is an upper limit. Whether the transfer is stable couples to orbital evolution, which nothing in this package models.
- **An eccentric orbit is our own extension.** The source has no eccentric treatment. This module evaluates the circular formula at each separation and averages in time over the arc where the transfer applies, so an eccentric result omits the wind the planet drives on the rest of the orbit and inherits a Roche geometry that is synchronous at no single phase.
- **The applicability criteria are this module's.** The source draws the comparison between the sonic radius and L1 qualitatively, and states no fluid condition at L1. The dispatched rate jumps across both criteria, by a measured factor of about $4 \times 10^{3}$ on one family for the first.
- **Free-molecular leakage through L1 is not modelled.** Where the gas at L1 is collisionless, molecules still leave through the saddle, but that is exospheric escape in the tidal potential, which no branch computes. The bound branch's rate stands, and on a planet close to that threshold it misses the tidal part.
- **A renamed state still reports a bound-flow rate.** Where the Roche screen renames a bound branch, the reported rate is that branch's own. It is a lower limit on the tidal transfer only where the transfer through L1 did not apply. Where it applied and lost, the module's own estimate of the transfer is `diagnostics['nozzle']['rate_kg_s']`, below the dispatched rate. The subflag and the reported extent of the atmosphere separate an atmosphere that overflows from one deep inside its lobe with only its sonic surface outside.

### Scope

- **Impacts are not dispatched.** The giant-impact channel has a reserved label but is invoked directly by the caller.

## The energy-limited entry point (`EL_escape`)

`EL_escape` applies one prescription to every state, so every limitation of that prescription passes through to any run that calls it, coupled or standalone:

- **No regime awareness.** `EL_escape` returns an energy-limited rate whether or not the state sustains a collisional XUV wind. Outside that regime (weak XUV flux, compact atmosphere, boil-off conditions, Roche-filling geometries) the rate can be wrong by orders of magnitude in either direction. The [dispatcher](regimes.md) decides the regime first.
- **Radiative cooling is not accounted for.** Line emission, molecular bands, and recombination divert absorbed XUV power away from driving the outflow, which is why the effective efficiency collapses for strongly bound planets (of order $10^{-2}$ above the threshold potential of Caldiroli et al. 2022 [^caldiroli], see the [energy-limited page](energy_limited.md)). `EL_escape` treats the efficiency $\epsilon$ as a constant input, so $\epsilon = 1$ is an unphysical upper limit and even canonical values overestimate the loss for compact planets. The dispatcher improves on this in two ways (a radiatively cooled wind temperature, and the recombination-limited and photon-limited rates), but its energy-limited efficiency is an input as well.
- **Bulk removal.** The rate is split over species by reservoir mass fractions, with no preferential loss of light species. The [fractionation closure](fractionation.md) resolves the partition when the dispatcher confirms a wind, which `EL_escape` alone cannot do. For close-in planets where fractionation matters, bulk-removal rates are a lower bound on how fast the atmospheric mean molecular weight grows.
- **$\epsilon$ is constant in time.** The efficiency in reality evolves with mass, radius, and flux, and fixed-$\epsilon$ histories can overestimate late-time loss. The fitted-efficiency option of the dispatcher captures the dependence on the potential, not a separate time dependence.

## The impact channel (`collision.impact_loss`, `collision.mass_loss`, `collision.mass_loss_roche2026`)

The channel evaluates empirical scaling laws (Kegerreis et al. 2020 for thin atmospheres, and Roche et al. 2026 for envelope mass fractions from 1% to 20%), rather than simulating collisions directly, and inherits the scope of the simulation suites behind them (see the [giant impacts page](impacts.md) for the fitted domains):

- **Atmospheric mass limits.** The `kegerreis2020` fit applies to thin atmospheres of order 1% of the target mass. The `roche2026` fit extends calibration over $f_\mathrm{atm} \in [0.01, 0.20]$, but envelopes outside that range require extrapolation. An airless target ($f_\mathrm{atm} = 0$) returns $X = 0$ immediately, while any non-zero $f_\mathrm{atm} < 10^{-6}$ is evaluated at the $10^{-6}$ stability clamp (only the fit arguments use the clamped value, while $M_\mathrm{t}^\mathrm{tot}$, $v_\mathrm{esc}$ and $Q'_\mathrm{R}$ use the raw $f_\mathrm{atm}$, so $X$ changes only slowly below $10^{-6}$ and drops to 0 at $f_\mathrm{atm} = 0$).
- **Target-side loss only.** Any atmosphere the impactor carries, and any volatile delivery into the merged body, is outside the function. The underlying simulations of Kegerreis et al. (2020) show a slow grazing collision with an atmosphere-hosting impactor can leave the target with about 85% of the two bodies' combined atmospheres, so treating the impactor as bare is a caller-side assumption.
- **Simulation-to-law scatter.** For Kegerreis et al. (2020), the scatter between simulations and the scaling law is a relative deviation: the median relative deviation is 9%, rising to about 20% relative deviation for slow, head-on collisions where outcomes are chaotic. For Roche et al. (2026), misfits are reported as absolute differences on the loss fraction ($|X_\mathrm{sim} - X_\mathrm{calc}|$): over their 790-row fitting set (296 of the 300 new simulations plus 494 runs at $f_\mathrm{atm} = 0.05$ from Roche et al. 2025, Fig. 3), the scaling law has mean and maximum absolute misfits of 0.0396 and 0.2224 respectively.
- **Physical simplifications.** The Roche et al. (2026, Sect. 4.1) scaling law does not account for pre-impact planetary rotation, surface liquid water oceans, thermal evolution, or a core mass fraction other than about 0.3 (Sect. 2.1). Magma-envelope miscibility at high pressures without a sharp boundary can cause the scaling law to overestimate the loss of massive envelopes on young planets (Sect. 4.1). The law covers the immediate shock- and vapour-plume-driven loss only. It neglects later thermally driven loss (an outflow driven by heat from the post-impact interior, Biersteker & Schlichting 2021), so for primordial H2-He envelopes the total loss can be higher (Roche et al. 2026, Sect. 4.1). That later loss becomes negligible for envelopes of higher mean molecular weight. All simulations assume H2-He envelopes. Heavier atmospheres ($\mathrm{CO}, \mathrm{CO}_2$) are less susceptible to shock-driven removal, so shock-driven loss is an upper limit for a given atmosphere mass (Sect. 4.3). For PROTEUS atmospheres of higher mean molecular weight ($\mathrm{H}_2\mathrm{O}, \mathrm{CO}_2, \mathrm{O}_2$), the H2-He loss fractions likewise serve as upper limits.
- **Fit-domain flags.** `collision.impact_loss` records out-of-range parameters in `flags` for the `roche2026` law (`'f_atm'`, `'M_t_earth'`, `'gamma'`, `'b'`, `'R_ratio'`, `'v_ratio'`). It also records diagnostic flags `'v_sub_escape'` for speeds below 0.99 mutual escape speed and `'X_FF_zero_energy'` when the zero-energy far-field loss exceeds zero (which can occur inside the fitted range). The direct functions `collision.mass_loss` and `collision.mass_loss_roche2026` evaluate without warning flags. Staying inside fitted bounds is the caller's responsibility.

## Not modelled by any entry point

- **Nonthermal escape.** Ion pickup, charge exchange, photochemical escape, sputtering, and unmagnetised ion outflow are absent everywhere. For present-day Earth and Venus these dominate the total loss, at rates of order $10^3$ g s$^{-1}$, many orders of magnitude below the rates ZEPHYRUS produces during the early high-XUV phase it targets. They matter for evolved, weakly irradiated states.
- **Magnetic fields.** No entry point knows about planetary or stellar magnetic fields, which can channel, throttle, or enhance the loss.
- **Photochemistry.** Hazes, aerosols, and photochemically produced species are not tracked. The composition the escape sees is the one the atmosphere model supplies.

## Upstream uncertainties

The XUV flux $F_\mathrm{XUV}$ entering every XUV-driven rate carries large intrinsic uncertainty from the stellar model: saturation timescales vary from about 10 to 300 Myr for Sun-like stars and up to a Gyr for fully convective M dwarfs depending on initial rotation, integrated XUV histories differ by factors of a few to ten between standard prescriptions, and the observational anchors are themselves absorbed by the interstellar medium. Integrated mass-loss histories inherit these factors on top of everything above.

## Practical implications

- Avoid $\epsilon > 0.3$ for rocky planets without a specific reason. $\epsilon \approx 0.15$ is the conservative baseline, and for strongly bound planets consider the fitted-efficiency option of the regime framework.
- For close-in planets around M dwarfs, where fractionation is expected, bulk-removal rates bound the growth of the atmospheric mean molecular weight from below. Use the dispatcher with fractionation on when the composition history matters.
- Treat single-prescription mass-loss histories as scenario calculations rather than predictions: the regime diagnostics reported beside every dispatcher verdict are the tool for judging how sensitive a given history is to the boundary placements.

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763

[^owenalvarez]: Owen, J. E., & Alvarez, M. A. (2016). UV Driven Evaporation of Close-in Planets: Energy-limited, Recombination-limited, and Photon-limited Flows. *The Astrophysical Journal, 816*(1), 34. https://doi.org/10.3847/0004-637X/816/1/34

[^debrecht]: Debrecht, A., Carroll-Nellenback, J., Frank, A., McCann, J., Murray-Clay, R., & Blackman, E. G. (2019). Photoevaporative flows from exoplanet atmospheres: a 3D radiative hydrodynamic parameter study. *Monthly Notices of the Royal Astronomical Society, 483*(2), 1481–1495. https://doi.org/10.1093/mnras/sty3212

[^jackson17]: Jackson, B., Arras, P., Penev, K., Peacock, S., & Marchant, P. (2017). A new model of Roche lobe overflow for short-period gaseous planets and binary stars. *The Astrophysical Journal, 835*(2), 145. https://doi.org/10.3847/1538-4357/835/2/145

[^owensch]: Owen, J. E., & Schlichting, H. E. (2024). Mapping out the parameter space for photoevaporation and core-powered mass-loss. *Monthly Notices of the Royal Astronomical Society, 528*(2), 1615–1629.

[^volkova]: Volkov, A. N., et al. (2011). Thermally driven atmospheric escape: transition from hydrodynamic to Jeans escape. *The Astrophysical Journal Letters, 729*(2), L24.

[^hunten]: Hunten, D. M., Pepin, R. O., & Walker, J. C. G. (1987). Mass Fractionation in Hydrodynamic Escape. *Icarus, 69*, 532–549.

[^z90]: Zahnle, K., Kasting, J. F., & Pollack, J. B. (1990). Mass Fractionation of Noble Gases in Diffusion-Limited Hydrodynamic Hydrogen Escape. *Icarus, 84*(2), 502–527.

[^zk86]: Zahnle, K. J., & Kasting, J. F. (1986). Mass Fractionation during Transonic Escape and Implications for Loss of Water from Mars and Venus. *Icarus, 68*(3), 462–480.

[^cp26]: Chatterjee, R. D., & Pierrehumbert, R. T. (2026). Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline. *The Astrophysical Journal, 998*(2), 236. https://doi.org/10.3847/1538-4357/ae2ffa
