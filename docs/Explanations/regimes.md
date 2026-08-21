# Escape regimes

This page defines the escape-regime framework behind `zephyrus.dispatch`: the quantities it evaluates, the order it evaluates them in, and the rate physics of every branch. The [model overview](model.md) gives the short version with a flowchart; here every threshold and equation is spelled out. The [energy-limited escape](energy_limited.md) and [fractionation](fractionation.md) pages define the two pieces that have their own pages.

One call takes one planetary state and returns one verdict. The inputs are the planet mass $M_\mathrm{p}$ and interior radius $R_\mathrm{p}$, the stellar mass $M_\star$, the orbit ($a$, $e$), the equilibrium temperature $T_\mathrm{eq}$, the XUV and interior heat fluxes $F_\mathrm{XUV}$ and $F_\mathrm{int}$, the photospheric opacity $\kappa$, and an atmosphere profile (pressure, radius, temperature, composition, and mean molecular mass per level, from the base to the top of the modeled atmosphere). The output is one of the five regime labels, a bulk mass-loss rate, per-species rates that sum to it, flags recording every clamp and fallback, and a diagnostics container that reports how close the state sat to each boundary. Every physically posed input returns a result; exceptions are reserved for malformed input.

## The evaluation order

1. The bolometric (boil-off) candidate is computed at every call. If the restricted Jeans parameter sits below its threshold, the atmosphere is boiling off and that candidate is the rate; XUV-driven escape needs a stable base to launch from, and a bolometrically boiling atmosphere has not built one yet, which is why this test precedes everything else [^owensch].
2. Otherwise the hydrodynamic candidate is assembled: the wind base is located on the profile, a thermostat sets the wind temperature, and the candidate is the smaller of the energy-limited and radiation-recombination-limited rates, with the winner naming the sub-label.
3. The sonic-point Knudsen number decides whether that wind is collisional enough to exist. If it is, the hydrodynamic label stands and the [fractionation closure](fractionation.md) partitions the rate over species; if not, the state re-routes to the hydrostatic branch.
4. The hydrostatic branch evaluates per-species Jeans escape on an extended upper structure. Its stability gate can send a thermally unstable exosphere back to the hydrodynamic rate.
5. Before the label is finalized, the Roche screen tests the active flow radius against the Hill sphere; an overflowing state is relabeled `roche_overflow`.
6. The final rate is the larger of the surviving branch rate and the luminosity-capped bolometric residual, labeled by the winner.

## Boil-off

A freshly formed or strongly heated planet can hold an atmosphere so distended that its outer layers sit beyond the sonic point of a thermal wind: the gas then flows out on the planet's own thermal energy alone. The activation criterion is the restricted Jeans parameter [^fossati]

$$\Lambda \;=\; \frac{G\,M_\mathrm{p}\,\mu}{k_\mathrm{B}\,T_\mathrm{eq}\,R_\mathrm{p}} \tag{1}$$

the ratio of a particle's gravitational binding energy at the surface to its thermal energy, built with the mean molecular mass $\mu$ of the atmosphere at the photospheric level and the Boltzmann constant $k_\mathrm{B}$. For isothermal gas $\Lambda = 2 R_\mathrm{B} / R_\mathrm{p}$ identically, where $R_\mathrm{B} = G M_\mathrm{p} / (2 c_\mathrm{s}^2)$ is the sonic (Bondi) radius at isothermal sound speed $c_\mathrm{s}$, so the shutoff Owen & Wu (2016) find at $R_\mathrm{p}/R_\mathrm{B} = 0.1$ is $\Lambda = 20$ for every composition [^owenwu]. That threshold is the default, with the literature spread of 15 to 35 reported as its band; its calibration on hydrogen-rich envelopes is an assumption the diagnostics keep visible.

While $\Lambda < 20$ the state is labeled `boiloff` and the rate is the closed-form transonic Parker wind of Owen & Wu (2016), evaluated at wind temperature $T_\mathrm{eq}/2^{1/4}$ (the recommendation of Misener et al. 2025 for the isothermal formulas [^misener]):

$$\dot{M}_\mathrm{Parker} \;=\; \frac{4\pi\,G\,M_\mathrm{p}\,\mathcal{M}}{\kappa\,c_\mathrm{s}}, \qquad \mathcal{M} = \sqrt{-W_0\!\left(-x^{-4}\,e^{\,3 - 4/x}\right)}, \quad x = \frac{R_\mathrm{launch}}{R_\mathrm{B}} \tag{2}$$

where $\mathcal{M}$ is the Mach number at the launch level (the photospheric level, radius $R_\mathrm{launch}$), $W_0$ is the principal branch of the Lambert function, and $\kappa$ is the photospheric opacity, which the rate scales inversely with. At $x = 1$ the launch level is sonic and $\mathcal{M} = 1$; for small $x$ the rate shuts off exponentially, which is the physical end of boil-off. The rate is capped by the Bondi-limited supply of Gupta & Schlichting (2020) [^gs20],

$$\dot{M}_\mathrm{B} \;=\; 4\pi R_\mathrm{B}^2\, c_\mathrm{s}\, \rho_\mathrm{launch}\, \exp\!\left(-\frac{G M_\mathrm{p}}{c_\mathrm{s}^2 R_\mathrm{launch}}\right) \tag{3}$$

with $\rho_\mathrm{launch}$ the mass density at the launch level. Past the $\Lambda$ gate the same machinery survives as a residual, additionally capped by the interior luminosity, $\dot{M}_\mathrm{E} = L / (g R_\mathrm{p})$ with $L = 4\pi R_\mathrm{p}^2 F_\mathrm{int}$ and $g$ the surface gravity [^gs19]. Keeping the residual luminosity-capped represents the slow late tail of core-powered mass loss without adjudicating the open dispute over how long it survives (Tang et al. 2024 argue it is brief [^tang]; Gupta & Schlichting argue it lasts); the Tang et al. timescale comparison runs as a diagnostic beside the rate, never as a gate.

## The hydrodynamic wind

Past the boil-off gate, stellar XUV heating can drive a fluid wind. Three pieces are assembled.

The wind base. XUV photons deposit their energy where the atmosphere first reaches unit optical depth to them, at the pressure level $P_\mathrm{base} = \mu g / \sigma_{\nu_0}$, with $\sigma_{\nu_0}$ the photoionization cross section at the representative photon energy (the level Lopez 2017 builds the wind on, about a nanobar [^lopez2017]; the cross section follows Murray-Clay et al. 2009 [^mc09]). The base is located on the supplied profile by fixed-point iteration; when the profile is too shallow to reach it, the level clamps to the profile top with the clamp distance recorded, or is evaluated on the extended upper structure (the `extend` option), and either way the choice is flagged.

The wind temperature. Rather than assuming the canonical $10^4$ K, a thermostat balances local photoionization heating against radiative cooling at the base. Heating follows the monochromatic-front approximation, $Q_\mathrm{heat} = n_0\, \sigma\, (F_\mathrm{XUV}/h\nu)\,(h\nu - E_\mathrm{ion})$, with $n_0$ the neutral density, $h\nu$ the representative photon energy, and $E_\mathrm{ion}$ the ionization potential. Four cooling channels oppose it: atomic line cooling by H, C, C$^+$, N, N$^+$, O, and O$^+$ in three-level statistical equilibrium (the machinery of Chatterjee & Pierrehumbert 2026 [^cp26] on the atomic data of Nakayama et al. 2022 [^nakayama]), the CO$_2$ 15 micron band and the atomic oxygen fine structure (Johnstone et al. 2018 [^johnstone]), and recombination cooling (with the radiative recombination fit of Badnell 2006 [^badnell]). The balance is bracketed between $T_\mathrm{eq}$ and $5 \times 10^4$ K; a balance without a root clamps to the nearer edge, flagged. A high clamp is the expected outcome at dense bases, where electron densities far above the forbidden-line critical densities quench the line coolants collisionally and the wind runs hot.

The two rate limits. The energy-limited rate is Eq. (1) of the [energy-limited page](energy_limited.md) in the Erkaev form (`scaling=2`, $\xi = R_\mathrm{Hill}/R_\mathrm{p}$), with the efficiency either fixed or taken from the fitted efficiency of Caldiroli et al. (2022) converted to that geometry [^caldiroli]. The radiation-recombination-limited (RR) rate follows the analytic chain of Murray-Clay et al. (2009) [^mc09]: at high flux the base ionization reaches equilibrium between photoionization and recombination, which fixes the base ion density to

$$n_+ \;=\; \sqrt{\frac{F_\mathrm{XUV}\, G M_\mathrm{p}}{h\nu_0\, \alpha_\mathrm{B}\, c_\mathrm{s}^2\, R_\mathrm{base}^2}} \tag{4}$$

with $\alpha_\mathrm{B}$ the case B recombination coefficient of the composition (carrying its $T^{-0.9}$ temperature dependence) and $R_\mathrm{base}$ the base radius. An isothermal wind then carries that density to the sonic radius $R_\mathrm{s} = G M_\mathrm{p} / (2 c_\mathrm{s}^2)$ with the barometric factor $e^{\,3/2 - \lambda_\mathrm{b}}$, where $\lambda_\mathrm{b} = G M_\mathrm{p} / (R_\mathrm{base} c_\mathrm{s}^2)$ is the Jeans parameter at the base, giving

$$\dot{M}_\mathrm{RR} \;=\; 4\pi\, \rho_\mathrm{s}\, c_\mathrm{s}\, R_\mathrm{s}^2, \qquad \rho_\mathrm{s} = \rho_\mathrm{base}\, e^{\,3/2 - \lambda_\mathrm{b}} \tag{5}$$

The energy in Eq. (4) has been spent ionizing and is re-radiated on recombination, which is why this limit can undercut the energy-limited one. The candidate is $\min(\dot{M}_\mathrm{EL}, \dot{M}_\mathrm{RR})$ and the winner names the sub-label, `hydrodynamic:EL` or `hydrodynamic:RR`. One caution travels with the label: the minimum selects RR through two physically different mechanisms, genuine recombination saturation at modest $\lambda_\mathrm{b}$ (the $\sqrt{F_\mathrm{XUV}}$ regime of Eq. 4) and plain barometric suppression at large $\lambda_\mathrm{b}$, where calling the result recombination-limited would be a category error; the diagnostics report which mechanism acted. When the computed sonic radius falls below the base, the sonic radius is floored at the base with the base density (no barometric factor below the base) and the state is flagged subcritical. The crossover flux between the two sub-labels is sensitive to the wind temperature the thermostat returns, which enters the RR chain through the sound speed, the barometric exponent, and the recombination coefficient.

## The collisionality switch

A fluid wind only exists if the gas is still collisional where it goes sonic. The switch compares the mean free path against the density scale height at the sonic point, following the construction of Chatterjee & Pierrehumbert (2026), their Eqs. 17 and 18 [^cp26]:

$$\mathrm{Kn}_\mathrm{sc} \;=\; \frac{\ell}{H_\mathrm{sc}}, \qquad \ell = \frac{1}{\sqrt{2}\,\sigma_\mathrm{C}\, n_\mathrm{sc}}, \qquad H_\mathrm{sc} = \frac{(1+\gamma)\, R_\mathrm{s}}{4 + \sqrt{2}\sqrt{5 - 3\gamma}} \tag{6}$$

where $\ell$ is the Maxwell mean free path, $n_\mathrm{sc}$ the particle density at the sonic point (taken from the isothermal wind of Eq. 5), $\gamma$ the polytropic index (1 for an isothermal wind), and $\sigma_\mathrm{C}$ the density-weighted collision cross section of the mixture. Cross sections come from a provenance-classed ladder: tabulated collision integrals where they exist (Laricchiuta et al. 2009 [^laricchiuta], validated against measured viscosities), a diffusion-coefficient inversion for hydrogen (Zahnle et al. 1990 [^z90] on the compilation of Zahnle & Kasting 1986 [^zk86]), and a geometric hard sphere as last resort, whose bias is documented and flagged.

A state with $\mathrm{Kn}_\mathrm{sc}$ at or below the threshold sustains the wind and keeps the hydrodynamic label; above it, the gas decouples before reaching sonic conditions and the state re-routes to the hydrostatic branch. The default threshold is 1, and its physical band is 0.1 to 3: kinetic simulations place the transition near 0.1 when the heating is deposited in a sharp layer and near 1 when it is distributed (Johnson et al. 2013 [^johnson]), so the band is heating-geometry physics rather than tuning freedom, and the diagnostics report the counterfactual labels at both band edges beside every verdict. For evolutionary use, a supplied previous regime label activates a hysteresis window (factor 1.5) around the threshold so a time-stepping track cannot chatter between branches on numerical noise.

## Hydrostatic escape

Where no wind exists, escape proceeds particle by particle from the exobase, the level where the mean free path first reaches the local scale height. All exobase quantities are evaluated on an extended upper structure: a Bates temperature profile $T(\zeta) = T_\mathrm{exo} - (T_\mathrm{exo} - T_\mathrm{top})\, e^{-\gamma_\mathrm{B} \zeta}$ integrated hydrostatically above the supplied profile top (in $\zeta = \ln(p_\mathrm{top}/p)$, shape parameter $\gamma_\mathrm{B}$; the form Yelle 2024 uses [^yelle]), with the exobase temperature $T_\mathrm{exo}$ prescribed by the caller. Extending the structure is not a refinement but a requirement: the Jeans parameter at the true exobase can differ from its photospheric value by an order of magnitude, and evaluating the escape on photospheric values biases rates toward false retention by up to three orders of magnitude (Johnson et al. 2013 [^johnson]). The prescribed $T_\mathrm{exo}$ (default 1000 K) is the branch's dominant sensitivity, because the rate depends on it exponentially; an optional estimator balances local heating against cooling at the profile top, but a conduction-free local balance is biased high by construction and is deliberately not the default.

Each species $i$ escapes with the Jeans effusion flux at the exobase (radius $r_\mathrm{exo}$, temperature $T_\mathrm{exo}$),

$$w_\mathrm{J} \;=\; \sqrt{\frac{k_\mathrm{B} T_\mathrm{exo}}{2\pi m_i}}\,(1 + \lambda_i)\, e^{-\lambda_i}, \qquad \lambda_i = \frac{G M_\mathrm{p} m_i}{k_\mathrm{B} T_\mathrm{exo}\, r_\mathrm{exo}} \tag{7}$$

where $m_i$ is the particle mass and $\lambda_i$ the species Jeans parameter, multiplied by the kinetic enhancement factor $C(\lambda)$ that direct simulation Monte Carlo runs find above the equilibrium Jeans flux: about 1.7 at $\lambda = 6$, falling to about 1.4 at $\lambda = 15$ (Volkov et al. 2011 [^volkova][^volkovb]), held constant beyond 15 as a flagged extrapolation. The escape of a minor species is additionally capped by how fast diffusion can resupply it through the background gas: the diffusion-limited flux $\Phi_\mathrm{l}$ follows the formulation of Yelle (2024) [^yelle] on binary diffusion coefficients that each carry a provenance class, and the two limits combine as the harmonic mean, $\Phi = \Phi_\mathrm{J}\,\Phi_\mathrm{l} / (\Phi_\mathrm{J} + \Phi_\mathrm{l})$, their Eq. 14. The dominant species supplies itself and takes the Jeans flux alone.

Two escape temperatures gate the branch's validity. The neutral escape temperature $T_\mathrm{esc} = G M_\mathrm{p} m / (2 k_\mathrm{B} r_\mathrm{exo})$ marks where thermal energy rivals binding energy; the plasma escape temperature is half of it, because in an ionized exosphere the ambipolar electric field shares each ion's binding energy with its electron (Chatterjee & Pierrehumbert 2026, their Eq. 34 [^cp26]). An exobase hotter than half the gating escape temperature cannot remain hydrostatic, and such states re-route to the hydrodynamic rate. Both temperatures are always computed; states where the two conventions disagree are flagged as contested, with both branch rates recorded, because the ion physics that would decide them is not modeled in this version. For the same reason, hydrostatic heavy-element rates are lower limits (the non-thermal channels that dominate heavy-species loss from real exospheres are absent), and every hydrostatic result carries a flag saying so.

## The Roche screen and overflow

Everything above assumes the flow is bound to the planet. Before the label is finalized, the active flow radius of the winning branch (the Bondi radius on the bolometric branch, the larger of the XUV and sonic radii on the hydrodynamic branch, the exobase radius on the hydrostatic branch) is tested against the periapsis Hill radius of Eq. (3) on the [energy-limited page](energy_limited.md). A flow that reaches the Hill sphere is not described by any of the four regimes above: the state is labeled `roche_overflow` and carries the Bondi-capped bolometric rate at the overflow geometry, with a subflag separating geometries whose Hill sphere sits inside the photosphere itself from those where only the flow reaches it. Near misses (flow radius within 1.5 Hill radii) raise a `near_roche` flag, because the tidal factor inflates the energy-limited rate steeply there; the flag reports, and never modifies, the rate.

## Boundaries are bands

Every threshold above carries a stated physical width, and the framework reports the width instead of hiding it behind a sharp switch. Beside every verdict, the diagnostics container carries: the counterfactual labels at the Knudsen band edges 0.1 and 3; the boil-off activation band 15 to 35; the transonic energy criterion of Johnson et al. (2013) [^johnson] (can the absorbed power drive the flow sonic at all); the Jeans-parameter triple of Guo (2024) [^guo], which translates the verdict into that taxonomy; both escape temperatures with the local ionization fraction; the tidally corrected critical exobase temperature of Erkaev et al. (2007) [^erkaev]; the fluid condition checked level by level below the sonic radius, after Owen & Jackson (2012) [^oj12]; the threshold-potential screens (the efficiency-collapse band of Caldiroli et al. 2022 [^caldiroli], and a wind-versus-hydrostatic screen commonly attributed to Salz et al. 2016, quoted from secondary literature and marked as such); the boil-off termination timescales of Tang et al. (2024) [^tang]; a snapshot self-consistency screen (would the dispatched rate have destroyed the atmosphere within the system age); and the coefficient provenance class of every species. The container is reporting only: nothing in the dispatch control flow reads it, and it has no off switch.

## Configuration

All knobs, their defaults, and their meanings are tabulated in the [parameter reference](../Reference/parameters.md); the defaults are the documented reference choices used throughout this page. Every field of the result, every flag, and every diagnostics group is tabulated in the [dispatch results reference](../Reference/results.md). The assumptions that remain on every result, whatever the knobs, are collected on the [limitations page](limitations.md).

For the framework in use rather than in principle, the [dispatcher tutorial](../Tutorials/dispatch.md) crosses two of the boundaries above on one planet, measures how far one of them moves across the width of its own criterion, and dispatches an atmosphere along a stellar history; the [troubleshooting guide](../How-to/troubleshooting.md) starts from a flag or an unexpected verdict instead.

---

[^owenwu]: Owen, J. E., & Wu, Y. (2016). Atmospheres of low-mass planets: the "boil-off". *The Astrophysical Journal, 817*(2), 107.

[^owensch]: Owen, J. E., & Schlichting, H. E. (2024). Mapping out the parameter space for photoevaporation and core-powered mass-loss. *Monthly Notices of the Royal Astronomical Society, 528*(2), 1615–1629.

[^fossati]: Fossati, L., et al. (2017). Aeronomical constraints to the minimum mass and maximum radius of hot low-mass planets. *Astronomy & Astrophysics, 598*, A90.

[^misener]: Misener, W., et al. (2025). Blowin' in the Nonisothermal Wind: Core-powered Mass Loss with Hydrodynamic Radiative Transfer. *The Astrophysical Journal, 980*(1), 152.

[^gs19]: Gupta, A., & Schlichting, H. E. (2019). Sculpting the valley in the radius distribution of small exoplanets as a by-product of planet formation: the core-powered mass-loss mechanism. *Monthly Notices of the Royal Astronomical Society, 487*(1), 24–33.

[^gs20]: Gupta, A., & Schlichting, H. E. (2020). Signatures of the core-powered mass-loss mechanism in the exoplanet population: dependence on stellar properties and observational predictions. *Monthly Notices of the Royal Astronomical Society, 493*(1), 792–806.

[^tang]: Tang, Y., et al. (2024). Assessing Core-powered Mass Loss in the Context of Early Boil-off: Minimal Long-lived Mass Loss for the Sub-Neptune Population. *The Astrophysical Journal, 976*(2), 221.

[^mc09]: Murray-Clay, R. A., Chiang, E. I., & Murray, N. (2009). Atmospheric Escape From Hot Jupiters. *The Astrophysical Journal, 693*(1), 23–42. https://doi.org/10.1088/0004-637X/693/1/23

[^lopez2017]: Lopez, E. D. (2017). Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor. *Monthly Notices of the Royal Astronomical Society, 472*(1), 245–253.

[^erkaev]: Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007). Roche lobe effects on the atmospheric loss from "Hot Jupiters". *Astronomy & Astrophysics, 472*(1), 329–334. https://doi.org/10.1051/0004-6361:20066929

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763

[^cp26]: Chatterjee, R., & Pierrehumbert, R. T. (2026). Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline. arXiv:2412.05188.

[^nakayama]: Nakayama, A., Ikoma, M., & Terada, N. (2022). Survival of Terrestrial N$_2$-O$_2$ Atmospheres in Violent XUV Environments through Efficient Atomic Line Radiative Cooling. *The Astrophysical Journal, 937*(2), 72. https://doi.org/10.3847/1538-4357/ac86ca

[^johnstone]: Johnstone, C. P., Güdel, M., Lammer, H., & Kislyakova, K. G. (2018). The Upper Atmospheres of Terrestrial Planets: Carbon Dioxide Cooling and the Earth's Thermospheric Evolution. *Astronomy & Astrophysics, 617*, A107. https://doi.org/10.1051/0004-6361/201832776

[^badnell]: Badnell, N. R. (2006). Radiative recombination data for modelling dynamic finite-density plasmas. *The Astrophysical Journal Supplement Series, 167*, 334. arXiv:astro-ph/0604144.

[^laricchiuta]: Laricchiuta, A., Bruno, D., Capitelli, M., et al. (2009). High temperature Mars atmosphere. Part I: transport cross sections. *The European Physical Journal D, 54*(3), 607–612. https://doi.org/10.1140/epjd/e2009-00192-7

[^z90]: Zahnle, K., Kasting, J. F., & Pollack, J. B. (1990). Mass Fractionation of Noble Gases in Diffusion-Limited Hydrodynamic Hydrogen Escape. *Icarus, 84*(2), 502–527.

[^zk86]: Zahnle, K. J., & Kasting, J. F. (1986). Mass Fractionation during Transonic Escape and Implications for Loss of Water from Mars and Venus. *Icarus, 68*(3), 462–480.

[^johnson]: Johnson, R. E., Volkov, A. N., & Erwin, J. T. (2013). Molecular-Kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow. *The Astrophysical Journal Letters, 768*(1), L4. https://doi.org/10.1088/2041-8205/768/1/L4

[^volkova]: Volkov, A. N., et al. (2011). Thermally driven atmospheric escape: transition from hydrodynamic to Jeans escape. *The Astrophysical Journal Letters, 729*(2), L24.

[^volkovb]: Volkov, A. N., Tucker, O. J., Erwin, J. T., & Johnson, R. E. (2011). Kinetic simulations of thermal escape from a single component atmosphere. *Physics of Fluids, 23*(6), 066601. https://doi.org/10.1063/1.3592253

[^yelle]: Yelle, R. V. (2024). Diffusion limited escape of hydrogen from Mars. *Icarus, 416*, 116099.

[^guo]: Guo (2024). Characterizing regimes of hydrodynamic escape of close-in low mass exoplanets. arXiv:2405.13283.

[^oj12]: Owen, J. E., & Jackson, A. P. (2012). Planetary evaporation by UV and X-ray radiation: basic hydrodynamics. *Monthly Notices of the Royal Astronomical Society, 425*(4), 2931. https://doi.org/10.1111/j.1365-2966.2012.21481.x
