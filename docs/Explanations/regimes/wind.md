# The hydrodynamic wind

Stellar XUV photons heat the upper atmosphere, and where the gas is collisional enough the heating drives a fluid wind. The dispatcher computes this wind on both sides of the [boil-off criterion](boiloff.md): while the criterion holds the wind competes with the bolometric outflow, and otherwise the [collisionality switch](collisionality.md) decides whether it exists at all. Its rate is the smallest of three limits built on one wind base and one wind temperature. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The wind base

XUV photons deposit their energy where the atmosphere first reaches unit optical depth to them. Over one scale height that is the pressure level of Lopez (2017) [^lopez2017],

$$P_\mathrm{base} \;=\; \frac{\mu_\mathrm{base}\, g_\mathrm{base}}{\sigma_{\nu_0,\mathrm{H}}}, \qquad g_\mathrm{base} = \frac{G M_\mathrm{p}}{R_\mathrm{base}^2} \tag{1}$$

about a nanobar, with $\mu_\mathrm{base}$ the profile's mean molecular mass at that level. The cross section is the hydrogen photoionisation cross section of Murray-Clay et al. (2009) at their representative photon energy $h\nu_0 = 20$ eV [^mc09], $\sigma_{\nu_0,\mathrm{H}} = 6 \times 10^{-18} (h\nu_0 / 13.6\ \mathrm{eV})^{-3}$ cm², used for every composition. Because $\mu_\mathrm{base}$ and $R_\mathrm{base}$ depend on the level the pressure selects, the base is found by fixed-point iteration on the profile.

Two settings change the base:

- `base_out_of_range`: when the profile top is deeper than $P_\mathrm{base}$, the default clamps the base to the profile top and records the clamp distance in decades of pressure (`base_clamped`). The value `'extend'` evaluates the base on the extended upper structure of the [hydrostatic branch](hydrostatic.md) instead (`base_extended`).
- `base_method`: a fixed pressure, or the XUV radius of the BOREAS solver, in place of Eq. (1).

## The wind temperature

The wind temperature $T_\mathrm{w}$ is not assumed to be the canonical $10^4$ K. A thermostat balances local photoionisation heating against radiative cooling at the base density $n_\mathrm{base}$ and the atomised base composition (element mole fractions $y_j$). $T_\mathrm{w}$ is the root of

$$(1 - f_+)\, n_\mathrm{base}\, \sigma_{\nu_0}\, \frac{F_\mathrm{XUV}}{h\nu_0}\, \left(h\nu_0 - E_\mathrm{ion}\right) \;=\; Q_\mathrm{lines} + Q_{\mathrm{CO_2}} + Q_\mathrm{O} + \tfrac{3}{2}\, k_\mathrm{B} T_\mathrm{w}\, \alpha_\mathrm{rec}(T_\mathrm{w})\, f_+^2\, n_\mathrm{base}^2 \tag{2}$$

**Heating.** The left side is the photoionisation heating rate per unit volume in the monochromatic-front approximation. Each ionisation leaves $h\nu_0 - E_\mathrm{ion}$ in the gas, with $h\nu_0$ the front's photon energy, $E_\mathrm{ion}$ its ionisation potential, and $\sigma_{\nu_0}$ its photoionisation cross section. $f_+$ is the ionisation fraction from local photoionisation-recombination balance, $f_+^2 / (1 - f_+) = \sigma_{\nu_0} F_\mathrm{XUV} / (h\nu_0\, \alpha_\mathrm{rec}\, n_\mathrm{base})$. The front follows the composition:

- the hydrogen front ($h\nu_0 = 20$ eV, $E_\mathrm{ion} = 13.6$ eV, $\sigma_{\nu_0} = \sigma_{\nu_0,\mathrm{H}}$) when the atomised hydrogen fraction is one half or more,
- the nitrogen-like front of Chatterjee & Pierrehumbert (2026) ($h\nu_0 = 33.6$ eV, $E_\mathrm{ion} = 14.53$ eV, $\sigma_{\nu_0} = 10^{-17}$ cm²) otherwise [^cp26].

**Cooling.** Four terms stand on the right, each with its own setting, and at least one must stay on:

- $Q_\mathrm{lines}$ (`cool_atomic`): atomic line cooling by H, C, C$^+$, N, N$^+$, O, and O$^+$, in three-level statistical equilibrium under electron impact, with every emitted photon escaping. This is the treatment of Chatterjee & Pierrehumbert (2026), their Eqs. (26) to (30), on the atomic data of Nakayama et al. (2022) [^nakayama]. The hydrogen system carries Lyman alpha.
- $Q_{\mathrm{CO_2}}$ (`cool_co2_band`): the CO$_2$ 15 micron band (Johnstone et al. 2018 [^johnstone]). It is a coolant of the region near the base, evaluated on molecular abundances.
- $Q_\mathrm{O}$ (`cool_o_finestructure`): the atomic oxygen fine-structure lines at 63 and 147 micron, from the same source.
- The last term (`cool_recombination`): the continuum part of recombination cooling, with $\alpha_\mathrm{rec}$ the radiative recombination fit of Badnell (2006) for nitrogen-dominated gas [^badnell] and the composition-weighted case B coefficient of Eq. (5) otherwise.

**The root.** The balance is scanned upward from $T_\mathrm{eq}$ to $5 \times 10^4$ K, and the first downward crossing of heating through cooling is taken, which selects the lowest stable root. A balance with no root in that range is clamped to the nearer edge and flagged `thermostat_clamped`. A high clamp is expected at dense bases, where electron densities far above the critical densities of the forbidden lines quench the line coolants and the wind runs hot. The upper edge is a validity limit: raising it finds a root near $10^5$ K, but the coolants are neutral three-level systems and gas at that temperature is fully ionised, so that root is outside the model. A clamped temperature is therefore the edge of the range, not a solution, and every rate below that is built on $T_\mathrm{w}$ inherits that.

## The energy-limited rate

If a fraction $\epsilon$ of the intercepted XUV power lifts gas out of the tidally reduced well, the rate is that power divided by the escape energy per unit mass. The dispatcher uses the form of Erkaev et al. (2007), their Eq. (21) [^erkaev], which is `scaling=2` of the [energy-limited entry point](../energy_limited.md):

$$\dot{M}_\mathrm{EL} \;=\; \frac{\epsilon\,\pi\,F_\mathrm{XUV}\,R_\mathrm{p}\,R_\mathrm{XUV}^2}{G\,M_\mathrm{p}\,K(\xi)}, \qquad R_\mathrm{XUV} = R_\mathrm{ph} \tag{3}$$

with $K$ and $\xi = R_\mathrm{Hill}/R_\mathrm{p}$ from [Eq. (2) of the overview](../regimes.md#tidal-geometry). The absorbing disk $\pi R_\mathrm{XUV}^2$ is taken at the photospheric working level and the barrier $G M_\mathrm{p} K / R_\mathrm{p}$ at the interior radius. The factor $\pi$ spreads the intercepted power over the full surface and is not a dayside cross section.

The efficiency $\epsilon$ is the setting `efficiency` (default 0.1). Under `efficiency_mode = 'caldiroli'` it is instead the fitted evaporation efficiency $\eta$ of Caldiroli et al. (2022, their Appendix A.1) [^caldiroli]. Their $\eta$ is defined against an $R_\mathrm{p}^3$ geometry, so it is converted to $\epsilon = \eta\, (R_\mathrm{p}/R_\mathrm{XUV})^2$ before use. Below the fit's flux bound the fixed efficiency is used instead, and a flag says so.

## The recombination-limited rate

At high flux the energy absorbed at the base is spent ionising and is re-radiated on recombination, so the base ionisation, not the energy budget, sets the rate. The calculation is the analytic one of Murray-Clay et al. (2009), their Section 3.2 [^mc09], generalised to any atomised composition with the mean-mass rule of Lopez (2017) [^lopez2017]. With hydrogen fully ionised, the heavier atoms singly ionised, and the electrons counted among the particles, the wind's mean mass per ion and per particle, its sound speed, sonic radius, and base Jeans parameter are

$$\mu_+ = \sum_j y_j A_j, \qquad \mu_\mathrm{w} = \frac{\mu_+}{2}, \qquad c_\mathrm{s} = \sqrt{\frac{k_\mathrm{B} T_\mathrm{w}}{\mu_\mathrm{w} m_\mathrm{p}}}, \qquad R_\mathrm{s} = \frac{G M_\mathrm{p}}{2 c_\mathrm{s}^2}, \qquad \lambda_\mathrm{b} = \frac{G M_\mathrm{p}}{R_\mathrm{base}\, c_\mathrm{s}^2} \tag{4}$$

with $A_j$ the atomic mass of element $j$ in atomic mass units. The code multiplies by the proton mass where Lopez writes the hydrogen atom mass. The recombination coefficient is the mole-fraction-weighted case B set,

$$\alpha_\mathrm{B}(T_\mathrm{w}) \;=\; \sum_j y_j\, \alpha_{\mathrm{B},j}(10^4\ \mathrm{K}) \left(\frac{T_\mathrm{w}}{10^4\ \mathrm{K}}\right)^{-0.9} \tag{5}$$

with hydrogen's $2.7 \times 10^{-13}$ cm³ s⁻¹ and its $T^{-0.9}$ scaling from Murray-Clay et al. (2009), their Eq. (7) [^mc09], and archival values at $10^4$ K for He, C, N, and O, which take hydrogen's exponent as an approximation (other elements take oxygen's value).

Balancing photoionisation of the unit-optical-depth column against recombination removes the cross section and fixes the base ion density,

$$n_+ \;=\; \sqrt{\frac{F_\mathrm{XUV}\, G M_\mathrm{p}}{h\nu_0\, \alpha_\mathrm{B}\, c_\mathrm{s}^2\, R_\mathrm{base}^2}} \tag{6}$$

with $h\nu_0$ the front energy of Eq. (2). An isothermal Parker wind carries the base mass density $n_+\, \mu_+\, m_\mathrm{p}$ to the sonic point, and the rate is the mass flux there over the sonic sphere:

$$\dot{M}_\mathrm{RR} \;=\; 4\pi\, \rho_\mathrm{s}\, c_\mathrm{s}\, R_\mathrm{s}^2, \qquad \rho_\mathrm{s} = n_+\, \mu_+\, m_\mathrm{p}\; e^{\,3/2 - \lambda_\mathrm{b}} \tag{7}$$

The barometric factor $e^{\,3/2 - \lambda_\mathrm{b}}$ is the exact isothermal Bernoulli value, including the kinetic energy at the sonic point. The analytic estimate of Murray-Clay et al. omits that term and carries $e^{\,2 - \lambda_\mathrm{b}}$. When the sonic radius falls below the base, $R_\mathrm{s}$ is floored at $R_\mathrm{base}$ with $\rho_\mathrm{s}$ the base density (no barometric factor), and the state is flagged `subcritical_sonic`. The rate rises as $F_\mathrm{XUV}^{1/2}$, the analytic exponent, where the numerical models of Murray-Clay et al. give 0.6.

## The photon-limited rate

A third limit counts photons rather than energy. Where recombination is slow and every escaping particle is ionised once, each intercepted ionising photon removes at most one particle. This is Eq. (10) of Owen & Alvarez (2016), generalised from pure hydrogen [^owenalvarez], with the photon energy of Eq. (2) and the mass per ion of Eq. (4):

$$\dot{M}_\mathrm{PL} \;=\; \pi R_\mathrm{XUV}^2\, \frac{F_\mathrm{XUV}}{h\nu_0}\, \mu_+ m_\mathrm{p} \tag{8}$$

It is evaluated on the disk of the energy-limited rate, so the two share one photon budget and their ratio $\dot{M}_\mathrm{EL} / \dot{M}_\mathrm{PL} = \epsilon / \epsilon_\mathrm{PL}$ does not depend on the absorbing radius. Here

$$\epsilon_\mathrm{PL} \;=\; \frac{G M_\mathrm{p}\, K\, \mu_+ m_\mathrm{p}}{h\nu_0\, R_\mathrm{p}} \tag{9}$$

is the fraction of a photon's energy that lifting one particle out of the well costs, reported as `efficiency_photon_limit`. The photon limit binds when the efficiency exceeds it, which happens on shallow wells ([measured below](#measured-behaviour)). The [limitations page](../limitations.md) explains why the count is an estimate rather than a bound.

## The hydrodynamic candidate

The three limits combine as

$$\dot{M}_\mathrm{hyd} \;=\; \min\left(\dot{M}_\mathrm{EL},\ \dot{M}_\mathrm{RR},\ \dot{M}_\mathrm{PL}\right) \tag{10}$$

and the smallest names the sub-label, `hydrodynamic:EL`, `hydrodynamic:RR`, or `hydrodynamic:PL`. `photon_limit = False` removes $\dot{M}_\mathrm{PL}$ from the minimum and `recombination_limit = False` removes $\dot{M}_\mathrm{RR}$. With both off the wind is the energy-limited rate alone, the configuration most of the literature computes. The recombination-limited rate is still evaluated when it is off, because the collisionality switch reads its sonic point, so the switch's verdict does not move with that setting.

**What an RR label means.** The minimum can select the recombination-limited rate for two physically different reasons:

- the recombination-limited base ionisation of Eq. (6) sets the rate, which is what the label names,
- the barometric factor of Eq. (7) is tiny at large $\lambda_\mathrm{b}$, so the wind cannot carry material to the sonic point, and calling the result recombination-limited would be a category error.

The barometric factor, reported as `rr_chain['barometric_factor']`, separates them: near 1 the sonic-point density is the base density, while several decades below 1 the rate is small because the wind is throttled. The flux scaling does not separate them, since the base ion density follows $\sqrt{F_\mathrm{XUV}}$ at every $\lambda_\mathrm{b}$. The flux at which the sub-labels cross depends on the wind temperature, which enters through the sound speed, the barometric exponent, and the recombination coefficient.

## Measured behaviour

Numbers measured with this code, on the states named:

- $\epsilon_\mathrm{PL}$ is 0.034 for a 1 Earth-mass, 1.05 Earth-radius hydrogen and helium envelope, and 0.055 at 2 Earth masses and 1.3 Earth radii, so at the default efficiency of 0.1 both winds are photon-limited. A 1 Earth-mass carbon dioxide atmosphere reaches it at 0.26.
- The CO$_2$ band carries a few per cent of the cooling at the temperatures the thermostat selects, and most of it at a base of 1000 to 3000 K.

[^lopez2017]: Lopez, E. D. (2017). Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor. *Monthly Notices of the Royal Astronomical Society, 472*(1), 245–253.

[^mc09]: Murray-Clay, R. A., Chiang, E. I., & Murray, N. (2009). Atmospheric Escape From Hot Jupiters. *The Astrophysical Journal, 693*(1), 23–42. https://doi.org/10.1088/0004-637X/693/1/23

[^cp26]: Chatterjee, R. D., & Pierrehumbert, R. T. (2026). Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline. *The Astrophysical Journal, 998*(2), 236. https://doi.org/10.3847/1538-4357/ae2ffa

[^nakayama]: Nakayama, A., Ikoma, M., & Terada, N. (2022). Survival of Terrestrial N$_2$-O$_2$ Atmospheres in Violent XUV Environments through Efficient Atomic Line Radiative Cooling. *The Astrophysical Journal, 937*(2), 72. https://doi.org/10.3847/1538-4357/ac86ca

[^johnstone]: Johnstone, C. P., Güdel, M., Lammer, H., & Kislyakova, K. G. (2018). The Upper Atmospheres of Terrestrial Planets: Carbon Dioxide Cooling and the Earth's Thermospheric Evolution. *Astronomy & Astrophysics, 617*, A107. https://doi.org/10.1051/0004-6361/201832776

[^badnell]: Badnell, N. R. (2006). Radiative recombination data for modelling dynamic finite-density plasmas. *The Astrophysical Journal Supplement Series, 167*, 334. arXiv:astro-ph/0604144.

[^erkaev]: Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007). Roche lobe effects on the atmospheric loss from "Hot Jupiters". *Astronomy & Astrophysics, 472*(1), 329–334. https://doi.org/10.1051/0004-6361:20066929

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763

[^owenalvarez]: Owen, J. E., & Alvarez, M. A. (2016). UV Driven Evaporation of Close-in Planets: Energy-limited, Recombination-limited, and Photon-limited Flows. *The Astrophysical Journal, 816*(1), 34. https://doi.org/10.3847/0004-637X/816/1/34
