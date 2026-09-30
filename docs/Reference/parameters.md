# ZEPHYRUS parameter reference

This is a reference page for all parameters and constants used in ZEPHYRUS. For the physical model, see the [model overview](../Explanations/model.md). For what a dispatch call returns, see [dispatch results](results.md).

---

## Physical constants (`constants.py`)

| Name | Symbol | Value | Units |
|---|---|---|---|
| `kb` | $k_B$ | $1.380649 \times 10^{-23}$ (exact) | J K⁻¹ |
| `kb_cgs` | $k_B$ | $1.380649 \times 10^{-16}$ (exact) | erg K⁻¹ |
| `G` | $G$ | $6.6743 \times 10^{-11}$ | m³ kg⁻¹ s⁻² |
| `G_cgs` | $G_\mathrm{cgs}$ | $6.6743 \times 10^{-8}$ | cm³ g⁻¹ s⁻² |
| `c` | $c$ | $2.99792458 \times 10^{8}$ (exact) | m s⁻¹ |
| `h_planck` | $h$ | $6.62607015 \times 10^{-34}$ (exact) | J s |
| `m_p` | $m_p$ | $1.67262192369 \times 10^{-27}$ | kg |
| `amu` | $u$ | $1.66053906660 \times 10^{-27}$ | kg |
| `rate_floor` | none | $5.300219 \times 10^{-35}$ | kg s⁻¹ |

## Unit conversions (`constants.py`)

| Name | Factor | Conversion |
|---|---|---|
| `s2yr` | $1/(3600 \cdot 24 \cdot 365)$ | seconds → years |
| `erg2joule` | $10^{-7}$ | erg → J |
| `au2m` | $1.496 \times 10^{11}$ | au → m |
| `au2cm` | $1.496 \times 10^{13}$ | au → cm |
| `ergpersecondtowatt` | $10^{-7}$ | erg s⁻¹ → W |
| `ergcm2stoWm2` | $10^{-3}$ | erg s⁻¹ cm⁻² → W m⁻² |
| `ev2joule` | $1.602176634 \times 10^{-19}$ (exact) | eV → J |

---

## Sun and Earth reference values (`planets_parameters.py`)

### Sun

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Rs` | $R_\odot$ | $6.957 \times 10^{8}$ | m |
| `Ms` | $M_\odot$ | $1.98847 \times 10^{30}$ | kg |
| `Ls` | $L_\odot$ | $3.828 \times 10^{26}$ | W |
| `age_sun` | none | $4.603 \times 10^{9}$ | yr |


### Earth

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Re` | $R_\oplus$ | $6.378 \times 10^{6}$ | m |
| `Me` | $M_\oplus$ | $5.9722 \times 10^{24}$ | kg |
| `Me_atm` | none | $5.15 \times 10^{18}$ | kg |
| `Fxuv_earth_10Myr` | $F_\mathrm{XUV,\oplus}(10\,\mathrm{Myr})$ | $14.67$ | W m⁻² |
| `Fxuv_earth_today` | $F_\mathrm{XUV,\oplus}$ | $4.64 \times 10^{-3}$ | W m⁻² |
| `age_earth` | none | $4.543 \times 10^{9}$ | yr |
| `e_earth` | none | $0.017$ | dimensionless |
| `a_earth` | none | $1$ | au |

`Fxuv_earth_10Myr` is taken from Fig. 9 of Wordsworth et al. (2018) [^wordsworth2018].

### Jupiter

From the IAU 2015 nominal values (Prša et al. 2016 [^prsa]).

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Rjup` | $R_\mathrm{Jup}$ | $7.1492 \times 10^{7}$ | m |
| `Mjup` | $M_\mathrm{Jup}$ | $1.8982 \times 10^{27}$ | kg |

---

## TOI-561 reference values (`planets_parameters.py`)

### TOI-561 (star)

From Weiss et al. (2021) [^weiss].

| Name | Value | Errorbar | Units |
|---|---|---|---|
| `R_TOI561` | $0.832\,R_\odot$ | $0.019\,R_\odot$ | m |
| `M_TOI561` | $0.805\,M_\odot$ | $0.030\,M_\odot$ | kg |
| `L_TOI561` | $0.522\,L_\odot$ | $0.017\,L_\odot$ | W |
| `age_TOI561` | $10 \times 10^{9}$ | $3 \times 10^{9}$ | yr |


### TOI-561 b (planet)

From Brinkman et al. (2023) [^brinkman].

| Name | Value | Errorbar | Units |
|---|---|---|---|
| `R_TOI561b` | $1.37\,R_\oplus$ | $0.04\,R_\oplus$ | m |
| `M_TOI561b` | $2.24\,M_\oplus$ | $0.20\,M_\oplus$ | kg |
| `a_TOI561b` | $0.0106$ | $0.0004$ | au |


---

## Dispatcher settings (`dispatcher.DispatchSettings`)

The knobs of the [escape-regime framework](../Explanations/regimes.md). Every default is the documented reference choice; the criteria thresholds carry the physical bands stated in the dispatcher page, which the diagnostics report beside every verdict.

| Name | Default | Options / units | Meaning |
|---|---|---|---|
| `base_method` | `'lopez'` | `'lopez'`, `'fixed_pressure'`, `'boreas'` | How the XUV wind base is located on the profile. The Lopez (2017) [^lopez2017] level is $P_\mathrm{base} = \mu g / \sigma_{\nu_0}$, about a nanobar; `'boreas'` uses the optional BOREAS solver and falls back to `'lopez'` with a flag when it is absent or does not converge. |
| `base_out_of_range` | `'clamp'` | `'clamp'`, `'extend'` | What happens when the profile top is deeper than the physical base level: clamp to the top level (flagged, distance recorded) or evaluate the base on the extended upper structure. Whether it engages depends on the state: the Lopez base is $\mu g / \sigma_{\nu_0}$, tens of nanobars on an Earth-mass carbon dioxide planet but below a nanobar on a low-gravity hydrogen envelope, so no single profile top clears it everywhere. Read `base_clamped` rather than assuming. |
| `P_photo` | 2000 | Pa | Photospheric-type level for the energy-limited geometric factor (20 mbar, after Baumeister et al. 2023 [^baumeister]). |
| `P_base_fixed` | 5.0 | Pa | Base pressure for the `'fixed_pressure'` method only. |
| `kn_crit` | 1.0 | $> 0$ | Sonic-point Knudsen threshold of the fluid-to-kinetic switch; the physical band 0.1 to 3 is a diagnostic constant, not a knob. |
| `kn_hysteresis` | 1.5 | $\geq 1$ (1 disables the window) | Hysteresis window factor around `kn_crit`, consumed only when a previous regime label is supplied. |
| `gate` | `'neutral'` | `'neutral'`, `'plasma'` | Which escape temperature gates the hydrostatic branch; both are always computed and disagreements are flagged as contested. |
| `efficiency` | 0.1 | $0 < \epsilon \leq 1$ | Energy-limited heating efficiency $\epsilon$. |
| `efficiency_mode` | `'fixed'` | `'fixed'`, `'caldiroli'` | Fixed $\epsilon$, or the Caldiroli et al. (2022) [^caldiroli] fitted efficiency converted to the Erkaev geometry, with a guarded fallback below its validity bound. |
| `T_exo_mode` | `'prescribed'` | `'prescribed'`, `'thermostat'` | Exobase temperature source. The prescribed value is the hydrostatic branch's dominant sensitivity; the thermostat estimator is biased high by construction. |
| `T_exo_value` | 1000 | K | The prescribed exobase temperature. |
| `cool_atomic` | `True` | `True`, `False` | Atomic line cooling (H, C, C+, N, N+, O, O+) in the wind thermostat. |
| `cool_co2_band` | `True` | `True`, `False` | CO2 15 micron band cooling (deexcitation rates measured over roughly 150 to 500 K). |
| `cool_o_finestructure` | `True` | `True`, `False` | Atomic O fine-structure cooling at 63 and 147 micron. |
| `cool_recombination` | `True` | `True`, `False` | Recombination (continuum) cooling. Disabling all four channels at once is rejected. |
| `fractionate` | `True` | `True`, `False` | Apply the N-species closure on confirmed hydrodynamic verdicts; otherwise split by reservoir mass fractions. |
| `tidal` | `True` | `True`, `False` | Apply the Erkaev et al. (2007) [^erkaev] tidal factor. It divides the energy-limited rate and the interior-luminosity cap alike, so that when the residual is admitted (`residual_mode`) the two candidates it compares measure one barrier; `False` sets $K_\mathrm{tide} = 1$ in both. |
| `photon_limit` | `True` | `True`, `False` | Cap the XUV wind at one particle per intercepted ionizing photon, the photon-limited rate of Owen & Alvarez (2016) [^owenalvarez] Eq. (10), on the energy-limited rate's disk. It binds where the efficiency exceeds `efficiency_photon_limit`, on shallow wells. `False` restores min(EL, RR), which is what most energy-limited codes in the literature compute. |
| `recombination_limit` | `True` | `True`, `False` | Admit the recombination-limited rate of Murray-Clay et al. (2009) [^mc09] as a candidate for the XUV wind. `False` leaves min(EL, PL), or the energy-limited rate alone together with `photon_limit = False`. The chain is still evaluated either way, because the sonic-point Knudsen switch reads its sonic density, so the collisionality verdict does not change with this setting. |
| `nozzle_temperature` | `'photospheric'` | `'photospheric'`, `'wind'` | Which level the L1 nozzle candidate is launched from: the photospheric level at the profile's own temperature (the construction of Jackson et al. 2017 [^jackson17]) or the wind base at the thermostat's wind temperature and mean mass (the upper envelope their Figure 9 explores). Both settings take the density, temperature, and mean mass from one level, which is what the Bernoulli cancellation behind the launch-level convention requires; the wind setting places the launch level on the wind's own isothermal column, anchored at the wind base with the ideal-gas density there for the wind's temperature and mean mass, rather than carrying a cold density into a hot sound speed. That column is a device for placing the level consistently with the sound speed evaluating the barrier, not a claim about structure below the anchor, which is far hotter than the atmosphere really is there. Where the anchor itself sits is a separate and physical question, and `base_method` moves it: the two are worth keeping apart, because the level convention cancels and the anchor does not. One consequence to watch under the `lopez` default: on an inflated envelope the wind base can sit outside the planet's Roche lobe, measured at 1.29 lobe radii on a 3 Earth-mass, 2.2 Earth-radius H/He case, in which case the exponent is clamped and the candidate reports the lobe-filling boundary value with `nozzle_saturated` raised. The temperature is the model's dominant uncertainty by its authors' own statement, and it also moves the applicability criterion through the sonic radius, so the setting can decide whether the branch competes at all; the spread across the two settings is analysis, not a module output. |
| `residual_mode` | `'off'` | `'off'`, `'luminosity_capped'` | Whether the bolometric candidate competes for the rate past the activation gate. It is computed and reported on every call either way. `'luminosity_capped'` admits it, capped by the interior luminosity, which is the core-powered rate of Gupta & Schlichting (2019) [^gs19] whose persistence Tang et al. (2024) [^tang] dispute. Off by default, following Tang et al.; the price is a jump in the dispatched rate at the activation gate, measured on the [escape regimes](../Explanations/regimes.md) page together with the dispute and the band past the gate where the candidate would win. |
| `lambda_crit` | 20.0 | $> 0$ | Boil-off activation threshold on the restricted Jeans parameter (literature band 15 to 35). |
| `gamma_bates` | 0.75 | $> 0$ | Shape parameter of the Bates temperature profile of the extended upper structure. |
| `kzz` | 300 | m² s⁻¹ | Eddy diffusion coefficient when the profile carries no `kzz` column. |
| `gamma_wind` | 1.0 | $1 \leq \gamma \leq 5/3$ | Polytropic index at the sonic point (1 for an isothermal wind). |
| `hydrostatic_levels_min` | 200 | integer $\geq 2$ | First quadrature grid of the diffusion-limited supply integrals. |
| `hydrostatic_levels_max` | 3200 | integer $\geq$ `hydrostatic_levels_min` | Refinement ceiling. Reaching it without meeting the target is reported, not raised. |
| `hydrostatic_rtol` | 0.01 | $> 0$ | Target relative change in the bulk hydrostatic rate between a grid and its refinement. The integrals are first order in the log-pressure step, so that change also estimates what is left to converge. |

## Dispatcher inputs (`dispatcher.EscapeInputs`)

| Name | Units | Meaning |
|---|---|---|
| `M_p`, `R_p` | kg, m | Planet (interior) mass and radius. |
| `M_star`, `a`, `e` | kg, m, dimensionless | Stellar mass, semi-major axis, eccentricity (the Hill radius is evaluated at periapsis). |
| `T_eq` | K | Equilibrium temperature; the boil-off wind runs at $T_\mathrm{eq} / 2^{1/4}$. |
| `F_xuv`, `F_bol`, `F_int` | W m⁻² | XUV flux, bolometric instellation, and interior heat flux (the luminosity cap). `F_bol` is validated and carried with the state but no branch reads it: the bolometric branch takes its temperature from `T_eq`. It is required so that a caller assembling a state cannot omit it and then find a later version silently reading zero. |
| `kappa_photo` | m² kg⁻¹ | Photospheric opacity; the boil-off rate scales as its inverse. |
| `profile` | not applicable | `profiles.Profile`: pressure, radius, temperature, per-species mixing ratios, and mean molecular mass per level, base to top. |
| `settings` | not applicable | The `DispatchSettings` block above. |
| `prev_regime` | not applicable | Optional previous regime label; activates the hysteresis window. |
| `atm_converged` | not applicable | Optional data-quality flag, surfaced as `stale_input`. |
| `age`, `reservoirs` | s, kg | Optional; consumed only by the snapshot self-consistency screen and the unfractionated split. |
| `dt` | s | Optional; carried for the caller's supply cap, never used by the dispatcher itself. |

---

[^wordsworth2018]: Wordsworth, R. D., Schaefer, L. K., & Fischer, R. A. (2018). Redox evolution via gravitational differentiation on low-mass planets: implications for abiotic oxygen, water loss, and habitability. *The Astronomical Journal, 155*(5), 195. https://doi.org/10.3847/1538-3881/aab608

[^prsa]: Prša, A., Harmanec, P., Torres, G., et al. (2016). Nominal Values for Selected Solar and Planetary Quantities: IAU 2015 Resolution B3. *The Astronomical Journal, 152*(2), 41. https://doi.org/10.3847/0004-6256/152/2/41

[^weiss]: Weiss, L. M., Dai, F., Huber, D., et al. (2021). The TESS-Keck Survey. II. An Ultra-short-period Rocky Planet and Its Siblings Transiting the Galactic Thick-disk Star TOI-561. *The Astronomical Journal, 161*(2), 56. https://doi.org/10.3847/1538-3881/abd409

[^brinkman]: Brinkman, C. L., Weiss, L. M., Dai, F., et al. (2023). TOI-561 b: A Low-density Ultra-short-period "Rocky" Planet around a Metal-poor Star. *The Astronomical Journal, 165*(3), 88. https://doi.org/10.3847/1538-3881/acad83

[^lopez2017]: Lopez, E. D. (2017). Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor. *Monthly Notices of the Royal Astronomical Society, 472*(1), 245–253.

[^baumeister]: Baumeister, P., Tosi, N., Brachmann, C., Grenfell, J. L., & Noack, L. (2023). Redox state and interior structure control on the long-term habitability of stagnant-lid planets. *Astronomy & Astrophysics, 675*, A122. https://doi.org/10.1051/0004-6361/202245791

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763

[^erkaev]: Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007). Roche lobe effects on the atmospheric loss from "Hot Jupiters". *Astronomy & Astrophysics, 472*(1), 329–334. https://doi.org/10.1051/0004-6361:20066929

[^owenalvarez]: Owen, J. E., & Alvarez, M. A. (2016). UV Driven Evaporation of Close-in Planets: Energy-limited, Recombination-limited, and Photon-limited Flows. *The Astrophysical Journal, 816*(1), 34. https://doi.org/10.3847/0004-637X/816/1/34

[^mc09]: Murray-Clay, R. A., Chiang, E. I., & Murray, N. (2009). Atmospheric Escape From Hot Jupiters. *The Astrophysical Journal, 693*(1), 23–42. https://doi.org/10.1088/0004-637X/693/1/23

[^jackson17]: Jackson, B., Arras, P., Penev, K., Peacock, S., & Marchant, P. (2017). A new model of Roche lobe overflow for short-period gaseous planets and binary stars. *The Astrophysical Journal, 835*(2), 145. https://doi.org/10.3847/1538-4357/835/2/145

[^gs19]: Gupta, A., & Schlichting, H. E. (2019). Sculpting the valley in the radius distribution of small exoplanets as a by-product of planet formation: the core-powered mass-loss mechanism. *Monthly Notices of the Royal Astronomical Society, 487*(1), 24–33.

[^tang]: Tang, Y., et al. (2024). Assessing Core-powered Mass Loss in the Context of Early Boil-off: Minimal Long-lived Mass Loss for the Sub-Neptune Population. *The Astrophysical Journal, 976*(2), 221.
