# ZEPHYRUS parameter reference

This is a reference page for all parameters and constants used in ZEPHYRUS. For the physical model, see the [model overview](../Explanations/model.md).

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

## Sun–Earth reference values (`planets_parameters.py`)

### Sun

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Rs` | $R_\odot$ | $6.957 \times 10^{8}$ | m |
| `Ms` | $M_\odot$ | $1.98847 \times 10^{30}$ | kg |
| `Ls` | $L_\odot$ | $3.828 \times 10^{26}$ | W |
| `age_sun` | – | $4.603 \times 10^{9}$ | yr |


### Earth

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Re` | $R_\oplus$ | $6.378 \times 10^{6}$ | m |
| `Me` | $M_\oplus$ | $5.9722 \times 10^{24}$ | kg |
| `Me_atm` | – | $5.15 \times 10^{18}$ | kg |
| `Fxuv_earth_10Myr` | $F_\mathrm{XUV,\oplus}(10\,\mathrm{Myr})$ | $14.67$ | W m⁻² |
| `Fxuv_earth_today` | $F_\mathrm{XUV,\oplus}$ | $4.64 \times 10^{-3}$ | W m⁻² |
| `age_earth` | – | $4.543 \times 10^{9}$ | yr |
| `e_earth` | – | $0.017$ | – |
| `a_earth` | – | $1$ | au |

`Fxuv_earth_10Myr` is taken from Fig. 9 of Wordsworth et al. (2018).

### Jupiter (IAU 2015 nominal values)

| Name | Symbol | Value | Units |
|---|---|---|---|
| `Rjup` | $R_\mathrm{Jup}$ | $7.1492 \times 10^{7}$ | m |
| `Mjup` | $M_\mathrm{Jup}$ | $1.8982 \times 10^{27}$ | kg |

---

## TOI-561 reference values (`planets_parameters.py`)

### TOI-561 (star, Weiss et al. 2021)

| Name | Value | Errorbar | Units |
|---|---|---|---|
| `R_TOI561` | $0.832\,R_\odot$ | $0.019\,R_\odot$ | m |
| `M_TOI561` | $0.805\,M_\odot$ | $0.030\,M_\odot$ | kg |
| `L_TOI561` | $0.522\,L_\odot$ | $0.017\,L_\odot$ | W |
| `age_TOI561` | $10 \times 10^{9}$ | $3 \times 10^{9}$ | yr |


### TOI-561 b (planet, Brinkman et al. 2023)

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
| `base_method` | `'lopez'` | `'lopez'`, `'fixed_pressure'`, `'boreas'` | How the XUV wind base is located on the profile. The Lopez (2017) level is $P_\mathrm{base} = \mu g / \sigma_{\nu_0}$, about a nanobar; `'boreas'` uses the optional BOREAS solver and falls back to `'lopez'` with a flag when it is absent or does not converge. |
| `base_out_of_range` | `'clamp'` | `'clamp'`, `'extend'` | Policy when the profile top is deeper than the physical base level: clamp to the top level (flagged, distance recorded) or evaluate the base on the extended upper structure. Profiles reaching below 1 nanobar never engage it. |
| `P_photo` | 2000 | Pa | Photospheric-type level for the energy-limited geometric factor (20 mbar, after Baumeister et al. 2023). |
| `P_base_fixed` | 5.0 | Pa | Base pressure for the `'fixed_pressure'` method only. |
| `kn_crit` | 1.0 | – | Sonic-point Knudsen threshold of the fluid-to-kinetic switch; the physical band 0.1 to 3 is a diagnostic constant, not a knob. |
| `kn_hysteresis` | 1.5 | – | Hysteresis window factor around `kn_crit`, consumed only when a previous regime label is supplied. |
| `gate` | `'neutral'` | `'neutral'`, `'plasma'` | Which escape temperature gates the hydrostatic branch; both are always computed and disagreements are flagged as contested. |
| `efficiency` | 0.1 | – | Energy-limited heating efficiency $\epsilon$. |
| `efficiency_mode` | `'fixed'` | `'fixed'`, `'caldiroli'` | Fixed $\epsilon$, or the Caldiroli et al. (2022) fitted efficiency converted to the Erkaev geometry, with a guarded fallback below its validity bound. |
| `T_exo_mode` | `'prescribed'` | `'prescribed'`, `'thermostat'` | Exobase temperature source. The prescribed value is the hydrostatic branch's dominant sensitivity; the thermostat estimator is biased high by construction. |
| `T_exo_value` | 1000 | K | The prescribed exobase temperature. |
| `cool_atomic` | `True` | – | Atomic line cooling (H, C, C+, N, N+, O, O+) in the wind thermostat. |
| `cool_co2_band` | `True` | – | CO2 15 micron band cooling (deexcitation rates measured over roughly 150 to 500 K). |
| `cool_o_finestructure` | `True` | – | Atomic O fine-structure cooling at 63 and 147 micron. |
| `cool_recombination` | `True` | – | Recombination (continuum) cooling. Disabling all four channels at once is rejected. |
| `fractionate` | `True` | – | Apply the N-species closure on confirmed hydrodynamic verdicts; otherwise split by reservoir mass fractions. |
| `tidal` | `True` | – | Apply the Erkaev et al. (2007) tidal factor to the energy-limited candidate. |
| `lambda_crit` | 20.0 | – | Boil-off activation threshold on the restricted Jeans parameter (literature band 15 to 35). |
| `gamma_bates` | 0.75 | – | Shape parameter of the Bates temperature profile of the extended upper structure. |
| `kzz` | 300 | m² s⁻¹ | Eddy diffusion coefficient when the profile carries no `kzz` column. |
| `gamma_wind` | 1.0 | – | Polytropic index at the sonic point (1 for an isothermal wind). |

## Dispatcher inputs (`dispatcher.EscapeInputs`)

| Name | Units | Meaning |
|---|---|---|
| `M_p`, `R_p` | kg, m | Planet (interior) mass and radius. |
| `M_star`, `a`, `e` | kg, m, – | Stellar mass, semi-major axis, eccentricity (the Hill radius is evaluated at periapsis). |
| `T_eq` | K | Equilibrium temperature; the boil-off wind runs at $T_\mathrm{eq} / 2^{1/4}$. |
| `F_xuv`, `F_bol`, `F_int` | W m⁻² | XUV flux, bolometric instellation, and interior heat flux (the luminosity cap). |
| `kappa_photo` | m² kg⁻¹ | Photospheric opacity; the boil-off rate scales as its inverse. |
| `profile` | – | `profiles.Profile`: pressure, radius, temperature, per-species mixing ratios, and mean molecular mass per level, base to top. |
| `settings` | – | The `DispatchSettings` block above. |
| `prev_regime` | – | Optional previous regime label; activates the hysteresis window. |
| `atm_converged` | – | Optional data-quality flag, surfaced as `stale_input`. |
| `age`, `reservoirs` | s, kg | Optional; consumed only by the snapshot self-consistency screen and the unfractionated split. |
| `dt` | s | Optional; carried for the caller's supply cap, never used by the dispatcher itself. |
