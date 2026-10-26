# Boil-off

A freshly formed or strongly heated planet can hold an atmosphere so distended that its outer layers approach the sonic point of a thermal wind. The gas then flows out on the planet's own thermal energy, with the stellar continuum keeping it near isothermal, before XUV heating matters. This page defines the criterion that recognises such a state, the rate of the resulting bolometric wind, its two caps, and what happens to the same rate once the criterion no longer holds. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The boil-off criterion

The criterion is the restricted Jeans parameter of Fossati et al. (2017) [^fossati],

$$\Lambda \;=\; \frac{G\,M_\mathrm{p}\,\mu_\mathrm{ph}}{k_\mathrm{B}\,T_\mathrm{eq}\,R_\mathrm{ph}} \tag{1}$$

the ratio of a particle's gravitational binding energy to its thermal energy, evaluated at the photospheric working level with the mean molecular mass there. The boil-off candidate competes while $\Lambda < \Lambda_\mathrm{c}$.

The level matters. The threshold below is calibrated at the photosphere. In a coupled run $R_\mathrm{p}$ is the interior radius, which on an inflated envelope sits well below the photosphere, and a parameter built there would end boil-off while the wind is still running.

The mean mass differs from the source. Fossati et al. write their parameter with the atomic hydrogen mass in place of $\mu_\mathrm{ph}$, so their $\Lambda$ and this one differ by the factor $\mu_\mathrm{ph} / m_\mathrm{H}$, and their numbers cannot be compared with these without that conversion. The composition mean mass makes the threshold one number for every composition. For isothermal gas at $T_\mathrm{eq}$, $\Lambda = 2 R_\mathrm{B}' / R_\mathrm{ph}$ identically, where $R_\mathrm{B}' = G M_\mathrm{p} \mu_\mathrm{ph} / (2 k_\mathrm{B} T_\mathrm{eq})$ is the Bondi radius at that temperature. Owen & Wu (2016) find that boil-off ends when the photosphere sits at a tenth of the Bondi radius [^owenwu], which is therefore $\Lambda = 20$ whatever the composition. That threshold $\Lambda_\mathrm{c}$ (`lambda_crit`) is the default, and the literature spread of 15 to 35 is reported as its band. Its calibration on hydrogen-rich envelopes is an assumption, which the [limitations page](../limitations.md) states.

## The bolometric wind

The wind runs at the temperature that Misener et al. (2025) recommend for the isothermal formulas [^misener], cooler than $T_\mathrm{eq}$. That temperature sets the sound speed and the sonic radius:

$$T_\mathrm{B} = \frac{T_\mathrm{eq}}{2^{1/4}}, \qquad c_\mathrm{B} = \sqrt{\frac{k_\mathrm{B} T_\mathrm{B}}{\mu_\mathrm{ph}}}, \qquad R_\mathrm{B} = \frac{G M_\mathrm{p}}{2\,c_\mathrm{B}^2} \tag{2}$$

$R_\mathrm{B}$ is larger than $R_\mathrm{B}'$ by $2^{1/4}$, so in the variable $x = R_\mathrm{ph}/R_\mathrm{B}$ the threshold $\Lambda = 20$ sits at $x = 2^{3/4} / 20 = 0.084$.

The rate is the closed-form transonic Parker wind of Owen & Wu (2016), their Eq. (6), with the Mach number at the launch level in the exact Lambert-function form of their Eq. (7) [^owenwu]:

$$\dot{M}_\mathrm{P} \;=\; \frac{4\pi\,G\,M_\mathrm{p}\,\mathcal{M}}{\kappa\,c_\mathrm{B}}, \qquad \mathcal{M} = \sqrt{-W_0\!\left(-x^{-4}\,e^{\,3 - 4/x}\right)}, \qquad x = \frac{R_\mathrm{ph}}{R_\mathrm{B}} \tag{3}$$

where $\mathcal{M}$ is the Mach number at the launch level and $W_0$ the principal branch of the Lambert function. The rate is $4\pi R_\mathrm{ph}^2 \rho u$ at a launch level of unit optical depth, where the pressure is $g/\kappa$, which is why it scales as $1/\kappa$. In a coupled run $\kappa$ should be the mean opacity the atmosphere model finds at the photospheric level. The exact form is used rather than their small-$x$ approximation, which absorbs an order-unity prefactor.

At $x = 1$ the launch level is sonic and $\mathcal{M} = 1$. A launch level outside the Bondi radius is clamped to $x = 1$ and flagged `bondi_inflated`. For small $x$ the Mach number, and with it the rate, falls exponentially, which is the physical end of boil-off.

## The Bondi cap

The rate is capped by the supply the Bondi sphere admits, in the form of Misener et al. (2025), their Eq. (10) [^misener]:

$$\dot{M}_\mathrm{B} \;=\; 4\pi R_\mathrm{B}^2\, c_\mathrm{B}\, \rho_\mathrm{ph}\, \exp\!\left(2 - \frac{2 R_\mathrm{B}}{R_\mathrm{ph}}\right) \tag{4}$$

This is the launch density carried hydrostatically to the sonic point and multiplied by the sonic speed over the sonic sphere, at the same wind temperature as Eq. (3), the temperature Misener et al. recommend for it. The launch level stands in for their radiative-convective boundary.

Gupta & Schlichting (2020), their Eq. (10), write the cap without the factor $e^2$ and at $T_\mathrm{eq}$ [^gs20]. Their form is larger than Eq. (4) by $e^{-2}\, 2^{-3/8} \exp[\Lambda (2^{1/4} - 1)]$, a factor of 4.6 at $\Lambda = 20$.

Against the Parker rate the cap is

$$\frac{\dot{M}_\mathrm{B}}{\dot{M}_\mathrm{P}} = \tau_\mathrm{ph}\, \frac{T_\mathrm{B}}{T_\mathrm{ph}}\, \exp\!\left(\frac{1}{2} - \frac{\mathcal{M}^2}{2}\right), \qquad \tau_\mathrm{ph} = \frac{\kappa P_\mathrm{ph} R_\mathrm{ph}^2}{G M_\mathrm{p}}$$

with $\tau_\mathrm{ph}$ the plane-parallel optical depth of the launch level to the supplied opacity. The cap therefore binds only on a launch level optically thin to that opacity, below an optical depth of about 0.6 on a level at the wind temperature. `tau_launch` and `binding_cap` in `diagnostics['bolometric']` say which case a state is in.

## After boil-off: the residual

Once $\Lambda \ge \Lambda_\mathrm{c}$ the atmosphere has contracted. The same wind is still computed, with one more cap: the outflow can no longer draw on the envelope's own inflation and is limited by the heat the interior supplies. Dividing the interior luminosity by the work per unit mass needed to lift gas out of the well gives the cap of Gupta & Schlichting (2019), the cooling-luminosity term of their Eq. (8) [^gs19]:

$$\dot{M}_\mathrm{E} \;=\; \frac{L_\mathrm{int}}{g_\mathrm{p}\,R_\mathrm{p}\,K}, \qquad L_\mathrm{int} = 4\pi R_\mathrm{p}^2 F_\mathrm{int}, \qquad g_\mathrm{p} = \frac{G M_\mathrm{p}}{R_\mathrm{p}^2} \tag{5}$$

with the tidal factor $K$ of [Eq. (2) of the overview](../regimes.md#tidal-geometry). The denominator is the barrier the energy-limited rate divides by, measured from the same radius. A planet close to filling its Roche lobe has a shallower barrier, so the same interior heat lifts more gas, and an admitted residual and the XUV rate it competes with measure one barrier between them.

The bolometric candidate is the minimum of the caps in force,

$$\dot{M}_\mathrm{bol} \;=\; \begin{cases} \min\left(\dot{M}_\mathrm{P},\ \dot{M}_\mathrm{B}\right), & \Lambda < \Lambda_\mathrm{c} \\ \min\left(\dot{M}_\mathrm{P},\ \dot{M}_\mathrm{B},\ \dot{M}_\mathrm{E}\right), & \Lambda \ge \Lambda_\mathrm{c} \end{cases} \tag{6}$$

and is reported in `diagnostics['bolometric']` on every call. Past the threshold it is called the residual. Whether the residual competes for the rate is the `residual_mode` setting:

- `'off'` (default): it is reported and does not compete.
- `'luminosity_capped'`: it enters the final comparison, and if it is the largest admitted rate it takes the rate and the label `boiloff`, with the flag `bolometric_residual`.

`competes` in the same diagnostics group says which setting was in force.

### Why the residual is off by default

Whether a bolometric wind persists after contraction is disputed, and the default takes the side that dispatches less.

- Gupta & Schlichting (2019) find that a bolometric wind fed by the cooling interior, core-powered mass loss, persists for gigayears, at the smaller of the Bondi-limited rate and the cooling-luminosity rate of Eq. (5) [^gs19].
- Tang et al. (2024) find that once boil-off is initialised self-consistently, the same wind removes at most a tenth of a percent of the envelope over gigayears [^tang]. They attribute the persistent rate to coupling the wind to the core luminosity, to a missing energy-loss term in the wind, and to initial conditions set before boil-off.

Here the dispute matters only in a narrow range past the threshold, where the closed-form rate of Eq. (3) has not yet shut off ([measured below](#measured-behaviour)). The default follows Tang et al.: it reports the residual and does not dispatch it. The boil-off criterion is then a jump in the dispatched rate wherever the bolometric rate is still the larger one at the threshold, and continuous where the XUV wind has already overtaken it. A run that wants the rate continuous sets `residual_mode = 'luminosity_capped'`, and the two settings differ by the candidate rate that the diagnostics print either way.

The termination timescale of Tang et al. (their Eq. 8) is computed on the closed-form wind $\min(\dot{M}_\mathrm{P}, \dot{M}_\mathrm{B})$, before the luminosity cap, since that cap and their cooling time are set by the same luminosity. It is reported beside the rate and never used as a criterion, so whichever setting is in force, the reader can see how long the wind would survive under their criterion.

## Competition with the XUV wind

While $\Lambda < \Lambda_\mathrm{c}$, the [hydrodynamic wind](wind.md) competes with the bolometric wind and the larger rate is dispatched with its own label. Owen & Schlichting (2024) give the mass loss to photoevaporation when XUV photons are absorbed inside the Bondi radius and drive the stronger outflow, and to the bolometric wind otherwise [^owensch]. The dispatcher applies the rate comparison and does not test where the photons are absorbed, a simplification the [limitations page](../limitations.md) states. The collisionality switch is not applied while the criterion holds, because a bolometric outflow is penetrated by XUV photons before it turns collisionless (their Eq. 7).

## Measured behaviour

Numbers measured with this code, on the states named:

- On a 3 Earth-mass hydrogen envelope at 0.1 au under $F_\mathrm{XUV} = 10$ W m⁻², the residual is set by the Parker rate of Eq. (3) almost everywhere past the threshold, since the closed form is still shutting off. It exceeds the XUV rate only while $\Lambda$ stays below about 21 to 23 (lower for a cooler envelope), by up to a factor of about 9 just past the threshold at 1000 to 1500 K, and has fallen to 3% to 7% of it at $\Lambda = 25$.
- At the threshold itself, on a 3 Earth-mass hydrogen and helium envelope at 1000 K and 0.0775 au, the default setting makes the dispatched rate drop by a factor of 5.2 under $F_\mathrm{XUV} = 10$ W m⁻² and 518 under 0.1 W m⁻², where the XUV rate is small (3.0 and 303 with `photon_limit = False`).

[^fossati]: Fossati, L., et al. (2017). Aeronomical constraints to the minimum mass and maximum radius of hot low-mass planets. *Astronomy & Astrophysics, 598*, A90.

[^owenwu]: Owen, J. E., & Wu, Y. (2016). Atmospheres of low-mass planets: the "boil-off". *The Astrophysical Journal, 817*(2), 107.

[^misener]: Misener, W., et al. (2025). Blowin' in the Nonisothermal Wind: Core-powered Mass Loss with Hydrodynamic Radiative Transfer. *The Astrophysical Journal, 980*(1), 152.

[^gs19]: Gupta, A., & Schlichting, H. E. (2019). Sculpting the valley in the radius distribution of small exoplanets as a by-product of planet formation: the core-powered mass-loss mechanism. *Monthly Notices of the Royal Astronomical Society, 487*(1), 24–33.

[^gs20]: Gupta, A., & Schlichting, H. E. (2020). Signatures of the core-powered mass-loss mechanism in the exoplanet population: dependence on stellar properties and observational predictions. *Monthly Notices of the Royal Astronomical Society, 493*(1), 792–806.

[^tang]: Tang, Y., et al. (2024). Assessing Core-powered Mass Loss in the Context of Early Boil-off: Minimal Long-lived Mass Loss for the Sub-Neptune Population. *The Astrophysical Journal, 976*(2), 221.

[^owensch]: Owen, J. E., & Schlichting, H. E. (2024). Mapping out the parameter space for photoevaporation and core-powered mass-loss. *Monthly Notices of the Royal Astronomical Society, 528*(2), 1615–1629.
