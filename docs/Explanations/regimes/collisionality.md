# The collisionality switch

A fluid wind exists only if the gas is still collisional where it becomes supersonic. Where it is not, molecules leave one by one from an exobase and the [hydrostatic branch](hydrostatic.md) gives the rate. The switch between the two compares the mean free path with the scale height at the sonic point of the [hydrodynamic wind](wind.md). It is not applied while the [boil-off criterion](boiloff.md) holds. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The Knudsen number at the sonic point

Following Chatterjee & Pierrehumbert (2026), their Eqs. (17) and (18) [^cp26],

$$\mathrm{Kn}_\mathrm{s} \;=\; \frac{\ell_\mathrm{s}}{H_\mathrm{s}}, \qquad \ell_\mathrm{s} = \frac{1}{\sqrt{2}\,\sigma_\mathrm{C}\, n_\mathrm{s}}, \qquad H_\mathrm{s} = \frac{(1+\gamma)\, R_\mathrm{s}}{4 + \sqrt{2}\sqrt{5 - 3\gamma}}, \qquad n_\mathrm{s} = \frac{\rho_\mathrm{s}}{\mu_+ m_\mathrm{p}} \tag{1}$$

where $\ell_\mathrm{s}$ is the Maxwell mean free path, $H_\mathrm{s}$ the scale height at the sonic radius $R_\mathrm{s}$, $n_\mathrm{s}$ the heavy-particle density at the sonic point from the sonic-point density $\rho_\mathrm{s}$ of the [recombination-limited rate](wind.md#the-recombination-limited-rate), and $\gamma$ the polytropic index (`gamma_wind`, 1 for an isothermal wind). $\sigma_\mathrm{C}$ is the collision cross section of the atomised mixture at the wind temperature, their Eq. (25),

$$\sigma_\mathrm{C} \;=\; \sum_j y_j\, \sigma_j(T_\mathrm{w}) \tag{2}$$

where $y_j$ is the mole fraction of element $j$ in the atomised composition of the wind base (every molecule split into its atoms), and $\sigma_j(T_\mathrm{w})$ is the collision cross section of element $j$ at the wind temperature.

Each element's cross section comes from the first available source in a fixed order:

1. tabulated collision integrals, as the momentum-transfer cross section $\pi \sigma^2 \Omega^{(1,1)*}$ (Laricchiuta et al. 2009, their Eqs. (2) to (4) [^laricchiuta]), validated against measured viscosities
2. for hydrogen, an inversion of the diffusion coefficient (Zahnle et al. 1990, their Eq. (30) [^z90], on the compilation of Zahnle & Kasting 1986 [^zk86])
3. a geometric hard sphere, whose bias is documented and flagged.

Six elements (Al, P, Cl, K, Ca, Ti) have no tabulated radius yet and fall back on an assumed one, flagged with their own provenance class. The sources of every radius are listed on the [knudsen validation page](../../Validation/knudsen.md).

The cross sections are neutral-neutral, while the wind can be substantially ionised. This is the neutral onset, which Chatterjee & Pierrehumbert choose deliberately. Including the ion channels would shorten the mean free path and lower $\mathrm{Kn}_\mathrm{s}$, so the switch errs toward hydrostatic verdicts on an ionised wind ([limitations](../limitations.md)).

## The threshold and its band

A state with $\mathrm{Kn}_\mathrm{s}$ at or below the threshold $\mathrm{Kn}_\mathrm{c}$ (`kn_crit`) sustains the wind, and the hydrodynamic candidate is admitted. Above it, the gas decouples before reaching sonic conditions, and only the hydrostatic candidate is admitted (with the wind joining it if the exosphere is unstable, see [hydrostatic escape](hydrostatic.md#the-exosphere-stability-test)).

The default threshold is 1, and its physical band is 0.1 to 3. The two ends come from different arguments:

- Kinetic simulations place the transition near 0.1 when the heating is deposited in a sharp layer and near 1 when it is distributed (Johnson et al. 2013 [^johnson]). The lower part of the band is the physics of the heating geometry, not tuning freedom.
- The upper edge is not Johnson's number. Chatterjee & Pierrehumbert (2026) [^cp26] argue, citing the same work, that the energy limit may survive to a sonic Knudsen number of 1 to 3 or beyond, and call how far it survives an open question.

The diagnostics report the label at both band edges beside every verdict (`diagnostics['knudsen']['counterfactual_labels']`). Each is the label a dispatch with the threshold at that edge returns, so a state the switch does not decide, one with the boil-off gate open for instance, reports the same label at both edges.

## Hysteresis

For a time-stepping caller, supplying the previous result's `rate_branch` as `prev_regime` opens a hysteresis window. The threshold becomes $h\,\mathrm{Kn}_\mathrm{c}$ while leaving a hydrodynamic state and $\mathrm{Kn}_\mathrm{c}/h$ while leaving a hydrostatic one, with $h$ the setting `kn_hysteresis` (default 1.5). An evolutionary track then cannot switch back and forth between branches on numerical noise. The branch, not the label, carries the memory, because a wind renamed `roche_overflow` by the Roche screen is still a wind. After a boil-off or an L1 transfer the switch stays sharp. `hysteresis_active` and `diagnostics['knudsen']['threshold_applied']` report it.

## Measured behaviour

- The recombination-limited calculation reports base ionisation fractions reaching 0.86 on heavy compositions, which is where the neutral cross sections matter most.

[^cp26]: Chatterjee, R. D., & Pierrehumbert, R. T. (2026). Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline. *The Astrophysical Journal, 998*(2), 236. https://doi.org/10.3847/1538-4357/ae2ffa

[^laricchiuta]: Laricchiuta, A., Bruno, D., Capitelli, M., et al. (2009). High temperature Mars atmosphere. Part I: transport cross sections. *The European Physical Journal D, 54*(3), 607–612. https://doi.org/10.1140/epjd/e2009-00192-7

[^z90]: Zahnle, K., Kasting, J. F., & Pollack, J. B. (1990). Mass Fractionation of Noble Gases in Diffusion-Limited Hydrodynamic Hydrogen Escape. *Icarus, 84*(2), 502–527.

[^zk86]: Zahnle, K. J., & Kasting, J. F. (1986). Mass Fractionation during Transonic Escape and Implications for Loss of Water from Mars and Venus. *Icarus, 68*(3), 462–480.

[^johnson]: Johnson, R. E., Volkov, A. N., & Erwin, J. T. (2013). Molecular-Kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow. *The Astrophysical Journal Letters, 768*(1), L4. https://doi.org/10.1088/2041-8205/768/1/L4
