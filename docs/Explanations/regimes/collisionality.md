# The collisionality switch

A fluid wind exists only if the gas is still collisional where it becomes supersonic. Where it is not, molecules leave one by one from an exobase and the [hydrostatic branch](hydrostatic.md) gives the rate. The switch between the two compares the mean free path with the scale height at the sonic point of the [hydrodynamic wind](wind.md). It is not applied while the [boil-off criterion](boiloff.md) holds. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The Knudsen number at the sonic point

Following Chatterjee & Pierrehumbert (2026), their Eqs. (17) and (18)[^cite-chatterjee2026],

$$\mathrm{Kn}_\mathrm{s} \;=\; \frac{\ell_\mathrm{s}}{H_\mathrm{s}}, \qquad \ell_\mathrm{s} = \frac{1}{\sqrt{2}\,\sigma_\mathrm{C}\, n_\mathrm{s}}, \qquad H_\mathrm{s} = \frac{(1+\gamma)\, R_\mathrm{s}}{4 + \sqrt{2}\sqrt{5 - 3\gamma}}, \qquad n_\mathrm{s} = \frac{\rho_\mathrm{s}}{\mu_+ m_\mathrm{p}} \tag{1}$$

where $\ell_\mathrm{s}$ is the Maxwell mean free path, $H_\mathrm{s}$ the scale height at the sonic radius $R_\mathrm{s}$, $n_\mathrm{s}$ the heavy-particle density at the sonic point from the sonic-point density $\rho_\mathrm{s}$ of the [recombination-limited rate](wind.md#the-recombination-limited-rate), and $\gamma$ the polytropic index (`gamma_wind`, 1 for an isothermal wind). $\sigma_\mathrm{C}$ is the collision cross section of the atomised mixture at the wind temperature, their Eq. (25),

$$\sigma_\mathrm{C} \;=\; \sum_j y_j\, \sigma_j(T_\mathrm{w}) \tag{2}$$

where $y_j$ is the mole fraction of element $j$ in the atomised composition of the wind base (every molecule split into its atoms), and $\sigma_j(T_\mathrm{w})$ is the collision cross section of element $j$ at the wind temperature.

Each element's cross section comes from the first available source in a fixed order:

1. tabulated collision integrals, as the momentum-transfer cross section $\pi \sigma^2 \Omega^{(1,1)*}$ (Laricchiuta et al. 2009, their Eqs. (2) to (4)[^cite-laricchiuta2009]), validated against measured viscosities
2. for hydrogen, an inversion of the diffusion coefficient (Zahnle et al. 1990, their Eq. (30)[^cite-zahnle1990], on the compilation of Zahnle & Kasting 1986[^cite-zahnle1986])
3. a geometric hard sphere, whose bias is documented and flagged.

Six elements (Al, P, Cl, K, Ca, Ti) have no tabulated radius yet and fall back on an assumed one, flagged with their own provenance class. The sources of every radius are listed on the [knudsen validation page](../../Validation/knudsen.md).

The cross sections are neutral-neutral, while the wind can be substantially ionised. This is the neutral onset, which Chatterjee & Pierrehumbert choose deliberately. Including the ion channels would shorten the mean free path and lower $\mathrm{Kn}_\mathrm{s}$, so the switch errs toward hydrostatic verdicts on an ionised wind ([limitations](../limitations.md)).

## The threshold and its band

A state with $\mathrm{Kn}_\mathrm{s}$ at or below the threshold $\mathrm{Kn}_\mathrm{c}$ (`kn_crit`) sustains the wind, and the hydrodynamic candidate is admitted. Above it, the gas decouples before reaching sonic conditions, and only the hydrostatic candidate is admitted (with the wind joining it if the exosphere is unstable, see [hydrostatic escape](hydrostatic.md#the-exosphere-stability-test)).

The default threshold is 1, and its physical band is 0.1 to 3. The two ends come from different arguments:

- Kinetic simulations place the transition near 0.1 when the heating is deposited in a sharp layer and near 1 when it is distributed (Johnson et al. 2013[^cite-johnson2013]). The lower part of the band is the physics of the heating geometry, not tuning freedom.
- The upper edge is not Johnson's number. Chatterjee & Pierrehumbert (2026)[^cite-chatterjee2026] argue, citing the same work, that the energy limit may survive to a sonic Knudsen number of 1 to 3 or beyond, and call how far it survives an open question.

The diagnostics report the label at both band edges beside every verdict (`diagnostics['knudsen']['counterfactual_labels']`). Each is the label a dispatch with the threshold at that edge returns when no previous branch is supplied (the L1 transfer's own fluid condition moving with it), so a state the switch does not decide, one with the boil-off gate open for instance, reports the same label at both edges.

## Hysteresis

For a time-stepping caller, supplying the previous result's `rate_branch` as `prev_regime` opens a hysteresis window. The threshold becomes $h\,\mathrm{Kn}_\mathrm{c}$ while leaving a hydrodynamic state and $\mathrm{Kn}_\mathrm{c}/h$ while leaving a hydrostatic one, with $h$ the setting `kn_hysteresis` (default 1.5). An evolutionary track then cannot switch back and forth between branches on numerical noise. The branch, not the label, carries the memory, because a wind renamed `roche_overflow` by the Roche screen is still a wind. After a boil-off or an L1 transfer the switch stays sharp. `hysteresis_active` and `diagnostics['knudsen']['threshold_applied']` report it.

## Measured behaviour

- The recombination-limited calculation reports base ionisation fractions above 0.99 on heavy compositions at high flux (a 0.1 Earth-mass nitrogen and oxygen atmosphere at $10^5$ W m⁻²), which is where the neutral cross sections matter most.

[^cite-chatterjee2026]: R. D. Chatterjee, R. T. Pierrehumbert, *[Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline](https://doi.org/10.3847/1538-4357/ae2ffa)*, The Astrophysical Journal, 998, 236, 2026. [SciX](https://scixplorer.org/abs/2026ApJ...998..236C/abstract).

[^cite-laricchiuta2009]: A. Laricchiuta, D. Bruno, M. Capitelli, C. Catalfamo, R. Celiberto, G. Colonna, P. Diomede, D. Giordano, C. Gorse, S. Longo, D. Pagano, F. Pirani, *[High temperature Mars atmosphere. Part I: transport cross sections](https://doi.org/10.1140/epjd/e2009-00192-7)*, The European Physical Journal D, 54, 607-612, 2009. [SciX](https://scixplorer.org/abs/2009EPJD...54..607L/abstract).

[^cite-zahnle1990]: K. Zahnle, J. F. Kasting, J. B. Pollack, *[Mass fractionation of noble gases in diffusion-limited hydrodynamic hydrogen escape](https://doi.org/10.1016/0019-1035(90)90050-J)*, Icarus, 84, 502-527, 1990. [SciX](https://scixplorer.org/abs/1990Icar...84..502Z/abstract).

[^cite-zahnle1986]: K. J. Zahnle, J. F. Kasting, *[Mass fractionation during transonic escape and implications for loss of water from Mars and Venus](https://doi.org/10.1016/0019-1035(86)90051-5)*, Icarus, 68, 462-480, 1986. [SciX](https://scixplorer.org/abs/1986Icar...68..462Z/abstract).

[^cite-johnson2013]: R. E. Johnson, A. N. Volkov, J. T. Erwin, *[Molecular-kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow](https://doi.org/10.1088/2041-8205/768/1/L4)*, The Astrophysical Journal Letters, 768, L4, 2013. [SciX](https://scixplorer.org/abs/2013ApJ...768L...4J/abstract).
