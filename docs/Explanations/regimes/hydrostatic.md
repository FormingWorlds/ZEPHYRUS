# Hydrostatic escape

Where the gas is too rarefied for a wind, escape proceeds particle by particle from the exobase, the level where the mean free path first reaches the local scale height. Each species escapes at the smaller of two bottlenecks: how fast it effuses from the exobase, and how fast diffusion resupplies it from below. The branch is natively per-species. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The upper structure

The atmosphere profile usually stops well below the exobase, so the branch builds an extended upper structure above the profile top, with composition and mean molecular mass $\mu_0$ frozen at their values there. The temperature follows the Bates profile in the form Yelle (2024) uses, their Eq. (19)[^cite-yelle2024], and the radius follows from hydrostatic balance, both in the log-pressure coordinate $\zeta = \ln(P_\mathrm{top}/P)$:

$$T(\zeta) \;=\; T_\infty - \left(T_\infty - T_\mathrm{top}\right) e^{-\gamma_\mathrm{B} \zeta}, \qquad \frac{\mathrm{d}r}{\mathrm{d}\zeta} = \frac{k_\mathrm{B}\, T\, r^2}{G M_\mathrm{p}\, \mu_0}, \qquad n = \frac{P_\mathrm{top}\, e^{-\zeta}}{k_\mathrm{B} T} \tag{1}$$

with $r(\zeta = 0) = R_0$, the shape parameter $\gamma_\mathrm{B}$ (`gamma_bates`, default 0.75), and the exospheric temperature $T_\infty$ that the profile approaches (`T_exo_value`, default 1000 K). Because the Bates profile integrates analytically, so does the radius,

$$\frac{1}{r(\zeta)} \;=\; \frac{1}{R_0} - \frac{k_\mathrm{B}}{G M_\mathrm{p}\, \mu_0} \left[T_\infty\, \zeta - \left(T_\infty - T_\mathrm{top}\right) \frac{1 - e^{-\gamma_\mathrm{B} \zeta}}{\gamma_\mathrm{B}}\right],$$

so the structure carries no step error. It ends where the local Jeans parameter reaches 2, since the structure is unbound beyond that point, and the flag `extension_unbound` is raised. The exobase is the lowest point of this structure where

$$\frac{1}{\sqrt{2}\,\sigma_\mathrm{C}\, n} \;\ge\; \frac{k_\mathrm{B}\, T\, r^2}{G M_\mathrm{p}\, \mu_0} \tag{2}$$

with $\sigma_\mathrm{C}$ the mixture cross section of [the collisionality switch](collisionality.md), taken over the species at the profile top at the local temperature. This defines $R_\mathrm{exo}$, $n_\mathrm{exo}$, and $T_\mathrm{exo} = T(R_\mathrm{exo})$.

The extension is a requirement, not a refinement. The Jeans parameter at the true exobase can differ from its photospheric value by an order of magnitude, and evaluating the escape on photospheric values biases rates toward false retention by up to three orders of magnitude (Johnson et al. 2013[^cite-johnson2013]).

**The exospheric temperature.** $T_\infty$ is prescribed, and it is the branch's dominant sensitivity, because the rate depends on it exponentially. A value below $T_\mathrm{top}$ would build a thermosphere that cools with height, whose exobase is more strongly bound than its anchor, so it is raised to $T_\mathrm{top}$ and flagged. An optional estimator (`T_exo_mode = 'thermostat'`) solves the heating and cooling balance of the [wind temperature](wind.md#the-wind-temperature) at the profile top instead. A local balance without conduction is biased high by construction, so it is not the default.

## Jeans effusion

Each species $i$, of particle mass $m_i$, leaves the exobase with the Jeans effusion velocity of Yelle (2024), their Eqs. (20) and (21)[^cite-yelle2024],

$$w_{\mathrm{J},i} \;=\; \sqrt{\frac{k_\mathrm{B} T_\mathrm{exo}}{2\pi m_i}}\,(1 + \lambda_i)\, e^{-\lambda_i}, \qquad \lambda_i = \frac{G M_\mathrm{p}\, m_i}{k_\mathrm{B}\, T_\mathrm{exo}\, R_\mathrm{exo}} \tag{3}$$

where $\lambda_i$ is the species Jeans parameter. Direct simulation Monte Carlo runs find the escape rate above the equilibrium Jeans rate, because collisions above the exobase and a non-Maxwellian distribution there feed the escaping tail: about 1.7 times at $\lambda = 6$, falling to about 1.4 at $\lambda = 15$ (Volkov et al. 2011a[^cite-volkov2011a], 2011b[^cite-volkov2011b]). The code interpolates that factor linearly,

$$C(\lambda) \;=\; \begin{cases} 1.7, & \lambda \le 6 \\ 1.7 - 0.3\,(\lambda - 6)/9, & 6 < \lambda < 15 \\ 1.4, & \lambda \ge 15 \end{cases} \tag{4}$$

holding the end values outside the simulated range and flagging it (`volkov_extrapolated`). The effusion flux, per unit area of the sphere of radius $R_0$ (Yelle 2024, their Eq. (15)), is

$$\tilde{\Phi}_{\mathrm{J},i} \;=\; \left(\frac{R_\mathrm{exo}}{R_0}\right)^2 C(\lambda_i)\; w_{\mathrm{J},i}\; \tilde{X}_i(\zeta_\mathrm{exo})\; n_\mathrm{exo} \tag{5}$$

with $\tilde{X}_i(\zeta_\mathrm{exo})$ the mole fraction at the exobase from Eq. (6).

## Diffusive supply

Below the exobase each species diffuses through the frozen background with molecular diffusion coefficient $D_i = b_i / n$ and eddy diffusion coefficient $K_{zz}$ (the profile's top value, or the setting `kzz`, default 300 m² s⁻¹). The mixture diffusion parameter follows Blanc's law, $b_i = (1 - X_i) / \sum_{k \ne i} (X_k / b_{ik})$, over binary parameters $b_{ik}$ that each carry a provenance class. Thermal diffusion enters through an effective mass, Yelle (2024) Eq. (4), with the thermal diffusion factor $\alpha_\mathrm{T} = -0.25$ that Yelle adopts for light species. The mole fraction in diffusive equilibrium, their Eq. (9), is

$$\begin{gathered} \tilde{X}_i(\zeta) \;=\; X_i(\zeta = 0)\, \exp\!\left[\int_0^{\zeta} \left(1 - \frac{\tilde{m}_i}{\mu_0}\right) \frac{D_i}{D_i + K_{zz}}\, \mathrm{d}\zeta'\right], \\ \tilde{m}_i = m_i + \alpha_\mathrm{T}\, \mu_0\, \frac{\mathrm{d}\ln T}{\mathrm{d}\zeta} \end{gathered} \tag{6}$$

with $X_i(\zeta = 0)$ the species mole fraction at the profile top. The limiting flux of their Eqs. (10) and (11), the largest flux diffusion can deliver to the exobase, is

$$\Phi_{\mathrm{l},i} \;=\; \left[\int_0^{\zeta_\mathrm{exo}} \frac{k_\mathrm{B}\, T\, R_0^2}{G M_\mathrm{p}\, \mu_0\, \tilde{X}_i\, n\, (D_i + K_{zz})}\, \mathrm{d}\zeta\right]^{-1} \tag{7}$$

The exobase is located by root finding on the closed-form structure, so it does not move with the quadrature grid of these integrals. Both integrals are first order in the log-pressure step. The quadrature grid is therefore doubled from `hydrostatic_levels_min` until two successive refinements each change the bulk rate, and every element whose rate clears the rate floor, by less than `hydrostatic_rtol` (default 1%). Two, so that one chance agreement between two grids does not end the refinement, and every element, since a trace whose supply integral is still moving can sit under a bulk rate that the dominant species sets alone. Refinement also stops at `hydrostatic_levels_max`, which must be at least four times `hydrostatic_levels_min` for two changes to fit. The levels used and the last change are reported in the diagnostics, and the finest grid's values are returned, never an extrapolation.

## The combined rate

The two bottlenecks combine by the harmonic mean of Yelle (2024), their Eq. (14)[^cite-yelle2024], which lies below both and tends to the smaller:

$$\Phi_i \;=\; \frac{\tilde{\Phi}_{\mathrm{J},i}\; \Phi_{\mathrm{l},i}}{\tilde{\Phi}_{\mathrm{J},i} + \Phi_{\mathrm{l},i}} \tag{8}$$

The most abundant species at the profile top supplies itself: it takes $\tilde{X}_i = X_i(\zeta = 0)$ and the Jeans flux alone, $\Phi_i = \tilde{\Phi}_{\mathrm{J},i}$. So does any species whose Jeans rate is already below the [rate floor](../regimes.md#the-rate-floor), to save the cost of its integrals. Its rate is then an estimate, not a bound, since diffusive enrichment can raise a light species' full rate above it (2 to 5 times on the trace hydrogen and helium profiles tested), enough to carry a rate just under the floor above it. Such species are marked `pruned` in the per-species detail. Every flux on this page is per unit area of one reference sphere, the profile top at radius $R_0$, following Yelle (2024). The two bottlenecks live at different heights (effusion at the exobase, diffusion between the profile top and the exobase), and Eq. (8) can only combine them once both are per unit area of the same sphere, which is why Eq. (5) rescales the exobase flux by $(R_\mathrm{exo}/R_0)^2$. The choice of sphere does not change the result: in a steady outflow the same number of particles crosses every sphere, so the flux times the sphere's area, which is the rate, is the same at $R_0$ as at $R_\mathrm{exo}$. The hydrostatic rate is

$$\dot{M}_\mathrm{hs} \;=\; \sum_i 4\pi R_0^2\, m_i\, \Phi_i \tag{9}$$

with each species rate mapped onto elements by the mass fractions of its formula. Heavy-element rates are lower limits, since the nonthermal channels that dominate heavy-species loss from real exospheres are absent, and every hydrostatic result carries the flag `hydrostatic_lower_limit`.

## The exosphere stability test

Two escape temperatures decide whether the exobase can stay hydrostatic at all:

$$T_\mathrm{esc} \;=\; \frac{G M_\mathrm{p}\, \mu_0}{2\, k_\mathrm{B}\, R_\mathrm{exo}}, \qquad T_\mathrm{esc,+} = \frac{T_\mathrm{esc}}{2} \tag{10}$$

The neutral escape temperature $T_\mathrm{esc}$ is where thermal energy rivals binding energy. The plasma escape temperature $T_\mathrm{esc,+}$ is half of it, because in an ionised exosphere the ambipolar electric field shares each ion's binding energy with its electron (Chatterjee & Pierrehumbert 2026, their Eq. (34)[^cite-chatterjee2026]). An exobase with $T_\mathrm{exo}$ above half the selected escape temperature is unstable by their Figure 10 criterion.

In such states the [hydrodynamic candidate](wind.md#the-hydrodynamic-candidate) competes with the Jeans rate, the larger one is dispatched with its own label, and the state is flagged `gate_rerouted`. Taking the larger keeps the dispatched rate from falling as the exosphere heats past the threshold. Neither candidate describes the outflow of an unstable exosphere, as the [limitations page](../limitations.md) explains.

The setting `gate` chooses the neutral (default) or the plasma escape temperature. Both are always computed, with the ionisation fraction at the exobase under the surface-averaged flux $F_\mathrm{XUV}/4$. States where the two choices disagree are flagged `contested_ion`, with both branch rates recorded, because the ion physics that would decide them is not modelled.

[^cite-yelle2024]: R. V. Yelle, *[Diffusion limited escape of hydrogen from Mars](https://doi.org/10.1016/j.icarus.2024.116099)*, Icarus, 416, 116099, 2024. [SciX](https://scixplorer.org/abs/2024Icar..41616099Y/abstract).

[^cite-johnson2013]: R. E. Johnson, A. N. Volkov, J. T. Erwin, *[Molecular-kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow](https://doi.org/10.1088/2041-8205/768/1/L4)*, The Astrophysical Journal Letters, 768, L4, 2013. [SciX](https://scixplorer.org/abs/2013ApJ...768L...4J/abstract).

[^cite-volkov2011a]: A. N. Volkov, R. E. Johnson, O. J. Tucker, J. T. Erwin, *[Thermally Driven Atmospheric Escape: Transition from Hydrodynamic to Jeans Escape](https://doi.org/10.1088/2041-8205/729/2/l24)*, The Astrophysical Journal Letters, 729, L24, 2011. [SciX](https://scixplorer.org/abs/2011ApJ...729L..24V/abstract).

[^cite-volkov2011b]: A. N. Volkov, O. J. Tucker, J. T. Erwin, R. E. Johnson, *[Kinetic simulations of thermal escape from a single component atmosphere](https://doi.org/10.1063/1.3592253)*, Physics of Fluids, 23, 066601, 2011. [SciX](https://scixplorer.org/abs/2011PhFl...23f6601V/abstract).

[^cite-chatterjee2026]: R. D. Chatterjee, R. T. Pierrehumbert, *[Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline](https://doi.org/10.3847/1538-4357/ae2ffa)*, The Astrophysical Journal, 998, 236, 2026. [SciX](https://scixplorer.org/abs/2026ApJ...998..236C/abstract).
