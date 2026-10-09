# Giant impacts

A giant collision removes part of the target planet's atmosphere in a single event: a shock launched by the impact accelerates atmosphere past the escape velocity, locally near the impact site and globally through ground motion. This is the second mass-loss channel of ZEPHYRUS, physically and numerically separate from the continuous escape of the [regime framework](regimes.md): it is applied per collision rather than per time step, and its rate question ("what fraction is lost in this event") replaces the continuous channel's ("how fast is mass leaving").

ZEPHYRUS provides two scaling prescriptions through the unified entry point `zephyrus.collision.impact_loss`:

- `kegerreis2020`: the thin-atmosphere scaling law of Kegerreis et al. (2020) [^kegerreis], also evaluated directly through `zephyrus.collision.mass_loss`.
- `roche2026`: the generalized envelope-dependent scaling law of Roche et al. (2026) [^roche2026], also evaluated directly through `zephyrus.collision.mass_loss_roche2026`.

Both return the target atmospheric mass loss fraction $X \in [0, 1]$. The unified dispatcher `impact_loss` returns an `ImpactLossResult` dataclass holding the loss fraction along with diagnostic flags and active stability clamps.

## The Kegerreis et al. (2020) law

The prescription of Kegerreis et al. (2020), their Eq. 1 [^kegerreis], evaluates the eroded fraction as:

$$X \;=\; \min\!\left(0.64 \left( \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^2 \left(\frac{M_\mathrm{i}}{M_\mathrm{tot}}\right)^{1/2} \left(\frac{\rho_\mathrm{i}}{\rho_\mathrm{t}}\right)^{1/2} f_M(b) \right)^{0.65},\; 1\right) \tag{1}$$

where $X$ is the fraction of the target atmosphere removed, subscript $\mathrm{i}$ denotes the impactor and $\mathrm{t}$ the target, $v_\mathrm{c}$ is the speed at first contact, $M_\mathrm{tot} = M_\mathrm{i} + M_\mathrm{t}$ is the total mass, $\rho_\mathrm{i}$ and $\rho_\mathrm{t}$ are the bulk densities of the atmosphere-free bodies, and $b \equiv \sin\beta$ is the dimensionless impact parameter for impact angle $\beta$ (0 head-on, 1 fully grazing). The prefactor and exponent are least-squares fits to the paper's suite of 259 smoothed-particle-hydrodynamics (SPH) simulations, each fitted with an uncertainty of 0.01. The mutual escape speed of the pair at contact is

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t} + M_\mathrm{i})}{R_\mathrm{t} + R_\mathrm{i}}} \tag{2}$$

with $R_\mathrm{t}$ and $R_\mathrm{i}$ the body radii at the base of any atmosphere, and $f_M(b)$ is the fractional interacting mass of the pair (their Eq. B1), built from density-weighted spherical caps of common height $d = (R_\mathrm{t} + R_\mathrm{i})(1 - b)$:

$$f_M \;=\; \frac{\rho_\mathrm{t}\, V^\mathrm{cap}_\mathrm{t} + \rho_\mathrm{i}\, V^\mathrm{cap}_\mathrm{i}}{\rho_\mathrm{t}\, V_\mathrm{t} + \rho_\mathrm{i}\, V_\mathrm{i}}, \qquad V^\mathrm{cap}_\mathrm{t,i} = \frac{\pi}{3}\, d^2 \left(3 R_\mathrm{t,i} - d\right) \tag{3}$$

where $V_\mathrm{t,i}$ are the full body volumes. At equal bulk densities $f_M$ reduces to the fractional interacting volume of their Eq. B2. The common-height caps are a linearized bookkeeping: outside the fitted geometry, for a much denser and much smaller impactor near head-on, the raw $f_M$ can leave $[0, 1]$ and vary non-monotonically with $b$, so ZEPHYRUS clamps $f_M$ to $[0, 1]$. Within the fitted domain the clamp never engages.

Three input conventions follow the paper and must be honored by the caller: $v_\mathrm{c}$ is the speed at first contact, not the relative speed at infinity; the masses and radii exclude any atmosphere, with radii taken at the base of the atmosphere; and the densities are bulk values of the atmosphere-free bodies.

### Fitted domain and accuracy

The Kegerreis et al. (2020) law is constrained by simulations spanning target masses of roughly 0.3 to 3 Earth masses, impactor masses down to about 0.05 Earth masses, bulk densities from about half to double Earth's, contact speeds of 1 to 3 $v_\mathrm{esc}$, all impact angles, and thin atmospheres of order 1% of the planet mass. The median deviation of the simulations from the law is 9%, rising to about 20% for slow, head-on impacts, whose outcomes are chaotic. The loss depends only mildly on the atmosphere mass in this thin-atmosphere regime, with a factor of 10 less atmosphere increasing the eroded fraction by roughly 10%; substantially thicker atmospheres, which can cushion the impactor, fall outside the law's domain.

## The Roche et al. (2026) law

Roche et al. (2026) [^roche2026] generalize atmospheric erosion scaling to account explicitly for envelope mass fraction $f_\mathrm{atm} \equiv M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$. The scaling law de-convolves atmospheric loss into near-field erosion ($X_\mathrm{NF}$, driven by ejecta plumes and atmospheric shock expansion near the impact site) and far-field erosion ($X_\mathrm{FF}$, driven by shock acceleration propagating through the planetary mantle):

$$X_\mathrm{atm} \;=\; f_\mathrm{NF}\, X_\mathrm{NF} + (1 - f_\mathrm{NF})\, X_\mathrm{FF} \tag{4}$$

where $X_\mathrm{atm}$ is clamped to $[0, 1]$. The fraction of the atmosphere located in the near-field region, $f_\mathrm{NF}$, is parameterized by a logistic sigmoid:

$$f_\mathrm{NF} \;=\; \frac{1}{1 + \exp\!\left(-\left(\zeta_1 + \zeta_2\,\gamma^{\zeta_3} + \zeta_4\left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^{\zeta_5} + \zeta_6\left(\frac{R_\mathrm{i}^\mathrm{r}}{R_\mathrm{t}^\mathrm{r}}\right)\right)\right)} \tag{5}$$

where $\gamma \equiv M_\mathrm{i}^\mathrm{tot} / M_\mathrm{t}^\mathrm{tot}$ is the impactor-to-target mass ratio, and $R_\mathrm{t}^\mathrm{r}$ and $R_\mathrm{i}^\mathrm{r}$ are the refractory (atmosphere-free) core-plus-mantle radii.

Near-field loss efficiency $X_\mathrm{NF}$ is described by:

$$X_\mathrm{NF} \;=\; \xi_1 + \xi_2\, b^{\xi_3} + \xi_4\, \gamma^{\xi_5} + \xi_6 \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^{\xi_7} \tag{6}$$

where the velocity ratio is subject to a physical near-field floor: when $v_\mathrm{c} / v_\mathrm{esc} < 1.0$, the velocity factor evaluates at $v_\mathrm{c} / v_\mathrm{esc} = 1.0$.

Each coefficient vector ($\boldsymbol{\zeta}$, $\boldsymbol{\xi}$, $\boldsymbol{\psi}$) is evaluated as a polynomial function of impact parameter $b$, envelope mass fraction $f_\mathrm{atm}$, and refractory target mass $M_\mathrm{t}^\mathrm{r} / M_\oplus$:

$$\zeta_i \;=\; q_{i1} + q_{i2}\, b + q_{i3}\, b^{q_{i4}} + q_{i5}\, f_\mathrm{atm} + q_{i6}\, f_\mathrm{atm}^2 + q_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{q_{i8}} \tag{7}$$

$$\xi_i \;=\; k_{i1} + k_{i2}\, b + k_{i3}\, b^{k_{i4}} + k_{i5}\, f_\mathrm{atm} + k_{i6}\, f_\mathrm{atm}^2 + k_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{k_{i8}} \tag{8}$$

$$\psi_i \;=\; s_{i1} + s_{i2}\, b + s_{i3}\, b^{s_{i4}} + s_{i5}\, f_\mathrm{atm} + s_{i6}\, f_\mathrm{atm}^2 + s_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{s_{i8}} \tag{9}$$

The far-field loss efficiency $X_\mathrm{FF}$ accounts for ground-shock transmission and scales with the modified specific impact energy $Q'_\mathrm{R}$:

$$X_\mathrm{FF} \;=\; \psi_1 + \psi_2\, b^{\psi_3} + \psi_4\, (\gamma + 0.05)^{\psi_5} + \psi_6 \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 0.05\right)^{\psi_7} + \psi_8 \left(\frac{Q'_\mathrm{R}}{Q''_\mathrm{R}}\right)^{\psi_9} \tag{10}$$

The energy denominator $Q''_\mathrm{R}$ normalizes the shock energy:

$$\frac{Q''_\mathrm{R}}{\mathrm{MJ\,kg^{-1}}} \;=\; p_1 \left(1 + p_2 \left(\frac{M_\mathrm{t}^\mathrm{tot}}{M_\oplus}\right)^{p_3}\right) \left(1 + Q'_\mathrm{R}\,(1 + m_\mathrm{ratio})\,(1 - b)^{p_4}\right)^{p_5} \tag{11}$$

where $m_\mathrm{ratio} \equiv M_\mathrm{i}^\mathrm{tot} / M_\mathrm{t}^\mathrm{tot}$. For grazing trajectories ($b \ge 1.0$), the geometry factor evaluates continuously to zero ($(1 - b)^{p_4} \to 0$).

### Impact energy and escape speed definitions

The mutual escape speed at contact uses total planetary masses and refractory radii:

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{tot})}{R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}}} \tag{12}$$

The modified specific impact energy $Q'_\mathrm{R}$ follows the interacting-mass formulation of Leinhardt & Stewart (2012) [^leinhardt2012] and Roche et al. (2025) [^roche2025]:

$$Q'_\mathrm{R} \;=\; \frac{1}{2}\,\frac{M_\mathrm{i}^\mathrm{tot}\, M_\mathrm{t}^\mathrm{tot}}{(M_\mathrm{i}^\mathrm{tot} + M_\mathrm{t}^\mathrm{tot})^2}\,\frac{M_\mathrm{i,interact}}{M_\mathrm{i}^\mathrm{tot}}\, v_\mathrm{c}^2 \tag{13}$$

For head-on and small-angle collisions ($b < b_\mathrm{crit} \equiv (R_\mathrm{t}^\mathrm{r} - R_\mathrm{i}^\mathrm{r}) / (R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r})$), the entire impactor enters the interacting volume so $M_\mathrm{i,interact} = M_\mathrm{i}^\mathrm{tot}$. For oblique impacts ($b \ge b_\mathrm{crit}$), the intersecting mass is calculated from the spherical cap intersection:

$$\frac{M_\mathrm{i,interact}}{M_\mathrm{i}^\mathrm{tot}} \;=\; \frac{3 R_\mathrm{i}^\mathrm{r} l^2 - l^3}{4\,(R_\mathrm{i}^\mathrm{r})^3}, \qquad l \equiv (R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r})(1 - b) \tag{14}$$

When bodies pass without physical contact ($b \ge 1.0$), $M_\mathrm{i,interact} = 0$ and $Q'_\mathrm{R} = 0$. In Eq. (11), $Q'_\mathrm{R}$ enters in units of $\mathrm{MJ\,kg^{-1}}$ ($10^{-6}\,\mathrm{J\,kg^{-1}}$).

### Symbols and units

| Symbol | Quantity | Units | Notes |
|---|---|---|---|
| $M_\mathrm{t}^\mathrm{tot}, M_\mathrm{i}^\mathrm{tot}$ | Total target and impactor masses | kg | Includes envelope mass |
| $M_\mathrm{t}^\mathrm{r}, M_\mathrm{i}^\mathrm{r}$ | Refractory target and impactor masses | kg | Core plus mantle mass, excluding atmosphere |
| $M_\oplus$ | Earth mass constant | kg | $5.972 \times 10^{24}$ kg (`zephyrus.constants.M_earth`) |
| $R_\mathrm{t}^\mathrm{r}, R_\mathrm{i}^\mathrm{r}$ | Refractory radii | m | Radius at base of atmosphere |
| $v_\mathrm{c}$ | Contact velocity | $\mathrm{m\,s^{-1}}$ | Speed at moment of contact |
| $v_\mathrm{esc}$ | Mutual escape speed | $\mathrm{m\,s^{-1}}$ | Calculated from total masses and refractory radii |
| $b$ | Dimensionless impact parameter | - | $\sin\beta \in [0, 1]$ |
| $\gamma$ | Mass ratio | - | $M_\mathrm{i}^\mathrm{tot} / M_\mathrm{t}^\mathrm{tot}$ |
| $f_\mathrm{atm}$ | Target atmosphere mass fraction | - | $M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$ |
| $Q'_\mathrm{R}$ | Modified specific impact energy | $\mathrm{J\,kg^{-1}}$ | Converted to $\mathrm{MJ\,kg^{-1}}$ in Eq. (11) |
| $X_\mathrm{atm}$ | Atmospheric loss fraction | - | Fractional loss bounded in $[0, 1]$ |

### Coefficient origin

The scaling coefficients comprise 61 fitted values derived from 300 three-dimensional SPH simulations across $f_\mathrm{atm} \in \{0.01, 0.1, 0.2\}$ combined with the 5% atmosphere simulation suite of Roche et al. (2025) [^roche2025]. All 61 values are transcribed at full precision from the published dataset [^roche2026]. Unlisted parameters from Table C of Roche et al. (2026) are fixed analytical constants: $q_{37} = 1$, $q_{48} = 1$, and all remaining unlisted terms are 0.

### Calibrated domain and stability policy

The simulation suite constrains the law over the following parameter space:

- Target refractory mass: $M_\mathrm{t}^\mathrm{r} \in [0.35, 5.0]\,M_\oplus$
- Atmosphere mass fraction: $f_\mathrm{atm} \in [0.01, 0.2]$
- Mass ratio: $\gamma \in [0.05, 0.5]$
- Impact parameter: $b \in [0.0, 0.9]$
- Contact velocity: $v_\mathrm{c} \in [1.0, 3.0]\,v_\mathrm{esc}$

To support planetary evolution and accretion calculations where conditions cross these empirical boundaries, `zephyrus.collision.impact_loss` applies a structured evaluation policy:

1. Diagnostic range flags. When inputs fall outside the calibrated range ($f_\mathrm{atm} \notin [0.01, 0.2]$, $M_\mathrm{t}^\mathrm{r} \notin [0.35, 5.0]\,M_\oplus$, $\gamma \notin [0.1, 0.5]$, $b > 0.9$, or $v_\mathrm{c} / v_\mathrm{esc} > 3.0$), the loss fraction is computed and the out-of-range parameter names are recorded in `ImpactLossResult.flags` (`'f_atm'`, `'M_t_earth'`, `'gamma'`, `'b'`, `'v_ratio'`).
2. Stability clamping. Outside empirical stability limits, arguments to the empirical fit are evaluated at the nearest bound ($f_\mathrm{atm} \in [10^{-6}, 0.4]$, $M_\mathrm{t}^\mathrm{r} \in [10^{-3}, 10]\,M_\oplus$, $\gamma \ge 10^{-3}$) to prevent unphysical numerical divergence. Applied bounds are recorded in `diagnostics['clamped']`. Physical quantities ($v_\mathrm{esc}$, $v_\mathrm{c} / v_\mathrm{esc}$, $Q'_\mathrm{R}$, $\gamma$) and the far-field mass ratio use the physical input masses and radii.
3. Zero atmosphere. When $f_\mathrm{atm} = 0.0$, the function returns $X = 0.0$ immediately.
4. Grazing continuity. For grazing collisions ($b \ge 1.0$), the mathematical formulation evaluates continuously to $X \approx 0.0315$ without singularity or divergence.

### Physical caveats

The underlying hydrodynamical simulations assume hydrogen-helium gas mixtures modeled with ideal gas equations of state. Heavy secondary atmospheres (steam, carbon dioxide, nitrogen) feature higher mean molecular weights and different shock impedance contrasts with the underlying mantle, which can alter ground-motion coupling and breakout velocities.

Colliding bodies are simulated without initial spin. Rapid retrograde or prograde planetary rotation modifies the effective surface gravity and shock breakout geometry.

The law predicts prompt hydrodynamic atmospheric ejection within hours after the collision. Prolonged thermal escape driven by magma ocean outgassing and high post-impact surface temperatures over kiloyear timescales is not included.

## Choosing between scaling laws

The `kegerreis2020` prescription is suitable for collisions into planets with thin atmospheres ($f_\mathrm{atm} \sim 0.01$) where density contrasts between impactor and target bodies are significant, or where an analytical single-bracket closed form is preferred.

The `roche2026` prescription is suitable for accretion and planetary evolution models where atmospheric envelope mass fractions vary dynamically between $0.001$ and $0.20$. Because thinner envelopes are stripped more efficiently by mantle shocks, using `roche2026` captures the enhanced vulnerability of small primordial atmospheres to complete erosion during giant impacts.

---

[^kegerreis]: Kegerreis, J. A., Eke, V. R., Catling, D. C., Massey, R. J., Teodoro, L. F. A., & Zahnle, K. J. (2020). Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density. *The Astrophysical Journal Letters, 901*(2), L31. https://doi.org/10.3847/2041-8213/abb5fb

[^roche2026]: Roche, M. J., Lock, S. J., Carter, P. J., & Leinhardt, Z. M. (2026). Giant impacts preferentially remove low-mass atmospheres: a generalised scaling law for impact-driven atmospheric loss. *The Astrophysical Journal Letters* (accepted), arXiv:2610.06077. https://doi.org/10.48550/arXiv.2610.06077. Dataset: Zenodo, https://doi.org/10.5281/zenodo.23192423.

[^roche2025]: Roche, M. J., Lock, S. J., Dou, J., Carter, P. J., & Leinhardt, Z. M. (2025). Giant impacts preferentially remove low-mass atmospheres: scaling laws for 5% atmospheres. *The Planetary Science Journal, 6*(6), 149. https://doi.org/10.3847/PSJ/add929

[^leinhardt2012]: Leinhardt, Z. M., & Stewart, S. T. (2012). Collisions between gravity-dominated bodies. I. Outcome regimes and scaling laws. *The Astrophysical Journal, 745*(1), 79. https://doi.org/10.1088/0004-637X/745/1/79
