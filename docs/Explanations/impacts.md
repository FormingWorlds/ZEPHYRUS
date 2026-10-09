# Giant impacts

A giant collision removes part of the target planet's atmosphere in a single event: a shock launched by the impact accelerates atmosphere past the escape velocity, locally near the impact site and globally through ground motion. This is the second mass-loss channel of ZEPHYRUS, physically and numerically separate from the continuous escape of the [regime framework](regimes.md): it is applied per collision rather than per time step, and its rate question ("what fraction is lost in this event") replaces the continuous channel's ("how fast is mass leaving"). A regime label is reserved for routing impacts through the same interface in the future; today the caller invokes the channel directly.

ZEPHYRUS provides two scaling prescriptions through the unified entry point `zephyrus.collision.impact_loss`:

- `kegerreis2020`: the thin-atmosphere scaling law of Kegerreis et al. (2020) [^kegerreis], also evaluated directly through `zephyrus.collision.mass_loss`.
- `roche2026`: the generalized envelope-dependent scaling law of Roche et al. (2026) [^roche2026], also evaluated directly through `zephyrus.collision.mass_loss_roche2026`.

Both return the target atmospheric mass loss fraction $X \in [0, 1]$. The unified dispatcher `impact_loss` returns an `ImpactLossResult` dataclass holding the loss fraction along with diagnostic flags and active stability clamps.

## The Kegerreis et al. (2020) law

The prescription of Kegerreis et al. (2020), their Eq. 1 [^kegerreis], evaluates the eroded fraction as:

$$X \;=\; \min\!\left(0.64 \left( \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^2 \left(\frac{M_\mathrm{i}}{M_\mathrm{tot}}\right)^{1/2} \left(\frac{\rho_\mathrm{i}}{\rho_\mathrm{t}}\right)^{1/2} f_M(b) \right)^{0.65},\; 1\right) \tag{K1}$$

where $X$ is the fraction of the target's atmosphere removed (capped at 1 for total erosion), subscript $\mathrm{i}$ denotes the impactor and $\mathrm{t}$ the target, $v_\mathrm{c}$ is the speed at first contact, $M_\mathrm{tot} = M_\mathrm{i} + M_\mathrm{t}$ is the total mass, $\rho_\mathrm{i}$ and $\rho_\mathrm{t}$ are the bulk densities of the atmosphere-free bodies, and $b \equiv \sin\beta$ is the dimensionless impact parameter for impact angle $\beta$ (0 head-on, 1 fully grazing). The prefactor and exponent are least-squares fits to the paper's suite of 259 smoothed-particle-hydrodynamics (SPH) simulations, each fitted with an uncertainty of 0.01. The mutual escape speed of the pair at contact is

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t} + M_\mathrm{i})}{R_\mathrm{t} + R_\mathrm{i}}} \tag{K2}$$

with $R_\mathrm{t}$ and $R_\mathrm{i}$ the body radii at the base of any atmosphere, and $f_M(b)$ is the fractional interacting mass of the pair (their Eq. B1), built from density-weighted spherical caps of common height $d = (R_\mathrm{t} + R_\mathrm{i})(1 - b)$:

$$f_M \;=\; \frac{\rho_\mathrm{t}\, V^\mathrm{cap}_\mathrm{t} + \rho_\mathrm{i}\, V^\mathrm{cap}_\mathrm{i}}{\rho_\mathrm{t}\, V_\mathrm{t} + \rho_\mathrm{i}\, V_\mathrm{i}}, \qquad V^\mathrm{cap}_\mathrm{t,i} = \frac{\pi}{3}\, d^2 \left(3 R_\mathrm{t,i} - d\right) \tag{K3}$$

where $V_\mathrm{t,i}$ are the full body volumes. At equal bulk densities $f_M$ reduces to the fractional interacting volume of their Eq. B2. The common-height caps are a linearized bookkeeping: outside the fitted geometry, for a much denser and much smaller impactor near head-on, the raw $f_M$ can leave $[0, 1]$ and vary non-monotonically with $b$, so ZEPHYRUS clamps $f_M$ to $[0, 1]$. Within the fitted domain the clamp never engages.

Three input conventions follow the paper and must be honored by the caller: $v_\mathrm{c}$ is the speed at first contact, not the relative speed at infinity; the masses and radii exclude any atmosphere, with radii taken at the base of the atmosphere; and the densities are bulk values of the atmosphere-free bodies.

### Fitted domain and accuracy

The Kegerreis et al. (2020) law is constrained by simulations spanning target masses of roughly 0.3 to 3 Earth masses, impactor masses down to about 0.05 Earth masses, bulk densities from about half to double Earth's, contact speeds of 1 to 3 $v_\mathrm{esc}$, all impact angles, and thin atmospheres of order 1% of the planet mass. The median deviation of the simulations from the law is 9%, rising to about 20% for slow, head-on impacts, whose outcomes are chaotic. The loss depends only mildly on the atmosphere mass in this thin-atmosphere regime, with a factor of 10 less atmosphere increasing the eroded fraction by roughly 10%; substantially thicker atmospheres, which can cushion the impactor, fall outside the law's domain. The function evaluates the law for any physically valid inputs and does not warn when masses, densities, or speeds leave the fitted ranges; staying inside them is the caller's responsibility, and the [limitations page](limitations.md) lists what the channel leaves out (impactor-side atmosphere, volatile delivery, and mantle stripping among them).

The returned fraction applies to the target's atmosphere as a whole. Consistent with the bulk-removal treatment of the continuous channel's unfractionated splits, the caller partitions the lost mass across atmospheric species without elemental fractionation.

## The Roche et al. (2026) law

Roche et al. (2026) [^roche2026] generalize atmospheric erosion scaling to account explicitly for envelope mass fraction $f_\mathrm{atm} \equiv M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$. The scaling law splits atmospheric loss into near-field erosion ($X_\mathrm{NF}$, loss near the impact site) and far-field erosion ($X_\mathrm{FF}$, loss through ground motion):

$$X_\mathrm{atm} \;=\; X_\mathrm{NF} + X_\mathrm{FF} \tag{R4}$$

where $X_\mathrm{atm}$ is clamped to $[0, 1]$. The fraction of the atmosphere located in the near-field region, $f_\mathrm{NF}$, is parameterized by a generalized logistic function:

$$f_\mathrm{NF} \;=\; \frac{\zeta_4}{\left(1 + \zeta_6 \exp\!\left(\zeta_3 \left(\frac{R_\mathrm{i}^\mathrm{r}}{R_\mathrm{t}^\mathrm{r}} - \zeta_2\right)\right)\right)^{\zeta_1}} + \zeta_5 \tag{R6}$$

clamped to $[0, 1]$, and the far-field envelope fraction is $f_\mathrm{FF} = 1 - f_\mathrm{NF}$. Here $R_\mathrm{t}^\mathrm{r}$ and $R_\mathrm{i}^\mathrm{r}$ are the refractory (atmosphere-free) core-plus-mantle radii.

Near-field loss $X_\mathrm{NF}$ is described by:

$$X_\mathrm{NF} \;=\; f_\mathrm{NF} \left(\xi_1 - \xi_2\, (b + \xi_3)^2\right) \tag{R8}$$

with a velocity floor at $v_\mathrm{c} / v_\mathrm{esc} = 1.0$, clamped to $[\max(0, X_\mathrm{NF}(v_\mathrm{c} = v_\mathrm{esc})), f_\mathrm{NF}]$.

Each coefficient vector ($\boldsymbol{\zeta}$, $\boldsymbol{\xi}$, $\boldsymbol{\psi}$) is evaluated as a polynomial function of impact parameter $b$, envelope mass fraction $f_\mathrm{atm}$, refractory target mass $M_\mathrm{t}^\mathrm{r} / M_\oplus$, and refractory impactor mass fraction $\gamma \equiv M_\mathrm{i}^\mathrm{r} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{r})$:

$$\zeta_i \;=\; q_{i1} + q_{i2}\, b + q_{i3}\, b^{q_{i4}} + q_{i5}\, f_\mathrm{atm} + q_{i6}\, f_\mathrm{atm}^2 + q_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{q_{i8}} \tag{R7}$$

$$\xi_i \;=\; k_{i1} + k_{i2}\, \gamma + k_{i3}\, (\gamma + 0.05)^2 + k_{i4} \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^{k_{i5}} + k_{i6} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 1.0\right) + k_{i7} \log_{10}(f_\mathrm{atm}) \tag{R9}$$

$$\psi_i \;=\; s_{i1} + s_{i2}\, (\gamma + 0.05)^{s_{i3}} + s_{i4} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 0.05\right)^{s_{i5}} + s_{i6} \log_{10}(f_\mathrm{atm}) \tag{R11}$$

The far-field loss $X_\mathrm{FF}$ accounts for ground motion and scales with the modified specific impact energy $Q'_\mathrm{R}$:

$$X_\mathrm{FF} \;=\; f_\mathrm{FF} \left(\psi_1 \exp\!\left(-\psi_2\, Q'_\mathrm{R} \left(1 + \frac{M_\mathrm{i}^\mathrm{r}}{M_\mathrm{t}^\mathrm{r}}\right) (1 - b)^{\psi_4}\right) + \psi_3\right) \tag{R10}$$

clamped to $[0, f_\mathrm{FF}]$, with $Q'_\mathrm{R}$ expressed in $\mathrm{MJ\,kg^{-1}}$.

The forms on this page follow the authors' published code, which reproduces their Fig. 3. The printed Eq. 6 has $\zeta_3 R - \zeta_2$ in the exponent, the printed Eq. 10 has $M_\mathrm{t}^\mathrm{tot}$ in the mass ratio, and the $X_\mathrm{NF}$ floor at $v_\mathrm{c} = v_\mathrm{esc}$ is in the code and not in the paper; the printed Table C values have 4 significant digits, which is not enough, because the constant term $q_{41}$ of $\zeta_4$ (161.057) and $\zeta_5$ (-161.022) almost cancel.

### Impact energy and escape speed definitions

The mutual escape speed at contact uses the total target mass and refractory impactor mass with refractory radii:

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{r})}{R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}}} \tag{R1}$$

The modified specific impact energy $Q'_\mathrm{R}$ follows the interacting-mass formulation of Leinhardt & Stewart (2012) [^leinhardt2012] and Roche et al. (2025) [^roche2025], Eqns. 13 to 18:

$$Q'_\mathrm{R} \;=\; \frac{\mu_\alpha}{\mu}\, Q_\mathrm{R}, \qquad Q_\mathrm{R} \;=\; \frac{\mu\, v_\mathrm{c}^2}{2\,(M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{tot})}$$

where $\mu \equiv M_\mathrm{i}^\mathrm{r} M_\mathrm{t}^\mathrm{tot} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{tot})$ and $\mu_\alpha \equiv \alpha M_\mathrm{i}^\mathrm{r} M_\mathrm{t}^\mathrm{tot} / (\alpha M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{tot})$. Here $\alpha$ is the interacting mass fraction of the impactor:

$$\alpha \;=\; \frac{3 R_\mathrm{i}^\mathrm{r} l^2 - l^3}{4\,(R_\mathrm{i}^\mathrm{r})^3}$$

with $l \equiv R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r} - (R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}) b$ if $(R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}) b + R_\mathrm{i}^\mathrm{r} > R_\mathrm{t}^\mathrm{r}$, and $\alpha = 1$ when the impactor is completely intercepted.

### Symbols and units

| Symbol | Quantity | Units | Notes |
|---|---|---|---|
| $M_\mathrm{t}^\mathrm{tot}$ | Total target mass | kg | Includes envelope mass |
| $M_\mathrm{t}^\mathrm{r}, M_\mathrm{i}^\mathrm{r}$ | Refractory target and impactor masses | kg | Core plus mantle mass, excluding atmosphere |
| $M_\oplus$ | Earth mass constant | kg | $5.9722 \times 10^{24}$ kg (`zephyrus.planets_parameters.Me`) |
| $R_\mathrm{t}^\mathrm{r}, R_\mathrm{i}^\mathrm{r}$ | Refractory radii | m | Radius at base of atmosphere |
| $v_\mathrm{c}$ | Contact velocity | $\mathrm{m\,s^{-1}}$ | Speed at moment of contact |
| $v_\mathrm{esc}$ | Mutual escape speed | $\mathrm{m\,s^{-1}}$ | Calculated from total target mass and refractory impactor mass |
| $b$ | Dimensionless impact parameter | - | $\sin\beta \in [0, 1]$ |
| $\gamma$ | Impactor mass fraction | - | $M_\mathrm{i}^\mathrm{r} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{r})$ |
| $f_\mathrm{atm}$ | Target atmosphere mass fraction | - | $M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$ |
| $Q'_\mathrm{R}$ | Modified specific impact energy | $\mathrm{MJ\,kg^{-1}}$ | Interacting impact energy per unit total mass |
| $X_\mathrm{atm}$ | Atmospheric loss fraction | - | Fractional loss bounded in $[0, 1]$ |

### Coefficient origin

The scaling coefficients comprise 61 numbers (Roche et al. 2026, pp. 5-6):

- The 23 near-field mass coefficients ($q_{ij}$, $\zeta_5$, $\zeta_6$) are fitted to a separate suite of initialised SPH planets without impacts ($R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \in [0.001, 1.0]$, $M_\mathrm{t}^\mathrm{r} \in [0.01, 5.0]\,M_\oplus$, $f_\mathrm{atm} \in [0.01, 0.2]$).
- The 38 loss coefficients ($k_{ij}$, $s_{ij}$) are fitted to 300 new impact simulations (with a single target mass of about 1 $M_\oplus$, and $f_\mathrm{atm} \in \{0.01, 0.1, 0.2\}$) and the $f_\mathrm{atm} = 0.05$ simulations of Roche et al. (2025) [^roche2025], which span multiple target masses. Away from 1 $M_\oplus$, the impact data cover only $f_\mathrm{atm} = 0.05$.

All 61 values are transcribed at full precision from the published dataset [^roche2026]. Terms fixed in Table C1-C3 ($q_{37} = q_{48} = 1$ and other unlisted parameters fixed to 0) are held in internal constants.

### Calibrated domain and stability policy

The simulation suite constrains the law over the following parameter space:

- Target refractory mass: $M_\mathrm{t}^\mathrm{r} \in [0.35, 5.0]\,M_\oplus$
- Atmosphere mass fraction: $f_\mathrm{atm} \in [0.01, 0.2]$
- Impactor mass fraction: $\gamma \in [0.1, 0.5]$
- Impact parameter: $b \in [0.0, 0.9]$
- Radius ratio: $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \in [0.001, 1.0]$ (Roche et al. 2026, Sect. 3.2, p. 5, from the suite of initialised SPH planets)
- Contact velocity: $v_\mathrm{c} \in [1.0, 3.0]\,v_\mathrm{esc}$

To support planetary evolution and accretion calculations where conditions cross these empirical boundaries, `zephyrus.collision.impact_loss` applies a structured evaluation policy:

1. Diagnostic range flags. When inputs fall outside the calibrated range with a 1% relative tolerance ($f_\mathrm{atm} \notin [0.01, 0.2]$, $M_\mathrm{t}^\mathrm{r} \notin [0.35, 5.0]\,M_\oplus$, $\gamma \notin [0.1, 0.5]$, $b > 0.9$, $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \notin [0.001, 1.0]$, or $v_\mathrm{c} / v_\mathrm{esc} > 3.0$), the loss fraction is computed and the out-of-range parameter names are recorded in `ImpactLossResult.flags` (`'f_atm'`, `'M_t_earth'`, `'gamma'`, `'b'`, `'R_ratio'`, `'v_ratio'`). Speeds below the escape speed ($v_\mathrm{c} < v_\mathrm{esc}$) are governed by the velocity floor and are not flagged. An airless target ($f_\mathrm{atm} = 0$) returns no flags.
2. Stability clamping. Outside empirical stability limits, arguments to the empirical fit are evaluated at the nearest bound ($f_\mathrm{atm} \in [10^{-6}, 0.4]$, $M_\mathrm{t}^\mathrm{r} \in [10^{-3}, 10]\,M_\oplus$, $\gamma \ge 10^{-3}$) to prevent unphysical numerical divergence. Applied bounds are recorded in `diagnostics['clamped']`. Physical quantities ($v_\mathrm{esc}$, $v_\mathrm{c} / v_\mathrm{esc}$, $Q'_\mathrm{R}$, $\gamma$) and the far-field mass ratio use the physical input masses and radii.
3. Zero atmosphere. When $f_\mathrm{atm} = 0.0$, the function returns $X = 0.0$ immediately.
4. Grazing collisions. The scaling law is calibrated for impact parameters $b \in [0.0, 0.9]$. For grazing collisions at $b = 1$, the formulation evaluates continuously to $X \approx 0.0315$ for the Earth-mass benchmark case ($M_\mathrm{t}^\mathrm{r} = 1.0\,M_\oplus$, $\gamma = 0.3$, $f_\mathrm{atm} = 0.01$, $v_\mathrm{c} / v_\mathrm{esc} = 1.5$), while oblique collisions with $b > 0.9$ trigger the `'b'` diagnostic flag.

### Physical caveats

The paper states that the simulations do not model several physical processes (Sect. 4.1, p. 7): rotation, atmosphere composition (all H2-He, with the Hubbard & MacFarlane 1980 [^hubbard1980] equation of state, not an ideal gas), temperature, oceans, miscible envelopes, and thermally driven loss after the impact (the simulations follow the first tens of hours after the impact, p. 7).

## Fitted data of each law

Both laws evaluate empirical fits and extrapolate outside their calibration data:

- The `kegerreis2020` law was fitted on simulations with thin atmospheres of order 1% of planet mass, target masses of roughly 0.3 to 3 $M_\oplus$, impactor masses down to about 0.05 $M_\oplus$, bulk densities from about half to double Earth's, and contact speeds of 1 to 3 $v_\mathrm{esc}$.
- The `roche2026` law was fitted on simulations with envelope mass fractions $f_\mathrm{atm} \in [0.01, 0.2]$, target refractory masses $M_\mathrm{t}^\mathrm{r} \in [0.35, 5.0]\,M_\oplus$, impactor mass fractions $\gamma \in [0.1, 0.5]$, impact parameters $b \in [0.0, 0.9]$, radius ratios $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \in [0.001, 1.0]$, and contact speeds $v_\mathrm{c} \in [1.0, 3.0]\,v_\mathrm{esc}$. Roche et al. (2026) find that less massive atmospheres are easier to remove.

---

[^kegerreis]: Kegerreis, J. A., Eke, V. R., Catling, D. C., Massey, R. J., Teodoro, L. F. A., & Zahnle, K. J. (2020). Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density. *The Astrophysical Journal Letters, 901*(2), L31. https://doi.org/10.3847/2041-8213/abb5fb

[^roche2026]: Roche, M. J., Lock, S. J., Carter, P. J., & Leinhardt, Z. M. (2026). Giant impacts preferentially remove low-mass atmospheres: a generalised scaling law for impact-driven atmospheric loss. *The Astrophysical Journal Letters* (accepted), arXiv:2610.06077. https://doi.org/10.48550/arXiv.2610.06077. Dataset: Zenodo, https://doi.org/10.5281/zenodo.23192423.

[^roche2025]: Roche, M. J., Lock, S. J., Dou, J., Carter, P. J., Kegerreis, J. A., & Leinhardt, Z. M. (2025). Atmospheric Loss during Giant Impacts: Mechanisms and Scaling of Near- and Far-field Loss. *The Planetary Science Journal, 6*, 149. https://doi.org/10.3847/PSJ/add929

[^leinhardt2012]: Leinhardt, Z. M., & Stewart, S. T. (2012). Collisions between gravity-dominated bodies. I. Outcome regimes and scaling laws. *The Astrophysical Journal, 745*(1), 79. https://doi.org/10.1088/0004-637X/745/1/79

[^hubbard1980]: Hubbard, W. B., & MacFarlane, J. J. (1980). Structure and evolution of Uranus and Neptune. *Journal of Geophysical Research: Solid Earth, 85*(B1), 225-234. https://doi.org/10.1029/JB085iB01p00225
