# Giant impacts

A giant collision removes part of the target planet's atmosphere in a single event: a shock launched by the impact accelerates atmosphere past the escape velocity, locally near the impact site and globally through ground motion. This is the second mass-loss channel of ZEPHYRUS, physically and numerically separate from the continuous escape of the [regime framework](regimes.md): it is applied per collision rather than per time step, and its rate question ("what fraction is lost in this event") replaces the continuous channel's ("how fast is mass leaving"). A regime label is reserved for routing impacts through the same interface in the future; today the caller invokes the channel directly.

ZEPHYRUS provides two scaling laws for the eroded fraction, selected with `zephyrus.collision.impact_loss(law=...)`, where `law` is `'kegerreis2020'` (thin atmospheres) or `'roche2026'` (envelopes of 1% to 20% of the planet mass). The function returns an `ImpactLossResult` dataclass with `fraction`, `flags`, and `diagnostics`; each law is also callable directly, as `mass_loss` and `mass_loss_roche2026`.

## The Kegerreis et al. (2020) law

`zephyrus.collision.mass_loss` implements the scaling law of Kegerreis et al. (2020), their Eq. 1 [^kegerreis]:

$$X \;=\; \min\!\left\{0.64 \left[ \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^2 \left(\frac{M_\mathrm{i}}{M_\mathrm{tot}}\right)^{1/2} \left(\frac{\rho_\mathrm{i}}{\rho_\mathrm{t}}\right)^{1/2} f_M(b) \right]^{0.65},\; 1\right\} \tag{1}$$

where $X$ is the fraction of the target's atmosphere removed (capped at 1 for total erosion), subscript $\mathrm{i}$ denotes the impactor and $\mathrm{t}$ the target, $v_\mathrm{c}$ is the speed at first contact, $M_\mathrm{tot} = M_\mathrm{i} + M_\mathrm{t}$ is the total mass, $\rho_\mathrm{i}$ and $\rho_\mathrm{t}$ are the bulk densities of the atmosphere-free bodies, and $b \equiv \sin\beta$ is the dimensionless impact parameter for impact angle $\beta$ (0 head-on, 1 fully grazing). The prefactor and exponent are least-squares fits to the paper's suite of 259 smoothed-particle-hydrodynamics simulations, each fitted with an uncertainty of 0.01. The mutual escape speed of the pair at contact is

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t} + M_\mathrm{i})}{R_\mathrm{t} + R_\mathrm{i}}} \tag{2}$$

with $R_\mathrm{t}$ and $R_\mathrm{i}$ the body radii at the base of any atmosphere, and $f_M(b)$ is the fractional interacting mass of the pair (their Eq. B1), built from density-weighted spherical caps of common height $d = (R_\mathrm{t} + R_\mathrm{i})(1 - b)$:

$$f_M \;=\; \frac{\rho_\mathrm{t}\, V^\mathrm{cap}_\mathrm{t} + \rho_\mathrm{i}\, V^\mathrm{cap}_\mathrm{i}}{\rho_\mathrm{t}\, V_\mathrm{t} + \rho_\mathrm{i}\, V_\mathrm{i}}, \qquad V^\mathrm{cap}_\mathrm{t,i} = \frac{\pi}{3}\, d^2 \left(3 R_\mathrm{t,i} - d\right) \tag{3}$$

where $V_\mathrm{t,i}$ are the full body volumes. At equal bulk densities $f_M$ reduces to the fractional interacting volume of their Eq. B2. The common-height caps are a linearized bookkeeping: outside the fitted geometry, for a much denser and much smaller impactor near head-on, the raw $f_M$ can leave $[0, 1]$ and vary non-monotonically with $b$, so ZEPHYRUS clamps $f_M$ to $[0, 1]$. Within the fitted domain the clamp never engages.

Three input conventions follow the paper and must be honored by the caller: $v_\mathrm{c}$ is the speed at first contact, not the relative speed at infinity; the masses and radii exclude any atmosphere, with radii taken at the base of the atmosphere; and the densities are bulk values of the atmosphere-free bodies.

### Fitted domain and accuracy

The law is constrained by simulations spanning target masses of roughly 0.3 to 3 Earth masses, impactor masses down to about 0.05 Earth masses, bulk densities from about half to double Earth's, contact speeds of 1 to 3 $v_\mathrm{esc}$, all impact angles, and thin atmospheres of order 1% of the planet mass. The median deviation of the simulations from the law is 9%, rising to about 20% for slow, head-on impacts, whose outcomes are chaotic. The loss depends only mildly on the atmosphere mass in this thin-atmosphere regime, with a factor of 10 less atmosphere increasing the eroded fraction by roughly 10%; substantially thicker atmospheres, which can cushion the impactor, fall outside the law's domain. The function evaluates the law for any physically valid inputs and does not warn when masses, densities, or speeds leave the fitted ranges; staying inside them is the caller's responsibility, and the [limitations page](limitations.md) lists what the channel leaves out (impactor-side atmosphere, volatile delivery, and mantle stripping among them).

## The Roche et al. (2026) law

Roche et al. (2026) [^roche2026] generalize atmospheric erosion scaling to account explicitly for envelope mass fraction $f_\mathrm{atm} \equiv M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$. The scaling law splits atmospheric loss into near-field erosion ($X_\mathrm{NF}$, loss near the impact site) and far-field erosion ($X_\mathrm{FF}$, loss through ground motion):

$$X_\mathrm{atm} \;=\; X_\mathrm{NF} + X_\mathrm{FF} \tag{4}$$

where $X_\mathrm{atm}$ is clamped to $[0, 1]$. The fraction of the atmosphere located in the near-field region, $f_\mathrm{NF}$, is parameterized by a generalized logistic function:

$$f_\mathrm{NF} \;=\; \frac{\zeta_4}{\left(1 + \zeta_6 \exp\!\left(\zeta_3 \left(\frac{R_\mathrm{i}^\mathrm{r}}{R_\mathrm{t}^\mathrm{r}} - \zeta_2\right)\right)\right)^{\zeta_1}} + \zeta_5 \tag{5}$$

clamped to $[0, 1]$, and the far-field envelope fraction is $f_\mathrm{FF} = 1 - f_\mathrm{NF}$. Here $R_\mathrm{t}^\mathrm{r}$ and $R_\mathrm{i}^\mathrm{r}$ are the refractory (atmosphere-free) core-plus-mantle radii.

Near-field loss $X_\mathrm{NF}$ is described by:

$$X_\mathrm{NF} \;=\; f_\mathrm{NF} \left(\xi_1 - \xi_2\, (b + \xi_3)^2\right) \tag{6}$$

with a velocity floor at $v_\mathrm{c} / v_\mathrm{esc} = 1.0$, clamped to $[\max(0, X_\mathrm{NF}(v_\mathrm{c} = v_\mathrm{esc})), f_\mathrm{NF}]$.

Each coefficient vector ($\boldsymbol{\zeta}$, $\boldsymbol{\xi}$, $\boldsymbol{\psi}$) is evaluated as an empirical function combining power-law, polynomial, and logarithmic dependencies on impact parameter $b$, envelope mass fraction $f_\mathrm{atm}$, refractory target mass $M_\mathrm{t}^\mathrm{r} / M_\oplus$, and refractory impactor mass fraction $\gamma \equiv M_\mathrm{i}^\mathrm{r} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{r})$:

$$\zeta_i \;=\; q_{i1} + q_{i2}\, b + q_{i3}\, b^{q_{i4}} + q_{i5}\, f_\mathrm{atm} + q_{i6}\, f_\mathrm{atm}^2 + q_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{q_{i8}} \tag{7}$$

$$\xi_i \;=\; k_{i1} + k_{i2}\, \gamma + k_{i3}\, (\gamma + 0.05)^2 + k_{i4} \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^{k_{i5}} + k_{i6} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 1.0\right) + k_{i7} \log_{10}(f_\mathrm{atm}) \tag{8}$$

$$\psi_i \;=\; s_{i1} + s_{i2}\, (\gamma + 0.05)^{s_{i3}} + s_{i4} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 0.05\right)^{s_{i5}} + s_{i6} \log_{10}(f_\mathrm{atm}) \tag{9}$$

The far-field loss $X_\mathrm{FF}$ accounts for ground motion and scales with the modified specific impact energy $Q'_\mathrm{R}$:

$$X_\mathrm{FF} \;=\; f_\mathrm{FF} \operatorname{clip}\!\left(\psi_1 \exp\!\left(-\psi_2\, Q'_\mathrm{R} \left(1 + \frac{M_\mathrm{i}^\mathrm{r}}{M_\mathrm{t}^\mathrm{r}}\right) (1 - b)^{\psi_4}\right) + \psi_3, 0, 1\right) \tag{10}$$

with $Q'_\mathrm{R}$ expressed in $\mathrm{MJ\,kg^{-1}}$. Part of the far-field loss does not depend on impact energy whenever $\psi_1 + \psi_3 > 0$, driven by the $s_{16} \log_{10}(f_\mathrm{atm})$ term of $\psi_1$. At zero impact energy the empirical fit yields $X_\mathrm{FF} = f_\mathrm{FF}\operatorname{clip}(\psi_1 + \psi_3, 0, 1)$, which is strictly positive when $\psi_1 + \psi_3 > 0$. Within the fitted parameter box, $\psi_1 + \psi_3 > 0$ occurs above threshold values of $\gamma$ that depend on target mass:

| $f_\mathrm{atm}$ | $M_\mathrm{t}^\mathrm{r} = 0.35\,M_\oplus$ | $M_\mathrm{t}^\mathrm{r} = 1.0\,M_\oplus$ | $M_\mathrm{t}^\mathrm{r} = 2.0\,M_\oplus$ | $M_\mathrm{t}^\mathrm{r} = 5.0\,M_\oplus$ |
|---|---|---|---|---|
| 0.01 | $\gamma \ge 0.227$ | $\gamma \ge 0.375$ | $\gamma \ge 0.466$ | never |
| 0.03 | $\gamma \ge 0.340$ | $\gamma \ge 0.472$ | never | never |
| 0.05 | $\gamma \ge 0.386$ | never | never | never |
| 0.10 | $\gamma \ge 0.447$ | never | never | never |
| 0.20 | never | never | never | never |

Below $f_\mathrm{atm} = 0.01$ this positive zero-energy region expands. Consequently, the `'X_FF_zero_energy'` diagnostic flag can be raised for scenarios within the fitted range, as in two of the authors' reference impacts: one with $\gamma = 0.5$ at $f_\mathrm{atm} = 0.01$ on a $1\,M_\oplus$ target and one with $\gamma = 0.4$ at $f_\mathrm{atm} = 0.05$ on a $0.35\,M_\oplus$ target. ZEPHYRUS evaluates the authors' fit and reports this zero-energy far-field contribution in `diagnostics['X_FF_zero_energy']` and the `'X_FF_zero_energy'` flag. For the authors' own Moon-forming impact scenarios (Table D1 of Roche et al. 2026, evaluated with radii scaling as $M^{1/4}\,R_\oplus$), the law gives the following loss fractions at $f_\mathrm{atm} = 10^{-4}$: CA01 0.208, R12 0.310, CS12 0.467, C12 0.466, LS18a 0.613, and LS18b 0.363. For $f_\mathrm{atm}$ from $10^{-4}$ to $10^{-6}$ the law gives CA01 0.208 to 0.296 and LS18a 0.613 to 0.764; Roche et al. (2026, Sect. 4.3) estimate about 20% to 30% for the canonical impact (CA01) and about 70% to 80% for the most energetic synestia-forming scenario (Lock et al. 2018) for an atmosphere of 100 bar or less. The law reaches 0.70 for LS18a between $f_\mathrm{atm} = 10^{-5}$ and $10^{-6}$ (0.688, 0.764), yielding 0.613 at $10^{-4}$. For CA01 at $f_\mathrm{atm} = 10^{-4}$, the zero-energy far-field part is 0.057 of 0.208 (27.5% of the loss, 75% of the far-field part). Over the 790 fitting rows, the mean and maximum absolute misfits against simulation data are 0.0396 and 0.2224.

The forms on this page follow the authors' published code, which reproduces their Fig. 3. The printed Roche et al. (2026), Eq. 6 has $\zeta_3 R - \zeta_2$ in the exponent, the printed Roche et al. (2026), Eq. 10 has $M_\mathrm{t}^\mathrm{tot}$ in the mass ratio, and the $X_\mathrm{NF}$ floor at $v_\mathrm{c} = v_\mathrm{esc}$ is in the code and not in the paper; Table C prints some values with only 2 significant digits ($s_{15}, k_{15}, s_{25}, s_{33}$), and the constant term $q_{41}$ of $\zeta_4$ (161.057) and $\zeta_5$ (-161.022) almost cancel, so values from the published dataset are needed.

### Impact energy and escape speed definitions

The mutual escape speed at contact uses the total target mass and refractory impactor mass with refractory radii:

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{r})}{R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}}} \tag{11}$$

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
| $R_\mathrm{t}^\mathrm{r}, R_\mathrm{i}^\mathrm{r}$ | Refractory radii | m | Radius at base of atmosphere (mantle contact surface; explicitly distinct from inner core-mantle boundary contact) |
| $v_\mathrm{c}$ | Contact velocity | $\mathrm{m\,s^{-1}}$ | Speed at moment of first surface mantle contact |
| $v_\mathrm{esc}$ | Mutual escape speed | $\mathrm{m\,s^{-1}}$ | Calculated from total target mass and refractory impactor mass |
| $b$ | Dimensionless impact parameter | - | $\sin\beta \in [0, 1]$ |
| $\gamma$ | Impactor mass fraction | - | $M_\mathrm{i}^\mathrm{r} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{r})$ |
| $f_\mathrm{atm}$ | Target atmosphere mass fraction | - | $M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$ |
| $Q'_\mathrm{R}$ | Modified specific impact energy | $\mathrm{MJ\,kg^{-1}}$ | Interacting impact energy per unit total mass |
| $X_\mathrm{atm}$ | Atmospheric loss fraction | - | Fractional loss bounded in $[0, 1]$ |

### Coefficient origin

The scaling coefficients comprise 61 stored values (Roche et al. 2026, Appendix C):

- The 23 near-field mass coefficients ($q_{ij}$, $\zeta_5$, $\zeta_6$) are fitted to a separate suite of initialised SPH planets without impacts (radius ratios 0.1 to 1.0 in steps of 0.1 plus 0.001, $M_\mathrm{t}^\mathrm{r} \in [0.01, 5.0]\,M_\oplus$, $f_\mathrm{atm} \in [0.01, 0.2]$); $q_{14} = q_{34} = 2$ are fixed exponents. The impact dataset reaches a radius ratio of 1.015.
- The 38 loss coefficients ($k_{ij}$, $s_{ij}$) are fitted to 296 of the 300 new impact simulations (96 at $f_\mathrm{atm} = 0.01$, 100 at 0.1, and 100 at 0.2, all on a single target mass of about $1\,M_\oplus$) and the 494 simulations at $f_\mathrm{atm} = 0.05$ of Roche et al. (2025) [^roche2025] (spanning multiple target masses), giving 790 rows total. The loss coefficients ($k, s$) were fitted on $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \ge 0.498$ only; the 0.001 lower bound of the radius-ratio range comes from the separate $f_\mathrm{NF}$ planet-initialisation suite. The velocity exponent is shared across near-field and far-field regimes ($k_{15} = k_{25} = k_{35}$). Away from $1\,M_\oplus$, the impact data cover only $f_\mathrm{atm} = 0.05$.

All 61 values are transcribed at full precision from the published dataset [^roche2026]. In Table C of Roche et al. (2026), some coefficients are printed with only 2 significant digits ($s_{15}, k_{15}, s_{25}, s_{33}$), making the dataset values necessary for reproduction. Terms fixed in Tables C1 to C3 ($q_{14} = q_{34} = 2$, $q_{37} = q_{48} = 1$, and unlisted parameters fixed to 0) are held in internal constants.

### Calibrated domain and stability policy

The simulation suite constrains the law over the following parameter space:

- Target refractory mass: $M_\mathrm{t}^\mathrm{r} \in [0.35, 5.0]\,M_\oplus$
- Atmosphere mass fraction: $f_\mathrm{atm} \in [0.01, 0.2]$
- Impactor mass fraction: $\gamma \in [0.1, 0.5]$
- Impact parameter: $b \in [0.0, 0.9]$
- Radius ratio: $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \in [0.001, 1.015]$ (from the 0.001 lower bound of the initialised suite to the 1.015 maximum of the impact data)
- Contact velocity: $v_\mathrm{c} \in [1.0, 3.0]\,v_\mathrm{esc}$

To support planetary evolution and accretion calculations where conditions cross these empirical boundaries, `zephyrus.collision.impact_loss` applies a structured evaluation policy:

1. Diagnostic range flags. When inputs fall outside the calibrated range with a 1% relative tolerance ($f_\mathrm{atm} \notin [0.01, 0.2]$, $M_\mathrm{t}^\mathrm{r} \notin [0.35, 5.0]\,M_\oplus$, $\gamma \notin [0.1, 0.5]$, $b > 0.9$, $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} \notin [0.001, 1.015]$, or $v_\mathrm{c} / v_\mathrm{esc} > 3.0$), the loss fraction is computed and the out-of-range parameter names are recorded in `ImpactLossResult.flags` (`'f_atm'`, `'M_t_earth'`, `'gamma'`, `'b'`, `'R_ratio'`, `'v_ratio'`). Where a stability bound equals a fitted bound (such as $\gamma = 0.5$), the 1% relative tolerance does not apply: any clamped value is flagged. Additional diagnostic flags record physical regimes: `'v_sub_escape'` triggers when contact speed falls below 0.99 mutual escape speed ($v_\mathrm{c} < 0.99\,v_\mathrm{esc}$, where the near-field velocity floor evaluates at $v_\mathrm{esc}$), and `'X_FF_zero_energy'` triggers when the zero-energy far-field loss exceeds zero ($X_\mathrm{FF,zero} > 0$, which can occur inside the fitted range). In PROTEUS the atmosphere fraction is typically $10^{-5}$ to $10^{-3}$, 1 to 3 decades below the calibrated 0.01, so `'f_atm'` is flagged. Over a 400-impact grid in the fitted box ($M_\mathrm{t}^\mathrm{r} \in \{0.35, 1, 2, 5\}\,M_\oplus$, $\gamma \in \{0.1, 0.2, 0.3, 0.4, 0.5\}$, $b \in \{0, 0.2, 0.4, 0.6, 0.8\}$, $v_\mathrm{c} / v_\mathrm{esc} \in \{1, 1.5, 2, 3\}$, radii from equal bulk density), $X_\mathrm{NF}$ reaches $f_\mathrm{NF}$ in 172/400 impacts at $f_\mathrm{atm} = 10^{-2}$, 300/400 at $10^{-3}$, and 400/400 at $10^{-4}$, while $X$ keeps rising through the far-field term (for CA01, $X$ evaluates to 0.056 at $10^{-2}$, 0.131 at $10^{-3}$, 0.208 at $10^{-4}$, 0.252 at $10^{-5}$, and 0.296 at $10^{-6}$). An airless target ($f_\mathrm{atm} = 0$) returns no flags.
2. Stability clamping. The bounds match the parameter range where Roche et al. (2026, Sect. 4.1) verified numerical stability ($f_\mathrm{atm} \in [10^{-6}, 0.4]$, $M_\mathrm{t}^\mathrm{r} \in [10^{-3}, 10]\,M_\oplus$, $\gamma \in [10^{-3}, 0.5]$). Clamped parameter values enter only the empirical fit arguments ($\boldsymbol{\zeta}, \boldsymbol{\xi}, \boldsymbol{\psi}$); physical quantities ($v_\mathrm{esc}$, $v_\mathrm{c} / v_\mathrm{esc}$, $Q'_\mathrm{R}$) and the far-field mass ratio use the physical input masses and radii. `diagnostics['gamma']` reports the raw impactor mass fraction. Applied bounds are recorded in `diagnostics['clamped']`.
3. Zero atmosphere. When $f_\mathrm{atm} = 0.0$, the function returns $X = 0.0$ immediately.
4. Grazing collisions. The scaling law is calibrated for impact parameters $b \in [0.0, 0.9]$. At $b = 1$ the geometry factor is 0, so the bracket of Eq. (10) is $\psi_1 + \psi_3$ and $X_\mathrm{FF}$ evaluates to the zero-energy value $f_\mathrm{FF}\operatorname{clip}(\psi_1 + \psi_3, 0, 1)$ (0 in the benchmark case $M_\mathrm{t}^\mathrm{r} = 1.0\,M_\oplus$, $\gamma = 0.3$, $f_\mathrm{atm} = 0.01$, $v_\mathrm{c} / v_\mathrm{esc} = 1.5$, and 0.080 at $f_\mathrm{atm} = 10^{-4}$ with $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} = (\gamma / (1 - \gamma))^{1/3}$); $X$ at $b = 1$ is $X_\mathrm{NF} + X_\mathrm{FF,zero}$ (0.0315 for the benchmark case), while collisions with $b > 0.9$ trigger the `'b'` diagnostic flag.

### Physical caveats

The Roche et al. (2026, Sect. 4.1) scaling law does not account for pre-impact planetary rotation, surface liquid water oceans, thermal evolution, or a core mass fraction other than about 0.3 (Sect. 2.1). All simulations assume H2-He envelopes governed by the Hubbard & MacFarlane 1980 [^hubbard1980] equation of state.

Physical effects operate in both directions:

- **Miscibility (overestimate):** Magma-envelope miscibility at high pressures without a sharp boundary can cause the scaling law to overestimate the loss of massive envelopes on young planets (Roche et al. 2026, Sect. 4.1).
- **Post-impact thermal loss (underestimate):** The law covers the immediate shock- and vapour-plume-driven loss only; it neglects later thermally driven loss (an outflow driven by heat from the post-impact interior, Biersteker & Schlichting 2021), so for primordial H2-He envelopes the total loss can be higher (Roche et al. 2026, Sect. 4.1); that later loss becomes negligible for envelopes of higher mean molecular weight.
- **Atmospheric composition (upper limit):** Heavier atmospheres ($\mathrm{CO}, \mathrm{CO}_2$) are less susceptible to shock-driven removal, so shock-driven loss is an upper limit for a given atmosphere mass (Roche et al. 2026, Sect. 4.3). For PROTEUS atmospheres of higher mean molecular weight ($\mathrm{H}_2\mathrm{O}, \mathrm{CO}_2, \mathrm{O}_2$), the H2-He loss fractions likewise serve as upper limits.

## Partition across species

The returned fraction applies to the target's atmosphere as a whole. Consistent with the bulk-removal treatment of the continuous channel's unfractionated splits, the caller partitions the lost mass across atmospheric species without elemental fractionation.

---

[^kegerreis]: Kegerreis, J. A., Eke, V. R., Catling, D. C., Massey, R. J., Teodoro, L. F. A., & Zahnle, K. J. (2020). Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density. *The Astrophysical Journal Letters, 901*(2), L31. https://doi.org/10.3847/2041-8213/abb5fb

[^roche2026]: Roche, M. J., Lock, S. J., Carter, P. J., & Leinhardt, Z. M. (2026). Giant impacts preferentially remove low-mass atmospheres: a generalised scaling law for impact-driven atmospheric loss. *The Astrophysical Journal Letters* (accepted), arXiv:2610.06077. https://doi.org/10.48550/arXiv.2610.06077. Dataset: Zenodo, https://doi.org/10.5281/zenodo.23192423.

[^roche2025]: Roche, M. J., Lock, S. J., Dou, J., Carter, P. J., Kegerreis, J. A., & Leinhardt, Z. M. (2025). Atmospheric Loss during Giant Impacts: Mechanisms and Scaling of Near- and Far-field Loss. *The Planetary Science Journal, 6*, 149. https://doi.org/10.3847/PSJ/add929

[^leinhardt2012]: Leinhardt, Z. M., & Stewart, S. T. (2012). Collisions between gravity-dominated bodies. I. Outcome regimes and scaling laws. *The Astrophysical Journal, 745*(1), 79. https://doi.org/10.1088/0004-637X/745/1/79

[^hubbard1980]: Hubbard, W. B., & MacFarlane, J. J. (1980). Structure and evolution of Uranus and Neptune. *Journal of Geophysical Research: Solid Earth, 85*(B1), 225-234. https://doi.org/10.1029/JB085iB01p00225
