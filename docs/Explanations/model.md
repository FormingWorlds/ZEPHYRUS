# ZEPHYRUS model overview

ZEPHYRUS models two channels of atmospheric mass loss for rocky exoplanets coupled to the [PROTEUS](https://proteus-framework.org) interior-atmosphere framework: the continuous, bulk hydrodynamic escape driven by stellar XUV irradiation, and the impulsive erosion caused by giant impacts during accretion. The continuous channel implements an energy-limited (EL) formalism following Watson et al. (1981) [^watson] and Lopez & Fortney (2013) [^lopez]; it is called at each PROTEUS time step with the current planetary radius and mass, the stellar XUV flux supplied by [MORS](https://proteus-framework.org/MORS), and the escape radius computed from the atmospheric structure produced by AGNI or JANUS. The mass-loss rate it returns is distributed across atmospheric species according to their elemental mass mixing ratios, so the atmosphere is depleted in bulk without elemental fractionation. The impulsive channel implements the giant-impact erosion scaling laws of Kegerreis et al. (2020) [^kegerreis] and Roche et al. (2026) [^roche2026], which return the fraction of the target's atmosphere removed by a single collision.

A model parameter reference can be found [here](../Reference/parameters.md).

## Energy-limited escape

The mass-loss rate is computed by `escape.EL_escape` as

$$\dot{M}_\mathrm{EL} = \frac{\epsilon\,\pi\,R^3_\mathrm{XUV}\,F_\mathrm{XUV}}{G\,M_p\,K_\mathrm{tide}} \tag{1}$$

where $\epsilon$ is the escape efficiency factor (`epsilon`), $R_\mathrm{XUV}$ is the planetary radius at which the atmosphere becomes optically thick to stellar XUV photons (`Rxuv`), $F_\mathrm{XUV}$ is the XUV flux received at the planet (`Fxuv`) supplied by MORS, $M_p$ is the planetary mass (`Mp`), $G$ is the gravitational constant, and $K_\mathrm{tide}$ is the tidal correction factor described below. The efficiency $\epsilon$ quantifies the fraction of incident XUV energy that is converted into work against gravity to drive the outflow; canonical values for rocky planets lie in the range $0.1 \leq \epsilon \leq 0.3$, although ZEPHYRUS accepts any $\epsilon \in (0, 1]$.

### Radius scaling

The cubic radius factor in the numerator of Eq. (1) is selected at runtime by the `scaling` argument of `escape.EL_escape`:

| `scaling` | Expression | Description |
|---|---|---|
| `2` | $R_p\,R^2_\mathrm{XUV}$ | Default; XUV-absorbing cross-section weighted by surface radius |
| `3` | $R^3_\mathrm{XUV}$ | All three powers taken at the XUV radius |

Both forms reduce to $R_p^3$ when $R_\mathrm{XUV} = R_p$, which is the conservative lower bound on the mass-loss rate adopted by Luger & Barnes (2015) [^luger] and Moore et al. (2023) [^moore]. Allowing $R_\mathrm{XUV} > R_p$ increases the effective XUV-absorbing area and therefore the escape rate. In PROTEUS, $R_\mathrm{XUV}$ is recomputed at each time step from the atmospheric pressure–temperature profile at a user-specified reference pressure $P_\mathrm{XUV}$.

### Tidal correction $K_\mathrm{tide}$

When the `tidal_contribution` flag is `True`, the effective gravitational potential is reduced by the host star's tidal field following the tidal reduction factor of Erkaev et al. (2007), eq. 17 [^erkaev]:

$$K_\mathrm{tide} = 1 - \frac{3}{2\xi} + \frac{1}{2\xi^3}, \qquad \xi = \frac{R_\mathrm{Hill}}{R_\mathrm{XUV}} \tag{2}$$

with the Hill radius

$$R_\mathrm{Hill} = a\,(1-e)\,\left(\frac{M_p}{3\,M_\star}\right)^{1/3} \tag{3}$$

where $a$ is the planetary semi-major axis, $e$ is the orbital eccentricity, and $M_\star$ is the stellar mass. Factoring the numerator gives $K_\mathrm{tide} = (\xi - 1)^2\,(2\xi + 1) / (2\xi^3)$, which is non-negative for every $\xi > 0$ with a double root at $\xi = 1$. In the physical regime $\xi > 1$ it lies in $(0, 1)$, rising toward 1 for $\xi \gg 1$ (the XUV radius well inside the Hill sphere) and falling toward 0 as the atmosphere expands toward the Roche lobe at $\xi = 1$; because the escape rate divides by $K_\mathrm{tide}$, the rate is enhanced by the tidal correction and diverges as $\xi \to 1$. The tidally corrected rate is therefore defined only for $\xi > 1$: ZEPHYRUS raises a `ValueError` for $\xi \le 1$, where the atmosphere reaches the Roche lobe and the energy-limited approximation no longer applies. When `tidal_contribution` is `False`, $K_\mathrm{tide} = 1$ is enforced.

---

## Giant-impact atmospheric erosion

A giant impact removes part of the target planet's atmosphere in a single event. ZEPHYRUS computes the eroded fraction with `collision.mass_loss`, which implements the scaling law of Kegerreis et al. (2020), their Eq. 1 [^kegerreis]:

$$X \approx 0.64 \left[ \left(\frac{v_c}{v_\mathrm{esc}}\right)^2 \left(\frac{M_i}{M_\mathrm{tot}}\right)^{1/2} \left(\frac{\rho_i}{\rho_t}\right)^{1/2} f_M(b) \right]^{0.65} \tag{4}$$

capped at 1 for total erosion, where subscript $i$ denotes the impactor, $t$ the target, $M_\mathrm{tot} = M_i + M_t$, and $b \equiv \sin\beta$ is the dimensionless impact parameter for impact angle $\beta$ (0 head-on, 1 fully grazing). The prefactor and exponent are least-squares fits to the paper's suite of 259 SPH simulations, each with an uncertainty of 0.01. The mutual escape speed of the pair at contact is

$$v_\mathrm{esc} = \sqrt{\frac{2\,G\,(M_t + M_i)}{R_t + R_i}} \tag{5}$$

and $f_M(b)$ is the fractional interacting mass of the pair (their Eq. B1), built from density-weighted spherical caps of common height $d = (R_t + R_i)(1 - b)$:

$$f_M = \frac{\rho_t V^\mathrm{cap}_t + \rho_i V^\mathrm{cap}_i}{\rho_t V_t + \rho_i V_i}, \qquad V^\mathrm{cap}_{t,i} = \frac{\pi}{3} d^2 \left(3 R_{t,i} - d\right) \tag{6}$$

where $V_{t,i}$ are the full body volumes. At equal bulk densities $f_M$ reduces exactly to the fractional interacting volume of their Eq. B2. The common-height caps are a linearised bookkeeping: outside the fitted geometry, for a much denser and much smaller impactor near head-on, the raw $f_M$ can leave $[0, 1]$ and vary non-monotonically with $b$, so ZEPHYRUS clamps $f_M$ to $[0, 1]$. Within the fitted domain the clamp never engages.

Three input conventions follow the paper and must be honoured by the caller: $v_c$ is the speed at first contact, not the relative speed at infinity; the masses and radii exclude any atmosphere, with radii taken at the base of the atmosphere; and the densities are bulk values of the atmosphere-free bodies.

The returned fraction applies to the target's atmosphere as a whole. Consistent with the bulk-removal treatment of the continuous channel, the caller partitions the lost mass across atmospheric species without elemental fractionation.

### The Roche et al. (2026) law

Roche et al. (2026) [^roche2026] generalize atmospheric erosion scaling to account explicitly for envelope mass fraction $f_\mathrm{atm} \equiv M_\mathrm{atm} / M_\mathrm{t}^\mathrm{tot}$. The scaling law splits atmospheric loss into near-field erosion ($X_\mathrm{NF}$, loss near the impact site) and far-field erosion ($X_\mathrm{FF}$, loss through ground motion):

$$X_\mathrm{atm} \;=\; X_\mathrm{NF} + X_\mathrm{FF} \tag{7}$$

where $X_\mathrm{atm}$ is clamped to $[0, 1]$. The fraction of the atmosphere located in the near-field region, $f_\mathrm{NF}$, is parameterized by a generalized logistic function:

$$f_\mathrm{NF} \;=\; \frac{\zeta_4}{\left(1 + \zeta_6 \exp\!\left(\zeta_3 \left(\frac{R_\mathrm{i}^\mathrm{r}}{R_\mathrm{t}^\mathrm{r}} - \zeta_2\right)\right)\right)^{\zeta_1}} + \zeta_5 \tag{8}$$

clamped to $[0, 1]$, and the far-field envelope fraction is $f_\mathrm{FF} = 1 - f_\mathrm{NF}$. Here $R_\mathrm{t}^\mathrm{r}$ and $R_\mathrm{i}^\mathrm{r}$ are the refractory (atmosphere-free) core-plus-mantle radii.

Near-field loss $X_\mathrm{NF}$ is described by:

$$X_\mathrm{NF} \;=\; f_\mathrm{NF} \left(\xi_1 - \xi_2\, (b + \xi_3)^2\right) \tag{9}$$

with a velocity floor at $v_\mathrm{c} / v_\mathrm{esc} = 1.0$, clamped to $[\max(0, X_\mathrm{NF}(v_\mathrm{c} = v_\mathrm{esc})), f_\mathrm{NF}]$.

Each coefficient vector ($\boldsymbol{\zeta}$, $\boldsymbol{\xi}$, $\boldsymbol{\psi}$) is evaluated as an empirical function combining power-law, polynomial, and logarithmic dependencies on impact parameter $b$, envelope mass fraction $f_\mathrm{atm}$, refractory target mass $M_\mathrm{t}^\mathrm{r} / M_\oplus$, and refractory impactor mass fraction $\gamma \equiv M_\mathrm{i}^\mathrm{r} / (M_\mathrm{i}^\mathrm{r} + M_\mathrm{t}^\mathrm{r})$:

$$\zeta_i \;=\; q_{i1} + q_{i2}\, b + q_{i3}\, b^{q_{i4}} + q_{i5}\, f_\mathrm{atm} + q_{i6}\, f_\mathrm{atm}^2 + q_{i7}\left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus}\right)^{q_{i8}} \tag{10}$$

$$\xi_i \;=\; k_{i1} + k_{i2}\, \gamma + k_{i3}\, (\gamma + 0.05)^2 + k_{i4} \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^{k_{i5}} + k_{i6} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 1.0\right) + k_{i7} \log_{10}(f_\mathrm{atm}) \tag{11}$$

$$\psi_i \;=\; s_{i1} + s_{i2}\, (\gamma + 0.05)^{s_{i3}} + s_{i4} \left(\frac{M_\mathrm{t}^\mathrm{r}}{M_\oplus} + 0.05\right)^{s_{i5}} + s_{i6} \log_{10}(f_\mathrm{atm}) \tag{12}$$

The far-field loss $X_\mathrm{FF}$ accounts for ground motion and scales with the modified specific impact energy $Q'_\mathrm{R}$:

$$X_\mathrm{FF} \;=\; f_\mathrm{FF} \operatorname{clip}\!\left(\psi_1 \exp\!\left(-\psi_2\, Q'_\mathrm{R} \left(1 + \frac{M_\mathrm{i}^\mathrm{r}}{M_\mathrm{t}^\mathrm{r}}\right) (1 - b)^{\psi_4}\right) + \psi_3, 0, 1\right) \tag{13}$$

with $Q'_\mathrm{R}$ expressed in $\mathrm{MJ\,kg^{-1}}$. Part of the far-field loss does not depend on impact energy whenever $\psi_1 + \psi_3 > 0$, driven by the $s_{16} \log_{10}(f_\mathrm{atm})$ term of $\psi_1$. At zero impact energy the empirical fit yields $X_\mathrm{FF} = f_\mathrm{FF}\operatorname{clip}(\psi_1 + \psi_3, 0, 1)$, which is strictly positive when $\psi_1 + \psi_3 > 0$. In the fitted parameter box, $\psi_1 + \psi_3 > 0$ in 3.5% of cases, occurring at $f_\mathrm{atm} = 0.01$ for $\gamma \gtrsim 0.23$, at $0.03$ for $\gamma \gtrsim 0.34$, at $0.05$ for $\gamma \gtrsim 0.39$, and at $0.1$ for $\gamma \gtrsim 0.45$ (never at $0.2$). Below $f_\mathrm{atm} = 0.01$ this positive zero-energy region expands. Consequently, the `'X_FF_zero_energy'` diagnostic flag can be raised for scenarios within the fitted range, as in Set A row 1 and Set B row 0. ZEPHYRUS evaluates the authors' fit and reports this zero-energy far-field contribution in `diagnostics['X_FF_zero_energy']` and the `'X_FF_zero_energy'` flag. For the authors' own Moon-forming impact scenarios (Table D1 of Roche et al. 2026, evaluated with radii scaling as $M^{1/4}\,R_\oplus$), the law gives the following loss fractions at $f_\mathrm{atm} = 10^{-4}$: CA01 0.208, R12 0.310, CS12 0.467, C12 0.466, LS18a 0.613, and LS18b 0.363. For $f_\mathrm{atm}$ from $10^{-4}$ to $10^{-6}$ the law gives CA01 0.208 to 0.296 and LS18a 0.613 to 0.764; Roche et al. (2026, Sect. 4.3) estimate about 20% to 30% for the canonical impact (CA01) and about 70% to 80% for the most energetic synestia-forming scenario (Lock et al. 2018) for an atmosphere of 100 bar or less. The law reaches the 70% to 80% range for LS18a only at $f_\mathrm{atm} \approx 10^{-5}$ to $10^{-6}$ (0.689 and 0.764), yielding 0.613 at $10^{-4}$. For CA01 at $f_\mathrm{atm} = 10^{-4}$, the zero-energy far-field component is 0.0572 (0.057 of 0.208, representing 27.5% of total loss and 75.0% of the far-field component). Over the 790 fitting rows, the mean and maximum absolute misfits against simulation data are 0.0396 and 0.2224.

The forms on this page follow the authors' published code, which reproduces their Fig. 3. The printed Eq. 6 has $\zeta_3 R - \zeta_2$ in the exponent, the printed Eq. 10 has $M_\mathrm{t}^\mathrm{tot}$ in the mass ratio, and the $X_\mathrm{NF}$ floor at $v_\mathrm{c} = v_\mathrm{esc}$ is in the code and not in the paper; the printed Table C values have 4 significant digits, which is not enough, because the constant term $q_{41}$ of $\zeta_4$ (161.057) and $\zeta_5$ (-161.022) almost cancel.

### Impact energy and escape speed definitions

The mutual escape speed at contact uses the total target mass and refractory impactor mass with refractory radii:

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{r})}{R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}}} \tag{14}$$

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

The scaling coefficients comprise 61 numbers (Roche et al. 2026, Appendix C):

- The 23 near-field mass coefficients ($q_{ij}$, $\zeta_5$, $\zeta_6$) are fitted to a separate suite of initialised SPH planets without impacts (radius ratios 0.1 to 1.0 in steps of 0.1 plus 0.001, $M_\mathrm{t}^\mathrm{r} \in [0.01, 5.0]\,M_\oplus$, $f_\mathrm{atm} \in [0.01, 0.2]$). The impact dataset reaches a radius ratio of 1.015.
- The 38 loss coefficients ($k_{ij}$, $s_{ij}$) are fitted to 300 new impact simulations (with a single target mass of about 1 $M_\oplus$, and $f_\mathrm{atm} \in \{0.01, 0.1, 0.2\}$) and the $f_\mathrm{atm} = 0.05$ simulations of Roche et al. (2025) [^roche2025], which span multiple target masses. Away from 1 $M_\oplus$, the impact data cover only $f_\mathrm{atm} = 0.05$.

All 61 values are transcribed at full precision from the published dataset [^roche2026]. Terms fixed in Table C1-C3 ($q_{37} = q_{48} = 1$ and other unlisted parameters fixed to 0) are held in internal constants.

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
4. Grazing collisions. The scaling law is calibrated for impact parameters $b \in [0.0, 0.9]$. At $b = 1$ the geometry factor is 0, so the bracket of Eq. (13) is $\psi_1 + \psi_3$ and $X_\mathrm{FF}$ evaluates to the zero-energy value $f_\mathrm{FF}\operatorname{clip}(\psi_1 + \psi_3, 0, 1)$ (0 in the benchmark case $M_\mathrm{t}^\mathrm{r} = 1.0\,M_\oplus$, $\gamma = 0.3$, $f_\mathrm{atm} = 0.01$, $v_\mathrm{c} / v_\mathrm{esc} = 1.5$, and 0.080 at $f_\mathrm{atm} = 10^{-4}$ with $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} = (\gamma / (1 - \gamma))^{1/3}$); $X$ at $b = 1$ is $X_\mathrm{NF} + X_\mathrm{FF,zero}$ (0.0315 for the benchmark case), while collisions with $b > 0.9$ trigger the `'b'` diagnostic flag.

### Physical caveats

The Roche et al. (2026, Sect. 4.1) scaling law does not account for pre-impact planetary rotation, surface liquid water oceans, thermal evolution, miscible envelopes, or a core mass fraction other than about 0.3 (Sect. 2.1). All simulations assume H2-He envelopes governed by the Hubbard & MacFarlane 1980 [^hubbard1980] equation of state. Heavier (CO or CO2) atmospheres are less easily removed (Roche et al. 2026, Sect. 4.3); we expect the same for the O2-rich atmospheres of PROTEUS, so H2-He-based loss fractions serve as upper limits. Hydrodynamic simulations track only the first tens of hours after the event (Roche et al. 2026, Sect. 4.1).

### Law selection with `impact_loss`

Callers select the desired erosion law with `zephyrus.collision.impact_loss(law=...)`, where `law` is `'kegerreis2020'` or `'roche2026'`. The function returns an `ImpactLossResult` dataclass with `fraction`, `flags`, and `diagnostics`.

---

## Coupling to PROTEUS

ZEPHYRUS treats atmospheric escape as a bulk process: at each PROTEUS time step the total mass-loss rate from Eq. (1) is partitioned across atmospheric species in proportion to their elemental mass mixing ratios as computed by CALLIOPE. No elemental fractionation between light and heavy species is imposed in the outflow itself; however, because only outgassed volatiles are subject to escape while dissolved species remain in the magma ocean reservoir, escape fractionates the planet's *total* (interior + atmosphere) volatile budget over time, preferentially retaining species that are highly soluble in silicate melts (e.g. H$_2$O, S$_2$). 

More about this in its dedicated [page](proteus.md).


## Regime of validity

The EL formalism is appropriate in the high-irradiation, hydrodynamic regime that dominates atmospheric loss during the first $\sim 10^6$–$10^8$ yr of evolution for close-in rocky planets [^watson][^lammer2003]. Outside this regime—at lower XUV fluxes or for less extended atmospheres—non-thermal escape (Jeans escape, ion pickup, charge exchange) becomes comparable to or exceeds the hydrodynamic rate, and the bulk EL prescription no longer applies. ZEPHYRUS does not currently include these processes; users should verify that the integrated XUV-driven loss exceeds non-thermal estimates (e.g. $\sim 10^7$–$10^8$ g s$^{-1}$ for an Earth-mass planet; Kislyakova et al. 2014 [^kislyakova]) before interpreting model outputs.

Similarly, the bulk-removal assumption breaks down when the hydrodynamic particle flux drops below the critical flux required to drag heavy species against gravity, at which point compositional fractionation in the outflow becomes significant [^wordsworth2018][^cherubim2024]. Following Yoshida et al. (2022) [^yoshida], the critical flux for H$_2$O in an H$_2$ background is $\approx 1.9 \times 10^{8}$ g s$^{-1}$.

The giant-impact erosion law (Eq. 4) is constrained by simulations spanning target masses of roughly 0.3 to 3 $M_\oplus$, impactor masses down to about 0.05 $M_\oplus$, bulk densities from about half to double Earth's, contact speeds of 1 to 3 $v_\mathrm{esc}$, all impact angles, and thin atmospheres of order 1 percent of the planet mass. The median deviation of the simulations from the law is 9 percent, rising to about 20 percent for slow, head-on impacts, whose outcomes are chaotic. The loss depends only mildly on the atmosphere mass in this thin-atmosphere regime, with a factor of 10 less atmosphere increasing the eroded fraction by roughly 10 percent; substantially thicker atmospheres, which can cushion the impactor, fall outside the law's regime.

---

[^watson]: Watson, A. J., Donahue, T. M., & Walker, J. C. G. (1981). The dynamics of a rapidly escaping atmosphere: applications to the evolution of Earth and Venus. *Icarus, 48*(2), 150–166. https://doi.org/10.1016/0019-1035(81)90101-9


[^lopez]: Lopez, E. D., & Fortney, J. J. (2013). The role of core mass in controlling evaporation: the Kepler radius distribution and the Kepler-36 density dichotomy. *The Astrophysical Journal, 776*(1), 2. https://doi.org/10.1088/0004-637X/776/1/2

[^erkaev]: Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007). Roche lobe effects on the atmospheric loss from "Hot Jupiters". *Astronomy & Astrophysics, 472*(1), 329–334. https://doi.org/10.1051/0004-6361:20066929

[^luger]: Luger, R., & Barnes, R. (2015). Extreme water loss and abiotic O$_2$ buildup on planets throughout the habitable zones of M dwarfs. *Astrobiology, 15*(2), 119–143. https://doi.org/10.1089/ast.2014.1231

[^moore]: Moore, K., Cowan, N. B., & Boukaré, C.-É. (2023). The role of magma oceans in maintaining surface water on rocky planets orbiting M-dwarfs. *Monthly Notices of the Royal Astronomical Society, 526*(4), 6235–6249. https://doi.org/10.1093/mnras/stad3138

[^lammer2003]: Lammer, H., Selsis, F., Ribas, I., et al. (2003). Atmospheric loss of exoplanets resulting from stellar X-ray and extreme-ultraviolet heating. *The Astrophysical Journal, 598*(2), L121–L124. https://doi.org/10.1086/380815

[^kislyakova]: Kislyakova, K. G., Johnstone, C. P., Odert, P., et al. (2014). Stellar wind interaction and pick-up ion escape of the Kepler-11 "super-Earths". *Astronomy & Astrophysics, 562*, A116. https://doi.org/10.1051/0004-6361/201322933

[^wordsworth2018]: Wordsworth, R. D., Schaefer, L. K., & Fischer, R. A. (2018). Redox evolution via gravitational differentiation on low-mass planets: implications for abiotic oxygen, water loss, and habitability. *The Astronomical Journal, 155*(5), 195. https://doi.org/10.3847/1538-3881/aab608

[^cherubim2024]: Cherubim, C., Wordsworth, R., Hu, R., & Shkolnik, E. (2024). Strong Fractionation of Deuterium and Helium in Sub-Neptune Atmospheres along the Radius Valley. *The Astrophysical Journal, 967*(2), 139. https://doi.org/10.3847/1538-4357/ad3e77

[^yoshida]: Yoshida, T., Terada, N., Ikoma, M., & Kuramoto, K. (2022). Less Effective Hydrodynamic Escape of H$_2$–H$_2$O Atmospheres on Terrestrial Planets Orbiting Pre-main-sequence M Dwarfs. *The Astrophysical Journal, 934*(2), 137. https://doi.org/10.3847/1538-4357/ac7be7

[^kegerreis]: Kegerreis, J. A., Eke, V. R., Catling, D. C., Massey, R. J., Teodoro, L. F. A., & Zahnle, K. J. (2020). Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density. *The Astrophysical Journal Letters, 901*(2), L31. https://doi.org/10.3847/2041-8213/abb5fb

[^roche2026]: Roche, M. J., Lock, S. J., Carter, P. J., & Leinhardt, Z. M. (2026). Giant impacts preferentially remove low-mass atmospheres: a generalised scaling law for impact-driven atmospheric loss. *The Astrophysical Journal Letters* (accepted), arXiv:2610.06077. https://doi.org/10.48550/arXiv.2610.06077. Dataset: Zenodo, https://doi.org/10.5281/zenodo.23192423.

[^roche2025]: Roche, M. J., Lock, S. J., Dou, J., Carter, P. J., Kegerreis, J. A., & Leinhardt, Z. M. (2025). Atmospheric Loss during Giant Impacts: Mechanisms and Scaling of Near- and Far-field Loss. *The Planetary Science Journal, 6*, 149. https://doi.org/10.3847/PSJ/add929

[^leinhardt2012]: Leinhardt, Z. M., & Stewart, S. T. (2012). Collisions between gravity-dominated bodies. I. Outcome regimes and scaling laws. *The Astrophysical Journal, 745*(1), 79. https://doi.org/10.1088/0004-637X/745/1/79

[^hubbard1980]: Hubbard, W. B., & MacFarlane, J. J. (1980). Structure and evolution of Uranus and Neptune. *Journal of Geophysical Research: Solid Earth, 85*(B1), 225-234. https://doi.org/10.1029/JB085iB01p00225