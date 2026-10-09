# Validation: `src/zephyrus/collision.py`

This page tracks the tests that anchor the behaviour of `zephyrus.collision` against published giant-impact atmospheric erosion scaling laws.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_collision.py::test_scaling_law_pins_the_kegerreis_closed_form` | Kegerreis et al. (2020), ApJL 901, L31, Eq. 1 (erosion scaling law, closed form) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the loss fraction for two identical Earth-like bodies head-on at 1.0 and 1.5 mutual escape speeds, where the law collapses to `X = 0.64 * (v_ratio^2 / sqrt(2))^0.65` with no geometry left. |
| `tests/test_collision.py::test_scaling_law_reproduces_kegerreis_table2_simulations` | Kegerreis et al. (2020), ApJL 901, L31, Tables 1 and 2 (SPH simulation suite) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the law against three simulated loss fractions from the paper's first suite (`b = 0.7`, `v_c = 3 v_esc`, impactor:target mass ratio `10^-0.5`) at 20% tolerance, the paper's stated simulation-to-law scatter. |
| `tests/test_collision.py::test_roche2026_oracle_reproduction` | Roche et al. (2026), arXiv:2610.06077, Eqs. 4-11, Table C1-C3, Zenodo doi:10.5281/zenodo.23192423 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins 12 reference rows of the scaling law evaluated by the authors (Zenodo doi:10.5281/zenodo.23192423, `scaling_law.csv`, column `X_atm_calc`): eight reference rows on $M_\mathrm{t} \approx 1\,M_\oplus$ targets ($f_\mathrm{atm} \in [0.01, 0.2]$) to abs tolerance 1e-10 on $X_\mathrm{atm}$; four reference rows on $M_\mathrm{t} \in [0.35, 4.98]\,M_\oplus$ targets to abs tolerance 2e-4. |
| `tests/test_collision.py::test_roche2026_oracle_set_a_impact_loss` | Roche et al. (2026), Zenodo doi:10.5281/zenodo.23192423 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins 8 reference rows on $M_\mathrm{t} \approx 1\,M_\oplus$ targets ($f_\mathrm{atm} \in [0.01, 0.2]$) evaluated through public `impact_loss` against published `X_atm_calc` within $2 \times 10^{-4}$ tolerance. |
| `tests/test_collision.py::test_roche2026_oracle_set_b_impact_loss` | Roche et al. (2026), Zenodo doi:10.5281/zenodo.23192423 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins 4 reference rows on $M_\mathrm{t} \in [0.35, 4.98]\,M_\oplus$ targets evaluated through public `impact_loss` against published `X_atm_calc` within $2 \times 10^{-4}$ tolerance. |
| `tests/test_collision.py::test_roche2026_specific_impact_energy_calculation` | Roche et al. (2025), PSJ 6, 149, Eqs. 13-18; Leinhardt & Stewart (2012), ApJ 745, 79 | [ADS 2025PSJ.....6..149R](https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..149R) | Pins the modified specific impact energy $Q'_\mathrm{R}$ computed from SI masses, radii, contact velocity, and impact parameter against the eight reference rows on $M_\mathrm{t} \approx 1\,M_\oplus$ targets to rel tolerance 1e-12 (matches to $7 \times 10^{-16}$). |
| `tests/test_collision.py::test_roche2026_mutual_escape_speed_calculation` | Roche et al. (2026), Eq. 1 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins the mutual escape speed $v_\mathrm{esc}$ using total target mass and refractory impactor mass ($M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{r}$) and refractory radii against contact velocity ratios for the eight reference rows on $M_\mathrm{t} \approx 1\,M_\oplus$ targets to rel tolerance 1e-3. |

## Self-consistency checks

| Test id | Reference property | Scope |
|---|---|---|
| `tests/test_collision.py::test_roche2026_velocity_floor_property` | Roche et al. (2026), Zenodo doi:10.5281/zenodo.23192423 (authors' code) | Verifies the near-field velocity floor where $v_\mathrm{c} / v_\mathrm{esc} = 0.5$ evaluates identically to $1.0$. |
| `tests/test_collision.py::test_roche2026_grazing_continuity_and_value` | Roche et al. (2026), Eq. 10 | At $b = 1$ the geometry factor is set to 0 to avoid $0 \times \infty$, and $X_\mathrm{FF}$ evaluates to the zero-energy value (0 in this benchmark case), giving continuous evaluation ($X \approx 0.031470$) without numerical instability. |
| `tests/test_collision.py::test_roche2026_table_d1_moon_forming_scenarios` | Roche et al. (2026), Table D1 | Regression pin for atmospheric erosion fractions across six Moon-forming scenarios computed by this code for the published scenario parameters. |

## Re-derivation notes

### The Kegerreis et al. (2020) law

`collision.mass_loss` returns the fractional atmospheric mass loss of the target,

```
X = min(0.64 * [ (v_c/v_esc)^2 * (M_i/M_tot)^(1/2) * (rho_i/rho_t)^(1/2) * f_M(b) ]^0.65, 1)
```

with the mutual escape speed `v_esc = sqrt(2 G (M_t + M_i) / (R_t + R_i))` and the fractional interacting mass `f_M` of the paper's Eqn. B1: density-weighted spherical caps of common height `d = (R_t + R_i)(1 - b)`, normalised by the density-weighted body volumes. At equal bulk densities `f_M` reduces exactly to the interacting volume `f_V` of Eqn. B2, which the test suite verifies as a property; a dedicated pin separates the two forms for an iron-rich impactor, where they differ by 0.0056 in `X` against a pin tolerance two orders of magnitude tighter.

The closed-form pin uses identical twin bodies so every bracket ratio except `M_i/M_tot = 1/2` is unity, giving `X = 0.64 * 0.5^0.325 = 0.510911` at contact speed equal to the mutual escape speed. The discrimination guards re-evaluate the law with the historical wrong mass-ratio denominator (`M_i/M_t`, giving 0.640000), a wrong outer exponent (0.5, giving 0.807261 at 1.5 `v_esc`), and a wrong velocity exponent (1, giving 0.664974); each sits far outside the `rel = 1e-4` pin.

The Table 2 pins reproduce the paper's own SPH results to 4% in the fast grazing regime the authors report fits tightest; the 20% tolerance covers their stated scatter (9% median, about 20% for slow head-on impacts). The three scenarios share one mass ratio across a factor of 5.6 in total mass, so the pinned cluster also exercises the paper's finding that the loss is independent of the system mass at fixed impactor:target ratio.

### The Roche et al. (2026) law

`collision.mass_loss_roche2026` and `collision.impact_loss` return the fractional atmospheric mass loss from:

```
X_atm = X_NF + X_FF
```

where `f_NF` is the near-field envelope mass fraction, `X_NF` is the near-field erosion fraction with a velocity floor at `v_c/v_esc = 1.0`, and `X_FF` is the far-field ground-shock erosion fraction scaling with modified specific impact energy $Q'_\mathrm{R}$ (in $\mathrm{MJ\,kg^{-1}}$).

Part of the far-field loss does not depend on impact energy when $\psi_1 + \psi_3 > 0$, reported in `diagnostics['X_FF_zero_energy']` and the `'X_FF_zero_energy'` flag; the zero-energy behaviour and Table D1 benchmark values are detailed in the [model overview](../Explanations/model.md#the-roche-et-al-2026-law).

The 12 reference oracle rows are drawn from the authors' published dataset (Zenodo doi:10.5281/zenodo.23192423, `scaling_law.csv`, column `X_atm_calc`): eight reference impacts onto $M_\mathrm{t} \approx 1\,M_\oplus$ targets with $f_\mathrm{atm} \in [0.01, 0.2]$ and four reference impacts onto planets with $M_\mathrm{t} \in [0.35, 4.98]\,M_\oplus$. Evaluated through `impact_loss`, all 12 rows match the published benchmark within $2 \times 10^{-4}$ tolerance. The residuals arise from $v_\mathrm{c}/v_\mathrm{esc}$ (the published contact velocity has two-digit rounding) and nominal $\gamma$ or mass ratio, whereas modified specific impact energy $Q'_\mathrm{R}$ matches the dataset value to $7 \times 10^{-16}$. Diagnostic range flags are checked with a 1% relative tolerance (`_ROCHE2026_RANGE_RTOL = 0.01`) on parameter bounds. Speeds below 0.99 mutual escape speed trigger the `'v_sub_escape'` diagnostic flag.

The mutual escape speed $v_\mathrm{esc}$ in the Roche scaling framework adopts the total target mass and refractory impactor mass, $M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{r}$, divided by the sum of refractory radii $R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}$. Using only refractory masses shifts $v_\mathrm{esc}$ and fails the $10^{-3}$ tolerance guard for every reference row on $M_\mathrm{t} \approx 1\,M_\oplus$ targets ($f_\mathrm{atm} \in [0.01, 0.2]$).

In the grazing test case ($M_\mathrm{t}^\mathrm{r} = 1\,M_\oplus$, $\gamma = 0.3$, $R_\mathrm{i}^\mathrm{r} / R_\mathrm{t}^\mathrm{r} = (\gamma / (1 - \gamma))^{1/3}$, $f_\mathrm{atm} = 0.01$, $v_\mathrm{c}/v_\mathrm{esc} = 1.5$), $X$ evaluates continuously to $X \approx 0.031470$ as $b \to 1$. Over the stability box ($\gamma \le 0.5$), $\psi_4$ spans $-0.53$ to $1.32$, and $X_\mathrm{FF}$ at $b = 1$ equals the zero-energy value ($X_\mathrm{FF,zero} = 0$ in this benchmark case, $0.080$ at $f_\mathrm{atm} = 10^{-4}$). The factor is set to 0 at $b = 1$ to avoid $0 \times \infty$; the evaluation stays continuous as $b \to 1$.

## Anchor type

Published benchmarks (SPH simulation suites and fitted scaling laws from the primary literature), plus closed-form analytical limits and continuity anchors.
