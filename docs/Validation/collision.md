# Validation: `src/zephyrus/collision.py`

This page tracks the tests that anchor the behaviour of `zephyrus.collision` against published giant-impact atmospheric erosion scaling laws.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_collision.py::test_scaling_law_pins_the_kegerreis_closed_form` | Kegerreis et al. (2020), ApJL 901, L31, Eqn. 1 (erosion scaling law, closed form) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the loss fraction for two identical Earth-like bodies head-on at 1.0 and 1.5 mutual escape speeds, where the law collapses to `X = 0.64 * (v_ratio^2 / sqrt(2))^0.65` with no geometry left. |
| `tests/test_collision.py::test_scaling_law_reproduces_kegerreis_table2_simulations` | Kegerreis et al. (2020), ApJL 901, L31, Tables 1 and 2 (SPH simulation suite) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the law against three simulated loss fractions from the paper's first suite (`b = 0.7`, `v_c = 3 v_esc`, impactor:target mass ratio `10^-0.5`) at 20% tolerance, the paper's stated simulation-to-law scatter. |
| `tests/test_collision.py::test_roche2026_oracle_reproduction` | Roche et al. (2026), arXiv:2610.06077, Eqns. 4-11, Table C1-C3, Zenodo doi:10.5281/zenodo.23192423 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins 12 reference simulation rows from the authors' published benchmark: Set A (1 Earth mass targets) to abs tolerance 1e-10 on $X_\mathrm{atm}$; Set B (0.35 to 4.98 Earth mass targets) to abs tolerance 1e-3. |
| `tests/test_collision.py::test_roche2026_specific_impact_energy_calculation` | Roche et al. (2025), PSJ 6, 149; Leinhardt & Stewart (2012), ApJ 745, 79 | [ADS 2025PSJ.....6..149R](https://ui.adsabs.harvard.edu/abs/2025PSJ.....6..149R) | Pins the modified specific impact energy $Q'_\mathrm{R}$ computed from SI masses, radii, contact velocity, and impact parameter against the 12 reference rows to rel tolerance 1e-12. |
| `tests/test_collision.py::test_roche2026_mutual_escape_speed_calculation` | Roche et al. (2026), arXiv:2610.06077, Eqn. 1 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins the mutual escape speed $v_\mathrm{esc}$ using total target-plus-impactor mass and refractory radii against contact velocity ratios to rel tolerance 1e-3. |
| `tests/test_collision.py::test_roche2026_velocity_floor_property` | Roche et al. (2026), arXiv:2610.06077, Section 3.2 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Verifies the physical near-field velocity floor where $v_\mathrm{c} / v_\mathrm{esc} = 0.5$ evaluates identically to $1.0$. Discrimination against uncorrected printed formulas (Eq. 6 and Eq. 10) is verified by mutant checks. |
| `tests/test_collision.py::test_roche2026_grazing_continuity_and_value` | Roche et al. (2026), arXiv:2610.06077, Eqn. 10 | [arXiv:2610.06077](https://arxiv.org/abs/2610.06077) | Pins the continuous mathematical evaluation of the Roche scaling law at grazing impact parameter $b = 1.0$ to the exact reference value $X \approx 0.031470$ without numerical instability. |

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
X_atm = f_NF * X_NF + (1 - f_NF) * X_FF
```

where `f_NF` is the near-field envelope mass fraction, `X_NF` is the near-field erosion fraction with a velocity floor at `v_c/v_esc = 1.0`, and `X_FF` is the far-field ground-shock erosion fraction scaling with `Q'_R / Q''_R`.

The 12 reference oracle rows are drawn from the authors' published SPH simulation dataset (Zenodo doi:10.5281/zenodo.23192423). In Set A (eight impacts onto 1 Earth mass targets with 1% to 20% atmospheres), the evaluated law matches the benchmark to 3.1e-13 in $X_\mathrm{FF}$ and $X_\mathrm{atm}$, and to 5.6e-17 in $f_\mathrm{NF}$ and $X_\mathrm{NF}$. In Set B (four impacts onto planets from 0.35 to 4.98 Earth masses), the maximum difference is 1.02e-4, below the 1e-3 benchmark tolerance.

The mutual escape speed $v_\mathrm{esc}$ in the Roche scaling framework adopts the total mass of target and impactor, $M_\mathrm{t}^\mathrm{tot} + M_\mathrm{i}^\mathrm{tot}$, divided by the sum of refractory radii $R_\mathrm{t}^\mathrm{r} + R_\mathrm{i}^\mathrm{r}$. Using only refractory masses shifts $v_\mathrm{esc}$ and fails the tolerance guard with relative error exceeding 1e-3 for substantial envelopes.

For grazing collisions where $b \ge 1.0$, the geometry factor $(1 - b)^{p_4}$ in $Q''_R$ is evaluated continuously so that $X_\mathrm{atm}$ approaches 0.031470 without numerical divergence or division by zero.

## Anchor type

Published benchmarks (SPH simulation suites and fitted scaling laws from the primary literature), plus closed-form analytical limits and continuity anchors.
