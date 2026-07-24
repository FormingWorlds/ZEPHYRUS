# Validation: `src/zephyrus/collision.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the behaviour of `zephyrus.collision` against the published giant-impact atmospheric erosion scaling law.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_collision.py::test_scaling_law_pins_the_kegerreis_closed_form` | Kegerreis et al. (2020), ApJL 901, L31, Eqn. 1 (erosion scaling law, closed form) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the loss fraction for two identical Earth-like bodies head-on at 1.0 and 1.5 mutual escape speeds, where the law collapses to `X = 0.64 * (v_ratio^2 / sqrt(2))^0.65` with no geometry left. |
| `tests/test_collision.py::test_scaling_law_reproduces_kegerreis_table2_simulations` | Kegerreis et al. (2020), ApJL 901, L31, Tables 1 and 2 (SPH simulation suite) | [ADS 2020ApJ...901L..31K](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..31K) | Pins the law against three simulated loss fractions from the paper's first suite (`b = 0.7`, `v_c = 3 v_esc`, impactor:target mass ratio `10^-0.5`) at 20% tolerance, the paper's stated simulation-to-law scatter. |

## Re-derivation note

`collision.mass_loss` returns the fractional atmospheric mass loss of the target,

```
X = min(0.64 * [ (v_c/v_esc)^2 * (M_i/M_tot)^(1/2) * (rho_i/rho_t)^(1/2) * f_M(b) ]^0.65, 1)
```

with the mutual escape speed `v_esc = sqrt(2 G (M_t + M_i) / (R_t + R_i))` and the fractional interacting mass `f_M` of the paper's Eqn. B1: density-weighted spherical caps of common height `d = (R_t + R_i)(1 - b)`, normalised by the density-weighted body volumes. At equal bulk densities `f_M` reduces exactly to the interacting volume `f_V` of Eqn. B2, which the test suite verifies as a property; a dedicated pin separates the two forms for an iron-rich impactor, where they differ by 0.0056 in `X` against a pin tolerance two orders of magnitude tighter.

The closed-form pin uses identical twin bodies so every bracket ratio except `M_i/M_tot = 1/2` is unity, giving `X = 0.64 * 0.5^0.325 = 0.510911` at contact speed equal to the mutual escape speed. The discrimination guards re-evaluate the law with the historical wrong mass-ratio denominator (`M_i/M_t`, giving 0.640000), a wrong outer exponent (0.5, giving 0.807261 at 1.5 `v_esc`), and a wrong velocity exponent (1, giving 0.664974); each sits far outside the `rel = 1e-4` pin.

The Table 2 pins reproduce the paper's own SPH results to within 4% in the fast grazing regime the authors report fits tightest; the 20% tolerance covers their stated scatter (9% median, about 20% for slow head-on impacts). The three scenarios share one mass ratio across a factor of 5.6 in total mass, so the pinned cluster also exercises the paper's finding that the loss is independent of the system mass at fixed impactor:target ratio.

## Anchor type

Published benchmark (SPH simulation suite and fitted scaling law from the same paper), plus the closed-form analytical limit of the fitted law.
