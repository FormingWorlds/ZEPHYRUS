# Validation: `src/zephyrus/boiloff.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the behaviour of `zephyrus.boiloff` against the published criteria for the boil-off regime.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_boiloff.py::test_published_threshold_equivalence_holds_for_atomic_hydrogen` | Fossati et al. (2017), A&A 598:A90, Eq. 10 (restricted Jeans parameter); Affolter et al. (2023), A&A 676:A119, Eq. 7 and the accompanying statement that `Lambda = 20` matches the Owen & Wu (2016) radius condition | [ADS 2017A&A...598A..90F](https://ui.adsabs.harvard.edu/abs/2017A%26A...598A..90F), [ADS 2023A&A...676A.119A](https://ui.adsabs.harvard.edu/abs/2023A%26A...676A.119A) | Pins the published equality of the two independent boil-off thresholds, and the mean-molecular-weight convention it depends on. |

## Re-derivation note

Two thresholds for the same regime appear in the literature. Owen & Wu (2016), ApJ 817:107, place the onset of boil-off where the planet radius reaches a tenth of the Bondi radius, `Rp = 0.1 R_B`. Fossati et al. (2017) instead use the restricted Jeans parameter `Lambda = G Mp m_H / (kb Teq Rp)`, evaluated for atomic hydrogen, and place the regime below `Lambda` of 15 to 35. Affolter et al. (2023) state that a `Lambda` of 20 is identical to the Owen & Wu radius condition.

Combining the two definitions gives

```
Rp / R_B = 2 / (mu * Lambda)
```

so the published equality holds exactly at `mu = 1`, and `Lambda = 20` then corresponds to `Rp / R_B = 0.1`. The identity is therefore specific to evaluating the Bondi radius for atomic hydrogen. Evaluated instead for a hydrogen-helium envelope (`mu = 2.35`, the default of `bondi_radius`) the same radius sits at `0.1 / 2.35 = 0.0426` of the Bondi radius, a factor 2.35 away.

The test pins both statements: the equality at `mu = 1`, and the offset at `mu = 2.35`, with a guard asserting the two differ by more than 0.05. This is why the surviving-mass calculation in `boiloff_mass_factor` is expressed in the radius ratio `alpha = Rp / R_B` rather than in `Lambda`: the radius form is the one Owen & Wu (2016), Eq. 16, is derived in, so it carries no conversion and cannot silently inherit the wrong convention.

Scale: `Lambda` is dimensionless and of order unity to a few hundred; the pinned case sits at exactly 20 by construction of the test mass. A unit slip in the Boltzmann constant or the hydrogen mass moves it by many orders of magnitude and fails the `rel=1e-12` pin.

## Anchor type

Published criterion and its algebraic equivalence, cross-checked between two independent formulations (Fossati et al. 2017 / Affolter et al. 2023 versus Owen & Wu 2016).

## Cross-references

- `src/zephyrus/boiloff.py`, `restricted_jeans` and `bondi_radius` docstring Notes sections.
- `src/zephyrus/escape.py`, `EL_escape` `K_tide_floor` parameter: the clamp that keeps the energy-limited rate finite as an atmosphere approaches its Roche lobe, which is the geometry this regime test identifies.
