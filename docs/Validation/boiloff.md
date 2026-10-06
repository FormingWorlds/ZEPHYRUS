# Validation: `src/zephyrus/boiloff.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the boil-off branch of `zephyrus.boiloff`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_boiloff.py::test_parker_mach_sonic_limit_and_shutoff` | Owen & Wu (2016), ApJ 817, 107 (isothermal transonic Parker wind, Lambert-W form) | The photospheric Mach number is 1 identically with the launch level at the Bondi radius (the analytical sonic-point limit), falls monotonically as the level retreats inward, and has collapsed by more than six decades at their published shutoff R_p/R_B = 0.1. |
| `tests/test_boiloff.py::test_lambda_equals_two_bondi_radii_over_rp` | Fossati et al. (2017), A&A 598, A90 (restricted Jeans parameter); Owen & Wu (2016), ApJ 817, 107 | The identity Lambda = 2 R_B / R_p holds for every mean molecular mass, which is what makes the Owen & Wu shutoff equal to Lambda = 20 for every composition; the literature band 15 to 35 brackets the default threshold. |
| `tests/test_boiloff.py::test_bondi_cap_is_misener_eq10_at_the_wind_temperature` | Misener et al. (2025), ApJ 980, 152, Eq. 10; Gupta & Schlichting (2020), MNRAS 493, 792, Eq. 10 | The Bondi cap equals the closed form at the wind temperature T_eq / 2^(1/4) on a hydrostatic launch level, differs from the same form at T_eq by 2^(-3/8) exp(Lambda_launch (2^(1/4) - 1)), and stands to the Parker rate as tau_launch (T_wind / T_launch) exp(1/2 - Mach^2 / 2), so it binds on an optically thin launch level and not on a thick one. |
| `tests/test_boiloff.py::test_luminosity_cap_carries_the_tidal_barrier_reduction` | Erkaev et al. (2007), A&A 472, 329, Eq. 17 and Table 1 | The luminosity cap divides by the tidally reduced escape barrier, so it rises by the enhancement factor 1 / K, whose value at xi = 3 is the 1.92 printed in their Table 1; K = 1 reproduces the untidal cap identically. |

## Notes

The Bondi cap follows Misener et al. (2025, ApJ 980, 152, Eq. 10) at their recommended wind temperature, which carries a factor e^2 the Gupta & Schlichting (2020, MNRAS 493, 792, Eq. 10) original lacks, and the luminosity cap Gupta & Schlichting (2019, MNRAS 487, 24, Eq. 8), whose barrier carries the Erkaev tidal factor at xi = R_Hill / R_p so that an admitted residual and the energy-limited rate it competes against measure the same barrier; the cap ordering and the gate semantics are asserted as physics invariants in the same file. The Tang et al. (2024, ApJ 976, 221, Eq. 8) timescale comparison runs as a diagnostic only.

## Anchor type

Analytical limit plus published benchmark (the printed shutoff).

Date of last comparison against the sources: 2026-09-28.
