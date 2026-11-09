# Validation: `src/zephyrus/boiloff.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the boil-off branch of `zephyrus.boiloff`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_boiloff.py::test_parker_mach_sonic_limit_and_shutoff` | Owen & Wu (2016)[^cite-owenwu2016] (isothermal transonic Parker wind, Lambert-W form) | The photospheric Mach number is 1 identically with the launch level at the Bondi radius (the analytical sonic-point limit), falls monotonically as the level retreats inward, and has collapsed by more than six decades at their published shutoff R_p/R_B = 0.1. |
| `tests/test_boiloff.py::test_lambda_equals_two_bondi_radii_over_rp` | Fossati et al. (2017)[^cite-fossati2017] (restricted Jeans parameter); Owen & Wu (2016) | The identity Lambda = 2 R_B / R_p holds for every mean molecular mass, which is what makes the Owen & Wu shutoff equal to Lambda = 20 for every composition. The literature band 15 to 35 brackets the default threshold. |
| `tests/test_boiloff.py::test_bondi_cap_is_misener_eq10_at_the_wind_temperature` | Misener et al. (2025)[^cite-misener2025], Eq. (10); Gupta & Schlichting (2020)[^cite-gupta2020], Eq. (10) | The Bondi cap equals the closed form at the wind temperature T_eq / 2^(1/4) on a hydrostatic launch level, differs from the same form at T_eq by 2^(-3/8) exp(Lambda_launch (2^(1/4) - 1)), and stands to the Parker rate as tau_launch (T_wind / T_launch) exp(1/2 - Mach^2 / 2), so it binds on an optically thin launch level and not on a thick one. |
| `tests/test_boiloff.py::test_luminosity_cap_carries_the_tidal_barrier_reduction` | Erkaev et al. (2007)[^cite-erkaev2007], Eq. (17) and Table 1 | The luminosity cap divides by the tidally reduced escape barrier, so it rises by the enhancement factor 1 / K, whose value at xi = 3 is the 1.92 printed in their Table 1. K = 1 reproduces the untidal cap identically. |

## Notes

The Bondi cap follows Misener et al. (2025, Eq. (10)) at their recommended wind temperature, which carries a factor e^2 the Gupta & Schlichting (2020, Eq. (10)) original lacks, and the luminosity cap Gupta & Schlichting (2019)[^cite-gupta2019] Eq. (8), whose barrier carries the Erkaev tidal factor at xi = R_Hill / R_p so that an admitted residual and the energy-limited rate it competes against measure the same barrier. The cap ordering and the gate semantics are asserted as physics invariants in the same file. The Tang et al. (2024)[^cite-tang2024] Eq. (8) timescale comparison runs as a diagnostic only.

## Anchor type

Analytical limit plus published benchmark (the printed shutoff).

Date of last comparison against the sources: 2026-09-28.

[^cite-owenwu2016]: J. E. Owen, Y. Wu, *[Atmospheres of Low-mass Planets: The "Boil-off"](https://doi.org/10.3847/0004-637X/817/2/107)*, The Astrophysical Journal, 817, 107, 2016. [SciX](https://scixplorer.org/abs/2016ApJ...817..107O/abstract).
[^cite-fossati2017]: L. Fossati, N. V. Erkaev, H. Lammer, P. E. Cubillos, P. Odert, I. Juvan, K. G. Kislyakova, M. Lendl, D. Kubyshkina, S. J. Bauer, *[Aeronomical constraints to the minimum mass and maximum radius of hot low-mass planets](https://doi.org/10.1051/0004-6361/201629716)*, Astronomy & Astrophysics, 598, A90, 2017. [SciX](https://scixplorer.org/abs/2017A%26A...598A..90F/abstract).
[^cite-misener2025]: W. Misener, M. Schulik, H. E. Schlichting, J. E. Owen, *[Blowin’ in the Nonisothermal Wind: Core-powered Mass Loss with Hydrodynamic Radiative Transfer](https://doi.org/10.3847/1538-4357/ada777)*, The Astrophysical Journal, 980, 152, 2025. [SciX](https://scixplorer.org/abs/2025ApJ...980..152M/abstract).
[^cite-gupta2020]: A. Gupta, H. E. Schlichting, *[Signatures of the core-powered mass-loss mechanism in the exoplanet population: dependence on stellar properties and observational predictions](https://doi.org/10.1093/mnras/staa315)*, Monthly Notices of the Royal Astronomical Society, 493, 792-806, 2020. [SciX](https://scixplorer.org/abs/2020MNRAS.493..792G/abstract).
[^cite-erkaev2007]: N. V. Erkaev, Yu. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).
[^cite-gupta2019]: A. Gupta, H. E. Schlichting, *[Sculpting the valley in the radius distribution of small exoplanets as a by-product of planet formation: the core-powered mass-loss mechanism](https://doi.org/10.1093/mnras/stz1230)*, Monthly Notices of the Royal Astronomical Society, 487, 24-33, 2019. [SciX](https://scixplorer.org/abs/2019MNRAS.487...24G/abstract).
[^cite-tang2024]: Y. Tang, J. J. Fortney, R. Murray-Clay, *[Assessing Core-powered Mass Loss in the Context of Early Boil-off: Minimal Long-lived Mass Loss for the Sub-Neptune Population](https://doi.org/10.3847/1538-4357/ad8567)*, The Astrophysical Journal, 976, 221, 2024. [SciX](https://scixplorer.org/abs/2024ApJ...976..221T/abstract).
