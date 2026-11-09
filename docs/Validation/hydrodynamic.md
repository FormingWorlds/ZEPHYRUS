# Validation: `src/zephyrus/hydrodynamic.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the hydrodynamic branch of `zephyrus.hydrodynamic`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_hydrodynamic.py::test_erkaev_table1_enhancement_factors` | Erkaev et al. (2007)[^cite-erkaev2007], Table 1 and Eq. (17) | The tidal factor reproduces all seven printed enhancement factors 1/K within 1%, falling monotonically toward 1 with xi. |
| `tests/test_hydrodynamic.py::test_caldiroli_fit_spot_values_and_flux_guard` | Caldiroli et al. (2022)[^cite-caldiroli2022], Appendix A.1 | Six spot evaluations across the validity box reproduce the fitted efficiency to 5% across three decades of collapse, and the complex-valued region below F_XUV/rho_p = $10^{2}$ (cgs) is rejected with a flag rather than evaluated. |
| `tests/test_hydrodynamic.py::test_wind_mean_masses_against_lopez_pairs` | Lopez (2017)[^cite-lopez2017] (printed wind mean-mass pairs) | The generalised ionised-wind rule reproduces the printed steam pair (3, 6) in atomic mass units and the per-ion entry 1.3 of the printed H/He pair. The per-particle entry 0.62 of that pair is not reachable: the rule makes the per-ion mass twice the per-particle mass, so 1.3 halves to 0.65, and the 4.6% gap is an inconsistency between the two printed entries rather than an error in the rule. The test pins 0.65 and asserts 0.62 is excluded, so a future tolerance cannot absorb the difference. |
| `tests/test_hydrodynamic.py::test_murray_clay_fiducial_hot_jupiter_anchors` | Murray-Clay et al. (2009)[^cite-murrayclay2009] (fiducial hot Jupiter) | The chain reproduces their base Jeans parameter 5.49 (5%), places the sonic point at lambda_b/2 planetary radii inside their stated 2 to 4, and lands their printed sonic-point Knudsen numbers ($10^{-4}$ and $10^{-5}$ at 450 and $5 \times 10^{5}$ erg cm^-2 s^-1) within a factor 3 using their Coulomb cross section. |
| `tests/test_hydrodynamic.py::test_photon_limited_rate_and_its_crossover_with_el` | Owen & Alvarez (2016)[^cite-owenalvarez2016], Eq. (10) (analytical limit) | The photon-limited rate is one particle per intercepted ionising photon, linear in the flux and free of the planet mass, and on the energy-limited rate's own disk it crosses the energy-limited rate at the efficiency G M_p K m_ion / (e_ion R_p), on asymmetric inputs that move the crossover under a swapped radius or a missing tidal factor. |

## Notes

Murray-Clay et al. fit their numerical models with flux exponents 0.6 (radiation-recombination limited) and 0.9 (energy limited) where the analytic chains carry 0.5 and 1.0. Rate-level comparisons against their models must budget that difference, which the diagnostics module documents as reporting constants.

## Anchor type

Published benchmarks (printed tables, fits, and worked values) and one analytical limit (the photon-limited crossover).

Date of last comparison against the sources: 2026-09-28.

[^cite-erkaev2007]: N. V. Erkaev, Yu. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).
[^cite-caldiroli2022]: A. Caldiroli, F. Haardt, E. Gallo, R. Spinelli, I. Malsky, E. Rauscher, *[Irradiation-driven escape of primordial planetary atmospheres. II. Evaporation efficiency of sub-Neptunes through hot Jupiters](https://doi.org/10.1051/0004-6361/202142763)*, Astronomy & Astrophysics, 663, A122, 2022. [SciX](https://scixplorer.org/abs/2022A%26A...663A.122C/abstract).
[^cite-lopez2017]: E. D. Lopez, *[Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor](https://doi.org/10.1093/mnras/stx1558)*, Monthly Notices of the Royal Astronomical Society, 472, 245-253, 2017. [SciX](https://scixplorer.org/abs/2017MNRAS.472..245L/abstract).
[^cite-murrayclay2009]: R. A. Murray-Clay, E. I. Chiang, N. Murray, *[Atmospheric Escape From Hot Jupiters](https://doi.org/10.1088/0004-637X/693/1/23)*, The Astrophysical Journal, 693, 23-42, 2009. [SciX](https://scixplorer.org/abs/2009ApJ...693...23M/abstract).
[^cite-owenalvarez2016]: J. E. Owen, M. A. Alvarez, *[UV Driven Evaporation of Close-in Planets: Energy-limited, Recombination-limited, and Photon-limited Flows](https://doi.org/10.3847/0004-637X/816/1/34)*, The Astrophysical Journal, 816, 34, 2016. [SciX](https://scixplorer.org/abs/2016ApJ...816...34O/abstract).
