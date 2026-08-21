# Validation: `src/zephyrus/hydrodynamic.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the hydrodynamic branch of `zephyrus.hydrodynamic`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_hydrodynamic.py::test_erkaev_table1_enhancement_factors` | Erkaev et al. (2007), A&A 472, 329, Table 1 and Eq. 17 | The tidal factor reproduces all seven printed enhancement factors 1/K within 1%, falling monotonically toward 1 with xi. |
| `tests/test_hydrodynamic.py::test_caldiroli_fit_spot_values_and_flux_guard` | Caldiroli et al. (2022), A&A 663, A122, Appendix A.1 | Six spot evaluations across the validity box reproduce the fitted efficiency to 5% across three decades of collapse, and the complex-valued region below F_XUV/rho_p = 1e2 (cgs) is rejected with a flag rather than evaluated. |
| `tests/test_hydrodynamic.py::test_wind_mean_masses_reproduce_lopez_pairs` | Lopez (2017), MNRAS 472, 245 (printed wind mean-mass pairs) | The generalized ionized-wind rule reproduces the printed H/He pair (0.62, 1.3) and steam pair (3, 6) in proton masses. |
| `tests/test_hydrodynamic.py::test_murray_clay_fiducial_hot_jupiter_anchors` | Murray-Clay et al. (2009), ApJ 693, 23 (fiducial hot Jupiter) | The chain reproduces their base Jeans parameter 5.49 (5%), places the sonic point at lambda_b/2 planetary radii inside their stated 2 to 4, and lands their printed sonic-point Knudsen numbers (1e-4 and 1e-5 at 450 and 5e5 erg cm^-2 s^-1) within a factor 3 using their Coulomb cross section. |

## Notes

Murray-Clay et al. fit their numerical models with flux exponents 0.6 (radiation-recombination limited) and 0.9 (energy limited) where the analytic chains carry 0.5 and 1.0; rate-level comparisons against their models must budget that difference, which the diagnostics module documents as reporting constants.

## Anchor type

Published benchmarks (printed tables, fits, and worked values).

Date of last comparison against the sources: 2026-08-20.
