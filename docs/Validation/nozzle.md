# Validation: `src/zephyrus/nozzle.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the Roche-lobe overflow nozzle rate of `zephyrus.nozzle` against the published model.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_nozzle.py::test_table2_planets_reproduce_published_rates` | Jackson et al. (2017), ApJ 835, 145, Table 2 (Kepler-21 b and CoRoT-24 b) with their Section 3 input prescriptions | The full closed form (their Eqs. 3, 10, 13, and 14 plus the Eggleton lobe radius) reproduces the two printed rates whose donors sit well inside their lobes, where the potential approximation is good, to 3%. |
| `tests/test_nozzle.py::test_lobe_filling_binaries_land_on_figure5` | Jackson et al. (2017) Figure 5 and Table 1 (Ritter 1988 Table A1 binaries, 0.8 solar-mass accretor) | The saturated lobe-filling limit lands on the printed solid curve at two donor masses within figure-reading tolerance, anchoring the boundary value used at and beyond lobe contact. |
| `tests/test_nozzle.py::test_curvature_equal_mass_value` | Jackson et al. (2017) Section 2.1, the printed equal-mass curvature A(1) = 8 and the small-mass-ratio expansion of the L1 position | The Eq. (10) fit returns the printed equal-mass value without fit error and carries the second-order small-q coefficient $2 \cdot 3^{2/3}$. |
| `tests/test_nozzle.py::test_potential_at_lobe_matches_numerical_l1` | The exact corotating Roche potential, L1 solved numerically on the star-planet axis | The Eq. (14) volume-averaged potential evaluated at the Eggleton lobe radius matches the exact L1 potential to better than $10^{-4}$ relative at a planetary mass ratio, and the escape barrier built from it to 0.5%. |

## Notes

The two near-lobe hot Jupiters of their Table 2 are deliberately not pinned: reproducing them requires the photospheric-radius distortion conversion of their Appendix, which this implementation omits as a stated convention, and the residual factor of 1.4 to 1.7 sits inside the paper's own quoted factor-of-two error from the approximate potentials. The rate is exponentially sensitive to the barrier there ($d\ln\dot{M}/d\ln r_\mathrm{ph} \approx G M_\mathrm{p}/(r_\mathrm{ph} v_\mathrm{th}^2)$, about 50 for a hot Jupiter), which is also why the two pinned planets carry a 3% tolerance: an exponent near 16 amplifies physical-constant conventions tenfold, while the transcription errors the pin exists to catch move the rate by factors of several to decades.

## Anchor type

Published benchmark (model rates and printed coefficients) plus an exact-geometry cross-check.

Date of last comparison against the source: 2026-08-31.
