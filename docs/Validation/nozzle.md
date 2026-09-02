# Validation: `src/zephyrus/nozzle.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the Roche-lobe overflow nozzle rate of `zephyrus.nozzle` against the published model.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_nozzle.py::test_table2_deep_launch_planets_reproduce_published_rates` | Jackson et al. (2017), ApJ 835, 145, Table 2 (Kepler-21 b, CoRoT-24 b, WASP-47 e, 55 Cnc e) with their Section 3 input prescriptions | The full closed form (their Eqs. 3, 10, 13, and 14 plus the Eggleton lobe radius) reproduces the four printed rates whose launch level sits below 0.4 Roche lobe radii, to 0.4, 2.1, 3.5, and 5.0 percent, asserted at 6 percent. |
| `tests/test_nozzle.py::test_lobe_filling_binaries_land_on_figure5` | Jackson et al. (2017) Figure 5 and Table 1 (Ritter 1988 Table A1 binaries, 0.8 solar-mass accretor) | The saturated lobe-filling limit lands 4.1 and 8.0 percent below values read off the printed solid curve at two donor masses, asserted at 30 percent to cover the figure read. This anchors the prefactor of the boundary value used at and beyond lobe contact; the barrier is switched off by construction there and is not exercised. |
| `tests/test_nozzle.py::test_curvature_equal_mass_value` | Jackson et al. (2017) Eq. (10) with $b_1 = 2 \cdot 3^{2/3}$, and the printed equal-mass value A(1) = 8 | The fit returns the printed equal-mass value without fit error and carries the second-order small-$q$ coefficient $2 \cdot 3^{2/3}$, which the exact root of their Eq. (7) confirms. |
| `tests/test_nozzle.py::test_potential_at_lobe_matches_numerical_l1` | The exact corotating Roche potential, L1 solved numerically on the star-planet axis | The escape barrier built from the Eq. (14) volume-averaged potential at the Eggleton lobe radius agrees with the exact one to 0.17 percent at a planetary mass ratio and a deep launch level, asserted at 0.5 percent. |

## Notes

**What the Table 2 pin covers, and what it does not.** Their Table 2 lists twenty-one objects. Four are pinned, and they are every row whose launch level sits below 0.4 Roche lobe radii. The other seventeen are not reproduced by this implementation, running from a factor 1.04 to 3.9e6 below the published rates, and the shortfall is monotone in launch depth rather than scattered. Two effects account for it, neither of which is a transcription question:

- The photospheric radius convention. Their Appendix converts a measured transit radius to the volume-equivalent radius of the distorted equipotential, and this implementation omits that conversion because its input is a one-dimensional profile radius rather than an observed one. The rate is exponential in the launch radius, with $\mathrm{d}\ln\dot{M}/\mathrm{d}\ln r_\mathrm{ph} = G M_\mathrm{p}/(r_\mathrm{ph} v_\mathrm{th}^2)$ running from 9.7 to 276 across their table, so the shift needed to close each gap is only 2.2 to 5.2 percent in radius. On the four confirmed near-lobe hot Jupiters that requirement matches the 2.6 to 3.7 percent the conversion produces, and their published post-conversion residuals of 1.43 to 1.69 are consistent with it.
- Assumed masses. Their six Kepler candidate rows carry $M_\mathrm{p} = 1\,M_\mathrm{Jup}$ by assumption rather than measurement, and need 3.8 to 5.2 percent in radius, which no radius conversion supplies. PTFO8-8695 is a third case: the paper states they assumed it just fills its lobe, so their density sits at the lobe radius rather than at the planet radius.

Replacing the Eq. (14) potential with the exact numerically solved L1 potential moves five of the twenty-one rows inside a factor of two to seven, so the potential approximation is a minor part of the shortfall and the radius convention is the dominant one.

**Where the potential approximation holds.** The 0.17 percent barrier agreement above is measured at a launch level of 0.21 lobe radii. Holding that depth, the agreement stays under 0.5 percent across every mass ratio from $10^{-7}$ to $10^{-2}$. Holding the mass ratio and varying the depth instead, it passes 0.5 percent near 0.42 lobe radii and reaches 16 percent at 0.95, which is why the Table 2 pin stops where it does. A quantity worth not quoting: the same comparison expressed against the total potential agrees to $4 \times 10^{-6}$, but the star's constant term is 99.94 percent of that total and cancels in the barrier, so a tolerance on it constrains the part Eq. (14) actually computes only to about 16 percent of itself.

**A correction to the primary, recorded where it is used.** Their Section 2.1 prints the small-mass-ratio asymptotic $A \approx 4 + 3 (M_\mathrm{d}/3M_\mathrm{a})^{1/3}$, whose coefficient is $3^{2/3}$. The Eq. (10) fit uses $b_1 = 2 \cdot 3^{2/3}$, twice that. Solving their Eq. (7) at the exact numerical L1 root gives $(A-4)/q^{1/3} \to 4.1602 = 2 \cdot 3^{2/3}$ as $q \to 0$, so the fit is right and the printed in-text asymptotic understates the coefficient by half. The implementation follows the fit.

**A caveat on the Figure 5 anchor.** Its two states are constructed by inverting the Eggleton fit at the printed launch radius, which places them on the saturation boundary by construction: perturbing the separation by 0.1 percent drops the rate by factors of 71 and 2786. The 30 percent window is therefore slack on the prefactor (density, sound speed, nozzle area) and tests the barrier not at all. One of the two rows returns `saturated = False`, because the round-trip lands the launch level one unit in the last place inside the lobe; the rate agrees with the boundary value to three parts in $10^{12}$ either way.

**A radius convention the tolerance does not cover.** The pinned rows use the equatorial Jupiter radius, 7.1492e7 m. Switching to the volumetric mean, 6.9911e7 m, moves CoRoT-24 b from 0.979 to 0.631 of the published rate. The convention is stated here rather than absorbed into the tolerance.

## Anchor type

Published benchmark (model rates and printed coefficients) plus an exact-geometry cross-check.

Date of last comparison against the source: 2026-09-01.
