# Validation: `src/zephyrus/hydrostatic.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the hydrostatic branch of `zephyrus.hydrostatic`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_hydrostatic.py::test_volkov_flat_factor_against_published_correction` | Volkov et al. (2011), Phys. Fluids 23, 066601, Eq. 9 and their printed c(lambda) table | The drifting-Maxwellian flux ratio reduces to unity at zero bulk velocity and reproduces the printed linear coefficient across lambda 1 to 106; the test also demonstrates that this correction and the flat kinetic factor the branch applies have opposite slopes in lambda, so applying both would double-count. |
| `tests/test_hydrostatic.py::test_yelle_figure1_mars_hydrogen_flux` | Yelle (2024), Icarus 416, 116099, Figure 1 (fully specified Mars model) | The branch reproduces the diffusion-limited hydrogen plateau at $2.4 \times 10^{8}$ cm^-2 s^-1 (40% tolerance: the binary H-CO2 coefficient source of the original calculation is not pinned in the paper), the saturation above 200 K, and the Jeans-limited collapse at 100 K checked as a regime. |

## Notes

The escape-temperature identities (the lambda = 2 neutral criterion and the plasma value at half of it, after Chatterjee & Pierrehumbert 2026, their Eq. 34) and the stoichiometric mass conservation of the element mapping are asserted as physics invariants in the same file.

## Anchor type

Published benchmark (printed coefficient table and a published model figure) plus the analytical rest limit.

Date of last comparison against the sources: 2026-08-20.
