# Validation: `src/zephyrus/knudsen.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the collision cross-section ladder of `zephyrus.knudsen` against laboratory data.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_knudsen.py::test_viscosities_reproduce_measurements_within_7_percent` | Laricchiuta et al. (2009), EPJD 54, 607 (phenomenological collision integrals, Eqs. 2 to 4 with their appendix Tables 3 to 5); measured 300 K dynamic viscosities of N2, O2, CO, and CO2 (CRC Handbook values) | The first Chapman-Enskog approximation built on the transcribed Omega(2,2)* integrals reproduces all four measured viscosities within 7 percent, anchoring the whole coefficient transcription on laboratory measurements. |

## Notes

The companion transcription-pin test freezes the momentum-transfer cross sections evaluated from the fit at transcription time, so any later coefficient corruption fails even where no measurement exists. The hydrogen route is pinned to its construction values (sigma(H-H) = 6.4e-20 m^2 at 1e4 K from Zahnle et al. 1990, Eq. 30, on the Zahnle & Kasting 1986 Table I diffusion parameter).

## Anchor type

Published benchmark (laboratory viscosity measurements).

Date of last comparison against the sources: 2026-08-20.
