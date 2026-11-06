# Validation: `src/zephyrus/knudsen.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the collision cross-section fallback order of `zephyrus.knudsen` against laboratory data.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_knudsen.py::test_viscosities_reproduce_measurements_within_7_percent` | Laricchiuta et al. (2009), EPJD 54, 607 (phenomenological collision integrals, Eqs. (2) to (4) with their appendix Tables 3 to 5); measured 300 K dynamic viscosities of N2, O2, CO, and CO2 (CRC Handbook values) | The first Chapman-Enskog approximation built on the transcribed Omega(2,2)* integrals reproduces all four measured viscosities within 7%, anchoring the whole coefficient transcription on laboratory measurements. |

## Notes

The companion transcription-pin test freezes the momentum-transfer cross sections evaluated from the fit at transcription time, so any later coefficient corruption fails even where no measurement exists. The hydrogen route is pinned to its construction values (sigma(H-H) = $6.4 \times 10^{-20}$ m^2 at $10^{4}$ K from Zahnle et al. 1990, Eq. (30), on the Zahnle & Kasting 1986 Table I diffusion parameter).

**Assumed radii.** The geometric hard-sphere fallback takes van der Waals radii from the package's table. Six elements carry no radius there, so aluminium, phosphorus, chlorine, potassium, calcium, and titanium reach the fallback on an assumed radius, with a provenance class of their own so an assumed number cannot be read as a published one. For two of them the gap is in the source: Bondi (1964) [^bondi] gives no radius for aluminium or calcium, among the 16 of 44 main-group elements missing from that table, which Mantina et al. (2009) later supplied on the same scale [^mantina]. Phosphorus, chlorine, and potassium are main-group elements Bondi does tabulate, and titanium is a transition metal. The table does not carry them yet.

## Anchor type

Published benchmark (laboratory viscosity measurements).

Date of last comparison against the sources: 2026-08-20.

[^bondi]: Bondi, A. (1964). van der Waals Volumes and Radii. *The Journal of Physical Chemistry, 68*(3), 441–451. https://doi.org/10.1021/j100785a001

[^mantina]: Mantina, M., Chamberlin, A. C., Valero, R., Cramer, C. J., & Truhlar, D. G. (2009). Consistent van der Waals Radii for the Whole Main Group. *The Journal of Physical Chemistry A, 113*(19), 5806–5812. https://doi.org/10.1021/jp8111556
