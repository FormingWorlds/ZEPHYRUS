# Validation: `src/zephyrus/knudsen.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the collision cross-section fallback order of `zephyrus.knudsen` against laboratory data.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_knudsen.py::test_viscosities_reproduce_measurements_within_7_percent` | Laricchiuta et al. (2009)[^cite-laricchiuta2009], EPJD 54, 607 (phenomenological collision integrals, Eqs. (2) to (4) with their appendix Tables 3 to 5); measured 300 K dynamic viscosities of N2, O2, CO, and CO2 (CRC Handbook values) | The first Chapman-Enskog approximation built on the transcribed Omega(2,2)* integrals reproduces all four measured viscosities within 7%, anchoring the whole coefficient transcription on laboratory measurements. |

## Notes

The companion transcription-pin test freezes the momentum-transfer cross sections evaluated from the fit at transcription time, so any later coefficient corruption fails even where no measurement exists. The hydrogen route is pinned to its construction values (sigma(H-H) = $6.4 \times 10^{-20}$ m^2 at $10^{4}$ K from Zahnle et al. 1990, Eq. (30), on the Zahnle & Kasting 1986 Table I diffusion parameter).

**Assumed radii.** The geometric hard-sphere fallback takes van der Waals radii from the package's table. Six elements carry no radius there, so aluminium, phosphorus, chlorine, potassium, calcium, and titanium reach the fallback on an assumed radius, with a provenance class of their own so an assumed number cannot be read as a published one. For two of them the gap is in the source: Bondi (1964) [^cite-bondi1964] gives no radius for aluminium or calcium, among the 16 of 44 main-group elements missing from that table, which Mantina et al. (2009) later supplied on the same scale [^cite-mantina2009]. Phosphorus, chlorine, and potassium are main-group elements Bondi does tabulate, and titanium is a transition metal. The table does not carry them yet.

## Anchor type

Published benchmark (laboratory viscosity measurements).

Date of last comparison against the sources: 2026-08-20.

[^cite-bondi1964]: A. Bondi, *[van der Waals Volumes and Radii](https://doi.org/10.1021/j100785a001)*, The Journal of Physical Chemistry, 68, 441-451, 1964. [SciX](https://scixplorer.org/abs/1964JPhCh..68..441B/abstract).

[^cite-mantina2009]: M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, D. G. Truhlar, *[Consistent van der Waals Radii for the Whole Main Group](https://doi.org/10.1021/jp8111556)*, The Journal of Physical Chemistry A, 113, 5806-5812, 2009. [SciX](https://scixplorer.org/abs/2009JPCA..113.5806M/abstract).
[^cite-laricchiuta2009]: A. Laricchiuta, D. Bruno, M. Capitelli, C. Catalfamo, R. Celiberto, G. Colonna, P. Diomede, D. Giordano, C. Gorse, S. Longo, D. Pagano, F. Pirani, *[High temperature Mars atmosphere. Part I: transport cross sections](https://doi.org/10.1140/epjd/e2009-00192-7)*, The European Physical Journal D, 54, 607-612, 2009. [SciX](https://scixplorer.org/abs/2009EPJD...54..607L/abstract).
