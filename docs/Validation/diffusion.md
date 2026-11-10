# Validation: `src/zephyrus/diffusion.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the binary-diffusion library of `zephyrus.diffusion` against its independent published compilations.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_diffusion.py::test_sn88_unit_reading_matches_zk86_table` | Sasaki & Nakazawa (1988)[^cite-sasaki1988], Table 1; Zahnle & Kasting (1986)[^cite-zahnle1986], Table I | The unit reading of the Sasaki & Nakazawa table (which prints no units) reproduces five in-H2 entries of the independent 1986 compilation to 3%. No other unit choice comes within two orders of magnitude. |
| `tests/test_diffusion.py::test_zk23_zk86_cross_compilation_agreement` | Zahnle & Kasting (2023)[^cite-zahnle2023], Table 2; Zahnle & Kasting (1986), Table I | The shared rows of the two compilations agree to 4% where measured (Marrero & Mason 1972[^cite-marrero1972] lineage) and within 35% where the 2023 rows are estimates, with the measured rows agreeing more tightly than the estimated ones. |

## Notes

The Eq. (10) scaling rule is validated in sample on the three atomic rows the 2023 authors themselves obtained by scaling, and out of sample on the Kr and Xe rows of the 1986 table (which are not sources of this library), landing inside the 30% scaled provenance class in natural log.

## Anchor type

Cross-compilation cross-check (two independent compilations of the Marrero & Mason 1972 measurements, 35 years apart).

Date of last comparison against the sources: 2026-08-20.

[^cite-sasaki1988]: S. Sasaki, K. Nakazawa, *[Origin of isotopic fractionation of terrestrial Xe: hydrodynamic fractionation during escape of the primordial H₂-He atmosphere](https://doi.org/10.1016/0012-821x(88)90120-3)*, Earth and Planetary Science Letters, 89, 323-334, 1988. [SciX](https://scixplorer.org/abs/1988E%26PSL..89..323S/abstract).
[^cite-zahnle1986]: K. J. Zahnle, J. F. Kasting, *[Mass fractionation during transonic escape and implications for loss of water from Mars and Venus](https://doi.org/10.1016/0019-1035(86)90051-5)*, Icarus, 68, 462-480, 1986. [SciX](https://scixplorer.org/abs/1986Icar...68..462Z/abstract).
[^cite-zahnle2023]: K. J. Zahnle, J. F. Kasting, *[Elemental and isotopic fractionation as fossils of water escape from Venus](https://doi.org/10.1016/j.gca.2023.09.023)*, Geochimica et Cosmochimica Acta, 361, 228-244, 2023. [SciX](https://scixplorer.org/abs/2023GeCoA.361..228Z/abstract).
[^cite-marrero1972]: T. R. Marrero, E. A. Mason, *[Gaseous Diffusion Coefficients](https://doi.org/10.1063/1.3253094)*, Journal of Physical and Chemical Reference Data, 1, 3-118, 1972. [SciX](https://scixplorer.org/abs/1972JPCRD...1....3M/abstract).
