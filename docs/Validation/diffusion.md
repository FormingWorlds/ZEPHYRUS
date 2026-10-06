# Validation: `src/zephyrus/diffusion.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the binary-diffusion library of `zephyrus.diffusion` against its independent published compilations.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_diffusion.py::test_sn88_unit_reading_matches_zk86_table` | Sasaki & Nakazawa (1988), EPSL 89, 323, Table 1; Zahnle & Kasting (1986), Icarus 68, 462, Table I | The unit reading of the Sasaki & Nakazawa table (which prints no units) reproduces five in-H2 entries of the independent 1986 compilation to 3%; no other unit choice comes within two orders of magnitude. |
| `tests/test_diffusion.py::test_zk23_zk86_cross_compilation_agreement` | Zahnle & Kasting (2023), GeCoA 361, 228, Table 2; Zahnle & Kasting (1986), Icarus 68, 462, Table I | The shared rows of the two compilations agree to 4% where measured (Marrero & Mason 1972 lineage) and within 35% where the 2023 rows are estimates, with the measured rows agreeing more tightly than the estimated ones. |

## Notes

The Eq. (10) scaling rule is validated in sample on the three atomic rows the 2023 authors themselves obtained by scaling, and out of sample on the Kr and Xe rows of the 1986 table (which are not sources of this library), landing inside the 30% scaled provenance class in natural log.

## Anchor type

Cross-compilation cross-check (two independent compilations of the Marrero & Mason 1972 measurements, 35 years apart).

Date of last comparison against the sources: 2026-08-20.
