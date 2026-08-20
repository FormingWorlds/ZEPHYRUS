# Validation: `src/zephyrus/diagnostics.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the reporting quantities of `zephyrus.diagnostics`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_diagnostics.py::test_erkaev_critical_temperature_normalization` | Erkaev et al. (2007), A&A 472, 329, Eq. 23 | The tidally corrected critical exobase temperature recovers the printed Jupiter normalization 1.45e5 K in the wide-orbit, exobase-at-radius limit, vanishes at the Roche lobe, and scales linearly with planet mass. |

## Notes

The Johnson et al. (2013, ApJL 768, L4, Eq. 10) transonic energy criterion and the Guo (2024, arXiv:2405.13283) regime triple are asserted through their published scalings and limits as physics invariants in the same file. The threshold-potential screen attributed to Salz et al. (2016) is quoted from secondary literature, which the module states beside the values.

## Anchor type

Published benchmark (printed normalization and closed form).

Date of last comparison against the sources: 2026-08-20.
