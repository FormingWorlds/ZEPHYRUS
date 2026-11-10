# Validation: `src/zephyrus/diagnostics.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the reporting quantities of `zephyrus.diagnostics`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_diagnostics.py::test_erkaev_critical_temperature_normalization` | Erkaev et al. (2007)[^cite-erkaev2007], Eq. (23) | The tidally corrected critical exobase temperature recovers the printed Jupiter normalisation $1.45 \times 10^{5}$ K in the wide-orbit, exobase-at-radius limit, vanishes at the Roche lobe, and scales linearly with planet mass. |

## Notes

The Johnson et al. (2013)[^cite-johnson2013] Eq. (10) transonic energy criterion and the Guo (2024)[^cite-guo2024] regime triple are asserted through their published scalings and limits as physics invariants in the same file. The threshold-potential screen of Salz et al. (2016)[^cite-salz2016] carries the two values of their abstract, 13.11 for the validity of energy-limited escape and about 13.6 for hydrodynamically stable thermospheres, read from the paper. The module reports the hydrogen-dominated scope of their grid beside the verdict.

## Anchor type

Published benchmark (printed normalisation and closed form).

Date of last comparison against the sources: 2026-08-20.

[^cite-erkaev2007]: N. V. Erkaev, Yu. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).
[^cite-johnson2013]: R. E. Johnson, A. N. Volkov, J. T. Erwin, *[Molecular-kinetic Simulations of Escape from the Ex-planet and Exoplanets: Criterion for Transonic Flow](https://doi.org/10.1088/2041-8205/768/1/L4)*, The Astrophysical Journal Letters, 768, L4, 2013. [SciX](https://scixplorer.org/abs/2013ApJ...768L...4J/abstract).
[^cite-guo2024]: J. H. Guo, *[Characterization of the regimes of hydrodynamic escape from low-mass exoplanets](https://doi.org/10.1038/s41550-024-02269-w)*, Nature Astronomy, 8, 920-928, 2024. [SciX](https://scixplorer.org/abs/2024NatAs...8..920G/abstract).
[^cite-salz2016]: M. Salz, P. C. Schneider, S. Czesla, J. H. M. M. Schmitt, *[Energy-limited escape revised. The transition from strong planetary winds to stable thermospheres](https://doi.org/10.1051/0004-6361/201527042)*, Astronomy & Astrophysics, 585, L2, 2016. [SciX](https://scixplorer.org/abs/2016A%26A...585L...2S/abstract).
