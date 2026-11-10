# Validation: `src/zephyrus/atomic_data.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the transcribed cooling data of `zephyrus.atomic_data`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_atomic_data.py::test_three_level_transcription_spot_values` | Nakayama et al. (2022)[^cite-nakayama2022], Appendix C Tables 2 to 5 | Spot Einstein coefficients (the N transauroral and auroral lines, hydrogen Lyman-alpha) and the corrected O+ level set (4S-2D-2P with the printed weights 4, 10, 6) pin the transcription against the printed tables. |
| `tests/test_atomic_data.py::test_badnell_fit_magnitude_slope_and_misprint_guard` | Badnell (2006)[^cite-badnell2006], Eqs. (1) and (2) | The radiative recombination fit lands in the published low $10^{-13}$ cm^3 s^-1 decade at $10^{4}$ K with the published falling slope, its six coefficients are pinned against the $Z = 7$, $N = 6$ row of the table, and it is discriminated against the garbled rendering of the same fit printed by Chatterjee & Pierrehumbert (2026)[^cite-chatterjee2026], their Eq. (35), which disagrees by more than a factor 2 and whose $T_2$ differs from the table in the order of two digits. |

## Notes

The CO2 15 micron band and the O fine-structure channels follow Johnstone et al. (2018)[^cite-johnstone2018] Eqs. (34) to (38) and (41) to (43). Their coronal (collision-limited) reductions and detailed-balance structure are asserted as physics invariants in the same test file. Two departures from the printed source are pinned rather than inherited. The detailed-balance exponential uses the band quantum in kelvin, $h\nu / k_\mathrm{B} = 959.7$ K, where the source prints 667: that 667 is the bending-mode wavenumber in cm$^{-1}$, and $hc$ times it is the same $1.325 \times 10^{-13}$ erg quantum the source prints two equations earlier, so the printed formula has dropped the second radiation constant of 1.4388 K cm. The excitation rate differs by a factor 2.7 at 300 K and 3.5% at the wind temperatures the thermostat selects. And the escape probability is capped at the non-LTE ceiling of 0.5, which the fitted thin-column branch would otherwise exceed below a column parameter of $3.4 \times 10^{-4}$ and diverge through as the column vanishes. The band deexcitation coefficients are measured only over roughly 150 to 500 K, and the band is evaluated on molecular abundances because it is a base-region coolant. Both limitations the module documents.

## Anchor type

Published benchmark (printed atomic tables and fit coefficients).

Date of last comparison against the sources: 2026-08-20.

[^cite-nakayama2022]: A. Nakayama, M. Ikoma, N. Terada, *[Survival of Terrestrial N₂-O₂ Atmospheres in Violent XUV Environments through Efficient Atomic Line Radiative Cooling](https://doi.org/10.3847/1538-4357/ac86ca)*, The Astrophysical Journal, 937, 72, 2022. [SciX](https://scixplorer.org/abs/2022ApJ...937...72N/abstract).
[^cite-badnell2006]: N. R. Badnell, *[Radiative Recombination Data for Modeling Dynamic Finite-Density Plasmas](https://doi.org/10.1086/508465)*, The Astrophysical Journal Supplement Series, 167, 334-342, 2006. [SciX](https://scixplorer.org/abs/2006ApJS..167..334B/abstract).
[^cite-chatterjee2026]: R. D. Chatterjee, R. T. Pierrehumbert, *[Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline](https://doi.org/10.3847/1538-4357/ae2ffa)*, The Astrophysical Journal, 998, 236, 2026. [SciX](https://scixplorer.org/abs/2026ApJ...998..236C/abstract).
[^cite-johnstone2018]: C. P. Johnstone, M. Güdel, H. Lammer, K. G. Kislyakova, *[Upper atmospheres of terrestrial planets: Carbon dioxide cooling and the Earth’s thermospheric evolution](https://doi.org/10.1051/0004-6361/201832776)*, Astronomy & Astrophysics, 617, A107, 2018. [SciX](https://scixplorer.org/abs/2018A%26A...617A.107J/abstract).
