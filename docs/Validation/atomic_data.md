# Validation: `src/zephyrus/atomic_data.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the transcribed cooling data of `zephyrus.atomic_data`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_atomic_data.py::test_three_level_transcription_spot_values` | Nakayama, Ikoma & Terada (2022), ApJ 937, 72 ([ADS 2022ApJ...937...72N](https://ui.adsabs.harvard.edu/abs/2022ApJ...937...72N)), Appendix C Tables 2 to 5 | Spot Einstein coefficients (the N transauroral and auroral lines, hydrogen Lyman-alpha) and the corrected O+ level set (4S-2D-2P with the printed weights 4, 10, 6) pin the transcription against the printed tables. |
| `tests/test_atomic_data.py::test_badnell_fit_magnitude_slope_and_misprint_guard` | Badnell (2006), ApJS 167, 334 (arXiv astro-ph/0604144), Eqs. 1 and 2 | The radiative recombination fit lands in the published low $10^{-13}$ cm^3 s^-1 decade at $10^{4}$ K with the published falling slope, its six coefficients are pinned against the $Z = 7$, $N = 6$ row of the table, and it is discriminated against the garbled rendering of the same fit printed by Chatterjee & Pierrehumbert (2026, ApJ 998, 236, their Eq. 35), which disagrees by more than a factor 2 and whose $T_2$ differs from the table in the order of two digits. |

## Notes

The CO2 15 micron band and the O fine-structure channels follow Johnstone et al. (2018, A&A 617, A107) Eqs. (34) to (38) and (41) to (43); their coronal (collision-limited) reductions and detailed-balance structure are asserted as physics invariants in the same test file. The band deexcitation coefficients are measured only over roughly 150 to 500 K, a limitation the module documents.

## Anchor type

Published benchmark (printed atomic tables and fit coefficients).

Date of last comparison against the sources: 2026-08-20.
