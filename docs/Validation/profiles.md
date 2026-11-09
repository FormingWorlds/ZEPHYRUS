# Validation: `src/zephyrus/profiles.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the escape working levels of `zephyrus.profiles` against the published literature.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_profiles.py::test_lopez_base_pressure_is_the_nanobar_level` | Lopez (2017)[^cite-lopez2017], MNRAS 472, 245, Section 2 (wind-base prescription); Murray-Clay et al. (2009)[^cite-murrayclay2009], ApJ 693, 23, Section 2.1 (the nanobar tau = 1 level) | Pins the Lopez base pressure `P_base = mu g / sigma_nu0` to the published nanobar scale (within a factor of a few of $10^{-4}$ Pa) for the Murray-Clay fiducial hot Jupiter, with exact linearity in gravity and insensitivity to the proton-mass versus atomic-weight convention below 1%. |

## Notes

Both papers quote the scale of the level rather than a precise value, so the pin is an order-of-magnitude anchor with sign and scale guards. The fixed-point iteration of the base pressure on a profile is exercised separately with clamping semantics in the same test file.

## Anchor type

Published benchmark (order-of-magnitude level).

Date of last comparison against the sources: 2026-08-20.

[^cite-lopez2017]: E. D. Lopez, *[Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor](https://doi.org/10.1093/mnras/stx1558)*, Monthly Notices of the Royal Astronomical Society, 472, 245-253, 2017. [SciX](https://scixplorer.org/abs/2017MNRAS.472..245L/abstract).
[^cite-murrayclay2009]: R. A. Murray-Clay, E. I. Chiang, N. Murray, *[Atmospheric Escape From Hot Jupiters](https://doi.org/10.1088/0004-637X/693/1/23)*, The Astrophysical Journal, 693, 23-42, 2009. [SciX](https://scixplorer.org/abs/2009ApJ...693...23M/abstract).
