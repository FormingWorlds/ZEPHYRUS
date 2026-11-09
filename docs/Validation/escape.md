# Validation: `zephyrus.escape.EL_escape`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the behaviour of `zephyrus.escape` against the published energy-limited mass-loss formulation and its closed-form analytical limit.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_escape.py::test_el_escape_scaling2_matches_erkaev2007_closed_form` | Erkaev et al. (2007)[^cite-erkaev2007], A&A 472:329, Eq. (21) (energy-limited mass-loss rate, analytical closed form) | [ADS 2007A&A...472..329E](https://ui.adsabs.harvard.edu/abs/2007A%26A...472..329E) | Pins the EL mass-loss rate for the default `scaling=2` radius term (`R^3 = Rp * Rxuv^2`) at the Earth-like reference geometry with no tidal correction. |
| `tests/test_escape.py::test_el_escape_scaling3_matches_lehmer_catling_closed_form` | Lehmer & Catling (2017)[^cite-lehmer2017], ApJ 845:130, Eq. (1) (analytical closed form) | [ADS 2017ApJ...845..130L](https://ui.adsabs.harvard.edu/abs/2017ApJ...845..130L) | Pins the EL mass-loss rate for the alternative `scaling=3` radius term (`R^3 = Rxuv^3`) at the same reference geometry. |

## Re-derivation note

`escape.EL_escape` returns the energy-limited mass-loss rate

```
escape_EL = epsilon * pi * R^3 * Fxuv / (G * Mp * K_tide)
```

where the radius term `R^3` is selected by the `scaling` argument: `scaling=2` (default) uses `Rp * Rxuv^2` following Erkaev et al. (2007); `scaling=3` uses `Rxuv^3` following Lehmer & Catling (2017). `K_tide` is unity in the no-tidal branch, so both pinned tests reduce the formula to the pure radius-scaling closed form.

Both tests evaluate the rate at a shared Earth-like reference geometry with `Rxuv` chosen distinct from `Rp` (`Rxuv = 1.2 * Rp`), so the two scaling branches differ by 20%. This is the point of the geometry: a degenerate `Rp = Rxuv` choice would make `Rp * Rxuv^2` and `Rxuv^3` identical and a regression that swapped the default scaling would pass silently. With the branches 20% apart each test carries a wrong-scaling discrimination guard that re-evaluates the sibling branch and asserts the difference exceeds 1% of the pinned value.

Scale: the pinned rates are order `1e6` to `1e7` kg s-1. Each test brackets the magnitude with `1e6 < val < 1e8`. A gross unit slip, such as radii in centimetres rather than metres, lands near `1e12` and fails this band; a dropped `pi` or `epsilon` factor stays within the band and is instead caught by the `rel=1e-9` pin. A sign flip is caught separately by an explicit `val > 0` assertion, since the closed-form rate is always non-negative for physically valid inputs.

## Anchor type

Analytical limit (closed-form energy-limited rate at a fixed geometry), one per radius-scaling branch. The published formulations of Erkaev et al. (2007) and Lehmer & Catling (2017) supply the exact algebraic form; the tests pin the hand-evaluated rate at the reference geometry to float rounding (`rel=1e-9`) because the source uses the same closed form and constants.

## Cross-references

- `src/zephyrus/hydrodynamic.py`, `EL_escape` docstring References section (`zephyrus.escape` re-exports it, so the released import path is unchanged). It cites Watson et al. (1981)[^cite-watson1981], Lammer et al. (2003)[^cite-lammer2003], Eq. (6), and Erkaev et al. (2007), Eq. (21), for the default `scaling=2` radius term. It cites Lopez et al. (2012)[^cite-lopez2012], Eq. (2), Lopez & Fortney (2013)[^cite-lopez2013], Eq. (1), and Lehmer & Catling (2017), Eq. (1), for the `scaling=3` term, and Erkaev et al. (2007), Eq. (17), for the tidal reduction factor.
- `docs/Explanations/energy_limited.md`: user-facing account of the energy-limited rate, its two radius scalings, and the tidal correction, with the equations numbered. `docs/Explanations/model.md` is the hub that links to it.

## Last comparison

2026-08-28, against `src/zephyrus/hydrodynamic.py` at the branch head.

[^cite-erkaev2007]: N. V. Erkaev, Y. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).
[^cite-lehmer2017]: O. R. Lehmer, D. C. Catling, *[Rocky Worlds Limited to ∼1.8 Earth Radii by Atmospheric Escape during a Star’s Extreme UV Saturation](https://doi.org/10.3847/1538-4357/aa8137)*, The Astrophysical Journal, 845, 130, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...845..130L/abstract).
[^cite-watson1981]: A. J. Watson, T. M. Donahue, J. C. G. Walker, *[The dynamics of a rapidly escaping atmosphere: Applications to the evolution of Earth and Venus](https://doi.org/10.1016/0019-1035(81)90101-9)*, Icarus, 48, 150-166, 1981. [SciX](https://scixplorer.org/abs/1981Icar...48..150W/abstract).
[^cite-lammer2003]: H. Lammer, F. Selsis, I. Ribas, E. F. Guinan, S. J. Bauer, W. W. Weiss, *[Atmospheric Loss of Exoplanets Resulting from Stellar X-Ray and Extreme-Ultraviolet Heating](https://doi.org/10.1086/380815)*, The Astrophysical Journal Letters, 598, L121-L124, 2003. [SciX](https://scixplorer.org/abs/2003ApJ...598L.121L/abstract).
[^cite-lopez2012]: E. D. Lopez, J. J. Fortney, N. Miller, *[How Thermal Evolution and Mass-loss Sculpt Populations of Super-Earths and Sub-Neptunes: Application to the Kepler-11 System and Beyond](https://doi.org/10.1088/0004-637X/761/1/59)*, The Astrophysical Journal, 761, 59, 2012. [SciX](https://scixplorer.org/abs/2012ApJ...761...59L/abstract).
[^cite-lopez2013]: E. D. Lopez, J. J. Fortney, *[The Role of Core Mass in Controlling Evaporation: The Kepler Radius Distribution and the Kepler-36 Density Dichotomy](https://doi.org/10.1088/0004-637X/776/1/2)*, The Astrophysical Journal, 776, 2, 2013. [SciX](https://scixplorer.org/abs/2013ApJ...776....2L/abstract).
