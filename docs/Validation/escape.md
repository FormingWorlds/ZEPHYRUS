# Validation: `src/zephyrus/escape.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the behaviour of `zephyrus.escape` against the published energy-limited mass-loss formulation and its closed-form analytical limit.

| Test id | Reference | Source page | Scope |
|---|---|---|---|
| `tests/test_escape.py::test_el_escape_scaling2_matches_lopez2012_closed_form` | Lopez, Fortney & Miller (2012), ApJ 761:59, Eqs. 2-4 (analytical closed form) | [ADS 2012ApJ...761...59L](https://ui.adsabs.harvard.edu/abs/2012ApJ...761...59L) | Pins the EL mass-loss rate for the default `scaling=2` radius term (`R^3 = Rp * Rxuv^2`) at the Earth-like reference geometry with no tidal correction. |
| `tests/test_escape.py::test_el_escape_scaling3_matches_lehmer_catling_closed_form` | Lehmer & Catling (2017), ApJ 845:130, Eq. 1 (analytical closed form) | [ADS 2017ApJ...845..130L](https://ui.adsabs.harvard.edu/abs/2017ApJ...845..130L) | Pins the EL mass-loss rate for the alternative `scaling=3` radius term (`R^3 = Rxuv^3`) at the same reference geometry. |

## Re-derivation note

`escape.EL_escape` returns the energy-limited mass-loss rate

```
escape_EL = epsilon * pi * R^3 * Fxuv / (G * Mp * K_tide)
```

where the radius term `R^3` is selected by the `scaling` argument: `scaling=2` (default) uses `Rp * Rxuv^2` following Lopez, Fortney & Miller (2012); `scaling=3` uses `Rxuv^3` following Lehmer & Catling (2017). `K_tide` is unity in the no-tidal branch, so both pinned tests reduce the formula to the pure radius-scaling closed form.

Both tests evaluate the rate at a shared Earth-like reference geometry with `Rxuv` chosen distinct from `Rp` (`Rxuv = 1.2 * Rp`), so the two scaling branches differ by 20%. This is the point of the geometry: a degenerate `Rp = Rxuv` choice would make `Rp * Rxuv^2` and `Rxuv^3` identical and a regression that swapped the default scaling would pass silently. With the branches 20% apart each test carries a wrong-scaling discrimination guard that re-evaluates the sibling branch and asserts the difference exceeds 1% of the pinned value.

Scale: the pinned rates are order `1e6` to `1e7` kg s-1. Each test brackets the magnitude with `1e6 < val < 1e8`. A gross unit slip, such as radii in centimetres rather than metres, lands near `1e12` and fails this band; a dropped `pi` or `epsilon` factor stays within the band and is instead caught by the `rel=1e-9` pin. A sign flip is caught separately by an explicit `val > 0` assertion, since the closed-form rate is always non-negative for physically valid inputs.

## Anchor type

Analytical limit (closed-form energy-limited rate at a fixed geometry), one per radius-scaling branch. The published formulations of Lopez, Fortney & Miller (2012) and Lehmer & Catling (2017) supply the exact algebraic form; the tests pin the hand-evaluated rate at the reference geometry to float rounding (`rel=1e-9`) because the source uses the same closed form and constants.

## Cross-references

- `src/zephyrus/escape.py`, `EL_escape` docstring References section: cites Lopez, Fortney & Miller (2012), Eqs. 2-4 and Lehmer & Catling (2017), Eq. 1 for the two radius-scaling variants, and Erkaev et al. (2007), Eq. 21 for the tidal correction.
- `docs/Explanations/model.md`: user-facing overview of the energy-limited escape model and its tidal correction.

## Last comparison

2026-07-10, against `src/zephyrus/escape.py` at the branch head.
