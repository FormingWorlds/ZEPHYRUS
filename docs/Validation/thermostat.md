# Validation: `src/zephyrus/thermostat.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the statistical-equilibrium cooling of `zephyrus.thermostat`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_thermostat.py::test_hydrogen_system_brackets_the_black_lyalpha_rate` | Black (1981) Lyman-alpha cooling rate as printed by Murray-Clay et al. (2009), ApJ 693, 23, their Eq. 6 | The hydrogen three-level system agrees with the Black rate at order unity, with the ratio declining from about 0.5 at 1e4 K to about 0.3 at 2e4 K because the effective collision strengths are frozen at 1e4 K and the three-level system carries no cascades; the bracket catches transcription errors in the constant and in the exp(-118348 K / T) activation. |

## Notes

Identity with the Black fit is not expected and not asserted: the two treatments differ in their collision-strength temperature dependence and cascade content, and the documented drift is part of the anchor. The detailed-balance (LTE) limit and the coronal linearity of the level populations are asserted as physics invariants in the same file.

## Anchor type

Published benchmark (order-unity consistency with a printed rate).

Date of last comparison against the sources: 2026-08-20.
